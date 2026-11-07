///|
/// AMBER NetCDF (classic CDF-1) trajectory reader/writer.
///
/// Classic NetCDF binary format, all big-endian: magic "CDF\001", int32
/// numrecs, a dimension list (tag 10), a global attribute list (tag 12 or
/// empty tag 0) and a variable list (tag 11). Names are pascal strings
/// (int32 length + bytes + zero padding to a 4-byte boundary). CDF-1 uses
/// int32 `begin` offsets, so files are limited to 2 GiB.
///
/// Record variables (first dimension = unlimited `frame`) are interleaved:
/// record `r` of a variable lives at `begin + r * recsize`, where
/// `recsize` is the sum of the (4-byte aligned) `vsize` of all record
/// variables.
///
/// AMBER conventions written/read here:
/// - dimensions: frame (unlimited), atom, spatial=3, and when a box is
///   present cell_spatial=3, cell_angular=3, angular=3, label=5
/// - variables: coordinates (float [frame,atom,spatial]), time (double
///   [frame]), spatial (char [spatial] = "xyz"), and when a box is present
///   cell_lengths (float [frame,cell_spatial]), cell_angles (float
///   [frame,cell_angular]), cell_spatial (char [cell_spatial] = "abc"),
///   cell_angular (char [cell_angular,label] = "alphabeta gamma")

///|
/// Error type for NetCDF operations.
pub suberror NcdfError {
  NcdfError(String)
}

///|
/// A single frame (snapshot) from an AMBER NetCDF trajectory.
/// `coords` is a flat 3*natoms array; `box_lengths` holds a,b,c and
/// `box_angles` alpha,beta,gamma in degrees (both empty when no box is
/// stored); `time` is 0.0 when the file has no time variable.
pub(all) struct NcdfFrame {
  coords : Array[Double]
  box_lengths : Array[Double]
  box_angles : Array[Double]
  time : Double
}

///|
/// Internal variable description used while building the header.
priv struct NcVarDef {
  name : String
  nctype : Int32
  dimids : Array[Int32]
  nelems : Int32
}

// ---------------------------------------------------------------------------
// Layout helpers
// ---------------------------------------------------------------------------

///|
/// Round `n` up to the next multiple of 4.
fn nc_pad4(n : Int32) -> Int32 {
  (n + 3) / 4 * 4
}

///|
/// Padded on-disk size (bytes) of `nelems` values of type `nctype`:
/// raw byte count rounded up to a 4-byte boundary (NC_BYTE=1, NC_CHAR=2,
/// NC_SHORT=3, NC_INT=4, NC_FLOAT=5, NC_DOUBLE=6).
fn nc_vsize(nctype : Int32, nelems : Int32) -> Int32 {
  let raw = if nctype == 1 || nctype == 2 {
    nelems
  } else if nctype == 3 {
    nelems * 2
  } else if nctype == 6 {
    nelems * 8
  } else {
    nelems * 4
  }
  nc_pad4(raw)
}

///|
/// Look up a variable by name, returning its index or -1.
fn nc_find_var(names : Array[String], target : String) -> Int32 {
  for i = 0; i < names.length(); i = i + 1 {
    if names[i] == target {
      return i
    }
  }
  -1
}

// ---------------------------------------------------------------------------
// IEEE 754 bit conversion (math decomposition, see dcd_f32_bits)
// ---------------------------------------------------------------------------

///|
fn nc_pow2(e : Int32) -> Double {
  if e >= 0 {
    let mut result = 1.0
    for i = 0; i < e; i = i + 1 {
      result = result * 2.0
    }
    result
  } else {
    let mut result = 1.0
    for i = 0; i < -e; i = i + 1 {
      result = result / 2.0
    }
    result
  }
}

///|
fn nc_f32_bits(v : Double) -> Int32 {
  if v != v {
    return 0x7FC00000
  }
  if v == 0.0 {
    return if 1.0 / v < 0.0 { 0x80000000 } else { 0 }
  }
  let sign = if v < 0.0 { 1 } else { 0 }
  let a = if v < 0.0 { 0.0 - v } else { v }
  if a > 3.4028234663852886e38 {
    return (sign << 31) | 0x7F800000
  }
  let mut e = 0
  let mut m = a
  if m >= 1.0 {
    while m >= 2.0 && e < 128 {
      m = m / 2.0
      e = e + 1
    }
  } else {
    while m < 1.0 && e > -127 {
      m = m * 2.0
      e = e - 1
    }
  }
  let mut exp = e + 127
  if exp <= 0 {
    return sign << 31
  }
  let mut mant = ((m - 1.0) * 8388608.0 + 0.5).to_int()
  if mant >= 8388608 {
    mant = 0
    exp = exp + 1
  }
  if exp >= 255 {
    return (sign << 31) | 0x7F800000
  }
  (sign << 31) | (exp << 23) | (mant & 0x7FFFFF)
}

///|
fn nc_f32_from_bits(bits : Int32) -> Double {
  let sign = if ((bits >> 31) & 1) == 1 { -1.0 } else { 1.0 }
  let exp = (bits >> 23) & 0xFF
  let mant = bits & 0x7FFFFF
  if exp == 0 {
    if mant == 0 {
      return sign * 0.0
    }
    return sign * mant.to_double() * nc_pow2(-149)
  }
  if exp == 255 {
    if mant == 0 {
      return sign * 1.0e308 * 10.0
    }
    return 0.0 / 0.0
  }
  sign * (1.0 + mant.to_double() / 8388608.0) * nc_pow2(exp - 127)
}

///|
fn nc_f64_bits(v : Double) -> UInt64 {
  if v != v {
    return 0x7FF8000000000000UL
  }
  if v == 0.0 {
    return if 1.0 / v < 0.0 { 0x8000000000000000UL } else { 0UL }
  }
  let sign = if v < 0.0 { 1UL } else { 0UL }
  let a = if v < 0.0 { 0.0 - v } else { v }
  if a > 1.7976931348623157e308 {
    return (sign << 63) | 0x7FF0000000000000UL
  }
  let mut e = 0
  let mut m = a
  if m >= 1.0 {
    while m >= 2.0 && e < 1024 {
      m = m / 2.0
      e = e + 1
    }
  } else {
    while m < 1.0 && e > -1023 {
      m = m * 2.0
      e = e - 1
    }
  }
  let mut exp = e + 1023
  if exp <= 0 {
    return sign << 63
  }
  let mut mant = ((m - 1.0) * 4503599627370496.0).to_uint64()
  if mant >= 4503599627370496UL {
    mant = 0UL
    exp = exp + 1
  }
  if exp >= 2047 {
    return (sign << 63) | 0x7FF0000000000000UL
  }
  (sign << 63) | (exp.to_uint64() << 52) | (mant & 0xFFFFFFFFFFFFFUL)
}

///|
fn nc_f64_from_bits(bits : UInt64) -> Double {
  let sign = if ((bits >> 63) & 1UL) == 1UL { -1.0 } else { 1.0 }
  let exp = ((bits >> 52) & 0x7FFUL).to_int()
  let mant = bits & 0xFFFFFFFFFFFFFUL
  if exp == 0 {
    if mant == 0UL {
      return sign * 0.0
    }
    return sign * mant.to_double() * nc_pow2(-1074)
  }
  if exp == 2047 {
    if mant == 0UL {
      return sign * 1.0e308 * 10.0
    }
    return 0.0 / 0.0
  }
  sign * (1.0 + mant.to_double() / 4503599627370496.0) * nc_pow2(exp - 1023)
}

// ---------------------------------------------------------------------------
// Big-endian read helpers (byte buffer = Array[Int32], values 0-255)
// ---------------------------------------------------------------------------

///|
fn nc_read_i32(
  data : Array[Int32],
  pos : Int32,
) -> (Int32, Int32) raise NcdfError {
  if pos + 4 > data.length() {
    raise NcdfError("unexpected end of data reading int at pos \{pos}")
  }
  let b0 = data[pos]
  let b1 = data[pos + 1]
  let b2 = data[pos + 2]
  let b3 = data[pos + 3]
  ((b0 << 24) | (b1 << 16) | (b2 << 8) | b3, pos + 4)
}

///|
fn nc_read_f32(
  data : Array[Int32],
  pos : Int32,
) -> (Double, Int32) raise NcdfError {
  if pos + 4 > data.length() {
    raise NcdfError("unexpected end of data reading float32 at pos \{pos}")
  }
  let b0 = data[pos]
  let b1 = data[pos + 1]
  let b2 = data[pos + 2]
  let b3 = data[pos + 3]
  let bits = (b0 << 24) | (b1 << 16) | (b2 << 8) | b3
  (nc_f32_from_bits(bits), pos + 4)
}

///|
fn nc_read_f64(
  data : Array[Int32],
  pos : Int32,
) -> (Double, Int32) raise NcdfError {
  if pos + 8 > data.length() {
    raise NcdfError("unexpected end of data reading float64 at pos \{pos}")
  }
  let b0 = data[pos].to_uint64()
  let b1 = data[pos + 1].to_uint64()
  let b2 = data[pos + 2].to_uint64()
  let b3 = data[pos + 3].to_uint64()
  let b4 = data[pos + 4].to_uint64()
  let b5 = data[pos + 5].to_uint64()
  let b6 = data[pos + 6].to_uint64()
  let b7 = data[pos + 7].to_uint64()
  let bits = (b0 << 56) |
    (b1 << 48) |
    (b2 << 40) |
    (b3 << 32) |
    (b4 << 24) |
    (b5 << 16) |
    (b6 << 8) |
    b7
  (nc_f64_from_bits(bits), pos + 8)
}

///|
/// Read a pascal-style name: int32 length + bytes + zero padding to a
/// 4-byte boundary (padding computed on the byte payload only).
fn nc_read_name(
  data : Array[Int32],
  pos : Int32,
) -> (String, Int32) raise NcdfError {
  let (n, p1) = nc_read_i32(data, pos)
  if n < 0 {
    raise NcdfError("invalid name length \{n} at pos \{pos}")
  }
  if p1 + n > data.length() {
    raise NcdfError("unexpected end of data reading name at pos \{pos}")
  }
  let sb = StringBuilder()
  for i = 0; i < n; i = i + 1 {
    sb.write_char(data[p1 + i].to_uint16().unsafe_to_char())
  }
  let end = p1 + n + (nc_pad4(n) - n)
  (sb.to_string(), end)
}

// ---------------------------------------------------------------------------
// Big-endian write helpers
// ---------------------------------------------------------------------------

///|
fn nc_write_i32(buf : Array[Int32], v : Int32) -> Unit {
  buf.push((v >> 24) & 0xFF)
  buf.push((v >> 16) & 0xFF)
  buf.push((v >> 8) & 0xFF)
  buf.push(v & 0xFF)
}

///|
fn nc_write_f32(buf : Array[Int32], v : Double) -> Unit {
  nc_write_i32(buf, nc_f32_bits(v))
}

///|
fn nc_write_f64(buf : Array[Int32], v : Double) -> Unit {
  let bits = nc_f64_bits(v)
  buf.push(((bits >> 56) & 0xFFUL).to_int())
  buf.push(((bits >> 48) & 0xFFUL).to_int())
  buf.push(((bits >> 40) & 0xFFUL).to_int())
  buf.push(((bits >> 32) & 0xFFUL).to_int())
  buf.push(((bits >> 24) & 0xFFUL).to_int())
  buf.push(((bits >> 16) & 0xFFUL).to_int())
  buf.push(((bits >> 8) & 0xFFUL).to_int())
  buf.push((bits & 0xFFUL).to_int())
}

///|
/// Write a pascal-style name: int32 length + bytes + zero padding to a
/// 4-byte boundary (padding computed on the byte payload only).
fn nc_write_name(buf : Array[Int32], s : String) -> Unit {
  let n = s.length()
  nc_write_i32(buf, n)
  for i = 0; i < n; i = i + 1 {
    buf.push(s.unsafe_get(i).to_int() & 0xFF)
  }
  let pad = nc_pad4(n) - n
  for i = 0; i < pad; i = i + 1 {
    buf.push(0)
  }
}

///|
/// Write char variable/attribute payload: bytes + zero padding to a
/// 4-byte boundary.
fn nc_write_char_data(buf : Array[Int32], s : String) -> Unit {
  let n = s.length()
  for i = 0; i < n; i = i + 1 {
    buf.push(s.unsafe_get(i).to_int() & 0xFF)
  }
  let pad = nc_pad4(n) - n
  for i = 0; i < pad; i = i + 1 {
    buf.push(0)
  }
}

// ---------------------------------------------------------------------------
// Write / read
// ---------------------------------------------------------------------------

///|
/// Serialize AMBER-style frames into a classic NetCDF CDF-1 byte stream.
/// Always writes dimensions frame (unlimited), atom, spatial plus the
/// variables coordinates (float) and time (double) and the char variable
/// spatial ("xyz"). If any frame carries a box, also writes dimensions
/// cell_spatial/cell_angular/angular/label and variables cell_lengths,
/// cell_angles, cell_spatial ("abc") and cell_angular ("alphabeta gamma").
pub fn ncdf_write_amber(
  frames : Array[NcdfFrame],
  natoms : Int32,
) -> Array[Int32] {
  let nframes = frames.length()
  let mut has_box = false
  for f in frames {
    if f.box_lengths.length() >= 3 {
      has_box = true
    }
  }
  let dim_names : Array[String] = ["frame", "atom", "spatial"]
  let dim_sizes : Array[Int32] = [0, natoms, 3]
  if has_box {
    dim_names.push("cell_spatial")
    dim_sizes.push(3)
    dim_names.push("cell_angular")
    dim_sizes.push(3)
    dim_names.push("angular")
    dim_sizes.push(3)
    dim_names.push("label")
    dim_sizes.push(5)
  }
  // Variables in data order: non-record variables first, then record ones.
  let vars : Array[NcVarDef] = Array(capacity=0)
  vars.push(NcVarDef::{ name: "spatial", nctype: 2, dimids: [2], nelems: 3, })
  if has_box {
    vars.push(NcVarDef::{
      name: "cell_spatial",
      nctype: 2,
      dimids: [3],
      nelems: 3,
    })
    vars.push(NcVarDef::{
      name: "cell_angular",
      nctype: 2,
      dimids: [4, 6],
      nelems: 15,
    })
  }
  vars.push(NcVarDef::{
    name: "coordinates",
    nctype: 5,
    dimids: [0, 1, 2],
    nelems: natoms * 3,
  })
  if has_box {
    vars.push(NcVarDef::{
      name: "cell_lengths",
      nctype: 5,
      dimids: [0, 3],
      nelems: 3,
    })
    vars.push(NcVarDef::{
      name: "cell_angles",
      nctype: 5,
      dimids: [0, 4],
      nelems: 3,
    })
  }
  vars.push(NcVarDef::{ name: "time", nctype: 6, dimids: [0], nelems: 1, })
  // Exact header length in bytes.
  let mut hlen = 8
  hlen = hlen + 8
  for i = 0; i < dim_names.length(); i = i + 1 {
    hlen = hlen + 4 + nc_pad4(dim_names[i].length()) + 4
  }
  hlen = hlen + 4
  hlen = hlen + 8
  for v in vars {
    hlen = hlen + 4 + nc_pad4(v.name.length())
    hlen = hlen + 4
    hlen = hlen + 4 * v.dimids.length()
    hlen = hlen + 4
    hlen = hlen + 12
  }
  // begin offsets: non-record data right after the header, then records.
  let mut nonrec_total = 0
  for v in vars {
    if v.dimids.length() > 0 && v.dimids[0] != 0 {
      nonrec_total = nonrec_total + nc_vsize(v.nctype, v.nelems)
    }
  }
  let recstart = hlen + nonrec_total
  let begins : Array[Int32] = Array(capacity=vars.length())
  let mut c1 = hlen
  let mut c2 = recstart
  for v in vars {
    let vs = nc_vsize(v.nctype, v.nelems)
    if v.dimids.length() > 0 && v.dimids[0] == 0 {
      begins.push(c2)
      c2 = c2 + vs
    } else {
      begins.push(c1)
      c1 = c1 + vs
    }
  }
  let buf : Array[Int32] = Array(capacity=0)
  buf.push(67) // 'C'
  buf.push(68) // 'D'
  buf.push(70) // 'F'
  buf.push(1) // CDF-1
  nc_write_i32(buf, nframes)
  nc_write_i32(buf, 10) // NC_DIMENSION
  nc_write_i32(buf, dim_names.length())
  for i = 0; i < dim_names.length(); i = i + 1 {
    nc_write_name(buf, dim_names[i])
    nc_write_i32(buf, dim_sizes[i])
  }
  nc_write_i32(buf, 0) // empty global attribute list
  nc_write_i32(buf, 11) // NC_VARIABLE
  nc_write_i32(buf, vars.length())
  for i = 0; i < vars.length(); i = i + 1 {
    let v = vars[i]
    nc_write_name(buf, v.name)
    nc_write_i32(buf, v.dimids.length())
    for d in v.dimids {
      nc_write_i32(buf, d)
    }
    nc_write_i32(buf, 0) // empty variable attribute list
    nc_write_i32(buf, v.nctype)
    nc_write_i32(buf, nc_vsize(v.nctype, v.nelems))
    nc_write_i32(buf, begins[i])
  }
  // Non-record variable payloads.
  nc_write_char_data(buf, "xyz")
  if has_box {
    nc_write_char_data(buf, "abc")
    nc_write_char_data(buf, "alphabeta gamma")
  }
  // Interleaved record payloads.
  for f in frames {
    for a = 0; a < natoms; a = a + 1 {
      for d = 0; d < 3; d = d + 1 {
        let k = a * 3 + d
        let val = if k < f.coords.length() { f.coords[k] } else { 0.0 }
        nc_write_f32(buf, val)
      }
    }
    if has_box {
      for d = 0; d < 3; d = d + 1 {
        let val = if d < f.box_lengths.length() {
          f.box_lengths[d]
        } else {
          0.0
        }
        nc_write_f32(buf, val)
      }
      for d = 0; d < 3; d = d + 1 {
        let val = if d < f.box_angles.length() { f.box_angles[d] } else { 0.0 }
        nc_write_f32(buf, val)
      }
    }
    nc_write_f64(buf, f.time)
  }
  buf
}

///|
/// Parse a classic NetCDF CDF-1 byte stream written in AMBER conventions
/// and return all record frames. Raises NcdfError on malformed input (bad
/// magic, invalid list tags, missing coordinates variable, truncated data).
pub fn ncdf_read_amber(data : Array[Int32]) -> Array[NcdfFrame] raise NcdfError {
  let n = data.length()
  if n < 8 {
    raise NcdfError("data too short for a NetCDF file")
  }
  if data[0] != 67 || data[1] != 68 || data[2] != 70 || data[3] != 1 {
    raise NcdfError("invalid NetCDF magic, expected CDF-1")
  }
  let (numrecs, p1) = nc_read_i32(data, 4)
  let mut pos = p1
  let dim_names : Array[String] = Array(capacity=0)
  let dim_sizes : Array[Int32] = Array(capacity=0)
  let (dtag, p2) = nc_read_i32(data, pos)
  pos = p2
  if dtag == 10 {
    let (ndims, p3) = nc_read_i32(data, pos)
    pos = p3
    for i = 0; i < ndims; i = i + 1 {
      let (nm, p4) = nc_read_name(data, pos)
      pos = p4
      let (sz, p5) = nc_read_i32(data, pos)
      pos = p5
      dim_names.push(nm)
      dim_sizes.push(sz)
    }
  } else if dtag != 0 {
    raise NcdfError("invalid dimension list tag \{dtag}")
  }
  let (gtag, p6) = nc_read_i32(data, pos)
  pos = p6
  if gtag == 12 {
    let (ng, p7) = nc_read_i32(data, pos)
    pos = p7
    for i = 0; i < ng; i = i + 1 {
      let (_, p8) = nc_read_name(data, pos)
      pos = p8
      let (t, p9) = nc_read_i32(data, pos)
      pos = p9
      let (ne, p10) = nc_read_i32(data, pos)
      pos = p10
      pos = pos + nc_vsize(t, ne)
    }
  } else if gtag != 0 {
    raise NcdfError("invalid global attribute list tag \{gtag}")
  }
  let (vtag, p11) = nc_read_i32(data, pos)
  pos = p11
  if vtag != 11 {
    raise NcdfError("missing variable list (tag \{vtag})")
  }
  let (nvars, p12) = nc_read_i32(data, pos)
  pos = p12
  let var_names : Array[String] = Array(capacity=0)
  let var_types : Array[Int32] = Array(capacity=0)
  let var_vsizes : Array[Int32] = Array(capacity=0)
  let var_begins : Array[Int32] = Array(capacity=0)
  let var_dimids : Array[Array[Int32]] = Array(capacity=0)
  for i = 0; i < nvars; i = i + 1 {
    let (nm, pa) = nc_read_name(data, pos)
    pos = pa
    let (nd, pb) = nc_read_i32(data, pos)
    pos = pb
    let ids : Array[Int32] = Array(capacity=0)
    for k = 0; k < nd; k = k + 1 {
      let (id, pc) = nc_read_i32(data, pos)
      pos = pc
      ids.push(id)
    }
    let (atag, pd) = nc_read_i32(data, pos)
    pos = pd
    if atag == 12 {
      let (na, pe) = nc_read_i32(data, pos)
      pos = pe
      for k = 0; k < na; k = k + 1 {
        let (_, pf) = nc_read_name(data, pos)
        pos = pf
        let (t, pg) = nc_read_i32(data, pos)
        pos = pg
        let (ne, ph) = nc_read_i32(data, pos)
        pos = ph
        pos = pos + nc_vsize(t, ne)
      }
    } else if atag != 0 {
      raise NcdfError("invalid attribute list tag \{atag} in variable \{nm}")
    }
    let (t, pi) = nc_read_i32(data, pos)
    pos = pi
    let (vs, pj) = nc_read_i32(data, pos)
    pos = pj
    let (bg, pk) = nc_read_i32(data, pos)
    pos = pk
    var_names.push(nm)
    var_types.push(t)
    var_vsizes.push(vs)
    var_begins.push(bg)
    var_dimids.push(ids)
  }
  let coords_idx = nc_find_var(var_names, "coordinates")
  if coords_idx < 0 {
    raise NcdfError("coordinates variable not found")
  }
  if var_types[coords_idx] != 5 {
    raise NcdfError("coordinates must be NC_FLOAT")
  }
  let dimids = var_dimids[coords_idx]
  if dimids.length() < 3 {
    raise NcdfError("coordinates must have at least 3 dimensions")
  }
  if dimids[0] >= dim_sizes.length() || dim_sizes[dimids[0]] != 0 {
    raise NcdfError("coordinates first dimension must be unlimited")
  }
  if dimids[1] >= dim_sizes.length() {
    raise NcdfError("coordinates atom dimension out of range")
  }
  let natoms = dim_sizes[dimids[1]]
  let mut recsize = 0
  for i = 0; i < nvars; i = i + 1 {
    if var_dimids[i].length() > 0 && var_dimids[i][0] == 0 {
      recsize = recsize + var_vsizes[i]
    }
  }
  let time_idx = nc_find_var(var_names, "time")
  let cl_idx = nc_find_var(var_names, "cell_lengths")
  let ca_idx = nc_find_var(var_names, "cell_angles")
  let frames : Array[NcdfFrame] = Array(capacity=0)
  for r = 0; r < numrecs; r = r + 1 {
    let base = var_begins[coords_idx] + r * recsize
    let mut p = base
    let coords : Array[Double] = Array(capacity=natoms * 3)
    for k = 0; k < natoms * 3; k = k + 1 {
      let (v, pn) = nc_read_f32(data, p)
      p = pn
      coords.push(v)
    }
    let mut time = 0.0
    if time_idx >= 0 && var_types[time_idx] == 6 {
      let (tv, _) = nc_read_f64(data, var_begins[time_idx] + r * recsize)
      time = tv
    }
    let bl : Array[Double] = Array(capacity=0)
    let ba : Array[Double] = Array(capacity=0)
    if cl_idx >= 0 &&
      ca_idx >= 0 &&
      var_types[cl_idx] == 5 &&
      var_types[ca_idx] == 5 {
      let mut q = var_begins[cl_idx] + r * recsize
      for k = 0; k < 3; k = k + 1 {
        let (v, qn) = nc_read_f32(data, q)
        q = qn
        bl.push(v)
      }
      let mut q2 = var_begins[ca_idx] + r * recsize
      for k = 0; k < 3; k = k + 1 {
        let (v, qn2) = nc_read_f32(data, q2)
        q2 = qn2
        ba.push(v)
      }
    }
    frames.push(NcdfFrame::{ coords, box_lengths: bl, box_angles: ba, time, })
  }
  frames
}

///|
/// Convert NetCDF frames into a Trajectory. Frames whose `coords` array
/// is empty are skipped.
pub fn ncdf_to_trajectory(frames : Array[NcdfFrame]) -> Trajectory {
  let out : Array[Array[Array[Double]]] = Array(capacity=0)
  for f in frames {
    if f.coords.length() >= 3 {
      let natoms = f.coords.length() / 3
      let frame : Array[Array[Double]] = Array(capacity=natoms)
      for a = 0; a < natoms; a = a + 1 {
        frame.push([f.coords[a * 3], f.coords[a * 3 + 1], f.coords[a * 3 + 2]])
      }
      out.push(frame)
    }
  }
  Trajectory::new(out)
}
