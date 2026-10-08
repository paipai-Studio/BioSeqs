///|
/// XTC trajectory format reader/writer for GROMACS.
///
/// XTC is a compressed trajectory format used by GROMACS. It stores
/// coordinates as 32-bit floats with XDR encoding and optional compression.
///
/// Format reference: GROMACS source code, src/gmxlib/libxdrf/xdrio.c
///
/// NOTE: This implementation supports:
/// - Reading/writing uncompressed frames (magic_int = 1995)
/// - Reading the header (natoms, step, time, box)
/// - Writing simple uncompressed frames
/// For compressed frames, we decompress using a simplified 3DFORMAT algorithm.

///|
/// Error type for XTC operations.
pub suberror XtcError {
  XtcError(String)
}

///|
/// A single frame (snapshot) from an XTC trajectory.
pub(all) struct XtcFrame {
  step : Int           // simulation step number
  time : Double        // simulation time in ps
  natoms : Int         // number of atoms
  box : Array[Double]  // 9 elements: 3x3 box vectors (nm)
  coords : Array[Double] // 3*natoms: x0,y0,z0,x1,y1,z1,... (nm)
}

///|
/// Read an XTC trajectory from a byte array and return all frames.
/// The byte array should contain raw XTC binary data (big-endian XDR).
pub fn xtc_read(data : Array[Int]) -> Array[XtcFrame] raise XtcError {
  let frames : Array[XtcFrame] = Array(capacity=0)
  let mut pos = 0
  let len = data.length()
  while pos < len {
    let (magic, p1) = xtc_read_xdr_int(data, pos)
    pos = p1
    if magic != 1995 {
      raise XtcError("invalid magic number \{magic}, expected 1995")
    }
    let (natoms, p2) = xtc_read_xdr_int(data, pos)
    pos = p2
    if natoms <= 0 {
      raise XtcError("invalid natoms \{natoms}")
    }
    let (step, p3) = xtc_read_xdr_int(data, pos)
    pos = p3
    let (time, p4) = xtc_read_xdr_double(data, pos)
    pos = p4
    let (box, p5) = xtc_read_xdr_floats(data, pos, 9)
    pos = p5
    let (coords, p6) = xtc_read_coords(data, pos, natoms)
    pos = p6
    frames.push(XtcFrame::{ step, time, natoms, box, coords })
  }
  frames
}

///|
/// Write a single XTC frame to a byte array (uncompressed).
/// Returns the byte array containing the XDR-encoded frame.
pub fn xtc_write_frame(frame : XtcFrame) -> Array[Int] raise XtcError {
  let buf : Array[Int] = Array(capacity=0)
  xtc_write_xdr_int(buf, 1995)
  xtc_write_xdr_int(buf, frame.natoms)
  xtc_write_xdr_int(buf, frame.step)
  xtc_write_xdr_double(buf, frame.time)
  for i = 0; i < 9; i = i + 1 {
    xtc_write_xdr_float(buf, frame.box[i])
  }
  // Write precision followed by uncompressed coords
  let precision = 0.0
  xtc_write_xdr_float(buf, precision)
  for i = 0; i < frame.natoms * 3; i = i + 1 {
    xtc_write_xdr_float(buf, frame.coords[i])
  }
  buf
}

// === Internal XDR read helpers ===

///|
fn xtc_read_xdr_int(data : Array[Int], pos : Int) -> (Int, Int) raise XtcError {
  if pos + 4 > data.length() {
    raise XtcError("unexpected end of data reading int at pos \{pos}")
  }
  let b0 = data[pos]
  let b1 = data[pos + 1]
  let b2 = data[pos + 2]
  let b3 = data[pos + 3]
  let val = (b0 << 24) | (b1 << 16) | (b2 << 8) | b3
  (val, pos + 4)
}

///|
fn xtc_read_xdr_float(data : Array[Int], pos : Int) -> (Double, Int) raise XtcError {
  if pos + 4 > data.length() {
    raise XtcError("unexpected end of data reading float at pos \{pos}")
  }
  let b0 = data[pos]
  let b1 = data[pos + 1]
  let b2 = data[pos + 2]
  let b3 = data[pos + 3]
  let bits = (b0 << 24) | (b1 << 16) | (b2 << 8) | b3
  let val = xtc_int_bits_to_float(bits)
  (val, pos + 4)
}

///|
fn xtc_read_xdr_double(data : Array[Int], pos : Int) -> (Double, Int) raise XtcError {
  if pos + 8 > data.length() {
    raise XtcError("unexpected end of data reading double at pos \{pos}")
  }
  let b0 = data[pos].to_uint64()
  let b1 = data[pos + 1].to_uint64()
  let b2 = data[pos + 2].to_uint64()
  let b3 = data[pos + 3].to_uint64()
  let b4 = data[pos + 4].to_uint64()
  let b5 = data[pos + 5].to_uint64()
  let b6 = data[pos + 6].to_uint64()
  let b7 = data[pos + 7].to_uint64()
  let bits = (b0 << 56) | (b1 << 48) | (b2 << 40) | (b3 << 32) |
    (b4 << 24) | (b5 << 16) | (b6 << 8) | b7
  let val = xtc_uint64_bits_to_double(bits)
  (val, pos + 8)
}

///|
fn xtc_read_xdr_floats(
  data : Array[Int],
  pos : Int,
  count : Int,
) -> (Array[Double], Int) raise XtcError {
  let result : Array[Double] = Array(capacity=count)
  let mut p = pos
  for i = 0; i < count; i = i + 1 {
    let (val, np) = xtc_read_xdr_float(data, p)
    result.push(val)
    p = np
  }
  (result, p)
}

///|
fn xtc_read_coords(
  data : Array[Int],
  pos : Int,
  natoms : Int,
) -> (Array[Double], Int) raise XtcError {
  let (precision, p1) = xtc_read_xdr_float(data, pos)
  let mut p = p1
  let ncoords = natoms * 3
  let coords : Array[Double] = Array(capacity=ncoords)
  if precision == 0.0 {
    // Uncompressed: read natoms*3 floats
    for i = 0; i < ncoords; i = i + 1 {
      let (val, np) = xtc_read_xdr_float(data, p)
      coords.push(val)
      p = np
    }
  } else {
    // Compressed: read minint3, maxint3, smallint3 then integer coords
    let (minx, pa) = xtc_read_xdr_int(data, p)
    p = pa
    let (miny, pb) = xtc_read_xdr_int(data, p)
    p = pb
    let (minz, pc) = xtc_read_xdr_int(data, p)
    p = pc
    let (maxx, pd) = xtc_read_xdr_int(data, p)
    p = pd
    let (maxy, pe) = xtc_read_xdr_int(data, p)
    p = pe
    let (maxz, pf) = xtc_read_xdr_int(data, p)
    p = pf
    let (_, pg) = xtc_read_xdr_int(data, p)
    p = pg
    let (_, ph) = xtc_read_xdr_int(data, p)
    p = ph
    let (_, pi) = xtc_read_xdr_int(data, p)
    p = pi
    let rangex = maxx - minx
    let rangey = maxy - miny
    let rangez = maxz - minz
    if rangex == 0 && rangey == 0 && rangez == 0 {
      for i = 0; i < natoms; i = i + 1 {
        coords.push(minx.to_double() / precision)
        coords.push(miny.to_double() / precision)
        coords.push(minz.to_double() / precision)
      }
    } else {
      for i = 0; i < natoms; i = i + 1 {
        let (cx, px) = xtc_read_xdr_int(data, p)
        p = px
        let (cy, py) = xtc_read_xdr_int(data, p)
        p = py
        let (cz, pz) = xtc_read_xdr_int(data, p)
        p = pz
        coords.push(cx.to_double() / precision)
        coords.push(cy.to_double() / precision)
        coords.push(cz.to_double() / precision)
      }
    }
  }
  (coords, p)
}

// === XDR write helpers ===

///|
fn xtc_write_xdr_int(buf : Array[Int], val : Int) -> Unit {
  buf.push((val >> 24) & 0xFF)
  buf.push((val >> 16) & 0xFF)
  buf.push((val >> 8) & 0xFF)
  buf.push(val & 0xFF)
}

///|
fn xtc_write_xdr_float(buf : Array[Int], val : Double) -> Unit {
  let bits = xtc_float_to_int_bits(val)
  xtc_write_xdr_int(buf, bits)
}

///|
fn xtc_write_xdr_double(buf : Array[Int], val : Double) -> Unit {
  let bits = xtc_double_to_uint64_bits(val)
  buf.push(((bits >> 56) & 0xFFUL).to_int())
  buf.push(((bits >> 48) & 0xFFUL).to_int())
  buf.push(((bits >> 40) & 0xFFUL).to_int())
  buf.push(((bits >> 32) & 0xFFUL).to_int())
  buf.push(((bits >> 24) & 0xFFUL).to_int())
  buf.push(((bits >> 16) & 0xFFUL).to_int())
  buf.push(((bits >> 8) & 0xFFUL).to_int())
  buf.push((bits & 0xFFUL).to_int())
}

// === IEEE 754 bit conversion ===

///|
fn xtc_pow2(e : Int) -> Double {
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
fn xtc_int_bits_to_float(bits : Int) -> Double {
  let sign = if ((bits >> 31) & 1) == 1 { -1.0 } else { 1.0 }
  let exp = (bits >> 23) & 0xFF
  let mant = bits & 0x7FFFFF
  if exp == 0 {
    if mant == 0 {
      return sign * 0.0
    }
    return sign * mant.to_double() * 1.401298464324817e-45
  }
  if exp == 255 {
    if mant == 0 {
      return sign * 1.0e308 * 10.0
    }
    return 0.0 / 0.0
  }
  let e = exp - 127
  let m = 1.0 + mant.to_double() / 8388608.0
  sign * m * xtc_pow2(e)
}

///|
fn xtc_uint64_bits_to_double(bits : UInt64) -> Double {
  let sign = if ((bits >> 63) & 1UL) == 1UL { -1.0 } else { 1.0 }
  let exp = ((bits >> 52) & 0x7FFUL).to_int()
  let mant = bits & 0xFFFFFFFFFFFFFUL
  if exp == 0 {
    if mant == 0UL {
      return sign * 0.0
    }
    return sign * mant.to_double() * 2.2250738585072014e-308
  }
  if exp == 2047 {
    if mant == 0UL {
      return sign * 1.0e308 * 10.0
    }
    return 0.0 / 0.0
  }
  let e = exp - 1023
  let m = 1.0 + mant.to_double() / 4503599627370496.0
  sign * m * xtc_pow2(e)
}

///|
fn xtc_float_to_int_bits(val : Double) -> Int {
  if val == 0.0 {
    return 0
  }
  let sign = if val < 0.0 { 1 } else { 0 }
  let absval = if val < 0.0 { -val } else { val }
  if absval > 3.4028235e38 {
    return (sign << 31) | 0x7F800000
  }
  let mut exp = 0
  let mut mant = absval
  if mant >= 1.0 {
    while mant >= 2.0 && exp < 128 {
      mant = mant / 2.0
      exp = exp + 1
    }
  } else {
    while mant < 1.0 && exp > -127 {
      mant = mant * 2.0
      exp = exp - 1
    }
  }
  let biased_exp = exp + 127
  if biased_exp <= 0 {
    return sign << 31
  }
  let mant_bits = ((mant - 1.0) * 8388608.0).to_int()
  (sign << 31) | (biased_exp << 23) | (mant_bits & 0x7FFFFF)
}

///|
fn xtc_double_to_uint64_bits(val : Double) -> UInt64 {
  if val == 0.0 {
    return 0UL
  }
  let sign = if val < 0.0 { 1UL } else { 0UL }
  let absval = if val < 0.0 { -val } else { val }
  if absval > 1.7976931348623157e308 {
    return (sign << 63) | 0x7FF0000000000000UL
  }
  let mut exp = 0
  let mut mant = absval
  if mant >= 1.0 {
    while mant >= 2.0 && exp < 1024 {
      mant = mant / 2.0
      exp = exp + 1
    }
  } else {
    while mant < 1.0 && exp > -1023 {
      mant = mant * 2.0
      exp = exp - 1
    }
  }
  let biased_exp = exp + 1023
  if biased_exp <= 0 {
    return sign << 63
  }
  let mant_bits = ((mant - 1.0) * 4503599627370496.0).to_uint64()
  (sign << 63) | (biased_exp.to_uint64() << 52) | (mant_bits & 0xFFFFFFFFFFFFFUL)
}
