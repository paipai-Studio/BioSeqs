///|
/// illuminaio-style IDAT raw file parsing (Illumina BeadArray).
///
/// Port of Bioconductor illuminaio's `readIDAT_nonenc` — the non-encrypted
/// IDAT format (version >= 3) produced by Illumina Infinium methylation
/// arrays (450K / EPIC). The file is a set of typed fields located by a
/// directory of `(fieldCode, byteOffset)` pairs: summary intensities
/// (Mean/SD/NBeads) and the per-probe Illumina bead ID, plus string metadata.
///
/// Encrypted expression-array IDAT files (version 1) are out of scope.

///|
pub suberror IdatError {
  IdatError(String)
}

///|
fn[T] idat_fail(msg : String) -> T raise IdatError {
  raise IdatError(msg)
}

///|
pub struct IdatResult {
  version : Int
  n_fields : Int
  n_snps_read : Int
  illumina_id : Array[Int]
  mean : Array[Int]
  sd : Array[Int]
  n_beads : Array[Int]
  barcode : String
  chip_type : String
}

///|
fn idat_u8(data : String, pos : Int) -> Int {
  data.unsafe_get(pos).to_int() & 0xFF
}

///|
fn idat_u16(data : String, pos : Int) -> Int {
  let b0 = data.unsafe_get(pos).to_int() & 0xFF
  let b1 = data.unsafe_get(pos + 1).to_int() & 0xFF
  b0 | (b1 << 8)
}

///|
fn idat_i32(data : String, pos : Int) -> Int {
  let b0 = data.unsafe_get(pos).to_int() & 0xFF
  let b1 = data.unsafe_get(pos + 1).to_int() & 0xFF
  let b2 = data.unsafe_get(pos + 2).to_int() & 0xFF
  let b3 = data.unsafe_get(pos + 3).to_int() & 0xFF
  b0 | (b1 << 8) | (b2 << 16) | (b3 << 24)
}

///|
/// Read a little-endian int64 (unsigned; positive values in practice).
fn idat_i64(data : String, pos : Int) -> Int {
  let mut v = 0
  let mut k = 0
  while k < 8 {
    v = v | ((data.unsafe_get(pos + k).to_int() & 0xFF) << (8 * k))
    k = k + 1
  }
  v
}

///|
/// Read the variable-length byte count (7-bit groups, MSB = continuation),
/// matching illuminaio's `readBytesToRead`.
fn idat_read_var_length(data : String, pos : Int) -> (Int, Int) raise IdatError {
  if pos >= data.length() {
    idat_fail("Truncated IDAT length prefix")
  }
  let mut m = idat_u8(data, pos)
  let mut p = pos + 1
  let mut n = m % 128
  let mut shift = 0
  while m / 128 == 1 {
    if p >= data.length() {
      idat_fail("Truncated IDAT length prefix")
    }
    m = idat_u8(data, p)
    p = p + 1
    shift = shift + 7
    let k = (m % 128) * (1 << shift)
    n = n + k
  }
  (n, p)
}

///|
/// Read a length-prefixed string field.
fn idat_read_string(data : String, pos : Int) -> (String, Int) raise IdatError {
  let (n, p) = idat_read_var_length(data, pos)
  if p + n > data.length() {
    idat_fail("Truncated IDAT string field")
  }
  let buf = StringBuilder::new()
  let mut k = 0
  while k < n {
    buf.write_char(data.unsafe_get(p + k).unsafe_to_char())
    k = k + 1
  }
  (buf.to_string(), p + n)
}

///|
fn idat_field_offset(
  offsets : Map[Int, Int],
  code : Int,
) -> Int? {
  offsets.get(code)
}

///|
/// Parse a non-encrypted IDAT file (the whole buffer).
pub fn parse_idat(data : String) -> IdatResult raise IdatError {
  if data.length() < 16 {
    idat_fail("IDAT file too short")
  }
  // magic
  if data.unsafe_get(0).to_int() != 'I'.to_int() ||
  data.unsafe_get(1).to_int() != 'D'.to_int() ||
  data.unsafe_get(2).to_int() != 'A'.to_int() ||
  data.unsafe_get(3).to_int() != 'T'.to_int() {
    idat_fail("Not an IDAT file (bad magic)")
  }
  let version = idat_i64(data, 4)
  let n_fields = idat_i32(data, 12)

  let offsets : Map[Int, Int] = Map([], capacity=n_fields)
  let mut p = 16
  for _i = 0; _i < n_fields; _i = _i + 1 {
    if p + 10 > data.length() {
      idat_fail("Truncated IDAT field directory")
    }
    let code = idat_u16(data, p)
    let offset = idat_i64(data, p + 2)
    offsets.set(code, offset)
    p = p + 10
  }

  // nSNPsRead
  let ns_off = idat_field_offset(offsets, 1000).unwrap_or(0)
  let n_snps_read = idat_i32(data, ns_off)

  let read_i32_array = fn(off : Int) -> Array[Int] raise IdatError {
    let arr : Array[Int] = Array(capacity=0)
    let mut k = 0
    while k < n_snps_read {
      arr.push(idat_i32(data, off + 4 * k))
      k = k + 1
    }
    arr
  }

  let read_u16_array = fn(off : Int) -> Array[Int] raise IdatError {
    let arr : Array[Int] = Array(capacity=0)
    let mut k = 0
    while k < n_snps_read {
      arr.push(idat_u16(data, off + 2 * k))
      k = k + 1
    }
    arr
  }

  let read_u8_array = fn(off : Int) -> Array[Int] raise IdatError {
    let arr : Array[Int] = Array(capacity=0)
    let mut k = 0
    while k < n_snps_read {
      arr.push(idat_u8(data, off + k))
      k = k + 1
    }
    arr
  }

  let illumina_id = match idat_field_offset(offsets, 102) {
    Some(off) => read_i32_array(off)
    None => []
  }
  let sd = match idat_field_offset(offsets, 103) {
    Some(off) => read_u16_array(off)
    None => []
  }
  let mean = match idat_field_offset(offsets, 104) {
    Some(off) => read_u16_array(off)
    None => []
  }
  let n_beads = match idat_field_offset(offsets, 107) {
    Some(off) => read_u8_array(off)
    None => []
  }
  let barcode = match idat_field_offset(offsets, 402) {
    Some(off) => idat_read_string(data, off).0
    None => ""
  }
  let chip_type = match idat_field_offset(offsets, 403) {
    Some(off) => idat_read_string(data, off).0
    None => ""
  }

  IdatResult::{
    version,
    n_fields,
    n_snps_read,
    illumina_id,
    mean,
    sd,
    n_beads,
    barcode,
    chip_type,
  }
}

///|
/// A sample IDAT fixture (3 probes) for tests and demos, base64-encoded.
pub fn idat_sample_base64() -> String {
  "SURBVAMAAAAAAAAABwAAAOgDVgAAAAAAAABmAFoAAAAAAAAAZwBmAAAAAAAAAGgAbAAAAAAAAABrAHIAAAAAAAAAkgF1AAAAAAAAAJMBgAAAAAAAAAADAAAAZQAAAMoAAAAvAQAACgAUAB4A6APQB7gLBQYHCjEyMzQ1Njc4OTATSHVtYW5NZXRoeWxhdGlvbjQ1MA=="
}