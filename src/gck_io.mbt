// Bio.SeqIO.GckIO - Parser and writer for the Gene Construction Kit (GCK)
// binary format, mirroring BioPython's Bio.SeqIO.GckIO.
//
// GCK is a binary format used by Gene Construction Kit software (Textco
// BioSoftware) for storing plasmid and sequence data with features.
//
// All multi-byte integers are BIG-ENDIAN.
//
// File layout:
//   Bytes 0-7   : Magic "GCK_FILE"
//   Byte  8     : Format version (1)
//   Byte  9     : Sequence type (0=DNA, 1=RNA, 2=Protein)
//   Byte  10    : Topology (0=linear, 1=circular)
//   Byte  11    : Reserved (0)
//   Bytes 12-13 : Feature count (uint16)
//   Then       : Feature data (variable, see below)
//   Then       : Sequence length (uint32)
//   Then       : Sequence bytes (ASCII, 1 byte per base)
//
// Feature layout (repeated feature_count times):
//   Bytes 0-1   : Name length (uint16)
//   Then        : Name bytes (ASCII)
//   Bytes 0-1   : Type length (uint16)
//   Then        : Type bytes (ASCII, e.g. "CDS", "promoter")
//   Byte  0     : Direction (0=none, 1=forward, 2=reverse)
//   Bytes 0-1   : Segment count (uint16)
//   For each segment:
//     Bytes 0-3 : Start (uint32, 0-based)
//     Bytes 0-3 : End (uint32, 0-based, inclusive)

// ============================================================================
// Data structures
// ============================================================================

///|
/// The molecular type stored in a GCK file.
pub enum GckSeqType {
  GckDna
  GckRna
  GckProtein
  GckUnknown
} derive(Eq)

///|
/// Convert a GCK sequence type to its integer code.
pub fn GckSeqType::to_int(self : GckSeqType) -> Int {
  match self {
    GckDna => 0
    GckRna => 1
    GckProtein => 2
    GckUnknown => 3
  }
}

///|
/// Convert an integer code to a GCK sequence type.
pub fn GckSeqType::from_int(code : Int) -> GckSeqType {
  match code {
    0 => GckDna
    1 => GckRna
    2 => GckProtein
    _ => GckUnknown
  }
}

///|
/// Convert a GCK sequence type to a human-readable string.
pub fn GckSeqType::to_string(self : GckSeqType) -> String {
  match self {
    GckDna => "DNA"
    GckRna => "RNA"
    GckProtein => "Protein"
    GckUnknown => "Unknown"
  }
}

///|
/// A feature annotation in a GCK file.
pub struct GckFeature {
  name : String
  type_ : String // e.g. "CDS", "promoter", "terminator"
  direction : Int // 0=none, 1=forward, 2=reverse
  segments : Array[(Int, Int)] // (start, end) 0-based inclusive
}

///|
/// Construct a GckFeature.
pub fn GckFeature::new(
  name~ : String,
  type_~ : String,
  direction~ : Int,
  segments~ : Array[(Int, Int)],
) -> GckFeature {
  GckFeature::{ name, type_, direction, segments }
}

///|
/// A parsed GCK file.
pub struct GckFile {
  pub mut version : Int
  pub mut seq_type : GckSeqType
  pub mut is_circular : Bool
  pub mut sequence : String
  features : Array[GckFeature]
}

///|
/// Construct an empty GckFile.
pub fn GckFile::new() -> GckFile {
  GckFile::{
    version: 1,
    seq_type: GckDna,
    is_circular: false,
    sequence: "",
    features: [],
  }
}

///|
/// Set the version.
pub fn GckFile::set_version(self : GckFile, v : Int) -> Unit {
  self.version = v
}

///|
/// Set the sequence type.
pub fn GckFile::set_seq_type(self : GckFile, t : GckSeqType) -> Unit {
  self.seq_type = t
}

///|
/// Set whether the sequence is circular.
pub fn GckFile::set_circular(self : GckFile, b : Bool) -> Unit {
  self.is_circular = b
}

///|
/// Set the sequence.
pub fn GckFile::set_sequence(self : GckFile, s : String) -> Unit {
  self.sequence = s
}

///|
/// Add a feature to the file.
pub fn GckFile::add_feature(self : GckFile, f : GckFeature) -> Unit {
  self.features.push(f)
}

// ============================================================================
// Binary helpers
// ============================================================================

///|
/// Convert an integer 0-15 to a lowercase hex digit character.
fn gck_hex_digit(n : Int) -> Char {
  if n < 10 {
    Char::from_int(48 + n)
  } else {
    Char::from_int(87 + n)
  }
}

///|
/// Convert a hex digit character's ASCII code to its integer value 0-15.
fn gck_hex_value(c : Int) -> Int {
  if c >= 48 && c <= 57 {
    c - 48
  } else if c >= 65 && c <= 70 {
    c - 55
  } else if c >= 97 && c <= 102 {
    c - 87
  } else {
    0
  }
}

///|
/// Read a uint16 big-endian from a byte array.
fn gck_read_u16_be(data : Array[Int], offset : Int) -> Int {
  (data[offset] & 0xFF) * 256 + (data[offset + 1] & 0xFF)
}

///|
/// Read a uint32 big-endian from a byte array.
fn gck_read_u32_be(data : Array[Int], offset : Int) -> Int {
  (data[offset] & 0xFF) * 16777216 +
  (data[offset + 1] & 0xFF) * 65536 +
  (data[offset + 2] & 0xFF) * 256 +
  (data[offset + 3] & 0xFF)
}

///|
/// Write a uint16 big-endian into a byte array at the given offset.
fn gck_write_u16_be(data : Array[Int], offset : Int, value : Int) -> Unit {
  data[offset] = (value >> 8) & 0xFF
  data[offset + 1] = value & 0xFF
}

///|
/// Write a uint32 big-endian into a byte array at the given offset.
fn gck_write_u32_be(data : Array[Int], offset : Int, value : Int) -> Unit {
  data[offset] = (value >> 24) & 0xFF
  data[offset + 1] = (value >> 16) & 0xFF
  data[offset + 2] = (value >> 8) & 0xFF
  data[offset + 3] = value & 0xFF
}

///|
/// Convert a hex string to a byte array.
fn gck_hex_to_bytes(hex : String) -> Array[Int] {
  let n = hex.length()
  let bytes = Array::make(n / 2, 0)
  for i in 0..<bytes.length() {
    let high = gck_hex_value(hex.unsafe_get(i * 2).to_int())
    let low = gck_hex_value(hex.unsafe_get(i * 2 + 1).to_int())
    bytes[i] = high * 16 + low
  }
  bytes
}

///|
/// Convert a byte array to a lowercase hex string.
fn gck_bytes_to_hex(data : Array[Int]) -> String {
  let sb = StringBuilder::new()
  for i in 0..<data.length() {
    let b = data[i] & 0xFF
    sb.write_char(gck_hex_digit((b >> 4) & 0xF))
    sb.write_char(gck_hex_digit(b & 0xF))
  }
  sb.to_string()
}

///|
/// Convert ASCII byte array slice to a string.
fn gck_ascii_to_string(data : Array[Int], offset : Int, length : Int) -> String {
  let sb = StringBuilder::new()
  for i in 0..<length {
    let b = data[offset + i] & 0xFF
    sb.write_char(Char::from_int(b))
  }
  sb.to_string()
}

///|
/// Convert a string to an ASCII byte array.
fn gck_string_to_ascii(s : String) -> Array[Int] {
  let n = s.length()
  let out = Array::make(n, 0)
  for i in 0..<n {
    out[i] = s.unsafe_get(i).to_int() & 0xFF
  }
  out
}

// ============================================================================
// Parsing
// ============================================================================

///|
/// Parse a GCK file from a hex string.
pub fn gck_parse(hex : String) -> GckFile {
  let bytes = gck_hex_to_bytes(hex)
  gck_from_bytes(bytes)
}

///|
/// Parse a GCK file from a byte array.
pub fn gck_from_bytes(data : Array[Int]) -> GckFile {
  let file = GckFile::new()
  let n = data.length()
  if n < 14 {
    return file
  }
  // Verify magic
  let magic = gck_ascii_to_string(data, 0, 8)
  if magic != "GCK_FILE" {
    return file
  }
  file.version = data[8] & 0xFF
  file.seq_type = GckSeqType::from_int(data[9] & 0xFF)
  file.is_circular = (data[10] & 0xFF) != 0
  // data[11] is reserved
  let feature_count = gck_read_u16_be(data, 12)
  let mut offset = 14
  // Parse features
  for _ in 0..<feature_count {
    if offset + 2 > n {
      break
    }
    let name_len = gck_read_u16_be(data, offset)
    offset = offset + 2
    if offset + name_len > n {
      break
    }
    let name = gck_ascii_to_string(data, offset, name_len)
    offset = offset + name_len
    if offset + 2 > n {
      break
    }
    let type_len = gck_read_u16_be(data, offset)
    offset = offset + 2
    if offset + type_len > n {
      break
    }
    let type_ = gck_ascii_to_string(data, offset, type_len)
    offset = offset + type_len
    if offset + 1 > n {
      break
    }
    let direction = data[offset] & 0xFF
    offset = offset + 1
    if offset + 2 > n {
      break
    }
    let seg_count = gck_read_u16_be(data, offset)
    offset = offset + 2
    let segs = Array::new()
    for _ in 0..<seg_count {
      if offset + 8 > n {
        break
      }
      let start = gck_read_u32_be(data, offset)
      let end = gck_read_u32_be(data, offset + 4)
      segs.push((start, end))
      offset = offset + 8
    }
    file.features.push(
      GckFeature::new(
        name=name,
        type_=type_,
        direction=direction,
        segments=segs,
      ),
    )
  }
  // Parse sequence
  if offset + 4 <= n {
    let seq_len = gck_read_u32_be(data, offset)
    offset = offset + 4
    if offset + seq_len <= n {
      file.sequence = gck_ascii_to_string(data, offset, seq_len)
    }
  }
  file
}

// ============================================================================
// Serialization
// ============================================================================

///|
/// Serialize a GckFile to a hex string.
pub fn gck_write(file : GckFile) -> String {
  let bytes = gck_to_bytes(file)
  gck_bytes_to_hex(bytes)
}

///|
/// Serialize a GckFile to a byte array.
pub fn gck_to_bytes(file : GckFile) -> Array[Int] {
  let out = Array::new()
  // Magic
  let magic = gck_string_to_ascii("GCK_FILE")
  for b in magic {
    out.push(b)
  }
  // Version, seq_type, topology, reserved
  out.push(file.version & 0xFF)
  out.push(file.seq_type.to_int() & 0xFF)
  out.push(if file.is_circular { 1 } else { 0 })
  out.push(0) // reserved
  // Feature count
  let fc = Array::make(2, 0)
  gck_write_u16_be(fc, 0, file.features.length())
  for b in fc {
    out.push(b)
  }
  // Features
  for f in file.features {
    // Name
    let name_bytes = gck_string_to_ascii(f.name)
    let nl = Array::make(2, 0)
    gck_write_u16_be(nl, 0, name_bytes.length())
    for b in nl {
      out.push(b)
    }
    for b in name_bytes {
      out.push(b)
    }
    // Type
    let type_bytes = gck_string_to_ascii(f.type_)
    let tl = Array::make(2, 0)
    gck_write_u16_be(tl, 0, type_bytes.length())
    for b in tl {
      out.push(b)
    }
    for b in type_bytes {
      out.push(b)
    }
    // Direction
    out.push(f.direction & 0xFF)
    // Segments
    let sc = Array::make(2, 0)
    gck_write_u16_be(sc, 0, f.segments.length())
    for b in sc {
      out.push(b)
    }
    for seg in f.segments {
      let (start, end) = seg
      let seg_bytes = Array::make(8, 0)
      gck_write_u32_be(seg_bytes, 0, start)
      gck_write_u32_be(seg_bytes, 4, end)
      for b in seg_bytes {
        out.push(b)
      }
    }
  }
  // Sequence length
  let seq_bytes = gck_string_to_ascii(file.sequence)
  let sl = Array::make(4, 0)
  gck_write_u32_be(sl, 0, seq_bytes.length())
  for b in sl {
    out.push(b)
  }
  // Sequence
  for b in seq_bytes {
    out.push(b)
  }
  out
}

// ============================================================================
// SeqRecord conversion
// ============================================================================

///|
/// Convert a GckFile to a SeqRecord.
pub fn gck_to_seqrecord(file : GckFile) -> SeqRecord {
  let id = "gck"
  let topo = if file.is_circular { " circular" } else { " linear" }
  let desc = file.seq_type.to_string() + topo
  SeqRecord::new(Seq::new(file.sequence), id=id, name=id, description=desc)
}

///|
/// Convert a SeqRecord to a GckFile (with default type DNA, linear).
pub fn gck_from_seqrecord(rec : SeqRecord) -> GckFile {
  let file = GckFile::new()
  file.seq_type = GckDna
  file.sequence = rec.seq.data
  file
}

// ============================================================================
// Utility functions
// ============================================================================

///|
/// Count features by type.
pub fn gck_count_features_by_type(file : GckFile, type_name : String) -> Int {
  let mut count = 0
  for f in file.features {
    if f.type_ == type_name {
      count = count + 1
    }
  }
  count
}

///|
/// Get the sequence length.
pub fn gck_seq_length(file : GckFile) -> Int {
  file.sequence.length()
}

///|
/// Generate a summary string.
pub fn gck_summary(file : GckFile) -> String {
  let sb = StringBuilder::new()
  sb.write_string("GCK File Summary:\n")
  sb.write_string("  Version: " + file.version.to_string() + "\n")
  sb.write_string("  Type: " + file.seq_type.to_string() + "\n")
  let topo = if file.is_circular { "circular" } else { "linear" }
  sb.write_string("  Topology: " + topo + "\n")
  sb.write_string("  Sequence Length: " + gck_seq_length(file).to_string() + " bp\n")
  sb.write_string("  Features: " + file.features.length().to_string() + "\n")
  if file.features.length() > 0 {
    sb.write_string("  Feature List:\n")
    for f in file.features {
      let dir_str = if f.direction == 1 {
        "forward"
      } else if f.direction == 2 {
        "reverse"
      } else {
        "none"
      }
      sb.write_string(
        "    " + f.name + " (" + f.type_ + ", " + dir_str + ") - " +
        f.segments.length().to_string() + " segment(s)\n",
      )
    }
  }
  sb.to_string()
}

///|
/// Create a sample GckFile for testing/demonstration.
pub fn gck_sample() -> GckFile {
  let file = GckFile::new()
  file.version = 1
  file.seq_type = GckDna
  file.is_circular = true
  file.sequence = "ATGGCTAGCAATGGGTCGCGATAGCGATAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAA"
  file.features.push(
    GckFeature::new(
      name="AmpR",
      type_="CDS",
      direction=1,
      segments=[(0, 19)],
    ),
  )
  file.features.push(
    GckFeature::new(
      name="lacZ",
      type_="CDS",
      direction=1,
      segments=[(30, 49)],
    ),
  )
  file.features.push(
    GckFeature::new(
      name="MCS",
      type_="misc",
      direction=0,
      segments=[(50, 63)],
    ),
  )
  file
}
