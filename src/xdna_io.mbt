// Bio.SeqIO.XdnaIO - Parser and writer for the XDNA binary format used by
// Geneious, mirroring BioPython's Bio.SeqIO.XdnaIO.
//
// XDNA is a binary format. All multi-byte integers are BIG-ENDIAN.
//
// Per-record layout (each record is self-contained):
//   Byte  0    : 0x00 (signature)
//   Byte  1    : 0x00 (signature)
//   Byte  2    : format version (usually 1)
//   Byte  3    : sequence type (0=DNA, 1=RNA, 2=Protein, 3=Unknown)
//   Bytes 4-7  : name length        (uint32 big-endian)
//   Bytes 8-11 : sequence length    (uint32 big-endian, in bases/chars)
//   Bytes 12-15: annotation length  (uint32 big-endian)
//   Bytes 16-19: reserved           (uint32 big-endian, 0)
//   Then       : name bytes         (ASCII, name_length bytes)
//   Then       : sequence bytes
//                - DNA: 2 bits per base (A=00, C=01, G=10, T=11),
//                       4 bases per byte, MSB first -> ceil(seq_len/4) bytes
//                - RNA / Protein / Unknown: raw ASCII bytes (1 byte per char)
//   Then       : annotation bytes   (ASCII, ann_length bytes)
//
// A complete file is one or more records concatenated, followed by a 4-byte
// big-endian checksum (sum of all preceding bytes mod 2^32).

// ============================================================================
// XdnaSeqType - the molecular type of an XDNA record
// ============================================================================

///|
/// The molecular type stored in an XDNA record's header byte 3.
pub enum XdnaSeqType {
  DnaType
  RnaType
  ProteinType
  UnknownType
} derive(Eq)

///|
/// Convert a sequence type to its integer code: DNA=0, RNA=1, Protein=2,
/// Unknown=3.
pub fn XdnaSeqType::to_int(self : XdnaSeqType) -> Int {
  match self {
    DnaType => 0
    RnaType => 1
    ProteinType => 2
    UnknownType => 3
  }
}

///|
/// Convert an integer code back to an XdnaSeqType. Values outside 0-3 map to
/// UnknownType.
pub fn XdnaSeqType::from_int(code : Int) -> XdnaSeqType {
  match code {
    0 => DnaType
    1 => RnaType
    2 => ProteinType
    _ => UnknownType
  }
}

///|
/// Short name of the type ("DnaType", "RnaType", ...).
pub fn XdnaSeqType::to_string(self : XdnaSeqType) -> String {
  match self {
    DnaType => "DnaType"
    RnaType => "RnaType"
    ProteinType => "ProteinType"
    UnknownType => "UnknownType"
  }
}

///|
/// Human-readable description of the type ("DNA", "RNA", "Protein", "Unknown").
pub fn XdnaSeqType::description(self : XdnaSeqType) -> String {
  match self {
    DnaType => "DNA"
    RnaType => "RNA"
    ProteinType => "Protein"
    UnknownType => "Unknown"
  }
}

// ============================================================================
// XdnaRecord - a single sequence record
// ============================================================================

///|
/// A single XDNA record: a named sequence of a known molecular type, with
/// optional annotation text and a checksum.
pub struct XdnaRecord {
  name : String
  sequence : String
  seq_type : XdnaSeqType
  annotations : String
  checksum : Int
}

///|
/// Construct an XdnaRecord. The checksum is computed from the record's bytes
/// when serialized; pass 0 here and it will be filled in by `xdna_to_bytes`.
pub fn XdnaRecord::new(
  name : String,
  sequence : String,
  seq_type : XdnaSeqType,
  annotations : String,
  checksum? : Int = 0,
) -> XdnaRecord {
  XdnaRecord::{ name, sequence, seq_type, annotations, checksum }
}

///|
/// Record name.
pub fn XdnaRecord::name(self : XdnaRecord) -> String {
  self.name
}

///|
/// Sequence string (decoded).
pub fn XdnaRecord::sequence(self : XdnaRecord) -> String {
  self.sequence
}

///|
/// Molecular type of the sequence.
pub fn XdnaRecord::seq_type(self : XdnaRecord) -> XdnaSeqType {
  self.seq_type
}

///|
/// Annotation text (raw).
pub fn XdnaRecord::annotations(self : XdnaRecord) -> String {
  self.annotations
}

///|
/// Stored checksum value.
pub fn XdnaRecord::checksum(self : XdnaRecord) -> Int {
  self.checksum
}

///|
/// Human-readable representation of this record.
pub fn XdnaRecord::to_string(self : XdnaRecord) -> String {
  let buf = StringBuilder::new()
  buf.write_string("XdnaRecord(name=")
  buf.write_string(self.name)
  buf.write_string(", seq_type=")
  buf.write_string(self.seq_type.to_string())
  buf.write_string(", length=")
  buf.write_string(self.sequence.length().to_string())
  buf.write_string(", annotations_length=")
  buf.write_string(self.annotations.length().to_string())
  buf.write_string(", checksum=")
  buf.write_string(self.checksum.to_string())
  buf.write_string(")")
  buf.to_string()
}

// ============================================================================
// XdnaFile - a complete XDNA document
// ============================================================================

///|
/// A complete XDNA file: one or more records and a format version.
pub struct XdnaFile {
  records : Array[XdnaRecord]
  mut version : Int
}

///|
/// Construct an empty XdnaFile with the default version (1).
pub fn XdnaFile::new() -> XdnaFile {
  XdnaFile::{ records: Array::new(), version: 1 }
}

///|
/// All records in the file.
pub fn XdnaFile::records(self : XdnaFile) -> Array[XdnaRecord] {
  self.records.copy()
}

///|
/// Format version (usually 1).
pub fn XdnaFile::version(self : XdnaFile) -> Int {
  self.version
}

///|
/// Number of records.
pub fn XdnaFile::n_records(self : XdnaFile) -> Int {
  self.records.length()
}

///|
/// Append a record to the file.
pub fn XdnaFile::add_record(self : XdnaFile, record : XdnaRecord) -> Unit {
  self.records.push(record)
}

// ============================================================================
// 2-bit DNA encoding (A=00, C=01, G=10, T=11)
// ============================================================================

///|
/// Convert a nucleotide character to its 2-bit code:
/// A -> 0, C -> 1, G -> 2, T -> 3. Returns -1 for any other character.
pub fn xdna_base_to_code(base : Char) -> Int {
  let c = base.to_int()
  if c == 'A'.to_int() {
    0
  } else if c == 'C'.to_int() {
    1
  } else if c == 'G'.to_int() {
    2
  } else if c == 'T'.to_int() {
    3
  } else {
    -1
  }
}

///|
/// Convert a 2-bit code back to a nucleotide character:
/// 0 -> A, 1 -> C, 2 -> G, 3 -> T. Any other code maps to 'N'.
pub fn xdna_code_to_base(code : Int) -> Char {
  if code == 0 {
    'A'
  } else if code == 1 {
    'C'
  } else if code == 2 {
    'G'
  } else if code == 3 {
    'T'
  } else {
    'N'
  }
}

///|
/// Pack a DNA string into 2-bit bytes (4 bases per byte, MSB first).
/// Unknown bases are packed as A (code 0). Returns ceil(length/4) bytes.
pub fn xdna_pack_sequence(sequence : String) -> Array[Int] {
  let n = sequence.length()
  let nbytes = (n + 3) / 4
  let packed = Array::make(nbytes, 0)
  for i in 0..<n {
    let code = xdna_base_to_code(sequence.unsafe_get(i).unsafe_to_char())
    let safe_code = if code < 0 { 0 } else { code }
    let byte_idx = i / 4
    let shift = 6 - (i % 4) * 2
    packed[byte_idx] = packed[byte_idx] | (safe_code << shift)
  }
  packed
}

///|
/// Unpack 2-bit bytes back into a DNA string of length `seq_len`.
pub fn xdna_unpack_sequence(packed : Array[Int], seq_len : Int) -> String {
  let sb = StringBuilder::new()
  for i in 0..<seq_len {
    let byte_idx = i / 4
    let shift = 6 - (i % 4) * 2
    let code = if byte_idx < packed.length() {
      (packed[byte_idx] >> shift) & 3
    } else {
      0
    }
    sb.write_char(xdna_code_to_base(code))
  }
  sb.to_string()
}

// ============================================================================
// Big-endian serialization helpers
// ============================================================================

///|
/// Read a 32-bit big-endian unsigned integer from `data` at `offset`.
pub fn xdna_read_u32_be(data : Array[Int], offset : Int) -> Int {
  (data[offset] << 24) |
  (data[offset + 1] << 16) |
  (data[offset + 2] << 8) |
  data[offset + 3]
}

///|
/// Write a 32-bit unsigned integer as 4 big-endian bytes.
pub fn xdna_write_u32_be(value : Int) -> Array[Int] {
  [
    (value >> 24) & 0xFF,
    (value >> 16) & 0xFF,
    (value >> 8) & 0xFF,
    value & 0xFF,
  ]
}

// ============================================================================
// Serialization (to / from byte arrays)
// ============================================================================

///|
/// Serialize an XdnaFile to a byte array. Each record is written with its
/// 20-byte header followed by name, sequence (2-bit packed for DNA, raw ASCII
/// otherwise), and annotation bytes. A 4-byte big-endian checksum
/// (sum of all preceding bytes) is appended at the end.
pub fn xdna_to_bytes(file : XdnaFile) -> Array[Int] {
  let bytes : Array[Int] = Array::new()
  for ri in 0..<file.records.length() {
    let rec = file.records[ri]
    // 20-byte record header.
    bytes.push(0)
    bytes.push(0)
    bytes.push(file.version & 0xFF)
    bytes.push(rec.seq_type.to_int() & 0xFF)
    let name_bytes = string_to_ascii_bytes(rec.name)
    let ann_bytes = string_to_ascii_bytes(rec.annotations)
    let seq_bytes : Array[Int] = match rec.seq_type {
      DnaType => xdna_pack_sequence(rec.sequence)
      _ => string_to_ascii_bytes(rec.sequence)
    }
    let name_len = name_bytes.length()
    let seq_len = rec.sequence.length()
    let ann_len = ann_bytes.length()
    for b in xdna_write_u32_be(name_len) {
      bytes.push(b)
    }
    for b in xdna_write_u32_be(seq_len) {
      bytes.push(b)
    }
    for b in xdna_write_u32_be(ann_len) {
      bytes.push(b)
    }
    for b in xdna_write_u32_be(0) {
      bytes.push(b)
    }
    for b in name_bytes {
      bytes.push(b)
    }
    for b in seq_bytes {
      bytes.push(b)
    }
    for b in ann_bytes {
      bytes.push(b)
    }
  }
  // Checksum: sum of all bytes mod 2^32, as big-endian 4 bytes.
  let mut sum = 0
  for i in 0..<bytes.length() {
    sum = (sum + bytes[i]) & 0xFFFFFFFF
  }
  for b in xdna_write_u32_be(sum) {
    bytes.push(b)
  }
  bytes
}

///|
/// Deserialize an XdnaFile from a byte array produced by `xdna_to_bytes`.
/// Trailing 4 bytes are read as the checksum (and stored on each record).
pub fn xdna_from_bytes(data : Array[Int]) -> XdnaFile {
  let file = XdnaFile::new()
  let n = data.length()
  let mut offset = 0
  // Determine checksum (last 4 bytes) if present.
  let stored_checksum = if n >= 4 {
    xdna_read_u32_be(data, n - 4)
  } else {
    0
  }
  // Parse records while at least one record header (20 bytes) fits before the
  // trailing checksum.
  while offset + 20 <= n - 4 || (offset + 20 <= n && n < 4) {
    // Header
    let _sig0 = data[offset]
    let _sig1 = data[offset + 1]
    let version = data[offset + 2] & 0xFF
    let type_code = data[offset + 3] & 0xFF
    let name_len = xdna_read_u32_be(data, offset + 4)
    let seq_len = xdna_read_u32_be(data, offset + 8)
    let ann_len = xdna_read_u32_be(data, offset + 12)
    // reserved at offset+16..offset+19
    offset = offset + 20
    // Name
    let name = ascii_bytes_to_string(data, offset, name_len)
    offset = offset + name_len
    // Sequence
    let seq_type = XdnaSeqType::from_int(type_code)
    let sequence : String = match seq_type {
      DnaType => {
        let nbytes = (seq_len + 3) / 4
        let packed = Array::make(nbytes, 0)
        for i in 0..<nbytes {
          if offset + i < n {
            packed[i] = data[offset + i]
          }
        }
        xdna_unpack_sequence(packed, seq_len)
      }
      _ => ascii_bytes_to_string(data, offset, seq_len)
    }
    let seq_byte_count : Int = match seq_type {
      DnaType => (seq_len + 3) / 4
      _ => seq_len
    }
    offset = offset + seq_byte_count
    // Annotations
    let annotations = ascii_bytes_to_string(data, offset, ann_len)
    offset = offset + ann_len
    // Set version from the first record.
    if file.records.length() == 0 {
      file.version = version
    }
    file.records.push(
      XdnaRecord::{ name, sequence, seq_type, annotations, checksum: stored_checksum },
    )
    // Stop if we have consumed everything up to the checksum.
    if offset >= n - 4 && n >= 4 {
      break
    }
    if offset >= n {
      break
    }
  }
  file
}

// ============================================================================
// I/O functions (hex string)
// ============================================================================

///|
/// Serialize an XdnaFile to a lowercase hex string (two hex digits per byte),
/// suitable for text storage.
pub fn xdna_write(file : XdnaFile) -> String {
  let bytes = xdna_to_bytes(file)
  let sb = StringBuilder::new()
  for i in 0..<bytes.length() {
    let b = bytes[i]
    sb.write_char(xdna_hex_digit((b >> 4) & 0xF))
    sb.write_char(xdna_hex_digit(b & 0xF))
  }
  sb.to_string()
}

///|
/// Parse a hex string produced by `xdna_write` back into an XdnaFile.
pub fn xdna_read(data : String) -> XdnaFile {
  let bytes : Array[Int] = Array::new()
  let len = data.length()
  let mut i = 0
  while i + 1 < len {
    let high = xdna_hex_value(data[i].to_int())
    let low = xdna_hex_value(data[i + 1].to_int())
    bytes.push(high * 16 + low)
    i = i + 2
  }
  xdna_from_bytes(bytes)
}

///|
/// Return the sequence string of a record (convenience accessor).
pub fn xdna_get_sequence_string(record : XdnaRecord) -> String {
  record.sequence
}

// ============================================================================
// SeqRecord integration
// ============================================================================

///|
/// Convert all XdnaRecords in a file to SeqRecord objects. The record name
/// becomes the id, name, and description.
pub fn xdna_to_seq_records(file : XdnaFile) -> Array[SeqRecord] {
  let records = Array::new()
  for i in 0..<file.records.length() {
    let rec = file.records[i]
    records.push(
      SeqRecord::new(
        Seq::new(rec.sequence),
        id=rec.name,
        name=rec.name,
        description=rec.name,
      ),
    )
  }
  records
}

///|
/// Build an XdnaFile from an array of SeqRecord objects, tagging every
/// record with the given sequence type.
pub fn xdna_from_seq_records(
  records : Array[SeqRecord],
  seq_type : XdnaSeqType,
) -> XdnaFile {
  let file = XdnaFile::new()
  for i in 0..<records.length() {
    let rec = records[i]
    file.add_record(
      XdnaRecord::new(rec.id, rec.seq.to_string(), seq_type, "", checksum=0),
    )
  }
  file
}

// ============================================================================
// Sample data
// ============================================================================

///|
/// Return a sample XdnaFile with 2 DNA records: "seq1" (ATGCATGC) and
/// "seq2" (GATTACA).
pub fn xdna_sample_file() -> XdnaFile {
  let file = XdnaFile::new()
  file.add_record(XdnaRecord::new("seq1", "ATGCATGC", DnaType, ""))
  file.add_record(XdnaRecord::new("seq2", "GATTACA", DnaType, ""))
  file
}

// ============================================================================
// Private helpers
// ============================================================================

///|
/// Convert a String to an array of ASCII byte values (one byte per char,
/// taking the low 8 bits of each code unit).
fn string_to_ascii_bytes(s : String) -> Array[Int] {
  let out = Array::new()
  let n = s.length()
  for i in 0..<n {
    out.push(s.unsafe_get(i).to_int() & 0xFF)
  }
  out
}

///|
/// Read `len` bytes from `data` starting at `offset` and return them as an
/// ASCII string.
fn ascii_bytes_to_string(data : Array[Int], offset : Int, len : Int) -> String {
  let sb = StringBuilder::new()
  for i in 0..<len {
    if offset + i < data.length() {
      sb.write_char(Int::unsafe_to_char(data[offset + i] & 0xFF))
    }
  }
  sb.to_string()
}

///|
/// Convert an integer 0-15 to a lowercase hex digit character.
fn xdna_hex_digit(n : Int) -> Char {
  if n < 10 {
    Int::unsafe_to_char('0'.to_int() + n)
  } else {
    Int::unsafe_to_char('a'.to_int() + (n - 10))
  }
}

///|
/// Convert a hex digit character's ASCII code to its integer value 0-15.
fn xdna_hex_value(c : Int) -> Int {
  let code_0 = '0'.to_int()
  let code_9 = '9'.to_int()
  let code_a = 'a'.to_int()
  let code_f = 'f'.to_int()
  let code_A = 'A'.to_int()
  let code_F = 'F'.to_int()
  if c >= code_0 && c <= code_9 {
    c - code_0
  } else if c >= code_a && c <= code_f {
    c - code_a + 10
  } else if c >= code_A && c <= code_F {
    c - code_A + 10
  } else {
    0
  }
}
