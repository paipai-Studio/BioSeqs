///|
/// BAM format parsing module.
/// BAM is the binary version of SAM format.

///|
/// BAM parsing errors.
pub suberror BamError {
  BamError(String)
}

///|
/// BAM magic number (BAM\1).
let bam_magic : Array[Int] = [66, 65, 77, 1]

///|
/// BGZF magic number (1f 8b 08 04).
let bgzf_magic : Array[Int] = [0x1f, 0x8b, 0x08, 0x04]

///|
/// Check if data starts with BAM magic number.
pub fn is_bam_magic(data : String) -> Bool {
  if data.length() < 4 {
    return false
  }
  for i = 0; i < 4; i = i + 1 {
    if data.unsafe_get(i).to_int() != bam_magic[i] {
      return false
    }
  }
  true
}

///|
/// Check if data starts with BGZF magic number.
pub fn is_bgzf(data : String) -> Bool {
  if data.length() < 4 {
    return false
  }
  for i = 0; i < 4; i = i + 1 {
    if data.unsafe_get(i).to_int() != bgzf_magic[i] {
      return false
    }
  }
  true
}

///|
/// Reference sequence entry from BAM header.
pub struct BamReference {
  name : String
  length : Int
} derive(Eq)

///|
/// BAM header structure.
pub struct BamHeader {
  text : String
  references : FixedArray[BamReference]
}

///|
/// CIGAR operation types for BAM (same as SAM).
pub enum BamCigarOp {
  Match
  Insert
  Delete
  Skip
  SoftClip
  HardClip
  Pad
  Equal
  Diff
} derive(Eq)

///|
/// CIGAR element: operation + length.
pub struct BamCigarElement {
  op : BamCigarOp
  length : Int
} derive(Eq)

///|
/// BAM alignment record.
pub struct BamRecord {
  ref_id : Int
  pos : Int
  bin : Int
  mapq : Int
  flag : Int
  l_read_name : Int
  cigar_len : Int
  l_seq : Int
  next_ref_id : Int
  next_pos : Int
  tlen : Int
  qname : String
  cigar : FixedArray[BamCigarElement]
  seq : String
  qual : String
}

///|
/// BAM file representation.
pub struct BamFile {
  header : BamHeader
  records : Array[BamRecord]
}

///|
/// Convert BAM CIGAR operation code to enum.
fn cigar_op_from_code(code : Int) -> BamCigarOp {
  match code {
    0 => BamCigarOp::Match
    1 => BamCigarOp::Insert
    2 => BamCigarOp::Delete
    3 => BamCigarOp::Skip
    4 => BamCigarOp::SoftClip
    5 => BamCigarOp::HardClip
    6 => BamCigarOp::Pad
    7 => BamCigarOp::Equal
    8 => BamCigarOp::Diff
    _ => BamCigarOp::Match
  }
}

///|
/// Read an Int from data at position, advance position.
fn read_int8(data : String, pos : Int) -> (Int, Int) {
  let val = data.unsafe_get(pos).to_int() & 0xFF
  (val, pos + 1)
}

///|
/// Read a UInt16 (little-endian) from data at position, advance position.
fn read_uint16_le(data : String, pos : Int) -> (Int, Int) {
  let b0 = data.unsafe_get(pos).to_int() & 0xFF
  let b1 = data.unsafe_get(pos + 1).to_int() & 0xFF
  let val = ((b1 << 8) | b0) & 0xFFFF
  (val, pos + 2)
}

///|
/// Read an Int32 (little-endian) from data at position, advance position.
fn read_int32_le(data : String, pos : Int) -> (Int, Int) {
  let b0 = data.unsafe_get(pos).to_int() & 0xFF
  let b1 = data.unsafe_get(pos + 1).to_int() & 0xFF
  let b2 = data.unsafe_get(pos + 2).to_int() & 0xFF
  let b3 = data.unsafe_get(pos + 3).to_int() & 0xFF
  let val = b0 | (b1 << 8) | (b2 << 16) | (b3 << 24)
  (val, pos + 4)
}

///|
/// Read a UInt32 (little-endian) from data at position, advance position.
fn read_uint32_le(data : String, pos : Int) -> (Int, Int) {
  let b0 = data.unsafe_get(pos).to_int() & 0xFF
  let b1 = data.unsafe_get(pos + 1).to_int() & 0xFF
  let b2 = data.unsafe_get(pos + 2).to_int() & 0xFF
  let b3 = data.unsafe_get(pos + 3).to_int() & 0xFF
  let val = b0 | (b1 << 8) | (b2 << 16) | (b3 << 24)
  (val & 0xFFFFFFFF, pos + 4)
}

///|
/// Read a null-terminated string from data at position, advance position.
fn read_cstring(data : String, pos : Int) -> (String, Int) {
  let mut end = pos
  while end < data.length() && data.unsafe_get(end) != 0 {
    end = end + 1
  }
  let s = data[pos:end].to_owned()
  (s, end + 1)
}

///|
/// Decode BAM sequence from 4-bit encoding.
/// Each byte encodes 2 nucleotides: high nibble first.
fn decode_bam_seq(encoded : String, length : Int) -> String {
  let buf : FixedArray[UInt16] = FixedArray::make(length, 0)
  let mut i = 0
  while i < length {
    let byte_pos = i / 2
    let byte = encoded.unsafe_get(byte_pos).to_int()
    let offset = (i % 2) * 4
    let code = (byte >> offset) & 0xF
    let c = match code {
      0 => 'A'.to_int().to_uint16()
      1 => 'C'.to_int().to_uint16()
      2 => 'G'.to_int().to_uint16()
      4 => 'T'.to_int().to_uint16()
      8 => 'N'.to_int().to_uint16()
      _ => 'N'.to_int().to_uint16()
    }
    buf[i] = c
    i = i + 1
  }
  unsafe_from_uint16_fixedarray(buf)
}

///|
/// Decode BAM quality from Phred+33 encoding.
fn decode_bam_qual(encoded : String, length : Int) -> String {
  let buf : FixedArray[UInt16] = FixedArray::make(length, 0)
  let mut i = 0
  while i < length {
    let q = (encoded.unsafe_get(i).to_int() + 33).to_uint16()
    buf[i] = q
    i = i + 1
  }
  unsafe_from_uint16_fixedarray(buf)
}

///|
/// Parse the BAM header section.
fn parse_bam_header(data : String, pos : Int) -> (BamHeader, Int) raise BamError {
  let (l_text, pos1) = read_int32_le(data, pos)
  
  if pos1 + l_text > data.length() {
    raise BamError("Truncated BAM header text")
  }
  let text = data[pos1:pos1 + l_text].to_owned()
  let pos2 = pos1 + l_text
  
  let (n_ref, pos3) = read_int32_le(data, pos2)
  
  let references : FixedArray[BamReference] = FixedArray::make(n_ref, BamReference::{ name: "", length: 0 })
  let mut ref_idx = 0
  let mut current_pos = pos3
  while ref_idx < n_ref {
    let (l_name, pos4) = read_int32_le(data, current_pos)
    
    if pos4 + l_name > data.length() {
      raise BamError("Truncated BAM reference name")
    }
    let name = data[pos4:pos4 + l_name - 1].to_owned()
    let pos5 = pos4 + l_name
    
    let (length, pos6) = read_int32_le(data, pos5)
    
    references[ref_idx] = BamReference::{ name, length }
    ref_idx = ref_idx + 1
    current_pos = pos6
  }
  
  (BamHeader::{ text, references }, current_pos)
}

///|
/// Parse a single BAM alignment record.
fn parse_bam_record(data : String, pos : Int) -> (BamRecord, Int) raise BamError {
  let (ref_id, pos1) = read_int32_le(data, pos)
  let (pos_val, pos2) = read_int32_le(data, pos1)
  let (bin, pos3) = read_uint32_le(data, pos2)
  let (mapq, pos4) = read_int8(data, pos3)
  let (flag, pos5) = read_uint16_le(data, pos4)
  let (l_read_name, pos6) = read_uint16_le(data, pos5)
  let (cigar_len, pos7) = read_uint16_le(data, pos6)
  let (seq_len, pos8) = read_uint32_le(data, pos7)
  let (next_ref_id, pos9) = read_int32_le(data, pos8)
  let (next_pos, pos10) = read_int32_le(data, pos9)
  let (tlen, pos11) = read_int32_le(data, pos10)
  
  let (qname, pos12) = read_cstring(data, pos11)
  
  let cigar : FixedArray[BamCigarElement] = FixedArray::make(cigar_len, BamCigarElement::{ op: BamCigarOp::Match, length: 0 })
  let mut cigar_idx = 0
  let mut current_pos = pos12
  while cigar_idx < cigar_len {
    let (cigar_op_val, pos_new) = read_uint32_le(data, current_pos)
    let length = cigar_op_val >> 4
    let op_code = cigar_op_val & 0xF
    cigar[cigar_idx] = BamCigarElement::{
      op: cigar_op_from_code(op_code),
      length,
    }
    cigar_idx = cigar_idx + 1
    current_pos = pos_new
  }
  
  let seq_bytes = (seq_len + 1) / 2
  if current_pos + seq_bytes > data.length() {
    raise BamError("Truncated BAM sequence")
  }
  let encoded_seq = data[current_pos:current_pos + seq_bytes].to_owned()
  let pos13 = current_pos + seq_bytes
  
  let seq = decode_bam_seq(encoded_seq, seq_len)
  
  let qual_bytes = seq_len
  if pos13 + qual_bytes > data.length() {
    raise BamError("Truncated BAM quality")
  }
  let encoded_qual = data[pos13:pos13 + qual_bytes].to_owned()
  let pos14 = pos13 + qual_bytes
  
  let qual = decode_bam_qual(encoded_qual, seq_len)
  
  (BamRecord::{
    ref_id,
    pos: pos_val,
    bin,
    mapq,
    flag,
    l_read_name,
    cigar_len,
    l_seq: seq_len,
    next_ref_id,
    next_pos,
    tlen,
    qname,
    cigar,
    seq,
    qual,
  }, pos14)
}

///|
/// Parse raw BAM data (uncompressed).
/// The input should be the raw BAM binary content after BGZF decompression.
pub fn parse_bam(data : String) -> BamFile raise BamError {
  let n = data.length()
  if n < 4 {
    raise BamError("BAM file too short")
  }
  
  for i = 0; i < 4; i = i + 1 {
    if data.unsafe_get(i).to_int() != bam_magic[i] {
      raise BamError("Invalid BAM magic number")
    }
  }
  
  let (header, pos) = parse_bam_header(data, 4)
  
  let records = Array::new()
  let mut current_pos = pos
  while current_pos < n {
    let (block_size, pos_new) = read_int32_le(data, current_pos)
    if pos_new + block_size > n {
      break
    }
    
    let (record, _) = parse_bam_record(data, pos_new)
    
    records.push(record)
    current_pos = pos_new + block_size + 4
  }
  
  BamFile::{ header, records }
}

///|
/// Get the reference name for a record.
pub fn BamRecord::get_reference_name(self : BamRecord, header : BamHeader) -> String {
  if self.ref_id >= 0 && self.ref_id < header.references.length() {
    header.references[self.ref_id].name
  } else {
    "*"
  }
}

///|
/// Check if the read is paired.
pub fn BamRecord::is_paired(self : BamRecord) -> Bool {
  (self.flag & 0x0001) != 0
}

///|
/// Check if the read is properly paired.
pub fn BamRecord::is_proper_pair(self : BamRecord) -> Bool {
  (self.flag & 0x0002) != 0
}

///|
/// Check if the read is unmapped.
pub fn BamRecord::is_unmapped(self : BamRecord) -> Bool {
  (self.flag & 0x0004) != 0
}

///|
/// Check if the mate is unmapped.
pub fn BamRecord::is_mate_unmapped(self : BamRecord) -> Bool {
  (self.flag & 0x0008) != 0
}

///|
/// Check if this is the reverse complement.
pub fn BamRecord::is_reverse(self : BamRecord) -> Bool {
  (self.flag & 0x010) != 0
}

///|
/// Check if mate is reverse complement.
pub fn BamRecord::is_mate_reverse(self : BamRecord) -> Bool {
  (self.flag & 0x020) != 0
}

///|
/// Check if this is the first read in pair.
pub fn BamRecord::is_read1(self : BamRecord) -> Bool {
  (self.flag & 0x040) != 0
}

///|
/// Check if this is the second read in pair.
pub fn BamRecord::is_read2(self : BamRecord) -> Bool {
  (self.flag & 0x080) != 0
}

///|
/// Check if this is a secondary alignment.
pub fn BamRecord::is_secondary(self : BamRecord) -> Bool {
  (self.flag & 0x100) != 0
}

///|
/// Check if the read is filtered out.
pub fn BamRecord::is_supplementary(self : BamRecord) -> Bool {
  (self.flag & 0x800) != 0
}

///|
/// Get number of records in BAM file.
pub fn BamFile::len(self : BamFile) -> Int {
  self.records.length()
}

///|
/// Get number of reference sequences.
pub fn BamFile::num_references(self : BamFile) -> Int {
  self.header.references.length()
}

///|
/// Get all reference sequence names.
pub fn BamFile::get_reference_names(self : BamFile) -> Array[String] {
  let result = Array::new()
  for ref_ in self.header.references {
    result.push(ref_.name)
  }
  result
}

///|
/// Get reference length by name.
pub fn BamFile::get_reference_length(self : BamFile, name : String) -> Int? {
  for ref_ in self.header.references {
    if ref_.name == name {
      return Some(ref_.length)
    }
  }
  None
}

///|
/// Parse BAM from base64-encoded raw (uncompressed) BAM data.
/// This is useful for testing and cases where the data is already decompressed.
pub fn parse_bam_base64(base64_str : String) -> BamFile raise BamError {
  let raw_data = base64_decode(base64_str)
  parse_bam(raw_data)
}

///|
/// Parse BAM from base64-encoded raw (uncompressed) BAM data.
/// Returns None if parsing fails (for use in main functions).
pub fn try_parse_bam_base64(base64_str : String) -> BamFile? {
  let result = parse_bam_base64(base64_str) catch { _ => return None }
  Some(result)
}

///|
/// Base64 decoding table.
let base64_table : FixedArray[Int] = {
  let t = FixedArray::make(256, -1)
  let chars = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/"
  for i = 0; i < chars.length(); i = i + 1 {
    t[chars.unsafe_get(i).to_int()] = i
  }
  t
}

///|
/// Decode base64 string to raw binary string.
pub fn base64_decode(s : String) -> String {
  let n = s.length()
  let eq_char = '='.to_int().to_uint16()
  let padding = if n >= 2 && s.unsafe_get(n - 1) == eq_char {
    if n >= 3 && s.unsafe_get(n - 2) == eq_char { 2 } else { 1 }
  } else {
    0
  }
  let out_len = (n * 3) / 4 - padding
  let buf : FixedArray[UInt16] = FixedArray::make(out_len, 0)
  let mut i = 0
  let mut j = 0
  while i < n {
    let c0 = base64_table[s.unsafe_get(i).to_int()]
    let c1 = if i + 1 < n { base64_table[s.unsafe_get(i + 1).to_int()] } else { 0 }
    let c2 = if i + 2 < n { base64_table[s.unsafe_get(i + 2).to_int()] } else { 0 }
    let c3 = if i + 3 < n { base64_table[s.unsafe_get(i + 3).to_int()] } else { 0 }
    
    buf[j] = ((c0 << 2) | ((c1 >> 4) & 0x03)).to_uint16()
    if j + 1 < out_len {
      buf[j + 1] = (((c1 & 0x0F) << 4) | ((c2 >> 2) & 0x0F)).to_uint16()
    }
    if j + 2 < out_len {
      buf[j + 2] = (((c2 & 0x03) << 6) | (c3 & 0x3F)).to_uint16()
    }
    
    i = i + 4
    j = j + 3
  }
  unsafe_from_uint16_fixedarray(buf)
}

///|
/// Get the end position of the alignment (similar to SAM).
pub fn BamRecord::get_end(self : BamRecord) -> Int {
  if self.pos == 0 || self.cigar.length() == 0 {
    return 0
  }
  
  let mut end = self.pos - 1
  let mut i = 0
  while i < self.cigar.length() {
    let element = self.cigar[i]
    let op = element.op
    if op == BamCigarOp::Match ||
      op == BamCigarOp::Delete ||
      op == BamCigarOp::Skip ||
      op == BamCigarOp::Equal ||
      op == BamCigarOp::Diff {
      end = end + element.length
    }
    i = i + 1
  }
  end
}

///|
/// Calculate the alignment length (similar to SAM).
pub fn BamRecord::get_alignment_length(self : BamRecord) -> Int {
  let mut length = 0
  let mut i = 0
  while i < self.cigar.length() {
    let element = self.cigar[i]
    let op = element.op
    if op == BamCigarOp::Match ||
      op == BamCigarOp::Insert ||
      op == BamCigarOp::Equal ||
      op == BamCigarOp::Diff {
      length = length + element.length
    }
    i = i + 1
  }
  length
}

///|
/// Convert BAM record to SAM-like string representation.
pub fn BamRecord::to_sam_string(self : BamRecord, header : BamHeader) -> String {
  let rname = self.get_reference_name(header)
  let cigar_str = StringBuilder::new()
  for elem in self.cigar {
    let op_char = match elem.op {
      BamCigarOp::Match => 'M'
      BamCigarOp::Insert => 'I'
      BamCigarOp::Delete => 'D'
      BamCigarOp::Skip => 'N'
      BamCigarOp::SoftClip => 'S'
      BamCigarOp::HardClip => 'H'
      BamCigarOp::Pad => 'P'
      BamCigarOp::Equal => '='
      BamCigarOp::Diff => 'X'
    }
    cigar_str.write_string(elem.length.to_string())
    cigar_str.write_char(op_char)
  }
  let rnext = if self.next_ref_id >= 0 && self.next_ref_id < header.references.length() {
    header.references[self.next_ref_id].name
  } else {
    "*"
  }
  "\{self.qname}\t\{self.flag}\t\{rname}\t\{self.pos}\t\{self.mapq}\t\{cigar_str.to_string()}\t\{rnext}\t\{self.next_pos}\t\{self.tlen}\t\{self.seq}\t\{self.qual}"
}

///|
/// Convert BAM file to SAM format string.
pub fn BamFile::to_sam(self : BamFile) -> String {
  let sb = StringBuilder::new()
  sb.write_string(self.header.text)
  for record in self.records {
    sb.write_string(record.to_sam_string(self.header))
    sb.write_char('\n')
  }
  sb.to_string()
}