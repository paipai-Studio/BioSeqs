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
  data.unsafe_get(0).to_int() == 0x1f &&
  data.unsafe_get(1).to_int() == 0x8b &&
  data.unsafe_get(2).to_int() == 0x08 &&
  data.unsafe_get(3).to_int() == 0x04
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
/// A typed value of a BAM auxiliary tag.
pub enum BamTagValue {
  Char(Char)
  Int(Int)
  Float(Double)
  Str(String)
  IntArray(Array[Int])
  FloatArray(Array[Double])
}

///|
/// A BAM auxiliary tag (optional field after the quality string).
pub struct BamTag {
  name : String
  type_code : Char
  value : BamTagValue
}

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
  tags : Array[BamTag]
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
/// Read a signed int8.
fn read_int8_signed(data : String, pos : Int) -> (Int, Int) {
  let b = data.unsafe_get(pos).to_int() & 0xFF
  let v = if b >= 128 { b - 256 } else { b }
  (v, pos + 1)
}

///|
/// Read a signed int16 (little-endian).
fn read_int16_le(data : String, pos : Int) -> (Int, Int) {
  let (u, np) = read_uint16_le(data, pos)
  let v = if u >= 32768 { u - 65536 } else { u }
  (v, np)
}

///|
/// Read an IEEE 754 float32 (4 little-endian bytes) as a Double.
fn read_float32_le(data : String, pos : Int) -> (Double, Int) {
  let b0 = data.unsafe_get(pos).to_int() & 0xFF
  let b1 = data.unsafe_get(pos + 1).to_int() & 0xFF
  let b2 = data.unsafe_get(pos + 2).to_int() & 0xFF
  let b3 = data.unsafe_get(pos + 3).to_int() & 0xFF
  let sign = if b3 >= 128 { -1.0 } else { 1.0 }
  let exponent = (b3 & 0x7F) * 2 + (b2 >> 7)
  let mantissa = (b2 & 0x7F) * 65536 + b1 * 256 + b0
  let v = if exponent == 255 {
    0.0 / 0.0
  } else if exponent == 0 {
    sign * mantissa.to_double() * @math.pow(2.0, -149.0)
  } else {
    sign * (1.0 + mantissa.to_double() / 8388608.0) * @math.pow(2.0, (exponent - 127).to_double())
  }
  (v, pos + 4)
}

///|
/// Parse BAM auxiliary tags from `pos` until `end`.
fn parse_bam_tags(
  data : String,
  pos : Int,
  end : Int,
) -> (Array[BamTag], Int) raise BamError {
  let tags : Array[BamTag] = Array::new()
  let mut p = pos
  while p < end {
    if p + 3 > end {
      raise BamError("Truncated BAM tag")
    }
    let name = data[p:p + 2].to_owned()
    let type_code = data.unsafe_get(p + 2).unsafe_to_char()
    p = p + 3
    if type_code == 'A' {
      if p + 1 > end {
        raise BamError("Truncated A tag")
      }
      tags.push(BamTag::{ name, type_code, value: BamTagValue::Char(data.unsafe_get(p).unsafe_to_char()) })
      p = p + 1
    } else if type_code == 'c' {
      let (v, np) = read_int8_signed(data, p)
      tags.push(BamTag::{ name, type_code, value: BamTagValue::Int(v) })
      p = np
    } else if type_code == 'C' {
      let (v, np) = read_int8(data, p)
      tags.push(BamTag::{ name, type_code, value: BamTagValue::Int(v) })
      p = np
    } else if type_code == 's' {
      let (v, np) = read_int16_le(data, p)
      tags.push(BamTag::{ name, type_code, value: BamTagValue::Int(v) })
      p = np
    } else if type_code == 'S' {
      let (v, np) = read_uint16_le(data, p)
      tags.push(BamTag::{ name, type_code, value: BamTagValue::Int(v) })
      p = np
    } else if type_code == 'i' {
      let (v, np) = read_int32_le(data, p)
      tags.push(BamTag::{ name, type_code, value: BamTagValue::Int(v) })
      p = np
    } else if type_code == 'I' {
      let (v, np) = read_uint32_le(data, p)
      tags.push(BamTag::{ name, type_code, value: BamTagValue::Int(v) })
      p = np
    } else if type_code == 'f' {
      let (v, np) = read_float32_le(data, p)
      tags.push(BamTag::{ name, type_code, value: BamTagValue::Float(v) })
      p = np
    } else if type_code == 'Z' || type_code == 'H' {
      let (v, np) = read_cstring(data, p)
      tags.push(BamTag::{ name, type_code, value: BamTagValue::Str(v) })
      p = np
    } else if type_code == 'B' {
      if p + 1 > end {
        raise BamError("Truncated B tag")
      }
      let subtype = data.unsafe_get(p).unsafe_to_char()
      p = p + 1
      let (count, np) = read_int32_le(data, p)
      p = np
      if subtype == 'f' {
        let arr : Array[Double] = Array::new()
        for _j = 0; _j < count; _j = _j + 1 {
          let (v, np2) = read_float32_le(data, p)
          arr.push(v)
          p = np2
        }
        tags.push(BamTag::{ name, type_code, value: BamTagValue::FloatArray(arr) })
      } else {
        let arr : Array[Int] = Array::new()
        for _j = 0; _j < count; _j = _j + 1 {
          let pair = if subtype == 'c' {
            read_int8_signed(data, p)
          } else if subtype == 'C' {
            read_int8(data, p)
          } else if subtype == 's' {
            read_int16_le(data, p)
          } else if subtype == 'S' {
            read_uint16_le(data, p)
          } else if subtype == 'i' {
            read_int32_le(data, p)
          } else if subtype == 'I' {
            read_uint32_le(data, p)
          } else {
            raise BamError("Unknown BAM array subtype")
          }
          arr.push(pair.0)
          p = pair.1
        }
        tags.push(BamTag::{ name, type_code, value: BamTagValue::IntArray(arr) })
      }
    } else {
      raise BamError("Unknown BAM tag type")
    }
  }
  (tags, p)
}

///|
/// Decode BAM sequence from 4-bit encoding.
/// Each byte encodes 2 nucleotides: high nibble first.
fn decode_bam_seq(encoded : String, length : Int) -> String {
  let lookup = "=ACMGRSVTWYHKDBN"
  let buf : FixedArray[UInt16] = FixedArray::make(length, 0)
  let mut i = 0
  while i < length {
    let byte_pos = i / 2
    let byte = encoded.unsafe_get(byte_pos).to_int()
    let code = if i % 2 == 0 { (byte >> 4) & 0xF } else { byte & 0xF }
    buf[i] = lookup.unsafe_get(code)
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
fn parse_bam_header(
  data : String,
  pos : Int,
) -> (BamHeader, Int) raise BamError {
  let (l_text, pos1) = read_int32_le(data, pos)

  if pos1 + l_text > data.length() {
    raise BamError("Truncated BAM header text")
  }
  let text = data[pos1:pos1 + l_text].to_owned()
  let pos2 = pos1 + l_text

  let (n_ref, pos3) = read_int32_le(data, pos2)

  let references : FixedArray[BamReference] = FixedArray::make(n_ref, BamReference::{
    name: "",
    length: 0,
  })
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

    references[ref_idx] = BamReference::{ name, length, }
    ref_idx = ref_idx + 1
    current_pos = pos6
  }

  (BamHeader::{ text, references, }, current_pos)
}

///|
/// Parse a single BAM alignment record (standard htslib layout).
fn parse_bam_record(
  data : String,
  pos : Int,
  end : Int,
) -> (BamRecord, Int) raise BamError {
  let (ref_id, pos1) = read_int32_le(data, pos)
  let (pos_val, pos2) = read_int32_le(data, pos1)
  let (bin_mq_nl, pos3) = read_uint32_le(data, pos2)
  let (flag_nc, pos4) = read_uint32_le(data, pos3)
  let (seq_len, pos5) = read_int32_le(data, pos4)
  let (next_ref_id, pos6) = read_int32_le(data, pos5)
  let (next_pos, pos7) = read_int32_le(data, pos6)
  let (tlen, pos8) = read_int32_le(data, pos7)

  let bin = bin_mq_nl >> 16
  let mapq = (bin_mq_nl >> 8) & 0xFF
  let l_read_name = bin_mq_nl & 0xFF
  let flag = flag_nc >> 16
  let cigar_len = flag_nc & 0xFFFF

  let (qname, pos9) = read_cstring(data, pos8)

  let cigar : FixedArray[BamCigarElement] = FixedArray::make(cigar_len, BamCigarElement::{
    op: BamCigarOp::Match,
    length: 0,
  })
  let mut cigar_idx = 0
  let mut current_pos = pos9
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
  if current_pos + seq_bytes > end {
    raise BamError("Truncated BAM sequence")
  }
  let encoded_seq = data[current_pos:current_pos + seq_bytes].to_owned()
  let pos10 = current_pos + seq_bytes

  let seq = decode_bam_seq(encoded_seq, seq_len)

  let qual_bytes = seq_len
  if pos10 + qual_bytes > end {
    raise BamError("Truncated BAM quality")
  }
  let encoded_qual = data[pos10:pos10 + qual_bytes].to_owned()
  let pos11 = pos10 + qual_bytes

  let qual = decode_bam_qual(encoded_qual, seq_len)

  let (tags, _) = parse_bam_tags(data, pos11, end)

  (
    BamRecord::{
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
      tags,
    },
    pos11,
  )
}
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

    let (record, _) = parse_bam_record(data, pos_new, pos_new + block_size)

    records.push(record)
    current_pos = pos_new + block_size
  }

  BamFile::{ header, records, }
}

///|
/// Get the reference name for a record.
pub fn BamRecord::get_reference_name(
  self : BamRecord,
  header : BamHeader,
) -> String {
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
/// Check if the read fails vendor QC.
pub fn BamRecord::is_qcfail(self : BamRecord) -> Bool {
  (self.flag & 0x200) != 0
}

///|
/// Check if the read is a PCR/optical duplicate.
pub fn BamRecord::is_duplicate(self : BamRecord) -> Bool {
  (self.flag & 0x400) != 0
}

///|
/// Get a CIGAR string representation of this record (e.g. "8M1I5M").
pub fn BamRecord::cigarstring(self : BamRecord) -> String {
  let buf = StringBuilder::new()
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
    buf.write_string(elem.length.to_string())
    buf.write_char(op_char)
  }
  buf.to_string()
}

///|
/// Look up an auxiliary tag by name (e.g. "NM", "MD").
pub fn BamRecord::get_tag(self : BamRecord, tag : String) -> BamTagValue? {
  for t in self.tags {
    if t.name == tag {
      return Some(t.value)
    }
  }
  None
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
    if n >= 3 && s.unsafe_get(n - 2) == eq_char {
      2
    } else {
      1
    }
  } else {
    0
  }
  let out_len = n * 3 / 4 - padding
  let buf : FixedArray[UInt16] = FixedArray::make(out_len, 0)
  let mut i = 0
  let mut j = 0
  while i < n {
    let c0 = base64_table[s.unsafe_get(i).to_int()]
    let c1 = if i + 1 < n {
      base64_table[s.unsafe_get(i + 1).to_int()]
    } else {
      0
    }
    let c2 = if i + 2 < n {
      base64_table[s.unsafe_get(i + 2).to_int()]
    } else {
      0
    }
    let c3 = if i + 3 < n {
      base64_table[s.unsafe_get(i + 3).to_int()]
    } else {
      0
    }

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
  let rnext = if self.next_ref_id >= 0 &&
    self.next_ref_id < header.references.length() {
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

///|
/// CRAM format parsing module.
/// CRAM is a compressed binary format for sequencing reads, similar to BAM but with better compression.

///|
/// CRAM parsing errors.
pub suberror CramError {
  CramError(String)
}

///|
/// CRAM magic number (CRAM\0).
let cram_magic : Array[Int] = [67, 82, 65, 77]

///|
/// CRAM format version.
pub struct CramVersion {
  major : Int
  minor : Int
}

///|
/// Check if data starts with CRAM magic number.
pub fn is_cram_magic(data : String) -> Bool {
  if data.length() < 4 {
    return false
  }
  let mut i = 0
  while i < 4 {
    if data.unsafe_get(i).to_int() != cram_magic[i] {
      return false
    }
    i = i + 1
  }
  true
}

///|
/// CRAM compression type.
pub enum CramCompression {
  Gzip
  Bzip2
  Lzma
  None
  Unknown(Int)
}

///|
/// CRAM reference sequence entry.
pub struct CramReference {
  name : String
  length : Int
  md5 : String
  uri : String
}

///|
/// CRAM file header structure.
pub struct CramHeader {
  version : CramVersion
  reference_count : Int
  references : Array[CramReference]
  compression : CramCompression
  stream_count : Int
  stream_info : Array[String]
  sam_header : String
}

///|
/// CRAM tag (key-value metadata).
pub struct CramTag {
  key : String
  value : String
}

///|
/// CRAM read features (encoding information).
pub struct CramReadFeatures {
  read_group : Int
  flags : Int
  mapping_quality : Int
  read_name : String
  cigar : Array[BamCigarElement]
  sequence : String
  quality : String
  tags : Array[CramTag]
}

///|
/// CRAM alignment record.
pub struct CramRecord {
  ref_id : Int
  pos : Int
  read_features : CramReadFeatures
}

///|
/// CRAM file representation.
pub struct CramFile {
  header : CramHeader
  records : Array[CramRecord]
}

///|
/// Parse CRAM version from data.
fn cram_parse_version(
  data : String,
  pos : Int,
) -> (CramVersion, Int) raise CramError {
  let (major, pos1) = read_uint32_le(data, pos)
  let (minor, pos2) = read_uint32_le(data, pos1)
  (CramVersion::{ major, minor, }, pos2)
}

///|
/// Parse CRAM compression type.
fn cram_parse_compression(code : Int) -> CramCompression {
  match code {
    0 => CramCompression::None
    1 => CramCompression::Gzip
    2 => CramCompression::Bzip2
    3 => CramCompression::Lzma
    _ => CramCompression::Unknown(code)
  }
}

///|
/// Parse CRAM header.
fn cram_parse_header(
  data : String,
  pos : Int,
) -> (CramHeader, Int) raise CramError {
  let (version, pos1) = cram_parse_version(data, pos)
  let (ref_count, pos2) = read_uint32_le(data, pos1)
  let (compression_code, pos3) = read_uint32_le(data, pos2)
  let (stream_count, pos4) = read_uint32_le(data, pos3)

  let compression = cram_parse_compression(compression_code)

  let references : Array[CramReference] = Array::new()
  let mut current_pos = pos4
  let mut ref_idx = 0
  while ref_idx < ref_count {
    let (name, pos5) = read_cstring(data, current_pos)
    let (length, pos6) = read_uint32_le(data, pos5)
    let (md5, pos7) = read_cstring(data, pos6)
    let (uri, pos8) = read_cstring(data, pos7)
    references.push(CramReference::{ name, length, md5, uri, })
    current_pos = pos8
    ref_idx = ref_idx + 1
  }

  let stream_info : Array[String] = Array::new()
  let mut stream_idx = 0
  while stream_idx < stream_count {
    let (info, pos_new) = read_cstring(data, current_pos)
    stream_info.push(info)
    current_pos = pos_new
    stream_idx = stream_idx + 1
  }

  let (sam_header_len, pos9) = read_uint32_le(data, current_pos)
  if pos9 + sam_header_len > data.length() {
    raise CramError("Truncated CRAM SAM header")
  }
  let sam_header = data[pos9:pos9 + sam_header_len].to_owned()
  let pos10 = pos9 + sam_header_len

  (
    CramHeader::{
      version,
      reference_count: ref_count,
      references,
      compression,
      stream_count,
      stream_info,
      sam_header,
    },
    pos10,
  )
}

///|
/// Parse CRAM read features.
fn cram_parse_read_features(
  data : String,
  pos : Int,
) -> (CramReadFeatures, Int) raise CramError {
  let (read_group, pos1) = read_uint32_le(data, pos)
  let (flags, pos2) = read_uint32_le(data, pos1)
  let (mapq, pos3) = read_int8(data, pos2)
  let (name_len, pos4) = read_uint32_le(data, pos3)

  if pos4 + name_len > data.length() {
    raise CramError("Truncated CRAM read name")
  }
  let read_name = data[pos4:pos4 + name_len].to_owned()
  let pos5 = pos4 + name_len

  let (cigar_len, pos6) = read_uint32_le(data, pos5)
  let cigar : Array[BamCigarElement] = Array::new()
  let mut current_pos = pos6
  let mut cigar_idx = 0
  while cigar_idx < cigar_len {
    let (cigar_op_val, pos_new) = read_uint32_le(data, current_pos)
    let length = cigar_op_val >> 4
    let op_code = cigar_op_val & 0xF
    let op = cigar_op_from_code(op_code)
    cigar.push(BamCigarElement::{ op, length, })
    current_pos = pos_new
    cigar_idx = cigar_idx + 1
  }

  let (seq_len, pos7) = read_uint32_le(data, current_pos)
  if pos7 + seq_len > data.length() {
    raise CramError("Truncated CRAM sequence")
  }
  let sequence = data[pos7:pos7 + seq_len].to_owned()
  let pos8 = pos7 + seq_len

  let (qual_len, pos9) = read_uint32_le(data, pos8)
  if pos9 + qual_len > data.length() {
    raise CramError("Truncated CRAM quality")
  }
  let quality = data[pos9:pos9 + qual_len].to_owned()
  let pos10 = pos9 + qual_len

  let (tag_count, pos11) = read_uint32_le(data, pos10)
  let tags : Array[CramTag] = Array::new()
  current_pos = pos11
  let mut tag_idx = 0
  while tag_idx < tag_count {
    let (key, pos_new) = read_cstring(data, current_pos)
    let (value, pos_new2) = read_cstring(data, pos_new)
    tags.push(CramTag::{ key, value, })
    current_pos = pos_new2
    tag_idx = tag_idx + 1
  }

  (
    CramReadFeatures::{
      read_group,
      flags,
      mapping_quality: mapq,
      read_name,
      cigar,
      sequence,
      quality,
      tags,
    },
    current_pos,
  )
}

///|
/// Parse CRAM record.
fn cram_parse_record(
  data : String,
  pos : Int,
) -> (CramRecord, Int) raise CramError {
  let (ref_id, pos1) = read_int32_le(data, pos)
  let (pos_val, pos2) = read_int32_le(data, pos1)
  let (read_features, pos3) = cram_parse_read_features(data, pos2)
  (CramRecord::{ ref_id, pos: pos_val, read_features, }, pos3)
}

///|
/// Parse raw CRAM data.
pub fn parse_cram(data : String) -> CramFile raise CramError {
  let n = data.length()
  if n < 4 {
    raise CramError("CRAM file too short")
  }

  let mut i = 0
  while i < 4 {
    if data.unsafe_get(i).to_int() != cram_magic[i] {
      raise CramError("Invalid CRAM magic number")
    }
    i = i + 1
  }

  let (header, pos) = cram_parse_header(data, 4)

  let records : Array[CramRecord] = Array::new()
  let mut current_pos = pos
  while current_pos < n - 8 {
    let (record_size, pos_new) = read_uint32_le(data, current_pos)
    if pos_new + record_size > n {
      break
    }
    let (record, _) = cram_parse_record(data, pos_new)
    records.push(record)
    current_pos = pos_new + record_size + 4
  }

  CramFile::{ header, records, }
}

///|
/// Parse CRAM from base64-encoded raw data.
pub fn parse_cram_base64(base64_str : String) -> CramFile raise CramError {
  let raw_data = base64_decode(base64_str)
  parse_cram(raw_data)
}

///|
/// Parse CRAM from base64-encoded raw data, returns None on failure.
pub fn try_parse_cram_base64(base64_str : String) -> CramFile? {
  let result = parse_cram_base64(base64_str) catch { _ => return None }
  Some(result)
}

///|
/// Convert CRAM record to BAM record (simplified).
pub fn CramRecord::to_bam_record(self : CramRecord) -> BamRecord {
  let cigar_fixed : FixedArray[BamCigarElement] = FixedArray::make(
    self.read_features.cigar.length(),
    BamCigarElement::{ op: BamCigarOp::Match, length: 0, },
  )
  let mut i = 0
  while i < self.read_features.cigar.length() {
    cigar_fixed[i] = self.read_features.cigar[i]
    i = i + 1
  }
  BamRecord::{
    ref_id: self.ref_id,
    pos: self.pos,
    bin: 0,
    mapq: self.read_features.mapping_quality,
    flag: self.read_features.flags,
    l_read_name: self.read_features.read_name.length(),
    cigar_len: self.read_features.cigar.length(),
    l_seq: self.read_features.sequence.length(),
    next_ref_id: -1,
    next_pos: 0,
    tlen: 0,
    qname: self.read_features.read_name,
    cigar: cigar_fixed,
    seq: self.read_features.sequence,
    qual: self.read_features.quality,
    tags: [],
  }
}

///|
/// Convert CRAM file to BAM file (simplified).
pub fn CramFile::to_bam(self : CramFile) -> BamFile {
  let refs_fixed : FixedArray[BamReference] = FixedArray::make(
    self.header.references.length(),
    BamReference::{ name: "", length: 0, },
  )
  let mut i = 0
  while i < self.header.references.length() {
    refs_fixed[i] = BamReference::{
      name: self.header.references[i].name,
      length: self.header.references[i].length,
    }
    i = i + 1
  }
  let header = BamHeader::{
    text: self.header.sam_header,
    references: refs_fixed,
  }
  let records : Array[BamRecord] = Array::new()
  let mut i = 0
  while i < self.records.length() {
    records.push(self.records[i].to_bam_record())
    i = i + 1
  }
  BamFile::{ header, records, }
}

///|
/// Get number of records in CRAM file.
pub fn CramFile::len(self : CramFile) -> Int {
  self.records.length()
}

///|
/// Get number of reference sequences.
pub fn CramFile::num_references(self : CramFile) -> Int {
  self.header.references.length()
}

///|
/// Get all reference sequence names.
pub fn CramFile::get_reference_names(self : CramFile) -> Array[String] {
  let result = Array::new()
  let mut i = 0
  while i < self.header.references.length() {
    result.push(self.header.references[i].name)
    i = i + 1
  }
  result
}

///|
/// Get reference length by name.
pub fn CramFile::get_reference_length(self : CramFile, name : String) -> Int? {
  let mut i = 0
  while i < self.header.references.length() {
    if self.header.references[i].name == name {
      return Some(self.header.references[i].length)
    }
    i = i + 1
  }
  None
}

///|
/// Check if the read is paired.
pub fn CramRecord::is_paired(self : CramRecord) -> Bool {
  (self.read_features.flags & 0x0001) != 0
}

///|
/// Check if the read is properly paired.
pub fn CramRecord::is_proper_pair(self : CramRecord) -> Bool {
  (self.read_features.flags & 0x0002) != 0
}

///|
/// Check if the read is unmapped.
pub fn CramRecord::is_unmapped(self : CramRecord) -> Bool {
  (self.read_features.flags & 0x0004) != 0
}

///|
/// Check if the read is reverse complement.
pub fn CramRecord::is_reverse(self : CramRecord) -> Bool {
  (self.read_features.flags & 0x010) != 0
}

///|
/// Check if this is the first read in pair.
pub fn CramRecord::is_read1(self : CramRecord) -> Bool {
  (self.read_features.flags & 0x040) != 0
}

///|
/// Check if this is the second read in pair.
pub fn CramRecord::is_read2(self : CramRecord) -> Bool {
  (self.read_features.flags & 0x080) != 0
}

///|
/// Get the reference name for a record.
pub fn CramRecord::get_reference_name(
  self : CramRecord,
  header : CramHeader,
) -> String {
  if self.ref_id >= 0 && self.ref_id < header.references.length() {
    header.references[self.ref_id].name
  } else {
    "*"
  }
}

///|
/// Sample CRAM content for testing.
pub fn sample_cram_content() -> String {
  let magic = "CRAM"
  let version = "\u{0002}\u{0000}\u{0000}\u{0000}"
  let ref_count = "\u{0001}\u{0000}\u{0000}\u{0000}"
  let compression = "\u{0000}\u{0000}\u{0000}\u{0000}"
  let stream_count = "\u{0000}\u{0000}\u{0000}\u{0000}"
  let ref_name = "chr1\u{0000}"
  let ref_length = "\u{0064}\u{0000}\u{0000}\u{0000}"
  let ref_md5 = "\u{0000}"
  let ref_uri = "\u{0000}"
  let sam_header_len = "\u{0014}\u{0000}\u{0000}\u{0000}"
  let sam_header = "@SQ SN:chr1 LN:100\u{0000}"
  let header = magic +
    version +
    ref_count +
    compression +
    stream_count +
    ref_name +
    ref_length +
    ref_md5 +
    ref_uri +
    sam_header_len +
    sam_header
  let ref_id = "\u{0000}\u{0000}\u{0000}\u{0000}"
  let pos = "\u{0001}\u{0000}\u{0000}\u{0000}"
  let read_group = "\u{0000}\u{0000}\u{0000}\u{0000}"
  let flags = "\u{0000}\u{0000}\u{0000}\u{0000}"
  let mapq = "\u{003c}"
  let name_len = "\u{0005}\u{0000}\u{0000}\u{0000}"
  let read_name = "read1\u{0000}"
  let cigar_len = "\u{0001}\u{0000}\u{0000}\u{0000}"
  let cigar_val = "\u{0008}\u{0000}\u{0000}\u{0000}"
  let seq_len = "\u{0008}\u{0000}\u{0000}\u{0000}"
  let sequence = "ATGCATGC"
  let qual_len = "\u{0008}\u{0000}\u{0000}\u{0000}"
  let quality = "!!!!!!!!"
  let tags_len = "\u{0000}\u{0000}\u{0000}\u{0000}"
  let record = ref_id +
    pos +
    read_group +
    flags +
    mapq +
    name_len +
    read_name +
    cigar_len +
    cigar_val +
    seq_len +
    sequence +
    qual_len +
    quality +
    tags_len
  header + record
}

///|
/// Sample CRAM file for testing.
pub fn sample_cram_file() -> CramFile? {
  let result = parse_cram(sample_cram_content()) catch { _ => return None }
  Some(result)
}

///|
/// 0-based (inclusive) reference start position.
pub fn BamRecord::reference_start(self : BamRecord) -> Int {
  self.pos
}

///|
/// 0-based (exclusive) reference end position.
pub fn BamRecord::reference_end(self : BamRecord) -> Int {
  let mut end = self.pos
  for i = 0; i < self.cigar.length(); i = i + 1 {
    let el = self.cigar[i]
    match el.op {
      BamCigarOp::Match | BamCigarOp::Delete | BamCigarOp::Skip | BamCigarOp::Equal | BamCigarOp::Diff =>
        end = end + el.length
      _ => ()
    }
  }
  end
}

///|
/// Reference positions covered by the read. Aligned bases (M/=/X) contribute a
/// position; with `full_length`, insertions/soft clips contribute `None`
/// placeholders. Deletions/skips are skipped entirely.
pub fn BamRecord::get_reference_positions(self : BamRecord, full_length : Bool) -> Array[Int?] {
  let out : Array[Int?] = Array::new()
  let mut rpos = self.pos
  for i = 0; i < self.cigar.length(); i = i + 1 {
    let el = self.cigar[i]
    let n = el.length
    match el.op {
      BamCigarOp::Match | BamCigarOp::Equal | BamCigarOp::Diff => {
        let mut k = 0
        while k < n {
          out.push(Some(rpos + k))
          k = k + 1
        }
        rpos = rpos + n
      }
      BamCigarOp::Insert | BamCigarOp::SoftClip => {
        if full_length {
          let mut k = 0
          while k < n {
            out.push(None)
            k = k + 1
          }
        }
      }
      BamCigarOp::Delete | BamCigarOp::Skip => rpos = rpos + n
      BamCigarOp::HardClip | BamCigarOp::Pad => ()
    }
  }
  out
}

///|
/// Number of aligned (M/=/X) reference bases of the read overlapping [beg, end).
pub fn BamRecord::get_overlap(self : BamRecord, beg : Int, end : Int) -> Int {
  let mut overlap = 0
  let mut rpos = self.pos
  for i = 0; i < self.cigar.length(); i = i + 1 {
    let el = self.cigar[i]
    let n = el.length
    match el.op {
      BamCigarOp::Match | BamCigarOp::Equal | BamCigarOp::Diff => {
        let mut k = 0
        while k < n {
          let p = rpos + k
          if p >= beg && p < end {
            overlap = overlap + 1
          }
          k = k + 1
        }
        rpos = rpos + n
      }
      BamCigarOp::Delete | BamCigarOp::Skip => rpos = rpos + n
      _ => ()
    }
  }
  overlap
}

///|
/// Query↔reference coordinate pairs. Insertions/soft clips yield
/// `(query, None)`, deletions/skips `(None, ref)`. With `matches_only`, only
/// aligned pairs (M/=/X) are kept.
pub fn BamRecord::get_aligned_pairs(self : BamRecord, matches_only : Bool) -> Array[(Int?, Int?)] {
  let out : Array[(Int?, Int?)] = Array::new()
  let mut qpos = 0
  let mut rpos = self.pos
  for i = 0; i < self.cigar.length(); i = i + 1 {
    let el = self.cigar[i]
    let n = el.length
    match el.op {
      BamCigarOp::Match | BamCigarOp::Equal | BamCigarOp::Diff => {
        let mut k = 0
        while k < n {
          out.push((Some(qpos + k), Some(rpos + k)))
          k = k + 1
        }
        qpos = qpos + n
        rpos = rpos + n
      }
      BamCigarOp::Insert | BamCigarOp::SoftClip => {
        if not(matches_only) {
          let mut k = 0
          while k < n {
            out.push((Some(qpos + k), None))
            k = k + 1
          }
        }
        qpos = qpos + n
      }
      BamCigarOp::Delete | BamCigarOp::Skip => {
        if not(matches_only) {
          let mut k = 0
          while k < n {
            out.push((None, Some(rpos + k)))
            k = k + 1
          }
        }
        rpos = rpos + n
      }
      BamCigarOp::HardClip | BamCigarOp::Pad => ()
    }
  }
  out
}
