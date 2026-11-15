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
let bam_magic : Array[Int32] = [66, 65, 77, 1]

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
#derive(Eq)
pub struct BamReference {
  name : String
  length : Int32
}

///|
pub extend BamReference with Eq::{equal, not_equal}

///|
/// BAM header structure.
pub struct BamHeader {
  text : String
  references : FixedArray[BamReference]
}

///|
/// CIGAR operation types for BAM (same as SAM).
#derive(Eq)
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
}

///|
pub extend BamCigarOp with Eq::{equal, not_equal}

///|
/// CIGAR element: operation + length.
#derive(Eq)
pub struct BamCigarElement {
  op : BamCigarOp
  length : Int32
}

///|
pub extend BamCigarElement with Eq::{equal, not_equal}

///|
/// A typed value of a BAM auxiliary tag.
pub enum BamTagValue {
  Char(Char)
  Int(Int32)
  Float(Double)
  Str(String)
  IntArray(Array[Int32])
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
  ref_id : Int32
  pos : Int32
  bin : Int32
  mapq : Int32
  flag : Int32
  l_read_name : Int32
  cigar_len : Int32
  l_seq : Int32
  next_ref_id : Int32
  next_pos : Int32
  tlen : Int32
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
fn cigar_op_from_code(code : Int32) -> BamCigarOp {
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
fn read_int8(data : String, pos : Int32) -> (Int32, Int32) {
  let val = data.unsafe_get(pos).to_int() & 0xFF
  (val, pos + 1)
}

///|
/// Read a UInt16 (little-endian) from data at position, advance position.
fn read_uint16_le(data : String, pos : Int32) -> (Int32, Int32) {
  let b0 = data.unsafe_get(pos).to_int() & 0xFF
  let b1 = data.unsafe_get(pos + 1).to_int() & 0xFF
  let val = ((b1 << 8) | b0) & 0xFFFF
  (val, pos + 2)
}

///|
/// Read an Int32 (little-endian) from data at position, advance position.
fn read_int32_le(data : String, pos : Int32) -> (Int32, Int32) {
  let b0 = data.unsafe_get(pos).to_int() & 0xFF
  let b1 = data.unsafe_get(pos + 1).to_int() & 0xFF
  let b2 = data.unsafe_get(pos + 2).to_int() & 0xFF
  let b3 = data.unsafe_get(pos + 3).to_int() & 0xFF
  let val = b0 | (b1 << 8) | (b2 << 16) | (b3 << 24)
  (val, pos + 4)
}

///|
/// Read a UInt32 (little-endian) from data at position, advance position.
fn read_uint32_le(data : String, pos : Int32) -> (Int32, Int32) {
  let b0 = data.unsafe_get(pos).to_int() & 0xFF
  let b1 = data.unsafe_get(pos + 1).to_int() & 0xFF
  let b2 = data.unsafe_get(pos + 2).to_int() & 0xFF
  let b3 = data.unsafe_get(pos + 3).to_int() & 0xFF
  let val = b0 | (b1 << 8) | (b2 << 16) | (b3 << 24)
  (val & 0xFFFFFFFF, pos + 4)
}

///|
/// Read a null-terminated string from data at position, advance position.
fn read_cstring(data : String, pos : Int32) -> (String, Int32) {
  let mut end = pos
  while end < data.length() && data.unsafe_get(end) != 0 {
    end = end + 1
  }
  let s = data[pos:end].to_owned()
  (s, end + 1)
}

///|
/// Read a signed int8.
fn read_int8_signed(data : String, pos : Int32) -> (Int32, Int32) {
  let b = data.unsafe_get(pos).to_int() & 0xFF
  let v = if b >= 128 { b - 256 } else { b }
  (v, pos + 1)
}

///|
/// Read a signed int16 (little-endian).
fn read_int16_le(data : String, pos : Int32) -> (Int32, Int32) {
  let (u, np) = read_uint16_le(data, pos)
  let v = if u >= 32768 { u - 65536 } else { u }
  (v, np)
}

///|
/// Read an IEEE 754 float32 (4 little-endian bytes) as a Double.
fn read_float32_le(data : String, pos : Int32) -> (Double, Int32) {
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
    sign *
    (1.0 + mantissa.to_double() / 8388608.0) *
    @math.pow(2.0, (exponent - 127).to_double())
  }
  (v, pos + 4)
}

///|
/// Parse BAM auxiliary tags from `pos` until `end`.
fn parse_bam_tags(
  data : String,
  pos : Int32,
  end : Int32,
) -> (Array[BamTag], Int32) raise BamError {
  let tags : Array[BamTag] = Array(capacity=0)
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
      tags.push(BamTag::{
        name,
        type_code,
        value: BamTagValue::Char(data.unsafe_get(p).unsafe_to_char()),
      })
      p = p + 1
    } else if type_code == 'c' {
      let (v, np) = read_int8_signed(data, p)
      tags.push(BamTag::{ name, type_code, value: BamTagValue::Int(v), })
      p = np
    } else if type_code == 'C' {
      let (v, np) = read_int8(data, p)
      tags.push(BamTag::{ name, type_code, value: BamTagValue::Int(v), })
      p = np
    } else if type_code == 's' {
      let (v, np) = read_int16_le(data, p)
      tags.push(BamTag::{ name, type_code, value: BamTagValue::Int(v), })
      p = np
    } else if type_code == 'S' {
      let (v, np) = read_uint16_le(data, p)
      tags.push(BamTag::{ name, type_code, value: BamTagValue::Int(v), })
      p = np
    } else if type_code == 'i' {
      let (v, np) = read_int32_le(data, p)
      tags.push(BamTag::{ name, type_code, value: BamTagValue::Int(v), })
      p = np
    } else if type_code == 'I' {
      let (v, np) = read_uint32_le(data, p)
      tags.push(BamTag::{ name, type_code, value: BamTagValue::Int(v), })
      p = np
    } else if type_code == 'f' {
      let (v, np) = read_float32_le(data, p)
      tags.push(BamTag::{ name, type_code, value: BamTagValue::Float(v), })
      p = np
    } else if type_code == 'Z' || type_code == 'H' {
      let (v, np) = read_cstring(data, p)
      tags.push(BamTag::{ name, type_code, value: BamTagValue::Str(v), })
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
        let arr : Array[Double] = Array(capacity=0)
        for _j = 0; _j < count; _j = _j + 1 {
          let (v, np2) = read_float32_le(data, p)
          arr.push(v)
          p = np2
        }
        tags.push(BamTag::{
          name,
          type_code,
          value: BamTagValue::FloatArray(arr),
        })
      } else {
        let arr : Array[Int32] = Array(capacity=0)
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
        tags.push(BamTag::{
          name,
          type_code,
          value: BamTagValue::IntArray(arr),
        })
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
fn decode_bam_seq(encoded : String, length : Int32) -> String {
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
fn decode_bam_qual(encoded : String, length : Int32) -> String {
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
  pos : Int32,
) -> (BamHeader, Int32) raise BamError {
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
  pos : Int32,
  end : Int32,
) -> (BamRecord, Int32) raise BamError {
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

///|
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

  let records = Array(capacity=0)
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
  let buf = StringBuilder()
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
pub fn BamFile::len(self : BamFile) -> Int32 {
  self.records.length()
}

///|
/// Get number of reference sequences.
pub fn BamFile::num_references(self : BamFile) -> Int32 {
  self.header.references.length()
}

///|
/// Get all reference sequence names.
pub fn BamFile::get_reference_names(self : BamFile) -> Array[String] {
  let result = Array(capacity=0)
  for ref_ in self.header.references {
    result.push(ref_.name)
  }
  result
}

///|
/// Get reference length by name.
pub fn BamFile::get_reference_length(self : BamFile, name : String) -> Int32? {
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
let base64_table : FixedArray[Int32] = {
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
pub fn BamRecord::get_end(self : BamRecord) -> Int32 {
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
pub fn BamRecord::get_alignment_length(self : BamRecord) -> Int32 {
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
  let cigar_str = StringBuilder()
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
  let sb = StringBuilder()
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
let cram_magic : Array[Int32] = [67, 82, 65, 77]

///|
/// CRAM format version.
pub struct CramVersion {
  major : Int32
  minor : Int32
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
  Unknown(Int32)
}

///|
/// CRAM reference sequence entry.
pub struct CramReference {
  name : String
  length : Int32
  md5 : String
  uri : String
}

///|
/// CRAM file header structure.
pub struct CramHeader {
  version : CramVersion
  reference_count : Int32
  references : Array[CramReference]
  compression : CramCompression
  stream_count : Int32
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
  read_group : Int32
  flags : Int32
  mapping_quality : Int32
  read_name : String
  cigar : Array[BamCigarElement]
  sequence : String
  quality : String
  tags : Array[CramTag]
}

///|
/// CRAM alignment record.
pub struct CramRecord {
  ref_id : Int32
  pos : Int32
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
fn cram_parse_version(data : String, pos : Int32) -> (CramVersion, Int32) {
  let (major, pos1) = read_uint32_le(data, pos)
  let (minor, pos2) = read_uint32_le(data, pos1)
  (CramVersion::{ major, minor, }, pos2)
}

///|
/// Parse CRAM compression type.
fn cram_parse_compression(code : Int32) -> CramCompression {
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
  pos : Int32,
) -> (CramHeader, Int32) raise CramError {
  let (version, pos1) = cram_parse_version(data, pos)
  let (ref_count, pos2) = read_uint32_le(data, pos1)
  let (compression_code, pos3) = read_uint32_le(data, pos2)
  let (stream_count, pos4) = read_uint32_le(data, pos3)

  let compression = cram_parse_compression(compression_code)

  let references : Array[CramReference] = Array(capacity=0)
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

  let stream_info : Array[String] = Array(capacity=0)
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
  pos : Int32,
) -> (CramReadFeatures, Int32) raise CramError {
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
  let cigar : Array[BamCigarElement] = Array(capacity=0)
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
  let tags : Array[CramTag] = Array(capacity=0)
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
  pos : Int32,
) -> (CramRecord, Int32) raise CramError {
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

  let records : Array[CramRecord] = Array(capacity=0)
  let mut current_pos = pos
  while current_pos < n - 8 {
    let (record_size, pos_new) = read_uint32_le(data, current_pos)
    if pos_new + record_size > n {
      break
    }
    let (record, _) = cram_parse_record(data, pos_new)
    records.push(record)
    current_pos = pos_new + record_size
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
  let records : Array[BamRecord] = Array(capacity=0)
  let mut i = 0
  while i < self.records.length() {
    records.push(self.records[i].to_bam_record())
    i = i + 1
  }
  BamFile::{ header, records, }
}

///|
/// Get number of records in CRAM file.
pub fn CramFile::len(self : CramFile) -> Int32 {
  self.records.length()
}

///|
/// Get number of reference sequences.
pub fn CramFile::num_references(self : CramFile) -> Int32 {
  self.header.references.length()
}

///|
/// Get all reference sequence names.
pub fn CramFile::get_reference_names(self : CramFile) -> Array[String] {
  let result = Array(capacity=0)
  let mut i = 0
  while i < self.header.references.length() {
    result.push(self.header.references[i].name)
    i = i + 1
  }
  result
}

///|
/// Get reference length by name.
pub fn CramFile::get_reference_length(self : CramFile, name : String) -> Int32? {
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
  // Build CRAM binary programmatically (avoids manual escape-sequence errors).
  let buf : Array[Int32] = []
  // Magic "CRAM"
  buf.push(67)
  buf.push(82)
  buf.push(65)
  buf.push(77)
  // Version: major=2 (uint32_le) + minor=0 (uint32_le)
  buf.push(2)
  buf.push(0)
  buf.push(0)
  buf.push(0)
  buf.push(0)
  buf.push(0)
  buf.push(0)
  buf.push(0)
  // ref_count=1
  buf.push(1)
  buf.push(0)
  buf.push(0)
  buf.push(0)
  // compression=0 (None)
  buf.push(0)
  buf.push(0)
  buf.push(0)
  buf.push(0)
  // stream_count=0
  buf.push(0)
  buf.push(0)
  buf.push(0)
  buf.push(0)
  // ref: name="chr1" (cstring), length=100, md5="" (cstring), uri="" (cstring)
  buf.push(99) // c
  buf.push(104) // h
  buf.push(114) // r
  buf.push(49) // 1
  buf.push(0) // null
  // ref_length=100 = 0x64
  buf.push(100)
  buf.push(0)
  buf.push(0)
  buf.push(0)
  // md5="" (empty cstring)
  buf.push(0)
  // uri="" (empty cstring)
  buf.push(0)
  // sam_header = "@SQ\tSN:chr1\tLN:100" (18 bytes)
  let hdr = "@SQ\tSN:chr1\tLN:100"
  buf.push(hdr.length())
  buf.push(0)
  buf.push(0)
  buf.push(0)
  for i = 0; i < hdr.length(); i = i + 1 {
    buf.push(hdr.unsafe_get(i).to_int() & 0xFF)
  }
  // Record body (built separately to know its size)
  let body : Array[Int32] = []
  // ref_id=0 (int32_le)
  body.push(0)
  body.push(0)
  body.push(0)
  body.push(0)
  // pos=1 (int32_le)
  body.push(1)
  body.push(0)
  body.push(0)
  body.push(0)
  // read_group=0
  body.push(0)
  body.push(0)
  body.push(0)
  body.push(0)
  // flags=0
  body.push(0)
  body.push(0)
  body.push(0)
  body.push(0)
  // mapq=60 (int8)
  body.push(60)
  // read_name: len=5 + "read1"
  body.push(5)
  body.push(0)
  body.push(0)
  body.push(0)
  body.push(114) // r
  body.push(101) // e
  body.push(97) // a
  body.push(100) // d
  body.push(49) // 1
  // cigar: count=1, val=100M = (100<<4)|0 = 1600 = 0x640
  body.push(1)
  body.push(0)
  body.push(0)
  body.push(0)
  body.push(0x40)
  body.push(0x06)
  body.push(0)
  body.push(0)
  // seq: len=8, "ATGCATGC"
  body.push(8)
  body.push(0)
  body.push(0)
  body.push(0)
  body.push(65) // A
  body.push(84) // T
  body.push(71) // G
  body.push(67) // C
  body.push(65) // A
  body.push(84) // T
  body.push(71) // G
  body.push(67) // C
  // qual: len=8, "!!!!!!!!"
  body.push(8)
  body.push(0)
  body.push(0)
  body.push(0)
  for _ in 0..<8 {
    body.push(33) // !
  }
  // tags: count=0
  body.push(0)
  body.push(0)
  body.push(0)
  body.push(0)
  // Record size prefix + body
  buf.push(body.length())
  buf.push(0)
  buf.push(0)
  buf.push(0)
  for i = 0; i < body.length(); i = i + 1 {
    buf.push(body[i])
  }
  // Convert to String
  let sb = StringBuilder()
  for i = 0; i < buf.length(); i = i + 1 {
    sb.write_char(buf[i].to_uint16().unsafe_to_char())
  }
  sb.to_string()
}

///|
/// Sample CRAM file for testing.
pub fn sample_cram_file() -> CramFile? {
  let result = parse_cram(sample_cram_content()) catch { _ => return None }
  Some(result)
}

///|
/// 0-based (inclusive) reference start position.
pub fn BamRecord::reference_start(self : BamRecord) -> Int32 {
  self.pos
}

///|
/// 0-based (exclusive) reference end position.
pub fn BamRecord::reference_end(self : BamRecord) -> Int32 {
  let mut end = self.pos
  for i = 0; i < self.cigar.length(); i = i + 1 {
    let el = self.cigar[i]
    match el.op {
      BamCigarOp::Match
      | BamCigarOp::Delete
      | BamCigarOp::Skip
      | BamCigarOp::Equal
      | BamCigarOp::Diff => end = end + el.length
      _ => ()
    }
  }
  end
}

///|
/// Reference positions covered by the read. Aligned bases (M/=/X) contribute a
/// position; with `full_length`, insertions/soft clips contribute `None`
/// placeholders. Deletions/skips are skipped entirely.
pub fn BamRecord::get_reference_positions(
  self : BamRecord,
  full_length : Bool,
) -> Array[Int32?] {
  let out : Array[Int32?] = Array(capacity=0)
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
      BamCigarOp::Insert | BamCigarOp::SoftClip =>
        if full_length {
          let mut k = 0
          while k < n {
            out.push(None)
            k = k + 1
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
pub fn BamRecord::get_overlap(
  self : BamRecord,
  beg : Int32,
  end : Int32,
) -> Int32 {
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
pub fn BamRecord::get_aligned_pairs(
  self : BamRecord,
  matches_only : Bool,
) -> Array[(Int32?, Int32?)] {
  let out : Array[(Int32?, Int32?)] = Array(capacity=0)
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
        if !matches_only {
          let mut k = 0
          while k < n {
            out.push((Some(qpos + k), None))
            k = k + 1
          }
        }
        qpos = qpos + n
      }
      BamCigarOp::Delete | BamCigarOp::Skip => {
        if !matches_only {
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

///|
pub struct PileupRead {
  record_index : Int32
  query_position : Int32?
  is_del : Bool
  is_refskip : Bool
  indel : Int32
}

///|
pub struct PileupColumn {
  reference_pos : Int32
  nsegments : Int32
  pileups : Array[PileupRead]
}

///|
/// The state of a record at reference position `p`, or None if it does not
/// cover `p`. Returns (query_position, is_del, is_refskip, indel).
fn bam_pileup_state(r : BamRecord, p : Int32) -> (Int32?, Bool, Bool, Int32)? {
  let mut qpos = 0
  let mut rpos = r.pos
  let mut i = 0
  while i < r.cigar.length() {
    let el = r.cigar[i]
    let n = el.length
    match el.op {
      BamCigarOp::Match | BamCigarOp::Equal | BamCigarOp::Diff => {
        if p >= rpos && p < rpos + n {
          let q = qpos + (p - rpos)
          // indel is reported on the base immediately preceding an I/D
          let ind = if p == rpos + n - 1 {
            if i + 1 < r.cigar.length() {
              match r.cigar[i + 1].op {
                BamCigarOp::Insert => r.cigar[i + 1].length
                BamCigarOp::Delete => -r.cigar[i + 1].length
                _ => 0
              }
            } else {
              0
            }
          } else {
            0
          }
          return Some((Some(q), false, false, ind))
        }
        qpos = qpos + n
        rpos = rpos + n
      }
      BamCigarOp::Insert | BamCigarOp::SoftClip => qpos = qpos + n
      BamCigarOp::Delete => {
        if p >= rpos && p < rpos + n {
          return Some((None, true, false, 0))
        }
        rpos = rpos + n
      }
      BamCigarOp::Skip => {
        if p >= rpos && p < rpos + n {
          return Some((None, true, true, 0))
        }
        rpos = rpos + n
      }
      BamCigarOp::HardClip | BamCigarOp::Pad => ()
    }
    i = i + 1
  }
  None
}

///|
/// Pileup of all records on `contig` over the half-open region [start, stop).
/// Returns one column per covered reference position, preserving record order
/// (file order), matching pysam's `stepper="all"` no-filtering behavior.
pub fn BamFile::pileup(
  self : BamFile,
  contig : String,
  start : Int32,
  stop : Int32,
) -> Array[PileupColumn] {
  let mut ref_id = -1
  for j = 0; j < self.header.references.length(); j = j + 1 {
    if self.header.references[j].name == contig {
      ref_id = j
      break
    }
  }
  let columns : Array[PileupColumn] = Array(capacity=0)
  for p = start; p < stop; p = p + 1 {
    let pileups : Array[PileupRead] = Array(capacity=0)
    for i = 0; i < self.records.length(); i = i + 1 {
      let r = self.records[i]
      if r.ref_id != ref_id {
        continue
      }
      match bam_pileup_state(r, p) {
        Some((qp, is_del, is_refskip, ind)) =>
          pileups.push(PileupRead::{
            record_index: i,
            query_position: qp,
            is_del,
            is_refskip,
            indel: ind,
          })
        None => ()
      }
    }
    if pileups.length() > 0 {
      columns.push(PileupColumn::{
        reference_pos: p,
        nsegments: pileups.length(),
        pileups,
      })
    }
  }
  columns
}

///|
/// The full query sequence (includes soft-clipped bases, no hard clips).
pub fn BamRecord::query_sequence(self : BamRecord) -> String {
  self.seq
}

///|
/// The full query base qualities (aligned with `query_sequence`).
pub fn BamRecord::query_qualities(self : BamRecord) -> String {
  self.qual
}

///|
/// Stored query length (`l_seq`).
pub fn BamRecord::query_length(self : BamRecord) -> Int32 {
  self.l_seq
}

///|
/// Query length inferred from the CIGAR (M + I + S + = + X), matching htslib's
/// `bam_cigar2qlen`.
pub fn BamRecord::infer_query_length(self : BamRecord) -> Int32 {
  let mut length = 0
  for i = 0; i < self.cigar.length(); i = i + 1 {
    let el = self.cigar[i]
    match el.op {
      BamCigarOp::Match
      | BamCigarOp::Insert
      | BamCigarOp::SoftClip
      | BamCigarOp::Equal
      | BamCigarOp::Diff => length = length + el.length
      _ => ()
    }
  }
  length
}

///|
pub struct FlagstatResult {
  total : Int32
  primary : Int32
  secondary : Int32
  supplementary : Int32
  duplicates : Int32
  primary_duplicates : Int32
  mapped : Int32
  primary_mapped : Int32
  paired : Int32
  read1 : Int32
  read2 : Int32
  properly_paired : Int32
  both_mapped : Int32
  singletons : Int32
  mate_diff_chr : Int32
  mate_diff_chr_mapq5 : Int32
}

///|
/// Flag statistics over all records, matching `samtools flagstat` (bam_stat.c
/// `flagstat_loop`). Counts are combined QC-pass + QC-fail.
pub fn BamFile::flagstat(self : BamFile) -> FlagstatResult {
  let mut total = 0
  let mut primary = 0
  let mut secondary = 0
  let mut supplementary = 0
  let mut duplicates = 0
  let mut primary_duplicates = 0
  let mut mapped = 0
  let mut primary_mapped = 0
  let mut paired = 0
  let mut read1 = 0
  let mut read2 = 0
  let mut properly_paired = 0
  let mut both_mapped = 0
  let mut singletons = 0
  let mut mate_diff_chr = 0
  let mut mate_diff_chr_mapq5 = 0
  for r in self.records {
    let flag = r.flag
    total = total + 1
    if (flag & 0x100) != 0 {
      secondary = secondary + 1
    } else if (flag & 0x800) != 0 {
      supplementary = supplementary + 1
    } else {
      primary = primary + 1
      if (flag & 0x1) != 0 {
        paired = paired + 1
        if (flag & 0x2) != 0 && (flag & 0x4) == 0 {
          properly_paired = properly_paired + 1
        }
        if (flag & 0x40) != 0 {
          read1 = read1 + 1
        }
        if (flag & 0x80) != 0 {
          read2 = read2 + 1
        }
        if (flag & 0x8) != 0 && (flag & 0x4) == 0 {
          singletons = singletons + 1
        }
        if (flag & 0x4) == 0 && (flag & 0x8) == 0 {
          both_mapped = both_mapped + 1
          if r.ref_id != r.next_ref_id {
            mate_diff_chr = mate_diff_chr + 1
            if r.mapq >= 5 {
              mate_diff_chr_mapq5 = mate_diff_chr_mapq5 + 1
            }
          }
        }
      }
      if (flag & 0x4) == 0 {
        primary_mapped = primary_mapped + 1
      }
      if (flag & 0x400) != 0 {
        primary_duplicates = primary_duplicates + 1
      }
    }
    if (flag & 0x4) == 0 {
      mapped = mapped + 1
    }
    if (flag & 0x400) != 0 {
      duplicates = duplicates + 1
    }
  }
  FlagstatResult::{
    total,
    primary,
    secondary,
    supplementary,
    duplicates,
    primary_duplicates,
    mapped,
    primary_mapped,
    paired,
    read1,
    read2,
    properly_paired,
    both_mapped,
    singletons,
    mate_diff_chr,
    mate_diff_chr_mapq5,
  }
}

///|
pub struct IdxstatsResult {
  ref_name : String
  ref_length : Int32
  mapped : Int32
  unmapped : Int32
}

///|
/// Per-reference read counts (mapped/unmapped), matching `samtools idxstats`
/// (bam_index.c `slow_idxstats`). A trailing "*" row carries unplaced reads.
pub fn BamFile::idxstats(self : BamFile) -> Array[IdxstatsResult] {
  let nref = self.header.references.length()
  let mapped_cnt = Array::make(nref + 1, 0)
  let unmapped_cnt = Array::make(nref + 1, 0)
  for r in self.records {
    let b = if r.ref_id < 0 { nref } else { r.ref_id }
    if (r.flag & 0x4) != 0 {
      unmapped_cnt[b] = unmapped_cnt[b] + 1
    } else {
      mapped_cnt[b] = mapped_cnt[b] + 1
    }
  }
  let out : Array[IdxstatsResult] = Array(capacity=0)
  for i = 0; i < nref; i = i + 1 {
    out.push(IdxstatsResult::{
      ref_name: self.header.references[i].name,
      ref_length: self.header.references[i].length,
      mapped: mapped_cnt[i],
      unmapped: unmapped_cnt[i],
    })
  }
  out.push(IdxstatsResult::{
    ref_name: "*",
    ref_length: 0,
    mapped: mapped_cnt[nref],
    unmapped: unmapped_cnt[nref],
  })
  out
}

///|
pub struct MpileupLine {
  chrom : String
  pos : Int32
  ref_base : String
  depth : Int32
  bases : String
  quals : String
}

///|
/// Build a single-character string from a code point. Used for the mapping
/// quality byte that follows the `^` read-start marker (and any other raw
/// base-quality byte), which `char_at` cannot produce for non-alphabetic
/// codes.
fn char_from_code(c : Int32) -> String {
  let buf : FixedArray[UInt16] = FixedArray::make(1, 0)
  buf[0] = c.to_uint16()
  unsafe_from_uint16_fixedarray(buf)
}

///|
/// Reference-position column pileup matching `samtools mpileup` / `bcftools
/// mpileup` tabular output (columns 1-6). Builds on `BamFile::pileup`
/// (`stepper="all"` semantics), and for each covered position emits:
///   - `chrom`, 1-based `pos`, uppercase `ref_base` (from `ref_seq`, 0-based),
///   - `depth` = number of reads covering (including deletions/refskips),
///   - `bases`: `.`/`,` forward/reverse match, `A..T`/`a..t` mismatch,
///     `N`/`n` unknown base, `>`/`<` forward/reverse refskip, `*` deletion,
///     `^`+mapq read start, `$` read end, `+<len><seq>` insertion and
///     `-<len><ref>` deletion markers (matching htslib `mplp_func`),
///   - `quals`: one Phred+33 byte per emitted base (no qual for
///     deletion/refskip/indel/start/end markers).
/// `ref_seq` is the reference sequence of `contig`; it must cover `[start, stop)`.
pub fn BamFile::mpileup(
  self : BamFile,
  contig : String,
  start : Int32,
  stop : Int32,
  ref_seq : String,
) -> Array[MpileupLine] {
  let columns = self.pileup(contig, start, stop)
  let lines : Array[MpileupLine] = Array(capacity=0)
  for col in columns {
    let p = col.reference_pos
    let ref_base = if p < ref_seq.length() {
      to_uppercase(ref_seq[p:p + 1].to_owned())
    } else {
      "N"
    }
    let bases = StringBuilder()
    let quals = StringBuilder()
    for pr in col.pileups {
      let r = self.records[pr.record_index]
      let rev = r.is_reverse()
      if p == r.pos {
        bases.write_char('^')
        bases.write_string(char_from_code(r.mapq + 33))
      }
      if pr.is_del {
        bases.write_char('*')
      } else if pr.is_refskip {
        bases.write_char(if rev { '<' } else { '>' })
      } else {
        match pr.query_position {
          None => ()
          Some(qpos) => {
            let read_base = if qpos < r.seq.length() {
              r.seq[qpos:qpos + 1].to_owned()
            } else {
              "N"
            }
            let qual_char = if qpos < r.qual.length() {
              r.qual[qpos:qpos + 1].to_owned()
            } else {
              "!"
            }
            if read_base == ref_base || read_base == "=" {
              bases.write_char(if rev { ',' } else { '.' })
              quals.write_string(qual_char)
            } else {
              bases.write_string(
                if rev {
                  to_lowercase(read_base)
                } else {
                  read_base
                },
              )
              quals.write_string(qual_char)
            }
          }
        }
      }
      if pr.indel > 0 {
        match pr.query_position {
          None => ()
          Some(qpos) => {
            bases.write_char('+')
            bases.write_string(pr.indel.to_string())
            bases.write_string(r.seq[qpos + 1:qpos + 1 + pr.indel].to_owned())
          }
        }
      } else if pr.indel < 0 {
        let n = -pr.indel
        bases.write_char('-')
        bases.write_string(n.to_string())
        bases.write_string(ref_seq[p + 1:p + 1 + n].to_owned())
      }
      if p == r.get_end() {
        bases.write_char('$')
      }
    }
    lines.push(MpileupLine::{
      chrom: contig,
      pos: p + 1,
      ref_base,
      depth: col.pileups.length(),
      bases: bases.to_string(),
      quals: quals.to_string(),
    })
  }
  lines
}

// ============================================================================
// samtools 子命令：sort / merge / depth / fastq
// ============================================================================

///|
/// Construct a BamFile from a header and records.
pub fn BamFile::new(header : BamHeader, records : Array[BamRecord]) -> BamFile {
  BamFile::{ header, records, }
}

///|
/// Coordinate-sort a BAM file by `(ref_id, pos)` — `samtools sort` (header
/// preserved, records reordered in place).
pub fn bam_sort(bam : BamFile) -> BamFile {
  let records = bam.records.copy()
  let mut i = 1
  while i < records.length() {
    let key = records[i]
    let mut j = i - 1
    while j >= 0 && bam_record_less(key, records[j]) {
      records[j + 1] = records[j]
      j = j - 1
    }
    records[j + 1] = key
    i = i + 1
  }
  BamFile::{ header: bam.header, records, }
}

///|
fn bam_record_less(a : BamRecord, b : BamRecord) -> Bool {
  if a.ref_id != b.ref_id {
    a.ref_id < b.ref_id
  } else {
    a.pos < b.pos
  }
}

///|
/// Merge multiple BAM files into one (`samtools merge`): records concatenated,
/// header taken from the first file.
pub fn bam_merge(files : Array[BamFile]) -> BamFile {
  if files.length() == 0 {
    return BamFile::{
      header: BamHeader::{
        text: "",
        references: FixedArray::make(0, BamReference::{ name: "", length: 0, }),
      },
      records: [],
    }
  }
  let header = files[0].header
  let records : Array[BamRecord] = Array(capacity=0)
  for f in files {
    for r in f.records {
      records.push(r)
    }
  }
  BamFile::{ header, records, }
}

///|
/// Per-base coverage depth on one contig — `samtools depth`. Returns
/// `(1-based position, depth)` for positions with depth > 0, via a difference
/// array over each record's reference span.
pub fn bam_depth(bam : BamFile, contig : String) -> Array[(Int32, Int32)] {
  let mut ref_id = -1
  let mut ref_len = 0
  for j = 0; j < bam.header.references.length(); j = j + 1 {
    if bam.header.references[j].name == contig {
      ref_id = j
      ref_len = bam.header.references[j].length
      break
    }
  }
  if ref_id < 0 {
    return []
  }
  let diff : Array[Int32] = Array::make(ref_len + 1, 0)
  for r in bam.records {
    if r.ref_id != ref_id {
      continue
    }
    let s = r.pos
    let e = r.get_end() + 1
    if e > s {
      diff[s] = diff[s] + 1
      diff[e] = diff[e] - 1
    }
  }
  let out : Array[(Int32, Int32)] = Array(capacity=0)
  let mut depth = 0
  let mut p = 0
  while p < ref_len {
    depth = depth + diff[p]
    if depth > 0 {
      out.push((p + 1, depth))
    }
    p = p + 1
  }
  out
}

///|
/// Write reads from a BAM file as FASTQ (`samtools fastq`); unmapped reads are
/// skipped. Sequences carry their phred qualities.
pub fn bam_fastq(bam : BamFile) -> String raise SeqIOError {
  let records : Array[SeqRecord] = Array(capacity=0)
  for rec in bam.records {
    if rec.is_unmapped() {
      continue
    }
    let qual_arr : Array[Int32] = Array(capacity=0)
    let qn = rec.qual.length()
    let mut qi = 0
    while qi < qn {
      qual_arr.push(rec.qual.unsafe_get(qi).to_int() - 33)
      qi = qi + 1
    }
    let letter_ann : Map[String, Array[Int32]] = Map([], capacity=0)
    letter_ann.set("phred_quality", qual_arr)
    records.push(SeqRecord::{
      seq: Seq::new(rec.seq),
      id: rec.qname,
      name: rec.qname,
      description: "",
      letter_annotations: letter_ann,
      annotations: Map([], capacity=0),
      features: [],
    })
  }
  write_fastq(records)
}

// ============================================================================
// samtools 子命令：markdup / stats
// ============================================================================

///|
fn bam_record_with_flag(rec : BamRecord, flag : Int32) -> BamRecord {
  BamRecord::{
    ref_id: rec.ref_id,
    pos: rec.pos,
    bin: rec.bin,
    mapq: rec.mapq,
    flag,
    l_read_name: rec.l_read_name,
    cigar_len: rec.cigar_len,
    l_seq: rec.l_seq,
    next_ref_id: rec.next_ref_id,
    next_pos: rec.next_pos,
    tlen: rec.tlen,
    qname: rec.qname,
    cigar: rec.cigar,
    seq: rec.seq,
    qual: rec.qual,
    tags: rec.tags,
  }
}

///|
/// Mark duplicate reads (`samtools markdup`, simplified): after coordinate
/// sorting, a read identical in `(ref_id, pos, seq)` to its predecessor is
/// flagged as duplicate (`0x400`).
pub fn bam_markdup(bam : BamFile) -> BamFile {
  let sorted = bam_sort(bam)
  let records = sorted.records.copy()
  let out : Array[BamRecord] = Array(capacity=0)
  let mut i = 0
  while i < records.length() {
    let r = records[i]
    let is_dup = if i > 0 {
      let prev = records[i - 1]
      prev.ref_id == r.ref_id && prev.pos == r.pos && prev.seq == r.seq
    } else {
      false
    }
    if is_dup {
      out.push(bam_record_with_flag(r, r.flag | 0x400))
    } else {
      out.push(r)
    }
    i = i + 1
  }
  BamFile::{ header: bam.header, records: out, }
}

///|
/// Summary statistics.
pub struct BamStatsResult {
  total_reads : Int32
  mapped_reads : Int32
  duplicate_reads : Int32
  avg_read_length : Double
}

///|
/// Basic BAM summary (`samtools stats` subset): total / mapped / duplicate
/// counts plus average read length.
pub fn bam_stats(bam : BamFile) -> BamStatsResult {
  let mut total = 0
  let mut mapped = 0
  let mut dup = 0
  let mut total_len = 0
  for r in bam.records {
    total = total + 1
    if !r.is_unmapped() {
      mapped = mapped + 1
    }
    if r.is_duplicate() {
      dup = dup + 1
    }
    total_len = total_len + r.seq.length()
  }
  let avg = if total == 0 {
    0.0
  } else {
    total_len.to_double() / total.to_double()
  }
  BamStatsResult::{
    total_reads: total,
    mapped_reads: mapped,
    duplicate_reads: dup,
    avg_read_length: avg,
  }
}

// ============================================================================
// pysam 高频访问：fetch / count_coverage / get_base_counts
// ============================================================================

///|
/// Fetch reads overlapping a region (`AlignmentFile.fetch`): 0-based half-open
/// `[start, stop)`. Returns records in file order whose reference span
/// intersects the region.
pub fn bam_fetch(
  bam : BamFile,
  contig : String,
  start : Int32,
  stop : Int32,
) -> Array[BamRecord] {
  let mut ref_id = -1
  for j = 0; j < bam.header.references.length(); j = j + 1 {
    if bam.header.references[j].name == contig {
      ref_id = j
      break
    }
  }
  let out : Array[BamRecord] = Array(capacity=0)
  for r in bam.records {
    if r.ref_id == ref_id && r.pos < stop && r.get_end() + 1 > start {
      out.push(r)
    }
  }
  out
}

///|
/// Per-base read coverage split by base (`AlignmentFile.count_coverage`):
/// returns four arrays (A, C, G, T) of length `stop - start`, each entry the
/// number of reads with that base at the position (deletions/ref-skips do not
/// contribute a base).
pub fn bam_count_coverage(
  bam : BamFile,
  contig : String,
  start : Int32,
  stop : Int32,
) -> (Array[Int32], Array[Int32], Array[Int32], Array[Int32]) {
  let n = stop - start
  let a = Array::make(n, 0)
  let c = Array::make(n, 0)
  let g = Array::make(n, 0)
  let t = Array::make(n, 0)
  let columns = bam.pileup(contig, start, stop)
  for col in columns {
    let idx = col.reference_pos - start
    if idx < 0 || idx >= n {
      continue
    }
    for pu in col.pileups {
      match pu.query_position {
        Some(qp) => {
          let read = bam.records[pu.record_index]
          let seq = read.query_sequence()
          if qp >= 0 && qp < seq.length() {
            let base = seq.unsafe_get(qp)
            let upper = if base >= 97 && base <= 122 { base - 32 } else { base }
            if upper == 65 {
              a[idx] = a[idx] + 1
            } else if upper == 67 {
              c[idx] = c[idx] + 1
            } else if upper == 71 {
              g[idx] = g[idx] + 1
            } else if upper == 84 {
              t[idx] = t[idx] + 1
            }
          }
        }
        None => ()
      }
    }
  }
  (a, c, g, t)
}

///|
/// Base counts at a single position (`(A, C, G, T)`), the per-position slice of
/// `bam_count_coverage`.
pub fn bam_get_base_counts(
  bam : BamFile,
  contig : String,
  pos : Int32,
) -> (Int32, Int32, Int32, Int32) {
  let (a, c, g, t) = bam_count_coverage(bam, contig, pos, pos + 1)
  let va = if a.length() > 0 { a[0] } else { 0 }
  let vc = if c.length() > 0 { c[0] } else { 0 }
  let vg = if g.length() > 0 { g[0] } else { 0 }
  let vt = if t.length() > 0 { t[0] } else { 0 }
  (va, vc, vg, vt)
}

// ============================================================================
// get_alignment_quality（NM tag 校准）
// ============================================================================

///|
/// Alignment quality (`AlignedSegment.get_alignment_quality`): `mapq` penalized
/// by 10 per edit (the NM tag), clamped to `[0, 100]`. Returns the raw mapq
/// when no NM tag is present.
pub fn bam_alignment_quality(rec : BamRecord) -> Double {
  let nm = match rec.get_tag("NM") {
    Some(BamTagValue::Int(v)) => v
    _ => 0
  }
  let q = rec.mapq.to_double() - 10.0 * nm.to_double()
  if q < 0.0 {
    0.0
  } else if q > 100.0 {
    100.0
  } else {
    q
  }
}

// ============================================================================
// P0 核心 API 补全：AlignedSegment 高频方法（query_alignment_* / cigar_stats /
// blocks / infer_read_length / get_reference_sequence / get_forward_* /
// has_tag / set_tag）
// ============================================================================

///|
/// Total query-consuming length of a CIGAR op (M/I/S/=/X; excludes D/N/H/P).
fn bam_cigar_op_qlen(op : BamCigarOp) -> Int32 {
  match op {
    BamCigarOp::Match
    | BamCigarOp::Insert
    | BamCigarOp::SoftClip
    | BamCigarOp::Equal
    | BamCigarOp::Diff => 1
    _ => 0
  }
}

///|
/// Total reference-consuming length of a CIGAR op (M/D/N/=/X; excludes I/S/H/P).
fn bam_cigar_op_rlen(op : BamCigarOp) -> Int32 {
  match op {
    BamCigarOp::Match
    | BamCigarOp::Delete
    | BamCigarOp::Skip
    | BamCigarOp::Equal
    | BamCigarOp::Diff => 1
    _ => 0
  }
}

///|
/// Read length inferred from CIGAR including H clips (`bam_cigar2rlen`):
/// sum of M/I/S/H/=/X lengths (Pysam `infer_read_length`).
pub fn BamRecord::infer_read_length(self : BamRecord) -> Int32 {
  let mut len = 0
  for el in self.cigar {
    match el.op {
      BamCigarOp::Match
      | BamCigarOp::Insert
      | BamCigarOp::SoftClip
      | BamCigarOp::HardClip
      | BamCigarOp::Equal
      | BamCigarOp::Diff => len = len + el.length
      _ => ()
    }
  }
  len
}

///|
/// First aligned query position (0-based), i.e. query offset of the first
/// non-softclip op. Returns None when the alignment contains no aligned bases
/// (only H/S clips). Matches Pysam `AlignedSegment.query_alignment_start`.
pub fn BamRecord::query_alignment_start(self : BamRecord) -> Int32? {
  let mut qpos = 0
  for el in self.cigar {
    match el.op {
      BamCigarOp::Match
      | BamCigarOp::Equal
      | BamCigarOp::Diff => return Some(qpos)
      BamCigarOp::Insert
      | BamCigarOp::SoftClip => qpos = qpos + el.length
      BamCigarOp::HardClip | BamCigarOp::Pad => ()
      _ => ()
    }
  }
  None
}

///|
/// Last aligned query position (exclusive). Returns None when the alignment
/// contains no aligned bases. Matches Pysam `AlignedSegment.query_alignment_end`.
pub fn BamRecord::query_alignment_end(self : BamRecord) -> Int32? {
  let start = match self.query_alignment_start() {
    None => return None
    Some(s) => s
  }
  let mut qpos = 0
  let mut first_aligned_seen = false
  let mut last_aligned_qend = start
  for el in self.cigar {
    match el.op {
      BamCigarOp::Match
      | BamCigarOp::Equal
      | BamCigarOp::Diff => {
        qpos = qpos + el.length
        first_aligned_seen = true
        last_aligned_qend = qpos
      }
      BamCigarOp::Insert
      | BamCigarOp::SoftClip => qpos = qpos + el.length
      BamCigarOp::HardClip | BamCigarOp::Pad => ()
      _ => ()
    }
  }
  if first_aligned_seen { Some(last_aligned_qend) } else { None }
}

///|
/// Length of the aligned portion of the query (end - start). Returns 0 when
/// no aligned bases. Matches Pysam `AlignedSegment.query_alignment_length`.
pub fn BamRecord::query_alignment_length(self : BamRecord) -> Int32 {
  let start = self.query_alignment_start()
  let end = self.query_alignment_end()
  match (start, end) {
    (Some(s), Some(e)) => e - s
    _ => 0
  }
}

///|
/// Query sequence with soft clips removed. Equivalent to
/// `query_sequence()[start..end]`. Matches Pysam
/// `AlignedSegment.query_alignment_sequence`.
pub fn BamRecord::query_alignment_sequence(self : BamRecord) -> String {
  let start = self.query_alignment_start()
  let end = self.query_alignment_end()
  match (start, end) {
    (Some(s), Some(e)) =>
      if e > s && s >= 0 && e <= self.seq.length() {
        self.seq[s:e].to_owned()
      } else {
        ""
      }
    _ => ""
  }
}

///|
/// Query base qualities with soft clips removed. Matches Pysam
/// `AlignedSegment.query_alignment_qualities_str`.
pub fn BamRecord::query_alignment_qualities(self : BamRecord) -> String {
  let start = self.query_alignment_start()
  let end = self.query_alignment_end()
  match (start, end) {
    (Some(s), Some(e)) =>
      if e > s && s >= 0 && e <= self.qual.length() {
        self.qual[s:e].to_owned()
      } else {
        ""
      }
    _ => ""
  }
}

///|
/// CIGAR operation counts in order `[M, I, D, N, S, H, P, =, X]`. Missing ops
/// produce a 0 entry. Matches Pysam
/// `AlignedSegment.get_cigar_stats()[0]`.
pub fn BamRecord::get_cigar_stats(self : BamRecord) -> Array[Int32] {
  let stats = Array::make(9, 0)
  for el in self.cigar {
    match el.op {
      BamCigarOp::Match => stats[0] = stats[0] + el.length
      BamCigarOp::Insert => stats[1] = stats[1] + el.length
      BamCigarOp::Delete => stats[2] = stats[2] + el.length
      BamCigarOp::Skip => stats[3] = stats[3] + el.length
      BamCigarOp::SoftClip => stats[4] = stats[4] + el.length
      BamCigarOp::HardClip => stats[5] = stats[5] + el.length
      BamCigarOp::Pad => stats[6] = stats[6] + el.length
      BamCigarOp::Equal => stats[7] = stats[7] + el.length
      BamCigarOp::Diff => stats[8] = stats[8] + el.length
    }
  }
  stats
}

///|
/// Contiguous reference-coordinate blocks covered by aligned bases (M/=/X),
/// split at D/N. Each block is `(ref_start, ref_end)` half-open. Soft clips
/// and insertions do not produce blocks. Matches Pysam
/// `AlignedSegment.blocks` / `get_blocks()`.
pub fn BamRecord::blocks(self : BamRecord) -> Array[(Int32, Int32)] {
  let out : Array[(Int32, Int32)] = Array(capacity=0)
  if self.cigar.length() == 0 {
    return out
  }
  let mut rpos = self.pos
  let mut block_start : Int32? = None
  let mut block_end : Int32? = None
  for el in self.cigar {
    match el.op {
      BamCigarOp::Match
      | BamCigarOp::Equal
      | BamCigarOp::Diff => {
        if block_start is None {
          block_start = Some(rpos)
          block_end = Some(rpos + el.length)
        } else {
          block_end = Some(block_end.unwrap() + el.length)
        }
        rpos = rpos + el.length
      }
      BamCigarOp::Delete
      | BamCigarOp::Skip => {
        if block_start is Some(_) {
          out.push((block_start.unwrap(), block_end.unwrap()))
          block_start = None
          block_end = None
        }
        rpos = rpos + el.length
      }
      _ => ()
    }
  }
  if block_start is Some(_) {
    out.push((block_start.unwrap(), block_end.unwrap()))
  }
  out
}

///|
/// Reference sequence covered by this alignment given the external contig
/// FASTA `ref_seq`. Returns `ref_seq[reference_start:reference_end]`. Caller
/// must ensure `ref_seq` covers the full reference contig. Matches Pysam
/// `AlignedSegment.get_reference_sequence()`.
pub fn BamRecord::get_reference_sequence(self : BamRecord, ref_seq : String) -> String {
  let s = self.reference_start()
  let e = self.reference_end()
  if s < 0 || e < s || e > ref_seq.length() {
    ""
  } else {
    ref_seq[s:e].to_owned()
  }
}

///|
/// Query sequence in original sequencing orientation (reverse-complemented
/// when the read maps to the minus strand). Matches Pysam
/// `AlignedSegment.get_forward_sequence()`.
pub fn BamRecord::get_forward_sequence(self : BamRecord) -> String {
  if self.is_reverse() {
    reverse_complement(self.seq)
  } else {
    self.seq
  }
}

///|
/// Query qualities in original sequencing orientation (reversed when minus
/// strand). Matches Pysam `AlignedSegment.get_forward_qualities()`.
pub fn BamRecord::get_forward_qualities(self : BamRecord) -> String {
  if self.is_reverse() {
    let n = self.qual.length()
    let buf : FixedArray[UInt16] = FixedArray::make(n, 0)
    let mut i = 0
    while i < n {
      buf[i] = self.qual.unsafe_get(n - 1 - i)
      i = i + 1
    }
    unsafe_from_uint16_fixedarray(buf)
  } else {
    self.qual
  }
}

///|
/// Check whether a tag is present on this record. Matches Pysam
/// `AlignedSegment.has_tag()`.
pub fn BamRecord::has_tag(self : BamRecord, tag : String) -> Bool {
  for t in self.tags {
    if t.name == tag {
      return true
    }
  }
  false
}

///|
/// Create a new BamRecord with `tag` set to `value`. If the tag already exists,
/// it is updated; otherwise appended. Existing `type_code` is preserved for an
/// update. Matches Pysam `AlignedSegment.set_tag()` (in-place mutation on
/// mutable records; here we follow the project's immutable-record convention).
pub fn BamRecord::set_tag(
  self : BamRecord,
  tag : String,
  value : BamTagValue,
) -> BamRecord {
  let new_tags = self.tags.copy()
  let mut found = false
  let mut i = 0
  while i < new_tags.length() {
    if new_tags[i].name == tag {
      let old_type = new_tags[i].type_code
      new_tags[i] = BamTag::{ name: tag, type_code: old_type, value, }
      found = true
      break
    }
    i = i + 1
  }
  if !found {
    // Default type codes per Pysam convention
    let type_code : Char = match value {
      BamTagValue::Char(_) => 'A'
      BamTagValue::Int(_) => 'i'
      BamTagValue::Float(_) => 'f'
      BamTagValue::Str(_) => 'Z'
      BamTagValue::IntArray(_) => 'B'
      BamTagValue::FloatArray(_) => 'B'
    }
    new_tags.push(BamTag::{ name: tag, type_code, value, })
  }
  BamRecord::{
    ref_id: self.ref_id,
    pos: self.pos,
    bin: self.bin,
    mapq: self.mapq,
    flag: self.flag,
    l_read_name: self.l_read_name,
    cigar_len: self.cigar_len,
    l_seq: self.l_seq,
    next_ref_id: self.next_ref_id,
    next_pos: self.next_pos,
    tlen: self.tlen,
    qname: self.qname,
    cigar: self.cigar,
    seq: self.seq,
    qual: self.qual,
    tags: new_tags,
  }
}

// ============================================================================
// Immutable-style BamRecord builder methods (for fixmate / merge workflows)
// ============================================================================

///|
/// Produce a new BamRecord with `next_ref_id` and `next_pos` updated.
/// Used by `samtools merge` + `samtools fixmate` pair repair.
pub fn BamRecord::with_mate(self : BamRecord, next_ref_id : Int32, next_pos : Int32) -> BamRecord {
  BamRecord::{
    ref_id: self.ref_id,
    pos: self.pos,
    bin: self.bin,
    mapq: self.mapq,
    flag: self.flag,
    l_read_name: self.l_read_name,
    cigar_len: self.cigar_len,
    l_seq: self.l_seq,
    next_ref_id,
    next_pos,
    tlen: self.tlen,
    qname: self.qname,
    cigar: self.cigar,
    seq: self.seq,
    qual: self.qual,
    tags: self.tags.copy(),
  }
}

///|
/// Produce a new BamRecord with `tlen` (insert size) updated.
pub fn BamRecord::with_tlen(self : BamRecord, tlen : Int32) -> BamRecord {
  BamRecord::{
    ref_id: self.ref_id,
    pos: self.pos,
    bin: self.bin,
    mapq: self.mapq,
    flag: self.flag,
    l_read_name: self.l_read_name,
    cigar_len: self.cigar_len,
    l_seq: self.l_seq,
    next_ref_id: self.next_ref_id,
    next_pos: self.next_pos,
    tlen,
    qname: self.qname,
    cigar: self.cigar,
    seq: self.seq,
    qual: self.qual,
    tags: self.tags.copy(),
  }
}

///|
/// Produce a new BamRecord with the FLAG field replaced.
pub fn BamRecord::with_flag(self : BamRecord, flag : Int32) -> BamRecord {
  BamRecord::{
    ref_id: self.ref_id,
    pos: self.pos,
    bin: self.bin,
    mapq: self.mapq,
    flag,
    l_read_name: self.l_read_name,
    cigar_len: self.cigar_len,
    l_seq: self.l_seq,
    next_ref_id: self.next_ref_id,
    next_pos: self.next_pos,
    tlen: self.tlen,
    qname: self.qname,
    cigar: self.cigar,
    seq: self.seq,
    qual: self.qual,
    tags: self.tags.copy(),
  }
}

///|
/// Per-base coverage split by base on `[start, stop)`. Returns four equal-length
/// arrays `(A, C, G, T)`. `quality_threshold` filters bases whose Phred score
/// is below the cutoff (default 0 = no filter). `read_filter` supports
/// `"nofilter"` (count every read) and `"all"` (skip UNMAP / SECONDARY /
/// QCFAIL / DUP). Matches Pysam `AlignmentFile.count_coverage`.
pub fn BamFile::count_coverage(
  self : BamFile,
  contig : String,
  start : Int32,
  stop : Int32,
  quality_threshold~ : Int32 = 0,
  read_filter~ : String = "nofilter",
) -> (Array[Int32], Array[Int32], Array[Int32], Array[Int32]) {
  let n = stop - start
  if n <= 0 {
    return ([], [], [], [])
  }
  let mut ref_id = -1
  for j = 0; j < self.header.references.length(); j = j + 1 {
    if self.header.references[j].name == contig {
      ref_id = j
      break
    }
  }
  if ref_id < 0 {
    return (Array::make(n, 0), Array::make(n, 0), Array::make(n, 0), Array::make(n, 0))
  }
  let a = Array::make(n, 0)
  let c = Array::make(n, 0)
  let g = Array::make(n, 0)
  let t = Array::make(n, 0)
  let columns = self.pileup(contig, start, stop)
  for col in columns {
    let idx = col.reference_pos - start
    if idx < 0 || idx >= n {
      continue
    }
    for pu in col.pileups {
      let read = self.records[pu.record_index]
      // read filter
      if read_filter == "all" {
        if read.is_unmapped() ||
          read.is_secondary() ||
          read.is_qcfail() ||
          read.is_duplicate() {
          continue
        }
      }
      // skip deletions / refskips
      let qp = match pu.query_position {
        Some(q) => q
        None => continue
      }
      if qp < 0 || qp >= read.seq.length() {
        continue
      }
      // quality filter (Phred score = qual - 33)
      if quality_threshold > 0 {
        let phred = read.qual.unsafe_get(qp).to_int() - 33
        if phred < quality_threshold {
          continue
        }
      }
      let base = read.seq.unsafe_get(qp)
      let upper = if base >= 97 && base <= 122 { base - 32 } else { base }
      if upper == 65 {
        a[idx] = a[idx] + 1
      } else if upper == 67 {
        c[idx] = c[idx] + 1
      } else if upper == 71 {
        g[idx] = g[idx] + 1
      } else if upper == 84 {
        t[idx] = t[idx] + 1
      }
    }
  }
  (a, c, g, t)
}

///|
/// Number of reads overlapping `[start, stop)` on `contig`. Matches Pysam
/// `AlignmentFile.count`.
pub fn BamFile::count(
  self : BamFile,
  contig : String,
  start : Int32,
  stop : Int32,
) -> Int32 {
  let mut ref_id = -1
  for j = 0; j < self.header.references.length(); j = j + 1 {
    if self.header.references[j].name == contig {
      ref_id = j
      break
    }
  }
  if ref_id < 0 {
    return 0
  }
  let mut total = 0
  for r in self.records {
    if r.ref_id == ref_id && !r.is_unmapped() &&
      r.pos < stop && r.get_end() + 1 > start {
      total = total + 1
    }
  }
  total
}

///|
/// Intronic sites inferred from CIGAR `N` operations within `[start, stop)`.
/// Returns a map of `(intron_start, intron_end)` to the number of reads
/// supporting each intron. Matches Pysam `AlignmentFile.find_introns`.
pub fn BamFile::find_introns(
  self : BamFile,
  contig : String,
  start : Int32,
  stop : Int32,
) -> Map[(Int32, Int32), Int32] {
  let introns : Map[(Int32, Int32), Int32] = Map([], capacity=0)
  let mut ref_id = -1
  for j = 0; j < self.header.references.length(); j = j + 1 {
    if self.header.references[j].name == contig {
      ref_id = j
      break
    }
  }
  if ref_id < 0 {
    return introns
  }
  for r in self.records {
    if r.ref_id != ref_id {
      continue
    }
    let mut rpos = r.pos
    for el in r.cigar {
      match el.op {
        BamCigarOp::Match
        | BamCigarOp::Equal
        | BamCigarOp::Diff => rpos = rpos + el.length
        BamCigarOp::Delete => rpos = rpos + el.length
        BamCigarOp::Skip => {
          let intron_start = rpos
          let intron_end = rpos + el.length
          // overlap check
          if intron_end > start && intron_start < stop {
            let key = (intron_start, intron_end)
            let cur = match introns.get(key) {
              Some(v) => v
              None => 0
            }
            introns.set(key, cur + 1)
          }
          rpos = rpos + el.length
        }
        _ => ()
      }
    }
  }
  introns
}

// ============================================================================
// P1-1: PileupColumn / PileupRead 扩展（零 struct 变更；Pysam 对齐）
// ============================================================================

///|
/// Return the full `BamRecord` this `PileupRead` refers to, looked up via the
/// given `BamFile`. Matches Pysam `PileupRead.alignment`.
pub fn PileupRead::alignment(self : PileupRead, bam : BamFile) -> BamRecord {
  bam.records[self.record_index]
}

///|
/// When `query_position` is None (deletion / refskip), return the query
/// position immediately following the skipped region; otherwise return the
/// stored value. Matches Pysam `PileupRead.query_position_or_next`.
pub fn PileupRead::query_position_or_next(
  self : PileupRead,
  bam : BamFile,
) -> Int32? {
  match self.query_position {
    Some(q) => Some(q)
    None => {
      let r = bam.records[self.record_index]
      // Walk the CIGAR from the pileup's reference pos to find the next
      // query-consuming op (M/I/S/=/X) that follows this D/N region.
      let p = self.alignment(bam).pos + 0 // reference pos of this column
      let mut rpos = r.pos
      let mut qpos = 0
      for el in r.cigar {
        match el.op {
          BamCigarOp::Match | BamCigarOp::Equal | BamCigarOp::Diff => {
            // first aligned base strictly after the skipped ref position
            // belongs to this op if rpos > p
            if rpos + el.length > p && rpos >= p {
              return Some(qpos)
            }
            if rpos + el.length > p {
              // D covers the p position; next aligned base is at qpos
              return Some(qpos)
            }
            qpos = qpos + el.length
            rpos = rpos + el.length
          }
          BamCigarOp::Insert | BamCigarOp::SoftClip => qpos = qpos + el.length
          BamCigarOp::Delete | BamCigarOp::Skip => rpos = rpos + el.length
          _ => ()
        }
      }
      None
    }
  }
}

///|
/// Query sequences of all aligned reads in this column. Deletions / ref-skips
/// are skipped. Matches Pysam `PileupColumn.get_query_sequences()`.
pub fn BamFile::pileup_get_query_sequences(
  self : BamFile,
  col : PileupColumn,
) -> Array[String] {
  let out : Array[String] = Array(capacity=0)
  for pu in col.pileups {
    let r = self.records[pu.record_index]
    match pu.query_position {
      Some(q) =>
        if q >= 0 && q < r.seq.length() {
          out.push(r.seq[q:q + 1].to_owned())
        } else {
          out.push("")
        }
      None => ()
    }
  }
  out
}

///|
/// Query base qualities (Phred+33 char) of all aligned reads in this column.
/// Deletions / ref-skips are skipped. Matches Pysam
/// `PileupColumn.get_query_qualities()`.
pub fn BamFile::pileup_get_query_qualities(
  self : BamFile,
  col : PileupColumn,
) -> Array[String] {
  let out : Array[String] = Array(capacity=0)
  for pu in col.pileups {
    let r = self.records[pu.record_index]
    match pu.query_position {
      Some(q) =>
        if q >= 0 && q < r.qual.length() {
          out.push(r.qual[q:q + 1].to_owned())
        } else {
          out.push("")
        }
      None => ()
    }
  }
  out
}

///|
/// Query positions of all aligned reads in this column. Matches Pysam
/// `PileupColumn.get_query_positions()`.
pub fn BamFile::pileup_get_query_positions(
  self : BamFile,
  col : PileupColumn,
) -> Array[Int32] {
  let out : Array[Int32] = Array(capacity=0)
  for pu in col.pileups {
    match pu.query_position {
      Some(q) => out.push(q)
      None => ()
    }
  }
  out
}

///|
/// Mapping quality (`mapq`) of every read covering this column, in the same
/// order as `col.pileups`. Matches Pysam `PileupColumn.get_mapping_qualities()`.
pub fn BamFile::pileup_get_mapping_qualities(
  self : BamFile,
  col : PileupColumn,
) -> Array[Int32] {
  let out : Array[Int32] = Array(capacity=0)
  for pu in col.pileups {
    out.push(self.records[pu.record_index].mapq)
  }
  out
}

///|
/// Number of aligned (non-deletion, non-ref-skip) reads in this column.
/// Matches Pysam `PileupColumn.get_num_aligned()`.
pub fn BamFile::pileup_get_num_aligned(
  self : BamFile,
  col : PileupColumn,
) -> Int32 {
  let mut n = 0
  for pu in col.pileups {
    if pu.query_position is Some(_) {
      n = n + 1
    }
  }
  n
}

///|
/// Reads covering this column whose base quality (Phred score, `qual - 33`) is
/// at least `min_base_quality`. Returns `(filtered column, original count)`.
/// Deletions / ref-skips are filtered out. Matches Pysam
/// `PileupColumn.set_min_base_quality(min_base_quality)` filtering semantics.
pub fn BamFile::pileup_set_min_base_quality(
  self : BamFile,
  col : PileupColumn,
  min_base_quality : Int32,
) -> PileupColumn {
  let filtered : Array[PileupRead] = Array(capacity=0)
  for pu in col.pileups {
    match pu.query_position {
      None => ()
      Some(q) => {
        let r = self.records[pu.record_index]
        if q >= 0 && q < r.qual.length() {
          let phred = r.qual.unsafe_get(q).to_int() - 33
          if phred >= min_base_quality {
            filtered.push(pu)
          }
        }
      }
    }
  }
  PileupColumn::{
    reference_pos: col.reference_pos,
    nsegments: filtered.length(),
    pileups: filtered,
  }
}

// ============================================================================
// P1-2: IndexedReads — 按 read name 建索引，O(1) 快速查找
// ============================================================================

///|
/// Build a read-name → record-array index over a BAM file. Matches Pysam
/// `IndexedReads` (`build` + `find`).
pub fn bam_build_indexed_reads(
  bam : BamFile,
) -> Map[String, Array[BamRecord]] {
  let idx : Map[String, Array[BamRecord]] = Map([], capacity=0)
  for r in bam.records {
    let arr = match idx.get(r.qname) {
      Some(existing) => {
        let copy = existing.copy()
        copy.push(r)
        copy
      }
      None => [r]
    }
    idx.set(r.qname, arr)
  }
  idx
}

///|
/// Find all reads (typically read1 + read2) sharing the same `qname` in
/// `index`. Returns an empty array when not found.
pub fn bam_indexed_find(
  index : Map[String, Array[BamRecord]],
  qname : String,
) -> Array[BamRecord] {
  match index.get(qname) {
    Some(arr) => arr.copy()
    None => []
  }
}

// ============================================================================
// P1-2: parse_region — UCSC/samtools 风格区域字符串解析
// ============================================================================

///|
/// Parse a samtools-style region string into `(contig, start, stop)` (0-based
/// half-open). Supported forms:
///
/// - `chr1`                       → `("chr1", 0, None)`
/// - `chr1:100`                   → `("chr1", 99, 100)`
/// - `chr1:100-200`               → `("chr1", 99, 200)`
/// - `chr1:100..200`              → `("chr1", 99, 200)`
/// - `chr1:100+50`                → `("chr1", 99, 149)`   (length, 1-based)
///
/// None values mean "to the contig boundary". Matches Pysam
/// `HTSFile.parse_region`.
pub fn parse_region(region : String) -> (String, Int32?, Int32?) {
  let colon = bam_index_of_char(region, ':'.to_int())
  if colon < 0 {
    return (region, None, None)
  }
  let contig = region[0:colon].to_owned()
  let interval = region[colon + 1:region.length()].to_owned()
  let dash = bam_index_of_char(interval, '-'.to_int())
  let dotdot = bam_index_of_char(interval, ':'.to_int()) // for chr1:100..200
  let plus = bam_index_of_char(interval, '+'.to_int())
  if dash >= 0 {
    let s = parse_int(interval[0:dash].to_owned()) - 1
    let e = parse_int(interval[dash + 1:interval.length()].to_owned())
    (contig, Some(s), Some(e))
  } else if dotdot >= 0 {
    let s = parse_int(interval[0:dotdot].to_owned()) - 1
    let e = parse_int(interval[dotdot + 1:interval.length()].to_owned())
    (contig, Some(s), Some(e))
  } else if plus >= 0 {
    let s = parse_int(interval[0:plus].to_owned()) - 1
    let len = parse_int(interval[plus + 1:interval.length()].to_owned())
    (contig, Some(s), Some(s + len))
  } else {
    let s = parse_int(interval) - 1
    (contig, Some(s), Some(s + 1))
  }
}

///|
fn bam_index_of_char(s : String, target : Int32) -> Int32 {
  let mut i = 0
  while i < s.length() {
    if s.unsafe_get(i).to_int() == target {
      return i
    }
    i = i + 1
  }
  -1
}

///|
/// Fetch by samtools-style region string. Convenience wrapper combining
/// `parse_region` + `BamFile::fetch`.
pub fn BamFile::fetch_region(self : BamFile, region : String) -> Array[BamRecord] {
  let (contig, start, stop) = parse_region(region)
  match (start, stop) {
    (Some(s), Some(e)) => bam_fetch(self, contig, s, e)
    (Some(s), None) => bam_fetch(self, contig, s, self.get_reference_length_or_zero(contig))
    (None, None) => {
      let out : Array[BamRecord] = Array(capacity=0)
      for r in self.records {
        if r.ref_id >= 0 {
          out.push(r)
        }
      }
      out
    }
    _ => bam_fetch(self, contig, 0, self.get_reference_length_or_zero(contig))
  }
}

///|
fn BamFile::get_reference_length_or_zero(self : BamFile, name : String) -> Int32 {
  match self.get_reference_length(name) {
    Some(l) => l
    None => 0
  }
}

// ============================================================================
// P2-3c: BamFile::build_index — samtools index 的内存实现
// ============================================================================

///|
/// Compute which BAI bins a given half-open `[beg, end)` interval overlaps.
/// Identical to `bai_reg2bins` in bam_index.mbt but exposed here for record
/// aggregation.
fn bam_bins_for_range(beg : Int32, end : Int32) -> Array[Int32] {
  let bins : Array[Int32] = []
  if beg >= end {
    return bins
  }
  let e = end - 1
  bins.push(0)
  let mut k = 1 + (beg >> 26)
  while k <= 1 + (e >> 26) {
    bins.push(k)
    k = k + 1
  }
  k = 9 + (beg >> 23)
  while k <= 9 + (e >> 23) {
    bins.push(k)
    k = k + 1
  }
  k = 73 + (beg >> 20)
  while k <= 73 + (e >> 20) {
    bins.push(k)
    k = k + 1
  }
  k = 585 + (beg >> 17)
  while k <= 585 + (e >> 17) {
    bins.push(k)
    k = k + 1
  }
  k = 4681 + (beg >> 14)
  while k <= 4681 + (e >> 14) {
    bins.push(k)
    k = k + 1
  }
  bins
}

///|
/// Build an in-memory BAI index from the BamFile's records. This mirrors
/// `samtools index` but works entirely on parsed records (no file I/O).
///
/// Since BamFile does not retain raw file offsets, we synthesize virtual
/// offsets linearly per record — sufficient for binning correctness and
/// query-result verification, but not for actual BGZF seek. The produced
/// index can be serialized with `serialize_bai` and parsed back by
/// `parse_bai` (round-trip testable).
pub fn BamFile::build_index(self : BamFile) -> BaiIndex {
  let n_ref = self.header.references.length()
  let references : Array[BaiReference] = Array(capacity=n_ref)
  let n_no_coor : Int64 = 0L
  // Per-ref accumulators
  let ref_bins : Array[Map[Int32, Array[BaiChunk]]] = Array(capacity=n_ref)
  let ref_linear_last : Array[Int64] = Array(capacity=n_ref)
  let ref_linear_count : Array[Int32] = Array(capacity=n_ref)
  // Init
  let mut i = 0
  while i < n_ref {
    ref_bins.push(Map([], capacity=0))
    ref_linear_last.push(0L)
    ref_linear_count.push(0)
    i = i + 1
  }
  // Linear "virtual offset" synthesizer
  let mut next_voff : Int64 = 65536L // reserve 1 virtual offset for header
  for r in self.records {
    if r.is_unmapped() || r.ref_id < 0 {
      continue
    }
    let beg = r.pos
    let end = beg + r.get_end()
    if end <= beg {
      continue
    }
    let ref_ = r.ref_id
    if ref_ >= n_ref {
      continue
    }
    let voff = next_voff
    let voff_end = voff + 256L // synthetic chunk end (256 bytes per record)
    next_voff = voff_end
    let chunk = BaiChunk::{ chunk_beg: voff, chunk_end: voff_end, }
    // Bins
    let bins = bam_bins_for_range(beg, end)
    for bin in bins {
      let cur = ref_bins[ref_].get(bin)
      match cur {
        Some(arr) => {
          let copy = arr.copy()
          copy.push(chunk)
          ref_bins[ref_].set(bin, copy)
        }
        None => {
          ref_bins[ref_].set(bin, [chunk])
        }
      }
    }
    // Linear index: every 16384 bp, record minimum voff
    let iv = beg >> 14
    // Expand ref_linear_last to fit iv+1
    while ref_linear_last.length() <= iv {
      ref_linear_last.push(9223372036854775807L) // sentinel "no record"
      ref_linear_count.push(0)
    }
    if voff < ref_linear_last[iv] {
      ref_linear_last[iv] = voff
    }
  }
  // Collapse sentinels back to the last known good offset
  i = 0
  while i < n_ref {
    let linear : Array[Int64] = []
    let ref_len = if i < self.header.references.length() {
      self.header.references[i].length
    } else {
      0
    }
    let n_ivs = (ref_len >> 14) + 1
    let mut last_good : Int64 = 0L
    let mut k = 0
    while k < n_ivs {
      if k < ref_linear_last.length() {
        if ref_linear_last[k] != 9223372036854775807L {
          last_good = ref_linear_last[k]
        }
      }
      linear.push(last_good)
      k = k + 1
    }
    references.push(BaiReference::{ bins: ref_bins[i], linear, })
    i = i + 1
  }
  BaiIndex::{
    n_ref,
    references,
    n_no_coor: Some(n_no_coor),
  }
}
