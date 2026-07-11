///|
/// SAM format parsing module, matching pysam SAM specification.

///|
/// SAM header line types.
pub enum SamHeaderType {
  HD
  SQ
  RG
  PG
  CO
} derive(Eq)

///|
/// SAM header record.
pub struct SamHeader {
  header_type : SamHeaderType
  tags : @hashmap.HashMap[String, String]
}

///|
/// CIGAR operation types.
pub enum CigarOp {
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
pub struct CigarElement {
  op : CigarOp
  length : Int
} derive(Eq)

///|
/// SAM alignment record (11 mandatory fields + optional tags).
pub struct SamRecord {
  qname : String
  flag : Int
  rname : String
  pos : Int
  mapq : Int
  cigar : Array[CigarElement]
  rnext : String
  pnext : Int
  tlen : Int
  seq : String
  qual : String
  tags : @hashmap.HashMap[String, String]
}

///|
/// SAM file representation.
pub struct SamFile {
  headers : Array[SamHeader]
  records : Array[SamRecord]
}

///|
/// Parse SAM header type from string.
fn parse_header_type(s : String) -> SamHeaderType {
  if s == "HD" {
    return SamHeaderType::HD
  }
  if s == "SQ" {
    return SamHeaderType::SQ
  }
  if s == "RG" {
    return SamHeaderType::RG
  }
  if s == "PG" {
    return SamHeaderType::PG
  }
  if s == "CO" {
    return SamHeaderType::CO
  }
  SamHeaderType::CO
}

///|
/// Parse CIGAR operation from char.
fn parse_cigar_op(c : Char) -> CigarOp {
  if c == 'M' {
    return CigarOp::Match
  }
  if c == 'I' {
    return CigarOp::Insert
  }
  if c == 'D' {
    return CigarOp::Delete
  }
  if c == 'N' {
    return CigarOp::Skip
  }
  if c == 'S' {
    return CigarOp::SoftClip
  }
  if c == 'H' {
    return CigarOp::HardClip
  }
  if c == 'P' {
    return CigarOp::Pad
  }
  if c == '=' {
    return CigarOp::Equal
  }
  if c == 'X' {
    return CigarOp::Diff
  }
  CigarOp::Match
}

///|
/// Parse CIGAR string into CigarElement array.
pub fn parse_cigar(cigar_str : String) -> Array[CigarElement] {
  let result = Array::new()
  if cigar_str == "*" {
    return result
  }

  let mut i = 0
  let n = cigar_str.length()
  while i < n {
    let c = cigar_str.unsafe_get(i)
    let ci = c.to_int()
    if ci >= '0'.to_int() && ci <= '9'.to_int() {
      let num_str = StringBuilder::new()
      while i < n {
        let d = cigar_str.unsafe_get(i)
        let di = d.to_int()
        if di >= '0'.to_int() && di <= '9'.to_int() {
          num_str.write_char(d.unsafe_to_char())
          i = i + 1
        } else {
          break
        }
      }
      let length = parse_int(num_str.to_string())
      if i < n {
        let op_char = cigar_str.unsafe_get(i).unsafe_to_char()
        let op = parse_cigar_op(op_char)
        result.push(CigarElement::{ op, length })
        i = i + 1
      }
    } else {
      i = i + 1
    }
  }
  result
}

///|
/// Parse SAM header line.
fn parse_sam_header(line : String) -> SamHeader? {
  if line.length() < 2 {
    return None
  }
  let first_char = line.unsafe_get(0).to_int()
  if first_char != '@'.to_int() {
    return None
  }

  let parts = split_by_tab(line)
  if parts.length() == 0 {
    return None
  }

  let header_part = parts[0]
  if header_part.length() < 2 {
    return None
  }
  let header_type_str = header_part[1:header_part.length()].to_owned()
  let header_type = parse_header_type(header_type_str)

  let tags = @hashmap.HashMap([])
  let mut i = 1
  while i < parts.length() {
    let tag = parts[i]
    let colon_pos = tag.find(":")
    if colon_pos.is_some() {
      let pos = colon_pos.unwrap()
      let key = tag[0:pos].to_owned()
      let value = if pos + 1 < tag.length() {
        tag[pos + 1:tag.length()].to_owned()
      } else {
        ""
      }
      tags.set(key, value)
    }
    i = i + 1
  }

  Some(SamHeader::{ header_type, tags })
}

///|
/// Parse SAM alignment record line.
fn parse_sam_record(line : String) -> SamRecord? {
  let parts = split_by_tab(line)
  if parts.length() < 11 {
    return None
  }

  let qname = parts[0]
  let flag = parse_int(parts[1])
  let rname = parts[2]
  let pos = parse_int(parts[3])
  let mapq = parse_int(parts[4])
  let cigar = parse_cigar(parts[5])
  let rnext = parts[6]
  let pnext = parse_int(parts[7])
  let tlen = parse_int(parts[8])
  let seq = parts[9]
  let qual = parts[10]

  let tags = @hashmap.HashMap([])
  let mut i = 11
  while i < parts.length() {
    let tag = parts[i]
    let colon_pos = tag.find(":")
    if colon_pos.is_some() {
      let pos1 = colon_pos.unwrap()
      let sub = tag[pos1 + 1:tag.length()]
      let colon_pos2 = sub.find(":")
      if colon_pos2.is_some() {
        let pos2 = colon_pos2.unwrap()
        let key = tag[0:pos1].to_owned()
        let value = sub[pos2 + 1:sub.length()].to_owned()
        tags.set(key, value)
      }
    }
    i = i + 1
  }

  Some(SamRecord::{
    qname,
    flag,
    rname,
    pos,
    mapq,
    cigar,
    rnext,
    pnext,
    tlen,
    seq,
    qual,
    tags,
  })
}

///|
/// Parse SAM content from string.
pub fn parse_sam(content : String) -> SamFile {
  let headers = Array::new()
  let records = Array::new()

  let lines = content.split("\n").to_array()
  let mut i = 0
  let n = lines.length()
  while i < n {
    let line = lines[i].to_owned()
    let line_len = line.length()
    if line_len == 0 {
      i = i + 1
      continue
    }

    if line.unsafe_get(0).to_int() == '@'.to_int() {
      let header_opt = parse_sam_header(line)
      if header_opt.is_some() {
        headers.push(header_opt.unwrap())
      }
    } else {
      let record_opt = parse_sam_record(line)
      if record_opt.is_some() {
        records.push(record_opt.unwrap())
      }
    }
    i = i + 1
  }

  SamFile::{ headers, records }
}

///|
/// Check if the read is paired.
pub fn SamRecord::is_paired(self : SamRecord) -> Bool {
  (self.flag & 0x0001) != 0
}

///|
/// Check if the read is properly paired.
pub fn SamRecord::is_proper_pair(self : SamRecord) -> Bool {
  (self.flag & 0x0002) != 0
}

///|
/// Check if the read is unmapped.
pub fn SamRecord::is_unmapped(self : SamRecord) -> Bool {
  (self.flag & 0x0004) != 0
}

///|
/// Check if the mate is unmapped.
pub fn SamRecord::is_mate_unmapped(self : SamRecord) -> Bool {
  (self.flag & 0x0008) != 0
}

///|
/// Check if this is the reverse complement.
pub fn SamRecord::is_reverse(self : SamRecord) -> Bool {
  (self.flag & 0x010) != 0
}

///|
/// Check if mate is reverse complement.
pub fn SamRecord::is_mate_reverse(self : SamRecord) -> Bool {
  (self.flag & 0x020) != 0
}

///|
/// Check if this is the first read in pair.
pub fn SamRecord::is_read1(self : SamRecord) -> Bool {
  (self.flag & 0x040) != 0
}

///|
/// Check if this is the second read in pair.
pub fn SamRecord::is_read2(self : SamRecord) -> Bool {
  (self.flag & 0x080) != 0
}

///|
/// Check if this is a secondary alignment.
pub fn SamRecord::is_secondary(self : SamRecord) -> Bool {
  (self.flag & 0x100) != 0
}

///|
/// Check if the read is filtered out.
pub fn SamRecord::is_supplementary(self : SamRecord) -> Bool {
  (self.flag & 0x800) != 0
}

///|
/// Calculate the end position of the alignment.
pub fn SamRecord::get_end(self : SamRecord) -> Int {
  if self.pos == 0 || self.cigar.length() == 0 {
    return 0
  }

  let mut end = self.pos - 1
  let mut i = 0
  while i < self.cigar.length() {
    let element = self.cigar[i]
    let op = element.op
    if op == CigarOp::Match ||
      op == CigarOp::Delete ||
      op == CigarOp::Skip ||
      op == CigarOp::Equal ||
      op == CigarOp::Diff {
      end = end + element.length
    }
    i = i + 1
  }
  end
}

///|
/// Calculate the alignment length.
pub fn SamRecord::get_alignment_length(self : SamRecord) -> Int {
  let mut length = 0
  let mut i = 0
  while i < self.cigar.length() {
    let element = self.cigar[i]
    let op = element.op
    if op == CigarOp::Match ||
      op == CigarOp::Insert ||
      op == CigarOp::Equal ||
      op == CigarOp::Diff {
      length = length + element.length
    }
    i = i + 1
  }
  length
}

///|
/// Get number of records in SAM file.
pub fn SamFile::len(self : SamFile) -> Int {
  self.records.length()
}

///|
/// Get number of headers in SAM file.
pub fn SamFile::num_headers(self : SamFile) -> Int {
  self.headers.length()
}

///|
/// Get all sequence names from @SQ headers.
pub fn SamFile::get_sequence_names(self : SamFile) -> Array[String] {
  let result = Array::new()
  let mut i = 0
  while i < self.headers.length() {
    let header = self.headers[i]
    if header.header_type == SamHeaderType::SQ {
      let sn_opt = header.tags.get("SN")
      if sn_opt.is_some() {
        result.push(sn_opt.unwrap())
      }
    }
    i = i + 1
  }
  result
}

///|
/// Get sequence length from @SQ header.
pub fn SamFile::get_sequence_length(self : SamFile, name : String) -> Int? {
  let mut i = 0
  while i < self.headers.length() {
    let header = self.headers[i]
    if header.header_type == SamHeaderType::SQ {
      let sn_opt = header.tags.get("SN")
      if sn_opt.is_some() && sn_opt.unwrap() == name {
        let ln_opt = header.tags.get("LN")
        if ln_opt.is_some() {
          return Some(parse_int(ln_opt.unwrap()))
        }
      }
    }
    i = i + 1
  }
  None
}
