///|
/// CRAM format — ITF-8/LTF-8 variable-length integer encoding, binary write
/// (`CramFile::to_cram_bytes`), and SAM text conversion.
///
/// The parse side lives in `bam.mbt` (`parse_cram`, `CramRecord`, `CramHeader`,
/// `CramFile`). This module adds the symmetric write + SAM export + proper
/// CRAM 3.0 variable-length integer codecs.

// ============================================================
// ITF-8 — CRAM 3.0 variable-length 32-bit integer encoding
// ============================================================

///|
/// Encode a 32-bit integer as ITF-8 (1–5 bytes).
///
/// ITF-8 layout:
/// - `0xxxxxxx`                         → 1 byte  (0 .. 127)
/// - `10xxxxxx xxxxxxxx`                → 2 bytes (0 .. 16383)
/// - `110xxxxx xxxxxxxx xxxxxxxx`       → 3 bytes
/// - `1110xxxx xxxxxxxx xxxxxxxx xxxxxxxx` → 4 bytes
/// - `11110xxx xxxxxxxx xxxxxxxx xxxxxxxx xxxxxxxx` → 5 bytes
pub fn itf8_encode(value : Int, buf : Array[Int]) -> Unit {
  let v = value & 0xFFFFFFFF
  if v <= 0x7F {
    buf.push(v)
  } else if v <= 0x3FFF {
    buf.push(0x80 | ((v >> 8) & 0x3F))
    buf.push(v & 0xFF)
  } else if v <= 0x1FFFFF {
    buf.push(0xC0 | ((v >> 16) & 0x1F))
    buf.push((v >> 8) & 0xFF)
    buf.push(v & 0xFF)
  } else if v <= 0xFFFFFFF {
    buf.push(0xE0 | ((v >> 24) & 0x0F))
    buf.push((v >> 16) & 0xFF)
    buf.push((v >> 8) & 0xFF)
    buf.push(v & 0xFF)
  } else {
    buf.push(0xF0 | ((v >> 28) & 0x07))
    buf.push((v >> 20) & 0xFF)
    buf.push((v >> 12) & 0xFF)
    buf.push((v >> 4) & 0xFF)
    buf.push(v & 0x0F)
  }
}

///|
/// Decode an ITF-8 integer from a string at the given position.
/// Returns `(value, new_position)` or `(-1, pos)` on error.
pub fn itf8_decode(data : String, pos : Int) -> (Int, Int) {
  if pos >= data.length() {
    return (-1, pos)
  }
  let b0 = data.unsafe_get(pos).to_int() & 0xFF
  if b0 <= 0x7F {
    (b0, pos + 1)
  } else if (b0 & 0xC0) == 0x80 {
    if pos + 1 >= data.length() {
      return (-1, pos)
    }
    let b1 = data.unsafe_get(pos + 1).to_int() & 0xFF
    (((b0 & 0x3F) << 8) | b1, pos + 2)
  } else if (b0 & 0xE0) == 0xC0 {
    if pos + 2 >= data.length() {
      return (-1, pos)
    }
    let b1 = data.unsafe_get(pos + 1).to_int() & 0xFF
    let b2 = data.unsafe_get(pos + 2).to_int() & 0xFF
    (((b0 & 0x1F) << 16) | (b1 << 8) | b2, pos + 3)
  } else if (b0 & 0xF0) == 0xE0 {
    if pos + 3 >= data.length() {
      return (-1, pos)
    }
    let b1 = data.unsafe_get(pos + 1).to_int() & 0xFF
    let b2 = data.unsafe_get(pos + 2).to_int() & 0xFF
    let b3 = data.unsafe_get(pos + 3).to_int() & 0xFF
    (((b0 & 0x0F) << 24) | (b1 << 16) | (b2 << 8) | b3, pos + 4)
  } else {
    // 0xF0 prefix — 5 bytes
    if pos + 4 >= data.length() {
      return (-1, pos)
    }
    let b1 = data.unsafe_get(pos + 1).to_int() & 0xFF
    let b2 = data.unsafe_get(pos + 2).to_int() & 0xFF
    let b3 = data.unsafe_get(pos + 3).to_int() & 0xFF
    let b4 = data.unsafe_get(pos + 4).to_int() & 0xFF
    (
      ((b0 & 0x07) << 28) | (b1 << 20) | (b2 << 12) | (b3 << 4) | (b4 & 0x0F),
      pos + 5,
    )
  }
}

// ============================================================
// LTF-8 — CRAM 3.0 variable-length 64-bit integer encoding
// ============================================================

///|
/// Encode a 64-bit integer as LTF-8 (1–9 bytes).
/// Uses the same high-bit-prefix scheme as ITF-8 but extends to 9 bytes.
pub fn ltf8_encode(value : Int, buf : Array[Int]) -> Unit {
  // For MoonBit Int (32-bit on this platform), we only need up to 5 bytes.
  // The full 9-byte LTF-8 is for 64-bit values beyond Int range.
  itf8_encode(value, buf)
}

///|
/// Decode an LTF-8 integer. Same as ITF-8 for 32-bit range.
pub fn ltf8_decode(data : String, pos : Int) -> (Int, Int) {
  itf8_decode(data, pos)
}

// ============================================================
// CRAM binary write — inverse of parse_cram
// ============================================================

///|
fn cram_write_u32_le(buf : Array[Int], v : Int) -> Unit {
  buf.push(v & 0xFF)
  buf.push((v >> 8) & 0xFF)
  buf.push((v >> 16) & 0xFF)
  buf.push((v >> 24) & 0xFF)
}

///|
fn cram_write_cstring(buf : Array[Int], s : String) -> Unit {
  for i = 0; i < s.length(); i = i + 1 {
    buf.push(s.unsafe_get(i).to_int() & 0xFF)
  }
  buf.push(0) // null terminator
}

///|
fn cram_cigar_to_code(op : BamCigarOp, length : Int) -> Int {
  let op_code = match op {
    BamCigarOp::Match => 0
    BamCigarOp::Insert => 1
    BamCigarOp::Delete => 2
    BamCigarOp::Skip => 3
    BamCigarOp::SoftClip => 4
    BamCigarOp::HardClip => 5
    BamCigarOp::Pad => 6
    BamCigarOp::Equal => 7
    BamCigarOp::Diff => 8
  }
  (length << 4) | op_code
}

///|
/// Serialize a `CramReadFeatures` to a byte array (matching the parse layout).
fn cram_write_read_features(buf : Array[Int], rf : CramReadFeatures) -> Unit {
  cram_write_u32_le(buf, rf.read_group)
  cram_write_u32_le(buf, rf.flags)
  buf.push(rf.mapping_quality & 0xFF) // int8
  // read name: length-prefixed raw bytes (no null terminator)
  cram_write_u32_le(buf, rf.read_name.length())
  for i = 0; i < rf.read_name.length(); i = i + 1 {
    buf.push(rf.read_name.unsafe_get(i).to_int() & 0xFF)
  }
  // cigar
  cram_write_u32_le(buf, rf.cigar.length())
  for i = 0; i < rf.cigar.length(); i = i + 1 {
    let el = rf.cigar[i]
    cram_write_u32_le(buf, cram_cigar_to_code(el.op, el.length))
  }
  // sequence
  cram_write_u32_le(buf, rf.sequence.length())
  for i = 0; i < rf.sequence.length(); i = i + 1 {
    buf.push(rf.sequence.unsafe_get(i).to_int() & 0xFF)
  }
  // quality
  cram_write_u32_le(buf, rf.quality.length())
  for i = 0; i < rf.quality.length(); i = i + 1 {
    buf.push(rf.quality.unsafe_get(i).to_int() & 0xFF)
  }
  // tags
  cram_write_u32_le(buf, rf.tags.length())
  for i = 0; i < rf.tags.length(); i = i + 1 {
    cram_write_cstring(buf, rf.tags[i].key)
    cram_write_cstring(buf, rf.tags[i].value)
  }
}

///|
/// Serialize a single `CramRecord` to a byte array (without the size prefix).
fn cram_write_record_body(buf : Array[Int], rec : CramRecord) -> Unit {
  cram_write_u32_le(buf, rec.ref_id)
  cram_write_u32_le(buf, rec.pos)
  cram_write_read_features(buf, rec.read_features)
}

///|
/// Serialize a `CramFile` to CRAM binary bytes (the inverse of `parse_cram`).
///
/// The output can be re-parsed by `parse_cram` to recover the same data —
/// a round-trip.
pub fn CramFile::to_cram_bytes(self : CramFile) -> String {
  let buf : Array[Int] = []

  // Magic
  buf.push(67) // C
  buf.push(82) // R
  buf.push(65) // A
  buf.push(77) // M

  // Version
  cram_write_u32_le(buf, self.header.version.major)
  cram_write_u32_le(buf, self.header.version.minor)

  // Reference count
  cram_write_u32_le(buf, self.header.reference_count)

  // Compression code
  let comp_code = match self.header.compression {
    CramCompression::None => 0
    CramCompression::Gzip => 1
    CramCompression::Bzip2 => 2
    CramCompression::Lzma => 3
    CramCompression::Unknown(code) => code
  }
  cram_write_u32_le(buf, comp_code)

  // Stream count
  cram_write_u32_le(buf, self.header.stream_count)

  // References
  for i = 0; i < self.header.references.length(); i = i + 1 {
    let ref_ = self.header.references[i]
    cram_write_cstring(buf, ref_.name)
    cram_write_u32_le(buf, ref_.length)
    cram_write_cstring(buf, ref_.md5)
    cram_write_cstring(buf, ref_.uri)
  }

  // Stream info
  for i = 0; i < self.header.stream_info.length(); i = i + 1 {
    cram_write_cstring(buf, self.header.stream_info[i])
  }

  // SAM header
  cram_write_u32_le(buf, self.header.sam_header.length())
  for i = 0; i < self.header.sam_header.length(); i = i + 1 {
    buf.push(self.header.sam_header.unsafe_get(i).to_int() & 0xFF)
  }

  // Records (each prefixed with its body size)
  for i = 0; i < self.records.length(); i = i + 1 {
    let body : Array[Int] = []
    cram_write_record_body(body, self.records[i])
    cram_write_u32_le(buf, body.length())
    for j = 0; j < body.length(); j = j + 1 {
      buf.push(body[j])
    }
  }

  // Convert byte array to String
  let sb = StringBuilder()
  for i = 0; i < buf.length(); i = i + 1 {
    sb.write_char(buf[i].to_uint16().unsafe_to_char())
  }
  sb.to_string()
}

// ============================================================
// CRAM → SAM text conversion
// ============================================================

///|
/// Convert a `CramRecord` to a SAM-format line (1-based position).
pub fn CramRecord::to_sam_line(
  self : CramRecord,
  header : CramHeader,
) -> String {
  let rname = self.get_reference_name(header)
  let flag = self.read_features.flags
  let pos_1based = self.pos + 1
  let mapq = self.read_features.mapping_quality
  // CIGAR string
  let cigar_sb = StringBuilder()
  for i = 0; i < self.read_features.cigar.length(); i = i + 1 {
    let el = self.read_features.cigar[i]
    let op_char = match el.op {
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
    cigar_sb.write_string(el.length.to_string())
    cigar_sb.write_char(op_char)
  }
  let cigar_str = cigar_sb.to_string()
  let cigar = if cigar_str.length() > 0 { cigar_str } else { "*" }
  let seq = self.read_features.sequence
  let qual = self.read_features.quality
  // Optional tags
  let tag_sb = StringBuilder()
  for i = 0; i < self.read_features.tags.length(); i = i + 1 {
    let tag = self.read_features.tags[i]
    tag_sb.write_char('\t')
    tag_sb.write_string(tag.key)
    tag_sb.write_char(':')
    tag_sb.write_string(tag.value)
  }
  "\{self.read_features.read_name}\t\{flag}\t\{rname}\t\{pos_1based}\t\{mapq}\t\{cigar}\t*\t0\t0\t\{seq}\t\{qual}\{tag_sb.to_string()}"
}

///|
/// Convert a `CramFile` to SAM text (header + records).
pub fn CramFile::to_sam(self : CramFile) -> String {
  let sb = StringBuilder()
  // SAM header
  if self.header.sam_header.length() > 0 {
    sb.write_string(self.header.sam_header)
    if !(self.header.sam_header.has_suffix("\n")) {
      sb.write_char('\n')
    }
  }
  // Records
  for i = 0; i < self.records.length(); i = i + 1 {
    sb.write_string(self.records[i].to_sam_line(self.header))
    sb.write_char('\n')
  }
  sb.to_string()
}