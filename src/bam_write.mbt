///|
/// BAM binary writer — the reverse of `parse_bam`. Encodes a `BamFile`
/// (header + records) into uncompressed BAM bytes using the htslib layout
/// documented on the read side (`bam.mbt`), and optionally BGZF-compresses via
/// `bgzf_compress`.
///
/// Note: the reader's `read_int32_le` does not sign-extend, so the signed BAM
/// fields (pos/next_pos/tlen) round-trip using that same unsigned convention.

///|
fn bam_write_u8(buf : FixedArray[UInt16], pos : Int, v : Int) -> Int {
  buf[pos] = (v & 0xFF).to_uint16()
  pos + 1
}

///|
fn bam_write_u16_le(buf : FixedArray[UInt16], pos : Int, v : Int) -> Int {
  buf[pos] = (v & 0xFF).to_uint16()
  buf[pos + 1] = ((v >> 8) & 0xFF).to_uint16()
  pos + 2
}

///|
fn bam_write_i32_le(buf : FixedArray[UInt16], pos : Int, v : Int) -> Int {
  buf[pos] = (v & 0xFF).to_uint16()
  buf[pos + 1] = ((v >> 8) & 0xFF).to_uint16()
  buf[pos + 2] = ((v >> 16) & 0xFF).to_uint16()
  buf[pos + 3] = ((v >> 24) & 0xFF).to_uint16()
  pos + 4
}

///|
let bam_base_lookup = "=ACMGRSVTWYHKDBN"

///|
fn bam_base_code(seq : String, i : Int) -> Int {
  let c = seq.unsafe_get(i)
  let code = 15
  for k = 0; k < 16; k = k + 1 {
    if bam_base_lookup.unsafe_get(k) == c {
      return k
    }
  }
  code
}

///|
fn cigar_op_code(op : BamCigarOp) -> Int {
  match op {
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
}

///|
fn bam_tag_size(tag : BamTag) -> Int {
  let v = match tag.value {
    BamTagValue::Char(_) => 1
    BamTagValue::Int(_) => match tag.type_code {
      'c' => 1
      'C' => 1
      's' => 2
      'S' => 2
      _ => 4
    }
    BamTagValue::Float(_) => 4
    BamTagValue::Str(s) => s.length() + 1
    BamTagValue::IntArray(arr) => 1 + 4 + arr.length() * 4
    BamTagValue::FloatArray(arr) => 1 + 4 + arr.length() * 4
  }
  3 + v
}

///|
/// IEEE-754 binary32 bit pattern (round-to-nearest) for a Double.
fn f32_bits(v : Double) -> Int {
  if v != v {
    return 0x7FC00000 // NaN
  }
  if v == 0.0 {
    return if (1.0 / v) < 0.0 { 0x80000000 } else { 0 }
  }
  if v > 3.4028234663852886e38 {
    return 0x7F800000 // +Inf
  }
  if v < -3.4028234663852886e38 {
    return 0xFF800000 // -Inf
  }
  let sign = if v < 0.0 { 1 } else { 0 }
  let a = if v < 0.0 { 0.0 - v } else { v }
  let mut e = 0
  let mut m = a
  if m >= 1.0 {
    while m >= 2.0 {
      m = m / 2.0
      e = e + 1
    }
  } else {
    while m < 1.0 {
      m = m * 2.0
      e = e - 1
    }
  }
  let mut exp = e + 127
  if exp <= 0 {
    return 0 // subnormal: approximated as zero (kept for simplicity)
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
fn bam_write_f32(buf : FixedArray[UInt16], pos : Int, v : Double) -> Int {
  bam_write_i32_le(buf, pos, f32_bits(v))
}

///|
fn bam_write_tag(tag : BamTag, buf : FixedArray[UInt16], pos : Int) -> Int {
  buf[pos] = tag.name.unsafe_get(0)
  buf[pos + 1] = tag.name.unsafe_get(1)
  buf[pos + 2] = tag.type_code.to_int().to_uint16()
  let mut p = pos + 3
  match tag.value {
    BamTagValue::Char(c) => {
      buf[p] = c.to_int().to_uint16()
      p = p + 1
    }
    BamTagValue::Int(v) => {
      if tag.type_code == 'c' || tag.type_code == 'C' {
        p = bam_write_u8(buf, p, v)
      } else if tag.type_code == 's' || tag.type_code == 'S' {
        p = bam_write_u16_le(buf, p, v)
      } else {
        p = bam_write_i32_le(buf, p, v)
      }
    }
    BamTagValue::Float(v) => {
      p = bam_write_f32(buf, p, v)
    }
    BamTagValue::Str(s) => {
      let k = s.length()
      let mut i = 0
      while i < k {
        buf[p + i] = s.unsafe_get(i)
        i = i + 1
      }
      buf[p + k] = 0
      p = p + k + 1
    }
    BamTagValue::IntArray(arr) => {
      buf[p] = ('i'.to_int()).to_uint16()
      p = p + 1
      p = bam_write_i32_le(buf, p, arr.length())
      for x in arr {
        p = bam_write_i32_le(buf, p, x)
      }
    }
    BamTagValue::FloatArray(arr) => {
      buf[p] = ('f'.to_int()).to_uint16()
      p = p + 1
      p = bam_write_i32_le(buf, p, arr.length())
      for x in arr {
        p = bam_write_f32(buf, p, x)
      }
    }
  }
  p
}

///|
fn bam_write_seq(
  seq : String,
  l_seq : Int,
  buf : FixedArray[UInt16],
  pos : Int,
) -> Int {
  let seq_bytes = (l_seq + 1) / 2
  let mut i = 0
  while i < l_seq {
    let code = bam_base_code(seq, i)
    let byte_pos = pos + i / 2
    let cur = buf[byte_pos].to_int()
    if i % 2 == 0 {
      buf[byte_pos] = (code << 4).to_uint16()
    } else {
      buf[byte_pos] = (cur | code).to_uint16()
    }
    i = i + 1
  }
  pos + seq_bytes
}

///|
fn encode_bam_record(rec : BamRecord) -> String {
  let l_read_name = rec.qname.length() + 1
  let cigar_len = rec.cigar.length()
  let l_seq = rec.l_seq
  let seq_bytes = (l_seq + 1) / 2
  let mut tags_size = 0
  for t in rec.tags {
    tags_size = tags_size + bam_tag_size(t)
  }
  let body = 32 + l_read_name + cigar_len * 4 + seq_bytes + l_seq + tags_size
  let total = 4 + body
  let buf : FixedArray[UInt16] = FixedArray::make(total, 0)

  let mut p = bam_write_i32_le(buf, 0, body)
  p = bam_write_i32_le(buf, p, rec.ref_id)
  p = bam_write_i32_le(buf, p, rec.pos)
  p = bam_write_i32_le(buf, p, (rec.bin << 16) | (rec.mapq << 8) | l_read_name)
  p = bam_write_i32_le(buf, p, (rec.flag << 16) | cigar_len)
  p = bam_write_i32_le(buf, p, l_seq)
  p = bam_write_i32_le(buf, p, rec.next_ref_id)
  p = bam_write_i32_le(buf, p, rec.next_pos)
  p = bam_write_i32_le(buf, p, rec.tlen)

  // qname + null terminator
  let mut i = 0
  while i < l_read_name - 1 {
    buf[p + i] = rec.qname.unsafe_get(i)
    i = i + 1
  }
  buf[p + l_read_name - 1] = 0
  p = p + l_read_name

  // CIGAR
  for elem in rec.cigar {
    p = bam_write_i32_le(buf, p, (elem.length << 4) | cigar_op_code(elem.op))
  }

  // SEQ
  p = bam_write_seq(rec.seq, l_seq, buf, p)

  // QUAL
  i = 0
  while i < l_seq {
    buf[p + i] = (rec.qual.unsafe_get(i).to_int() - 33).to_uint16()
    i = i + 1
  }
  p = p + l_seq

  // AUX tags
  for t in rec.tags {
    p = bam_write_tag(t, buf, p)
  }

  unsafe_from_uint16_fixedarray(buf)
}

///|
fn encode_bam_header(header : BamHeader) -> String {
  let l_text = header.text.length()
  let mut ref_bytes = 0
  for i = 0; i < header.references.length(); i = i + 1 {
    ref_bytes = ref_bytes + 4 + header.references[i].name.length() + 1 + 4
  }
  let total = 4 + 4 + l_text + 4 + ref_bytes
  let buf : FixedArray[UInt16] = FixedArray::make(total, 0)
  buf[0] = (0x42).to_uint16()
  buf[1] = (0x41).to_uint16()
  buf[2] = (0x4D).to_uint16()
  buf[3] = (0x01).to_uint16()
  let mut p = 4
  p = bam_write_i32_le(buf, p, l_text)
  let mut i = 0
  while i < l_text {
    buf[p + i] = header.text.unsafe_get(i)
    i = i + 1
  }
  p = p + l_text
  p = bam_write_i32_le(buf, p, header.references.length())
  for j = 0; j < header.references.length(); j = j + 1 {
    let r = header.references[j]
    let name_len = r.name.length()
    p = bam_write_i32_le(buf, p, name_len + 1)
    let mut k = 0
    while k < name_len {
      buf[p + k] = r.name.unsafe_get(k)
      k = k + 1
    }
    p = p + name_len
    buf[p] = 0
    p = p + 1
    p = bam_write_i32_le(buf, p, r.length)
  }
  unsafe_from_uint16_fixedarray(buf)
}

///|
/// Serialize a `BamFile` into uncompressed BAM bytes (magic `BAM\1` + header +
/// records). The inverse of `parse_bam`.
pub fn BamFile::to_bam_bytes(self : BamFile) -> String {
  let buf = StringBuilder()
  buf.write_string(encode_bam_header(self.header))
  for rec in self.records {
    buf.write_string(encode_bam_record(rec))
  }
  buf.to_string()
}

///|
/// Serialize a `BamFile` into BGZF-compressed BAM bytes (valid `.bam` on-disk
/// representation). The inverse of `parse_bam_from_bgzf`.
pub fn BamFile::to_bgzf_bytes(self : BamFile) -> String raise BgzfError {
  bgzf_compress(self.to_bam_bytes())
}