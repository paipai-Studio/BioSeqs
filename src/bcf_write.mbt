///|
/// BCF 2.2 binary writer — the reverse of `parse_bcf`. Encodes a `BcfFile`
/// (header text + records) into raw BCF bytes using the htslib layout
/// documented on the read side (`bcf.mbt`). Byte helpers (`bam_write_*`,
/// `f32_bits`) are shared with the BAM writer.

///|
/// Minimal int sub-type (1=i8, 2=i16, 3=i32) able to hold every value in `vs`.
fn bcf_int_type(vs : Array[Int]) -> Int {
  let mut m = 0
  for v in vs {
    if v > m {
      m = v
    }
  }
  if m < 128 {
    1
  } else if m < 32768 {
    2
  } else {
    3
  }
}

///|
fn bcf_tint_len(v : Int) -> Int {
  if v < 128 {
    2
  } else if v < 32768 {
    3
  } else {
    5
  }
}

///|
fn bcf_size_head_len(count : Int) -> Int {
  if count < 15 {
    1
  } else {
    1 + bcf_tint_len(count)
  }
}

///|
fn bcf_str_len(s : String) -> Int {
  bcf_size_head_len(s.length()) + s.length()
}

///|
fn bcf_ivec_len(vs : Array[Int]) -> Int {
  if vs.length() == 0 {
    1
  } else {
    let t = bcf_int_type(vs)
    let w = if t == 1 { 1 } else if t == 2 { 2 } else { 4 }
    bcf_size_head_len(vs.length()) + vs.length() * w
  }
}

///|
fn bcf_values_to_ints(values : Array[BcfValue]) -> Array[Int] {
  let out : Array[Int] = []
  for v in values {
    match v {
      BcfInt(x) => out.push(x)
      _ => out.push(0)
    }
  }
  out
}

///|
fn bcf_val_len(values : Array[BcfValue]) -> Int {
  match values[0] {
    BcfFlag => 1
    BcfStr(s) => bcf_str_len(s)
    BcfInt(_) => bcf_ivec_len(bcf_values_to_ints(values))
    BcfFloat(_) => bcf_size_head_len(values.length()) + values.length() * 4
  }
}

///|
/// Write a typed integer (type byte 1/2/3 + value), used for ids and large
/// counts.
fn bcf_w_tint(buf : FixedArray[UInt16], pos : Int, v : Int) -> Int {
  if v < 128 {
    buf[pos] = 0x01
    bam_write_u8(buf, pos + 1, v)
  } else if v < 32768 {
    buf[pos] = 0x02
    bam_write_u16_le(buf, pos + 1, v)
  } else {
    buf[pos] = 0x03
    bam_write_i32_le(buf, pos + 1, v)
  }
}

///|
/// Write the typed-value size header (`count<<4|type`, with a typed-int count
/// when `count >= 15`).
fn bcf_w_size(buf : FixedArray[UInt16], pos : Int, count : Int, type_ : Int) -> Int {
  if count < 15 {
    buf[pos] = ((count << 4) | type_).to_uint16()
    pos + 1
  } else {
    buf[pos] = (0xF0 | type_).to_uint16()
    bcf_w_tint(buf, pos + 1, count)
  }
}

///|
fn bcf_w_string(buf : FixedArray[UInt16], pos : Int, s : String) -> Int {
  let p = bcf_w_size(buf, pos, s.length(), 7)
  let mut i = 0
  while i < s.length() {
    buf[p + i] = s.unsafe_get(i)
    i = i + 1
  }
  p + s.length()
}

///|
fn bcf_w_ivec(buf : FixedArray[UInt16], pos : Int, vs : Array[Int]) -> Int {
  if vs.length() == 0 {
    buf[pos] = 0x00
    return pos + 1
  }
  let t = bcf_int_type(vs)
  let mut p = bcf_w_size(buf, pos, vs.length(), t)
  for v in vs {
    if t == 1 {
      p = bam_write_u8(buf, p, v)
    } else if t == 2 {
      p = bam_write_u16_le(buf, p, v)
    } else {
      p = bam_write_i32_le(buf, p, v)
    }
  }
  p
}

///|
fn bcf_w_value(buf : FixedArray[UInt16], pos : Int, values : Array[BcfValue]) -> Int {
  match values[0] {
    BcfFlag => {
      buf[pos] = 0x00
      pos + 1
    }
    BcfStr(s) => bcf_w_string(buf, pos, s)
    BcfInt(_) => bcf_w_ivec(buf, pos, bcf_values_to_ints(values))
    BcfFloat(_) => {
      let mut p = bcf_w_size(buf, pos, values.length(), 5)
      for v in values {
        match v {
          BcfFloat(f) => p = bam_write_f32(buf, p, f)
          _ => p = bam_write_f32(buf, p, 0.0)
        }
      }
      p
    }
  }
}

///|
fn bcf_encode_record(rec : BcfRecord) -> String {
  // shared / individual byte lengths
  let mut shared_len = (bcf_str_len(rec.id) + bcf_str_len(rec.ref_seq) + bcf_ivec_len(rec.filter_ids))
  for a in rec.alt {
    shared_len = shared_len + bcf_str_len(a)
  }
  for info in rec.info {
    shared_len = shared_len + bcf_tint_len(info.key_id) + bcf_val_len(info.values)
  }
  let mut indiv_len = 0
  for fmt in rec.formats {
    indiv_len = indiv_len + bcf_tint_len(fmt.key_id) + bcf_val_len(fmt.values)
  }

  let l_shared = 24 + shared_len
  let total = 32 + shared_len + indiv_len
  let buf : FixedArray[UInt16] = FixedArray::make(total, 0)

  let mut p = bam_write_i32_le(buf, 0, l_shared)
  p = bam_write_i32_le(buf, p, indiv_len)
  p = bam_write_i32_le(buf, p, rec.rid)
  p = bam_write_i32_le(buf, p, rec.pos)
  p = bam_write_i32_le(buf, p, rec.rlen)
  p = match rec.qual {
    None => bam_write_f32(buf, p, 0.0 / 0.0)
    Some(q) => bam_write_f32(buf, p, q)
  }
  p = bam_write_u16_le(buf, p, rec.info.length())
  p = bam_write_u16_le(buf, p, rec.alt.length() + 1)
  // n_sample (24-bit) + n_fmt (8-bit)
  let ns = rec.n_sample & 0xFFFFFF
  buf[p] = (ns & 0xFF).to_uint16()
  buf[p + 1] = ((ns >> 8) & 0xFF).to_uint16()
  buf[p + 2] = ((ns >> 16) & 0xFF).to_uint16()
  buf[p + 3] = (rec.n_fmt & 0xFF).to_uint16()
  p = p + 4

  // shared block
  p = bcf_w_string(buf, p, rec.id)
  p = bcf_w_string(buf, p, rec.ref_seq)
  for a in rec.alt {
    p = bcf_w_string(buf, p, a)
  }
  p = bcf_w_ivec(buf, p, rec.filter_ids)
  for info in rec.info {
    p = bcf_w_tint(buf, p, info.key_id)
    p = bcf_w_value(buf, p, info.values)
  }
  // individual block
  for fmt in rec.formats {
    p = bcf_w_tint(buf, p, fmt.key_id)
    p = bcf_w_value(buf, p, fmt.values)
  }

  unsafe_from_uint16_fixedarray(buf)
}

///|
/// Serialize a parsed `BcfFile` back to raw BCF 2.2 bytes (magic + header text
/// + records).
pub fn BcfFile::to_bcf_bytes(self : BcfFile) -> String {
  let text = self.header.text
  let l_text = text.length()
  let mut records_len = 0
  let encoded : Array[String] = []
  for rec in self.records {
    let e = bcf_encode_record(rec)
    records_len = records_len + e.length()
    encoded.push(e)
  }

  let total = 5 + 4 + l_text + records_len
  let buf : FixedArray[UInt16] = FixedArray::make(total, 0)
  // magic "BCF\x02\x02"
  buf[0] = 66
  buf[1] = 67
  buf[2] = 70
  buf[3] = 2
  buf[4] = 2
  let mut p = 5
  p = bam_write_i32_le(buf, p, l_text)
  // header text
  let mut i = 0
  while i < l_text {
    buf[p + i] = text.unsafe_get(i)
    i = i + 1
  }
  p = p + l_text
  // records
  for e in encoded {
    let mut j = 0
    while j < e.length() {
      buf[p + j] = e.unsafe_get(j)
      j = j + 1
    }
    p = p + e.length()
  }
  unsafe_from_uint16_fixedarray(buf)
}