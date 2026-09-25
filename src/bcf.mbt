///|
/// BCF (binary VCF) parser — the binary companion to text VCF, as written by
/// htslib / bcftools (`BCFv2.2`).
///
/// Port of the htslib BCF2 format. A BCF file is a single BGZF stream whose
/// decompressed contents are: the magic `BCF\2\2`, a little-endian text
/// length, the VCF header text (the header dictionary is implicit in the
/// parsed `##` lines), then one record per variant. Each record is a 32-byte
/// fixed header followed by a shared block (ID, REF/ALT, FILTER, INFO) and an
/// individual block (FORMAT + sample values), all using htslib's typed
/// encoding (a size byte `count<<4 | type`, plus int8/int16/int32, float32,
/// char or null values).
///
/// The compressed BGZF layer is handled by `src/bgzf.mbt`; `parse_bcf`
/// expects the already-decompressed BCF byte stream.

///|
pub suberror BcfError {
  BcfError(String)
}

///|
fn[T] bcf_fail(msg : String) -> T raise BcfError {
  raise BcfError(msg)
}

///|
pub enum BcfValue {
  BcfInt(Int)
  BcfFloat(Double)
  BcfStr(String)
  BcfFlag
} derive(Show)

///|
/// A parsed BCF header: the raw text plus the name dictionaries derived from
/// it. FILTER/INFO/FORMAT ids share one integer space (`id_names`).
pub struct BcfHeader {
  text : String
  contig_names : Array[String]
  contig_lengths : Array[Int?]
  id_names : Array[String]
  filter_ids : Array[Int]
  info_ids : Array[Int]
  format_ids : Array[Int]
  sample_names : Array[String]
}

///|
pub struct BcfInfo {
  key_id : Int
  values_type : Int
  values : Array[BcfValue]
}

///|
pub struct BcfFormat {
  key_id : Int
  values_type : Int
  values_per_sample : Int
  values : Array[BcfValue]
}

///|
pub struct BcfRecord {
  rid : Int
  chrom : String
  pos : Int
  rlen : Int
  qual : Double?
  id : String
  ref_seq : String
  alt : Array[String]
  filter_ids : Array[Int]
  info : Array[BcfInfo]
  formats : Array[BcfFormat]
  n_sample : Int
  n_fmt : Int
}

///|
pub struct BcfFile {
  header : BcfHeader
  records : Array[BcfRecord]
}

// ---------------------------------------------------------------------------
// Little-endian readers
// ---------------------------------------------------------------------------

///|
fn bcf_u8(data : String, pos : Int) -> Int {
  data.unsafe_get(pos).to_int() & 0xFF
}

///|
fn bcf_check(data : String, pos : Int, len : Int) -> Unit raise BcfError {
  if pos < 0 || pos + len > data.length() {
    bcf_fail("Premature end of BCF data at offset \{pos} (length \{data.length()})")
  }
}

///|
fn bcf_read_u16(data : String, pos : Int) -> Int raise BcfError {
  bcf_check(data, pos, 2)
  bcf_u8(data, pos) | (bcf_u8(data, pos + 1) << 8)
}

///|
fn bcf_read_u32(data : String, pos : Int) -> Int raise BcfError {
  bcf_check(data, pos, 4)
  bcf_u8(data, pos) |
    (bcf_u8(data, pos + 1) << 8) |
    (bcf_u8(data, pos + 2) << 16) |
    (bcf_u8(data, pos + 3) << 24)
}

///|
fn bcf_read_i8(data : String, pos : Int) -> Int {
  let b = bcf_u8(data, pos)
  if b >= 128 {
    b - 256
  } else {
    b
  }
}

///|
fn bcf_read_i16(data : String, pos : Int) -> Int raise BcfError {
  let v = bcf_read_u16(data, pos)
  if v >= 32768 {
    v - 65536
  } else {
    v
  }
}

///|
/// Decode an IEEE 754 float32 (4 little-endian bytes) to a Double. Non-finite
/// values (exponent all ones, e.g. the missing sentinel) become NaN.
fn bcf_f32(b0 : Int, b1 : Int, b2 : Int, b3 : Int) -> Double {
  let sign = if b3 >= 128 { -1.0 } else { 1.0 }
  let exponent = (b3 & 0x7F) * 2 + (b2 >> 7)
  let mantissa = (b2 & 0x7F) * 65536 + b1 * 256 + b0
  if exponent == 255 {
    return 0.0 / 0.0
  }
  if exponent == 0 {
    sign * mantissa.to_double() * @math.pow(2.0, -149.0)
  } else {
    sign * (1.0 + mantissa.to_double() / 8388608.0) * @math.pow(
      2.0,
      (exponent - 127).to_double(),
    )
  }
}

///|
fn bcf_read_f32(data : String, pos : Int) -> Double raise BcfError {
  bcf_check(data, pos, 4)
  bcf_f32(bcf_u8(data, pos), bcf_u8(data, pos + 1), bcf_u8(data, pos + 2), bcf_u8(data, pos + 3))
}

// ---------------------------------------------------------------------------
// Typed-value decoding
// ---------------------------------------------------------------------------

///|
/// Decode a size byte (`count<<4 | type`). When the high nibble is 15, the
/// count is a typed integer that follows. Returns (count, type, next).
fn bcf_read_size(data : String, pos : Int) -> (Int, Int, Int) raise BcfError {
  bcf_check(data, pos, 1)
  let b = bcf_u8(data, pos)
  let type_ = b & 0xF
  let high = b >> 4
  if high != 15 {
    (high, type_, pos + 1)
  } else {
    let (count, next) = bcf_read_typed_int(data, pos + 1)
    (count, type_, next)
  }
}

///|
/// Decode a single typed integer (a size byte `1<<4|type` + the value).
fn bcf_read_typed_int(data : String, pos : Int) -> (Int, Int) raise BcfError {
  bcf_check(data, pos, 2)
  let type_ = bcf_u8(data, pos) & 0xF
  if type_ == 1 {
    let v = bcf_read_i8(data, pos + 1)
    (v, pos + 2)
  } else if type_ == 2 {
    let v = bcf_read_i16(data, pos + 1)
    (v, pos + 3)
  } else if type_ == 3 {
    let v = bcf_read_u32(data, pos + 1)
    (v, pos + 5)
  } else {
    bcf_fail("Invalid typed integer type \{type_}")
  }
}

///|
/// Decode `count` typed values of the given type. Returns (values, next).
fn bcf_decode_values(
  data : String,
  pos : Int,
  count : Int,
  type_ : Int,
) -> (Array[BcfValue], Int) raise BcfError {
  let values : Array[BcfValue] = []
  let mut p = pos
  if type_ == 0 {
    // BCF_BT_NULL (flag): no data
    values.push(BcfFlag)
  } else if type_ == 7 {
    // BCF_BT_CHAR: a single string of `count` bytes
    bcf_check(data, p, count)
    values.push(BcfStr(data[p:p + count].to_owned()))
    p = p + count
  } else if type_ == 1 || type_ == 2 || type_ == 3 {
    let width = if type_ == 1 { 1 } else if type_ == 2 { 2 } else { 4 }
    bcf_check(data, p, count * width)
    for _i = 0; _i < count; _i = _i + 1 {
      let v = if type_ == 1 {
        bcf_read_i8(data, p)
      } else if type_ == 2 {
        bcf_read_i16(data, p)
      } else {
        bcf_read_u32(data, p)
      }
      values.push(BcfInt(v))
      p = p + width
    }
  } else if type_ == 5 {
    // BCF_BT_FLOAT
    bcf_check(data, p, count * 4)
    for _i = 0; _i < count; _i = _i + 1 {
      values.push(BcfFloat(bcf_read_f32(data, p)))
      p = p + 4
    }
  } else {
    bcf_fail("Unknown BCF value type \{type_}")
  }
  (values, p)
}

///|
/// Decode a typed string (a size byte `count<<4 | CHAR` + `count` bytes).
fn bcf_read_string(data : String, pos : Int) -> (String, Int) raise BcfError {
  let (count, type_, next) = bcf_read_size(data, pos)
  if type_ != 7 {
    bcf_fail("Expected string type, got \{type_}")
  }
  if count == 0 {
    ("", next)
  } else {
    bcf_check(data, next, count)
    (data[next:next + count].to_owned(), next + count)
  }
}

///|
/// Decode a typed vector of integers (FILTER). A null vector or a single
/// integer are both handled by the generic size byte.
fn bcf_read_int_vector(data : String, pos : Int) -> (Array[Int], Int) raise BcfError {
  let (count, type_, next) = bcf_read_size(data, pos)
  if type_ == 0 {
    return ([], next)
  }
  if type_ != 1 && type_ != 2 && type_ != 3 {
    bcf_fail("Expected integer vector, got type \{type_}")
  }
  let width = if type_ == 1 { 1 } else if type_ == 2 { 2 } else { 4 }
  let out : Array[Int] = []
  let mut p = next
  for _i = 0; _i < count; _i = _i + 1 {
    let v = if type_ == 1 {
      bcf_read_i8(data, p)
    } else if type_ == 2 {
      bcf_read_i16(data, p)
    } else {
      bcf_read_u32(data, p)
    }
    out.push(v)
    p = p + width
  }
  (out, p)
}

// ---------------------------------------------------------------------------
// Header text parsing
// ---------------------------------------------------------------------------

///|
/// Split the header text into lines (NUL / CRLF tolerant).
fn bcf_lines(text : String) -> Array[String] {
  let raw = split_by_char(text, '\n'.to_int())
  let out : Array[String] = []
  for line in raw {
    let mut s = line
    if s.length() > 0 && s.unsafe_get(s.length() - 1).to_int() == 13 {
      s = s[0:s.length() - 1].to_owned()
    }
    // strip a trailing NUL (htslib appends one after the final newline)
    if s.length() > 0 && s.unsafe_get(s.length() - 1).to_int() == 0 {
      s = s[0:s.length() - 1].to_owned()
    }
    if s.length() > 0 {
      out.push(s)
    }
  }
  out
}

///|
/// Extract the value of a `TAG=` field inside a `##...=<...>` line.
fn bcf_tag_value(line : String, tag : String) -> String {
  let needle = tag + "="
  let mut i = 0
  while i + needle.length() <= line.length() {
    if starts_with(line, i, needle) {
      let boundary = i == 0 ||
        line.unsafe_get(i - 1).to_int() == 44 ||
        line.unsafe_get(i - 1).to_int() == 60
      if boundary {
        let mut j = i + needle.length()
        let buf = StringBuilder()
        while j < line.length() {
          let c = line.unsafe_get(j).to_int()
          if c == 44 || c == 62 {
            break
          }
          buf.write_char(line.unsafe_get(j).unsafe_to_char())
          j = j + 1
        }
        return buf.to_string()
      }
    }
    i = i + 1
  }
  ""
}

///|
fn bcf_parse_header_text(text : String) -> BcfHeader {
  let contig_names : Array[String] = []
  let contig_lengths : Array[Int?] = []
  let id_names : Array[String] = []
  let filter_ids : Array[Int] = []
  let info_ids : Array[Int] = []
  let format_ids : Array[Int] = []
  let sample_names : Array[String] = []
  let mut id_counter = 0
  for line in bcf_lines(text) {
    if starts_with(line, 0, "##contig=") {
      contig_names.push(bcf_tag_value(line, "ID"))
      let len = bcf_tag_value(line, "length")
      contig_lengths.push(if len.length() == 0 { None } else { Some(parse_int(len)) })
    } else if starts_with(line, 0, "##FILTER=") {
      id_names.push(bcf_tag_value(line, "ID"))
      filter_ids.push(id_counter)
      id_counter = id_counter + 1
    } else if starts_with(line, 0, "##INFO=") {
      id_names.push(bcf_tag_value(line, "ID"))
      info_ids.push(id_counter)
      id_counter = id_counter + 1
    } else if starts_with(line, 0, "##FORMAT=") {
      id_names.push(bcf_tag_value(line, "ID"))
      format_ids.push(id_counter)
      id_counter = id_counter + 1
    } else if starts_with(line, 0, "#CHROM") {
      let cols = split_by_tab(line)
      let mut j = 9
      while j < cols.length() {
        sample_names.push(cols[j])
        j = j + 1
      }
    } else {
      // ##fileformat and other structured lines carry no id space of interest
    }
  }
  BcfHeader::{
    text,
    contig_names,
    contig_lengths,
    id_names,
    filter_ids,
    info_ids,
    format_ids,
    sample_names,
  }
}

// ---------------------------------------------------------------------------
// Parsing
// ---------------------------------------------------------------------------

///|
fn bcf_id_name(header : BcfHeader, id : Int) -> String {
  if id >= 0 && id < header.id_names.length() {
    header.id_names[id]
  } else {
    id.to_string()
  }
}

///|
fn bcf_contig_name(header : BcfHeader, rid : Int) -> String {
  if rid >= 0 && rid < header.contig_names.length() {
    header.contig_names[rid]
  } else {
    rid.to_string()
  }
}

///|
/// Parse a single BCF record starting at `off`, using `header` for name
/// resolution. Returns (record, next offset).
fn bcf_parse_record(
  data : String,
  off : Int,
  header : BcfHeader,
) -> (BcfRecord, Int) raise BcfError {
  bcf_check(data, off, 32)
  let l_shared = bcf_read_u32(data, off)
  let l_indiv = bcf_read_u32(data, off + 4)
  let rid = bcf_read_u32(data, off + 8)
  let pos = bcf_read_u32(data, off + 12)
  let rlen = bcf_read_u32(data, off + 16)
  let qual_f = bcf_read_f32(data, off + 20)
  let n_info = bcf_read_u16(data, off + 24)
  let n_allele = bcf_read_u16(data, off + 26)
  let n_sample = bcf_read_u32(data, off + 28) & 0xffffff
  let n_fmt = bcf_u8(data, off + 31)

  let qual : Double? = if qual_f.is_nan() { None } else { Some(qual_f) }

  // shared block: 24 bytes of the fixed header are counted inside l_shared
  let shared_len = l_shared - 24
  if shared_len < 0 {
    bcf_fail("Invalid BCF record: l_shared < 24")
  }
  let mut p = off + 32
  let shared_end = p + shared_len
  let indiv_end = shared_end + l_indiv

  // ID
  let (id, p2) = bcf_read_string(data, p)
  p = p2
  if p > shared_end {
    bcf_fail("BCF ID field overruns the shared block")
  }

  // REF + ALT
  let (ref_seq, p3) = bcf_read_string(data, p)
  p = p3
  let alt : Array[String] = []
  for i = 1; i < n_allele; i = i + 1 {
    let (a, p4) = bcf_read_string(data, p)
    p = p4
    alt.push(a)
  }

  // FILTER
  let (filter_ids, p5) = bcf_read_int_vector(data, p)
  p = p5

  // INFO
  let info : Array[BcfInfo] = []
  for i = 0; i < n_info; i = i + 1 {
    let (key_id, p6) = bcf_read_typed_int(data, p)
    p = p6
    let (count, type_, p7) = bcf_read_size(data, p)
    p = p7
    let (values, p8) = bcf_decode_values(data, p, count, type_)
    p = p8
    info.push(BcfInfo::{ key_id, values_type: type_, values })
  }

  // individual block: FORMAT fields
  let formats : Array[BcfFormat] = []
  let mut q = shared_end
  for i = 0; i < n_fmt; i = i + 1 {
    let (key_id, q2) = bcf_read_typed_int(data, q)
    q = q2
    let (n_per, type_, q3) = bcf_read_size(data, q)
    q = q3
    let total = n_per * n_sample
    let (values, q4) = bcf_decode_values(data, q, total, type_)
    q = q4
    formats.push(BcfFormat::{
      key_id,
      values_type: type_,
      values_per_sample: n_per,
      values,
    })
  }

  let record = BcfRecord::{
    rid,
    chrom: bcf_contig_name(header, rid),
    pos,
    rlen,
    qual,
    id,
    ref_seq,
    alt,
    filter_ids,
    info,
    formats,
    n_sample,
    n_fmt,
  }
  (record, indiv_end)
}

///|
/// Parse a BCF header from the raw (decompressed) byte stream.
pub fn parse_bcf_header(data : String) -> BcfHeader raise BcfError {
  bcf_check(data, 0, 5)
  if data.unsafe_get(0).to_int() != 66 ||
    data.unsafe_get(1).to_int() != 67 ||
    data.unsafe_get(2).to_int() != 70 ||
    data.unsafe_get(3).to_int() != 2 ||
    data.unsafe_get(4).to_int() != 2 {
    bcf_fail("Invalid BCF2 magic string")
  }
  let l_text = bcf_read_u32(data, 5)
  bcf_check(data, 9, l_text)
  let text = data[9:9 + l_text].to_owned()
  bcf_parse_header_text(text)
}

///|
/// Parse a whole BCF byte stream (header + all records).
pub fn parse_bcf(data : String) -> BcfFile raise BcfError {
  let header = parse_bcf_header(data)
  let mut off = 9 + bcf_read_u32(data, 5)
  let records : Array[BcfRecord] = []
  while off < data.length() {
    let (record, next) = bcf_parse_record(data, off, header)
    records.push(record)
    off = next
  }
  BcfFile::{ header, records }
}

///|
/// Decode a packed GT integer to its allele index (None for missing).
pub fn bcf_gt_allele(encoded : Int) -> Int? {
  if encoded <= 0 {
    None
  } else {
    Some((encoded >> 1) - 1)
  }
}

///|
/// Whether a packed GT integer marks a phased allele.
pub fn bcf_gt_phased(encoded : Int) -> Bool {
  (encoded & 1) == 1
}

///|
/// The name for a FILTER/INFO/FORMAT id, falling back to the numeric id.
pub fn BcfHeader::id_name(self : BcfHeader, id : Int) -> String {
  bcf_id_name(self, id)
}

///|
/// The contig name for a reference id, falling back to the numeric id.
pub fn BcfHeader::contig_name(self : BcfHeader, rid : Int) -> String {
  bcf_contig_name(self, rid)
}

///|
/// The number of records in a parsed BCF file.
pub fn BcfFile::len(self : BcfFile) -> Int {
  self.records.length()
}