// rtracklayer-style UCSC BigWig binary format.
//
// BigWig is the density-track sibling of BigBed: both share the bbi container
// (64-byte header, chromosome B+ tree, interval R-tree, zlib blocks, total
// summary and zoom levels).  The low-level container machinery lives in
// bigbed.mbt; this module owns the bigWig-specific on-disk parts:
//
//   * magic 0x888FFC26 (version 3/4), field count 0 and no AutoSQL
//   * a u64 section count before the data blocks
//   * 24-byte section headers followed by bedGraph (type 1), variableStep
//     (type 2) or fixedStep (type 3) items
//   * 40-byte total summary (u64 valid count + four f64) and 24-byte zoom
//     headers pointing at R-trees of 32-byte f32 summary elements
//
// The reference implementation is libBigWig (as used by pyBigWig); writer
// output is cross-validated against pyBigWig and the reader accepts both
// libBigWig and Kent tool output conventions.

///|
pub suberror BigWigError {
  BigWigError(String)
}

///|
/// One dense interval: chromosome name, 0-based half-open range, value.

///|
pub struct BigWigInterval {
  chrom : String
  start : Int
  end : Int
  value : Double
} derive(Eq, Debug)

///|
pub fn BigWigInterval::create(
  chrom : String,
  start : Int,
  end : Int,
  value : Double,
) -> BigWigInterval {
  BigWigInterval::{ chrom, start, end, value, }
}

///|
/// Aggregate statistics over a range (total summary or a zoom element).
pub struct BigWigSummary {
  valid_count : Int
  min_val : Double
  max_val : Double
  sum_data : Double
  sum_squares : Double
} derive(Eq, Debug)

///|
/// One zoom (reduction) level: reduction width and file offsets.
pub struct BigWigZoomHeader {
  reduction_level : Int
  data_offset : Int
  index_offset : Int
} derive(Eq, Debug)

///|
pub struct BigWigHeader {
  little_endian : Bool
  version : Int
  zoom_levels : Int
  chromosome_tree_offset : Int
  full_data_offset : Int
  full_index_offset : Int
  total_summary_offset : Int
  uncompress_buffer_size : Int
} derive(Eq, Debug)

///|
pub struct BigWigWriteConfig {
  block_size : Int
  compress : Bool
} derive(Eq, Debug)

///|
pub struct BigWigFile {
  header : BigWigHeader
  targets : Array[BigBedTarget]
  source : Array[Int]
  blocks : Array[BigBedBlockIndex]
  zoom_headers : Array[BigWigZoomHeader]
  total_summary : BigWigSummary
  section_count : Int
}

///|
let bigwig_magic : Int = 0x888FFC26

///|
/// libBigWig on-disk section type codes (Kent/Kent-compatible 24-byte header).
let bigwig_bed_graph : Int = 1

///|
let bigwig_variable_step : Int = 2

///|
let bigwig_fixed_step : Int = 3

///|
/// Sections hold at most 32 KiB uncompressed; bedGraph items are 12 bytes.
let bigwig_max_bedgraph_items : Int = 2728

///|
fn[T] bigwig_wrap(operation : () -> T raise BigBedError) -> T raise BigWigError {
  operation() catch {
    BigBedError(message) => raise BigWigError(message)
  }
}

///|
pub fn BigWigWriteConfig::default() -> BigWigWriteConfig {
  BigWigWriteConfig::{ block_size: 64, compress: true, }
}

///|
pub fn BigWigWriteConfig::create(
  compress? : Bool = true,
  block_size? : Int = 64,
) -> BigWigWriteConfig raise BigWigError {
  bigwig_wrap(fn() -> BigWigWriteConfig raise BigBedError {
    if block_size < 1 || block_size > 65535 {
      bigbed_fail("BigWig block size must be between 1 and 65535")
    }
    BigWigWriteConfig::{ block_size, compress, }
  })
}

// ---------------------------------------------------------------------------
// IEEE 754 float codecs (little-endian on disk)
//
// The local Int type is 32 bit, so f64 patterns are carried as a
// (low32, high32) pair.
// ---------------------------------------------------------------------------

///|
/// Interpret a signed i32 read as its unsigned f64 magnitude.
fn bigwig_u32_to_double(value : Int) -> Double {
  if value < 0 {
    value.to_double() + 4294967296.0
  } else {
    value.to_double()
  }
}

///|
/// Split a non-negative integer-valued Double in [0, 2^53) into its low and
/// high 32-bit halves, rounding to nearest with ties to even.
fn bigwig_split_52(value : Double) -> (Int, Int) {
  let scaled = value.round()
  let high = (scaled / 4294967296.0).floor().to_int()
  let low_full = scaled - high.to_double() * 4294967296.0
  let low_u16 = (low_full / 65536.0).floor().to_int()
  let low_l = (low_full - low_u16.to_double() * 65536.0).round().to_int()
  let low = (low_u16 << 16) | low_l
  (low, high)
}

///|
/// Encode an f64 as (little-endian low 32 bits, high 32 bits).
fn bigwig_f64_words(value : Double) -> (Int, Int) {
  if value.is_nan() {
    return (0, 0x7FF80000)
  }
  let sign_bit : Int = if value < 0.0 { 1 } else { 0 }
  let abs = if value < 0.0 { -value } else { value }
  if abs == 0.0 {
    return (0, sign_bit << 31)
  }
  if abs > 1.7976931348623157e308 {
    return (0, (sign_bit << 31) | 0x7FF00000)
  }
  let mut normalized = abs
  let mut exponent = 0
  while normalized >= 2.0 {
    normalized = normalized / 2.0
    exponent = exponent + 1
  }
  while normalized < 1.0 {
    normalized = normalized * 2.0
    exponent = exponent - 1
  }
  let biased = exponent + 1023
  if biased >= 0x7FF {
    return (0, (sign_bit << 31) | 0x7FF00000)
  }
  if biased <= 0 {
    // Subnormal f64: significand is abs * 2^1074 in [0, 2^52).
    let mut scaled = abs
    for _ in 0..<1074 {
      scaled = scaled * 2.0
    }
    let (low, high) = bigwig_split_52(scaled)
    return (low, (sign_bit << 31) | (high & 0x000FFFFF))
  }
  let (low, raw_high) = bigwig_split_52(normalized * 4503599627370496.0)
  let mut high = raw_high
  let mut exp_field = biased
  if high == 0x00200000 {
    // Rounding carried into the implicit bit.
    high = 0
    exp_field = exp_field + 1
  }
  (low, (sign_bit << 31) | ((exp_field & 0x7FF) << 20) | (high & 0x000FFFFF))
}

///|
fn bigwig_f64_from_words(low : Int, high : Int) -> Double {
  let sign = if high < 0 { -1.0 } else { 1.0 }
  let exponent = (high >> 20) & 0x7FF
  let mantissa_high = (high & 0x000FFFFF).to_double()
  let mantissa_low = bigwig_u32_to_double(low)
  let mantissa = mantissa_high * 4294967296.0 + mantissa_low
  if exponent == 0 {
    if mantissa == 0.0 {
      return sign * 0.0
    }
    return sign * mantissa * 4.450147717014403e-308
  }
  if exponent == 0x7FF {
    if mantissa == 0.0 {
      return sign * (1.0 / 0.0)
    } else {
      return 0.0 / 0.0
    }
  }
  let significand = mantissa + 4503599627370496.0
  let mut scale = 1.0
  let mut power = exponent - 1075
  while power > 0 {
    scale = scale * 2.0
    power = power - 1
  }
  while power < 0 {
    scale = scale / 2.0
    power = power + 1
  }
  sign * significand * scale
}

///|
fn bigwig_push_f64(output : Array[Int], value : Double) -> Unit {
  let (low, high) = bigwig_f64_words(value)
  bigbed_push_u32(output, low)
  bigbed_push_u32(output, high)
}

///|
fn bigwig_read_f64(
  data : Array[Int],
  position : Int,
  endian : BigBedEndian,
) -> Double raise BigBedError {
  if endian == BigBedLittleEndian {
    let low = bigbed_read_u32(data, position, endian)
    let high = bigbed_read_u32(data, position + 4, endian)
    bigwig_f64_from_words(low, high)
  } else {
    let high = bigbed_read_u32(data, position, endian)
    let low = bigbed_read_u32(data, position + 4, endian)
    bigwig_f64_from_words(low, high)
  }
}

///|
/// Encode an f64 value as IEEE 754 f32 bits (round to nearest, ties to even).
fn bigwig_f32_bits(value : Double) -> Int {
  if value.is_nan() {
    return 0x7FC00000
  }
  let sign_bit : Int = if value < 0.0 { 1 } else { 0 }
  let mut abs = if value < 0.0 { -value } else { value }
  if abs == 0.0 {
    return sign_bit << 31
  }
  if abs > 3.4028235e38 {
    return (sign_bit << 31) | 0x7F800000
  }
  let mut exponent = 0
  while abs >= 2.0 {
    abs = abs / 2.0
    exponent = exponent + 1
  }
  while abs < 1.0 {
    abs = abs * 2.0
    exponent = exponent - 1
  }
  let biased = exponent + 127
  if biased >= 0xFF {
    return (sign_bit << 31) | 0x7F800000
  }
  if biased <= 0 {
    let mantissa = (abs * 7.0368744177664e44).round().to_int()
    return (sign_bit << 31) | (mantissa & 0x007FFFFF)
  }
  let mut significand = (abs * 8388608.0).round().to_int()
  let mut exp_field = biased
  if significand == 16777216 {
    significand = 8388608
    exp_field = exp_field + 1
  }
  (sign_bit << 31) | (exp_field << 23) | (significand & 0x007FFFFF)
}

///|
fn bigwig_f32_from_bits(bits : Int) -> Double {
  let sign = if bits < 0 { -1.0 } else { 1.0 }
  let exponent = (bits >> 23) & 0xFF
  let mantissa = (bits & 0x007FFFFF).to_double()
  if exponent == 0 {
    if mantissa == 0.0 {
      return sign * 0.0
    }
    return sign * mantissa * 1.401298464324817e-45
  }
  if exponent == 0xFF {
    if mantissa == 0.0 {
      return sign * (1.0 / 0.0)
    } else {
      return 0.0 / 0.0
    }
  }
  let mut scale = 1.0
  let mut power = exponent - 150
  while power > 0 {
    scale = scale * 2.0
    power = power - 1
  }
  while power < 0 {
    scale = scale / 2.0
    power = power + 1
  }
  sign * (mantissa + 8388608.0) * scale
}

///|
fn bigwig_push_f32(output : Array[Int], value : Double) -> Unit {
  bigbed_push_u32(output, bigwig_f32_bits(value))
}

///|
fn bigwig_read_f32(
  data : Array[Int],
  position : Int,
  endian : BigBedEndian,
) -> Double raise BigBedError {
  bigwig_f32_from_bits(bigbed_read_u32(data, position, endian))
}

// ---------------------------------------------------------------------------
// Header / summary parsing
// ---------------------------------------------------------------------------

///|
fn bigwig_read_header(
  data : Array[Int],
) -> (BigWigHeader, BigBedEndian) raise BigBedError {
  bigbed_check_range(data, 0, 64, "header")
  // bigWig shares the bbi container but uses magic 0x888FFC26, unlike BigBed.
  let little = data[0] | (data[1] << 8) | (data[2] << 16) | (data[3] << 24)
  let endian = if little == bigwig_magic {
    BigBedLittleEndian
  } else {
    let big = (data[0] << 24) | (data[1] << 16) | (data[2] << 8) | data[3]
    if big == bigwig_magic {
      BigBedBigEndian
    } else {
      bigbed_fail("Input does not contain the BigWig magic number")
      BigBedLittleEndian
    }
  }
  let version = bigbed_read_u16(data, 4, endian)
  if version != 3 && version != 4 {
    bigbed_fail(
      "Unsupported BigWig version " +
      version.to_string() +
      "; expected version 3 or 4",
    )
  }
  let field_count = bigbed_read_u16(data, 32, endian)
  let defined = bigbed_read_u16(data, 34, endian)
  if field_count != 0 || defined != 0 {
    bigbed_fail("BigWig files must not declare BED fields")
  }
  let zoom_levels = bigbed_read_u16(data, 6, endian)
  let chrom_tree_offset = bigbed_read_u64_local(data, 8, endian)
  let full_data_offset = bigbed_read_u64_local(data, 16, endian)
  let full_index_offset = bigbed_read_u64_local(data, 24, endian)
  let total_summary_offset = bigbed_read_u64_local(data, 44, endian)
  let header = BigWigHeader::{
    little_endian: endian == BigBedLittleEndian,
    version,
    zoom_levels,
    chromosome_tree_offset: chrom_tree_offset,
    full_data_offset,
    full_index_offset,
    total_summary_offset,
    uncompress_buffer_size: bigbed_read_u32(data, 52, endian),
  }
  let zoom_table_end = 64 + zoom_levels * 24
  if chrom_tree_offset < zoom_table_end ||
    full_data_offset <= chrom_tree_offset ||
    full_index_offset < full_data_offset ||
    full_index_offset >= data.length() {
    bigbed_fail("BigWig header contains invalid section offsets")
  }
  if total_summary_offset != 0 &&
    (
      total_summary_offset < zoom_table_end ||
      total_summary_offset + 40 > chrom_tree_offset
    ) {
    bigbed_fail("BigWig total summary offset is outside the metadata section")
  }
  (header, endian)
}

///|
fn bigwig_read_zoom_headers(
  data : Array[Int],
  count : Int,
  endian : BigBedEndian,
) -> Array[BigWigZoomHeader] raise BigBedError {
  let headers : Array[BigWigZoomHeader] = []
  for i in 0..<count {
    let position = 64 + i * 24
    let reduction = bigbed_read_u32(data, position, endian)
    let data_offset = bigbed_read_u64_local(data, position + 8, endian)
    let index_offset = bigbed_read_u64_local(data, position + 16, endian)
    if reduction == 0 || data_offset == 0 || index_offset == 0 {
      bigbed_fail("BigWig zoom header contains a null offset")
    }
    headers.push(BigWigZoomHeader::{
      reduction_level: reduction,
      data_offset,
      index_offset,
    })
  }
  headers
}

///|
fn bigwig_read_summary(
  data : Array[Int],
  offset : Int,
  endian : BigBedEndian,
) -> BigWigSummary raise BigBedError {
  if offset == 0 {
    return BigWigSummary::{
      valid_count: 0,
      min_val: 0.0,
      max_val: 0.0,
      sum_data: 0.0,
      sum_squares: 0.0,
    }
  }
  bigbed_check_range(data, offset, 40, "total summary")
  BigWigSummary::{
    valid_count: bigbed_read_u64_local(data, offset, endian),
    min_val: bigwig_read_f64(data, offset + 8, endian),
    max_val: bigwig_read_f64(data, offset + 16, endian),
    sum_data: bigwig_read_f64(data, offset + 24, endian),
    sum_squares: bigwig_read_f64(data, offset + 32, endian),
  }
}

///|
pub fn bigwig_parse(data : Array[Int]) -> BigWigFile raise BigWigError {
  bigwig_wrap(fn() -> BigWigFile raise BigBedError {
    let (header, endian) = bigwig_read_header(data)
    let trailing = bigbed_read_u32(data, data.length() - 4, endian)
    if trailing != bigwig_magic {
      bigbed_fail("BigWig trailing magic number is missing")
    }
    let zoom_headers = bigwig_read_zoom_headers(
      data,
      header.zoom_levels,
      endian,
    )
    let summary = bigwig_read_summary(data, header.total_summary_offset, endian)
    let targets = bigbed_read_chrom_tree(
      data,
      header.chromosome_tree_offset,
      endian,
    )
    let section_count = bigbed_read_u64_local(
      data,
      header.full_data_offset,
      endian,
    )
    let (blocks, indexed_count) = bigbed_read_rtree(
      data,
      header.full_index_offset,
      endian,
    )
    if indexed_count != section_count {
      bigbed_fail("BigWig data and R-tree section counts disagree")
    }
    BigWigFile::{
      header,
      targets,
      source: data.copy(),
      blocks,
      zoom_headers,
      total_summary: summary,
      section_count,
    }
  })
}

// ---------------------------------------------------------------------------
// Section decoding
// ---------------------------------------------------------------------------

///|
/// Decode one (possibly compressed) data section into intervals.
fn bigwig_decode_section(
  file : BigWigFile,
  block : BigBedBlockIndex,
) -> Array[BigWigInterval] raise BigBedError {
  let endian = if file.header.little_endian {
    BigBedLittleEndian
  } else {
    BigBedBigEndian
  }
  let payload : Array[Int] = []
  for i in block.data_offset..<(block.data_offset + block.data_size) {
    payload.push(file.source[i])
  }
  let raw = if file.header.uncompress_buffer_size > 0 {
    bigbed_zlib_unstore(payload, file.header.uncompress_buffer_size)
  } else if payload.length() >= 2 && payload[0] == 0x78 {
    // Some writers compress blocks while leaving the buffer size at zero.
    bigbed_zlib_unstore(payload, block.data_size * 256 + 65536)
  } else {
    payload
  }
  bigbed_check_range(raw, 0, 24, "section header")
  let chrom_id = bigbed_read_u32(raw, 0, endian)
  let chrom_start = bigbed_read_u32(raw, 4, endian)
  let chrom_end = bigbed_read_u32(raw, 8, endian)
  let item_step = bigbed_read_u32(raw, 12, endian)
  let item_span = bigbed_read_u32(raw, 16, endian)
  let section_type = raw[20]
  let item_count = bigbed_read_u16(raw, 22, endian)
  if chrom_id < 0 || chrom_id >= file.targets.length() {
    bigbed_fail("BigWig section references an unknown chromosome id")
  }
  let chrom_name = file.targets[chrom_id].name
  let chrom_length = file.targets[chrom_id].length
  if chrom_start < 0 || chrom_end > chrom_length || chrom_start > chrom_end {
    bigbed_fail("BigWig section bounds exceed the chromosome length")
  }
  let intervals : Array[BigWigInterval] = []
  if section_type == bigwig_bed_graph {
    let need = item_count * 12
    bigbed_check_range(raw, 24, need, "bedGraph section")
    for i in 0..<item_count {
      let position = 24 + i * 12
      let start = bigbed_read_u32(raw, position, endian)
      let end = bigbed_read_u32(raw, position + 4, endian)
      let value = bigwig_read_f32(raw, position + 8, endian)
      if start < chrom_start || end > chrom_end || start >= end {
        bigbed_fail("BigWig bedGraph item is outside its section bounds")
      }
      intervals.push(BigWigInterval::{ chrom: chrom_name, start, end, value, })
    }
  } else if section_type == bigwig_variable_step {
    if item_span == 0 {
      bigbed_fail("BigWig variableStep section is missing its span")
    }
    let need = item_count * 8
    bigbed_check_range(raw, 24, need, "variableStep section")
    for i in 0..<item_count {
      let position = 24 + i * 8
      let start = bigbed_read_u32(raw, position, endian)
      let value = bigwig_read_f32(raw, position + 4, endian)
      let end = start + item_span
      if start < chrom_start || end > chrom_end {
        bigbed_fail("BigWig variableStep item is outside its section bounds")
      }
      intervals.push(BigWigInterval::{ chrom: chrom_name, start, end, value, })
    }
  } else if section_type == bigwig_fixed_step {
    if item_step == 0 || item_span == 0 {
      bigbed_fail("BigWig fixedStep section is missing its step or span")
    }
    let need = item_count * 4
    bigbed_check_range(raw, 24, need, "fixedStep section")
    for i in 0..<item_count {
      let value = bigwig_read_f32(raw, 24 + i * 4, endian)
      let start = chrom_start + i * item_step
      let end = start + item_span
      if end > chrom_end {
        bigbed_fail("BigWig fixedStep item is outside its section bounds")
      }
      intervals.push(BigWigInterval::{ chrom: chrom_name, start, end, value, })
    }
  } else {
    bigbed_fail("Unsupported BigWig section type " + section_type.to_string())
  }
  intervals
}

///|
fn bigwig_target_id(file : BigWigFile, chrom : String) -> Int raise BigBedError {
  for i in 0..<file.targets.length() {
    if file.targets[i].name == chrom {
      return i
    }
  }
  raise BigBedError("Unknown BigWig chromosome '" + chrom + "'")
}

///|
fn bigwig_query_bounds(
  file : BigWigFile,
  chrom : String,
  start : Int,
  end_ : Int?,
) -> (Int, Int, Int) raise BigBedError {
  let chrom_id = bigwig_target_id(file, chrom)
  let length = file.targets[chrom_id].length
  let qend = match end_ {
    Some(value) => value
    None => length
  }
  if start < 0 || qend > length || start >= qend {
    bigbed_fail(
      "BigWig query range [" +
      start.to_string() +
      "," +
      qend.to_string() +
      ") is outside chromosome '" +
      chrom +
      "' (length " +
      length.to_string() +
      ")",
    )
  }
  (chrom_id, start, qend)
}

///|
/// All intervals overlapping the query, in section/item order. The intervals
/// are returned unclipped, matching libBigWig/pyBigWig `intervals()`.
pub fn BigWigFile::intervals(
  self : BigWigFile,
  chrom : String,
  start? : Int = 0,
  end_? : Int? = None,
) -> Array[BigWigInterval] raise BigWigError {
  bigwig_wrap(fn() -> Array[BigWigInterval] raise BigBedError {
    let (chrom_id, qstart, qend) = bigwig_query_bounds(self, chrom, start, end_)
    let result : Array[BigWigInterval] = []
    for block in self.blocks {
      if bigbed_block_overlaps(block, chrom_id, qstart, qend) {
        for interval in bigwig_decode_section(self, block) {
          if interval.start < qend && qstart < interval.end {
            result.push(interval)
          }
        }
      }
    }
    result
  })
}

///|
/// Per-base values over [start, end); uncovered bases are NaN. When intervals
/// overlap, the later section/item wins (libBigWig semantics).
pub fn BigWigFile::values(
  self : BigWigFile,
  chrom : String,
  start? : Int = 0,
  end_? : Int? = None,
) -> Array[Double] raise BigWigError {
  bigwig_wrap(fn() -> Array[Double] raise BigBedError {
    let (chrom_id, qstart, qend) = bigwig_query_bounds(self, chrom, start, end_)
    let nan = 0.0 / 0.0
    let values = Array::make(qend - qstart, nan)
    for block in self.blocks {
      if bigbed_block_overlaps(block, chrom_id, qstart, qend) {
        for interval in bigwig_decode_section(self, block) {
          let lo = if interval.start > qstart { interval.start } else { qstart }
          let hi = if interval.end < qend { interval.end } else { qend }
          for base in lo..<hi {
            values[base - qstart] = interval.value
          }
        }
      }
    }
    values
  })
}

// ---------------------------------------------------------------------------
// Statistics
// ---------------------------------------------------------------------------

///|
priv struct BigWigBin {
  mut count : Double
  mut min_val : Double
  mut max_val : Double
  mut sum_data : Double
  mut sum_squares : Double
}

///|
fn bigwig_empty_bins(n_bins : Int) -> Array[BigWigBin] {
  let bins : Array[BigWigBin] = []
  for _ in 0..<n_bins {
    bins.push(BigWigBin::{
      count: 0,
      min_val: 0.0,
      max_val: 0.0,
      sum_data: 0.0,
      sum_squares: 0.0,
    })
  }
  bins
}

///|
fn bigwig_bin_add(bin : BigWigBin, value : Double, bases : Int) -> Unit {
  let weighted = value * bases.to_double()
  bin.sum_data = bin.sum_data + weighted
  bin.sum_squares = bin.sum_squares + value * weighted
  if bin.count == 0 {
    bin.min_val = value
    bin.max_val = value
  } else {
    if value < bin.min_val {
      bin.min_val = value
    }
    if value > bin.max_val {
      bin.max_val = value
    }
  }
  bin.count = bin.count + bases.to_double()
}

///|
/// Cumulative-floor bin boundary: bin i spans [bound(i), bound(i+1)) within
/// a query of `total` bases split into `n_bins` bins (pyBigWig semantics).
fn bigwig_bin_bound(total : Int, n_bins : Int, index : Int) -> Int {
  let scaled = index.to_double() * total.to_double() / n_bins.to_double()
  scaled.floor().to_int()
}

///|
/// Inverse boundary: the bin index containing a position `relative` bases
/// into the query.
fn bigwig_bin_index(total : Int, n_bins : Int, relative : Int) -> Int {
  let scaled = relative.to_double() * n_bins.to_double() / total.to_double()
  scaled.floor().to_int()
}

///|
fn bigwig_stat_value(bin : BigWigBin, bin_span : Int, stat : String) -> Double? {
  if stat == "cov" {
    return Some(bin.count / bin_span.to_double())
  }
  if bin.count == 0.0 {
    return None
  }
  let mean = bin.sum_data / bin.count
  match stat {
    "mean" => Some(mean)
    "sum" => Some(bin.sum_data)
    "min" => Some(bin.min_val)
    "max" => Some(bin.max_val)
    "std" =>
      // Sample standard deviation with covered bases as the sample size,
      // matching libBigWig/pyBigWig; single-base samples return 0.
      if bin.count < 2.0 {
        Some(0.0)
      } else {
        let variance = (bin.sum_squares - bin.sum_data * mean) /
          (bin.count - 1.0)
        Some(if variance < 0.0 { 0.0 } else { variance.sqrt() })
      }
    _ => None
  }
}

///|
fn bigwig_raw_stats(
  file : BigWigFile,
  chrom_id : Int,
  qstart : Int,
  qend : Int,
  n_bins : Int,
  stat : String,
) -> Array[Double?] raise BigBedError {
  let bins = bigwig_empty_bins(n_bins)
  let total = qend - qstart
  for block in file.blocks {
    if bigbed_block_overlaps(block, chrom_id, qstart, qend) {
      for interval in bigwig_decode_section(file, block) {
        if interval.start >= qend || qstart >= interval.end {
          continue
        }
        let lo = if interval.start > qstart { interval.start } else { qstart }
        let hi = if interval.end < qend { interval.end } else { qend }
        let mut first_bin = bigwig_bin_index(total, n_bins, lo - qstart)
        let last_bin = bigwig_bin_index(total, n_bins, hi - 1 - qstart)
        if first_bin >= n_bins {
          first_bin = n_bins - 1
        }
        let bin_limit = last_bin + 1
        for bin_index in first_bin..<bin_limit {
          let bstart = qstart + bigwig_bin_bound(total, n_bins, bin_index)
          let bend = qstart + bigwig_bin_bound(total, n_bins, bin_index + 1)
          let overlap_start = if lo > bstart { lo } else { bstart }
          let overlap_end = if hi < bend { hi } else { bend }
          if overlap_end > overlap_start {
            bigwig_bin_add(
              bins[bin_index],
              interval.value,
              overlap_end - overlap_start,
            )
          }
        }
      }
    }
  }
  let result : Array[Double?] = []
  for i in 0..<n_bins {
    let span = bigwig_bin_bound(total, n_bins, i + 1) -
      bigwig_bin_bound(total, n_bins, i)
    result.push(bigwig_stat_value(bins[i], span, stat))
  }
  result
}

///|
priv struct BigWigZoomElement {
  chrom_id : Int
  start : Int
  end : Int
  valid_count : Int
  min_val : Double
  max_val : Double
  sum_data : Double
  sum_squares : Double
}

///|
fn bigwig_decode_zoom_block(
  file : BigWigFile,
  data_offset : Int,
  data_size : Int,
) -> Array[BigWigZoomElement] raise BigBedError {
  let endian = if file.header.little_endian {
    BigBedLittleEndian
  } else {
    BigBedBigEndian
  }
  let payload : Array[Int] = []
  for i in data_offset..<(data_offset + data_size) {
    payload.push(file.source[i])
  }
  let cap = if file.header.uncompress_buffer_size > 0 {
    file.header.uncompress_buffer_size
  } else {
    data_size * 256 + 65536
  }
  let raw = bigbed_zlib_unstore(payload, cap)
  let elements : Array[BigWigZoomElement] = []
  let mut cursor = 0
  while cursor + 32 <= raw.length() {
    let chrom_id = bigbed_read_u32(raw, cursor, endian)
    let start = bigbed_read_u32(raw, cursor + 4, endian)
    let end = bigbed_read_u32(raw, cursor + 8, endian)
    let valid_count = bigbed_read_u32(raw, cursor + 12, endian)
    elements.push(BigWigZoomElement::{
      chrom_id,
      start,
      end,
      valid_count,
      min_val: bigwig_read_f32(raw, cursor + 16, endian),
      max_val: bigwig_read_f32(raw, cursor + 20, endian),
      sum_data: bigwig_read_f32(raw, cursor + 24, endian),
      sum_squares: bigwig_read_f32(raw, cursor + 28, endian),
    })
    cursor = cursor + 32
  }
  if cursor != raw.length() {
    bigbed_fail("BigWig zoom block is not a whole number of elements")
  }
  elements
}

///|
/// Read the R-tree leaves of a zoom index and decode every zoom payload.
fn bigwig_zoom_elements(
  file : BigWigFile,
  zoom : BigWigZoomHeader,
) -> Array[BigWigZoomElement] raise BigBedError {
  let endian = if file.header.little_endian {
    BigBedLittleEndian
  } else {
    BigBedBigEndian
  }
  let leaves : Array[BigBedBlockIndex] = []
  bigbed_read_rtree_node(file.source, zoom.index_offset + 48, endian, leaves, 0)
  let elements : Array[BigWigZoomElement] = []
  // Zoom data starts with a u32 counting the compressed data blocks (i.e. the
  // R-tree leaves); each block may hold several 32-byte zoom records.
  for leaf in leaves {
    let decoded = bigwig_decode_zoom_block(
      file,
      leaf.data_offset,
      leaf.data_size,
    )
    for element in decoded {
      elements.push(element)
    }
  }
  let prefix_count = bigbed_read_u32(file.source, zoom.data_offset, endian)
  if prefix_count != leaves.length() {
    bigbed_fail("BigWig zoom data block count does not match its index")
  }
  elements
}

///|
fn bigwig_choose_zoom(file : BigWigFile, bin_width : Int) -> BigWigZoomHeader? {
  let mut selected : BigWigZoomHeader? = None
  for zoom in file.zoom_headers {
    if zoom.reduction_level <= bin_width {
      selected = Some(zoom)
    }
  }
  selected
}

///|
fn bigwig_zoom_stats(
  file : BigWigFile,
  zoom : BigWigZoomHeader,
  chrom_id : Int,
  qstart : Int,
  qend : Int,
  n_bins : Int,
  stat : String,
) -> Array[Double?] raise BigBedError {
  let elements = bigwig_zoom_elements(file, zoom)
  let bins = bigwig_empty_bins(n_bins)
  let total = qend - qstart
  for element in elements {
    if element.chrom_id != chrom_id ||
      element.start >= qend ||
      qstart >= element.end {
      continue
    }
    let lo = if element.start > qstart { element.start } else { qstart }
    let hi = if element.end < qend { element.end } else { qend }
    if lo >= hi {
      continue
    }
    let mut first_bin = bigwig_bin_index(total, n_bins, lo - qstart)
    let last_bin = bigwig_bin_index(total, n_bins, hi - 1 - qstart)
    if first_bin >= n_bins {
      first_bin = n_bins - 1
    }
    let bin_limit = last_bin + 1
    for bin_index in first_bin..<bin_limit {
      let bstart = qstart + bigwig_bin_bound(total, n_bins, bin_index)
      let bend = qstart + bigwig_bin_bound(total, n_bins, bin_index + 1)
      let overlap_start = if lo > bstart { lo } else { bstart }
      let overlap_end = if hi < bend { hi } else { bend }
      if overlap_end <= overlap_start {
        continue
      }
      // libBigWig weights each zoom element by the fraction overlapping the
      // bin, including covered-base counts.
      let fraction = (overlap_end - overlap_start).to_double() /
        (element.end - element.start).to_double()
      let bin = bins[bin_index]
      let was_empty = bin.count == 0.0
      bin.count = bin.count + element.valid_count.to_double() * fraction
      bin.sum_data = bin.sum_data + element.sum_data * fraction
      bin.sum_squares = bin.sum_squares + element.sum_squares * fraction
      if was_empty {
        bin.min_val = element.min_val
        bin.max_val = element.max_val
      } else {
        if element.min_val < bin.min_val {
          bin.min_val = element.min_val
        }
        if element.max_val > bin.max_val {
          bin.max_val = element.max_val
        }
      }
    }
  }
  let result : Array[Double?] = []
  for i in 0..<n_bins {
    let span = bigwig_bin_bound(total, n_bins, i + 1) -
      bigwig_bin_bound(total, n_bins, i)
    result.push(bigwig_stat_value(bins[i], span, stat))
  }
  result
}

///|
/// Binned statistics. `stat` is one of mean/min/max/cov/sum/std (libBigWig
/// names). Zoom data is used when a zoom level exists whose reduction is no
/// wider than a bin; otherwise statistics are computed exactly from raw
/// sections. Empty bins return None (coverage returns 0.0).
pub fn BigWigFile::stats(
  self : BigWigFile,
  chrom : String,
  start? : Int = 0,
  end_? : Int? = None,
  n_bins? : Int = 1,
  stat? : String = "mean",
) -> Array[Double?] raise BigWigError {
  bigwig_wrap(fn() -> Array[Double?] raise BigBedError {
    let (chrom_id, qstart, qend) = bigwig_query_bounds(self, chrom, start, end_)
    if n_bins < 1 {
      bigbed_fail("BigWig stats require at least one bin")
    }
    if stat != "mean" &&
      stat != "min" &&
      stat != "max" &&
      stat != "cov" &&
      stat != "sum" &&
      stat != "std" {
      bigbed_fail("Unknown BigWig statistic '" + stat + "'")
    }
    let bin_width = (qend - qstart) / n_bins
    match bigwig_choose_zoom(self, bin_width) {
      Some(zoom) =>
        bigwig_zoom_stats(self, zoom, chrom_id, qstart, qend, n_bins, stat)
      None => bigwig_raw_stats(self, chrom_id, qstart, qend, n_bins, stat)
    }
  })
}

// ---------------------------------------------------------------------------
// Accessors
// ---------------------------------------------------------------------------

///|
pub fn BigWigFile::chroms(self : BigWigFile) -> Array[(String, Int)] {
  let result : Array[(String, Int)] = []
  for target in self.targets {
    result.push((target.name, target.length))
  }
  result
}

///|
pub fn BigWigFile::header(self : BigWigFile) -> BigWigHeader {
  self.header
}

///|
pub fn BigWigFile::total_summary(self : BigWigFile) -> BigWigSummary {
  self.total_summary
}

///|
pub fn BigWigFile::zoom_levels(self : BigWigFile) -> Array[BigWigZoomHeader] {
  self.zoom_headers.copy()
}

///|
pub fn BigWigFile::section_count(self : BigWigFile) -> Int {
  self.section_count
}

// ---------------------------------------------------------------------------
// Writer
// ---------------------------------------------------------------------------

///|
fn bigwig_validate_inputs(
  targets : Array[BigBedTarget],
  intervals : Array[BigWigInterval],
) -> Array[Int] raise BigBedError {
  if targets.length() == 0 {
    bigbed_fail("BigWig output requires at least one chromosome")
  }
  for i in 0..<targets.length() {
    if targets[i].name.length() == 0 {
      bigbed_fail("BigWig chromosome names must be non-empty")
    }
    if targets[i].length <= 0 {
      bigbed_fail("BigWig chromosome lengths must be positive")
    }
    for j in (i + 1)..<targets.length() {
      if targets[i].name == targets[j].name {
        bigbed_fail(
          "Duplicate BigWig chromosome name '" + targets[i].name + "'",
        )
      }
    }
  }
  for interval in intervals {
    if bigbed_target_index(targets, interval.chrom) < 0 {
      bigbed_fail("Unknown BigWig chromosome '" + interval.chrom + "'")
    }
  }
  // Resolve and sort every interval by (chrom id, start).
  let order : Array[Int] = []
  for i in 0..<intervals.length() {
    order.push(i)
  }
  order.sort_by(fn(left : Int, right : Int) -> Int {
    let a = intervals[left]
    let b = intervals[right]
    let ia = bigbed_target_index(targets, a.chrom)
    let ib = bigbed_target_index(targets, b.chrom)
    if ia != ib {
      ia - ib
    } else if a.start != b.start {
      a.start - b.start
    } else {
      a.end - b.end
    }
  })
  let mut previous_chrom = -1
  let mut previous_end = 0
  for index in order {
    let interval = intervals[index]
    let chrom_id = bigbed_target_index(targets, interval.chrom)
    let length = targets[chrom_id].length
    if interval.start < 0 ||
      interval.end > length ||
      interval.start >= interval.end {
      bigbed_fail(
        "BigWig interval " +
        interval.chrom +
        ":[" +
        interval.start.to_string() +
        "," +
        interval.end.to_string() +
        ") is outside the chromosome (length " +
        length.to_string() +
        ")",
      )
    }
    if interval.value.is_nan() || interval.value.abs() > 3.4028235e38 {
      bigbed_fail(
        "BigWig value at " +
        interval.chrom +
        ":" +
        interval.start.to_string() +
        " is not representable as f32",
      )
    }
    if chrom_id != previous_chrom {
      previous_chrom = chrom_id
      previous_end = interval.end
    } else {
      if interval.start < previous_end {
        bigbed_fail(
          "Overlapping BigWig intervals on chromosome '" +
          interval.chrom +
          "' at " +
          interval.start.to_string(),
        )
      }
      previous_end = interval.end
    }
  }
  order
}

///|
fn bigwig_section_bytes(
  intervals : Array[BigWigInterval],
  start_index : Int,
  end_index : Int,
  chrom_id : Int,
) -> Array[Int] {
  let first = intervals[start_index]
  let mut section_end = first.end
  for i in (start_index + 1)..<end_index {
    if intervals[i].end > section_end {
      section_end = intervals[i].end
    }
  }
  let raw : Array[Int] = []
  bigbed_push_u32(raw, chrom_id)
  bigbed_push_u32(raw, first.start)
  bigbed_push_u32(raw, section_end)
  bigbed_push_u32(raw, 0)
  bigbed_push_u32(raw, 0)
  raw.push(bigwig_bed_graph)
  raw.push(0)
  bigbed_push_u16(raw, end_index - start_index)
  for i in start_index..<end_index {
    let interval = intervals[i]
    bigbed_push_u32(raw, interval.start)
    bigbed_push_u32(raw, interval.end)
    bigwig_push_f32(raw, interval.value)
  }
  raw
}

///|
fn bigwig_patch_header(
  output : Array[Int],
  chromosome_tree_offset : Int,
  full_data_offset : Int,
  full_index_offset : Int,
  total_summary_offset : Int,
  uncompress_buffer_size : Int,
) -> Unit {
  bigbed_patch_u32(output, 0, bigwig_magic)
  bigbed_patch_u16(output, 4, 4)
  bigbed_patch_u16(output, 6, 0)
  bigbed_patch_u64_local(output, 8, chromosome_tree_offset)
  bigbed_patch_u64_local(output, 16, full_data_offset)
  bigbed_patch_u64_local(output, 24, full_index_offset)
  bigbed_patch_u16(output, 32, 0)
  bigbed_patch_u16(output, 34, 0)
  bigbed_patch_u64_local(output, 36, 0)
  bigbed_patch_u64_local(output, 44, total_summary_offset)
  bigbed_patch_u32(output, 52, uncompress_buffer_size)
  bigbed_patch_u64_local(output, 56, 0)
}

///|
pub fn bigwig_write(
  targets : Array[BigBedTarget],
  intervals : Array[BigWigInterval],
  config? : BigWigWriteConfig = BigWigWriteConfig::default(),
) -> Array[Int] raise BigWigError {
  bigwig_wrap(fn() -> Array[Int] raise BigBedError {
    let order = bigwig_validate_inputs(targets, intervals)
    let sorted : Array[BigWigInterval] = []
    for index in order {
      sorted.push(intervals[index])
    }
    let output : Array[Int] = Array::make(64, 0)
    let total_summary_offset = output.length()
    for _ in 0..<40 {
      output.push(0)
    }
    let chrom_block_size = if targets.length() < config.block_size {
      targets.length()
    } else {
      config.block_size
    }
    let chromosome_tree_offset = output.length()
    bigbed_write_chrom_tree(output, targets, chrom_block_size)
    // Build every data section first so the u64 section count can precede the
    // payloads in the output stream.
    // (chrom id, start base, end base, encoded payload) per section.
    let sections : Array[(Int, Int, Int, Array[Int])] = []
    let mut maximum_uncompressed = 0
    let mut index = 0
    while index < sorted.length() {
      let chrom_id = bigbed_target_index(targets, sorted[index].chrom)
      let mut end_index = index
      while end_index < sorted.length() &&
            bigbed_target_index(targets, sorted[end_index].chrom) == chrom_id &&
            end_index - index < bigwig_max_bedgraph_items {
        end_index = end_index + 1
      }
      let raw = bigwig_section_bytes(sorted, index, end_index, chrom_id)
      if raw.length() > maximum_uncompressed {
        maximum_uncompressed = raw.length()
      }
      let payload = if config.compress { bigbed_zlib_store(raw) } else { raw }
      let mut section_end = sorted[index].end
      for i in (index + 1)..<end_index {
        if sorted[i].end > section_end {
          section_end = sorted[i].end
        }
      }
      sections.push((chrom_id, sorted[index].start, section_end, payload))
      index = end_index
    }
    let blocks : Array[BigBedBlockIndex] = []
    let full_data_offset = output.length()
    bigbed_push_u64_local(output, sections.length())
    for section in sections {
      let (section_chrom, section_start, section_end, payload) = section
      let data_offset = output.length()
      bigbed_append(output, payload)
      blocks.push(BigBedBlockIndex::{
        start_chrom: section_chrom,
        start_base: section_start,
        end_chrom: section_chrom,
        end_base: section_end,
        data_offset,
        data_size: payload.length(),
      })
    }
    // Total summary (u64 count + four f64 values).
    let mut valid_count = 0
    let mut min_val = 0.0
    let mut max_val = 0.0
    let mut has_value = false
    let mut sum_data = 0.0
    let mut sum_squares = 0.0
    for interval in sorted {
      let width = interval.end - interval.start
      valid_count = valid_count + width
      if !has_value {
        min_val = interval.value
        max_val = interval.value
        has_value = true
      } else {
        if interval.value < min_val {
          min_val = interval.value
        }
        if interval.value > max_val {
          max_val = interval.value
        }
      }
      sum_data = sum_data + interval.value * width.to_double()
      sum_squares = sum_squares +
        interval.value * interval.value * width.to_double()
    }
    bigbed_patch_u64_local(output, total_summary_offset, valid_count)
    let summary_bytes : Array[Int] = []
    bigwig_push_f64(summary_bytes, min_val)
    bigwig_push_f64(summary_bytes, max_val)
    bigwig_push_f64(summary_bytes, sum_data)
    bigwig_push_f64(summary_bytes, sum_squares)
    for i in 0..<32 {
      output[total_summary_offset + 8 + i] = summary_bytes[i]
    }
    let full_index_offset = output.length()
    bigbed_write_rtree(
      output,
      blocks,
      config.block_size,
      blocks.length(),
      1,
      52,
    )
    bigbed_push_u32(output, bigwig_magic)
    bigwig_patch_header(
      output,
      chromosome_tree_offset,
      full_data_offset,
      full_index_offset,
      total_summary_offset,
      if config.compress && maximum_uncompressed > 0 {
        maximum_uncompressed
      } else {
        0
      },
    )
    output
  })
}
