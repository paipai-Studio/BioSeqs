///|
/// Tabix (.tbi) index parsing and region querying - the random-access index
/// for coordinate-sorted, bgzip-compressed tab-delimited files (VCF, GFF,
/// BED, ...), as produced by `tabix`.
///
/// Port of the Rsamtools / htslib tabix reader. Tabix uses the same binning
/// scheme as the BAM index (BAI): a hierarchy of fixed-size bins plus a
/// linear index every 16384 bp. The per-reference binning structures are
/// therefore shared with `src/bam_index.mbt`.

///|
pub suberror TabixError {
  TabixError(String)
}

///|
fn[T] tbx_fail(msg : String) -> T raise TabixError {
  raise TabixError(msg)
}

///|
/// A parsed Tabix index.
///
/// `format` is 0 (generic), 1 (SAM) or 2 (VCF); `col_seq` / `col_beg` /
/// `col_end` are the 1-based file columns for sequence, start and end
/// (`col_end` is 0 when it defaults to `col_beg`); `meta` is the comment
/// character; `skip` the number of header lines; `names` the reference
/// sequence names. The per-reference binning is identical to BAI.
pub struct TabixIndex {
  n_ref : Int
  format : Int
  col_seq : Int
  col_beg : Int
  col_end : Int
  meta : Int
  skip : Int
  names : Array[String]
  references : Array[BaiReference]
}

// ---------------------------------------------------------------------------
// Little-endian binary readers (mirror src/bam_index.mbt, raising TabixError)
// ---------------------------------------------------------------------------

///|
fn tbx_byte(data : String, pos : Int) -> Int64 {
  (data.unsafe_get(pos).to_int() & 0xFF).to_int64()
}

///|
fn tbx_check_range(
  data : String,
  pos : Int,
  len : Int,
  what : String,
) -> Unit raise TabixError {
  if pos < 0 || pos + len > data.length() {
    tbx_fail(
      "Cannot read \{what} at offset \{pos}: premature end of data (length \{data.length()})",
    )
  }
}

///|
fn tbx_read_u32(data : String, pos : Int) -> Int raise TabixError {
  tbx_check_range(data, pos, 4, "32-bit integer")
  (data.unsafe_get(pos).to_int() & 0xFF) |
    ((data.unsafe_get(pos + 1).to_int() & 0xFF) << 8) |
    ((data.unsafe_get(pos + 2).to_int() & 0xFF) << 16) |
    ((data.unsafe_get(pos + 3).to_int() & 0xFF) << 24)
}

///|
fn tbx_read_u64(data : String, pos : Int) -> Int64 raise TabixError {
  tbx_check_range(data, pos, 8, "64-bit integer")
  tbx_byte(data, pos) |
    (tbx_byte(data, pos + 1) << 8) |
    (tbx_byte(data, pos + 2) << 16) |
    (tbx_byte(data, pos + 3) << 24) |
    (tbx_byte(data, pos + 4) << 32) |
    (tbx_byte(data, pos + 5) << 40) |
    (tbx_byte(data, pos + 6) << 48) |
    (tbx_byte(data, pos + 7) << 56)
}

///|
fn tbx_check_magic(data : String) -> Bool {
  if data.length() < 4 {
    return false
  }
  data.unsafe_get(0).to_int() == 0x54 &&
    data.unsafe_get(1).to_int() == 0x42 &&
    data.unsafe_get(2).to_int() == 0x49 &&
    data.unsafe_get(3).to_int() == 0x01
}

///|
/// Split the NUL-terminated sequence names (total `l_nm` bytes).
fn tbx_parse_names(data : String, pos : Int, l_nm : Int) -> Array[String] {
  let names : Array[String] = []
  let mut cur = StringBuilder::new()
  for i = 0; i < l_nm; i = i + 1 {
    let byte = data.unsafe_get(pos + i).to_int() & 0xFF
    if byte == 0 {
      names.push(cur.to_string())
      cur = StringBuilder::new()
    } else {
      cur.write_char(data.unsafe_get(pos + i).unsafe_to_char())
    }
  }
  names
}

// ---------------------------------------------------------------------------
// Parsing and querying
// ---------------------------------------------------------------------------

///|
/// Parse a Tabix (.tbi) index from its raw bytes.
pub fn parse_tabix(data : String) -> TabixIndex raise TabixError {
  if !tbx_check_magic(data) {
    tbx_fail("Not a Tabix index: bad magic bytes")
  }
  let n_ref = tbx_read_u32(data, 4)
  let format = tbx_read_u32(data, 8)
  let col_seq = tbx_read_u32(data, 12)
  let col_beg = tbx_read_u32(data, 16)
  let col_end = tbx_read_u32(data, 20)
  let meta = tbx_read_u32(data, 24)
  let skip = tbx_read_u32(data, 28)
  let l_nm = tbx_read_u32(data, 32)
  let names = tbx_parse_names(data, 36, l_nm)
  let references : Array[BaiReference] = []
  let mut pos = 36 + l_nm
  for r = 0; r < n_ref; r = r + 1 {
    let n_bin = tbx_read_u32(data, pos)
    pos = pos + 4
    let bins : Map[Int, Array[BaiChunk]] = Map([], capacity=n_bin)
    for b = 0; b < n_bin; b = b + 1 {
      let bin = tbx_read_u32(data, pos)
      pos = pos + 4
      let n_chunk = tbx_read_u32(data, pos)
      pos = pos + 4
      let chunks : Array[BaiChunk] = []
      for c = 0; c < n_chunk; c = c + 1 {
        let beg = tbx_read_u64(data, pos)
        let end = tbx_read_u64(data, pos + 8)
        pos = pos + 16
        chunks.push(BaiChunk::{ chunk_beg: beg, chunk_end: end })
      }
      bins.set(bin, chunks)
    }
    let n_intv = tbx_read_u32(data, pos)
    pos = pos + 4
    let linear : Array[Int64] = []
    for iv = 0; iv < n_intv; iv = iv + 1 {
      linear.push(tbx_read_u64(data, pos))
      pos = pos + 8
    }
    references.push(BaiReference::{ bins, linear })
  }
  TabixIndex::{
    n_ref,
    format,
    col_seq,
    col_beg,
    col_end,
    meta,
    skip,
    names,
    references,
  }
}

///|
/// Find the chunks of a Tabix index that overlap a half-open `[beg, end)`
/// interval on reference `ref_id`.
///
/// Returns the linear-index minimum offset together with every overlapping
/// chunk (BGZF virtual offsets), ready for `BgzfReader::seek`.
pub fn tabix_query(
  index : TabixIndex,
  ref_id : Int,
  beg : Int,
  end : Int,
) -> BaiQueryResult raise TabixError {
  if ref_id < 0 || ref_id >= index.references.length() {
    tbx_fail("Reference id \{ref_id} is out of range (n_ref = \{index.n_ref})")
  }
  if beg < 0 || end < beg {
    tbx_fail("Invalid query interval [\{beg}, \{end})")
  }
  let reference = index.references[ref_id]
  let min_off = if reference.linear.length() > 0 {
    let iv = beg >> 14
    if iv < reference.linear.length() {
      reference.linear[iv]
    } else {
      reference.linear[reference.linear.length() - 1]
    }
  } else {
    0L
  }
  let chunks : Array[BaiChunk] = []
  let bins = bai_reg2bins(beg, end)
  for bin in bins {
    match reference.bins.get(bin) {
      Some(cs) =>
        for c in cs {
          chunks.push(c)
        }
      None => ()
    }
  }
  BaiQueryResult::{ min_off, chunks }
}