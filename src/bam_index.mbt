///|
/// BAI (BAM index) parsing and region querying - the random-access index for
/// BAM files, as produced by `samtools index`.
///
/// Port of the Rsamtools / htslib binning index. BAI shares its binning
/// scheme with tabix: a hierarchy of fixed-size bins (top bin 0 down to
/// 16 kb bins) plus a linear index sampled every 16384 bp. `parse_bai`
/// decodes the on-disk structure; `bai_query` returns the BGZF virtual
/// offset chunks that overlap a half-open `[beg, end)` interval, ready to be
/// fed to `BgzfReader::seek` (src/bgzf.mbt) for real random access.

///|
pub suberror BaiIndexError {
  BaiIndexError(String)
}

///|
fn[T] bai_fail(msg : String) -> T raise BaiIndexError {
  raise BaiIndexError(msg)
}

///|
/// A contiguous chunk of the BAM file, given as a pair of BGZF virtual
/// offsets (block offset in the upper 48 bits, within-block offset in the
/// lower 16 bits).
pub struct BaiChunk {
  chunk_beg : Int64
  chunk_end : Int64
}

///|
/// Per-reference binning data: a map of bin number to its chunk list, plus
/// the linear index (smallest virtual offset every 16 kb).
pub struct BaiReference {
  bins : Map[Int, Array[BaiChunk]]
  linear : Array[Int64]
}

///|
/// A parsed BAI index.
///
/// `n_no_coor` is the number of records without coordinates (the trailing
/// `uint64`), when present.
pub struct BaiIndex {
  n_ref : Int
  references : Array[BaiReference]
  n_no_coor : Int64?
}

///|
/// The result of a region query: the linear-index minimum offset plus every
/// overlapping chunk.
pub struct BaiQueryResult {
  min_off : Int64
  chunks : Array[BaiChunk]
}

///|
/// The list of bin numbers that may contain a half-open interval
/// `[beg, end)`. Mirrors htslib's `reg2bins` (Heng Li): bin 0 is always
/// included, followed by the five fixed-size levels.
pub fn bai_reg2bins(beg : Int, end : Int) -> Array[Int] {
  let bins : Array[Int] = []
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

// ---------------------------------------------------------------------------
// Little-endian binary readers
// ---------------------------------------------------------------------------

///|
fn bai_byte(data : String, pos : Int) -> Int64 {
  (data.unsafe_get(pos).to_int() & 0xFF).to_int64()
}

///|
fn bai_check_range(
  data : String,
  pos : Int,
  len : Int,
  what : String,
) -> Unit raise BaiIndexError {
  if pos < 0 || pos + len > data.length() {
    bai_fail(
      "Cannot read \{what} at offset \{pos}: premature end of data (length \{data.length()})",
    )
  }
}

///|
fn bai_read_u32(data : String, pos : Int) -> Int raise BaiIndexError {
  bai_check_range(data, pos, 4, "32-bit integer")
  (data.unsafe_get(pos).to_int() & 0xFF) |
    ((data.unsafe_get(pos + 1).to_int() & 0xFF) << 8) |
    ((data.unsafe_get(pos + 2).to_int() & 0xFF) << 16) |
    ((data.unsafe_get(pos + 3).to_int() & 0xFF) << 24)
}

///|
fn bai_read_u64(data : String, pos : Int) -> Int64 raise BaiIndexError {
  bai_check_range(data, pos, 8, "64-bit integer")
  bai_byte(data, pos) |
    (bai_byte(data, pos + 1) << 8) |
    (bai_byte(data, pos + 2) << 16) |
    (bai_byte(data, pos + 3) << 24) |
    (bai_byte(data, pos + 4) << 32) |
    (bai_byte(data, pos + 5) << 40) |
    (bai_byte(data, pos + 6) << 48) |
    (bai_byte(data, pos + 7) << 56)
}

///|
fn bai_check_magic(data : String) -> Bool {
  if data.length() < 4 {
    return false
  }
  data.unsafe_get(0).to_int() == 0x42 &&
    data.unsafe_get(1).to_int() == 0x41 &&
    data.unsafe_get(2).to_int() == 0x49 &&
    data.unsafe_get(3).to_int() == 0x01
}

// ---------------------------------------------------------------------------
// Parsing and querying
// ---------------------------------------------------------------------------

///|
/// Parse a BAI index from its raw bytes.
pub fn parse_bai(data : String) -> BaiIndex raise BaiIndexError {
  if !bai_check_magic(data) {
    bai_fail("Not a BAI index: bad magic bytes")
  }
  let n_ref = bai_read_u32(data, 4)
  let references : Array[BaiReference] = []
  let mut pos = 8
  for r = 0; r < n_ref; r = r + 1 {
    let n_bin = bai_read_u32(data, pos)
    pos = pos + 4
    let bins : Map[Int, Array[BaiChunk]] = Map([], capacity=n_bin)
    for b = 0; b < n_bin; b = b + 1 {
      let bin = bai_read_u32(data, pos)
      pos = pos + 4
      let n_chunk = bai_read_u32(data, pos)
      pos = pos + 4
      let chunks : Array[BaiChunk] = []
      for c = 0; c < n_chunk; c = c + 1 {
        let beg = bai_read_u64(data, pos)
        let end = bai_read_u64(data, pos + 8)
        pos = pos + 16
        chunks.push(BaiChunk::{ chunk_beg: beg, chunk_end: end })
      }
      bins.set(bin, chunks)
    }
    let n_intv = bai_read_u32(data, pos)
    pos = pos + 4
    let linear : Array[Int64] = []
    for iv = 0; iv < n_intv; iv = iv + 1 {
      linear.push(bai_read_u64(data, pos))
      pos = pos + 8
    }
    references.push(BaiReference::{ bins, linear })
  }
  let n_no_coor : Int64? = if pos + 8 <= data.length() {
    Some(bai_read_u64(data, pos))
  } else {
    None
  }
  BaiIndex::{ n_ref, references, n_no_coor }
}

///|
/// Find the chunks of a BAI index that overlap a half-open `[beg, end)`
/// interval on reference `ref_id`.
///
/// Returns the linear-index minimum offset (the virtual offset to start
/// scanning from) together with every chunk whose bin covers the interval.
pub fn bai_query(
  index : BaiIndex,
  ref_id : Int,
  beg : Int,
  end : Int,
) -> BaiQueryResult raise BaiIndexError {
  if ref_id < 0 || ref_id >= index.references.length() {
    bai_fail("Reference id \{ref_id} is out of range (n_ref = \{index.n_ref})")
  }
  if beg < 0 || end < beg {
    bai_fail("Invalid query interval [\{beg}, \{end})")
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