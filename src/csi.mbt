///|
/// CSI (Coordinate Sorted Index) — htslib's index for VCF/BCF (and generic
/// tab-delimited files), the companion to BCF like BAI is to BAM.
///
/// Port of htslib's `hts_idx_load_core` for `HTS_FMT_CSI`. Key differences from
/// BAI/Tabix: each bin carries a `loffset` (virtual offset of the first
/// overlapping record), there is NO linear index, and `n_no_coor` (unplaced
/// records) is always written at the end.

/// Error raised by CSI parsing.
pub suberror CsiError {
  CsiError(String)
}

///|
fn[T] csi_fail(msg : String) -> T raise CsiError {
  raise CsiError(msg)
}

///|
fn csi_byte(data : String, pos : Int) -> Int64 {
  (data.unsafe_get(pos).to_int() & 0xFF).to_int64()
}

///|
fn csi_check_range(data : String, pos : Int, len : Int, what : String) -> Unit raise CsiError {
  if pos < 0 || pos + len > data.length() {
    csi_fail(
      "Cannot read \{what} at offset \{pos}: premature end (length \{data.length()})",
    )
  }
}

///|
fn csi_read_u32(data : String, pos : Int) -> Int raise CsiError {
  csi_check_range(data, pos, 4, "u32")
  (data.unsafe_get(pos).to_int() & 0xFF) |
    ((data.unsafe_get(pos + 1).to_int() & 0xFF) << 8) |
    ((data.unsafe_get(pos + 2).to_int() & 0xFF) << 16) |
    ((data.unsafe_get(pos + 3).to_int() & 0xFF) << 24)
}

///|
fn csi_read_u64(data : String, pos : Int) -> Int64 raise CsiError {
  csi_check_range(data, pos, 8, "u64")
  csi_byte(data, pos) |
    (csi_byte(data, pos + 1) << 8) |
    (csi_byte(data, pos + 2) << 16) |
    (csi_byte(data, pos + 3) << 24) |
    (csi_byte(data, pos + 4) << 32) |
    (csi_byte(data, pos + 5) << 40) |
    (csi_byte(data, pos + 6) << 48) |
    (csi_byte(data, pos + 7) << 56)
}

///|
pub struct CsiChunk {
  chunk_beg : Int64
  chunk_end : Int64
}

///|
pub struct CsiBin {
  bin : Int
  loffset : Int64
  chunks : Array[CsiChunk]
}

///|
pub struct CsiReference {
  bins : Array[CsiBin]
}

///|
pub struct CsiIndex {
  min_shift : Int
  n_lvls : Int
  n_no_coor : Int64
  references : Array[CsiReference]
}

///|
/// Check the "CSI\1" magic.
pub fn is_csi_magic(data : String) -> Bool {
  if data.length() < 4 {
    return false
  }
  data.unsafe_get(0).to_int() == 0x43 && // 'C'
  data.unsafe_get(1).to_int() == 0x53 && // 'S'
  data.unsafe_get(2).to_int() == 0x49 && // 'I'
  data.unsafe_get(3).to_int() == 0x01
}

///|
/// Parse a raw (decompressed) CSI index.
pub fn parse_csi(data : String) -> CsiIndex raise CsiError {
  if not(is_csi_magic(data)) {
    csi_fail("Invalid CSI magic")
  }
  let min_shift = csi_read_u32(data, 4)
  let n_lvls = csi_read_u32(data, 8)
  let l_meta = csi_read_u32(data, 12)
  let mut pos = 16 + l_meta
  let n_ref = csi_read_u32(data, pos)
  pos = pos + 4
  let references : Array[CsiReference] = Array(capacity=0)
  for _r = 0; _r < n_ref; _r = _r + 1 {
    let n_bin = csi_read_u32(data, pos)
    pos = pos + 4
    let bins : Array[CsiBin] = Array(capacity=0)
    for _b = 0; _b < n_bin; _b = _b + 1 {
      let bin = csi_read_u32(data, pos)
      pos = pos + 4
      let loffset = csi_read_u64(data, pos)
      pos = pos + 8
      let n_chunk = csi_read_u32(data, pos)
      pos = pos + 4
      let chunks : Array[CsiChunk] = Array(capacity=0)
      for _c = 0; _c < n_chunk; _c = _c + 1 {
        let beg = csi_read_u64(data, pos)
        pos = pos + 8
        let end = csi_read_u64(data, pos)
        pos = pos + 8
        chunks.push(CsiChunk::{ chunk_beg: beg, chunk_end: end })
      }
      bins.push(CsiBin::{ bin, loffset, chunks })
    }
    references.push(CsiReference::{ bins })
  }
  let n_no_coor = if pos + 8 <= data.length() {
    csi_read_u64(data, pos)
  } else {
    0L
  }
  CsiIndex::{ min_shift, n_lvls, n_no_coor, references }
}

///|
/// The list of bin numbers that may contain a half-open interval `[beg, end)`,
/// generalized by `min_shift` and `n_lvls` (mirrors htslib `hts_reg2bins`).
pub fn csi_reg2bins(beg : Int, end : Int, min_shift : Int, n_lvls : Int) -> Array[Int] {
  let bins : Array[Int] = Array(capacity=0)
  if beg >= end {
    return bins
  }
  let e = end - 1
  bins.push(0)
  let mut l = 1
  while l <= n_lvls {
    let s = min_shift + 3 * (n_lvls - l)
    let offset = ((1 << (3 * l)) - 1) / 7
    let mut k = (beg >> s) + offset
    let upper = (e >> s) + offset
    while k <= upper {
      bins.push(k)
      k = k + 1
    }
    l = l + 1
  }
  bins
}
///|
pub struct CsiQueryResult {
  min_off : Int64
  chunks : Array[CsiChunk]
}

///|
fn csi_bin_first(l : Int) -> Int {
  ((1 << (3 * l)) - 1) / 7
}

///|
fn csi_find_bin(reference : CsiReference, bin : Int) -> CsiBin? {
  for b in reference.bins {
    if b.bin == bin {
      return Some(b)
    }
  }
  None
}

///|
/// Query a half-open interval `[beg, end)` against a CSI reference, returning
/// the chunks (virtual-offset pairs) covering it plus the minimum starting
/// offset (htslib's `loff` walk: the first overlapping record's offset).
pub fn csi_query(
  index : CsiIndex,
  ref_id : Int,
  beg : Int,
  end : Int,
) -> CsiQueryResult raise CsiError {
  if ref_id < 0 || ref_id >= index.references.length() {
    csi_fail("Reference id \{ref_id} out of range (n_ref = \{index.references.length()})")
  }
  if beg < 0 || end < beg {
    csi_fail("Invalid query interval [\{beg}, \{end})")
  }
  let reference = index.references[ref_id]
  // min_off: walk up the deepest bin covering `beg` until an existing bin (htslib)
  let rel_off = beg >> index.min_shift
  let mut bin = csi_bin_first(index.n_lvls) + rel_off
  let mut min_off = 0L
  let mut found = false
  while not(found) && bin != 0 {
    match csi_find_bin(reference, bin) {
      Some(b) => {
        min_off = b.loffset
        found = true
      }
      None => {
        let parent = (bin - 1) >> 3
        let first = (parent << 3) + 1
        if bin > first {
          bin = bin - 1
        } else {
          bin = parent
        }
      }
    }
  }
  if not(found) {
    match csi_find_bin(reference, 0) {
      Some(b) => min_off = b.loffset
      None => ()
    }
  }
  // chunks from covering bins
  let chunks : Array[CsiChunk] = Array(capacity=0)
  let bins = csi_reg2bins(beg, end, index.min_shift, index.n_lvls)
  for b in bins {
    match csi_find_bin(reference, b) {
      Some(cb) =>
        for c in cb.chunks {
          chunks.push(c)
        }
      None => ()
    }
  }
  CsiQueryResult::{ min_off, chunks }
}
