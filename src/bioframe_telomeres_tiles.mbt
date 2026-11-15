/// bioframe_telomeres_tiles.mbt — Special genomic intervals B-functions.
///
/// - bf_telomeres:       build telomeric intervals per chromosome
/// - bf_tiled_regions:   build variable-size tiled regions over a genome
///
/// Bioframe 0.8.0  provides these in `bioframe.spec` / `bioframe.arrops`.
/// MoonBit target here: keep types simple (BedInterval-array-in /
/// BedInterval-array-out), compose easily with existing bf_merge /
/// bf_complement / bf_expand.

///|
/// A single tiled region row.
pub struct BfTile {
  chrom : String
  start : Int
  end : Int
  index : Int32   /// which tile along the chromosome (0-based)
  is_last : Bool  /// true if this tile is the last on its chrom
}

// ================================================================
// 1. Telomeres
// ================================================================

///|
/// Build telomeric intervals at each end of every chromosome.
///
/// Telomeres are stretches of repetitive DNA at chromosomal termini.
/// This function doesn't know the *actual* telomere repeat tract length
/// (that requires specialised assembly or long-read data), but produces
/// *placeholder intervals* of a chosen `pad` bp from each chromosome end
/// — exactly like Bioframe's viewframes.
///
/// - `chromsizes`: [(name, length)]  (use bf_assembly_chromsizes /
///   bf_fetch_chromsizes_from_ucsc)
/// - `pad`: telomere length per side.  Defaults follow typical vertebrate
///   estimates but can be overridden for any genome.
///   * 10_000 bp is a common default for hg38/mm39.
///   * 500 bp suffices for yeast / fly.
///
/// Returns Array[BedInterval] — two intervals per linear chromosome
/// (p-telomere at start, q-telomere at end).  Telomere intervals never
/// overlap the chromosome end itself: the interval [0, pad) sits at the
/// start, and [length - pad, length) sits at the end.
pub fn bf_telomeres(
  chromsizes : Array[(String, Int)],
  pad~ : Int = 10000,
) -> Array[BedInterval] {
  let out : Array[BedInterval] = Array(capacity=0)
  for cs in chromsizes {
    let chrom = cs.0
    let length = cs.1
    if length <= 2 * pad {
      // Too short — skip (can't meaningfully define two telomeres)
      continue
    }
    // p-arm telomere
    out.push(BedInterval::new(chrom, 0, pad))
    // q-arm telomere
    out.push(BedInterval::new(chrom, length - pad, length))
  }
  out
}

///|
/// Like bf_telomeres but returns rows with explicit side label ("p"/"q").
pub struct BfTelomereRow {
  chrom : String
  start : Int
  end : Int
  side : Char  /// 'p' or 'q'
}

///|
/// bf_telomeres variant with explicit p/q side labels.
pub fn bf_telomeres_labelled(
  chromsizes : Array[(String, Int)],
  pad~ : Int = 10000,
) -> Array[BfTelomereRow] {
  let out : Array[BfTelomereRow] = Array(capacity=0)
  for cs in chromsizes {
    let chrom = cs.0
    let length = cs.1
    if length <= 2 * pad {
      continue
    }
    out.push(BfTelomereRow::{
      chrom,
      start: 0,
      end: pad,
      side: 'p',
    })
    out.push(BfTelomereRow::{
      chrom,
      start: length - pad,
      end: length,
      side: 'q',
    })
  }
  out
}

// ================================================================
// 2. Tiled regions (variable-size)
// ================================================================

///|
/// Build a tiling of variable-width regions over every chromosome.
///
/// - `chromsizes`: as above
/// - `tile_size`: width of each tile in bp  (default 1_000_000)
/// - `stride`:    step between tile starts.  If `stride == tile_size`
///                tiles are adjacent; if `stride < tile_size` they
///                overlap.  Default `tile_size` (adjacent, no gap).
/// - `start_offset`:  skip this many bp at the start of each chrom
///                    (e.g. to exclude the telomere).  Default 0.
///
/// Returns Array[BfTile].  Last tile on each chromosome is clipped to
/// the chromosome length and marked `is_last = true`.
pub fn bf_tiled_regions(
  chromsizes : Array[(String, Int)],
  tile_size~ : Int = 1000000,
  stride~ : Int = 1000000,
  start_offset~ : Int = 0,
) -> Array[BfTile] {
  let s = if stride <= 0 { tile_size } else { stride }
  let ts = if tile_size <= 0 { 1 } else { tile_size }
  let out : Array[BfTile] = Array(capacity=0)
  for cs in chromsizes {
    let chrom = cs.0
    let length = cs.1
    let mut idx : Int32 = 0
    let mut cur = start_offset
    while cur + ts <= length {
      out.push(BfTile::{
        chrom,
        start: cur,
        end: cur + ts,
        index: idx,
        is_last: false,
      })
      idx = idx + 1
      cur = cur + s
    }
    if cur < length {
      out.push(BfTile::{
        chrom,
        start: cur,
        end: length,
        index: idx,
        is_last: true,
      })
    }
  }
  out
}

///|
/// Convenience: return tiled regions as BedInterval array (no index info).
pub fn bf_tiled_as_bed(
  chromsizes : Array[(String, Int)],
  tile_size~ : Int = 1000000,
  stride~ : Int = 1000000,
) -> Array[BedInterval] {
  let tiles = bf_tiled_regions(chromsizes, tile_size=1000000, stride=1000000)
  let out : Array[BedInterval] = []
  for t in tiles {
    out.push(BedInterval::new(t.chrom, t.start, t.end))
  }
  out
}
