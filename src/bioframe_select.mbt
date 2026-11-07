/// Chromosome selection, dropping, and cutting operations.
///
/// Bioframe equivalent: `select`, `drop`, `cut`, `expand`

///|
/// Select intervals on specific chromosomes.
pub fn bf_select(
  intervals : Array[BedInterval],
  chroms : Array[String],
) -> Array[BedInterval] {
  let chrom_set : Set[String] = Set::from_array(chroms)
  let result : Array[BedInterval] = Array(capacity=0)
  for iv in intervals {
    if chrom_set.contains(iv.chrom) {
      result.push(iv)
    }
  }
  result
}

///|
/// Drop intervals on specific chromosomes.
pub fn bf_drop(
  intervals : Array[BedInterval],
  chroms : Array[String],
) -> Array[BedInterval] {
  let chrom_set : Set[String] = Set::from_array(chroms)
  let result : Array[BedInterval] = Array(capacity=0)
  for iv in intervals {
    if !chrom_set.contains(iv.chrom) {
      result.push(iv)
    }
  }
  result
}

///|
/// Cut intervals at specified positions.
/// For each cut position on a chromosome, any interval that straddles
/// the position is split into two.
///
/// cuts: array of (chrom, position) pairs
pub fn bf_cut(
  intervals : Array[BedInterval],
  cuts : Array[(String, Int32)],
) -> Array[BedInterval] {
  // Group cut positions by chromosome
  let cut_map : Map[String, Array[Int32]] = Map([])
  for c in cuts {
    match cut_map.get(c.0) {
      Some(arr) => arr.push(c.1)
      None => cut_map[c.0] = [c.1]
    }
  }
  let result : Array[BedInterval] = Array(capacity=0)
  for iv in intervals {
    match cut_map.get(iv.chrom) {
      Some(positions) => {
        // Sort positions and split interval at each position within (start, end)
        let sorted_pos = positions.copy()
        sorted_pos.sort()
        let mut cur_start = iv.start
        for p in sorted_pos {
          if p > cur_start && p < iv.end {
            result.push(BedInterval::{
              chrom: iv.chrom,
              start: cur_start,
              end: p,
            })
            cur_start = p
          }
        }
        result.push(BedInterval::{
          chrom: iv.chrom,
          start: cur_start,
          end: iv.end,
        })
      }
      None => result.push(iv)
    }
  }
  result
}

///|
/// Expand intervals by a fixed amount on each side.
/// Equivalent to bedtools slop but as a standalone function.
pub fn bf_expand(
  intervals : Array[BedInterval],
  pad : Int32,
  chromsizes~ : Array[(String, Int32)] = [],
) -> Array[BedInterval] {
  let bound_map : Map[String, Int32] = Map([])
  for cs in chromsizes {
    bound_map[cs.0] = cs.1
  }
  let result : Array[BedInterval] = Array(capacity=0)
  for iv in intervals {
    let new_start = if iv.start - pad < 0 { 0 } else { iv.start - pad }
    let new_end = match bound_map.get(iv.chrom) {
      Some(length) => if iv.end + pad > length { length } else { iv.end + pad }
      None => iv.end + pad
    }
    result.push(BedInterval::{
      chrom: iv.chrom,
      start: new_start,
      end: new_end,
    })
  }
  result
}

///|
/// Select intervals by start position range (across all chromosomes).
pub fn bf_select_by_start(
  intervals : Array[BedInterval],
  min_start : Int32,
  max_start : Int32,
) -> Array[BedInterval] {
  let result : Array[BedInterval] = Array(capacity=0)
  for iv in intervals {
    if iv.start >= min_start && iv.start < max_start {
      result.push(iv)
    }
  }
  result
}

///|
/// Select intervals by length range.
pub fn bf_select_by_length(
  intervals : Array[BedInterval],
  min_len : Int32,
  max_len : Int32,
) -> Array[BedInterval] {
  let result : Array[BedInterval] = Array(capacity=0)
  for iv in intervals {
    let len = iv.end - iv.start
    if len >= min_len && len <= max_len {
      result.push(iv)
    }
  }
  result
}
