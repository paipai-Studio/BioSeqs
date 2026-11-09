/// bioframe_transform.mbt — DataFrame-style tile, bisect, stack, flatten
/// Items 7, 8, 9, 10: bf_tile, bf_bisect, bf_stack, bf_flatten

///|
/// Tile a genome (or a subset of chromsizes) into fixed-size non-overlapping windows.
/// Returns windows of size `tile_size`, with optional step offset.
pub fn bf_tile(
  chromsizes : Array[BfChromSize],
  tile_size~ : Int32 = 1000,
  step~ : Int32 = 0, // 0 means step == tile_size (non-overlapping)
) -> Array[BedInterval] {
  let results : Array[BedInterval] = Array(capacity=0)
  let effective_step = if step <= 0 { tile_size } else { step }
  for ci = 0; ci < chromsizes.length(); ci = ci + 1 {
    let cs = chromsizes[ci]
    let len = cs.length
    let mut pos = 0
    while pos < len {
      let end = if pos + tile_size > len { len } else { pos + tile_size }
      results.push(BedInterval::new(cs.name, pos, end))
      pos = pos + effective_step
    }
  }
  results
}

///|
/// Bisect intervals at specified breakpoint positions.
/// Each interval that spans a breakpoint is split into two at that position.
/// Breakpoints are specified as (chrom, position) pairs.
pub fn bf_bisect(
  intervals : Array[BedInterval],
  breakpoints : Array[(String, Int32)],
) -> Array[BedInterval] {
  let results : Array[BedInterval] = Array(capacity=0)
  for i = 0; i < intervals.length(); i = i + 1 {
    let iv = intervals[i]
    // Find breakpoints within this interval
    let bps : Array[Int32] = Array(capacity=0)
    for j = 0; j < breakpoints.length(); j = j + 1 {
      let (bp_chrom, bp_pos) = breakpoints[j]
      if bp_chrom == iv.chrom && bp_pos > iv.start && bp_pos < iv.end {
        bps.push(bp_pos)
      }
    }
    // Sort breakpoints
    bps.sort_by(fn(a, b) { if a < b { -1 } else if a > b { 1 } else { 0 } })
    // Split interval at each breakpoint
    let mut prev = iv.start
    for j = 0; j < bps.length(); j = j + 1 {
      results.push(BedInterval::new(iv.chrom, prev, bps[j]))
      prev = bps[j]
    }
    results.push(BedInterval::new(iv.chrom, prev, iv.end))
  }
  results
}

///|
/// Stack overlapping intervals into non-overlapping layers.
/// Returns (interval, layer_number) pairs where intervals on the same layer
/// do not overlap each other.
pub struct BfStackResult {
  chrom : String
  start : Int32
  end : Int32
  layer : Int32 // 0-based layer index
  source_idx : Int32 // index into original array
}

///|
pub fn bf_stack(intervals : Array[BedInterval]) -> Array[BfStackResult] {
  // Sort by (chrom, start, end)
  let sorted = bf_stack_sort(intervals)
  let results : Array[BfStackResult] = Array(capacity=0)
  // Track end position of last interval on each layer per chrom
  let layer_ends : Array[(String, Int32)] = Array(capacity=0)
  for i = 0; i < sorted.length(); i = i + 1 {
    let (iv, src_idx) = sorted[i]
    // Find first layer where this interval doesn't overlap
    let mut assigned_layer = -1
    for l = 0; l < layer_ends.length(); l = l + 1 {
      let (le_chrom, le_end) = layer_ends[l]
      if le_chrom == iv.chrom && le_end <= iv.start {
        assigned_layer = l
        layer_ends[l] = (iv.chrom, iv.end)
        break
      }
    }
    if assigned_layer < 0 {
      assigned_layer = layer_ends.length()
      layer_ends.push((iv.chrom, iv.end))
    }
    results.push(BfStackResult::{
      chrom: iv.chrom,
      start: iv.start,
      end: iv.end,
      layer: assigned_layer,
      source_idx: src_idx,
    })
  }
  results
}

///|
/// Flatten overlapping intervals into a disjoint set.
/// Returns non-overlapping intervals that cover exactly the same genomic space
/// as the input. Each output interval is tagged with the count of input
/// intervals that cover it.
pub struct BfFlattenResult {
  chrom : String
  start : Int32
  end : Int32
  depth : Int32 // how many input intervals overlap this segment
}

///|
pub fn bf_flatten(intervals : Array[BedInterval]) -> Array[BfFlattenResult] {
  if intervals.length() == 0 {
    return Array(capacity=0)
  }
  // Collect all boundary events: (position, +1 for start, -1 for end)
  // Group by chrom
  let chroms : Array[String] = Array(capacity=0)
  let events : Array[Array[(Int32, Int32)]] = Array(capacity=0)
  for i = 0; i < intervals.length(); i = i + 1 {
    let iv = intervals[i]
    let mut found = false
    for j = 0; j < chroms.length(); j = j + 1 {
      if chroms[j] == iv.chrom {
        events[j].push((iv.start, 1))
        events[j].push((iv.end, -1))
        found = true
      }
    }
    if !found {
      chroms.push(iv.chrom)
      let arr : Array[(Int32, Int32)] = Array(capacity=0)
      arr.push((iv.start, 1))
      arr.push((iv.end, -1))
      events.push(arr)
    }
  }
  // Process each chrom
  let results : Array[BfFlattenResult] = Array(capacity=0)
  for c = 0; c < chroms.length(); c = c + 1 {
    let evs = events[c]
    // Sort by position, with starts before ends at same position
    evs.sort_by(fn(a, b) {
      if a.0 < b.0 {
        -1
      } else if a.0 > b.0 {
        1
      } else if a.1 > b.1 {
        -1 // +1 (start) before -1 (end)
      } else {
        1
      }
    })
    let mut depth = 0
    let mut prev_pos : Int32 = 0
    let mut has_prev = false
    for e = 0; e < evs.length(); e = e + 1 {
      let (pos, delta) = evs[e]
      if has_prev && depth > 0 && pos > prev_pos {
        results.push(BfFlattenResult::{
          chrom: chroms[c],
          start: prev_pos,
          end: pos,
          depth,
        })
      }
      depth = depth + delta
      prev_pos = pos
      has_prev = true
    }
  }
  results
}

///|
/// Internal: sort intervals and track original indices.
fn bf_stack_sort(intervals : Array[BedInterval]) -> Array[(BedInterval, Int32)] {
  let arr : Array[(BedInterval, Int32)] = Array(capacity=intervals.length())
  for i = 0; i < intervals.length(); i = i + 1 {
    arr.push((intervals[i], i))
  }
  arr.sort_by(fn(a, b) {
    let ai = a.0
    let bi = b.0
    if ai.chrom < bi.chrom {
      -1
    } else if ai.chrom > bi.chrom {
      1
    } else if ai.start < bi.start {
      -1
    } else if ai.start > bi.start {
      1
    } else if ai.end < bi.end {
      -1
    } else if ai.end > bi.end {
      1
    } else {
      0
    }
  })
  arr
}
