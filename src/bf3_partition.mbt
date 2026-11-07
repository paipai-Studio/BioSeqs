// Bioframe third-round: partition and sliding (items 10-13)
// genomic_partition, slide_windows, cut_at_positions, split_at_boundaries

///|
pub fn bf3_genomic_partition(
  chrom : String,
  chrom_size : Int32,
  n_parts : Int32,
) -> Array[(String, Int32, Int32, Int32)] {
  // Partition a chromosome into n_parts roughly equal non-overlapping regions
  // Returns (chrom, start, end, part_index)
  let result : Array[(String, Int32, Int32, Int32)] = []
  let part_size = chrom_size / n_parts
  let remainder = chrom_size % n_parts
  let mut current = 0
  for i in 0..<n_parts {
    let size = part_size + (if i < remainder { 1 } else { 0 })
    result.push((chrom, current, current + size, i))
    current = current + size
  }
  result
}

///|
pub fn bf3_genomic_partition_multi(
  chroms : Array[(String, Int32)],
  n_parts : Int32,
) -> Array[(String, Int32, Int32, Int32)] {
  // Partition multiple chromosomes into n_parts total regions
  let mut total_size = 0
  for i in 0..<chroms.length() {
    let (_, sz) = chroms[i]
    total_size = total_size + sz
  }
  let target_size = total_size / n_parts
  let result : Array[(String, Int32, Int32, Int32)] = []
  let mut part_idx = 0
  let mut accumulated = 0
  let mut part_start = 0
  let mut current_chrom = ""
  for ci in 0..<chroms.length() {
    let (chrom, chrom_size) = chroms[ci]
    current_chrom = chrom
    part_start = 0
    accumulated = 0
    let mut pos = 0
    while pos < chrom_size {
      let remaining = chrom_size - pos
      let needed = target_size - accumulated
      if remaining >= needed && needed > 0 {
        pos = pos + needed
        result.push((current_chrom, part_start, pos, part_idx))
        part_idx = part_idx + 1
        accumulated = 0
        part_start = pos
      } else if remaining < needed {
        pos = chrom_size
        accumulated = accumulated + remaining
        if accumulated >= target_size {
          result.push((current_chrom, part_start, pos, part_idx))
          part_idx = part_idx + 1
          accumulated = 0
          part_start = pos
        }
      } else {
        pos = chrom_size
      }
    }
    if part_start < chrom_size && accumulated > 0 {
      result.push((current_chrom, part_start, chrom_size, part_idx))
      part_idx = part_idx + 1
    }
  }
  result
}

///|
pub struct Bf3Window {
  chrom : String
  start : Int32
  end_ : Int32
  window_index : Int32
}

///|
pub fn bf3_slide_windows(
  chrom : String,
  chrom_size : Int32,
  window_size : Int32,
  step : Int32,
) -> Array[Bf3Window] {
  // Generate sliding windows across a chromosome
  let windows : Array[Bf3Window] = []
  let mut idx = 0
  let mut start = 0
  while start + window_size <= chrom_size {
    windows.push({ chrom, start, end_: start + window_size, window_index: idx, })
    idx = idx + 1
    start = start + step
  }
  // Handle last partial window
  if start < chrom_size && start + window_size > chrom_size {
    windows.push({ chrom, start, end_: chrom_size, window_index: idx, })
  }
  windows
}

///|
pub fn bf3_cut_at_positions(
  intervals : Array[(String, Int32, Int32)],
  cut_positions : Array[(String, Int32)],
) -> Array[(String, Int32, Int32)] {
  // Cut intervals at specified positions
  // Each interval that spans a cut position is split into two
  let result : Array[(String, Int32, Int32)] = []
  for i in 0..<intervals.length() {
    let (chrom, start, end_) = intervals[i]
    // Find all cut positions within this interval
    let cuts : Array[Int32] = []
    for j in 0..<cut_positions.length() {
      let (c_chrom, c_pos) = cut_positions[j]
      if c_chrom == chrom && c_pos > start && c_pos < end_ {
        cuts.push(c_pos)
      }
    }
    cuts.sort_by(fn(a, b) { if a < b { -1 } else if a > b { 1 } else { 0 } })
    let mut current_start = start
    for k in 0..<cuts.length() {
      result.push((chrom, current_start, cuts[k]))
      current_start = cuts[k]
    }
    result.push((chrom, current_start, end_))
  }
  result
}

///|
pub fn bf3_split_at_boundaries(
  intervals : Array[(String, Int32, Int32)],
  boundaries : Array[(String, Int32, Int32)],
) -> Array[(String, Int32, Int32)] {
  // Split intervals at the boundaries of another interval set
  // Collect all boundary positions from the boundary set
  let positions : Array[(String, Int32)] = []
  for i in 0..<boundaries.length() {
    let (chrom, start, end_) = boundaries[i]
    positions.push((chrom, start))
    positions.push((chrom, end_))
  }
  bf3_cut_at_positions(intervals, positions)
}
