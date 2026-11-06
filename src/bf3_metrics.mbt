// Bioframe third-round: interval metrics and statistics (items 14-18)
// frac_genome_covered, interval_density, bounds, gap_lengths, length_distribution

///|
pub fn bf3_frac_genome_covered(
  intervals : Array[(String, Int32, Int32)],
  chrom_sizes : Array[(String, Int32)],
) -> Double {
  // Calculate fraction of genome covered by intervals
  // First merge overlapping intervals per chromosome, then sum lengths
  let mut total_genome = 0
  for i in 0..<chrom_sizes.length() {
    let (_, sz) = chrom_sizes[i]
    total_genome = total_genome + sz
  }
  if total_genome == 0 {
    return 0.0
  }
  // Group by chrom and merge
  let chroms : Array[String] = []
  for i in 0..<intervals.length() {
    let (chrom, _, _) = intervals[i]
    let mut found = false
    for j in 0..<chroms.length() {
      if chroms[j] == chrom {
        found = true
        break
      }
    }
    if not(found) {
      chroms.push(chrom)
    }
  }
  let mut total_covered = 0
  for ci in 0..<chroms.length() {
    let chrom = chroms[ci]
    let segs : Array[(Int32, Int32)] = []
    for i in 0..<intervals.length() {
      let (c, s, e) = intervals[i]
      if c == chrom {
        segs.push((s, e))
      }
    }
    segs.sort_by(fn(a, b) {
      let (sa, _) = a
      let (sb, _) = b
      if sa < sb {
        -1
      } else if sa > sb {
        1
      } else {
        0
      }
    })
    // Merge
    let mut covered = 0
    let mut cur_start = -1
    let mut cur_end = -1
    for i in 0..<segs.length() {
      let (s, e) = segs[i]
      if cur_start < 0 {
        cur_start = s
        cur_end = e
      } else if s <= cur_end {
        if e > cur_end {
          cur_end = e
        }
      } else {
        covered = covered + (cur_end - cur_start)
        cur_start = s
        cur_end = e
      }
    }
    if cur_start >= 0 {
      covered = covered + (cur_end - cur_start)
    }
    total_covered = total_covered + covered
  }
  total_covered.to_double() / total_genome.to_double()
}

///|
pub fn bf3_interval_density(
  intervals : Array[(Int32, Int32)],
  region_start : Int32,
  region_end : Int32,
) -> Double {
  // Calculate interval density within a region (covered bp / region size)
  let region_size = region_end - region_start
  if region_size <= 0 {
    return 0.0
  }
  let mut total_covered = 0
  for i in 0..<intervals.length() {
    let (start, end_) = intervals[i]
    let lo = if start > region_start { start } else { region_start }
    let hi = if end_ < region_end { end_ } else { region_end }
    if hi > lo {
      total_covered = total_covered + (hi - lo)
    }
  }
  total_covered.to_double() / region_size.to_double()
}

///|
pub fn bf3_bounds(
  intervals : Array[(String, Int32, Int32)],
) -> (String, Int32, Int32)? {
  // Get overall bounds of interval set (min start, max end)
  if intervals.is_empty() {
    return None
  }
  let (first_chrom, first_start, first_end) = intervals[0]
  let mut min_start = first_start
  let mut max_end = first_end
  let mut chrom = first_chrom
  for i in 1..<intervals.length() {
    let (c, s, e) = intervals[i]
    if s < min_start {
      min_start = s
      chrom = c
    }
    if e > max_end {
      max_end = e
    }
  }
  Some((chrom, min_start, max_end))
}

///|
pub fn bf3_gap_lengths(
  intervals : Array[(String, Int32, Int32)],
) -> Array[(String, Int32)] {
  // Calculate gap lengths between adjacent intervals (sorted by chrom, start)
  let sorted = intervals.copy()
  sorted.sort_by(fn(a, b) {
    let (ca, sa, _) = a
    let (cb, sb, _) = b
    if ca < cb {
      -1
    } else if ca > cb {
      1
    } else if sa < sb {
      -1
    } else if sa > sb {
      1
    } else {
      0
    }
  })
  let gaps : Array[(String, Int32)] = []
  for i in 1..<sorted.length() {
    let (prev_chrom, _, prev_end) = sorted[i - 1]
    let (cur_chrom, cur_start, _) = sorted[i]
    if cur_chrom == prev_chrom && cur_start > prev_end {
      gaps.push((cur_chrom, cur_start - prev_end))
    }
  }
  gaps
}

///|
pub struct Bf3LengthStats {
  count : Int32
  total_length : Int32
  mean_length : Double
  min_length : Int32
  max_length : Int32
  median_length : Double
  q25 : Double
  q75 : Double
  std_length : Double
}

///|
pub fn bf3_length_distribution(
  intervals : Array[(String, Int32, Int32)],
) -> Bf3LengthStats {
  // Compute length distribution statistics
  let lengths : Array[Int32] = []
  for i in 0..<intervals.length() {
    let (_, start, end_) = intervals[i]
    lengths.push(end_ - start)
  }
  if lengths.is_empty() {
    return {
      count: 0,
      total_length: 0,
      mean_length: 0.0,
      min_length: 0,
      max_length: 0,
      median_length: 0.0,
      q25: 0.0,
      q75: 0.0,
      std_length: 0.0,
    }
  }
  lengths.sort_by(fn(a, b) { if a < b { -1 } else if a > b { 1 } else { 0 } })
  let n = lengths.length()
  let mut total = 0
  let mut min_val = lengths[0]
  let mut max_val = lengths[0]
  for i in 0..<n {
    total = total + lengths[i]
    if lengths[i] < min_val {
      min_val = lengths[i]
    }
    if lengths[i] > max_val {
      max_val = lengths[i]
    }
  }
  let mean = total.to_double() / n.to_double()
  let mut variance = 0.0
  for i in 0..<n {
    let diff = lengths[i].to_double() - mean
    variance = variance + diff * diff
  }
  variance = variance / n.to_double()
  let std = variance.sqrt()
  // Quantiles
  let q25_idx = (n.to_double() * 0.25).floor().to_int()
  let q50_idx = (n.to_double() * 0.50).floor().to_int()
  let q75_idx = (n.to_double() * 0.75).floor().to_int()
  let q25i = if q25_idx >= n { n - 1 } else { q25_idx }
  let q75i = if q75_idx >= n { n - 1 } else { q75_idx }
  {
    count: n,
    total_length: total,
    mean_length: mean,
    min_length: min_val,
    max_length: max_val,
    median_length: lengths[q50_idx].to_double(),
    q25: lengths[q25i].to_double(),
    q75: lengths[q75i].to_double(),
    std_length: std,
  }
}
