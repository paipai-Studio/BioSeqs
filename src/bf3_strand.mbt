// Bioframe third-round: strand-aware and weighted operations (items 33-34)
// strand_aware_overlap, score_weighted_coverage

///|
pub fn bf3_strand_aware_overlap(
  intervals_a : Array[(String, Int32, Int32, String)],
  intervals_b : Array[(String, Int32, Int32, String)],
) -> Array[(Int32, Int32, Int32)] {
  // Find overlapping pairs where strands match (+/+ or -/-)
  // Returns (index_a, index_b, overlap_length)
  let result : Array[(Int32, Int32, Int32)] = []
  for i in 0..<intervals_a.length() {
    let (chrom_a, start_a, end_a, strand_a) = intervals_a[i]
    for j in 0..<intervals_b.length() {
      let (chrom_b, start_b, end_b, strand_b) = intervals_b[j]
      if chrom_a == chrom_b && strand_a == strand_b && strand_a != "." {
        let ov_start = if start_a > start_b { start_a } else { start_b }
        let ov_end = if end_a < end_b { end_a } else { end_b }
        if ov_end > ov_start {
          result.push((i, j, ov_end - ov_start))
        }
      }
    }
  }
  result
}

///|
pub fn bf3_strand_aware_overlap_count(
  intervals_a : Array[(String, Int32, Int32, String)],
  intervals_b : Array[(String, Int32, Int32, String)],
) -> Array[(Int32, Int32)] {
  // For each interval in A, count strand-matching overlaps in B
  // Returns (index_a, count)
  let result : Array[(Int32, Int32)] = []
  for i in 0..<intervals_a.length() {
    let (chrom_a, start_a, end_a, strand_a) = intervals_a[i]
    let mut count = 0
    for j in 0..<intervals_b.length() {
      let (chrom_b, start_b, end_b, strand_b) = intervals_b[j]
      if chrom_a == chrom_b && strand_a == strand_b && strand_a != "." {
        if start_b < end_a && end_b > start_a {
          count = count + 1
        }
      }
    }
    result.push((i, count))
  }
  result
}

///|
pub fn bf3_score_weighted_coverage(
  intervals : Array[(String, Int32, Int32, Double)],
  _chrom_size : Int32,
) -> Array[(Int32, Double)] {
  // Calculate score-weighted coverage at each position
  // Each interval contributes its score to the coverage signal
  let events : Array[(Int32, Double)] = [] // (pos, score_delta)
  for i in 0..<intervals.length() {
    let (_, start, end_, score) = intervals[i]
    events.push((start, score))
    events.push((end_, 0.0 - score))
  }
  events.sort_by(fn(a, b) {
    let (pa, _) = a
    let (pb, _) = b
    if pa < pb {
      -1
    } else if pa > pb {
      1
    } else {
      0
    }
  })
  // Build signal as step function
  let signal : Array[(Int32, Double)] = []
  let mut current_signal = 0.0
  for i in 0..<events.length() {
    let (pos, delta) = events[i]
    current_signal = current_signal + delta
    signal.push((pos, current_signal))
  }
  signal
}

///|
pub fn bf3_score_weighted_coverage_regions(
  intervals : Array[(String, Int32, Int32, Double)],
  chrom : String,
  _chrom_size : Int32,
) -> Array[(String, Int32, Int32, Double)] {
  // Score-weighted coverage as bedGraph-like regions
  let events : Array[(Int32, Double)] = []
  for i in 0..<intervals.length() {
    let (_, start, end_, score) = intervals[i]
    events.push((start, score))
    events.push((end_, 0.0 - score))
  }
  events.sort_by(fn(a, b) {
    let (pa, _) = a
    let (pb, _) = b
    if pa < pb {
      -1
    } else if pa > pb {
      1
    } else {
      0
    }
  })
  let result : Array[(String, Int32, Int32, Double)] = []
  let mut current_signal = 0.0
  let mut prev_pos = 0
  for i in 0..<events.length() {
    let (pos, delta) = events[i]
    if pos > prev_pos && current_signal != 0.0 {
      result.push((chrom, prev_pos, pos, current_signal))
    }
    current_signal = current_signal + delta
    prev_pos = pos
  }
  result
}
