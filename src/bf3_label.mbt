// Bioframe third-round: labeling and filtering (items 19-22)
// label_overlap, select_by_overlap_count, mask/unmask, rank_by_score

///|
pub struct Bf3LabeledInterval {
  chrom : String
  start : Int32
  end_ : Int32
  name : String
  score : Double
  is_overlapping : Bool
}

///|
pub fn bf3_label_overlap(
  intervals : Array[(String, Int32, Int32, String, Double)],
  others : Array[(String, Int32, Int32)],
) -> Array[Bf3LabeledInterval] {
  // Label each interval with whether it overlaps any interval in 'others'
  let result : Array[Bf3LabeledInterval] = []
  for i in 0..<intervals.length() {
    let (chrom, start, end_, name, score) = intervals[i]
    let mut has_overlap = false
    for j in 0..<others.length() {
      let (o_chrom, o_start, o_end) = others[j]
      if o_chrom == chrom && o_start < end_ && o_end > start {
        has_overlap = true
        break
      }
    }
    result.push({
      chrom,
      start,
      end_,
      name,
      score,
      is_overlapping: has_overlap,
    })
  }
  result
}

///|
pub fn bf3_select_by_overlap_count(
  intervals : Array[(String, Int32, Int32)],
  others : Array[(String, Int32, Int32)],
  min_count : Int32,
  max_count : Int32,
) -> Array[(String, Int32, Int32, Int32)] {
  // Select intervals that overlap with [min_count, max_count] intervals from 'others'
  let result : Array[(String, Int32, Int32, Int32)] = []
  for i in 0..<intervals.length() {
    let (chrom, start, end_) = intervals[i]
    let mut count = 0
    for j in 0..<others.length() {
      let (o_chrom, o_start, o_end) = others[j]
      if o_chrom == chrom && o_start < end_ && o_end > start {
        count = count + 1
      }
    }
    if count >= min_count && count <= max_count {
      result.push((chrom, start, end_, count))
    }
  }
  result
}

///|
pub struct Bf3MaskedInterval {
  chrom : String
  start : Int32
  end_ : Int32
  is_masked : Bool
  mask_reason : String
}

///|
pub fn bf3_mask(
  intervals : Array[(String, Int32, Int32)],
  mask_regions : Array[(String, Int32, Int32)],
  reason : String,
) -> Array[Bf3MaskedInterval] {
  // Mask intervals that overlap with mask regions (e.g. blacklist)
  let result : Array[Bf3MaskedInterval] = []
  for i in 0..<intervals.length() {
    let (chrom, start, end_) = intervals[i]
    let mut masked = false
    let mut mask_reason_ = ""
    for j in 0..<mask_regions.length() {
      let (m_chrom, m_start, m_end) = mask_regions[j]
      if m_chrom == chrom && m_start < end_ && m_end > start {
        masked = true
        mask_reason_ = reason
        break
      }
    }
    result.push({
      chrom,
      start,
      end_,
      is_masked: masked,
      mask_reason: mask_reason_,
    })
  }
  result
}

///|
pub fn bf3_unmask(
  masked_intervals : Array[Bf3MaskedInterval],
) -> Array[(String, Int32, Int32)] {
  // Return only unmasked intervals
  let result : Array[(String, Int32, Int32)] = []
  for i in 0..<masked_intervals.length() {
    let m = masked_intervals[i]
    if not(m.is_masked) {
      result.push((m.chrom, m.start, m.end_))
    }
  }
  result
}

///|
pub fn bf3_filter_masked(
  masked_intervals : Array[Bf3MaskedInterval],
) -> Array[(String, Int32, Int32)] {
  // Alias for unmask — return only unmasked
  bf3_unmask(masked_intervals)
}

///|
pub struct Bf3RankedInterval {
  chrom : String
  start : Int32
  end_ : Int32
  score : Double
  rank : Int32 // 1-based rank (highest score = rank 1)
  percentile : Double // 0.0-100.0
}

///|
pub fn bf3_rank_by_score(
  intervals : Array[(String, Int32, Int32, Double)],
  ascending : Bool,
) -> Array[Bf3RankedInterval] {
  // Rank intervals by score (descending by default, ascending if specified)
  let n = intervals.length()
  let indexed : Array[(Int32, String, Int32, Int32, Double)] = []
  for i in 0..<n {
    let (chrom, start, end_, score) = intervals[i]
    indexed.push((i, chrom, start, end_, score))
  }
  indexed.sort_by(fn(a, b) {
    let (_, _, _, _, sa) = a
    let (_, _, _, _, sb) = b
    if ascending {
      if sa < sb {
        -1
      } else if sa > sb {
        1
      } else {
        0
      }
    } else if sa > sb {
      -1
    } else if sa < sb {
      1
    } else {
      0
    }
  })
  let result : Array[Bf3RankedInterval] = []
  for r in 0..<n {
    let (_, chrom, start, end_, score) = indexed[r]
    let percentile = (n - r - 1).to_double() / n.to_double() * 100.0
    result.push({ chrom, start, end_, score, rank: r + 1, percentile, })
  }
  result
}
