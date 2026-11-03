/// Overlap mode classification.
///
/// Bioframe equivalent: different overlap semantics used in `overlap()`
/// with `how` parameter.
///
/// Four overlap modes:
/// - any: at least 1bp shared (default, most common)
/// - contained: A is fully contained within B (B.start <= A.start && A.end <= B.end)
/// - within: B is fully contained within A (A.start <= B.start && B.end <= A.end)
/// - strict: overlap but neither contains the other

///|
/// Check if two intervals share at least 1bp (any overlap)
pub fn bf_overlap_any(a : BedInterval, b : BedInterval) -> Bool {
  a.chrom == b.chrom && a.start < b.end && b.start < a.end
}

///|
/// Check if interval a is fully contained within interval b
pub fn bf_overlap_contained(a : BedInterval, b : BedInterval) -> Bool {
  a.chrom == b.chrom && b.start <= a.start && a.end <= b.end
}

///|
/// Check if interval b is fully contained within interval a
pub fn bf_overlap_within(a : BedInterval, b : BedInterval) -> Bool {
  a.chrom == b.chrom && a.start <= b.start && b.end <= a.end
}

///|
/// Check if intervals overlap but neither contains the other
pub fn bf_overlap_strict(a : BedInterval, b : BedInterval) -> Bool {
  bf_overlap_any(a, b) &&
  !bf_overlap_contained(a, b) &&
  !bf_overlap_within(a, b)
}

///|
/// Find all pairs (a, b) where a and b overlap with the given mode.
/// Returns array of (index_a, index_b) pairs.
pub fn bf_overlap_pairs(
  a : Array[BedInterval],
  b : Array[BedInterval],
  mode~ : String = "any",
) -> Array[(Int32, Int32)] {
  let result : Array[(Int32, Int32)] = Array(capacity=0)
  for i = 0; i < a.length(); i = i + 1 {
    for j = 0; j < b.length(); j = j + 1 {
      let overlaps = match mode {
        "any" => bf_overlap_any(a[i], b[j])
        "contained" => bf_overlap_contained(a[i], b[j])
        "within" => bf_overlap_within(a[i], b[j])
        "strict" => bf_overlap_strict(a[i], b[j])
        _ => bf_overlap_any(a[i], b[j])
      }
      if overlaps {
        result.push((i, j))
      }
    }
  }
  result
}

///|
/// Filter intervals in `a` that overlap with at least one in `b` using the given mode.
pub fn bf_overlap_filter(
  a : Array[BedInterval],
  b : Array[BedInterval],
  mode~ : String = "any",
) -> Array[BedInterval] {
  let result : Array[BedInterval] = Array(capacity=0)
  for ai in a {
    let mut found = false
    for bi in b {
      let overlaps = match mode {
        "any" => bf_overlap_any(ai, bi)
        "contained" => bf_overlap_contained(ai, bi)
        "within" => bf_overlap_within(ai, bi)
        "strict" => bf_overlap_strict(ai, bi)
        _ => bf_overlap_any(ai, bi)
      }
      if overlaps {
        found = true
        break
      }
    }
    if found {
      result.push(ai)
    }
  }
  result
}

///|
/// Count overlaps for each interval in `a` with the given mode.
pub fn bf_count_overlaps(
  a : Array[BedInterval],
  b : Array[BedInterval],
  mode~ : String = "any",
) -> Array[Int32] {
  let result : Array[Int32] = Array(capacity=0)
  for ai in a {
    let mut count : Int32 = 0
    for bi in b {
      let overlaps = match mode {
        "any" => bf_overlap_any(ai, bi)
        "contained" => bf_overlap_contained(ai, bi)
        "within" => bf_overlap_within(ai, bi)
        "strict" => bf_overlap_strict(ai, bi)
        _ => bf_overlap_any(ai, bi)
      }
      if overlaps {
        count = count + 1
      }
    }
    result.push(count)
  }
  result
}
