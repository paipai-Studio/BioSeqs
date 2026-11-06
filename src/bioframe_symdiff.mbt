/// Symmetric difference of interval sets.
///
/// Bioframe equivalent: `symmetric_diff` (A ∪ B) - (A ∩ B)
///
/// Returns intervals that are in exactly one of the two sets,
/// not in both. This is the XOR of interval sets.

///|
/// Compute symmetric difference of two interval sets.
/// Result: (A - B) ∪ (B - A), where subtraction removes overlapping portions.
pub fn bf_symmetric_diff(
  a : Array[BedInterval],
  b : Array[BedInterval],
) -> Array[BedInterval] {
  // A - B: portions of A not covered by any B interval
  let a_minus_b = bf_subtract_intervals(a, b)
  // B - A: portions of B not covered by any A interval
  let b_minus_a = bf_subtract_intervals(b, a)
  // Merge the two results
  bedtools_merge(a_minus_b + b_minus_a)
}

///|
/// Subtract intervals in b from a, returning the portions of a
/// that don't overlap with any b interval.
/// This is more precise than bedtools_subtract when we need
/// the actual remaining fragments.
fn bf_subtract_intervals(
  a : Array[BedInterval],
  b : Array[BedInterval],
) -> Array[BedInterval] {
  let result : Array[BedInterval] = Array(capacity=0)
  for ai in a {
    // Find all b intervals on the same chrom that overlap with ai
    let overlapping : Array[BedInterval] = Array(capacity=0)
    for bi in b {
      if bi.chrom == ai.chrom && bi.start < ai.end && ai.start < bi.end {
        overlapping.push(bi)
      }
    }
    if overlapping.length() == 0 {
      result.push(ai)
    } else {
      // Sort overlapping by start
      overlapping.sort_by(fn(x, y) {
        if x.start < y.start {
          -1
        } else if x.start > y.start {
          1
        } else {
          0
        }
      })
      // Walk through and emit gaps
      let mut cur = ai.start
      for oi in overlapping {
        if oi.start > cur {
          result.push(BedInterval::{
            chrom: ai.chrom,
            start: cur,
            end: oi.start,
          })
        }
        if oi.end > cur {
          cur = oi.end
        }
      }
      if cur < ai.end {
        result.push(BedInterval::{ chrom: ai.chrom, start: cur, end: ai.end, })
      }
    }
  }
  result
}

///|
/// Check if two interval sets are equal (same merged regions).
pub fn bf_sets_equal(a : Array[BedInterval], b : Array[BedInterval]) -> Bool {
  let merged_a = bedtools_merge(a)
  let merged_b = bedtools_merge(b)
  if merged_a.length() != merged_b.length() {
    return false
  }
  // Sort both for comparison
  let sorted_a = merged_a.copy()
  sorted_a.sort_by(fn(x, y) {
    if x.chrom < y.chrom {
      -1
    } else if x.chrom > y.chrom {
      1
    } else if x.start < y.start {
      -1
    } else if x.start > y.start {
      1
    } else if x.end < y.end {
      -1
    } else if x.end > y.end {
      1
    } else {
      0
    }
  })
  let sorted_b = merged_b.copy()
  sorted_b.sort_by(fn(x, y) {
    if x.chrom < y.chrom {
      -1
    } else if x.chrom > y.chrom {
      1
    } else if x.start < y.start {
      -1
    } else if x.start > y.start {
      1
    } else if x.end < y.end {
      -1
    } else if x.end > y.end {
      1
    } else {
      0
    }
  })
  for i = 0; i < sorted_a.length(); i = i + 1 {
    if sorted_a[i].chrom != sorted_b[i].chrom ||
      sorted_a[i].start != sorted_b[i].start ||
      sorted_a[i].end != sorted_b[i].end {
      return false
    }
  }
  true
}

///|
/// Compute Jaccard similarity between two interval sets.
/// Jaccard = |A ∩ B| / |A ∪ B|
pub fn bf_jaccard(a : Array[BedInterval], b : Array[BedInterval]) -> Double {
  let merged_a = bedtools_merge(a)
  let merged_b = bedtools_merge(b)
  let union = bedtools_merge(merged_a + merged_b)
  let intersect = bf_intersect_sets(merged_a, merged_b)
  let union_size = bf_total_length(union)
  let intersect_size = bf_total_length(intersect)
  if union_size == 0 {
    return 0.0
  }
  intersect_size.to_double() / union_size.to_double()
}

///|
fn bf_intersect_sets(
  a : Array[BedInterval],
  b : Array[BedInterval],
) -> Array[BedInterval] {
  let result : Array[BedInterval] = Array(capacity=0)
  for ai in a {
    for bi in b {
      if ai.chrom == bi.chrom && ai.start < bi.end && bi.start < ai.end {
        let start = if ai.start > bi.start { ai.start } else { bi.start }
        let end = if ai.end < bi.end { ai.end } else { bi.end }
        result.push(BedInterval::{ chrom: ai.chrom, start, end, })
      }
    }
  }
  bedtools_merge(result)
}

///|
fn bf_total_length(intervals : Array[BedInterval]) -> Int32 {
  let mut total : Int32 = 0
  for iv in intervals {
    total = total + (iv.end - iv.start)
  }
  total
}
