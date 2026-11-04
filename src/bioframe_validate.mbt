/// Interval validation functions.
///
/// Bioframe equivalent: `is_sorted`, `is_contained`, `is_overlapping`,
/// `is_tight`, `is_disjoint`
///
/// These predicates check whether a set of genomic intervals satisfies
/// specific structural properties, useful for debugging and as preconditions
/// for algorithms that assume sorted/disjoint input.

///|
/// Check if intervals are sorted by (chrom, start, end).
/// Within each chromosome, intervals must have non-decreasing start,
/// and ties broken by non-decreasing end.
pub fn bf_is_sorted(intervals : Array[BedInterval]) -> Bool {
  for i = 1; i < intervals.length(); i = i + 1 {
    let prev = intervals[i - 1]
    let cur = intervals[i]
    if prev.chrom > cur.chrom {
      return false
    }
    if prev.chrom == cur.chrom {
      if prev.start > cur.start {
        return false
      }
      if prev.start == cur.start && prev.end > cur.end {
        return false
      }
    }
  }
  true
}

///|
/// Check if all intervals are contained within given bounds (chromsizes).
/// bounds: array of (chrom, length) pairs.
pub fn bf_is_contained(
  intervals : Array[BedInterval],
  bounds : Array[(String, Int32)],
) -> Bool {
  let bound_map : Map[String, Int32] = Map([])
  for b in bounds {
    bound_map[b.0] = b.1
  }
  for iv in intervals {
    match bound_map.get(iv.chrom) {
      Some(length) => if iv.start < 0 || iv.end > length { return false }
      None => return false
    }
  }
  true
}

///|
/// Check if any pair of intervals on the same chromosome overlaps.
/// Overlap: start1 < end2 && start2 < end1 (at least 1bp shared).
pub fn bf_is_overlapping(intervals : Array[BedInterval]) -> Bool {
  for i = 0; i < intervals.length(); i = i + 1 {
    for j = i + 1; j < intervals.length(); j = j + 1 {
      let a = intervals[i]
      let b = intervals[j]
      if a.chrom == b.chrom && a.start < b.end && b.start < a.end {
        return true
      }
    }
  }
  false
}

///|
/// Check if intervals are "tight": no interval has zero or negative length,
/// and no two adjacent intervals on the same chrom are book-ended
/// (i.e., end_i == start_{i+1} would mean they could be merged).
pub fn bf_is_tight(intervals : Array[BedInterval]) -> Bool {
  for i = 0; i < intervals.length(); i = i + 1 {
    if intervals[i].end <= intervals[i].start {
      return false
    }
  }
  for i = 1; i < intervals.length(); i = i + 1 {
    let prev = intervals[i - 1]
    let cur = intervals[i]
    if prev.chrom == cur.chrom && prev.end == cur.start {
      return false
    }
  }
  true
}

///|
/// Check if all intervals are pairwise disjoint (no overlaps).
pub fn bf_is_disjoint(intervals : Array[BedInterval]) -> Bool {
  for i = 0; i < intervals.length(); i = i + 1 {
    for j = i + 1; j < intervals.length(); j = j + 1 {
      let a = intervals[i]
      let b = intervals[j]
      if a.chrom == b.chrom && a.start < b.end && b.start < a.end {
        return false
      }
    }
  }
  true
}

///|
/// Comprehensive validation: returns a report of all failed checks.
pub struct BfValidationReport {
  is_sorted : Bool
  is_contained : Bool
  is_overlapping : Bool
  is_tight : Bool
  is_disjoint : Bool
}

///|
pub fn BfValidationReport::is_sorted(self : BfValidationReport) -> Bool {
  self.is_sorted
}

///|
pub fn BfValidationReport::is_contained(self : BfValidationReport) -> Bool {
  self.is_contained
}

///|
pub fn BfValidationReport::is_overlapping(self : BfValidationReport) -> Bool {
  self.is_overlapping
}

///|
pub fn BfValidationReport::is_tight(self : BfValidationReport) -> Bool {
  self.is_tight
}

///|
pub fn BfValidationReport::is_disjoint(self : BfValidationReport) -> Bool {
  self.is_disjoint
}

///|
pub fn bf_validate(
  intervals : Array[BedInterval],
  bounds~ : Array[(String, Int32)] = [],
) -> BfValidationReport {
  {
    is_sorted: bf_is_sorted(intervals),
    is_contained: if bounds.length() > 0 {
      bf_is_contained(intervals, bounds)
    } else {
      true
    },
    is_overlapping: bf_is_overlapping(intervals),
    is_tight: bf_is_tight(intervals),
    is_disjoint: bf_is_disjoint(intervals),
  }
}
