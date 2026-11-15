/// bioframe_stats.mbt — DataFrame-style coverage signal and length statistics
/// Items 11, 14: bf_coverage_signal, bf_length_stats

///|
/// Compute per-base coverage signal track for a set of intervals.
/// Returns a signal array per chromosome: (chrom, [coverage at each base]).
/// Note: for large genomes this is memory-intensive; use for small regions.
pub struct BfSignalTrack {
  chrom : String
  start : Int32 // start coordinate of the signal array
  signal : Array[Int32] // coverage at each base position
}

///|
/// Compute per-base coverage signal for intervals on a specific chromosome.
/// Only computes signal for the range [region_start, region_end).
/// For whole-chromosome, pass region_start=0, region_end=chrom_length.
pub fn bf_coverage_signal(
  intervals : Array[BedInterval],
  chrom~ : String,
  region_start~ : Int32 = 0,
  region_end~ : Int32 = 0, // 0 means auto-detect from intervals
) -> BfSignalTrack {
  // Determine region end
  let mut max_end : Int32 = region_end
  if region_end == 0 {
    max_end = 0
    for i = 0; i < intervals.length(); i = i + 1 {
      if intervals[i].chrom == chrom && intervals[i].end > max_end {
        max_end = intervals[i].end
      }
    }
  }
  let track_len = max_end - region_start
  let signal = Array(capacity=track_len)
  for i = 0; i < track_len; i = i + 1 {
    signal.push(0)
  }
  // Add coverage from each interval
  for i = 0; i < intervals.length(); i = i + 1 {
    let iv = intervals[i]
    if iv.chrom != chrom {
      continue
    }
    let s = if iv.start > region_start { iv.start } else { region_start }
    let e = if iv.end > max_end { max_end } else { iv.end }
    for pos = s; pos < e; pos = pos + 1 {
      let idx = pos - region_start
      if idx >= 0 && idx < track_len {
        signal[idx] = signal[idx] + 1
      }
    }
  }
  BfSignalTrack::{ chrom, start: region_start, signal, }
}

///|
/// Statistics over interval lengths.
pub struct BfLengthStats {
  n : Int32 // number of intervals
  min_len : Int32
  max_len : Int32
  mean_len : Double
  median_len : Double
  total_bp : Int64 // total base pairs covered
  n50 : Int32 // N50 statistic (contig-style)
  q1 : Double // 25th percentile
  q3 : Double // 75th percentile
  stdev : Double // standard deviation
}

///|
/// Compute statistics over interval lengths.
/// Returns min, max, mean, median, N50, quartiles, and stdev.
pub fn bf_length_stats(intervals : Array[BedInterval]) -> BfLengthStats {
  let n = intervals.length()
  if n == 0 {
    return BfLengthStats::{
      n: 0,
      min_len: 0,
      max_len: 0,
      mean_len: 0.0,
      median_len: 0.0,
      total_bp: 0,
      n50: 0,
      q1: 0.0,
      q3: 0.0,
      stdev: 0.0,
    }
  }
  // Compute lengths
  let lengths : Array[Int32] = Array(capacity=n)
  let mut total_bp : Int64 = 0
  for i = 0; i < n; i = i + 1 {
    let len = intervals[i].end - intervals[i].start
    lengths.push(len)
    total_bp = total_bp + len.to_int64()
  }
  // Sort lengths
  lengths.sort_by(fn(a, b) { if a < b { -1 } else if a > b { 1 } else { 0 } })
  let min_len = lengths[0]
  let max_len = lengths[n - 1]
  // Mean
  let sum = bf_sum_lengths(lengths)
  let mean_len = sum.to_double() / n.to_double()
  // Median
  let median_len = if n % 2 == 1 {
    lengths[n / 2].to_double()
  } else {
    (lengths[n / 2 - 1] + lengths[n / 2]).to_double() / 2.0
  }
  // Quartiles
  let q1 = bf_percentile(lengths, 25)
  let q3 = bf_percentile(lengths, 75)
  // Stdev
  let mut var_sum = 0.0
  for i = 0; i < n; i = i + 1 {
    let diff = lengths[i].to_double() - mean_len
    var_sum = var_sum + diff * diff
  }
  let stdev = @math.pow(var_sum / n.to_double(), 0.5)
  // N50: sort descending, accumulate until >= 50% of total
  let desc : Array[Int32] = Array(capacity=n)
  for i = 0; i < n; i = i + 1 {
    desc.push(lengths[n - 1 - i])
  }
  let half = sum / 2
  let mut cum = 0
  let mut n50 = 0
  for i = 0; i < n; i = i + 1 {
    cum = cum + desc[i]
    if cum >= half {
      n50 = desc[i]
      break
    }
  }
  BfLengthStats::{
    n,
    min_len,
    max_len,
    mean_len,
    median_len,
    total_bp,
    n50,
    q1,
    q3,
    stdev,
  }
}

///|
fn bf_sum_lengths(lengths : Array[Int32]) -> Int32 {
  let mut sum = 0
  for i = 0; i < lengths.length(); i = i + 1 {
    sum = sum + lengths[i]
  }
  sum
}

///|
/// Linear interpolation percentile from sorted array.
fn bf_percentile(sorted : Array[Int32], pct : Int32) -> Double {
  let n = sorted.length()
  if n == 0 {
    return 0.0
  }
  if n == 1 {
    return sorted[0].to_double()
  }
  let rank = pct.to_double() / 100.0 * (n - 1).to_double()
  let lo = rank.floor().to_int()
  let hi = if lo + 1 < n { lo + 1 } else { lo }
  let frac = rank - lo.to_double()
  sorted[lo].to_double() + frac * (sorted[hi] - sorted[lo]).to_double()
}

// ================================================================
// New: per-interval coverage (Bioframe-style)
// ================================================================

///|
/// Coverage of each interval in A by intervals in B.
///
/// For every interval in `a`, compute how many base pairs of it are
/// covered by any interval in `b`.  Overlaps in B are merged *first*
/// so that overlapping B regions don't double-count.
///
/// Returns an Int32 array with exactly the same length as `a`.
///
///     bf_coverage_per_interval(a, b)
///
/// Example (no merge step shown for brevity):
///     A:  [0, 100), [200, 300)
///     B:  [50, 150)
///     → coverage: [50, 0]
///
pub fn bf_coverage_per_interval(
  a : Array[BedInterval],
  b : Array[BedInterval],
) -> Array[Int32] {
  // Merge b so we count each covered bp only once
  let merged_b = bf_coverage_merge(b)
  let out : Array[Int32] = Array(capacity=a.length())
  for ai in a {
    let mut covered : Int32 = 0
    // Only scan B on the same chromosome that overlaps this A interval
    for bi in merged_b {
      if bi.chrom != ai.chrom {
        continue
      }
      if bi.end <= ai.start || bi.start >= ai.end {
        continue
      }
      let overlap_start = if bi.start > ai.start { bi.start } else { ai.start }
      let overlap_end = if bi.end < ai.end { bi.end } else { ai.end }
      covered = covered + (overlap_end - overlap_start)
    }
    out.push(covered)
  }
  out
}

///|
/// Same as bf_coverage_per_interval but also returns the fraction of
/// each A interval that is covered (Double in [0, 1]).
pub fn bf_coverage_per_interval_fraction(
  a : Array[BedInterval],
  b : Array[BedInterval],
) -> Array[Double] {
  let cov = bf_coverage_per_interval(a, b)
  let out : Array[Double] = Array(capacity=a.length())
  let mut i = 0
  while i < a.length() {
    let len = a[i].end - a[i].start
    if len <= 0 {
      out.push(0.0)
    } else {
      out.push(cov[i].to_double() / len.to_double())
    }
    i = i + 1
  }
  out
}

///|
/// Merge helper for coverage — simplified, no min_dist, same-chrom only.
fn bf_coverage_merge(intervals : Array[BedInterval]) -> Array[BedInterval] {
  if intervals.length() == 0 {
    return Array(capacity=0)
  }
  let sorted = intervals.copy()
  sorted.sort_by(fn(x, y) {
    if x.chrom < y.chrom { -1 }
    else if x.chrom > y.chrom { 1 }
    else if x.start < y.start { -1 }
    else if x.start > y.start { 1 }
    else if x.end < y.end { -1 }
    else if x.end > y.end { 1 }
    else { 0 }
  })
  let out : Array[BedInterval] = Array(capacity=0)
  let mut cur = sorted[0]
  for i = 1; i < sorted.length(); i = i + 1 {
    let iv = sorted[i]
    if iv.chrom == cur.chrom && iv.start <= cur.end {
      if iv.end > cur.end {
        cur = BedInterval::new(cur.chrom, cur.start, iv.end)
      }
    } else {
      out.push(cur)
      cur = iv
    }
  }
  out.push(cur)
  out
}
