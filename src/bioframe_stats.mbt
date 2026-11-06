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
