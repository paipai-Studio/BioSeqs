/// pysam_bam_insert_size.mbt — Insert size distribution statistics
/// Inspired by pysam AlignmentFile: template length analysis

///|
/// Insert size distribution statistics.
pub struct PysamInsertSizeStats {
  n_pairs : Int32
  mean : Double
  stdev : Double
  median : Double
  min : Int32
  max : Int32
  q1 : Double // 25th percentile
  q3 : Double // 75th percentile
  histogram : Array[(Int32, Int32)] // (insert_size, count)
}

///|
/// Calculate insert size distribution from SAM.
pub fn pysam_insert_size_stats(
  sam_text : String,
  bin_size~ : Int32 = 1,
) -> PysamInsertSizeStats {
  let lines = sam_text.split("\n").to_array()
  let sizes : Array[Int32] = Array(capacity=0)
  for i = 0; i < lines.length(); i = i + 1 {
    let line = lines[i].to_owned()
    if line.length() == 0 || line[0].to_int() == 64 {
      continue
    }
    let fields = line.split("\t").to_array()
    if fields.length() < 9 {
      continue
    }
    let flag = parse_int(fields[1].to_owned())
    // Only count properly paired, first-in-pair, non-duplicate
    if (flag & 1) == 0 ||
      (flag & 2) == 0 ||
      (flag & 64) == 0 ||
      (flag & 1024) != 0 {
      continue
    }
    let tlen = parse_int(fields[8].to_owned())
    if tlen > 0 {
      sizes.push(tlen)
    }
  }
  let n = sizes.length()
  if n == 0 {
    return {
      n_pairs: 0,
      mean: 0.0,
      stdev: 0.0,
      median: 0.0,
      min: 0,
      max: 0,
      q1: 0.0,
      q3: 0.0,
      histogram: [],
    }
  }
  // Sort for percentiles
  sizes.sort_by(fn(a, b) { if a < b { -1 } else if a > b { 1 } else { 0 } })
  let mut total = 0
  for i = 0; i < n; i = i + 1 {
    total = total + sizes[i]
  }
  let mean = total.to_double() / n.to_double()
  let mut sq_sum = 0.0
  for i = 0; i < n; i = i + 1 {
    let diff = sizes[i].to_double() - mean
    sq_sum = sq_sum + diff * diff
  }
  let stdev = @math.pow(sq_sum / n.to_double(), 0.5)
  let median = pysam_is_percentile(sizes, 0.5)
  let q1 = pysam_is_percentile(sizes, 0.25)
  let q3 = pysam_is_percentile(sizes, 0.75)
  let min_val = sizes[0]
  let max_val = sizes[n - 1]
  // Build histogram
  let hist : Array[(Int32, Int32)] = Array(capacity=0)
  let mut i = 0
  while i < n {
    let bin = sizes[i] / bin_size * bin_size
    let mut count = 0
    while i < n && sizes[i] / bin_size * bin_size == bin {
      count = count + 1
      i = i + 1
    }
    hist.push((bin, count))
  }
  {
    n_pairs: n,
    mean,
    stdev,
    median,
    min: min_val,
    max: max_val,
    q1,
    q3,
    histogram: hist,
  }
}

///|
/// Calculate percentile from sorted array.
fn pysam_is_percentile(sorted : Array[Int32], p : Double) -> Double {
  let n = sorted.length()
  if n == 0 {
    return 0.0
  }
  let idx = p * (n - 1).to_double()
  let lo = idx.floor().to_int()
  let hi = if lo + 1 < n { lo + 1 } else { lo }
  let frac = idx - lo.to_double()
  sorted[lo].to_double() * (1.0 - frac) + sorted[hi].to_double() * frac
}

///|
/// Filter reads by insert size range.
pub fn pysam_filter_insert_size(
  sam_text : String,
  min_size : Int32,
  max_size : Int32,
) -> String {
  let lines = sam_text.split("\n").to_array()
  let sb = StringBuilder()
  for i = 0; i < lines.length(); i = i + 1 {
    let line = lines[i].to_owned()
    if line.length() == 0 {
      continue
    }
    if line[0].to_int() == 64 {
      let _ = sb.write(line + "\n")
      continue
    }
    let fields = line.split("\t").to_array()
    if fields.length() < 9 {
      let _ = sb.write(line + "\n")
      continue
    }
    let tlen = parse_int(fields[8].to_owned())
    let abs_tlen = if tlen < 0 { -tlen } else { tlen }
    if abs_tlen >= min_size && abs_tlen <= max_size {
      let _ = sb.write(line + "\n")
    }
  }
  sb.to_string()
}
