/// pysam_bam_qual_stats.mbt — Quality score statistics and trimming
/// Inspired by pysam quality analysis

///|
/// Per-position quality statistics.
pub struct PysamQualPosStats {
  pos : Int32
  mean : Double
  median : Double
  min : Int32
  max : Int32
  q25 : Double
  q75 : Double
}

///|
/// Calculate per-position quality statistics across reads.
pub fn pysam_qual_by_position(
  sam_text : String,
  max_pos~ : Int32 = 200,
) -> Array[PysamQualPosStats] {
  let lines = sam_text.split("\n").to_array()
  // Collect qualities per position
  let quals : Array[Array[Int32]] = Array::make(max_pos, Array(capacity=0))
  for i = 0; i < lines.length(); i = i + 1 {
    let line = lines[i].to_owned()
    if line.length() == 0 || line[0].to_int() == 64 {
      continue
    }
    let fields = line.split("\t").to_array()
    if fields.length() < 11 {
      continue
    }
    let qual_str = fields[10].to_owned()
    if qual_str == "*" {
      continue
    }
    for q = 0; q < qual_str.length() && q < max_pos; q = q + 1 {
      quals[q].push(qual_str[q].to_int() - 33)
    }
  }
  // Compute stats
  let result : Array[PysamQualPosStats] = Array(capacity=0)
  for p = 0; p < max_pos; p = p + 1 {
    if quals[p].length() == 0 {
      continue
    }
    let sorted = quals[p].copy()
    sorted.sort_by(fn(a, b) { if a < b { -1 } else if a > b { 1 } else { 0 } })
    let n = sorted.length()
    let mut total = 0
    for j = 0; j < n; j = j + 1 {
      total = total + sorted[j]
    }
    let mean = total.to_double() / n.to_double()
    let median = sorted[n / 2].to_double()
    let min_val = sorted[0]
    let max_val = sorted[n - 1]
    let q25 = sorted[n / 4].to_double()
    let q75 = sorted[3 * n / 4].to_double()
    result.push({ pos: p, mean, median, min: min_val, max: max_val, q25, q75, })
  }
  result
}

///|
/// Quality trimming by sliding window (like Trimmomatic SLIDINGWINDOW).
pub fn pysam_trim_qual_sliding(
  seq : String,
  qual : String,
  window~ : Int32 = 4,
  threshold~ : Int32 = 20,
  phred~ : Int32 = 33,
) -> (String, String) {
  let n = seq.length()
  if n < window || qual.length() < n {
    return (seq, qual)
  }
  // Scan from left, find first position where window avg < threshold
  let mut trim_pos = n
  for i = 0; i <= n - window; i = i + 1 {
    let mut sum = 0
    for j = 0; j < window; j = j + 1 {
      sum = sum + (qual[i + j].to_int() - phred)
    }
    if sum / window < threshold {
      trim_pos = i
      break
    }
  }
  // Also scan from right
  let mut trim_end = 0
  for i = trim_pos - 1; i >= window - 1; i = i - 1 {
    let mut sum = 0
    for j = 0; j < window; j = j + 1 {
      sum = sum + (qual[i - j].to_int() - phred)
    }
    if sum / window >= threshold {
      trim_end = i + 1
      break
    }
  }
  if trim_end <= trim_pos {
    return ("", "")
  }
  (seq[trim_pos:trim_end].to_owned(), qual[trim_pos:trim_end].to_owned())
}

///|
/// Mean quality of a read.
pub fn pysam_read_mean_qual(qual : String, phred~ : Int32 = 33) -> Double {
  if qual.length() == 0 {
    return 0.0
  }
  let mut total = 0
  for i = 0; i < qual.length(); i = i + 1 {
    total = total + (qual[i].to_int() - phred)
  }
  total.to_double() / qual.length().to_double()
}

///|
/// Filter reads by mean quality.
pub fn pysam_filter_by_mean_qual(
  sam_text : String,
  min_mean_qual~ : Int32 = 20,
  phred~ : Int32 = 33,
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
    if fields.length() < 11 {
      let _ = sb.write(line + "\n")
      continue
    }
    let qual = fields[10].to_owned()
    let mean_q = pysam_read_mean_qual(qual, phred~)
    if mean_q >= min_mean_qual.to_double() {
      let _ = sb.write(line + "\n")
    }
  }
  sb.to_string()
}
