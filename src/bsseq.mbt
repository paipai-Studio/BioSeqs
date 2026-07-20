///|
/// bsseq module - Bisulfite sequencing data analysis.
/// Inspired by Bioconductor bsseq package.
///
/// Provides core utilities for BS-seq data:
/// - Parse BED and Bismark methylation output
/// - Compute coverage and methylation rates
/// - Filter records by coverage / context (CpG, CHG, CHH)
/// - BSmooth local smoothing of methylation
/// - Differentially methylated region (DMR) detection via t-statistics
/// - Combine multiple samples into coverage matrices

///|
/// Natural logarithm of 2.
let bsseq_ln2 : Double = 0.69314718055994530942

///|
/// Newton's method for square root. Returns 0 for negative input.
pub fn bsseq_sqrt(x : Double) -> Double {
  if x <= 0.0 {
    return 0.0
  }
  let mut guess = if x > 1.0 { x } else { 1.0 }
  let mut i = 0
  while i < 60 {
    let next = (guess + x / guess) / 2.0
    let diff = (next - guess).abs()
    guess = next
    if diff < 0.0000000000001 {
      return guess
    }
    i = i + 1
  }
  guess
}

///|
/// Natural logarithm. Uses sqrt-based range reduction plus Taylor series.
pub fn bsseq_ln(x : Double) -> Double {
  if x <= 0.0 {
    return 0.0
  }
  if x == 1.0 {
    return 0.0
  }
  let mut v = x
  let mut k = 0
  while (v > 2.0 || v < 0.5) && k < 40 {
    v = bsseq_sqrt(v)
    k = k + 1
  }
  let y = (v - 1.0) / (v + 1.0)
  let y2 = y * y
  let mut sum = y
  let mut term = y
  let mut i = 1
  while i <= 60 {
    term = term * y2
    let delta = term / (2 * i + 1).to_double()
    sum = sum + delta
    if delta.abs() < 0.0000000000001 {
      break
    }
    i = i + 1
  }
  sum * 2.0 * (1 << k).to_double()
}

///|
/// Exponential function. Uses range reduction plus Taylor series.
pub fn bsseq_exp(x : Double) -> Double {
  if x < -50.0 {
    return 0.0
  }
  if x > 50.0 {
    return 1000000000000000000000000000000.0
  }
  let mut v = x
  let mut k = 0
  while v.abs() > 0.5 && k < 40 {
    v = v / 2.0
    k = k + 1
  }
  let mut sum = 1.0
  let mut term = 1.0
  let mut i = 1
  while i <= 40 {
    term = term * v / i.to_double()
    sum = sum + term
    if term.abs() < 0.0000000000001 {
      break
    }
    i = i + 1
  }
  let mut result = sum
  while k > 0 {
    result = result * result
    k = k - 1
  }
  result
}

///|
/// Base-2 logarithm.
pub fn bsseq_log2(x : Double) -> Double {
  bsseq_ln(x) / bsseq_ln2
}

///|
/// Power via exp/ln.
pub fn bsseq_pow(base : Double, exp : Double) -> Double {
  if base <= 0.0 {
    return 0.0
  }
  bsseq_exp(exp * bsseq_ln(base))
}

///|
/// Lanczos approximation for ln(gamma(x)).
pub fn bsseq_ln_gamma(x : Double) -> Double {
  let cof = [
    76.18009172947146, -86.50532032941677, 24.01409824083091, -1.231739572450155,
    0.001208650973866179, -0.000005395239384953,
  ]
  let mut x1 = x
  let mut tmp = x1 + 5.5
  tmp = (x1 + 0.5) * bsseq_ln(tmp) - tmp
  let mut ser = 1.000000000190015
  let mut j = 0
  while j < cof.length() {
    x1 = x1 + 1.0
    ser = ser + cof[j] / x1
    j = j + 1
  }
  tmp + bsseq_ln(2.5066282746310005 * ser)
}

///|
/// Continued fraction for incomplete beta function.
fn bsseq_betacf(x : Double, a : Double, b : Double) -> Double {
  let max_iter = 100
  let eps = 0.0000003
  let mut c = 1.0
  let mut d = 1.0 - (a + b) * x / (a + 1.0)
  if d.abs() < eps {
    d = eps
  }
  d = 1.0 / d
  let mut h = d
  let mut m = 1
  while m <= max_iter {
    let m2 = (2 * m).to_double()
    let mut aa = m.to_double() *
      (b - m.to_double()) *
      x /
      ((a + m2 - 1.0) * (a + m2))
    d = 1.0 + aa * d
    if d.abs() < eps {
      d = eps
    }
    c = 1.0 + aa / c
    if c.abs() < eps {
      c = eps
    }
    d = 1.0 / d
    h = h * d * c
    aa = -(a + m.to_double()) *
      (a + b + m.to_double()) *
      x /
      ((a + m2) * (a + m2 + 1.0))
    d = 1.0 + aa * d
    if d.abs() < eps {
      d = eps
    }
    c = 1.0 + aa / c
    if c.abs() < eps {
      c = eps
    }
    d = 1.0 / d
    let del = d * c
    h = h * del
    if del.abs() < eps {
      break
    }
    m = m + 1
  }
  h
}

///|
/// Regularized incomplete beta function I_x(a, b).
pub fn bsseq_betainc(x : Double, a : Double, b : Double) -> Double {
  if x <= 0.0 {
    return 0.0
  }
  if x >= 1.0 {
    return 1.0
  }
  let lbeta_val = bsseq_ln_gamma(a) + bsseq_ln_gamma(b) - bsseq_ln_gamma(a + b)
  let mut bt = bsseq_ln(x) * a + bsseq_ln(1.0 - x) * b - lbeta_val
  bt = bsseq_exp(bt)
  if x < (a + 1.0) / (a + b + 2.0) {
    bt * bsseq_betacf(x, a, b) / a
  } else {
    1.0 - bt * bsseq_betacf(1.0 - x, b, a) / b
  }
}

///|
/// Standard normal CDF using Abramowitz and Stegun erf approximation.
pub fn bsseq_normal_cdf(z : Double) -> Double {
  let sign = if z < 0.0 { -1.0 } else { 1.0 }
  let x = z.abs() / bsseq_sqrt(2.0)
  let t = 1.0 / (1.0 + 0.3275911 * x)
  let erf = 1.0 -
    (
      0.254829592 * t -
      0.284496736 * t * t +
      1.421413741 * t * t * t -
      1.453152027 * t * t * t * t +
      1.061405429 * t * t * t * t * t
    ) *
    bsseq_exp(-x * x)
  0.5 * (1.0 + sign * erf)
}

///|
/// Two-tailed p-value for Student's t distribution.
/// p = P(|T| > |t|) for T ~ t(df).
pub fn bsseq_t_test_pval(t : Double, df : Double) -> Double {
  if df <= 0.0 {
    return 1.0
  }
  if t == 0.0 {
    return 1.0
  }
  let x = df / (df + t * t)
  let p = bsseq_betainc(x, df / 2.0, 0.5)
  if p < 0.0 {
    0.0
  } else if p > 1.0 {
    1.0
  } else {
    p
  }
}

///|
/// Copy and sort an array of doubles (ascending) using selection sort.
pub fn bsseq_sort_doubles(arr : Array[Double]) -> Array[Double] {
  let sorted = arr.copy()
  let n = sorted.length()
  let mut i = 0
  while i < n - 1 {
    let mut min_idx = i
    let mut j = i + 1
    while j < n {
      if sorted[j] < sorted[min_idx] {
        min_idx = j
      }
      j = j + 1
    }
    if min_idx != i {
      let tmp = sorted[i]
      sorted[i] = sorted[min_idx]
      sorted[min_idx] = tmp
    }
    i = i + 1
  }
  sorted
}

///|
/// Mean of an array of doubles.
pub fn bsseq_mean(arr : Array[Double]) -> Double {
  let n = arr.length()
  if n == 0 {
    return 0.0
  }
  let mut sum = 0.0
  let mut i = 0
  while i < n {
    sum = sum + arr[i]
    i = i + 1
  }
  sum / n.to_double()
}

///|
/// Median of an array of doubles.
pub fn bsseq_median(arr : Array[Double]) -> Double {
  let n = arr.length()
  if n == 0 {
    return 0.0
  }
  let sorted = bsseq_sort_doubles(arr)
  if n % 2 == 0 {
    (sorted[n / 2 - 1] + sorted[n / 2]) / 2.0
  } else {
    sorted[n / 2]
  }
}

///|
/// Sample variance (denominator n-1).
pub fn bsseq_variance(arr : Array[Double]) -> Double {
  let n = arr.length()
  if n < 2 {
    return 0.0
  }
  let m = bsseq_mean(arr)
  let mut ss = 0.0
  let mut i = 0
  while i < n {
    let d = arr[i] - m
    ss = ss + d * d
    i = i + 1
  }
  ss / (n - 1).to_double()
}

///|
/// Standard deviation.
pub fn bsseq_sd(arr : Array[Double]) -> Double {
  bsseq_sqrt(bsseq_variance(arr))
}

///|
/// Two-sample Welch's t-test.
/// Returns (t_statistic, degrees_of_freedom).
pub fn bsseq_t_test(x : Array[Double], y : Array[Double]) -> (Double, Double) {
  let nx = x.length()
  let ny = y.length()
  if nx < 2 || ny < 2 {
    return (0.0, 0.0)
  }
  let mx = bsseq_mean(x)
  let my = bsseq_mean(y)
  let vx = bsseq_variance(x)
  let vy = bsseq_variance(y)
  let se2 = vx / nx.to_double() + vy / ny.to_double()
  if se2 <= 0.0 {
    return (0.0, (nx + ny - 2).to_double())
  }
  let se = bsseq_sqrt(se2)
  let t = (mx - my) / se
  let num = se2 * se2
  let den = vx / nx.to_double() * (vx / nx.to_double()) / (nx - 1).to_double() +
    vy / ny.to_double() * (vy / ny.to_double()) / (ny - 1).to_double()
  let df = if den > 0.0 { num / den } else { (nx + ny - 2).to_double() }
  (t, df)
}

///|
/// Benjamini-Hochberg FDR adjustment.
pub fn bsseq_bh_fdr(pvalues : Array[Double]) -> Array[Double] {
  let n = pvalues.length()
  if n == 0 {
    return []
  }
  let indices = Array::make(n, 0)
  let mut i = 0
  while i < n {
    indices[i] = i
    i = i + 1
  }
  // Sort indices by p-value ascending (bubble sort).
  let mut swapped = true
  let mut k = n - 1
  while k > 0 && swapped {
    swapped = false
    let mut j = 0
    while j < k {
      if pvalues[indices[j]] > pvalues[indices[j + 1]] {
        let tmp = indices[j]
        indices[j] = indices[j + 1]
        indices[j + 1] = tmp
        swapped = true
      }
      j = j + 1
    }
    k = k - 1
  }
  let adj = Array::make(n, 0.0)
  let mut min_adj = 1.0
  let mut idx = n - 1
  while idx >= 0 {
    let rank = (idx + 1).to_double()
    let adj_val = pvalues[indices[idx]] * n.to_double() / rank
    if adj_val < min_adj {
      min_adj = adj_val
    }
    adj[indices[idx]] = min_adj
    idx = idx - 1
  }
  i = 0
  while i < n {
    if adj[i] > 1.0 {
      adj[i] = 1.0
    }
    i = i + 1
  }
  adj
}

///|
/// Trim leading/trailing whitespace from a string.
fn bsseq_trim(s : String) -> String {
  let n = s.length()
  if n == 0 {
    return ""
  }
  let mut start = 0
  let mut end = n
  while start < end {
    let c = s.unsafe_get(start).to_int()
    if c == 32 || c == 9 || c == 13 || c == 10 {
      start = start + 1
    } else {
      break
    }
  }
  while end > start {
    let c = s.unsafe_get(end - 1).to_int()
    if c == 32 || c == 9 || c == 13 || c == 10 {
      end = end - 1
    } else {
      break
    }
  }
  if start >= end {
    ""
  } else {
    s[start:end].to_string()
  }
}

///|
/// Check whether a line is a header / comment line.
fn bsseq_is_header_line(line : String) -> Bool {
  if line.length() == 0 {
    return true
  }
  let c = line.unsafe_get(0).to_int()
  if c == '#'.to_int() {
    return true
  }
  // Lines beginning with "chr" but containing "chrom" as a header word.
  if line.length() >= 5 {
    let prefix = line[0:5].to_string()
    if prefix == "chrom" || prefix == "Chrom" {
      return true
    }
  }
  false
}

///|
/// Compute methylation rate meth / (meth + unmeth). Returns 0 when total is 0.
pub fn bsseq_get_methylation_rate(meth : Int, unmeth : Int) -> Double {
  let total = meth + unmeth
  if total <= 0 {
    return 0.0
  }
  meth.to_double() / total.to_double()
}

///|
/// A single bisulfite sequencing methylation record.
pub struct BSRecord {
  chr : String
  pos : Int
  strand : String
  meth : Int
  unmeth : Int
  context : String
}

///|
pub fn BSRecord::new(
  chr : String,
  pos : Int,
  strand : String,
  meth : Int,
  unmeth : Int,
) -> BSRecord {
  BSRecord::{ chr, pos, strand, meth, unmeth, context: "CpG" }
}

///|
pub fn BSRecord::new_with_context(
  chr : String,
  pos : Int,
  strand : String,
  meth : Int,
  unmeth : Int,
  context : String,
) -> BSRecord {
  BSRecord::{ chr, pos, strand, meth, unmeth, context }
}

///|
/// Total coverage for this record.
pub fn BSRecord::coverage(self : BSRecord) -> Int {
  self.meth + self.unmeth
}

///|
/// Methylation rate in [0, 1].
pub fn BSRecord::rate(self : BSRecord) -> Double {
  bsseq_get_methylation_rate(self.meth, self.unmeth)
}

///|
pub fn BSRecord::to_string(self : BSRecord) -> String {
  "BSRecord(chr=" +
  self.chr +
  ", pos=" +
  self.pos.to_string() +
  ", strand=" +
  self.strand +
  ", meth=" +
  self.meth.to_string() +
  ", unmeth=" +
  self.unmeth.to_string() +
  ", context=" +
  self.context +
  ")"
}

///|
/// Container for a single-sample BS-seq dataset.
pub struct BSData {
  records : Array[BSRecord]
  sample_names : Array[String]
  chr_names : Array[String]
  coverage : Array[Int]
  methylation : Array[Double]
}

///|
pub fn BSData::new(
  records : Array[BSRecord],
  sample_names : Array[String],
) -> BSData {
  let coverage : Array[Int] = Array::new()
  let methylation : Array[Double] = Array::new()
  let chr_set : Map[String, Bool] = Map([], capacity=32)
  let mut i = 0
  while i < records.length() {
    let cov = records[i].meth + records[i].unmeth
    coverage.push(cov)
    methylation.push(
      bsseq_get_methylation_rate(records[i].meth, records[i].unmeth),
    )
    chr_set[records[i].chr] = true
    i = i + 1
  }
  let chr_names = chr_set.keys().collect()
  BSData::{ records, sample_names, chr_names, coverage, methylation }
}

///|
pub fn BSData::n_records(self : BSData) -> Int {
  self.records.length()
}

///|
pub fn BSData::n_samples(self : BSData) -> Int {
  self.sample_names.length()
}

///|
pub fn BSData::to_string(self : BSData) -> String {
  let mut result = "BSData\n"
  result = result + "------\n"
  result = result + "Records: " + self.n_records().to_string() + "\n"
  result = result + "Samples: " + self.n_samples().to_string() + "\n"
  result = result + "Sample names:\n"
  let mut i = 0
  while i < self.sample_names.length() && i < 5 {
    result = result + "  - " + self.sample_names[i] + "\n"
    i = i + 1
  }
  if self.sample_names.length() > 5 {
    result = result +
      "  ... and " +
      (self.sample_names.length() - 5).to_string() +
      " more\n"
  }
  result = result + "Chromosomes: "
  i = 0
  while i < self.chr_names.length() && i < 10 {
    if i > 0 {
      result = result + ", "
    }
    result = result + self.chr_names[i]
    i = i + 1
  }
  result = result + "\n"
  let mean_cov = if self.coverage.length() > 0 {
    let mut s = 0.0
    let mut j = 0
    while j < self.coverage.length() {
      s = s + self.coverage[j].to_double()
      j = j + 1
    }
    s / self.coverage.length().to_double()
  } else {
    0.0
  }
  let mean_meth = bsseq_mean(self.methylation)
  result = result + "Mean coverage: " + mean_cov.to_string() + "\n"
  result = result + "Mean methylation: " + mean_meth.to_string() + "\n"
  result
}

///|
/// Result of BSmooth local smoothing.
pub struct BSmoothResult {
  chr : String
  pos : Array[Int]
  smoothed_meth : Array[Double]
  raw_meth : Array[Double]
  coverage : Array[Int]
  window_size : Int
}

///|
pub fn BSmoothResult::to_string(self : BSmoothResult) -> String {
  let mut result = "BSmoothResult\n"
  result = result + "-------------\n"
  result = result + "Chromosome: " + self.chr + "\n"
  result = result + "CpGs: " + self.pos.length().to_string() + "\n"
  result = result + "Window size: " + self.window_size.to_string() + "\n"
  let mean_smoothed = bsseq_mean(self.smoothed_meth)
  let mean_raw = bsseq_mean(self.raw_meth)
  result = result +
    "Mean smoothed methylation: " +
    mean_smoothed.to_string() +
    "\n"
  result = result + "Mean raw methylation: " + mean_raw.to_string() + "\n"
  result
}

///|
/// A single differentially methylated region.
pub struct DMRegion {
  chr : String
  start : Int
  end : Int
  n_cpg : Int
  mean_meth_diff : Double
  pvalue : Double
  fdr : Double
  direction : String
}

///|
pub fn DMRegion::new(
  chr : String,
  start : Int,
  end : Int,
  n_cpg : Int,
  mean_meth_diff : Double,
  pvalue : Double,
  fdr : Double,
  direction : String,
) -> DMRegion {
  DMRegion::{ chr, start, end, n_cpg, mean_meth_diff, pvalue, fdr, direction }
}

///|
pub fn DMRegion::to_string(self : DMRegion) -> String {
  "DMRegion(chr=" +
  self.chr +
  ", start=" +
  self.start.to_string() +
  ", end=" +
  self.end.to_string() +
  ", n_cpg=" +
  self.n_cpg.to_string() +
  ", mean_diff=" +
  self.mean_meth_diff.to_string() +
  ", pvalue=" +
  self.pvalue.to_string() +
  ", fdr=" +
  self.fdr.to_string() +
  ", direction=" +
  self.direction +
  ")"
}

///|
/// Result of DMR detection.
/// Named BsseqDMRResult to avoid clashing with minfi's DMRResult in the same package.
pub struct BsseqDMRResult {
  regions : Array[DMRegion]
  n_dmr : Int
  n_significant : Int
  threshold : Double
}

///|
pub fn BsseqDMRResult::new(
  regions : Array[DMRegion],
  threshold : Double,
) -> BsseqDMRResult {
  let n_dmr = regions.length()
  let mut n_sig = 0
  let mut i = 0
  while i < regions.length() {
    if regions[i].fdr < threshold {
      n_sig = n_sig + 1
    }
    i = i + 1
  }
  BsseqDMRResult::{ regions, n_dmr, n_significant: n_sig, threshold }
}

///|
pub fn BsseqDMRResult::to_string(self : BsseqDMRResult) -> String {
  let mut result = "BsseqDMRResult\n"
  result = result + "-------------\n"
  result = result + "Total DMRs: " + self.n_dmr.to_string() + "\n"
  result = result +
    "Significant (FDR < " +
    self.threshold.to_string() +
    "): " +
    self.n_significant.to_string() +
    "\n"
  let mut i = 0
  while i < self.regions.length() && i < 10 {
    let r = self.regions[i]
    result = result +
      "  " +
      r.chr +
      ":" +
      r.start.to_string() +
      "-" +
      r.end.to_string() +
      " (n_cpg=" +
      r.n_cpg.to_string() +
      ", diff=" +
      r.mean_meth_diff.to_string() +
      ", p=" +
      r.pvalue.to_string() +
      ", fdr=" +
      r.fdr.to_string() +
      ", " +
      r.direction +
      ")\n"
    i = i + 1
  }
  if self.regions.length() > 10 {
    result = result +
      "  ... and " +
      (self.regions.length() - 10).to_string() +
      " more\n"
  }
  result
}

///|
/// Result of a per-position t-test.
pub struct TStatResult {
  t_stat : Array[Double]
  pvalue : Array[Double]
  fdr : Array[Double]
  mean_diff : Array[Double]
}

///|
pub fn TStatResult::to_string(self : TStatResult) -> String {
  let mut result = "TStatResult\n"
  result = result + "-----------\n"
  result = result + "CpGs: " + self.t_stat.length().to_string() + "\n"
  let mean_t = bsseq_mean(self.t_stat)
  let mean_diff = bsseq_mean(self.mean_diff)
  result = result + "Mean t-stat: " + mean_t.to_string() + "\n"
  result = result + "Mean diff: " + mean_diff.to_string() + "\n"
  let mut n_sig = 0
  let mut i = 0
  while i < self.pvalue.length() {
    if self.pvalue[i] < 0.05 {
      n_sig = n_sig + 1
    }
    i = i + 1
  }
  result = result + "Significant (p<0.05): " + n_sig.to_string() + "\n"
  result
}

///|
/// Coverage matrix across multiple samples.
pub struct CovMatrix {
  chr : Array[String]
  pos : Array[Int]
  cov : Array[Array[Int]]
  meth : Array[Array[Int]]
}

///|
pub fn CovMatrix::to_string(self : CovMatrix) -> String {
  let mut result = "CovMatrix\n"
  result = result + "---------\n"
  result = result + "Positions: " + self.pos.length().to_string() + "\n"
  let n_samples = if self.cov.length() > 0 { self.cov[0].length() } else { 0 }
  result = result + "Samples: " + n_samples.to_string() + "\n"
  let mut total_cov = 0
  let mut total_meth = 0
  let mut i = 0
  while i < self.cov.length() {
    let mut j = 0
    while j < self.cov[i].length() {
      total_cov = total_cov + self.cov[i][j]
      total_meth = total_meth + self.meth[i][j]
      j = j + 1
    }
    i = i + 1
  }
  result = result + "Total coverage: " + total_cov.to_string() + "\n"
  result = result + "Total methylated: " + total_meth.to_string() + "\n"
  result
}

///|
/// Parse a BED-like methylation file.
/// Each line: chrom, start, end, name, meth, unmeth, [strand], [context].
/// Columns are tab- or space-separated. Lines starting with '#' or "chrom" are skipped.
pub fn bsseq_parse_bed(content : String) -> BSData {
  let lines = split_by_char(content, '\n'.to_int())
  let records : Array[BSRecord] = Array::new()
  let mut i = 0
  while i < lines.length() {
    let raw = lines[i]
    let line = bsseq_trim(raw)
    if line.length() == 0 {
      i = i + 1
      continue
    }
    if bsseq_is_header_line(line) {
      i = i + 1
      continue
    }
    // Split by tab first; if only one column results, fall back to space.
    let mut parts = split_by_char(line, 9) // tab
    if parts.length() < 2 {
      parts = split_by_char(line, ' '.to_int())
    }
    if parts.length() < 6 {
      i = i + 1
      continue
    }
    let chr = parts[0]
    let start = parse_int(parts[1])
    let meth = parse_int(parts[4])
    let unmeth = parse_int(parts[5])
    let strand = if parts.length() >= 7 {
      let s = parts[6]
      if s == "+" || s == "-" {
        s
      } else {
        "+"
      }
    } else {
      "+"
    }
    let context = if parts.length() >= 8 {
      let c = parts[7]
      if c == "CpG" || c == "CHG" || c == "CHH" {
        c
      } else {
        "CpG"
      }
    } else {
      "CpG"
    }
    // Position: prefer start; if start is 0 and we have end, use end.
    let pos = if start > 0 {
      start
    } else if parts.length() >= 3 {
      parse_int(parts[2])
    } else {
      0
    }
    records.push(
      BSRecord::new_with_context(chr, pos, strand, meth, unmeth, context),
    )
    i = i + 1
  }
  BSData::new(records, ["sample1"])
}

///|
/// Parse Bismark methylation extractor output.
/// Each line: chrom, pos, strand, meth_count, unmeth_count, context, [trinucleotide].
/// Lines starting with '#' or "Bismark" are skipped.
pub fn bsseq_parse_bismark(content : String) -> BSData {
  let lines = split_by_char(content, '\n'.to_int())
  let records : Array[BSRecord] = Array::new()
  let mut i = 0
  while i < lines.length() {
    let raw = lines[i]
    let line = bsseq_trim(raw)
    if line.length() == 0 {
      i = i + 1
      continue
    }
    if bsseq_is_header_line(line) {
      i = i + 1
      continue
    }
    // Skip the Bismark methylation extractor preamble line.
    if line.length() >= 6 {
      let prefix = line[0:6].to_string()
      if prefix == "Bismar" {
        i = i + 1
        continue
      }
    }
    let mut parts = split_by_char(line, 9)
    if parts.length() < 2 {
      parts = split_by_char(line, ' '.to_int())
    }
    if parts.length() < 5 {
      i = i + 1
      continue
    }
    let chr = parts[0]
    let pos = parse_int(parts[1])
    let strand = if parts.length() >= 3 {
      let s = parts[2]
      if s == "+" || s == "-" {
        s
      } else {
        "+"
      }
    } else {
      "+"
    }
    let meth = parse_int(parts[3])
    let unmeth = parse_int(parts[4])
    let context = if parts.length() >= 6 {
      let c = parts[5]
      if c == "CpG" || c == "CHG" || c == "CHH" {
        c
      } else {
        "CpG"
      }
    } else {
      "CpG"
    }
    records.push(
      BSRecord::new_with_context(chr, pos, strand, meth, unmeth, context),
    )
    i = i + 1
  }
  BSData::new(records, ["sample1"])
}

///|
/// Return the per-record coverage array.
pub fn bsseq_get_coverage(bsdata : BSData) -> Array[Int] {
  let cov : Array[Int] = Array::new()
  let mut i = 0
  while i < bsdata.records.length() {
    let r = bsdata.records[i]
    cov.push(r.meth + r.unmeth)
    i = i + 1
  }
  cov
}

///|
/// Return the per-record methylation rate array (values in [0,1]).
pub fn bsseq_get_methylation(bsdata : BSData) -> Array[Double] {
  let meth : Array[Double] = Array::new()
  let mut i = 0
  while i < bsdata.records.length() {
    let r = bsdata.records[i]
    meth.push(bsseq_get_methylation_rate(r.meth, r.unmeth))
    i = i + 1
  }
  meth
}

///|
/// Filter records by total coverage (inclusive bounds).
pub fn bsseq_filter_by_coverage(
  bsdata : BSData,
  min_cov : Int,
  max_cov : Int,
) -> BSData {
  let filtered : Array[BSRecord] = Array::new()
  let mut i = 0
  while i < bsdata.records.length() {
    let r = bsdata.records[i]
    let cov = r.meth + r.unmeth
    if cov >= min_cov && cov <= max_cov {
      filtered.push(r)
    }
    i = i + 1
  }
  BSData::new(filtered, bsdata.sample_names.copy())
}

///|
/// Filter records by methylation context ("CpG", "CHG", or "CHH").
pub fn bsseq_filter_by_context(bsdata : BSData, context : String) -> BSData {
  let filtered : Array[BSRecord] = Array::new()
  let mut i = 0
  while i < bsdata.records.length() {
    let r = bsdata.records[i]
    if r.context == context {
      filtered.push(r)
    }
    i = i + 1
  }
  BSData::new(filtered, bsdata.sample_names.copy())
}

///|
/// Sliding-window smoothing of methylation rates.
/// `window` is the number of CpGs to include on each side of the center.
/// The methylation rate at position i is the coverage-weighted mean of
/// rates in [i-window, i+window].
pub fn bsseq_smooth_methylation(
  meth : Array[Int],
  cov : Array[Int],
  window : Int,
) -> Array[Double] {
  let n = meth.length()
  let result : Array[Double] = Array::make(n, 0.0)
  if n == 0 {
    return result
  }
  let w = if window < 0 { 0 } else { window }
  let mut i = 0
  while i < n {
    let mut sum_meth = 0.0
    let mut sum_cov = 0.0
    let lo = if i - w < 0 { 0 } else { i - w }
    let hi = if i + w >= n { n - 1 } else { i + w }
    let mut j = lo
    while j <= hi {
      let c = cov[j].to_double()
      if c > 0.0 {
        let rate = meth[j].to_double() / c
        sum_meth = sum_meth + rate * c
        sum_cov = sum_cov + c
      }
      j = j + 1
    }
    result[i] = if sum_cov > 0.0 { sum_meth / sum_cov } else { 0.0 }
    i = i + 1
  }
  result
}

///|
/// BSmooth smoothing for a single-sample BSData.
/// Groups records by chromosome, then applies coverage-weighted local smoothing
/// using a window of ±window_size/2 CpGs around each site.
pub fn bsseq_bsmooth(bsdata : BSData, window_size : Int) -> BSmoothResult {
  // If there is only one chromosome, smooth directly.
  if bsdata.chr_names.length() <= 1 {
    let chr = if bsdata.chr_names.length() == 1 {
      bsdata.chr_names[0]
    } else {
      ""
    }
    let n = bsdata.records.length()
    let pos : Array[Int] = Array::make(n, 0)
    let meth_arr : Array[Int] = Array::make(n, 0)
    let cov_arr : Array[Int] = Array::make(n, 0)
    let raw : Array[Double] = Array::make(n, 0.0)
    let mut i = 0
    while i < n {
      let r = bsdata.records[i]
      pos[i] = r.pos
      meth_arr[i] = r.meth
      cov_arr[i] = r.meth + r.unmeth
      raw[i] = bsseq_get_methylation_rate(r.meth, r.unmeth)
      i = i + 1
    }
    let smoothed = bsseq_smooth_methylation(meth_arr, cov_arr, window_size)
    return BSmoothResult::{
      chr,
      pos,
      smoothed_meth: smoothed,
      raw_meth: raw,
      coverage: cov_arr,
      window_size,
    }
  }
  // Multiple chromosomes: smooth the first chromosome for simplicity.
  let target_chr = bsdata.chr_names[0]
  let pos : Array[Int] = Array::new()
  let meth_arr : Array[Int] = Array::new()
  let cov_arr : Array[Int] = Array::new()
  let raw : Array[Double] = Array::new()
  let mut i = 0
  while i < bsdata.records.length() {
    let r = bsdata.records[i]
    if r.chr == target_chr {
      pos.push(r.pos)
      meth_arr.push(r.meth)
      cov_arr.push(r.meth + r.unmeth)
      raw.push(bsseq_get_methylation_rate(r.meth, r.unmeth))
    }
    i = i + 1
  }
  let smoothed = bsseq_smooth_methylation(meth_arr, cov_arr, window_size)
  BSmoothResult::{
    chr: target_chr,
    pos,
    smoothed_meth: smoothed,
    raw_meth: raw,
    coverage: cov_arr,
    window_size,
  }
}

///|
/// Compute per-position t-statistics for two groups.
/// `group1` and `group2` are parallel arrays of per-position mean methylation
/// values (e.g., smoothed methylation per CpG, averaged within each group).
/// Returns per-position t-statistic, p-value, BH-FDR, and mean difference.
pub fn bsseq_compute_tstat(
  group1 : Array[Double],
  group2 : Array[Double],
) -> TStatResult {
  let n = group1.length()
  if n == 0 || group2.length() != n {
    return TStatResult::{ t_stat: [], pvalue: [], fdr: [], mean_diff: [] }
  }
  let mean_diff : Array[Double] = Array::make(n, 0.0)
  let t_stat : Array[Double] = Array::make(n, 0.0)
  let pvalue : Array[Double] = Array::make(n, 1.0)
  let mut i = 0
  while i < n {
    mean_diff[i] = group1[i] - group2[i]
    i = i + 1
  }
  // Pooled SD across all positions: use variance of per-position differences.
  let var_diff = bsseq_variance(mean_diff)
  let sd_pooled = bsseq_sqrt(var_diff + 0.0000001)
  let df = if n > 1 { (n - 1).to_double() } else { 1.0 }
  i = 0
  while i < n {
    t_stat[i] = mean_diff[i] / sd_pooled
    pvalue[i] = bsseq_t_test_pval(t_stat[i], df)
    i = i + 1
  }
  let fdr = bsseq_bh_fdr(pvalue)
  TStatResult::{ t_stat, pvalue, fdr, mean_diff }
}

///|
/// Detect differentially methylated regions from per-position t-statistics.
/// Consecutive positions whose |t-stat| exceeds `threshold` and whose genomic
/// gap is <= `max_gap` are merged into a single DMR.
pub fn bsseq_find_dmr(
  tstat_result : TStatResult,
  positions : Array[Int],
  chr : String,
  threshold : Double,
  max_gap : Int,
) -> BsseqDMRResult {
  let n = tstat_result.t_stat.length()
  if n == 0 || positions.length() != n {
    return BsseqDMRResult::new([], threshold)
  }
  // Collect indices of significant positions, in input order.
  let sig_idx : Array[Int] = Array::new()
  let mut i = 0
  while i < n {
    if tstat_result.t_stat[i].abs() >= threshold {
      sig_idx.push(i)
    }
    i = i + 1
  }
  if sig_idx.length() == 0 {
    return BsseqDMRResult::new([], threshold)
  }
  // Sort significant indices by position (selection sort).
  let m = sig_idx.length()
  let mut a = 0
  while a < m - 1 {
    let mut min_idx = a
    let mut b = a + 1
    while b < m {
      if positions[sig_idx[b]] < positions[sig_idx[min_idx]] {
        min_idx = b
      }
      b = b + 1
    }
    if min_idx != a {
      let tmp = sig_idx[a]
      sig_idx[a] = sig_idx[min_idx]
      sig_idx[min_idx] = tmp
    }
    a = a + 1
  }
  // Merge consecutive significant positions into regions.
  let regions : Array[DMRegion] = Array::new()
  let mut start_idx = 0
  let mut k = 1
  while k <= m {
    let should_break = if k == m {
      true
    } else {
      let prev_pos = positions[sig_idx[k - 1]]
      let cur_pos = positions[sig_idx[k]]
      cur_pos - prev_pos > max_gap
    }
    if should_break {
      // Region spans sig_idx[start_idx .. k-1].
      let first = sig_idx[start_idx]
      let last = sig_idx[k - 1]
      let region_start = positions[first]
      let region_end = positions[last]
      let n_cpg = k - start_idx
      // Mean methylation difference over the region.
      let diffs : Array[Double] = Array::new()
      let pvals : Array[Double] = Array::new()
      let fdrs : Array[Double] = Array::new()
      let mut q = start_idx
      while q < k {
        diffs.push(tstat_result.mean_diff[sig_idx[q]])
        pvals.push(tstat_result.pvalue[sig_idx[q]])
        fdrs.push(tstat_result.fdr[sig_idx[q]])
        q = q + 1
      }
      let mean_diff = bsseq_mean(diffs)
      let min_p = {
        let mut v = 1.0
        let mut r = 0
        while r < pvals.length() {
          if pvals[r] < v {
            v = pvals[r]
          }
          r = r + 1
        }
        v
      }
      let min_fdr = {
        let mut v = 1.0
        let mut r = 0
        while r < fdrs.length() {
          if fdrs[r] < v {
            v = fdrs[r]
          }
          r = r + 1
        }
        v
      }
      let direction = if mean_diff > 0.0 {
        "hyper"
      } else if mean_diff < 0.0 {
        "hypo"
      } else {
        "none"
      }
      regions.push(
        DMRegion::new(
          chr, region_start, region_end, n_cpg, mean_diff, min_p, min_fdr, direction,
        ),
      )
      start_idx = k
    }
    k = k + 1
  }
  BsseqDMRResult::new(regions, threshold)
}

///|
/// Statistical test for differential methylation between two groups at
/// a single region. Inputs are parallel arrays of methylated read counts
/// and total coverage for each CpG in the region, for each group.
/// Returns (t_statistic, p_value).
pub fn bsseq_dmr_test(
  meth1 : Array[Int],
  cov1 : Array[Int],
  meth2 : Array[Int],
  cov2 : Array[Int],
) -> (Double, Double) {
  let n = meth1.length()
  if n == 0 || cov1.length() != n || meth2.length() != n || cov2.length() != n {
    return (0.0, 1.0)
  }
  let rates1 : Array[Double] = Array::make(n, 0.0)
  let rates2 : Array[Double] = Array::make(n, 0.0)
  let mut i = 0
  while i < n {
    rates1[i] = if cov1[i] > 0 {
      meth1[i].to_double() / cov1[i].to_double()
    } else {
      0.0
    }
    rates2[i] = if cov2[i] > 0 {
      meth2[i].to_double() / cov2[i].to_double()
    } else {
      0.0
    }
    i = i + 1
  }
  let (t, df) = bsseq_t_test(rates1, rates2)
  let p = bsseq_t_test_pval(t, df)
  (t, p)
}

///|
/// Return DMRs whose direction is "hyper" (increased methylation in group 1).
pub fn bsseq_find_hyper_dmr(dmr_result : BsseqDMRResult) -> Array[DMRegion] {
  let result : Array[DMRegion] = Array::new()
  let mut i = 0
  while i < dmr_result.regions.length() {
    if dmr_result.regions[i].direction == "hyper" {
      result.push(dmr_result.regions[i])
    }
    i = i + 1
  }
  result
}

///|
/// Return DMRs whose direction is "hypo" (decreased methylation in group 1).
pub fn bsseq_find_hypo_dmr(dmr_result : BsseqDMRResult) -> Array[DMRegion] {
  let result : Array[DMRegion] = Array::new()
  let mut i = 0
  while i < dmr_result.regions.length() {
    if dmr_result.regions[i].direction == "hypo" {
      result.push(dmr_result.regions[i])
    }
    i = i + 1
  }
  result
}

///|
/// Compute the difference in methylation between two rates (group1 - group2).
pub fn bsseq_compute_methylation_diff(meth1 : Double, meth2 : Double) -> Double {
  meth1 - meth2
}

///|
/// Merge adjacent CpG records on the same chromosome whose positions are
/// within `max_gap` base pairs. Methylated and unmethylated counts are summed.
/// Records on different chromosomes are never merged.
pub fn bsseq_merge_cpgs(
  records : Array[BSRecord],
  max_gap : Int,
) -> Array[BSRecord] {
  let n = records.length()
  if n == 0 {
    return []
  }
  // Copy and sort by (chr, pos) using selection sort.
  let sorted = records.copy()
  let mut i = 0
  while i < n - 1 {
    let mut min_idx = i
    let mut j = i + 1
    while j < n {
      let a_chr_lt = sorted[j].chr < sorted[min_idx].chr
      let a_chr_eq = sorted[j].chr == sorted[min_idx].chr
      let a_pos_lt = sorted[j].pos < sorted[min_idx].pos
      if a_chr_lt || (a_chr_eq && a_pos_lt) {
        min_idx = j
      }
      j = j + 1
    }
    if min_idx != i {
      let tmp = sorted[i]
      sorted[i] = sorted[min_idx]
      sorted[min_idx] = tmp
    }
    i = i + 1
  }
  let merged : Array[BSRecord] = Array::new()
  let mut cur = sorted[0]
  let mut k = 1
  while k < n {
    let r = sorted[k]
    if r.chr == cur.chr && r.pos - cur.pos <= max_gap {
      // Merge: sum counts. Keep earlier position; preserve strand/context.
      cur = BSRecord::new_with_context(
        cur.chr,
        cur.pos,
        cur.strand,
        cur.meth + r.meth,
        cur.unmeth + r.unmeth,
        cur.context,
      )
    } else {
      merged.push(cur)
      cur = r
    }
    k = k + 1
  }
  merged.push(cur)
  merged
}

///|
/// Compute the coverage-weighted mean methylation rate over a genomic region.
pub fn bsseq_region_methylation(
  bsdata : BSData,
  chr : String,
  start : Int,
  end : Int,
) -> Double {
  let mut sum_meth = 0.0
  let mut sum_cov = 0.0
  let mut i = 0
  while i < bsdata.records.length() {
    let r = bsdata.records[i]
    if r.chr == chr && r.pos >= start && r.pos <= end {
      let c = (r.meth + r.unmeth).to_double()
      if c > 0.0 {
        sum_meth = sum_meth + r.meth.to_double()
        sum_cov = sum_cov + c
      }
    }
    i = i + 1
  }
  if sum_cov > 0.0 {
    sum_meth / sum_cov
  } else {
    0.0
  }
}

///|
/// Normalize methylation values.
/// Supported methods:
/// - "none": no normalization (returns a copy).
/// - "global": rescale per-record methylation rates so the global mean
///   becomes 0.5, then recompute meth/unmeth counts preserving coverage.
pub fn bsseq_normalize(bsdata : BSData, method : String) -> BSData {
  if method == "none" {
    let copied : Array[BSRecord] = Array::new()
    let mut i = 0
    while i < bsdata.records.length() {
      let r = bsdata.records[i]
      copied.push(
        BSRecord::new_with_context(
          r.chr,
          r.pos,
          r.strand,
          r.meth,
          r.unmeth,
          r.context,
        ),
      )
      i = i + 1
    }
    return BSData::new(copied, bsdata.sample_names.copy())
  }
  if method == "global" {
    let mean_meth = bsseq_mean(bsdata.methylation)
    if mean_meth <= 0.0 || mean_meth >= 1.0 {
      let copied : Array[BSRecord] = Array::new()
      let mut i = 0
      while i < bsdata.records.length() {
        let r = bsdata.records[i]
        copied.push(
          BSRecord::new_with_context(
            r.chr,
            r.pos,
            r.strand,
            r.meth,
            r.unmeth,
            r.context,
          ),
        )
        i = i + 1
      }
      return BSData::new(copied, bsdata.sample_names.copy())
    }
    let scale = 0.5 / mean_meth
    let normalized : Array[BSRecord] = Array::new()
    let mut i = 0
    while i < bsdata.records.length() {
      let r = bsdata.records[i]
      let cov = r.meth + r.unmeth
      let raw_rate = bsseq_get_methylation_rate(r.meth, r.unmeth)
      let mut new_rate = raw_rate * scale
      if new_rate < 0.0 {
        new_rate = 0.0
      }
      if new_rate > 1.0 {
        new_rate = 1.0
      }
      let new_meth = (new_rate * cov.to_double() + 0.5).to_int()
      let new_unmeth = cov - new_meth
      normalized.push(
        BSRecord::new_with_context(
          r.chr,
          r.pos,
          r.strand,
          new_meth,
          new_unmeth,
          r.context,
        ),
      )
      i = i + 1
    }
    return BSData::new(normalized, bsdata.sample_names.copy())
  }
  // Unknown method: return a copy.
  let copied : Array[BSRecord] = Array::new()
  let mut i = 0
  while i < bsdata.records.length() {
    let r = bsdata.records[i]
    copied.push(
      BSRecord::new_with_context(
        r.chr,
        r.pos,
        r.strand,
        r.meth,
        r.unmeth,
        r.context,
      ),
    )
    i = i + 1
  }
  BSData::new(copied, bsdata.sample_names.copy())
}

///|
/// Build a human-readable summary of a BSData object.
pub fn bsseq_summary(bsdata : BSData) -> String {
  let mut result = "BS-seq Data Summary\n"
  result = result + "===================\n"
  result = result + "Records: " + bsdata.n_records().to_string() + "\n"
  result = result + "Samples: " + bsdata.n_samples().to_string() + "\n"
  result = result +
    "Chromosomes: " +
    bsdata.chr_names.length().to_string() +
    "\n"
  let mean_cov = bsseq_mean(
    {
      let arr : Array[Double] = Array::new()
      let mut i = 0
      while i < bsdata.coverage.length() {
        arr.push(bsdata.coverage[i].to_double())
        i = i + 1
      }
      arr
    },
  )
  let median_cov = {
    let arr : Array[Double] = Array::new()
    let mut i = 0
    while i < bsdata.coverage.length() {
      arr.push(bsdata.coverage[i].to_double())
      i = i + 1
    }
    bsseq_median(arr)
  }
  let mean_meth = bsseq_mean(bsdata.methylation)
  let median_meth = bsseq_median(bsdata.methylation)
  result = result + "Mean coverage: " + mean_cov.to_string() + "\n"
  result = result + "Median coverage: " + median_cov.to_string() + "\n"
  result = result + "Mean methylation: " + mean_meth.to_string() + "\n"
  result = result + "Median methylation: " + median_meth.to_string() + "\n"
  // Context breakdown.
  let mut n_cpg = 0
  let mut n_chg = 0
  let mut n_chh = 0
  let mut i = 0
  while i < bsdata.records.length() {
    match bsdata.records[i].context {
      "CpG" => n_cpg = n_cpg + 1
      "CHG" => n_chg = n_chg + 1
      "CHH" => n_chh = n_chh + 1
      _ => ()
    }
    i = i + 1
  }
  result = result + "CpG sites: " + n_cpg.to_string() + "\n"
  result = result + "CHG sites: " + n_chg.to_string() + "\n"
  result = result + "CHH sites: " + n_chh.to_string() + "\n"
  result
}

///|
/// Combine multiple single-sample BSData objects into a coverage matrix.
/// The output positions are the union of (chr, pos) pairs across all samples,
/// sorted by chromosome then position. Missing values are filled with 0.
pub fn bsseq_combine_samples(samples : Array[BSData]) -> CovMatrix {
  let n_samples = samples.length()
  if n_samples == 0 {
    return CovMatrix::{ chr: [], pos: [], cov: [], meth: [] }
  }
  // Build a map from "chr:pos" to (chr, pos) and per-sample lookup.
  let key_to_chr : Map[String, String] = Map([], capacity=1024)
  let key_to_pos : Map[String, Int] = Map([], capacity=1024)
  let keys : Array[String] = Array::new()
  let mut s = 0
  while s < n_samples {
    let sample = samples[s]
    let mut r = 0
    while r < sample.records.length() {
      let rec = sample.records[r]
      let key = rec.chr + ":" + rec.pos.to_string()
      if !key_to_chr.contains(key) {
        key_to_chr[key] = rec.chr
        key_to_pos[key] = rec.pos
        keys.push(key)
      }
      r = r + 1
    }
    s = s + 1
  }
  // Sort keys by (chr, pos). Selection sort on (chr_string, pos_int).
  let n_keys = keys.length()
  let mut i = 0
  while i < n_keys - 1 {
    let mut min_idx = i
    let mut j = i + 1
    while j < n_keys {
      let chr_j = key_to_chr[keys[j]]
      let chr_m = key_to_chr[keys[min_idx]]
      let pos_j = key_to_pos[keys[j]]
      let pos_m = key_to_pos[keys[min_idx]]
      let lt = chr_j < chr_m || (chr_j == chr_m && pos_j < pos_m)
      if lt {
        min_idx = j
      }
      j = j + 1
    }
    if min_idx != i {
      let tmp = keys[i]
      keys[i] = keys[min_idx]
      keys[min_idx] = tmp
    }
    i = i + 1
  }
  // Build output arrays.
  let chr_arr : Array[String] = Array::new()
  let pos_arr : Array[Int] = Array::new()
  let cov_arr : Array[Array[Int]] = Array::new()
  let meth_arr : Array[Array[Int]] = Array::new()
  let mut k = 0
  while k < n_keys {
    let key = keys[k]
    chr_arr.push(key_to_chr[key])
    pos_arr.push(key_to_pos[key])
    let target_chr = key_to_chr[key]
    let target_pos = key_to_pos[key]
    let cov_row : Array[Int] = Array::make(n_samples, 0)
    let meth_row : Array[Int] = Array::make(n_samples, 0)
    let mut s2 = 0
    while s2 < n_samples {
      let sample = samples[s2]
      let mut r2 = 0
      while r2 < sample.records.length() {
        let rec = sample.records[r2]
        if rec.chr == target_chr && rec.pos == target_pos {
          cov_row[s2] = rec.meth + rec.unmeth
          meth_row[s2] = rec.meth
          break
        }
        r2 = r2 + 1
      }
      s2 = s2 + 1
    }
    cov_arr.push(cov_row)
    meth_arr.push(meth_row)
    k = k + 1
  }
  CovMatrix::{ chr: chr_arr, pos: pos_arr, cov: cov_arr, meth: meth_arr }
}

///|
/// Create a small set of example BSRecords for testing.
pub fn bsseq_create_example_records() -> Array[BSRecord] {
  [
    BSRecord::new_with_context("chr1", 1000, "+", 8, 2, "CpG"),
    BSRecord::new_with_context("chr1", 1100, "+", 9, 1, "CpG"),
    BSRecord::new_with_context("chr1", 1200, "+", 5, 5, "CpG"),
    BSRecord::new_with_context("chr1", 1300, "+", 3, 7, "CpG"),
    BSRecord::new_with_context("chr1", 1400, "+", 7, 3, "CpG"),
    BSRecord::new_with_context("chr1", 1500, "+", 6, 4, "CpG"),
    BSRecord::new_with_context("chr1", 1600, "+", 4, 6, "CpG"),
    BSRecord::new_with_context("chr1", 1700, "+", 2, 8, "CpG"),
    BSRecord::new_with_context("chr1", 1800, "+", 8, 2, "CpG"),
    BSRecord::new_with_context("chr1", 1900, "+", 9, 1, "CpG"),
  ]
}

///|
/// Create an example BSData for testing.
pub fn bsseq_create_example_bsdata() -> BSData {
  BSData::new(bsseq_create_example_records(), ["sample1"])
}

///|
/// Create an example CovMatrix for testing.
pub fn bsseq_create_example_cov_matrix() -> CovMatrix {
  let s1 = bsseq_create_example_bsdata()
  let recs2 : Array[BSRecord] = [
    BSRecord::new_with_context("chr1", 1000, "+", 4, 6, "CpG"),
    BSRecord::new_with_context("chr1", 1100, "+", 3, 7, "CpG"),
    BSRecord::new_with_context("chr1", 1200, "+", 5, 5, "CpG"),
    BSRecord::new_with_context("chr1", 1300, "+", 8, 2, "CpG"),
    BSRecord::new_with_context("chr1", 1400, "+", 7, 3, "CpG"),
    BSRecord::new_with_context("chr1", 1500, "+", 2, 8, "CpG"),
    BSRecord::new_with_context("chr1", 1600, "+", 1, 9, "CpG"),
    BSRecord::new_with_context("chr1", 1700, "+", 9, 1, "CpG"),
    BSRecord::new_with_context("chr1", 1800, "+", 5, 5, "CpG"),
    BSRecord::new_with_context("chr1", 1900, "+", 6, 4, "CpG"),
  ]
  let s2 = BSData::new(recs2, ["sample2"])
  bsseq_combine_samples([s1, s2])
}
