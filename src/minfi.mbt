///|
/// minfi module - Illumina DNA methylation array preprocessing.
/// Inspired by Bioconductor minfi package.
///
/// Provides core utilities for 450K/EPIC array data:
/// - Red/Green channel container (RGChannelSet)
/// - Methylated/Unmethylated signal container (MethylSet)
/// - Beta/M-value container (GenomicRatioSet)
/// - Background correction (NOOB, Illumina)
/// - Normalization (quantile, functional)
/// - Detection p-values, QC reports, sex prediction
/// - Differentially methylated positions/regions (DMP/DMR) analysis
/// - Probe filtering and SNP probe removal

///|
/// Natural logarithm of 2.
let minfi_ln2 : Double = 0.69314718055994530942

///|
/// Newton's method for square root. Returns 0 for negative input.
pub fn minfi_sqrt(x : Double) -> Double {
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
pub fn minfi_ln(x : Double) -> Double {
  if x <= 0.0 {
    return 0.0
  }
  if x == 1.0 {
    return 0.0
  }
  // Bring x close to 1 via repeated sqrt: ln(x) = 2^k * ln(sqrt^k(x))
  let mut v = x
  let mut k = 0
  while (v > 2.0 || v < 0.5) && k < 40 {
    v = minfi_sqrt(v)
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
pub fn minfi_exp(x : Double) -> Double {
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
pub fn minfi_log2(x : Double) -> Double {
  minfi_ln(x) / minfi_ln2
}

///|
/// Power via exp/ln.
pub fn minfi_pow(base : Double, exp : Double) -> Double {
  if base <= 0.0 {
    return 0.0
  }
  minfi_exp(exp * minfi_ln(base))
}

///|
/// Lanczos approximation for ln(gamma(x)).
fn minfi_lngamma(x : Double) -> Double {
  let cof = [
    76.18009172947146, -86.50532032941677, 24.01409824083091, -1.231739572450155,
    0.001208650973866179, -0.000005395239384953,
  ]
  let mut x1 = x
  let mut tmp = x1 + 5.5
  tmp = (x1 + 0.5) * minfi_ln(tmp) - tmp
  let mut ser = 1.000000000190015
  let mut j = 0
  while j < cof.length() {
    x1 = x1 + 1.0
    ser = ser + cof[j] / x1
    j = j + 1
  }
  tmp + minfi_ln(2.5066282746310005 * ser)
}

///|
/// Continued fraction for incomplete beta function.
fn minfi_betacf(x : Double, a : Double, b : Double) -> Double {
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
fn minfi_betainc(x : Double, a : Double, b : Double) -> Double {
  if x <= 0.0 {
    return 0.0
  }
  if x >= 1.0 {
    return 1.0
  }
  let lbeta_val = minfi_lngamma(a) + minfi_lngamma(b) - minfi_lngamma(a + b)
  let mut bt = minfi_ln(x) * a + minfi_ln(1.0 - x) * b - lbeta_val
  bt = minfi_exp(bt)
  if x < (a + 1.0) / (a + b + 2.0) {
    bt * minfi_betacf(x, a, b) / a
  } else {
    1.0 - bt * minfi_betacf(1.0 - x, b, a) / b
  }
}

///|
/// Standard normal CDF using Abramowitz and Stegun erf approximation.
fn minfi_normal_cdf(z : Double) -> Double {
  let sign = if z < 0.0 { -1.0 } else { 1.0 }
  let x = z.abs() / minfi_sqrt(2.0)
  let t = 1.0 / (1.0 + 0.3275911 * x)
  let erf = 1.0 -
    (
      0.254829592 * t -
      0.284496736 * t * t +
      1.421413741 * t * t * t -
      1.453152027 * t * t * t * t +
      1.061405429 * t * t * t * t * t
    ) *
    minfi_exp(-x * x)
  0.5 * (1.0 + sign * erf)
}

///|
/// Two-tailed p-value for Student's t distribution.
/// p = P(|T| > |t|) for T ~ t(df).
pub fn minfi_t_test_pval(t : Double, df : Double) -> Double {
  if df <= 0.0 {
    return 1.0
  }
  if t == 0.0 {
    return 1.0
  }
  let x = df / (df + t * t)
  let p = minfi_betainc(x, df / 2.0, 0.5)
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
pub fn minfi_sort_doubles(arr : Array[Double]) -> Array[Double] {
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
/// Mean of an array.
pub fn minfi_mean(arr : Array[Double]) -> Double {
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
/// Median of an array.
pub fn minfi_median(arr : Array[Double]) -> Double {
  let n = arr.length()
  if n == 0 {
    return 0.0
  }
  let sorted = minfi_sort_doubles(arr)
  if n % 2 == 0 {
    (sorted[n / 2 - 1] + sorted[n / 2]) / 2.0
  } else {
    sorted[n / 2]
  }
}

///|
/// Sample variance (denominator n-1).
pub fn minfi_variance(arr : Array[Double]) -> Double {
  let n = arr.length()
  if n < 2 {
    return 0.0
  }
  let m = minfi_mean(arr)
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
pub fn minfi_sd(arr : Array[Double]) -> Double {
  minfi_sqrt(minfi_variance(arr))
}

///|
/// Return the q-th quantile (q in [0,1]) of an array using linear interpolation.
pub fn minfi_quantile_value(arr : Array[Double], q : Double) -> Double {
  let n = arr.length()
  if n == 0 {
    return 0.0
  }
  let sorted = minfi_sort_doubles(arr)
  let pos = q * (n - 1).to_double()
  let lo = pos.to_int()
  let hi = lo + 1
  if hi >= n {
    return sorted[n - 1]
  }
  let frac = pos - lo.to_double()
  sorted[lo] * (1.0 - frac) + sorted[hi] * frac
}

///|
/// Two-sample Welch's t-test.
/// Returns (t_statistic, degrees_of_freedom).
pub fn minfi_t_test(x : Array[Double], y : Array[Double]) -> (Double, Double) {
  let nx = x.length()
  let ny = y.length()
  if nx < 2 || ny < 2 {
    return (0.0, 0.0)
  }
  let mx = minfi_mean(x)
  let my = minfi_mean(y)
  let vx = minfi_variance(x)
  let vy = minfi_variance(y)
  let se2 = vx / nx.to_double() + vy / ny.to_double()
  if se2 <= 0.0 {
    return (0.0, (nx + ny - 2).to_double())
  }
  let se = minfi_sqrt(se2)
  let t = (mx - my) / se
  let num = se2 * se2
  let den = vx / nx.to_double() * (vx / nx.to_double()) / (nx - 1).to_double() +
    vy / ny.to_double() * (vy / ny.to_double()) / (ny - 1).to_double()
  let df = if den > 0.0 { num / den } else { (nx + ny - 2).to_double() }
  (t, df)
}

///|
/// Benjamini-Hochberg FDR adjustment.
pub fn minfi_bh_fdr(pvalues : Array[Double]) -> Array[Double] {
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
  let mut i = n - 1
  while i >= 0 {
    let rank = (i + 1).to_double()
    let adj_val = pvalues[indices[i]] * n.to_double() / rank
    if adj_val < min_adj {
      min_adj = adj_val
    }
    adj[indices[i]] = min_adj
    i = i - 1
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
/// Quantile normalize a probe x sample matrix across samples.
/// Each row is a probe; columns are samples. Returns normalized matrix
/// where each sample has the same per-column distribution.
pub fn minfi_quantile_normalize(
  data : Array[Array[Double]],
) -> Array[Array[Double]] {
  let n_probes = data.length()
  if n_probes == 0 {
    return data
  }
  let n_samples = data[0].length()
  if n_samples == 0 {
    return data
  }
  // Sort each sample column independently (column-major).
  let sorted_cols : Array[Array[Double]] = Array::new()
  let mut s = 0
  while s < n_samples {
    let col : Array[Double] = Array::make(n_probes, 0.0)
    let mut p = 0
    while p < n_probes {
      col[p] = data[p][s]
      p = p + 1
    }
    sorted_cols.push(minfi_sort_doubles(col))
    s = s + 1
  }
  // Compute the mean of each rank across samples.
  let rank_means : Array[Double] = Array::make(n_probes, 0.0)
  let mut r = 0
  while r < n_probes {
    let mut sum = 0.0
    let mut s = 0
    while s < n_samples {
      sum = sum + sorted_cols[s][r]
      s = s + 1
    }
    rank_means[r] = sum / n_samples.to_double()
    r = r + 1
  }
  // Assign each value its rank within its sample, then map to rank_means.
  let result : Array[Array[Double]] = Array::new()
  let mut p = 0
  while p < n_probes {
    let out_row : Array[Double] = Array::make(n_samples, 0.0)
    let mut s = 0
    while s < n_samples {
      // Compute rank of data[p][s] within sample s (0-indexed).
      let val = data[p][s]
      let mut lower = 0
      let mut equal = 0
      let mut q = 0
      while q < n_probes {
        if data[q][s] < val {
          lower = lower + 1
        } else if data[q][s] == val {
          equal = equal + 1
        }
        q = q + 1
      }
      // Average rank for ties.
      let rank = lower.to_double() + (equal - 1).to_double() / 2.0
      let lo = rank.to_int()
      let hi = lo + 1
      let frac = rank - lo.to_double()
      let v = if hi >= n_probes {
        rank_means[n_probes - 1]
      } else {
        rank_means[lo] * (1.0 - frac) + rank_means[hi] * frac
      }
      out_row[s] = v
      s = s + 1
    }
    result.push(out_row)
    p = p + 1
  }
  result
}

///|
pub struct RGChannelSet {
  red : Array[Array[Double]]
  green : Array[Array[Double]]
  sample_names : Array[String]
  probe_names : Array[String]
  pheno_data : Map[String, String]
}

///|
pub fn RGChannelSet::new(
  red : Array[Array[Double]],
  green : Array[Array[Double]],
  sample_names : Array[String],
  probe_names : Array[String],
) -> RGChannelSet {
  RGChannelSet::{
    red,
    green,
    sample_names,
    probe_names,
    pheno_data: Map([], capacity=8),
  }
}

///|
pub fn RGChannelSet::n_probes(self : RGChannelSet) -> Int {
  self.red.length()
}

///|
pub fn RGChannelSet::n_samples(self : RGChannelSet) -> Int {
  if self.red.length() > 0 {
    self.red[0].length()
  } else {
    0
  }
}

///|
pub fn RGChannelSet::to_string(self : RGChannelSet) -> String {
  let mut result = "RGChannelSet\n"
  result = result + "------------\n"
  result = result + "Probes: " + self.n_probes().to_string() + "\n"
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
  result
}

///|
pub struct MethylSet {
  methylated : Array[Array[Double]]
  unmethylated : Array[Array[Double]]
  sample_names : Array[String]
  probe_names : Array[String]
  pheno_data : Map[String, String]
}

///|
pub fn MethylSet::new(
  methylated : Array[Array[Double]],
  unmethylated : Array[Array[Double]],
  sample_names : Array[String],
  probe_names : Array[String],
) -> MethylSet {
  MethylSet::{
    methylated,
    unmethylated,
    sample_names,
    probe_names,
    pheno_data: Map([], capacity=8),
  }
}

///|
pub fn MethylSet::n_probes(self : MethylSet) -> Int {
  self.methylated.length()
}

///|
pub fn MethylSet::n_samples(self : MethylSet) -> Int {
  if self.methylated.length() > 0 {
    self.methylated[0].length()
  } else {
    0
  }
}

///|
pub fn MethylSet::to_string(self : MethylSet) -> String {
  let mut result = "MethylSet\n"
  result = result + "---------\n"
  result = result + "Probes: " + self.n_probes().to_string() + "\n"
  result = result + "Samples: " + self.n_samples().to_string() + "\n"
  // Median intensities summary.
  let mut s = 0
  while s < self.n_samples() && s < 5 {
    let m_vals : Array[Double] = Array::new()
    let u_vals : Array[Double] = Array::new()
    let mut p = 0
    while p < self.n_probes() {
      m_vals.push(self.methylated[p][s])
      u_vals.push(self.unmethylated[p][s])
      p = p + 1
    }
    result = result +
      "  " +
      self.sample_names[s] +
      ": medianM=" +
      minfi_median(m_vals).to_string() +
      ", medianU=" +
      minfi_median(u_vals).to_string() +
      "\n"
    s = s + 1
  }
  result
}

///|
pub struct GenomicRatioSet {
  mut beta : Array[Array[Double]]
  mut mvalues : Array[Array[Double]]
  sample_names : Array[String]
  probe_names : Array[String]
  mut annotation : Map[String, String]
  mut pheno_data : Map[String, String]
}

///|
pub fn GenomicRatioSet::new(
  beta : Array[Array[Double]],
  mvalues : Array[Array[Double]],
  sample_names : Array[String],
  probe_names : Array[String],
) -> GenomicRatioSet {
  GenomicRatioSet::{
    beta,
    mvalues,
    sample_names,
    probe_names,
    annotation: Map([], capacity=8),
    pheno_data: Map([], capacity=8),
  }
}

///|
pub fn GenomicRatioSet::n_probes(self : GenomicRatioSet) -> Int {
  self.beta.length()
}

///|
pub fn GenomicRatioSet::n_samples(self : GenomicRatioSet) -> Int {
  if self.beta.length() > 0 {
    self.beta[0].length()
  } else {
    0
  }
}

///|
pub fn GenomicRatioSet::to_string(self : GenomicRatioSet) -> String {
  let mut result = "GenomicRatioSet\n"
  result = result + "---------------\n"
  result = result + "Probes: " + self.n_probes().to_string() + "\n"
  result = result + "Samples: " + self.n_samples().to_string() + "\n"
  // Overall beta summary.
  let all_betas : Array[Double] = Array::new()
  let mut p = 0
  while p < self.n_probes() {
    let mut s = 0
    while s < self.n_samples() {
      all_betas.push(self.beta[p][s])
      s = s + 1
    }
    p = p + 1
  }
  if all_betas.length() > 0 {
    result = result +
      "Median beta: " +
      minfi_median(all_betas).to_string() +
      "\n"
    result = result + "Mean beta: " + minfi_mean(all_betas).to_string() + "\n"
  }
  result
}

///|
pub struct DMRRegion {
  chr : String
  start : Int
  end : Int
  n_probes : Int
  mean_beta_diff : Double
  pvalue : Double
  fdr : Double
}

///|
pub fn DMRRegion::new(
  chr : String,
  start : Int,
  end : Int,
  n_probes : Int,
  mean_beta_diff : Double,
  pvalue : Double,
  fdr : Double,
) -> DMRRegion {
  DMRRegion::{ chr, start, end, n_probes, mean_beta_diff, pvalue, fdr }
}

///|
pub struct DMRResult {
  regions : Array[DMRRegion]
  n_dmr : Int
  n_significant : Int
  threshold : Double
}

///|
pub fn DMRResult::new(
  regions : Array[DMRRegion],
  threshold : Double,
) -> DMRResult {
  let n_dmr = regions.length()
  let mut n_sig = 0
  let mut i = 0
  while i < regions.length() {
    if regions[i].fdr < threshold {
      n_sig = n_sig + 1
    }
    i = i + 1
  }
  DMRResult::{ regions, n_dmr, n_significant: n_sig, threshold }
}

///|
pub fn DMRResult::to_string(self : DMRResult) -> String {
  let mut result = "DMRResult\n"
  result = result + "--------\n"
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
      " (n=" +
      r.n_probes.to_string() +
      ", diff=" +
      r.mean_beta_diff.to_string() +
      ", fdr=" +
      r.fdr.to_string() +
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
pub struct DMPResult {
  probe_names : Array[String]
  chr : Array[String]
  pos : Array[Int]
  pvalues : Array[Double]
  fdrs : Array[Double]
  beta_diffs : Array[Double]
  n_dmp : Int
}

///|
pub fn DMPResult::new(
  probe_names : Array[String],
  chr : Array[String],
  pos : Array[Int],
  pvalues : Array[Double],
  fdrs : Array[Double],
  beta_diffs : Array[Double],
) -> DMPResult {
  let mut n_dmp = 0
  let mut i = 0
  while i < fdrs.length() {
    if fdrs[i] < 0.05 {
      n_dmp = n_dmp + 1
    }
    i = i + 1
  }
  DMPResult::{ probe_names, chr, pos, pvalues, fdrs, beta_diffs, n_dmp }
}

///|
pub fn DMPResult::empty() -> DMPResult {
  DMPResult::{
    probe_names: [],
    chr: [],
    pos: [],
    pvalues: [],
    fdrs: [],
    beta_diffs: [],
    n_dmp: 0,
  }
}

///|
pub fn DMPResult::to_string(self : DMPResult) -> String {
  let mut result = "DMPResult\n"
  result = result + "---------\n"
  result = result +
    "Total probes tested: " +
    self.probe_names.length().to_string() +
    "\n"
  result = result +
    "Significant DMPs (FDR < 0.05): " +
    self.n_dmp.to_string() +
    "\n"
  // Show top 10 by FDR.
  let n = self.probe_names.length()
  let indices = Array::make(n, 0)
  let mut i = 0
  while i < n {
    indices[i] = i
    i = i + 1
  }
  // Sort indices by FDR ascending (bubble sort).
  let mut swapped = true
  let mut k = n - 1
  while k > 0 && swapped {
    swapped = false
    let mut j = 0
    while j < k {
      if self.fdrs[indices[j]] > self.fdrs[indices[j + 1]] {
        let tmp = indices[j]
        indices[j] = indices[j + 1]
        indices[j + 1] = tmp
        swapped = true
      }
      j = j + 1
    }
    k = k - 1
  }
  let mut shown = 0
  i = 0
  while i < n && shown < 10 {
    let idx = indices[i]
    result = result +
      "  " +
      self.probe_names[idx] +
      " (" +
      self.chr[idx] +
      ":" +
      self.pos[idx].to_string() +
      ") p=" +
      self.pvalues[idx].to_string() +
      ", fdr=" +
      self.fdrs[idx].to_string() +
      ", diff=" +
      self.beta_diffs[idx].to_string() +
      "\n"
    shown = shown + 1
    i = i + 1
  }
  result
}

///|
pub struct QCReport {
  sample_names : Array[String]
  median_meth : Array[Double]
  median_unmeth : Array[Double]
  median_total : Array[Double]
  detection_pvalues : Array[Double]
  failed_probes : Array[Int]
  bisulfite_conversion : Array[Double]
  sex_predictions : Array[String]
}

///|
pub fn QCReport::new(
  sample_names : Array[String],
  median_meth : Array[Double],
  median_unmeth : Array[Double],
  median_total : Array[Double],
  detection_pvalues : Array[Double],
  failed_probes : Array[Int],
  bisulfite_conversion : Array[Double],
  sex_predictions : Array[String],
) -> QCReport {
  QCReport::{
    sample_names,
    median_meth,
    median_unmeth,
    median_total,
    detection_pvalues,
    failed_probes,
    bisulfite_conversion,
    sex_predictions,
  }
}

///|
pub fn QCReport::to_string(self : QCReport) -> String {
  let mut result = "QCReport\n"
  result = result + "--------\n"
  result = result + "Samples: " + self.sample_names.length().to_string() + "\n"
  let mut i = 0
  while i < self.sample_names.length() {
    result = result +
      "  " +
      self.sample_names[i] +
      ": medM=" +
      self.median_meth[i].to_string() +
      ", medU=" +
      self.median_unmeth[i].to_string() +
      ", medTotal=" +
      self.median_total[i].to_string() +
      ", detP=" +
      self.detection_pvalues[i].to_string() +
      ", failed=" +
      self.failed_probes[i].to_string() +
      ", bisulfite=" +
      self.bisulfite_conversion[i].to_string() +
      ", sex=" +
      self.sex_predictions[i] +
      "\n"
    i = i + 1
  }
  result
}

///|
/// Compute beta value = M / (M + U + offset).
pub fn minfi_compute_beta(
  methylated : Double,
  unmethylated : Double,
  offset : Double,
) -> Double {
  let denom = methylated + unmethylated + offset
  if denom <= 0.0 {
    return 0.0
  }
  let b = methylated / denom
  if b < 0.0 {
    0.0
  } else if b > 1.0 {
    1.0
  } else {
    b
  }
}

///|
/// Compute M value = log2(M / U). Uses offset of 1 to avoid log(0).
pub fn minfi_compute_mvalue(
  methylated : Double,
  unmethylated : Double,
) -> Double {
  let m = if methylated < 1.0 { 1.0 } else { methylated }
  let u = if unmethylated < 1.0 { 1.0 } else { unmethylated }
  minfi_log2(m / u)
}

///|
/// NOOB (Normal-Exponential Out-of-Band) background correction.
/// For each sample, estimates background as the 5th percentile of channel
/// intensities, subtracts it from each probe, and floors at 1.0.
pub fn minfi_preprocess_noob(rgset : RGChannelSet) -> MethylSet {
  let n_probes = rgset.red.length()
  let n_samples = if n_probes > 0 { rgset.red[0].length() } else { 0 }
  // Per-sample background.
  let red_bg : Array[Double] = Array::make(n_samples, 0.0)
  let green_bg : Array[Double] = Array::make(n_samples, 0.0)
  let mut s = 0
  while s < n_samples {
    let red_vals : Array[Double] = Array::make(n_probes, 0.0)
    let green_vals : Array[Double] = Array::make(n_probes, 0.0)
    let mut p = 0
    while p < n_probes {
      red_vals[p] = rgset.red[p][s]
      green_vals[p] = rgset.green[p][s]
      p = p + 1
    }
    red_bg[s] = minfi_quantile_value(red_vals, 0.05)
    green_bg[s] = minfi_quantile_value(green_vals, 0.05)
    s = s + 1
  }
  let methylated : Array[Array[Double]] = Array::new()
  let unmethylated : Array[Array[Double]] = Array::new()
  let mut p = 0
  while p < n_probes {
    let m_row : Array[Double] = Array::make(n_samples, 0.0)
    let u_row : Array[Double] = Array::make(n_samples, 0.0)
    let mut s = 0
    while s < n_samples {
      let m_val = rgset.green[p][s] - green_bg[s]
      let u_val = rgset.red[p][s] - red_bg[s]
      m_row[s] = if m_val < 1.0 { 1.0 } else { m_val }
      u_row[s] = if u_val < 1.0 { 1.0 } else { u_val }
      s = s + 1
    }
    methylated.push(m_row)
    unmethylated.push(u_row)
    p = p + 1
  }
  MethylSet::new(
    methylated,
    unmethylated,
    rgset.sample_names,
    rgset.probe_names,
  )
}

///|
/// Illumina background correction.
/// Subtracts a per-sample mean-background estimated from low-intensity probes.
pub fn minfi_preprocess_illumina(rgset : RGChannelSet) -> MethylSet {
  let n_probes = rgset.red.length()
  let n_samples = if n_probes > 0 { rgset.red[0].length() } else { 0 }
  let red_bg : Array[Double] = Array::make(n_samples, 0.0)
  let green_bg : Array[Double] = Array::make(n_samples, 0.0)
  let mut s = 0
  while s < n_samples {
    let red_low : Array[Double] = Array::new()
    let green_low : Array[Double] = Array::new()
    let mut p = 0
    while p < n_probes {
      red_low.push(rgset.red[p][s])
      green_low.push(rgset.green[p][s])
      p = p + 1
    }
    // Background mean estimated from bottom 10% of signals.
    let red_q = minfi_quantile_value(red_low, 0.10)
    let green_q = minfi_quantile_value(green_low, 0.10)
    let mut red_sum = 0.0
    let mut red_n = 0
    let mut green_sum = 0.0
    let mut green_n = 0
    p = 0
    while p < n_probes {
      if rgset.red[p][s] <= red_q {
        red_sum = red_sum + rgset.red[p][s]
        red_n = red_n + 1
      }
      if rgset.green[p][s] <= green_q {
        green_sum = green_sum + rgset.green[p][s]
        green_n = green_n + 1
      }
      p = p + 1
    }
    red_bg[s] = if red_n > 0 { red_sum / red_n.to_double() } else { 0.0 }
    green_bg[s] = if green_n > 0 {
      green_sum / green_n.to_double()
    } else {
      0.0
    }
    s = s + 1
  }
  let methylated : Array[Array[Double]] = Array::new()
  let unmethylated : Array[Array[Double]] = Array::new()
  let mut p = 0
  while p < n_probes {
    let m_row : Array[Double] = Array::make(n_samples, 0.0)
    let u_row : Array[Double] = Array::make(n_samples, 0.0)
    let mut s = 0
    while s < n_samples {
      let m_val = rgset.green[p][s] - green_bg[s]
      let u_val = rgset.red[p][s] - red_bg[s]
      m_row[s] = if m_val < 1.0 { 1.0 } else { m_val }
      u_row[s] = if u_val < 1.0 { 1.0 } else { u_val }
      s = s + 1
    }
    methylated.push(m_row)
    unmethylated.push(u_row)
    p = p + 1
  }
  MethylSet::new(
    methylated,
    unmethylated,
    rgset.sample_names,
    rgset.probe_names,
  )
}

///|
/// Quantile normalization preprocessing.
/// Builds a MethylSet directly from channels (Type II convention:
/// methylated = green, unmethylated = red) then quantile-normalizes M and U
/// independently across samples.
pub fn minfi_preprocess_quantile(rgset : RGChannelSet) -> MethylSet {
  let n_probes = rgset.red.length()
  let n_samples = if n_probes > 0 { rgset.red[0].length() } else { 0 }
  let methylated : Array[Array[Double]] = Array::new()
  let unmethylated : Array[Array[Double]] = Array::new()
  let mut p = 0
  while p < n_probes {
    let m_row : Array[Double] = Array::make(n_samples, 0.0)
    let u_row : Array[Double] = Array::make(n_samples, 0.0)
    let mut s = 0
    while s < n_samples {
      let m_val = rgset.green[p][s]
      let u_val = rgset.red[p][s]
      m_row[s] = if m_val < 1.0 { 1.0 } else { m_val }
      u_row[s] = if u_val < 1.0 { 1.0 } else { u_val }
      s = s + 1
    }
    methylated.push(m_row)
    unmethylated.push(u_row)
    p = p + 1
  }
  let norm_m = minfi_quantile_normalize(methylated)
  let norm_u = minfi_quantile_normalize(unmethylated)
  MethylSet::new(norm_m, norm_u, rgset.sample_names, rgset.probe_names)
}

///|
/// Functional normalization.
/// Performs NOOB background correction, converts to ratio set, then
/// quantile-normalizes the M-values across samples.
pub fn minfi_preprocess_funnorm(rgset : RGChannelSet) -> GenomicRatioSet {
  let mset = minfi_preprocess_noob(rgset)
  let grset = minfi_ratio_convert(mset)
  // Quantile-normalize M-values across samples.
  grset.mvalues = minfi_quantile_normalize(grset.mvalues)
  // Recompute beta from normalized M-values.
  let n_probes = grset.mvalues.length()
  let n_samples = if n_probes > 0 { grset.mvalues[0].length() } else { 0 }
  let new_beta : Array[Array[Double]] = Array::new()
  let mut p = 0
  while p < n_probes {
    let row : Array[Double] = Array::make(n_samples, 0.0)
    let mut s = 0
    while s < n_samples {
      // beta = 2^M / (1 + 2^M)
      let m = grset.mvalues[p][s]
      let two_m = minfi_pow(2.0, m)
      let b = two_m / (1.0 + two_m)
      row[s] = if b < 0.0 { 0.0 } else if b > 1.0 { 1.0 } else { b }
      s = s + 1
    }
    new_beta.push(row)
    p = p + 1
  }
  grset.beta = new_beta
  grset
}

///|
/// Convert a MethylSet to a GenomicRatioSet by computing beta and M values.
pub fn minfi_ratio_convert(mset : MethylSet) -> GenomicRatioSet {
  let n_probes = mset.methylated.length()
  let n_samples = if n_probes > 0 { mset.methylated[0].length() } else { 0 }
  let beta : Array[Array[Double]] = Array::new()
  let mvalues : Array[Array[Double]] = Array::new()
  let offset = 100.0
  let mut p = 0
  while p < n_probes {
    let beta_row : Array[Double] = Array::make(n_samples, 0.0)
    let m_row : Array[Double] = Array::make(n_samples, 0.0)
    let mut s = 0
    while s < n_samples {
      let m_val = mset.methylated[p][s]
      let u_val = mset.unmethylated[p][s]
      beta_row[s] = minfi_compute_beta(m_val, u_val, offset)
      m_row[s] = minfi_compute_mvalue(m_val, u_val)
      s = s + 1
    }
    beta.push(beta_row)
    mvalues.push(m_row)
    p = p + 1
  }
  let grset = GenomicRatioSet::new(
    beta,
    mvalues,
    mset.sample_names,
    mset.probe_names,
  )
  grset
}

///|
/// Get the beta value matrix from a GenomicRatioSet.
pub fn minfi_get_beta(grset : GenomicRatioSet) -> Array[Array[Double]] {
  grset.beta
}

///|
/// Get the M-value matrix from a GenomicRatioSet.
pub fn minfi_get_mvalues(grset : GenomicRatioSet) -> Array[Array[Double]] {
  grset.mvalues
}

///|
/// Compute detection p-values for each probe.
/// Returns one p-value per probe = mean across samples of the right-tailed
/// normal p-value of (max(red, green) - bg_mean) / bg_sd, where the
/// background is estimated per-sample from the 5th percentile intensities.
pub fn minfi_detection_pvalue(rgset : RGChannelSet) -> Array[Double] {
  let n_probes = rgset.red.length()
  let n_samples = if n_probes > 0 { rgset.red[0].length() } else { 0 }
  let bg_mean : Array[Double] = Array::make(n_samples, 0.0)
  let bg_sd : Array[Double] = Array::make(n_samples, 0.0)
  let mut s = 0
  while s < n_samples {
    let low : Array[Double] = Array::new()
    let mut p = 0
    while p < n_probes {
      let m = if rgset.red[p][s] < rgset.green[p][s] {
        rgset.red[p][s]
      } else {
        rgset.green[p][s]
      }
      low.push(m)
      p = p + 1
    }
    let q = minfi_quantile_value(low, 0.05)
    let bg_vals : Array[Double] = Array::new()
    p = 0
    while p < n_probes {
      let m = if rgset.red[p][s] < rgset.green[p][s] {
        rgset.red[p][s]
      } else {
        rgset.green[p][s]
      }
      if m <= q {
        bg_vals.push(m)
      }
      p = p + 1
    }
    bg_mean[s] = minfi_mean(bg_vals)
    let sd = minfi_sd(bg_vals)
    bg_sd[s] = if sd < 1.0 { 1.0 } else { sd }
    s = s + 1
  }
  let result : Array[Double] = Array::new()
  let mut p = 0
  while p < n_probes {
    let mut sum_p = 0.0
    let mut s = 0
    while s < n_samples {
      let sig = if rgset.red[p][s] > rgset.green[p][s] {
        rgset.red[p][s]
      } else {
        rgset.green[p][s]
      }
      let z = (sig - bg_mean[s]) / bg_sd[s]
      // Right-tailed p-value: P(Z > z) = 1 - Phi(z).
      let pv = 1.0 - minfi_normal_cdf(z)
      sum_p = sum_p + pv
      s = s + 1
    }
    let mean_p = if n_samples > 0 { sum_p / n_samples.to_double() } else { 1.0 }
    result.push(
      if mean_p < 0.0 {
        0.0
      } else if mean_p > 1.0 {
        1.0
      } else {
        mean_p
      },
    )
    p = p + 1
  }
  result
}

///|
/// Generate a QC report for an RGChannelSet.
pub fn minfi_get_qc(rgset : RGChannelSet) -> QCReport {
  let n_probes = rgset.red.length()
  let n_samples = if n_probes > 0 { rgset.red[0].length() } else { 0 }
  let det_p = minfi_detection_pvalue(rgset)
  // Build a temporary MethylSet for median M/U.
  let mset = minfi_preprocess_noob(rgset)
  let median_meth : Array[Double] = Array::make(n_samples, 0.0)
  let median_unmeth : Array[Double] = Array::make(n_samples, 0.0)
  let median_total : Array[Double] = Array::make(n_samples, 0.0)
  let detection_pvalues : Array[Double] = Array::make(n_samples, 0.0)
  let failed_probes : Array[Int] = Array::make(n_samples, 0)
  let bisulfite_conversion : Array[Double] = Array::make(n_samples, 0.0)
  let mut s = 0
  while s < n_samples {
    let m_vals : Array[Double] = Array::new()
    let u_vals : Array[Double] = Array::new()
    let t_vals : Array[Double] = Array::new()
    let mut p = 0
    while p < n_probes {
      m_vals.push(mset.methylated[p][s])
      u_vals.push(mset.unmethylated[p][s])
      t_vals.push(mset.methylated[p][s] + mset.unmethylated[p][s])
      p = p + 1
    }
    median_meth[s] = minfi_median(m_vals)
    median_unmeth[s] = minfi_median(u_vals)
    median_total[s] = minfi_median(t_vals)
    // Per-sample mean detection p-value.
    let mut psum = 0.0
    let mut failed = 0
    p = 0
    while p < n_probes {
      psum = psum + det_p[p]
      if det_p[p] > 0.01 {
        failed = failed + 1
      }
      p = p + 1
    }
    detection_pvalues[s] = if n_probes > 0 {
      psum / n_probes.to_double()
    } else {
      1.0
    }
    failed_probes[s] = failed
    // Bisulfite conversion heuristic: ratio of median M to (M + U).
    let denom = median_meth[s] + median_unmeth[s]
    bisulfite_conversion[s] = if denom > 0.0 {
      median_meth[s] / denom
    } else {
      0.0
    }
    s = s + 1
  }
  // Build a quick GenomicRatioSet for sex prediction.
  let grset = minfi_ratio_convert(mset)
  let sex_predictions = minfi_estimate_sex(grset)
  QCReport::new(
    rgset.sample_names,
    median_meth,
    median_unmeth,
    median_total,
    detection_pvalues,
    failed_probes,
    bisulfite_conversion,
    sex_predictions,
  )
}

///|
/// Estimate sample sex from chrX-like probe methylation.
/// Uses probes whose name contains "chrX" if present; otherwise falls back
/// to the lowest-beta probes (assumed to be on chrX). Female: mean beta
/// near 0.5 (X-inactivation). Male: mean beta skewed (one X).
pub fn minfi_estimate_sex(grset : GenomicRatioSet) -> Array[String] {
  let n_probes = grset.beta.length()
  let n_samples = grset.n_samples()
  // Collect chrX probes by name.
  let chrx_indices : Array[Int] = Array::new()
  let mut p = 0
  while p < n_probes {
    let name = grset.probe_names[p]
    if name.contains("chrX") || name.contains("X_q") || name.contains("_X_") {
      chrx_indices.push(p)
    }
    p = p + 1
  }
  let result : Array[String] = Array::new()
  let mut s = 0
  while s < n_samples {
    let vals : Array[Double] = Array::new()
    if chrx_indices.length() > 0 {
      let mut i = 0
      while i < chrx_indices.length() {
        vals.push(grset.beta[chrx_indices[i]][s])
        i = i + 1
      }
    } else {
      // Fallback: use the lowest 5% beta probes for this sample.
      let all : Array[Double] = Array::new()
      let mut p = 0
      while p < n_probes {
        all.push(grset.beta[p][s])
        p = p + 1
      }
      let sorted = minfi_sort_doubles(all)
      let cutoff = n_probes / 20
      let mut i = 0
      while i < cutoff && i < sorted.length() {
        vals.push(sorted[i])
        i = i + 1
      }
    }
    let mean_beta = if vals.length() > 0 { minfi_mean(vals) } else { 0.5 }
    let threshold = 0.4
    let sex = if mean_beta < threshold {
      "M"
    } else if mean_beta > threshold && mean_beta < 0.6 {
      "F"
    } else {
      "F"
    }
    result.push(sex)
    s = s + 1
  }
  result
}

///|
/// Differentially methylated positions analysis.
/// For each probe, runs a two-sample Welch t-test between group 0 and group 1
/// (as defined in the groups array, indexed by sample). Returns p-values,
/// BH-FDR adjusted p-values, and per-probe beta differences (mean1 - mean0).
pub fn minfi_find_dmps(
  grset : GenomicRatioSet,
  groups : Array[Int],
) -> DMPResult {
  let n_probes = grset.beta.length()
  let n_samples = grset.n_samples()
  if n_samples != groups.length() {
    return DMPResult::empty()
  }
  // Find the two groups actually present.
  let group_a = 0
  let group_b = 1
  let probe_names : Array[String] = Array::new()
  let chr_arr : Array[String] = Array::new()
  let pos_arr : Array[Int] = Array::new()
  let pvalues : Array[Double] = Array::new()
  let beta_diffs : Array[Double] = Array::new()
  let mut p = 0
  while p < n_probes {
    let g0_vals : Array[Double] = Array::new()
    let g1_vals : Array[Double] = Array::new()
    let mut s = 0
    while s < n_samples {
      if groups[s] == group_a {
        g0_vals.push(grset.beta[p][s])
      } else if groups[s] == group_b {
        g1_vals.push(grset.beta[p][s])
      }
      s = s + 1
    }
    let (pval, diff) = if g0_vals.length() >= 2 && g1_vals.length() >= 2 {
      let (t, df) = minfi_t_test(g0_vals, g1_vals)
      let pv = minfi_t_test_pval(t, df)
      let m0 = minfi_mean(g0_vals)
      let m1 = minfi_mean(g1_vals)
      (pv, m1 - m0)
    } else {
      (1.0, 0.0)
    }
    pvalues.push(pval)
    beta_diffs.push(diff)
    probe_names.push(grset.probe_names[p])
    // Derive a pseudo chromosome and position from the probe name when possible.
    let (chr_s, pos_i) = minfi_parse_probe_loc(grset.probe_names[p], p)
    chr_arr.push(chr_s)
    pos_arr.push(pos_i)
    p = p + 1
  }
  let fdrs = minfi_bh_fdr(pvalues)
  DMPResult::new(probe_names, chr_arr, pos_arr, pvalues, fdrs, beta_diffs)
}

///|
/// Parse a pseudo genomic location from a probe name.
/// Names like "cg00000001_chr1_1000" yield (chr1, 1000); otherwise a
/// synthetic location is generated so downstream DMR merging still works.
fn minfi_parse_probe_loc(
  probe_name : String,
  fallback_idx : Int,
) -> (String, Int) {
  let parts = split_by_char(probe_name, '_'.to_int())
  if parts.length() >= 3 {
    let chr_part = parts[parts.length() - 2]
    let pos_part = parts[parts.length() - 1]
    let pos = match parse_double(pos_part) {
      Some(v) => v.to_int()
      None => fallback_idx * 1000 + 1
    }
    (chr_part, pos)
  } else {
    ("chr1", fallback_idx * 1000 + 1)
  }
}

///|
/// Differentially methylated regions analysis.
/// Groups significant DMPs (FDR < threshold) by chromosome, sorts by position,
/// then merges DMPs whose adjacent positions are within max_gap base pairs.
pub fn minfi_find_dmrs(dmp_result : DMPResult, max_gap : Int) -> DMRResult {
  let n = dmp_result.probe_names.length()
  // Group significant probe indices by chromosome.
  let chr_to_indices : Map[String, Array[Int]] = Map([], capacity=24)
  let mut i = 0
  while i < n {
    if dmp_result.fdrs[i] < 0.05 {
      let c = dmp_result.chr[i]
      match chr_to_indices.get(c) {
        Some(a) => a.push(i)
        None => {
          let a : Array[Int] = [i]
          chr_to_indices[c] = a
        }
      }
    }
    i = i + 1
  }
  let regions : Array[DMRRegion] = Array::new()
  let chr_keys = chr_to_indices.keys().collect()
  let mut ci = 0
  while ci < chr_keys.length() {
    let c = chr_keys[ci]
    let indices = match chr_to_indices.get(c) {
      Some(a) => a
      None => {
        ci = ci + 1
        continue
      }
    }
    // Sort indices by position (selection sort).
    let m = indices.length()
    let mut a = 0
    while a < m - 1 {
      let mut min_idx = a
      let mut b = a + 1
      while b < m {
        if dmp_result.pos[indices[b]] < dmp_result.pos[indices[min_idx]] {
          min_idx = b
        }
        b = b + 1
      }
      if min_idx != a {
        let tmp = indices[a]
        indices[a] = indices[min_idx]
        indices[min_idx] = tmp
      }
      a = a + 1
    }
    // Walk through and merge adjacent within max_gap.
    let mut idx = 0
    while idx < m {
      let start_idx = indices[idx]
      let mut end_idx = start_idx
      let mut n_probes_in_region = 1
      let mut sum_diff = dmp_result.beta_diffs[start_idx]
      let mut min_p = dmp_result.pvalues[start_idx]
      let mut min_fdr = dmp_result.fdrs[start_idx]
      while idx + 1 < m {
        let curr_pos = dmp_result.pos[indices[idx]]
        let next_pos = dmp_result.pos[indices[idx + 1]]
        if next_pos - curr_pos <= max_gap {
          idx = idx + 1
          end_idx = indices[idx]
          n_probes_in_region = n_probes_in_region + 1
          sum_diff = sum_diff + dmp_result.beta_diffs[end_idx]
          if dmp_result.pvalues[end_idx] < min_p {
            min_p = dmp_result.pvalues[end_idx]
          }
          if dmp_result.fdrs[end_idx] < min_fdr {
            min_fdr = dmp_result.fdrs[end_idx]
          }
        } else {
          break
        }
      }
      regions.push(
        DMRRegion::new(
          c,
          dmp_result.pos[start_idx],
          dmp_result.pos[end_idx],
          n_probes_in_region,
          sum_diff / n_probes_in_region.to_double(),
          min_p,
          min_fdr,
        ),
      )
      idx = idx + 1
    }
    ci = ci + 1
  }
  DMRResult::new(regions, 0.05)
}

///|
/// Filter low-quality probes.
/// Drops probes whose beta-value standard deviation across samples is below
/// the given threshold (low-variance probes are uninformative for
/// differential analysis).
pub fn minfi_filter_probes(
  grset : GenomicRatioSet,
  threshold : Double,
) -> GenomicRatioSet {
  let n_probes = grset.beta.length()
  let new_beta : Array[Array[Double]] = Array::new()
  let new_mvalues : Array[Array[Double]] = Array::new()
  let new_probe_names : Array[String] = Array::new()
  let mut p = 0
  while p < n_probes {
    let row = grset.beta[p]
    let sd = minfi_sd(row)
    if sd >= threshold {
      new_beta.push(grset.beta[p])
      new_mvalues.push(grset.mvalues[p])
      new_probe_names.push(grset.probe_names[p])
    }
    p = p + 1
  }
  let result = GenomicRatioSet::new(
    new_beta,
    new_mvalues,
    grset.sample_names,
    new_probe_names,
  )
  result.annotation = grset.annotation
  result.pheno_data = grset.pheno_data
  result
}

///|
/// Remove SNP-related probes.
/// Drops probes whose names start with "rs" (the dbSNP identifier prefix).
pub fn minfi_drop_snp_probes(grset : GenomicRatioSet) -> GenomicRatioSet {
  let n_probes = grset.beta.length()
  let new_beta : Array[Array[Double]] = Array::new()
  let new_mvalues : Array[Array[Double]] = Array::new()
  let new_probe_names : Array[String] = Array::new()
  let mut p = 0
  while p < n_probes {
    let name = grset.probe_names[p]
    if !name.has_prefix("rs") {
      new_beta.push(grset.beta[p])
      new_mvalues.push(grset.mvalues[p])
      new_probe_names.push(name)
    }
    p = p + 1
  }
  let result = GenomicRatioSet::new(
    new_beta,
    new_mvalues,
    grset.sample_names,
    new_probe_names,
  )
  result.annotation = grset.annotation
  result.pheno_data = grset.pheno_data
  result
}

///|
/// Map probes to genome positions.
/// Returns a Map from probe name to (chromosome, position). Uses positions
/// encoded in probe names when available (e.g. "cg1234_chr1_5000");
/// otherwise assigns synthetic positions on chr1.
pub fn minfi_map_to_genome(
  grset : GenomicRatioSet,
) -> Map[String, (String, Int)] {
  let result : Map[String, (String, Int)] = Map(
    [],
    capacity=grset.probe_names.length(),
  )
  let mut p = 0
  while p < grset.probe_names.length() {
    let name = grset.probe_names[p]
    let (chr_s, pos_i) = minfi_parse_probe_loc(name, p)
    result[name] = (chr_s, pos_i)
    p = p + 1
  }
  result
}

///|
/// Create an example RGChannelSet with synthetic 450K-like data.
pub fn minfi_create_example_rgset() -> RGChannelSet {
  let sample_names = ["sample1", "sample2", "sample3", "sample4"]
  let probe_names = [
    "cg00000001_chr1_1000", "cg00000002_chr1_2000", "cg00000003_chr1_3000", "cg00000004_chr1_4000",
    "cg00000005_chr1_5000", "cg00000006_chr1_6000", "cg00000007_chr1_7000", "cg00000008_chr1_8000",
    "cg00000009_chr1_9000", "cg00000010_chr1_10000", "cg00000011_chrX_1000", "cg00000012_chrX_2000",
    "cg00000013_chrX_3000", "cg00000014_chrX_4000", "cg00000015_chrX_5000",
  ]
  // Red (unmethylated) intensities: 4 samples x 15 probes.
  let red = [
    [200.0, 180.0, 220.0, 190.0],
    [300.0, 310.0, 290.0, 305.0],
    [150.0, 160.0, 155.0, 165.0],
    [400.0, 410.0, 395.0, 420.0],
    [250.0, 240.0, 260.0, 255.0],
    [350.0, 360.0, 345.0, 365.0],
    [180.0, 190.0, 175.0, 195.0],
    [420.0, 430.0, 410.0, 440.0],
    [270.0, 280.0, 265.0, 285.0],
    [320.0, 330.0, 315.0, 335.0],
    [200.0, 700.0, 210.0, 720.0],
    [220.0, 680.0, 230.0, 690.0],
    [180.0, 650.0, 190.0, 660.0],
    [240.0, 720.0, 250.0, 730.0],
    [210.0, 700.0, 220.0, 710.0],
  ]
  // Green (methylated) intensities.
  let green = [
    [800.0, 820.0, 780.0, 810.0],
    [600.0, 590.0, 610.0, 595.0],
    [900.0, 910.0, 895.0, 905.0],
    [500.0, 490.0, 510.0, 485.0],
    [750.0, 740.0, 760.0, 745.0],
    [550.0, 540.0, 560.0, 545.0],
    [850.0, 860.0, 845.0, 865.0],
    [480.0, 470.0, 490.0, 475.0],
    [730.0, 720.0, 740.0, 725.0],
    [680.0, 670.0, 690.0, 675.0],
    [780.0, 300.0, 770.0, 290.0],
    [760.0, 320.0, 750.0, 310.0],
    [820.0, 350.0, 810.0, 340.0],
    [740.0, 280.0, 730.0, 270.0],
    [790.0, 310.0, 780.0, 300.0],
  ]
  RGChannelSet::new(red, green, sample_names, probe_names)
}

///|
/// Create an example GenomicRatioSet with synthetic beta values.
pub fn minfi_create_example_grset() -> GenomicRatioSet {
  let sample_names = ["sample1", "sample2", "sample3", "sample4"]
  let probe_names = [
    "cg00000001_chr1_1000", "cg00000002_chr1_2000", "cg00000003_chr1_3000", "cg00000004_chr1_4000",
    "cg00000005_chr1_5000", "cg00000006_chr1_6000", "cg00000007_chr1_7000", "cg00000008_chr1_8000",
    "cg00000009_chr1_9000", "cg00000010_chr1_10000", "cg00000011_chrX_1000", "cg00000012_chrX_2000",
    "cg00000013_chrX_3000", "cg00000014_chrX_4000", "cg00000015_chrX_5000",
  ]
  // Beta values: probes x samples. samples 1-2 = group 0, 3-4 = group 1.
  let beta = [
    [0.8, 0.82, 0.5, 0.52],
    [0.7, 0.72, 0.4, 0.42],
    [0.9, 0.91, 0.85, 0.86],
    [0.6, 0.61, 0.55, 0.56],
    [0.85, 0.86, 0.5, 0.51],
    [0.75, 0.76, 0.7, 0.71],
    [0.92, 0.93, 0.6, 0.61],
    [0.5, 0.51, 0.45, 0.46],
    [0.83, 0.84, 0.55, 0.56],
    [0.78, 0.79, 0.77, 0.78],
    [0.45, 0.46, 0.48, 0.47],
    [0.5, 0.51, 0.49, 0.5],
    [0.55, 0.56, 0.54, 0.55],
    [0.42, 0.43, 0.45, 0.44],
    [0.48, 0.49, 0.5, 0.49],
  ]
  let mvalues : Array[Array[Double]] = Array::new()
  let mut p = 0
  while p < beta.length() {
    let row : Array[Double] = Array::new()
    let mut s = 0
    while s < beta[p].length() {
      let b = beta[p][s]
      // M = log2(b / (1 - b)), with clamping.
      let bb = if b < 0.001 { 0.001 } else if b > 0.999 { 0.999 } else { b }
      row.push(minfi_log2(bb / (1.0 - bb)))
      s = s + 1
    }
    mvalues.push(row)
    p = p + 1
  }
  GenomicRatioSet::new(beta, mvalues, sample_names, probe_names)
}

///|
/// Create an example DMPResult for testing DMR analysis.
pub fn minfi_create_example_dmp_result() -> DMPResult {
  let probe_names = [
    "cg00000001_chr1_1000", "cg00000002_chr1_1100", "cg00000003_chr1_1500", "cg00000004_chr1_5000",
    "cg00000005_chr1_5100", "cg00000006_chr1_5200", "cg00000007_chr1_9000", "cg00000008_chr1_9100",
    "cg00000009_chr2_1000", "cg00000010_chr2_1100",
  ]
  let chr = [
    "chr1", "chr1", "chr1", "chr1", "chr1", "chr1", "chr1", "chr1", "chr2", "chr2",
  ]
  let pos = [1000, 1100, 1500, 5000, 5100, 5200, 9000, 9100, 1000, 1100]
  let pvalues = [
    0.0001, 0.0002, 0.5, 0.001, 0.002, 0.4, 0.0005, 0.6, 0.0003, 0.7,
  ]
  let fdrs = [0.001, 0.0015, 0.04, 0.005, 0.006, 0.03, 0.003, 0.65, 0.002, 0.7]
  let beta_diffs = [0.3, 0.25, 0.05, 0.2, 0.18, 0.02, 0.28, 0.01, 0.22, 0.04]
  DMPResult::new(probe_names, chr, pos, pvalues, fdrs, beta_diffs)
}

///|
/// Demo function exercising the major minfi entry points.
/// Not named `main` because `src/` is a library package; the project entry
/// point lives in `cmd/main/`. Call `minfi_demo()` from a main package
/// (or examples/minfi_demo) to run this.
pub fn minfi_demo() -> Unit {
  println("=== minfi module demo ===")
  println("")
  // 1. Build example RGChannelSet and show summary.
  let rgset = minfi_create_example_rgset()
  println(rgset.to_string())
  // 2. NOOB preprocessing.
  let mset_noob = minfi_preprocess_noob(rgset)
  println(mset_noob.to_string())
  // 3. Illumina preprocessing.
  let mset_illumina = minfi_preprocess_illumina(rgset)
  println(mset_illumina.to_string())
  // 4. Quantile preprocessing.
  let mset_quantile = minfi_preprocess_quantile(rgset)
  println(mset_quantile.to_string())
  // 5. Functional normalization (returns GenomicRatioSet directly).
  let grset_funnorm = minfi_preprocess_funnorm(rgset)
  println(grset_funnorm.to_string())
  // 6. Convert MethylSet to GenomicRatioSet.
  let grset = minfi_ratio_convert(mset_noob)
  println(grset.to_string())
  // 7. Get beta and M values.
  let beta = minfi_get_beta(grset)
  let mvalues = minfi_get_mvalues(grset)
  println(
    "Beta matrix: " +
    beta.length().to_string() +
    " probes x " +
    (if beta.length() > 0 { beta[0].length() } else { 0 }).to_string() +
    " samples",
  )
  println(
    "M-value matrix: " +
    mvalues.length().to_string() +
    " probes x " +
    (if mvalues.length() > 0 { mvalues[0].length() } else { 0 }).to_string() +
    " samples",
  )
  // 8. Detection p-values.
  let det_p = minfi_detection_pvalue(rgset)
  println("Detection p-values (first 5): ")
  let mut i = 0
  while i < det_p.length() && i < 5 {
    println("  " + rgset.probe_names[i] + ": " + det_p[i].to_string())
    i = i + 1
  }
  println("")
  // 9. QC report.
  let qc = minfi_get_qc(rgset)
  println(qc.to_string())
  // 10. Sex prediction.
  let sexes = minfi_estimate_sex(grset)
  println("Sex predictions:")
  i = 0
  while i < sexes.length() {
    println("  " + grset.sample_names[i] + ": " + sexes[i])
    i = i + 1
  }
  println("")
  // 11. DMP analysis.
  let groups = [0, 0, 1, 1]
  let dmp_result = minfi_find_dmps(grset, groups)
  println(dmp_result.to_string())
  // 12. DMR analysis from example DMP result (which has clear clusters).
  let example_dmp = minfi_create_example_dmp_result()
  let dmr_result = minfi_find_dmrs(example_dmp, 500)
  println(dmr_result.to_string())
  // 13. Probe filtering.
  let filtered = minfi_filter_probes(grset, 0.01)
  println(
    "After filtering (sd >= 0.01): " +
    filtered.n_probes().to_string() +
    " probes remain",
  )
  // 14. SNP probe removal.
  let no_snp = minfi_drop_snp_probes(grset)
  println(
    "After SNP removal: " + no_snp.n_probes().to_string() + " probes remain",
  )
  // 15. Map probes to genome.
  let loc_map = minfi_map_to_genome(grset)
  let keys = loc_map.keys().collect()
  println("Mapped " + keys.length().to_string() + " probes to genome positions")
  // 16. Math helpers sanity check.
  println("")
  println("Math helpers:")
  println("  minfi_sqrt(2.0) = " + minfi_sqrt(2.0).to_string())
  println("  minfi_ln(2.0) = " + minfi_ln(2.0).to_string())
  println("  minfi_exp(1.0) = " + minfi_exp(1.0).to_string())
  println("  minfi_log2(8.0) = " + minfi_log2(8.0).to_string())
  println(
    "  minfi_compute_beta(800, 200, 100) = " +
    minfi_compute_beta(800.0, 200.0, 100.0).to_string(),
  )
  println(
    "  minfi_compute_mvalue(800, 200) = " +
    minfi_compute_mvalue(800.0, 200.0).to_string(),
  )
}
