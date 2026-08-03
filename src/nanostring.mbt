///|
/// NanoString nCounter module - Digital molecular barcode counting analysis.
/// Inspired by the Bioconductor NanoStringNCTools package.
///
/// NanoString nCounter measures gene expression directly without amplification
/// using color-coded molecular barcodes. This module provides:
/// - Sample and data containers with control/housekeeping/endogenous splits
/// - Positive control, negative control, housekeeping and content normalization
/// - Quality control (imaging, linearity, limit of quantification)
/// - Differential expression analysis with t-tests and BH-FDR

// ============================================================================
// Private statistical helpers
// ============================================================================

///|
/// Arithmetic mean of an array of doubles. Returns 0.0 when empty.
fn ns_mean(arr : Array[Double]) -> Double {
  let n = arr.length()
  if n == 0 {
    return 0.0
  }
  let mut s = 0.0
  for v in arr {
    s = s + v
  }
  s / n.to_double()
}

///|
/// Sample variance (n-1 denominator) given the precomputed mean.
fn ns_variance(arr : Array[Double], mean : Double) -> Double {
  let n = arr.length()
  if n < 2 {
    return 0.0
  }
  let mut s = 0.0
  for v in arr {
    let d = v - mean
    s = s + d * d
  }
  s / (n - 1).to_double()
}

///|
/// Sample standard deviation given the precomputed mean.
fn ns_std(arr : Array[Double], mean : Double) -> Double {
  ns_variance(arr, mean).sqrt()
}

///|
/// Lanczos approximation for ln(gamma(x)).
fn ns_lgamma(x : Double) -> Double {
  let cof = [
    76.18009172947146,
    -86.50532032941677,
    24.01409824083091,
    -1.231739572450155,
    0.001208650973866179,
    -0.000005395239384953,
  ]
  let mut x1 = x
  let mut tmp = x1 + 5.5
  tmp = (x1 + 0.5) * @math.ln(tmp) - tmp
  let mut ser = 1.000000000190015
  let mut j = 0
  while j < cof.length() {
    x1 = x1 + 1.0
    ser = ser + cof[j] / x1
    j = j + 1
  }
  tmp + @math.ln(2.5066282746310005 * ser)
}

///|
/// Continued fraction for the incomplete beta function.
fn ns_betacf(x : Double, a : Double, b : Double) -> Double {
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
    h = h * d * c
    if (h - 1.0).abs() < eps {
      return h
    }
    m = m + 1
  }
  h
}

///|
/// Regularized incomplete beta function I_x(a, b).
fn ns_pbeta(x : Double, a : Double, b : Double) -> Double {
  if x <= 0.0 {
    return 0.0
  }
  if x >= 1.0 {
    return 1.0
  }
  let lbeta_val = ns_lgamma(a) + ns_lgamma(b) - ns_lgamma(a + b)
  let mut bt = @math.ln(x) * a + @math.ln(1.0 - x) * b - lbeta_val
  bt = @math.exp(bt)
  let mut betaf = 0.0
  if x < (a + 1.0) / (a + b + 2.0) {
    betaf = bt * ns_betacf(x, a, b) / a
  } else {
    betaf = 1.0 - bt * ns_betacf(1.0 - x, b, a) / b
  }
  if betaf < 0.0 {
    0.0
  } else if betaf > 1.0 {
    1.0
  } else {
    betaf
  }
}

///|
/// Student's t-distribution CDF.
fn ns_t_cdf(t : Double, df : Double) -> Double {
  if df <= 0.0 {
    return 0.5
  }
  let x = df / (df + t * t)
  let p = 0.5 * ns_pbeta(x, df / 2.0, 0.5)
  if t >= 0.0 {
    1.0 - p
  } else {
    p
  }
}

///|
/// Simple linear regression of y on x. Returns (slope, intercept, r_squared).
fn ns_linear_regression(
  x : Array[Double],
  y : Array[Double],
) -> (Double, Double, Double) {
  let n = x.length()
  if n < 2 || n != y.length() {
    return (0.0, 0.0, 0.0)
  }
  let nd = n.to_double()
  let mut sum_x = 0.0
  let mut sum_y = 0.0
  let mut sum_xy = 0.0
  let mut sum_xx = 0.0
  let mut sum_yy = 0.0
  for i in 0..<n {
    sum_x = sum_x + x[i]
    sum_y = sum_y + y[i]
    sum_xy = sum_xy + x[i] * y[i]
    sum_xx = sum_xx + x[i] * x[i]
    sum_yy = sum_yy + y[i] * y[i]
  }
  let denom = nd * sum_xx - sum_x * sum_x
  if denom == 0.0 {
    return (0.0, 0.0, 0.0)
  }
  let slope = (nd * sum_xy - sum_x * sum_y) / denom
  let intercept = (sum_y - slope * sum_x) / nd
  let ss_tot = sum_yy - sum_y * sum_y / nd
  let ss_res = sum_yy - intercept * sum_y - slope * sum_xy
  let r_sq = if ss_tot == 0.0 {
    0.0
  } else {
    1.0 - ss_res / ss_tot
  }
  let r_clamped = if r_sq < 0.0 {
    0.0
  } else if r_sq > 1.0 {
    1.0
  } else {
    r_sq
  }
  (slope, intercept, r_clamped)
}

// ============================================================================
// Data structures
// ============================================================================

///|
/// A single NanoString sample with its control and endogenous counts.
pub struct NsSample {
  sample_id : String
  raw_counts : Array[Int]
  positive_controls : Array[Int]
  negative_controls : Array[Int]
  housekeeping : Array[Int]
  endogenous : Array[Int]
}

///|
/// Construct an NsSample, splitting raw_counts into control and endogenous
/// subsets using the provided index arrays.
pub fn NsSample::new(
  sample_id : String,
  raw_counts : Array[Int],
  positive_control_indices : Array[Int],
  negative_control_indices : Array[Int],
  housekeeping_indices : Array[Int],
  endogenous_indices : Array[Int],
) -> NsSample {
  let positive_controls : Array[Int] = Array::new()
  for idx in positive_control_indices {
    if idx >= 0 && idx < raw_counts.length() {
      positive_controls.push(raw_counts[idx])
    }
  }
  let negative_controls : Array[Int] = Array::new()
  for idx in negative_control_indices {
    if idx >= 0 && idx < raw_counts.length() {
      negative_controls.push(raw_counts[idx])
    }
  }
  let housekeeping : Array[Int] = Array::new()
  for idx in housekeeping_indices {
    if idx >= 0 && idx < raw_counts.length() {
      housekeeping.push(raw_counts[idx])
    }
  }
  let endogenous : Array[Int] = Array::new()
  for idx in endogenous_indices {
    if idx >= 0 && idx < raw_counts.length() {
      endogenous.push(raw_counts[idx])
    }
  }
  NsSample::{
    sample_id,
    raw_counts,
    positive_controls,
    negative_controls,
    housekeeping,
    endogenous,
  }
}

///|
/// Sample identifier.
pub fn NsSample::sample_id(self : NsSample) -> String {
  self.sample_id
}

///|
/// Raw count vector for the whole sample.
pub fn NsSample::raw_counts(self : NsSample) -> Array[Int] {
  self.raw_counts
}

///|
/// Positive control counts.
pub fn NsSample::positive_controls(self : NsSample) -> Array[Int] {
  self.positive_controls
}

///|
/// Negative control counts.
pub fn NsSample::negative_controls(self : NsSample) -> Array[Int] {
  self.negative_controls
}

///|
/// Housekeeping gene counts.
pub fn NsSample::housekeeping(self : NsSample) -> Array[Int] {
  self.housekeeping
}

///|
/// Endogenous gene counts.
pub fn NsSample::endogenous(self : NsSample) -> Array[Int] {
  self.endogenous
}

///|
/// Container for a NanoString experiment with multiple samples.
pub struct NsNanostringData {
  samples : Array[NsSample]
  gene_names : Array[String]
  mut n_genes : Int
  mut n_samples : Int
}

///|
/// Construct an empty NsNanostringData with the given endogenous gene names.
pub fn NsNanostringData::new(
  gene_names : Array[String],
) -> NsNanostringData {
  NsNanostringData::{
    samples: Array::new(),
    gene_names,
    n_genes: gene_names.length(),
    n_samples: 0,
  }
}

///|
/// Append a sample to the experiment.
pub fn NsNanostringData::add_sample(
  self : NsNanostringData,
  sample : NsSample,
) -> Unit {
  self.samples.push(sample)
  self.n_samples = self.samples.length()
  self.n_genes = self.gene_names.length()
}

///|
/// Number of samples in the experiment.
pub fn NsNanostringData::n_samples(self : NsNanostringData) -> Int {
  self.n_samples
}

///|
/// Number of endogenous genes.
pub fn NsNanostringData::n_genes(self : NsNanostringData) -> Int {
  self.n_genes
}

///|
/// Gene names.
pub fn NsNanostringData::gene_names(self : NsNanostringData) -> Array[String] {
  self.gene_names
}

///|
/// Sample array.
pub fn NsNanostringData::samples(self : NsNanostringData) -> Array[NsSample] {
  self.samples
}

///|
/// Result of normalizing a single sample.
pub struct NsNormResult {
  sample_id : String
  norm_factors : Array[Double]
  normalized_counts : Array[Double]
  method : String
}

///|
/// Construct an NsNormResult.
pub fn NsNormResult::new(
  sample_id : String,
  norm_factors : Array[Double],
  normalized_counts : Array[Double],
  method : String,
) -> NsNormResult {
  NsNormResult::{
    sample_id,
    norm_factors,
    normalized_counts,
    method,
  }
}

///|
/// Sample identifier.
pub fn NsNormResult::sample_id(self : NsNormResult) -> String {
  self.sample_id
}

///|
/// Normalization factors applied to the sample.
pub fn NsNormResult::norm_factors(self : NsNormResult) -> Array[Double] {
  self.norm_factors
}

///|
/// Normalized endogenous counts.
pub fn NsNormResult::normalized_counts(self : NsNormResult) -> Array[Double] {
  self.normalized_counts
}

///|
/// Normalization method name.
pub fn NsNormResult::method(self : NsNormResult) -> String {
  self.method
}

///|
/// Differential expression result for a single gene.
pub struct NsDeResult {
  gene : String
  log_fc : Double
  p_value : Double
  adj_p_value : Double
  significant : Bool
}

///|
/// Construct an NsDeResult.
pub fn NsDeResult::new(
  gene : String,
  log_fc : Double,
  p_value : Double,
  adj_p_value : Double,
  significant : Bool,
) -> NsDeResult {
  NsDeResult::{ gene, log_fc, p_value, adj_p_value, significant }
}

///|
/// Gene name.
pub fn NsDeResult::gene(self : NsDeResult) -> String {
  self.gene
}

///|
/// Log2 fold change (treatment - control).
pub fn NsDeResult::log_fc(self : NsDeResult) -> Double {
  self.log_fc
}

///|
/// Raw p-value.
pub fn NsDeResult::p_value(self : NsDeResult) -> Double {
  self.p_value
}

///|
/// BH-adjusted p-value.
pub fn NsDeResult::adj_p_value(self : NsDeResult) -> Double {
  self.adj_p_value
}

///|
/// Significance flag (adj p < 0.05).
pub fn NsDeResult::significant(self : NsDeResult) -> Bool {
  self.significant
}

// ============================================================================
// Utility functions
// ============================================================================

///|
/// Geometric mean of an array (only positive values contribute).
/// Returns 0.0 when no positive values are present.
pub fn ns_geometric_mean(values : Array[Double]) -> Double {
  let mut sum_log = 0.0
  let mut n_valid = 0
  for v in values {
    if v > 0.0 {
      sum_log = sum_log + @math.ln(v)
      n_valid = n_valid + 1
    }
  }
  if n_valid > 0 {
    @math.exp(sum_log / n_valid.to_double())
  } else {
    0.0
  }
}

///|
/// Log2 transform counts with a 0.5 pseudocount: log2(x + 0.5).
pub fn ns_log2_counts(counts : Array[Double]) -> Array[Double] {
  let result : Array[Double] = Array::new()
  for c in counts {
    result.push(@math.log2(c + 0.5))
  }
  result
}

// ============================================================================
// Normalization
// ============================================================================

///|
/// Positive-control normalization.
/// Each sample is scaled so that the geometric mean of its positive controls
/// matches the geometric mean of those means across all samples.
pub fn ns_positive_control_norm(
  data : NsNanostringData,
) -> Array[NsNormResult] {
  let samples = data.samples()
  let n_genes = data.n_genes()
  let sample_geo_means : Array[Double] = Array::new()
  for sample in samples {
    let pos_doubles : Array[Double] = Array::new()
    for c in sample.positive_controls() {
      pos_doubles.push(c.to_double())
    }
    sample_geo_means.push(ns_geometric_mean(pos_doubles))
  }
  let target = ns_geometric_mean(sample_geo_means)
  let results : Array[NsNormResult] = Array::new()
  for i in 0..<samples.length() {
    let sample = samples[i]
    let gm = sample_geo_means[i]
    let factor = if gm > 0.0 {
      target / gm
    } else {
      1.0
    }
    let norm_counts : Array[Double] = Array::new()
    let endo = sample.endogenous()
    for g in 0..<n_genes {
      if g < endo.length() {
        norm_counts.push(endo[g].to_double() * factor)
      } else {
        norm_counts.push(0.0)
      }
    }
    results.push(
      NsNormResult::new(sample.sample_id(), [factor], norm_counts, "positive_control"),
    )
  }
  results
}

///|
/// Negative-control normalization.
/// A background threshold (mean + 2*SD of negative controls) is subtracted
/// from every count; values falling at or below zero are floored to 0.5.
pub fn ns_negative_control_norm(
  data : NsNanostringData,
) -> Array[NsNormResult] {
  let samples = data.samples()
  let n_genes = data.n_genes()
  let results : Array[NsNormResult] = Array::new()
  for sample in samples {
    let neg_doubles : Array[Double] = Array::new()
    for c in sample.negative_controls() {
      neg_doubles.push(c.to_double())
    }
    let mean = ns_mean(neg_doubles)
    let sd = ns_std(neg_doubles, mean)
    let threshold = mean + 2.0 * sd
    let norm_counts : Array[Double] = Array::new()
    let endo = sample.endogenous()
    for g in 0..<n_genes {
      if g < endo.length() {
        let v = endo[g].to_double() - threshold
        norm_counts.push(if v < 0.5 { 0.5 } else { v })
      } else {
        norm_counts.push(0.5)
      }
    }
    results.push(
      NsNormResult::new(
        sample.sample_id(),
        [1.0],
        norm_counts,
        "negative_control",
      ),
    )
  }
  results
}

///|
/// Housekeeping-gene normalization.
/// Each sample is scaled so that the geometric mean of its housekeeping gene
/// counts matches the geometric mean of those means across all samples.
pub fn ns_housekeeping_norm(
  data : NsNanostringData,
) -> Array[NsNormResult] {
  let samples = data.samples()
  let n_genes = data.n_genes()
  let sample_geo_means : Array[Double] = Array::new()
  for sample in samples {
    let hk_doubles : Array[Double] = Array::new()
    for c in sample.housekeeping() {
      hk_doubles.push(c.to_double())
    }
    sample_geo_means.push(ns_geometric_mean(hk_doubles))
  }
  let target = ns_geometric_mean(sample_geo_means)
  let results : Array[NsNormResult] = Array::new()
  for i in 0..<samples.length() {
    let sample = samples[i]
    let gm = sample_geo_means[i]
    let factor = if gm > 0.0 {
      target / gm
    } else {
      1.0
    }
    let norm_counts : Array[Double] = Array::new()
    let endo = sample.endogenous()
    for g in 0..<n_genes {
      if g < endo.length() {
        norm_counts.push(endo[g].to_double() * factor)
      } else {
        norm_counts.push(0.0)
      }
    }
    results.push(
      NsNormResult::new(
        sample.sample_id(),
        [factor],
        norm_counts,
        "housekeeping",
      ),
    )
  }
  results
}

///|
/// Content normalization: positive control normalization followed by
/// housekeeping normalization on the positive-control-scaled counts.
pub fn ns_content_norm(data : NsNanostringData) -> Array[NsNormResult] {
  let pos_results = ns_positive_control_norm(data)
  let samples = data.samples()
  let n_genes = data.n_genes()
  // Build housekeeping geometric means from the positive-control-scaled
  // housekeeping counts.
  let hk_geo_means : Array[Double] = Array::new()
  for i in 0..<samples.length() {
    let sample = samples[i]
    let pos_factor = if pos_results[i].norm_factors().length() > 0 {
      pos_results[i].norm_factors()[0]
    } else {
      1.0
    }
    let hk_doubles : Array[Double] = Array::new()
    for c in sample.housekeeping() {
      hk_doubles.push(c.to_double() * pos_factor)
    }
    hk_geo_means.push(ns_geometric_mean(hk_doubles))
  }
  let target = ns_geometric_mean(hk_geo_means)
  let results : Array[NsNormResult] = Array::new()
  for i in 0..<samples.length() {
    let pos_factor = if pos_results[i].norm_factors().length() > 0 {
      pos_results[i].norm_factors()[0]
    } else {
      1.0
    }
    let gm = hk_geo_means[i]
    let hk_factor = if gm > 0.0 {
      target / gm
    } else {
      1.0
    }
    let norm_counts : Array[Double] = Array::new()
    let pos_counts = pos_results[i].normalized_counts()
    for g in 0..<n_genes {
      if g < pos_counts.length() {
        norm_counts.push(pos_counts[g] * hk_factor)
      } else {
        norm_counts.push(0.0)
      }
    }
    results.push(
      NsNormResult::new(
        samples[i].sample_id(),
        [pos_factor, hk_factor],
        norm_counts,
        "content",
      ),
    )
  }
  results
}

// ============================================================================
// Quality control
// ============================================================================

///|
/// Imaging QC: flag samples whose average positive-control binding density
/// (used here as a proxy) is at least `threshold`.
pub fn ns_imaging_qc(data : NsNanostringData, threshold : Int) -> Array[Bool] {
  let results : Array[Bool] = Array::new()
  for sample in data.samples() {
    let pos = sample.positive_controls()
    let mut total = 0
    for c in pos {
      total = total + c
    }
    let avg = if pos.length() > 0 {
      total / pos.length()
    } else {
      0
    }
    results.push(avg >= threshold)
  }
  results
}

///|
/// Positive control linearity QC.
/// Fits a linear model of log2(positive control counts) against log2 of the
/// known dilution-series concentrations. Returns (r_squared, slope).
pub fn ns_pos_linearity_qc(
  data : NsNanostringData,
) -> (Double, Double) {
  // Standard NanoString positive control dilution series (fM).
  let concentrations = [128.0, 32.0, 8.0, 2.0, 0.5, 0.125]
  let log_conc : Array[Double] = Array::new()
  for c in concentrations {
    log_conc.push(@math.log2(c))
  }
  let xs : Array[Double] = Array::new()
  let ys : Array[Double] = Array::new()
  for sample in data.samples() {
    let pos = sample.positive_controls()
    for i in 0..<pos.length() {
      if i < log_conc.length() {
        let v = pos[i].to_double()
        if v > 0.0 {
          xs.push(log_conc[i])
          ys.push(@math.log2(v))
        }
      }
    }
  }
  let (slope, _intercept, r_sq) = ns_linear_regression(xs, ys)
  (r_sq, slope)
}

///|
/// Limit of quantification QC.
/// For each sample, count endogenous genes whose raw count is below `loq`.
pub fn ns_loq_qc(data : NsNanostringData, loq : Int) -> Array[Int] {
  let results : Array[Int] = Array::new()
  for sample in data.samples() {
    let mut count = 0
    for c in sample.endogenous() {
      if c < loq {
        count = count + 1
      }
    }
    results.push(count)
  }
  results
}

///|
/// Build a textual QC summary for a NanoString experiment.
pub fn ns_qc_summary(data : NsNanostringData) -> String {
  let samples = data.samples()
  let imaging = ns_imaging_qc(data, 100)
  let (r_sq, slope) = ns_pos_linearity_qc(data)
  let loq = ns_loq_qc(data, 50)
  let buf = StringBuilder::new()
  buf.write_string("NanoString QC Summary\n")
  buf.write_string("======================\n")
  buf.write_string("Samples: ")
  buf.write_string(samples.length().to_string())
  buf.write_string("\nPositive control linearity R^2: ")
  buf.write_string(r_sq.to_string())
  buf.write_string("\nPositive control slope: ")
  buf.write_string(slope.to_string())
  buf.write_string("\nImaging QC (threshold=100):\n")
  for i in 0..<samples.length() {
    buf.write_string("  ")
    buf.write_string(samples[i].sample_id())
    buf.write_string(": pass=")
    buf.write_string(imaging[i].to_string())
    buf.write_string(", below LOQ=")
    buf.write_string(loq[i].to_string())
    buf.write_string("\n")
  }
  buf.to_string()
}

// ============================================================================
// Differential expression
// ============================================================================

///|
/// Welch's two-sample t-test. Returns (t_statistic, p_value).
pub fn ns_t_test(
  g1 : Array[Double],
  g2 : Array[Double],
) -> (Double, Double) {
  let n1 = g1.length()
  let n2 = g2.length()
  if n1 < 2 || n2 < 2 {
    return (0.0, 1.0)
  }
  let m1 = ns_mean(g1)
  let m2 = ns_mean(g2)
  let v1 = ns_variance(g1, m1)
  let v2 = ns_variance(g2, m2)
  let n1d = n1.to_double()
  let n2d = n2.to_double()
  let se = (v1 / n1d + v2 / n2d).sqrt()
  if se <= 0.0 {
    return (0.0, 1.0)
  }
  let t = (m1 - m2) / se
  let term1 = v1 / n1d
  let term2 = v2 / n2d
  let num = (term1 + term2) * (term1 + term2)
  let den = term1 * term1 / (n1 - 1).to_double() +
    term2 * term2 / (n2 - 1).to_double()
  let df = if den > 0.0 {
    num / den
  } else {
    (n1 + n2 - 2).to_double()
  }
  let p = 2.0 * (1.0 - ns_t_cdf(t.abs(), df))
  let pc = if p < 0.0 {
    0.0
  } else if p > 1.0 {
    1.0
  } else {
    p
  }
  (t, pc)
}

///|
/// Benjamini-Hochberg FDR correction.
pub fn ns_bh_correct(p_values : Array[Double]) -> Array[Double] {
  let n = p_values.length()
  if n == 0 {
    return Array::new()
  }
  let indices = Array::make(n, 0)
  for i in 0..<n {
    indices[i] = i
  }
  let mut swapped = true
  let mut k = n - 1
  while k > 0 && swapped {
    swapped = false
    let mut i = 0
    while i < k {
      if p_values[indices[i]] > p_values[indices[i + 1]] {
        let tmp = indices[i]
        indices[i] = indices[i + 1]
        indices[i + 1] = tmp
        swapped = true
      }
      i = i + 1
    }
    k = k - 1
  }
  let adj = Array::make(n, 1.0)
  let mut min_adj = 1.0
  let mut i = n - 1
  while i >= 0 {
    let rank = (i + 1).to_double()
    let val = p_values[indices[i]] * n.to_double() / rank
    if val < min_adj {
      min_adj = val
    }
    adj[indices[i]] = if min_adj > 1.0 {
      1.0
    } else {
      min_adj
    }
    i = i - 1
  }
  adj
}

///|
/// Differential expression analysis from normalized counts.
/// `norm_results` holds one NsNormResult per sample, `group` flags treatment
/// (true) vs control (false). Counts are log2-transformed and tested per gene
/// with a Welch t-test, then BH-FDR corrected.
pub fn ns_differential_expression(
  norm_results : Array[NsNormResult],
  group : Array[Bool],
  gene_names : Array[String],
) -> Array[NsDeResult] {
  let n_samples = norm_results.length()
  let n_genes = gene_names.length()
  let p_values = Array::make(n_genes, 1.0)
  let log_fcs = Array::make(n_genes, 0.0)
  let t_stats = Array::make(n_genes, 0.0)
  for g in 0..<n_genes {
    let treat : Array[Double] = Array::new()
    let ctrl : Array[Double] = Array::new()
    for s in 0..<n_samples {
      let counts = norm_results[s].normalized_counts()
      if g < counts.length() {
        let v = @math.log2(counts[g] + 0.5)
        if s < group.length() && group[s] {
          treat.push(v)
        } else {
          ctrl.push(v)
        }
      }
    }
    let (t_stat, p_val) = ns_t_test(treat, ctrl)
    t_stats[g] = t_stat
    p_values[g] = p_val
    log_fcs[g] = ns_mean(treat) - ns_mean(ctrl)
  }
  let adj_p = ns_bh_correct(p_values)
  let results : Array[NsDeResult] = Array::new()
  for g in 0..<n_genes {
    results.push(
      NsDeResult::new(
        gene_names[g],
        log_fcs[g],
        p_values[g],
        adj_p[g],
        adj_p[g] < 0.05,
      ),
    )
  }
  results
}

// ============================================================================
// Summaries
// ============================================================================

///|
/// Build a textual summary of a NanoString experiment.
pub fn ns_summary(data : NsNanostringData) -> String {
  let buf = StringBuilder::new()
  buf.write_string("NanoString Data Summary\n")
  buf.write_string("=======================\n")
  buf.write_string("Samples: ")
  buf.write_string(data.n_samples().to_string())
  buf.write_string("\nGenes: ")
  buf.write_string(data.n_genes().to_string())
  buf.write_string("\n")
  for sample in data.samples() {
    let mut total = 0
    for c in sample.raw_counts() {
      total = total + c
    }
    buf.write_string("  ")
    buf.write_string(sample.sample_id())
    buf.write_string(": total counts = ")
    buf.write_string(total.to_string())
    buf.write_string("\n")
  }
  buf.to_string()
}

///|
/// Build a textual summary of a differential expression analysis.
pub fn ns_de_summary(results : Array[NsDeResult]) -> String {
  let n = results.length()
  let mut n_up = 0
  let mut n_down = 0
  let mut n_sig = 0
  for r in results {
    if r.significant() {
      n_sig = n_sig + 1
      if r.log_fc() > 0.0 {
        n_up = n_up + 1
      } else {
        n_down = n_down + 1
      }
    }
  }
  let buf = StringBuilder::new()
  buf.write_string("Differential Expression Summary\n")
  buf.write_string("================================\n")
  buf.write_string("Genes tested: ")
  buf.write_string(n.to_string())
  buf.write_string("\nSignificant (adj p < 0.05): ")
  buf.write_string(n_sig.to_string())
  buf.write_string("\n  Up-regulated: ")
  buf.write_string(n_up.to_string())
  buf.write_string("\n  Down-regulated: ")
  buf.write_string(n_down.to_string())
  buf.write_string("\n")
  buf.to_string()
}

// ============================================================================
// Sample data
// ============================================================================

///|
/// A sample NanoString experiment: 20 endogenous genes x 6 samples
/// (3 control + 3 treatment). Each sample's raw_counts vector has 32 entries:
/// 20 endogenous genes followed by 6 positive controls and 6 negative controls.
/// Genes 0-4 are up-regulated and genes 5-9 are down-regulated in treatment.
/// The first positive control indices double as housekeeping gene positions.
pub fn ns_sample_data() -> NsNanostringData {
  let gene_names = [
    "GENE01",
    "GENE02",
    "GENE03",
    "GENE04",
    "GENE05",
    "GENE06",
    "GENE07",
    "GENE08",
    "GENE09",
    "GENE10",
    "GENE11",
    "GENE12",
    "GENE13",
    "GENE14",
    "GENE15",
    "GENE16",
    "GENE17",
    "GENE18",
    "GENE19",
    "GENE20",
  ]
  let data = NsNanostringData::new(gene_names)
  let pos_indices = [20, 21, 22, 23, 24, 25]
  let neg_indices = [26, 27, 28, 29, 30, 31]
  let hk_indices = [0, 5, 10, 15]
  let endo_indices = [
    0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19,
  ]
  // Positive control dilution series counts (proportional to 128,32,8,2,0.5,0.125)
  // with mild per-sample variation.
  let pos_control_profiles = [
    [6400, 1600, 420, 110, 35, 12],
    [6350, 1580, 415, 108, 33, 11],
    [6420, 1610, 425, 112, 36, 12],
    [6380, 1595, 418, 109, 34, 11],
    [6410, 1605, 422, 111, 35, 12],
    [6395, 1590, 416, 107, 32, 10],
  ]
  // Negative control counts (low background with noise).
  let neg_control_profiles = [
    [14, 11, 16, 9, 12, 13],
    [15, 10, 14, 11, 13, 12],
    [13, 12, 15, 10, 11, 14],
    [16, 13, 17, 12, 14, 15],
    [14, 11, 15, 10, 12, 13],
    [15, 12, 16, 11, 13, 12],
  ]
  let sample_ids = ["C1", "C2", "C3", "T1", "T2", "T3"]
  for s in 0..<6 {
    let raw : Array[Int] = Array::new()
    // Endogenous genes 0-4: up-regulated in treatment.
    for g in 0..<5 {
      if s < 3 {
        raw.push(40 + g * 2 + s)
      } else {
        raw.push(320 + g * 5 + (s - 3) * 4)
      }
    }
    // Endogenous genes 5-9: down-regulated in treatment.
    for g in 0..<5 {
      if s < 3 {
        raw.push(280 + g * 4 + s)
      } else {
        raw.push(35 + g * 2 + (s - 3))
      }
    }
    // Endogenous genes 10-19: stable across groups.
    for g in 0..<10 {
      raw.push(120 + g * 6 + s)
    }
    // Positive controls.
    for v in pos_control_profiles[s] {
      raw.push(v)
    }
    // Negative controls.
    for v in neg_control_profiles[s] {
      raw.push(v)
    }
    data.add_sample(NsSample::new(
      sample_ids[s],
      raw,
      pos_indices,
      neg_indices,
      hk_indices,
      endo_indices,
    ))
  }
  data
}
