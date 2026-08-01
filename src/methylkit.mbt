///|
/// methylKit module - Bisulfite sequencing methylation analysis.
/// Inspired by the Bioconductor methylKit package.
///
/// Provides core utilities for BS-seq methylation data:
/// - MethylCytosine / MethylSample containers for per-base counts
/// - Coverage and context based filtering
/// - Coverage normalization (median scaling)
/// - Differential methylation via Fisher's exact test with BH-FDR
/// - Genome tiling windows and DMR detection
/// - Sample correlation and k-means clustering
/// - BED export helpers and sample data

///|
/// Round a non-negative / negative double to the nearest integer.
fn methylkit_round(x : Double) -> Int {
  if x >= 0.0 {
    (x + 0.5).to_int()
  } else {
    -((-x) + 0.5).to_int()
  }
}

///|
/// Mean of an array of doubles. Returns 0.0 for an empty array.
fn methylkit_mean(arr : Array[Double]) -> Double {
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
/// Copy and sort (ascending) an array of doubles using selection sort.
fn methylkit_sort_doubles(arr : Array[Double]) -> Array[Double] {
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
/// Median of an array of doubles. Returns 0.0 for an empty array.
fn methylkit_median(arr : Array[Double]) -> Double {
  let n = arr.length()
  if n == 0 {
    return 0.0
  }
  let sorted = methylkit_sort_doubles(arr)
  if n % 2 == 1 {
    sorted[n / 2]
  } else {
    (sorted[n / 2 - 1] + sorted[n / 2]) / 2.0
  }
}

///|
/// Logarithm of a binomial coefficient using a precomputed log-factorial table.
fn methylkit_log_choose_tbl(
  n : Int,
  k : Int,
  lf : Array[Double],
) -> Double {
  if k < 0 || k > n {
    return -700.0
  }
  if k == 0 || k == n {
    return 0.0
  }
  lf[n] - lf[k] - lf[n - k]
}

///|
/// A single methylated cytosine record from bisulfite sequencing.
pub struct MethylCytosine {
  chr : String
  position : Int
  strand : String
  context : String
  methylated : Int
  unmethylated : Int
}

///|
/// Construct a MethylCytosine.
pub fn MethylCytosine::new(
  chr : String,
  position : Int,
  strand : String,
  context : String,
  methylated : Int,
  unmethylated : Int,
) -> MethylCytosine {
  MethylCytosine::{ chr, position, strand, context, methylated, unmethylated }
}

///|
/// Chromosome name.
pub fn MethylCytosine::chr(self : MethylCytosine) -> String {
  self.chr
}

///|
/// Genomic position.
pub fn MethylCytosine::position(self : MethylCytosine) -> Int {
  self.position
}

///|
/// Strand ("+" or "-").
pub fn MethylCytosine::strand(self : MethylCytosine) -> String {
  self.strand
}

///|
/// Sequence context ("CG", "CHG" or "CHH").
pub fn MethylCytosine::context(self : MethylCytosine) -> String {
  self.context
}

///|
/// Methylated read count.
pub fn MethylCytosine::methylated(self : MethylCytosine) -> Int {
  self.methylated
}

///|
/// Unmethylated read count.
pub fn MethylCytosine::unmethylated(self : MethylCytosine) -> Int {
  self.unmethylated
}

///|
/// Total coverage (methylated + unmethylated).
pub fn MethylCytosine::coverage(self : MethylCytosine) -> Int {
  self.methylated + self.unmethylated
}

///|
/// Methylation percentage in [0, 100]. Returns 0.0 when coverage is zero.
pub fn MethylCytosine::methylation_pct(self : MethylCytosine) -> Double {
  let cov = self.coverage()
  if cov == 0 {
    return 0.0
  }
  self.methylated.to_double() / cov.to_double() * 100.0
}

///|
/// Human-readable representation.
pub fn MethylCytosine::to_string(self : MethylCytosine) -> String {
  "MethylCytosine(chr=" +
  self.chr +
  ", position=" +
  self.position.to_string() +
  ", strand=" +
  self.strand +
  ", context=" +
  self.context +
  ", methylated=" +
  self.methylated.to_string() +
  ", unmethylated=" +
  self.unmethylated.to_string() +
  ", coverage=" +
  self.coverage().to_string() +
  ", meth_pct=" +
  self.methylation_pct().to_string() +
  ")"
}

///|
/// A sample of bisulfite sequencing methylation calls.
pub struct MethylSample {
  sample_id : String
  coverage : Array[MethylCytosine]
  mut n_cpgs : Int
}

///|
/// Construct an empty MethylSample.
pub fn MethylSample::new(sample_id : String) -> MethylSample {
  MethylSample::{
    sample_id,
    coverage: Array::new(),
    n_cpgs: 0,
  }
}

///|
/// Construct a MethylSample from an existing array of cytosines.
pub fn MethylSample::from_coverage(
  sample_id : String,
  coverage : Array[MethylCytosine],
) -> MethylSample {
  MethylSample::{ sample_id, coverage, n_cpgs: coverage.length() }
}

///|
/// Sample identifier.
pub fn MethylSample::sample_id(self : MethylSample) -> String {
  self.sample_id
}

///|
/// Array of cytosine records.
pub fn MethylSample::coverage(self : MethylSample) -> Array[MethylCytosine] {
  self.coverage
}

///|
/// Number of cytosines stored in the sample.
pub fn MethylSample::n_cpgs(self : MethylSample) -> Int {
  self.n_cpgs
}

///|
/// Append a cytosine record to the sample.
pub fn MethylSample::add_cytosine(
  self : MethylSample,
  mc : MethylCytosine,
) -> Unit {
  self.coverage.push(mc)
  self.n_cpgs = self.n_cpgs + 1
}

///|
/// Look up a cytosine by chromosome and position. Returns None if absent.
pub fn MethylSample::get_cytosine(
  self : MethylSample,
  chr : String,
  position : Int,
) -> MethylCytosine? {
  let mut i = 0
  while i < self.coverage.length() {
    let mc = self.coverage[i]
    if mc.chr == chr && mc.position == position {
      return Some(mc)
    }
    i = i + 1
  }
  None
}

///|
/// All methylation percentages for cytosines in the sample.
pub fn MethylSample::methylation_pcts(self : MethylSample) -> Array[Double] {
  let pcts : Array[Double] = Array::new()
  let mut i = 0
  while i < self.coverage.length() {
    pcts.push(self.coverage[i].methylation_pct())
    i = i + 1
  }
  pcts
}

///|
/// Mean methylation percentage across all cytosines.
pub fn MethylSample::mean_methylation(self : MethylSample) -> Double {
  methylkit_mean(self.methylation_pcts())
}

///|
/// Return a new sample keeping only cytosines within [min_cov, max_cov].
pub fn MethylSample::filter_by_coverage(
  self : MethylSample,
  min_cov : Int,
  max_cov : Int,
) -> MethylSample {
  let kept : Array[MethylCytosine] = Array::new()
  let mut i = 0
  while i < self.coverage.length() {
    let mc = self.coverage[i]
    let cov = mc.coverage()
    if cov >= min_cov && cov <= max_cov {
      kept.push(mc)
    }
    i = i + 1
  }
  MethylSample::from_coverage(self.sample_id, kept)
}

///|
/// A single differential methylation test result.
pub struct MethylDiffResult {
  chr : String
  position : Int
  p_value : Double
  q_value : Double
  meth_diff : Double
  significant : Bool
}

///|
/// Chromosome of the tested position.
pub fn MethylDiffResult::chr(self : MethylDiffResult) -> String {
  self.chr
}

///|
/// Genomic position of the tested position.
pub fn MethylDiffResult::position(self : MethylDiffResult) -> Int {
  self.position
}

///|
/// Raw p-value from the differential test.
pub fn MethylDiffResult::p_value(self : MethylDiffResult) -> Double {
  self.p_value
}

///|
/// BH-adjusted q-value.
pub fn MethylDiffResult::q_value(self : MethylDiffResult) -> Double {
  self.q_value
}

///|
/// Methylation difference (treatment pct - control pct).
pub fn MethylDiffResult::meth_diff(self : MethylDiffResult) -> Double {
  self.meth_diff
}

///|
/// Whether this position is a significant DMP.
pub fn MethylDiffResult::significant(self : MethylDiffResult) -> Bool {
  self.significant
}

///|
/// Human-readable representation.
pub fn MethylDiffResult::to_string(self : MethylDiffResult) -> String {
  "MethylDiffResult(chr=" +
  self.chr +
  ", position=" +
  self.position.to_string() +
  ", p_value=" +
  self.p_value.to_string() +
  ", q_value=" +
  self.q_value.to_string() +
  ", meth_diff=" +
  self.meth_diff.to_string() +
  ", significant=" +
  self.significant.to_string() +
  ")"
}

///|
/// Configuration parameters for methylKit analyses.
pub struct MethylKitConfig {
  min_coverage : Int
  max_coverage : Int
  context_filter : String
  min_perc_samples : Double
  min_diff : Double
  q_threshold : Double
}

///|
/// Default configuration:
/// min_coverage=10, max_coverage=999, context_filter="CG",
/// min_perc_samples=0.6, min_diff=25.0, q_threshold=0.05.
pub fn MethylKitConfig::new() -> MethylKitConfig {
  MethylKitConfig::{
    min_coverage: 10,
    max_coverage: 999,
    context_filter: "CG",
    min_perc_samples: 0.6,
    min_diff: 25.0,
    q_threshold: 0.05,
  }
}

///|
/// Construct a configuration with explicit parameters.
pub fn MethylKitConfig::with_params(
  min_coverage : Int,
  max_coverage : Int,
  context : String,
  min_perc : Double,
  min_diff : Double,
  q_threshold : Double,
) -> MethylKitConfig {
  MethylKitConfig::{
    min_coverage,
    max_coverage,
    context_filter: context,
    min_perc_samples: min_perc,
    min_diff,
    q_threshold,
  }
}

///|
pub fn MethylKitConfig::min_coverage(self : MethylKitConfig) -> Int {
  self.min_coverage
}

///|
pub fn MethylKitConfig::max_coverage(self : MethylKitConfig) -> Int {
  self.max_coverage
}

///|
pub fn MethylKitConfig::context_filter(self : MethylKitConfig) -> String {
  self.context_filter
}

///|
pub fn MethylKitConfig::min_perc_samples(self : MethylKitConfig) -> Double {
  self.min_perc_samples
}

///|
pub fn MethylKitConfig::min_diff(self : MethylKitConfig) -> Double {
  self.min_diff
}

///|
pub fn MethylKitConfig::q_threshold(self : MethylKitConfig) -> Double {
  self.q_threshold
}

///|
/// Filter a single sample by coverage and sequence context.
pub fn methylkit_filter_sample(
  sample : MethylSample,
  config : MethylKitConfig,
) -> MethylSample {
  let kept : Array[MethylCytosine] = Array::new()
  let mut i = 0
  while i < sample.coverage.length() {
    let mc = sample.coverage[i]
    let cov = mc.coverage()
    if cov >= config.min_coverage &&
      cov <= config.max_coverage &&
      mc.context == config.context_filter {
      kept.push(mc)
    }
    i = i + 1
  }
  MethylSample::from_coverage(sample.sample_id, kept)
}

///|
/// Filter all samples by coverage/context and require each retained CpG to be
/// present in at least `min_perc_samples` fraction of the (already filtered)
/// samples.
pub fn methylkit_filter_samples(
  samples : Array[MethylSample],
  config : MethylKitConfig,
) -> Array[MethylSample] {
  let filtered : Array[MethylSample] = Array::new()
  let mut i = 0
  while i < samples.length() {
    filtered.push(methylkit_filter_sample(samples[i], config))
    i = i + 1
  }
  let n_samples = filtered.length().to_double()
  if n_samples == 0.0 {
    return filtered
  }
  // Count how many samples contain each (chr, position) key.
  let counts : Map[String, Int] = Map([], capacity=256)
  let mut s = 0
  while s < filtered.length() {
    let seen : Map[String, Bool] = Map([], capacity=128)
    let mut k = 0
    while k < filtered[s].coverage.length() {
      let mc = filtered[s].coverage[k]
      let key = mc.chr + ":" + mc.position.to_string()
      seen[key] = true
      k = k + 1
    }
    let keys = seen.keys().collect()
    let mut ki = 0
    while ki < keys.length() {
      let key = keys[ki]
      counts[key] = counts.get(key).unwrap_or(0) + 1
      ki = ki + 1
    }
    s = s + 1
  }
  let threshold = config.min_perc_samples * n_samples
  let keep : Map[String, Bool] = Map([], capacity=256)
  let keys = counts.keys().collect()
  let mut ki = 0
  while ki < keys.length() {
    let key = keys[ki]
    if counts.get(key).unwrap_or(0).to_double() >= threshold {
      keep[key] = true
    }
    ki = ki + 1
  }
  let result : Array[MethylSample] = Array::new()
  let mut s = 0
  while s < filtered.length() {
    let new_cov : Array[MethylCytosine] = Array::new()
    let mut k = 0
    while k < filtered[s].coverage.length() {
      let mc = filtered[s].coverage[k]
      let key = mc.chr + ":" + mc.position.to_string()
      match keep.get(key) {
        Some(_) => new_cov.push(mc)
        None => ()
      }
      k = k + 1
    }
    result.push(MethylSample::from_coverage(filtered[s].sample_id, new_cov))
    s = s + 1
  }
  result
}

///|
/// Coverage normalization: scale each cytosine's counts so that its total
/// coverage equals the median per-cytosine coverage of the sample.
pub fn methylkit_normalize(sample : MethylSample) -> MethylSample {
  let covs : Array[Double] = Array::new()
  let mut i = 0
  while i < sample.coverage.length() {
    covs.push(sample.coverage[i].coverage().to_double())
    i = i + 1
  }
  let med = methylkit_median(covs)
  let new_cov : Array[MethylCytosine] = Array::new()
  let mut i = 0
  while i < sample.coverage.length() {
    let mc = sample.coverage[i]
    let cov = mc.coverage().to_double()
    if cov <= 0.0 {
      new_cov.push(
        MethylCytosine::new(mc.chr, mc.position, mc.strand, mc.context, 0, 0),
      )
    } else {
      let scale = med / cov
      let nm = methylkit_round(mc.methylated.to_double() * scale)
      let nu = methylkit_round(mc.unmethylated.to_double() * scale)
      new_cov.push(
        MethylCytosine::new(mc.chr, mc.position, mc.strand, mc.context, nm, nu),
      )
    }
    i = i + 1
  }
  MethylSample::from_coverage(sample.sample_id, new_cov)
}

///|
/// Apply per-sample normalization to every sample.
pub fn methylkit_normalize_samples(
  samples : Array[MethylSample],
) -> Array[MethylSample] {
  let result : Array[MethylSample] = Array::new()
  let mut i = 0
  while i < samples.length() {
    result.push(methylkit_normalize(samples[i]))
    i = i + 1
  }
  result
}

///|
/// Two-sided Fisher's exact test p-value for the 2x2 table
/// [[mc1.methylated, mc1.unmethylated], [mc2.methylated, mc2.unmethylated]].
/// Uses the hypergeometric distribution evaluated via log-factorials.
pub fn methylkit_fisher_test(
  mc1 : MethylCytosine,
  mc2 : MethylCytosine,
) -> Double {
  let a = mc1.methylated
  let b = mc1.unmethylated
  let c = mc2.methylated
  let d = mc2.unmethylated
  let r1 = a + b
  let r2 = c + d
  let c1 = a + c
  let c2 = b + d
  let n = a + b + c + d
  if n == 0 {
    return 1.0
  }
  if r1 == 0 || r2 == 0 || c1 == 0 || c2 == 0 {
    return 1.0
  }
  // Precompute log-factorial table up to n.
  let log_fact : Array[Double] = Array::make(n + 1, 0.0)
  for i in 2..=n {
    log_fact[i] = log_fact[i - 1] + @math.ln(i.to_double())
  }
  let log_pmf = fn(aa : Int) -> Double {
    let bb = r1 - aa
    let cc = c1 - aa
    let dd = r2 - cc
    if bb < 0 || cc < 0 || dd < 0 {
      return -700.0
    }
    methylkit_log_choose_tbl(r1, aa, log_fact) +
    methylkit_log_choose_tbl(r2, cc, log_fact) -
    methylkit_log_choose_tbl(n, c1, log_fact)
  }
  let log_p_obs = log_pmf(a)
  let p_obs = @math.exp(log_p_obs)
  // Range of valid table entries for the first cell.
  let lo = if c1 - r2 > 0 { c1 - r2 } else { 0 }
  let hi = if r1 < c1 { r1 } else { c1 }
  let mut p_sum = 0.0
  for aa in lo..=hi {
    let p = @math.exp(log_pmf(aa))
    if p <= p_obs + 1.0e-12 {
      p_sum = p_sum + p
    }
  }
  if p_sum > 1.0 {
    1.0
  } else if p_sum < 0.0 {
    0.0
  } else {
    p_sum
  }
}

///|
/// Benjamini-Hochberg FDR correction of a vector of p-values.
pub fn methylkit_bh_correct(p_values : Array[Double]) -> Array[Double] {
  let n = p_values.length()
  if n == 0 {
    return Array::new()
  }
  // Indexed sort by ascending p-value (selection sort).
  let indices : Array[Int] = Array::new()
  let sorted_p : Array[Double] = Array::new()
  let mut i = 0
  while i < n {
    indices.push(i)
    sorted_p.push(p_values[i])
    i = i + 1
  }
  let mut i = 0
  while i < n - 1 {
    let mut min_idx = i
    let mut j = i + 1
    while j < n {
      if sorted_p[j] < sorted_p[min_idx] {
        min_idx = j
      }
      j = j + 1
    }
    if min_idx != i {
      let tp = sorted_p[i]
      sorted_p[i] = sorted_p[min_idx]
      sorted_p[min_idx] = tp
      let ti = indices[i]
      indices[i] = indices[min_idx]
      indices[min_idx] = ti
    }
    i = i + 1
  }
  // Walk from largest to smallest rank, enforcing monotonicity.
  let adj_sorted : Array[Double] = Array::make(n, 1.0)
  let mut min_adj = 1.0
  let mut i = n - 1
  while i >= 0 {
    let rank = i + 1
    let adj = sorted_p[i] * n.to_double() / rank.to_double()
    if adj < min_adj {
      min_adj = adj
    }
    let clamped = if min_adj > 1.0 { 1.0 } else { min_adj }
    adj_sorted[i] = clamped
    i = i - 1
  }
  // Map back to original order.
  let result : Array[Double] = Array::make(n, 1.0)
  let mut i = 0
  while i < n {
    result[indices[i]] = adj_sorted[i]
    i = i + 1
  }
  result
}

///|
/// Test every CpG shared between two samples for differential methylation.
/// Applies BH-FDR correction and marks significant positions where
/// q_value < q_threshold AND abs(meth_diff) >= min_diff.
pub fn methylkit_diff_meth(
  sample1 : MethylSample,
  sample2 : MethylSample,
  config : MethylKitConfig,
) -> Array[MethylDiffResult] {
  let p_values : Array[Double] = Array::new()
  let info : Array[(String, Int, Double)] = Array::new()
  let mut i = 0
  while i < sample1.coverage.length() {
    let mc1 = sample1.coverage[i]
    match sample2.get_cytosine(mc1.chr, mc1.position) {
      Some(mc2) => {
        let p = methylkit_fisher_test(mc1, mc2)
        let diff = mc1.methylation_pct() - mc2.methylation_pct()
        p_values.push(p)
        info.push((mc1.chr, mc1.position, diff))
      }
      None => ()
    }
    i = i + 1
  }
  let q_values = methylkit_bh_correct(p_values)
  let results : Array[MethylDiffResult] = Array::new()
  let mut i = 0
  while i < info.length() {
    let (chr, pos, diff) = info[i]
    let q = q_values[i]
    let sig = q < config.q_threshold && diff.abs() >= config.min_diff
    results.push(
      MethylDiffResult::{
        chr,
        position: pos,
        p_value: p_values[i],
        q_value: q,
        meth_diff: diff,
        significant: sig,
      },
    )
    i = i + 1
  }
  results
}

///|
/// Tile the genome into windows and test each window for differential
/// methylation between the first two samples (control vs treatment).
/// Requires at least two samples; otherwise returns an empty array.
pub fn methylkit_tiling(
  samples : Array[MethylSample],
  window_size : Int,
  step_size : Int,
) -> Array[MethylDiffResult] {
  if samples.length() < 2 || window_size <= 0 || step_size <= 0 {
    return Array::new()
  }
  let s1 = samples[0]
  let s2 = samples[1]
  // Determine per-chromosome min/max positions across both samples.
  let chr_min : Map[String, Int] = Map([], capacity=16)
  let chr_max : Map[String, Int] = Map([], capacity=16)
  let collect = fn(sample : MethylSample) -> Unit {
    let mut k = 0
    while k < sample.coverage.length() {
      let mc = sample.coverage[k]
      let new_mn = match chr_min.get(mc.chr) {
        Some(v) => if mc.position < v { mc.position } else { v }
        None => mc.position
      }
      let new_mx = match chr_max.get(mc.chr) {
        Some(v) => if mc.position > v { mc.position } else { v }
        None => mc.position
      }
      chr_min[mc.chr] = new_mn
      chr_max[mc.chr] = new_mx
      k = k + 1
    }
  }
  collect(s1)
  collect(s2)
  let chrs = chr_min.keys().collect()
  let p_values : Array[Double] = Array::new()
  let info : Array[(String, Int, Double)] = Array::new()
  let mut ci = 0
  while ci < chrs.length() {
    let chr = chrs[ci]
    let mn = chr_min.get(chr).unwrap_or(0)
    let mx = chr_max.get(chr).unwrap_or(0)
    let mut start = mn
    while start <= mx {
      let end = start + window_size
      let mut m1 = 0
      let mut u1 = 0
      let mut m2 = 0
      let mut u2 = 0
      let mut k = 0
      while k < s1.coverage.length() {
        let mc = s1.coverage[k]
        if mc.chr == chr && mc.position >= start && mc.position < end {
          m1 = m1 + mc.methylated
          u1 = u1 + mc.unmethylated
        }
        k = k + 1
      }
      let mut k = 0
      while k < s2.coverage.length() {
        let mc = s2.coverage[k]
        if mc.chr == chr && mc.position >= start && mc.position < end {
          m2 = m2 + mc.methylated
          u2 = u2 + mc.unmethylated
        }
        k = k + 1
      }
      if m1 + u1 > 0 && m2 + u2 > 0 {
        let mc1 = MethylCytosine::new(chr, start, "+", "CG", m1, u1)
        let mc2 = MethylCytosine::new(chr, start, "+", "CG", m2, u2)
        let p = methylkit_fisher_test(mc1, mc2)
        let diff = mc1.methylation_pct() - mc2.methylation_pct()
        p_values.push(p)
        info.push((chr, start, diff))
      }
      start = start + step_size
    }
    ci = ci + 1
  }
  let q_values = methylkit_bh_correct(p_values)
  let results : Array[MethylDiffResult] = Array::new()
  let mut i = 0
  while i < info.length() {
    let (chr, pos, diff) = info[i]
    results.push(
      MethylDiffResult::{
        chr,
        position: pos,
        p_value: p_values[i],
        q_value: q_values[i],
        meth_diff: diff,
        significant: q_values[i] < 0.05 && diff.abs() >= 25.0,
      },
    )
    i = i + 1
  }
  results
}

///|
/// Identify differentially methylated regions: groups of consecutive
/// significant CpGs (within 1000 bp of each other on the same chromosome).
/// Returns (chr, start, end) tuples.
pub fn methylkit_get_dmr(
  diff_results : Array[MethylDiffResult],
) -> Array[(String, Int, Int)] {
  let n = diff_results.length()
  if n == 0 {
    return Array::new()
  }
  // Copy and sort by (chr, position) using selection sort.
  let sorted = diff_results.copy()
  let mut i = 0
  while i < n - 1 {
    let mut min_idx = i
    let mut j = i + 1
    while j < n {
      if sorted[j].chr < sorted[min_idx].chr ||
        (sorted[j].chr == sorted[min_idx].chr &&
          sorted[j].position < sorted[min_idx].position) {
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
  let dmrs : Array[(String, Int, Int)] = Array::new()
  let mut i = 0
  while i < n {
    if sorted[i].significant {
      let chr = sorted[i].chr
      let start_pos = sorted[i].position
      let mut end_pos = sorted[i].position
      let mut j = i + 1
      while j < n &&
        sorted[j].significant &&
        sorted[j].chr == chr &&
        sorted[j].position - end_pos <= 1000 {
        end_pos = sorted[j].position
        j = j + 1
      }
      dmrs.push((chr, start_pos, end_pos))
      i = j
    } else {
      i = i + 1
    }
  }
  dmrs
}

///|
/// Pearson correlation of methylation percentages at CpGs shared between
/// two samples. Returns 0.0 when fewer than two shared CpGs are available.
pub fn methylkit_sample_correlation(
  sample1 : MethylSample,
  sample2 : MethylSample,
) -> Double {
  let xs : Array[Double] = Array::new()
  let ys : Array[Double] = Array::new()
  let mut i = 0
  while i < sample1.coverage.length() {
    let mc1 = sample1.coverage[i]
    match sample2.get_cytosine(mc1.chr, mc1.position) {
      Some(mc2) => {
        xs.push(mc1.methylation_pct())
        ys.push(mc2.methylation_pct())
      }
      None => ()
    }
    i = i + 1
  }
  let m = xs.length()
  if m < 2 {
    return 0.0
  }
  let nd = m.to_double()
  let mut sx = 0.0
  let mut sy = 0.0
  let mut sxy = 0.0
  let mut sx2 = 0.0
  let mut sy2 = 0.0
  let mut i = 0
  while i < m {
    let x = xs[i]
    let y = ys[i]
    sx = sx + x
    sy = sy + y
    sxy = sxy + x * y
    sx2 = sx2 + x * x
    sy2 = sy2 + y * y
    i = i + 1
  }
  let num = nd * sxy - sx * sy
  let den = ((nd * sx2 - sx * sx) * (nd * sy2 - sy * sy)).sqrt()
  if den == 0.0 {
    0.0
  } else {
    num / den
  }
}

///|
/// NxN Pearson correlation matrix for N samples.
pub fn methylkit_correlation_matrix(
  samples : Array[MethylSample],
) -> Array[Array[Double]] {
  let n = samples.length()
  let mat : Array[Array[Double]] = Array::new()
  let mut i = 0
  while i < n {
    let row : Array[Double] = Array::new()
    let mut j = 0
    while j < n {
      row.push(methylkit_sample_correlation(samples[i], samples[j]))
      j = j + 1
    }
    mat.push(row)
    i = i + 1
  }
  mat
}

///|
/// Simple k-means clustering of samples based on their mean methylation
/// profile. Returns a vector of cluster labels (one per sample).
pub fn methylkit_cluster_samples(
  samples : Array[MethylSample],
  n_clusters : Int,
) -> Array[Int] {
  let n = samples.length()
  if n == 0 {
    return Array::new()
  }
  let k = n_clusters
  if k <= 0 {
    return Array::new()
  }
  if k == 1 {
    return Array::make(n, 0)
  }
  if k >= n {
    let labels : Array[Int] = Array::new()
    let mut i = 0
    while i < n {
      labels.push(i)
      i = i + 1
    }
    return labels
  }
  // 1D feature: mean methylation.
  let features : Array[Array[Double]] = Array::new()
  let mut i = 0
  while i < n {
    features.push([samples[i].mean_methylation()])
    i = i + 1
  }
  // Initialize centroids from evenly spaced samples.
  let centroids : Array[Array[Double]] = Array::new()
  let mut i = 0
  while i < k {
    let idx = i * (n - 1) / (k - 1)
    centroids.push(features[idx].copy())
    i = i + 1
  }
  let labels = Array::make(n, 0)
  let mut changed = true
  let mut iter = 0
  while changed && iter < 100 {
    changed = false
    // Assignment step.
    let mut i = 0
    while i < n {
      let mut best = 0
      let mut best_dist = 1000000000000.0
      let mut c = 0
      while c < k {
        let d0 = features[i][0] - centroids[c][0]
        let dist = d0 * d0
        if dist < best_dist {
          best_dist = dist
          best = c
        }
        c = c + 1
      }
      if labels[i] != best {
        changed = true
        labels[i] = best
      }
      i = i + 1
    }
    // Update step.
    let sums = Array::make(k, 0.0)
    let counts = Array::make(k, 0)
    let mut i = 0
    while i < n {
      let c = labels[i]
      sums[c] = sums[c] + features[i][0]
      counts[c] = counts[c] + 1
      i = i + 1
    }
    let mut c = 0
    while c < k {
      if counts[c] > 0 {
        centroids[c][0] = sums[c] / counts[c].to_double()
      }
      c = c + 1
    }
    iter = iter + 1
  }
  labels
}

///|
/// Human-readable summary of a sample.
pub fn methylkit_sample_summary(sample : MethylSample) -> String {
  let covs : Array[Double] = Array::new()
  let mut i = 0
  while i < sample.coverage.length() {
    covs.push(sample.coverage[i].coverage().to_double())
    i = i + 1
  }
  let mean_cov = methylkit_mean(covs)
  let mut s = "MethylSample Summary\n"
  s = s + "====================\n"
  s = s + "Sample ID: " + sample.sample_id + "\n"
  s = s + "Number of CpGs: " + sample.n_cpgs().to_string() + "\n"
  s = s + "Mean methylation: " + sample.mean_methylation().to_string() + "%\n"
  s = s + "Mean coverage: " + mean_cov.to_string() + "\n"
  s
}

///|
/// Human-readable summary of a differential methylation result set.
pub fn methylkit_diff_summary(results : Array[MethylDiffResult]) -> String {
  let n = results.length()
  let mut n_sig = 0
  let mut n_hyper = 0
  let mut n_hypo = 0
  let mut i = 0
  while i < n {
    if results[i].significant {
      n_sig = n_sig + 1
      if results[i].meth_diff > 0.0 {
        n_hyper = n_hyper + 1
      } else {
        n_hypo = n_hypo + 1
      }
    }
    i = i + 1
  }
  let mut s = "Differential Methylation Summary\n"
  s = s + "================================\n"
  s = s + "Total CpGs tested: " + n.to_string() + "\n"
  s = s + "Significant DMPs: " + n_sig.to_string() + "\n"
  s = s + "  Hyper-methylated: " + n_hyper.to_string() + "\n"
  s = s + "  Hypo-methylated: " + n_hypo.to_string() + "\n"
  s
}

///|
/// Convert a cytosine to a BED line (chrom, start, end, meth_pct).
pub fn methylkit_cpg_to_bed(mc : MethylCytosine) -> String {
  let start = mc.position
  let end = mc.position + 1
  mc.chr +
  "\t" +
  start.to_string() +
  "\t" +
  end.to_string() +
  "\t" +
  mc.methylation_pct().to_string()
}

///|
/// Convert a DMR tuple (chr, start, end) to a BED line.
pub fn methylkit_dmr_to_bed(dmr : (String, Int, Int)) -> String {
  let (chr, start, end) = dmr
  chr + "\t" + start.to_string() + "\t" + end.to_string()
}

///|
/// Build two example samples (control and treatment), each with 20 CpGs.
/// The treatment sample has 5 hyper-methylated and 5 hypo-methylated CpGs
/// relative to the control.
pub fn methylkit_sample_data() -> (MethylSample, MethylSample) {
  let control_cov : Array[MethylCytosine] = Array::new()
  let treatment_cov : Array[MethylCytosine] = Array::new()
  let mut i = 0
  while i < 20 {
    let pos = 1000 * (i + 1)
    // Control: moderate ~50% methylation, 20x coverage.
    control_cov.push(MethylCytosine::new("chr1", pos, "+", "CG", 10, 10))
    let (tm, tu) = if i < 5 {
      // Hyper-methylated in treatment (~90%).
      (18, 2)
    } else if i < 10 {
      // Hypo-methylated in treatment (~10%).
      (2, 18)
    } else {
      // Similar to control (~50%).
      (10, 10)
    }
    treatment_cov.push(MethylCytosine::new("chr1", pos, "+", "CG", tm, tu))
    i = i + 1
  }
  (
    MethylSample::from_coverage("control", control_cov),
    MethylSample::from_coverage("treatment", treatment_cov),
  )
}
