///|
/// csaw module, inspired by Bioconductor csaw package.
/// Provides window-based ChIP-seq differential binding analysis.
/// Supports sliding window counting, normalization, filtering,
/// and negative binomial GLM-based differential testing.

/// Normalization method for window counts
pub enum CswNormMethod {
  /// TMM (trimmed mean of M-values)
  TMM
  /// Loess normalization based on input
  Loess
  /// Simple library size normalization
  Library
  /// No normalization
  None
}

/// Filtering status for a window
pub enum CswFilterStatus {
  /// Window passes filter
  Pass
  /// Window fails abundance filter
  LowAbundance
  /// Window fails spatial filter (merged with neighbor)
  Merged
  /// Window excluded
  Excluded
} derive(Eq)

/// A single genomic window
pub struct CswWindow {
  /// Window ID
  id : String
  /// Chromosome
  chrom : String
  /// Start position (0-based)
  start : Int
  /// End position (0-based, exclusive)
  end : Int
  /// Counts per sample
  counts : Array[Double]
  /// Filter status
  filter_status : CswFilterStatus
  /// Filter abundance (log2 CPB)
  filter_abundance : Double
}

/// Collection of windows across samples
pub struct CswDataSet {
  /// Windows
  windows : Array[CswWindow]
  /// Sample names
  samples : Array[String]
  /// Library sizes
  lib_sizes : Array[Double]
  /// Normalization factors
  norm_factors : Array[Double]
  /// Filtered windows (after filtering)
  filtered_windows : Array[CswWindow]
  /// Binned window indices (for region finding)
  bin_indices : Array[Int]
}

/// Normalization result
pub struct CswNormResult {
  /// Normalization factors per sample
  norm_factors : Array[Double]
  /// Normalization method used
  norm_method : CswNormMethod
  /// Reference sample index (for TMM)
  ref_sample : Int
  /// Log-ratios used for TMM
  log_ratios : Array[Array[Double]]
}

/// Differential binding result
pub struct CswResult {
  /// Window indices
  window_indices : Array[Int]
  /// Log2 fold changes
  log_fc : Array[Double]
  /// Standard errors
  std_err : Array[Double]
  /// Wald statistics
  wald_stat : Array[Double]
  /// P-values
  p_values : Array[Double]
  /// Adjusted p-values (BH)
  adj_p_values : Array[Double]
  /// Significant flags
  significant : Array[Bool]
}

/// Region result (merged adjacent significant windows)
pub struct CswRegion {
  /// Region ID
  id : String
  /// Chromosome
  chrom : String
  /// Start position
  start : Int
  /// End position
  end : Int
  /// Combined p-value (Fisher's method)
  combined_p_value : Double
  /// Maximum logFC
  max_log_fc : Double
  /// Number of merged windows
  n_windows : Int
  /// Window indices
  window_indices : Array[Int]
}

// ===== Window methods =====

/// Create a new window
pub fn CswWindow::new(id : String, chrom : String, start : Int, end : Int, counts : Array[Double]) -> CswWindow {
  CswWindow::{
    id,
    chrom,
    start,
    end,
    counts,
    filter_status: CswFilterStatus::Pass,
    filter_abundance: 0.0
  }
}

/// Get window width
pub fn CswWindow::csw_width(self : CswWindow) -> Int {
  self.end - self.start
}

/// Get count for a sample index
pub fn CswWindow::csw_count(self : CswWindow, sample_idx : Int) -> Double {
  self.counts[sample_idx]
}

/// Get mean count across samples
pub fn CswWindow::csw_mean_count(self : CswWindow) -> Double {
  let n = self.counts.length()
  if n == 0 { return 0.0 }
  let mut sum = 0.0
  let mut i = 0
  while i < n {
    sum = sum + self.counts[i]
    i = i + 1
  }
  sum / n.to_double()
}

// ===== DataSet methods =====

/// Create a new CswDataSet
pub fn CswDataSet::new(windows : Array[CswWindow], samples : Array[String], lib_sizes : Array[Double]) -> CswDataSet {
  let n_samples = samples.length()
  let norm_factors : Array[Double] = Array::new()
  let mut i = 0
  while i < n_samples {
    norm_factors.push(1.0)
    i = i + 1
  }
  CswDataSet::{
    windows,
    samples,
    lib_sizes,
    norm_factors,
    filtered_windows: windows.copy(),
    bin_indices: Array::new()
  }
}

/// Get number of windows
pub fn CswDataSet::csw_n_windows(self : CswDataSet) -> Int {
  self.windows.length()
}

/// Get number of samples
pub fn CswDataSet::csw_n_samples(self : CswDataSet) -> Int {
  self.samples.length()
}

/// Get a window by index
pub fn CswDataSet::csw_get_window(self : CswDataSet, idx : Int) -> CswWindow {
  self.windows[idx]
}

/// Get filtered windows
pub fn CswDataSet::csw_get_filtered(self : CswDataSet) -> Array[CswWindow] {
  self.filtered_windows.copy()
}

// ===== Sliding window counting =====

///|
/// Generate sliding windows across a chromosome.
pub fn csw_make_windows(chrom : String, start : Int, end : Int, width : Int, shift : Int) -> Array[CswWindow] {
  let windows : Array[CswWindow] = Array::new()
  let mut pos = start
  while pos + width <= end {
    let w = CswWindow::new(
      "win_" + chrom + "_" + pos.to_string(),
      chrom,
      pos,
      pos + width,
      Array::new()
    )
    windows.push(w)
    pos = pos + shift
  }
  windows
}

///|
/// Assign counts to windows from read positions.
pub fn csw_assign_counts(windows : Array[CswWindow], read_positions : Array[Int]) -> Array[CswWindow] {
  let result : Array[CswWindow] = Array::new()
  let mut i = 0
  while i < windows.length() {
    let w = windows[i]
    let mut count = 0.0
    let mut j = 0
    while j < read_positions.length() {
      if read_positions[j] >= w.start && read_positions[j] < w.end {
        count = count + 1.0
      }
      j = j + 1
    }
    w.counts.push(count)
    result.push(w)
    i = i + 1
  }
  result
}

// ===== Normalization =====

///|
/// Compute library size normalization factors.
pub fn csw_norm_library(lib_sizes : Array[Double]) -> Array[Double] {
  let n = lib_sizes.length()
  let factors : Array[Double] = Array::new()
  let mean = csw_mean_double(lib_sizes)
  let mut i = 0
  while i < n {
    let f = if lib_sizes[i] > 0.0 { mean / lib_sizes[i] } else { 1.0 }
    factors.push(f)
    i = i + 1
  }
  factors
}

///|
/// Compute TMM normalization factors.
/// Simplified version: finds a reference sample and computes M-values.
pub fn csw_norm_tmm(data : CswDataSet) -> CswNormResult {
  let n_samples = data.samples.length()
  let n_windows = data.windows.length()
  let lib_sizes = data.lib_sizes
  let factors : Array[Double] = Array::new()
  let mut i = 0

  // Find reference sample (closest to mean)
  let mean_lib = csw_mean_double(lib_sizes)
  let mut ref_idx = 0
  let mut min_diff = 1.0e30
  let mut s = 0
  while s < n_samples {
    let diff_val = if lib_sizes[s] > 0.0 { (lib_sizes[s] - mean_lib).abs() } else { 1.0e30 }
    if diff_val < min_diff {
      min_diff = diff_val
      ref_idx = s
    }
    s = s + 1
  }

  // Compute log-ratios for normalization
  let log_ratios : Array[Array[Double]] = Array::new()
  while i < n_samples {
    let sample_log_ratios : Array[Double] = Array::new()
    if i == ref_idx {
      factors.push(1.0)
      log_ratios.push(sample_log_ratios)
      i = i + 1
      continue
    }
    let mut k = 0
    let mut m_vals_sum = 0.0
    let mut m_count = 0
    while k < n_windows {
      let ref_count = data.windows[k].counts[ref_idx]
      let obs_count = data.windows[k].counts[i]
      let ref_norm = if lib_sizes[ref_idx] > 0.0 { ref_count / lib_sizes[ref_idx] } else { 0.0 }
      let obs_norm = if lib_sizes[i] > 0.0 { obs_count / lib_sizes[i] } else { 0.0 }
      if ref_norm > 0.0 && obs_norm > 0.0 {
        let m_val = @math.ln(obs_norm / ref_norm) / @math.ln(2.0)
        sample_log_ratios.push(m_val)
        m_vals_sum = m_vals_sum + m_val
        m_count = m_count + 1
      }
      k = k + 1
    }
    let trimmed_m = if m_count > 0 { m_vals_sum / m_count.to_double() } else { 0.0 }
    let factor = @math.exp(-trimmed_m * @math.ln(2.0))
    factors.push(factor)
    log_ratios.push(sample_log_ratios)
    i = i + 1
  }
  CswNormResult::{
    norm_factors: factors,
    norm_method: CswNormMethod::TMM,
    ref_sample: ref_idx,
    log_ratios
  }
}

///|
/// Apply normalization factors to window counts.
pub fn csw_apply_normalization(data : CswDataSet, norm_factors : Array[Double]) -> CswDataSet {
  let n_samples = data.samples.length()
  let n_windows = data.windows.length()
  let new_windows : Array[CswWindow] = Array::new()
  let mut i = 0
  while i < n_windows {
    let w = data.windows[i]
    let new_counts : Array[Double] = Array::new()
    let mut j = 0
    while j < n_samples {
      let lib_adj = if data.lib_sizes[j] > 0.0 { 1.0 / data.lib_sizes[j] } else { 1.0 }
      let norm_count = w.counts[j] * lib_adj * norm_factors[j]
      new_counts.push(norm_count)
      j = j + 1
    }
    new_windows.push(CswWindow::{
      id: w.id,
      chrom: w.chrom,
      start: w.start,
      end: w.end,
      counts: new_counts,
      filter_status: w.filter_status,
      filter_abundance: w.filter_abundance
    })
    i = i + 1
  }
  CswDataSet::{
    windows: new_windows,
    samples: data.samples,
    lib_sizes: data.lib_sizes,
    norm_factors,
    filtered_windows: new_windows.copy(),
    bin_indices: Array::new()
  }
}

// ===== Filtering =====

///|
/// Filter windows by abundance (log2 CPB threshold).
pub fn csw_filter_abundance(data : CswDataSet, log2_cpb_threshold : Double) -> CswDataSet {
  let n_windows = data.windows.length()
  let filtered : Array[CswWindow] = Array::new()
  let mut i = 0
  while i < n_windows {
    let w = data.windows[i]
    let total_cpb = csw_total_cpb(w.counts, data.lib_sizes)
    let log2_cpb = if total_cpb > 0.0 { @math.ln(total_cpb) / @math.ln(2.0) } else { -1.0e9 }
    let status = if log2_cpb >= log2_cpb_threshold { CswFilterStatus::Pass } else { CswFilterStatus::LowAbundance }
    let fw = CswWindow::{
      id: w.id,
      chrom: w.chrom,
      start: w.start,
      end: w.end,
      counts: w.counts,
      filter_status: status,
      filter_abundance: log2_cpb
    }
    filtered.push(fw)
    i = i + 1
  }
  CswDataSet::{
    windows: data.windows,
    samples: data.samples,
    lib_sizes: data.lib_sizes,
    norm_factors: data.norm_factors,
    filtered_windows: filtered,
    bin_indices: data.bin_indices
  }
}

///|
/// Filter windows by spatial merging (merge adjacent windows with similar counts).
pub fn csw_filter_merged(data : CswDataSet, tol : Double) -> CswDataSet {
  let filtered : Array[CswWindow] = data.filtered_windows.copy()
  let n = filtered.length()
  let mut i = 0
  while i < n {
    if filtered[i].filter_status == CswFilterStatus::Pass {
      let mut j = i + 1
      while j < n &&
            filtered[j].filter_status == CswFilterStatus::Pass &&
            filtered[i].chrom == filtered[j].chrom &&
            csw_adjacent_window(filtered[i], filtered[j]) {
        let mean_i = filtered[i].csw_mean_count()
        let mean_j = filtered[j].csw_mean_count()
        let ratio = if mean_i > 0.0 { mean_j / mean_i } else { 1.0 }
        let log_ratio_val = (@math.ln(ratio) / @math.ln(2.0)).abs()
        if log_ratio_val < tol {
          filtered[j] = CswWindow::{
            id: filtered[j].id,
            chrom: filtered[j].chrom,
            start: filtered[j].start,
            end: filtered[j].end,
            counts: filtered[j].counts,
            filter_status: CswFilterStatus::Merged,
            filter_abundance: filtered[j].filter_abundance
          }
        }
        j = j + 1
      }
    }
    i = i + 1
  }
  CswDataSet::{
    windows: data.windows,
    samples: data.samples,
    lib_sizes: data.lib_sizes,
    norm_factors: data.norm_factors,
    filtered_windows: filtered,
    bin_indices: data.bin_indices
  }
}

///|
/// Check if two windows are adjacent (on same chromosome, within gap).
fn csw_adjacent_window(a : CswWindow, b : CswWindow) -> Bool {
  let gap = b.start - a.end
  gap >= 0 && gap <= 500
}

// ===== Differential binding =====

///|
/// Perform differential binding testing using negative binomial GLM.
/// Simplified Wald test for two-group comparison.
pub fn csw_test_differential(
  data : CswDataSet,
  group1_indices : Array[Int],
  group2_indices : Array[Int],
  dispersion : Double
) -> CswResult {
  let filtered = data.filtered_windows
  let n_windows = filtered.length()
  let log_fc : Array[Double] = Array::new()
  let std_err : Array[Double] = Array::new()
  let wald_stat : Array[Double] = Array::new()
  let p_values : Array[Double] = Array::new()
  let window_indices : Array[Int] = Array::new()

  let mut i = 0
  while i < n_windows {
    let w = filtered[i]
    if w.filter_status != CswFilterStatus::Pass {
      i = i + 1
      continue
    }
    let (mean1, variance1) = csw_mean_var_group(w.counts, group1_indices)
    let (mean2, variance2) = csw_mean_var_group(w.counts, group2_indices)
    let logfc = if mean1 > 0.0 && mean2 > 0.0 {
      @math.ln(mean2 / mean1) / @math.ln(2.0)
    } else {
      0.0
    }
    let disp = if dispersion > 0.0 { dispersion } else { 0.1 }
    let se = if variance1 > 0.0 || variance2 > 0.0 {
      (disp * (mean1 / csw_mean_int(group1_indices) + mean2 / csw_mean_int(group2_indices))).sqrt()
    } else {
      1.0
    }
    let wald = if se > 0.0 { logfc / se } else { 0.0 }
    let pval = csw_norm_p_value(wald)
    log_fc.push(logfc)
    std_err.push(se)
    wald_stat.push(wald)
    p_values.push(pval)
    window_indices.push(i)
    i = i + 1
  }

  // BH correction
  let adj_p = csw_bh_correction(p_values)
  let significant : Array[Bool] = Array::new()
  let mut j = 0
  while j < adj_p.length() {
    significant.push(adj_p[j] < 0.05)
    j = j + 1
  }
  CswResult::{
    window_indices,
    log_fc,
    std_err,
    wald_stat,
    p_values,
    adj_p_values: adj_p,
    significant
  }
}

///|
/// Find regions by merging adjacent significant windows.
pub fn csw_find_regions(data : CswDataSet, result : CswResult, merge_tol : Int) -> Array[CswRegion] {
  let regions : Array[CswRegion] = Array::new()
  let n = result.window_indices.length()
  let mut i = 0
  while i < n {
    if !result.significant[i] {
      i = i + 1
      continue
    }
    let start_idx = i
    let mut end_idx = i
    while end_idx + 1 < n &&
          result.significant[end_idx + 1] &&
          csw_can_merge_windows(data, result.window_indices[end_idx], result.window_indices[end_idx + 1], merge_tol) {
      end_idx = end_idx + 1
    }
    let region = csw_build_region(data, result, start_idx, end_idx)
    regions.push(region)
    i = end_idx + 1
  }
  regions
}

///|
/// Check if two windows can be merged into a region.
fn csw_can_merge_windows(data : CswDataSet, idx1 : Int, idx2 : Int, tol : Int) -> Bool {
  let w1 = data.filtered_windows[idx1]
  let w2 = data.filtered_windows[idx2]
  if w1.chrom != w2.chrom { return false }
  let gap = w2.start - w1.end
  gap >= 0 && gap <= tol
}

///|
/// Build a region from merged windows.
fn csw_build_region(data : CswDataSet, result : CswResult, start : Int, end : Int) -> CswRegion {
  let win_indices : Array[Int] = Array::new()
  let mut i = start
  let mut combined_chi2 = 0.0
  let mut max_lfc = -1.0e30
  let mut chrom = ""
  let mut rstart = 1.0e30
  let mut rend = 0.0
  while i <= end {
    let wi = result.window_indices[i]
    win_indices.push(wi)
    combined_chi2 = combined_chi2 + result.wald_stat[i] * result.wald_stat[i]
    if result.log_fc[i] > max_lfc {
      max_lfc = result.log_fc[i]
    }
    let w = data.filtered_windows[wi]
    chrom = w.chrom
    if w.start.to_double() < rstart { rstart = w.start.to_double() }
    if w.end.to_double() > rend { rend = w.end.to_double() }
    i = i + 1
  }
  let n_windows = end - start + 1
  let combined_p = if combined_chi2 > 0.0 {
    csw_chi2_sf(combined_chi2, (n_windows * 2).to_double())
  } else {
    1.0
  }
  CswRegion::{
    id: "region_" + start.to_string() + "_" + end.to_string(),
    chrom,
    start: rstart.to_int(),
    end: rend.to_int(),
    combined_p_value: combined_p,
    max_log_fc: max_lfc,
    n_windows,
    window_indices: win_indices
  }
}

// ===== Helper functions =====

fn csw_mean_double(arr : Array[Double]) -> Double {
  let n = arr.length()
  if n == 0 { return 0.0 }
  let mut sum = 0.0
  let mut i = 0
  while i < n {
    sum = sum + arr[i]
    i = i + 1
  }
  sum / n.to_double()
}

fn csw_mean_int(arr : Array[Int]) -> Double {
  let n = arr.length()
  if n == 0 { return 0.0 }
  let mut sum = 0.0
  let mut i = 0
  while i < n {
    sum = sum + arr[i].to_double()
    i = i + 1
  }
  sum / n.to_double()
}

fn csw_total_cpb(counts : Array[Double], lib_sizes : Array[Double]) -> Double {
  let mut total = 0.0
  let mut i = 0
  while i < counts.length() {
    let lib = if lib_sizes[i] > 0.0 { lib_sizes[i] } else { 1.0 }
    total = total + counts[i] / lib
    i = i + 1
  }
  total
}

fn csw_mean_var_group(counts : Array[Double], indices : Array[Int]) -> (Double, Double) {
  let n = indices.length()
  if n == 0 { return (0.0, 0.0) }
  let mut sum = 0.0
  let mut i = 0
  while i < n {
    sum = sum + counts[indices[i]]
    i = i + 1
  }
  let mean = sum / n.to_double()
  let mut variance = 0.0
  let mut j = 0
  while j < n {
    let diff = counts[indices[j]] - mean
    variance = variance + diff * diff
    j = j + 1
  }
  variance = if n > 1 { variance / (n.to_double() - 1.0) } else { 0.0 }
  (mean, variance)
}

///|
/// Approximate two-sided p-value from normal distribution.
fn csw_norm_p_value(z : Double) -> Double {
  let abs_z = z.abs()
  let p_approx = 2.0 * (1.0 - csw_norm_cdf(abs_z))
  if p_approx < 1.0e-300 { 1.0e-300 } else { p_approx }
}

///|
/// Normal CDF approximation.
fn csw_norm_cdf(x : Double) -> Double {
  let t = 1.0 / (1.0 + 0.2316419 * x)
  let d = 0.3989422804014327 * @math.exp(-0.5 * x * x)
  let p = d * t * (0.3193815 + t * (-0.3565638 + t * (1.781478 + t * (-1.821256 + t * 1.3302744))))
  1.0 - p
}

///|
/// Benjamini-Hochberg correction.
fn csw_bh_correction(p_values : Array[Double]) -> Array[Double] {
  let n = p_values.length()
  if n == 0 { return Array::new() }
  let indexed : Array[(Double, Int)] = Array::new()
  let mut i = 0
  while i < n {
    indexed.push((p_values[i], i))
    i = i + 1
  }
  // Sort by p-value
  let sorted = indexed.copy()
  let mut a = 0
  while a < sorted.length() {
    let mut b = a + 1
    while b < sorted.length() {
      if sorted[b].0 < sorted[a].0 {
        let tmp = sorted[a]
        sorted[a] = sorted[b]
        sorted[b] = tmp
      }
      b = b + 1
    }
    a = a + 1
  }
  // Compute adjusted p-values
  let adj : Array[Double] = Array::new()
  let mut j = 0
  while j < n {
    let rank = j + 1
    let val = sorted[j].0 * n.to_double() / rank.to_double()
    adj.push(val)
    j = j + 1
  }
  // Enforce monotonicity
  let mut k = n - 2
  while k >= 0 {
    if adj[k] > adj[k + 1] {
      adj[k] = adj[k + 1]
    }
    k = k - 1
  }
  // Map back to original order
  let result : Array[Double] = Array::new()
  let mut m = 0
  while m < n {
    result.push(0.0)
    m = m + 1
  }
  let mut p = 0
  while p < n {
    result[sorted[p].1] = if adj[p] < 1.0 { adj[p] } else { 1.0 }
    p = p + 1
  }
  result
}

///|
/// Chi-squared survival function approximation.
fn csw_chi2_sf(x : Double, df : Double) -> Double {
  if x <= 0.0 { return 1.0 }
  if df <= 0.0 { return 1.0 }
  // Wilson-Hilferty approximation
  let z = @math.pow(x / df, 1.0 / 3.0)
  let mu = 1.0 - 2.0 / (9.0 * df)
  let sigma_val = (2.0 / (9.0 * df)).sqrt()
  let z_approx = (z - mu) / sigma_val
  csw_norm_cdf(z_approx)
}
