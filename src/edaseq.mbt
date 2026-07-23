///|
/// EDASeq module, inspired by Bioconductor EDASeq package.
/// Provides RNA-seq exploratory data analysis with GC-content and
/// gene length normalization, plus between/within-lane normalization.

/// Normalization type
pub enum EDASeqNormType {
  /// Within-lane normalization (GC content or gene length)
  WithinLane
  /// Between-lane normalization
  BetweenLane
  /// Full normalization (within + between)
  Full
}

/// Within-lane normalization method
pub enum EDASeqWithinMethod {
  /// GC content normalization
  GCContent
  /// Gene length normalization
  GeneLength
  /// Both GC and length (loess)
  Both
}

/// Between-lane normalization method
pub enum EDASeqBetweenMethod {
  /// Median ratio normalization
  MedianRatio
  /// Quantile normalization
  Quantile
  /// Library size normalization
  Library
}

/// Gene annotation for normalization
pub struct EDASeqGeneAnno {
  /// Gene ID
  gene_id : String
  /// GC content (0-1)
  gc_content : Double
  /// Gene length (bp)
  length : Double
  /// Effective length (adjusted for fragment size)
  effective_length : Double
}

/// Expression set for EDASeq
pub struct EDASeqDataSet {
  /// Gene IDs
  gene_ids : Array[String]
  /// Sample names
  sample_ids : Array[String]
  /// Count matrix (genes x samples)
  counts : Array[Array[Double]]
  /// Gene annotations
  annotations : Array[EDASeqGeneAnno]
  /// Library sizes
  lib_sizes : Array[Double]
}

/// Result of within-lane normalization
pub struct EDASeqWithinResult {
  /// Normalized counts
  norm_counts : Array[Array[Double]]
  /// GC loess fits per sample
  gc_fits : Array[EDASeqFit]
  /// Length loess fits per sample
  length_fits : Array[EDASeqFit]
  /// Method used
  norm_method : EDASeqWithinMethod
}

/// A loess fit result
pub struct EDASeqFit {
  /// Smoothing span
  span : Double
  /// Fitted values
  fitted : Array[Double]
  /// Residuals
  residuals : Array[Double]
  /// Whether fit was successful
  converged : Bool
}

/// Result of between-lane normalization
pub struct EDASeqBetweenResult {
  /// Normalized counts
  norm_counts : Array[Array[Double]]
  /// Size factors per sample
  size_factors : Array[Double]
  /// Method used
  norm_method : EDASeqBetweenMethod
}

/// Complete EDASeq normalization result
pub struct EDASeqResult {
  /// Within-lane normalized counts (if applied)
  within_norm_counts : Array[Array[Double]]
  /// Between-lane normalized counts (if applied)
  between_norm_counts : Array[Array[Double]]
  /// Final normalized counts
  norm_counts : Array[Array[Double]]
  /// Method type
  norm_type : EDASeqNormType
  /// Gene IDs
  gene_ids : Array[String]
  /// Sample IDs
  sample_ids : Array[String]
}

/// Parameters for EDASeq normalization
pub struct EDASeqParams {
  /// Normalization type
  norm_type : EDASeqNormType
  /// Within-lane method
  within_method : EDASeqWithinMethod
  /// Between-lane method
  between_method : EDASeqBetweenMethod
  /// Loess span for GC/length normalization
  loess_span : Double
  /// Minimum count threshold
  min_count : Double
  /// Whether to log-transform
  log_transform : Bool
}

/// Create default parameters
pub fn EDASeqParams::new() -> EDASeqParams {
  EDASeqParams::{
    norm_type: EDASeqNormType::Full,
    within_method: EDASeqWithinMethod::GCContent,
    between_method: EDASeqBetweenMethod::MedianRatio,
    loess_span: 0.75,
    min_count: 0.5,
    log_transform: true
  }
}

// ===== Gene annotation =====

/// Create a new gene annotation
pub fn EDASeqGeneAnno::new(gene_id : String, gc_content : Double, length : Double) -> EDASeqGeneAnno {
  let eff_len = if length > 200.0 { length } else { 200.0 }
  EDASeqGeneAnno::{
    gene_id,
    gc_content,
    length,
    effective_length: eff_len
  }
}

/// Compute RPKM scaling factor (1000/gene_length)
pub fn EDASeqGeneAnno::eda_rpkm_scale(self : EDASeqGeneAnno) -> Double {
  if self.effective_length > 0.0 {
    1000.0 / self.effective_length
  } else {
    1.0
  }
}

// ===== DataSet methods =====

/// Create a new EDASeqDataSet
pub fn EDASeqDataSet::new(
  gene_ids : Array[String],
  sample_ids : Array[String],
  counts : Array[Array[Double]],
  annotations : Array[EDASeqGeneAnno]
) -> EDASeqDataSet {
  let n_samples = sample_ids.length()
  let lib_sizes : Array[Double] = Array::new()
  let mut s = 0
  while s < n_samples {
    let mut sum = 0.0
    let mut g = 0
    while g < counts.length() {
      sum = sum + counts[g][s]
      g = g + 1
    }
    lib_sizes.push(sum)
    s = s + 1
  }
  EDASeqDataSet::{
    gene_ids,
    sample_ids,
    counts,
    annotations,
    lib_sizes
  }
}

/// Get number of genes
pub fn EDASeqDataSet::eda_n_genes(self : EDASeqDataSet) -> Int {
  self.gene_ids.length()
}

/// Get number of samples
pub fn EDASeqDataSet::eda_n_samples(self : EDASeqDataSet) -> Int {
  self.sample_ids.length()
}

/// Get a gene's counts across samples
pub fn EDASeqDataSet::eda_gene_counts(self : EDASeqDataSet, gene_idx : Int) -> Array[Double] {
  self.counts[gene_idx].copy()
}

/// Get a sample's counts across genes
pub fn EDASeqDataSet::eda_sample_counts(self : EDASeqDataSet, sample_idx : Int) -> Array[Double] {
  let result : Array[Double] = Array::new()
  let mut i = 0
  while i < self.counts.length() {
    result.push(self.counts[i][sample_idx])
    i = i + 1
  }
  result
}

// ===== Within-lane normalization =====

///|
/// Compute GC content-based within-lane normalization.
pub fn eda_normalize_within_lane(data : EDASeqDataSet, params : EDASeqParams) -> EDASeqWithinResult {
  let n_genes = data.gene_ids.length()
  let n_samples = data.sample_ids.length()
  let gc_fits : Array[EDASeqFit] = Array::new()
  let length_fits : Array[EDASeqFit] = Array::new()
  let norm_counts : Array[Array[Double]] = Array::new()

  let mut s = 0
  while s < n_samples {
    let sample_counts = data.eda_sample_counts(s)
    let log_counts : Array[Double] = Array::new()
    let gc_values : Array[Double] = Array::new()
    let length_values : Array[Double] = Array::new()
    let mut g = 0
    while g < n_genes {
      let c = if sample_counts[g] > 0.0 { sample_counts[g] } else { params.min_count }
      log_counts.push(@math.ln(c))
      gc_values.push(data.annotations[g].gc_content)
      length_values.push(data.annotations[g].effective_length)
      g = g + 1
    }
    // Apply within-lane normalization
    let (norm_sample_log, gc_fit, len_fit) = eda_apply_within_norm(
      log_counts, gc_values, length_values, data.lib_sizes[s], params
    )
    gc_fits.push(gc_fit)
    length_fits.push(len_fit)
    let col_norm : Array[Double] = Array::new()
    let mut ni = 0
    while ni < n_genes {
      col_norm.push(@math.exp(norm_sample_log[ni]))
      ni = ni + 1
    }
    norm_counts.push(col_norm)
    s = s + 1
  }

  // Transpose norm_counts from samples x genes to genes x samples
  let result : Array[Array[Double]] = Array::new()
  let mut g2 = 0
  while g2 < n_genes {
    let row : Array[Double] = Array::new()
    let mut s2 = 0
    while s2 < n_samples {
      row.push(norm_counts[s2][g2])
      s2 = s2 + 1
    }
    result.push(row)
    g2 = g2 + 1
  }

  EDASeqWithinResult::{
    norm_counts: result,
    gc_fits,
    length_fits,
    norm_method: params.within_method
  }
}

///|
/// Apply within-lane normalization for a single sample.
fn eda_apply_within_norm(
  log_counts : Array[Double],
  gc_values : Array[Double],
  length_values : Array[Double],
  lib_size : Double,
  params : EDASeqParams
) -> (Array[Double], EDASeqFit, EDASeqFit) {
  let n = log_counts.length()
  let gc_fit_default = EDASeqFit::{
    span: params.loess_span,
    fitted: Array::new(),
    residuals: Array::new(),
    converged: false
  }
  let len_fit_default = EDASeqFit::{
    span: params.loess_span,
    fitted: Array::new(),
    residuals: Array::new(),
    converged: false
  }
  let lib_adj = if lib_size > 0.0 { lib_size } else { 1.0 }
  let log_lib = @math.ln(lib_adj)
  // Normalize by library size first
  let base_norm : Array[Double] = Array::new()
  let mut i = 0
  while i < n {
    base_norm.push(log_counts[i] - log_lib)
    i = i + 1
  }
  // Apply GC correction
  let (gc_corrected, gc_fit) = eda_gc_correction(base_norm, gc_values, params.loess_span)
  // Apply length correction
  let (final_norm, len_fit) = eda_length_correction(gc_corrected, length_values, params.loess_span)
  (final_norm, gc_fit, len_fit)
}

///|
/// GC content-based loess correction.
fn eda_gc_correction(log_counts : Array[Double], gc_values : Array[Double], span : Double) -> (Array[Double], EDASeqFit) {
  let n = log_counts.length()
  if n < 4 {
    let result = log_counts.copy()
    let fit = EDASeqFit::{ span, fitted: Array::new(), residuals: Array::new(), converged: false }
    return (result, fit)
  }
  // Apply loess smoothing
  let fitted = eda_loess_smooth(gc_values, log_counts, span)
  let residuals : Array[Double] = Array::new()
  let corrected : Array[Double] = Array::new()
  let mut i = 0
  while i < n {
    let res = log_counts[i] - fitted[i]
    residuals.push(res)
    // Add back the overall mean to keep the data centered
    corrected.push(res)
    i = i + 1
  }
  // Add mean back to preserve overall expression levels
  let mean_res = eda_mean(corrected)
  let final_corrected : Array[Double] = Array::new()
  let mut j = 0
  while j < n {
    final_corrected.push(corrected[j] + mean_res)
    j = j + 1
  }
  (final_corrected, EDASeqFit::{ span, fitted, residuals, converged: true })
}

///|
/// Gene length-based loess correction.
fn eda_length_correction(log_counts : Array[Double], length_values : Array[Double], span : Double) -> (Array[Double], EDASeqFit) {
  let n = log_counts.length()
  if n < 4 {
    let result = log_counts.copy()
    let fit = EDASeqFit::{ span, fitted: Array::new(), residuals: Array::new(), converged: false }
    return (result, fit)
  }
  // Use log-length for smoothing
  let log_lengths : Array[Double] = Array::new()
  let mut i = 0
  while i < n {
    let l = if length_values[i] > 1.0 { length_values[i] } else { 1.0 }
    log_lengths.push(@math.ln(l))
    i = i + 1
  }
  let fitted = eda_loess_smooth(log_lengths, log_counts, span)
  let residuals : Array[Double] = Array::new()
  let corrected : Array[Double] = Array::new()
  let mut j = 0
  while j < n {
    let res = log_counts[j] - fitted[j]
    residuals.push(res)
    corrected.push(res)
    j = j + 1
  }
  let mean_res = eda_mean(corrected)
  let final_corrected : Array[Double] = Array::new()
  let mut k = 0
  while k < n {
    final_corrected.push(corrected[k] + mean_res)
    k = k + 1
  }
  (final_corrected, EDASeqFit::{ span, fitted, residuals, converged: true })
}

///|
/// Simple loess smoothing (locally weighted regression).
fn eda_loess_smooth(x : Array[Double], y : Array[Double], span : Double) -> Array[Double] {
  let n = x.length()
  if n < 3 { return y.copy() }
  let result : Array[Double] = Array::new()
  let half_span = span * 0.5
  let mut i = 0
  while i < n {
    let xi = x[i]
    // Collect neighbor indices
    let neighbors : Array[(Double, Double)] = Array::new()
    let mut j = 0
    while j < n {
      let dist = (x[j] - xi).abs()
      if dist <= half_span {
        neighbors.push((x[j], y[j]))
      }
      j = j + 1
    }
    let n_neighbors = neighbors.length()
    let effective_n = if n_neighbors >= 2 { n_neighbors } else { n }
    // Compute weighted regression
    let mut sum_wx = 0.0
    let mut sum_wy = 0.0
    let mut sum_wxx = 0.0
    let mut sum_wxy = 0.0
    let max_dist = if half_span > 0.0 { half_span } else { 1.0 }
    let mut k = 0
    while k < effective_n {
      let (xn, yn) = if k < n_neighbors { neighbors[k] } else { (x[k], y[k]) }
      let dist = (xn - xi).abs()
      let u = if max_dist > 0.0 { dist / max_dist } else { 0.0 }
      let w = if u <= 1.0 { (1.0 - u * u) * (1.0 - u * u) } else { 0.0 }
      sum_wx = sum_wx + w * xn
      sum_wy = sum_wy + w * yn
      sum_wxx = sum_wxx + w * xn * xn
      sum_wxy = sum_wxy + w * xn * yn
      k = k + 1
    }
    let denom = effective_n.to_double() * sum_wxx - sum_wx * sum_wx
    let intercept = if denom.abs() > 1.0e-10 {
      (sum_wy * sum_wxx - sum_wx * sum_wxy) / denom
    } else {
      eda_mean(y)
    }
    let slope = if denom.abs() > 1.0e-10 {
      (effective_n.to_double() * sum_wxy - sum_wx * sum_wy) / denom
    } else {
      0.0
    }
    let fitted = slope * xi + intercept
    result.push(fitted)
    i = i + 1
  }
  result
}

// ===== Between-lane normalization =====

///|
/// Compute between-lane normalization.
pub fn eda_normalize_between_lane(data : EDASeqDataSet, params : EDASeqParams) -> EDASeqBetweenResult {
  let n_genes = data.gene_ids.length()
  let n_samples = data.sample_ids.length()
  // Compute size factors
  let size_factors = eda_compute_size_factors(data.counts, params)
  // Apply normalization
  let norm_counts : Array[Array[Double]] = Array::new()
  let mut g = 0
  while g < n_genes {
    let row : Array[Double] = Array::new()
    let mut s = 0
    while s < n_samples {
      let sf = if size_factors[s] > 0.0 { size_factors[s] } else { 1.0 }
      row.push(data.counts[g][s] / sf)
      s = s + 1
    }
    norm_counts.push(row)
    g = g + 1
  }
  EDASeqBetweenResult::{
    norm_counts,
    size_factors,
    norm_method: params.between_method
  }
}

///|
/// Compute size factors for between-lane normalization.
fn eda_compute_size_factors(counts : Array[Array[Double]], params : EDASeqParams) -> Array[Double] {
  let n_genes = counts.length()
  let n_samples = if n_genes > 0 { counts[0].length() } else { 0 }
  let size_factors : Array[Double] = Array::new()
  if n_samples == 0 { return size_factors }

  match params.between_method {
    EDASeqBetweenMethod::MedianRatio => {
      // DESeq2-style median ratio normalization
      // Compute geometric means per gene
      let geo_means : Array[Double] = Array::new()
      let mut g = 0
      while g < n_genes {
        let row = counts[g]
        let mut log_sum = 0.0
        let mut s = 0
        while s < n_samples {
          let c = if row[s] > 0.0 { row[s] } else { params.min_count }
          log_sum = log_sum + @math.ln(c)
          s = s + 1
        }
        geo_means.push(@math.exp(log_sum / n_samples.to_double()))
        g = g + 1
      }
      // Compute size factors as median of ratios
      let mut s2 = 0
      while s2 < n_samples {
        let ratios : Array[Double] = Array::new()
        let mut g2 = 0
        while g2 < n_genes {
          let c = if counts[g2][s2] > 0.0 { counts[g2][s2] } else { params.min_count }
          if geo_means[g2] > 0.0 {
            ratios.push(c / geo_means[g2])
          }
          g2 = g2 + 1
        }
        size_factors.push(eda_median(ratios))
        s2 = s2 + 1
      }
    }
    EDASeqBetweenMethod::Quantile => {
      // Quantile normalization
      let mut s = 0
      while s < n_samples {
        let col : Array[Double] = Array::new()
        let mut g = 0
        while g < n_genes {
          col.push(counts[g][s])
          g = g + 1
        }
        size_factors.push(eda_mean(col))
        s = s + 1
      }
    }
    EDASeqBetweenMethod::Library => {
      // Library size normalization
      let total : Array[Double] = Array::new()
      let mut s = 0
      while s < n_samples {
        let mut sum = 0.0
        let mut g = 0
        while g < n_genes {
          sum = sum + counts[g][s]
          g = g + 1
        }
        total.push(sum)
        s = s + 1
      }
      let mean_total = eda_mean(total)
      let mut s2 = 0
      while s2 < n_samples {
        let sf = if total[s2] > 0.0 { total[s2] / mean_total } else { 1.0 }
        size_factors.push(sf)
        s2 = s2 + 1
      }
    }
  }
  size_factors
}

/// Compute median
fn eda_median(arr : Array[Double]) -> Double {
  let n = arr.length()
  if n == 0 { return 1.0 }
  let sorted = eda_sort(arr)
  if n % 2 == 1 { sorted[n / 2] } else { (sorted[n / 2 - 1] + sorted[n / 2]) / 2.0 }
}

/// Compute mean
fn eda_mean(arr : Array[Double]) -> Double {
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

/// Sort array
fn eda_sort(arr : Array[Double]) -> Array[Double] {
  let result = arr.copy()
  let n = result.length()
  let mut i = 1
  while i < n {
    let key = result[i]
    let mut j = i - 1
    while j >= 0 && result[j] > key {
      result[j + 1] = result[j]
      j = j - 1
    }
    result[j + 1] = key
    i = i + 1
  }
  result
}

// ===== Full normalization pipeline =====

///|
/// Run the complete EDASeq normalization pipeline.
pub fn eda_run(data : EDASeqDataSet, params : EDASeqParams) -> EDASeqResult {
  let n_genes = data.gene_ids.length()
  let n_samples = data.sample_ids.length()
  let within_counts : Array[Array[Double]] = Array::new()
  let between_counts : Array[Array[Double]] = Array::new()

  match params.norm_type {
    EDASeqNormType::WithinLane => {
      let within = eda_normalize_within_lane(data, params)
      EDASeqResult::{
        within_norm_counts: within.norm_counts,
        between_norm_counts: Array::new(),
        norm_counts: within.norm_counts,
        norm_type: params.norm_type,
        gene_ids: data.gene_ids,
        sample_ids: data.sample_ids
      }
    }
    EDASeqNormType::BetweenLane => {
      let between = eda_normalize_between_lane(data, params)
      EDASeqResult::{
        within_norm_counts: Array::new(),
        between_norm_counts: between.norm_counts,
        norm_counts: between.norm_counts,
        norm_type: params.norm_type,
        gene_ids: data.gene_ids,
        sample_ids: data.sample_ids
      }
    }
    EDASeqNormType::Full => {
      // Apply within-lane first, then between-lane
      let within = eda_normalize_within_lane(data, params)
      let within_data = EDASeqDataSet::new(
        data.gene_ids,
        data.sample_ids,
        within.norm_counts,
        data.annotations
      )
      let between = eda_normalize_between_lane(within_data, params)
      EDASeqResult::{
        within_norm_counts: within.norm_counts,
        between_norm_counts: between.norm_counts,
        norm_counts: between.norm_counts,
        norm_type: params.norm_type,
        gene_ids: data.gene_ids,
        sample_ids: data.sample_ids
      }
    }
  }
}

// ===== RPKM/TPM calculation =====

///|
/// Compute RPKM (Reads Per Kilobase per Million) normalization.
pub fn eda_rpkm(counts : Array[Array[Double]], annotations : Array[EDASeqGeneAnno], lib_sizes : Array[Double]) -> Array[Array[Double]] {
  let n_genes = counts.length()
  let n_samples = if n_genes > 0 { counts[0].length() } else { 0 }
  let rpkm : Array[Array[Double]] = Array::new()
  let mut g = 0
  while g < n_genes {
    let row : Array[Double] = Array::new()
    let anno = annotations[g]
    let scale = anno.eda_rpkm_scale()
    let mut s = 0
    while s < n_samples {
      let lib_scale = if lib_sizes[s] > 0.0 { 1.0e6 / lib_sizes[s] } else { 0.0 }
      let val = counts[g][s] * scale * lib_scale
      row.push(val)
      s = s + 1
    }
    rpkm.push(row)
    g = g + 1
  }
  rpkm
}

///|
/// Compute TPM (Transcripts Per Million) normalization.
pub fn eda_tpm(counts : Array[Array[Double]], annotations : Array[EDASeqGeneAnno], lib_sizes : Array[Double]) -> Array[Array[Double]] {
  let n_genes = counts.length()
  let n_samples = if n_genes > 0 { counts[0].length() } else { 0 }
  let tpm : Array[Array[Double]] = Array::new()
  let mut g = 0
  while g < n_genes {
    let row : Array[Double] = Array::new()
    let anno = annotations[g]
    let rpk_scale = anno.eda_rpkm_scale()
    let mut s = 0
    while s < n_samples {
      let rpkm_val = counts[g][s] * rpk_scale
      row.push(rpkm_val)
      s = s + 1
    }
    tpm.push(row)
    g = g + 1
  }
  // Convert RPK to TPM by normalizing to million
  let result : Array[Array[Double]] = Array::new()
  let mut s2 = 0
  while s2 < n_samples {
    let col : Array[Double] = Array::new()
    let mut total_rpk = 0.0
    let mut g2 = 0
    while g2 < n_genes {
      total_rpk = total_rpk + tpm[g2][s2]
      g2 = g2 + 1
    }
    let scale = if total_rpk > 0.0 { 1.0e6 / total_rpk } else { 0.0 }
    let mut g3 = 0
    while g3 < n_genes {
      col.push(tpm[g3][s2] * scale)
      g3 = g3 + 1
    }
    result.push(col)
    s2 = s2 + 1
  }
  // Transpose back to genes x samples
  let final_result : Array[Array[Double]] = Array::new()
  let mut g4 = 0
  while g4 < n_genes {
    let row : Array[Double] = Array::new()
    let mut s3 = 0
    while s3 < n_samples {
      row.push(result[s3][g4])
      s3 = s3 + 1
    }
    final_result.push(row)
    g4 = g4 + 1
  }
  final_result
}

// ===== Accessor methods =====

/// Get a normalized count
pub fn EDASeqResult::eda_get_count(self : EDASeqResult, gene_idx : Int, sample_idx : Int) -> Double {
  self.norm_counts[gene_idx][sample_idx]
}

/// Get normalized counts for a gene
pub fn EDASeqResult::eda_get_gene(self : EDASeqResult, gene_idx : Int) -> Array[Double] {
  self.norm_counts[gene_idx].copy()
}

/// Get normalized counts for a sample
pub fn EDASeqResult::eda_get_sample(self : EDASeqResult, sample_idx : Int) -> Array[Double] {
  let result : Array[Double] = Array::new()
  let mut i = 0
  while i < self.norm_counts.length() {
    result.push(self.norm_counts[i][sample_idx])
    i = i + 1
  }
  result
}

/// Get number of genes
pub fn EDASeqResult::eda_n_genes(self : EDASeqResult) -> Int {
  self.gene_ids.length()
}

/// Get number of samples
pub fn EDASeqResult::eda_n_samples(self : EDASeqResult) -> Int {
  self.sample_ids.length()
}

/// Get gene ID
pub fn EDASeqResult::eda_get_gene_id(self : EDASeqResult, idx : Int) -> String {
  self.gene_ids[idx]
}

/// Get sample ID
pub fn EDASeqResult::eda_get_sample_id(self : EDASeqResult, idx : Int) -> String {
  self.sample_ids[idx]
}
