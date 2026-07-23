///|
/// DSS module, inspired by Bioconductor DSS package.
/// Provides dispersion shrinkage and differential analysis for sequencing data.
/// Supports RNA-seq differential expression analysis with dispersion estimation,
/// and differential methylation analysis for bisulfite sequencing data.

/// Dispersion estimation methods
pub enum DSSDispMethod {
  /// Maximum likelihood estimation
  ML
  /// Shrinkage estimator
  Shrink
}

/// Result of dispersion estimation
pub struct DSSDispResult {
  /// Dispersion values for each gene/probe
  dispersions : Array[Double]
  /// Shrunken dispersions
  shrunken_dispersions : Array[Double]
  /// Prior degrees of freedom
  prior_df : Double
  /// Prior variance
  prior_var : Double
}

/// Gene count data structure
pub struct DSSCountData {
  /// Gene IDs
  gene_ids : Array[String]
  /// Count matrix (genes x samples)
  counts : Array[Array[Double]]
  /// Sample groups
  groups : Array[String]
  /// Library sizes
  lib_sizes : Array[Double]
}

/// Result of differential expression analysis
pub struct DSSDEResult {
  /// Gene IDs
  gene_ids : Array[String]
  /// Mean expression in group 1
  mean1 : Array[Double]
  /// Mean expression in group 2
  mean2 : Array[Double]
  /// Dispersion
  dispersion : Array[Double]
  /// Wald statistics
  wald_stat : Array[Double]
  /// P-values
  p_values : Array[Double]
  /// Adjusted p-values
  adj_p_values : Array[Double]
  /// Significant gene indices
  significant : Array[Bool]
}

/// CpG site methylation data
pub struct DSSMethylData {
  /// CpG site IDs
  cpg_ids : Array[String]
  /// Chromosome
  chrom : Array[String]
  /// Position
  pos : Array[Int]
  /// Methylation counts
  meth_counts : Array[Array[Double]]
  /// Total counts
  total_counts : Array[Array[Double]]
  /// Sample groups
  groups : Array[String]
}

/// Result of differential methylation analysis
pub struct DSSDMResult {
  /// CpG site IDs
  cpg_ids : Array[String]
  /// Mean methylation in group 1
  meth_mean1 : Array[Double]
  /// Mean methylation in group 2
  meth_mean2 : Array[Double]
  /// Methylation difference
  delta_meth : Array[Double]
  /// Test statistics
  stats : Array[Double]
  /// P-values
  p_values : Array[Double]
  /// Adjusted p-values
  adj_p_values : Array[Double]
  /// Significant site indices
  significant : Array[Bool]
}

/// Result of DMR detection
pub struct DSSDMRResult {
  /// DMR IDs
  dmr_ids : Array[String]
  /// Chromosome
  chrom : Array[String]
  /// Start position
  start : Array[Int]
  /// End position
  end : Array[Int]
  /// Number of CpGs in DMR
  n_cpgs : Array[Int]
  /// Mean methylation difference
  mean_delta : Array[Double]
  /// DMR statistics
  stats : Array[Double]
  /// P-values
  p_values : Array[Double]
}

// ===== DSSDispResult methods =====

/// Create new dispersion result
pub fn DSSDispResult::new(dispersions : Array[Double], shrunken_dispersions : Array[Double], prior_df : Double, prior_var : Double) -> DSSDispResult {
  DSSDispResult::{ dispersions, shrunken_dispersions, prior_df, prior_var }
}

/// Get dispersion values
pub fn DSSDispResult::dss_get_dispersions(self : DSSDispResult) -> Array[Double] {
  self.dispersions
}

/// Get shrunken dispersion values
pub fn DSSDispResult::dss_get_shrunken_dispersions(self : DSSDispResult) -> Array[Double] {
  self.shrunken_dispersions
}

// ===== DSSCountData methods =====

/// Create new count data
pub fn DSSCountData::new(gene_ids : Array[String], counts : Array[Array[Double]], groups : Array[String], lib_sizes : Array[Double]) -> DSSCountData {
  DSSCountData::{ gene_ids, counts, groups, lib_sizes }
}

/// Get gene count for specific gene and sample
pub fn DSSCountData::dss_get_count(self : DSSCountData, gene_idx : Int, sample_idx : Int) -> Double {
  self.counts[gene_idx][sample_idx]
}

/// Get number of genes
pub fn DSSCountData::dss_n_genes(self : DSSCountData) -> Int {
  self.gene_ids.length()
}

/// Get number of samples
pub fn DSSCountData::dss_n_samples(self : DSSCountData) -> Int {
  self.groups.length()
}

/// Normalize counts by library size
pub fn DSSCountData::dss_normalize(self : DSSCountData) -> Array[Array[Double]] {
  let n_genes = self.counts.length()
  let n_samples = if n_genes > 0 { self.counts[0].length() } else { 0 }
  let normalized : Array[Array[Double]] = Array::new()
  
  let mut i = 0
  while i < n_genes {
    let row : Array[Double] = Array::new()
    let mut j = 0
    while j < n_samples {
      let lib_size = if j < self.lib_sizes.length() { self.lib_sizes[j] } else { 1.0 }
      row.push(self.counts[i][j] / lib_size * 1000000.0)
      j = j + 1
    }
    normalized.push(row)
    i = i + 1
  }
  
  normalized
}

/// Filter low count genes
pub fn DSSCountData::dss_filter_low_counts(self : DSSCountData, min_count : Double, min_samples : Int) -> DSSCountData {
  let new_gene_ids : Array[String] = Array::new()
  let new_counts : Array[Array[Double]] = Array::new()
  
  let mut i = 0
  while i < self.gene_ids.length() {
    let count_row = self.counts[i]
    let mut passing_samples = 0
    
    let mut j = 0
    while j < count_row.length() {
      if count_row[j] >= min_count {
        passing_samples = passing_samples + 1
      }
      j = j + 1
    }
    
    if passing_samples >= min_samples {
      new_gene_ids.push(self.gene_ids[i])
      new_counts.push(count_row)
    }
    
    i = i + 1
  }
  
  DSSCountData::new(new_gene_ids, new_counts, self.groups, self.lib_sizes)
}

// ===== DSSDEResult methods =====

/// Create new DSR result
pub fn DSSDEResult::new(gene_ids : Array[String], mean1 : Array[Double], mean2 : Array[Double], 
             dispersion : Array[Double], wald_stat : Array[Double], 
             p_values : Array[Double], adj_p_values : Array[Double],
             significant : Array[Bool]) -> DSSDEResult {
  DSSDEResult::{ gene_ids, mean1, mean2, dispersion, wald_stat, p_values, adj_p_values, significant }
}

/// Get number of genes
pub fn DSSDEResult::dss_n_genes(self : DSSDEResult) -> Int {
  self.gene_ids.length()
}

/// Get significant gene count
pub fn DSSDEResult::dss_n_significant(self : DSSDEResult) -> Int {
  let mut count = 0
  let mut i = 0
  while i < self.significant.length() {
    if self.significant[i] {
      count = count + 1
    }
    i = i + 1
  }
  count
}

/// Get gene ID at index
pub fn DSSDEResult::dss_get_gene_id(self : DSSDEResult, idx : Int) -> String {
  self.gene_ids[idx]
}

/// Get p-value at index
pub fn DSSDEResult::dss_get_pvalue(self : DSSDEResult, idx : Int) -> Double {
  self.p_values[idx]
}

/// Get adjusted p-value at index
pub fn DSSDEResult::dss_get_adj_pvalue(self : DSSDEResult, idx : Int) -> Double {
  self.adj_p_values[idx]
}

/// Get fold change at index (log2 scale)
pub fn DSSDEResult::dss_get_logfc(self : DSSDEResult, idx : Int) -> Double {
  if self.mean1[idx] > 0.0 && self.mean2[idx] > 0.0 {
    @math.log2(self.mean2[idx] / self.mean1[idx])
  } else {
    0.0
  }
}

// ===== DSSMethylData methods =====

/// Create new methylation data
pub fn DSSMethylData::new(cpg_ids : Array[String], chrom : Array[String], pos : Array[Int],
             meth_counts : Array[Array[Double]], total_counts : Array[Array[Double]],
             groups : Array[String]) -> DSSMethylData {
  DSSMethylData::{ cpg_ids, chrom, pos, meth_counts, total_counts, groups }
}

/// Get number of CpG sites
pub fn DSSMethylData::dss_n_cpgs(self : DSSMethylData) -> Int {
  self.cpg_ids.length()
}

/// Get number of samples
pub fn DSSMethylData::dss_n_samples(self : DSSMethylData) -> Int {
  self.groups.length()
}

/// Calculate methylation proportion for a site and sample
pub fn DSSMethylData::dss_get_meth_prop(self : DSSMethylData, cpg_idx : Int, sample_idx : Int) -> Double {
  let total = self.total_counts[cpg_idx][sample_idx]
  if total > 0.0 {
    self.meth_counts[cpg_idx][sample_idx] / total
  } else {
    0.0
  }
}

// ===== DSSDMResult methods =====

/// Create new DM result
pub fn DSSDMResult::new(cpg_ids : Array[String], meth_mean1 : Array[Double], meth_mean2 : Array[Double],
             delta_meth : Array[Double], stats : Array[Double],
             p_values : Array[Double], adj_p_values : Array[Double],
             significant : Array[Bool]) -> DSSDMResult {
  DSSDMResult::{ cpg_ids, meth_mean1, meth_mean2, delta_meth, stats, p_values, adj_p_values, significant }
}

/// Get number of sites
pub fn DSSDMResult::dss_n_sites(self : DSSDMResult) -> Int {
  self.cpg_ids.length()
}

/// Get significant site count
pub fn DSSDMResult::dss_n_significant(self : DSSDMResult) -> Int {
  let mut count = 0
  let mut i = 0
  while i < self.significant.length() {
    if self.significant[i] {
      count = count + 1
    }
    i = i + 1
  }
  count
}

/// Get delta methylation at index
pub fn DSSDMResult::dss_get_delta_meth(self : DSSDMResult, idx : Int) -> Double {
  self.delta_meth[idx]
}

// ===== DSSDMRResult methods =====

/// Create new DMR result
pub fn DSSDMRResult::new(dmr_ids : Array[String], chrom : Array[String], start : Array[Int],
             end : Array[Int], n_cpgs : Array[Int], mean_delta : Array[Double],
             stats : Array[Double], p_values : Array[Double]) -> DSSDMRResult {
  DSSDMRResult::{ dmr_ids, chrom, start, end, n_cpgs, mean_delta, stats, p_values }
}

/// Get number of DMRs
pub fn DSSDMRResult::dss_n_dmrs(self : DSSDMRResult) -> Int {
  self.dmr_ids.length()
}

/// Get DMR at index
pub fn DSSDMRResult::dss_get_dmr(self : DSSDMRResult, idx : Int) -> (String, String, Int, Int, Int, Double) {
  (self.dmr_ids[idx], self.chrom[idx], self.start[idx], self.end[idx], 
   self.n_cpgs[idx], self.mean_delta[idx])
}

// ===== Helper functions =====

/// Calculate mean of array
fn dss_mean(arr : Array[Double]) -> Double {
  if arr.length() == 0 {
    return 0.0
  }
  let mut sum = 0.0
  let mut i = 0
  while i < arr.length() {
    sum = sum + arr[i]
    i = i + 1
  }
  sum / arr.length().to_double()
}

/// Calculate variance of array
fn dss_var(arr : Array[Double]) -> Double {
  if arr.length() < 2 {
    return 0.0
  }
  let m = dss_mean(arr)
  let mut sum_sq = 0.0
  let mut i = 0
  while i < arr.length() {
    sum_sq = sum_sq + (arr[i] - m).pow(2.0)
    i = i + 1
  }
  sum_sq / (arr.length() - 1).to_double()
}

// ===== Public API functions =====

/// Estimate per-gene dispersion using method of moments
pub fn dss_estimate_dispersion(counts_data : DSSCountData) -> DSSDispResult {
  let n_genes = counts_data.dss_n_genes()
  let n_samples = counts_data.dss_n_samples()
  
  let dispersions : Array[Double] = Array::new()
  let norm_counts = counts_data.dss_normalize()
  
  let mut i = 0
  while i < n_genes {
    let expr = norm_counts[i]
    let m = dss_mean(expr)
    let v = dss_var(expr)
    
    let disp = if m > 0.0 {
      v / m.pow(2.0)
    } else {
      0.0000000001
    }
    
    dispersions.push(disp)
    i = i + 1
  }
  
  // Shrinkage estimation
  let mean_disp = dss_mean(dispersions)
  let var_disp = dss_var(dispersions)
  
  let prior_df = if var_disp > 0.0 {
    2.0 * mean_disp.pow(2.0) / var_disp
  } else {
    10.0
  }
  let prior_var = mean_disp
  
  let shrunken : Array[Double] = Array::new()
  i = 0
  while i < n_genes {
    let disp = dispersions[i]
    // Empirical Bayes shrinkage
    let post_df = prior_df + n_samples.to_double()
    let post_var = (prior_df * prior_var + (n_samples.to_double() - 1.0) * disp) / post_df
    shrunken.push(post_var)
    i = i + 1
  }
  
  DSSDispResult::new(dispersions, shrunken, prior_df, prior_var)
}

/// Perform differential expression analysis using Wald test
pub fn dss_deseq(counts_data : DSSCountData, disp_result : DSSDispResult) -> DSSDEResult {
  let n_genes = counts_data.dss_n_genes()
  let n_samples = counts_data.dss_n_samples()
  let norm_counts = counts_data.dss_normalize()
  
  let group1_idx : Array[Int] = Array::new()
  let group2_idx : Array[Int] = Array::new()
  
  let mut j = 0
  while j < n_samples {
    if counts_data.groups[j] == "group1" {
      group1_idx.push(j)
    } else {
      group2_idx.push(j)
    }
    j = j + 1
  }
  
  let n1 = group1_idx.length().to_double()
  let n2 = group2_idx.length().to_double()
  
  let mean1_arr : Array[Double] = Array::new()
  let mean2_arr : Array[Double] = Array::new()
  let disp_arr : Array[Double] = Array::new()
  let wald_stat_arr : Array[Double] = Array::new()
  let p_values_arr : Array[Double] = Array::new()
  let gene_ids_arr : Array[String] = Array::new()
  
  let mut i = 0
  while i < n_genes {
    let expr = norm_counts[i]
    
    let sum1 = if n1 > 0.0 {
      let mut s = 0.0
      let mut k = 0
      while k < group1_idx.length() {
        s = s + expr[group1_idx[k]]
        k = k + 1
      }
      s
    } else {
      0.0
    }
    
    let sum2 = if n2 > 0.0 {
      let mut s = 0.0
      let mut k = 0
      while k < group2_idx.length() {
        s = s + expr[group2_idx[k]]
        k = k + 1
      }
      s
    } else {
      0.0
    }
    
    let m1 = if n1 > 0.0 { sum1 / n1 } else { 0.0 }
    let m2 = if n2 > 0.0 { sum2 / n2 } else { 0.0 }
    
    let disp = if i < disp_result.shrunken_dispersions.length() {
      disp_result.shrunken_dispersions[i]
    } else {
      0.1
    }
    
    // Wald statistic
    let se = (disp * (m1 / n1 + m2 / n2)).sqrt()
    let wald = if se > 0.0000000001 {
      (m2 - m1) / se
    } else {
      0.0
    }
    
    // Two-sided p-value approximation
    let p_val = if wald.abs() < 37.0 {
      2.0 * (1.0 - dss_normal_cdf(wald.abs()))
    } else {
      0.0
    }
    
    mean1_arr.push(m1)
    mean2_arr.push(m2)
    disp_arr.push(disp)
    wald_stat_arr.push(wald)
    p_values_arr.push(p_val)
    gene_ids_arr.push(counts_data.gene_ids[i])
    
    i = i + 1
  }
  
  // Multiple testing correction (Benjamini-Hochberg)
  let adj_p_values = dss_benjamini_hochberg(p_values_arr)
  
  // Determine significance (FDR < 0.05)
  let significant : Array[Bool] = Array::new()
  let mut k = 0
  while k < n_genes {
    significant.push(adj_p_values[k] < 0.05)
    k = k + 1
  }
  
  DSSDEResult::new(gene_ids_arr, mean1_arr, mean2_arr, disp_arr, wald_stat_arr,
                   p_values_arr, adj_p_values, significant)
}

/// Perform differential methylation analysis
pub fn dss_dml(methyl_data : DSSMethylData) -> DSSDMResult {
  let n_sites = methyl_data.dss_n_cpgs()
  let n_samples = methyl_data.dss_n_samples()
  
  let group1_idx : Array[Int] = Array::new()
  let group2_idx : Array[Int] = Array::new()
  
  let mut j = 0
  while j < n_samples {
    if methyl_data.groups[j] == "group1" {
      group1_idx.push(j)
    } else {
      group2_idx.push(j)
    }
    j = j + 1
  }
  
  let n1 = group1_idx.length().to_double()
  let n2 = group2_idx.length().to_double()
  
  let meth_mean1_arr : Array[Double] = Array::new()
  let meth_mean2_arr : Array[Double] = Array::new()
  let delta_meth_arr : Array[Double] = Array::new()
  let stats_arr : Array[Double] = Array::new()
  let p_values_arr : Array[Double] = Array::new()
  let cpg_ids_arr : Array[String] = Array::new()
  
  let mut i = 0
  while i < n_sites {
    let mut sum1 = 0.0
    let mut total1 = 0.0
    let mut sum2 = 0.0
    let mut total2 = 0.0
    
    let mut k = 0
    while k < group1_idx.length() {
      let idx = group1_idx[k]
      sum1 = sum1 + methyl_data.meth_counts[i][idx]
      total1 = total1 + methyl_data.total_counts[i][idx]
      k = k + 1
    }
    
    k = 0
    while k < group2_idx.length() {
      let idx = group2_idx[k]
      sum2 = sum2 + methyl_data.meth_counts[i][idx]
      total2 = total2 + methyl_data.total_counts[i][idx]
      k = k + 1
    }
    
    let m1 = if total1 > 0.0 { sum1 / total1 } else { 0.0 }
    let m2 = if total2 > 0.0 { sum2 / total2 } else { 0.0 }
    let delta = m2 - m1
    
    // Test statistic (approximation)
    let se = if total1 + total2 > 0.0 {
      (m1 * (1.0 - m1) / total1 + m2 * (1.0 - m2) / total2).sqrt()
    } else {
      1.0
    }
    
    let stat = if se > 0.0000000001 { delta / se } else { 0.0 }
    let p_val = if stat.abs() < 37.0 {
      2.0 * (1.0 - dss_normal_cdf(stat.abs()))
    } else {
      0.0
    }
    
    meth_mean1_arr.push(m1)
    meth_mean2_arr.push(m2)
    delta_meth_arr.push(delta)
    stats_arr.push(stat)
    p_values_arr.push(p_val)
    cpg_ids_arr.push(methyl_data.cpg_ids[i])
    
    i = i + 1
  }
  
  // Multiple testing correction
  let adj_p_values = dss_benjamini_hochberg(p_values_arr)
  
  let significant : Array[Bool] = Array::new()
  let mut s = 0
  while s < n_sites {
    significant.push(adj_p_values[s] < 0.05)
    s = s + 1
  }
  
  DSSDMResult::new(cpg_ids_arr, meth_mean1_arr, meth_mean2_arr, delta_meth_arr,
                   stats_arr, p_values_arr, adj_p_values, significant)
}

/// Detect differentially methylated regions
pub fn dss_dmr(methyl_data : DSSMethylData, dm_result : DSSDMResult, min_cpgs : Int, min_delta : Double) -> DSSDMRResult {
  let dmr_ids : Array[String] = Array::new()
  let chrom_arr : Array[String] = Array::new()
  let start_arr : Array[Int] = Array::new()
  let end_arr : Array[Int] = Array::new()
  let n_cpgs_arr : Array[Int] = Array::new()
  let mean_delta_arr : Array[Double] = Array::new()
  let stats_arr : Array[Double] = Array::new()
  let p_values_arr : Array[Double] = Array::new()
  
  let mut i = 0
  while i < dm_result.cpg_ids.length() {
    if dm_result.significant[i] && dm_result.delta_meth[i].abs() >= min_delta {
      // Start of potential DMR
      let current_chrom = if i < methyl_data.chrom.length() { methyl_data.chrom[i] } else { "chr1" }
      let current_pos = if i < methyl_data.pos.length() { methyl_data.pos[i] } else { i }
      
      let mut j = i + 1
      let mut cpgs = 1
      let mut total_delta = dm_result.delta_meth[i]
      let mut total_stat = dm_result.stats[i]
      let max_gap = 10000
      let mut last_pos = current_pos
      
      while j < dm_result.cpg_ids.length() && j < i + 100 {
        if dm_result.significant[j] && dm_result.delta_meth[j].abs() >= min_delta {
          let next_chrom = if j < methyl_data.chrom.length() { methyl_data.chrom[j] } else { current_chrom }
          let next_pos = if j < methyl_data.pos.length() { methyl_data.pos[j] } else { last_pos + 100 }
          
          if next_chrom == current_chrom && next_pos - last_pos <= max_gap {
            cpgs = cpgs + 1
            total_delta = total_delta + dm_result.delta_meth[j]
            total_stat = total_stat + dm_result.stats[j]
            last_pos = next_pos
            j = j + 1
          } else {
            break
          }
        } else {
          break
        }
      }
      
      if cpgs >= min_cpgs {
        let dmr_id = "DMR_" + (dmr_ids.length() + 1).to_string()
        dmr_ids.push(dmr_id)
        chrom_arr.push(current_chrom)
        start_arr.push(current_pos)
        end_arr.push(last_pos)
        n_cpgs_arr.push(cpgs)
        mean_delta_arr.push(total_delta / cpgs.to_double())
        stats_arr.push(total_stat)
        p_values_arr.push(0.0)
      }
      
      i = j
    } else {
      i = i + 1
    }
  }
  
  DSSDMRResult::new(dmr_ids, chrom_arr, start_arr, end_arr, n_cpgs_arr,
                    mean_delta_arr, stats_arr, p_values_arr)
}

/// Standard normal CDF approximation (renamed to avoid conflict with ballgown)
fn dss_normal_cdf(x : Double) -> Double {
  // Abramowitz and Stegun approximation
  let t = 1.0 / (1.0 + 0.2316419 * x.abs())
  let d = 0.3989422804014327 * @math.exp(-x * x / 2.0)
  let p = d * t * (0.3193815 + t * (-0.3565638 + t * (1.781478 + t * (-1.821256 + t * 1.330274))))
  if x > 0.0 { 1.0 - p } else { p }
}

/// Benjamini-Hochberg multiple testing correction (renamed to avoid conflict with ballgown)
fn dss_benjamini_hochberg(p_values : Array[Double]) -> Array[Double] {
  let n = p_values.length()
  let result = Array::make(n, 1.0)
  
  if n == 0 {
    return result
  }
  
  // Create index array and sort by p-value
  let indices : Array[Int] = Array::new()
  let sorted_pvals : Array[Double] = Array::new()
  
  let mut i = 0
  while i < n {
    indices.push(i)
    sorted_pvals.push(p_values[i])
    i = i + 1
  }
  
  // Simple bubble sort
  let mut j = 0
  while j < n {
    let mut k = j + 1
    while k < n {
      if sorted_pvals[k] < sorted_pvals[j] {
        let tmp_p = sorted_pvals[j]
        sorted_pvals[j] = sorted_pvals[k]
        sorted_pvals[k] = tmp_p
        
        let tmp_i = indices[j]
        indices[j] = indices[k]
        indices[k] = tmp_i
      }
      k = k + 1
    }
    j = j + 1
  }
  
  // Apply BH correction
  let adj = Array::make(n, 1.0)
  let mut last = 1.0
  
  i = n - 1
  while i >= 0 {
    let rank = (i + 1).to_double()
    let val = sorted_pvals[i] * n.to_double() / rank
    last = if val < last { val } else { last }
    adj[i] = if last <= 1.0 { last } else { 1.0 }
    i = i - 1
  }
  
  // Map back to original order
  i = 0
  while i < n {
    result[indices[i]] = adj[i]
    i = i + 1
  }
  
  result
}

/// Create example count data for testing
pub fn dss_create_example_counts() -> DSSCountData {
  let gene_ids = ["GeneA", "GeneB", "GeneC", "GeneD", "GeneE"]
  let counts = [
    [100.0, 120.0, 110.0, 200.0, 220.0, 210.0],
    [50.0, 48.0, 52.0, 55.0, 53.0, 51.0],
    [200.0, 180.0, 220.0, 600.0, 620.0, 580.0],
    [30.0, 32.0, 28.0, 29.0, 31.0, 30.0],
    [80.0, 85.0, 75.0, 82.0, 78.0, 84.0],
  ]
  let groups = ["group1", "group1", "group1", "group2", "group2", "group2"]
  let lib_sizes = [1.0, 1.0, 1.0, 1.0, 1.0, 1.0]
  
  DSSCountData::new(gene_ids, counts, groups, lib_sizes)
}

/// Create example methylation data for testing
pub fn dss_create_example_methylation() -> DSSMethylData {
  let cpg_ids = ["CpG1", "CpG2", "CpG3", "CpG4", "CpG5"]
  let chrom = ["chr1", "chr1", "chr1", "chr1", "chr1"]
  let pos = [1000, 2000, 3000, 4000, 5000]
  let meth_counts = [
    [18.0, 20.0, 19.0, 5.0, 4.0, 6.0],
    [15.0, 14.0, 16.0, 30.0, 32.0, 28.0],
    [8.0, 9.0, 7.0, 25.0, 24.0, 26.0],
    [45.0, 42.0, 48.0, 10.0, 12.0, 8.0],
    [22.0, 24.0, 20.0, 23.0, 21.0, 22.0],
  ]
  let total_counts = [
    [20.0, 22.0, 21.0, 20.0, 21.0, 19.0],
    [20.0, 19.0, 21.0, 35.0, 34.0, 36.0],
    [10.0, 11.0, 9.0, 30.0, 29.0, 31.0],
    [50.0, 48.0, 52.0, 15.0, 16.0, 14.0],
    [25.0, 26.0, 24.0, 25.0, 24.0, 26.0],
  ]
  let groups = ["group1", "group1", "group1", "group2", "group2", "group2"]
  
  DSSMethylData::new(cpg_ids, chrom, pos, meth_counts, total_counts, groups)
}
