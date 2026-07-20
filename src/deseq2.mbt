///|
/// DESeq2 - Differential Expression Analysis for RNA-seq data.
/// Inspired by Bioconductor DESeq2 package.

///|
pub struct DESeqDataSet {
  counts : Array[Array[Int]]
  row_names : Array[String]
  col_names : Array[String]
  design : Array[Array[Double]]
  size_factors : Array[Double]
  dispersions : Array[Double]
}

///|
pub fn DESeqDataSet::new(
  counts : Array[Array[Int]],
  row_names : Array[String],
  col_names : Array[String],
  design : Array[Array[Double]],
) -> DESeqDataSet {
  let n_genes = counts.length()
  let n_samples = if counts.length() > 0 { counts[0].length() } else { 0 }
  
  DESeqDataSet::{
    counts,
    row_names,
    col_names,
    design,
    size_factors: Array::make(n_samples, 1.0),
    dispersions: Array::make(n_genes, 0.1),
  }
}

///|
pub struct DESeqResults {
  row_names : Array[String]
  base_mean : Array[Double]
  log2_fold_change : Array[Double]
  lfc_se : Array[Double]
  stat : Array[Double]
  p_value : Array[Double]
  padj : Array[Double]
}

///|
pub fn DESeqResults::new(
  row_names : Array[String],
  base_mean : Array[Double],
  log2_fold_change : Array[Double],
  lfc_se : Array[Double],
  stat : Array[Double],
  p_value : Array[Double],
  padj : Array[Double],
) -> DESeqResults {
  DESeqResults::{
    row_names,
    base_mean,
    log2_fold_change,
    lfc_se,
    stat,
    p_value,
    padj,
  }
}

///|
fn geometric_mean(x : Array[Double]) -> Double {
  let mut sum_log = 0.0
  let mut n_valid = 0
  for val in x {
    if val > 0.0 {
      sum_log = sum_log + log_nat(val)
      n_valid = n_valid + 1
    }
  }
  if n_valid > 0 {
    exp(sum_log / n_valid.to_double())
  } else {
    0.0
  }
}

///|
pub fn estimate_size_factors(dds : DESeqDataSet) -> DESeqDataSet {
  let n_genes = dds.counts.length()
  let n_samples = dds.col_names.length()
  
  let geo_means : Array[Double] = Array::new()
  for i = 0; i < n_genes; i = i + 1 {
    let gene_counts : Array[Double] = Array::new()
    for j = 0; j < n_samples; j = j + 1 {
      gene_counts.push(dds.counts[i][j].to_double())
    }
    geo_means.push(geometric_mean(gene_counts))
  }
  
  let size_factors : Array[Double] = Array::new()
  for j = 0; j < n_samples; j = j + 1 {
    let ratios : Array[Double] = Array::new()
    for i = 0; i < n_genes; i = i + 1 {
      if geo_means[i] > 0.0 && dds.counts[i][j] > 0 {
        ratios.push(dds.counts[i][j].to_double() / geo_means[i])
      }
    }
    
    let n_ratios = ratios.length()
    let sorted_ratios : Array[Double] = Array::new()
    for r in ratios {
      sorted_ratios.push(r)
    }
    for i = 0; i < n_ratios - 1; i = i + 1 {
      for k = 0; k < n_ratios - 1 - i; k = k + 1 {
        if sorted_ratios[k] > sorted_ratios[k + 1] {
          let tmp = sorted_ratios[k]
          sorted_ratios[k] = sorted_ratios[k + 1]
          sorted_ratios[k + 1] = tmp
        }
      }
    }
    
    if n_ratios > 0 {
      let median = if n_ratios % 2 == 0 {
        (sorted_ratios[n_ratios / 2 - 1] + sorted_ratios[n_ratios / 2]) / 2.0
      } else {
        sorted_ratios[n_ratios / 2]
      }
      size_factors.push(median)
    } else {
      size_factors.push(1.0)
    }
  }
  
  let mut mean_sf = 0.0
  for sf in size_factors {
    mean_sf = mean_sf + sf
  }
  mean_sf = mean_sf / n_samples.to_double()
  
  let norm_sf : Array[Double] = Array::new()
  for sf in size_factors {
    norm_sf.push(sf / mean_sf)
  }
  
  DESeqDataSet::{
    counts: dds.counts,
    row_names: dds.row_names,
    col_names: dds.col_names,
    design: dds.design,
    size_factors: norm_sf,
    dispersions: dds.dispersions,
  }
}

///|
pub fn normalize_counts(dds : DESeqDataSet) -> Array[Array[Double]] {
  let n_genes = dds.counts.length()
  let n_samples = dds.col_names.length()
  
  let normalized : Array[Array[Double]] = Array::new()
  for i = 0; i < n_genes; i = i + 1 {
    let row : Array[Double] = Array::new()
    for j = 0; j < n_samples; j = j + 1 {
      if dds.size_factors[j] > 0.0 {
        row.push(dds.counts[i][j].to_double() / dds.size_factors[j])
      } else {
        row.push(dds.counts[i][j].to_double())
      }
    }
    normalized.push(row)
  }
  normalized
}

///|
pub fn log2_cpm(dds : DESeqDataSet) -> Array[Array[Double]] {
  let n_genes = dds.counts.length()
  let n_samples = dds.col_names.length()
  
  let normalized = normalize_counts(dds)
  
  let lib_sizes : Array[Double] = Array::new()
  for j = 0; j < n_samples; j = j + 1 {
    let mut sum = 0.0
    for i = 0; i < n_genes; i = i + 1 {
      sum = sum + normalized[i][j]
    }
    lib_sizes.push(sum)
  }
  
  let log2_cpm : Array[Array[Double]] = Array::new()
  for i = 0; i < n_genes; i = i + 1 {
    let row : Array[Double] = Array::new()
    for j = 0; j < n_samples; j = j + 1 {
      let cpm = if lib_sizes[j] > 0.0 {
        normalized[i][j] * 1000000.0 / lib_sizes[j]
      } else {
        0.0
      }
      row.push(if cpm > 0.0 { log_nat(cpm) / log_nat(2.0) } else { 0.0 })
    }
    log2_cpm.push(row)
  }
  log2_cpm
}

///|
pub fn estimate_dispersions(dds : DESeqDataSet) -> DESeqDataSet {
  let n_genes = dds.counts.length()
  let n_samples = dds.col_names.length()
  
  let normalized = normalize_counts(dds)
  
  let log_means : Array[Double] = Array::new()
  for i = 0; i < n_genes; i = i + 1 {
    let mut sum = 0.0
    for j = 0; j < n_samples; j = j + 1 {
      sum = sum + normalized[i][j]
    }
    log_means.push(if sum > 0.0 { log_nat(sum / n_samples.to_double()) } else { -10.0 })
  }
  
  let dispersions : Array[Double] = Array::new()
  for i = 0; i < n_genes; i = i + 1 {
    let mut sum_counts = 0.0
    for j = 0; j < n_samples; j = j + 1 {
      sum_counts = sum_counts + dds.counts[i][j].to_double()
    }
    
    if sum_counts > 0.0 {
      let mu_hat = sum_counts / n_samples.to_double()
      let mut var_obs = 0.0
      for j = 0; j < n_samples; j = j + 1 {
        let diff = dds.counts[i][j].to_double() - mu_hat
        var_obs = var_obs + diff * diff
      }
      var_obs = var_obs / (n_samples - 1).to_double()
      
      if mu_hat > 0.0 {
        let disp_est = (var_obs / mu_hat - 1.0) / mu_hat
        dispersions.push(if disp_est > 0.001 { disp_est } else { 0.001 })
      } else {
        dispersions.push(0.1)
      }
    } else {
      dispersions.push(0.1)
    }
  }
  
  let alpha = 0.1
  let beta = 0.5
  
  let final_dispersions : Array[Double] = Array::new()
  for i = 0; i < n_genes; i = i + 1 {
    let raw_disp = dispersions[i]
    let prior_disp = alpha * exp(-beta * log_means[i])
    let weight = 1.0 / (1.0 / n_samples.to_double() + 1.0)
    final_dispersions.push((1.0 - weight) * prior_disp + weight * raw_disp)
  }
  
  DESeqDataSet::{
    counts: dds.counts,
    row_names: dds.row_names,
    col_names: dds.col_names,
    design: dds.design,
    size_factors: dds.size_factors,
    dispersions: final_dispersions,
  }
}

///|
fn pnorm(x : Double) -> Double {
  let t = 1.0 / (1.0 + 0.2316419 * x.abs())
  let d = 0.3989423 * exp(-x * x / 2.0)
  let prob = d * t * (0.3193815 + t * (-0.3565638 + t * (1.781478 + t * (-1.821256 + t * 1.330274))))
  if x >= 0.0 {
    1.0 - prob
  } else {
    prob
  }
}

///|
pub fn deseq(dds : DESeqDataSet) -> DESeqResults {
  let n_genes = dds.counts.length()
  let n_samples = dds.col_names.length()
  let n_coefs = if dds.design.length() > 0 { dds.design[0].length() } else { 0 }
  
  let base_mean : Array[Double] = Array::new()
  let log2_fold_change : Array[Double] = Array::new()
  let lfc_se : Array[Double] = Array::new()
  let stat : Array[Double] = Array::new()
  let p_value : Array[Double] = Array::new()
  
  let normalized = normalize_counts(dds)
  
  let intercept : Array[Double] = Array::new()
  let slope : Array[Double] = Array::new()
  let slope_se : Array[Double] = Array::new()
  
  for i = 0; i < n_genes; i = i + 1 {
    let mut sum_x = 0.0
    let mut sum_y = 0.0
    let mut sum_xy = 0.0
    let mut sum_x2 = 0.0
    let mut n = 0.0
    
    for j = 0; j < n_samples; j = j + 1 {
      let x = dds.design[j][1]
      let y = if normalized[i][j] > 0.0 { log_nat(normalized[i][j]) } else { 0.0 }
      sum_x = sum_x + x
      sum_y = sum_y + y
      sum_xy = sum_xy + x * y
      sum_x2 = sum_x2 + x * x
      n = n + 1.0
    }
    
    let denom = n * sum_x2 - sum_x * sum_x
    let b1 = if denom != 0.0 { (n * sum_xy - sum_x * sum_y) / denom } else { 0.0 }
    let b0 = (sum_y - b1 * sum_x) / n
    
    let mut sse = 0.0
    for j = 0; j < n_samples; j = j + 1 {
      let x = dds.design[j][1]
      let y = if normalized[i][j] > 0.0 { log_nat(normalized[i][j]) } else { 0.0 }
      let pred = b0 + b1 * x
      sse = sse + (y - pred) * (y - pred)
    }
    
    let se_b1 = if denom != 0.0 {
      ((sse / (n - 2.0)) / denom).sqrt()
    } else {
      1.0
    }
    
    intercept.push(b0)
    slope.push(b1)
    slope_se.push(se_b1)
    
    let mut sum_fitted = 0.0
    for j = 0; j < n_samples; j = j + 1 {
      let x = dds.design[j][1]
      sum_fitted = sum_fitted + exp(b0 + b1 * x)
    }
    base_mean.push(sum_fitted / n)
    
    log2_fold_change.push(b1 * log_nat(2.0))
    lfc_se.push(se_b1 * log_nat(2.0))
    
    let t_stat = if se_b1 > 0.0 { b1 / se_b1 } else { 0.0 }
    stat.push(t_stat)
    
    p_value.push(2.0 * (1.0 - pnorm(t_stat.abs())))
  }
  
  let padj = deseq2_adjust_p_values(p_value)
  
  DESeqResults::new(
    dds.row_names,
    base_mean,
    log2_fold_change,
    lfc_se,
    stat,
    p_value,
    padj,
  )
}

///|
pub fn results(dds : DESeqDataSet, contrast? : Int = -1) -> DESeqResults {
  deseq(dds)
}

///|
fn deseq2_adjust_p_values(p_values : Array[Double]) -> Array[Double] {
  let n = p_values.length()
  if n == 0 {
    return []
  }
  
  let indices : Array[Int] = Array::new()
  for i = 0; i < n; i = i + 1 {
    indices.push(i)
  }
  
  let sorted_indices = Array::new()
  for idx in indices {
    sorted_indices.push(idx)
  }
  for i = 0; i < n - 1; i = i + 1 {
    for j = 0; j < n - 1 - i; j = j + 1 {
      if p_values[sorted_indices[j]] > p_values[sorted_indices[j + 1]] {
        let tmp = sorted_indices[j]
        sorted_indices[j] = sorted_indices[j + 1]
        sorted_indices[j + 1] = tmp
      }
    }
  }
  
  let padj : Array[Double] = Array::make(n, 1.0)
  
  let mut min_padj = 1.0
  let mut i = n - 1
  while i >= 0 {
    let rank = (i + 1).to_double()
    let adj = p_values[sorted_indices[i]] * n.to_double() / rank
    if adj < min_padj {
      min_padj = adj
    }
    padj[sorted_indices[i]] = min_padj
    i = i - 1
  }
  
  let final_padj : Array[Double] = Array::new()
  for p in padj {
    final_padj.push(if p > 1.0 { 1.0 } else { p })
  }
  final_padj
}

///|
pub fn lfc_shrink(dds : DESeqDataSet, res : DESeqResults, coef? : Int = -1) -> DESeqResults {
  let n_genes = dds.counts.length()
  
  let shrunk_lfc : Array[Double] = Array::new()
  let shrunk_se : Array[Double] = Array::new()
  
  let prior_var = 0.5
  
  for i = 0; i < n_genes; i = i + 1 {
    let raw_lfc = res.log2_fold_change[i]
    let raw_se = res.lfc_se[i]
    
    let se2 = raw_se * raw_se
    let weight = se2 / (se2 + prior_var)
    
    shrunk_lfc.push(weight * raw_lfc)
    shrunk_se.push((se2 * prior_var / (se2 + prior_var)).sqrt())
  }
  
  DESeqResults::new(
    res.row_names,
    res.base_mean,
    shrunk_lfc,
    shrunk_se,
    res.stat,
    res.p_value,
    res.padj,
  )
}

///|
pub fn significant_genes(
  results : DESeqResults,
  alpha? : Double = 0.05,
  lfc_threshold? : Double = 0.0,
) -> Array[String] {
  let genes : Array[String] = Array::new()
  for i = 0; i < results.row_names.length(); i = i + 1 {
    if results.padj[i] < alpha && results.log2_fold_change[i].abs() > lfc_threshold {
      genes.push(results.row_names[i])
    }
  }
  genes
}

///|
pub fn top_genes(results : DESeqResults, n? : Int = 10) -> Array[(String, Double, Double)] {
  let n_genes = results.row_names.length()
  let actual_n = if n < n_genes { n } else { n_genes }
  
  let indices : Array[Int] = Array::new()
  for i = 0; i < n_genes; i = i + 1 {
    indices.push(i)
  }
  
  let sorted_indices = Array::new()
  for idx in indices {
    sorted_indices.push(idx)
  }
  for i = 0; i < n_genes - 1; i = i + 1 {
    for j = 0; j < n_genes - 1 - i; j = j + 1 {
      if results.padj[sorted_indices[j]] > results.padj[sorted_indices[j + 1]] {
        let tmp = sorted_indices[j]
        sorted_indices[j] = sorted_indices[j + 1]
        sorted_indices[j + 1] = tmp
      }
    }
  }
  
  let top : Array[(String, Double, Double)] = Array::new()
  for i = 0; i < actual_n; i = i + 1 {
    let idx = sorted_indices[i]
    top.push((results.row_names[idx], results.log2_fold_change[idx], results.padj[idx]))
  }
  top
}

///|
pub fn sample_deseq_dataset() -> DESeqDataSet {
  let genes = [
    "GeneA", "GeneB", "GeneC", "GeneD", "GeneE", "GeneF", "GeneG", "GeneH", "GeneI", "GeneJ",
    "GeneK", "GeneL", "GeneM", "GeneN", "GeneO", "GeneP", "GeneQ", "GeneR", "GeneS", "GeneT",
  ]
  let samples = ["Ctrl1", "Ctrl2", "Ctrl3", "Treat1", "Treat2", "Treat3"]
  
  let counts : Array[Array[Int]] = Array::new()
  for i = 0; i < 20; i = i + 1 {
    let row : Array[Int] = Array::new()
    let base = 100.0 + i.to_double() * 50.0
    
    for j = 0; j < 3; j = j + 1 {
      let noise = ((i * 7 + j * 13) % 20 - 10).to_double() * 0.1
      row.push((base * (1.0 + noise)).to_int())
    }
    
    for j = 3; j < 6; j = j + 1 {
      let fold_change = if i < 5 { 2.5 } else if i >= 15 { 0.4 } else { 1.0 }
      let noise = ((i * 7 + j * 13) % 20 - 10).to_double() * 0.1
      row.push((base * fold_change * (1.0 + noise)).to_int())
    }
    
    counts.push(row)
  }
  
  let design : Array[Array[Double]] = [
    [1.0, 0.0], [1.0, 0.0], [1.0, 0.0],
    [1.0, 1.0], [1.0, 1.0], [1.0, 1.0],
  ]
  
  DESeqDataSet::new(counts, genes, samples, design)
}
