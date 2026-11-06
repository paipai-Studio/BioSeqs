/// Scanpy second-round: Experimental preprocessing extensions
/// Inspired by scanpy.experimental: HVG pearson residuals, VST recipe

///|
/// Select highly variable genes using Pearson residuals.
/// Genes with high variance in Pearson residuals are considered HVG.
pub fn scanpy_highly_variable_genes_pearson_residuals(
  x : Array[Array[Double]],
  n_top_genes~ : Int32 = 2000,
  min_mean~ : Double = 0.0,
  max_mean~ : Double = 100.0,
  theta~ : Double = 100.0,
) -> Array[(Int32, Bool, Double)] {
  let n_cells = x.length()
  let n_genes = if n_cells > 0 { x[0].length() } else { 0 }
  // Compute gene means and variances of Pearson residuals
  let stats = Array(capacity=n_genes)
  for g in 0..<n_genes {
    let mut sum = 0.0
    for i in 0..<n_cells {
      sum = sum + x[i][g]
    }
    let mean = if n_cells > 0 { sum / n_cells.to_double() } else { 0.0 }
    // Pearson residual: (x - mean) / sqrt(mean + mean^2/theta)
    let mut var_resid = 0.0
    for i in 0..<n_cells {
      let denom = @math.pow(mean + mean * mean / theta, 0.5)
      let resid = if denom > 0.0 { (x[i][g] - mean) / denom } else { 0.0 }
      var_resid = var_resid + resid * resid
    }
    if n_cells > 0 {
      var_resid = var_resid / n_cells.to_double()
    }
    stats.push((g, mean, var_resid))
  }
  // Filter by mean range
  let filtered = Array(capacity=n_genes)
  for s in stats {
    if s.1 >= min_mean && s.1 <= max_mean {
      filtered.push(s)
    }
  }
  // Sort by variance descending, select top n
  filtered.sort_by(fn(a, b) {
    if a.2 > b.2 {
      -1
    } else if a.2 < b.2 {
      1
    } else {
      0
    }
  })
  let hvg_set = Set([])
  let k = if n_top_genes < filtered.length() {
    n_top_genes
  } else {
    filtered.length()
  }
  for i in 0..<k {
    hvg_set.add(filtered[i].0)
  }
  // Build result: (gene_index, is_hvg, variance)
  let result = Array(capacity=n_genes)
  for s in stats {
    result.push((s.0, hvg_set.contains(s.0), s.2))
  }
  result
}

///|
/// Variance Stabilizing Transform (VST) recipe.
/// Applies a VST based on the negative binomial parameters.
pub fn scanpy_pp_recipe_vst(
  x : Array[Array[Double]],
  theta~ : Double = 100.0,
  min_cells~ : Int32 = 5,
  n_top_genes~ : Int32 = 2000,
) -> Array[Array[Double]] {
  let n_cells = x.length()
  let n_genes = if n_cells > 0 { x[0].length() } else { 0 }
  // Step 1: Filter genes expressed in < min_cells
  let gene_counts = Array(capacity=n_genes)
  for g in 0..<n_genes {
    let mut count = 0
    for i in 0..<n_cells {
      if x[i][g] > 0.0 {
        count = count + 1
      }
    }
    gene_counts.push(count)
  }
  // Step 2: Compute gene means
  let gene_means = Array(capacity=n_genes)
  for g in 0..<n_genes {
    let mut sum = 0.0
    for i in 0..<n_cells {
      sum = sum + x[i][g]
    }
    gene_means.push(if n_cells > 0 { sum / n_cells.to_double() } else { 0.0 })
  }
  // Step 3: Apply, VST: asinh(sqrt(x)) - asinh(sqrt(mean))
  // This stabilizes variance for NB distribution
  let result = Array(capacity=n_cells)
  for i in 0..<n_cells {
    let row = Array(capacity=n_genes)
    for g in 0..<n_genes {
      if gene_counts[g] >= min_cells && gene_means[g] > 0.0 {
        let v = x[i][g]
        let transformed = @math.ln(v + @math.pow(v * v + theta * theta, 0.5)) -
          @math.ln(
            gene_means[g] +
            @math.pow(gene_means[g] * gene_means[g] + theta * theta, 0.5),
          )
        row.push(transformed)
      } else {
        row.push(0.0)
      }
    }
    result.push(row)
  }
  ignore(n_top_genes)
  result
}

///|
/// Normalize using Pearson residuals with theta and clip parameters.
/// resid = (x - mean) / sqrt(mean + mean^2/theta), clipped to [-clip, clip]
pub fn scanpy_normalize_pearson_residuals_ext(
  x : Array[Array[Double]],
  theta~ : Double = 100.0,
  clip~ : Double = 30.0,
) -> Array[Array[Double]] {
  let n_cells = x.length()
  let n_genes = if n_cells > 0 { x[0].length() } else { 0 }
  // Compute gene means
  let gene_means = Array(capacity=n_genes)
  for g in 0..<n_genes {
    let mut sum = 0.0
    for i in 0..<n_cells {
      sum = sum + x[i][g]
    }
    gene_means.push(if n_cells > 0 { sum / n_cells.to_double() } else { 0.0 })
  }
  // Apply Pearson residual normalization
  let result = Array(capacity=n_cells)
  for i in 0..<n_cells {
    let row = Array(capacity=n_genes)
    for g in 0..<n_genes {
      let mean = gene_means[g]
      if mean > 0.0 {
        let denom = @math.pow(mean + mean * mean / theta, 0.5)
        let resid = (x[i][g] - mean) / denom
        // Clip to [-clip, clip]
        let clipped = if resid > clip {
          clip
        } else if resid < -clip {
          -clip
        } else {
          resid
        }
        row.push(clipped)
      } else {
        row.push(0.0)
      }
    }
    result.push(row)
  }
  result
}
