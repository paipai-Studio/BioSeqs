/// Scanpy-inspired: pp.recipe_seurat, pp.recipe_zheng17 — standard preprocessing recipes

///|
/// Seurat-style preprocessing recipe:
/// 1. Normalize total counts to target_sum
/// 2. Log1p transform
/// 3. Select highly variable genes (by dispersion)
/// 4. Scale to unit variance (per gene)
///
/// Parameters:
/// - `x`: cells × genes count matrix
/// - `target_sum`: target total counts per cell (default 10000)
/// - `n_top_genes`: number of HVGs to select (default 2000)
pub fn scanpy_recipe_seurat(
  x : Array[Array[Double]],
  target_sum~ : Double = 10000.0,
  n_top_genes~ : Int32 = 2000,
) -> Array[Array[Double]] {
  let n_cells = x.length()
  if n_cells == 0 {
    return []
  }
  let n_genes = x[0].length()
  // Step 1: Normalize total
  let normalized = Array(capacity=n_cells)
  for i in 0..<n_cells {
    let mut total = 0.0
    for j in 0..<n_genes {
      total += x[i][j]
    }
    let row = Array(capacity=n_genes)
    if total > 0.0 {
      let scale = target_sum / total
      for j in 0..<n_genes {
        row.push(x[i][j] * scale)
      }
    } else {
      for _ in 0..<n_genes {
        row.push(0.0)
      }
    }
    normalized.push(row)
  }
  // Step 2: Log1p
  for i in 0..<n_cells {
    for j in 0..<n_genes {
      normalized[i][j] = @math.ln(1.0 + normalized[i][j])
    }
  }
  // Step 3: Select HVGs by dispersion
  let gene_means = Array(capacity=n_genes)
  let gene_vars = Array(capacity=n_genes)
  for j in 0..<n_genes {
    let mut s = 0.0
    for i in 0..<n_cells {
      s += normalized[i][j]
    }
    let mean = s / n_cells.to_double()
    gene_means.push(mean)
    let mut v = 0.0
    for i in 0..<n_cells {
      let diff = normalized[i][j] - mean
      v += diff * diff
    }
    gene_vars.push(v / n_cells.to_double())
  }
  // Dispersion = variance / mean (skip zero-mean genes)
  let dispersions = Array(capacity=n_genes)
  for j in 0..<n_genes {
    if gene_means[j] > 0.0 {
      dispersions.push((gene_vars[j], j))
    } else {
      dispersions.push((0.0, j))
    }
  }
  // Sort by dispersion descending
  dispersions.sort_by(fn(a, b) {
    if a.0 > b.0 {
      -1
    } else if a.0 < b.0 {
      1
    } else {
      0
    }
  })
  let n_select = if n_top_genes < n_genes { n_top_genes } else { n_genes }
  let hvg_set : Set[Int32] = Set([], capacity=n_select)
  for i in 0..<n_select {
    hvg_set.add(dispersions[i].1)
  }
  // Step 4: Scale to unit variance, keep only HVGs
  let result = Array(capacity=n_cells)
  for i in 0..<n_cells {
    let row = Array(capacity=n_select)
    for j in 0..<n_genes {
      if hvg_set.contains(j) {
        let std = gene_vars[j].sqrt()
        if std > 0.0 {
          row.push((normalized[i][j] - gene_means[j]) / std)
        } else {
          row.push(0.0)
        }
      }
    }
    result.push(row)
  }
  result
}

///|
/// Zheng 2017 recipe (smart-seq2):
/// 1. Filter genes expressed in < min_cells
/// 2. Normalize to target_sum
/// 3. Log1p
/// 4. Select HVGs by z-scored dispersion
/// 5. Scale to unit variance, zero mean
///
/// Parameters:
/// - `x`: cells × genes count matrix
/// - `target_sum`: target total (default 10000)
/// - `min_cells`: minimum cells expressing a gene (default 1)
/// - `n_top_genes`: number of HVGs (default 2000)
pub fn scanpy_recipe_zheng17(
  x : Array[Array[Double]],
  target_sum~ : Double = 10000.0,
  min_cells~ : Int32 = 1,
  n_top_genes~ : Int32 = 2000,
) -> Array[Array[Double]] {
  let n_cells = x.length()
  if n_cells == 0 {
    return []
  }
  let n_genes = x[0].length()
  // Step 1: Filter genes
  let gene_expressed = Array(capacity=n_genes)
  for j in 0..<n_genes {
    let mut cnt = 0
    for i in 0..<n_cells {
      if x[i][j] > 0.0 {
        cnt += 1
      }
    }
    gene_expressed.push(cnt >= min_cells)
  }
  // Collect kept gene indices
  let kept_genes = Array(capacity=n_genes)
  for j in 0..<n_genes {
    if gene_expressed[j] {
      kept_genes.push(j)
    }
  }
  let n_kept = kept_genes.length()
  if n_kept == 0 {
    return []
  }
  // Step 2: Normalize total
  let normalized = Array(capacity=n_cells)
  for i in 0..<n_cells {
    let mut total = 0.0
    for j in 0..<n_genes {
      total += x[i][j]
    }
    let row = Array(capacity=n_kept)
    if total > 0.0 {
      let scale = target_sum / total
      for k in 0..<n_kept {
        row.push(x[i][kept_genes[k]] * scale)
      }
    } else {
      for _ in 0..<n_kept {
        row.push(0.0)
      }
    }
    normalized.push(row)
  }
  // Step 3: Log1p
  for i in 0..<n_cells {
    for j in 0..<n_kept {
      normalized[i][j] = @math.ln(1.0 + normalized[i][j])
    }
  }
  // Step 4: Compute means and dispersions
  let means = Array(capacity=n_kept)
  let vars = Array(capacity=n_kept)
  for j in 0..<n_kept {
    let mut s = 0.0
    for i in 0..<n_cells {
      s += normalized[i][j]
    }
    let mean = s / n_cells.to_double()
    means.push(mean)
    let mut v = 0.0
    for i in 0..<n_cells {
      let diff = normalized[i][j] - mean
      v += diff * diff
    }
    vars.push(v / n_cells.to_double())
  }
  // Dispersion = var / mean
  let disps = Array(capacity=n_kept)
  for j in 0..<n_kept {
    if means[j] > 0.0 {
      disps.push(vars[j] / means[j])
    } else {
      disps.push(0.0)
    }
  }
  // Z-score dispersions
  let mut disp_mean = 0.0
  for d in disps {
    disp_mean += d
  }
  disp_mean /= n_kept.to_double()
  let mut disp_var = 0.0
  for d in disps {
    let diff = d - disp_mean
    disp_var += diff * diff
  }
  disp_var /= n_kept.to_double()
  let disp_std = disp_var.sqrt()
  let z_disps = Array(capacity=n_kept)
  for j in 0..<n_kept {
    if disp_std > 0.0 {
      z_disps.push(((disps[j] - disp_mean) / disp_std, j))
    } else {
      z_disps.push((0.0, j))
    }
  }
  // Sort by z-dispersion descending
  z_disps.sort_by(fn(a, b) {
    if a.0 > b.0 {
      -1
    } else if a.0 < b.0 {
      1
    } else {
      0
    }
  })
  let n_select = if n_top_genes < n_kept { n_top_genes } else { n_kept }
  // Step 5: Select HVGs and scale
  let result = Array(capacity=n_cells)
  for i in 0..<n_cells {
    let row = Array(capacity=n_select)
    for k in 0..<n_select {
      let j = z_disps[k].1
      let std = vars[j].sqrt()
      if std > 0.0 {
        row.push((normalized[i][j] - means[j]) / std)
      } else {
        row.push(0.0)
      }
    }
    result.push(row)
  }
  result
}
