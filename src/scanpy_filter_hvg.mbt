/// Scanpy-inspired: pp.filter_highly_variable_genes — filter to only HVGs

///|
/// Filter a matrix to keep only highly variable genes.
///
/// Parameters:
/// - `x`: cells × genes matrix
/// - `var_names`: gene names
/// - `hvg_mask`: boolean mask (true = keep)
pub fn scanpy_filter_hvg(
  x : Array[Array[Double]],
  var_names : Array[String],
  hvg_mask : Array[Bool],
) -> (Array[Array[Double]], Array[String]) {
  let n_cells = x.length()
  if n_cells == 0 {
    return ([], [])
  }
  let n_genes = x[0].length()
  // Collect indices to keep
  let keep = Array(capacity=n_genes)
  for j in 0..<n_genes {
    if j < hvg_mask.length() && hvg_mask[j] {
      keep.push(j)
    }
  }
  if keep.length() == 0 {
    return ([], [])
  }
  // Filter matrix
  let result = Array(capacity=n_cells)
  for i in 0..<n_cells {
    let row = Array(capacity=keep.length())
    for k in 0..<keep.length() {
      row.push(x[i][keep[k]])
    }
    result.push(row)
  }
  // Filter var names
  let new_names = Array(capacity=keep.length())
  for k in 0..<keep.length() {
    if keep[k] < var_names.length() {
      new_names.push(var_names[keep[k]])
    } else {
      new_names.push("")
    }
  }
  (result, new_names)
}

///|
/// Compute highly variable genes using dispersion-based method.
/// Returns a boolean mask and dispersion scores.
///
/// Parameters:
/// - `x`: cells × genes log-normalized matrix
/// - `n_top_genes`: number of HVGs to select (default 2000)
/// - `flavor`: "seurat" or "cell_ranger" (default "seurat")
pub fn scanpy_hvg_dispersion(
  x : Array[Array[Double]],
  n_top_genes~ : Int32 = 2000,
  flavor~ : String = "seurat",
) -> (Array[Bool], Array[Double]) {
  let n_cells = x.length()
  let n_genes = if n_cells > 0 { x[0].length() } else { 0 }
  if n_genes == 0 {
    return ([], [])
  }
  // Compute means and variances
  let means = Array(capacity=n_genes)
  let vars = Array(capacity=n_genes)
  for j in 0..<n_genes {
    let mut s = 0.0
    for i in 0..<n_cells {
      s += x[i][j]
    }
    let mean = s / n_cells.to_double()
    means.push(mean)
    let mut v = 0.0
    for i in 0..<n_cells {
      let diff = x[i][j] - mean
      v += diff * diff
    }
    if flavor == "cell_ranger" {
      // Cell Ranger: use binomial variance
      vars.push(v / n_cells.to_double())
      // Seurat: use sample variance
    } else if n_cells > 1 {
      vars.push(v / (n_cells - 1).to_double())
    } else {
      vars.push(0.0)
    }
  }
  // Compute normalized dispersion
  let disps = Array(capacity=n_genes)
  for j in 0..<n_genes {
    if means[j] > 0.0 && vars[j] > 0.0 {
      disps.push(vars[j] / means[j])
    } else {
      disps.push(0.0)
    }
  }
  // Log-transform means and disps for binning
  let log_means = Array(capacity=n_genes)
  let log_disps = Array(capacity=n_genes)
  for j in 0..<n_genes {
    log_means.push(@math.ln(means[j] + 1.0))
    log_disps.push(@math.ln(disps[j] + 1.0))
  }
  // Bin by log_mean and z-score log_disp within each bin
  let n_bins = 20
  let mut min_lm = 1.0e100
  let mut max_lm = -1.0e100
  for v in log_means {
    if v < min_lm {
      min_lm = v
    }
    if v > max_lm {
      max_lm = v
    }
  }
  let bin_width = if max_lm > min_lm {
    (max_lm - min_lm) / n_bins.to_double()
  } else {
    1.0
  }
  // Compute z-scored dispersion per bin
  let z_disps = Array(capacity=n_genes)
  for b in 0..<n_bins {
    let lo = min_lm + b.to_double() * bin_width
    let hi = lo + bin_width
    // Collect genes in this bin
    let bin_genes = Array(capacity=0)
    for j in 0..<n_genes {
      if log_means[j] >= lo && (log_means[j] < hi || b == n_bins - 1) {
        bin_genes.push(j)
      }
    }
    if bin_genes.length() == 0 {
      continue
    }
    // Mean and std of log_disp in bin
    let mut s = 0.0
    for j in bin_genes {
      s += log_disps[j]
    }
    let bin_mean = s / bin_genes.length().to_double()
    let mut v = 0.0
    for j in bin_genes {
      let diff = log_disps[j] - bin_mean
      v += diff * diff
    }
    let bin_std = (v / bin_genes.length().to_double()).sqrt()
    for j in bin_genes {
      if bin_std > 0.0 {
        z_disps.push(((log_disps[j] - bin_mean) / bin_std, j))
      } else {
        z_disps.push((0.0, j))
      }
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
  // Select top N
  let mask = Array(capacity=n_genes)
  for _ in 0..<n_genes {
    mask.push(false)
  }
  let n_select = if n_top_genes < z_disps.length() {
    n_top_genes
  } else {
    z_disps.length()
  }
  for k in 0..<n_select {
    mask[z_disps[k].1] = true
  }
  (mask, disps)
}
