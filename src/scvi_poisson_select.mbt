// scvi_poisson_select.mbt - Poisson gene selection
// Inspired by scvi-tools poisson_gene_selection
// Feature selection via Poisson deviance, robust to dropout

///|
/// Compute Poisson deviance for each gene
pub fn scvi_poisson_deviance(
  x : Array[Array[Double]],
  n_top_genes : Int32,
) -> (Array[Int32], Array[Double]) {
  let n_cells = x.length()
  let n_genes = x[0].length()
  // Compute row sums (library sizes) and column sums
  let lib_sizes = Array(capacity=n_cells)
    |> a => {
      for i in 0..<n_cells {
        let mut s = 0.0
        for g in 0..<n_genes {
          s += x[i][g]
        }
        a.push(s)
      }
      a
    }
  let gene_sums = Array(capacity=n_genes)
    |> a => {
      for g in 0..<n_genes {
        let mut s = 0.0
        for i in 0..<n_cells {
          s += x[i][g]
        }
        a.push(s)
      }
      a
    }
  let total = scvi_array_sum(gene_sums)
  // Compute Poisson deviance per gene
  let deviances = Array(capacity=n_genes)
  for g in 0..<n_genes {
    let mut deviance = 0.0
    for i in 0..<n_cells {
      let y = x[i][g]
      if y > 0.0 {
        // Expected under independence: mu = lib_size * gene_sum / total
        let mu = lib_sizes[i] * gene_sums[g] / (total + 0.0000000001)
        if mu > 0.0 {
          deviance += 2.0 * (y * @math.ln(y / mu) - (y - mu))
        }
      }
    }
    deviances.push(deviance)
  }
  // Rank genes by deviance (higher = more deviant = more informative)
  let indexed = Array(capacity=n_genes)
    |> a => {
      for g in 0..<n_genes {
        a.push((g, deviances[g]))
      }
      a
    }
  indexed.sort_by(fn(a, b) {
    if a.1 > b.1 {
      -1
    } else if a.1 < b.1 {
      1
    } else {
      0
    }
  })
  let selected = Array(capacity=n_top_genes)
  let selected_deviances = Array(capacity=n_top_genes)
  for i in 0..<n_top_genes {
    if i < indexed.length() {
      selected.push(indexed[i].0)
      selected_deviances.push(indexed[i].1)
    }
  }
  (selected, selected_deviances)
}

///|
/// Filter genes by Poisson deviance threshold
pub fn scvi_poisson_filter(
  x : Array[Array[Double]],
  threshold : Double,
) -> Array[Bool] {
  let n_genes = x[0].length()
  let (_, deviances) = scvi_poisson_deviance(x, n_genes)
  let result = Array(capacity=n_genes)
  for g in 0..<n_genes {
    result.push(deviances[g] > threshold)
  }
  result
}

///|
/// Compute gene selection statistics
pub fn scvi_gene_selection_stats(
  x : Array[Array[Double]],
) -> (Array[Double], Array[Double], Array[Double]) {
  let n_cells = x.length()
  let n_genes = x[0].length()
  let means = Array(capacity=n_genes)
  let variances = Array(capacity=n_genes)
  let zero_fracs = Array(capacity=n_genes)
  for g in 0..<n_genes {
    let mut sum = 0.0
    let mut sum_sq = 0.0
    let mut n_zero = 0
    for i in 0..<n_cells {
      let val = x[i][g]
      sum += val
      sum_sq += val * val
      if val < 0.5 {
        n_zero += 1
      }
    }
    let mean = sum / n_cells.to_double()
    let var_ = sum_sq / n_cells.to_double() - mean * mean
    means.push(mean)
    variances.push(var_)
    zero_fracs.push(n_zero.to_double() / n_cells.to_double())
  }
  (means, variances, zero_fracs)
}
