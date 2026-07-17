///|
/// SVA - Surrogate Variable Analysis.
/// Inspired by Bioconductor sva package.
/// Implements ComBat batch effect correction and SVA for hidden variation.

///|
/// Result of ComBat batch effect correction.
pub struct ComBatResult {
  adjusted : Array[Array[Double]]
  gamma_hat : Array[Array[Double]]
  delta_hat : Array[Array[Double]]
}

///|
/// Create a new ComBat result.
pub fn ComBatResult::new(
  adjusted : Array[Array[Double]],
  gamma_hat : Array[Array[Double]],
  delta_hat : Array[Array[Double]],
) -> ComBatResult {
  ComBatResult::{ adjusted, gamma_hat, delta_hat }
}

///|
/// Calculate column means of a matrix.
fn col_means(mat : Array[Array[Double]]) -> Array[Double] {
  let nrow = mat.length()
  if nrow == 0 {
    return []
  }
  let ncol = mat[0].length()
  let means : Array[Double] = []
  for j = 0; j < ncol; j = j + 1 {
    let mut sum = 0.0
    for i = 0; i < nrow; i = i + 1 {
      sum = sum + mat[i][j]
    }
    means.push(sum / nrow.to_double())
  }
  means
}

///|
/// Calculate row means of a matrix.
fn row_means(mat : Array[Array[Double]]) -> Array[Double] {
  let nrow = mat.length()
  let means : Array[Double] = []
  for i = 0; i < nrow; i = i + 1 {
    let row = mat[i]
    let ncol = row.length()
    let mut sum = 0.0
    for j = 0; j < ncol; j = j + 1 {
      sum = sum + row[j]
    }
    means.push(sum / ncol.to_double())
  }
  means
}

///|
/// Calculate row variances of a matrix.
fn row_vars(mat : Array[Array[Double]]) -> Array[Double] {
  let nrow = mat.length()
  let vars : Array[Double] = []
  for i = 0; i < nrow; i = i + 1 {
    let row = mat[i]
    let ncol = row.length()
    let mean = {
      let mut sum = 0.0
      for j = 0; j < ncol; j = j + 1 {
        sum = sum + row[j]
      }
      sum / ncol.to_double()
    }
    let mut var_sum = 0.0
    for j = 0; j < ncol; j = j + 1 {
      let diff = row[j] - mean
      var_sum = var_sum + diff * diff
    }
    vars.push(var_sum / (ncol - 1).to_double())
  }
  vars
}

///|
/// Standardize a matrix to have mean 0 and variance 1 by rows.
fn standardize_rows(
  mat : Array[Array[Double]],
) -> (Array[Array[Double]], Array[Double], Array[Double]) {
  let nrow = mat.length()
  let means = row_means(mat)
  let vars = row_vars(mat)
  let result : Array[Array[Double]] = []
  for i = 0; i < nrow; i = i + 1 {
    let row = mat[i]
    let ncol = row.length()
    let mean = means[i]
    let sd = if vars[i] > 0.0 { vars[i].sqrt() } else { 1.0 }
    let new_row : Array[Double] = []
    for j = 0; j < ncol; j = j + 1 {
      new_row.push((row[j] - mean) / sd)
    }
    result.push(new_row)
  }
  (result, means, vars)
}

///|
/// ComBat batch effect correction using empirical Bayes approach.
/// Simplified implementation for batch adjustment.
pub fn combat(
  dat : Array[Array[Double]],
  batch : Array[Int],
  mod? : Array[Array[Double]] = [],
) -> ComBatResult {
  let nrow = dat.length()
  let ncol = dat[0].length()

  let batches = {
    let unique : Array[Int] = []
    for b in batch {
      let mut found = false
      for u in unique {
        if u == b {
          found = true
          break
        }
      }
      if !found {
        unique.push(b)
      }
    }
    unique
  }
  let n_batch = batches.length()

  let (s_data, means, vars) = standardize_rows(dat)

  let batch_design : Array[Array[Double]] = []
  for j = 0; j < ncol; j = j + 1 {
    let row : Array[Double] = [1.0]
    for k = 1; k < n_batch; k = k + 1 {
      if batch[j] == batches[k] {
        row.push(1.0)
      } else {
        row.push(0.0)
      }
    }
    batch_design.push(row)
  }

  let gamma_hat : Array[Array[Double]] = []
  let delta_hat : Array[Array[Double]] = []

  for i = 0; i < nrow; i = i + 1 {
    let row = s_data[i]
    let gene_gamma : Array[Double] = []
    let gene_delta : Array[Double] = []

    for k = 0; k < n_batch; k = k + 1 {
      let batch_vals : Array[Double] = []
      for j = 0; j < ncol; j = j + 1 {
        if batch[j] == batches[k] {
          batch_vals.push(row[j])
        }
      }
      let n_b = batch_vals.length()
      if n_b > 0 {
        let mut mean_b = 0.0
        for v in batch_vals {
          mean_b = mean_b + v
        }
        mean_b = mean_b / n_b.to_double()
        gene_gamma.push(mean_b)

        if n_b > 1 {
          let mut var_b = 0.0
          for v in batch_vals {
            let diff = v - mean_b
            var_b = var_b + diff * diff
          }
          gene_delta.push(var_b / (n_b - 1).to_double())
        } else {
          gene_delta.push(1.0)
        }
      } else {
        gene_gamma.push(0.0)
        gene_delta.push(1.0)
      }
    }
    gamma_hat.push(gene_gamma)
    delta_hat.push(gene_delta)
  }

  let gamma_bar : Array[Double] = []
  let tau2 : Array[Double] = []
  for k = 0; k < n_batch; k = k + 1 {
    let mut sum_g = 0.0
    for i = 0; i < nrow; i = i + 1 {
      sum_g = sum_g + gamma_hat[i][k]
    }
    gamma_bar.push(sum_g / nrow.to_double())

    let mut var_sum = 0.0
    let gb = gamma_bar[k]
    for i = 0; i < nrow; i = i + 1 {
      let diff = gamma_hat[i][k] - gb
      var_sum = var_sum + diff * diff
    }
    tau2.push(var_sum / (nrow - 1).to_double())
  }

  let lambda : Array[Double] = []
  let theta : Array[Double] = []
  for k = 0; k < n_batch; k = k + 1 {
    let mut sum_d = 0.0
    for i = 0; i < nrow; i = i + 1 {
      sum_d = sum_d + delta_hat[i][k]
    }
    let mean_d = sum_d / nrow.to_double()
    let mut sum_sq = 0.0
    for i = 0; i < nrow; i = i + 1 {
      let diff = delta_hat[i][k] - mean_d
      sum_sq = sum_sq + diff * diff
    }
    let var_d = sum_sq / (nrow - 1).to_double()
    lambda.push(mean_d * mean_d / var_d + 2.0)
    theta.push(var_d / mean_d)
  }

  let gamma_star : Array[Array[Double]] = []
  let delta_star : Array[Array[Double]] = []
  for i = 0; i < nrow; i = i + 1 {
    let gs : Array[Double] = []
    let ds : Array[Double] = []
    for k = 0; k < n_batch; k = k + 1 {
      let g_hat = gamma_hat[i][k]
      let d_hat = delta_hat[i][k]
      let t2 = tau2[k]
      let gb = gamma_bar[k]

      let gamma_post_mean = if t2 > 0.0 {
        (g_hat / d_hat + gb / t2) / (1.0 / d_hat + 1.0 / t2)
      } else {
        g_hat
      }

      let d_post = if lambda[k] > 0.0 {
        (lambda[k] * theta[k] + d_hat * 3.0) / (lambda[k] + 2.0)
      } else {
        d_hat
      }

      gs.push(gamma_post_mean)
      ds.push(d_post)
    }
    gamma_star.push(gs)
    delta_star.push(ds)
  }

  let adjusted : Array[Array[Double]] = []
  for i = 0; i < nrow; i = i + 1 {
    let new_row : Array[Double] = []
    let mean_i = means[i]
    let var_i = vars[i]
    let sd_i = if var_i > 0.0 { var_i.sqrt() } else { 1.0 }
    for j = 0; j < ncol; j = j + 1 {
      let batch_idx = {
        let mut idx = 0
        for k = 0; k < n_batch; k = k + 1 {
          if batch[j] == batches[k] {
            idx = k
            break
          }
        }
        idx
      }
      let gs = gamma_star[i][batch_idx]
      let ds = delta_star[i][batch_idx]
      let sd_ds = if ds > 0.0 { ds.sqrt() } else { 1.0 }
      let adj_std = (s_data[i][j] - gs) / sd_ds
      new_row.push(adj_std * sd_i + mean_i)
    }
    adjusted.push(new_row)
  }

  ComBatResult::new(adjusted, gamma_hat, delta_hat)
}

///|
/// Result of Surrogate Variable Analysis.
pub struct SvaResult {
  sv : Array[Array[Double]]
  pprob_gam : Array[Double]
  n_sv : Int
}

///|
/// Create a new SVA result.
pub fn SvaResult::new(
  sv : Array[Array[Double]],
  pprob_gam : Array[Double],
  n_sv : Int,
) -> SvaResult {
  SvaResult::{ sv, pprob_gam, n_sv }
}

///|
/// Simple PCA using power iteration to get top components.
fn pca_simple(
  mat : Array[Array[Double]],
  n_components : Int,
) -> (Array[Array[Double]], Array[Double]) {
  let nrow = mat.length()
  let ncol = mat[0].length()

  let means = col_means(mat)
  let centered : Array[Array[Double]] = []
  for i = 0; i < nrow; i = i + 1 {
    let row : Array[Double] = []
    for j = 0; j < ncol; j = j + 1 {
      row.push(mat[i][j] - means[j])
    }
    centered.push(row)
  }

  let eigenvectors : Array[Array[Double]] = []
  let eigenvalues : Array[Double] = []
  let mut residual = centered

  for c = 0; c < n_components; c = c + 1 {
    let v : Array[Double] = []
    for j = 0; j < ncol; j = j + 1 {
      v.push(1.0 / ncol.to_double().sqrt())
    }

    for _iter = 0; _iter < 100; _iter = _iter + 1 {
      let v_new : Array[Double] = []
      for j = 0; j < ncol; j = j + 1 {
        v_new.push(0.0)
      }

      for i = 0; i < nrow; i = i + 1 {
        let mut dot = 0.0
        for j = 0; j < ncol; j = j + 1 {
          dot = dot + residual[i][j] * v[j]
        }
        for j = 0; j < ncol; j = j + 1 {
          v_new[j] = v_new[j] + residual[i][j] * dot
        }
      }

      let mut norm = 0.0
      for j = 0; j < ncol; j = j + 1 {
        norm = norm + v_new[j] * v_new[j]
      }
      norm = norm.sqrt()
      if norm > 0.0 {
        for j = 0; j < ncol; j = j + 1 {
          v[j] = v_new[j] / norm
        }
      }
    }

    let mut eigenvalue = 0.0
    for i = 0; i < nrow; i = i + 1 {
      let mut dot = 0.0
      for j = 0; j < ncol; j = j + 1 {
        dot = dot + residual[i][j] * v[j]
      }
      eigenvalue = eigenvalue + dot * dot
    }
    eigenvalue = eigenvalue / (nrow - 1).to_double()

    eigenvectors.push(v)
    eigenvalues.push(eigenvalue)

    let new_residual : Array[Array[Double]] = []
    for i = 0; i < nrow; i = i + 1 {
      let mut proj = 0.0
      for j = 0; j < ncol; j = j + 1 {
        proj = proj + residual[i][j] * v[j]
      }
      let row : Array[Double] = []
      for j = 0; j < ncol; j = j + 1 {
        row.push(residual[i][j] - proj * v[j])
      }
      new_residual.push(row)
    }
    residual = new_residual
  }

  (eigenvectors, eigenvalues)
}

///|
/// Estimate number of surrogate variables.
fn num_sv(dat : Array[Array[Double]], mod : Array[Array[Double]]) -> Int {
  let nrow = dat.length()
  let ncol = dat[0].length()
  let (s_data, _, _) = standardize_rows(dat)

  let resid : Array[Array[Double]] = []
  for i = 0; i < nrow; i = i + 1 {
    let row : Array[Double] = []
    for j = 0; j < ncol; j = j + 1 {
      row.push(s_data[i][j])
    }
    resid.push(row)
  }

  let max_comp = if 10 < ncol - 1 { 10 } else { ncol - 1 }
  let (_, eigenvalues) = pca_simple(resid, max_comp)

  let mut n_sv = 0
  let threshold = 1.0
  for i = 0; i < eigenvalues.length(); i = i + 1 {
    if eigenvalues[i] > threshold {
      n_sv = n_sv + 1
    }
  }
  if n_sv == 0 { n_sv = 1 }
  let max_sv = ncol - 2
  if n_sv > max_sv { max_sv } else { n_sv }
}

///|
/// Perform Surrogate Variable Analysis (simplified implementation).
pub fn sva(
  dat : Array[Array[Double]],
  mod : Array[Array[Double]],
  mod0? : Array[Array[Double]] = [],
  n_sv? : Int = 0,
) -> SvaResult {
  let nrow = dat.length()
  let ncol = dat[0].length()
  let (s_data, _, _) = standardize_rows(dat)

  let n = if n_sv > 0 { n_sv } else { num_sv(dat, mod) }

  let resid : Array[Array[Double]] = []
  for i = 0; i < nrow; i = i + 1 {
    let row : Array[Double] = []
    for j = 0; j < ncol; j = j + 1 {
      row.push(s_data[i][j])
    }
    resid.push(row)
  }

  let transposed : Array[Array[Double]] = []
  for j = 0; j < ncol; j = j + 1 {
    let row : Array[Double] = []
    for i = 0; i < nrow; i = i + 1 {
      row.push(resid[i][j])
    }
    transposed.push(row)
  }

  let (eigenvectors, _eigenvalues) = pca_simple(transposed, n)

  let sv = eigenvectors

  let pprob_gam : Array[Double] = []
  for i = 0; i < nrow; i = i + 1 {
    let mut rss0 = 0.0
    for j = 0; j < ncol; j = j + 1 {
      rss0 = rss0 + s_data[i][j] * s_data[i][j]
    }
    let mut rss1 = 0.0
    for j = 0; j < ncol; j = j + 1 {
      let mut fitted = 0.0
      for k = 0; k < n; k = k + 1 {
        fitted = fitted + sv[k][j] * 0.0
      }
      let res = s_data[i][j] - fitted
      rss1 = rss1 + res * res
    }
    let f_stat = if rss1 > 0.0 {
      (rss0 - rss1) / n.to_double() / (rss1 / (ncol - n - 1).to_double())
    } else {
      1.0
    }
    pprob_gam.push(if f_stat > 1.0 { 0.8 } else { 0.2 })
  }

  SvaResult::new(sv, pprob_gam, n)
}

///|
/// Remove batch effects using simple mean-centering per batch.
pub fn combat_simple(
  dat : Array[Array[Double]],
  batch : Array[Int],
) -> Array[Array[Double]] {
  let nrow = dat.length()
  let ncol = dat[0].length()

  let batches = {
    let unique : Array[Int] = []
    for b in batch {
      let mut found = false
      for u in unique {
        if u == b {
          found = true
          break
        }
      }
      if !found {
        unique.push(b)
      }
    }
    unique
  }
  let n_batch = batches.length()

  let result : Array[Array[Double]] = []
  for i = 0; i < nrow; i = i + 1 {
    let row = dat[i]

    let batch_means_arr : Array[Double] = []
    let batch_counts_arr : Array[Int] = []
    for k = 0; k < n_batch; k = k + 1 {
      batch_means_arr.push(0.0)
      batch_counts_arr.push(0)
    }

    for j = 0; j < ncol; j = j + 1 {
      let b = batch[j]
      let mut batch_idx = 0
      for k = 0; k < n_batch; k = k + 1 {
        if batches[k] == b {
          batch_idx = k
          break
        }
      }
      batch_counts_arr[batch_idx] = batch_counts_arr[batch_idx] + 1
      batch_means_arr[batch_idx] = batch_means_arr[batch_idx] + row[j]
    }

    let mut grand_mean = 0.0
    for j = 0; j < ncol; j = j + 1 {
      grand_mean = grand_mean + row[j]
    }
    grand_mean = grand_mean / ncol.to_double()

    let new_row : Array[Double] = []
    for j = 0; j < ncol; j = j + 1 {
      let b = batch[j]
      let mut batch_idx = 0
      for k = 0; k < n_batch; k = k + 1 {
        if batches[k] == b {
          batch_idx = k
          break
        }
      }
      let count = batch_counts_arr[batch_idx]
      let b_mean = if count > 0 {
        batch_means_arr[batch_idx] / count.to_double()
      } else {
        grand_mean
      }
      new_row.push(row[j] - b_mean + grand_mean)
    }
    result.push(new_row)
  }

  result
}

///|
/// Sample data matrix for testing/demo purposes.
pub fn sample_expression_matrix() -> (Array[Array[Double]], Array[String], Array[String], Array[Int]) {
  let genes = [
    "GeneA", "GeneB", "GeneC", "GeneD", "GeneE",
    "GeneF", "GeneG", "GeneH", "GeneI", "GeneJ",
    "GeneK", "GeneL", "GeneM", "GeneN", "GeneO",
    "GeneP", "GeneQ", "GeneR", "GeneS", "GeneT",
  ]
  let samples = [
    "s1", "s2", "s3", "s4", "s5", "s6",
    "s7", "s8", "s9", "s10", "s11", "s12",
  ]
  let batch = [0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1]

  let dat : Array[Array[Double]] = []
  for i = 0; i < 20; i = i + 1 {
    let row : Array[Double] = []
    for j = 0; j < 12; j = j + 1 {
      let base = 5.0 + i.to_double() * 0.2
      let batch_effect = if batch[j] == 1 { 1.5 + i.to_double() * 0.1 } else { 0.0 }
      let noise = (i * 7 + j * 13).to_double() % 10.0 / 10.0 - 0.5
      row.push(base + batch_effect + noise)
    }
    dat.push(row)
  }

  (dat, genes, samples, batch)
}
