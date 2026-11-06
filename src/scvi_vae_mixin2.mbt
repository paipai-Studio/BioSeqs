///|
/// VAE Mixin methods: aggregated posterior, differential abundance, importance weights,
/// feature correlation, likelihood parameters, latent library size, posterior predictive,
/// reconstruction error.
pub fn scvi_get_aggregated_posterior(
  latent : Array[Array[Double]],
  component_type : String,
) -> (Array[Double], Array[Double]) {
  let n = latent.length()
  if n == 0 {
    return ([], [])
  }
  let d = latent[0].length()
  let means = Array(capacity=d)
  let vars_ = Array(capacity=d)
  for j in 0..<d {
    let mut sum = 0.0
    for i in 0..<n {
      if j < latent[i].length() {
        sum = sum + latent[i][j]
      }
    }
    means.push(sum / n.to_double())
  }
  for j in 0..<d {
    let mut sum = 0.0
    for i in 0..<n {
      if j < latent[i].length() {
        let diff = latent[i][j] - means[j]
        sum = sum + diff * diff
      }
    }
    let variance = sum / n.to_double()
    let adjusted = match component_type {
      "student_t" => variance * (n.to_double() / (n.to_double() - 2.0))
      _ => variance
    }
    vars_.push(adjusted)
  }
  (means, vars_)
}

///|
pub fn scvi_differential_abundance(
  _latent : Array[Array[Double]],
  sample_labels : Array[Int32],
  n_samples : Int32,
) -> Array[Double] {
  let sample_log_probs = Array(capacity=n_samples)
  for s in 0..<n_samples {
    let mut count = 0
    for i in 0..<sample_labels.length() {
      if sample_labels[i] == s {
        count = count + 1
      }
    }
    let prob = count.to_double() / sample_labels.length().to_double()
    if prob > 0.0 {
      sample_log_probs.push(@math.ln(prob))
    } else {
      sample_log_probs.push(-999.0)
    }
  }
  sample_log_probs
}

///|
pub fn scvi_get_importance_weights(
  latent : Array[Array[Double]],
  x : Array[Array[Double]],
  decoder : Array[Array[Double]],
) -> Array[Double] {
  let weights = Array(capacity=latent.length())
  for i in 0..<latent.length() {
    let mut log_w = 0.0
    for g in 0..<x[i].length() {
      let mut rate = 0.0
      for j in 0..<latent[i].length() {
        if j < decoder[g].length() {
          rate = rate + latent[i][j] * decoder[g][j]
        }
      }
      rate = @math.exp(rate)
      if rate > 0.0 {
        log_w = log_w + x[i][g] * @math.ln(rate) - rate
      }
    }
    weights.push(@math.exp(log_w))
  }
  let mut max_w = weights[0]
  for i in 1..<weights.length() {
    if weights[i] > max_w {
      max_w = weights[i]
    }
  }
  let mut sum_w = 0.0
  for i in 0..<weights.length() {
    weights[i] = weights[i] / max_w
    sum_w = sum_w + weights[i]
  }
  for i in 0..<weights.length() {
    weights[i] = weights[i] / sum_w * weights.length().to_double()
  }
  weights
}

///|
pub fn scvi_get_feature_correlation_matrix(
  latent : Array[Array[Double]],
  decoder : Array[Array[Double]],
  n_samples : Int32,
) -> Array[Array[Double]] {
  let n_genes = decoder.length()
  let n_latent = if latent.length() > 0 { latent[0].length() } else { 0 }
  let gene_expr = Array(capacity=n_samples)
  for s in 0..<n_samples {
    let z = Array(capacity=n_latent)
    for j in 0..<n_latent {
      z.push(latent[s % latent.length()][j])
    }
    let expr = Array(capacity=n_genes)
    for g in 0..<n_genes {
      let mut sum = 0.0
      for j in 0..<z.length() {
        if j < decoder[g].length() {
          sum = sum + z[j] * decoder[g][j]
        }
      }
      expr.push(@math.exp(sum))
    }
    gene_expr.push(expr)
  }
  let corr = Array(capacity=n_genes)
  for g1 in 0..<n_genes {
    let row = Array(capacity=n_genes)
    let mut mean1 = 0.0
    for s in 0..<n_samples {
      if g1 < gene_expr[s].length() {
        mean1 = mean1 + gene_expr[s][g1]
      }
    }
    mean1 = mean1 / n_samples.to_double()
    for g2 in 0..<n_genes {
      let mut mean2 = 0.0
      for s in 0..<n_samples {
        if g2 < gene_expr[s].length() {
          mean2 = mean2 + gene_expr[s][g2]
        }
      }
      mean2 = mean2 / n_samples.to_double()
      let mut cov = 0.0
      let mut var1 = 0.0
      let mut var2 = 0.0
      for s in 0..<n_samples {
        if g1 < gene_expr[s].length() && g2 < gene_expr[s].length() {
          let d1 = gene_expr[s][g1] - mean1
          let d2 = gene_expr[s][g2] - mean2
          cov = cov + d1 * d2
          var1 = var1 + d1 * d1
          var2 = var2 + d2 * d2
        }
      }
      let denom = @math.pow(var1 * var2, 0.5)
      if denom > 0.0 {
        row.push(cov / denom)
      } else {
        row.push(0.0)
      }
    }
    corr.push(row)
  }
  corr
}

///|
pub fn scvi_get_likelihood_parameters(
  latent : Array[Array[Double]],
  decoder : Array[Array[Double]],
  dispersion : Array[Double],
) -> Array[(Double, Double)] {
  let n = latent.length()
  let params = Array(capacity=n)
  for i in 0..<n {
    let mut mean_rate = 0.0
    for g in 0..<decoder.length() {
      let mut rate = 0.0
      for j in 0..<latent[i].length() {
        if j < decoder[g].length() {
          rate = rate + latent[i][j] * decoder[g][j]
        }
      }
      mean_rate = mean_rate + @math.exp(rate)
    }
    mean_rate = mean_rate / decoder.length().to_double()
    let mean_disp = if dispersion.length() > 0 { dispersion[0] } else { 1.0 }
    params.push((mean_rate, mean_disp))
  }
  params
}

///|
pub fn scvi_get_latent_library_size(
  latent : Array[Array[Double]],
  library_decoder : Array[Double],
) -> Array[Double] {
  let lib_sizes = Array(capacity=latent.length())
  for i in 0..<latent.length() {
    let mut sum = 0.0
    for j in 0..<latent[i].length() {
      if j < library_decoder.length() {
        sum = sum + latent[i][j] * library_decoder[j]
      }
    }
    lib_sizes.push(@math.exp(sum))
  }
  lib_sizes
}

///|
pub fn scvi_posterior_predictive_sample(
  latent : Array[Array[Double]],
  decoder : Array[Array[Double]],
  dispersion : Double,
  n_samples : Int32,
) -> Array[Array[Double]] {
  let samples = Array(capacity=n_samples)
  for s in 0..<n_samples {
    let idx = s % latent.length()
    let sample = Array(capacity=decoder.length())
    for g in 0..<decoder.length() {
      let mut rate = 0.0
      for j in 0..<latent[idx].length() {
        if j < decoder[g].length() {
          rate = rate + latent[idx][j] * decoder[g][j]
        }
      }
      rate = @math.exp(rate)
      let mean = rate
      let variance = mean + mean * mean / dispersion
      let sample_val = if variance > 0.0 {
        mean + (variance - mean) * 0.5
      } else {
        mean
      }
      sample.push(if sample_val < 0.0 { 0.0 } else { sample_val })
    }
    samples.push(sample)
  }
  samples
}

///|
pub fn scvi_get_reconstruction_error(
  x : Array[Double],
  x_hat : Array[Double],
) -> Double {
  let mut error = 0.0
  for i in 0..<x.length() {
    if i < x_hat.length() && x_hat[i] > 0.0 {
      error = error + x[i] * @math.ln(x[i] / x_hat[i]) - x[i] + x_hat[i]
    }
  }
  error
}

///|
pub fn scvi_get_ranked_features(
  _latent : Array[Array[Double]],
  decoder : Array[Array[Double]],
  target_label : Int32,
) -> Array[(Int32, Double)] {
  let n_genes = decoder.length()
  let scores = Array(capacity=n_genes)
  for g in 0..<n_genes {
    let mut score = 0.0
    for j in 0..<decoder[g].length() {
      score = score + decoder[g][j] * decoder[g][j]
    }
    if target_label < n_genes {
      score = score * (1.0 + target_label.to_double() * 0.01)
    }
    scores.push((g, score))
  }
  scores.sort_by(fn(a, b) {
    if a.1 > b.1 {
      -1
    } else if a.1 < b.1 {
      1
    } else {
      0
    }
  })
  scores
}

///|
pub fn scvi_get_batch_representation(
  batch_embeddings : Array[Array[Double]],
  batch_indices : Array[Int32],
) -> Array[Array[Double]] {
  let result = Array(capacity=batch_indices.length())
  for i in 0..<batch_indices.length() {
    let b = batch_indices[i]
    if b >= 0 && b < batch_embeddings.length() {
      result.push(batch_embeddings[b])
    } else {
      result.push([])
    }
  }
  result
}
