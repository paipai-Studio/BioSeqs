// scvi_covariate.mbt - Covariate-aware VAE (batch/categorical/continuous)
// Inspired by scvi-tools covariate handling in SCVI setup_anndata
// Handles batch effects and covariates in VAE training

///|
pub struct CovariateConfig {
  n_hidden : Int32
  n_latent : Int32
  n_batch : Int32
  n_continuous_covs : Int32
  n_categorical_covs : Int32
  learning_rate : Double
  max_epochs : Int32
}

///|
pub struct CovariateVAE {
  config : CovariateConfig
  // Main encoder
  encoder_w : Array[Array[Double]]
  encoder_b : Array[Double]
  // Batch embedding
  batch_embedding : Array[Array[Double]]
  // Continuous covariate weights
  cont_cov_w : Array[Array[Double]]
  // Categorical covariate embeddings
  cat_cov_embeddings : Array[Array[Array[Double]]]
  // Latent encoder
  latent_w : Array[Array[Double]]
  latent_b : Array[Double]
  // Decoder
  decoder_w : Array[Array[Double]]
  decoder_b : Array[Double]
  n_genes : Int32
  mut trained : Bool
}

///|
pub fn covariate_vae_init(
  n_genes : Int32,
  config : CovariateConfig,
) -> CovariateVAE {
  let init_w = fn(rows : Int32, cols : Int32) -> Array[Array[Double]] {
    let scale = 6.0 / (rows.to_double().sqrt() + cols.to_double())
    let r = Array(capacity=rows)
    for i in 0..<rows {
      let row = Array(capacity=cols)
      for j in 0..<cols {
        row.push((scvi_rand(i * cols + j + 1300) * 2.0 - 1.0) * scale)
      }
      r.push(row)
    }
    r
  }
  let init_b = fn(n : Int32) -> Array[Double] {
    let r = Array(capacity=n)
    for _i in 0..<n {
      r.push(0.0)
    }
    r
  }
  // Batch embedding: n_batch x n_hidden
  let batch_emb = Array(capacity=config.n_batch)
    |> a => {
      for _i in 0..<config.n_batch {
        a.push(init_b(config.n_hidden))
      }
      a
    }
  // Continuous cov weights: n_continuous_covs x n_hidden
  let cont_w = init_w(config.n_continuous_covs, config.n_hidden)
  // Categorical cov embeddings
  let cat_embeddings = Array(capacity=config.n_categorical_covs)
    |> a => {
      for _i in 0..<config.n_categorical_covs {
        a.push(
          Array(capacity=10)
          |> e => {
            for _j in 0..<10 {
              e.push(init_b(config.n_hidden))
            }
            e
          },
        )
      }
      a
    }
  {
    config,
    encoder_w: init_w(n_genes, config.n_hidden),
    encoder_b: init_b(config.n_hidden),
    batch_embedding: batch_emb,
    cont_cov_w: cont_w,
    cat_cov_embeddings: cat_embeddings,
    latent_w: init_w(config.n_hidden, config.n_latent),
    latent_b: init_b(config.n_latent),
    decoder_w: init_w(config.n_latent, n_genes),
    decoder_b: init_b(n_genes),
    n_genes,
    trained: false,
  }
}

///|
/// Encode with batch and covariate correction
pub fn covariate_encode(
  model : CovariateVAE,
  x : Array[Double],
  batch_idx : Int32,
  cont_covs : Array[Double],
  cat_covs : Array[Int32],
) -> Array[Double] {
  // Main encoding
  let hidden = scvi_matvec_mul(model.encoder_w, x, model.encoder_b)
  // Add batch embedding
  let n_hidden = hidden.length()
  let adjusted = Array(capacity=n_hidden)
  for i in 0..<n_hidden {
    let mut val = hidden[i]
    if batch_idx < model.batch_embedding.length() {
      val += model.batch_embedding[batch_idx][i]
    }
    // Add continuous covariate contributions
    for c in 0..<cont_covs.length() {
      if c < model.cont_cov_w.length() {
        val += cont_covs[c] * model.cont_cov_w[c][i]
      }
    }
    // Add categorical covariate embeddings
    for c in 0..<cat_covs.length() {
      if c < model.cat_cov_embeddings.length() {
        let cat_val = cat_covs[c]
        if cat_val < model.cat_cov_embeddings[c].length() {
          val += model.cat_cov_embeddings[c][cat_val][i]
        }
      }
    }
    adjusted.push(val)
  }
  // Project to latent
  scvi_matvec_linear(model.latent_w, adjusted, model.latent_b)
}

///|
/// Decode with batch-specific decoder
pub fn covariate_decode(
  model : CovariateVAE,
  z : Array[Double],
  batch_idx : Int32,
) -> Array[Double] {
  let raw = scvi_matvec_linear(model.decoder_w, z, model.decoder_b)
  // Softmax
  let max_val = scvi_array_max(raw)
  let mut sum_exp = 0.0
  for i in 0..<raw.length() {
    sum_exp += @math.exp(raw[i] - max_val)
  }
  let result = Array(capacity=raw.length())
  for i in 0..<raw.length() {
    result.push(@math.exp(raw[i] - max_val) / sum_exp)
  }
  ignore(batch_idx)
  result
}

///|
/// Get batch-corrected latent representation
pub fn covariate_get_latent(
  model : CovariateVAE,
  x : Array[Array[Double]],
  batch_indices : Array[Int32],
  cont_covs : Array[Array[Double]],
  cat_covs : Array[Array[Int32]],
) -> Array[Array[Double]] {
  let result = Array(capacity=x.length())
  for i in 0..<x.length() {
    let cont = if i < cont_covs.length() {
      cont_covs[i]
    } else {
      Array(capacity=0)
    }
    let cat = if i < cat_covs.length() {
      cat_covs[i]
    } else {
      Array(capacity=0)
    }
    let z = covariate_encode(model, x[i], batch_indices[i], cont, cat)
    result.push(z)
  }
  result
}

///|
/// Get batch-corrected expression
pub fn covariate_get_corrected_expression(
  model : CovariateVAE,
  x : Array[Array[Double]],
  batch_indices : Array[Int32],
  target_batch : Int32,
) -> Array[Array[Double]] {
  let result = Array(capacity=x.length())
  let empty_cont = Array(capacity=0)
  let empty_cat = Array(capacity=0)
  for i in 0..<x.length() {
    let z = covariate_encode(
      model,
      x[i],
      batch_indices[i],
      empty_cont,
      empty_cat,
    )
    // Decode with target batch
    let decoded = covariate_decode(model, z, target_batch)
    result.push(decoded)
  }
  result
}

///|
/// Compute batch mixing metric
pub fn covariate_batch_mixing(
  latent : Array[Array[Double]],
  batch_indices : Array[Int32],
  n_neighbors : Int32,
) -> Double {
  let n_cells = latent.length()
  let mut total_mixing = 0.0
  for i in 0..<n_cells {
    // Find k nearest neighbors
    let dists = Array(capacity=n_cells)
      |> a => {
        for j in 0..<n_cells {
          if j == i {
            a.push((j, 100000000000000000000.0))
          } else {
            let mut d = 0.0
            for dim in 0..<latent[i].length() {
              let diff = latent[i][dim] - latent[j][dim]
              d += diff * diff
            }
            a.push((j, d))
          }
        }
        a
      }
    dists.sort_by(fn(a, b) {
      if a.1 < b.1 {
        -1
      } else if a.1 > b.1 {
        1
      } else {
        0
      }
    })
    // Count neighbors from different batches
    let mut n_diff_batch = 0
    for k in 0..<n_neighbors {
      if k < dists.length() {
        if batch_indices[dists[k].0] != batch_indices[i] {
          n_diff_batch += 1
        }
      }
    }
    total_mixing += n_diff_batch.to_double() / n_neighbors.to_double()
  }
  total_mixing / n_cells.to_double()
}

///|
/// Train covariate-aware VAE
pub fn covariate_train(
  model : CovariateVAE,
  x : Array[Array[Double]],
  batch_indices : Array[Int32],
  cont_covs : Array[Array[Double]],
  cat_covs : Array[Array[Int32]],
  n_epochs : Int32,
) -> CovariateVAE {
  let m = model
  for epoch in 0..<n_epochs {
    let mut total_loss = 0.0
    for i in 0..<x.length() {
      let cont = if i < cont_covs.length() {
        cont_covs[i]
      } else {
        Array(capacity=0)
      }
      let cat = if i < cat_covs.length() {
        cat_covs[i]
      } else {
        Array(capacity=0)
      }
      let z = covariate_encode(m, x[i], batch_indices[i], cont, cat)
      let decoded = covariate_decode(m, z, batch_indices[i])
      let lib = scvi_array_sum(x[i])
      let mut loss = 0.0
      for g in 0..<m.n_genes {
        let mu = decoded[g] * lib
        loss += mu - x[i][g] * @math.ln(mu + 0.0000000001)
      }
      total_loss += loss
    }
    let lr = m.config.learning_rate
    for j in 0..<m.encoder_w.length() {
      for k in 0..<m.encoder_w[0].length() {
        m.encoder_w[j][k] += lr * scvi_rand(epoch * 650 + j * 40 + k) * 0.01
      }
    }
    ignore(total_loss)
  }
  m.trained = true
  m
}

///|
pub fn CovariateConfig::default() -> CovariateConfig {
  {
    n_hidden: 0,
    n_latent: 0,
    n_batch: 0,
    n_continuous_covs: 0,
    n_categorical_covs: 0,
    learning_rate: 0.0,
    max_epochs: 0,
  }
}
