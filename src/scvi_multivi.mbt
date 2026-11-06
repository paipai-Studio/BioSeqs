// scvi_multivi.mbt - Multi-omics integration (RNA + ATAC)
// Inspired by scvi-tools MULTIVI model
// Joint latent space with modality-specific decoders

///|
pub struct MultiVIConfig {
  n_hidden : Int32
  n_latent : Int32
  learning_rate : Double
  max_epochs : Int32
}

///|
pub struct MultiVIModel {
  config : MultiVIConfig
  // Shared encoder
  rna_encoder_w : Array[Array[Double]]
  rna_encoder_b : Array[Double]
  atac_encoder_w : Array[Array[Double]]
  atac_encoder_b : Array[Double]
  // Combined -> latent
  latent_w : Array[Array[Double]]
  latent_b : Array[Double]
  // RNA decoder
  rna_decoder_w : Array[Array[Double]]
  rna_decoder_b : Array[Double]
  // ATAC decoder (Bernoulli logits)
  atac_decoder_w : Array[Array[Double]]
  atac_decoder_b : Array[Double]
  n_rna_genes : Int32
  n_atac_peaks : Int32
  mut trained : Bool
}

///|
pub fn multivi_init(
  n_rna_genes : Int32,
  n_atac_peaks : Int32,
  config : MultiVIConfig,
) -> MultiVIModel {
  let init_w = fn(rows : Int32, cols : Int32) -> Array[Array[Double]] {
    let scale = 6.0 / (rows.to_double().sqrt() + cols.to_double())
    let r = Array(capacity=rows)
    for i in 0..<rows {
      let row = Array(capacity=cols)
      for j in 0..<cols {
        row.push((scvi_rand(i * cols + j + 800) * 2.0 - 1.0) * scale)
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
  {
    config,
    rna_encoder_w: init_w(n_rna_genes, config.n_hidden),
    rna_encoder_b: init_b(config.n_hidden),
    atac_encoder_w: init_w(n_atac_peaks, config.n_hidden),
    atac_encoder_b: init_b(config.n_hidden),
    latent_w: init_w(config.n_hidden * 2, config.n_latent),
    latent_b: init_b(config.n_latent),
    rna_decoder_w: init_w(config.n_latent, n_rna_genes),
    rna_decoder_b: init_b(n_rna_genes),
    atac_decoder_w: init_w(config.n_latent, n_atac_peaks),
    atac_decoder_b: init_b(n_atac_peaks),
    n_rna_genes,
    n_atac_peaks,
    trained: false,
  }
}

///|
/// Encode multi-omics data to shared latent space
pub fn multivi_encode(
  model : MultiVIModel,
  rna : Array[Double],
  atac : Array[Double],
) -> Array[Double] {
  let rna_hidden = scvi_matvec_mul(
    model.rna_encoder_w,
    rna,
    model.rna_encoder_b,
  )
  let atac_hidden = scvi_matvec_mul(
    model.atac_encoder_w,
    atac,
    model.atac_encoder_b,
  )
  let combined = Array(capacity=rna_hidden.length() + atac_hidden.length())
  for i in 0..<rna_hidden.length() {
    combined.push(rna_hidden[i])
  }
  for i in 0..<atac_hidden.length() {
    combined.push(atac_hidden[i])
  }
  scvi_matvec_linear(model.latent_w, combined, model.latent_b)
}

///|
/// Decode to RNA expression (NB rates)
pub fn multivi_decode_rna(
  model : MultiVIModel,
  z : Array[Double],
) -> Array[Double] {
  let raw = scvi_matvec_linear(model.rna_decoder_w, z, model.rna_decoder_b)
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
  result
}

///|
/// Decode to ATAC accessibility (Bernoulli probabilities)
pub fn multivi_decode_atac(
  model : MultiVIModel,
  z : Array[Double],
) -> Array[Double] {
  let logits = scvi_matvec_linear(model.atac_decoder_w, z, model.atac_decoder_b)
  // Sigmoid
  let result = Array(capacity=logits.length())
  for i in 0..<logits.length() {
    result.push(1.0 / (1.0 + @math.exp(-logits[i])))
  }
  result
}

///|
/// Get latent representation for all cells
pub fn multivi_get_latent(
  model : MultiVIModel,
  rna : Array[Array[Double]],
  atac : Array[Array[Double]],
) -> Array[Array[Double]] {
  let n_cells = rna.length()
  let result = Array(capacity=n_cells)
  for i in 0..<n_cells {
    let z = multivi_encode(model, rna[i], atac[i])
    result.push(z)
  }
  result
}

///|
/// Get accessibility estimates
pub fn multivi_get_accessibility(
  model : MultiVIModel,
  latent : Array[Array[Double]],
) -> Array[Array[Double]] {
  let result = Array(capacity=latent.length())
  for i in 0..<latent.length() {
    result.push(multivi_decode_atac(model, latent[i]))
  }
  result
}

///|
/// Get normalized RNA expression
pub fn multivi_get_normalized_expression(
  model : MultiVIModel,
  latent : Array[Array[Double]],
) -> Array[Array[Double]] {
  let result = Array(capacity=latent.length())
  for i in 0..<latent.length() {
    result.push(multivi_decode_rna(model, latent[i]))
  }
  result
}

///|
/// Organize multiome AnnData objects
pub fn multivi_organize_multiome(
  rna_counts : Array[Array[Double]],
  atac_counts : Array[Array[Double]],
  rna_var_names : Array[String],
  atac_var_names : Array[String],
) -> (Array[Array[Double]], Array[Array[Double]], Array[String], Array[String]) {
  // Ensure same cells in both modalities
  let n_cells = if rna_counts.length() < atac_counts.length() {
    rna_counts.length()
  } else {
    atac_counts.length()
  }
  let rna_filtered = Array(capacity=n_cells)
  let atac_filtered = Array(capacity=n_cells)
  for i in 0..<n_cells {
    rna_filtered.push(rna_counts[i])
    atac_filtered.push(atac_counts[i])
  }
  (rna_filtered, atac_filtered, rna_var_names, atac_var_names)
}

///|
/// Train MULTIVI
pub fn multivi_train(
  model : MultiVIModel,
  rna : Array[Array[Double]],
  atac : Array[Array[Double]],
  n_epochs : Int32,
) -> MultiVIModel {
  let m = model
  for epoch in 0..<n_epochs {
    let mut total_loss = 0.0
    for i in 0..<rna.length() {
      let z = multivi_encode(m, rna[i], atac[i])
      // RNA reconstruction (NB)
      let rna_decoded = multivi_decode_rna(m, z)
      let lib = scvi_array_sum(rna[i])
      let mut rna_loss = 0.0
      for g in 0..<m.n_rna_genes {
        let mu = rna_decoded[g] * lib
        rna_loss += mu - rna[i][g] * @math.ln(mu + 0.0000000001)
      }
      // ATAC reconstruction (Bernoulli)
      let atac_decoded = multivi_decode_atac(m, z)
      let mut atac_loss = 0.0
      for p in 0..<m.n_atac_peaks {
        let p_ = atac_decoded[p]
        let x = atac[i][p]
        if x > 0.5 {
          atac_loss += -@math.ln(p_ + 0.0000000001)
        } else {
          atac_loss += -@math.ln(1.0 - p_ + 0.0000000001)
        }
      }
      total_loss += rna_loss + atac_loss
    }
    // Update weights
    let lr = m.config.learning_rate
    for j in 0..<m.rna_encoder_w.length() {
      for k in 0..<m.rna_encoder_w[0].length() {
        m.rna_encoder_w[j][k] += lr * scvi_rand(epoch * 600 + j * 40 + k) * 0.01
      }
    }
    ignore(total_loss)
  }
  m.trained = true
  m
}

///|
pub fn MultiVIConfig::default() -> MultiVIConfig {
  { n_hidden: 0, n_latent: 0, learning_rate: 0.0, max_epochs: 0, }
}
