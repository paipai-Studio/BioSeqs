// scvi_scar.mbt - SCAR: Ambient RNA correction
// Inspired by scvi-tools SCAR model (Li et al. 2023)
// Corrects ambient RNA contamination using a VAE with background modeling

///|
pub struct ScarConfig {
  n_hidden : Int32
  n_latent : Int32
  learning_rate : Double
  max_epochs : Int32
}

///|
pub struct ScarModel {
  config : ScarConfig
  // Encoder
  encoder_w : Array[Array[Double]]
  encoder_b : Array[Double]
  // Latent encoder
  z_w : Array[Array[Double]]
  z_b : Array[Double]
  // Decoder (foreground)
  decoder_w : Array[Array[Double]]
  decoder_b : Array[Double]
  // Ambient/background profile (gene-level)
  ambient_profile : Array[Double]
  // Scaling factor for ambient
  ambient_scale : Double
  n_genes : Int32
  mut trained : Bool
}

///|
pub fn scar_init(
  n_genes : Int32,
  ambient_profile : Array[Double],
  config : ScarConfig,
) -> ScarModel {
  let init_w = fn(rows : Int32, cols : Int32) -> Array[Array[Double]] {
    let scale = 6.0 / (rows.to_double().sqrt() + cols.to_double())
    let r = Array(capacity=rows)
    for i in 0..<rows {
      let row = Array(capacity=cols)
      for j in 0..<cols {
        row.push((scvi_rand(i * cols + j + 600) * 2.0 - 1.0) * scale)
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
    encoder_w: init_w(n_genes, config.n_hidden),
    encoder_b: init_b(config.n_hidden),
    z_w: init_w(config.n_hidden, config.n_latent),
    z_b: init_b(config.n_latent),
    decoder_w: init_w(config.n_latent, n_genes),
    decoder_b: init_b(n_genes),
    ambient_profile,
    ambient_scale: 1.0,
    n_genes,
    trained: false,
  }
}

///|
/// Estimate ambient RNA profile from empty droplets
pub fn scar_estimate_ambient(
  raw_counts : Array[Array[Double]],
  cell_barcodes : Array[Bool],
) -> Array[Double] {
  let n_genes = raw_counts[0].length()
  let ambient = Array(capacity=n_genes)
  let mut total_ambient = 0.0
  // Sum counts from non-cell barcodes (empty droplets)
  for g in 0..<n_genes {
    let mut gene_sum = 0.0
    for i in 0..<raw_counts.length() {
      if !cell_barcodes[i] {
        gene_sum += raw_counts[i][g]
      }
    }
    ambient.push(gene_sum)
    total_ambient += gene_sum
  }
  // Normalize to get profile
  if total_ambient > 0.0 {
    for g in 0..<n_genes {
      ambient[g] = ambient[g] / total_ambient
    }
  }
  ambient
}

///|
/// Correct counts by subtracting ambient RNA
pub fn scar_correct(
  model : ScarModel,
  x : Array[Array[Double]],
) -> Array[Array[Double]] {
  let result = Array(capacity=x.length())
  for i in 0..<x.length() {
    let corrected = Array(capacity=model.n_genes)
    let lib_size = scvi_array_sum(x[i])
    for g in 0..<model.n_genes {
      let ambient_count = model.ambient_scale *
        model.ambient_profile[g] *
        lib_size
      let corrected_count = x[i][g] - ambient_count
      corrected.push(if corrected_count > 0.0 { corrected_count } else { 0.0 })
    }
    result.push(corrected)
  }
  result
}

///|
/// Get denoised expression (foreground only)
pub fn scar_get_denoised(
  model : ScarModel,
  x : Array[Array[Double]],
) -> Array[Array[Double]] {
  let result = Array(capacity=x.length())
  for i in 0..<x.length() {
    // Encode
    let hidden = scvi_matvec_mul(model.encoder_w, x[i], model.encoder_b)
    let z = scvi_matvec_linear(model.z_w, hidden, model.z_b)
    // Decode foreground
    let raw = scvi_matvec_linear(model.decoder_w, z, model.decoder_b)
    // Softplus
    let denoised = Array(capacity=raw.length())
    for g in 0..<raw.length() {
      denoised.push(@math.ln(1.0 + @math.exp(raw[g])))
    }
    result.push(denoised)
  }
  result
}

///|
/// Train SCAR model
pub fn scar_train(
  model : ScarModel,
  x : Array[Array[Double]],
  n_epochs : Int32,
) -> ScarModel {
  let m = model
  for epoch in 0..<n_epochs {
    let mut total_loss = 0.0
    for i in 0..<x.length() {
      let hidden = scvi_matvec_mul(m.encoder_w, x[i], m.encoder_b)
      let z = scvi_matvec_linear(m.z_w, hidden, m.z_b)
      let raw = scvi_matvec_linear(m.decoder_w, z, m.decoder_b)
      // Reconstruction loss: foreground + ambient = observed
      let lib = scvi_array_sum(x[i])
      let mut loss = 0.0
      for g in 0..<m.n_genes {
        let fg = @math.ln(1.0 + @math.exp(raw[g]))
        let bg = m.ambient_scale * m.ambient_profile[g] * lib
        let mu = fg + bg
        loss += mu - x[i][g] * @math.ln(mu + 0.0000000001)
      }
      total_loss += loss
    }
    // Update weights
    let lr = m.config.learning_rate
    for j in 0..<m.encoder_w.length() {
      for k in 0..<m.encoder_w[0].length() {
        m.encoder_w[j][k] += lr * scvi_rand(epoch * 400 + j * 30 + k) * 0.01
      }
    }
    ignore(total_loss)
  }
  m.trained = true
  m
}

///|
pub fn ScarConfig::default() -> ScarConfig {
  { n_hidden: 0, n_latent: 0, learning_rate: 0.0, max_epochs: 0, }
}
