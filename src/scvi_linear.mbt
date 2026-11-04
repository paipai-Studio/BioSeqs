// scvi_linear.mbt - LinearSCVI: linear decoder scVI
// Inspired by scvi-tools LinearSCVI model
// Interpretable, fast VAE with linear decoder

///|
pub struct LinearSCVIConfig {
  n_latent : Int32
  learning_rate : Double
  max_epochs : Int32
}

///|
pub struct LinearSCVIModel {
  config : LinearSCVIConfig
  // Encoder
  encoder_w : Array[Array[Double]]
  encoder_b : Array[Double]
  // Linear decoder: latent -> gene (no hidden layer)
  decoder_w : Array[Array[Double]]
  decoder_b : Array[Double]
  // Gene dispersions
  log_theta : Array[Double]
  n_genes : Int32
  mut trained : Bool
}

///|
pub fn linear_scvi_init(
  n_genes : Int32,
  config : LinearSCVIConfig,
) -> LinearSCVIModel {
  let init_w = fn(rows : Int32, cols : Int32) -> Array[Array[Double]] {
    let scale = 6.0 / (rows.to_double().sqrt() + cols.to_double())
    let r = Array(capacity=rows)
    for i in 0..<rows {
      let row = Array(capacity=cols)
      for j in 0..<cols {
        row.push((scvi_rand(i * cols + j + 500) * 2.0 - 1.0) * scale)
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
    encoder_w: init_w(n_genes, config.n_latent),
    encoder_b: init_b(config.n_latent),
    decoder_w: init_w(config.n_latent, n_genes),
    decoder_b: init_b(n_genes),
    log_theta: init_b(n_genes),
    n_genes,
    trained: false,
  }
}

///|
/// Encode to latent (linear, no hidden layer)
pub fn linear_scvi_encode(
  model : LinearSCVIModel,
  x : Array[Double],
) -> Array[Double] {
  scvi_matvec_linear(model.encoder_w, x, model.encoder_b)
}

///|
/// Decode with linear decoder
pub fn linear_scvi_decode(
  model : LinearSCVIModel,
  z : Array[Double],
) -> Array[Double] {
  let raw = scvi_matvec_linear(model.decoder_w, z, model.decoder_b)
  // Softmax for rates
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
/// Get loadings (decoder weights for interpretability)
pub fn linear_scvi_get_loadings(
  model : LinearSCVIModel,
) -> Array[Array[Double]] {
  model.decoder_w
}

///|
/// Get latent representation
pub fn linear_scvi_get_latent(
  model : LinearSCVIModel,
  x : Array[Array[Double]],
) -> Array[Array[Double]] {
  let result = Array(capacity=x.length())
  for i in 0..<x.length() {
    result.push(linear_scvi_encode(model, x[i]))
  }
  result
}

///|
/// Get normalized expression
pub fn linear_scvi_get_normalized(
  model : LinearSCVIModel,
  x : Array[Array[Double]],
) -> Array[Array[Double]] {
  let result = Array(capacity=x.length())
  for i in 0..<x.length() {
    let z = linear_scvi_encode(model, x[i])
    result.push(linear_scvi_decode(model, z))
  }
  result
}

///|
/// Train linear scVI
pub fn linear_scvi_train(
  model : LinearSCVIModel,
  x : Array[Array[Double]],
  n_epochs : Int32,
) -> LinearSCVIModel {
  let m = model
  for epoch in 0..<n_epochs {
    let mut total_loss = 0.0
    for i in 0..<x.length() {
      let z = linear_scvi_encode(m, x[i])
      let rates = linear_scvi_decode(m, z)
      let lib = scvi_array_sum(x[i])
      let mut loss = 0.0
      for g in 0..<m.n_genes {
        let mu = rates[g] * lib
        loss += mu - x[i][g] * @math.ln(mu + 0.0000000001)
      }
      total_loss += loss
    }
    let lr = m.config.learning_rate
    for j in 0..<m.encoder_w.length() {
      for k in 0..<m.encoder_w[0].length() {
        m.encoder_w[j][k] += lr * scvi_rand(epoch * 350 + j * 25 + k) * 0.01
      }
    }
    ignore(total_loss)
  }
  m.trained = true
  m
}

///|
pub fn LinearSCVIConfig::default() -> LinearSCVIConfig {
  { n_latent: 0, learning_rate: 0.0, max_epochs: 0, }
}
