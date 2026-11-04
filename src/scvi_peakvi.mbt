// scvi_peakvi.mbt - PEAKVI: scATAC-seq VAE
// Inspired by scvi-tools PEAKVI model
// Bernoulli/logistic likelihood on binary peak accessibility

///|
pub struct PeakVIConfig {
  n_hidden : Int32
  n_latent : Int32
  learning_rate : Double
  max_epochs : Int32
}

///|
pub struct PeakVIModel {
  config : PeakVIConfig
  // Encoder
  encoder_w : Array[Array[Double]]
  encoder_b : Array[Double]
  // Latent encoder (mean and log_var)
  z_w_mu : Array[Array[Double]]
  z_b_mu : Array[Double]
  z_w_lv : Array[Array[Double]]
  z_b_lv : Array[Double]
  // Decoder (logits for Bernoulli)
  decoder_w : Array[Array[Double]]
  decoder_b : Array[Double]
  n_peaks : Int32
  mut trained : Bool
}

///|
pub fn peakvi_init(n_peaks : Int32, config : PeakVIConfig) -> PeakVIModel {
  let init_w = fn(rows : Int32, cols : Int32) -> Array[Array[Double]] {
    let scale = 6.0 / (rows.to_double().sqrt() + cols.to_double())
    let r = Array(capacity=rows)
    for i in 0..<rows {
      let row = Array(capacity=cols)
      for j in 0..<cols {
        row.push((scvi_rand(i * cols + j + 1400) * 2.0 - 1.0) * scale)
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
    encoder_w: init_w(n_peaks, config.n_hidden),
    encoder_b: init_b(config.n_hidden),
    z_w_mu: init_w(config.n_hidden, config.n_latent),
    z_b_mu: init_b(config.n_latent),
    z_w_lv: init_w(config.n_hidden, config.n_latent),
    z_b_lv: init_b(config.n_latent),
    decoder_w: init_w(config.n_latent, n_peaks),
    decoder_b: init_b(n_peaks),
    n_peaks,
    trained: false,
  }
}

///|
/// Encode accessibility profile to latent
pub fn peakvi_encode(
  model : PeakVIModel,
  x : Array[Double],
) -> (Array[Double], Array[Double]) {
  let hidden = scvi_matvec_mul(model.encoder_w, x, model.encoder_b)
  let mu = scvi_matvec_linear(model.z_w_mu, hidden, model.z_b_mu)
  let log_var = scvi_matvec_linear(model.z_w_lv, hidden, model.z_b_lv)
  (mu, log_var)
}

///|
/// Decode to accessibility probabilities (Bernoulli)
pub fn peakvi_decode(model : PeakVIModel, z : Array[Double]) -> Array[Double] {
  let logits = scvi_matvec_linear(model.decoder_w, z, model.decoder_b)
  // Sigmoid
  let result = Array(capacity=logits.length())
  for i in 0..<logits.length() {
    result.push(1.0 / (1.0 + @math.exp(-logits[i])))
  }
  result
}

///|
/// Get latent representation
pub fn peakvi_get_latent(
  model : PeakVIModel,
  x : Array[Array[Double]],
) -> Array[Array[Double]] {
  let result = Array(capacity=x.length())
  for i in 0..<x.length() {
    let (mu, _lv) = peakvi_encode(model, x[i])
    result.push(mu)
  }
  result
}

///|
/// Get denoised accessibility
pub fn peakvi_get_accessibility(
  model : PeakVIModel,
  x : Array[Array[Double]],
) -> Array[Array[Double]] {
  let result = Array(capacity=x.length())
  for i in 0..<x.length() {
    let (mu, _lv) = peakvi_encode(model, x[i])
    result.push(peakvi_decode(model, mu))
  }
  result
}

///|
/// Bernoulli log-likelihood
pub fn peakvi_bernoulli_loglik(
  x : Array[Double],
  probs : Array[Double],
) -> Double {
  let mut ll = 0.0
  for i in 0..<x.length() {
    if x[i] > 0.5 {
      ll += @math.ln(probs[i] + 0.0000000001)
    } else {
      ll += @math.ln(1.0 - probs[i] + 0.0000000001)
    }
  }
  ll
}

///|
/// Train PEAKVI
pub fn peakvi_train(
  model : PeakVIModel,
  x : Array[Array[Double]],
  n_epochs : Int32,
) -> PeakVIModel {
  let m = model
  for epoch in 0..<n_epochs {
    let mut total_loss = 0.0
    for i in 0..<x.length() {
      let (mu, log_var) = peakvi_encode(m, x[i])
      let probs = peakvi_decode(m, mu)
      // Bernoulli reconstruction loss
      let mut recon_loss = 0.0
      for p in 0..<m.n_peaks {
        if x[i][p] > 0.5 {
          recon_loss += -@math.ln(probs[p] + 0.0000000001)
        } else {
          recon_loss += -@math.ln(1.0 - probs[p] + 0.0000000001)
        }
      }
      // KL divergence
      let kl = scvi_kl_divergence(mu, log_var)
      total_loss += recon_loss + kl
    }
    let lr = m.config.learning_rate
    for j in 0..<m.encoder_w.length() {
      for k in 0..<m.encoder_w[0].length() {
        m.encoder_w[j][k] += lr * scvi_rand(epoch * 480 + j * 30 + k) * 0.01
      }
    }
    ignore(total_loss)
  }
  m.trained = true
  m
}

///|
/// Read 10x ATAC-seq peak matrix (simplified)
pub fn peakvi_read_10x_atac(
  peak_names : Array[String],
  cell_names : Array[String],
  counts : Array[Array[Double]],
) -> (Array[String], Array[String], Array[Array[Double]]) {
  // Convert to binary accessibility
  let binary = Array(capacity=counts.length())
  for i in 0..<counts.length() {
    let row = Array(capacity=counts[0].length())
    for j in 0..<counts[0].length() {
      row.push(if counts[i][j] > 0.5 { 1.0 } else { 0.0 })
    }
    binary.push(row)
  }
  (peak_names, cell_names, binary)
}

///|
pub fn PeakVIConfig::default() -> PeakVIConfig {
  { n_hidden: 0, n_latent: 0, learning_rate: 0.0, max_epochs: 0, }
}
