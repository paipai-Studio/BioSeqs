// scvi_vae.mbt - scVI core: Variational Autoencoder for single-cell RNA-seq
// Inspired by scvi-tools SCVI model (Lopez et al. 2018)
// Provides dimensionality reduction, batch correction, and denoising

///|
pub struct ScVIConfig {
  n_hidden : Int32
  n_latent : Int32
  n_layers : Int32
  dropout_rate : Double
  dispersion : String
  gene_likelihood : String
  latent_distribution : String
  learning_rate : Double
  max_epochs : Int32
}

///|
pub struct ScVIModel {
  config : ScVIConfig
  // Encoder weights: input -> hidden
  encoder_w1 : Array[Array[Double]]
  encoder_b1 : Array[Double]
  // Encoder -> latent mean
  encoder_w_mu : Array[Array[Double]]
  encoder_b_mu : Array[Double]
  // Encoder -> latent log_var
  encoder_w_lv : Array[Array[Double]]
  encoder_b_lv : Array[Double]
  // Decoder weights: latent -> hidden
  decoder_w1 : Array[Array[Double]]
  decoder_b1 : Array[Double]
  // Decoder -> output (gene dispersion rates)
  decoder_w_out : Array[Array[Double]]
  decoder_b_out : Array[Double]
  // Gene-specific dispersion (log theta)
  log_theta : Array[Double]
  // Library size log mean/let mut per batch
  library_log_mean : Array[Double]
  library_log_var : Array[Double]
  n_batch : Int32
  n_genes : Int32
  mut trained : Bool
}

///|
pub fn ScVIConfig::default() -> ScVIConfig {
  {
    n_hidden: 128,
    n_latent: 10,
    n_layers: 1,
    dropout_rate: 0.1,
    dispersion: "gene",
    gene_likelihood: "zinb",
    latent_distribution: "normal",
    learning_rate: 0.001,
    max_epochs: 100,
  }
}

///|
pub fn scvi_init_model(
  n_genes : Int32,
  n_batch : Int32,
  config : ScVIConfig,
) -> ScVIModel {
  let rng_seed = 42
  // Xavier initialization
  let init_w = fn(rows : Int32, cols : Int32) -> Array[Array[Double]] {
    let scale = 6.0 / (rows.to_double().sqrt() + cols.to_double())
    let result = Array(capacity=rows)
    for i in 0..<rows {
      let row = Array(capacity=cols)
      for j in 0..<cols {
        let val = scvi_rand(rng_seed + i * cols + j) * 2.0 * scale - scale
        row.push(val)
      }
      result.push(row)
    }
    result
  }
  let init_b = fn(n : Int32) -> Array[Double] {
    let result = Array(capacity=n)
    for _i in 0..<n {
      result.push(0.0)
    }
    result
  }
  {
    config,
    encoder_w1: init_w(config.n_hidden, n_genes),
    encoder_b1: init_b(config.n_hidden),
    encoder_w_mu: init_w(config.n_latent, config.n_hidden),
    encoder_b_mu: init_b(config.n_latent),
    encoder_w_lv: init_w(config.n_latent, config.n_hidden),
    encoder_b_lv: init_b(config.n_latent),
    decoder_w1: init_w(config.n_hidden, config.n_latent),
    decoder_b1: init_b(config.n_hidden),
    decoder_w_out: init_w(n_genes, config.n_hidden),
    decoder_b_out: init_b(n_genes),
    log_theta: init_b(n_genes),
    library_log_mean: init_b(n_batch),
    library_log_var: init_b(n_batch),
    n_batch,
    n_genes,
    trained: false,
  }
}

///|
fn scvi_rand(seed : Int32) -> Double {
  let s = (seed * 1103515245 + 12345).reinterpret_as_uint()
  (s % 10000).to_double() / 10000.0
}

///|
fn scvi_relu(x : Double) -> Double {
  if x > 0.0 {
    x
  } else {
    0.0
  }
}

///|
fn scvi_matvec_mul(
  weights : Array[Array[Double]],
  input : Array[Double],
  bias : Array[Double],
) -> Array[Double] {
  let n_out = weights.length()
  let result = Array(capacity=n_out)
  for i in 0..<n_out {
    let mut sum = bias[i]
    let row = weights[i]
    for j in 0..<row.length() {
      sum += row[j] * input[j]
    }
    result.push(scvi_relu(sum))
  }
  result
}

///|
fn scvi_matvec_linear(
  weights : Array[Array[Double]],
  input : Array[Double],
  bias : Array[Double],
) -> Array[Double] {
  let n_out = weights.length()
  let result = Array(capacity=n_out)
  for i in 0..<n_out {
    let mut sum = bias[i]
    let row = weights[i]
    for j in 0..<row.length() {
      sum += row[j] * input[j]
    }
    result.push(sum)
  }
  result
}

///|
/// Encode a single cell to latent space (mean and log_var)
pub fn scvi_encode(
  model : ScVIModel,
  x : Array[Double],
) -> (Array[Double], Array[Double]) {
  let hidden = scvi_matvec_mul(model.encoder_w1, x, model.encoder_b1)
  let mu = scvi_matvec_linear(model.encoder_w_mu, hidden, model.encoder_b_mu)
  let log_var = scvi_matvec_linear(
    model.encoder_w_lv,
    hidden,
    model.encoder_b_lv,
  )
  (mu, log_var)
}

///|
/// Decode latent to gene expression rates
pub fn scvi_decode(model : ScVIModel, z : Array[Double]) -> Array[Double] {
  let hidden = scvi_matvec_mul(model.decoder_w1, z, model.decoder_b1)
  let raw = scvi_matvec_linear(model.decoder_w_out, hidden, model.decoder_b_out)
  // Softmax to get normalized rates
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
fn scvi_array_max(arr : Array[Double]) -> Double {
  let mut m = arr[0]
  for i in 1..<arr.length() {
    if arr[i] > m {
      m = arr[i]
    }
  }
  m
}

///|
/// Get latent representation for all cells
pub fn scvi_get_latent_representation(
  model : ScVIModel,
  x : Array[Array[Double]],
) -> Array[Array[Double]] {
  let n_cells = x.length()
  let result = Array(capacity=n_cells)
  for i in 0..<n_cells {
    let (mu, _log_var) = scvi_encode(model, x[i])
    result.push(mu)
  }
  result
}

///|
/// Get normalized (denoised) expression
pub fn scvi_get_normalized_expression(
  model : ScVIModel,
  x : Array[Array[Double]],
) -> Array[Array[Double]] {
  let n_cells = x.length()
  let result = Array(capacity=n_cells)
  for i in 0..<n_cells {
    let (mu, log_var) = scvi_encode(model, x[i])
    // Sample z from N(mu, exp(log_var))
    let z = Array(capacity=mu.length())
    for j in 0..<mu.length() {
      let std = @math.exp(log_var[j] * 0.5)
      z.push(mu[j] + std * scvi_rand(i * 100 + j))
    }
    let rates = scvi_decode(model, z)
    result.push(rates)
  }
  result
}

///|
/// Compute reconstruction loss (negative log-likelihood)
pub fn scvi_reconstruction_loss(
  model : ScVIModel,
  x : Array[Double],
  z : Array[Double],
) -> Double {
  let rates = scvi_decode(model, z)
  let library = scvi_array_sum(x)
  let mut loss = 0.0
  for g in 0..<model.n_genes {
    let mu = rates[g] * library
    let theta = @math.exp(model.log_theta[g])
    // Negative binomial log-likelihood (simplified)
    let x_g = x[g]
    if x_g > 0.0 {
      loss += -(x_g * @math.ln(mu + 0.0000000001) -
        mu -
        theta * @math.ln(theta + mu + 0.0000000001) +
        x_g * @math.ln(theta + mu + 0.0000000001))
    } else {
      loss += -(-mu - theta * @math.ln(theta + mu + 0.0000000001))
    }
  }
  loss
}

///|
fn scvi_array_sum(arr : Array[Double]) -> Double {
  let mut s = 0.0
  for i in 0..<arr.length() {
    s += arr[i]
  }
  s
}

///|
/// KL divergence between latent posterior and prior N(0, I)
pub fn scvi_kl_divergence(
  mu : Array[Double],
  log_var : Array[Double],
) -> Double {
  let mut kl = 0.0
  for i in 0..<mu.length() {
    kl += 0.5 * (log_var[i] - @math.exp(log_var[i]) - mu[i] * mu[i] + 1.0)
  }
  kl
}

///|
/// Compute ELBO (Evidence Lower Bound) for a cell
pub fn scvi_compute_elbo(model : ScVIModel, x : Array[Double]) -> Double {
  let (mu, log_var) = scvi_encode(model, x)
  let z = mu // Use mean for deterministic estimate
  let recon = scvi_reconstruction_loss(model, x, z)
  let kl = scvi_kl_divergence(mu, log_var)
  -(recon + kl)
}

///|
/// Train the VAE on count data (simplified gradient-free optimization)
pub fn scvi_train(
  model : ScVIModel,
  x : Array[Array[Double]],
  n_epochs : Int32,
) -> ScVIModel {
  let m = model
  if x.length() == 0 {
    m.trained = true
    return m
  }
  // Gene-wise inverse dispersion is kept fixed during this training pass.
  let theta = Array(capacity=m.n_genes)
  for g in 0..<m.n_genes {
    theta.push(@math.exp(m.log_theta[g]))
  }
  // Allocate Adam moments with shapes matching the model parameters.
  let template = ScviAgTape::new()
  ignore(template.matrix_param(m.encoder_w1))
  ignore(template.vector_param(m.encoder_b1))
  ignore(template.matrix_param(m.encoder_w_mu))
  ignore(template.vector_param(m.encoder_b_mu))
  ignore(template.matrix_param(m.encoder_w_lv))
  ignore(template.vector_param(m.encoder_b_lv))
  ignore(template.matrix_param(m.decoder_w1))
  ignore(template.vector_param(m.decoder_b1))
  ignore(template.matrix_param(m.decoder_w_out))
  ignore(template.vector_param(m.decoder_b_out))
  let adam = scvi_adam_init(template, m.config.learning_rate)
  let n_cells = x.length().to_double()
  for _epoch in 0..<n_epochs {
    let t = ScviAgTape::new()
    let w1 = t.matrix_param(m.encoder_w1)
    let b1 = t.vector_param(m.encoder_b1)
    let w_mu = t.matrix_param(m.encoder_w_mu)
    let b_mu = t.vector_param(m.encoder_b_mu)
    let w_lv = t.matrix_param(m.encoder_w_lv)
    let b_lv = t.vector_param(m.encoder_b_lv)
    let wd1 = t.matrix_param(m.decoder_w1)
    let bd1 = t.vector_param(m.decoder_b1)
    let w_out = t.matrix_param(m.decoder_w_out)
    let b_out = t.vector_param(m.decoder_b_out)
    let mut total = t.leaf([0.0])
    for i in 0..<x.length() {
      let xi = x[i]
      let lib = scvi_array_sum(xi)
      let lib_size = if lib <= 0.0 { 1.0 } else { lib }
      let xin = t.leaf(xi)
      let hidden = t.relu(t.dense(xin, w1, b1))
      let mu = t.dense(hidden, w_mu, b_mu)
      let log_var = t.dense(hidden, w_lv, b_lv)
      let decoded = t.relu(t.dense(mu, wd1, bd1))
      let raw = t.dense(decoded, w_out, b_out)
      let rates = t.softmax(raw)
      let mean_counts = t.scale(rates, lib_size)
      let nll = t.nb_nll(mean_counts, xi, theta)
      let kl = t.kl_divergence(mu, log_var)
      total = t.add(total, t.add(nll, kl))
    }
    let loss = t.scale(total, 1.0 / n_cells)
    t.backward(loss)
    ignore(t.clip_grad_norm(5.0))
    scvi_adam_step(t, adam)
  }
  m.trained = true
  m
}

///|
/// Get marginal log-likelihood estimate
pub fn scvi_get_marginal_ll(
  model : ScVIModel,
  x : Array[Array[Double]],
  n_samples : Int32,
) -> Double {
  let mut total_ll = 0.0
  for i in 0..<x.length() {
    let mut sample_ll = 0.0
    for _s in 0..<n_samples {
      sample_ll += scvi_compute_elbo(model, x[i])
    }
    total_ll += sample_ll / n_samples.to_double()
  }
  total_ll
}
