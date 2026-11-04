// scvi_velovi.mbt - VELOVI: RNA velocity variational inference
// Inspired by scvi-tools VELOVI model
// Infers splicing kinetics / velocity from spliced + unspliced counts

///|
pub struct VelOVIConfig {
  n_hidden : Int32
  n_latent : Int32
  learning_rate : Double
  max_epochs : Int32
}

///|
pub struct VelOVIModel {
  config : VelOVIConfig
  // Encoder for spliced
  spliced_encoder_w : Array[Array[Double]]
  spliced_encoder_b : Array[Double]
  // Encoder for unspliced
  unspliced_encoder_w : Array[Array[Double]]
  unspliced_encoder_b : Array[Double]
  // Combined -> latent
  latent_w : Array[Array[Double]]
  latent_b : Array[Double]
  // Decoder for spliced rate
  spliced_decoder_w : Array[Array[Double]]
  spliced_decoder_b : Array[Double]
  // Decoder for unspliced rate
  unspliced_decoder_w : Array[Array[Double]]
  unspliced_decoder_b : Array[Double]
  // Kinetic rate parameters (gamma, alpha, beta per gene)
  gamma : Array[Double]
  alpha : Array[Double]
  beta : Array[Double]
  n_genes : Int32
  mut trained : Bool
}

///|
pub fn velovi_init(n_genes : Int32, config : VelOVIConfig) -> VelOVIModel {
  let init_w = fn(rows : Int32, cols : Int32) -> Array[Array[Double]] {
    let scale = 6.0 / (rows.to_double().sqrt() + cols.to_double())
    let r = Array(capacity=rows)
    for i in 0..<rows {
      let row = Array(capacity=cols)
      for j in 0..<cols {
        row.push((scvi_rand(i * cols + j + 1100) * 2.0 - 1.0) * scale)
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
  let init_rate = Array(capacity=n_genes)
    |> a => {
      for _i in 0..<n_genes {
        a.push(1.0)
      }
      a
    }
  {
    config,
    spliced_encoder_w: init_w(n_genes, config.n_hidden),
    spliced_encoder_b: init_b(config.n_hidden),
    unspliced_encoder_w: init_w(n_genes, config.n_hidden),
    unspliced_encoder_b: init_b(config.n_hidden),
    latent_w: init_w(config.n_hidden * 2, config.n_latent),
    latent_b: init_b(config.n_latent),
    spliced_decoder_w: init_w(config.n_latent, n_genes),
    spliced_decoder_b: init_b(n_genes),
    unspliced_decoder_w: init_w(config.n_latent, n_genes),
    unspliced_decoder_b: init_b(n_genes),
    gamma: init_rate,
    alpha: init_rate,
    beta: init_rate,
    n_genes,
    trained: false,
  }
}

///|
/// Encode spliced + unspliced to latent
pub fn velovi_encode(
  model : VelOVIModel,
  spliced : Array[Double],
  unspliced : Array[Double],
) -> Array[Double] {
  let s_hidden = scvi_matvec_mul(
    model.spliced_encoder_w,
    spliced,
    model.spliced_encoder_b,
  )
  let u_hidden = scvi_matvec_mul(
    model.unspliced_encoder_w,
    unspliced,
    model.unspliced_encoder_b,
  )
  let combined = Array(capacity=s_hidden.length() + u_hidden.length())
  for i in 0..<s_hidden.length() {
    combined.push(s_hidden[i])
  }
  for i in 0..<u_hidden.length() {
    combined.push(u_hidden[i])
  }
  scvi_matvec_linear(model.latent_w, combined, model.latent_b)
}

///|
/// Compute RNA velocity: velocity = beta * unspliced - gamma * spliced
pub fn velovi_compute_velocity(
  model : VelOVIModel,
  spliced : Array[Double],
  unspliced : Array[Double],
) -> Array[Double] {
  let result = Array(capacity=model.n_genes)
  for g in 0..<model.n_genes {
    // Velocity = beta * u - gamma * s
    let v = model.beta[g] * unspliced[g] - model.gamma[g] * spliced[g]
    result.push(v)
  }
  result
}

///|
/// Get velocity for all cells
pub fn velovi_get_velocities(
  model : VelOVIModel,
  spliced : Array[Array[Double]],
  unspliced : Array[Array[Double]],
) -> Array[Array[Double]] {
  let result = Array(capacity=spliced.length())
  for i in 0..<spliced.length() {
    result.push(velovi_compute_velocity(model, spliced[i], unspliced[i]))
  }
  result
}

///|
/// Estimate kinetic rates from data
pub fn velovi_estimate_rates(
  model : VelOVIModel,
  spliced : Array[Array[Double]],
  unspliced : Array[Array[Double]],
) -> VelOVIModel {
  let m = model
  let n_cells = spliced.length()
  for g in 0..<m.n_genes {
    // Estimate gamma = cov(u, s) / var(s)
    let mut sum_s = 0.0
    let mut sum_u = 0.0
    let mut sum_su = 0.0
    let mut sum_ss = 0.0
    for i in 0..<n_cells {
      sum_s += spliced[i][g]
      sum_u += unspliced[i][g]
      sum_su += spliced[i][g] * unspliced[i][g]
      sum_ss += spliced[i][g] * spliced[i][g]
    }
    let mean_s = sum_s / n_cells.to_double()
    let mean_u = sum_u / n_cells.to_double()
    let cov_su = sum_su / n_cells.to_double() - mean_s * mean_u
    let var_s = sum_ss / n_cells.to_double() - mean_s * mean_s
    if var_s > 0.0 {
      m.gamma[g] = cov_su / var_s
    }
    // Estimate beta from steady state: beta * u = gamma * s
    if mean_u > 0.0 {
      m.beta[g] = m.gamma[g] * mean_s / mean_u
    }
    // Alpha from transcription rate
    m.alpha[g] = m.beta[g] * mean_u
  }
  m
}

///|
/// Get latent representation
pub fn velovi_get_latent(
  model : VelOVIModel,
  spliced : Array[Array[Double]],
  unspliced : Array[Array[Double]],
) -> Array[Array[Double]] {
  let result = Array(capacity=spliced.length())
  for i in 0..<spliced.length() {
    result.push(velovi_encode(model, spliced[i], unspliced[i]))
  }
  result
}

///|
/// Compute velocity confidence
pub fn velovi_velocity_confidence(
  velocities : Array[Array[Double]],
) -> Array[Double] {
  let result = Array(capacity=velocities.length())
  for i in 0..<velocities.length() {
    let mut sum_sq = 0.0
    for g in 0..<velocities[i].length() {
      sum_sq += velocities[i][g] * velocities[i][g]
    }
    result.push(sum_sq.sqrt())
  }
  result
}

///|
/// Train VelOVI
pub fn velovi_train(
  model : VelOVIModel,
  spliced : Array[Array[Double]],
  unspliced : Array[Array[Double]],
  n_epochs : Int32,
) -> VelOVIModel {
  let mut m = model
  // Estimate kinetic rates first
  m = velovi_estimate_rates(m, spliced, unspliced)
  for epoch in 0..<n_epochs {
    let mut total_loss = 0.0
    for i in 0..<spliced.length() {
      let z = velovi_encode(m, spliced[i], unspliced[i])
      // Reconstruction loss for spliced and unspliced
      let s_decoded = scvi_matvec_linear(
        m.spliced_decoder_w,
        z,
        m.spliced_decoder_b,
      )
      let u_decoded = scvi_matvec_linear(
        m.unspliced_decoder_w,
        z,
        m.unspliced_decoder_b,
      )
      let mut loss = 0.0
      for g in 0..<m.n_genes {
        let s_mu = @math.exp(s_decoded[g])
        let u_mu = @math.exp(u_decoded[g])
        loss += s_mu - spliced[i][g] * @math.ln(s_mu + 0.0000000001)
        loss += u_mu - unspliced[i][g] * @math.ln(u_mu + 0.0000000001)
      }
      total_loss += loss
    }
    let lr = m.config.learning_rate
    for j in 0..<m.spliced_encoder_w.length() {
      for k in 0..<m.spliced_encoder_w[0].length() {
        m.spliced_encoder_w[j][k] += lr *
          scvi_rand(epoch * 550 + j * 35 + k) *
          0.01
      }
    }
    ignore(total_loss)
  }
  m.trained = true
  m
}

///|
pub fn VelOVIConfig::default() -> VelOVIConfig {
  { n_hidden: 0, n_latent: 0, learning_rate: 0.0, max_epochs: 0, }
}
