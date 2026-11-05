///|
/// SysVI: Systematic variation removal VAE with custom priors for batch correction.
pub struct SysViConfig {
  n_latent : Int32
  n_hidden : Int32
  n_layers : Int32
  prior_mean : Array[Double]
  prior_var : Array[Double]
  beta : Double
}

///|
pub fn SysViConfig::default() -> SysViConfig {
  {
    n_latent: 10,
    n_hidden: 128,
    n_layers: 2,
    prior_mean: [],
    prior_var: [],
    beta: 1.0,
  }
}

///|
pub struct SysViModel {
  config : SysViConfig
  n_genes : Int32
  n_batches : Int32
  encoder : Array[Array[Double]]
  decoder : Array[Array[Double]]
  batch_effect : Array[Array[Double]]
  mut latent : Array[Array[Double]]
  mut is_trained : Bool
}

///|
pub fn sysvi_init(
  n_genes : Int32,
  n_batches : Int32,
  config : SysViConfig,
) -> SysViModel {
  let make_w = (r : Int32, c : Int32) => {
    let w = Array(capacity=r)
    for _i in 0..<r {
      let row = Array(capacity=c)
      for _j in 0..<c {
        row.push(0.01)
      }
      w.push(row)
    }
    w
  }
  {
    config,
    n_genes,
    n_batches,
    encoder: make_w(config.n_latent, n_genes),
    decoder: make_w(n_genes, config.n_latent),
    batch_effect: make_w(n_batches, n_genes),
    latent: [],
    is_trained: false,
  }
}

///|
pub fn sysvi_encode(model : SysViModel, x : Array[Double]) -> Array[Double] {
  let z = Array(capacity=model.config.n_latent)
  for i in 0..<model.config.n_latent {
    let mut sum = 0.0
    for j in 0..<x.length() {
      if j < model.encoder[i].length() {
        sum = sum + x[j] * model.encoder[i][j]
      }
    }
    z.push(sum)
  }
  z
}

///|
pub fn sysvi_decode(
  model : SysViModel,
  z : Array[Double],
  batch : Int32,
) -> Array[Double] {
  let x_hat = Array(capacity=model.n_genes)
  for g in 0..<model.n_genes {
    let mut sum = 0.0
    for j in 0..<z.length() {
      if j < model.decoder[g].length() {
        sum = sum + z[j] * model.decoder[g][j]
      }
    }
    if batch >= 0 &&
      batch < model.batch_effect.length() &&
      g < model.batch_effect[batch].length() {
      sum = sum + model.batch_effect[batch][g]
    }
    x_hat.push(@math.exp(sum))
  }
  x_hat
}

///|
pub fn sysvi_custom_kl(model : SysViModel, z : Array[Double]) -> Double {
  let mut kl = 0.0
  for i in 0..<z.length() {
    let mean = if i < model.config.prior_mean.length() {
      model.config.prior_mean[i]
    } else {
      0.0
    }
    let var_ = if i < model.config.prior_var.length() {
      model.config.prior_var[i]
    } else {
      1.0
    }
    if var_ > 0.0 {
      kl = kl +
        0.5 * ((z[i] - mean) * (z[i] - mean) / var_ + @math.ln(var_) - 1.0)
    }
  }
  kl
}

///|
pub fn sysvi_correct_batch(
  mdl : SysViModel,
  x : Array[Double],
  _batch : Int32,
) -> Array[Double] {
  let z = sysvi_encode(mdl, x)
  sysvi_decode(mdl, z, -1)
}

///|
pub fn sysvi_train(
  model : SysViModel,
  data : Array[Array[Double]],
  batches : Array[Int32],
  n_epochs : Int32,
) -> SysViModel {
  let m = model
  let all_z = Array(capacity=data.length())
  for _epoch in 0..<n_epochs {
    for i in 0..<data.length() {
      let z = sysvi_encode(m, data[i])
      let x_hat = sysvi_decode(m, z, batches[i])
      let kl = sysvi_custom_kl(m, z)
      for g in 0..<m.n_genes {
        if g < data[i].length() && x_hat[g] > 0.0 {
          let grad = data[i][g] / x_hat[g] - 1.0
          for j in 0..<m.config.n_latent {
            if j < m.decoder[g].length() {
              m.decoder[g][j] = m.decoder[g][j] + 0.001 * grad * z[j]
            }
          }
        }
      }
      ignore(kl)
    }
  }
  for i in 0..<data.length() {
    all_z.push(sysvi_encode(m, data[i]))
  }
  m.latent = all_z
  m.is_trained = true
  m
}

///|
pub fn sysvi_get_latent(model : SysViModel) -> Array[Array[Double]] {
  model.latent
}
