///|
/// POISSONVI: VAE for scATAC-seq data using Poisson observation model.
/// Alternative to PeakVI which uses Bernoulli.
pub struct PoissonViConfig {
  n_latent : Int32
  n_layers : Int32
  n_hidden : Int32
  dropout_rate : Double
  region_factors : Bool
}

///|
pub fn PoissonViConfig::default() -> PoissonViConfig {
  {
    n_latent: 10,
    n_layers: 2,
    n_hidden: 128,
    dropout_rate: 0.1,
    region_factors: true,
  }
}

///|
pub struct PoissonViModel {
  config : PoissonViConfig
  n_regions : Int32
  encoder_weights : Array[Array[Double]]
  decoder_weights : Array[Array[Double]]
  region_bias : Array[Double]
  mut latent : Array[Array[Double]]
  mut is_trained : Bool
}

///|
pub fn poissonvi_init(
  n_regions : Int32,
  config : PoissonViConfig,
) -> PoissonViModel {
  let enc = Array(capacity=config.n_hidden)
  for _i in 0..<config.n_hidden {
    let row = Array(capacity=n_regions)
    for _j in 0..<n_regions {
      row.push(0.01)
    }
    enc.push(row)
  }
  let dec = Array(capacity=n_regions)
  for _i in 0..<n_regions {
    let row = Array(capacity=config.n_latent)
    for _j in 0..<config.n_latent {
      row.push(0.01)
    }
    dec.push(row)
  }
  let bias = Array(capacity=n_regions)
  for _i in 0..<n_regions {
    bias.push(0.0)
  }
  {
    config,
    n_regions,
    encoder_weights: enc,
    decoder_weights: dec,
    region_bias: bias,
    latent: [],
    is_trained: false,
  }
}

///|
pub fn poissonvi_encode(
  model : PoissonViModel,
  x : Array[Double],
) -> Array[Double] {
  let z = Array(capacity=model.config.n_latent)
  for i in 0..<model.config.n_latent {
    let mut sum = 0.0
    for j in 0..<x.length() {
      if j < model.encoder_weights[i].length() {
        sum = sum + x[j] * model.encoder_weights[i][j]
      }
    }
    z.push(sum)
  }
  z
}

///|
pub fn poissonvi_decode(
  model : PoissonViModel,
  z : Array[Double],
) -> Array[Double] {
  let rate = Array(capacity=model.n_regions)
  for r in 0..<model.n_regions {
    let mut sum = model.region_bias[r]
    for j in 0..<z.length() {
      if j < model.decoder_weights[r].length() {
        sum = sum + z[j] * model.decoder_weights[r][j]
      }
    }
    rate.push(@math.exp(sum))
  }
  rate
}

///|
pub fn poissonvi_loglik(
  model : PoissonViModel,
  x : Array[Double],
  z : Array[Double],
) -> Double {
  let rates = poissonvi_decode(model, z)
  let mut ll = 0.0
  for r in 0..<x.length() {
    if r < rates.length() && rates[r] > 0.0 {
      ll = ll + x[r] * @math.ln(rates[r]) - rates[r]
    }
  }
  ll
}

///|
pub fn poissonvi_get_accessibility(
  model : PoissonViModel,
  z : Array[Double],
) -> Array[Double] {
  poissonvi_decode(model, z)
}

///|
pub fn poissonvi_train(
  model : PoissonViModel,
  data : Array[Array[Double]],
  n_epochs : Int32,
) -> PoissonViModel {
  let m = model
  let all_z = Array(capacity=data.length())
  for _epoch in 0..<n_epochs {
    for i in 0..<data.length() {
      let z = poissonvi_encode(m, data[i])
      let rates = poissonvi_decode(m, z)
      for r in 0..<m.n_regions {
        if r < data[i].length() && rates[r] > 0.0 {
          let grad = data[i][r] / rates[r] - 1.0
          for j in 0..<m.config.n_latent {
            if j < m.decoder_weights[r].length() {
              m.decoder_weights[r][j] = m.decoder_weights[r][j] +
                0.001 * grad * z[j] * rates[r]
            }
          }
        }
      }
    }
  }
  for i in 0..<data.length() {
    all_z.push(poissonvi_encode(m, data[i]))
  }
  m.latent = all_z
  m.is_trained = true
  m
}
