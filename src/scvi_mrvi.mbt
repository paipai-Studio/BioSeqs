///|
/// MRVI: Multi-resolution VAE for hierarchical modeling of single-cell data.
pub struct MrviConfig {
  n_latent : Int32
  n_hidden : Int32
  n_levels : Int32
  resolution_factors : Array[Double]
}

///|
pub fn MrviConfig::default() -> MrviConfig {
  {
    n_latent: 10,
    n_hidden: 128,
    n_levels: 3,
    resolution_factors: [1.0, 0.5, 0.25],
  }
}

///|
pub struct MrviModel {
  config : MrviConfig
  n_genes : Int32
  encoders : Array[Array[Array[Double]]]
  decoders : Array[Array[Array[Double]]]
  mut latent : Array[Array[Double]]
  mut is_trained : Bool
}

///|
pub fn mrvi_init(n_genes : Int32, config : MrviConfig) -> MrviModel {
  let encs = Array(capacity=config.n_levels)
  let decs = Array(capacity=config.n_levels)
  for level in 0..<config.n_levels {
    let lat_dim = config.n_latent
    let enc = Array(capacity=lat_dim)
    for _i in 0..<lat_dim {
      let row = Array(capacity=n_genes)
      for _j in 0..<n_genes {
        row.push(0.01 * (level + 1).to_double())
      }
      enc.push(row)
    }
    encs.push(enc)
    let dec = Array(capacity=n_genes)
    for _i in 0..<n_genes {
      let row = Array(capacity=lat_dim)
      for _j in 0..<lat_dim {
        row.push(0.01)
      }
      dec.push(row)
    }
    decs.push(dec)
  }
  {
    config,
    n_genes,
    encoders: encs,
    decoders: decs,
    latent: [],
    is_trained: false,
  }
}

///|
pub fn mrvi_encode_at_level(
  model : MrviModel,
  x : Array[Double],
  level : Int32,
) -> Array[Double] {
  let z = Array(capacity=model.config.n_latent)
  if level < 0 || level >= model.encoders.length() {
    return z
  }
  for i in 0..<model.config.n_latent {
    let mut sum = 0.0
    for j in 0..<x.length() {
      if j < model.encoders[level][i].length() {
        sum = sum + x[j] * model.encoders[level][i][j]
      }
    }
    z.push(sum)
  }
  z
}

///|
pub fn mrvi_decode_at_level(
  model : MrviModel,
  z : Array[Double],
  level : Int32,
) -> Array[Double] {
  let x_hat = Array(capacity=model.n_genes)
  if level < 0 || level >= model.decoders.length() {
    return x_hat
  }
  for g in 0..<model.n_genes {
    let mut sum = 0.0
    for j in 0..<z.length() {
      if j < model.decoders[level][g].length() {
        sum = sum + z[j] * model.decoders[level][g][j]
      }
    }
    x_hat.push(@math.exp(sum))
  }
  x_hat
}

///|
pub fn mrvi_encode_multi_resolution(
  model : MrviModel,
  x : Array[Double],
) -> Array[Array[Double]] {
  let result = Array(capacity=model.config.n_levels)
  for level in 0..<model.config.n_levels {
    result.push(mrvi_encode_at_level(model, x, level))
  }
  result
}

///|
pub fn mrvi_coarse_to_fine(
  model : MrviModel,
  z_coarse : Array[Double],
  level_coarse : Int32,
  level_fine : Int32,
) -> Array[Double] {
  let z_fine = Array(capacity=z_coarse.length())
  let scale = if level_coarse < model.config.resolution_factors.length() &&
    level_fine < model.config.resolution_factors.length() {
    model.config.resolution_factors[level_fine] /
    model.config.resolution_factors[level_coarse]
  } else {
    1.0
  }
  for i in 0..<z_coarse.length() {
    z_fine.push(z_coarse[i] * scale)
  }
  z_fine
}

///|
pub fn mrvi_train(
  model : MrviModel,
  data : Array[Array[Double]],
  n_epochs : Int32,
) -> MrviModel {
  let m = model
  let all_z = Array(capacity=data.length())
  for _epoch in 0..<n_epochs {
    for i in 0..<data.length() {
      for level in 0..<m.config.n_levels {
        let z = mrvi_encode_at_level(m, data[i], level)
        let x_hat = mrvi_decode_at_level(m, z, level)
        for g in 0..<m.n_genes {
          if g < data[i].length() && x_hat[g] > 0.0 {
            let grad = data[i][g] / x_hat[g] - 1.0
            for j in 0..<m.config.n_latent {
              if j < m.decoders[level][g].length() {
                m.decoders[level][g][j] = m.decoders[level][g][j] +
                  0.001 * grad * z[j]
              }
            }
          }
        }
      }
    }
  }
  for i in 0..<data.length() {
    all_z.push(mrvi_encode_at_level(m, data[i], 0))
  }
  m.latent = all_z
  m.is_trained = true
  m
}

///|
pub fn mrvi_get_latent(model : MrviModel) -> Array[Array[Double]] {
  model.latent
}
