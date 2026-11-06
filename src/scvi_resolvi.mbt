///|
/// RESOLVI: Resolution improvement VAE for enhancing spatial transcriptomics resolution.
pub struct ResolViConfig {
  n_latent : Int32
  n_hidden : Int32
  upsample_factor : Int32
  n_layers : Int32
}

///|
pub fn ResolViConfig::default() -> ResolViConfig {
  { n_latent: 10, n_hidden: 128, upsample_factor: 2, n_layers: 2, }
}

///|
pub struct ResolViModel {
  config : ResolViConfig
  n_genes : Int32
  n_spots : Int32
  encoder : Array[Array[Double]]
  decoder : Array[Array[Double]]
  spatial_coords : Array[(Double, Double)]
  mut latent : Array[Array[Double]]
  mut is_trained : Bool
}

///|
pub fn resolvi_init(
  n_genes : Int32,
  n_spots : Int32,
  spatial_coords : Array[(Double, Double)],
  config : ResolViConfig,
) -> ResolViModel {
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
    n_spots,
    spatial_coords,
    encoder: make_w(config.n_latent, n_genes),
    decoder: make_w(n_genes, config.n_latent),
    latent: [],
    is_trained: false,
  }
}

///|
pub fn resolvi_encode(model : ResolViModel, x : Array[Double]) -> Array[Double] {
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
pub fn resolvi_decode(model : ResolViModel, z : Array[Double]) -> Array[Double] {
  let x_hat = Array(capacity=model.n_genes)
  for g in 0..<model.n_genes {
    let mut sum = 0.0
    for j in 0..<z.length() {
      if j < model.decoder[g].length() {
        sum = sum + z[j] * model.decoder[g][j]
      }
    }
    x_hat.push(@math.exp(sum))
  }
  x_hat
}

///|
pub fn resolvi_upsample_spatial(
  model : ResolViModel,
) -> Array[(Double, Double)] {
  let upsampled = Array(
    capacity=model.n_spots *
      model.config.upsample_factor *
      model.config.upsample_factor,
  )
  let f = model.config.upsample_factor.to_double()
  for i in 0..<model.n_spots {
    let (x, y) = model.spatial_coords[i]
    for dx in 0..<model.config.upsample_factor {
      for dy in 0..<model.config.upsample_factor {
        let nx = x + dx.to_double() / f - 0.5 / f
        let ny = y + dy.to_double() / f - 0.5 / f
        upsampled.push((nx, ny))
      }
    }
  }
  upsampled
}

///|
pub fn resolvi_predict_subspot(
  model : ResolViModel,
  neighbors : Array[Array[Double]],
  weights : Array[Double],
) -> Array[Double] {
  let result = Array(capacity=model.n_genes)
  for g in 0..<model.n_genes {
    let mut sum = 0.0
    let mut wsum = 0.0
    for n in 0..<neighbors.length() {
      if g < neighbors[n].length() && n < weights.length() {
        sum = sum + neighbors[n][g] * weights[n]
        wsum = wsum + weights[n]
      }
    }
    if wsum > 0.0 {
      result.push(sum / wsum)
    } else {
      result.push(0.0)
    }
  }
  result
}

///|
pub fn resolvi_train(
  model : ResolViModel,
  data : Array[Array[Double]],
  n_epochs : Int32,
) -> ResolViModel {
  let m = model
  let all_z = Array(capacity=data.length())
  for _epoch in 0..<n_epochs {
    for i in 0..<data.length() {
      let z = resolvi_encode(m, data[i])
      let x_hat = resolvi_decode(m, z)
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
    }
  }
  for i in 0..<data.length() {
    all_z.push(resolvi_encode(m, data[i]))
  }
  m.latent = all_z
  m.is_trained = true
  m
}
