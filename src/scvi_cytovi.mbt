///|
/// CYTOVI: VAE for flow/mass cytometry data with preprocessing.
pub struct CytoViConfig {
  n_latent : Int32
  n_hidden : Int32
  n_markers : Int32
  cofactor : Double
  n_layers : Int32
}

///|
pub fn CytoViConfig::default() -> CytoViConfig {
  { n_latent: 10, n_hidden: 128, n_markers: 0, cofactor: 5.0, n_layers: 2, }
}

///|
pub struct CytoViModel {
  config : CytoViConfig
  n_markers : Int32
  encoder : Array[Array[Double]]
  decoder : Array[Array[Double]]
  mut latent : Array[Array[Double]]
  mut is_trained : Bool
}

///|
pub fn cytovi_init(config : CytoViConfig) -> CytoViModel {
  let n = if config.n_markers > 0 { config.n_markers } else { 10 }
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
    n_markers: n,
    encoder: make_w(config.n_latent, n),
    decoder: make_w(n, config.n_latent),
    latent: [],
    is_trained: false,
  }
}

///|
pub fn cytovi_arcsinh_transform(
  x : Array[Double],
  cofactor : Double,
) -> Array[Double] {
  let result = Array(capacity=x.length())
  for i in 0..<x.length() {
    let val = x[i] / cofactor
    let asinh = @math.ln(val + @math.pow(val * val + 1.0, 0.5))
    result.push(asinh)
  }
  result
}

///|
pub fn cytovi_preprocess(
  model : CytoViModel,
  data : Array[Array[Double]],
) -> Array[Array[Double]] {
  let result = Array(capacity=data.length())
  for i in 0..<data.length() {
    result.push(cytovi_arcsinh_transform(data[i], model.config.cofactor))
  }
  result
}

///|
pub fn cytovi_encode(model : CytoViModel, x : Array[Double]) -> Array[Double] {
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
pub fn cytovi_decode(model : CytoViModel, z : Array[Double]) -> Array[Double] {
  let x_hat = Array(capacity=model.n_markers)
  for m in 0..<model.n_markers {
    let mut sum = 0.0
    for j in 0..<z.length() {
      if j < model.decoder[m].length() {
        sum = sum + z[j] * model.decoder[m][j]
      }
    }
    x_hat.push(sum)
  }
  x_hat
}

///|
pub fn cytovi_train(
  model : CytoViModel,
  data : Array[Array[Double]],
  n_epochs : Int32,
) -> CytoViModel {
  let m = model
  let processed = cytovi_preprocess(m, data)
  let all_z = Array(capacity=processed.length())
  for _epoch in 0..<n_epochs {
    for i in 0..<processed.length() {
      let z = cytovi_encode(m, processed[i])
      let x_hat = cytovi_decode(m, z)
      for g in 0..<m.n_markers {
        if g < processed[i].length() {
          let grad = processed[i][g] - x_hat[g]
          for j in 0..<m.config.n_latent {
            if j < m.decoder[g].length() {
              m.decoder[g][j] = m.decoder[g][j] + 0.001 * grad * z[j]
            }
          }
        }
      }
    }
  }
  for i in 0..<processed.length() {
    all_z.push(cytovi_encode(m, processed[i]))
  }
  m.latent = all_z
  m.is_trained = true
  m
}

///|
pub fn cytovi_get_latent(model : CytoViModel) -> Array[Array[Double]] {
  model.latent
}
