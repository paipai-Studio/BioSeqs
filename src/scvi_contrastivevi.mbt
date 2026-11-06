///|
/// ContrastiveVI: Contrastive learning VAE for distinguishing foreground vs background.
pub struct ContrastiveViConfig {
  n_latent : Int32
  n_hidden : Int32
  n_layers : Int32
  beta : Double
  temperature : Double
}

///|
pub fn ContrastiveViConfig::default() -> ContrastiveViConfig {
  { n_latent: 10, n_hidden: 128, n_layers: 2, beta: 1.0, temperature: 0.5, }
}

///|
pub struct ContrastiveViModel {
  config : ContrastiveViConfig
  n_features : Int32
  encoder_fg : Array[Array[Double]]
  encoder_bg : Array[Array[Double]]
  decoder : Array[Array[Double]]
  mut latent : Array[Array[Double]]
  mut is_trained : Bool
}

///|
pub fn contrastivevi_init(
  n_features : Int32,
  config : ContrastiveViConfig,
) -> ContrastiveViModel {
  let make_weights = (rows : Int32, cols : Int32) => {
    let w = Array(capacity=rows)
    for _i in 0..<rows {
      let row = Array(capacity=cols)
      for _j in 0..<cols {
        row.push(0.01)
      }
      w.push(row)
    }
    w
  }
  {
    config,
    n_features,
    encoder_fg: make_weights(config.n_latent, n_features),
    encoder_bg: make_weights(config.n_latent, n_features),
    decoder: make_weights(n_features, config.n_latent),
    latent: [],
    is_trained: false,
  }
}

///|
pub fn contrastivevi_encode_fg(
  model : ContrastiveViModel,
  x : Array[Double],
) -> Array[Double] {
  let z = Array(capacity=model.config.n_latent)
  for i in 0..<model.config.n_latent {
    let mut sum = 0.0
    for j in 0..<x.length() {
      if j < model.encoder_fg[i].length() {
        sum = sum + x[j] * model.encoder_fg[i][j]
      }
    }
    z.push(sum)
  }
  z
}

///|
pub fn contrastivevi_encode_bg(
  model : ContrastiveViModel,
  x : Array[Double],
) -> Array[Double] {
  let z = Array(capacity=model.config.n_latent)
  for i in 0..<model.config.n_latent {
    let mut sum = 0.0
    for j in 0..<x.length() {
      if j < model.encoder_bg[i].length() {
        sum = sum + x[j] * model.encoder_bg[i][j]
      }
    }
    z.push(sum)
  }
  z
}

///|
pub fn contrastivevi_decode(
  model : ContrastiveViModel,
  z : Array[Double],
) -> Array[Double] {
  let x_hat = Array(capacity=model.n_features)
  for i in 0..<model.n_features {
    let mut sum = 0.0
    for j in 0..<z.length() {
      if j < model.decoder[i].length() {
        sum = sum + z[j] * model.decoder[i][j]
      }
    }
    x_hat.push(@math.exp(sum))
  }
  x_hat
}

///|
pub fn contrastivevi_contrastive_loss(
  model : ContrastiveViModel,
  z_fg : Array[Double],
  z_bg : Array[Double],
) -> Double {
  let mut sim = 0.0
  for i in 0..<z_fg.length() {
    if i < z_bg.length() {
      sim = sim + z_fg[i] * z_bg[i]
    }
  }
  sim = sim / model.config.temperature
  let mut norm_fg = 0.0
  for i in 0..<z_fg.length() {
    norm_fg = norm_fg + z_fg[i] * z_fg[i]
  }
  let mut norm_bg = 0.0
  for i in 0..<z_bg.length() {
    norm_bg = norm_bg + z_bg[i] * z_bg[i]
  }
  let denom = @math.pow(norm_fg * norm_bg, 0.5)
  if denom > 0.0 {
    sim = sim / denom
  }
  -sim
}

///|
pub fn contrastivevi_train(
  model : ContrastiveViModel,
  data : Array[Array[Double]],
  is_foreground : Array[Bool],
  n_epochs : Int32,
) -> ContrastiveViModel {
  let m = model
  let all_z = Array(capacity=data.length())
  for _epoch in 0..<n_epochs {
    for i in 0..<data.length() {
      let z_fg = contrastivevi_encode_fg(m, data[i])
      let _z_bg = contrastivevi_encode_bg(m, data[i])
      let x_hat = contrastivevi_decode(m, z_fg)
      for j in 0..<m.n_features {
        if j < data[i].length() && x_hat[j] > 0.0 {
          let grad = data[i][j] / x_hat[j] - 1.0
          for k in 0..<m.config.n_latent {
            if k < m.decoder[j].length() {
              m.decoder[j][k] = m.decoder[j][k] + 0.001 * grad * z_fg[k]
            }
          }
        }
      }
      if is_foreground[i] {
        for k in 0..<m.config.n_latent {
          m.encoder_fg[0][k % m.encoder_fg[0].length()] = m.encoder_fg[0][k %
            m.encoder_fg[0].length()] +
            0.0001
        }
      }
    }
  }
  for i in 0..<data.length() {
    all_z.push(contrastivevi_encode_fg(m, data[i]))
  }
  m.latent = all_z
  m.is_trained = true
  m
}

///|
pub fn contrastivevi_get_latent(
  model : ContrastiveViModel,
) -> Array[Array[Double]] {
  model.latent
}
