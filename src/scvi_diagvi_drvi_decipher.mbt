///|
/// DIAGVI: Diagonal-covariance VAE variant.
/// DRVI: Dimensionality reduction VAE with interpretability mixin.
/// Decipher: VAE with custom training plan, rotation utilities, and trajectory analysis.
pub struct DiagViConfig {
  n_latent : Int32
  n_hidden : Int32
  n_layers : Int32
  diag_prior : Bool
}

///|
pub fn DiagViConfig::default() -> DiagViConfig {
  { n_latent: 10, n_hidden: 128, n_layers: 2, diag_prior: true, }
}

///|
pub struct DiagViModel {
  config : DiagViConfig
  n_features : Int32
  encoder_mean : Array[Array[Double]]
  encoder_logvar : Array[Array[Double]]
  decoder : Array[Array[Double]]
  mut latent : Array[Array[Double]]
  mut is_trained : Bool
}

///|
pub fn diagvi_init(n_features : Int32, config : DiagViConfig) -> DiagViModel {
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
    n_features,
    encoder_mean: make_w(config.n_latent, n_features),
    encoder_logvar: make_w(config.n_latent, n_features),
    decoder: make_w(n_features, config.n_latent),
    latent: [],
    is_trained: false,
  }
}

///|
pub fn diagvi_encode(
  model : DiagViModel,
  x : Array[Double],
) -> (Array[Double], Array[Double]) {
  let mean = Array(capacity=model.config.n_latent)
  let logvar = Array(capacity=model.config.n_latent)
  for i in 0..<model.config.n_latent {
    let mut m_sum = 0.0
    let mut v_sum = 0.0
    for j in 0..<x.length() {
      if j < model.encoder_mean[i].length() {
        m_sum = m_sum + x[j] * model.encoder_mean[i][j]
        v_sum = v_sum + x[j] * model.encoder_logvar[i][j]
      }
    }
    mean.push(m_sum)
    logvar.push(v_sum)
  }
  (mean, logvar)
}

///|
pub fn diagvi_diag_kl(mean : Array[Double], logvar : Array[Double]) -> Double {
  let mut kl = 0.0
  for i in 0..<mean.length() {
    if i < logvar.length() {
      kl = kl +
        0.5 * (@math.exp(logvar[i]) + mean[i] * mean[i] - 1.0 - logvar[i])
    }
  }
  kl
}

///|
pub fn diagvi_train(
  model : DiagViModel,
  data : Array[Array[Double]],
  n_epochs : Int32,
) -> DiagViModel {
  let m = model
  let all_z = Array(capacity=data.length())
  for _epoch in 0..<n_epochs {
    for i in 0..<data.length() {
      let (mean, _logvar) = diagvi_encode(m, data[i])
      let z = Array(capacity=mean.length())
      for j in 0..<mean.length() {
        z.push(mean[j])
      }
      let x_hat = Array(capacity=m.n_features)
      for g in 0..<m.n_features {
        let mut sum = 0.0
        for j in 0..<z.length() {
          if j < m.decoder[g].length() {
            sum = sum + z[j] * m.decoder[g][j]
          }
        }
        x_hat.push(sum)
      }
      for g in 0..<m.n_features {
        if g < data[i].length() {
          let grad = data[i][g] - x_hat[g]
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
    let (mean, _) = diagvi_encode(m, data[i])
    all_z.push(mean)
  }
  m.latent = all_z
  m.is_trained = true
  m
}

// DRVI: Dimensionality reduction VAE with interpretability

///|
pub struct DrviConfig {
  n_latent : Int32
  n_hidden : Int32
  sparsity_lambda : Double
  interpretability_weight : Double
}

///|
pub fn DrviConfig::default() -> DrviConfig {
  {
    n_latent: 10,
    n_hidden: 128,
    sparsity_lambda: 0.01,
    interpretability_weight: 0.1,
  }
}

///|
pub struct DrviModel {
  config : DrviConfig
  n_features : Int32
  encoder : Array[Array[Double]]
  decoder : Array[Array[Double]]
  mut feature_importance : Array[Double]
  mut latent : Array[Array[Double]]
  mut is_trained : Bool
}

///|
pub fn drvi_init(n_features : Int32, config : DrviConfig) -> DrviModel {
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
  let imp = Array(capacity=n_features)
  for _i in 0..<n_features {
    imp.push(0.0)
  }
  {
    config,
    n_features,
    encoder: make_w(config.n_latent, n_features),
    decoder: make_w(n_features, config.n_latent),
    feature_importance: imp,
    latent: [],
    is_trained: false,
  }
}

///|
pub fn drvi_compute_importance(model : DrviModel) -> Array[Double] {
  let importance = Array(capacity=model.n_features)
  for g in 0..<model.n_features {
    let mut score = 0.0
    for j in 0..<model.config.n_latent {
      if j < model.decoder[g].length() {
        score = score + model.decoder[g][j] * model.decoder[g][j]
      }
    }
    importance.push(score)
  }
  importance
}

///|
pub fn drvi_train(
  model : DrviModel,
  data : Array[Array[Double]],
  n_epochs : Int32,
) -> DrviModel {
  let m = model
  let all_z = Array(capacity=data.length())
  for _epoch in 0..<n_epochs {
    for i in 0..<data.length() {
      let z = Array(capacity=m.config.n_latent)
      for k in 0..<m.config.n_latent {
        let mut sum = 0.0
        for j in 0..<data[i].length() {
          if j < m.encoder[k].length() {
            sum = sum + data[i][j] * m.encoder[k][j]
          }
        }
        z.push(sum)
      }
      let x_hat = Array(capacity=m.n_features)
      for g in 0..<m.n_features {
        let mut sum = 0.0
        for j in 0..<z.length() {
          if j < m.decoder[g].length() {
            sum = sum + z[j] * m.decoder[g][j]
          }
        }
        x_hat.push(sum)
      }
      for g in 0..<m.n_features {
        if g < data[i].length() {
          let grad = data[i][g] - x_hat[g]
          for j in 0..<m.config.n_latent {
            if j < m.decoder[g].length() {
              m.decoder[g][j] = m.decoder[g][j] +
                0.001 * grad * z[j] -
                m.config.sparsity_lambda * m.decoder[g][j]
            }
          }
        }
      }
    }
  }
  m.feature_importance = drvi_compute_importance(m)
  for i in 0..<data.length() {
    let z = Array(capacity=m.config.n_latent)
    for k in 0..<m.config.n_latent {
      let mut sum = 0.0
      for j in 0..<data[i].length() {
        if j < m.encoder[k].length() {
          sum = sum + data[i][j] * m.encoder[k][j]
        }
      }
      z.push(sum)
    }
    all_z.push(z)
  }
  m.latent = all_z
  m.is_trained = true
  m
}

// Decipher: VAE with trajectory analysis

///|
pub struct DecipherConfig {
  n_latent : Int32
  n_hidden : Int32
  rotation_angle : Double
}

///|
pub fn DecipherConfig::default() -> DecipherConfig {
  { n_latent: 10, n_hidden: 128, rotation_angle: 0.0, }
}

///|
pub struct DecipherModel {
  config : DecipherConfig
  n_features : Int32
  encoder : Array[Array[Double]]
  decoder : Array[Array[Double]]
  mut latent : Array[Array[Double]]
  mut pseudotime : Array[Double]
  mut is_trained : Bool
}

///|
pub fn decipher_init(
  n_features : Int32,
  config : DecipherConfig,
) -> DecipherModel {
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
    n_features,
    encoder: make_w(config.n_latent, n_features),
    decoder: make_w(n_features, config.n_latent),
    latent: [],
    pseudotime: [],
    is_trained: false,
  }
}

///|
pub fn decipher_rotate_latent(
  model : DecipherModel,
  z : Array[Double],
) -> Array[Double] {
  let angle = model.config.rotation_angle
  let cos_a = @math.cos(angle)
  let sin_a = @math.sin(angle)
  let rotated = Array(capacity=z.length())
  let n = z.length()
  for i in 0..<n {
    let j = (i + 1) % n
    rotated.push(z[i] * cos_a - z[j] * sin_a)
  }
  rotated
}

///|
pub fn decipher_compute_pseudotime(model : DecipherModel) -> Array[Double] {
  let n = model.latent.length()
  let pt = Array(capacity=n)
  for i in 0..<n {
    let mut sum = 0.0
    for j in 0..<model.latent[i].length() {
      sum = sum + model.latent[i][j] * model.latent[i][j]
    }
    pt.push(@math.pow(sum, 0.5))
  }
  pt
}

///|
pub fn decipher_train(
  model : DecipherModel,
  data : Array[Array[Double]],
  n_epochs : Int32,
) -> DecipherModel {
  let m = model
  let all_z = Array(capacity=data.length())
  for _epoch in 0..<n_epochs {
    for i in 0..<data.length() {
      let z = Array(capacity=m.config.n_latent)
      for k in 0..<m.config.n_latent {
        let mut sum = 0.0
        for j in 0..<data[i].length() {
          if j < m.encoder[k].length() {
            sum = sum + data[i][j] * m.encoder[k][j]
          }
        }
        z.push(sum)
      }
      let x_hat = Array(capacity=m.n_features)
      for g in 0..<m.n_features {
        let mut sum = 0.0
        for j in 0..<z.length() {
          if j < m.decoder[g].length() {
            sum = sum + z[j] * m.decoder[g][j]
          }
        }
        x_hat.push(sum)
      }
      for g in 0..<m.n_features {
        if g < data[i].length() {
          let grad = data[i][g] - x_hat[g]
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
    let z = Array(capacity=m.config.n_latent)
    for k in 0..<m.config.n_latent {
      let mut sum = 0.0
      for j in 0..<data[i].length() {
        if j < m.encoder[k].length() {
          sum = sum + data[i][j] * m.encoder[k][j]
        }
      }
      z.push(sum)
    }
    all_z.push(z)
  }
  m.latent = all_z
  m.pseudotime = decipher_compute_pseudotime(m)
  m.is_trained = true
  m
}

///|
pub fn _scvi_diagvi_drvi_decipher_dummy() -> Unit {
  let _ = DiagViConfig::{
    n_latent: 1,
    n_hidden: 1,
    n_layers: 1,
    diag_prior: false,
  }
  let _ = DrviConfig::{
    n_latent: 1,
    n_hidden: 1,
    sparsity_lambda: 0.0,
    interpretability_weight: 0.0,
  }
  let _ = DecipherConfig::{ n_latent: 1, n_hidden: 1, rotation_angle: 0.0, }
}
