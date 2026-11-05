///|
/// METHYLVI: VAE for DNA methylation data (methylcytosine counts).
/// METHYLANVI: Semi-supervised extension for methylation data annotation.
pub struct MethylViConfig {
  n_latent : Int32
  n_hidden : Int32
  n_layers : Int32
  dropout_rate : Double
}

///|
pub fn MethylViConfig::default() -> MethylViConfig {
  { n_latent: 10, n_hidden: 128, n_layers: 2, dropout_rate: 0.1, }
}

///|
pub struct MethylViModel {
  config : MethylViConfig
  n_cpgs : Int32
  encoder : Array[Array[Double]]
  decoder_meth : Array[Array[Double]]
  decoder_unmeth : Array[Array[Double]]
  mut latent : Array[Array[Double]]
  mut is_trained : Bool
}

///|
pub fn methylvi_init(n_cpgs : Int32, config : MethylViConfig) -> MethylViModel {
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
    n_cpgs,
    encoder: make_w(config.n_latent, n_cpgs),
    decoder_meth: make_w(n_cpgs, config.n_latent),
    decoder_unmeth: make_w(n_cpgs, config.n_latent),
    latent: [],
    is_trained: false,
  }
}

///|
pub fn methylvi_encode(
  model : MethylViModel,
  meth : Array[Double],
  unmeth : Array[Double],
) -> Array[Double] {
  let z = Array(capacity=model.config.n_latent)
  for i in 0..<model.config.n_latent {
    let mut sum = 0.0
    for j in 0..<meth.length() {
      if j < model.encoder[i].length() {
        sum = sum + meth[j] * model.encoder[i][j]
      }
    }
    for j in 0..<unmeth.length() {
      if j < model.encoder[i].length() {
        sum = sum + unmeth[j] * model.encoder[i][j] * 0.5
      }
    }
    z.push(sum)
  }
  z
}

///|
pub fn methylvi_decode(
  model : MethylViModel,
  z : Array[Double],
) -> (Array[Double], Array[Double]) {
  let meth_rate = Array(capacity=model.n_cpgs)
  let unmeth_rate = Array(capacity=model.n_cpgs)
  for c in 0..<model.n_cpgs {
    let mut sum_m = 0.0
    let mut sum_u = 0.0
    for j in 0..<z.length() {
      if j < model.decoder_meth[c].length() {
        sum_m = sum_m + z[j] * model.decoder_meth[c][j]
      }
      if j < model.decoder_unmeth[c].length() {
        sum_u = sum_u + z[j] * model.decoder_unmeth[c][j]
      }
    }
    meth_rate.push(@math.exp(sum_m))
    unmeth_rate.push(@math.exp(sum_u))
  }
  (meth_rate, unmeth_rate)
}

///|
pub fn methylvi_get_methylation_fraction(
  model : MethylViModel,
  z : Array[Double],
) -> Array[Double] {
  let (meth_rate, unmeth_rate) = methylvi_decode(model, z)
  let frac = Array(capacity=model.n_cpgs)
  for c in 0..<model.n_cpgs {
    let total = meth_rate[c] + unmeth_rate[c]
    if total > 0.0 {
      frac.push(meth_rate[c] / total)
    } else {
      frac.push(0.5)
    }
  }
  frac
}

///|
pub fn methylvi_train(
  model : MethylViModel,
  meth_data : Array[Array[Double]],
  unmeth_data : Array[Array[Double]],
  n_epochs : Int32,
) -> MethylViModel {
  let m = model
  let all_z = Array(capacity=meth_data.length())
  for _epoch in 0..<n_epochs {
    for i in 0..<meth_data.length() {
      let z = methylvi_encode(m, meth_data[i], unmeth_data[i])
      let (meth_hat, _unmeth_hat) = methylvi_decode(m, z)
      for c in 0..<m.n_cpgs {
        if c < meth_data[i].length() && meth_hat[c] > 0.0 {
          let grad = meth_data[i][c] / meth_hat[c] - 1.0
          for j in 0..<m.config.n_latent {
            if j < m.decoder_meth[c].length() {
              m.decoder_meth[c][j] = m.decoder_meth[c][j] + 0.001 * grad * z[j]
            }
          }
        }
      }
    }
  }
  for i in 0..<meth_data.length() {
    all_z.push(methylvi_encode(m, meth_data[i], unmeth_data[i]))
  }
  m.latent = all_z
  m.is_trained = true
  m
}

// METHYLANVI: Semi-supervised extension

///|
pub struct MethylanviModel {
  base : MethylViModel
  n_labels : Int32
  classifier : Array[Array[Double]]
  label_bias : Array[Double]
}

///|
pub fn methylanvi_init(
  base : MethylViModel,
  n_labels : Int32,
) -> MethylanviModel {
  let cls = Array(capacity=n_labels)
  for _i in 0..<n_labels {
    let row = Array(capacity=base.config.n_latent)
    for _j in 0..<base.config.n_latent {
      row.push(0.01)
    }
    cls.push(row)
  }
  let bias = Array(capacity=n_labels)
  for _i in 0..<n_labels {
    bias.push(0.0)
  }
  { base, n_labels, classifier: cls, label_bias: bias, }
}

///|
pub fn methylanvi_predict(model : MethylanviModel, z : Array[Double]) -> Int32 {
  let mut best = 0
  let mut best_score = -999999999999.0
  for l in 0..<model.n_labels {
    let mut score = model.label_bias[l]
    for j in 0..<z.length() {
      if j < model.classifier[l].length() {
        score = score + z[j] * model.classifier[l][j]
      }
    }
    if score > best_score {
      best_score = score
      best = l
    }
  }
  best
}

///|
pub fn methylanvi_get_latent(model : MethylanviModel) -> Array[Array[Double]] {
  model.base.latent
}
