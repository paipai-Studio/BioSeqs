///|
/// SCVIVA: VAE for spatial niche differential expression analysis.
pub struct ScvivaConfig {
  n_latent : Int32
  n_hidden : Int32
  n_niches : Int32
  n_layers : Int32
}

///|
pub fn ScvivaConfig::default() -> ScvivaConfig {
  { n_latent: 10, n_hidden: 128, n_niches: 5, n_layers: 2, }
}

///|
pub struct ScvivaModel {
  config : ScvivaConfig
  n_genes : Int32
  encoder : Array[Array[Double]]
  decoder : Array[Array[Double]]
  niche_classifier : Array[Array[Double]]
  niche_bias : Array[Double]
  mut latent : Array[Array[Double]]
  mut niche_labels : Array[Int32]
  mut is_trained : Bool
}

///|
pub fn scviva_init(n_genes : Int32, config : ScvivaConfig) -> ScvivaModel {
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
  let bias = Array(capacity=config.n_niches)
  for _i in 0..<config.n_niches {
    bias.push(0.0)
  }
  {
    config,
    n_genes,
    encoder: make_w(config.n_latent, n_genes),
    decoder: make_w(n_genes, config.n_latent),
    niche_classifier: make_w(config.n_niches, config.n_latent),
    niche_bias: bias,
    latent: [],
    niche_labels: [],
    is_trained: false,
  }
}

///|
pub fn scviva_encode(model : ScvivaModel, x : Array[Double]) -> Array[Double] {
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
pub fn scviva_assign_niche(model : ScvivaModel, z : Array[Double]) -> Int32 {
  let mut best = 0
  let mut best_score = -999999999999.0
  for n in 0..<model.config.n_niches {
    let mut score = model.niche_bias[n]
    for j in 0..<z.length() {
      if j < model.niche_classifier[n].length() {
        score = score + z[j] * model.niche_classifier[n][j]
      }
    }
    if score > best_score {
      best_score = score
      best = n
    }
  }
  best
}

///|
pub fn scviva_niche_de(
  model : ScvivaModel,
  data : Array[Array[Double]],
  niche1 : Int32,
  niche2 : Int32,
) -> Array[(Int32, Double)] {
  let n_genes = model.n_genes
  let mean1 = Array(capacity=n_genes)
  let mean2 = Array(capacity=n_genes)
  let mut count1 = 0
  let mut count2 = 0
  for _g in 0..<n_genes {
    mean1.push(0.0)
    mean2.push(0.0)
  }
  for i in 0..<data.length() {
    if i < model.niche_labels.length() {
      if model.niche_labels[i] == niche1 {
        for g in 0..<n_genes {
          if g < data[i].length() {
            mean1[g] = mean1[g] + data[i][g]
          }
        }
        count1 = count1 + 1
      } else if model.niche_labels[i] == niche2 {
        for g in 0..<n_genes {
          if g < data[i].length() {
            mean2[g] = mean2[g] + data[i][g]
          }
        }
        count2 = count2 + 1
      }
    }
  }
  let results = Array(capacity=n_genes)
  for g in 0..<n_genes {
    let m1 = if count1 > 0 { mean1[g] / count1.to_double() } else { 0.0 }
    let m2 = if count2 > 0 { mean2[g] / count2.to_double() } else { 0.0 }
    let lfc = if m2 > 0.0 { @math.ln(m1 / m2) } else { 0.0 }
    results.push((g, lfc))
  }
  results.sort_by(fn(a, b) {
    let abs_a = if a.1 < 0.0 { -a.1 } else { a.1 }
    let abs_b = if b.1 < 0.0 { -b.1 } else { b.1 }
    if abs_a > abs_b {
      -1
    } else if abs_a < abs_b {
      1
    } else {
      0
    }
  })
  results
}

///|
pub fn scviva_marker_classifier(
  model : ScvivaModel,
  gene_expr : Array[Double],
  niche : Int32,
) -> Double {
  let mut score = 0.0
  for g in 0..<gene_expr.length() {
    if niche < model.config.n_niches &&
      g < model.niche_classifier[niche].length() {
      score = score + gene_expr[g] * model.niche_classifier[niche][g]
    }
  }
  score
}

///|
pub fn scviva_train(
  model : ScvivaModel,
  data : Array[Array[Double]],
  n_epochs : Int32,
) -> ScvivaModel {
  let m = model
  let all_z = Array(capacity=data.length())
  let all_niches = Array(capacity=data.length())
  for _epoch in 0..<n_epochs {
    for i in 0..<data.length() {
      let z = scviva_encode(m, data[i])
      let niche = scviva_assign_niche(m, z)
      let x_hat = Array(capacity=m.n_genes)
      for g in 0..<m.n_genes {
        let mut sum = 0.0
        for j in 0..<z.length() {
          if j < m.decoder[g].length() {
            sum = sum + z[j] * m.decoder[g][j]
          }
        }
        x_hat.push(@math.exp(sum))
      }
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
      ignore(niche)
    }
  }
  for i in 0..<data.length() {
    let z = scviva_encode(m, data[i])
    all_z.push(z)
    all_niches.push(scviva_assign_niche(m, z))
  }
  m.latent = all_z
  m.niche_labels = all_niches
  m.is_trained = true
  m
}
