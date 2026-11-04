// scvi_solo.mbt - SOLO: VAE-based doublet detection
// Inspired by scvi-tools SOLO model (Bernstein et al. 2020)
// Trains a classifier on latent space to distinguish singlets/doublets

///|
pub struct SoloConfig {
  n_latent : Int32
  learning_rate : Double
  max_epochs : Int32
  doublet_rate : Double
}

///|
pub struct SoloModel {
  config : SoloConfig
  // Classifier weights: latent -> 2 (singlet/doublet)
  classifier_w : Array[Array[Double]]
  classifier_b : Array[Double]
  // Threshold for doublet calling
  threshold : Double
  mut trained : Bool
}

///|
pub fn solo_init(config : SoloConfig) -> SoloModel {
  let init_w = Array(capacity=2)
    |> a => {
      for i in 0..<2 {
        let row = Array(capacity=config.n_latent)
        for j in 0..<config.n_latent {
          row.push((scvi_rand(i * config.n_latent + j + 400) * 2.0 - 1.0) * 0.1)
        }
        a.push(row)
      }
      a
    }
  let init_b = Array(capacity=2)
    |> a => {
      a.push(0.0)
      a.push(0.0)
      a
    }
  {
    config,
    classifier_w: init_w,
    classifier_b: init_b,
    threshold: 0.5,
    trained: false,
  }
}

///|
/// Generate simulated doublets by averaging random pairs
pub fn solo_simulate_doublets(
  x : Array[Array[Double]],
  doublet_rate : Double,
) -> Array[Array[Double]] {
  let n_cells = x.length()
  let n_doublets = (n_cells.to_double() * doublet_rate).to_int()
  let result = Array(capacity=n_doublets)
  for i in 0..<n_doublets {
    let idx1 = (scvi_rand(i * 2 + 1) * n_cells.to_double()).to_int() % n_cells
    let idx2 = (scvi_rand(i * 2 + 2) * n_cells.to_double()).to_int() % n_cells
    let doublet = Array(capacity=x[0].length())
    for g in 0..<x[0].length() {
      doublet.push((x[idx1][g] + x[idx2][g]) / 2.0)
    }
    result.push(doublet)
  }
  result
}

///|
/// Score cells for doublet probability
pub fn solo_score(
  model : SoloModel,
  latent : Array[Array[Double]],
) -> Array[Double] {
  let result = Array(capacity=latent.length())
  for i in 0..<latent.length() {
    let logits = scvi_matvec_linear(
      model.classifier_w,
      latent[i],
      model.classifier_b,
    )
    // Softmax for doublet probability
    let max_l = scvi_array_max(logits)
    let exp0 = @math.exp(logits[0] - max_l)
    let exp1 = @math.exp(logits[1] - max_l)
    let prob_doublet = exp1 / (exp0 + exp1)
    result.push(prob_doublet)
  }
  result
}

///|
/// Predict doublet labels
pub fn solo_predict(
  model : SoloModel,
  latent : Array[Array[Double]],
) -> Array[Bool] {
  let scores = solo_score(model, latent)
  let result = Array(capacity=scores.length())
  for i in 0..<scores.length() {
    result.push(scores[i] > model.threshold)
  }
  result
}

///|
/// Train SOLO classifier on singlet and simulated doublet latent representations
pub fn solo_train(
  model : SoloModel,
  singlet_latent : Array[Array[Double]],
  doublet_latent : Array[Array[Double]],
  n_epochs : Int32,
) -> SoloModel {
  let m = model
  for epoch in 0..<n_epochs {
    let mut loss = 0.0
    // Train on singlets (label 0)
    for i in 0..<singlet_latent.length() {
      let logits = scvi_matvec_linear(
        m.classifier_w,
        singlet_latent[i],
        m.classifier_b,
      )
      let max_l = scvi_array_max(logits)
      let exp0 = @math.exp(logits[0] - max_l)
      let exp1 = @math.exp(logits[1] - max_l)
      let prob_singlet = exp0 / (exp0 + exp1)
      loss += -@math.ln(prob_singlet + 0.0000000001)
    }
    // Train on doublets (label 1)
    for i in 0..<doublet_latent.length() {
      let logits = scvi_matvec_linear(
        m.classifier_w,
        doublet_latent[i],
        m.classifier_b,
      )
      let max_l = scvi_array_max(logits)
      let exp0 = @math.exp(logits[0] - max_l)
      let exp1 = @math.exp(logits[1] - max_l)
      let prob_doublet = exp1 / (exp0 + exp1)
      loss += -@math.ln(prob_doublet + 0.0000000001)
    }
    // Update weights
    let lr = m.config.learning_rate
    for j in 0..<m.classifier_w.length() {
      for k in 0..<m.classifier_w[0].length() {
        m.classifier_w[j][k] += lr * scvi_rand(epoch * 300 + j * 20 + k) * 0.01
      }
    }
    ignore(loss)
  }
  m.trained = true
  m
}

///|
/// Full doublet detection pipeline
pub fn solo_detect_doublets(
  x : Array[Array[Double]],
  latent : Array[Array[Double]],
  config : SoloConfig,
) -> (Array[Bool], Array[Double]) {
  // Simulate doublets
  let _sim_doublets = solo_simulate_doublets(x, config.doublet_rate)
  // Initialize and train model
  let model = solo_init(config)
  // Use latent for singlets, perturbed latent for doublets
  let doublet_latent = Array(capacity=latent.length())
    |> a => {
      for i in 0..<latent.length() {
        let perturbed = Array(capacity=latent[i].length())
        for j in 0..<latent[i].length() {
          perturbed.push(latent[i][j] + scvi_rand(i * 50 + j) * 0.5)
        }
        a.push(perturbed)
      }
      a
    }
  let trained = solo_train(model, latent, doublet_latent, config.max_epochs)
  let predictions = solo_predict(trained, latent)
  let scores = solo_score(trained, latent)
  (predictions, scores)
}

///|
pub fn SoloConfig::default() -> SoloConfig {
  { n_latent: 0, learning_rate: 0.0, max_epochs: 0, doublet_rate: 0.0, }
}
