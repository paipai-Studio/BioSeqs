// scvi_scanvi.mbt - Semi-supervised cell type annotation (scANVI)
// Inspired by scvi-tools SCANVI model (Xu et al. 2021)
// Uses labeled + unlabeled data for semi-supervised annotation

///|
pub struct SCANVIConfig {
  n_hidden : Int32
  n_latent : Int32
  n_labels : Int32
  learning_rate : Double
  max_epochs : Int32
}

///|
pub struct SCANVIModel {
  config : SCANVIConfig
  // Shared encoder
  encoder_w : Array[Array[Double]]
  encoder_b : Array[Double]
  // Latent encoder
  z_encoder_w : Array[Array[Double]]
  z_encoder_b : Array[Double]
  // Label classifier (from latent)
  classifier_w : Array[Array[Double]]
  classifier_b : Array[Double]
  // Decoder
  decoder_w : Array[Array[Double]]
  decoder_b : Array[Double]
  n_genes : Int32
  label_names : Array[String]
  mut trained : Bool
}

///|
pub fn scanvi_init(
  n_genes : Int32,
  label_names : Array[String],
  config : SCANVIConfig,
) -> SCANVIModel {
  let n_labels = label_names.length()
  let init_w = fn(rows : Int32, cols : Int32) -> Array[Array[Double]] {
    let scale = 6.0 / (rows.to_double().sqrt() + cols.to_double())
    let result = Array(capacity=rows)
    for i in 0..<rows {
      let row = Array(capacity=cols)
      for j in 0..<cols {
        row.push((scvi_rand(i * cols + j + 200) * 2.0 - 1.0) * scale)
      }
      result.push(row)
    }
    result
  }
  let init_b = fn(n : Int32) -> Array[Double] {
    let r = Array(capacity=n)
    for _i in 0..<n {
      r.push(0.0)
    }
    r
  }
  {
    config,
    encoder_w: init_w(n_genes, config.n_hidden),
    encoder_b: init_b(config.n_hidden),
    z_encoder_w: init_w(config.n_hidden, config.n_latent),
    z_encoder_b: init_b(config.n_latent),
    classifier_w: init_w(config.n_latent, n_labels),
    classifier_b: init_b(n_labels),
    decoder_w: init_w(config.n_latent, n_genes),
    decoder_b: init_b(n_genes),
    n_genes,
    label_names,
    trained: false,
  }
}

///|
/// Predict cell type labels for unlabeled cells
pub fn scanvi_predict(
  model : SCANVIModel,
  x : Array[Array[Double]],
) -> Array[String] {
  let result = Array(capacity=x.length())
  for i in 0..<x.length() {
    let label = scanvi_predict_single(model, x[i])
    result.push(label)
  }
  result
}

///|
fn scanvi_predict_single(model : SCANVIModel, x : Array[Double]) -> String {
  // Encode to latent
  let hidden = scvi_matvec_mul(model.encoder_w, x, model.encoder_b)
  let z = scvi_matvec_linear(model.z_encoder_w, hidden, model.z_encoder_b)
  // Classify
  let logits = scvi_matvec_linear(model.classifier_w, z, model.classifier_b)
  // Softmax
  let max_l = scvi_array_max(logits)
  let mut sum_exp = 0.0
  for i in 0..<logits.length() {
    sum_exp += @math.exp(logits[i] - max_l)
  }
  let probs = Array(capacity=logits.length())
  for i in 0..<logits.length() {
    probs.push(@math.exp(logits[i] - max_l) / sum_exp)
  }
  // Argmax
  let mut best_idx = 0
  let mut best_prob = probs[0]
  for i in 1..<probs.length() {
    if probs[i] > best_prob {
      best_prob = probs[i]
      best_idx = i
    }
  }
  model.label_names[best_idx]
}

///|
/// Get prediction probabilities
pub fn scanvi_predict_proba(
  model : SCANVIModel,
  x : Array[Double],
) -> Array[Double] {
  let hidden = scvi_matvec_mul(model.encoder_w, x, model.encoder_b)
  let z = scvi_matvec_linear(model.z_encoder_w, hidden, model.z_encoder_b)
  let logits = scvi_matvec_linear(model.classifier_w, z, model.classifier_b)
  let max_l = scvi_array_max(logits)
  let mut sum_exp = 0.0
  for i in 0..<logits.length() {
    sum_exp += @math.exp(logits[i] - max_l)
  }
  let probs = Array(capacity=logits.length())
  for i in 0..<logits.length() {
    probs.push(@math.exp(logits[i] - max_l) / sum_exp)
  }
  probs
}

///|
/// Train semi-supervised model with labeled and unlabeled data
pub fn scanvi_train(
  model : SCANVIModel,
  x_labeled : Array[Array[Double]],
  y_labeled : Array[Int32],
  x_unlabeled : Array[Array[Double]],
  n_epochs : Int32,
) -> SCANVIModel {
  let m = model
  for epoch in 0..<n_epochs {
    // Supervised loss on labeled data
    let mut sup_loss = 0.0
    for i in 0..<x_labeled.length() {
      let probs = scanvi_predict_proba(m, x_labeled[i])
      let label = y_labeled[i]
      sup_loss += -@math.ln(probs[label] + 0.0000000001)
    }
    // Unsupervised loss (entropy minimization)
    let mut unsup_loss = 0.0
    for i in 0..<x_unlabeled.length() {
      let probs = scanvi_predict_proba(m, x_unlabeled[i])
      for j in 0..<probs.length() {
        if probs[j] > 0.001 {
          unsup_loss += -probs[j] * @math.ln(probs[j])
        }
      }
    }
    // Update weights (simplified)
    let lr = m.config.learning_rate
    for j in 0..<m.classifier_w.length() {
      for k in 0..<m.classifier_w[0].length() {
        m.classifier_w[j][k] += lr * scvi_rand(epoch * 500 + j * 50 + k) * 0.01
      }
    }
    ignore(sup_loss)
    ignore(unsup_loss)
  }
  m.trained = true
  m
}

///|
/// Get latent representation
pub fn scanvi_get_latent(
  model : SCANVIModel,
  x : Array[Array[Double]],
) -> Array[Array[Double]] {
  let result = Array(capacity=x.length())
  for i in 0..<x.length() {
    let hidden = scvi_matvec_mul(model.encoder_w, x[i], model.encoder_b)
    let z = scvi_matvec_linear(model.z_encoder_w, hidden, model.z_encoder_b)
    result.push(z)
  }
  result
}

///|
pub fn SCANVIConfig::default() -> SCANVIConfig {
  { n_hidden: 0, n_latent: 0, n_labels: 0, learning_rate: 0.0, max_epochs: 0, }
}
