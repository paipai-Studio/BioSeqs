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
    encoder_w: init_w(config.n_hidden, n_genes),
    encoder_b: init_b(config.n_hidden),
    z_encoder_w: init_w(config.n_latent, config.n_hidden),
    z_encoder_b: init_b(config.n_latent),
    classifier_w: init_w(n_labels, config.n_latent),
    classifier_b: init_b(n_labels),
    decoder_w: init_w(n_genes, config.n_latent),
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
  let n_labeled = x_labeled.length()
  let n_unlabeled = x_unlabeled.length()
  let n_total = n_labeled + n_unlabeled
  if n_total == 0 {
    m.trained = true
    return m
  }
  let theta = Array(capacity=m.n_genes)
  for _g in 0..<m.n_genes {
    theta.push(1.0)
  }
  let template = ScviAgTape::new()
  ignore(template.matrix_param(m.encoder_w))
  ignore(template.vector_param(m.encoder_b))
  ignore(template.matrix_param(m.z_encoder_w))
  ignore(template.vector_param(m.z_encoder_b))
  ignore(template.matrix_param(m.classifier_w))
  ignore(template.vector_param(m.classifier_b))
  ignore(template.matrix_param(m.decoder_w))
  ignore(template.vector_param(m.decoder_b))
  let adam = scvi_adam_init(template, m.config.learning_rate)
  let n_total_d = n_total.to_double()
  for _epoch in 0..<n_epochs {
    let t = ScviAgTape::new()
    let w_enc = t.matrix_param(m.encoder_w)
    let b_enc = t.vector_param(m.encoder_b)
    let w_z = t.matrix_param(m.z_encoder_w)
    let b_z = t.vector_param(m.z_encoder_b)
    let w_cls = t.matrix_param(m.classifier_w)
    let b_cls = t.vector_param(m.classifier_b)
    let w_dec = t.matrix_param(m.decoder_w)
    let b_dec = t.vector_param(m.decoder_b)
    let mut total = t.leaf([0.0])
    for i in 0..<n_labeled {
      let xi = x_labeled[i]
      let lib = scvi_array_sum(xi)
      let lib_size = if lib <= 0.0 { 1.0 } else { lib }
      let xin = t.leaf(xi)
      let hidden = t.relu(t.dense(xin, w_enc, b_enc))
      let z = t.dense(hidden, w_z, b_z)
      let probs = t.softmax(t.dense(z, w_cls, b_cls))
      let ce = t.neg(t.log(t.select(probs, y_labeled[i])))
      let rates = t.softmax(t.dense(z, w_dec, b_dec))
      let nll = t.nb_nll(t.scale(rates, lib_size), xi, theta)
      total = t.add(total, t.add(ce, nll))
    }
    for i in 0..<n_unlabeled {
      let xi = x_unlabeled[i]
      let lib = scvi_array_sum(xi)
      let lib_size = if lib <= 0.0 { 1.0 } else { lib }
      let xin = t.leaf(xi)
      let hidden = t.relu(t.dense(xin, w_enc, b_enc))
      let z = t.dense(hidden, w_z, b_z)
      let probs = t.softmax(t.dense(z, w_cls, b_cls))
      let ent = t.neg(t.sum(t.mul(probs, t.log(probs))))
      let rates = t.softmax(t.dense(z, w_dec, b_dec))
      let nll = t.nb_nll(t.scale(rates, lib_size), xi, theta)
      total = t.add(total, t.add(ent, nll))
    }
    let loss = t.scale(total, 1.0 / n_total_d)
    t.backward(loss)
    ignore(t.clip_grad_norm(5.0))
    scvi_adam_step(t, adam)
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
pub fn SCANVIConfig::new(
  n_hidden : Int32,
  n_latent : Int32,
  n_labels : Int32,
  learning_rate : Double,
  max_epochs : Int32,
) -> SCANVIConfig {
  { n_hidden, n_latent, n_labels, learning_rate, max_epochs, }
}

///|
pub fn SCANVIConfig::default() -> SCANVIConfig {
  { n_hidden: 0, n_latent: 0, n_labels: 0, learning_rate: 0.0, max_epochs: 0, }
}
