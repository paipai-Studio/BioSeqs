///|
/// TOTALANVI: Semi-supervised extension of totalVI for CITE-seq cell annotation.
/// Combines RNA + protein data with a classifier head for cell type prediction.
pub struct TotalAnviConfig {
  n_latent : Int32
  n_layers : Int32
  n_hidden : Int32
  dropout_rate : Double
  use_batch_norm : Bool
  use_layer_norm : Bool
  n_labels : Int32
  classifier_parameters : Array[Double]
}

///|
pub fn TotalAnviConfig::default() -> TotalAnviConfig {
  {
    n_latent: 10,
    n_layers: 2,
    n_hidden: 128,
    dropout_rate: 0.1,
    use_batch_norm: true,
    use_layer_norm: false,
    n_labels: 0,
    classifier_parameters: [],
  }
}

///|
pub struct TotalAnviModel {
  config : TotalAnviConfig
  n_genes : Int32
  n_proteins : Int32
  encoder_weights : Array[Array[Double]]
  decoder_rna_weights : Array[Array[Double]]
  decoder_protein_weights : Array[Array[Double]]
  classifier_weights : Array[Array[Double]]
  label_logits : Array[Double]
  mut is_trained : Bool
}

///|
pub fn totalanvi_init(
  n_genes : Int32,
  n_proteins : Int32,
  config : TotalAnviConfig,
) -> TotalAnviModel {
  let input_dim = n_genes + n_proteins
  let enc = Array(capacity=config.n_layers)
  let prev = input_dim
  for _i in 0..<config.n_layers {
    let row = Array(capacity=prev)
    for _j in 0..<prev {
      row.push(0.01)
    }
    enc.push(row)
  }
  let dec_rna = Array(capacity=config.n_latent)
  for _i in 0..<config.n_latent {
    let row = Array(capacity=n_genes)
    for _j in 0..<n_genes {
      row.push(0.01)
    }
    dec_rna.push(row)
  }
  let dec_prot = Array(capacity=config.n_latent)
  for _i in 0..<config.n_latent {
    let row = Array(capacity=n_proteins)
    for _j in 0..<n_proteins {
      row.push(0.01)
    }
    dec_prot.push(row)
  }
  let cls = Array(capacity=config.n_labels)
  for _i in 0..<config.n_labels {
    let row = Array(capacity=config.n_latent)
    for _j in 0..<config.n_latent {
      row.push(0.01)
    }
    cls.push(row)
  }
  let logits = Array(capacity=config.n_labels)
  for _i in 0..<config.n_labels {
    logits.push(0.0)
  }
  {
    config,
    n_genes,
    n_proteins,
    encoder_weights: enc,
    decoder_rna_weights: dec_rna,
    decoder_protein_weights: dec_prot,
    classifier_weights: cls,
    label_logits: logits,
    is_trained: false,
  }
}

///|
pub fn totalanvi_encode(
  model : TotalAnviModel,
  rna : Array[Double],
  protein : Array[Double],
) -> Array[Double] {
  let z = Array(capacity=model.config.n_latent)
  let combined = Array(capacity=rna.length() + protein.length())
  for i in 0..<rna.length() {
    combined.push(rna[i])
  }
  for i in 0..<protein.length() {
    combined.push(protein[i])
  }
  for _i in 0..<model.config.n_latent {
    let mut sum = 0.0
    for j in 0..<combined.length() {
      if j < model.encoder_weights[0].length() {
        sum = sum + combined[j] * model.encoder_weights[0][j]
      }
    }
    z.push(sum)
  }
  z
}

///|
pub fn totalanvi_predict(model : TotalAnviModel, z : Array[Double]) -> Int32 {
  let mut best = 0
  let mut best_score = -999999999999.0
  for l in 0..<model.config.n_labels {
    let mut score = model.label_logits[l]
    for j in 0..<z.length() {
      if j < model.classifier_weights[l].length() {
        score = score + z[j] * model.classifier_weights[l][j]
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
pub fn totalanvi_predict_proba(
  model : TotalAnviModel,
  z : Array[Double],
) -> Array[Double] {
  let scores = Array(capacity=model.config.n_labels)
  for l in 0..<model.config.n_labels {
    let mut score = model.label_logits[l]
    for j in 0..<z.length() {
      if j < model.classifier_weights[l].length() {
        score = score + z[j] * model.classifier_weights[l][j]
      }
    }
    scores.push(score)
  }
  let mut max_s = scores[0]
  for i in 1..<scores.length() {
    if scores[i] > max_s {
      max_s = scores[i]
    }
  }
  let mut sum_exp = 0.0
  for i in 0..<scores.length() {
    scores[i] = @math.exp(scores[i] - max_s)
    sum_exp = sum_exp + scores[i]
  }
  for i in 0..<scores.length() {
    scores[i] = scores[i] / sum_exp
  }
  scores
}

///|
pub fn totalanvi_train(
  model : TotalAnviModel,
  rna_data : Array[Array[Double]],
  protein_data : Array[Array[Double]],
  labels : Array[Int32],
  n_epochs : Int32,
) -> TotalAnviModel {
  let m = model
  for _epoch in 0..<n_epochs {
    for i in 0..<rna_data.length() {
      let _z = totalanvi_encode(m, rna_data[i], protein_data[i])
      let _pred = totalanvi_predict(m, _z)
      if _pred != labels[i] {
        for j in 0..<m.config.n_latent {
          if j < m.classifier_weights[labels[i]].length() {
            m.classifier_weights[labels[i]][j] = m.classifier_weights[labels[i]][j] +
              0.01 * _z[j]
          }
        }
      }
    }
  }
  m.is_trained = true
  m
}
