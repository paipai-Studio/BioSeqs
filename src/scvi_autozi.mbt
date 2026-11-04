// scvi_autozi.mbt - AUTOZI: automatic zero-inflation detection
// Inspired by scvi-tools AUTOZI model
// Learns per-gene zero-inflation probabilities for gene selection

///|
pub struct AutoZIConfig {
  learning_rate : Double
  max_epochs : Int32
  alpha_prior : Double
  beta_prior : Double
}

///|
pub struct AutoZIModel {
  config : AutoZIConfig
  // Per-gene Bernoulli dropout parameters (logit)
  gene_logits : Array[Double]
  // Decoder
  decoder_w : Array[Array[Double]]
  decoder_b : Array[Double]
  // Encoder
  encoder_w : Array[Array[Double]]
  encoder_b : Array[Double]
  n_genes : Int32
  mut trained : Bool
}

///|
pub fn autozi_init(n_genes : Int32, config : AutoZIConfig) -> AutoZIModel {
  let init_w = fn(rows : Int32, cols : Int32) -> Array[Array[Double]] {
    let scale = 6.0 / (rows.to_double().sqrt() + cols.to_double())
    let r = Array(capacity=rows)
    for i in 0..<rows {
      let row = Array(capacity=cols)
      for j in 0..<cols {
        row.push((scvi_rand(i * cols + j + 700) * 2.0 - 1.0) * scale)
      }
      r.push(row)
    }
    r
  }
  let init_b = fn(n : Int32) -> Array[Double] {
    let r = Array(capacity=n)
    for _i in 0..<n {
      r.push(0.0)
    }
    r
  }
  let gene_logits = Array(capacity=n_genes)
    |> a => {
      for _i in 0..<n_genes {
        a.push(0.0) // Start with 50% dropout
      }
      a
    }
  {
    config,
    gene_logits,
    decoder_w: init_w(config.max_epochs.min(10), n_genes),
    decoder_b: init_b(n_genes),
    encoder_w: init_w(n_genes, config.max_epochs.min(10)),
    encoder_b: init_b(config.max_epochs.min(10)),
    n_genes,
    trained: false,
  }
}

///|
/// Get alpha/beta parameters (zero-inflation probabilities)
pub fn autozi_get_alphas_betas(
  model : AutoZIModel,
) -> (Array[Double], Array[Double]) {
  let n_genes = model.n_genes
  let alphas = Array(capacity=n_genes)
  let betas = Array(capacity=n_genes)
  for g in 0..<n_genes {
    // Sigmoid to get dropout probability
    let logit = model.gene_logits[g]
    let prob = 1.0 / (1.0 + @math.exp(-logit))
    // Alpha = probability of gene being "active" (1 - dropout)
    alphas.push(1.0 - prob)
    // Beta = dropout probability
    betas.push(prob)
  }
  (alphas, betas)
}

///|
/// Get per-gene dropout probabilities
pub fn autozi_get_dropout_probs(model : AutoZIModel) -> Array[Double] {
  let result = Array(capacity=model.n_genes)
  for g in 0..<model.n_genes {
    let logit = model.gene_logits[g]
    result.push(1.0 / (1.0 + @math.exp(-logit)))
  }
  result
}

///|
/// Select genes based on zero-inflation (low dropout = informative)
pub fn autozi_select_genes(
  model : AutoZIModel,
  n_top_genes : Int32,
) -> Array[Int32] {
  let dropout_probs = autozi_get_dropout_probs(model)
  let indexed = Array(capacity=model.n_genes)
    |> a => {
      for g in 0..<model.n_genes {
        a.push((g, dropout_probs[g]))
      }
      a
    }
  // Sort by dropout probability (ascending = least dropout first)
  indexed.sort_by(fn(a, b) {
    if a.1 < b.1 {
      -1
    } else if a.1 > b.1 {
      1
    } else {
      0
    }
  })
  let result = Array(capacity=n_top_genes)
  for i in 0..<n_top_genes {
    if i < indexed.length() {
      result.push(indexed[i].0)
    }
  }
  result
}

///|
/// Train AUTOZI model
pub fn autozi_train(
  model : AutoZIModel,
  x : Array[Array[Double]],
  n_epochs : Int32,
) -> AutoZIModel {
  let m = model
  let n_cells = x.length()
  for epoch in 0..<n_epochs {
    // Update gene logits based on zero fraction
    for g in 0..<m.n_genes {
      let mut n_zero = 0
      let mut _n_nonzero = 0
      for i in 0..<n_cells {
        if x[i][g] < 0.5 {
          n_zero += 1
        } else {
          _n_nonzero += 1
        }
      }
      // Update logit: push towards 0 if many zeros, towards 1 if few zeros
      let zero_frac = n_zero.to_double() / n_cells.to_double()
      let target_logit = @math.ln(zero_frac + 0.0000000001) -
        @math.ln(1.0 - zero_frac + 0.0000000001)
      let lr = m.config.learning_rate
      m.gene_logits[g] += lr * (target_logit - m.gene_logits[g])
      // Regularization towards prior
      let prior_grad = (m.config.alpha_prior - m.config.beta_prior) * 0.01
      m.gene_logits[g] += lr * prior_grad
    }
    ignore(epoch)
  }
  m.trained = true
  m
}

///|
/// Classify genes as dropout-inflated vs truly zero
pub fn autozi_classify_genes(
  model : AutoZIModel,
  threshold : Double,
) -> Array[Bool] {
  let dropout_probs = autozi_get_dropout_probs(model)
  let result = Array(capacity=model.n_genes)
  for g in 0..<model.n_genes {
    // True = gene is significantly zero-inflated
    result.push(dropout_probs[g] > threshold)
  }
  result
}

///|
pub fn AutoZIConfig::default() -> AutoZIConfig {
  { learning_rate: 0.0, max_epochs: 0, alpha_prior: 0.0, beta_prior: 0.0, }
}
