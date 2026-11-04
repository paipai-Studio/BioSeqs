// scvi_lda.mbt - Amortized LDA: topic model for single-cell
// Inspired by scvi-tools AmortizedLDA model
// Discovers cell programs / topics over genes

///|
pub struct LDAConfig {
  n_topics : Int32
  learning_rate : Double
  max_epochs : Int32
}

///|
pub struct LDAModel {
  config : LDAConfig
  // Topic-gene distribution (beta): topic x gene
  topic_gene : Array[Array[Double]]
  // Cell-topic distribution (theta): cell x topic
  cell_topic : Array[Array[Double]]
  // Encoder for amortized inference
  encoder_w : Array[Array[Double]]
  encoder_b : Array[Double]
  n_genes : Int32
  n_cells : Int32
  mut trained : Bool
}

///|
pub fn lda_init(
  n_genes : Int32,
  n_cells : Int32,
  config : LDAConfig,
) -> LDAModel {
  let init_w = fn(rows : Int32, cols : Int32) -> Array[Array[Double]] {
    let scale = 6.0 / (rows.to_double().sqrt() + cols.to_double())
    let r = Array(capacity=rows)
    for i in 0..<rows {
      let row = Array(capacity=cols)
      for j in 0..<cols {
        row.push((scvi_rand(i * cols + j + 1200) * 2.0 - 1.0) * scale)
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
  // Initialize topic-gene uniformly
  let topic_gene = Array(capacity=config.n_topics)
    |> a => {
      for _t in 0..<config.n_topics {
        let row = Array(capacity=n_genes)
        for _g in 0..<n_genes {
          row.push(1.0 / n_genes.to_double())
        }
        a.push(row)
      }
      a
    }
  // Initialize cell-topic uniformly
  let cell_topic = Array(capacity=n_cells)
    |> a => {
      for _c in 0..<n_cells {
        let row = Array(capacity=config.n_topics)
        for _t in 0..<config.n_topics {
          row.push(1.0 / config.n_topics.to_double())
        }
        a.push(row)
      }
      a
    }
  {
    config,
    topic_gene,
    cell_topic,
    encoder_w: init_w(n_genes, config.n_topics),
    encoder_b: init_b(config.n_topics),
    n_genes,
    n_cells,
    trained: false,
  }
}

///|
/// Get topic proportions for a cell (amortized inference)
pub fn lda_get_topic_proportions(
  model : LDAModel,
  x : Array[Double],
  cell_idx : Int32,
) -> Array[Double] {
  if cell_idx < model.cell_topic.length() {
    return model.cell_topic[cell_idx]
  }
  // Encode
  let logits = scvi_matvec_linear(model.encoder_w, x, model.encoder_b)
  // Softmax
  let max_l = scvi_array_max(logits)
  let mut sum_exp = 0.0
  for i in 0..<logits.length() {
    sum_exp += @math.exp(logits[i] - max_l)
  }
  let result = Array(capacity=logits.length())
  for i in 0..<logits.length() {
    result.push(@math.exp(logits[i] - max_l) / sum_exp)
  }
  result
}

///|
/// Get topic-gene distribution
pub fn lda_get_topic_genes(model : LDAModel) -> Array[Array[Double]] {
  model.topic_gene
}

///|
/// Get dominant topic for each cell
pub fn lda_get_dominant_topics(model : LDAModel) -> Array[Int32] {
  let result = Array(capacity=model.n_cells)
  for c in 0..<model.n_cells {
    let mut best_topic = 0
    let mut best_prob = model.cell_topic[c][0]
    for t in 1..<model.config.n_topics {
      if model.cell_topic[c][t] > best_prob {
        best_prob = model.cell_topic[c][t]
        best_topic = t
      }
    }
    result.push(best_topic)
  }
  result
}

///|
/// Get top genes per topic
pub fn lda_get_top_genes_per_topic(
  model : LDAModel,
  n_top : Int32,
) -> Array[Array[Int32]] {
  let result = Array(capacity=model.config.n_topics)
  for t in 0..<model.config.n_topics {
    let indexed = Array(capacity=model.n_genes)
      |> a => {
        for g in 0..<model.n_genes {
          a.push((g, model.topic_gene[t][g]))
        }
        a
      }
    indexed.sort_by(fn(a, b) {
      if a.1 > b.1 {
        -1
      } else if a.1 < b.1 {
        1
      } else {
        0
      }
    })
    let top = Array(capacity=n_top)
    for i in 0..<n_top {
      if i < indexed.length() {
        top.push(indexed[i].0)
      }
    }
    result.push(top)
  }
  result
}

///|
/// Compute perplexity
pub fn lda_perplexity(model : LDAModel, x : Array[Array[Double]]) -> Double {
  let mut total_log_lik = 0.0
  let mut total_count = 0.0
  for c in 0..<x.length() {
    for g in 0..<model.n_genes {
      let mut expected = 0.0
      for t in 0..<model.config.n_topics {
        expected += model.cell_topic[c][t] * model.topic_gene[t][g]
      }
      if x[c][g] > 0.0 && expected > 0.0 {
        total_log_lik += x[c][g] * @math.ln(expected + 0.0000000001)
        total_count += x[c][g]
      }
    }
  }
  if total_count > 0.0 {
    @math.exp(-total_log_lik / total_count)
  } else {
    100000000000000000000.0
  }
}

///|
/// Train LDA model (variational EM, simplified)
pub fn lda_train(
  model : LDAModel,
  x : Array[Array[Double]],
  n_epochs : Int32,
) -> LDAModel {
  let m = model
  for epoch in 0..<n_epochs {
    // E-step: update cell-topic distributions
    for c in 0..<m.n_cells {
      let lib = scvi_array_sum(x[c])
      for t in 0..<m.config.n_topics {
        let mut score = 0.0
        for g in 0..<m.n_genes {
          if x[c][g] > 0.0 && m.topic_gene[t][g] > 0.0 {
            score += x[c][g] * m.topic_gene[t][g] / (lib + 0.0000000001)
          }
        }
        m.cell_topic[c][t] = score + 0.01 // Smoothing
      }
      // Normalize
      let mut total = 0.0
      for t in 0..<m.config.n_topics {
        total += m.cell_topic[c][t]
      }
      if total > 0.0 {
        for t in 0..<m.config.n_topics {
          m.cell_topic[c][t] /= total
        }
      }
    }
    // M-step: update topic-gene distributions
    for t in 0..<m.config.n_topics {
      for g in 0..<m.n_genes {
        let mut score = 0.0
        for c in 0..<m.n_cells {
          if x[c][g] > 0.0 {
            score += m.cell_topic[c][t] * x[c][g]
          }
        }
        m.topic_gene[t][g] = score + 0.01 // Smoothing
      }
      // Normalize
      let mut total = 0.0
      for g in 0..<m.n_genes {
        total += m.topic_gene[t][g]
      }
      if total > 0.0 {
        for g in 0..<m.n_genes {
          m.topic_gene[t][g] /= total
        }
      }
    }
    ignore(epoch)
  }
  m.trained = true
  m
}

///|
pub fn LDAConfig::default() -> LDAConfig {
  { n_topics: 0, learning_rate: 0.0, max_epochs: 0, }
}
