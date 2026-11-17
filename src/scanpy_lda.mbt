///|
/// LDA topic model over expression counts — `scanpy.external` style topic
/// modeling of cells ("documents" = cells, "words" = genes). Collapsed Gibbs
/// sampling (Griffiths & Steyvers 2004) with a deterministic LCG so results
/// are reproducible. Expression values are floored to integer counts
/// (continuous values contribute their integer part, as in Latent Dirichlet
/// Allocation on count matrices).

///|
/// LDA result: `doc_topic` is n_obs × n_topics (per-cell topic mixture),
/// `topic_term` is n_topics × n_vars (per-topic gene distribution). Both are
/// normalized probability rows.
pub struct ScanpyLda {
  doc_topic : Array[Array[Double]]
  topic_term : Array[Array[Double]]
  n_topics : Int32
  n_iter : Int32
}

///|
/// LDA hyperparameters: symmetric Dirichlet priors over document-topic and
/// topic-term distributions.
let scanpy_lda_alpha : Double = 0.1

///|
let scanpy_lda_beta : Double = 0.01

///|
/// Fit an LDA topic model on the AnnData count matrix (`layer` selects an
/// alternative layer; empty means `x`). Stores the cell-topic mixture as
/// `obsm["X_lda"]` and the gene-topic distribution (n_vars × n_topics) as
/// `varm["lda_topics"]`, and returns the full result.
pub fn scanpy_lda_topics(
  adata : AnnData,
  n_topics~ : Int32 = 5,
  n_iter~ : Int32 = 100,
  seed~ : Int32 = 1,
  layer~ : String = "",
) -> ScanpyLda {
  let x = if layer.length() > 0 && adata.layers.contains(layer) {
    adata.layers[layer]
  } else {
    adata.x_matrix()
  }
  let res = scanpy_lda_fit(x, n_topics~, n_iter~, seed~)
  let _ = adata.set_obsm("X_lda", res.doc_topic)
  // varm stores n_vars × n_topics (transpose of topic_term)
  let v = adata.n_vars()
  let kt = res.n_topics
  let varm_mat : Array[Array[Double]] = Array(capacity=v)
  for g = 0; g < v; g = g + 1 {
    let row : Array[Double] = Array(capacity=kt)
    for t = 0; t < kt; t = t + 1 {
      row.push(res.topic_term[t][g])
    }
    varm_mat.push(row)
  }
  let _ = adata.set_varm("lda_topics", varm_mat)
  res
}

///|
/// LDA (collapsed Gibbs) on a plain doc × term count matrix.
pub fn scanpy_lda_fit(
  x : Array[Array[Double]],
  n_topics~ : Int32 = 5,
  n_iter~ : Int32 = 100,
  seed~ : Int32 = 1,
) -> ScanpyLda {
  let d_docs = x.length()
  let v = if d_docs > 0 { x[0].length() } else { 0 }
  let k = if n_topics < 1 { 1 } else { n_topics }
  let iters = if n_iter < 1 { 1 } else { n_iter }
  let doc_topic : Array[Array[Double]] = Array(capacity=d_docs)
  let topic_term : Array[Array[Double]] = Array(capacity=k)
  for _d in 0..<d_docs {
    doc_topic.push(Array::make(k, 0.0))
  }
  for _t in 0..<k {
    topic_term.push(Array::make(v, 0.0))
  }
  if d_docs == 0 || v == 0 {
    return ScanpyLda::{ doc_topic, topic_term, n_topics: k, n_iter: iters, }
  }
  // token list per document: (gene, topic assignment)
  let doc_tokens : Array[Array[(Int32, Int32)]] = Array(capacity=d_docs)
  for d = 0; d < d_docs; d = d + 1 {
    let toks : Array[(Int32, Int32)] = Array(capacity=0)
    for g = 0; g < v; g = g + 1 {
      let c = x[d][g].floor().to_int()
      for _r in 0..<c {
        toks.push((g, 0))
      }
    }
    doc_tokens.push(toks)
  }
  // count tables
  let n_wt : Array[Array[Int32]] = Array(capacity=v)
  for _g in 0..<v {
    n_wt.push(Array::make(k, 0))
  }
  let n_dt : Array[Array[Int32]] = Array(capacity=d_docs)
  for _d in 0..<d_docs {
    n_dt.push(Array::make(k, 0))
  }
  let n_t = Array::make(k, 0)
  // initialize with the deterministic LCG
  let mut rng = if seed <= 0 { 1 } else { seed }
  for d = 0; d < d_docs; d = d + 1 {
    let toks = doc_tokens[d]
    for i = 0; i < toks.length(); i = i + 1 {
      rng = skbd_lcg(rng)
      let t = rng % k
      toks[i] = (toks[i].0, t)
      n_wt[toks[i].0][t] = n_wt[toks[i].0][t] + 1
      n_dt[d][t] = n_dt[d][t] + 1
      n_t[t] = n_t[t] + 1
    }
  }
  let v_beta = v.to_double() * scanpy_lda_beta
  // collapsed Gibbs iterations
  for _it in 0..<iters {
    for d = 0; d < d_docs; d = d + 1 {
      let toks = doc_tokens[d]
      for i = 0; i < toks.length(); i = i + 1 {
        let g = toks[i].0
        let old = toks[i].1
        n_wt[g][old] = n_wt[g][old] - 1
        n_dt[d][old] = n_dt[d][old] - 1
        n_t[old] = n_t[old] - 1
        // sample new topic ∝ (n_wt + β)(n_dt + α)/(n_t + Vβ)
        let mut total = 0.0
        let probs : Array[Double] = Array(capacity=k)
        for t = 0; t < k; t = t + 1 {
          let p = (n_wt[g][t].to_double() + scanpy_lda_beta) *
            (n_dt[d][t].to_double() + scanpy_lda_alpha) /
            (n_t[t].to_double() + v_beta)
          probs.push(p)
          total = total + p
        }
        rng = skbd_lcg(rng)
        let r = (rng % 1000000).to_double() / 1000000.0 * total
        let mut acc = 0.0
        let mut chosen = k - 1
        for t = 0; t < k; t = t + 1 {
          acc = acc + probs[t]
          if r < acc {
            chosen = t
            break
          }
        }
        toks[i] = (g, chosen)
        n_wt[g][chosen] = n_wt[g][chosen] + 1
        n_dt[d][chosen] = n_dt[d][chosen] + 1
        n_t[chosen] = n_t[chosen] + 1
      }
    }
  }
  // normalized posteriors
  for d = 0; d < d_docs; d = d + 1 {
    let mut s = 0.0
    for t = 0; t < k; t = t + 1 {
      s = s + n_dt[d][t].to_double() + scanpy_lda_alpha
    }
    for t = 0; t < k; t = t + 1 {
      doc_topic[d][t] = (n_dt[d][t].to_double() + scanpy_lda_alpha) / s
    }
  }
  for t = 0; t < k; t = t + 1 {
    let s = n_t[t].to_double() + v_beta
    for g = 0; g < v; g = g + 1 {
      topic_term[t][g] = (n_wt[g][t].to_double() + scanpy_lda_beta) / s
    }
  }
  ScanpyLda::{ doc_topic, topic_term, n_topics: k, n_iter: iters, }
}
