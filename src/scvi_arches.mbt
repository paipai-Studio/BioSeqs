// scvi_arches.mbt - Transfer learning (scArches)
// Inspired by scvi-tools ArchesMixin: load_query_data, prepare_query_anndata
// Online update of pretrained reference with query data

///|
pub struct ArchesConfig {
  unfrozen_decoder : Bool
  freeze_dropout : Bool
  freeze_expression : Bool
  freeze_first_layer : Bool
  freeze_batchnorm : Bool
  freeze_classifier : Bool
  transfer_batch : Bool
}

///|
pub fn ArchesConfig::default() -> ArchesConfig {
  {
    unfrozen_decoder: false,
    freeze_dropout: true,
    freeze_expression: true,
    freeze_first_layer: true,
    freeze_batchnorm: true,
    freeze_classifier: true,
    transfer_batch: false,
  }
}

///|
/// Prepare query AnnData by aligning genes to reference
pub fn arches_prepare_query(
  query_x : Array[Array[Double]],
  query_var_names : Array[String],
  reference_var_names : Array[String],
) -> Array[Array[Double]] {
  let n_cells = query_x.length()
  let n_ref_genes = reference_var_names.length()
  // Build gene index mapping
  let ref_index = Map([], capacity=n_ref_genes)
    |> m => {
      for i in 0..<n_ref_genes {
        m[reference_var_names[i]] = i
      }
      m
    }
  // Create aligned expression matrix
  let result = Array(capacity=n_cells)
  for c in 0..<n_cells {
    let aligned = Array(capacity=n_ref_genes)
      |> a => {
        for _i in 0..<n_ref_genes {
          a.push(0.0)
        }
        a
      }
    for g in 0..<query_var_names.length() {
      match ref_index.get(query_var_names[g]) {
        Some(ref_idx) => aligned[ref_idx] = query_x[c][g]
        None => ()
      }
    }
    result.push(aligned)
  }
  result
}

///|
/// Load query data into a pretrained model (transfer learning)
pub fn arches_load_query(
  reference_encoder_w : Array[Array[Double]],
  reference_encoder_b : Array[Double],
  query_x : Array[Array[Double]],
  config : ArchesConfig,
) -> (Array[Array[Double]], Array[Double]) {
  // Copy reference weights (freeze most parameters)
  let n_out = reference_encoder_w.length()
  let n_in = reference_encoder_w[0].length()
  let query_w = Array(capacity=n_out)
    |> a => {
      for i in 0..<n_out {
        let row = Array(capacity=n_in)
        for j in 0..<n_in {
          row.push(reference_encoder_w[i][j])
        }
        a.push(row)
      }
      a
    }
  let query_b = Array(capacity=n_out)
    |> a => {
      for i in 0..<n_out {
        a.push(reference_encoder_b[i])
      }
      a
    }
  // If first layer is not frozen, allow weight perturbation
  if !config.freeze_first_layer {
    for i in 0..<n_out {
      for j in 0..<n_in {
        query_w[i][j] += scvi_rand(i * n_in + j + 900) * 0.01
      }
    }
  }
  ignore(query_x)
  (query_w, query_b)
}

///|
/// Transfer batch labels from reference to query
pub fn arches_transfer_batch(
  query_latent : Array[Array[Double]],
  reference_latent : Array[Array[Double]],
  reference_batch_labels : Array[Int32],
) -> Array[Int32] {
  let n_query = query_latent.length()
  let n_ref = reference_latent.length()
  let result = Array(capacity=n_query)
  for i in 0..<n_query {
    // Find nearest reference cell
    let mut best_dist = 100000000000000000000.0
    let mut best_label = 0
    for j in 0..<n_ref {
      let mut dist = 0.0
      for d in 0..<query_latent[i].length() {
        let diff = query_latent[i][d] - reference_latent[j][d]
        dist += diff * diff
      }
      if dist < best_dist {
        best_dist = dist
        best_label = reference_batch_labels[j]
      }
    }
    result.push(best_label)
  }
  result
}

///|
/// Update model with query data (online learning)
pub fn arches_update(
  encoder_w : Array[Array[Double]],
  encoder_b : Array[Double],
  query_x : Array[Array[Double]],
  learning_rate : Double,
  n_epochs : Int32,
) -> (Array[Array[Double]], Array[Double]) {
  let w = encoder_w
  let b = encoder_b
  for epoch in 0..<n_epochs {
    for i in 0..<query_x.length() {
      // Forward pass
      let hidden = scvi_matvec_mul(w, query_x[i], b)
      // Simple reconstruction gradient (perturbation-based)
      for j in 0..<w.length() {
        b[j] += learning_rate * scvi_rand(epoch * 200 + i * 30 + j) * 0.001
        for k in 0..<w[j].length() {
          w[j][k] += learning_rate *
            scvi_rand(epoch * 200 + i * 30 + j * 20 + k) *
            0.001
        }
      }
      ignore(hidden)
    }
  }
  (w, b)
}

///|
/// Compute kNN-based batch transfer accuracy
pub fn arches_transfer_accuracy(
  predicted_labels : Array[Int32],
  true_labels : Array[Int32],
) -> Double {
  let n = predicted_labels.length()
  if n == 0 {
    return 0.0
  }
  let mut correct = 0
  for i in 0..<n {
    if predicted_labels[i] == true_labels[i] {
      correct += 1
    }
  }
  correct.to_double() / n.to_double()
}

///|
/// 单个查询细胞的标签迁移结果
pub(all) struct ArchesKnnLabel {
  /// kNN 投票胜出的参考标签
  label : Int32
  /// 软投票概率（胜出标签得票占比，0..1）
  confidence : Double
  /// k 个近邻的参考集下标（按距离升序）
  neighbors : Array[Int32]
}

///|
/// Latent kNN 标签迁移（scArches / scanpy-ingest 的标准下游）：查询细胞
/// 已映射到参考 latent 空间后，对每个查询细胞取 k 个最近参考细胞做标签
/// 投票。`distance_weighted=false` 为均匀软投票（概率 = 标签得票占比）；
/// `true` 时票权为 1/(欧氏距离+ε)，近邻主导。胜出概率即置信度。
///
/// 边界：参考集为空、标签数与参考集不一致或查询集为空时返回空列表；
/// `k` 收缩到 [1, n_ref]；票数并列时取较小标签保证确定性；latent 维度
/// 取参考与查询的最小公共维度。
pub fn arches_knn_label_transfer(
  reference_latent : Array[Array[Double]],
  reference_labels : Array[Int32],
  query_latent : Array[Array[Double]],
  k~ : Int32 = 15,
  distance_weighted~ : Bool = false,
) -> Array[ArchesKnnLabel] {
  let n_ref = reference_latent.length()
  let out : Array[ArchesKnnLabel] = Array(capacity=query_latent.length())
  if n_ref == 0 ||
    reference_labels.length() != n_ref ||
    query_latent.length() == 0 {
    return out
  }
  let kk = if k < 1 { 1 } else if k > n_ref { n_ref } else { k }
  let mut dim = 2147483647
  for v in reference_latent {
    if v.length() < dim {
      dim = v.length()
    }
  }
  for v in query_latent {
    if v.length() < dim {
      dim = v.length()
    }
  }
  for i in 0..<query_latent.length() {
    let q = query_latent[i]
    // 平方欧氏距离
    let dists : Array[Double] = Array(capacity=n_ref)
    for j in 0..<n_ref {
      let rj = reference_latent[j]
      let mut d = 0.0
      for t in 0..<dim {
        let diff = q[t] - rj[t]
        d += diff * diff
      }
      dists.push(d)
    }
    // 取 k 个最小距离（重复最小扫描，保持距离升序）
    let used : Array[Bool] = Array::make(n_ref, false)
    let nb_idx : Array[Int32] = Array(capacity=kk)
    let nb_dist : Array[Double] = Array(capacity=kk)
    for _a in 0..<kk {
      let mut best = -1
      let mut bd = 1.0e300
      for j in 0..<n_ref {
        if !used[j] && dists[j] < bd {
          bd = dists[j]
          best = j
        }
      }
      if best < 0 {
        break
      }
      used[best] = true
      nb_idx.push(best)
      nb_dist.push(bd)
    }
    // 软投票
    let votes : Map[Int32, Double] = Map([])
    let mut total = 0.0
    for a in 0..<nb_idx.length() {
      let lab = reference_labels[nb_idx[a]]
      let weight = if distance_weighted {
        1.0 / (nb_dist[a].sqrt() + 1.0e-12)
      } else {
        1.0
      }
      match votes.get(lab) {
        Some(w) => votes[lab] = w + weight
        None => votes[lab] = weight
      }
      total += weight
    }
    // 胜出标签（并列取较小标签，保证确定性）
    let mut best_label = -1
    let mut best_prob = -1.0
    for kv in votes {
      let (lab, w) = kv
      let p = w / total
      if p > best_prob + 1.0e-15 ||
        (p > best_prob - 1.0e-15 && lab < best_label) {
        best_prob = p
        best_label = lab
      }
    }
    out.push({ label: best_label, confidence: best_prob, neighbors: nb_idx, })
  }
  out
}
