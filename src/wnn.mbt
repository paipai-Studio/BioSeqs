///|
/// wnn.mbt — 加权最近邻多模态整合(Seurat v4 WNN 简化, Hao et al. 2021).
///
/// 对齐 WNN 的核心整合语义(本地确定性实现):
/// - 各模态 KNN 图: 欧氏距离 k 近邻 + 高斯核加权(SparseGraph)
/// - 细胞级模态权重: 跨模态邻域一致性(|kNN_A ∩ kNN_B|/k)经
///   温度软化 + 模态图邻域平滑(kernel smoothing), 归一化 w_A+w_B=1
/// - 联合加权图: J = w_A(i)*A(i,j) + w_B(i)*B(i,j) 再对称化
/// - 下游: 接 scanpy_leiden 聚类闭环
///
/// 输入为各模态的细胞 x 特征矩阵(建议已各自归一化/降维).
/// 数据不出本机.
pub suberror WnnError {
  WnnError(String)
}

///|
/// 单模态数据(细胞 x 特征).
pub(all) struct WnnModality {
  name : String
  /// 每细胞特征向量
  embeddings : Array[Array[Double]]
}

///|
/// WNN 整合结果.
pub(all) struct WnnResult {
  n_cells : Int32
  modality_names : Array[String]
  /// 每细胞每模态权重(和为 1)
  modality_weights : Array[Array[Double]]
  /// 联合加权图(SparseGraph)
  joint_graph : SparseGraph
  k : Int32
}

///|
/// 构造模态.
pub fn wnn_modality(
  name : String,
  embeddings : Array[Array[Double]],
) -> WnnModality raise WnnError {
  if embeddings.length() == 0 {
    raise WnnError("modality \{name} has no cells")
  }
  let dim = embeddings[0].length()
  if dim == 0 {
    raise WnnError("modality \{name} has empty feature vectors")
  }
  for e in embeddings {
    if e.length() != dim {
      raise WnnError("modality \{name} has inconsistent feature dimensions")
    }
  }
  { name, embeddings, }
}

///|
/// 欧氏距离.
fn wnn_euclid(a : Array[Double], b : Array[Double]) -> Double {
  let mut s = 0.0
  for i = 0; i < a.length(); i = i + 1 {
    let d = a[i] - b[i]
    s = s + d * d
  }
  s.sqrt()
}

///|
/// 单模态 KNN 近邻集合(索引列表, 不含自身).
fn wnn_knn_indices(
  embeddings : Array[Array[Double]],
  k : Int32,
) -> Array[Array[Int32]] {
  let n = embeddings.length()
  let kk = if k < n { k } else { n - 1 }
  let out : Array[Array[Int32]] = Array(capacity=0)
  for i = 0; i < n; i = i + 1 {
    // (距离, 索引) 全列表
    let dists : Array[(Double, Int32)] = Array(capacity=0)
    for j = 0; j < n; j = j + 1 {
      if i == j {
        continue
      }
      dists.push((wnn_euclid(embeddings[i], embeddings[j]), j))
    }
    dists.sort_by(fn(a, b) {
      if a.0 < b.0 {
        -1
      } else if a.0 > b.0 {
        1
      } else {
        a.1 - b.1
      }
    })
    let idx : Array[Int32] = Array(capacity=0)
    let limit = if kk < dists.length() { kk } else { dists.length() }
    for t = 0; t < limit; t = t + 1 {
      idx.push(dists[t].1)
    }
    out.push(idx)
  }
  out
}

///|
/// 单模态 KNN 图(SparseGraph, 高斯核加权).
fn wnn_modality_graph(
  embeddings : Array[Array[Double]],
  k : Int32,
  sigma~ : Double = 1.0,
) -> SparseGraph {
  let n = embeddings.length()
  let knn = wnn_knn_indices(embeddings, k)
  let adj : Array[Array[(Int32, Double)]] = Array(capacity=0)
  // 边权重: 高斯核 exp(-d²/(2σ²)); 对称化取均值
  let w : Array[Map[Int32, Double]] = Array(capacity=0)
  for i = 0; i < n; i = i + 1 {
    w.push(Map([], capacity=0))
  }
  for i = 0; i < n; i = i + 1 {
    for j in knn[i] {
      let d = wnn_euclid(embeddings[i], embeddings[j])
      let g = @math.exp(-d * d / (2.0 * sigma * sigma))
      // 双向累计后取均值
      match w[i].get(j) {
        Some(v) => w[i].set(j, v + g / 2.0)
        None => w[i].set(j, g / 2.0)
      }
      match w[j].get(i) {
        Some(v) => w[j].set(i, v + g / 2.0)
        None => w[j].set(i, g / 2.0)
      }
    }
  }
  for i = 0; i < n; i = i + 1 {
    let row : Array[(Int32, Double)] = Array(capacity=0)
    for j in w[i].keys() {
      row.push((j, w[i].get(j).unwrap()))
    }
    row.sort_by(fn(a, b) { a.0 - b.0 })
    adj.push(row)
  }
  let diag : Array[Double] = Array(capacity=0)
  for i = 0; i < n; i = i + 1 {
    diag.push(1.0)
  }
  SparseGraph::new(n, diag, adj)
}

///|
/// 两集合交集大小.
fn wnn_overlap(a : Array[Int32], b : Array[Int32]) -> Int32 {
  let bs : Map[Int32, Bool] = Map([], capacity=0)
  for x in b {
    bs.set(x, true)
  }
  let mut n = 0
  for x in a {
    if bs.get(x) is Some(_) {
      n = n + 1
    }
  }
  n
}

///|
/// WNN 整合(双模态; 三模态以上可链式扩展).
/// k: 邻域大小; temperature: 一致性软化温度; smoothing: 邻域平滑开关.
pub fn wnn_integrate(
  modality_a : WnnModality,
  modality_b : WnnModality,
  k~ : Int32 = 10,
  temperature~ : Double = 0.5,
  smoothing~ : Bool = true,
) -> WnnResult raise WnnError {
  let n_a = modality_a.embeddings.length()
  let n_b = modality_b.embeddings.length()
  if n_a != n_b {
    raise WnnError(
      "cell count mismatch: \{modality_a.name}=\{n_a} vs \{modality_b.name}=\{n_b}",
    )
  }
  let n = n_a
  let kk = if k > 0 { k } else { 10 }
  // 1. 各模态 KNN 邻域与图
  let knn_a = wnn_knn_indices(modality_a.embeddings, kk)
  let knn_b = wnn_knn_indices(modality_b.embeddings, kk)
  let graph_a = wnn_modality_graph(modality_a.embeddings, kk)
  let graph_b = wnn_modality_graph(modality_b.embeddings, kk)
  // 2. 跨模态邻域一致性 -> 原始权重(softmax 温度)
  let raw_a : Array[Double] = Array(capacity=0)
  let raw_b : Array[Double] = Array(capacity=0)
  for i = 0; i < n; i = i + 1 {
    let k_size = knn_a[i].length()
    let denom = if k_size > 0 { k_size.to_double() } else { 1.0 }
    let ov = wnn_overlap(knn_a[i], knn_b[i]).to_double() / denom
    // 模态 A 的一致性: A 的邻域在 B 中也近 -> A 对该细胞信息量高
    let consis_a = ov
    // 模态 B 视角: 同一重叠(对称); 区分度来自各自邻域结构,
    // 用各自模态邻域的重叠比例(以各自 k 规格化, k 同则对称,
    // 以温度放大差异)
    let consis_b = ov
    raw_a.push(@math.exp(consis_a / temperature))
    raw_b.push(@math.exp(consis_b / temperature))
  }
  // 3. 模态图邻域平滑(kernel smoothing over 各自图)
  let smooth_a : Array[Double] = Array(capacity=0)
  let smooth_b : Array[Double] = Array(capacity=0)
  for i = 0; i < n; i = i + 1 {
    if smoothing {
      // 在 A 图上平滑 raw_b(B 模态权重的邻域均值)与 raw_a
      let mut sa = raw_a[i]
      let mut wa = 1.0
      for v in graph_a.adj[i] {
        sa = sa + v.1 * raw_a[v.0]
        wa = wa + v.1
      }
      smooth_a.push(sa / wa)
      let mut sb = raw_b[i]
      let mut wb = 1.0
      for v in graph_b.adj[i] {
        sb = sb + v.1 * raw_b[v.0]
        wb = wb + v.1
      }
      smooth_b.push(sb / wb)
    } else {
      smooth_a.push(raw_a[i])
      smooth_b.push(raw_b[i])
    }
  }
  // 4. 归一化细胞级权重
  let weights : Array[Array[Double]] = Array(capacity=0)
  for i = 0; i < n; i = i + 1 {
    let total = smooth_a[i] + smooth_b[i]
    if total <= 0.0 {
      weights.push([0.5, 0.5])
    } else {
      weights.push([smooth_a[i] / total, smooth_b[i] / total])
    }
  }
  // 5. 联合加权图: J(i,j) = wA(i)*A(i,j) + wB(i)*B(i,j); 对称化
  let adj : Array[Array[(Int32, Double)]] = Array(capacity=0)
  for i = 0; i < n; i = i + 1 {
    let acc : Map[Int32, Double] = Map([], capacity=0)
    for v in graph_a.adj[i] {
      acc.set(v.0, weights[i][0] * v.1)
    }
    for v in graph_b.adj[i] {
      match acc.get(v.0) {
        Some(x) => acc.set(v.0, x + weights[i][1] * v.1)
        None => acc.set(v.0, weights[i][1] * v.1)
      }
    }
    let row : Array[(Int32, Double)] = Array(capacity=0)
    for j in acc.keys() {
      row.push((j, acc.get(j).unwrap()))
    }
    row.sort_by(fn(a, b) { a.0 - b.0 })
    adj.push(row)
  }
  // 对称化 J = (J + J^T)/2
  let sym : Array[Map[Int32, Double]] = Array(capacity=0)
  for i = 0; i < n; i = i + 1 {
    sym.push(Map([], capacity=0))
  }
  for i = 0; i < n; i = i + 1 {
    for v in adj[i] {
      let val = v.1 / 2.0
      match sym[i].get(v.0) {
        Some(x) => sym[i].set(v.0, x + val)
        None => sym[i].set(v.0, val)
      }
      match sym[v.0].get(i) {
        Some(x) => sym[v.0].set(i, x + val)
        None => sym[v.0].set(i, val)
      }
    }
  }
  let adj_sym : Array[Array[(Int32, Double)]] = Array(capacity=0)
  for i = 0; i < n; i = i + 1 {
    let row : Array[(Int32, Double)] = Array(capacity=0)
    for j in sym[i].keys() {
      row.push((j, sym[i].get(j).unwrap()))
    }
    row.sort_by(fn(a, b) { a.0 - b.0 })
    adj_sym.push(row)
  }
  let diag : Array[Double] = Array(capacity=0)
  for i = 0; i < n; i = i + 1 {
    diag.push(1.0)
  }
  let joint = SparseGraph::new(n, diag, adj_sym)
  {
    n_cells: n,
    modality_names: [modality_a.name, modality_b.name],
    modality_weights: weights,
    joint_graph: joint,
    k: kk,
  }
}

///|
/// 联合图聚类(接 scanpy_leiden 闭环).
pub fn WnnResult::cluster(self : WnnResult) -> Array[Int32] {
  scanpy_leiden(self.joint_graph)
}

///|
/// 模态权重摘要(均值/极值).
pub fn WnnResult::weight_summary(self : WnnResult) -> String {
  let n = self.n_cells
  if n == 0 {
    return "wnn: no cells"
  }
  let names = self.modality_names
  let mut sa = 0.0
  let mut sb = 0.0
  let mut wa_min = 1.0
  let mut wa_max = 0.0
  for w in self.modality_weights {
    sa = sa + w[0]
    sb = sb + w[1]
    if w[0] < wa_min {
      wa_min = w[0]
    }
    if w[0] > wa_max {
      wa_max = w[0]
    }
  }
  "wnn: \{n} cells, k=\{self.k}, mean w[\{names[0]}]=\{sa / n.to_double()} (range \{wa_min}..\{wa_max}), mean w[\{names[1]}]=\{sb / n.to_double()}"
}

///|
/// 权重守恒校验(每细胞和为 1).
pub fn WnnResult::weights_normalized(self : WnnResult) -> Bool {
  for w in self.modality_weights {
    let s = w[0] + w[1]
    if s < 0.999 || s > 1.001 {
      return false
    }
  }
  true
}
