///|
/// 空间邻域富集分析(squidpy nhood_enrichment 语义).
///
/// 基于空间坐标构建 kNN 邻域图, 统计各类别细胞邻域内出现其他类别的频数,
/// 并通过标签置换构建零分布, 输出类别 x 类别 z-score 矩阵.
/// 数据不出本机.

///|
/// 邻域富集分析错误.
pub suberror NhoodError {
  NhoodError(String)
}

///|
/// 邻域富集结果.
pub struct NhoodEnrichment {
  labels : Array[String]
  observed : Array[Array[Double]]
  perm_mean : Array[Array[Double]]
  perm_std : Array[Array[Double]]
  zscore : Array[Array[Double]]
  k : Int32
  n_perms : Int32
}

///|
/// 类别名(首次出现顺序).
pub fn NhoodEnrichment::labels(self : NhoodEnrichment) -> Array[String] {
  self.labels
}

///|
/// 观察计数矩阵.
pub fn NhoodEnrichment::observed(
  self : NhoodEnrichment,
) -> Array[Array[Double]] {
  self.observed
}

///|
/// 置换零分布均值.
pub fn NhoodEnrichment::perm_mean(
  self : NhoodEnrichment,
) -> Array[Array[Double]] {
  self.perm_mean
}

///|
/// 置换零分布标准差.
pub fn NhoodEnrichment::perm_std(
  self : NhoodEnrichment,
) -> Array[Array[Double]] {
  self.perm_std
}

///|
/// 富集 z-score 矩阵.
pub fn NhoodEnrichment::zscore(self : NhoodEnrichment) -> Array[Array[Double]] {
  self.zscore
}

///|
/// 邻域大小 k.
pub fn NhoodEnrichment::k(self : NhoodEnrichment) -> Int32 {
  self.k
}

///|
/// 置换次数.
pub fn NhoodEnrichment::n_perms(self : NhoodEnrichment) -> Int32 {
  self.n_perms
}

///|
/// 按类别名取索引, 不存在返回 -1.
pub fn NhoodEnrichment::label_index(
  self : NhoodEnrichment,
  name : String,
) -> Int32 {
  let mut i = 0
  while i < self.labels.length() {
    if self.labels[i] == name {
      return i
    }
    i = i + 1
  }
  -1
}

///|
/// 取类别 a 邻域内类别 b 的 z-score.
pub fn NhoodEnrichment::zscore_of(
  self : NhoodEnrichment,
  a : String,
  b : String,
) -> Double {
  let ia = self.label_index(a)
  let ib = self.label_index(b)
  if ia < 0 || ib < 0 {
    return 0.0
  }
  self.zscore[ia][ib]
}

// ---------------------------------------------------------------------------
// kNN 邻域图
// ---------------------------------------------------------------------------

///|
/// 按空间坐标构建 kNN 邻域图(不含自身).
/// 返回每个细胞的邻居索引数组, 邻居按距离升序.
pub fn spatial_knn_edges(
  coords : Array[(Double, Double)],
  k : Int32,
) -> Array[Array[Int32]] {
  let n = coords.length()
  let edges : Array[Array[Int32]] = Array(capacity=n)
  for i in 0..<n {
    let (xi, yi) = coords[i]
    let idx : Array[Int32] = Array(capacity=n)
    let dists : Array[Double] = Array(capacity=n)
    for j in 0..<n {
      if j != i {
        let (xj, yj) = coords[j]
        let dx = xi - xj
        let dy = yi - yj
        idx.push(j)
        dists.push(dx * dx + dy * dy)
      }
    }
    let kk : Int32 = k
    let m = if kk < idx.length() { kk } else { idx.length() }
    for a in 0..<m {
      let mut best = a
      let mut b = a + 1
      while b < idx.length() {
        if dists[b] < dists[best] {
          best = b
        }
        b = b + 1
      }
      if best != a {
        let t = idx[a]
        idx[a] = idx[best]
        idx[best] = t
        let td = dists[a]
        dists[a] = dists[best]
        dists[best] = td
      }
    }
    let out : Array[Int32] = Array(capacity=m)
    for a in 0..<m {
      out.push(idx[a])
    }
    edges.push(out)
  }
  edges
}

// ---------------------------------------------------------------------------
// 邻域富集统计
// ---------------------------------------------------------------------------

///|
/// 构造 n x n 的零矩阵(逐行新建, 避免共享引用).
fn nhood_zero_matrix(n : Int32) -> Array[Array[Double]] {
  let rows : Array[Array[Double]] = Array(capacity=n)
  for _ in 0..<n {
    let row : Array[Double] = Array(capacity=n)
    for _ in 0..<n {
      row.push(0.0)
    }
    rows.push(row)
  }
  rows
}

///|
/// 按 (源类别, 邻居类别) 累计邻域计数.
fn nhood_accumulate(
  edges : Array[Array[Int32]],
  label_idx : Array[Int32],
  mat : Array[Array[Double]],
) -> Unit {
  let mut i = 0
  while i < edges.length() {
    let ci = label_idx[i]
    for nb in edges[i] {
      let cj = label_idx[nb]
      mat[ci][cj] = mat[ci][cj] + 1.0
    }
    i = i + 1
  }
}

///|
/// 空间邻域富集(观察 + 标签置换零分布).
///
/// 参数:
/// - coords: 细胞空间坐标
/// - labels: 每个细胞的类别标签, 顺序须与 coords 一致
/// - k: 邻域大小
/// - n_perms: 置换次数
/// - seed: 置换随机种子(保证可复现)
pub fn nhood_enrichment(
  coords : Array[(Double, Double)],
  labels : Array[String],
  k~ : Int32 = 6,
  n_perms~ : Int32 = 100,
  seed~ : Int32 = 0,
) -> NhoodEnrichment raise NhoodError {
  if coords.length() != labels.length() {
    raise NhoodError("coords 与 labels 长度不一致")
  }
  let edges = spatial_knn_edges(coords, k)
  // 类别去重(保持首次出现顺序)
  let uniq : Array[String] = Array(capacity=0)
  let index : Map[String, Int32] = Map([], capacity=0)
  for l in labels {
    if index.get(l) is None {
      index[l] = uniq.length()
      uniq.push(l)
    }
  }
  let n_labels = uniq.length()
  let label_idx : Array[Int32] = Array(capacity=labels.length())
  for l in labels {
    match index.get(l) {
      Some(v) => label_idx.push(v)
      None => label_idx.push(-1)
    }
  }
  // 观察计数
  let observed = nhood_zero_matrix(n_labels)
  nhood_accumulate(edges, label_idx, observed)
  // 置换零分布
  let sum_mat = nhood_zero_matrix(n_labels)
  let sq_mat = nhood_zero_matrix(n_labels)
  let rng = Rng::from_int(seed)
  let perm : Array[String] = labels.copy()
  let mut p = 0
  while p < n_perms {
    rng.shuffle(perm)
    let pidx : Array[Int32] = Array(capacity=perm.length())
    for l in perm {
      match index.get(l) {
        Some(v) => pidx.push(v)
        None => pidx.push(-1)
      }
    }
    let mat = nhood_zero_matrix(n_labels)
    nhood_accumulate(edges, pidx, mat)
    for a in 0..<n_labels {
      for b in 0..<n_labels {
        let v = mat[a][b]
        sum_mat[a][b] = sum_mat[a][b] + v
        sq_mat[a][b] = sq_mat[a][b] + v * v
      }
    }
    p = p + 1
  }
  let mean_mat = nhood_zero_matrix(n_labels)
  let std_mat = nhood_zero_matrix(n_labels)
  let z_mat = nhood_zero_matrix(n_labels)
  if n_perms > 0 {
    let np = n_perms.to_double()
    for a in 0..<n_labels {
      for b in 0..<n_labels {
        let mean = sum_mat[a][b] / np
        let variance = sq_mat[a][b] / np - mean * mean
        let std = if variance > 0.0 { variance.sqrt() } else { 0.0 }
        mean_mat[a][b] = mean
        std_mat[a][b] = std
        z_mat[a][b] = if std > 0.0 {
          (observed[a][b] - mean) / std
        } else {
          0.0
        }
      }
    }
  }
  NhoodEnrichment::{
    labels: uniq,
    observed,
    perm_mean: mean_mat,
    perm_std: std_mat,
    zscore: z_mat,
    k,
    n_perms,
  }
}

///|
/// 类别 a 到类别 b 的邻域富集是否显著(|z| >= threshold).
pub fn nhood_is_enriched(
  result : NhoodEnrichment,
  a : String,
  b : String,
  threshold~ : Double = 2.0,
) -> Bool {
  let z = result.zscore_of(a, b)
  if z >= threshold {
    true
  } else if z <= -threshold {
    true
  } else {
    false
  }
}
