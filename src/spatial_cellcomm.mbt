///|
/// 空间感知细胞通讯(CellChat spatial 语义简化).
///
/// 在通用 CellChat 配体-受体打分基础上引入空间约束: 只有空间邻域内
/// 同时表达配体(sender)与受体(receiver)的细胞对计数, 通讯强度按
/// 类型平均表达与空间支持对数量计算.
/// 数据不出本机.

///|
/// 空间细胞通讯错误.
pub suberror SpatialCommError {
  SpatialCommError(String)
}

///|
/// 单条空间配体-受体通讯结果.
pub struct SpatialLRScore {
  ligand : String
  receptor : String
  source : String
  target : String
  n_pairs : Int32
  mean_ligand : Double
  mean_receptor : Double
  score : Double
}

///|
/// 配体.
pub fn SpatialLRScore::ligand(self : SpatialLRScore) -> String {
  self.ligand
}

///|
/// 受体.
pub fn SpatialLRScore::receptor(self : SpatialLRScore) -> String {
  self.receptor
}

///|
/// 发送方细胞类型.
pub fn SpatialLRScore::source(self : SpatialLRScore) -> String {
  self.source
}

///|
/// 接收方细胞类型.
pub fn SpatialLRScore::target(self : SpatialLRScore) -> String {
  self.target
}

///|
/// 空间支持细胞对数量.
pub fn SpatialLRScore::n_pairs(self : SpatialLRScore) -> Int32 {
  self.n_pairs
}

///|
/// 发送方配体平均表达.
pub fn SpatialLRScore::mean_ligand(self : SpatialLRScore) -> Double {
  self.mean_ligand
}

///|
/// 接收方受体平均表达.
pub fn SpatialLRScore::mean_receptor(self : SpatialLRScore) -> Double {
  self.mean_receptor
}

///|
/// 通讯强度得分.
pub fn SpatialLRScore::score(self : SpatialLRScore) -> Double {
  self.score
}

///|
/// 基因名到表达矩阵列索引.
fn spatial_gene_index(genes : Array[String]) -> Map[String, Int32] {
  let idx : Map[String, Int32] = Map([], capacity=0)
  let mut i = 0
  while i < genes.length() {
    idx[genes[i]] = i
    i = i + 1
  }
  idx
}

///|
/// 取某类细胞的配体/受体平均表达(未表达的细胞计 0).
fn spatial_type_mean(
  expr : Array[Array[Double]],
  cell_idx : Array[Int32],
  col : Int32,
) -> Double {
  if cell_idx.length() == 0 {
    return 0.0
  }
  let mut sum = 0.0
  for c in cell_idx {
    sum = sum + expr[c][col]
  }
  sum / cell_idx.length().to_double()
}

///|
/// 空间细胞通讯分析.
///
/// 参数:
/// - coords: 细胞空间坐标
/// - labels: 细胞类型标签, 顺序与 coords 一致
/// - expr: 表达矩阵(细胞 x 基因)
/// - genes: 基因名, 顺序与 expr 列一致
/// - lr_pairs: 配体-受体对(可复用 cellchat_lr_database)
/// - k: 邻域大小
/// - min_pairs: 最少空间支持细胞对(低于该值丢弃)
/// 结果按 score 降序.
pub fn spatial_cellchat(
  coords : Array[(Double, Double)],
  labels : Array[String],
  expr : Array[Array[Double]],
  genes : Array[String],
  lr_pairs : Array[LRPair],
  k~ : Int32 = 6,
  min_pairs~ : Int32 = 1,
) -> Array[SpatialLRScore] raise SpatialCommError {
  let n = coords.length()
  if labels.length() != n || expr.length() != n {
    raise SpatialCommError("coords/labels/expr 行数不一致")
  }
  let edges = spatial_knn_edges(coords, k)
  let gene_idx = spatial_gene_index(genes)
  // 细胞类型去重与分桶
  let types : Array[String] = Array(capacity=0)
  let type_map : Map[String, Int32] = Map([], capacity=0)
  for l in labels {
    if type_map.get(l) is None {
      type_map[l] = types.length()
      types.push(l)
    }
  }
  let members : Array[Array[Int32]] = Array(capacity=types.length())
  for _ in 0..<types.length() {
    let bucket : Array[Int32] = Array(capacity=0)
    members.push(bucket)
  }
  let mut ci = 0
  while ci < n {
    match type_map.get(labels[ci]) {
      Some(t) => members[t].push(ci)
      None => ()
    }
    ci = ci + 1
  }
  let out : Array[SpatialLRScore] = Array(capacity=0)
  for pair in lr_pairs {
    let li = match gene_idx.get(pair.ligand) {
      Some(v) => v
      None => -1
    }
    let ri = match gene_idx.get(pair.receptor) {
      Some(v) => v
      None => -1
    }
    if li < 0 || ri < 0 {
      continue
    }
    let mut s = 0
    while s < types.length() {
      let mut t = 0
      while t < types.length() {
        // 统计 s 类细胞 -> t 类邻居 的空间支持对
        let mut pairs = 0
        for i in members[s] {
          if expr[i][li] <= 0.0 {
            continue
          }
          for nb in edges[i] {
            if labels[nb] == types[t] && expr[nb][ri] > 0.0 {
              pairs = pairs + 1
            }
          }
        }
        if pairs >= min_pairs && pairs > 0 {
          let ml = spatial_type_mean(expr, members[s], li)
          let mr = spatial_type_mean(expr, members[t], ri)
          out.push(SpatialLRScore::{
            ligand: pair.ligand,
            receptor: pair.receptor,
            source: types[s],
            target: types[t],
            n_pairs: pairs,
            mean_ligand: ml,
            mean_receptor: mr,
            score: ml * mr,
          })
        }
        t = t + 1
      }
      s = s + 1
    }
  }
  // 按 score 降序
  let mut a = 0
  while a < out.length() {
    let mut b = a + 1
    while b < out.length() {
      if out[b].score > out[a].score {
        let tmp = out[a]
        out[a] = out[b]
        out[b] = tmp
      }
      b = b + 1
    }
    a = a + 1
  }
  out
}

///|
/// 汇总通讯结果到 (source, target) 的通讯总量.
pub fn spatial_cellchat_aggregate(
  scores : Array[SpatialLRScore],
) -> Map[String, Double] {
  let agg : Map[String, Double] = Map([], capacity=0)
  for s in scores {
    let key = s.source + "->" + s.target
    match agg.get(key) {
      Some(v) => agg[key] = v + s.score
      None => agg[key] = s.score
    }
  }
  agg
}
