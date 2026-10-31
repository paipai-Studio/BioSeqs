///|
/// 空间邻域富集与空间细胞通讯演示.
///
/// 1. 6x6 网格, 左半 A 类 / 右半 B 类 -> 邻域富集 z-score
/// 2. 2x2 布局 sender/receiver -> 空间配体-受体通讯
/// 3. 通讯结果按 source->target 聚合

///|
/// 构造 n x n 网格坐标.
fn demo_grid(n : Int32) -> Array[(Double, Double)] {
  let out : Array[(Double, Double)] = Array(capacity=0)
  let mut i = 0
  while i < n {
    let mut j = 0
    while j < n {
      out.push((i.to_double(), j.to_double()))
      j = j + 1
    }
    i = i + 1
  }
  out
}

///|
/// 左半 A / 右半 B 标签.
fn demo_split_labels(
  coords : Array[(Double, Double)],
  cut : Double,
) -> Array[String] {
  let out : Array[String] = Array(capacity=0)
  for c in coords {
    let (x, _) = c
    if x < cut {
      out.push("A")
    } else {
      out.push("B")
    }
  }
  out
}

///|
/// 主函数.
fn main raise {
  println("=== 空间邻域富集 ===")
  let coords = demo_grid(6)
  let labels = demo_split_labels(coords, 3.0)
  println(
    "spot 数: " +
    coords.length().to_string() +
    ", 类别: A(左半) / B(右半)",
  )
  let enr = @src.nhood_enrichment(coords, labels, k=6, n_perms=50, seed=1)
  println("类别顺序: " + repr(enr.labels()))
  println("z-score 矩阵 (行=源类别, 列=邻域类别):")
  for row in enr.zscore() {
    println("  " + repr(row))
  }
  println("A 邻域内 A 富集 z = " + enr.zscore_of("A", "A").to_string())
  println("A 邻域内 B 富集 z = " + enr.zscore_of("A", "B").to_string())

  println("\n=== 空间细胞通讯 ===")
  let comm_coords : Array[(Double, Double)] = [
    (0.0, 0.0),
    (1.0, 0.0),
    (0.0, 1.0),
    (1.0, 1.0),
  ]
  let comm_labels : Array[String] = ["S", "S", "R", "R"]
  let expr : Array[Array[Double]] = [
    [5.0, 0.0],
    [4.0, 0.0],
    [0.0, 3.0],
    [0.0, 2.0],
  ]
  let genes : Array[String] = ["LIG", "REC"]
  let pairs : Array[@src.LRPair] = [@src.lr_pair("LIG", "REC")]
  let scores = @src.spatial_cellchat(
    comm_coords,
    comm_labels,
    expr,
    genes,
    pairs,
    k=2,
  )
  println("检出通讯条目数: " + scores.length().to_string())
  for s in scores {
    println(
      "  " +
      s.source() +
      " -> " +
      s.target() +
      "  " +
      s.ligand() +
      "-" +
      s.receptor() +
      "  pairs=" +
      s.n_pairs().to_string() +
      "  score=" +
      s.score().to_string(),
    )
  }
  let agg = @src.spatial_cellchat_aggregate(scores)
  println("按 source->target 聚合:")
  for key in agg.keys() {
    println("  " + key + " : " + agg[key].to_string())
  }
  println("\n演示完成")
}
