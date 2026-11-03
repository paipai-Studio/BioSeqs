///|
/// WNN 多模态整合演示: RNA+ADT 双模态 -> 细胞级权重 -> 联合图 -> 聚类.
fn main raise {
  println("=== WNN 加权最近邻多模态整合 (Seurat v4 风格) ===")

  // 模拟: 两群细胞, RNA 与 ADT 各自编码群结构(坐标尺度不同)
  let rng = @src.Rng::from_int(7)
  let rna : Array[Array[Double]] = []
  let adt : Array[Array[Double]] = []
  // 群 A: RNA ~(0,0), ADT ~(1,5)
  for i = 0; i < 15; i = i + 1 {
    rna.push([rng.uniform(-0.3, 0.3), rng.uniform(-0.3, 0.3)])
    adt.push([1.0 + rng.uniform(-0.5, 0.5), 5.0 + rng.uniform(-0.5, 0.5)])
  }
  // 群 B: RNA ~(10,10), ADT ~(8,12)
  for i = 0; i < 15; i = i + 1 {
    rna.push([10.0 + rng.uniform(-0.3, 0.3), 10.0 + rng.uniform(-0.3, 0.3)])
    adt.push([8.0 + rng.uniform(-0.5, 0.5), 12.0 + rng.uniform(-0.5, 0.5)])
  }
  let m_rna = @src.wnn_modality("RNA", rna)
  let m_adt = @src.wnn_modality("ADT", adt)
  println("细胞 \{rna.length()} 个, RNA 与 ADT 双模态")

  // 1. WNN 整合
  println("\n1. WNN 整合(k=8)")
  let result = @src.wnn_integrate(m_rna, m_adt, k=8)
  println("  " + result.weight_summary())
  println("  权重守恒(每细胞和=1): \{result.weights_normalized()}")

  // 2. 细胞级权重分布
  println("\n2. 细胞级模态权重(前 6 个)")
  for i = 0; i < 6; i = i + 1 {
    println(
      "  cell\{i}: w[RNA]=\{result.modality_weights[i][0]} w[ADT]=\{result.modality_weights[i][1]}",
    )
  }

  // 3. 联合图
  println("\n3. 联合加权图(对称)")
  let g = result.joint_graph
  println("  顶点数: \{g.num_vertices()}")
  let mut n_edges = 0
  for i = 0; i < g.num_vertices(); i = i + 1 {
    n_edges = n_edges + g.adj[i].length()
  }
  println("  有向边数(无向 x2): \{n_edges}")
  println(
    "  cell0 邻居: \{g.adj[0].length()} 个, cell15 邻居: \{g.adj[15].length()} 个",
  )
  println("  cell0 度: \{g.degree(0)}")

  // 4. 联合图聚类
  println("\n4. 联合图聚类(scanpy_leiden 闭环)")
  let labels = result.cluster()
  let mut s0 = ""
  let mut s1 = ""
  for i = 0; i < 15; i = i + 1 {
    s0 = s0 + labels[i].to_string()
  }
  for i = 15; i < 30; i = i + 1 {
    s1 = s1 + labels[i].to_string()
  }
  println("  群 A 标签: \{s0}")
  println("  群 B 标签: \{s1}")
  println(
    "  两群分离: \{labels[0] != labels[15]} (群 A=\{labels[0]}, 群 B=\{labels[15]})",
  )

  // 5. 模态结构强弱对比
  println("\n5. 模态结构强弱对权重的影响")
  // RNA 强结构(紧簇), ADT 近随机 -> RNA 权重应高
  let rna2 : Array[Array[Double]] = []
  for i = 0; i < 10; i = i + 1 {
    rna2.push([0.0, 0.0])
  }
  for i = 0; i < 10; i = i + 1 {
    rna2.push([100.0, 100.0])
  }
  let adt2 : Array[Array[Double]] = []
  for i = 0; i < 20; i = i + 1 {
    adt2.push([i.to_double() * 0.1, i.to_double() % 2.0 * 0.05])
  }
  let a2 = @src.wnn_modality("RNA", rna2)
  let b2 = @src.wnn_modality("ADT", adt2)
  let r2 = @src.wnn_integrate(a2, b2, k=3)
  let mut sa = 0.0
  for w in r2.modality_weights {
    sa = sa + w[0]
  }
  println("  RNA 强结构时平均 w[RNA] = \{sa / 20.0} (应显著 > 0.5)")

  println(
    "\n完成: 双模态 KNN -> 跨模态一致性权重 -> 联合图 -> 聚类全链路本地完成.",
  )
}
