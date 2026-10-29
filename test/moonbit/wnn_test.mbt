// ===========================================================================
// wnn_test.mbt — 加权最近邻多模态整合测试
// ===========================================================================

///|
/// 构造双模态数据: 两群细胞, RNA 模态区分度高, ADT 模态噪声大.
fn wnn_test_data() -> (@src.WnnModality, @src.WnnModality) raise {
  // 两紧凑群(小抖动), 群间远离; 确定性抖动
  let rng = @src.Rng::from_int(7)
  let rna : Array[Array[Double]] = []
  for i = 0; i < 15; i = i + 1 {
    rna.push([rng.uniform(-0.3, 0.3), rng.uniform(-0.3, 0.3)])
  }
  for i = 0; i < 15; i = i + 1 {
    rna.push([10.0 + rng.uniform(-0.3, 0.3), 10.0 + rng.uniform(-0.3, 0.3)])
  }
  // ADT: 同样两群(坐标尺度不同)
  let adt : Array[Array[Double]] = []
  for i = 0; i < 15; i = i + 1 {
    adt.push([1.0 + rng.uniform(-0.5, 0.5), 5.0 + rng.uniform(-0.5, 0.5)])
  }
  for i = 0; i < 15; i = i + 1 {
    adt.push([8.0 + rng.uniform(-0.5, 0.5), 12.0 + rng.uniform(-0.5, 0.5)])
  }
  let m_rna = @src.wnn_modality("RNA", rna)
  let m_adt = @src.wnn_modality("ADT", adt)
  (m_rna, m_adt)
}

///|
test "wnn modality validation" {
  let mut raised = false
  try {
    let m = @src.wnn_modality("empty", [])
    ignore(m)
  } catch {
    _ => raised = true
  }
  assert_true(raised)
  raised = false
  try {
    let m = @src.wnn_modality("bad", [[1.0], [1.0, 2.0]])
    ignore(m)
  } catch {
    _ => raised = true
  }
  assert_true(raised)
}

///|
test "wnn integrate basic structure" {
  let (rna, adt) = wnn_test_data()
  let result = @src.wnn_integrate(rna, adt, k=5)
  assert_eq(result.n_cells, 30)
  assert_eq(result.modality_names.length(), 2)
  assert_eq(result.modality_names[0], "RNA")
  assert_eq(result.modality_names[1], "ADT")
  assert_eq(result.modality_weights.length(), 30)
  // 权重守恒
  assert_true(result.weights_normalized())
}

///|
test "wnn weights in valid range" {
  let (rna, adt) = wnn_test_data()
  let result = @src.wnn_integrate(rna, adt, k=5)
  for w in result.modality_weights {
    assert_true(w[0] >= 0.0 && w[0] <= 1.0)
    assert_true(w[1] >= 0.0 && w[1] <= 1.0)
  }
}

///|
test "wnn joint graph symmetric" {
  let (rna, adt) = wnn_test_data()
  let result = @src.wnn_integrate(rna, adt, k=5)
  // 对称性: adj[i] 含 (j,w) 则 adj[j] 含 (i,w)
  for i = 0; i < result.n_cells; i = i + 1 {
    for v in result.joint_graph.adj[i] {
      let j = v.0
      let mut found = false
      for u in result.joint_graph.adj[j] {
        if u.0 == i {
          found = true
          let d = if u.1 > v.1 { u.1 - v.1 } else { v.1 - u.1 }
          assert_true(d < 1.0e-9)
        }
      }
      assert_true(found)
    }
  }
}

///|
test "wnn joint graph preserves structure" {
  let (rna, adt) = wnn_test_data()
  let result = @src.wnn_integrate(rna, adt, k=5)
  // 同群细胞应有联合边
  let adj = result.joint_graph.adj
  let mut same_cluster_edges = 0
  for i = 0; i < 15; i = i + 1 {
    for v in adj[i] {
      if v.0 < 15 {
        same_cluster_edges = same_cluster_edges + 1
      }
    }
  }
  assert_true(same_cluster_edges > 20)
}

///|
test "wnn cluster recovers two populations" {
  let (rna, adt) = wnn_test_data()
  // k=8: 群内连通充分, leiden 干净分离两群
  let result = @src.wnn_integrate(rna, adt, k=8)
  let labels = result.cluster()
  assert_eq(labels.length(), 30)
  let l0 = labels[0]
  let l1 = labels[15]
  for i = 0; i < 15; i = i + 1 {
    assert_eq(labels[i], l0)
  }
  for i = 15; i < 30; i = i + 1 {
    assert_eq(labels[i], l1)
  }
  assert_true(l0 != l1)
}

///|
test "wnn identical modalities equal weights" {
  // 两模态完全相同 -> 一致性满 -> 权重接近 0.5/0.5
  let coords : Array[Array[Double]] = []
  for i = 0; i < 12; i = i + 1 {
    coords.push([i.to_double() * 0.7, i.to_double() % 3.0 * 1.1])
  }
  let a = @src.wnn_modality("A", coords)
  let b = @src.wnn_modality("B", coords)
  let result = @src.wnn_integrate(a, b, k=4)
  for w in result.modality_weights {
    assert_true(w[0] > 0.49 && w[0] < 0.51)
    assert_true(w[1] > 0.49 && w[1] < 0.51)
  }
}

///|
test "wnn cell count mismatch raises" {
  let a = @src.wnn_modality("A", [[1.0, 2.0], [3.0, 4.0]])
  let b = @src.wnn_modality("B", [[1.0, 2.0]])
  let mut raised = false
  try {
    let r = @src.wnn_integrate(a, b, k=1)
    ignore(r)
  } catch {
    _ => raised = true
  }
  assert_true(raised)
}

///|
test "wnn no smoothing option" {
  let (rna, adt) = wnn_test_data()
  let smoothed = @src.wnn_integrate(rna, adt, k=5, smoothing=true)
  let raw = @src.wnn_integrate(rna, adt, k=5, smoothing=false)
  // 两种模式都产出有效权重
  assert_true(smoothed.weights_normalized())
  assert_true(raw.weights_normalized())
}

///|
test "wnn weight summary" {
  let (rna, adt) = wnn_test_data()
  let result = @src.wnn_integrate(rna, adt, k=5)
  let s = result.weight_summary()
  assert_true(s.contains("wnn"))
  assert_true(s.contains("30 cells"))
  assert_true(s.contains("RNA"))
  assert_true(s.contains("ADT"))
}

///|
test "wnn modality-specific structure reflected in weights" {
  // RNA 模态强分离(簇内紧簇间远), ADT 模态近随机
  let rna : Array[Array[Double]] = []
  for i = 0; i < 10; i = i + 1 {
    rna.push([0.0, 0.0])
  }
  for i = 0; i < 10; i = i + 1 {
    rna.push([100.0, 100.0])
  }
  let adt : Array[Array[Double]] = []
  for i = 0; i < 20; i = i + 1 {
    adt.push([i.to_double() * 0.1, i.to_double() % 2.0 * 0.05])
  }
  let a = @src.wnn_modality("RNA", rna)
  let b = @src.wnn_modality("ADT", adt)
  let result = @src.wnn_integrate(a, b, k=3)
  // RNA 权重应显著高于 ADT(结构信号强)
  let mut sa = 0.0
  for w in result.modality_weights {
    sa = sa + w[0]
  }
  let mean_a = sa / 20.0
  assert_true(mean_a > 0.6)
}
