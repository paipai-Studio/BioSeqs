// scanpy_scib_metrics_test.mbt — scIB 整合评估指标补全（ARI/NMI/iLISI/kBET/
//                                 batch PCA entropy）测试
// ===========================================================================

///|
test "scanpy_scib_ari_nmi_agreement" {
  // 完全一致 → ARI / NMI = 1
  let a = [0, 0, 0, 1, 1, 1]
  assert_true((@src.scanpy_scib_ari(a, a) - 1.0).abs() < 1.0e-12)
  assert_true((@src.scanpy_scib_nmi(a, a) - 1.0).abs() < 1.0e-12)
  // 重命名类别不改分
  assert_true(
    (@src.scanpy_scib_ari(a, [2, 2, 2, 7, 7, 7]) - 1.0).abs() < 1.0e-12,
  )
  // 交叉错配 → 远低于 1
  let mixed_pred = [0, 1, 0, 1, 0, 1]
  let ari_cross = @src.scanpy_scib_ari(a, mixed_pred)
  assert_true(ari_cross < 0.1)
  assert_true(@src.scanpy_scib_nmi(a, mixed_pred) < 0.7)
  // 恒定标签对之间：信息一致（都无信息）→ NMI = 1；ARI 退化取 1（sklearn 约定）
  let flat = [3, 3, 3, 3, 3, 3]
  assert_true(
    (@src.scanpy_scib_nmi(flat, [9, 9, 9, 9, 9, 9]) - 1.0).abs() < 1.0e-12,
  )
  // 长度不匹配 → 0
  assert_eq(@src.scanpy_scib_ari(a, [0, 1]), 0.0)
  assert_eq(@src.scanpy_scib_nmi(a, [0, 1]), 0.0)
}

///|
test "scanpy_scib_ilisi_kbet_mixing" {
  // 沿线交替批次 → 邻域充分混合，iLISI 接近 2（两批各半），kBET 高
  let mixed = [
    [0.0, 0.0],
    [0.1, 0.0],
    [0.2, 0.0],
    [0.3, 0.0],
    [0.4, 0.0],
    [0.5, 0.0],
  ]
  let mixed_batch = [0, 1, 0, 1, 0, 1]
  let ilisi_mix = @src.scanpy_scib_ilisi(mixed, mixed_batch, k=5)
  let kbet_mix = @src.scanpy_scib_kbet(mixed, mixed_batch, k=5)
  assert_true(ilisi_mix > 1.9 && ilisi_mix < 2.1)
  assert_true(kbet_mix > 0.8)
  // 空间分离的批次 → 邻域批纯净，iLISI 接近 1，kBET 崩塌
  let split_batch = [0, 0, 0, 1, 1, 1]
  let split = [
    [0.0, 0.0],
    [0.1, 0.0],
    [0.2, 0.0],
    [5.0, 0.0],
    [5.1, 0.0],
    [5.2, 0.0],
  ]
  let ilisi_split = @src.scanpy_scib_ilisi(split, split_batch, k=3)
  let kbet_split = @src.scanpy_scib_kbet(split, split_batch, k=3)
  assert_true(ilisi_split < 1.5)
  assert_true(kbet_split < kbet_mix)
  // 同一函数可用于 label mixing（分组粒度泛化）
  let labels = [0, 1, 0, 1, 0, 1]
  assert_true(@src.scanpy_scib_ilisi(mixed, labels, k=5) > 1.9)
}

///|
test "scanpy_scib_pca_entropy_batch_structure" {
  // 主方向上交替批次：每个 PC 分箱内批次分布均匀 → 熵接近最大（1）
  let mixed = [
    [0.0, 0.0],
    [0.1, 0.0],
    [0.2, 0.0],
    [0.3, 0.0],
    [0.4, 0.0],
    [0.5, 0.0],
    [0.6, 0.0],
    [0.7, 0.0],
  ]
  let mixed_batch = [0, 1, 0, 1, 0, 1, 0, 1]
  let h_mix = @src.scanpy_scib_pca_entropy(mixed, mixed_batch)
  assert_true(h_mix > 0.9)
  // 批次沿主 PC 完全分块：低分箱全是批 0、高分箱全是批 1 → 熵崩塌
  let split = [
    [0.0, 0.0],
    [0.1, 0.0],
    [0.2, 0.0],
    [0.3, 0.0],
    [5.0, 0.0],
    [5.1, 0.0],
    [5.2, 0.0],
    [5.3, 0.0],
  ]
  let split_batch = [0, 0, 0, 0, 1, 1, 1, 1]
  let h_split = @src.scanpy_scib_pca_entropy(split, split_batch)
  assert_true(h_split < 0.6)
  assert_true(h_split < h_mix)
  // 单批次平凡满分；空输入为 0
  assert_eq(@src.scanpy_scib_pca_entropy(mixed, [7, 7, 7, 7, 7, 7, 7, 7]), 1.0)
  assert_eq(@src.scanpy_scib_pca_entropy([], []), 0.0)
}

///|
test "scanpy_scib_metrics_full_bundle" {
  // 混合批次 + 与标签一致的簇划分 → 全家桶字段齐备且方向正确
  // 8 个点配 4 个默认分箱：每箱 2 点、批次交替，pca_entropy 接近 1
  let mixed = [
    [0.0, 0.0],
    [0.1, 0.0],
    [0.2, 0.0],
    [0.3, 0.0],
    [0.4, 0.0],
    [0.5, 0.0],
    [0.6, 0.0],
    [0.7, 0.0],
  ]
  let mixed_batch = [0, 1, 0, 1, 0, 1, 0, 1]
  let labels = [1, 1, 0, 0, 1, 1, 0, 0]
  let clusters = [1, 1, 0, 0, 1, 1, 0, 0]
  let m = @src.scanpy_scib_metrics(mixed, mixed_batch, labels, clusters, k=5)
  assert_true((m.ari - 1.0).abs() < 1.0e-12)
  assert_true((m.nmi - 1.0).abs() < 1.0e-12)
  assert_true(m.pca_entropy > 0.9)
  // clusters 缺省（空数组）→ ARI/NMI 记 0
  let m0 = @src.scanpy_scib_metrics(mixed, mixed_batch, labels, [], k=5)
  assert_eq(m0.ari, 0.0)
  assert_eq(m0.nmi, 0.0)
  // 噪声簇划分 → ARI 掉到接近随机
  let m_noise = @src.scanpy_scib_metrics(
    mixed,
    mixed_batch,
    labels,
    [0, 1, 0, 1, 0, 1, 0, 1],
    k=5,
  )
  assert_true(m_noise.ari < 0.1)
}
