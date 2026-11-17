// scvi_arches_test.mbt — scArches 迁移学习测试（含 latent kNN 标签迁移）

///|
test "arches knn label transfer separates clusters" {
  let ref_latent = [
    [0.0, 0.0],
    [0.1, 0.0],
    [0.0, 0.1],
    [0.1, 0.1],
    [5.0, 5.0],
    [5.1, 5.0],
    [5.0, 5.1],
  ]
  let ref_labels : Array[Int32] = [0, 0, 0, 0, 1, 1, 1]
  let query = [[0.05, 0.05], [5.0, 5.0], [5.05, 5.0]]
  let res = @src.arches_knn_label_transfer(ref_latent, ref_labels, query, k=3)
  assert_eq(res.length(), 3)
  // 两簇内的查询细胞均转移到正确标签且高置信
  assert_eq(res[0].label, 0)
  assert_true(res[0].confidence > 0.999)
  assert_eq(res[0].neighbors.length(), 3)
  assert_eq(res[0].neighbors[0], 0)
  assert_eq(res[1].label, 1)
  assert_eq(res[2].label, 1)
  // 近邻按距离升序（[5,5] 自身距离 0）
  assert_eq(res[1].neighbors[0], 4)
}

///|
test "arches knn transfer tie and weighting" {
  // k=2 两邻各占一票 → 平局取较小标签，置信度 0.5
  let tie = @src.arches_knn_label_transfer(
    [[0.0, 0.0], [1.0, 1.0]],
    [0, 1],
    [[0.5, 0.5]],
    k=2,
  )
  assert_eq(tie[0].label, 0)
  assert_true(@moonbit.assert_num_close(tie[0].confidence, 0.5, 1.0e-9))
  // 均匀投票：一个近邻 label 1 + 一个远邻 label 0 → 平局取较小标签
  let ref_latent = [[0.0], [0.1], [2.0]]
  let ref_labels : Array[Int32] = [0, 0, 1]
  let uniform = @src.arches_knn_label_transfer(
    ref_latent,
    ref_labels,
    [[1.8]],
    k=2,
  )
  assert_eq(uniform[0].label, 0)
  assert_true(@moonbit.assert_num_close(uniform[0].confidence, 0.5, 1.0e-9))
  // 距离加权后近距离近邻主导：同一点翻转为 label 1
  let weighted = @src.arches_knn_label_transfer(
    ref_latent,
    ref_labels,
    [[1.8]],
    k=2,
    distance_weighted=true,
  )
  assert_eq(weighted[0].label, 1)
  assert_true(weighted[0].confidence > 0.8)
}

///|
test "arches knn transfer edge cases" {
  let ref_latent = [[0.0], [1.0]]
  let ref_labels : Array[Int32] = [0, 1]
  // k 大于参考集大小 → 收缩到 n_ref；两邻平局取较小标签
  let res = @src.arches_knn_label_transfer(
    ref_latent,
    ref_labels,
    [[0.9]],
    k=10,
  )
  assert_eq(res.length(), 1)
  assert_eq(res[0].neighbors.length(), 2)
  assert_eq(res[0].label, 0)
  assert_true(@moonbit.assert_num_close(res[0].confidence, 0.5, 1.0e-9))
  // k=1 → 最近邻标签，置信度 1
  let k1 = @src.arches_knn_label_transfer(ref_latent, ref_labels, [[0.9]], k=1)
  assert_eq(k1[0].label, 1)
  assert_true(k1[0].confidence > 0.999)
  // 标签数量与参考集不一致 → 空
  assert_eq(
    @src.arches_knn_label_transfer(ref_latent, [0], [[0.5]]).length(),
    0,
  )
  // 空参考 / 空查询 → 空
  assert_eq(@src.arches_knn_label_transfer([], [], [[0.5]]).length(), 0)
  assert_eq(
    @src.arches_knn_label_transfer(ref_latent, ref_labels, []).length(),
    0,
  )
}

///|
test "arches prepare_query_and_transfer_accuracy" {
  // 基因对齐：查询缺 G3，多 G4 → 按参考顺序对齐，缺失补 0
  let aligned = @src.arches_prepare_query(
    [[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]],
    ["G1", "G2", "G4"],
    ["G1", "G2", "G3"],
  )
  assert_eq(aligned.length(), 2)
  assert_true(@moonbit.assert_num_close(aligned[0][0], 1.0, 1.0e-12))
  assert_true(@moonbit.assert_num_close(aligned[0][2], 0.0, 1.0e-12))
  assert_true(@moonbit.assert_num_close(aligned[1][2], 0.0, 1.0e-12))
  // 迁移准确率
  assert_true(
    @moonbit.assert_num_close(
      @src.arches_transfer_accuracy([0, 1, 1], [0, 1, 0]),
      2.0 / 3.0,
      1.0e-12,
    ),
  )
  assert_true(@src.arches_transfer_accuracy([], []) == 0.0)
}
