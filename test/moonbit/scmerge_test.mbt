///|
test "scmerge_find_seg_basic" {
  // Create two small batches with 3 cells each, 5 genes
  let batch1 = [
    [1.0, 2.0, 3.0, 4.0, 5.0],
    [1.1, 2.1, 3.1, 4.1, 5.1],
    [0.9, 1.9, 2.9, 3.9, 4.9],
  ]
  let batch2 = [
    [2.0, 3.0, 4.0, 5.0, 6.0],
    [2.1, 3.1, 4.1, 5.1, 6.1],
    [1.9, 2.9, 3.9, 4.9, 5.9],
  ]
  let batches = [batch1, batch2]
  let seg = @src.scm_find_seg(batches, n_seg=2)
  assert_eq(seg.length(), 2)
  // All indices should be valid (0-4)
  let mut i = 0
  while i < seg.length() {
    assert_true(seg[i] >= 0 && seg[i] < 5)
    i = i + 1
  }
}

///|
test "scmerge_find_seg_empty" {
  let seg = @src.scm_find_seg([], n_seg=10)
  assert_eq(seg.length(), 0)
}

///|
test "scmerge_find_seg_all_stable" {
  // All genes have the same variance -> all are equally stable
  let batch1 = [[1.0, 2.0, 3.0]]
  let batch2 = [[4.0, 5.0, 6.0]]
  let batches = [batch1, batch2]
  // With only 1 gene, should return it
  let seg = @src.scm_find_seg(batches, n_seg=5)
  // With 3 genes and n_seg=5, should return all 3 genes
  assert_eq(seg.length(), 3)
}

///|
test "scmerge_zscore_matrix" {
  let data = [[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]]
  let z = @src.scm_zscore_matrix(data)
  assert_eq(z.length(), 2)
  assert_eq(z[0].length(), 3)
  // Each row should be approximately z-scored
  let mut i = 0
  while i < z.length() {
    let mut sum = 0.0
    let mut j = 0
    while j < z[i].length() {
      sum = sum + z[i][j]
      j = j + 1
    }
    // Mean should be ~0
    assert_true(sum.abs() < 1.0e-10)
    i = i + 1
  }
}

///|
test "scmerge_ruv4_basic" {
  // Create simple data: 3 genes, 6 cells (3 per batch)
  let data = [
    [1.0, 1.1, 0.9, 5.0, 5.1, 4.9],
    [2.0, 2.1, 1.9, 6.0, 6.1, 5.9],
    [3.0, 3.1, 2.9, 7.0, 7.1, 6.9],
  ]
  let batch_idx = [0, 0, 0, 1, 1, 1]
  let seg_genes = [0] // Gene 0 is the most stable
  // This is a structural test - just verify it runs without error
  let corrected = @src.scm_ruv4(data, batch_idx, seg_genes, 1)
  assert_eq(corrected.length(), 3)
  assert_eq(corrected[0].length(), 6)
}

///|
test "scmerge_ruv4_empty" {
  let corrected = @src.scm_ruv4([], [], [], 1)
  assert_eq(corrected.length(), 0)
}

///|
test "scmerge_invert_matrix" {
  // Identity matrix
  let identity = [[1.0, 0.0], [0.0, 1.0]]
  let inv = @src.scm_invert_matrix(identity, 2)
  // Inverse of identity is identity
  assert_true((inv[0][0] - 1.0).abs() < 1.0e-10)
  assert_true((inv[0][1] - 0.0).abs() < 1.0e-10)
  assert_true((inv[1][0] - 0.0).abs() < 1.0e-10)
  assert_true((inv[1][1] - 1.0).abs() < 1.0e-10)
}

///|
test "scmerge_invert_matrix_2x2" {
  // Simple 2x2 matrix
  let mat = [[4.0, 7.0], [2.0, 6.0]]
  let inv = @src.scm_invert_matrix(mat, 2)
  // A * A^{-1} should be identity
  let mut i = 0
  while i < 2 {
    let mut j = 0
    while j < 2 {
      let mut sum = 0.0
      let mut k = 0
      while k < 2 {
        sum = sum + mat[i][k] * inv[k][j]
        k = k + 1
      }
      let expected = if i == j { 1.0 } else { 0.0 }
      assert_true((sum - expected).abs() < 1.0e-10)
      j = j + 1
    }
    i = i + 1
  }
}

///|
test "scmerge_scmerge_basic" {
  // Create two small batches with 3 cells each, 4 genes
  let batch1 = [
    [1.0, 2.0, 3.0, 4.0],
    [1.1, 2.1, 3.1, 4.1],
    [0.9, 1.9, 2.9, 3.9],
  ]
  let batch2 = [
    [5.0, 6.0, 7.0, 8.0],
    [5.1, 6.1, 7.1, 8.1],
    [4.9, 5.9, 6.9, 7.9],
  ]
  let batches = [batch1, batch2]
  let result = @src.scmerge(batches, n_seg=2, n_factors=1)
  assert_eq(result.corrected_data.length(), 6)
  assert_eq(result.corrected_data[0].length(), 4)
  assert_eq(result.seg_genes.length(), 2)
}

///|
test "scmerge_scmerge_empty" {
  let result = @src.scmerge([], n_seg=10, n_factors=2)
  assert_eq(result.corrected_data.length(), 0)
}

///|
test "scmerge_scmerge_single_batch" {
  let batch = [[1.0, 2.0], [3.0, 4.0]]
  let result = @src.scmerge([batch], n_seg=1, n_factors=1)
  assert_eq(result.corrected_data.length(), 2)
  assert_eq(result.corrected_data[0].length(), 2)
}

///|
test "scmerge_scmerge_result_new" {
  let result = @src.ScMergeResult::new(
    [[1.0, 2.0], [3.0, 4.0]],
    [0, 1],
    [[0.0], [1.0]],
    [0, 0],
  )
  assert_eq(result.corrected_data.length(), 2)
  assert_eq(result.seg_genes.length(), 2)
  assert_eq(result.estimated_batch_factors.length(), 2)
  assert_eq(result.cell_type_labels.length(), 2)
}

///|
test "scmerge_randomized_svd_basic" {
  let matrix = [[1.0, 2.0, 3.0], [4.0, 5.0, 6.0], [7.0, 8.0, 9.0]]
  let v = @src.scm_randomized_svd(matrix, 2)
  // Result should be (n_cols x k) = (3 x 2)
  assert_eq(v.length(), 3)
  assert_eq(v[0].length(), 2)
}

///|
test "scmerge_randomized_svd_empty" {
  let v = @src.scm_randomized_svd([], 2)
  assert_eq(v.length(), 0)
}
