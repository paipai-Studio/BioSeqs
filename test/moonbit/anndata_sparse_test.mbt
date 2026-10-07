///|
/// Tests for AnnData sparse storage (P3 optimization).
///
/// Single-cell expression matrices are ~90% zeros.  Sparse storage reduces
/// memory from ~1.6 GB (10K × 20K dense) to ~50 MB (CSR/COO).

///|
/// Build a small sparse matrix: 3 cells × 4 genes, 5 non-zeros.
fn ads_test_sparse_matrix() -> @src.SparseArray {
  // [[0, 2, 0, 0],
  //  [0, 0, 3, 0],
  //  [1, 0, 0, 4]]
  let coords = [[0, 1], [1, 2], [2, 0], [2, 3], [0, 3]]
  let values = [2.0, 3.0, 1.0, 4.0, 5.0]
  @src.SparseArray::from_coo([3, 4], coords, values) catch {
    _ => abort("test matrix construction failed")
  }
}

///|
test "anndata_sparse_construct_shape" {
  let s = ads_test_sparse_matrix()
  let ad = @src.AnnData::new_sparse(s, ["c1", "c2", "c3"], ["g1", "g2", "g3", "g4"])
  assert_eq(ad.n_obs(), 3)
  assert_eq(ad.n_vars(), 4)
  let (n, m) = ad.shape()
  assert_eq(n, 3)
  assert_eq(m, 4)
  assert_true(ad.is_sparse())
}

///|
test "anndata_sparse_nnz_and_density" {
  let s = ads_test_sparse_matrix()
  let ad = @src.AnnData::new_sparse(s, ["c1", "c2", "c3"], ["g1", "g2", "g3", "g4"])
  // 5 non-zero entries in a 3×4 = 12 matrix
  assert_eq(ad.nnz(), 5)
  let d = ad.density()
  // 5/12 ≈ 0.4167
  assert_true(d > 0.41)
  assert_true(d < 0.42)
}

///|
test "anndata_sparse_x_matrix_materialize" {
  let s = ads_test_sparse_matrix()
  let ad = @src.AnnData::new_sparse(s, ["c1", "c2", "c3"], ["g1", "g2", "g3", "g4"])
  let dense = ad.x_matrix()
  assert_eq(dense.length(), 3)
  assert_eq(dense[0].length(), 4)
  // Check non-zero values
  assert_eq(dense[0][1], 2.0)
  assert_eq(dense[0][3], 5.0)
  assert_eq(dense[1][2], 3.0)
  assert_eq(dense[2][0], 1.0)
  assert_eq(dense[2][3], 4.0)
  // Check zeros
  assert_eq(dense[0][0], 0.0)
  assert_eq(dense[1][0], 0.0)
}

///|
test "anndata_sparse_get_x" {
  let s = ads_test_sparse_matrix()
  let ad = @src.AnnData::new_sparse(s, ["c1", "c2", "c3"], ["g1", "g2", "g3", "g4"])
  assert_eq(ad.get_x(0, 1), 2.0)
  assert_eq(ad.get_x(1, 2), 3.0)
  assert_eq(ad.get_x(2, 0), 1.0)
  assert_eq(ad.get_x(2, 3), 4.0)
  assert_eq(ad.get_x(0, 0), 0.0)
  assert_eq(ad.get_x(1, 0), 0.0)
}

///|
test "anndata_sparse_row_sums" {
  let s = ads_test_sparse_matrix()
  let ad = @src.AnnData::new_sparse(s, ["c1", "c2", "c3"], ["g1", "g2", "g3", "g4"])
  let rs = ad.row_sums()
  assert_eq(rs.length(), 3)
  // row 0: 2 + 5 = 7
  assert_eq(rs[0], 7.0)
  // row 1: 3
  assert_eq(rs[1], 3.0)
  // row 2: 1 + 4 = 5
  assert_eq(rs[2], 5.0)
}

///|
test "anndata_sparse_col_sums" {
  let s = ads_test_sparse_matrix()
  let ad = @src.AnnData::new_sparse(s, ["c1", "c2", "c3"], ["g1", "g2", "g3", "g4"])
  let cs = ad.col_sums()
  assert_eq(cs.length(), 4)
  // col 0: 1
  assert_eq(cs[0], 1.0)
  // col 1: 2
  assert_eq(cs[1], 2.0)
  // col 2: 3
  assert_eq(cs[2], 3.0)
  // col 3: 5 + 4 = 9
  assert_eq(cs[3], 9.0)
}

///|
test "anndata_sparse_row_means" {
  let s = ads_test_sparse_matrix()
  let ad = @src.AnnData::new_sparse(s, ["c1", "c2", "c3"], ["g1", "g2", "g3", "g4"])
  let rm = ad.row_means()
  // row 0: 7/4 = 1.75
  assert_eq(rm[0], 1.75)
  // row 1: 3/4 = 0.75
  assert_eq(rm[1], 0.75)
  // row 2: 5/4 = 1.25
  assert_eq(rm[2], 1.25)
}

///|
test "anndata_sparse_x_sum" {
  let s = ads_test_sparse_matrix()
  let ad = @src.AnnData::new_sparse(s, ["c1", "c2", "c3"], ["g1", "g2", "g3", "g4"])
  // 2 + 3 + 1 + 4 + 5 = 15
  assert_eq(ad.x_sum(), 15.0)
}

///|
test "anndata_dense_to_sparse_conversion" {
  let x = [
    [0.0, 2.0, 0.0, 0.0],
    [0.0, 0.0, 3.0, 0.0],
    [1.0, 0.0, 0.0, 4.0],
  ]
  let ad = @src.AnnData::new(x, ["c1", "c2", "c3"], ["g1", "g2", "g3", "g4"])
  assert_false(ad.is_sparse())
  let sp = ad.to_sparse()
  assert_true(sp.is_sparse())
  assert_eq(sp.nnz(), 4)
  assert_eq(sp.n_obs(), 3)
  assert_eq(sp.n_vars(), 4)
  // Values preserved
  assert_eq(sp.get_x(0, 1), 2.0)
  assert_eq(sp.get_x(1, 2), 3.0)
  assert_eq(sp.get_x(2, 0), 1.0)
  assert_eq(sp.get_x(2, 3), 4.0)
  assert_eq(sp.get_x(0, 0), 0.0)
}

///|
test "anndata_sparse_to_dense_conversion" {
  let s = ads_test_sparse_matrix()
  let ad = @src.AnnData::new_sparse(s, ["c1", "c2", "c3"], ["g1", "g2", "g3", "g4"])
  assert_true(ad.is_sparse())
  let dn = ad.to_dense()
  assert_false(dn.is_sparse())
  assert_eq(dn.n_obs(), 3)
  assert_eq(dn.n_vars(), 4)
  // Dense x is populated
  assert_eq(dn.x[0][1], 2.0)
  assert_eq(dn.x[1][2], 3.0)
  assert_eq(dn.x[2][0], 1.0)
  assert_eq(dn.x[2][3], 4.0)
}

///|
test "anndata_sparse_subset" {
  let s = ads_test_sparse_matrix()
  let ad = @src.AnnData::new_sparse(s, ["c1", "c2", "c3"], ["g1", "g2", "g3", "g4"])
  // Subset rows [2, 0] × cols [0, 3]
  let sub = ad.subset([2, 0], [0, 3])
  assert_true(sub.is_sparse())
  assert_eq(sub.n_obs(), 2)
  assert_eq(sub.n_vars(), 2)
  // Original [2,0]=1, [2,3]=4, [0,0]=0, [0,3]=5
  // Sub: row 0 = orig row 2, row 1 = orig row 0
  // col 0 = orig col 0, col 1 = orig col 3
  assert_eq(sub.get_x(0, 0), 1.0) // orig[2][0]
  assert_eq(sub.get_x(0, 1), 4.0) // orig[2][3]
  assert_eq(sub.get_x(1, 0), 0.0) // orig[0][0]
  assert_eq(sub.get_x(1, 1), 5.0) // orig[0][3]
}

///|
test "anndata_sparse_copy_preserves_sparse" {
  let s = ads_test_sparse_matrix()
  let ad = @src.AnnData::new_sparse(s, ["c1", "c2", "c3"], ["g1", "g2", "g3", "g4"])
    .set_obs("batch", ["A", "B", "A"])
  let cp = ad.copy()
  assert_true(cp.is_sparse())
  assert_eq(cp.n_obs(), 3)
  assert_eq(cp.nnz(), 5)
  // Metadata preserved
  let batch = cp.obs_column("batch")
  match batch {
    Some(arr) => assert_eq(arr[0], "A")
    None => assert_true(false)
  }
}

///|
test "anndata_dense_nnz" {
  let x = [
    [0.0, 2.0, 0.0],
    [0.0, 0.0, 3.0],
    [1.0, 0.0, 0.0],
  ]
  let ad = @src.AnnData::new(x, ["c1", "c2", "c3"], ["g1", "g2", "g3"])
  assert_eq(ad.nnz(), 3)
  assert_false(ad.is_sparse())
  let d = ad.density()
  // 3/9 ≈ 0.333
  assert_true(d > 0.33)
  assert_true(d < 0.34)
}

///|
test "anndata_dense_row_sums" {
  let x = [
    [1.0, 2.0, 3.0],
    [4.0, 5.0, 6.0],
  ]
  let ad = @src.AnnData::new(x, ["c1", "c2"], ["g1", "g2", "g3"])
  let rs = ad.row_sums()
  assert_eq(rs[0], 6.0)
  assert_eq(rs[1], 15.0)
  let cs = ad.col_sums()
  assert_eq(cs[0], 5.0)
  assert_eq(cs[1], 7.0)
  assert_eq(cs[2], 9.0)
}

///|
test "anndata_roundtrip_dense_sparse_dense" {
  let x = [
    [1.0, 0.0, 3.0, 0.0],
    [0.0, 5.0, 0.0, 7.0],
    [0.0, 0.0, 0.0, 0.0],
  ]
  let ad = @src.AnnData::new(x, ["c1", "c2", "c3"], ["g1", "g2", "g3", "g4"])
  let sp = ad.to_sparse()
  let dn = sp.to_dense()
  // Round-trip preserves values
  assert_eq(dn.x[0][0], 1.0)
  assert_eq(dn.x[0][2], 3.0)
  assert_eq(dn.x[1][1], 5.0)
  assert_eq(dn.x[1][3], 7.0)
  assert_eq(dn.x[2][0], 0.0)
  assert_eq(dn.x[2][3], 0.0)
  assert_eq(dn.nnz(), 4)
}

///|
test "anndata_sparse_empty_matrix" {
  // 2×3 all-zero matrix
  let s = @src.SparseArray::zeros([2, 3])
  let ad = @src.AnnData::new_sparse(s, ["c1", "c2"], ["g1", "g2", "g3"])
  assert_eq(ad.n_obs(), 2)
  assert_eq(ad.n_vars(), 3)
  assert_eq(ad.nnz(), 0)
  assert_eq(ad.density(), 0.0)
  assert_eq(ad.x_sum(), 0.0)
  let rs = ad.row_sums()
  assert_eq(rs[0], 0.0)
  assert_eq(rs[1], 0.0)
}
