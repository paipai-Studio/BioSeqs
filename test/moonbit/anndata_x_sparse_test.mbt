///|
/// Tests for first-class CSR/CSC access to AnnData X.

///|
fn xs_build_dense() -> @src.AnnData {
  @src.AnnData::new(
    [[1.0, 0.0, 2.0], [0.0, 0.0, 3.0], [4.0, 5.0, 0.0]],
    ["c1", "c2", "c3"],
    ["g1", "g2", "g3"],
  )
}

///|
test "x_sparse_dense_to_csr_csc" {
  let csr = @src.anndata_dense_to_csr([[1.0, 0.0, 2.0], [0.0, 3.0, 0.0]])
  assert_eq(@src.anndata_csr_nnz(csr), 3)
  assert_eq(csr.n_rows, 2)
  assert_eq(csr.n_cols, 3)
  let csc = @src.anndata_dense_to_csc([[1.0, 0.0], [0.0, 3.0]])
  let dense = @src.anndata_csc_to_dense(csc)
  assert_eq(dense[0][0], 1.0)
  assert_eq(dense[1][1], 3.0)
}

///|
test "x_sparse_new_csr_and_roundtrip" {
  let ad = xs_build_dense()
  let csr = @src.anndata_x_as_csr(ad)
  assert_eq(@src.anndata_csr_nnz(csr), 5)
  // construct straight from CSR — stays sparse
  let from_csr = match
    @src.AnnData::new_csr(csr, ["c1", "c2", "c3"], ["g1", "g2", "g3"]) {
    Some(a) => a
    None => fail("new_csr should succeed")
  }
  assert_eq(from_csr.x_layout(), "sparse")
  assert_eq(from_csr.n_obs(), 3)
  assert_eq(from_csr.n_vars(), 3)
  assert_eq(from_csr.get_x(0, 0), 1.0)
  assert_eq(from_csr.get_x(2, 1), 5.0)
  assert_eq(from_csr.get_x(1, 0), 0.0)
  // CSR request round-trips through the sparse storage
  let csr2 = @src.anndata_x_as_csr(from_csr)
  let dense2 = @src.anndata_csr_to_dense(csr2)
  assert_eq(dense2[0][2], 2.0)
  assert_eq(dense2[2][0], 4.0)
}

///|
test "x_sparse_new_csc_and_layout" {
  let ad = xs_build_dense()
  let csc = @src.anndata_x_as_csc(ad)
  let from_csc = match
    @src.AnnData::new_csc(csc, ["c1", "c2", "c3"], ["g1", "g2", "g3"]) {
    Some(a) => a
    None => fail("new_csc should succeed")
  }
  assert_eq(from_csc.x_layout(), "sparse")
  assert_eq(from_csc.get_x(1, 2), 3.0)
  // dense layout reported as dense
  assert_eq(ad.x_layout(), "dense")
}

///|
test "x_sparse_set_x_csr_inplace" {
  let ad = @src.AnnData::new([[1.0, 2.0], [3.0, 4.0]], ["c1", "c2"], [
    "g1", "g2",
  ])
  let csr = @src.anndata_dense_to_csr([[9.0, 0.0], [0.0, 7.0]])
  ignore(ad.set_x_csr(csr))
  assert_eq(ad.x_layout(), "sparse")
  assert_eq(ad.get_x(0, 0), 9.0)
  assert_eq(ad.get_x(0, 1), 0.0)
  assert_eq(ad.get_x(1, 1), 7.0)
  // malformed CSR (indptr shorter than n_rows + 1) keeps the original
  // matrix — detected inside anndata_csr_to_sparse
  let bad_dense = @src.anndata_dense_to_csr([[1.0]])
  let broken = @src.anndata_csr_make(
    bad_dense.data,
    bad_dense.indices,
    [0],
    1,
    1,
  )
  ignore(ad.set_x_csr(broken))
  assert_eq(ad.x_layout(), "sparse")
  assert_eq(ad.get_x(0, 0), 9.0)
}

///|
test "x_sparse_csr_slice_rows_and_cols" {
  let csr = @src.anndata_dense_to_csr([
    [1.0, 0.0, 2.0],
    [0.0, 3.0, 0.0],
    [4.0, 0.0, 5.0],
  ])
  // row slice keeps sparsity
  let rows = @src.anndata_csr_slice_rows(csr, [2, 0])
  assert_eq(rows.n_rows, 2)
  let rows_dense = @src.anndata_csr_to_dense(rows)
  assert_eq(rows_dense[0][0], 4.0)
  assert_eq(rows_dense[0][2], 5.0)
  assert_eq(rows_dense[1][0], 1.0)
  // col slice renumbers columns
  let cols = @src.anndata_csr_slice_cols(csr, [2, 0])
  assert_eq(cols.n_cols, 2)
  let cols_dense = @src.anndata_csr_to_dense(cols)
  assert_eq(cols_dense[0][0], 2.0) // old col 2 -> new col 0
  assert_eq(cols_dense[0][1], 1.0) // old col 0 -> new col 1
  // old col 1 was dropped, so row 1 has no values left
  assert_eq(cols_dense[1][0], 0.0)
  assert_eq(cols_dense[1][1], 0.0)
  // single row accessor
  let row1 = @src.anndata_csr_get_row(csr, 1)
  assert_eq(row1[1], 3.0)
  assert_eq(row1[0], 0.0)
}

///|
test "x_sparse_csc_to_csr" {
  let csc = @src.anndata_dense_to_csc([[1.0, 0.0], [2.0, 3.0]])
  let csr = @src.anndata_csc_to_csr(csc)
  let dense = @src.anndata_csr_to_dense(csr)
  assert_eq(dense[0][0], 1.0)
  assert_eq(dense[1][0], 2.0)
  assert_eq(dense[1][1], 3.0)
}
