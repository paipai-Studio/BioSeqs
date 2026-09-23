///|
/// Tests for the scanpy `read_10x_mtx` port (Matrix Market parsing + AnnData
/// assembly). Values are hand-computed from a tiny 10x fixture.

///|
test "scanpy_parse_mtx" {
  let mtx = "%%MatrixMarket matrix coordinate integer general\n% comment\n3 2 4\n1 1 10\n1 2 5\n2 1 3\n3 2 7\n"
  let (mat, n_rows, n_cols) = @src.scanpy_parse_mtx(mtx)
  assert_eq(n_rows, 3)
  assert_eq(n_cols, 2)
  assert_true(assert_close(mat[0][0], 10.0, 1.0e-9))
  assert_true(assert_close(mat[0][1], 5.0, 1.0e-9))
  assert_true(assert_close(mat[1][0], 3.0, 1.0e-9))
  assert_true(assert_close(mat[2][1], 7.0, 1.0e-9))
  assert_true(assert_close(mat[2][0], 0.0, 1.0e-12))
}

///|
test "scanpy_read_10x_mtx" {
  let mtx = "%%MatrixMarket matrix coordinate integer general\n3 2 4\n1 1 10\n1 2 5\n2 1 3\n3 2 7\n"
  let features = "g1\tsym1\tGene Expression\ng2\tsym2\tGene Expression\ng3\tsym3\tAntibody Capture\n"
  let barcodes = "cell1\ncell2\n"
  let adata = @src.scanpy_read_10x_mtx(mtx, features, barcodes, "gene_symbols", true, true)
  // 3 genes × 2 cells, transposed to 2×3; gex_only drops the Antibody gene
  assert_eq(adata.n_obs(), 2)
  assert_eq(adata.n_vars(), 2)
  assert_eq((adata.get_obs_names()[0], adata.get_obs_names()[1]), ("cell1", "cell2"))
  assert_eq((adata.get_var_names()[0], adata.get_var_names()[1]), ("sym1", "sym2"))
  // X = [[10, 3], [5, 0]]
  let x = adata.x_matrix()
  assert_true(assert_close(x[0][0], 10.0, 1.0e-9))
  assert_true(assert_close(x[0][1], 3.0, 1.0e-9))
  assert_true(assert_close(x[1][0], 5.0, 1.0e-9))
  assert_true(assert_close(x[1][1], 0.0, 1.0e-12))
  match adata.var_column("gene_ids") {
    Some(ids) => assert_eq((ids[0], ids[1]), ("g1", "g2"))
    None => assert_true(false)
  }
  match adata.var_column("feature_types") {
    Some(ft) => assert_eq((ft[0], ft[1]), ("Gene Expression", "Gene Expression"))
    None => assert_true(false)
  }
}

///|
test "scanpy_make_index_unique" {
  let u = @src.scanpy_make_index_unique(["a", "b", "a", "a", "c"])
  assert_eq((u[0], u[1], u[2], u[3], u[4]), ("a", "b", "a-1", "a-2", "c"))
}