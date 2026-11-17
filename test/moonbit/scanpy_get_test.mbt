///|
/// Tests for scanpy `get` layer: obs_df / var_df / rank_genes_groups_df.

///|
test "scanpy_get_obs_df" {
  let x = [[1.0, 2.0, 3.0], [4.0, 5.0, 6.0], [7.0, 8.0, 9.0]]
  let adata = @src.AnnData::new(x, ["c1", "c2", "c3"], ["G1", "G2", "G3"])
  let _ = adata.set_obs("ct", ["A", "A", "B"])
  // all genes
  let r1 = @src.scanpy_get_obs_df(adata)
  assert_eq(r1.df.n_cols(), 3)
  assert_eq(r1.df.n_rows(), 3)
  match r1.df.get_column("G2") {
    Some(v) => assert_eq(v[2], 8.0)
    None => assert_true(false)
  }
  // subset of genes
  let r2 = @src.scanpy_get_obs_df(adata, genes=["G3", "G1"])
  assert_eq(r2.df.n_cols(), 2)
  assert_eq(r2.df.col_names[0], "G3")
  match r2.df.get_row_by_name("c2") {
    Some(row) => {
      assert_eq(row[0], 6.0)
      assert_eq(row[1], 4.0)
    }
    None => assert_true(false)
  }
  // meta carries obs columns
  match r1.meta.get("ct") {
    Some(v) => assert_eq(v[1], "A")
    None => assert_true(false)
  }
  // missing gene -> zero column
  let r3 = @src.scanpy_get_obs_df(adata, genes=["NOPE"])
  match r3.df.get_column("NOPE") {
    Some(v) => {
      assert_eq(v.length(), 3)
      assert_eq(v[0], 0.0)
    }
    None => assert_true(false)
  }
}

///|
test "scanpy_get_var_df" {
  let x = [[1.0, 2.0], [3.0, 4.0], [5.0, 6.0]]
  let adata = @src.AnnData::new(x, ["c1", "c2", "c3"], ["G1", "G2"])
  let _ = adata.set_var("hvg", ["yes", "no"])
  let r = @src.scanpy_get_var_df(adata, cells=["c2", "c1"])
  assert_eq(r.df.n_cols(), 2)
  assert_eq(r.df.n_rows(), 2)
  assert_eq(r.df.col_names[0], "c2")
  match r.df.get_row_by_name("G1") {
    Some(row) => {
      assert_eq(row[0], 3.0)
      assert_eq(row[1], 1.0)
    }
    None => assert_true(false)
  }
  match r.meta.get("hvg") {
    Some(v) => assert_eq(v[0], "yes")
    None => assert_true(false)
  }
}

///|
test "scanpy_get_rank_genes_groups_df" {
  let x = [[10.0, 0.0], [12.0, 1.0], [1.0, 9.0], [0.0, 11.0]]
  let adata = @src.AnnData::new(x, ["c1", "c2", "c3", "c4"], ["GA", "GB"])
  let _ = adata.set_obs("ct", ["A", "A", "B", "B"])
  // rank GA vs GB with the t-test tool
  let res = @src.scanpy_rank_genes_groups(x, ["A", "A", "B", "B"], "B", "A")
  let payload = @src.scanpy_rank_genes_groups_json(
    ["A"],
    ["GA", "GB"],
    [res.map(fn(r) { r.t_stat })],
    [res.map(fn(r) { r.p_value })],
    [res.map(fn(r) { r.p_value })],
    [res.map(fn(r) { r.logfoldchange })],
  )
  let _ = adata.set_uns("rank_genes_groups", payload)
  match @src.scanpy_get_rank_genes_groups_df(adata, group="A") {
    Some(df) => {
      assert_eq(df.df.n_rows(), 2)
      assert_eq(df.df.n_cols(), 4)
      assert_eq(df.df.row_names[0], "GA")
      match df.df.get_row_by_name("GA") {
        Some(row) => assert_true(row[0] > 0.0) // GA enriched in A
        None => assert_true(false)
      }
    }
    None => assert_true(false)
  }
  // logFC cutoff keeps only GA (lfc=3 for GA, |lfc|≈2.87 for GB)
  match
    @src.scanpy_get_rank_genes_groups_df(adata, group="A", logfc_cutoff=3.0) {
    Some(df) => {
      assert_eq(df.df.n_rows(), 1)
      assert_eq(df.df.row_names[0], "GA")
    }
    None => assert_true(false)
  }
  // unknown group -> None
  assert_true(
    match @src.scanpy_get_rank_genes_groups_df(adata, group="ZZ") {
      Some(_) => false
      None => true
    },
  )
}
