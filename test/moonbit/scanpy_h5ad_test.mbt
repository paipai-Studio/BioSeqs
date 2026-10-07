///|
/// Tests for h5ad text-format IO (JSON-based AnnData serialization).

///|
test "h5ad_dense_roundtrip" {
  let x = [
    [1.0, 0.0, 3.0],
    [0.0, 5.0, 0.0],
    [7.0, 0.0, 9.0],
  ]
  let ad = @src.AnnData::new(x, ["c1", "c2", "c3"], ["g1", "g2", "g3"])
  let json = @src.scanpy_write_h5ad(ad)
  let ad2 = match @src.scanpy_read_h5ad(json) {
    Some(a) => a
    None => abort("h5ad parse failed")
  }
  assert_eq(ad2.n_obs(), 3)
  assert_eq(ad2.n_vars(), 3)
  // Check X values
  assert_eq(ad2.get_x(0, 0), 1.0)
  assert_eq(ad2.get_x(0, 2), 3.0)
  assert_eq(ad2.get_x(1, 1), 5.0)
  assert_eq(ad2.get_x(2, 0), 7.0)
  assert_eq(ad2.get_x(2, 2), 9.0)
  assert_eq(ad2.get_x(0, 1), 0.0)
  // Check names
  assert_eq(ad2.get_obs_names()[0], "c1")
  assert_eq(ad2.get_obs_names()[2], "c3")
  assert_eq(ad2.get_var_names()[0], "g1")
  assert_eq(ad2.get_var_names()[2], "g3")
}

///|
test "h5ad_dense_with_metadata_roundtrip" {
  let x = [[1.0, 2.0], [3.0, 4.0], [5.0, 6.0]]
  let ad = @src.AnnData::new(x, ["c1", "c2", "c3"], ["g1", "g2"])
    .set_obs("batch", ["A", "B", "A"])
    .set_var("gene_type", ["protein", "lncRNA"])
    .set_obsm("X_pca", [[1.0, 2.0], [3.0, 4.0], [5.0, 6.0]])
  let json = @src.scanpy_write_h5ad(ad)
  let ad2 = match @src.scanpy_read_h5ad(json) {
    Some(a) => a
    None => abort("h5ad parse failed")
  }
  // Check obs metadata
  let batch = ad2.obs_column("batch")
  match batch {
    Some(arr) => {
      assert_eq(arr[0], "A")
      assert_eq(arr[1], "B")
      assert_eq(arr[2], "A")
    }
    None => assert_true(false)
  }
  // Check var metadata
  let gt = ad2.var_column("gene_type")
  match gt {
    Some(arr) => {
      assert_eq(arr[0], "protein")
      assert_eq(arr[1], "lncRNA")
    }
    None => assert_true(false)
  }
  // Check obsm
  let pca = ad2.obsm["X_pca"]
  assert_eq(pca[0][0], 1.0)
  assert_eq(pca[1][1], 4.0)
  assert_eq(pca[2][0], 5.0)
}

///|
test "h5ad_sparse_roundtrip" {
  // Build a sparse AnnData
  let coords = [[0, 1], [1, 0], [2, 2], [0, 2]]
  let values = [2.0, 3.0, 5.0, 7.0]
  let s = @src.SparseArray::from_coo([3, 3], coords, values) catch {
    _ => abort("sparse construction failed")
  }
  let ad = @src.AnnData::new_sparse(s, ["c1", "c2", "c3"], ["g1", "g2", "g3"])
  let json = @src.scanpy_write_h5ad(ad)
  let ad2 = match @src.scanpy_read_h5ad(json) {
    Some(a) => a
    None => abort("h5ad parse failed")
  }
  assert_true(ad2.is_sparse())
  assert_eq(ad2.n_obs(), 3)
  assert_eq(ad2.n_vars(), 3)
  assert_eq(ad2.nnz(), 4)
  // Check values
  assert_eq(ad2.get_x(0, 1), 2.0)
  assert_eq(ad2.get_x(1, 0), 3.0)
  assert_eq(ad2.get_x(2, 2), 5.0)
  assert_eq(ad2.get_x(0, 2), 7.0)
  assert_eq(ad2.get_x(0, 0), 0.0)
  assert_eq(ad2.get_x(1, 1), 0.0)
}

///|
test "h5ad_sparse_with_metadata_roundtrip" {
  let coords = [[0, 0], [1, 1]]
  let values = [10.0, 20.0]
  let s = @src.SparseArray::from_coo([2, 2], coords, values) catch {
    _ => abort("sparse construction failed")
  }
  let ad = @src.AnnData::new_sparse(s, ["cell1", "cell2"], ["gene1", "gene2"])
    .set_obs("tissue", ["brain", "liver"])
    .set_var("mt", ["no", "yes"])
  let json = @src.scanpy_write_h5ad(ad)
  let ad2 = match @src.scanpy_read_h5ad(json) {
    Some(a) => a
    None => abort("h5ad parse failed")
  }
  assert_true(ad2.is_sparse())
  // Check obs
  let tissue = ad2.obs_column("tissue")
  match tissue {
    Some(arr) => {
      assert_eq(arr[0], "brain")
      assert_eq(arr[1], "liver")
    }
    None => assert_true(false)
  }
  // Check var
  let mt = ad2.var_column("mt")
  match mt {
    Some(arr) => {
      assert_eq(arr[0], "no")
      assert_eq(arr[1], "yes")
    }
    None => assert_true(false)
  }
  // Check sparse values
  assert_eq(ad2.get_x(0, 0), 10.0)
  assert_eq(ad2.get_x(1, 1), 20.0)
}

///|
test "h5ad_empty_matrix" {
  let ad = @src.AnnData::new([], [], [])
  let json = @src.scanpy_write_h5ad(ad)
  let ad2 = match @src.scanpy_read_h5ad(json) {
    Some(a) => a
    None => abort("h5ad parse failed")
  }
  assert_eq(ad2.n_obs(), 0)
  assert_eq(ad2.n_vars(), 0)
}

///|
test "h5ad_with_layers_roundtrip" {
  let x = [[1.0, 2.0], [3.0, 4.0]]
  let ad = @src.AnnData::new(x, ["c1", "c2"], ["g1", "g2"])
    .set_layer("counts", [[10.0, 20.0], [30.0, 40.0]])
    .set_layer("logcounts", [[1.0, 2.0], [3.0, 4.0]])
  let json = @src.scanpy_write_h5ad(ad)
  let ad2 = match @src.scanpy_read_h5ad(json) {
    Some(a) => a
    None => abort("h5ad parse failed")
  }
  let counts = ad2.layers["counts"]
  assert_eq(counts[0][0], 10.0)
  assert_eq(counts[1][1], 40.0)
  let logcounts = ad2.layers["logcounts"]
  assert_eq(logcounts[0][1], 2.0)
  assert_eq(logcounts[1][0], 3.0)
}

///|
test "h5ad_with_obsp_roundtrip" {
  let x = [[1.0, 2.0], [3.0, 4.0]]
  let ad = @src.AnnData::new(x, ["c1", "c2"], ["g1", "g2"])
    .set_obsp("connectivities", [[1.0, 0.5], [0.5, 1.0]])
  let json = @src.scanpy_write_h5ad(ad)
  let ad2 = match @src.scanpy_read_h5ad(json) {
    Some(a) => a
    None => abort("h5ad parse failed")
  }
  let conn = ad2.obsp["connectivities"]
  assert_eq(conn[0][0], 1.0)
  assert_eq(conn[0][1], 0.5)
  assert_eq(conn[1][0], 0.5)
  assert_eq(conn[1][1], 1.0)
}

///|
test "h5ad_double_roundtrip_stable" {
  // Write → read → write → read, second output should match first
  let x = [[1.0, 0.0, 3.0], [0.0, 5.0, 0.0]]
  let ad = @src.AnnData::new(x, ["c1", "c2"], ["g1", "g2", "g3"])
    .set_obs("batch", ["A", "B"])
  let json1 = @src.scanpy_write_h5ad(ad)
  let ad1 = match @src.scanpy_read_h5ad(json1) {
    Some(a) => a
    None => abort("first parse failed")
  }
  let json2 = @src.scanpy_write_h5ad(ad1)
  let ad2 = match @src.scanpy_read_h5ad(json2) {
    Some(a) => a
    None => abort("second parse failed")
  }
  // Values stable across double roundtrip
  assert_eq(ad2.n_obs(), 2)
  assert_eq(ad2.n_vars(), 3)
  assert_eq(ad2.get_x(0, 0), 1.0)
  assert_eq(ad2.get_x(0, 2), 3.0)
  assert_eq(ad2.get_x(1, 1), 5.0)
  let batch = ad2.obs_column("batch")
  match batch {
    Some(arr) => {
      assert_eq(arr[0], "A")
      assert_eq(arr[1], "B")
    }
    None => assert_true(false)
  }
}
