///|
/// Tests for AnnData setter/getter completeness (set_varm, set_varp, set_uns)
/// and binary h5ad format round-trip.

///|
test "anndata_set_varm" {
  let x = [[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]]
  let ad = @src.AnnData::new(x, ["c1", "c2"], ["g1", "g2", "g3"])
    .set_varm("VF", [[1.0, 0.0], [0.0, 1.0], [0.5, 0.5]])
  let vf = ad.get_varm("VF")
  match vf {
    Some(m) => {
      assert_eq(m.length(), 3)
      assert_eq(m[0][0], 1.0)
      assert_eq(m[1][1], 1.0)
      assert_eq(m[2][0], 0.5)
    }
    None => assert_true(false)
  }
  // Missing key
  let missing = ad.get_varm("nonexistent")
  match missing {
    Some(_) => assert_true(false)
    None => assert_true(true)
  }
}

///|
test "anndata_set_varp" {
  let x = [[1.0, 2.0], [3.0, 4.0]]
  let ad = @src.AnnData::new(x, ["c1", "c2"], ["g1", "g2"])
    .set_varp("correlation", [[1.0, 0.8], [0.8, 1.0]])
  let corr = ad.varp["correlation"]
  assert_eq(corr[0][0], 1.0)
  assert_eq(corr[0][1], 0.8)
  assert_eq(corr[1][0], 0.8)
  assert_eq(corr[1][1], 1.0)
}

///|
test "anndata_set_uns" {
  let x = [[1.0, 2.0], [3.0, 4.0]]
  let ad = @src.AnnData::new(x, ["c1", "c2"], ["g1", "g2"])
    .set_uns("title", "Sample experiment")
    .set_uns("date", "2024-01-15")
    .set_uns("n_clusters", "3")
  let title = ad.get_uns("title")
  match title {
    Some(s) => assert_eq(s, "Sample experiment")
    None => assert_true(false)
  }
  let date = ad.get_uns("date")
  match date {
    Some(s) => assert_eq(s, "2024-01-15")
    None => assert_true(false)
  }
  let missing = ad.get_uns("nonexistent")
  match missing {
    Some(_) => assert_true(false)
    None => assert_true(true)
  }
}

///|
test "anndata_get_layer_and_obsm" {
  let x = [[1.0, 2.0], [3.0, 4.0]]
  let ad = @src.AnnData::new(x, ["c1", "c2"], ["g1", "g2"])
    .set_layer("counts", [[10.0, 20.0], [30.0, 40.0]])
    .set_obsm("X_pca", [[1.0, 2.0], [3.0, 4.0]])
  let layer = ad.get_layer("counts")
  match layer {
    Some(m) => assert_eq(m[0][0], 10.0)
    None => assert_true(false)
  }
  let obsm = ad.get_obsm("X_pca")
  match obsm {
    Some(m) => assert_eq(m[1][1], 4.0)
    None => assert_true(false)
  }
  // Missing keys
  match ad.get_layer("nonexistent") {
    Some(_) => assert_true(false)
    None => assert_true(true)
  }
  match ad.get_obsm("nonexistent") {
    Some(_) => assert_true(false)
    None => assert_true(true)
  }
}

///|
test "h5ad_binary_dense_roundtrip" {
  let x = [
    [1.0, 0.0, 3.0],
    [0.0, 5.0, 0.0],
    [7.0, 0.0, 9.0],
  ]
  let ad = @src.AnnData::new(x, ["c1", "c2", "c3"], ["g1", "g2", "g3"])
    .set_obs("batch", ["A", "B", "A"])
    .set_var("type", ["G", "G", "G"])
  let bytes = @src.scanpy_write_h5ad_binary(ad)
  let ad2 = match @src.scanpy_read_h5ad_binary(bytes) {
    Some(a) => a
    None => abort("binary h5ad parse failed")
  }
  assert_eq(ad2.n_obs(), 3)
  assert_eq(ad2.n_vars(), 3)
  assert_eq(ad2.get_x(0, 0), 1.0)
  assert_eq(ad2.get_x(0, 2), 3.0)
  assert_eq(ad2.get_x(1, 1), 5.0)
  assert_eq(ad2.get_x(2, 0), 7.0)
  assert_eq(ad2.get_x(2, 2), 9.0)
  assert_eq(ad2.get_x(0, 1), 0.0)
  // Names
  assert_eq(ad2.get_obs_names()[0], "c1")
  assert_eq(ad2.get_var_names()[2], "g3")
  // Metadata
  let batch = ad2.obs_column("batch")
  match batch {
    Some(arr) => assert_eq(arr[1], "B")
    None => assert_true(false)
  }
}

///|
test "h5ad_binary_sparse_roundtrip" {
  let coords = [[0, 1], [1, 0], [2, 2]]
  let values = [2.0, 3.0, 5.0]
  let s = @src.SparseArray::from_coo([3, 3], coords, values) catch {
    _ => abort("sparse construction failed")
  }
  let ad = @src.AnnData::new_sparse(s, ["c1", "c2", "c3"], ["g1", "g2", "g3"])
    .set_obs("tissue", ["brain", "liver", "heart"])
  let bytes = @src.scanpy_write_h5ad_binary(ad)
  let ad2 = match @src.scanpy_read_h5ad_binary(bytes) {
    Some(a) => a
    None => abort("binary h5ad parse failed")
  }
  assert_true(ad2.is_sparse())
  assert_eq(ad2.n_obs(), 3)
  assert_eq(ad2.n_vars(), 3)
  assert_eq(ad2.nnz(), 3)
  assert_eq(ad2.get_x(0, 1), 2.0)
  assert_eq(ad2.get_x(1, 0), 3.0)
  assert_eq(ad2.get_x(2, 2), 5.0)
  assert_eq(ad2.get_x(0, 0), 0.0)
  // Metadata
  let tissue = ad2.obs_column("tissue")
  match tissue {
    Some(arr) => {
      assert_eq(arr[0], "brain")
      assert_eq(arr[1], "liver")
      assert_eq(arr[2], "heart")
    }
    None => assert_true(false)
  }
}

///|
test "h5ad_binary_with_all_fields" {
  let x = [[1.0, 2.0], [3.0, 4.0]]
  let ad = @src.AnnData::new(x, ["c1", "c2"], ["g1", "g2"])
    .set_obs("batch", ["A", "B"])
    .set_var("mt", ["no", "yes"])
    .set_layer("counts", [[10.0, 20.0], [30.0, 40.0]])
    .set_obsm("X_pca", [[1.0, 2.0], [3.0, 4.0]])
    .set_varm("VF", [[1.0, 0.0], [0.0, 1.0]])
    .set_obsp("conn", [[1.0, 0.5], [0.5, 1.0]])
    .set_varp("corr", [[1.0, 0.8], [0.8, 1.0]])
    .set_uns("title", "test experiment")
  let bytes = @src.scanpy_write_h5ad_binary(ad)
  let ad2 = match @src.scanpy_read_h5ad_binary(bytes) {
    Some(a) => a
    None => abort("binary h5ad parse failed")
  }
  // Check all fields
  assert_eq(ad2.get_x(0, 0), 1.0)
  assert_eq(ad2.get_x(1, 1), 4.0)
  // obs
  let batch = ad2.obs_column("batch")
  match batch {
    Some(arr) => assert_eq(arr[0], "A")
    None => assert_true(false)
  }
  // var
  let mt = ad2.var_column("mt")
  match mt {
    Some(arr) => assert_eq(arr[1], "yes")
    None => assert_true(false)
  }
  // layers
  let counts = ad2.layers["counts"]
  assert_eq(counts[0][0], 10.0)
  assert_eq(counts[1][1], 40.0)
  // obsm
  let pca = ad2.obsm["X_pca"]
  assert_eq(pca[0][1], 2.0)
  // varm
  let vf = ad2.varm["VF"]
  assert_eq(vf[0][0], 1.0)
  assert_eq(vf[1][1], 1.0)
  // obsp
  let conn = ad2.obsp["conn"]
  assert_eq(conn[0][1], 0.5)
  // varp
  let corr = ad2.varp["corr"]
  assert_eq(corr[1][0], 0.8)
  // uns
  let title = ad2.get_uns("title")
  match title {
    Some(s) => assert_eq(s, "test experiment")
    None => assert_true(false)
  }
}

///|
test "h5ad_binary_invalid_magic" {
  // Not a valid BH5AD file
  let bytes = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]
  let result = @src.scanpy_read_h5ad_binary(bytes)
  match result {
    Some(_) => assert_true(false)
    None => assert_true(true)
  }
}

///|
test "h5ad_binary_empty_matrix" {
  let ad = @src.AnnData::new([], [], [])
  let bytes = @src.scanpy_write_h5ad_binary(ad)
  let ad2 = match @src.scanpy_read_h5ad_binary(bytes) {
    Some(a) => a
    None => abort("binary h5ad parse failed")
  }
  assert_eq(ad2.n_obs(), 0)
  assert_eq(ad2.n_vars(), 0)
}

///|
test "h5ad_json_with_varm_varp_uns_roundtrip" {
  // Verify that the JSON h5ad format also preserves varm/varp/uns
  let x = [[1.0, 2.0], [3.0, 4.0]]
  let ad = @src.AnnData::new(x, ["c1", "c2"], ["g1", "g2"])
    .set_varm("VF", [[1.0, 0.0], [0.0, 1.0]])
    .set_varp("corr", [[1.0, 0.8], [0.8, 1.0]])
    .set_uns("title", "test")
  let json = @src.scanpy_write_h5ad(ad)
  let ad2 = match @src.scanpy_read_h5ad(json) {
    Some(a) => a
    None => abort("json h5ad parse failed")
  }
  // varm
  let vf = ad2.varm["VF"]
  assert_eq(vf[0][0], 1.0)
  assert_eq(vf[1][1], 1.0)
  // varp
  let corr = ad2.varp["corr"]
  assert_eq(corr[0][1], 0.8)
  // uns
  let title = ad2.get_uns("title")
  match title {
    Some(s) => assert_eq(s, "test")
    None => assert_true(false)
  }
}
