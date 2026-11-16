///|
/// Tests for 5 new AnnData features:
///   1. obs_vector(k) / var_vector(k)
///   2. rename_categories(key, mapping)
///   3. read_mtx / write_mtx
///   4. AnnCollectionView

// ============================================================================
// Fixture: small 3x2 AnnData with obs/var names + one layer
// ============================================================================

fn make_fixture() -> @src.AnnData {
  let x = [[1.0, 2.0], [3.0, 4.0], [5.0, 6.0]]
  let ad = @src.AnnData::new(x, ["Cell_1", "Cell_2", "Cell_3"], [
    "Gene_1", "Gene_2",
  ])
  let layer_mat = [[10.0, 20.0], [30.0, 40.0], [50.0, 60.0]]
  let _ = ad.set_layer("counts", layer_mat)
  ad
}

// ============================================================================
// 1. obs_vector(k) tests
// ============================================================================

///|
test "obs_vector_by_var_name" {
  let ad = make_fixture()
  let v = match ad.obs_vector("Gene_1") {
    Some(v) => v
    None => fail("expected Some")
  }
  assert_eq(v.length(), 3)
  assert_true(@moonbit.assert_num_close(v[0], 1.0, 1.0e-10))
  assert_true(@moonbit.assert_num_close(v[1], 3.0, 1.0e-10))
  assert_true(@moonbit.assert_num_close(v[2], 5.0, 1.0e-10))
}

///|
test "obs_vector_by_var_int_index_string" {
  let ad = make_fixture()
  let v0 = match ad.obs_vector("0") {
    Some(v) => v
    None => fail("expected Some")
  }
  let v1 = match ad.obs_vector("1") {
    Some(v) => v
    None => fail("expected Some")
  }
  assert_eq(v0.length(), 3)
  assert_true(@moonbit.assert_num_close(v0[0], 1.0, 1.0e-10))
  assert_true(@moonbit.assert_num_close(v1[0], 2.0, 1.0e-10))
  assert_true(@moonbit.assert_num_close(v1[2], 6.0, 1.0e-10))
}

///|
test "obs_vector_with_layer" {
  let ad = make_fixture()
  let v = match ad.obs_vector("Gene_2", layer="counts") {
    Some(v) => v
    None => fail("expected Some")
  }
  assert_eq(v.length(), 3)
  assert_true(@moonbit.assert_num_close(v[0], 20.0, 1.0e-10))
  assert_true(@moonbit.assert_num_close(v[1], 40.0, 1.0e-10))
  assert_true(@moonbit.assert_num_close(v[2], 60.0, 1.0e-10))
}

///|
test "obs_vector_not_found_returns_none" {
  let ad = make_fixture()
  match ad.obs_vector("Gene_XYZ") {
    Some(_) => fail("should be None")
    None => ()
  }
  match ad.obs_vector("Gene_1", layer="nonexistent") {
    Some(_) => fail("bad layer should be None")
    None => ()
  }
}

// ============================================================================
// 2. var_vector(k) tests
// ============================================================================

///|
test "var_vector_by_obs_name" {
  let ad = make_fixture()
  let v = match ad.var_vector("Cell_2") {
    Some(v) => v
    None => fail("expected Some")
  }
  assert_eq(v.length(), 2)
  assert_true(@moonbit.assert_num_close(v[0], 3.0, 1.0e-10))
  assert_true(@moonbit.assert_num_close(v[1], 4.0, 1.0e-10))
}

///|
test "var_vector_by_obs_int_index_string" {
  let ad = make_fixture()
  let v = match ad.var_vector("2") {
    Some(v) => v
    None => fail("expected Some")
  }
  assert_eq(v.length(), 2)
  assert_true(@moonbit.assert_num_close(v[0], 5.0, 1.0e-10))
  assert_true(@moonbit.assert_num_close(v[1], 6.0, 1.0e-10))
}

///|
test "var_vector_with_layer" {
  let ad = make_fixture()
  let v = match ad.var_vector("Cell_1", layer="counts") {
    Some(v) => v
    None => fail("expected Some")
  }
  assert_eq(v.length(), 2)
  assert_true(@moonbit.assert_num_close(v[0], 10.0, 1.0e-10))
  assert_true(@moonbit.assert_num_close(v[1], 20.0, 1.0e-10))
}

///|
test "var_vector_not_found_returns_none" {
  let ad = make_fixture()
  match ad.var_vector("Cell_Nope") {
    Some(_) => fail("should be None")
    None => ()
  }
}

// ============================================================================
// 3. rename_categories tests
// ============================================================================

///|
test "rename_categories_roundtrip" {
  let ad = make_fixture().set_obs("batch", ["WT", "KO", "WT"])
  @src.anndata_obs_factorize(ad, "batch")
  let before = match ad.obs_column("batch") {
    Some(c) => c
    None => fail("batch missing")
  }
  assert_eq(before[0], "0")
  assert_eq(before[1], "1")
  assert_eq(before[2], "0")
  let ok = ad.rename_categories("batch", [
    ("WT", "WildType"),
    ("KO", "Knockout"),
  ])
  assert_true(ok)
  let after = match ad.obs_column("batch") {
    Some(c) => c
    None => fail("batch missing after rename")
  }
  assert_eq(after[0], "WildType")
  assert_eq(after[1], "Knockout")
  assert_eq(after[2], "WildType")
  match ad.get_uns("batch_categories") {
    Some(@src.JsonArr(arr)) => {
      match arr[0] {
        @src.JsonStr(s) => assert_eq(s, "WildType")
        _ => fail("first cat wrong")
      }
      match arr[1] {
        @src.JsonStr(s) => assert_eq(s, "Knockout")
        _ => fail("second cat wrong")
      }
    }
    _ => fail("uns categories missing")
  }
}

///|
test "rename_categories_list_syntax" {
  let ad = make_fixture().set_obs("batch", ["A", "B", "A"])
  @src.anndata_obs_factorize(ad, "batch")
  let ok = ad.rename_categories_list("batch", ["Alpha", "Beta"])
  assert_true(ok)
  let col = match ad.obs_column("batch") {
    Some(c) => c
    None => fail("missing")
  }
  assert_eq(col[0], "Alpha")
  assert_eq(col[1], "Beta")
  assert_eq(col[2], "Alpha")
}

///|
test "AnnDataCategorical_rename_method" {
  let cat = @src.anndata_encode_categorical(["A", "B", "A", "C"])
  assert_eq(cat.n_categories(), 3)
  let renamed = cat.rename([("A", "Alpha"), ("C", "Charlie")])
  assert_eq(renamed.n_categories(), 3)
  assert_eq(renamed.category(0), "Alpha")
  assert_eq(renamed.category(1), "B")
  assert_eq(renamed.category(2), "Charlie")
  assert_eq(renamed.codes.length(), cat.codes.length())
}

// ============================================================================
// 4. read_mtx / write_mtx roundtrip
// ============================================================================

///|
test "write_mtx_header_format" {
  let ad = make_fixture()
  let mtx = @src.anndata_write_mtx(ad)
  let lines = mtx.split("\n").to_array()
  // 直接检查行内容
  let first_line = lines[0].to_string()
  assert_true(first_line.has_prefix("%%MatrixMarket"))
  let second_line = lines[1].to_string()
  assert_eq(second_line, "%")
  let header = lines[2].to_string().split(" ").to_array()
  assert_eq(header[0], "3")
  assert_eq(header[1], "2")
  assert_eq(header[2], "6")
}

///|
test "read_mtx_basic" {
  let mtx = "%%MatrixMarket matrix coordinate real general\n%\n3 2 6\n1 1 1\n1 2 2\n2 1 3\n2 2 4\n3 1 5\n3 2 6\n"
  let ad = match @src.anndata_read_mtx(mtx) {
    Some(a) => a
    None => fail("expected Some AnnData")
  }
  assert_eq(ad.n_obs(), 3)
  assert_eq(ad.n_vars(), 2)
  let x = ad.x_matrix()
  assert_true(@moonbit.assert_num_close(x[0][0], 1.0, 1.0e-10))
  assert_true(@moonbit.assert_num_close(x[0][1], 2.0, 1.0e-10))
  assert_true(@moonbit.assert_num_close(x[2][0], 5.0, 1.0e-10))
  assert_true(@moonbit.assert_num_close(x[2][1], 6.0, 1.0e-10))
  assert_eq(ad.obs_names[0], "Cell_1")
  assert_eq(ad.var_names[1], "Gene_2")
}

///|
test "read_mtx_sparse_zeros" {
  let mtx = "%%MatrixMarket matrix coordinate real general\n%\n3 3 4\n1 1 10\n2 2 20\n3 3 30\n1 3 5\n"
  let ad = match @src.anndata_read_mtx(mtx) {
    Some(a) => a
    None => fail("expected Some")
  }
  let x = ad.x_matrix()
  assert_true(@moonbit.assert_num_close(x[0][0], 10.0, 1.0e-10))
  assert_true(@moonbit.assert_num_close(x[0][2], 5.0, 1.0e-10))
  assert_true(@moonbit.assert_num_close(x[1][1], 20.0, 1.0e-10))
  assert_true(@moonbit.assert_num_close(x[2][2], 30.0, 1.0e-10))
  assert_true(@moonbit.assert_num_close(x[0][1], 0.0, 1.0e-10))
  assert_true(@moonbit.assert_num_close(x[1][0], 0.0, 1.0e-10))
}

///|
test "mtx_roundtrip_write_read" {
  let ad_orig = make_fixture()
  let mtx = @src.anndata_write_mtx(ad_orig)
  let ad_back = match @src.anndata_read_mtx(mtx) {
    Some(a) => a
    None => fail("roundtrip failed")
  }
  assert_eq(ad_back.n_obs(), ad_orig.n_obs())
  assert_eq(ad_back.n_vars(), ad_orig.n_vars())
  let xo = ad_orig.x_matrix()
  let xb = ad_back.x_matrix()
  for i = 0; i < xo.length(); i = i + 1 {
    for j = 0; j < xo[i].length(); j = j + 1 {
      assert_true(@moonbit.assert_num_close(xo[i][j], xb[i][j], 1.0e-10))
    }
  }
}

///|
test "read_mtx_empty_or_bad_returns_none" {
  match @src.anndata_read_mtx("") {
    None => ()
    Some(_) => fail("empty should be None")
  }
}

// ============================================================================
// 5. AnnCollectionView tests
// ============================================================================

///|
fn make_two_collections() -> (@src.AnnData, @src.AnnData, @src.AnnCollection) {
  let a = @src.AnnData::new([[1.0, 2.0], [3.0, 4.0]], ["c1", "c2"], ["g1", "g2"])
  let b = @src.AnnData::new([[5.0, 6.0]], ["c3"], ["g1", "g2"])
  let coll = @src.AnnCollection::from([a, b], label="batch")
  (a, b, coll)
}

///|
test "AnnCollectionView_full_no_filter" {
  let (_, _, coll) = make_two_collections()
  let view = @src.AnnCollectionView::from(coll)
  assert_eq(view.n_obs(), 3)
  assert_eq(view.n_vars(), 2)
  let v = match view.get_x(0, 0) {
    Some(v) => v
    None => fail("expected Some")
  }
  assert_true(@moonbit.assert_num_close(v, 1.0, 1.0e-10))
  let v2 = match view.get_x(2, 1) {
    Some(v) => v
    None => fail("expected Some")
  }
  assert_true(@moonbit.assert_num_close(v2, 6.0, 1.0e-10))
}

///|
test "AnnCollectionView_with_obs_filter" {
  let (_, _, coll) = make_two_collections()
  let view = @src.AnnCollectionView::with_obs(coll, [1, 2])
  assert_eq(view.n_obs(), 2)
  assert_eq(view.n_vars(), 2)
  let row0 = match view.get_obs(0) {
    Some(r) => r
    None => fail("expected Some")
  }
  assert_true(@moonbit.assert_num_close(row0[0], 3.0, 1.0e-10))
  assert_true(@moonbit.assert_num_close(row0[1], 4.0, 1.0e-10))
  let row1 = match view.get_obs(1) {
    Some(r) => r
    None => fail("expected Some")
  }
  assert_true(@moonbit.assert_num_close(row1[0], 5.0, 1.0e-10))
}

///|
test "AnnCollectionView_with_var_filter" {
  let (_, _, coll) = make_two_collections()
  let view = @src.AnnCollectionView::with_vars(coll, [1])
  assert_eq(view.n_obs(), 3)
  assert_eq(view.n_vars(), 1)
  let v0 = match view.get_x(0, 0) {
    Some(v) => v
    None => fail("expected Some")
  }
  assert_true(@moonbit.assert_num_close(v0, 2.0, 1.0e-10))
  let v2 = match view.get_x(2, 0) {
    Some(v) => v
    None => fail("expected Some")
  }
  assert_true(@moonbit.assert_num_close(v2, 6.0, 1.0e-10))
}

///|
test "AnnCollectionView_with_both_filters" {
  let (_, _, coll) = make_two_collections()
  let view = @src.AnnCollectionView::with_both(coll, [0, 2], [1])
  assert_eq(view.n_obs(), 2)
  assert_eq(view.n_vars(), 1)
  let v0 = match view.get_x(0, 0) {
    Some(v) => v
    None => fail("expected Some")
  }
  assert_true(@moonbit.assert_num_close(v0, 2.0, 1.0e-10))
  let v1 = match view.get_x(1, 0) {
    Some(v) => v
    None => fail("expected Some")
  }
  assert_true(@moonbit.assert_num_close(v1, 6.0, 1.0e-10))
}

///|
test "AnnCollectionView_to_adata" {
  let (_, _, coll) = make_two_collections()
  let view = @src.AnnCollectionView::with_both(coll, [1, 2], [0])
  let ad = view.to_adata()
  assert_eq(ad.n_obs(), 2)
  assert_eq(ad.n_vars(), 1)
  let x = ad.x_matrix()
  assert_true(@moonbit.assert_num_close(x[0][0], 3.0, 1.0e-10))
  assert_true(@moonbit.assert_num_close(x[1][0], 5.0, 1.0e-10))
  assert_eq(ad.var_names[0], "g1")
  assert_eq(ad.obs_names[0], "c2")
  assert_eq(ad.obs_names[1], "c3")
}

///|
test "AnnCollectionView_out_of_bounds_returns_none" {
  let (_, _, coll) = make_two_collections()
  let view = @src.AnnCollectionView::from(coll)
  match view.get_obs(99) {
    None => ()
    Some(_) => fail("out of bounds should be None")
  }
  match view.get_x(0, 99) {
    None => ()
    Some(_) => fail("out of bounds should be None")
  }
}
