///|
/// Tests for AnnData raw attachment (adata.raw as a first-class slot).

///|
test "anndata_raw_attach_and_read" {
  let x = [[1.0, 2.0], [3.0, 4.0]]
  let ad = @src.AnnData::new(x, ["c1", "c2"], ["g1", "g2"])
  assert_eq(ad.has_raw(), false)
  let raw = @src.anndata_raw_capture(ad)
  // mutate X afterwards; the raw snapshot keeps the original values
  ad.x[0][0] = 100.0
  let ad2 = ad.set_raw(raw)
  assert_eq(ad2.has_raw(), true)
  match ad2.get_raw() {
    Some(r) => {
      assert_eq(r.get_x(0, 0), 1.0)
      assert_eq(r.get_x(0, 1), 2.0)
      assert_eq(r.n_obs(), 2)
      assert_eq(r.n_vars(), 2)
      assert_eq(r.var_names()[1], "g2")
    }
    None => fail("raw should be attached")
  }
}

///|
test "anndata_raw_survives_subset" {
  let x = [[1.0, 2.0, 3.0], [4.0, 5.0, 6.0], [7.0, 8.0, 9.0]]
  let ad = @src.AnnData::new(x, ["c1", "c2", "c3"], ["g1", "g2", "g3"])
  let raw = @src.anndata_raw_capture(ad)
  let sub = ad.set_raw(raw).subset([2, 0], [0, 1])
  match sub.get_raw() {
    Some(r) => {
      // raw rows follow obs_idx [2, 0]; raw cols keep all genes
      assert_eq(r.n_obs(), 2)
      assert_eq(r.get_x(0, 0), 7.0)
      assert_eq(r.get_x(1, 0), 1.0)
    }
    None => fail("raw should follow subset")
  }
}

///|
test "anndata_raw_survives_copy_and_sparse_switch" {
  let x = [[1.0, 2.0], [3.0, 4.0]]
  let ad = @src.AnnData::new(x, ["c1", "c2"], ["g1", "g2"])
  let raw = @src.anndata_raw_capture(ad)
  let with_raw = ad.set_raw(raw)
  let cp = with_raw.copy()
  assert_eq(cp.has_raw(), true)
  let sp = with_raw.to_sparse()
  assert_eq(sp.has_raw(), true)
  let de = sp.to_dense()
  assert_eq(de.has_raw(), true)
}

///|
test "anndata_raw_subset_obs_api" {
  let x = [[1.0, 2.0], [3.0, 4.0], [5.0, 6.0]]
  let ad = @src.AnnData::new(x, ["c1", "c2", "c3"], ["g1", "g2"])
  let raw = @src.anndata_raw_capture(ad)
  let sub = raw.subset_obs([0, 2])
  assert_eq(sub.n_obs(), 2)
  assert_eq(sub.get_x(0, 0), 1.0)
  assert_eq(sub.get_x(1, 1), 6.0)
  // column access
  let col = sub.get_col(0)
  assert_eq(col[0], 1.0)
  assert_eq(col[1], 5.0)
  // var annotation carried over
  assert_eq(sub.var_column("missing"), None)
  assert_eq(sub.sum(), 1.0 + 2.0 + 5.0 + 6.0)
}
