///|
/// Tests for the AnnData single-cell container.

///|
fn ads_col(values : Array[String]?, i : Int) -> String {
  match values {
    Some(arr) => if i < arr.length() { arr[i] } else { "MISSING" }
    None => "MISSING"
  }
}

///|
test "anndata_construct_shape" {
  let x = [[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]]
  let ad = @src.AnnData::new(x, ["c1", "c2"], ["g1", "g2", "g3"])
  assert_eq(ad.n_obs(), 2)
  assert_eq(ad.n_vars(), 3)
  let (n, m) = ad.shape()
  assert_eq(n, 2)
  assert_eq(m, 3)
  assert_eq(ad.obs_names[0], "c1")
  assert_eq(ad.var_names[2], "g3")
  assert_eq(ad.x[1][2], 6.0)
}

///|
test "anndata_set_and_get_columns" {
  let x = [[1.0, 2.0], [3.0, 4.0]]
  let ad = @src.AnnData::new(x, ["c1", "c2"], ["g1", "g2"])
    .set_obs("batch", ["A", "B"])
    .set_var("gene", ["G1", "G2"])
  assert_eq(ads_col(ad.obs_column("batch"), 0), "A")
  assert_eq(ads_col(ad.obs_column("batch"), 1), "B")
  assert_eq(ads_col(ad.var_column("gene"), 1), "G2")
  assert_eq(ads_col(ad.obs_column("nonexistent"), 0), "MISSING")
}

///|
test "anndata_subset" {
  let x = [
    [1.0, 2.0, 3.0, 4.0],
    [5.0, 6.0, 7.0, 8.0],
    [9.0, 10.0, 11.0, 12.0],
  ]
  let ad = @src.AnnData::new(x, ["c1", "c2", "c3"], ["g1", "g2", "g3", "g4"])
    .set_obs("batch", ["A", "B", "A"])
    .set_var("gene", ["G1", "G2", "G3", "G4"])
    .set_obsm("X_pca", [[1.0, 2.0], [3.0, 4.0], [5.0, 6.0]])

  let sub = ad.subset([2, 0], [1, 3])
  assert_eq(sub.n_obs(), 2)
  assert_eq(sub.n_vars(), 2)
  assert_eq(sub.obs_names[0], "c3")
  assert_eq(sub.obs_names[1], "c1")
  assert_eq(sub.var_names[0], "g2")
  assert_eq(sub.var_names[1], "g4")
  // X reindexed: [2,0] rows × [1,3] cols
  assert_eq(sub.x[0][0], 10.0)
  assert_eq(sub.x[0][1], 12.0)
  assert_eq(sub.x[1][0], 2.0)
  assert_eq(sub.x[1][1], 4.0)
  // obs column subset
  assert_eq(ads_col(sub.obs_column("batch"), 0), "A")
  assert_eq(ads_col(sub.obs_column("batch"), 1), "A")
  // var column subset
  assert_eq(ads_col(sub.var_column("gene"), 0), "G2")
  assert_eq(ads_col(sub.var_column("gene"), 1), "G4")
  // obsm subset along obs
  let pca = sub.obsm["X_pca"]
  assert_eq(pca[0][0], 5.0)
  assert_eq(pca[0][1], 6.0)
  assert_eq(pca[1][0], 1.0)
}

///|
test "anndata_copy_is_deep" {
  let x = [[1.0, 2.0], [3.0, 4.0]]
  let ad = @src.AnnData::new(x, ["c1", "c2"], ["g1", "g2"]).set_obs("batch", ["A", "B"])
  let cp = ad.copy()
  let cp2 = cp.set_obs("batch", ["X", "Y"])
  // original obs unchanged by mutating the copy
  assert_eq(ads_col(ad.obs_column("batch"), 0), "A")
  assert_eq(ads_col(cp2.obs_column("batch"), 0), "X")
  // X deep-copied: mutating copy does not touch original
  cp.x[0][0] = 99.0
  assert_eq(ad.x[0][0], 1.0)
  assert_eq(cp.x[0][0], 99.0)
}