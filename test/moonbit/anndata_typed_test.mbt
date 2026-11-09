///|
/// Tests for typed / nullable obs & var columns.

///|
fn tc_col_str(ad : @src.AnnData, key : String, i : Int32) -> String {
  match ad.obs_column(key) {
    Some(arr) => if i < arr.length() { arr[i] } else { "MISSING" }
    None => "MISSING"
  }
}

///|
test "typed_col_construct_and_strings" {
  let c = @src.anndata_col_int("n_genes", [100, 200, 300])
  assert_eq(c.kind_name(), "integer")
  assert_eq(c.len(), 3)
  assert_eq(c.n_missing(), 0)
  let strs = c.to_strings()
  assert_eq(strs[0], "100")
  assert_eq(strs[2], "300")
  let d = @src.anndata_col_double("score", [0.5, 1.5])
  assert_eq(d.kind_name(), "double")
  assert_eq(d.to_strings()[1], "1.5")
  let b = @src.anndata_col_bool("is_tumor", [true, false])
  assert_eq(b.kind_name(), "boolean")
  assert_eq(b.to_strings()[0], "true")
  let cat = @src.anndata_col_cat("batch", ["A", "B"], [0, 1, 1])
  assert_eq(cat.kind_name(), "categorical")
  assert_eq(cat.cat_value(2), "B")
}

///|
test "typed_col_nullable_missing" {
  let na : Array[Bool] = [false, true, false]
  let c = @src.anndata_col_int_nullable("n_genes", [10, 0, 30], na)
  assert_eq(c.n_missing(), 1)
  assert_eq(c.is_na(1), true)
  assert_eq(c.is_na(0), false)
  // missing entries stringify to ""
  assert_eq(c.to_strings()[1], "")
  let db = @src.anndata_col_double_nullable("log1p_total", [1.5, 0.0], na)
  assert_eq(db.n_missing(), 1)
  let bo = @src.anndata_col_bool_nullable("qc_pass", [true, false, true], na)
  assert_eq(bo.n_missing(), 1)
}

///|
test "typed_obs_roundtrip_int_double_bool" {
  let ad = @src.AnnData::new([[1.0], [2.0], [3.0]], ["c1", "c2", "c3"], ["g1"])
  ignore(ad.set_obs_typed(@src.anndata_col_int("n_genes", [100, 200, 300])))
  ignore(ad.set_obs_typed(@src.anndata_col_double("score", [0.5, 1.5, 2.5])))
  ignore(ad.set_obs_typed(@src.anndata_col_bool("pass", [true, false, true])))
  // stringified values landed in obs
  assert_eq(tc_col_str(ad, "n_genes", 1), "200")
  // typed view recovers the dtype
  match ad.get_obs_typed("n_genes") {
    Some(c) => {
      assert_eq(c.kind_name(), "integer")
      match c.ints.get(2) {
        Some(v) => assert_eq(v, 300)
        None => fail("ints[2] missing")
      }
    }
    None => fail("typed n_genes missing")
  }
  match ad.get_obs_typed("score") {
    Some(c) =>
      match c.doubles.get(1) {
        Some(v) => assert_eq(v, 1.5)
        None => fail("doubles[1] missing")
      }
    None => fail("typed score missing")
  }
  match ad.get_obs_typed("pass") {
    Some(c) => {
      assert_eq(c.kind_name(), "boolean")
      match c.bools.get(1) {
        Some(v) => assert_eq(v, false)
        None => fail("bools[1] missing")
      }
    }
    None => fail("typed pass missing")
  }
}

///|
test "typed_obs_nullable_roundtrip" {
  let ad = @src.AnnData::new([[1.0], [2.0], [3.0]], ["c1", "c2", "c3"], ["g1"])
  let na : Array[Bool] = [false, true, false]
  ignore(
    ad.set_obs_typed(@src.anndata_col_int_nullable("n_genes", [10, 0, 30], na)),
  )
  match ad.get_obs_typed("n_genes") {
    Some(c) => {
      assert_eq(c.n_missing(), 1)
      assert_eq(c.is_na(1), true)
      match c.ints.get(2) {
        Some(v) => assert_eq(v, 30)
        None => fail("ints[2] missing")
      }
    }
    None => fail("typed nullable missing")
  }
}

///|
test "typed_var_roundtrip_and_categorical" {
  let ad = @src.AnnData::new([[1.0, 2.0]], ["c1"], ["g1", "g2"])
  ignore(ad.set_var_typed(@src.anndata_col_double("mean", [0.25, 0.75])))
  match ad.get_var_typed("mean") {
    Some(c) => {
      assert_eq(c.kind_name(), "double")
      match c.doubles.get(0) {
        Some(v) => assert_eq(v, 0.25)
        None => fail("doubles[0] missing")
      }
    }
    None => fail("typed var missing")
  }
  // categorical column round-trips through uns categories
  ignore(ad.set_obs_typed(@src.anndata_col_cat("batch", ["A", "B"], [0, 1, 1])))
  match ad.get_obs_typed("batch") {
    Some(c) => {
      assert_eq(c.kind_name(), "categorical")
      assert_eq(c.cat_value(0), "A")
      assert_eq(c.cat_value(2), "B")
    }
    None => fail("typed batch missing")
  }
}

///|
test "typed_take_subset" {
  let na : Array[Bool] = [false, true, false, false]
  let c = @src.anndata_col_double_nullable("score", [1.0, 2.0, 3.0, 4.0], na)
  let taken = c.take([2, 0])
  assert_eq(taken.len(), 2)
  match taken.doubles.get(0) {
    Some(v) => assert_eq(v, 3.0)
    None => fail("take[0] missing")
  }
  match taken.doubles.get(1) {
    Some(v) => assert_eq(v, 1.0)
    None => fail("take[1] missing")
  }
  assert_eq(taken.n_missing(), 0)
  assert_eq(taken.kind_name(), "double")
}

///|
test "typed_missing_column_returns_none" {
  let ad = @src.AnnData::new([[1.0]], ["c1"], ["g1"])
  match ad.get_obs_typed("nope") {
    None => ()
    Some(_) => fail("unknown obs column should be None")
  }
  match ad.get_var_typed("nope") {
    None => ()
    Some(_) => fail("unknown var column should be None")
  }
}
