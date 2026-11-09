///|
/// Tests for view implicit-modification warning semantics
/// (SettingWithCopyWarning / ImplicitModificationWarning analogue).

///|
fn vw_build() -> @src.AnnData {
  @src.AnnData::new(
    [[1.0, 2.0, 3.0], [4.0, 5.0, 6.0], [7.0, 8.0, 9.0]],
    ["c1", "c2", "c3"],
    ["g1", "g2", "g3"],
  )
}

///|
test "view_set_x_full_no_warning" {
  let ad = vw_build()
  let view = @src.anndata_view(ad)
  let warning = view.set_x(0, 0, 99.0)
  assert_eq(warning, None)
  assert_eq(view.is_modified(), true)
  // full view writes still reach the source
  assert_eq(ad.get_x(0, 0), 99.0)
  assert_eq(view.warning(), None)
}

///|
test "view_set_x_subset_warns" {
  let ad = vw_build()
  let view = @src.anndata_view_obs(ad, [0, 2])
  match view.set_x(1, 0, 70.0) {
    Some(w) => {
      // the pandas/anndata warning text
      assert_eq(w.has_prefix("ImplicitModificationWarning"), true)
      assert_eq(w.contains(".copy()"), true)
      assert_eq(w.contains("slice from an AnnData"), true)
    }
    None => fail("subset write should warn")
  }
  assert_eq(view.is_modified(), true)
  assert_eq(view.needs_copy(), true)
  // write went through to the underlying source row (obs index 2)
  assert_eq(ad.get_x(2, 0), 70.0)
  // the pending warning is retrievable
  match view.warning() {
    Some(_) => ()
    None => fail("warning should be pending")
  }
}

///|
test "view_set_x_out_of_range_no_write" {
  let ad = vw_build()
  let view = @src.anndata_view_obs(ad, [0, 2])
  // out of view range: no write, no modification
  assert_eq(view.set_x(5, 0, 1.0), None)
  assert_eq(view.is_modified(), false)
  assert_eq(view.needs_copy(), false)
  assert_eq(view.warning(), None)
}

///|
test "view_writeback_marks_modified" {
  let ad = vw_build()
  let view = @src.anndata_view_both(ad, [0, 1], [0, 1])
  assert_eq(view.is_modified(), false)
  view.writeback([[10.0, 20.0], [30.0, 40.0]])
  assert_eq(view.is_modified(), true)
  assert_eq(view.needs_copy(), true)
  assert_eq(ad.get_x(0, 0), 10.0)
  assert_eq(ad.get_x(1, 1), 40.0)
  // the materialized copy is independent of the source
  let cp = view.copy()
  assert_eq(cp.get_x(0, 0), 10.0)
}

///|
test "view_var_and_both_views_warn" {
  let ad = vw_build()
  let vvar = @src.anndata_view_var(ad, [1])
  match vvar.set_x(0, 0, 5.5) {
    Some(_) => ()
    None => fail("var subset write should warn")
  }
  assert_eq(ad.get_x(0, 1), 5.5)
  let vboth = @src.anndata_view_both(ad, [1], [1])
  match vboth.set_x(0, 0, 8.5) {
    Some(_) => ()
    None => fail("both-axes subset write should warn")
  }
  assert_eq(ad.get_x(1, 1), 8.5)
}
