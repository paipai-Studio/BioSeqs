///|
/// Tests for generic write/read dispatch by filename extension.

///|
fn dp_build() -> @src.AnnData {
  @src.AnnData::new([[1.0, 2.0], [3.0, 4.0]], ["c1", "c2"], ["g1", "g2"])
  .set_obs("batch", ["A", "B"])
  .set_var("hv", ["true", "false"])
}

///|
test "dispatch_format_detection" {
  assert_eq(@src.anndata_write_format("data.h5ad"), "h5ad")
  assert_eq(@src.anndata_write_format("store.zarr"), "zarr")
  assert_eq(@src.anndata_write_format("sample.loom"), "loom")
  assert_eq(@src.anndata_write_format("table.csv"), "csvs")
  assert_eq(@src.anndata_write_format("table.csvs"), "csvs")
  // unknown extension falls back to h5ad (like AnnData.write)
  assert_eq(@src.anndata_write_format("data.bin"), "h5ad")
  assert_eq(@src.anndata_write_format("data"), "h5ad")
}

///|
test "dispatch_h5ad_roundtrip" {
  let ad = dp_build()
  let content = @src.anndata_write(ad, "a.h5ad")
  match @src.anndata_read("a.h5ad", content) {
    Some(back) => {
      assert_eq(back.n_obs(), 2)
      assert_eq(back.get_x(0, 0), 1.0)
      match back.obs_column("batch") {
        Some(col) => assert_eq(col[1], "B")
        None => fail("obs missing")
      }
    }
    None => fail("h5ad dispatch roundtrip failed")
  }
}

///|
test "dispatch_zarr_roundtrip" {
  let ad = dp_build()
  let content = @src.anndata_write(ad, "a.zarr")
  match @src.anndata_read("a.zarr", content) {
    Some(back) => {
      assert_eq(back.n_obs(), 2)
      assert_eq(back.n_vars(), 2)
      assert_eq(back.get_x(1, 1), 4.0)
    }
    None => fail("zarr dispatch roundtrip failed")
  }
}

///|
test "dispatch_loom_roundtrip" {
  let ad = dp_build()
  let content = @src.anndata_write(ad, "a.loom")
  match @src.anndata_read("a.loom", content) {
    Some(back) => {
      assert_eq(back.n_obs(), 2)
      match back.var_column("hv") {
        Some(col) => assert_eq(col[0], "true")
        None => fail("var missing after loom roundtrip")
      }
    }
    None => fail("loom dispatch roundtrip failed")
  }
}

///|
test "dispatch_csvs_roundtrip" {
  let ad = dp_build()
  let content = @src.anndata_write(ad, "a.csv")
  // bundle carries the three sections
  assert_eq(content.contains("===OBS==="), true)
  assert_eq(content.contains("===VAR==="), true)
  match @src.anndata_read("a.csv", content) {
    Some(back) => {
      assert_eq(back.n_obs(), 2)
      assert_eq(back.get_x(0, 1), 2.0)
    }
    None => fail("csvs dispatch roundtrip failed")
  }
}

///|
test "dispatch_roundtrip_helper" {
  let ad = dp_build()
  match @src.anndata_write_read_roundtrip(ad, "rt.h5ad") {
    Some(back) => assert_eq(back.get_x(1, 0), 3.0)
    None => fail("roundtrip helper failed")
  }
  match @src.anndata_write_read_roundtrip(ad, "rt.zarr") {
    Some(back) => assert_eq(back.get_x(1, 0), 3.0)
    None => fail("roundtrip helper failed for zarr")
  }
}
