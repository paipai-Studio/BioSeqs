///|
/// Tests for the chunked, compressed Zarr v2 store.

///|
fn z2_build() -> @src.AnnData {
  @src.AnnData::new(
    [
      [1.0, 2.0, 3.0, 4.0],
      [1.0, 2.0, 3.0, 4.0],
      [5.0, 6.0, 7.0, 8.0],
      [9.0, 0.0, 1.0, 2.0],
    ],
    ["c1", "c2", "c3", "c4"],
    ["g1", "g2", "g3", "g4"],
  )
  .set_obs("batch", ["A", "A", "B", "B"])
  .set_var("hv", ["true", "true", "false", "false"])
  .set_uns("title", @src.JsonStr("z2"))
}

///|
test "zarr2_chunk_math" {
  assert_eq(@src.anndata_zarr_chunk_count(4, 2), 2)
  assert_eq(@src.anndata_zarr_chunk_count(5, 2), 3)
  assert_eq(@src.anndata_zarr_chunk_count(0, 2), 0)
  let store = @src.ZarrStore::new("/tmp/z", chunk_rows=2, chunk_cols=3)
  // ragged edge chunk: 4x4 grid with 2x3 chunks
  match @src.anndata_zarr_chunk_shape(store, 4, 4, 1, 1) {
    Some((rows, cols)) => {
      assert_eq(rows, 2)
      assert_eq(cols, 1)
    }
    None => fail("edge chunk shape missing")
  }
  assert_eq(@src.anndata_zarr_chunk_shape(store, 4, 4, 2, 0), None)
}

///|
test "zarr2_rle_roundtrip" {
  let lines = ["1,2,3", "1,2,3", "4,5,6", "4,5,6", "4,5,6"]
  let encoded = @src.zarr_rle_encode_lines(lines)
  // compressed form is shorter than raw concatenation
  assert_eq(encoded.has_prefix("1,2,3"), true)
  match @src.zarr_rle_decode_lines(encoded) {
    decoded => {
      assert_eq(decoded.length(), 5)
      assert_eq(decoded[0], "1,2,3")
      assert_eq(decoded[4], "4,5,6")
    }
  }
}

///|
test "zarr2_write_read_roundtrip" {
  let ad = z2_build()
  let store = @src.anndata_zarr_write(
    ad,
    "/demo.zarr",
    chunk_rows=2,
    chunk_cols=2,
  )
  // real zarr v2 structure files exist
  assert_eq(store.has("/demo.zarr/.zgroup"), true)
  assert_eq(store.has("/demo.zarr/X/.zarray"), true)
  assert_eq(store.has("/demo.zarr/X/0.0"), true)
  assert_eq(store.has("/demo.zarr/X/1.1"), true)
  assert_eq(store.has("/demo.zarr/.zattrs"), true)
  // chunk count: 2x2 chunks of 2x2
  assert_eq(store.n_files(), 3 + 4 + 2 + 2)
  match @src.anndata_zarr_read(store) {
    Some(back) => {
      assert_eq(back.n_obs(), 4)
      assert_eq(back.n_vars(), 4)
      assert_eq(back.get_x(0, 0), 1.0)
      assert_eq(back.get_x(2, 3), 8.0)
      assert_eq(back.get_x(3, 0), 9.0)
      assert_eq(back.obs_names[2], "c3")
      assert_eq(back.var_names[3], "g4")
      match back.obs_column("batch") {
        Some(col) => assert_eq(col[0], "A")
        None => fail("obs batch missing")
      }
      match back.var_column("hv") {
        Some(col) => assert_eq(col[3], "false")
        None => fail("var hv missing")
      }
      match back.get_uns("title") {
        Some(@src.JsonStr(s)) => assert_eq(s, "z2")
        _ => fail("uns title missing from zarr")
      }
    }
    None => fail("zarr read failed")
  }
}

///|
test "zarr2_write_read_roundtrip_compressed" {
  let ad = z2_build()
  let store = @src.anndata_zarr_write(
    ad,
    "/c.zarr",
    chunk_rows=2,
    chunk_cols=2,
    compressor="rle",
  )
  // compressor recorded in .zarray metadata
  match store.get("/c.zarr/X/.zarray") {
    Some(content) => assert_eq(content.contains("zlib-rle"), true)
    None => fail(".zarray missing")
  }
  match @src.anndata_zarr_read(store) {
    Some(back) => {
      assert_eq(back.get_x(0, 0), 1.0)
      assert_eq(back.get_x(2, 2), 7.0)
      assert_eq(back.get_x(3, 3), 2.0)
    }
    None => fail("compressed zarr read failed")
  }
}

///|
test "zarr2_hierarchy_sorted" {
  let ad = z2_build()
  let store = @src.anndata_zarr_write(ad, "/h.zarr", chunk_rows=2, chunk_cols=2)
  let paths = @src.anndata_zarr_hierarchy(store)
  // sorted output
  let mut sorted = true
  for i = 1; i < paths.length(); i = i + 1 {
    if paths[i - 1] > paths[i] {
      sorted = false
    }
  }
  assert_eq(sorted, true)
  assert_eq(paths.length(), store.n_files())
}
