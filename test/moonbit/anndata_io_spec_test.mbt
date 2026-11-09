///|
/// Tests for the full h5ad encoding-type / encoding-version spec.

///|
test "io_spec_registry_names_and_versions" {
  assert_eq(@src.anndata_enc_name(@src.SpecArray), "array")
  assert_eq(@src.anndata_enc_version(@src.SpecArray), "0.2.0")
  assert_eq(@src.anndata_enc_name(@src.SpecDataFrame), "dataframe")
  assert_eq(@src.anndata_enc_version(@src.SpecDataFrame), "0.2.0")
  assert_eq(@src.anndata_enc_name(@src.SpecCSR), "csr_matrix")
  assert_eq(@src.anndata_enc_version(@src.SpecCSR), "0.1.0")
  assert_eq(
    @src.anndata_enc_name(@src.SpecCategoricalMasked),
    "categorical-masked",
  )
  assert_eq(@src.anndata_enc_name(@src.SpecNullableInteger), "nullable-integer")
  assert_eq(@src.anndata_enc_name(@src.SpecNullableBoolean), "nullable-boolean")
  assert_eq(
    @src.anndata_enc_name(@src.SpecNullableStringArray),
    "nullable-string-array",
  )
  assert_eq(@src.anndata_enc_name(@src.SpecStringArray), "string-array")
  assert_eq(@src.anndata_enc_name(@src.SpecDict), "dict")
  // lookup by name round-trips
  match @src.anndata_enc_from_name("dataframe") {
    Some(e) => assert_eq(@src.anndata_enc_name(e), "dataframe")
    None => fail("dataframe lookup failed")
  }
  match @src.anndata_enc_from_name("bogus") {
    None => ()
    Some(_) => fail("bogus encoding should not resolve")
  }
}

///|
test "io_spec_attrs_and_check" {
  let attrs = @src.anndata_enc_attrs(@src.SpecArray)
  match attrs {
    @src.JsonObj(m) => {
      match m.get("encoding-type") {
        Some(@src.JsonStr(s)) => assert_eq(s, "array")
        _ => fail("encoding-type missing")
      }
      match m.get("encoding-version") {
        Some(@src.JsonStr(s)) => assert_eq(s, "0.2.0")
        _ => fail("encoding-version missing")
      }
    }
    _ => fail("attrs not an object")
  }
  // attrs -> encoding round-trip
  match @src.anndata_enc_from_attrs(attrs) {
    Some(e) => assert_eq(@src.anndata_enc_name(e), "array")
    None => fail("from_attrs failed")
  }
  // version validation
  assert_eq(@src.anndata_enc_check("array", "0.2.0"), None)
  match @src.anndata_enc_check("array", "9.9.9") {
    Some(_) => ()
    None => fail("bad version should fail")
  }
  match @src.anndata_enc_check("nope", "0.1.0") {
    Some(_) => ()
    None => fail("unknown encoding should fail")
  }
}

///|
test "io_spec_sparse_validation" {
  assert_eq(@src.anndata_sparse_spec_check(3, 3, 4, 3), None)
  match @src.anndata_sparse_spec_check(2, 3, 4, 3) {
    Some(_) => ()
    None => fail("data/indices mismatch should fail")
  }
  match @src.anndata_sparse_spec_check(3, 3, 3, 3) {
    Some(_) => ()
    None => fail("indptr length should fail")
  }
}

///|
test "io_spec_h5ad_roundtrip_dense" {
  let ad = @src.AnnData::new([[1.0, 2.0], [3.0, 4.0]], ["c1", "c2"], [
      "g1", "g2",
    ])
    .set_obs("batch", ["A", "B"])
    .set_var("hv", ["true", "false"])
    .set_layer("counts", [[1.0, 0.0], [0.0, 1.0]])
    .set_obsm("X_pca", [[1.0, 2.0], [3.0, 4.0]])
    .set_uns("title", @src.JsonStr("demo"))
  let encoded = @src.anndata_h5ad_encode_string(ad)
  match @src.anndata_h5ad_decode_string(encoded) {
    Some(back) => {
      assert_eq(back.n_obs(), 2)
      assert_eq(back.n_vars(), 2)
      assert_eq(back.get_x(0, 0), 1.0)
      assert_eq(back.get_x(1, 1), 4.0)
      assert_eq(back.obs_names[0], "c1")
      assert_eq(back.var_names[1], "g2")
      match back.obs_column("batch") {
        Some(col) => assert_eq(col[1], "B")
        None => fail("obs batch missing")
      }
      match back.var_column("hv") {
        Some(col) => assert_eq(col[0], "true")
        None => fail("var hv missing")
      }
      match back.get_layer("counts") {
        Some(l) => assert_eq(l[1][0], 0.0)
        None => fail("layer missing")
      }
      match back.get_obsm("X_pca") {
        Some(m) => assert_eq(m[1][1], 4.0)
        None => fail("obsm missing")
      }
      match back.get_uns("title") {
        Some(@src.JsonStr(s)) => assert_eq(s, "demo")
        _ => fail("uns title missing")
      }
    }
    None => fail("h5ad decode failed")
  }
}

///|
test "io_spec_h5ad_roundtrip_sparse" {
  let sparse = @src.SparseArray::from_coo([2, 2], [[0, 0], [1, 1]], [5.0, 7.0]) catch {
    _ => fail("valid coords should not fail")
  }
  let ad = @src.AnnData::new_sparse(sparse, ["c1", "c2"], ["g1", "g2"])
  let encoded = @src.anndata_h5ad_encode_string(ad)
  match @src.anndata_h5ad_decode_string(encoded) {
    Some(back) => {
      assert_eq(back.is_sparse(), true)
      assert_eq(back.get_x(0, 0), 5.0)
      assert_eq(back.get_x(0, 1), 0.0)
      assert_eq(back.get_x(1, 1), 7.0)
      assert_eq(back.nnz(), 2)
    }
    None => fail("sparse h5ad decode failed")
  }
}

///|
test "io_spec_h5ad_roundtrip_raw" {
  let ad = @src.AnnData::new([[1.0, 2.0], [3.0, 4.0]], ["c1", "c2"], [
    "g1", "g2",
  ])
  let raw = @src.anndata_raw_capture(ad)
  ignore(ad.set_raw(raw))
  let encoded = @src.anndata_h5ad_encode_string(ad)
  match @src.anndata_h5ad_decode_string(encoded) {
    Some(back) =>
      match back.get_raw() {
        Some(r) => {
          assert_eq(r.n_obs(), 2)
          assert_eq(r.get_x(1, 0), 3.0)
          assert_eq(r.var_names()[0], "g1")
        }
        None => fail("raw not restored")
      }
    None => fail("h5ad decode failed")
  }
}

///|
test "io_spec_h5ad_malformed_returns_none" {
  match @src.anndata_h5ad_decode_string("not json") {
    None => ()
    Some(_) => fail("malformed h5ad should not decode")
  }
  match @src.anndata_h5ad_decode_string("[1, 2, 3]") {
    None => ()
    Some(_) => fail("non-object h5ad should not decode")
  }
}
