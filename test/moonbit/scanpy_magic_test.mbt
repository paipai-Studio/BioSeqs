///|
/// Tests for MAGIC diffusion imputation.

///|
test "scanpy_magic_matrix" {
  // c1 and c2 identical, c3 different: imputation must keep c1 == c2 and
  // pull c3's gene 0 up toward the neighbors
  let x = [[10.0, 0.0], [10.0, 0.0], [0.0, 10.0]]
  let y = @src.scanpy_magic_matrix(x, t=1, k=2)
  assert_eq(y.length(), 3)
  assert_eq(y[0].length(), 2)
  for g in 0..<2 {
    assert_true(@moonbit.assert_num_close(y[0][g], y[1][g], 1.0e-10))
  }
  // diffusion moved c3's gene 0 off its original 0
  assert_true(y[2][0] > 0.0)
  // and c3's gene 1 decreased from 10 toward the neighbors' 0
  assert_true(y[2][1] < 10.0)
}

///|
test "scanpy_magic_adata_layer" {
  let x = [[4.0, 0.0], [4.0, 0.0], [0.0, 4.0]]
  let adata = @src.AnnData::new(x, ["c1", "c2", "c3"], ["A", "B"])
  let out = @src.scanpy_magic(adata, t=2, k=2, layer="magic")
  // stored as layer, original x untouched
  assert_eq(adata.n_obs(), 3)
  assert_eq(adata.x[0][0], 4.0)
  match adata.get_layer("magic") {
    Some(l) => {
      assert_eq(l.length(), 3)
      assert_true(l[2][0] > 0.0)
    }
    None => assert_true(false)
  }
  assert_eq(out.length(), 3)
  // in-place mode: empty layer name replaces x
  let out2 = @src.scanpy_magic(adata, t=1, k=2)
  assert_true(out2[2][0] > 0.0)
  assert_true(adata.x[2][0] > 0.0)
}
