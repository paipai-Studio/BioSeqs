///|
/// Test the Weir–Cockerham Fst port. Expected values are taken verbatim from
/// the scikit-allel docstring example for `weir_cockerham_fst`.

///|
test "skallel_weir_cockerham_fst" {
  let g = [
    [[0, 0], [0, 0], [1, 1], [1, 1]],
    [[0, 1], [0, 1], [0, 1], [0, 1]],
    [[0, 0], [0, 0], [0, 0], [0, 0]],
    [[0, 1], [1, 2], [1, 1], [2, 2]],
    [[0, 0], [1, 1], [0, 1], [-1, -1]],
  ]
  let subpops = [[0, 1], [2, 3]]
  let (a, b, c) = @src.skallel_weir_cockerham_fst(g, subpops, 2)

  let exp_a = [
    [0.5, 0.5, 0.0],
    [0.0, 0.0, 0.0],
    [0.0, 0.0, 0.0],
    [0.0, -0.125, -0.125],
    [-0.375, -0.375, 0.0],
  ]
  let exp_b = [
    [0.0, 0.0, 0.0],
    [-0.25, -0.25, 0.0],
    [0.0, 0.0, 0.0],
    [0.0, 0.125, 0.25],
    [0.41666667, 0.41666667, 0.0],
  ]
  let exp_c = [
    [0.0, 0.0, 0.0],
    [0.5, 0.5, 0.0],
    [0.0, 0.0, 0.0],
    [0.125, 0.25, 0.125],
    [0.16666667, 0.16666667, 0.0],
  ]
  for v = 0; v < 5; v = v + 1 {
    for d = 0; d < 3; d = d + 1 {
      assert_true(assert_close(a[v][d], exp_a[v][d], 1.0e-5))
      assert_true(assert_close(b[v][d], exp_b[v][d], 1.0e-5))
      assert_true(assert_close(c[v][d], exp_c[v][d], 1.0e-5))
    }
  }

  // per-variant Fst (averaging over alleles); nan -> 0 for the monomorphic variant
  let fst = @src.skallel_weir_cockerham_fst_variant(g, subpops, 2)
  let exp_fst = [1.0, 0.0, 0.0, -0.4, -1.8]
  for v = 0; v < 5; v = v + 1 {
    assert_true(assert_close(fst[v], exp_fst[v], 1.0e-5))
  }

  // block estimator: numerator sums to 0 here, so Fst ≈ 0
  let block = @src.skallel_weir_cockerham_fst_block(g, subpops, 2)
  assert_true(block.abs() < 1.0e-5)
}