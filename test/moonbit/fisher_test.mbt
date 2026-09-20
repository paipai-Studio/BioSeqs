///|
/// Tests for Fisher's exact test (bedtools fisher port). Expected values are
/// hand-computed from the hypergeometric distribution and the canonical
/// [[4,0],[0,4]] "lady tasting tea" table (two-sided p = 2/70).

///|
test "fisher_exact_hand" {
  // [[4,0],[0,4]]: N=8, n1_=4, n_1=4, P(k=4)=C(4,4)C(4,0)/70 = 1/70
  // two-tailed = 2/70 ≈ 0.0285714
  let r = @src.fisher_exact(4, 0, 0, 4)
  assert_true(assert_close(r.right, 1.0 / 70.0, 1.0e-9))
  assert_true(assert_close(r.two_tail, 2.0 / 70.0, 1.0e-9))
  // [[3,0],[0,3]]: two-tailed = 2/20 = 0.1
  let r2 = @src.fisher_exact(3, 0, 0, 3)
  assert_true(assert_close(r2.two_tail, 0.1, 1.0e-9))
}

///|
test "fisher_exact_ratio" {
  // [[1,9],[11,3]] — classic Fisher's-exact example; ratio = (1/9)/(11/3)
  let r = @src.fisher_exact(1, 9, 11, 3)
  assert_true(assert_close(r.two_tail, 0.002759, 1.0e-5))
  assert_true(assert_close(r.ratio, 0.0303030303, 1.0e-6))
  // no-overlap table -> ratio finite; symmetric 1×1
  let r2 = @src.fisher_exact(1, 1, 1, 1)
  assert_true(assert_close(r2.two_tail, 1.0, 1.0e-12))
  assert_true(assert_close(r2.ratio, 1.0, 1.0e-12))
}

///|
test "bedtools_fisher_smoke" {
  // one A overlapping one B; the "neither" cell is estimated from genome size.
  let r = @src.bedtools_fisher([@src.BedInterval::new("chr1", 0, 100)], [@src.BedInterval::new("chr1", 50, 150)], 1000)
  // ratio is NaN because the A-only / B-only cells are 0
  assert_true(r.ratio != r.ratio)
}