///|
/// Tests for the skallel selection (EHH / iHS / nSL / standardization) port.
/// Expected values are taken from scikit-allel's `test_selection.py` cases
/// (`test_ihh_scan_a`, `test_ihh_scan_c`) and hand-verification.

///|
test "skallel_ehh_decay" {
  // identical haplotypes -> EHH stays 1
  let e1 = @src.skallel_ehh_decay([[0, 0], [0, 0]])
  assert_true(assert_close(e1[0], 1.0, 1.0e-9))
  assert_true(assert_close(e1[1], 1.0, 1.0e-9))
  // opposite haplotypes -> EHH 0 everywhere
  let e2 = @src.skallel_ehh_decay([[0, 1], [0, 1]])
  assert_true(assert_close(e2[0], 0.0, 1.0e-12))
  assert_true(assert_close(e2[1], 0.0, 1.0e-12))
}

///|
test "skallel_ihh_scan" {
  // test_ihh_scan_a: 1 pair, homozygosity over all variants
  let h1 = [[0, 0], [0, 0], [0, 0]]
  let r1 = @src.skallel_ihh_scan(h1, [10.0, 10.0], 0.0, true)
  assert_eq(r1.length(), 3)
  assert_true(assert_close(r1[0], 0.0, 1.0e-9))
  assert_true(assert_close(r1[1], 10.0, 1.0e-9))
  assert_true(assert_close(r1[2], 20.0, 1.0e-9))
  // test_ihh_scan_c: 1 pair, homozygosity decays
  let h2 = [[0, 1], [0, 0], [0, 0]]
  let r2 = @src.skallel_ihh_scan(h2, [10.0, 10.0], 0.0, false)
  assert_true(assert_close(r2[0], 0.0, 1.0e-9))
  assert_true(assert_close(r2[1], 5.0, 1.0e-9))
  assert_true(assert_close(r2[2], 15.0, 1.0e-9))
}

///|
test "skallel_nsl_ihs_symmetric" {
  // 4 variants x 6 haplotypes, first 3 all-ref, last 3 all-alt
  let h = [
    [0, 0, 0, 1, 1, 1],
    [0, 0, 0, 1, 1, 1],
    [0, 0, 0, 1, 1, 1],
    [0, 0, 0, 1, 1, 1],
  ]
  let nsl = @src.skallel_nsl(h)
  assert_eq(nsl.length(), 4)
  for i = 0; i < 4; i = i + 1 {
    assert_true(assert_close(nsl[i], 0.0, 1.0e-9))
  }
  let ihs = @src.skallel_ihs(h, [1, 2, 3, 4], 0.0, false)
  assert_eq(ihs.length(), 4)
}

///|
test "skallel_standardize" {
  let z = @src.skallel_standardize([10.0, 20.0, 30.0, 40.0])
  // mean 0, std 1
  let mut mean = 0.0
  for x in z {
    mean = mean + x
  }
  mean = mean / 4.0
  let mut variance = 0.0
  for x in z {
    let d = x - mean
    variance = variance + d * d
  }
  let std = (variance / 4.0).sqrt()
  assert_true(mean.abs() < 1.0e-12)
  assert_true(assert_close(std, 1.0, 1.0e-9))
}

///|
test "skallel_standardize_by_allele_count" {
  // three two-item bins: aac [1,2] [3,4] [5,6], each standardized to ±1
  let z = @src.skallel_standardize_by_allele_count(
    [1.0, 2.0, 3.0, 4.0, 5.0, 6.0],
    [1, 2, 3, 4, 5, 6],
    3,
  )
  assert_eq(z.length(), 6)
  for i = 0; i < 6; i = i + 1 {
    let expect = if i % 2 == 0 { -1.0 } else { 1.0 }
    assert_true(assert_close(z[i], expect, 1.0e-9))
  }
}