///|
/// Tests for the skallel cross-population (XPEHH / XPNSL) and PBS statistics.

///|
test "skallel_xpehh_xpnsl_symmetric" {
  // 4 variants x 6 haplotypes, first 3 all-ref, last 3 all-alt
  let h = [
    [0, 0, 0, 1, 1, 1],
    [0, 0, 0, 1, 1, 1],
    [0, 0, 0, 1, 1, 1],
    [0, 0, 0, 1, 1, 1],
  ]
  // same population -> log(1) = 0
  let xpe = @src.skallel_xpehh(h, h, [1, 2, 3, 4], 0.0, false)
  assert_eq(xpe.length(), 4)
  for i = 0; i < 4; i = i + 1 {
    assert_true(assert_close(xpe[i], 0.0, 1.0e-9))
  }
  let xpn = @src.skallel_xpnsl(h, h)
  assert_eq(xpn.length(), 4)
  for i = 0; i < 4; i = i + 1 {
    assert_true(assert_close(xpn[i], 0.0, 1.0e-9))
  }
}

///|
test "skallel_pbs" {
  // identity -> Fst 0 -> PBS 0
  let ac = [[2, 2], [2, 2]]
  let p0 = @src.skallel_pbs(ac, ac, ac, 1, 1, false)
  assert_eq(p0.length(), 2)
  for i = 0; i < 2; i = i + 1 {
    assert_true(assert_close(p0[i], 0.0, 1.0e-9))
  }

  // pop1 (50/50) vs pop2/pop3 (all-alt): Fst(1,·)=1/3, Fst(2,3)=0
  let ac1 = [[2, 2], [2, 2]]
  let ac2 = [[0, 4], [0, 4]]
  let ac3 = [[0, 4], [0, 4]]
  let ln15 = @math.ln(1.5)
  let un = @src.skallel_pbs(ac1, ac2, ac3, 2, 2, false)
  assert_eq(un.length(), 1)
  assert_true(assert_close(un[0], ln15, 1.0e-6))
  let nm = @src.skallel_pbs(ac1, ac2, ac3, 2, 2, true)
  assert_true(assert_close(nm[0], ln15 / (1.0 + ln15), 1.0e-6))
}