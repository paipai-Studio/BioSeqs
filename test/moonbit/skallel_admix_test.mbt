///|
/// Tests for the skallel admixture / selection / windowed ports. Values are
/// hand-computed from the scikit-allel formulas.

///|
test "skallel_patterson" {
  // h_hat: [2,2] -> 2*2/(4*3) = 1/3 ; [4,0] -> 0
  let hh = @src.skallel_h_hat([[2, 2], [4, 0]])
  assert_true(assert_close(hh[0], 1.0 / 3.0, 1.0e-6))
  assert_true(assert_close(hh[1], 0.0, 1.0e-12))

  // f2: (a-b)^2 - ha/sa - hb/sb = 0.25 - 1/12
  let f2 = @src.skallel_patterson_f2([[2, 2], [2, 2]], [[0, 4], [0, 4]])
  for i = 0; i < 2; i = i + 1 {
    assert_true(assert_close(f2[i], 0.25 - 1.0 / 12.0, 1.0e-6))
  }

  // f3: T = (c-a)(c-b) - hc/sc ; B = 2*hc
  let (t, b) = @src.skallel_patterson_f3(
    [[2, 2], [2, 2]],
    [[0, 4], [0, 4]],
    [[0, 4], [0, 4]],
  )
  for i = 0; i < 2; i = i + 1 {
    assert_true(assert_close(t[i], 0.25 - 1.0 / 12.0, 1.0e-6))
    assert_true(assert_close(b[i], 2.0 / 3.0, 1.0e-6))
  }

  // D: num=(a-b)(c-d), den=(a+b-2ab)(c+d-2cd); average = -1
  let (num, den) = @src.skallel_patterson_d(
    [[2, 2], [2, 2]],
    [[0, 4], [0, 4]],
    [[0, 4], [0, 4]],
    [[4, 0], [4, 0]],
  )
  for i = 0; i < 2; i = i + 1 {
    assert_true(assert_close(num[i], -0.5, 1.0e-6))
    assert_true(assert_close(den[i], 0.5, 1.0e-6))
  }
  assert_true(assert_close(
    @src.skallel_patterson_d_average(
      [[2, 2], [2, 2]],
      [[0, 4], [0, 4]],
      [[0, 4], [0, 4]],
      [[4, 0], [4, 0]],
    ),
    -1.0,
    1.0e-6,
  ))
}

///|
test "skallel_garud_h" {
  let h = [[0, 0, 0, 1], [0, 0, 1, 1]]
  // distinct frequencies [0.5, 0.25, 0.25]
  let hd = @src.skallel_haplotype_diversity(h)
  assert_true(assert_close(hd, 0.625 * 4.0 / 3.0, 1.0e-6))
  let (h1, h12, h123, h2_h1) = @src.skallel_garud_h(h)
  assert_true(assert_close(h1, 0.375, 1.0e-6))
  assert_true(assert_close(h12, 0.625, 1.0e-6))
  assert_true(assert_close(h123, 1.0, 1.0e-6))
  assert_true(assert_close(h2_h1, 1.0 / 3.0, 1.0e-6))
}

///|
test "skallel_windowed_tajima_d" {
  let win = @src.skallel_index_windows(4, 2, 2)
  assert_eq(win.length(), 2)
  // monomorphic variants -> D = 0 in every window
  let d = @src.skallel_windowed_tajima_d(
    [[4, 0], [4, 0], [4, 0], [4, 0]],
    4,
    2,
    2,
  )
  assert_eq(d.length(), 2)
  for i = 0; i < 2; i = i + 1 {
    assert_true(assert_close(d[i], 0.0, 1.0e-12))
  }
}