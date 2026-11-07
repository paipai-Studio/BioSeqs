///|
/// Tests for mda_dihedral_acf.mbt (dihedral autocorrelation functions).

///|
test "mda_dihedral_acf_cos_constant" {
  // Constant dihedral: cos(phi(0)-phi(t)) = 1 for every lag.
  let acf = @src.mda_dihedral_autocorrelation_cos([0.5, 0.5, 0.5])
  assert_eq(acf.length(), 3)
  assert_true(@moonbit.assert_num_close(acf[0], 1.0, 1.0e-9))
  assert_true(@moonbit.assert_num_close(acf[1], 1.0, 1.0e-9))
  assert_true(@moonbit.assert_num_close(acf[2], 1.0, 1.0e-9))
}

///|
test "mda_dihedral_acf_cos_alternating" {
  // phi alternates 0, pi, 0, pi -> C(t) = cos(t*pi) = 1, -1, 1, -1.
  let pi = 3.141592653589793
  let acf = @src.mda_dihedral_autocorrelation_cos([0.0, pi, 0.0, pi])
  assert_eq(acf.length(), 4)
  assert_true(@moonbit.assert_num_close(acf[0], 1.0, 1.0e-9))
  assert_true(@moonbit.assert_num_close(acf[1], -1.0, 1.0e-9))
  assert_true(@moonbit.assert_num_close(acf[2], 1.0, 1.0e-9))
  assert_true(@moonbit.assert_num_close(acf[3], -1.0, 1.0e-9))
}

///|
test "mda_dihedral_acf_cos_empty" {
  let acf = @src.mda_dihedral_autocorrelation_cos([])
  assert_eq(acf.length(), 0)
}

///|
test "mda_dihedral_acf_product_linear_ramp" {
  // phi = 0,1,2,3; mean = 1.5; dev = -1.5,-0.5,0.5,1.5; var = 5.
  // C0 = 5/5 = 1; C1 = 1.25/5 = 0.25; C2 = -1.5/5 = -0.3; C3 = -2.25/5 = -0.45.
  let acf = @src.mda_dihedral_autocorrelation_product([0.0, 1.0, 2.0, 3.0])
  assert_eq(acf.length(), 4)
  assert_true(@moonbit.assert_num_close(acf[0], 1.0, 1.0e-9))
  assert_true(@moonbit.assert_num_close(acf[1], 0.25, 1.0e-9))
  assert_true(@moonbit.assert_num_close(acf[2], -0.3, 1.0e-9))
  assert_true(@moonbit.assert_num_close(acf[3], -0.45, 1.0e-9))
}

///|
test "mda_dihedral_acf_product_zero_variance" {
  // Constant series has zero variance -> all zeros.
  let acf = @src.mda_dihedral_autocorrelation_product([1.0, 1.0, 1.0])
  assert_eq(acf.length(), 3)
  assert_true(@moonbit.assert_num_close(acf[0], 0.0, 1.0e-12))
  assert_true(@moonbit.assert_num_close(acf[2], 0.0, 1.0e-12))
}

///|
test "mda_dihedral_acf_ensemble_average" {
  // res0 alternates 0,pi (ACF 1,-1,1,-1); res1 is constant (ACF all 1).
  // Ensemble average = 1, 0, 1, 0.
  let pi = 3.141592653589793
  let acf = @src.mda_dihedral_autocorrelation_ensemble([
    [0.0, pi, 0.0, pi],
    [0.3, 0.3, 0.3, 0.3],
  ])
  assert_eq(acf.length(), 4)
  assert_true(@moonbit.assert_num_close(acf[0], 1.0, 1.0e-9))
  assert_true(@moonbit.assert_num_close(acf[1], 0.0, 1.0e-9))
  assert_true(@moonbit.assert_num_close(acf[2], 1.0, 1.0e-9))
  assert_true(@moonbit.assert_num_close(acf[3], 0.0, 1.0e-9))
}

///|
test "mda_dihedral_acf_ensemble_empty" {
  assert_eq(@src.mda_dihedral_autocorrelation_ensemble([]).length(), 0)
}

///|
test "mda_correlation_time_piecewise" {
  // ACF 1,0.9,0.8,0.7,0.5,0.3,0.1; threshold = 1/e ~= 0.368.
  // Trapezoid over lags 0..4 (acf[5] < threshold stops): 0.95+0.85+0.75+0.6+0.4 = 3.55.
  let acf = [1.0, 0.9, 0.8, 0.7, 0.5, 0.3, 0.1]
  assert_true(
    @moonbit.assert_num_close(@src.mda_correlation_time(acf), 3.55, 1.0e-9),
  )
}

///|
test "mda_correlation_time_too_short" {
  assert_true(
    @moonbit.assert_num_close(@src.mda_correlation_time([1.0]), 0.0, 1.0e-12),
  )
}

///|
test "mda_relaxation_time_exponential_fit" {
  // C(t) = exp(-t/2) -> ln C linear in t with slope -0.5 -> tau = 2.
  // n = 8 -> max_fit = 4 -> regression over t = 1, 2, 3.
  let acf = [
    @math.exp(-0.0),
    @math.exp(-0.5),
    @math.exp(-1.0),
    @math.exp(-1.5),
    @math.exp(-2.0),
    @math.exp(-2.5),
    @math.exp(-3.0),
    @math.exp(-3.5),
  ]
  assert_true(
    @moonbit.assert_num_close(@src.mda_relaxation_time(acf), 2.0, 1.0e-6),
  )
}

///|
test "mda_relaxation_time_growing_acf" {
  // Increasing ACF gives a non-negative slope -> function returns 0.
  let acf = [1.0, 1.1, 1.2, 0.5, 0.5, 0.5]
  assert_true(
    @moonbit.assert_num_close(@src.mda_relaxation_time(acf), 0.0, 1.0e-12),
  )
}

///|
test "mda_relaxation_time_too_short" {
  assert_true(
    @moonbit.assert_num_close(
      @src.mda_relaxation_time([0.5, 0.4]),
      0.0,
      1.0e-12,
    ),
  )
}
