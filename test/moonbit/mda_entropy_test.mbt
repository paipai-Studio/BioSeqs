///|
/// Tests for mda_entropy.mbt (quasi-harmonic / Schlitter entropy).

///|
test "mda_quasi_harmonic_entropy_single_mode" {
  // S/kB = 0.5 * ln(2*pi*kB*T / lambda)
  let s = @src.mda_quasi_harmonic_entropy([1.0])
  let expected = 0.5 * @math.ln(2.0 * 3.141592653589793 * 300.0 / 1.0)
  assert_true(@moonbit.assert_num_close(s, expected, 1.0e-12))
}

///|
test "mda_quasi_harmonic_entropy_temperature_and_kB" {
  // S/kB = 0.5 * ln(2*pi*kB*T / lambda) with T = 100, kB = 2.
  let s = @src.mda_quasi_harmonic_entropy([4.0], temperature=100.0, kB=2.0)
  let expected = 0.5 * @math.ln(2.0 * 3.141592653589793 * 200.0 / 4.0)
  assert_true(@moonbit.assert_num_close(s, expected, 1.0e-12))
}

///|
test "mda_quasi_harmonic_entropy_skips_nonpositive" {
  // Non-positive eigenvalues are skipped: only lambda=2 contributes.
  let s = @src.mda_quasi_harmonic_entropy([0.0, -1.0, 2.0])
  let expected = 0.5 * @math.ln(2.0 * 3.141592653589793 * 300.0 / 2.0)
  assert_true(@moonbit.assert_num_close(s, expected, 1.0e-12))
}

///|
test "mda_quasi_harmonic_entropy_empty" {
  assert_true(
    @moonbit.assert_num_close(@src.mda_quasi_harmonic_entropy([]), 0.0, 1.0e-12),
  )
}

///|
test "mda_trajectory_entropy_static_trajectory" {
  // Identical frames -> all variances 0 -> entropy 0.
  let traj = [
    [[0.0, 0.0, 0.0], [1.0, 1.0, 1.0]],
    [[0.0, 0.0, 0.0], [1.0, 1.0, 1.0]],
  ]
  assert_true(
    @moonbit.assert_num_close(@src.mda_trajectory_entropy(traj), 0.0, 1.0e-12),
  )
}

///|
test "mda_trajectory_entropy_short" {
  // Fewer than 2 frames -> 0.
  assert_true(
    @moonbit.assert_num_close(
      @src.mda_trajectory_entropy([[[0.0, 0.0, 0.0]]]),
      0.0,
      1.0e-12,
    ),
  )
}

///|
test "mda_schlitter_entropy_single_mode" {
  // S = 0.5 + ln(sqrt(2*pi*e / lambda))
  let s = @src.mda_schlitter_entropy([1.0])
  let expected = 0.5 +
    @math.ln((2.0 * 3.141592653589793 * 2.718281828459045 / 1.0).sqrt())
  assert_true(@moonbit.assert_num_close(s, expected, 1.0e-12))
}

///|
test "mda_schlitter_entropy_additive_modes" {
  // Two modes contribute additively; non-positive eigenvalues are skipped.
  let s = @src.mda_schlitter_entropy([1.0, -3.0, 4.0], kB=2.0)
  let expected = 2.0 *
    (
      0.5 +
      @math.ln((2.0 * 3.141592653589793 * 2.718281828459045 / 1.0).sqrt()) +
      0.5 +
      @math.ln((2.0 * 3.141592653589793 * 2.718281828459045 / 4.0).sqrt())
    )
  assert_true(@moonbit.assert_num_close(s, expected, 1.0e-12))
}

///|
test "mda_entropy_difference" {
  // S([2]) - S([1]) = 0.5 * ln(2*pi*T/2) - 0.5 * ln(2*pi*T) = -0.5*ln(2).
  let d = @src.mda_entropy_difference([1.0], [2.0])
  assert_true(@moonbit.assert_num_close(d, -0.5 * @math.ln(2.0), 1.0e-12))
}

///|
test "mda_entropy_free_energy" {
  assert_true(
    @moonbit.assert_num_close(
      @src.mda_entropy_free_energy(2.0, temperature=300.0),
      -600.0,
      1.0e-12,
    ),
  )
  assert_true(
    @moonbit.assert_num_close(
      @src.mda_entropy_free_energy(-1.5),
      450.0,
      1.0e-12,
    ),
  )
}
