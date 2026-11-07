///|
/// Tests for mda_diffusion.mbt (diffusion coefficient, VACF, Green-Kubo).

///|
test "mda_diffusion_coefficient_linear_msd" {
  // MSD = 6*t (t = 0..7, dt = 1): slope = 6 -> D = 6/(2*3) = 1.
  let msd = [0.0, 6.0, 12.0, 18.0, 24.0, 30.0, 36.0, 42.0]
  assert_true(
    @moonbit.assert_num_close(@src.mda_diffusion_coefficient(msd), 1.0, 1.0e-9),
  )
}

///|
test "mda_diffusion_coefficient_dimensionality" {
  // Same MSD with dimension=2 -> D = 6/4 = 1.5.
  let msd = [0.0, 6.0, 12.0, 18.0, 24.0, 30.0, 36.0, 42.0]
  assert_true(
    @moonbit.assert_num_close(
      @src.mda_diffusion_coefficient(msd, dimension=2),
      1.5,
      1.0e-9,
    ),
  )
}

///|
test "mda_diffusion_coefficient_dt_scaling" {
  // MSD = 6*i for lag i with dt = 2 -> MSD = 3*t -> slope 3 -> D = 0.5.
  let msd = [0.0, 6.0, 12.0, 18.0, 24.0, 30.0, 36.0, 42.0]
  assert_true(
    @moonbit.assert_num_close(
      @src.mda_diffusion_coefficient(msd, dt=2.0),
      0.5,
      1.0e-9,
    ),
  )
}

///|
test "mda_diffusion_coefficient_degenerate_inputs" {
  assert_true(
    @moonbit.assert_num_close(@src.mda_diffusion_coefficient([]), 0.0, 1.0e-12),
  )
  assert_true(
    @moonbit.assert_num_close(
      @src.mda_diffusion_coefficient([1.0]),
      0.0,
      1.0e-12,
    ),
  )
  // Zero time base -> denominator 0 -> 0.
  assert_true(
    @moonbit.assert_num_close(
      @src.mda_diffusion_coefficient([0.0, 0.0], dt=0.0),
      0.0,
      1.0e-12,
    ),
  )
}

///|
test "mda_diffusion_from_trajectory_ballistic" {
  // Single atom moving with x(t) = t: MSD = t^2 for lags 0..3.
  // Regression on the second half (lags 1..3) gives slope 4 -> D = 4/6.
  let traj = [
    [[0.0, 0.0, 0.0]],
    [[1.0, 0.0, 0.0]],
    [[2.0, 0.0, 0.0]],
    [[3.0, 0.0, 0.0]],
  ]
  assert_true(
    @moonbit.assert_num_close(
      @src.mda_diffusion_from_trajectory(traj),
      4.0 / 6.0,
      1.0e-9,
    ),
  )
}

///|
test "mda_velocity_autocorrelation_constant" {
  // Constant velocity (1,0,0): C(t) = 1 for every lag.
  let v = [
    [[1.0, 0.0, 0.0]],
    [[1.0, 0.0, 0.0]],
    [[1.0, 0.0, 0.0]],
    [[1.0, 0.0, 0.0]],
  ]
  let vacf = @src.mda_velocity_autocorrelation(v)
  assert_eq(vacf.length(), 4)
  for c in vacf {
    assert_true(@moonbit.assert_num_close(c, 1.0, 1.0e-9))
  }
}

///|
test "mda_velocity_autocorrelation_alternating" {
  // v alternates +1/-1 along x: C(t) = (-1)^t.
  let v = [
    [[1.0, 0.0, 0.0]],
    [[-1.0, 0.0, 0.0]],
    [[1.0, 0.0, 0.0]],
    [[-1.0, 0.0, 0.0]],
  ]
  let vacf = @src.mda_velocity_autocorrelation(v)
  assert_true(@moonbit.assert_num_close(vacf[0], 1.0, 1.0e-9))
  assert_true(@moonbit.assert_num_close(vacf[1], -1.0, 1.0e-9))
  assert_true(@moonbit.assert_num_close(vacf[2], 1.0, 1.0e-9))
  assert_true(@moonbit.assert_num_close(vacf[3], -1.0, 1.0e-9))
}

///|
test "mda_velocity_autocorrelation_zero_and_empty" {
  // All-zero velocities (vacf0 = 0) -> all zeros; empty input -> [].
  let v = [[[0.0, 0.0, 0.0]], [[0.0, 0.0, 0.0]]]
  let vacf = @src.mda_velocity_autocorrelation(v)
  assert_eq(vacf.length(), 2)
  assert_true(@moonbit.assert_num_close(vacf[1], 0.0, 1.0e-12))
  assert_eq(@src.mda_velocity_autocorrelation([]).length(), 0)
}

///|
test "mda_green_kubo_diffusion_trapezoid" {
  // Integral of [1, 0.5, 0] with dt=1: 0.75 + 0.25 = 1 -> D = 1/3.
  let d = @src.mda_green_kubo_diffusion([1.0, 0.5, 0.0])
  assert_true(@moonbit.assert_num_close(d, 1.0 / 3.0, 1.0e-12))
  // Constant unit VACF over 3 points: integral = 2 -> D = 2/3.
  let d2 = @src.mda_green_kubo_diffusion([1.0, 1.0, 1.0])
  assert_true(@moonbit.assert_num_close(d2, 2.0 / 3.0, 1.0e-12))
  assert_true(
    @moonbit.assert_num_close(
      @src.mda_green_kubo_diffusion([1.0]),
      0.0,
      1.0e-12,
    ),
  )
}

///|
test "mda_stokes_einstein_radius" {
  // r = kT / (6*pi*eta*D), evaluated symbolically for D = 1e-9.
  let r = @src.mda_stokes_einstein_radius(1.0e-9)
  let expected = 1.380649e-23 *
    300.0 /
    (6.0 * 3.141592653589793 * 0.001 * 1.0e-9)
  assert_true(@moonbit.assert_num_close(r, expected, 1.0e-6))
  // Non-positive diffusion -> 0.
  assert_true(
    @moonbit.assert_num_close(
      @src.mda_stokes_einstein_radius(0.0),
      0.0,
      1.0e-30,
    ),
  )
}
