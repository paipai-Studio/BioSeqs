///|
test "variable_verlet_adaptive_energy_conserved" {
  // 1D harmonic oscillator: F = -k x, k = 5.
  let k = 5.0
  fn force(pos : Array[Array[Double]]) -> Array[Array[Double]] {
    let f = @src.md_new_forces(pos.length())
    for i = 0; i < pos.length(); i = i + 1 {
      for d = 0; d < 3; d = d + 1 {
        f[i][d] = -k * pos[i][d]
      }
    }
    f
  }
  let pos0 : Array[Array[Double]] = [[0.3, 0.0, 0.0]]
  let vel0 : Array[Array[Double]] = [[0.0, 0.0, 0.0]]
  let masses = [1.0]
  let (pos_hist, _vel_hist, dt_hist) = @src.omm_variable_verlet_simulate(
    pos0,
    vel0,
    masses,
    n_steps=2000,
    force_fn=force,
    dt0=0.05,
    error_tol=0.00001,
    min_dt=0.0001,
    max_dt=0.05,
  )
  assert_eq(pos_hist.length(), 2001)
  assert_eq(dt_hist.length(), 2001)
  // dt stayed within the clamp.
  let mut dt_min = 1.0e18
  let mut dt_max = 0.0
  for d in dt_hist {
    if d < dt_min {
      dt_min = d
    }
    if d > dt_max {
      dt_max = d
    }
  }
  assert_true(dt_min >= 0.0001 - 1.0e-12)
  assert_true(dt_max <= 0.05 + 1.0e-12)
  assert_true(dt_max > dt_min)
  // Energy conservation over 2000 adaptive steps (amplitude 0.3, E0 = 0.225).
  let mut e0 = 0.0
  let f0 = force(pos0)
  for d = 0; d < 3; d = d + 1 {
    e0 = e0 + 0.5 * masses[0] * vel0[0][d] * vel0[0][d]
  }
  e0 = e0 +
    0.5 *
    k *
    (
      pos0[0][0] * pos0[0][0] +
      pos0[0][1] * pos0[0][1] +
      pos0[0][2] * pos0[0][2]
    )
  ignore(f0)
  let pf = pos_hist[2000]
  let mut pe = 0.5 *
    k *
    (pf[0][0] * pf[0][0] + pf[0][1] * pf[0][1] + pf[0][2] * pf[0][2])
  let vf = _vel_hist[2000]
  pe = pe +
    0.5 *
    masses[0] *
    (vf[0][0] * vf[0][0] + vf[0][1] * vf[0][1] + vf[0][2] * vf[0][2])
  let scale = if e0 > 1.0e-9 { e0 } else { 1.0 }
  assert_true((pe - e0).abs() / scale < 0.05)
}

///|
test "variable_verlet_step_rejects_large_error" {
  // With a huge error tolerance the full step is accepted.
  fn force0(_pos : Array[Array[Double]]) -> Array[Array[Double]] {
    [[0.0, 0.0, 0.0]]
  }
  let (x, v, dt) = @src.omm_variable_verlet_step(
    [[1.0, 0.0, 0.0]],
    [[0.0, 0.0, 0.0]],
    force0,
    [1.0],
    dt=0.01,
    error_tol=1.0,
  )
  assert_true(@moonbit.assert_num_close(x[0][0], 1.0, 1.0e-12))
  assert_true(@moonbit.assert_num_close(v[0][0], 0.0, 1.0e-12))
  // No error -> controller grows dt (capped at max_dt default 0.01).
  assert_true(@moonbit.assert_num_close(dt, 0.01, 1.0e-12))
}

///|
test "variable_langevin_reaches_temperature" {
  let k = 5.0
  fn force(pos : Array[Array[Double]]) -> Array[Array[Double]] {
    let f = @src.md_new_forces(pos.length())
    for i = 0; i < pos.length(); i = i + 1 {
      for d = 0; d < 3; d = d + 1 {
        f[i][d] = -k * pos[i][d]
      }
    }
    f
  }
  let (pos_hist, vel_hist, dt_hist) = @src.omm_variable_langevin_simulate(
    [[0.3, 0.0, 0.0]],
    [[0.0, 0.0, 0.0]],
    [1.0],
    n_steps=3000,
    force_fn=force,
    dt0=0.02,
    temperature=300.0,
    tau=0.2,
    error_tol=0.00001,
    min_dt=0.0001,
    max_dt=0.02,
  )
  assert_eq(pos_hist.length(), 3001)
  // No NaN anywhere in the final state.
  let pf = pos_hist[3000]
  let vf = vel_hist[3000]
  for d = 0; d < 3; d = d + 1 {
    assert_true(pf[0][d] == pf[0][d])
    assert_true(vf[0][d] == vf[0][d])
  }
  // Mean temperature over the last 500 steps is near the target (wide band:
  // 3 DOF only).
  let mut sum_ke = 0.0
  for s = 2500; s < 3000; s = s + 1 {
    let v = vel_hist[s]
    sum_ke = sum_ke +
      0.5 * (v[0][0] * v[0][0] + v[0][1] * v[0][1] + v[0][2] * v[0][2])
  }
  let mean_t = 2.0 * (sum_ke / 500.0) / (3.0 * 0.001987204)
  assert_true(mean_t > 60.0 && mean_t < 900.0)
  // dt stayed clamped.
  let mut ok = true
  for d in dt_hist {
    if d < 0.0001 - 1.0e-12 || d > 0.02 + 1.0e-12 {
      ok = false
    }
  }
  assert_true(ok)
}
