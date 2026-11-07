///|
test "omm leapfrog step free particle uniform motion" {
  let positions = [[0.0, 0.0, 0.0]]
  let velocities_half = [[2.0, 0.0, 0.0]]
  let forces = [[0.0, 0.0, 0.0]]
  let masses = [1.0]
  let (x_new, v_new) = @src.omm_leapfrog_step(
    positions, velocities_half, forces, masses, 0.1,
  )
  assert_true(@moonbit.assert_num_close(x_new[0][0], 0.2, 1.0e-12))
  assert_true(@moonbit.assert_num_close(v_new[0][0], 2.0, 1.0e-12))
}

///|
test "omm leapfrog step constant force updates half step velocity" {
  let positions = [[0.0, 0.0, 0.0]]
  let velocities_half = [[0.0, 0.0, 0.0]]
  let forces = [[1.0, 0.0, 0.0]]
  let masses = [1.0]
  let (x_new, v_new) = @src.omm_leapfrog_step(
    positions, velocities_half, forces, masses, 0.1,
  )
  // v(dt/2) = 0 + 0.1*1/1 = 0.1 ; x(dt) = 0 + 0.1*0.1 = 0.01
  assert_true(@moonbit.assert_num_close(v_new[0][0], 0.1, 1.0e-12))
  assert_true(@moonbit.assert_num_close(x_new[0][0], 0.01, 1.0e-12))
}

///|
test "omm leapfrog to full velocity half kick" {
  let velocities_half = [[2.0, 0.0, 0.0]]
  let forces = [[1.0, 0.0, 0.0]]
  let masses = [1.0]
  let v_full = @src.omm_leapfrog_to_full_velocity(
    velocities_half, forces, masses, 0.1,
  )
  assert_true(@moonbit.assert_num_close(v_full[0][0], 2.05, 1.0e-12))
}

///|
test "omm leapfrog simulate free particle linear trajectory" {
  let force_fn = fn(_p : Array[Array[Double]]) -> Array[Array[Double]] {
    [[0.0, 0.0, 0.0]]
  }
  let (pos_hist, vel_hist) = @src.omm_leapfrog_simulate(
    [[0.0, 0.0, 0.0]],
    [[2.0, 0.0, 0.0]],
    [1.0],
    0.1,
    3,
    force_fn,
  )
  assert_eq(pos_hist.length(), 4)
  assert_eq(vel_hist.length(), 4)
  assert_true(@moonbit.assert_num_close(pos_hist[1][0][0], 0.2, 1.0e-12))
  assert_true(@moonbit.assert_num_close(pos_hist[2][0][0], 0.4, 1.0e-12))
  assert_true(@moonbit.assert_num_close(pos_hist[3][0][0], 0.6, 1.0e-12))
  for i = 0; i < 4; i = i + 1 {
    assert_true(@moonbit.assert_num_close(vel_hist[i][0][0], 2.0, 1.0e-12))
  }
}

///|
test "omm leapfrog simulate harmonic oscillator approximates energy" {
  // k=1, m=1, x0=1, v0=0, dt=0.001, 200 steps -> t=0.2
  // analytic: x = cos(0.2), E = 0.5 ; leapfrog conserves E within dt^2
  let force_fn = fn(p : Array[Array[Double]]) -> Array[Array[Double]] {
    [[-p[0][0], 0.0, 0.0]]
  }
  let (pos_hist, vel_hist) = @src.omm_leapfrog_simulate(
    [[1.0, 0.0, 0.0]],
    [[0.0, 0.0, 0.0]],
    [1.0],
    0.001,
    200,
    force_fn,
  )
  assert_eq(pos_hist.length(), 201)
  let x_final = pos_hist[200][0][0]
  assert_true(@moonbit.assert_num_close(x_final, @math.cos(0.2), 5.0e-3))
  // reconstruct full-step velocity and check total energy
  let f_final = force_fn(pos_hist[200])
  let v_full = @src.omm_leapfrog_to_full_velocity(
    vel_hist[200],
    f_final,
    [1.0],
    0.001,
  )
  let energy = 0.5 * x_final * x_final + 0.5 * v_full[0][0] * v_full[0][0]
  assert_true(@moonbit.assert_num_close(energy, 0.5, 1.0e-2))
}
