///|
test "omm_respa_harmonic_matches_velocity_verlet" {
  // fast forces zero: r-RESPA degenerates to standard velocity Verlet
  // analytic: x1 = x0 + dt*v0 - 0.5*dt^2*k*x0, v1 = v0 - dt*k*(x0+x1)/2
  let k = 1.0
  let masses = [1.0]
  let pos0 = [[1.0, 0.0, 0.0]]
  let vel0 = [[0.0, 0.0, 0.0]]
  let fast = fn(p : Array[Array[Double]]) -> Array[Array[Double]] {
    let f : Array[Array[Double]] = Array(capacity=p.length())
    for _i = 0; _i < p.length(); _i = _i + 1 {
      f.push([0.0, 0.0, 0.0])
    }
    f
  }
  let slow = fn(p : Array[Array[Double]]) -> Array[Array[Double]] {
    let f : Array[Array[Double]] = Array(capacity=p.length())
    for i = 0; i < p.length(); i = i + 1 {
      f.push([-k * p[i][0], 0.0, 0.0])
    }
    f
  }
  let dt = 0.01
  let (pos, vel) = @src.omm_respa_step(pos0, vel0, fast, slow, masses, dt, 1)
  assert_true(@moonbit.assert_num_close(pos[0][0], 0.99995, 1.0e-12))
  assert_true(@moonbit.assert_num_close(vel[0][0], -0.00999975, 1.0e-12))
}

///|
test "omm_respa_constant_fast_force_exact" {
  // constant force is integrated exactly: x = 0.5*a*dt^2, v = a*dt
  let masses = [1.0]
  let pos0 = [[0.0, 0.0, 0.0]]
  let vel0 = [[0.0, 0.0, 0.0]]
  let fast = fn(_p : Array[Array[Double]]) -> Array[Array[Double]] {
    let f : Array[Array[Double]] = Array(capacity=1)
    f.push([2.0, 0.0, 0.0])
    f
  }
  let slow = fn(_p : Array[Array[Double]]) -> Array[Array[Double]] {
    let f : Array[Array[Double]] = Array(capacity=1)
    f.push([0.0, 0.0, 0.0])
    f
  }
  let dt = 0.01
  let (pos, vel) = @src.omm_respa_step(pos0, vel0, fast, slow, masses, dt, 2)
  assert_true(
    @moonbit.assert_num_close(pos[0][0], 0.5 * 2.0 * dt * dt, 1.0e-12),
  )
  assert_true(@moonbit.assert_num_close(vel[0][0], 2.0 * dt, 1.0e-12))
  assert_true(@moonbit.assert_num_close(pos[0][1], 0.0, 1.0e-12))
  assert_true(@moonbit.assert_num_close(vel[0][2], 0.0, 1.0e-12))
}

///|
test "omm_respa_simulate_conserves_harmonic_energy" {
  let k = 1.0
  let masses = [1.0]
  let fast = fn(_p : Array[Array[Double]]) -> Array[Array[Double]] {
    let f : Array[Array[Double]] = Array(capacity=1)
    f.push([0.0, 0.0, 0.0])
    f
  }
  let slow = fn(p : Array[Array[Double]]) -> Array[Array[Double]] {
    let f : Array[Array[Double]] = Array(capacity=p.length())
    for i = 0; i < p.length(); i = i + 1 {
      f.push([-k * p[i][0], 0.0, 0.0])
    }
    f
  }
  let pos0 = [[1.0, 0.0, 0.0]]
  let vel0 = [[0.0, 0.0, 0.0]]
  let history = @src.omm_respa_simulate(
    pos0,
    vel0,
    fast,
    slow,
    masses,
    0.01,
    50,
    n_inner=2,
  )
  assert_eq(history.length(), 51)
  // first frame is the initial configuration
  assert_true(@moonbit.assert_num_close(history[0][0][0], 1.0, 1.0e-12))
  // step manually and check energy conservation
  let mut p = pos0
  let mut v = vel0
  for _i = 0; _i < 50; _i = _i + 1 {
    let (p2, v2) = @src.omm_respa_step(p, v, fast, slow, masses, 0.01, 2)
    p = p2
    v = v2
  }
  let energy = 0.5 * v[0][0] * v[0][0] + 0.5 * k * p[0][0] * p[0][0]
  assert_true((energy - 0.5).abs() < 1.0e-4)
}
