///|
test "omm custom integrator explicit euler free particle" {
  let integ = @src.OmmCustomIntegrator::new(0.1)
  integ.add_step(variable="v", expression="v + dt * f / m")
  integ.add_step(variable="x", expression="x + dt * v")
  let (pos, vel) = integ.step(
    [[0.0, 0.0, 0.0]],
    [[2.0, 0.0, 0.0]],
    [[0.0, 0.0, 0.0]],
    [1.0],
  )
  assert_true(@moonbit.assert_num_close(vel[0][0], 2.0, 1.0e-12))
  assert_true(@moonbit.assert_num_close(pos[0][0], 0.2, 1.0e-12))
}

///|
test "omm custom integrator constant force updates velocity first" {
  // semi-implicit Euler: v_new = v + dt*f/m, then x_new = x + dt*v_new
  let integ = @src.OmmCustomIntegrator::new(0.1)
  integ.add_step(variable="v", expression="v + dt * f / m")
  integ.add_step(variable="x", expression="x + dt * v")
  let (pos, vel) = integ.step(
    [[0.0, 0.0, 0.0]],
    [[0.0, 0.0, 0.0]],
    [[1.0, 0.0, 0.0]],
    [1.0],
  )
  assert_true(@moonbit.assert_num_close(vel[0][0], 0.1, 1.0e-12))
  assert_true(@moonbit.assert_num_close(pos[0][0], 0.01, 1.0e-12))
}

///|
test "omm custom integrator harmonic oscillator semi implicit step" {
  // force = -x ; x=1, v=0, dt=0.1 -> v=-0.1, x = 1 + 0.1*(-0.1) = 0.99
  let integ = @src.OmmCustomIntegrator::new(0.1)
  integ.add_step(variable="v", expression="v - dt * x * k / m")
  integ.add_step(variable="x", expression="x + dt * v")
  integ.add_global_variable("k", 1.0)
  let (pos, vel) = integ.step(
    [[1.0, 0.0, 0.0]],
    [[0.0, 0.0, 0.0]],
    [[0.0, 0.0, 0.0]],
    [1.0],
  )
  assert_true(@moonbit.assert_num_close(vel[0][0], -0.1, 1.0e-12))
  assert_true(@moonbit.assert_num_close(pos[0][0], 0.99, 1.0e-12))
}

///|
test "omm custom integrator global variable damping" {
  let integ = @src.OmmCustomIntegrator::new(0.01)
  integ.add_global_variable("damping", 0.5)
  integ.add_step(variable="v", expression="v * damping")
  let (_pos, vel) = integ.step(
    [[0.0, 0.0, 0.0]],
    [[4.0, 0.0, 0.0]],
    [[0.0, 0.0, 0.0]],
    [1.0],
  )
  assert_true(@moonbit.assert_num_close(vel[0][0], 2.0, 1.0e-12))
}

///|
test "omm custom integrator expression precedence and parentheses" {
  // x + 2.0 * (v + 1.0) with x=1, v=2 -> 1 + 6 = 7
  let integ = @src.OmmCustomIntegrator::new(0.01)
  integ.add_step(variable="x", expression="x + 2.0 * ( v + 1.0 )")
  let (pos, _vel) = integ.step(
    [[1.0, 0.0, 0.0]],
    [[2.0, 0.0, 0.0]],
    [[0.0, 0.0, 0.0]],
    [1.0],
  )
  assert_true(@moonbit.assert_num_close(pos[0][0], 7.0, 1.0e-12))
}

///|
test "omm custom integrator division by zero guarded" {
  let integ = @src.OmmCustomIntegrator::new(0.01)
  integ.add_step(variable="x", expression="x / 0.0")
  let (pos, _vel) = integ.step(
    [[5.0, 0.0, 0.0]],
    [[0.0, 0.0, 0.0]],
    [[0.0, 0.0, 0.0]],
    [1.0],
  )
  assert_true(@moonbit.assert_num_close(pos[0][0], 0.0, 1.0e-12))
}

///|
test "omm custom integrator updates all particles independently" {
  let integ = @src.OmmCustomIntegrator::new(0.1)
  integ.add_step(variable="x", expression="x + dt * v")
  let (pos, _vel) = integ.step(
    [[1.0, 0.0, 0.0], [2.0, 0.0, 0.0]],
    [[1.0, 0.0, 0.0], [-2.0, 0.0, 0.0]],
    [[0.0, 0.0, 0.0], [0.0, 0.0, 0.0]],
    [1.0, 1.0],
  )
  assert_true(@moonbit.assert_num_close(pos[0][0], 1.1, 1.0e-12))
  assert_true(@moonbit.assert_num_close(pos[1][0], 1.8, 1.0e-12))
}

///|
test "omm custom integrator mass scaling" {
  // f=2, m=4 -> a = 0.5 ; v = 0 + 0.1*2/4 = 0.05
  let integ = @src.OmmCustomIntegrator::new(0.1)
  integ.add_step(variable="v", expression="v + dt * f / m")
  let (_pos, vel) = integ.step(
    [[0.0, 0.0, 0.0]],
    [[0.0, 0.0, 0.0]],
    [[2.0, 0.0, 0.0]],
    [4.0],
  )
  assert_true(@moonbit.assert_num_close(vel[0][0], 0.05, 1.0e-12))
}
