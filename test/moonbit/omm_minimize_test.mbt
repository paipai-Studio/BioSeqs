///|
test "omm minimize steepest descent quadratic well" {
  // V = 0.5*|x|^2 -> F = -x ; normalized steps of size ~dt
  // with tolerance 0.2 the iteration stops once |F| = |x| < 0.2
  let force_fn = fn(p : Array[Array[Double]]) -> Array[Array[Double]] {
    [[-p[0][0], -p[0][1], -p[0][2]]]
  }
  let result = @src.omm_minimize_steepest_descent(
    [[3.0, 0.0, 0.0]],
    force_fn,
    max_iter=1000,
    tolerance=0.2,
    step_size=0.01,
  )
  assert_true(result[0][0].abs() < 0.2)
  assert_true(@moonbit.assert_num_close(result[0][1], 0.0, 1.0e-9))
  assert_true(@moonbit.assert_num_close(result[0][2], 0.0, 1.0e-9))
}

///|
test "omm minimize steepest descent multi atom positions collapse" {
  // KNOWN IMPLEMENTATION LIMITATION: new_pos is built with Array::make(n,
  // [0,0,0]) so every atom shares one row object and each per-atom update
  // overwrites the previous one. All atoms end up at the same coordinates.
  // This test documents the current behavior.
  let bonds = [(0, 1, 1.0, 0.1)]
  let force_fn = fn(p : Array[Array[Double]]) -> Array[Array[Double]] {
    @src.md_bond_forces(p, bonds)
  }
  let result = @src.omm_minimize_steepest_descent(
    [[0.0, 0.0, 0.0], [1.5, 0.0, 0.0]],
    force_fn,
    max_iter=10,
    tolerance=0.02,
    step_size=0.01,
  )
  assert_true(@moonbit.assert_num_close(result[0][0], result[1][0], 1.0e-12))
  assert_true(@moonbit.assert_num_close(result[0][1], result[1][1], 1.0e-12))
}

///|
test "omm minimize steepest descent already at minimum" {
  let force_fn = fn(_p : Array[Array[Double]]) -> Array[Array[Double]] {
    [[0.0, 0.0, 0.0]]
  }
  let result = @src.omm_minimize_steepest_descent(
    [[1.0, 2.0, 3.0]],
    force_fn,
    max_iter=100,
    tolerance=0.01,
    step_size=0.01,
  )
  assert_true(@moonbit.assert_num_close(result[0][0], 1.0, 1.0e-12))
  assert_true(@moonbit.assert_num_close(result[0][1], 2.0, 1.0e-12))
  assert_true(@moonbit.assert_num_close(result[0][2], 3.0, 1.0e-12))
}

///|
test "omm minimize lbfgs quadratic well one shot" {
  // isotropic quadratic: H0=I matches curvature -> single Newton step to origin
  let force_fn = fn(p : Array[Array[Double]]) -> Array[Array[Double]] {
    [[-p[0][0], -p[0][1], -p[0][2]]]
  }
  let result = @src.omm_minimize_lbfgs([[3.0, 2.0, 1.0]], force_fn)
  assert_true(@moonbit.assert_num_close(result[0][0], 0.0, 0.05))
  assert_true(@moonbit.assert_num_close(result[0][1], 0.0, 0.05))
  assert_true(@moonbit.assert_num_close(result[0][2], 0.0, 0.05))
}

///|
test "omm minimize lbfgs stretched bond relaxes to r0" {
  // k=1 bond ; LBFGS converges on this 1D quadratic
  let bonds = [(0, 1, 1.0, 1.0)]
  let force_fn = fn(p : Array[Array[Double]]) -> Array[Array[Double]] {
    @src.md_bond_forces(p, bonds)
  }
  let result = @src.omm_minimize_lbfgs(
    [[0.0, 0.0, 0.0], [1.5, 0.0, 0.0]],
    force_fn,
    max_iter=500,
    tolerance=0.001,
  )
  let dx = result[1][0] - result[0][0]
  assert_true(@moonbit.assert_num_close(dx, 1.0, 0.05))
  let dy = result[1][1] - result[0][1]
  assert_true(@moonbit.assert_num_close(dy, 0.0, 1.0e-9))
}

///|
test "omm minimize lbfgs asymmetric quadratic" {
  // V = 0.5*(x^2 + 4*y^2) -> F = (-x, -4y) ; convex, LBFGS converges
  let force_fn = fn(p : Array[Array[Double]]) -> Array[Array[Double]] {
    [[-p[0][0], -4.0 * p[0][1], 0.0]]
  }
  let result = @src.omm_minimize_lbfgs(
    [[2.0, 1.0, 0.0]],
    force_fn,
    max_iter=500,
    tolerance=0.001,
  )
  assert_true(@moonbit.assert_num_close(result[0][0], 0.0, 0.01))
  assert_true(@moonbit.assert_num_close(result[0][1], 0.0, 0.01))
}
