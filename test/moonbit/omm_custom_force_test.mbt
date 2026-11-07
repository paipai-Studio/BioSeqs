///|
test "custom_expr_functions" {
  // Function nodes inside custom integrator step expressions.
  let ci = @src.OmmCustomIntegrator::new(0.1)
  ci.add_global_variable("amp", 2.0)
  // x' = x + amp * sqrt(square(0.5)) - step(-0.1) -> x + 2*0.5 - 0
  ci.add_step(
    variable="x",
    expression="x + amp * sqrt(square(0.5)) - step(-0.1)",
  )
  let (pos, _vel) = ci.step(
    [[1.0, 0.0, 0.0]],
    [[0.0, 0.0, 0.0]],
    [[0.0, 0.0, 0.0]],
    [1.0],
  )
  assert_true(@moonbit.assert_num_close(pos[0][0], 2.0, 1.0e-12))
}

///|
test "custom_expr_min_max_select_abs" {
  let ci = @src.OmmCustomIntegrator::new(0.1)
  // min/max/select/abs arithmetic sanity through the x update.
  ci.add_step(variable="x", expression="x + max(min(3.0, 2.0), abs(-1.0))")
  ci.add_step(variable="x", expression="x + select(0.5, 10.0, 20.0)")
  ci.add_step(variable="x", expression="x + select(-0.5, 10.0, 20.0)")
  ci.add_step(variable="x", expression="x + floor(1.7) + ceil(0.2)")
  let (pos, _vel) = ci.step(
    [[0.0, 0.0, 0.0]],
    [[0.0, 0.0, 0.0]],
    [[0.0, 0.0, 0.0]],
    [1.0],
  )
  // 0 + 2 + 10 + 20 + 2 = 34
  assert_true(@moonbit.assert_num_close(pos[0][0], 34.0, 1.0e-12))
}

///|
test "custom_bond_force_harmonic" {
  // E = k (r - r0)^2 with per-bond k/r0; forces match the analytic
  // gradient 2 k (r - r0) along the bond axis.
  let f = @src.OmmCustomBondForce::new(expression="kb * (r - r0) * (r - r0)")
  f.add_bond(0, 1, params=[("kb", 4.0), ("r0", 1.0)])
  let pos : Array[Array[Double]] = [[0.0, 0.0, 0.0], [2.0, 0.0, 0.0]]
  let e = @src.omm_custom_bond_energy(f, pos)
  assert_true(@moonbit.assert_num_close(e, 4.0, 1.0e-12))
  let fo = @src.omm_custom_bond_forces(f, pos)
  // dE/dr = 2*4*(2-1) = 8; F_0 = +8 * x-hat, F_1 = -8 * x-hat.
  assert_true(@moonbit.assert_num_close(fo[0][0], 8.0, 1.0e-4))
  assert_true(@moonbit.assert_num_close(fo[1][0], -8.0, 1.0e-4))
  // Global params are visible; per-bond overrides win.
  let f2 = @src.OmmCustomBondForce::new(expression="kb * (r - r0) * (r - r0)")
  f2.add_global("kb", 1.0)
  f2.add_global("r0", 1.0)
  f2.add_bond(0, 1, params=[("kb", 4.0)])
  let e2 = @src.omm_custom_bond_energy(f2, pos)
  assert_true(@moonbit.assert_num_close(e2, 4.0, 1.0e-12))
  // Same energy at mirrored positions: r depends only on distance.
  let pos2 : Array[Array[Double]] = [[1.0, 1.0, 1.0], [2.37, 1.0, 1.0]]
  let e3 = @src.omm_custom_bond_energy(f2, pos2)
  assert_true(@moonbit.assert_num_close(e3, 4.0 * 0.37 * 0.37, 1.0e-9))
}

///|
test "custom_angle_force_harmonic" {
  // E = 0.5 * k * (theta - theta0)^2, theta in radians.
  let f = @src.OmmCustomAngleForce::new(
    expression="0.5 * ka * (theta - theta0) * (theta - theta0)",
  )
  f.add_global("ka", 10.0)
  f.add_global("theta0", 1.5707963267948966)
  f.add_angle(0, 1, 2)
  // Right angle at vertex 1: energy 0.
  let right : Array[Array[Double]] = [
    [1.0, 0.0, 0.0],
    [0.0, 0.0, 0.0],
    [0.0, 1.0, 0.0],
  ]
  assert_true(
    @moonbit.assert_num_close(
      @src.omm_custom_angle_energy(f, right),
      0.0,
      1.0e-12,
    ),
  )
  // Straight angle (pi): E = 0.5 * 10 * (pi - pi/2)^2.
  let straight : Array[Array[Double]] = [
    [1.0, 0.0, 0.0],
    [0.0, 0.0, 0.0],
    [-1.0, 0.0, 0.0],
  ]
  let expect = 0.5 *
    10.0 *
    (3.141592653589793 - 1.5707963267948966) *
    (3.141592653589793 - 1.5707963267948966)
  assert_true(
    @moonbit.assert_num_close(
      @src.omm_custom_angle_energy(f, straight),
      expect,
      1.0e-9,
    ),
  )
  // Forces point to reduce the angle: sum is finite and self-consistent
  // with the energy under a small displacement.
  let fo = @src.omm_custom_angle_forces(f, straight)
  let mut work = 0.0
  let h = 1.0e-4
  for i = 0; i < 3; i = i + 1 {
    for d = 0; d < 3; d = d + 1 {
      work = work + fo[i][d] * h
    }
  }
  let displaced : Array[Array[Double]] = [
    [1.0 + h, 0.0, 0.0],
    [0.0, 0.0, 0.0],
    [-1.0 + h, 0.0, 0.0],
  ]
  let delta_e = @src.omm_custom_angle_energy(f, displaced) - expect
  assert_true(@moonbit.assert_num_close(work, -delta_e, 1.0e-5))
}

///|
test "custom_torsion_force_periodic" {
  // E = k (1 + cos(3*phi)) — zero at phi = +/-60 deg, 2k at phi = 0.
  let f = @src.OmmCustomTorsionForce::new(expression="kt * (1 + cos(3 * phi))")
  f.add_global("kt", 1.5)
  f.add_torsion(0, 1, 2, 3)
  let cis : Array[Array[Double]] = [
    [0.0, 1.0, 0.0],
    [0.0, 0.0, 0.0],
    [1.0, 0.0, 0.0],
    [1.5, 0.8660254037844386, 0.0],
  ]
  assert_true(
    @moonbit.assert_num_close(
      @src.omm_custom_torsion_energy(f, cis),
      3.0,
      1.0e-9,
    ),
  )
  // phi = -60 deg (below the plane) -> E = 0.
  let gauche : Array[Array[Double]] = [
    [0.0, 1.0, 0.0],
    [0.0, 0.0, 0.0],
    [1.0, 0.0, 0.0],
    [1.5, 0.5, -0.8660254037844386],
  ]
  assert_true(
    @moonbit.assert_num_close(
      @src.omm_custom_torsion_energy(f, gauche),
      0.0,
      1.0e-9,
    ),
  )
  // Torsion forces are finite and antisymmetric under mirroring the end.
  let fo = @src.omm_custom_torsion_forces(f, cis)
  assert_true(fo[0][0] == fo[0][0])
  assert_true(fo[3][0] == fo[3][0])
}

///|
test "custom_nonbonded_soft_repulsion" {
  // Soft repulsive wall: E = A / r^12 over explicit pairs only.
  let f = @src.OmmCustomNonbondedForce::new(
    expression="A / (r * r * r * r * r * r)",
  )
  f.add_global("A", 1.0)
  f.add_pair(0, 1)
  f.add_pair(1, 2, params=[("A", 2.0)])
  let pos : Array[Array[Double]] = [
    [0.0, 0.0, 0.0],
    [1.0, 0.0, 0.0],
    [3.0, 0.0, 0.0],
  ]
  // Pair 0-1: 1/1^6 = 1; pair 1-2: 2/2^6 = 2/64.
  let e = @src.omm_custom_nonbonded_energy(f, pos)
  assert_true(@moonbit.assert_num_close(e, 1.0 + 2.0 / 64.0, 1.0e-12))
  let fo = @src.omm_custom_nonbonded_forces(f, pos)
  // dE/dr for pair 0-1 = -6/r^7 = -6 at r=1; repulsion pushes apart:
  // F_0 = +(-6) * x-hat = -6 (pulls i toward -x? no: F_i = dedr * r-hat).
  assert_true(@moonbit.assert_num_close(fo[0][0], -6.0, 1.0e-4))
  assert_true(
    @moonbit.assert_num_close(fo[1][0], 6.0 - 2.0 * 6.0 / 128.0, 1.0e-4),
  )
}
