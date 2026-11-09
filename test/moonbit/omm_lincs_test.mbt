///|
test "lincs_recovers_single_bond" {
  // One constrained bond pulled apart: positions and velocities restored.
  let constraints : Array[@src.OmmLincsConstraint] = [{ i: 0, j: 1, d: 1.0, }]
  let masses = [1.0, 1.0]
  let pos : Array[Array[Double]] = [[0.0, 0.0, 0.0], [1.4, 0.0, 0.0]]
  let vel : Array[Array[Double]] = [[0.0, 0.0, 0.0], [0.3, 0.0, 0.0]]
  let (p2, v2) = @src.omm_lincs_apply(pos, vel, constraints, masses)
  let dx = p2[1][0] - p2[0][0]
  let dy = p2[1][1] - p2[0][1]
  let dz = p2[1][2] - p2[0][2]
  let r = @math.pow(dx * dx + dy * dy + dz * dz, 0.5)
  assert_true(@moonbit.assert_num_close(r, 1.0, 1.0e-6))
  // Radial velocity component removed: (v_j - v_i) . rhat = 0.
  let rhat_x = dx / r
  let radial = (v2[1][0] - v2[0][0]) * rhat_x
  assert_true(radial.abs() < 1.0e-6)
  // COM velocity preserved by the projection.
  let com0 = (vel[0][0] + vel[1][0]) / 2.0
  let com1 = (v2[0][0] + v2[1][0]) / 2.0
  assert_true(@moonbit.assert_num_close(com0, com1, 1.0e-9))
}

///|
test "lincs_three_bond_chain" {
  // A 3-atom chain with two constraints, unequal masses.
  let constraints : Array[@src.OmmLincsConstraint] = [
    { i: 0, j: 1, d: 1.0, },
    { i: 1, j: 2, d: 1.0, },
  ]
  let masses = [12.0, 1.0, 12.0]
  let pos : Array[Array[Double]] = [
    [0.0, 0.0, 0.0],
    [1.3, 0.2, 0.0],
    [2.2, 0.4, 0.0],
  ]
  let vel : Array[Array[Double]] = [
    [0.1, 0.0, 0.0],
    [0.0, 0.2, 0.0],
    [0.3, 0.0, 0.1],
  ]
  let (p2, v2) = @src.omm_lincs_apply(pos, vel, constraints, masses)
  // Both bond lengths at their targets.
  let dx1 = p2[1][0] - p2[0][0]
  let dy1 = p2[1][1] - p2[0][1]
  let r1 = @math.pow(dx1 * dx1 + dy1 * dy1, 0.5)
  assert_true(@moonbit.assert_num_close(r1, 1.0, 1.0e-5))
  let dx2 = p2[2][0] - p2[1][0]
  let dy2 = p2[2][1] - p2[1][1]
  let r2 = @math.pow(dx2 * dx2 + dy2 * dy2, 0.5)
  assert_true(@moonbit.assert_num_close(r2, 1.0, 1.0e-5))
  // Velocities transverse to both bonds.
  let rad1 = (v2[1][0] - v2[0][0]) * (dx1 / r1) +
    (v2[1][1] - v2[0][1]) * (dy1 / r1)
  assert_true(rad1.abs() < 1.0e-5)
  let rad2 = (v2[2][0] - v2[1][0]) * (dx2 / r2) +
    (v2[2][1] - v2[1][1]) * (dy2 / r2)
  assert_true(rad2.abs() < 1.0e-5)
  // Mass-weighted COM displacement is preserved (projection property).
  let com0_x = (12.0 * pos[0][0] + 1.0 * pos[1][0] + 12.0 * pos[2][0]) / 25.0
  let com1_x = (12.0 * p2[0][0] + 1.0 * p2[1][0] + 12.0 * p2[2][0]) / 25.0
  assert_true(@moonbit.assert_num_close(com0_x, com1_x, 1.0e-6))
}

///|
test "lincs_water_geometry" {
  // Rigid water: two O-H constraints + no H-H constraint (LINCS handles
  // the coupled pair through the power series).
  let constraints : Array[@src.OmmLincsConstraint] = [
    { i: 0, j: 1, d: 0.9572, },
    { i: 0, j: 2, d: 0.9572, },
  ]
  let masses = [15.999, 1.008, 1.008]
  let pos : Array[Array[Double]] = [
    [0.0, 0.0, 0.0],
    [0.8, 0.5, 0.0],
    [-0.4, 0.9, 0.0],
  ]
  let vel : Array[Array[Double]] = [
    [0.0, 0.0, 0.0],
    [0.5, 0.0, 0.0],
    [0.0, 0.5, 0.0],
  ]
  let (p2, _v2) = @src.omm_lincs_apply(pos, vel, constraints, masses, n_iter=8)
  let d1 = @math.pow(
    (p2[1][0] - p2[0][0]) * (p2[1][0] - p2[0][0]) +
    (p2[1][1] - p2[0][1]) * (p2[1][1] - p2[0][1]),
    0.5,
  )
  let d2 = @math.pow(
    (p2[2][0] - p2[0][0]) * (p2[2][0] - p2[0][0]) +
    (p2[2][1] - p2[0][1]) * (p2[2][1] - p2[0][1]),
    0.5,
  )
  assert_true(@moonbit.assert_num_close(d1, 0.9572, 1.0e-5))
  assert_true(@moonbit.assert_num_close(d2, 0.9572, 1.0e-5))
  // All finite.
  for i = 0; i < 3; i = i + 1 {
    for d = 0; d < 3; d = d + 1 {
      assert_true(p2[i][d] == p2[i][d])
    }
  }
}
