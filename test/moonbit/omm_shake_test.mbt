///|
test "omm_shake_restores_bond_length" {
  // single O-H1 constraint taken from the TIP3P template
  let cs = @src.omm_tip3p_constraints(0, 1, 2)
  let constraints = [cs[0]]
  let old = [[0.0, 0.0, 0.0], [0.9572, 0.0, 0.0]]
  let unconstrained = [[0.0, 0.0, 0.0], [1.4358, 0.0, 0.0]]
  let masses = [1.0, 1.0]
  let result = @src.omm_shake(unconstrained, old, constraints, masses)
  let dx = result[1][0] - result[0][0]
  assert_true(@moonbit.assert_num_close(dx, 0.9572, 1.0e-9))
  // equal masses: midpoint preserved
  let mid = (result[0][0] + result[1][0]) * 0.5
  assert_true(@moonbit.assert_num_close(mid, 0.7179, 1.0e-9))
}

///|
test "omm_shake_no_constraints_passthrough" {
  let positions = [[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]]
  let old = [[0.0, 0.0, 0.0], [0.0, 0.0, 0.0]]
  let result = @src.omm_shake(positions, old, [], [1.0, 1.0])
  assert_eq(result.length(), 2)
  assert_true(@moonbit.assert_num_close(result[0][0], 1.0, 1.0e-12))
  assert_true(@moonbit.assert_num_close(result[1][2], 6.0, 1.0e-12))
}

///|
test "omm_shake_restores_tip3p_water_geometry" {
  let d_oh = 0.9572
  let d_hh = 1.5139
  let half = @math.asin(d_hh / (2.0 * d_oh))
  let hx = d_oh * @math.sin(half)
  let hz = d_oh * @math.cos(half)
  let old = [[0.0, 0.0, 0.0], [hx, 0.0, hz], [-hx, 0.0, hz]]
  // stretch H1 along the O-H1 direction by 20 percent
  let unconstrained = [
    [0.0, 0.0, 0.0],
    [hx * 1.2, 0.0, hz * 1.2],
    [-hx, 0.0, hz],
  ]
  let constraints = @src.omm_tip3p_constraints(0, 1, 2)
  let masses = [15.9994, 1.008, 1.008]
  let result = @src.omm_shake(unconstrained, old, constraints, masses)
  let r_oh1 = @math.pow(
    (result[1][0] - result[0][0]) * (result[1][0] - result[0][0]) +
    (result[1][1] - result[0][1]) * (result[1][1] - result[0][1]) +
    (result[1][2] - result[0][2]) * (result[1][2] - result[0][2]),
    0.5,
  )
  assert_true(@moonbit.assert_num_close(r_oh1, 0.9572, 1.0e-5))
  let r_oh2 = @math.pow(
    (result[2][0] - result[0][0]) * (result[2][0] - result[0][0]) +
    (result[2][1] - result[0][1]) * (result[2][1] - result[0][1]) +
    (result[2][2] - result[0][2]) * (result[2][2] - result[0][2]),
    0.5,
  )
  assert_true(@moonbit.assert_num_close(r_oh2, 0.9572, 1.0e-5))
  let r_hh = @math.pow(
    (result[2][0] - result[1][0]) * (result[2][0] - result[1][0]) +
    (result[2][1] - result[1][1]) * (result[2][1] - result[1][1]) +
    (result[2][2] - result[1][2]) * (result[2][2] - result[1][2]),
    0.5,
  )
  assert_true(@moonbit.assert_num_close(r_hh, 1.5139, 1.0e-5))
}

///|
test "omm_tip3p_and_tip4p_constraints" {
  let c = @src.omm_tip3p_constraints(0, 1, 2)
  assert_eq(c.length(), 3)
  assert_eq(c[0].i, 0)
  assert_eq(c[0].j, 1)
  assert_true(@moonbit.assert_num_close(c[0].distance, 0.9572, 1.0e-12))
  assert_eq(c[1].i, 0)
  assert_eq(c[1].j, 2)
  assert_eq(c[2].i, 1)
  assert_eq(c[2].j, 2)
  assert_true(@moonbit.assert_num_close(c[2].distance, 1.5139006545, 1.0e-12))
  let c4 = @src.omm_tip4p_constraints(3, 4, 5)
  assert_eq(c4.length(), 3)
  assert_eq(c4[0].i, 3)
  assert_eq(c4[0].j, 4)
  assert_true(@moonbit.assert_num_close(c4[1].distance, 0.9572, 1.0e-12))
}

///|
test "omm_rattle_removes_radial_relative_velocity" {
  let positions = [[0.0, 0.0, 0.0], [0.9572, 0.0, 0.0]]
  let cs = @src.omm_tip3p_constraints(0, 1, 2)
  let constraints = [cs[0]]
  let masses = [1.0, 1.0]
  // relative velocity along the bond is redistributed equally
  let velocities = [[1.0, 0.0, 0.0], [2.0, 0.0, 0.0]]
  let v_new = @src.omm_rattle(velocities, positions, constraints, masses)
  assert_true(@moonbit.assert_num_close(v_new[0][0], 1.5, 1.0e-9))
  assert_true(@moonbit.assert_num_close(v_new[1][0], 1.5, 1.0e-9))
  // purely tangential relative velocity is untouched
  let tangential = [[0.0, 1.0, 0.0], [0.0, 1.0, 0.0]]
  let v2 = @src.omm_rattle(tangential, positions, constraints, masses)
  assert_true(@moonbit.assert_num_close(v2[0][1], 1.0, 1.0e-9))
  assert_true(@moonbit.assert_num_close(v2[1][1], 1.0, 1.0e-9))
}
