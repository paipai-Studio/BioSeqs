///|
/// Tests for biotite_space_group.mbt — crystallographic symmetry.
test "biotite_space_group_p1_identity" {
  let ops = @src.biotite_space_group_transforms("P1")
  assert_eq(ops.length(), 1)
  let p = ops[0].apply([0.1, 0.2, 0.3])
  assert_true(@moonbit.assert_num_close(p[0], 0.1, 1.0e-12))
  assert_true(@moonbit.assert_num_close(p[1], 0.2, 1.0e-12))
  assert_true(@moonbit.assert_num_close(p[2], 0.3, 1.0e-12))
}

///|
test "biotite_space_group_inversion" {
  let ops = @src.biotite_space_group_transforms("P-1")
  assert_eq(ops.length(), 2)
  let p = ops[1].apply([0.2, 0.3, 0.4])
  assert_true(@moonbit.assert_num_close(p[0], -0.2, 1.0e-12))
  assert_true(@moonbit.assert_num_close(p[1], -0.3, 1.0e-12))
  assert_true(@moonbit.assert_num_close(p[2], -0.4, 1.0e-12))
}

///|
test "biotite_space_group_screw_axis" {
  let ops = @src.biotite_space_group_transforms("P21")
  assert_eq(ops.length(), 2)
  let p = ops[1].apply([0.1, 0.1, 0.1])
  assert_true(@moonbit.assert_num_close(p[0], -0.1, 1.0e-12))
  assert_true(@moonbit.assert_num_close(p[1], 0.6, 1.0e-12))
  assert_true(@moonbit.assert_num_close(p[2], -0.1, 1.0e-12))
}

///|
test "biotite_space_group_equivalent_positions" {
  let positions = @src.biotite_space_group_equivalent_positions("P-1", [
    0.2, 0.3, 0.4,
  ])
  assert_eq(positions.length(), 2)
  assert_true(@moonbit.assert_num_close(positions[1][0], 0.8, 1.0e-12))
  assert_true(@moonbit.assert_num_close(positions[1][1], 0.7, 1.0e-12))
  assert_true(@moonbit.assert_num_close(positions[1][2], 0.6, 1.0e-12))
}

///|
test "biotite_space_group_operation_inverse" {
  let ops = @src.biotite_space_group_transforms("P212121")
  let op = ops[1]
  let coord : Array[Double] = [0.25, 0.1, 0.3]
  let roundtrip = op.inverse().apply(op.apply(coord))
  assert_true(@moonbit.assert_num_close(roundtrip[0], 0.25, 1.0e-12))
  assert_true(@moonbit.assert_num_close(roundtrip[1], 0.1, 1.0e-12))
  assert_true(@moonbit.assert_num_close(roundtrip[2], 0.3, 1.0e-12))
}

///|
test "biotite_space_group_unknown_falls_back_to_identity" {
  let ops = @src.biotite_space_group_transforms("UNKNOWN")
  assert_eq(ops.length(), 1)
}
