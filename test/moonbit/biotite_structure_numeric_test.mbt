///|
/// Tests for biotite_structure_numeric.mbt.
test "biotite_vector_dot_and_angle" {
  let x = @src.Vector3::new(1.0, 0.0, 0.0)
  let y = @src.Vector3::new(0.0, 1.0, 0.0)
  assert_true(
    @moonbit.assert_num_close(@src.biotite_vector_dot(x, y), 0.0, 1.0e-12),
  )
  assert_true(
    @moonbit.assert_num_close(@src.biotite_vector_angle(x, y), 90.0, 1.0e-6),
  )
}

///|
test "biotite_vector_distance_and_norm" {
  let a = @src.Vector3::new(0.0, 0.0, 0.0)
  let b = @src.Vector3::new(3.0, 4.0, 0.0)
  assert_true(
    @moonbit.assert_num_close(@src.biotite_vector_distance(a, b), 5.0, 1.0e-12),
  )
  assert_true(
    @moonbit.assert_num_close(@src.biotite_vector_norm(b), 5.0, 1.0e-12),
  )
}

///|
test "biotite_vector_dihedral" {
  let p1 = @src.Vector3::new(1.0, 0.0, 0.0)
  let p2 = @src.Vector3::new(0.0, 0.0, 0.0)
  let p3 = @src.Vector3::new(0.0, 1.0, 0.0)
  let p4 = @src.Vector3::new(1.0, 1.0, 0.0)
  assert_true(
    @moonbit.assert_num_close(
      @src.biotite_vector_dihedral(p1, p2, p3, p4),
      0.0,
      1.0e-6,
    ),
  )
}
