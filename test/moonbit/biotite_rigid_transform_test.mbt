///|
/// Tests for biotite_rigid_transform.mbt — rigid transforms and fitting.
test "biotite_rigid_identity" {
  let t = @src.BiotiteRigidTransformation::identity()
  let p = t.apply([1.0, 2.0, 3.0])
  assert_true(@moonbit.assert_num_close(p[0], 1.0, 1.0e-12))
  assert_true(@moonbit.assert_num_close(p[1], 2.0, 1.0e-12))
  assert_true(@moonbit.assert_num_close(p[2], 3.0, 1.0e-12))
}

///|
test "biotite_rigid_translation_and_inverse" {
  let t = @src.BiotiteRigidTransformation::new(
    [1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0],
    [1.0, 2.0, 3.0],
  )
  let p = t.apply([1.0, 1.0, 1.0])
  assert_true(@moonbit.assert_num_close(p[0], 2.0, 1.0e-12))
  assert_true(@moonbit.assert_num_close(p[1], 3.0, 1.0e-12))
  assert_true(@moonbit.assert_num_close(p[2], 4.0, 1.0e-12))
  let roundtrip = t.inverse().apply(p)
  assert_true(@moonbit.assert_num_close(roundtrip[0], 1.0, 1.0e-12))
  assert_true(@moonbit.assert_num_close(roundtrip[1], 1.0, 1.0e-12))
  assert_true(@moonbit.assert_num_close(roundtrip[2], 1.0, 1.0e-12))
}

///|
test "biotite_rigid_compose" {
  let translate = @src.BiotiteRigidTransformation::new(
    [1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0],
    [1.0, 0.0, 0.0],
  )
  let composed = translate.compose(translate)
  let p = composed.apply([0.0, 0.0, 0.0])
  assert_true(@moonbit.assert_num_close(p[0], 2.0, 1.0e-12))
}

///|
test "biotite_rigid_from_alignment_translation" {
  // Three non-collinear points avoid the degenerate case of the quaternion fit.
  let fixed : Array[Double] = [1.0, 2.0, 3.0, 2.0, 2.0, 3.0, 1.0, 3.0, 3.0]
  let mobile : Array[Double] = [0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0, 0.0]
  let t = @src.BiotiteRigidTransformation::from_alignment(fixed, mobile)
  assert_true(t.rmsd(fixed, mobile) < 0.001)
}

///|
test "biotite_rigid_from_alignment_rotation" {
  let fixed : Array[Double] = [0.0, 1.0, 0.0, -1.0, 0.0, 0.0, 0.0, 0.0, 1.0]
  let mobile : Array[Double] = [1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0]
  let t = @src.BiotiteRigidTransformation::from_alignment(fixed, mobile)
  assert_true(t.rmsd(fixed, mobile) < 0.001)
}
