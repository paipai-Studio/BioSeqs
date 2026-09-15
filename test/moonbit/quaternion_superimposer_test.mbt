///|
/// Tests for quaternion_superimposer module (Bio.PDB.Superimposer - Quaternion method).

test "qs_point_norm" {
  let p = @src.QsPoint::new(3.0, 4.0, 0.0)
  let n = @src.qs_point_norm(p)
  assert_true((n - 5.0).abs() < 1.0e-9)
}

test "qs_point_dot" {
  let a = @src.QsPoint::new(1.0, 0.0, 0.0)
  let b = @src.QsPoint::new(0.0, 1.0, 0.0)
  let d = @src.qs_point_dot(a, b)
  assert_eq(d, 0.0)
}

test "qs_point_cross" {
  let a = @src.QsPoint::new(1.0, 0.0, 0.0)
  let b = @src.QsPoint::new(0.0, 1.0, 0.0)
  let c = @src.qs_point_cross(a, b)
  assert_true(c.x.abs() < 1.0e-9)
  assert_true(c.y.abs() < 1.0e-9)
  assert_true((c.z - 1.0).abs() < 1.0e-9)
}

test "quaternion_identity" {
  let q = @src.quaternion_identity()
  assert_true((q.w - 1.0).abs() < 1.0e-9)
  assert_true(q.x.abs() < 1.0e-9)
  assert_true(q.y.abs() < 1.0e-9)
  assert_true(q.z.abs() < 1.0e-9)
}

test "quaternion_multiply" {
  let q1 = @src.quaternion_identity()
  let q2 = @src.quaternion_identity()
  let q = @src.quaternion_multiply(q1, q2)
  assert_true((q.w - 1.0).abs() < 1.0e-9)
}

test "quaternion_conjugate" {
  let q = @src.Quaternion::new(1.0, 2.0, 3.0, 4.0)
  let qc = @src.quaternion_conjugate(q)
  assert_true((qc.w - 1.0).abs() < 1.0e-9)
  assert_true((qc.x + 2.0).abs() < 1.0e-9)
  assert_true((qc.y + 3.0).abs() < 1.0e-9)
  assert_true((qc.z + 4.0).abs() < 1.0e-9)
}

test "quaternion_to_matrix" {
  let q = @src.quaternion_identity()
  let mat = @src.quaternion_to_matrix(q)
  // Identity quaternion should give identity matrix
  assert_true((mat[0][0] - 1.0).abs() < 1.0e-9)
  assert_true((mat[1][1] - 1.0).abs() < 1.0e-9)
  assert_true((mat[2][2] - 1.0).abs() < 1.0e-9)
}

test "qs_centroid" {
  let points = [
    @src.QsPoint::new(0.0, 0.0, 0.0),
    @src.QsPoint::new(2.0, 0.0, 0.0),
    @src.QsPoint::new(0.0, 2.0, 0.0),
  ]
  let c = @src.qs_centroid(points)
  assert_true((c.x - 2.0/3.0).abs() < 1.0e-9)
  assert_true((c.y - 2.0/3.0).abs() < 1.0e-9)
}

test "qs_superimpose identical" {
  let a = [
    @src.QsPoint::new(1.0, 0.0, 0.0),
    @src.QsPoint::new(0.0, 1.0, 0.0),
    @src.QsPoint::new(0.0, 0.0, 1.0),
  ]
  let b = [
    @src.QsPoint::new(1.0, 0.0, 0.0),
    @src.QsPoint::new(0.0, 1.0, 0.0),
    @src.QsPoint::new(0.0, 0.0, 1.0),
  ]
  let result = @src.qs_superimpose(a, b)
  // Identical structures should have RMSD = 0
  assert_true(result.rmsd < 1.0e-9)
}

test "qs_superimpose rotated" {
  let a = [
    @src.QsPoint::new(1.0, 0.0, 0.0),
    @src.QsPoint::new(0.0, 1.0, 0.0),
    @src.QsPoint::new(0.0, 0.0, 1.0),
  ]
  let b = [
    @src.QsPoint::new(0.0, 1.0, 0.0),
    @src.QsPoint::new(-1.0, 0.0, 0.0),
    @src.QsPoint::new(0.0, 0.0, 1.0),
  ]
  let result = @src.qs_superimpose(a, b)
  // b is a rotated by +90 deg about z, so the fit should be exact
  assert_true(result.rmsd < 1.0e-9)
}


test "qs_superimpose translated" {
  let a = [
    @src.QsPoint::new(1.0, 0.0, 0.0),
    @src.QsPoint::new(0.0, 1.0, 0.0),
    @src.QsPoint::new(0.0, 0.0, 1.0),
  ]
  let b = [
    @src.QsPoint::new(2.0, 1.0, 1.0),
    @src.QsPoint::new(1.0, 2.0, 1.0),
    @src.QsPoint::new(1.0, 1.0, 2.0),
  ]
  let result = @src.qs_superimpose(a, b)
  // Translation is removed by centering, so the fit should be exact
  assert_true(result.rmsd < 1.0e-9)
}

test "qs_superimpose empty" {
  let a : Array[@src.QsPoint] = []
  let b : Array[@src.QsPoint] = []
  let result = @src.qs_superimpose(a, b)
  assert_eq(result.rmsd, 0.0)
}

test "qs_apply_transform" {
  let points = [
    @src.QsPoint::new(1.0, 0.0, 0.0),
    @src.QsPoint::new(0.0, 1.0, 0.0),
  ]
  let rot = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]
  let trans = [1.0, 1.0, 1.0]
  let transformed = @src.qs_apply_transform(points, rot, trans)
  assert_true((transformed[0].x - 2.0).abs() < 1.0e-9)
  assert_true((transformed[1].y - 2.0).abs() < 1.0e-9)
}

test "qs_rmsd" {
  let a = [
    @src.QsPoint::new(0.0, 0.0, 0.0),
    @src.QsPoint::new(1.0, 0.0, 0.0),
  ]
  let b = [
    @src.QsPoint::new(0.0, 0.0, 0.0),
    @src.QsPoint::new(1.0, 0.0, 0.0),
  ]
  let rmsd = @src.qs_rmsd(a, b)
  assert_eq(rmsd, 0.0)
}

test "qs_example_structures" {
  let (a, b) = @src.qs_example_structures()
  assert_eq(a.length(), 3)
  assert_eq(b.length(), 3)
  let result = @src.qs_superimpose(a, b)
  assert_true(result.rmsd < 1.0e-9)
}

test "qs_superimpose_apply_transform_roundtrip" {
  let (a, b) = @src.qs_example_structures()
  let result = @src.qs_superimpose(a, b)
  let fitted = @src.qs_apply_transform(b, result.rotation, result.translation)
  assert_true(@src.qs_rmsd(fitted, a) < 1.0e-9)
}

test "qs_superimpose_rotation_is_proper" {
  let (a, b) = @src.qs_example_structures()
  let r = @src.qs_superimpose(a, b).rotation
  let det = r[0][0] * (r[1][1] * r[2][2] - r[1][2] * r[2][1]) -
    r[0][1] * (r[1][0] * r[2][2] - r[1][2] * r[2][0]) +
    r[0][2] * (r[1][0] * r[2][1] - r[1][1] * r[2][0])
  assert_true((det - 1.0).abs() < 1.0e-9)
}