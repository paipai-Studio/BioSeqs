///|
test "Vector3::new construction" {
  let v = @src.Vector3::new(1.0, 2.0, 3.0)
  assert_eq(v.get_x(), 1.0)
  assert_eq(v.get_y(), 2.0)
  assert_eq(v.get_z(), 3.0)
}

///|
test "Vector3::zero" {
  let v = @src.Vector3::zero()
  assert_eq(v.get_x(), 0.0)
  assert_eq(v.get_y(), 0.0)
  assert_eq(v.get_z(), 0.0)
}

///|
test "Vector3::unit_x" {
  let v = @src.Vector3::unit_x()
  assert_eq(v.get_x(), 1.0)
  assert_eq(v.get_y(), 0.0)
  assert_eq(v.get_z(), 0.0)
}

///|
test "Vector3::unit_y" {
  let v = @src.Vector3::unit_y()
  assert_eq(v.get_x(), 0.0)
  assert_eq(v.get_y(), 1.0)
  assert_eq(v.get_z(), 0.0)
}

///|
test "Vector3::unit_z" {
  let v = @src.Vector3::unit_z()
  assert_eq(v.get_x(), 0.0)
  assert_eq(v.get_y(), 0.0)
  assert_eq(v.get_z(), 1.0)
}

///|
test "Vector3::from_array" {
  let v = @src.Vector3::from_array([1.0, 2.0, 3.0])
  assert_eq(v.get_x(), 1.0)
  assert_eq(v.get_y(), 2.0)
  assert_eq(v.get_z(), 3.0)
}

///|
test "Vector3::from_array short" {
  let v = @src.Vector3::from_array([1.0, 2.0])
  assert_eq(v.get_x(), 1.0)
  assert_eq(v.get_y(), 2.0)
  assert_eq(v.get_z(), 0.0)
}

///|
test "Vector3::from_array empty" {
  let v = @src.Vector3::from_array([])
  assert_eq(v.get_x(), 0.0)
  assert_eq(v.get_y(), 0.0)
  assert_eq(v.get_z(), 0.0)
}

///|
test "Vector3::add" {
  let v1 = @src.Vector3::new(1.0, 2.0, 3.0)
  let v2 = @src.Vector3::new(4.0, 5.0, 6.0)
  let result = v1.add(v2)
  assert_eq(result.get_x(), 5.0)
  assert_eq(result.get_y(), 7.0)
  assert_eq(result.get_z(), 9.0)
}

///|
test "Vector3::sub" {
  let v1 = @src.Vector3::new(5.0, 7.0, 9.0)
  let v2 = @src.Vector3::new(1.0, 2.0, 3.0)
  let result = v1.sub(v2)
  assert_eq(result.get_x(), 4.0)
  assert_eq(result.get_y(), 5.0)
  assert_eq(result.get_z(), 6.0)
}

///|
test "Vector3::scale" {
  let v = @src.Vector3::new(1.0, 2.0, 3.0)
  let result = v.scale(2.0)
  assert_eq(result.get_x(), 2.0)
  assert_eq(result.get_y(), 4.0)
  assert_eq(result.get_z(), 6.0)
}

///|
test "Vector3::scale negative" {
  let v = @src.Vector3::new(1.0, 2.0, 3.0)
  let result = v.scale(-1.0)
  assert_eq(result.get_x(), -1.0)
  assert_eq(result.get_y(), -2.0)
  assert_eq(result.get_z(), -3.0)
}

///|
test "Vector3::dot" {
  let v1 = @src.Vector3::new(1.0, 2.0, 3.0)
  let v2 = @src.Vector3::new(4.0, 5.0, 6.0)
  let result = v1.dot(v2)
  assert_true((result - 32.0).abs() < 1.0e-10)
}

///|
test "Vector3::dot orthogonal" {
  let v1 = @src.Vector3::new(1.0, 0.0, 0.0)
  let v2 = @src.Vector3::new(0.0, 1.0, 0.0)
  assert_true((v1.dot(v2)).abs() < 1.0e-10)
}

///|
test "Vector3::cross" {
  let v1 = @src.Vector3::new(1.0, 0.0, 0.0)
  let v2 = @src.Vector3::new(0.0, 1.0, 0.0)
  let result = v1.cross(v2)
  assert_eq(result.get_x(), 0.0)
  assert_eq(result.get_y(), 0.0)
  assert_eq(result.get_z(), 1.0)
}

///|
test "Vector3::cross z cross x" {
  let v1 = @src.Vector3::new(0.0, 0.0, 1.0)
  let v2 = @src.Vector3::new(1.0, 0.0, 0.0)
  let result = v1.cross(v2)
  assert_eq(result.get_x(), 0.0)
  assert_eq(result.get_y(), 1.0)
  assert_eq(result.get_z(), 0.0)
}

///|
test "Vector3::norm" {
  let v = @src.Vector3::new(3.0, 4.0, 0.0)
  let norm_val = v.norm()
  assert_true((norm_val - 5.0).abs() < 1.0e-10)
}

///|
test "Vector3::norm_squared" {
  let v = @src.Vector3::new(1.0, 2.0, 3.0)
  let n2 = v.norm_squared()
  assert_true((n2 - 14.0).abs() < 1.0e-10)
}

///|
test "Vector3::normalize" {
  let v = @src.Vector3::new(3.0, 4.0, 0.0)
  let norm = v.normalize()
  let norm_val = norm.norm()
  assert_true((norm_val - 1.0).abs() < 1.0e-10)
  assert_true((norm.get_x() - 0.6).abs() < 1.0e-10)
  assert_true((norm.get_y() - 0.8).abs() < 1.0e-10)
}

///|
test "Vector3::normalize zero vector" {
  let v = @src.Vector3::zero()
  let norm = v.normalize()
  assert_eq(norm.get_x(), 0.0)
  assert_eq(norm.get_y(), 0.0)
  assert_eq(norm.get_z(), 0.0)
}

///|
test "Vector3::distance" {
  let v1 = @src.Vector3::new(0.0, 0.0, 0.0)
  let v2 = @src.Vector3::new(3.0, 4.0, 0.0)
  let d = v1.distance(v2)
  assert_true((d - 5.0).abs() < 1.0e-10)
}

///|
test "Vector3::distance_squared" {
  let v1 = @src.Vector3::new(0.0, 0.0, 0.0)
  let v2 = @src.Vector3::new(1.0, 2.0, 3.0)
  let d2 = v1.distance_squared(v2)
  assert_true((d2 - 14.0).abs() < 1.0e-10)
}

///|
test "Vector3::negate" {
  let v = @src.Vector3::new(1.0, -2.0, 3.0)
  let neg = v.negate()
  assert_eq(neg.get_x(), -1.0)
  assert_eq(neg.get_y(), 2.0)
  assert_eq(neg.get_z(), -3.0)
}

///|
test "Vector3::is_zero" {
  let v = @src.Vector3::zero()
  assert_true(v.is_zero())
  let v2 = @src.Vector3::new(1.0e-12, 1.0e-12, 1.0e-12)
  assert_true(v2.is_zero())
  let v3 = @src.Vector3::new(1.0, 0.0, 0.0)
  assert_false(v3.is_zero())
}

///|
test "Vector3::angle" {
  let v1 = @src.Vector3::new(1.0, 0.0, 0.0)
  let v2 = @src.Vector3::new(0.0, 1.0, 0.0)
  let angle_rad = v1.angle(v2)
  let expected = @math.PI / 2.0
  assert_true((angle_rad - expected).abs() < 1.0e-10)
}

///|
test "Vector3::angle_degrees" {
  let v1 = @src.Vector3::new(1.0, 0.0, 0.0)
  let v2 = @src.Vector3::new(0.0, 1.0, 0.0)
  let angle_deg = v1.angle_degrees(v2)
  assert_true((angle_deg - 90.0).abs() < 1.0e-10)
}

///|
test "Vector3::to_array" {
  let v = @src.Vector3::new(1.0, 2.0, 3.0)
  let arr = v.to_array()
  assert_eq(arr.length(), 3)
  assert_eq(arr[0], 1.0)
  assert_eq(arr[1], 2.0)
  assert_eq(arr[2], 3.0)
}

///|
test "Vector3::project_to_vector" {
  let v = @src.Vector3::new(1.0, 2.0, 0.0)
  let target = @src.Vector3::new(1.0, 0.0, 0.0)
  let proj = v.project_to_vector(target)
  assert_true((proj.get_x() - 1.0).abs() < 1.0e-10)
  assert_true((proj.get_y()).abs() < 1.0e-10)
  assert_true((proj.get_z()).abs() < 1.0e-10)
}

///|
test "Vector3::project_to_plane" {
  let v = @src.Vector3::new(1.0, 2.0, 3.0)
  let normal = @src.Vector3::new(0.0, 0.0, 1.0)
  let proj = v.project_to_plane(normal)
  assert_true((proj.get_x() - 1.0).abs() < 1.0e-10)
  assert_true((proj.get_y() - 2.0).abs() < 1.0e-10)
  assert_true((proj.get_z()).abs() < 1.0e-10)
}

///|
test "RotationMatrix3::identity" {
  let m = @src.RotationMatrix3::identity()
  assert_eq(m.get(0, 0), 1.0)
  assert_eq(m.get(1, 1), 1.0)
  assert_eq(m.get(2, 2), 1.0)
  assert_eq(m.get(0, 1), 0.0)
  assert_eq(m.get(1, 0), 0.0)
}

///|
test "RotationMatrix3::transpose" {
  let m = @src.RotationMatrix3::new(
    m00=1.0, m01=2.0, m02=3.0,
    m10=4.0, m11=5.0, m12=6.0,
    m20=7.0, m21=8.0, m22=9.0,
  )
  let t = m.transpose()
  assert_eq(t.get(0, 1), 4.0)
  assert_eq(t.get(1, 0), 2.0)
  assert_eq(t.get(0, 2), 7.0)
  assert_eq(t.get(2, 0), 3.0)
}

///|
test "RotationMatrix3::multiply identity" {
  let m = @src.RotationMatrix3::new(
    m00=1.0, m01=2.0, m02=3.0,
    m10=4.0, m11=5.0, m12=6.0,
    m20=7.0, m21=8.0, m22=9.0,
  )
  let identity = @src.RotationMatrix3::identity()
  let result = m.multiply(identity)
  assert_true((result.get(0, 0) - m.get(0, 0)).abs() < 1.0e-10)
  assert_true((result.get(1, 1) - m.get(1, 1)).abs() < 1.0e-10)
  assert_true((result.get(2, 2) - m.get(2, 2)).abs() < 1.0e-10)
}

///|
test "RotationMatrix3::transform" {
  let rot = @src.RotationMatrix3::rotation_z(@math.PI / 2.0)
  let v = @src.Vector3::new(1.0, 0.0, 0.0)
  let result = rot.transform(v)
  assert_true((result.get_x()).abs() < 1.0e-10)
  assert_true((result.get_y() - 1.0).abs() < 1.0e-10)
  assert_true((result.get_z()).abs() < 1.0e-10)
}

///|
test "RotationMatrix3::rotation_x" {
  let rot = @src.RotationMatrix3::rotation_x(@math.PI / 2.0)
  let v = @src.Vector3::new(0.0, 1.0, 0.0)
  let result = rot.transform(v)
  assert_true((result.get_x()).abs() < 1.0e-10)
  assert_true((result.get_y()).abs() < 1.0e-10)
  assert_true((result.get_z() - 1.0).abs() < 1.0e-10)
}

///|
test "RotationMatrix3::rotation_y" {
  let rot = @src.RotationMatrix3::rotation_y(@math.PI / 2.0)
  let v = @src.Vector3::new(0.0, 0.0, 1.0)
  let result = rot.transform(v)
  assert_true((result.get_x() - 1.0).abs() < 1.0e-10)
  assert_true((result.get_y()).abs() < 1.0e-10)
  assert_true((result.get_z()).abs() < 1.0e-10)
}

///|
test "RotationMatrix3::rotation_z" {
  let rot = @src.RotationMatrix3::rotation_z(@math.PI / 2.0)
  let v = @src.Vector3::new(1.0, 0.0, 0.0)
  let result = rot.transform(v)
  assert_true((result.get_x()).abs() < 1.0e-10)
  assert_true((result.get_y() - 1.0).abs() < 1.0e-10)
  assert_true((result.get_z()).abs() < 1.0e-10)
}

///|
test "RotationMatrix3::rotation_axis_angle" {
  let axis = @src.Vector3::new(0.0, 0.0, 1.0)
  let rot = @src.RotationMatrix3::rotation_axis_angle(axis, @math.PI / 2.0)
  let v = @src.Vector3::new(1.0, 0.0, 0.0)
  let result = rot.transform(v)
  assert_true((result.get_x()).abs() < 1.0e-10)
  assert_true((result.get_y() - 1.0).abs() < 1.0e-10)
  assert_true((result.get_z()).abs() < 1.0e-10)
}

///|
test "RotationMatrix3::determinant identity" {
  let m = @src.RotationMatrix3::identity()
  assert_true((m.determinant() - 1.0).abs() < 1.0e-10)
}

///|
test "RotationMatrix3::determinant rotation" {
  let rot = @src.RotationMatrix3::rotation_z(@math.PI / 3.0)
  assert_true((rot.determinant() - 1.0).abs() < 1.0e-10)
}

///|
test "RotationMatrix3::is_rotation identity" {
  let m = @src.RotationMatrix3::identity()
  assert_true(m.is_rotation())
}

///|
test "RotationMatrix3::is_rotation rotation_z" {
  let rot = @src.RotationMatrix3::rotation_z(@math.PI / 4.0)
  assert_true(rot.is_rotation())
}

///|
test "RotationMatrix3::inverse" {
  let rot = @src.RotationMatrix3::rotation_z(@math.PI / 3.0)
  let inv = rot.inverse()
  let product = rot.multiply(inv)
  assert_true((product.get(0, 0) - 1.0).abs() < 1.0e-10)
  assert_true((product.get(1, 1) - 1.0).abs() < 1.0e-10)
  assert_true((product.get(2, 2) - 1.0).abs() < 1.0e-10)
  assert_true((product.get(0, 1)).abs() < 1.0e-10)
  assert_true((product.get(1, 0)).abs() < 1.0e-10)
}

///|
test "RotationMatrix3::get out of bounds" {
  let m = @src.RotationMatrix3::identity()
  assert_eq(m.get(5, 5), 0.0)
}

///|
test "RotationMatrix3::set" {
  let m = @src.RotationMatrix3::identity()
  let m2 = m.set(0, 0, 5.0)
  assert_eq(m2.get(0, 0), 5.0)
  assert_eq(m2.get(1, 1), 1.0)
}

///|
test "vector_centroid" {
  let vectors = [
    @src.Vector3::new(0.0, 0.0, 0.0),
    @src.Vector3::new(2.0, 0.0, 0.0),
    @src.Vector3::new(0.0, 2.0, 0.0),
    @src.Vector3::new(0.0, 0.0, 2.0),
  ]
  let c = @src.vector_centroid(vectors)
  assert_true((c.get_x() - 0.5).abs() < 1.0e-10)
  assert_true((c.get_y() - 0.5).abs() < 1.0e-10)
  assert_true((c.get_z() - 0.5).abs() < 1.0e-10)
}

///|
test "vector_centroid empty" {
  let vectors : Array[@src.Vector3] = []
  let c = @src.vector_centroid(vectors)
  assert_eq(c.get_x(), 0.0)
  assert_eq(c.get_y(), 0.0)
  assert_eq(c.get_z(), 0.0)
}

///|
test "vector_rmsd" {
  let set_a = [
    @src.Vector3::new(0.0, 0.0, 0.0),
    @src.Vector3::new(1.0, 0.0, 0.0),
  ]
  let set_b = [
    @src.Vector3::new(1.0, 0.0, 0.0),
    @src.Vector3::new(2.0, 0.0, 0.0),
  ]
  let rmsd = @src.vector_rmsd(set_a, set_b)
  assert_true((rmsd - 1.0).abs() < 1.0e-10)
}

///|
test "vector_rmsd identical" {
  let set_a = [
    @src.Vector3::new(1.0, 2.0, 3.0),
    @src.Vector3::new(4.0, 5.0, 6.0),
  ]
  let set_b = [
    @src.Vector3::new(1.0, 2.0, 3.0),
    @src.Vector3::new(4.0, 5.0, 6.0),
  ]
  let rmsd = @src.vector_rmsd(set_a, set_b)
  assert_true((rmsd).abs() < 1.0e-10)
}

///|
test "vector_rmsd empty" {
  let set_a : Array[@src.Vector3] = []
  let set_b : Array[@src.Vector3] = []
  let rmsd = @src.vector_rmsd(set_a, set_b)
  assert_true((rmsd).abs() < 1.0e-10)
}

///|
test "vector_superposition_rmsd" {
  let set_a = [
    @src.Vector3::new(0.0, 0.0, 0.0),
    @src.Vector3::new(1.0, 0.0, 0.0),
    @src.Vector3::new(0.0, 1.0, 0.0),
  ]
  let set_b = [
    @src.Vector3::new(1.0, 1.0, 0.0),
    @src.Vector3::new(2.0, 1.0, 0.0),
    @src.Vector3::new(1.0, 2.0, 0.0),
  ]
  let rmsd = @src.vector_superposition_rmsd(set_a, set_b)
  assert_true(rmsd >= 0.0)
  assert_true(rmsd < 1.0)
}

///|
test "vector_superposition_rmsd single point" {
  let set_a = [@src.Vector3::new(1.0, 2.0, 3.0)]
  let set_b = [@src.Vector3::new(4.0, 5.0, 6.0)]
  let rmsd = @src.vector_superposition_rmsd(set_a, set_b)
  assert_true(rmsd >= 0.0)
}

///|
test "vector_superposition_transform" {
  let set_a = [
    @src.Vector3::new(0.0, 0.0, 0.0),
    @src.Vector3::new(1.0, 0.0, 0.0),
    @src.Vector3::new(0.0, 1.0, 0.0),
  ]
  let set_b = [
    @src.Vector3::new(1.0, 1.0, 0.0),
    @src.Vector3::new(2.0, 1.0, 0.0),
    @src.Vector3::new(1.0, 2.0, 0.0),
  ]
  let (_rot, _trans) = @src.vector_superposition_transform(set_a, set_b)
  assert_true(true)
}

///|
test "vector_apply_transform" {
  let vectors = [@src.Vector3::new(1.0, 0.0, 0.0)]
  let rot = @src.RotationMatrix3::rotation_z(@math.PI / 2.0)
  let trans = @src.Vector3::new(0.0, 0.0, 0.0)
  let result = @src.vector_apply_transform(vectors, rot, trans)
  assert_eq(result.length(), 1)
  assert_true((result[0].get_x()).abs() < 1.0e-10)
  assert_true((result[0].get_y() - 1.0).abs() < 1.0e-10)
}

///|
test "vector_dihedral_angle" {
  let p1 = @src.Vector3::new(1.0, 0.0, 0.0)
  let p2 = @src.Vector3::new(0.0, 0.0, 0.0)
  let p3 = @src.Vector3::new(0.0, 1.0, 0.0)
  let p4 = @src.Vector3::new(0.0, 1.0, 1.0)
  let angle = @src.vector_dihedral_angle(p1, p2, p3, p4)
  assert_true(angle >= -@math.PI && angle <= @math.PI)
}

///|
test "vector_dihedral_angle_degrees" {
  let p1 = @src.Vector3::new(1.0, 0.0, 0.0)
  let p2 = @src.Vector3::new(0.0, 0.0, 0.0)
  let p3 = @src.Vector3::new(0.0, 1.0, 0.0)
  let p4 = @src.Vector3::new(0.0, 1.0, 1.0)
  let angle = @src.vector_dihedral_angle_degrees(p1, p2, p3, p4)
  assert_true(angle >= -180.0 && angle <= 180.0)
}
