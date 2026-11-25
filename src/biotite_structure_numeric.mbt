///|
/// biotite_structure_numeric.mbt — vector and dihedral helpers.

///|
/// Dot product of two vectors.
pub fn biotite_vector_dot(a : Vector3, b : Vector3) -> Double {
  a.dot(b)
}

///|
/// Angle between two vectors in degrees.
pub fn biotite_vector_angle(a : Vector3, b : Vector3) -> Double {
  a.angle_degrees(b)
}

///|
/// Euclidean distance between two points.
pub fn biotite_vector_distance(a : Vector3, b : Vector3) -> Double {
  a.distance(b)
}

///|
/// Euclidean norm of a vector.
pub fn biotite_vector_norm(a : Vector3) -> Double {
  a.norm()
}

///|
/// Dihedral angle in degrees for a four-point torsion.
pub fn biotite_vector_dihedral(
  p1 : Vector3,
  p2 : Vector3,
  p3 : Vector3,
  p4 : Vector3,
) -> Double {
  biotite_dihedral_degrees(p1, p2, p3, p4)
}
