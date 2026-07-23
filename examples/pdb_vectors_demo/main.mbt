///|
/// PDB Vectors Demo - Demonstrates 3D vector and rotation matrix operations.
/// Inspired by Biopython's Bio.PDB.vectors module.
///
/// Run with: moon run examples/pdb_vectors_demo/main.mbt
fn main {
  println("=== PDB Vectors Demo ===\n")

  // 1. Vector creation
  println("1. Vector creation:")
  let v1 = @bio.Vector3::new(1.0, 2.0, 3.0)
  let v2 = @bio.Vector3::new(4.0, 5.0, 6.0)
  println("   v1: (\{v1.get_x()}, \{v1.get_y()}, \{v1.get_z()})")
  println("   v2: (\{v2.get_x()}, \{v2.get_y()}, \{v2.get_z()})")

  let unit_x = @bio.Vector3::unit_x()
  let unit_y = @bio.Vector3::unit_y()
  let unit_z = @bio.Vector3::unit_z()
  println("   Unit vectors created successfully")

  // 2. Vector arithmetic
  println("\n2. Vector arithmetic:")
  let v_add = v1.add(v2)
  println("   v1 + v2: (\{v_add.get_x()}, \{v_add.get_y()}, \{v_add.get_z()})")

  let v_sub = v1.sub(v2)
  println("   v1 - v2: (\{v_sub.get_x()}, \{v_sub.get_y()}, \{v_sub.get_z()})")

  let v_scale = v1.scale(2.0)
  println("   v1 * 2: (\{v_scale.get_x()}, \{v_scale.get_y()}, \{v_scale.get_z()})")

  // 3. Vector products and norms
  println("\n3. Vector products and norms:")
  let dot = v1.dot(v2)
  println("   v1 . v2 (dot product): \{dot}")

  let cross = v1.cross(v2)
  println("   v1 x v2 (cross product): (\{cross.get_x()}, \{cross.get_y()}, \{cross.get_z()})")

  let norm = v1.norm()
  println("   |v1| (norm): \{norm}")

  let n1 = v1.normalize()
  println("   Normalized v1: (\{n1.get_x()}, \{n1.get_y()}, \{n1.get_z()})")

  // 4. Vector distance and angle
  println("\n4. Distance and angle:")
  let dist_sq = v1.distance_squared(v2)
  println("   Distance^2(v1, v2): \{dist_sq}")

  let angle = v1.angle(v2)
  println("   Angle(v1, v2) [rad]: \{angle}")

  let angle_deg = v1.angle_degrees(v2)
  println("   Angle(v1, v2) [deg]: \{angle_deg}")

  // 5. Rotation matrix
  println("\n5. Rotation matrices:")
  let ident = @bio.RotationMatrix3::identity()
  println("   Identity matrix det: \{ident.determinant()}")

  let rot_x = @bio.RotationMatrix3::rotation_x(1.5708)
  println("   Rotation by 90 around X - det: \{rot_x.determinant()}")
  println("   Is rotation: \{rot_x.is_rotation()}")

  // 6. Rotation axis-angle
  println("\n6. Axis-angle rotation:")
  let axis = @bio.Vector3::new(1.0, 1.0, 0.0)
  let rot_axis = @bio.RotationMatrix3::rotation_axis_angle(axis, 1.5708)
  println("   Rotation around (1,1,0) by 90 - det: \{rot_axis.determinant()}")

  // 7. Transform a vector
  println("\n7. Vector transformation:")
  let transformed = rot_x.transform(v1)
  println("   Rotated v1 90 around X: (\{transformed.get_x()}, \{transformed.get_y()}, \{transformed.get_z()})")

  // 8. Superposition (Kabsch algorithm)
  println("\n8. Structural superposition:")
  let set_a = [
    @bio.Vector3::new(1.0, 0.0, 0.0),
    @bio.Vector3::new(0.0, 1.0, 0.0),
    @bio.Vector3::new(0.0, 0.0, 1.0),
  ]
  let set_b = [
    @bio.Vector3::new(0.0, 1.0, 0.0),
    @bio.Vector3::new(-1.0, 0.0, 0.0),
    @bio.Vector3::new(0.0, 0.0, 1.0),
  ]
  let rmsd = @bio.vector_superposition_rmsd(set_a, set_b)
  println("   Superposition RMSD: \{rmsd}")

  let (rot, trans) = @bio.vector_superposition_transform(set_a, set_b)
  println("   Optimal rotation found, translation: (\{trans.get_x()}, \{trans.get_y()}, \{trans.get_z()})")

  let rotated_set = @bio.vector_apply_transform(set_a, rot, trans)
  println("   Applied transformation to \{rotated_set.length()} vectors")

  // 9. Dihedral angle
  println("\n9. Dihedral angle:")
  let p1 = @bio.Vector3::new(1.0, 0.0, 0.0)
  let p2 = @bio.Vector3::new(0.0, 1.0, 0.0)
  let p3 = @bio.Vector3::new(0.0, 0.0, 1.0)
  let p4 = @bio.Vector3::new(1.0, 1.0, 1.0)
  let dihedral = @bio.vector_dihedral_angle(p1, p2, p3, p4)
  println("   Dihedral angle [rad]: \{dihedral}")
  let dihedral_deg = @bio.vector_dihedral_angle_degrees(p1, p2, p3, p4)
  println("   Dihedral angle [deg]: \{dihedral_deg}")

  // 10. Centroid
  println("\n10. Centroid calculation:")
  let centroid = @bio.vector_centroid(set_a)
  println("   Centroid: (\{centroid.get_x()}, \{centroid.get_y()}, \{centroid.get_z()})")

  println("\n=== Demo Complete ===")
}