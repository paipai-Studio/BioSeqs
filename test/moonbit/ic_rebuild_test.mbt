///|
/// Tests for ic_rebuild module (Bio.PDB.ic_rebuild - NeRF algorithm).

test "ic_point_distance" {
  let a = @src.IcPoint::new(0.0, 0.0, 0.0)
  let b = @src.IcPoint::new(3.0, 4.0, 0.0)
  let d = @src.ic_point_distance(a, b)
  assert_true((d - 5.0).abs() < 1.0e-9)
}

test "ic_point_cross" {
  let a = @src.IcPoint::new(1.0, 0.0, 0.0)
  let b = @src.IcPoint::new(0.0, 1.0, 0.0)
  let c = @src.ic_point_cross(a, b)
  // (1,0,0) × (0,1,0) = (0,0,1)
  assert_true(c.z.abs() - 1.0 < 1.0e-9)
  assert_true(c.x.abs() < 1.0e-9)
  assert_true(c.y.abs() < 1.0e-9)
}

test "ic_point_normalize" {
  let a = @src.IcPoint::new(3.0, 0.0, 0.0)
  let n = @src.ic_point_normalize(a)
  assert_true((n.x - 1.0).abs() < 1.0e-9)
  assert_true(n.y.abs() < 1.0e-9)
  assert_true(n.z.abs() < 1.0e-9)
}

test "ic_nerf collinear" {
  // Non-collinear reference: A=(0,0,0), B=(1,0,0), C=(2,1,0)
  // bond_length=1, bond_angle=pi (180°, collinear), dihedral=0
  // D should extend the B-C bond
  let a = @src.IcPoint::new(0.0, 0.0, 0.0)
  let b = @src.IcPoint::new(1.0, 0.0, 0.0)
  let c = @src.IcPoint::new(2.0, 1.0, 0.0)
  let d = @src.ic_nerf(a, b, c, 1.0, 3.14159265358979, 0.0)
  // |CD| = 1
  let dist = @src.ic_point_distance(c, d)
  assert_true((dist - 1.0).abs() < 1.0e-6)
}

test "ic_nerf right_angle" {
  // Non-collinear reference: A=(0,0,0), B=(1,0,0), C=(1,1,0)
  // bond_length=1, bond_angle=pi/2 (90°), dihedral=0
  let a = @src.IcPoint::new(0.0, 0.0, 0.0)
  let b = @src.IcPoint::new(1.0, 0.0, 0.0)
  let c = @src.IcPoint::new(1.0, 1.0, 0.0)
  let d = @src.ic_nerf(a, b, c, 1.0, 1.57079632679490, 0.0)
  // |CD| should be 1 (bond_length)
  let dist = @src.ic_point_distance(c, d)
  assert_true((dist - 1.0).abs() < 1.0e-6)
}

test "ic_nerf right_angle" {
  // Non-collinear reference: A=(0,0,0), B=(1,0,0), C=(1,1,0)
  // bond_length=1, bond_angle=pi/2 (90°), dihedral=0
  let a = @src.IcPoint::new(0.0, 0.0, 0.0)
  let b = @src.IcPoint::new(1.0, 0.0, 0.0)
  let c = @src.IcPoint::new(1.0, 1.0, 0.0)
  let d = @src.ic_nerf(a, b, c, 1.0, 1.57079632679490, 0.0)
  // D = C + d2.x*bc_unit + d2.y*m + d2.z*n
  // bc_unit=(0,1,0), n=(0,0,1), m=(-1,0,0)
  // d2=(cos(pi/2), cos(0)*sin(pi/2), sin(0)*sin(pi/2)) ≈ (0, 1, 0)
  // D = (1,1,0) + 0*(0,1,0) + 1*(-1,0,0) + 0*(0,0,1) = (0, 1, 0)
  let dist = @src.ic_point_distance(c, d)
  // |CD| should be 1 (bond_length)
  assert_true((dist - 1.0).abs() < 1.0e-6)
}

test "ic_compute_bond_length" {
  let a = @src.IcPoint::new(0.0, 0.0, 0.0)
  let b = @src.IcPoint::new(0.0, 0.0, 1.5)
  let bl = @src.ic_compute_bond_length(a, b)
  assert_true((bl - 1.5).abs() < 1.0e-9)
}

test "ic_compute_bond_angle" {
  // a=(0,0,0), b=(1,0,0), c=(1,1,0): angle at b = 90°
  let a = @src.IcPoint::new(0.0, 0.0, 0.0)
  let b = @src.IcPoint::new(1.0, 0.0, 0.0)
  let c = @src.IcPoint::new(1.0, 1.0, 0.0)
  let ba = @src.ic_compute_bond_angle(a, b, c)
  let pi_half = 1.57079632679490
  assert_true((ba - pi_half).abs() < 1.0e-6)
}

test "ic_compute_dihedral planar" {
  // Four atoms in a zigzag planar arrangement
  let a = @src.IcPoint::new(0.0, 0.0, 0.0)
  let b = @src.IcPoint::new(1.0, 0.0, 0.0)
  let c = @src.IcPoint::new(1.5, 1.0, 0.0)
  let d = @src.IcPoint::new(2.5, 1.0, 0.0)
  let dh = @src.ic_compute_dihedral(a, b, c, d)
  // All in plane: dihedral = 0 or pi
  assert_true(dh.abs() < 1.0e-6 || (dh.abs() - 3.14159265358979).abs() < 1.0e-6)
}

test "ic_compute_dihedral perpendicular" {
  // a=(0,0,0), b=(1,0,0), c=(1,1,0), d=(1,1,1): dihedral = 90°
  let a = @src.IcPoint::new(0.0, 0.0, 0.0)
  let b = @src.IcPoint::new(1.0, 0.0, 0.0)
  let c = @src.IcPoint::new(1.0, 1.0, 0.0)
  let d = @src.IcPoint::new(1.0, 1.0, 1.0)
  let dh = @src.ic_compute_dihedral(a, b, c, d)
  let pi_half = 1.57079632679490
  assert_true((dh - pi_half).abs() < 1.0e-6 || (dh + pi_half).abs() < 1.0e-6)
}

test "ic_rebuild_chain basic" {
  let p0 = @src.IcPoint::new(0.0, 0.0, 0.0)
  let p1 = @src.IcPoint::new(1.5, 0.0, 0.0)
  let p2 = @src.IcPoint::new(2.5, 1.0, 0.0)
  let ics : Array[@src.IcInternalCoord] = [
    @src.IcInternalCoord::{ atom_name: "D", bond_length: 1.46, bond_angle: 1.91, dihedral: -0.52 },
  ]
  let result = @src.ic_rebuild_chain([p0, p1, p2], ics)
  assert_eq(result.length(), 4)
  assert_eq(result[3].name, "D")
}

test "ic_verify_roundtrip" {
  // Build a chain of 6 atoms with known coordinates (non-degenerate)
  let positions : Array[@src.IcPoint] = [
    @src.IcPoint::new(0.0, 0.0, 0.0),
    @src.IcPoint::new(1.5, 0.0, 0.0),
    @src.IcPoint::new(2.5, 1.0, 0.0),
    @src.IcPoint::new(3.0, 2.0, 0.5),
    @src.IcPoint::new(2.5, 3.0, 1.0),
    @src.IcPoint::new(1.5, 3.5, 0.8),
  ]
  // Verify each internal coord extraction + NeRF rebuild for atom 3
  let bl = @src.ic_compute_bond_length(positions[2], positions[3])
  let ba = @src.ic_compute_bond_angle(positions[1], positions[2], positions[3])
  let dh = @src.ic_compute_dihedral(positions[0], positions[1], positions[2], positions[3])
  let d = @src.ic_nerf(positions[0], positions[1], positions[2], bl, ba, dh)
  let dev = @src.ic_point_distance(positions[3], d)
  // Single atom rebuild should be near-exact
  assert_true(dev < 1.0e-6)
  // Full round-trip
  let full_dev = @src.ic_verify_roundtrip(positions)
  assert_true(full_dev < 1.0e-6)
}

test "ic_rebuild_example" {
  let result = @src.ic_rebuild_example()
  assert_eq(result.length(), 6)
  assert_eq(result[3].name, "CA")
  assert_eq(result[4].name, "CB")
  assert_eq(result[5].name, "CG")
}