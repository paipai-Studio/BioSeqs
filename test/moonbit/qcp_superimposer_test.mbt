///|
test "qcp_superimposer_atom_coordinate" {
  let coord = @src.QCPAtomCoordinate::new(1.0, 2.0, 3.0)

  assert_eq(coord.x, 1.0)
  assert_eq(coord.y, 2.0)
  assert_eq(coord.z, 3.0)
}

///|
test "qcp_superimposer_atom_distance" {
  let coord1 = @src.QCPAtomCoordinate::new(0.0, 0.0, 0.0)
  let coord2 = @src.QCPAtomCoordinate::new(3.0, 4.0, 0.0)

  let dist = coord1.distance(coord2)
  assert_eq(dist, 5.0)
}

///|
test "qcp_superimposer_create" {
  let superimposer = @src.QCPSuperimposer::new()

  assert_eq(superimposer.get_rmsd(), 0.0)
}

///|
test "qcp_superimposer_identical_coords" {
  let atoms1 = [
    @src.QCPAtomCoordinate::new(1.0, 0.0, 0.0),
    @src.QCPAtomCoordinate::new(0.0, 1.0, 0.0),
    @src.QCPAtomCoordinate::new(0.0, 0.0, 1.0),
  ]
  let atoms2 = atoms1.copy()

  let superimposer = @src.bio_qcp_superimpose(atoms1, atoms2)

  assert_true(superimposer.get_rmsd() < 0.01)
}

///|
test "qcp_superimposer_translated_coords" {
  let atoms1 = [
    @src.QCPAtomCoordinate::new(1.0, 0.0, 0.0),
    @src.QCPAtomCoordinate::new(0.0, 1.0, 0.0),
    @src.QCPAtomCoordinate::new(0.0, 0.0, 1.0),
  ]
  let atoms2 = [
    @src.QCPAtomCoordinate::new(2.0, 1.0, 1.0),
    @src.QCPAtomCoordinate::new(1.0, 2.0, 1.0),
    @src.QCPAtomCoordinate::new(1.0, 1.0, 2.0),
  ]

  let superimposer = @src.bio_qcp_superimpose(atoms1, atoms2)

  assert_true(superimposer.get_rmsd() < 0.01)
}

///|
test "qcp_superimposer_calculate_rmsd" {
  let atoms1 = [
    @src.QCPAtomCoordinate::new(0.0, 0.0, 0.0),
    @src.QCPAtomCoordinate::new(1.0, 0.0, 0.0),
    @src.QCPAtomCoordinate::new(0.0, 1.0, 0.0),
  ]
  let atoms2 = [
    @src.QCPAtomCoordinate::new(0.0, 0.0, 0.0),
    @src.QCPAtomCoordinate::new(2.0, 0.0, 0.0),
    @src.QCPAtomCoordinate::new(0.0, 2.0, 0.0),
  ]

  let rmsd = @src.bio_qcp_calculate_rmsd(atoms1, atoms2)

  assert_true(rmsd > 0.0)
}

///|
test "qcp_superimposer_apply_transform" {
  let atoms = [
    @src.QCPAtomCoordinate::new(1.0, 0.0, 0.0),
    @src.QCPAtomCoordinate::new(0.0, 1.0, 0.0),
  ]
  let rot = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]
  let tran = [1.0, 1.0, 0.0]

  let transformed = @src.bio_qcp_apply_transform(atoms, rot, tran)

  assert_eq(transformed.length(), 2)
  assert_true(transformed[0].x > 1.9)
}