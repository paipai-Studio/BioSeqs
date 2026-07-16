///|
test "svd_superimposer_atom_coordinate" {
  let coord = @src.AtomCoordinate::new(1.0, 2.0, 3.0)

  assert_eq(coord.x, 1.0)
  assert_eq(coord.y, 2.0)
  assert_eq(coord.z, 3.0)
}

///|
test "svd_superimposer_atom_distance" {
  let coord1 = @src.AtomCoordinate::new(0.0, 0.0, 0.0)
  let coord2 = @src.AtomCoordinate::new(3.0, 4.0, 0.0)

  let dist = coord1.distance(coord2)
  assert_eq(dist, 5.0)
}

///|
test "svd_superimposer_create" {
  let superimposer = @src.SVDSuperimposer::new()

  assert_eq(superimposer.get_rmsd(), 0.0)
}

///|
test "svd_superimposer_identical_coords" {
  let atoms1 = [
    @src.AtomCoordinate::new(1.0, 0.0, 0.0),
    @src.AtomCoordinate::new(0.0, 1.0, 0.0),
    @src.AtomCoordinate::new(0.0, 0.0, 1.0),
  ]
  let atoms2 = atoms1.copy()

  let (transformed, rmsd) = @src.superimpose(atoms1, atoms2)

  assert_eq(rmsd < 0.5, true)
}

///|
test "svd_superimposer_translated_coords" {
  let atoms1 = [
    @src.AtomCoordinate::new(1.0, 0.0, 0.0),
    @src.AtomCoordinate::new(0.0, 1.0, 0.0),
    @src.AtomCoordinate::new(0.0, 0.0, 1.0),
  ]
  let atoms2 = [
    @src.AtomCoordinate::new(3.0, 2.0, 1.0),
    @src.AtomCoordinate::new(2.0, 3.0, 1.0),
    @src.AtomCoordinate::new(2.0, 2.0, 2.0),
  ]

  let (transformed, rmsd) = @src.superimpose(atoms1, atoms2)

  assert_eq(rmsd < 0.5, true)
}

///|
test "svd_superimposer_rotated_coords" {
  let atoms1 = [
    @src.AtomCoordinate::new(1.0, 0.0, 0.0),
    @src.AtomCoordinate::new(0.0, 1.0, 0.0),
    @src.AtomCoordinate::new(0.0, 0.0, 1.0),
  ]
  let atoms2 = [
    @src.AtomCoordinate::new(0.0, 1.0, 0.0),
    @src.AtomCoordinate::new(-1.0, 0.0, 0.0),
    @src.AtomCoordinate::new(0.0, 0.0, 1.0),
  ]

  let (transformed, rmsd) = @src.superimpose(atoms1, atoms2)

  assert_eq(rmsd < 0.5, true)
}

///|
test "svd_superimposer_rmsd_calculation" {
  let atoms1 = [
    @src.AtomCoordinate::new(1.0, 0.0, 0.0),
    @src.AtomCoordinate::new(0.0, 1.0, 0.0),
    @src.AtomCoordinate::new(0.0, 0.0, 1.0),
  ]
  let atoms2 = [
    @src.AtomCoordinate::new(0.0, 1.0, 0.0),
    @src.AtomCoordinate::new(-1.0, 0.0, 0.0),
    @src.AtomCoordinate::new(0.0, 0.0, 1.0),
  ]

  let superimposer = @src.SVDSuperimposer::set(atoms1, atoms2)
  let rot = superimposer.get_rotation()

  assert_eq(rot.length(), 3)
  assert_eq(rot[0].length(), 3)
}

///|
test "svd_superimposer_get_translation" {
  let atoms1 = [
    @src.AtomCoordinate::new(1.0, 0.0, 0.0),
    @src.AtomCoordinate::new(0.0, 1.0, 0.0),
    @src.AtomCoordinate::new(0.0, 0.0, 1.0),
  ]
  let atoms2 = [
    @src.AtomCoordinate::new(3.0, 2.0, 1.0),
    @src.AtomCoordinate::new(2.0, 3.0, 1.0),
    @src.AtomCoordinate::new(2.0, 2.0, 2.0),
  ]

  let superimposer = @src.SVDSuperimposer::set(atoms1, atoms2)
  let tran = superimposer.get_translation()

  assert_eq(tran.length(), 3)
}

///|
test "svd_superimposer_calculate_rmsd" {
  let coords1 = [
    @src.AtomCoordinate::new(0.0, 0.0, 0.0),
    @src.AtomCoordinate::new(1.0, 0.0, 0.0),
    @src.AtomCoordinate::new(0.0, 1.0, 0.0),
  ]
  let coords2 = [
    @src.AtomCoordinate::new(0.0, 0.0, 0.0),
    @src.AtomCoordinate::new(1.0, 0.0, 0.0),
    @src.AtomCoordinate::new(0.0, 1.0, 0.0),
  ]

  let rmsd = @src.calculate_rmsd(coords1, coords2)
  assert_eq(rmsd, 0.0)
}

///|
test "svd_superimposer_apply" {
  let atoms = [
    @src.AtomCoordinate::new(1.0, 0.0, 0.0),
    @src.AtomCoordinate::new(0.0, 1.0, 0.0),
  ]
  let superimposer = @src.SVDSuperimposer::new()

  let transformed = superimposer.apply(atoms)
  assert_eq(transformed.length(), 2)
}
