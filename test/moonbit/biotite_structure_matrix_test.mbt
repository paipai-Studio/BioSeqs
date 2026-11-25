///|
/// Tests for biotite_structure_matrix.mbt.
fn bsm_atom(
  name : String,
  res_name : String,
  res_id : Int32,
  x : Double,
  y : Double,
  z : Double,
) -> @src.Atom {
  @src.Atom::new(
    name~,
    coord=@src.Vector3::new(x, y, z),
    resname=res_name,
    chainid='A',
    resseq=res_id,
    element="C",
  )
}

///|
test "biotite_distance_matrix" {
  let arr = @src.AtomArray::from_atoms([
    bsm_atom("C1", "LIG", 1, 0.0, 0.0, 0.0),
    bsm_atom("C2", "LIG", 1, 3.0, 4.0, 0.0),
    bsm_atom("C3", "LIG", 1, 0.0, 0.0, 5.0),
  ])
  let matrix = @src.biotite_get_distance_matrix(arr)
  assert_eq(matrix.length(), 3)
  assert_true(@moonbit.assert_num_close(matrix[0][0], 0.0, 1.0e-12))
  assert_true(@moonbit.assert_num_close(matrix[0][1], 5.0, 1.0e-12))
  assert_true(@moonbit.assert_num_close(matrix[0][2], 5.0, 1.0e-12))
}

///|
test "biotite_angle_and_dihedral" {
  let arr = @src.AtomArray::from_atoms([
    bsm_atom("A", "LIG", 1, 1.0, 0.0, 0.0),
    bsm_atom("B", "LIG", 1, 0.0, 0.0, 0.0),
    bsm_atom("C", "LIG", 1, 0.0, 1.0, 0.0),
    bsm_atom("D", "LIG", 1, 1.0, 1.0, 0.0),
  ])
  let angles = @src.biotite_get_angle_matrix(arr, [(0, 1, 2)])
  assert_eq(angles.length(), 1)
  assert_true(@moonbit.assert_num_close(angles[0], 90.0, 1.0e-6))
  let dihedrals = @src.biotite_get_dihedral_matrix(arr, [[0, 1, 2, 3]])
  assert_eq(dihedrals.length(), 1)
  assert_true(@moonbit.assert_num_close(dihedrals[0], 0.0, 1.0e-6))
}

///|
test "biotite_sequence_identity" {
  let a = @src.AtomArray::from_atoms([
    bsm_atom("CA", "ALA", 1, 0.0, 0.0, 0.0),
    bsm_atom("CA", "GLY", 2, 1.0, 0.0, 0.0),
  ])
  let b = @src.AtomArray::from_atoms([
    bsm_atom("CA", "ALA", 1, 0.0, 0.0, 0.0),
    bsm_atom("CA", "GLY", 2, 1.0, 0.0, 0.0),
  ])
  assert_true(
    @moonbit.assert_num_close(
      @src.biotite_get_sequence_identity(a, b),
      1.0,
      1.0e-12,
    ),
  )
  let c = @src.AtomArray::from_atoms([
    bsm_atom("CA", "ALA", 1, 0.0, 0.0, 0.0),
    bsm_atom("CA", "VAL", 2, 1.0, 0.0, 0.0),
  ])
  assert_true(
    @moonbit.assert_num_close(
      @src.biotite_get_sequence_identity(a, c),
      0.5,
      1.0e-12,
    ),
  )
}
