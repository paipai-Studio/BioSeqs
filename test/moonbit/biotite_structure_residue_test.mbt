///|
/// Tests for biotite_structure_residue.mbt.
fn bsres_atom(name : String, res_name : String, res_id : Int32) -> @src.Atom {
  @src.Atom::new(
    name~,
    coord=@src.Vector3::new(0.0, 0.0, 0.0),
    resname=res_name,
    chainid='A',
    resseq=res_id,
    element="C",
  )
}

///|
test "biotite_get_residues" {
  let arr = @src.AtomArray::from_atoms([
    bsres_atom("N", "ALA", 1),
    bsres_atom("CA", "ALA", 1),
    bsres_atom("CA", "GLY", 2),
  ])
  let residues = @src.biotite_get_residues(arr)
  assert_eq(residues.length(), 2)
  assert_eq(residues[0].0, "A")
  assert_eq(residues[0].1, 1)
  assert_eq(residues[0].2, "ALA")
  assert_eq(residues[0].3, 0)
  assert_eq(residues[1].3, 2)
}

///|
test "biotite_residue_iter_and_count" {
  let arr = @src.AtomArray::from_atoms([
    bsres_atom("N", "ALA", 1),
    bsres_atom("CA", "GLY", 2),
    bsres_atom("CA", "SER", 3),
  ])
  let starts = @src.biotite_residue_iter(arr)
  assert_eq(starts.length(), 3)
  assert_eq(starts[0], 0)
  assert_eq(starts[1], 1)
  assert_eq(@src.biotite_get_residue_count(arr), 3)
}

///|
test "biotite_residue_masses" {
  let arr = @src.AtomArray::from_atoms([
    bsres_atom("CA", "ALA", 1),
    bsres_atom("CA", "GLY", 2),
  ])
  let masses = @src.biotite_residue_masses(arr)
  assert_eq(masses.length(), 2)
  assert_true(@moonbit.assert_num_close(masses[0], 89.09, 1.0e-9))
  assert_true(@moonbit.assert_num_close(masses[1], 75.07, 1.0e-9))
}
