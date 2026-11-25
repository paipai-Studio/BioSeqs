///|
/// Tests for biotite_structure_info.mbt.
test "biotite_residue_code_conversion" {
  assert_eq(@src.biotite_one_letter_code("ALA"), "A")
  assert_eq(@src.biotite_one_letter_code("TRP"), "W")
  assert_eq(@src.biotite_one_letter_code("UNK"), "X")
  assert_eq(@src.biotite_three_letter_code("A"), "ALA")
  assert_eq(@src.biotite_three_letter_code("W"), "TRP")
  assert_eq(@src.biotite_three_letter_code("X"), "UNK")
}

///|
test "biotite_residue_and_atom_mass" {
  assert_true(
    @moonbit.assert_num_close(@src.biotite_residue_mass("GLY"), 75.07, 1.0e-9),
  )
  assert_true(
    @moonbit.assert_num_close(@src.biotite_residue_mass("TRP"), 204.23, 1.0e-9),
  )
  assert_eq(@src.biotite_residue_mass("UNK"), 0.0)
  assert_true(
    @moonbit.assert_num_close(@src.biotite_atom_mass("C"), 12.011, 1.0e-9),
  )
  assert_true(
    @moonbit.assert_num_close(@src.biotite_atom_mass("cl"), 35.45, 1.0e-9),
  )
  assert_eq(@src.biotite_atom_mass("XX"), 0.0)
}

///|
test "biotite_residue_charge" {
  assert_eq(@src.biotite_residue_charge("ASP"), -1)
  assert_eq(@src.biotite_residue_charge("GLU"), -1)
  assert_eq(@src.biotite_residue_charge("LYS"), 1)
  assert_eq(@src.biotite_residue_charge("ALA"), 0)
}

///|
test "biotite_structure_mass" {
  let arr = @src.AtomArray::from_atoms([
    @src.Atom::new(
      name="C1",
      coord=@src.Vector3::new(0.0, 0.0, 0.0),
      element="C",
    ),
    @src.Atom::new(
      name="O1",
      coord=@src.Vector3::new(1.0, 0.0, 0.0),
      element="O",
    ),
  ])
  let mass = @src.biotite_structure_mass(arr)
  assert_true(@moonbit.assert_num_close(mass, 28.01, 1.0e-6))
}
