///|
/// Tests for biotite_structure_integrity.mbt.
fn bsi_atom(
  name : String,
  res_name : String,
  res_id : Int32,
  x : Double,
) -> @src.Atom {
  @src.Atom::new(
    name~,
    coord=@src.Vector3::new(x, 0.0, 0.0),
    resname=res_name,
    chainid='A',
    resseq=res_id,
    element="C",
  )
}

///|
test "biotite_check_residues_complete" {
  let arr = @src.AtomArray::from_atoms([
    bsi_atom("N", "ALA", 1, 0.0),
    bsi_atom("CA", "ALA", 1, 1.0),
    bsi_atom("C", "ALA", 1, 2.0),
    bsi_atom("O", "ALA", 1, 3.0),
  ])
  assert_eq(@src.biotite_check_residues(arr, ["N", "CA", "C", "O"]).length(), 0)
}

///|
test "biotite_check_residues_missing" {
  let arr = @src.AtomArray::from_atoms([
    bsi_atom("N", "ALA", 1, 0.0),
    bsi_atom("CA", "ALA", 1, 1.0),
    bsi_atom("C", "ALA", 1, 2.0),
    bsi_atom("O", "ALA", 1, 3.0),
  ])
  let missing = @src.biotite_check_residues(arr, ["N", "CA", "C", "O", "CB"])
  assert_eq(missing.length(), 1)
  assert_true(missing[0].contains("ALA"))
  assert_true(missing[0].contains("CB"))
}

///|
test "biotite_check_double_bonds_default" {
  let arr = @src.AtomArray::from_atoms([
    bsi_atom("C1", "LIG", 1, 0.0),
    bsi_atom("C2", "LIG", 1, 1.5),
  ])
  assert_eq(@src.biotite_check_double_bonds(arr).length(), 0)
}

///|
test "biotite_check_charges" {
  let atoms : Array[@src.BiotiteFormalChargeAtom] = [
    { element: "C", bond_order_sum: 4.0, n_lone_pairs: -1, },
    { element: "O", bond_order_sum: 2.0, n_lone_pairs: 2, },
  ]
  assert_true(@src.biotite_check_charges(atoms, 0))
  assert_true(!@src.biotite_check_charges(atoms, 1))
}
