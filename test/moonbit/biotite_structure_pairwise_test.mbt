///|
/// Tests for biotite_structure_pairwise.mbt.
fn bsp_atom(
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
test "biotite_shared_residues_and_atoms" {
  let a = @src.AtomArray::from_atoms([
    bsp_atom("CA", "ALA", 1, 0.0),
    bsp_atom("CA", "GLY", 2, 1.0),
    bsp_atom("CA", "SER", 3, 2.0),
  ])
  let b = @src.AtomArray::from_atoms([
    bsp_atom("CA", "GLY", 2, 1.0),
    bsp_atom("CA", "SER", 3, 2.0),
    bsp_atom("CA", "THR", 4, 3.0),
  ])
  assert_eq(@src.biotite_find_shared_residues(a, b).length(), 2)
  assert_eq(@src.biotite_find_shared_atoms(a, b).length(), 2)
}

///|
test "biotite_shared_atoms_distinguish_atom_names" {
  let a = @src.AtomArray::from_atoms([
    bsp_atom("CA", "ALA", 1, 0.0),
    bsp_atom("CB", "ALA", 1, 0.0),
  ])
  let b = @src.AtomArray::from_atoms([bsp_atom("CB", "ALA", 1, 0.0)])
  assert_eq(@src.biotite_find_shared_atoms(a, b).length(), 1)
}

///|
test "biotite_connected_path" {
  let bonds : Array[@src.Biotite3Bond] = [
    { atom1: 0, atom2: 1, bond_type: 1, },
    { atom1: 1, atom2: 2, bond_type: 1, },
    { atom1: 2, atom2: 3, bond_type: 1, },
  ]
  let path = @src.biotite_find_connected_path(bonds, 4, 0, 3)
  assert_eq(path.length(), 4)
  assert_eq(path[0], 0)
  assert_eq(path[3], 3)
  assert_eq(@src.biotite_find_connected_path(bonds, 4, 0, 9).length(), 0)
}

///|
test "biotite_connected_path_disconnected" {
  let bonds : Array[@src.Biotite3Bond] = [
    { atom1: 0, atom2: 1, bond_type: 1, },
    { atom1: 2, atom2: 3, bond_type: 1, },
  ]
  assert_eq(@src.biotite_find_connected_path(bonds, 4, 0, 3).length(), 0)
}
