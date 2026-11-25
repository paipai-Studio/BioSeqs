///|
/// Tests for biotite_structure3.mbt — Atom views, bonds, AtomArray adapters.
fn b3_atom(
  name : String,
  res_name : String,
  res_id : Int32,
  chain_id : String,
  x : Double,
) -> @src.Biotite3Atom {
  {
    atom_name: name,
    element: "C",
    res_name,
    res_id,
    chain_id,
    coord: (x, 0.0, 0.0),
    b_factor: 0.0,
    occupancy: 1.0,
    hetero: false,
  }
}

///|
test "biotite3_connect_via_distances" {
  let atoms : Array[@src.Biotite3Atom] = [
    b3_atom("C1", "LIG", 1, "A", 0.0),
    b3_atom("C2", "LIG", 1, "A", 1.5),
    b3_atom("C3", "LIG", 1, "A", 3.0),
    b3_atom("C4", "LIG", 1, "A", 4.5),
  ]
  let bonds = @src.biotite3_connect_via_distances(atoms)
  assert_eq(bonds.length(), 3)
  let connected = @src.biotite3_find_connected(bonds, 4, 0)
  assert_eq(connected.length(), 4)
  let rotatable = @src.biotite3_find_rotatable_bonds(bonds, atoms)
  assert_eq(rotatable.length(), 1)
  let typed = @src.biotite3_infer_bond_types(bonds, atoms)
  assert_eq(typed.length(), 3)
  assert_eq(typed[0].bond_type, 1)
}

///|
test "biotite3_connect_via_residue_names" {
  let atoms : Array[@src.Biotite3Atom] = [
    b3_atom("N", "ALA", 1, "A", 0.0),
    b3_atom("CA", "ALA", 1, "A", 1.5),
    b3_atom("C", "ALA", 1, "A", 3.0),
    b3_atom("O", "ALA", 1, "A", 4.0),
  ]
  let bonds = @src.biotite3_connect_via_residue_names(atoms)
  assert_true(bonds.length() >= 3)
}

///|
test "biotite3_atom_array_adapters" {
  let atoms : Array[@src.Biotite3Atom] = [
    b3_atom("C1", "LIG", 1, "A", 0.0),
    b3_atom("C2", "LIG", 1, "A", 1.5),
  ]
  let array = @src.biotite3_atoms_to_atom_array(atoms)
  assert_eq(array.len(), 2)
  let views = @src.biotite3_atoms_from_atom_array(array)
  assert_eq(views.length(), 2)
  assert_eq(views[0].atom_name, "C1")
  assert_eq(views[0].chain_id, "A")
  let bonds = @src.biotite3_connect_via_distances_array(array)
  assert_eq(bonds.length(), 1)
}

///|
test "biotite3_atom_array_stack" {
  let atoms : Array[@src.Biotite3Atom] = [b3_atom("CA", "ALA", 1, "A", 0.0)]
  let stack = @src.biotite3_atom_array_stack(atoms, [
    [(0.0, 0.0, 0.0)],
    [(1.0, 0.0, 0.0)],
  ])
  assert_eq(stack.n_models, 2)
  assert_eq(@src.biotite3_atom_array_stack_get_model(stack, 0).length(), 1)
  assert_eq(@src.biotite3_atom_array_stack_get_model(stack, 5).length(), 0)
  let combined = @src.biotite3_atom_array_stack_concatenate([stack, stack])
  assert_eq(combined.n_models, 4)
}
