///|
test "calculate_hse_basic" {
  let ca = @src.PDBAtom::new("CA", 0.0, 0.0, 0.0)
  let cb = @src.PDBAtom::new("CB", 1.0, 0.0, 0.0)
  let n = @src.PDBAtom::new("N", 0.0, 1.0, 0.0)
  
  let all_ca : Array[@src.PDBAtom] = Array::new()
  all_ca.push(ca)
  all_ca.push(@src.PDBAtom::new("CA", 5.0, 0.0, 0.0))
  all_ca.push(@src.PDBAtom::new("CA", 0.0, 5.0, 0.0))
  all_ca.push(@src.PDBAtom::new("CA", 0.0, 0.0, 5.0))
  
  let result = @src.calculate_hse(ca, cb, n, all_ca)
  assert_eq(result.residue_name, "CA")
  assert_true(result.hse_up >= 0.0)
  assert_true(result.hse_down >= 0.0)
  assert_true(result.hse_total >= 0.0)
}

///|
test "calculate_hse_single_atom" {
  let ca = @src.PDBAtom::new("CA", 0.0, 0.0, 0.0)
  let cb = @src.PDBAtom::new("CB", 1.0, 0.0, 0.0)
  let n = @src.PDBAtom::new("N", 0.0, 1.0, 0.0)
  
  let all_ca : Array[@src.PDBAtom] = Array::new()
  all_ca.push(ca)
  
  let result = @src.calculate_hse(ca, cb, n, all_ca)
  assert_eq(result.residue_name, "CA")
}

///|
test "create_example_atoms" {
  let atoms = @src.create_example_atoms()
  assert_eq(atoms.length(), 9)
}

///|
test "PDBAtom_new" {
  let atom = @src.PDBAtom::new("CA", 1.0, 2.0, 3.0)
  assert_eq(atom.name, "CA")
  assert_eq(atom.x, 1.0)
  assert_eq(atom.y, 2.0)
  assert_eq(atom.z, 3.0)
}