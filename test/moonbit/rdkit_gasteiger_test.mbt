///|
/// Regression tests for rdkit_gasteiger.mbt. The adjacency list used
/// to be built with Array::make(n, []), which shares one inner array
/// across all atoms — every atom then saw every bond's neighbors and
/// the charges collapsed. These tests pin the fixed behavior.
test "gasteiger charges directional C-O" {
  let mol = @src.parse_smiles("CO")
  let q = @src.rdkit_gasteiger_charges(mol)
  assert_eq(q.length(), 2)
  // oxygen pulls density from carbon
  assert_true(q[1] < 0.0)
  assert_true(q[0] > 0.0)
}

///|
test "gasteiger charges distinguish atoms on a chain" {
  // With the shared-adjacency bug all atoms received the same
  // neighbor list and the same charge; on a C-C-O chain the terminal
  // carbons must differ clearly from the oxygen.
  let mol = @src.parse_smiles("CCO")
  let q = @src.rdkit_gasteiger_charges(mol)
  assert_eq(q.length(), 3)
  assert_true(q[2] < 0.0)
  assert_true(q[0] >= 0.0)
  assert_true(q[0] != q[2])
  assert_true(q[1] != q[2])
}

///|
test "gasteiger charges conserve neutral molecule" {
  let mol = @src.parse_smiles("CCO")
  let q = @src.rdkit_gasteiger_charges(mol)
  let mut total = 0.0
  for v in q {
    total = total + v
  }
  assert_true(@moonbit.assert_num_close(total, 0.0, 1.0e-9))
}

///|
test "gasteiger atom charge helper" {
  let mol = @src.parse_smiles("CO")
  let q_o = @src.rdkit_atom_charge(mol, 1)
  assert_true(q_o < 0.0)
}
