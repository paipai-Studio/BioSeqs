test "pdbqt_parse_basic" {
  let content = "ATOM      1  N   ALA A   1       0.000   0.000   0.000  0.000  N\nATOM      2  CA  ALA A   1       1.500   0.000   0.000  0.000  C\nEND\n"
  let atoms = @src.pdbqt_parse(content)
  assert_eq(atoms.length(), 2)
  assert_eq(atoms[0].atom_name, "N")
  assert_eq(atoms[0].res_name, "ALA")
  assert_eq(atoms[1].atom_name, "CA")
}

test "pdbqt_parse_charge_type" {
  let content = "ATOM      1  CA  ALA A   1       1.000   2.000   3.000 -0.100  C\n"
  let atoms = @src.pdbqt_parse(content)
  assert_eq(atoms.length(), 1)
  assert_eq(atoms[0].atom_type, "C")
  assert_eq(atoms[0].partial_charge, -0.1)
}

test "pdbqt_parse_hetatm" {
  let content = "HETATM    1  O   HOH A   1       5.000   5.000   5.000 -0.500  OA\n"
  let atoms = @src.pdbqt_parse(content)
  assert_eq(atoms.length(), 1)
  assert_true(atoms[0].is_hetatm)
  assert_eq(atoms[0].atom_type, "OA")
}

test "pdbqt_atom_type_mapping" {
  assert_eq(@src.pdbqt_atom_type("C", false), "C")
  assert_eq(@src.pdbqt_atom_type("C", true), "A")
  assert_eq(@src.pdbqt_atom_type("O", false), "OA")
  assert_eq(@src.pdbqt_atom_type("H", false), "HD")
  assert_eq(@src.pdbqt_atom_type("N", false), "NA")
  assert_eq(@src.pdbqt_atom_type("S", false), "SA")
}
