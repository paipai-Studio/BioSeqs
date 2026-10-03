///|
/// Tests for Mol2/SDF I/O.

test "parse_mol2_basic" {
  let content = "@<TRIPOS>MOLECULE\ntest_mol\n   3  2 0 0 0\nSMALL\nNO_CHARGES\n@<TRIPOS>ATOM\n   1 C1    0.0000    0.0000    0.0000 C.3       1 LIG     0.0000\n   2 C2    1.5400    0.0000    0.0000 C.3       1 LIG     0.0000\n   3 O1    2.5400    0.0000    0.0000 O.3       1 LIG     0.0000\n@<TRIPOS>BOND\n   1   1   2 1\n   2   2   3 1\n"
  let mols = @src.parse_mol2(content)
  assert_eq(mols.length(), 1)
  assert_eq(@src.Mol2Molecule::get_name(mols[0]), "test_mol")
  assert_eq(@src.Mol2Molecule::num_atoms(mols[0]), 3)
  assert_eq(@src.Mol2Molecule::num_bonds(mols[0]), 2)
}

test "parse_mol2_atoms" {
  let content = "@<TRIPOS>MOLECULE\nmol\n   2  0 0 0 0\nSMALL\nNO_CHARGES\n@<TRIPOS>ATOM\n   1 N1    1.0000    2.0000    3.0000 N.2       1 RES     0.5000\n   2 C1    4.0000    5.0000    6.0000 C.ar      1 RES     0.0000\n"
  let mols = @src.parse_mol2(content)
  let atoms = @src.Mol2Molecule::get_atoms(mols[0])
  assert_eq(@src.Mol2Atom::get_serial(atoms[0]), 1)
  assert_eq(@src.Mol2Atom::get_name(atoms[0]), "N1")
  assert_true((@src.Mol2Atom::get_x(atoms[0]) - 1.0).abs() < 0.001)
  assert_true((@src.Mol2Atom::get_y(atoms[0]) - 2.0).abs() < 0.001)
  assert_eq(@src.Mol2Atom::get_sybyl_type(atoms[0]), "N.2")
  assert_eq(@src.Mol2Atom::get_subst_name(atoms[0]), "RES")
  assert_true((@src.Mol2Atom::get_charge(atoms[0]) - 0.5).abs() < 0.001)
}

test "parse_mol2_bonds" {
  let content = "@<TRIPOS>MOLECULE\nmol\n   3  2 0 0 0\nSMALL\nNO_CHARGES\n@<TRIPOS>ATOM\n   1 C1    0.0000    0.0000    0.0000 C.3       1 LIG     0.0000\n   2 C2    1.5400    0.0000    0.0000 C.3       1 LIG     0.0000\n   3 O1    2.5400    0.0000    0.0000 O.3       1 LIG     0.0000\n@<TRIPOS>BOND\n   1   1   2 1\n   2   2   3 ar\n"
  let mols = @src.parse_mol2(content)
  let bonds = @src.Mol2Molecule::get_bonds(mols[0])
  assert_eq(@src.Mol2Bond::get_origin(bonds[0]), 1)
  assert_eq(@src.Mol2Bond::get_target(bonds[0]), 2)
  assert_eq(@src.Mol2Bond::get_bond_type(bonds[0]), "1")
  assert_eq(@src.Mol2Bond::get_bond_type(bonds[1]), "ar")
}

test "parse_mol2_multi" {
  let content = "@<TRIPOS>MOLECULE\nmol1\n   1  0 0 0 0\nSMALL\nNO_CHARGES\n@<TRIPOS>ATOM\n   1 C1    0.0000    0.0000    0.0000 C.3       1 LIG     0.0000\n@<TRIPOS>MOLECULE\nmol2\n   1  0 0 0 0\nSMALL\nNO_CHARGES\n@<TRIPOS>ATOM\n   1 N1    1.0000    0.0000    0.0000 N.3       1 LIG     0.0000\n"
  let mols = @src.parse_mol2(content)
  assert_eq(mols.length(), 2)
  assert_eq(@src.Mol2Molecule::get_name(mols[0]), "mol1")
  assert_eq(@src.Mol2Molecule::get_name(mols[1]), "mol2")
}

test "parse_mol2_empty" {
  assert_eq(@src.parse_mol2("").length(), 0)
}

test "mol2_roundtrip" {
  let content = "@<TRIPOS>MOLECULE\ntest\n   2  1 0 0 0\nSMALL\nNO_CHARGES\n@<TRIPOS>ATOM\n   1 C1    0.0000    0.0000    0.0000 C.3       1 LIG     0.0000\n   2 O1    1.5000    0.0000    0.0000 O.3       1 LIG     0.0000\n@<TRIPOS>BOND\n   1   1   2 2\n"
  let mols = @src.parse_mol2(content)
  let output = @src.mol2_to_string(mols[0])
  let reparsed = @src.parse_mol2(output)
  assert_eq(@src.Mol2Molecule::num_atoms(reparsed[0]), 2)
  assert_eq(@src.Mol2Molecule::num_bonds(reparsed[0]), 1)
  let bonds = @src.Mol2Molecule::get_bonds(reparsed[0])
  assert_eq(@src.Mol2Bond::get_bond_type(bonds[0]), "2")
}

test "mol2_to_chain" {
  let content = "@<TRIPOS>MOLECULE\nprotein\n   3  0 0 0 0\nSMALL\nNO_CHARGES\n@<TRIPOS>ATOM\n   1 CA    0.0000    0.0000    0.0000 C.3       1 ALA     0.0000\n   2 CA    3.8000    0.0000    0.0000 C.3       2 GLY     0.0000\n   3 CA    7.6000    0.0000    0.0000 C.3       3 VAL     0.0000\n"
  let mols = @src.parse_mol2(content)
  let chain = @src.mol2_to_chain(mols[0])
  let residues = chain.get_residues()
  assert_eq(residues.length(), 3)
  assert_eq(residues[0].resname, "ALA")
  assert_eq(residues[1].resname, "GLY")
  assert_eq(residues[2].resname, "VAL")
}

test "mol2_to_chain_grouped" {
  let content = "@<TRIPOS>MOLECULE\nprotein\n   4  0 0 0 0\nSMALL\nNO_CHARGES\n@<TRIPOS>ATOM\n   1 N     0.0000    0.0000    0.0000 N.3       1 ALA     0.0000\n   2 CA    1.5000    0.0000    0.0000 C.3       1 ALA     0.0000\n   3 N     3.8000    0.0000    0.0000 N.3       2 GLY     0.0000\n   4 CA    5.3000    0.0000    0.0000 C.3       2 GLY     0.0000\n"
  let mols = @src.parse_mol2(content)
  let chain = @src.mol2_to_chain(mols[0])
  let residues = chain.get_residues()
  assert_eq(residues.length(), 2)
  assert_eq(residues[0].get_num_atoms(), 2)
  assert_eq(residues[1].get_num_atoms(), 2)
}

test "sdf_to_string_basic" {
  let atoms = [@src.SdfAtom::new("C", 0.0, 0.0, 0.0, 0), @src.SdfAtom::new("O", 1.5, 0.0, 0.0, 0)]
  let bonds = [@src.SdfBond::new(0, 1, 2)]
  let mol = @src.SdfMolecule::new("test_sdf", atoms, bonds)
  let output = @src.sdf_to_string(mol)
  assert_true(output.contains("test_sdf"))
  assert_true(output.contains("M  END"))
  assert_true(output.contains("$$$$"))
}

test "sdf_to_string_charge" {
  let atoms = [@src.SdfAtom::new("N", 0.0, 0.0, 0.0, 1), @src.SdfAtom::new("O", 1.5, 0.0, 0.0, -1)]
  let mol = @src.SdfMolecule::new("charged", atoms, [])
  let output = @src.sdf_to_string(mol)
  assert_true(output.contains("charged"))
  assert_true(output.contains("M  END"))
}

test "sdf_to_string_props" {
  let atoms = [@src.SdfAtom::new("C", 0.0, 0.0, 0.0, 0)]
  let mol = @src.SdfMolecule::new("prop_test", atoms, [])
  mol.properties["logP"] = "2.5"
  mol.properties["SMILES"] = "CCO"
  let output = @src.sdf_to_string(mol)
  assert_true(output.contains("logP"))
  assert_true(output.contains("2.5"))
  assert_true(output.contains("SMILES"))
  assert_true(output.contains("CCO"))
}

test "sdf_roundtrip" {
  let content = "ethanol\n  BioSeqs\n\n  3  2  0  0  0  0  0  0  0  0999 V2000\n    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    1.5400    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n    2.5400    0.0000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0\n  1  2  1  0  0  0  0\n  2  3  1  0  0  0  0\nM  END\n$$$$\n"
  let mols = @src.chem_parse_sdf(content)
  assert_eq(mols.length(), 1)
  let output = @src.sdf_to_string(mols[0])
  let reparsed = @src.chem_parse_sdf(output)
  assert_eq(reparsed.length(), 1)
  assert_eq(reparsed[0].atoms.length(), 3)
  assert_eq(reparsed[0].bonds.length(), 2)
}

test "sdf_to_chain" {
  let atoms = [@src.SdfAtom::new("C", 0.0, 0.0, 0.0, 0), @src.SdfAtom::new("N", 1.5, 0.0, 0.0, 0), @src.SdfAtom::new("O", 3.0, 0.0, 0.0, 0)]
  let mol = @src.SdfMolecule::new("test", atoms, [])
  let chain = @src.sdf_to_chain(mol)
  let residues = chain.get_residues()
  assert_eq(residues.length(), 3)
  assert_eq(residues[0].resname, "C")
  assert_eq(residues[1].resname, "N")
}

test "chain_to_mol2" {
  let residues : Array[@src.Residue] = Array(capacity=0)
  for r = 0; r < 3; r = r + 1 {
    let x = r.to_double() * 3.8
    let ca = @src.Atom::new(name="CA", coord=@src.Vector3::new(x, 0.0, 0.0), resname="ALA", chainid='A', resseq=r + 1, element="C")
    residues.push(@src.Residue::new(resname="ALA", chainid='A', resseq=r + 1, atoms=[ca]))
  }
  let chain = @src.Chain::new(id='A', residues~)
  let mol = @src.chain_to_mol2(chain, "protein")
  assert_eq(@src.Mol2Molecule::get_name(mol), "protein")
  assert_eq(@src.Mol2Molecule::num_atoms(mol), 3)
  let atoms = @src.Mol2Molecule::get_atoms(mol)
  assert_eq(@src.Mol2Atom::get_subst_name(atoms[0]), "ALA")
}

test "chain_to_sdf" {
  let residues : Array[@src.Residue] = Array(capacity=0)
  for r = 0; r < 2; r = r + 1 {
    let x = r.to_double() * 3.8
    let ca = @src.Atom::new(name="CA", coord=@src.Vector3::new(x, 0.0, 0.0), resname="GLY", chainid='A', resseq=r + 1, element="C")
    residues.push(@src.Residue::new(resname="GLY", chainid='A', resseq=r + 1, atoms=[ca]))
  }
  let chain = @src.Chain::new(id='A', residues~)
  let mol = @src.chain_to_sdf(chain, "gly_chain")
  assert_eq(mol.name, "gly_chain")
  assert_eq(mol.atoms.length(), 2)
  assert_eq(mol.atoms[0].symbol, "C")
}

test "mol2_accessors" {
  let content = "@<TRIPOS>MOLECULE\nacc_test\n   2  1 0 0 0\nSMALL\nUSER_CHARGES\n@<TRIPOS>ATOM\n   1 C1    0.0000    0.0000    0.0000 C.3       1 LIG     0.0000\n   2 O1    1.5000    0.0000    0.0000 O.3       1 LIG     0.0000\n@<TRIPOS>BOND\n   1   1   2 1\n"
  let mols = @src.parse_mol2(content)
  let mol = mols[0]
  assert_eq(@src.Mol2Molecule::get_name(mol), "acc_test")
  assert_eq(@src.Mol2Molecule::num_atoms(mol), 2)
  assert_eq(@src.Mol2Molecule::num_bonds(mol), 1)
  assert_eq(@src.Mol2Molecule::get_mol_type(mol), "SMALL")
  assert_eq(@src.Mol2Molecule::get_charge_type(mol), "USER_CHARGES")
}

test "mol2_chain_roundtrip" {
  let residues : Array[@src.Residue] = Array(capacity=0)
  for r = 0; r < 3; r = r + 1 {
    let x = r.to_double() * 3.8
    let ca = @src.Atom::new(name="CA", coord=@src.Vector3::new(x, 0.0, 0.0), resname="ALA", chainid='A', resseq=r + 1, element="C")
    residues.push(@src.Residue::new(resname="ALA", chainid='A', resseq=r + 1, atoms=[ca]))
  }
  let chain = @src.Chain::new(id='A', residues~)
  let mol = @src.chain_to_mol2(chain, "protein")
  let mol2_str = @src.mol2_to_string(mol)
  let reparsed = @src.parse_mol2(mol2_str)
  assert_eq(@src.Mol2Molecule::num_atoms(reparsed[0]), 3)
  let chain2 = @src.mol2_to_chain(reparsed[0])
  assert_eq(chain2.get_residues().length(), 3)
}
