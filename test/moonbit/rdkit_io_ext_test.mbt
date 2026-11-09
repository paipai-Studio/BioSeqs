///|
/// Tests for rdkit_io_ext.mbt — molblock, SDF, RXN, .smi and InChI I/O.

///|
test "molblock_roundtrip_ethanol" {
  let mol = @src.parse_smiles("CCO")
  let block = @src.rdkit_mol_to_molblock(mol, name="ethanol")
  // counts line
  assert_true(block.contains("V2000"))
  assert_true(block.contains("M  END"))
  let back = @src.rdkit_molblock_to_mol(block)
  assert_eq(back.num_atoms(), 3)
  assert_eq(back.num_bonds(), 2)
  assert_eq(back.atom_symbols(), ["C", "C", "O"])
}

///|
test "molblock_roundtrip_kekule_benzene" {
  let mol = @src.parse_smiles("C1=CC=CC=C1")
  let block = @src.rdkit_mol_to_molblock(mol)
  let back = @src.rdkit_molblock_to_mol(block)
  assert_eq(back.num_atoms(), 6)
  assert_eq(back.num_bonds(), 6)
  // aromaticity is re-perceived on parse
  assert_eq(@src.mol_aromatic_ring_count(back), 1)
}

///|
test "molblock_charge_decode" {
  // atom charge codes: 3 -> +1, 5 -> -1 (V2000 column 37-39)
  let block = "mol\n  BioSeqs          2D\n\n" +
    "  2  1  0  0  0  0  0  0  0  0999 V2000\n" +
    "    0.0000    0.0000    0.0000 N   0  3  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0\n" +
    "    0.0000    0.0000    0.0000 O   0  5  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0\n" +
    "  1  2  1  0\nM  END\n"
  let mol = @src.rdkit_molblock_to_mol(block)
  assert_eq(mol.num_atoms(), 2)
  assert_eq(mol.atoms[0].charge, 1)
  assert_eq(mol.atoms[1].charge, -1)
}

///|
test "sdf_roundtrip" {
  let mols = [@src.parse_smiles("CCO"), @src.parse_smiles("CCN")]
  let sdf = @src.rdkit_write_sdf(mols, ["ethanol", "ethylamine"])
  assert_true(sdf.contains("$$$$"))
  let back = @src.rdkit_read_sdf(sdf)
  assert_eq(back.length(), 2)
  assert_eq(back[0].atom_symbols(), ["C", "C", "O"])
  assert_eq(back[1].atom_symbols(), ["C", "C", "N"])
}

///|
test "rxn_roundtrip" {
  let rxn = @src.rdkit_parse_reaction("CCO.CC(=O)O>>CCOC(C)=O")
  let text = @src.rdkit_write_rxn_file(rxn)
  assert_true(text.contains("$RXN"))
  assert_true(text.contains("$MOL"))
  let parsed = @src.rdkit_parse_rxn_file(text)
  assert_eq(parsed.reactants.length(), 2)
  assert_eq(parsed.products.length(), 1)
  // reactant SMILES are canonicalized through the molblock
  let canon_ethanol = @src.mol_to_smiles(@src.parse_smiles("CCO"))
  assert_true(parsed.reactants.contains(canon_ethanol))
}

///|
test "smi_library_roundtrip" {
  let text = "# comment line\nCCO\tethanol\nCCN ethylamine\n\nC"
  let records = @src.rdkit_read_smi(text)
  assert_eq(records.length(), 3)
  assert_eq(records[0].1, "ethanol")
  assert_eq(records[1].1, "ethylamine")
  assert_eq(records[2].2.num_atoms(), 1)
  let out = @src.rdkit_write_smi([("CCO", "ethanol"), ("CCN", "ethylamine")])
  let back = @src.rdkit_read_smi(out)
  assert_eq(back.length(), 2)
  assert_eq(back[0].0, "CCO")
  assert_eq(back[1].1, "ethylamine")
}

///|
test "inchi_to_smiles_ethanol" {
  let smi = @src.rdkit_inchi_to_smiles("InChI=1S/C2H6O/c1-2-3/h3H,2H2,1H3")
  match smi {
    Some(s) => assert_true(s == "CCO" || s == "OCC")
    None => fail("InChI ethanol should parse")
  }
}

///|
test "inchi_to_smiles_benzene_aromatic" {
  let smi = @src.rdkit_inchi_to_smiles("InChI=1S/C6H6/c1-2-4-6-3-5-1/h1-6H")
  match smi {
    Some(s) => {
      // six-membered ring backbone, round-trips through the parser
      let back = @src.parse_smiles(s)
      assert_eq(back.num_atoms(), 6)
      assert_eq(back.num_bonds(), 6)
      assert_eq(@src.mol_ring_count(back), 1)
    }
    None => fail("InChI benzene should parse")
  }
}

///|
test "inchi_to_smiles_branching" {
  // isobutane: central C bonded to three methyls via a branch
  let smi = @src.rdkit_inchi_to_smiles("InChI=1S/C4H10/c1-4(2)3/h4H,1-3H3")
  match smi {
    Some(s) => assert_eq(@src.parse_smiles(s).num_atoms(), 4)
    None => fail("InChI isobutane should parse")
  }
}

///|
test "inchi_to_smiles_rejects_bad_input" {
  // not an InChI
  assert_true(@src.rdkit_inchi_to_smiles("not-an-inchi") is None)
  // out-of-range atom reference in the connection layer
  assert_true(@src.rdkit_inchi_to_smiles("InChI=1S/CH4/c1-99/h1H4") is None)
}
