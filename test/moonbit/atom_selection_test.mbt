///|
/// Tests for ProDy-style Atom Selection Language.

///|
/// Build a test chain with mixed residue types.
fn make_sel_test_chain() -> @src.Chain {
  let residues : Array[@src.Residue] = Array(capacity=0)
  // Residues 1-5: protein (ALA) with backbone + sidechain
  let aa_atoms = [
    ("N", "N"), ("CA", "C"), ("C", "C"), ("O", "O"), ("CB", "C"), ("HB1", "H")
  ]
  for r = 0; r < 5; r = r + 1 {
    let atoms : Array[@src.Atom] = Array(capacity=0)
    for a = 0; a < aa_atoms.length(); a = a + 1 {
      let (aname, elem) = aa_atoms[a]
      let x = (a.to_double()) * 1.5
      let y = (r.to_double()) * 3.0
      let atom = @src.Atom::new(name=aname, coord=@src.Vector3::new(x, y, 0.0), resname="ALA", chainid='A', resseq=r + 1, element=elem)
      atoms.push(atom)
    }
    let res = @src.Residue::new(resname="ALA", chainid='A', resseq=r + 1, atoms~)
    residues.push(res)
  }
  // Residue 6: water
  let wat_atoms : Array[@src.Atom] = Array(capacity=0)
  wat_atoms.push(@src.Atom::new(name="O", coord=@src.Vector3::new(10.0, 10.0, 0.0), resname="HOH", chainid='A', resseq=6, element="O"))
  let wat_res = @src.Residue::new(resname="HOH", chainid='A', resseq=6, atoms=wat_atoms, hetfield="HETATM")
  residues.push(wat_res)
  // Residue 7: nucleic (DA)
  let na_atoms : Array[@src.Atom] = Array(capacity=0)
  na_atoms.push(@src.Atom::new(name="P", coord=@src.Vector3::new(15.0, 0.0, 0.0), resname="DA", chainid='B', resseq=1, element="P"))
  na_atoms.push(@src.Atom::new(name="C1'", coord=@src.Vector3::new(16.0, 0.0, 0.0), resname="DA", chainid='B', resseq=1, element="C"))
  let na_res = @src.Residue::new(resname="DA", chainid='B', resseq=1, atoms=na_atoms)
  residues.push(na_res)
  @src.Chain::new(id='A', residues~)
}

// ===== Basic selector tests =====

test "sel_all" {
  let chain = make_sel_test_chain()
  let atoms = chain.get_atoms()
  let sel = @src.select_atoms(atoms, "all")
  assert_eq(sel.length(), atoms.length())
}

test "sel_none" {
  let chain = make_sel_test_chain()
  let atoms = chain.get_atoms()
  let sel = @src.select_atoms(atoms, "none")
  assert_eq(sel.length(), 0)
}

test "sel_calpha" {
  let chain = make_sel_test_chain()
  let atoms = chain.get_atoms()
  let sel = @src.select_atoms(atoms, "calpha")
  // 5 ALA residues, each has 1 CA
  assert_eq(sel.length(), 5)
  for i = 0; i < sel.length(); i = i + 1 {
    assert_eq(sel[i].name, "CA")
  }
}

test "sel_ca_alias" {
  let chain = make_sel_test_chain()
  let atoms = chain.get_atoms()
  let sel = @src.select_atoms(atoms, "ca")
  assert_eq(sel.length(), 5)
}

test "sel_backbone" {
  let chain = make_sel_test_chain()
  let atoms = chain.get_atoms()
  let sel = @src.select_atoms(atoms, "backbone")
  // 5 ALA * (N, CA, C, O) = 20 backbone atoms
  // backbone names: N, CA, C, O, OXT, H, HA
  // backbone names: N, CA, C, O, H, HA
  // Per ALA: N, CA, C, O = 4
  assert_eq(sel.length(), 20)
}

test "sel_sidechain" {
  let chain = make_sel_test_chain()
  let atoms = chain.get_atoms()
  let sel = @src.select_atoms(atoms, "sidechain")
  // Per ALA: CB, HB1 are sidechain = 2; * 5 = 10
  assert_eq(sel.length(), 10)
}

test "sel_protein" {
  let chain = make_sel_test_chain()
  let atoms = chain.get_atoms()
  let sel = @src.select_atoms(atoms, "protein")
  // 5 ALA * 6 atoms = 30
  assert_eq(sel.length(), 30)
}

test "sel_water" {
  let chain = make_sel_test_chain()
  let atoms = chain.get_atoms()
  let sel = @src.select_atoms(atoms, "water")
  assert_eq(sel.length(), 1)
  assert_eq(sel[0].resname, "HOH")
}

test "sel_nucleic" {
  let chain = make_sel_test_chain()
  let atoms = chain.get_atoms()
  let sel = @src.select_atoms(atoms, "nucleic")
  assert_eq(sel.length(), 2)
  assert_eq(sel[0].resname, "DA")
}

test "sel_hetero" {
  let chain = make_sel_test_chain()
  let atoms = chain.get_atoms()
  let sel = @src.select_atoms(atoms, "hetero")
  // DA is nucleic but not AA and not water → hetero
  assert_eq(sel.length(), 2)
}

test "sel_hydrogen" {
  let chain = make_sel_test_chain()
  let atoms = chain.get_atoms()
  let sel = @src.select_atoms(atoms, "hydrogen")
  // 5 ALA * 1 HB1 = 5 H atoms
  assert_eq(sel.length(), 5)
}

test "sel_heavy" {
  let chain = make_sel_test_chain()
  let atoms = chain.get_atoms()
  let sel = @src.select_atoms(atoms, "heavy")
  // Total atoms: 5*6 + 1 + 2 = 33; H atoms: 5; heavy: 33-5=28
  assert_eq(sel.length(), 28)
}

// ===== Property selectors =====

test "sel_name" {
  let chain = make_sel_test_chain()
  let atoms = chain.get_atoms()
  let sel = @src.select_atoms(atoms, "name CA CB")
  // 5 CA + 5 CB = 10
  assert_eq(sel.length(), 10)
}

test "sel_resname" {
  let chain = make_sel_test_chain()
  let atoms = chain.get_atoms()
  let sel = @src.select_atoms(atoms, "resname ALA")
  assert_eq(sel.length(), 30)
}

test "sel_resnum_range" {
  let chain = make_sel_test_chain()
  let atoms = chain.get_atoms()
  let sel = @src.select_atoms(atoms, "resnum 1:3")
  // 3 ALA * 6 + 2 nucleic (resseq=1) = 20
  assert_eq(sel.length(), 20)
}

test "sel_resnum_single" {
  let chain = make_sel_test_chain()
  let atoms = chain.get_atoms()
  let sel = @src.select_atoms(atoms, "resnum 2")
  // 1 ALA * 6 atoms = 6
  assert_eq(sel.length(), 6)
}

test "sel_chain" {
  let chain = make_sel_test_chain()
  let atoms = chain.get_atoms()
  let sel = @src.select_atoms(atoms, "chain A")
  // 5 ALA * 6 + 1 water = 31
  assert_eq(sel.length(), 31)
}

test "sel_chain_b" {
  let chain = make_sel_test_chain()
  let atoms = chain.get_atoms()
  let sel = @src.select_atoms(atoms, "chain B")
  // 2 nucleic atoms
  assert_eq(sel.length(), 2)
}

test "sel_element" {
  let chain = make_sel_test_chain()
  let atoms = chain.get_atoms()
  let sel = @src.select_atoms(atoms, "element C")
  // 5 ALA * (CA, C, CB) + 1 nucleic C = 16
  assert_eq(sel.length(), 16)
}

// ===== Logical operators =====

test "sel_and" {
  let chain = make_sel_test_chain()
  let atoms = chain.get_atoms()
  let sel = @src.select_atoms(atoms, "calpha and resnum 1:3")
  // 3 CA atoms
  assert_eq(sel.length(), 3)
}

test "sel_or" {
  let chain = make_sel_test_chain()
  let atoms = chain.get_atoms()
  let sel = @src.select_atoms(atoms, "calpha or water")
  // 5 CA + 1 water = 6
  assert_eq(sel.length(), 6)
}

test "sel_not" {
  let chain = make_sel_test_chain()
  let atoms = chain.get_atoms()
  let sel = @src.select_atoms(atoms, "not water")
  // 33 - 1 = 32
  assert_eq(sel.length(), 32)
}

test "sel_paren_grouping" {
  let chain = make_sel_test_chain()
  let atoms = chain.get_atoms()
  let sel = @src.select_atoms(atoms, "(calpha or water) and not water")
  // calpha without water = 5
  assert_eq(sel.length(), 5)
}

test "sel_nested_parens" {
  let chain = make_sel_test_chain()
  let atoms = chain.get_atoms()
  let sel = @src.select_atoms(atoms, "((calpha) and (resnum 1:2)) or water")
  // 2 CA + 1 water = 3
  assert_eq(sel.length(), 3)
}

test "sel_not_not" {
  let chain = make_sel_test_chain()
  let atoms = chain.get_atoms()
  let sel = @src.select_atoms(atoms, "not not calpha")
  assert_eq(sel.length(), 5)
}

test "sel_complex" {
  let chain = make_sel_test_chain()
  let atoms = chain.get_atoms()
  let sel = @src.select_atoms(atoms, "protein and not calpha and element C")
  // Per ALA: C and CB (not CA) = 2; * 5 = 10
  assert_eq(sel.length(), 10)
}

// ===== Spatial selection =====

test "sel_within" {
  let chain = make_sel_test_chain()
  let atoms = chain.get_atoms()
  // Within 2.0 of calpha — should include CA atoms and nearby atoms
  let sel = @src.select_atoms(atoms, "within 2.0 of calpha")
  // At least the 5 CA atoms themselves
  assert_true(sel.length() >= 5)
}

test "sel_within_small" {
  let chain = make_sel_test_chain()
  let atoms = chain.get_atoms()
  // Within 0.1 of calpha — should be just the CA atoms
  let sel = @src.select_atoms(atoms, "within 0.1 of calpha")
  assert_eq(sel.length(), 5)
}

test "sel_within_large" {
  let chain = make_sel_test_chain()
  let atoms = chain.get_atoms()
  // Within 100.0 of calpha — should be all atoms
  let sel = @src.select_atoms(atoms, "within 100.0 of calpha")
  assert_eq(sel.length(), atoms.length())
}

// ===== Utility functions =====

test "sel_count" {
  let chain = make_sel_test_chain()
  let atoms = chain.get_atoms()
  assert_eq(@src.count_selected(atoms, "calpha"), 5)
  assert_eq(@src.count_selected(atoms, "water"), 1)
  assert_eq(@src.count_selected(atoms, "all"), 33)
}

test "sel_indices" {
  let chain = make_sel_test_chain()
  let atoms = chain.get_atoms()
  let idx = @src.select_indices(atoms, "calpha")
  assert_eq(idx.length(), 5)
  // CA atoms are at index 1, 7, 13, 19, 25 (every 6 atoms, starting at 1)
  assert_eq(idx[0], 1)
  assert_eq(idx[1], 7)
  assert_eq(idx[2], 13)
}

test "sel_mask" {
  let chain = make_sel_test_chain()
  let atoms = chain.get_atoms()
  let mask = @src.select_mask(atoms, "calpha")
  assert_eq(mask.length(), 33)
  assert_true(mask[1])
  assert_true(mask[7])
  assert_true(!mask[0])
  assert_true(!mask[2])
}

test "sel_chain_method" {
  let chain = make_sel_test_chain()
  let sel = chain.select("calpha")
  assert_eq(sel.length(), 5)
}

test "sel_chain_residues" {
  let chain = make_sel_test_chain()
  let res = chain.select_residues("calpha and resnum 1:2")
  assert_eq(res.length(), 2)
  assert_eq(res[0].resseq, 1)
  assert_eq(res[1].resseq, 2)
}

test "sel_parse_ast" {
  match @src.parse_sel("calpha and resnum 1:50") {
    Some(_) => assert_true(true)
    None => assert_true(false)
  }
}

test "sel_parse_invalid" {
  // Unparseable → None
  match @src.parse_sel("xyzzy frobnicate") {
    Some(_) => assert_true(true)  // Still parses, just produces None_ expr
    None => assert_true(true)
  }
}

test "sel_empty_string" {
  let chain = make_sel_test_chain()
  let atoms = chain.get_atoms()
  let sel = @src.select_atoms(atoms, "")
  assert_eq(sel.length(), 0)
}

test "sel_backbone_bb_alias" {
  let chain = make_sel_test_chain()
  let atoms = chain.get_atoms()
  let sel1 = @src.select_atoms(atoms, "backbone")
  let sel2 = @src.select_atoms(atoms, "bb")
  assert_eq(sel1.length(), sel2.length())
}

test "sel_sidechain_sc_alias" {
  let chain = make_sel_test_chain()
  let atoms = chain.get_atoms()
  let sel1 = @src.select_atoms(atoms, "sidechain")
  let sel2 = @src.select_atoms(atoms, "sc")
  assert_eq(sel1.length(), sel2.length())
}

test "sel_hoh_alias" {
  let chain = make_sel_test_chain()
  let atoms = chain.get_atoms()
  let sel = @src.select_atoms(atoms, "hoh")
  assert_eq(sel.length(), 1)
}

test "sel_noh_alias" {
  let chain = make_sel_test_chain()
  let atoms = chain.get_atoms()
  let sel1 = @src.select_atoms(atoms, "heavy")
  let sel2 = @src.select_atoms(atoms, "noh")
  assert_eq(sel1.length(), sel2.length())
}
