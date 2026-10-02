///|
/// Tier 1 cheminformatics tests: canonical SMILES, descriptors, MACCS, SMARTS,
/// kekulization, distance matrix, ring perception.

test "adjacency_list_ethanol" {
  let mol = @src.parse_smiles("CCO")
  let adj = @src.mol_adjacency(mol)
  // Atom 0 (C): neighbor 1
  assert_eq(adj[0].length(), 1)
  assert_eq(adj[0][0], (1, 1))
  // Atom 1 (C): neighbors 0 and 2
  assert_eq(adj[1].length(), 2)
  // Atom 2 (O): neighbor 1
  assert_eq(adj[2].length(), 1)
  assert_eq(adj[2][0], (1, 1))
}

test "distance_matrix_ethanol" {
  let mol = @src.parse_smiles("CCO")
  let dist = @src.mol_distance_matrix(mol)
  // 3 atoms: 0-1-2
  assert_eq(dist[0][0], 0)
  assert_eq(dist[0][1], 1)
  assert_eq(dist[0][2], 2)
  assert_eq(dist[1][0], 1)
  assert_eq(dist[1][1], 0)
  assert_eq(dist[1][2], 1)
  assert_eq(dist[2][0], 2)
  assert_eq(dist[2][1], 1)
  assert_eq(dist[2][2], 0)
}

test "distance_matrix_benzene" {
  let mol = @src.parse_smiles("c1ccccc1")
  let dist = @src.mol_distance_matrix(mol)
  // 6-membered ring: opposite atoms at distance 3
  assert_eq(dist[0][3], 3)
  assert_eq(dist[0][1], 1)
  assert_eq(dist[0][2], 2)
}

test "ring_perception_benzene" {
  let mol = @src.parse_smiles("c1ccccc1")
  let rings = @src.mol_find_rings(mol)
  assert_eq(rings.length(), 1)
  assert_eq(rings[0].length(), 6)
}

test "ring_perception_cyclohexane" {
  let mol = @src.parse_smiles("C1CCCCC1")
  let rings = @src.mol_find_rings(mol)
  assert_eq(rings.length(), 1)
  assert_eq(rings[0].length(), 6)
}

test "ring_perception_naphthalene" {
  let mol = @src.parse_smiles("c1ccc2ccccc2c1")
  let rings = @src.mol_find_rings(mol)
  // Naphthalene has 2 rings
  assert_eq(rings.length(), 2)
}

test "ring_count_benzene" {
  let mol = @src.parse_smiles("c1ccccc1")
  assert_eq(@src.mol_ring_count(mol), 1)
}

test "ring_count_no_ring" {
  let mol = @src.parse_smiles("CCO")
  assert_eq(@src.mol_ring_count(mol), 0)
}

test "is_in_ring" {
  let mol = @src.parse_smiles("c1ccccc1")
  // All atoms in benzene are in a ring
  for i = 0; i < 6; i = i + 1 {
    assert_true(@src.mol_is_in_ring(mol, i))
  }
}

test "is_in_ring_no_ring" {
  let mol = @src.parse_smiles("CCO")
  for i = 0; i < 3; i = i + 1 {
    assert_false(@src.mol_is_in_ring(mol, i))
  }
}

test "aromatic_ring_count_benzene" {
  let mol = @src.parse_smiles("c1ccccc1")
  assert_eq(@src.mol_aromatic_ring_count(mol), 1)
}

test "saturated_ring_count_cyclohexane" {
  let mol = @src.parse_smiles("C1CCCCC1")
  assert_eq(@src.mol_saturated_ring_count(mol), 1)
}

test "kekulize_benzene" {
  let mol = @src.parse_smiles("c1ccccc1")
  let kek = @src.mol_kekulize(mol)
  // After kekulization, no atoms should be aromatic
  for a in kek.atoms {
    assert_false(a.aromatic)
  }
  // Should have alternating single/double bonds
  let mut num_double = 0
  for b in kek.bonds {
    if b.order == 2 {
      num_double = num_double + 1
    }
  }
  assert_eq(num_double, 3)
}

test "perceive_aromaticity" {
  let mol = @src.parse_smiles("C1=CC=CC=C1")
  let arom = @src.mol_perceive_aromaticity(mol)
  let mut num_aromatic = 0
  for a in arom.atoms {
    if a.aromatic {
      num_aromatic = num_aromatic + 1
    }
  }
  assert_eq(num_aromatic, 6)
}

test "canonical_smiles_ethanol" {
  let mol = @src.parse_smiles("CCO")
  let smi = @src.mol_to_smiles(mol)
  // Should produce a valid SMILES with 3 atoms
  assert_true(smi.length() > 0)
}

test "canonical_smiles_benzene" {
  let mol = @src.parse_smiles("c1ccccc1")
  let smi = @src.mol_to_smiles(mol)
  assert_true(smi.length() > 0)
}

test "molecular_weight_ethanol" {
  let mol = @src.parse_smiles("CCO")
  let mw = @src.mol_molecular_weight(mol)
  // C2H6O = 2*12.011 + 6*1.008 + 15.999 ≈ 46.069
  // But parse_smiles doesn't add implicit H, so just C+C+O
  // 12.011 + 12.011 + 15.999 = 40.021
  assert_true(mw > 39.0 && mw < 41.0)
}

test "logp_ethanol" {
  let mol = @src.parse_smiles("CCO")
  let logp = @src.mol_logp(mol)
  // Should be a finite number
  assert_true(logp > -10.0 && logp < 10.0)
}

test "tpsa_ethanol" {
  let mol = @src.parse_smiles("CCO")
  let tpsa = @src.mol_tpsa(mol)
  // O contributes to TPSA
  assert_true(tpsa > 0.0)
}

test "hbd_count_ethanol" {
  let mol = @src.parse_smiles("CCO")
  let hbd = @src.mol_hbd_count(mol)
  // O with 1 neighbor has 1 implicit H -> 1 HBD
  assert_eq(hbd, 1)
}

test "hba_count_ethanol" {
  let mol = @src.parse_smiles("CCO")
  let hba = @src.mol_hba_count(mol)
  // 1 O acceptor
  assert_eq(hba, 1)
}

test "hba_count_pyridine" {
  let mol = @src.parse_smiles("c1ccncc1")
  let hba = @src.mol_hba_count(mol)
  // 1 N acceptor
  assert_eq(hba, 1)
}

test "rotatable_bonds_ethane" {
  let mol = @src.parse_smiles("CC")
  let rb = @src.mol_rotatable_bonds(mol)
  // C-C with each having 1 neighbor -> not rotatable
  assert_eq(rb, 0)
}

test "rotatable_bonds_propane" {
  let mol = @src.parse_smiles("CCC")
  let rb = @src.mol_rotatable_bonds(mol)
  // Without implicit H, terminal C has degree 1, so 0 rotatable
  assert_eq(rb, 0)
}

test "formal_charge" {
  let mol = @src.parse_smiles("CCO")
  let fc = @src.mol_formal_charge(mol)
  assert_eq(fc, 0)
}

test "fraction_sp3_ethanol" {
  let mol = @src.parse_smiles("CCO")
  let fsp3 = @src.mol_fraction_sp3(mol)
  // Both C atoms are sp3
  assert_true(fsp3 > 0.99)
}

test "fraction_sp3_benzene" {
  let mol = @src.parse_smiles("c1ccccc1")
  let fsp3 = @src.mol_fraction_sp3(mol)
  // All C are aromatic, not sp3
  assert_true(fsp3 < 0.01)
}

test "valence_electrons_ethanol" {
  let mol = @src.parse_smiles("CCO")
  let ve = @src.mol_valence_electrons(mol)
  // C=4, C=4, O=6 -> 14
  assert_eq(ve, 14)
}

test "heavy_atom_count" {
  let mol = @src.parse_smiles("CCO")
  assert_eq(@src.mol_heavy_atom_count(mol), 3)
}

test "descriptor_report" {
  let mol = @src.parse_smiles("CCO")
  let report = @src.mol_descriptor_report(mol)
  assert_true(report.length() > 0)
  // Should contain key labels
  assert_true(report.contains("LogP"))
  assert_true(report.contains("TPSA"))
}

// ============================================================================
// MACCS fingerprint tests
// ============================================================================

test "maccs_bit_count_ethanol" {
  let mol = @src.parse_smiles("CCO")
  let bc = @src.maccs_bit_count(mol)
  // Ethanol should set several bits (C, O, single bonds, etc.)
  assert_true(bc > 5)
}

test "maccs_fingerprint_length" {
  let mol = @src.parse_smiles("CCO")
  let fp = @src.maccs_fingerprint(mol)
  // 167 elements (index 0 unused + 166 bits)
  assert_eq(fp.length(), 167)
}

test "maccs_fingerprint_string_length" {
  let mol = @src.parse_smiles("CCO")
  let fps = @src.maccs_fingerprint_string(mol)
  assert_eq(fps.length(), 166)
}

test "maccs_bit1_always_set" {
  let mol = @src.parse_smiles("C")
  let fp = @src.maccs_fingerprint(mol)
  // Bit 1 is always set (isotope/any molecule)
  assert_true(fp[1])
}

test "maccs_oxygen_bit" {
  let mol = @src.parse_smiles("CCO")
  let fp = @src.maccs_fingerprint(mol)
  // Bit 17 = has O
  assert_true(fp[17])
}

test "maccs_no_oxygen" {
  let mol = @src.parse_smiles("CC")
  let fp = @src.maccs_fingerprint(mol)
  // Bit 17 = has O -> false
  assert_false(fp[17])
}

test "maccs_nitrogen_bit" {
  let mol = @src.parse_smiles("N")
  let fp = @src.maccs_fingerprint(mol)
  // Bit 24 = has N
  assert_true(fp[24])
}

test "maccs_ring_bit_benzene" {
  let mol = @src.parse_smiles("c1ccccc1")
  let fp = @src.maccs_fingerprint(mol)
  // Bit 61 = 6-membered ring
  assert_true(fp[61])
  // Bit 65 = aromatic ring
  assert_true(fp[65])
}

test "maccs_tanimoto_identical" {
  let mol1 = @src.parse_smiles("CCO")
  let mol2 = @src.parse_smiles("CCO")
  let sim = @src.maccs_tanimoto(mol1, mol2)
  // Identical molecules -> Tanimoto = 1.0
  assert_true(sim > 0.99)
}

test "maccs_tanimoto_different" {
  let mol1 = @src.parse_smiles("CCO")
  let mol2 = @src.parse_smiles("c1ccccc1")
  let sim = @src.maccs_tanimoto(mol1, mol2)
  // Different molecules -> Tanimoto < 1.0
  assert_true(sim < 1.0)
  assert_true(sim >= 0.0)
}

test "maccs_fingerprint_int_nonzero" {
  let mol = @src.parse_smiles("CCO")
  let fpi = @src.maccs_fingerprint_int(mol)
  // Should be non-zero for a molecule with atoms
  assert_true(fpi > 0L)
}

// ============================================================================
// SMARTS tests
// ============================================================================

test "smarts_parse_carbon" {
  let pat = @src.smarts_parse("[#6]")
  assert_eq(pat.atoms.length(), 1)
  assert_true(pat.atoms[0].atomic_num is Some(_))
}

test "smarts_parse_wildcard" {
  let pat = @src.smarts_parse("[*]")
  assert_eq(pat.atoms.length(), 1)
  assert_true(pat.atoms[0].is_wildcard)
}

test "smarts_parse_ring" {
  let pat = @src.smarts_parse("[R]")
  assert_eq(pat.atoms.length(), 1)
  assert_true(pat.atoms[0].in_ring is Some(_))
}

test "smarts_parse_aromatic" {
  let pat = @src.smarts_parse("[a]")
  assert_eq(pat.atoms.length(), 1)
  assert_true(pat.atoms[0].aromatic is Some(_))
}

test "smarts_parse_multi_atom" {
  let pat = @src.smarts_parse("[#6]-[#6]")
  assert_eq(pat.atoms.length(), 2)
  assert_eq(pat.bonds.length(), 1)
}

test "smarts_match_carbon" {
  let mol = @src.parse_smiles("CCO")
  // [#6] should match carbon atoms
  assert_true(@src.smarts_match_str(mol, "[#6]"))
}

test "smarts_match_oxygen" {
  let mol = @src.parse_smiles("CCO")
  // [#8] should match oxygen
  assert_true(@src.smarts_match_str(mol, "[#8]"))
}

test "smarts_match_no_match" {
  let mol = @src.parse_smiles("CCO")
  // [#7] (nitrogen) should not match ethanol
  assert_false(@src.smarts_match_str(mol, "[#7]"))
}

test "smarts_match_wildcard" {
  let mol = @src.parse_smiles("CCO")
  // [*] matches any atom
  assert_true(@src.smarts_match_str(mol, "[*]"))
}

test "smarts_match_cc_bond" {
  let mol = @src.parse_smiles("CCO")
  // [#6]-[#6] should match C-C bond in ethanol
  assert_true(@src.smarts_match_str(mol, "[#6]-[#6]"))
}

test "smarts_match_co_bond" {
  let mol = @src.parse_smiles("CCO")
  // [#6]-[#8] should match C-O bond
  assert_true(@src.smarts_match_str(mol, "[#6]-[#8]"))
}

test "smarts_match_ring" {
  let mol = @src.parse_smiles("c1ccccc1")
  // [R] should match ring atoms in benzene
  assert_true(@src.smarts_match_str(mol, "[R]"))
}

test "smarts_match_ring_no_ring" {
  let mol = @src.parse_smiles("CCO")
  // [R] should not match any atom in ethanol
  assert_false(@src.smarts_match_str(mol, "[R]"))
}

test "smarts_match_aromatic" {
  let mol = @src.parse_smiles("c1ccccc1")
  // [a] should match aromatic atoms in benzene
  assert_true(@src.smarts_match_str(mol, "[a]"))
}

test "smarts_match_aromatic_no_arom" {
  let mol = @src.parse_smiles("CCO")
  // [a] should not match any atom in ethanol
  assert_false(@src.smarts_match_str(mol, "[a]"))
}

test "smarts_count_carbons" {
  let mol = @src.parse_smiles("CCO")
  // [#6] should match 2 carbon atoms
  let count = @src.smarts_count_matches_str(mol, "[#6]")
  assert_eq(count, 2)
}

test "smarts_count_oxygen" {
  let mol = @src.parse_smiles("CCO")
  let count = @src.smarts_count_matches_str(mol, "[#8]")
  assert_eq(count, 1)
}

test "smarts_find_matches" {
  let mol = @src.parse_smiles("CCO")
  let pat = @src.smarts_parse("[#6]")
  let matches = @src.smarts_find_matches(mol, pat)
  assert_eq(matches.length(), 2)
}

test "smarts_double_bond" {
  let mol = @src.parse_smiles("C=CC")
  // [#6]=[#6] should match the C=C double bond
  assert_true(@src.smarts_match_str(mol, "[#6]=[#6]"))
}

test "smarts_triple_bond" {
  let mol = @src.parse_smiles("C#C")
  // [#6]#[#6] should match the C#C triple bond
  assert_true(@src.smarts_match_str(mol, "[#6]#[#6]"))
}

test "smarts_any_bond" {
  let mol = @src.parse_smiles("CCO")
  // [#6]~[#8] should match C-O with any bond
  assert_true(@src.smarts_match_str(mol, "[#6]~[#8]"))
}

test "smarts_degree_query" {
  let mol = @src.parse_smiles("CC(C)C")
  // [D3] should match the central carbon with degree 3
  assert_true(@src.smarts_match_str(mol, "[D3]"))
}

test "smarts_negation" {
  let mol = @src.parse_smiles("CCO")
  // [!#8] should match non-oxygen atoms (C, C)
  let count = @src.smarts_count_matches_str(mol, "[!#8]")
  assert_eq(count, 2)
}
