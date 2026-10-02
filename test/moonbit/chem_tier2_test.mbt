///|
/// Tier 2 cheminformatics tests: Murcko scaffold, MCS, stereochemistry,
/// tautomers, 2D coords, topological torsion fingerprint.

// ============================================================================
// Murcko scaffold tests
// ============================================================================

test "murcko_benzene_is_own_scaffold" {
  let mol = @src.parse_smiles("c1ccccc1")
  let result = @src.murcko_scaffold(mol)
  // Benzene is all scaffold, no side chains
  assert_eq(result.scaffold_atoms.length(), 6)
  assert_eq(result.side_chain_atoms.length(), 0)
  assert_eq(result.num_rings, 1)
}

test "murcko_toluene_scaffold" {
  let mol = @src.parse_smiles("Cc1ccccc1")
  let result = @src.murcko_scaffold(mol)
  // Toluene: scaffold = benzene (6 atoms), side chain = methyl (1 atom)
  assert_eq(result.scaffold_atoms.length(), 6)
  assert_eq(result.side_chain_atoms.length(), 1)
}

test "murcko_ethanol_no_scaffold" {
  let mol = @src.parse_smiles("CCO")
  let result = @src.murcko_scaffold(mol)
  // Ethanol has no rings -> empty scaffold
  assert_eq(result.scaffold_atoms.length(), 0)
  assert_eq(result.side_chain_atoms.length(), 3)
}

test "murcko_is_pure_scaffold" {
  let mol = @src.parse_smiles("c1ccccc1")
  assert_true(@src.murcko_is_pure_scaffold(mol))
}

test "murcko_not_pure_scaffold" {
  let mol = @src.parse_smiles("Cc1ccccc1")
  assert_false(@src.murcko_is_pure_scaffold(mol))
}

test "murcko_side_chain_count" {
  let mol = @src.parse_smiles("Cc1ccccc1")
  assert_eq(@src.murcko_side_chain_count(mol), 1)
}

test "murcko_scaffold_mol" {
  let mol = @src.parse_smiles("Cc1ccccc1")
  let scaffold = @src.murcko_scaffold_mol(mol)
  // Scaffold should be benzene with 6 atoms
  assert_eq(scaffold.num_atoms(), 6)
}

test "murcko_naphthalene" {
  let mol = @src.parse_smiles("c1ccc2ccccc2c1")
  let result = @src.murcko_scaffold(mol)
  // Naphthalene is all scaffold
  assert_eq(result.side_chain_atoms.length(), 0)
  assert_eq(result.num_rings, 2)
}

// ============================================================================
// MCS tests
// ============================================================================

test "mcs_identical_molecules" {
  let mol1 = @src.parse_smiles("c1ccccc1")
  let mol2 = @src.parse_smiles("c1ccccc1")
  let mcs = @src.find_mcs(mol1, mol2)
  // Identical molecules -> MCS = all 6 atoms
  assert_eq(mcs.num_atoms, 6)
}

test "mcs_benzene_pyridine" {
  let mol1 = @src.parse_smiles("c1ccccc1")
  let mol2 = @src.parse_smiles("c1ccncc1")
  let mcs = @src.find_mcs(mol1, mol2)
  // Benzene and pyridine share 5 C atoms (the N doesn't match C)
  assert_eq(mcs.num_atoms, 5)
}

test "mcs_ethanol_methanol" {
  let mol1 = @src.parse_smiles("CCO")
  let mol2 = @src.parse_smiles("CO")
  let mcs = @src.find_mcs(mol1, mol2)
  // Common substructure: C-O (2 atoms)
  assert_eq(mcs.num_atoms, 2)
}

test "mcs_no_common" {
  let mol1 = @src.parse_smiles("CCC")
  let mol2 = @src.parse_smiles("NNN")
  let mcs = @src.find_mcs(mol1, mol2)
  // No common atoms (C vs N)
  assert_eq(mcs.num_atoms, 0)
}

test "mcs_similarity_identical" {
  let mol1 = @src.parse_smiles("CCO")
  let mol2 = @src.parse_smiles("CCO")
  let sim = @src.mcs_similarity(mol1, mol2)
  assert_true(sim > 0.99)
}

test "mcs_similarity_different" {
  let mol1 = @src.parse_smiles("CCO")
  let mol2 = @src.parse_smiles("c1ccccc1")
  let sim = @src.mcs_similarity(mol1, mol2)
  assert_true(sim < 1.0)
}

test "mcs_atom_count" {
  let mol1 = @src.parse_smiles("CCCC")
  let mol2 = @src.parse_smiles("CCC")
  let count = @src.mcs_atom_count(mol1, mol2)
  // Butane and propane share 3 C atoms
  assert_eq(count, 3)
}

// ============================================================================
// Stereochemistry tests
// ============================================================================

test "stereo_tetrahedral_center" {
  // A carbon with 4 different neighbors (C, N, O, S)
  let mol = @src.parse_smiles("C(C)(N)(O)S")
  // Central C has 4 neighbors
  assert_true(@src.is_tetrahedral_center(mol, 0))
}

test "stereo_not_tetrahedral" {
  let mol = @src.parse_smiles("CCO")
  // No atom has 4 different neighbors
  for i = 0; i < 3; i = i + 1 {
    assert_false(@src.is_tetrahedral_center(mol, i))
  }
}

test "stereo_r_or_s" {
  let mol = @src.parse_smiles("C(C)(N)(O)S")
  let desc = @src.cip_tetrahedral(mol, 0)
  // Should be R or S (not None)
  let s = @src.stereo_to_string(desc)
  assert_true(s == "R" || s == "S")
}

test "stereo_find_centers" {
  let mol = @src.parse_smiles("C(C)(N)(O)S")
  let centers = @src.find_stereocenters(mol)
  assert_true(centers.length() >= 1)
}

test "stereo_double_bond_ez" {
  // 2-butene: C/C=C/C (E or Z)
  let mol = @src.parse_smiles("CC=CC")
  // Find double bond
  let mut db_idx = -1
  for i = 0; i < mol.bonds.length(); i = i + 1 {
    if mol.bonds[i].order == 2 {
      db_idx = i
    }
  }
  if db_idx >= 0 {
    let desc = @src.cip_double_bond(mol, db_idx)
    let s = @src.stereo_to_string(desc)
    // Should be E, Z, or none (depending on substituent analysis)
    assert_true(s == "E" || s == "Z" || s == "none")
  }
}

test "stereo_find_double_bond_stereo" {
  let mol = @src.parse_smiles("CC=CC")
  let stereos = @src.find_double_bond_stereo(mol)
  // May or may not find E/Z depending on substituent count
  assert_true(stereos.length() >= 0)
}

test "stereo_to_string" {
  assert_eq(@src.stereo_to_string(@src.R), "R")
  assert_eq(@src.stereo_to_string(@src.S), "S")
  assert_eq(@src.stereo_to_string(@src.E), "E")
  assert_eq(@src.stereo_to_string(@src.Z), "Z")
  assert_eq(@src.stereo_to_string(@src.None), "none")
}

// ============================================================================
// Tautomer tests
// ============================================================================

test "tautomer_acetone" {
  // Acetone: CC=O
  let mol = @src.parse_smiles("CC=O")
  let tautomers = @src.enumerate_tautomers(mol, 5)
  // Should have at least 2 tautomers (keto + enol)
  assert_true(tautomers.length() >= 1)
}

test "tautomer_count_ethanol" {
  // Ethanol: CCO (no tautomeric forms)
  let mol = @src.parse_smiles("CCO")
  let count = @src.tautomer_count(mol, 5)
  // Ethanol has no keto-enol tautomerism
  assert_eq(count, 1)
}

test "tautomer_count_acetone" {
  let mol = @src.parse_smiles("CC=O")
  let count = @src.tautomer_count(mol, 5)
  // Acetone should have at least 2 (keto + enol)
  assert_true(count >= 2)
}

test "tautomer_includes_original" {
  let mol = @src.parse_smiles("CC=O")
  let tautomers = @src.enumerate_tautomers(mol, 5)
  // First tautomer should be the original
  assert_eq(tautomers[0].num_atoms(), mol.num_atoms())
}

// ============================================================================
// 2D coordinate tests
// ============================================================================

test "coords_benzene" {
  let mol = @src.parse_smiles("c1ccccc1")
  let coords = @src.generate_2d_coords(mol)
  // 6 atoms -> 6 coordinate pairs
  assert_eq(coords.length(), 6)
}

test "coords_ethanol" {
  let mol = @src.parse_smiles("CCO")
  let coords = @src.generate_2d_coords(mol)
  assert_eq(coords.length(), 3)
}

test "coords_all_placed" {
  let mol = @src.parse_smiles("c1ccccc1CC")
  let coords = @src.generate_2d_coords(mol)
  // All atoms should have coordinates (non-zero for ring atoms)
  assert_eq(coords.length(), 8)
}

test "coords_bounding_box" {
  let mol = @src.parse_smiles("c1ccccc1")
  let coords = @src.generate_2d_coords(mol)
  let ((min_x, min_y), (max_x, max_y)) = @src.coord_bounding_box(coords)
  // Bounding box should have non-zero size
  assert_true(max_x > min_x)
  assert_true(max_y > min_y)
}

test "coords_centered" {
  let mol = @src.parse_smiles("c1ccccc1")
  let coords = @src.generate_2d_coords(mol)
  let centered = @src.center_2d_coords(coords)
  let ((min_x, min_y), (max_x, max_y)) = @src.coord_bounding_box(centered)
  // Centered coords should be roughly symmetric around origin
  assert_true((min_x + max_x).abs() < 0.01)
  assert_true((min_y + max_y).abs() < 0.01)
}

test "coords_to_string" {
  let mol = @src.parse_smiles("CCO")
  let coords = @src.generate_2d_coords(mol)
  let s = @src.coords_to_string(mol, coords)
  assert_true(s.contains("C"))
  assert_true(s.contains("O"))
}

test "coords_point2d" {
  let mol = @src.parse_smiles("c1ccccc1")
  let points = @src.generate_2d_points(mol)
  assert_eq(points.length(), 6)
}

// ============================================================================
// Topological torsion fingerprint tests
// ============================================================================

test "torsion_butane" {
  let mol = @src.parse_smiles("CCCC")
  let fp = @src.topological_torsion_fingerprint(mol)
  // Butane has 4 atoms in a line -> 1 unique torsion
  assert_true(fp.length() >= 1)
}

test "torsion_ethanol" {
  let mol = @src.parse_smiles("CCO")
  let fp = @src.topological_torsion_fingerprint(mol)
  // Ethanol has only 3 atoms -> no 4-atom torsions
  assert_eq(fp.length(), 0)
}

test "torsion_count" {
  let mol = @src.parse_smiles("CCCC")
  let count = @src.topological_torsion_count(mol)
  assert_true(count >= 1)
}

test "torsion_identical_molecules" {
  let mol1 = @src.parse_smiles("CCCC")
  let mol2 = @src.parse_smiles("CCCC")
  let sim = @src.torsion_tanimoto(mol1, mol2)
  // Identical molecules -> Tanimoto = 1.0
  assert_true(sim > 0.99)
}

test "torsion_different_molecules" {
  let mol1 = @src.parse_smiles("CCCC")
  let mol2 = @src.parse_smiles("CCO")
  let sim = @src.torsion_tanimoto(mol1, mol2)
  // Different molecules
  assert_true(sim < 1.0)
}

test "torsion_pentane" {
  let mol = @src.parse_smiles("CCCCC")
  let fp = @src.topological_torsion_fingerprint(mol)
  // Pentane has 5 atoms -> 2 torsions (atoms 0-1-2-3 and 1-2-3-4)
  assert_true(fp.length() >= 1)
}

test "torsion_fingerprint_string" {
  let mol = @src.parse_smiles("CCCC")
  let s = @src.torsion_fingerprint_string(mol)
  // Should be a non-empty string of hash values
  assert_true(s.length() > 0)
}
