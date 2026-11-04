///|
/// Tests for the ChemmineR module (Bioconductor ChemmineR equivalent).
///
/// Covers SDF atom/bond/molecule structures, molecular descriptors
/// (MW, formula, ring count, heavy atoms), atom-pair fingerprints,
/// Tanimoto/Dice similarity, substructure search, SDF parsing and
/// serialization, and sample data.

// ===========================================================================
// 1. SdfAtom
// ===========================================================================

///|
test "chem_atom_creation_and_accessors" {
  let a = @src.SdfAtom::new("C", 1.0, 2.0, 3.0, 0)
  assert_eq(a.symbol(), "C")
  assert_true(a.x() == 1.0)
  assert_true(a.y() == 2.0)
  assert_true(a.z() == 3.0)
  assert_eq(a.charge(), 0)
}

///|
test "chem_atom_distance" {
  let a = @src.SdfAtom::new("C", 0.0, 0.0, 0.0, 0)
  let b = @src.SdfAtom::new("C", 3.0, 4.0, 0.0, 0)
  let d = a.distance(b)
  assert_true(d > 4.99)
  assert_true(d < 5.01)
}

///|
test "chem_atom_distance_zero" {
  let a = @src.SdfAtom::new("O", 1.0, 1.0, 1.0, 0)
  assert_true(a.distance(a) == 0.0)
}

// ===========================================================================
// 2. SdfBond
// ===========================================================================

///|
test "chem_bond_creation_and_accessors" {
  let b = @src.SdfBond::new(0, 1, 2)
  assert_eq(b.atom1(), 0)
  assert_eq(b.atom2(), 1)
  assert_eq(b.bond_type(), 2)
  assert_eq(b.bond_label(), "double")
}

///|
test "chem_bond_labels" {
  assert_eq(@src.SdfBond::new(0, 1, 1).bond_label(), "single")
  assert_eq(@src.SdfBond::new(0, 1, 2).bond_label(), "double")
  assert_eq(@src.SdfBond::new(0, 1, 3).bond_label(), "triple")
  assert_eq(@src.SdfBond::new(0, 1, 4).bond_label(), "aromatic")
  assert_eq(@src.SdfBond::new(0, 1, 9).bond_label(), "unknown")
}

// ===========================================================================
// 3. SdfMolecule
// ===========================================================================

///|
test "chem_molecule_creation" {
  let mol = @src.chem_sample_water()
  assert_eq(mol.name(), "water")
  assert_eq(mol.n_atoms(), 3)
  assert_eq(mol.n_bonds(), 2)
}

///|
test "chem_molecule_properties" {
  let mol = @src.chem_sample_ethanol()
  assert_true(mol.has_property("MW"))
  assert_eq(mol.get_property("MW"), "46.069")
  assert_false(mol.has_property("XYZ"))
  assert_eq(mol.get_property("XYZ"), "")
}

///|
test "chem_molecule_set_property" {
  let mol = @src.chem_sample_water()
  mol.set_property("pKa", "15.7")
  assert_true(mol.has_property("pKa"))
  assert_eq(mol.get_property("pKa"), "15.7")
}

// ===========================================================================
// 4. Molecular descriptors
// ===========================================================================

///|
test "chem_molecular_weight_water" {
  let mol = @src.chem_sample_water()
  let mw = @src.chem_sdf_molecular_weight(mol)
  // H2O: 2*1.008 + 15.999 = 18.015
  assert_true(mw > 18.0)
  assert_true(mw < 18.1)
}

///|
test "chem_molecular_weight_methane" {
  let mol = @src.chem_sample_methane()
  let mw = @src.chem_sdf_molecular_weight(mol)
  // CH4: 12.011 + 4*1.008 = 16.043
  assert_true(mw > 16.0)
  assert_true(mw < 16.1)
}

///|
test "chem_molecular_formula_water" {
  let mol = @src.chem_sample_water()
  let formula = @src.chem_molecular_formula(mol)
  assert_eq(formula, "H2O")
}

///|
test "chem_molecular_formula_methane" {
  let mol = @src.chem_sample_methane()
  let formula = @src.chem_molecular_formula(mol)
  assert_eq(formula, "CH4")
}

///|
test "chem_molecular_formula_ethanol" {
  let mol = @src.chem_sample_ethanol()
  let formula = @src.chem_molecular_formula(mol)
  assert_eq(formula, "C2H6O")
}

///|
test "chem_molecular_formula_benzene" {
  let mol = @src.chem_sample_benzene()
  let formula = @src.chem_molecular_formula(mol)
  assert_eq(formula, "C6H6")
}

///|
test "chem_heavy_atom_count_water" {
  let mol = @src.chem_sample_water()
  assert_eq(@src.chem_heavy_atom_count(mol), 1)
}

///|
test "chem_heavy_atom_count_ethanol" {
  let mol = @src.chem_sample_ethanol()
  // C2H6O: 3 heavy atoms (2C + 1O)
  assert_eq(@src.chem_heavy_atom_count(mol), 3)
}

///|
test "chem_count_element_ethanol" {
  let mol = @src.chem_sample_ethanol()
  assert_eq(@src.chem_count_element(mol, "C"), 2)
  assert_eq(@src.chem_count_element(mol, "H"), 6)
  assert_eq(@src.chem_count_element(mol, "O"), 1)
  assert_eq(@src.chem_count_element(mol, "N"), 0)
}

///|
test "chem_ring_count_benzene" {
  let mol = @src.chem_sample_benzene()
  // Benzene has 1 ring.
  assert_eq(@src.chem_ring_count(mol), 1)
}

///|
test "chem_ring_count_acyclic" {
  let mol = @src.chem_sample_ethanol()
  // Ethanol has 0 rings.
  assert_eq(@src.chem_ring_count(mol), 0)
}

///|
test "chem_ring_count_aspirin" {
  let mol = @src.chem_sample_aspirin()
  // Aspirin has 1 ring (benzene).
  assert_eq(@src.chem_ring_count(mol), 1)
}

// ===========================================================================
// 5. Atom-pair fingerprints
// ===========================================================================

///|
test "chem_fingerprint_water" {
  let mol = @src.chem_sample_water()
  let fp = @src.chem_atom_pair_fingerprint(mol)
  // Water has 3 atoms but only 2 unique atom-pair hashes (O-H distance
  // appears twice), so the fingerprint set collapses to 2 entries.
  assert_eq(fp.length(), 2)
}

///|
test "chem_fingerprint_methane" {
  let mol = @src.chem_sample_methane()
  let fp = @src.chem_atom_pair_fingerprint(mol)
  // Methane has 5 atoms (10 pairs) but only 2 unique atom-pair hashes
  // (C-H at one distance, H-H at another), so the fingerprint set
  // collapses to 2 entries.
  assert_eq(fp.length(), 2)
}

///|
test "chem_fingerprint_sorted" {
  let mol = @src.chem_sample_ethanol()
  let fp = @src.chem_atom_pair_fingerprint(mol)
  // Verify sorted ascending.
  for i in 1..<fp.length() {
    assert_true(fp[i] >= fp[i - 1])
  }
}

///|
test "chem_fingerprint_empty_molecule" {
  let atoms : Array[@src.SdfAtom] = []
  let bonds : Array[@src.SdfBond] = []
  let mol = @src.SdfMolecule::new("empty", atoms, bonds)
  let fp = @src.chem_atom_pair_fingerprint(mol)
  assert_eq(fp.length(), 0)
}

// ===========================================================================
// 6. Similarity coefficients
// ===========================================================================

///|
test "chem_tanimoto_identical_fingerprints" {
  let fp = [1, 2, 3, 4, 5]
  let sim = @src.chem_tanimoto(fp, fp)
  assert_true(sim == 1.0)
}

///|
test "chem_tanimoto_disjoint_fingerprints" {
  let fp1 = [1, 2, 3]
  let fp2 = [4, 5, 6]
  let sim = @src.chem_tanimoto(fp1, fp2)
  assert_true(sim == 0.0)
}

///|
test "chem_tanimoto_partial_overlap" {
  let fp1 = [1, 2, 3, 4]
  let fp2 = [3, 4, 5, 6]
  let sim = @src.chem_tanimoto(fp1, fp2)
  // Intersection = 2, Union = 6, Tanimoto = 2/6 = 0.333...
  assert_true(sim > 0.33)
  assert_true(sim < 0.34)
}

///|
test "chem_dice_identical_fingerprints" {
  let fp = [1, 2, 3]
  let sim = @src.chem_dice(fp, fp)
  assert_true(sim == 1.0)
}

///|
test "chem_dice_partial_overlap" {
  let fp1 = [1, 2, 3, 4]
  let fp2 = [3, 4, 5, 6]
  let sim = @src.chem_dice(fp1, fp2)
  // Intersection = 2, total = 8, Dice = 4/8 = 0.5
  assert_true(sim == 0.5)
}

///|
test "chem_tanimoto_empty_fingerprints" {
  let sim = @src.chem_tanimoto([], [])
  assert_true(sim == 0.0)
}

///|
test "chem_similarity_methane_vs_methane" {
  let mol = @src.chem_sample_methane()
  let fp = @src.chem_atom_pair_fingerprint(mol)
  let sim = @src.chem_tanimoto(fp, fp)
  assert_true(sim == 1.0)
}

// ===========================================================================
// 7. Substructure search
// ===========================================================================

///|
test "chem_substructure_self_match" {
  let mol = @src.chem_sample_ethanol()
  assert_true(@src.chem_substructure_match(mol, mol))
}

///|
test "chem_substructure_water_in_ethanol" {
  // Water (O-H) should be a substructure of ethanol (has O-H bond).
  let query = @src.chem_sample_water()
  let target = @src.chem_sample_ethanol()
  assert_true(@src.chem_substructure_match(query, target))
}

///|
test "chem_substructure_benzene_in_aspirin" {
  // Benzene ring should be a substructure of aspirin.
  let query = @src.chem_sample_benzene()
  let target = @src.chem_sample_aspirin()
  assert_true(@src.chem_substructure_match(query, target))
}

///|
test "chem_substructure_aspirin_not_in_benzene" {
  // Aspirin is larger than benzene; should not match.
  let query = @src.chem_sample_aspirin()
  let target = @src.chem_sample_benzene()
  assert_false(@src.chem_substructure_match(query, target))
}

///|
test "chem_substructure_empty_query" {
  let target = @src.chem_sample_water()
  let atoms : Array[@src.SdfAtom] = []
  let bonds : Array[@src.SdfBond] = []
  let query = @src.SdfMolecule::new("empty", atoms, bonds)
  assert_true(@src.chem_substructure_match(query, target))
}

// ===========================================================================
// 8. SDF parsing
// ===========================================================================

///|
test "chem_parse_sdf_single_molecule" {
  let sdf = @src.chem_sample_sdf_text()
  let mols = @src.chem_parse_sdf(sdf)
  assert_eq(mols.length(), 2)
  // First molecule is ethanol.
  assert_eq(mols[0].name(), "ethanol")
  assert_eq(mols[0].n_atoms(), 9)
  assert_eq(mols[0].n_bonds(), 8)
}

///|
test "chem_parse_sdf_preserves_properties" {
  let sdf = @src.chem_sample_sdf_text()
  let mols = @src.chem_parse_sdf(sdf)
  // Ethanol has MW property.
  assert_true(mols[0].has_property("MW"))
  assert_eq(mols[0].get_property("MW"), "46.069")
}

///|
test "chem_parse_sdf_methane" {
  let sdf = @src.chem_sample_sdf_text()
  let mols = @src.chem_parse_sdf(sdf)
  // Second molecule is methane.
  assert_eq(mols[1].name(), "methane")
  assert_eq(mols[1].n_atoms(), 5)
  assert_eq(mols[1].n_bonds(), 4)
}

///|
test "chem_parse_sdf_empty" {
  let mols = @src.chem_parse_sdf("")
  assert_eq(mols.length(), 0)
}

// ===========================================================================
// 9. SDF serialization
// ===========================================================================

///|
test "chem_molecule_to_sdf_basic" {
  let mol = @src.chem_sample_water()
  let sdf = @src.chem_molecule_to_sdf(mol)
  assert_true(sdf.contains("water"))
  assert_true(sdf.contains("V2000"))
  assert_true(sdf.contains("$$$$"))
}

///|
test "chem_molecules_to_sdf_multiple" {
  let mols = [@src.chem_sample_water(), @src.chem_sample_methane()]
  let sdf = @src.chem_molecules_to_sdf(mols)
  // Should contain both molecule names and two $$$$ terminators.
  assert_true(sdf.contains("water"))
  assert_true(sdf.contains("methane"))
}

///|
test "chem_sdf_round_trip" {
  let original = @src.chem_sample_ethanol()
  let sdf = @src.chem_molecule_to_sdf(original)
  let parsed = @src.chem_parse_sdf(sdf)
  assert_eq(parsed.length(), 1)
  assert_eq(parsed[0].name(), "ethanol")
  assert_eq(parsed[0].n_atoms(), 9)
  assert_eq(parsed[0].n_bonds(), 8)
  assert_true(parsed[0].has_property("MW"))
}

// ===========================================================================
// 10. Summary
// ===========================================================================

///|
test "chem_molecule_summary_water" {
  let mol = @src.chem_sample_water()
  let summary = @src.chem_molecule_summary(mol)
  assert_true(summary.contains("water"))
  assert_true(summary.contains("H2O"))
  assert_true(summary.contains("Atoms: 3"))
  assert_true(summary.contains("Bonds: 2"))
}

///|
test "chem_molecule_summary_benzene" {
  let mol = @src.chem_sample_benzene()
  let summary = @src.chem_molecule_summary(mol)
  assert_true(summary.contains("benzene"))
  assert_true(summary.contains("C6H6"))
  assert_true(summary.contains("Rings: 1"))
}

// ===========================================================================
// 11. Sample data validation
// ===========================================================================

///|
test "chem_sample_aspirin_structure" {
  let mol = @src.chem_sample_aspirin()
  assert_eq(mol.name(), "aspirin")
  assert_eq(mol.n_atoms(), 21)
  assert_eq(mol.n_bonds(), 21)
  assert_true(mol.has_property("DRUG"))
  assert_eq(mol.get_property("DRUG"), "true")
}

///|
test "chem_sample_aspirin_formula" {
  let mol = @src.chem_sample_aspirin()
  let formula = @src.chem_molecular_formula(mol)
  assert_eq(formula, "C9H8O4")
}

///|
test "chem_sample_aspirin_ring_count" {
  let mol = @src.chem_sample_aspirin()
  // Aspirin has 1 ring (benzene).
  assert_eq(@src.chem_ring_count(mol), 1)
}

///|
test "chem_sample_aspirin_heavy_atoms" {
  let mol = @src.chem_sample_aspirin()
  // C9H8O4: 13 heavy atoms (9C + 4O).
  assert_eq(@src.chem_heavy_atom_count(mol), 13)
}

// ===========================================================================
// 12. Edge cases
// ===========================================================================

///|
test "chem_fingerprint_single_atom" {
  let atoms = [@src.SdfAtom::new("C", 0.0, 0.0, 0.0, 0)]
  let bonds : Array[@src.SdfBond] = []
  let mol = @src.SdfMolecule::new("single", atoms, bonds)
  let fp = @src.chem_atom_pair_fingerprint(mol)
  // Single atom has no pairs.
  assert_eq(fp.length(), 0)
}

///|
test "chem_tanimoto_different_sizes" {
  let fp1 = [1, 2, 3]
  let fp2 = [1, 2, 3, 4, 5, 6]
  let sim = @src.chem_tanimoto(fp1, fp2)
  // Intersection = 3, Union = 6, Tanimoto = 0.5
  assert_true(sim == 0.5)
}

///|
test "chem_count_element_case_insensitive" {
  let mol = @src.chem_sample_water()
  assert_eq(@src.chem_count_element(mol, "o"), 1)
  assert_eq(@src.chem_count_element(mol, "O"), 1)
  assert_eq(@src.chem_count_element(mol, "h"), 2)
  assert_eq(@src.chem_count_element(mol, "H"), 2)
}
