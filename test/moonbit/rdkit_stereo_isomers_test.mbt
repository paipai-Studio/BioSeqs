///|
/// Tests for rdkit_stereo_isomers.mbt — potential center detection,
/// isomer counting/enumeration and annotated stereo SMILES.

///|
test "potential_centers_lactic_acid" {
  // CC(O)C(=O)O: center atom 1 has neighbours C,O,C with distinct CIP
  // levels (CH3 vs COOH differ by degree)
  let mol = @src.parse_smiles("CC(O)C(=O)O")
  let centers = @src.rdkit_find_potential_stereocenters(mol)
  assert_eq(centers, [1])
}

///|
test "potential_centers_none_achiral" {
  // 2-butanol? no — use ethanol (no center) and dichloromethane-like
  // CC(Cl)C is chiral; propane is not
  let mol = @src.parse_smiles("CCC")
  assert_eq(@src.rdkit_find_potential_stereocenters(mol).length(), 0)
  // symmetric centre: C(Cl)(Cl)C — neighbours Cl,Cl,C with repeats
  let mol2 = @src.parse_smiles("CC(Cl)Cl")
  assert_eq(@src.rdkit_find_potential_stereocenters(mol2).length(), 0)
}

///|
test "count_stereoisomers_two_centers" {
  // 2,3-dichlorobutane: two potential centers -> 4 graph combinations
  let mol = @src.parse_smiles("CC(Cl)C(Cl)C")
  assert_eq(@src.rdkit_count_stereoisomers(mol), 4)
}

///|
test "enumerate_stereoisomers_lactic_full" {
  let mol = @src.parse_smiles("CC(O)C(=O)O")
  let isomers = @src.rdkit_enumerate_stereoisomers(mol)
  assert_eq(isomers.length(), 2)
  assert_eq(isomers[0].label, "R")
  assert_eq(isomers[1].label, "S")
  // assignments point at the center atom
  assert_eq(isomers[0].assignments[0].0, 1)
  assert_true(isomers[0].assignments[0].1)
  assert_true(!isomers[1].assignments[0].1)
}

///|
test "enumerate_stereoisomers_annotated_smiles" {
  let mol = @src.parse_smiles("CC(O)C(=O)O")
  let isomers = @src.rdkit_enumerate_stereoisomers(mol)
  // the two annotations must differ and carry the H marker
  assert_true(isomers[0].smiles.contains("@"))
  assert_true(isomers[0].smiles.contains("[C"))
  assert_true(isomers[0].smiles != isomers[1].smiles)
}

///|
test "enumerate_stereoisomers_max_isomers_cap" {
  let mol = @src.parse_smiles("CC(Cl)C(Cl)C")
  let capped = @src.rdkit_enumerate_stereoisomers(mol, max_isomers=3)
  assert_eq(capped.length(), 3)
  let full = @src.rdkit_enumerate_stereoisomers(mol)
  assert_eq(full.length(), 4)
  // labels are ordered R-first
  assert_eq(full[0].label, "R,R")
  assert_eq(full[3].label, "S,S")
}

///|
test "stereo_smiles_degree4_center" {
  // C(F)(Cl)(Br)I: degree-4 center -> [C@] / [C@@] without H
  let mol = @src.parse_smiles("C(F)(Cl)(Br)I")
  let centers = @src.rdkit_find_potential_stereocenters(mol)
  assert_eq(centers.length(), 1)
  let isomers = @src.rdkit_enumerate_stereoisomers(mol)
  assert_eq(isomers.length(), 2)
  assert_true(isomers[0].smiles.contains("[C@]"))
  assert_true(isomers[1].smiles.contains("[C@@]"))
}

///|
test "count_stereoisomers_achiral_returns_one" {
  // no potential centers -> the molecule itself (RDKit convention)
  assert_eq(@src.rdkit_count_stereoisomers(@src.parse_smiles("CCO")), 1)
}
