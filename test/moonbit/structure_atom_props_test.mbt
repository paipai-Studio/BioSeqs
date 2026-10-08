test "atom_mass_basic" {
  // Standard atomic masses
  assert_eq(@src.atom_mass("C"), 12.011)
  assert_eq(@src.atom_mass("N"), 14.007)
  assert_eq(@src.atom_mass("O"), 15.999)
  assert_eq(@src.atom_mass("H"), 1.008)
  assert_eq(@src.atom_mass("S"), 32.06)
  assert_eq(@src.atom_mass("P"), 30.974)
}

test "atom_mass_canonicalization" {
  // Should canonicalize case
  assert_eq(@src.atom_mass("c"), 12.011)
  assert_eq(@src.atom_mass("CL"), 35.45)
  assert_eq(@src.atom_mass("cl"), 35.45)
  assert_eq(@src.atom_mass("Fe"), 55.845)
  assert_eq(@src.atom_mass("FE"), 55.845)
}

test "atom_mass_unknown" {
  assert_eq(@src.atom_mass("Xx"), 0.0)
}

test "atom_vdw_radius_basic" {
  assert_eq(@src.atom_vdw_radius("C"), 1.7)
  assert_eq(@src.atom_vdw_radius("N"), 1.55)
  assert_eq(@src.atom_vdw_radius("O"), 1.52)
  assert_eq(@src.atom_vdw_radius("H"), 1.2)
  assert_eq(@src.atom_vdw_radius("S"), 1.8)
}

test "covalent_radius_basic" {
  assert_eq(@src.covalent_radius("C"), 0.76)
  assert_eq(@src.covalent_radius("N"), 0.71)
  assert_eq(@src.covalent_radius("O"), 0.66)
  assert_eq(@src.covalent_radius("H"), 0.31)
}

test "formal_charge_asp" {
  // Aspartate carboxylate oxygens carry -1
  assert_eq(@src.formal_charge("ASP", "OD1"), -1)
  assert_eq(@src.formal_charge("ASP", "OD2"), -1)
}

test "formal_charge_lys" {
  // Lysine NZ carries +1
  assert_eq(@src.formal_charge("LYS", "NZ"), 1)
}

test "formal_charge_glu" {
  assert_eq(@src.formal_charge("GLU", "OE1"), -1)
  assert_eq(@src.formal_charge("GLU", "OE2"), -1)
}

test "formal_charge_neutral" {
  // Backbone atoms are neutral
  assert_eq(@src.formal_charge("ALA", "CA"), 0)
  assert_eq(@src.formal_charge("GLY", "N"), 0)
}
