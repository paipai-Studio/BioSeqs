test "atom_mass_basic" {
  // Hydrogen
  let h = @src.atom_mass("H")
  assert_true(h > 1.0 && h < 1.1)
  // Carbon
  let c = @src.atom_mass("C")
  assert_true(c > 11.9 && c < 12.1)
  // Nitrogen
  let n = @src.atom_mass("N")
  assert_true(n > 13.9 && n < 14.1)
  // Oxygen
  let o = @src.atom_mass("O")
  assert_true(o > 15.9 && o < 16.1)
  // Unknown
  let unk = @src.atom_mass("XX")
  assert_eq(unk, 0.0)
  // Case insensitive
  assert_eq(@src.atom_mass("c"), @src.atom_mass("C"))
}

test "atom_vdw_radius_basic" {
  // Carbon VdW radius ≈ 1.70 Å
  let c = @src.atom_vdw_radius("C")
  assert_true(c > 1.6 && c < 1.8)
  // Hydrogen VdW radius ≈ 1.20 Å
  let h = @src.atom_vdw_radius("H")
  assert_true(h > 1.1 && h < 1.3)
  // Unknown
  assert_eq(@src.atom_vdw_radius("XX"), 0.0)
}

test "covalent_radius_basic" {
  // Carbon covalent radius ≈ 0.76 Å
  let c = @src.covalent_radius("C")
  assert_true(c > 0.7 && c < 0.8)
  // Hydrogen covalent radius ≈ 0.31 Å
  let h = @src.covalent_radius("H")
  assert_true(h > 0.25 && h < 0.35)
  // Unknown
  assert_eq(@src.covalent_radius("XX"), 0.0)
}

test "formal_charge_basic" {
  // ASP OD1 → -1
  assert_eq(@src.formal_charge("ASP", "OD1"), -1)
  // LYS NZ → +1
  assert_eq(@src.formal_charge("LYS", "NZ"), 1)
  // ALA CA → 0
  assert_eq(@src.formal_charge("ALA", "CA"), 0)
  // OXT → -1
  assert_eq(@src.formal_charge("GLY", "OXT"), -1)
}
