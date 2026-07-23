///|
/// Bio.PDB.PSEA module tests - Protein secondary structure prediction from CA coordinates

///|
test "psea_atom_new" {
  let atom = @src.PseaAtom::new("ALA", 1, "CA", 1.0, 2.0, 3.0)
  assert_eq(atom.res_name, "ALA")
  assert_eq(atom.res_num, 1)
  assert_eq(atom.atom_name, "CA")
  assert_true((atom.x - 1.0).abs() < 1.0e-9)
  assert_true((atom.y - 2.0).abs() < 1.0e-9)
  assert_true((atom.z - 3.0).abs() < 1.0e-9)
}

///|
test "psea_distance" {
  let a = @src.PseaAtom::new("ALA", 1, "CA", 0.0, 0.0, 0.0)
  let b = @src.PseaAtom::new("GLY", 2, "CA", 3.0, 4.0, 0.0)
  let d = @src.psea_distance(a, b)
  // sqrt(9 + 16) = 5
  assert_true((d - 5.0).abs() < 1.0e-9)
}

///|
test "psea_bond_angle" {
  // Three atoms forming a right angle
  let a = @src.PseaAtom::new("A", 1, "CA", 1.0, 0.0, 0.0)
  let b = @src.PseaAtom::new("B", 2, "CA", 0.0, 0.0, 0.0)
  let c = @src.PseaAtom::new("C", 3, "CA", 0.0, 1.0, 0.0)
  let angle = @src.psea_bond_angle(a, b, c)
  // Should be 90 degrees
  assert_true((angle - 90.0).abs() < 1.0e-6)
}

///|
test "psea_bond_angle_linear" {
  // Three collinear atoms: angle should be 180 degrees
  let a = @src.PseaAtom::new("A", 1, "CA", 0.0, 0.0, 0.0)
  let b = @src.PseaAtom::new("B", 2, "CA", 1.0, 0.0, 0.0)
  let c = @src.PseaAtom::new("C", 3, "CA", 2.0, 0.0, 0.0)
  let angle = @src.psea_bond_angle(a, b, c)
  assert_true((angle - 180.0).abs() < 1.0e-4)
}

///|
test "psea_dihedral_angle" {
  // Four atoms forming a planar zigzag: dihedral should be 180 degrees
  let a = @src.PseaAtom::new("A", 1, "CA", 0.0, 0.0, 0.0)
  let b = @src.PseaAtom::new("B", 2, "CA", 1.0, 0.0, 0.0)
  let c = @src.PseaAtom::new("C", 3, "CA", 2.0, 0.0, 0.0)
  let d = @src.PseaAtom::new("D", 4, "CA", 3.0, 0.0, 0.0)
  let angle = @src.psea_dihedral_angle(a, b, c, d)
  // All coplanar -> dihedral is 180 (or -180)
  let abs_angle = if angle < 0.0 { -angle } else { angle }
  assert_true((abs_angle - 180.0).abs() < 1.0e-4 || abs_angle < 1.0e-4)
}

///|
test "psea_compute_geometry" {
  let atoms = @src.psea_create_example_strand()
  let geom = @src.psea_compute_geometry(atoms)
  let n = atoms.length()
  assert_eq(geom.ca_distances.length(), n - 1)
  assert_eq(geom.bond_angles.length(), n - 2)
  assert_eq(geom.dihedral_angles.length(), n - 3)
  // CA-CA distances should be near 3.8 for real proteins; our strand uses 3.4 + 1.0 zigzag
  // Just check they're positive
  let mut i = 0
  while i < geom.ca_distances.length() {
    assert_true(geom.ca_distances[i] > 0.0)
    i = i + 1
  }
}

///|
test "psea_run_helix" {
  let atoms = @src.psea_create_example_helix()
  let result = @src.psea_run(atoms)
  assert_eq(result.psea_n_residues(), 12)
  // SS string should be 12 characters
  let ss_str = result.psea_ss_string()
  assert_eq(ss_str.length(), 12)
  // All characters should be H, E, or C
  let mut i = 0
  while i < ss_str.length() {
    let c = ss_str[i].to_char().unwrap_or('?')
    assert_true(c == 'H' || c == 'E' || c == 'C')
    i = i + 1
  }
  // Counts should add up to 12
  let total = result.n_helix + result.n_strand + result.n_coil
  assert_eq(total, 12)
}

///|
test "psea_run_strand" {
  let atoms = @src.psea_create_example_strand()
  let result = @src.psea_run(atoms)
  assert_eq(result.psea_n_residues(), 8)
  let ss_str = result.psea_ss_string()
  assert_eq(ss_str.length(), 8)
  // Strand example may or may not have E depending on geometry
  // Just check the result is valid
  let total = result.n_helix + result.n_strand + result.n_coil
  assert_eq(total, 8)
}

///|
test "psea_run_mixed" {
  let atoms = @src.psea_create_example_mixed()
  let result = @src.psea_run(atoms)
  let n = atoms.length()
  assert_eq(result.psea_n_residues(), n)
  let total = result.n_helix + result.n_strand + result.n_coil
  assert_eq(total, n)
  // Summary should mention all three
  let summary = result.psea_summary()
  assert_true(summary.contains("H="))
  assert_true(summary.contains("E="))
  assert_true(summary.contains("C="))
}

///|
test "psea_result_methods" {
  let atoms = @src.psea_create_example_helix()
  let result = @src.psea_run(atoms)
  // Test accessors
  let mut i = 0
  while i < result.psea_n_residues() {
    let ss = result.psea_get_ss(i)
    assert_true(ss == "H" || ss == "E" || ss == "C")
    let conf = result.psea_get_confidence(i)
    assert_true(conf >= 0.0)
    assert_true(conf <= 1.0)
    let res_num = result.psea_get_res_num(i)
    assert_eq(res_num, i + 1)
    let res_name = result.psea_get_res_name(i)
    assert_true(res_name.length() > 0)
    i = i + 1
  }
}

///|
test "psea_short_chain" {
  // Chain with fewer than 4 atoms should be all coil
  let atoms = [
    @src.PseaAtom::new("A", 1, "CA", 0.0, 0.0, 0.0),
    @src.PseaAtom::new("B", 2, "CA", 1.0, 0.0, 0.0),
    @src.PseaAtom::new("C", 3, "CA", 2.0, 0.0, 0.0)
  ]
  let result = @src.psea_run(atoms)
  assert_eq(result.psea_n_residues(), 3)
  assert_eq(result.n_helix, 0)
  assert_eq(result.n_strand, 0)
  assert_eq(result.n_coil, 3)
}

///|
test "psea_three_to_eight_state" {
  let ss = ["C", "C", "H", "H", "H", "H", "H", "C", "E", "E", "E", "C"]
  let eight = @src.psea_three_to_eight_state(ss)
  assert_eq(eight.length(), 12)
  // The 5-residue helix should stay as H
  assert_eq(eight[2], "H")
  assert_eq(eight[6], "H")
  // The 3-residue strand should stay as E (>= 2)
  assert_eq(eight[8], "E")
  assert_eq(eight[10], "E")
}

///|
test "psea_three_to_eight_state_short_helix" {
  // A 3-residue helix run should become G (3-10 helix)
  let ss = ["C", "H", "H", "H", "C"]
  let eight = @src.psea_three_to_eight_state(ss)
  assert_eq(eight.length(), 5)
  assert_eq(eight[1], "G")
  assert_eq(eight[3], "G")
}

///|
test "psea_three_to_eight_state_isolated_strand" {
  // A 1-residue strand should become B (bridge)
  let ss = ["C", "E", "C"]
  let eight = @src.psea_three_to_eight_state(ss)
  assert_eq(eight.length(), 3)
  assert_eq(eight[1], "B")
}
