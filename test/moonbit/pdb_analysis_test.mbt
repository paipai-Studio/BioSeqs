///|
test "pdb_analysis: calc_dihedral" {
  let a = @src.Vector3::new(0.0, 0.0, 0.0)
  let b = @src.Vector3::new(1.0, 0.0, 0.0)
  let c = @src.Vector3::new(1.0, 1.0, 0.0)
  let d = @src.Vector3::new(1.0, 1.0, 1.0)
  let angle = @src.calc_dihedral(a, b, c, d)
  let abs_angle = if angle < 0.0 { -angle } else { angle }
  assert_true(abs_angle > 80.0 && abs_angle < 100.0)
}

///|
test "pdb_analysis: SSType to_code and from_code" {
  assert_eq(@src.ss_alpha_helix().to_code(), 'H')
  assert_eq(@src.ss_extended_strand().to_code(), 'E')
  assert_eq(@src.ss_coil().to_code(), '-')
  assert_eq(@src.SSType::from_code('H'), @src.ss_alpha_helix())
  assert_eq(@src.SSType::from_code('E'), @src.ss_extended_strand())
  assert_eq(@src.SSType::from_code('-'), @src.ss_coil())
}

///|
test "pdb_analysis: create_demo_structure" {
  let s = @src.create_demo_structure()
  assert_eq(s.get_num_models(), 1)
  assert_eq(s.get_num_chains(), 1)
  assert_eq(s.get_num_residues(), 4)
  assert_eq(s.get_num_atoms(), 16)
}

///|
test "pdb_analysis: calc_chain_dihedrals" {
  let s = @src.create_demo_structure()
  let model = s.find_model(0)
  assert_true(model is Some(_))
  let chain = model.unwrap().find_chain('A')
  assert_true(chain is Some(_))
  let dihedrals = @src.calc_chain_dihedrals(chain.unwrap())
  assert_eq(dihedrals.length(), 4)
  assert_true(dihedrals[0].1 is None)
  assert_true(dihedrals[1].1 is Some(_))
  assert_true(dihedrals[3].2 is None)
}

///|
test "pdb_analysis: calc_ca_distance_matrix" {
  let s = @src.create_demo_structure()
  let chain = s.find_model(0).unwrap().find_chain('A').unwrap()
  let matrix = @src.calc_ca_distance_matrix(chain)
  assert_eq(matrix.length(), 4)
  assert_eq(matrix[0].length(), 4)
  assert_true(matrix[0][0] < 0.001)
  assert_true(matrix[1][1] < 0.001)
  assert_true(matrix[0][3] > 0.0)
}

///|
test "pdb_analysis: calc_contact_map" {
  let s = @src.create_demo_structure()
  let chain = s.find_model(0).unwrap().find_chain('A').unwrap()
  let contact = @src.calc_contact_map(chain, threshold=20.0)
  assert_eq(contact.length(), 4)
  assert_true(contact[0][1])
  assert_false(contact[0][0])
}

///|
test "pdb_analysis: detect_hydrogen_bonds" {
  let s = @src.create_demo_structure()
  let chain = s.find_model(0).unwrap().find_chain('A').unwrap()
  let hbonds = @src.detect_hydrogen_bonds(chain, max_dist=5.0, min_angle=0.0)
  assert_true(hbonds.length() > 0)
}

///|
test "pdb_analysis: assign_secondary_structure" {
  let s = @src.create_demo_structure()
  let chain = s.find_model(0).unwrap().find_chain('A').unwrap()
  let ss = @src.assign_secondary_structure(chain)
  assert_eq(ss.length(), 4)
}

///|
test "pdb_analysis: count_secondary_structure" {
  let ss : Array[@src.SSType] = [
    @src.ss_alpha_helix(),
    @src.ss_alpha_helix(),
    @src.ss_coil(),
    @src.ss_extended_strand(),
  ]
  let counts = @src.count_secondary_structure(ss.map(fn(s) { (1, s) }))
  assert_eq(counts.get("H"), Some(2))
  assert_eq(counts.get("-"), Some(1))
  assert_eq(counts.get("E"), Some(1))
}

///|
test "pdb_analysis: calc_radius_of_gyration" {
  let s = @src.create_demo_structure()
  let chain = s.find_model(0).unwrap().find_chain('A').unwrap()
  let rg = @src.calc_radius_of_gyration(chain)
  assert_true(rg is Some(_))
  assert_true(rg.unwrap() > 0.0)
}

///|
test "pdb_analysis: ramachandran_plot" {
  let s = @src.create_demo_structure()
  let chain = s.find_model(0).unwrap().find_chain('A').unwrap()
  let plot = @src.ramachandran_plot(chain)
  assert_true(plot.length() > 0)
}

///|
test "pdb_analysis: classify_ramachandran" {
  assert_eq(
    @src.classify_ramachandran(-60.0, -45.0),
    @src.ramachandran_favored(),
  )
  assert_eq(
    @src.classify_ramachandran(-120.0, 130.0),
    @src.ramachandran_favored(),
  )
  assert_eq(@src.classify_ramachandran(60.0, 30.0), @src.ramachandran_allowed())
}

///|
test "pdb_analysis: vdw_radius" {
  assert_true((@src.vdw_radius("C") - 1.7).abs() < 0.001)
  assert_true((@src.vdw_radius("CA") - 1.7).abs() < 0.001)
  assert_true((@src.vdw_radius("N") - 1.55).abs() < 0.001)
  assert_true((@src.vdw_radius("O") - 1.52).abs() < 0.001)
  assert_true((@src.vdw_radius("S") - 1.8).abs() < 0.001)
}

///|
test "pdb_analysis: calc_sasa" {
  let s = @src.create_demo_structure()
  let chain = s.find_model(0).unwrap().find_chain('A').unwrap()
  let sasa = @src.calc_sasa(chain, probe_radius=1.4, n_points=50)
  assert_true(sasa.length() > 0)
  // All SASA values should be non-negative
  for entry in sasa {
    assert_true(entry.1 >= 0.0)
  }
}

///|
test "pdb_analysis: calc_total_sasa" {
  let s = @src.create_demo_structure()
  let chain = s.find_model(0).unwrap().find_chain('A').unwrap()
  let total = @src.calc_total_sasa(chain)
  assert_true(total > 0.0)
}

///|
test "pdb_analysis: calc_residue_sasa" {
  let s = @src.create_demo_structure()
  let chain = s.find_model(0).unwrap().find_chain('A').unwrap()
  let res_sasa = @src.calc_residue_sasa(chain, probe_radius=1.4)
  assert_eq(res_sasa.length(), 4)
  // Residue 1 should be ALA
  assert_eq(res_sasa[0].1, "ALA")
  // All SASA values should be non-negative
  for entry in res_sasa {
    assert_true(entry.2 >= 0.0)
  }
}

///|
test "pdb_analysis: ramachandran_quality" {
  let s = @src.create_demo_structure()
  let chain = s.find_model(0).unwrap().find_chain('A').unwrap()
  let (favored, allowed, generous) = @src.ramachandran_quality(chain)
  // Fractions should sum to approximately 1
  let total = favored + allowed + generous
  assert_true(total > 0.99 && total < 1.01)
  assert_true(favored >= 0.0 && favored <= 1.0)
  assert_true(allowed >= 0.0 && allowed <= 1.0)
  assert_true(generous >= 0.0 && generous <= 1.0)
}

///|
test "pdb_analysis: get_hydrophobicity" {
  // Test some known values from Kyte-Doolittle scale
  assert_true((@src.get_hydrophobicity("ALA") - 1.8).abs() < 0.001)
  assert_true((@src.get_hydrophobicity("ARG") - (-4.5)).abs() < 0.001)
  assert_true((@src.get_hydrophobicity("ILE") - 4.5).abs() < 0.001)
  assert_true((@src.get_hydrophobicity("LYS") - (-3.9)).abs() < 0.001)
  assert_true((@src.get_hydrophobicity("PHE") - 2.8).abs() < 0.001)
  assert_true((@src.get_hydrophobicity("VAL") - 4.2).abs() < 0.001)
  // Unknown residue should return 0
  assert_eq(@src.get_hydrophobicity("XXX"), 0.0)
}

///|
test "pdb_analysis: calc_sequence_property_matrix" {
  let s = @src.create_demo_structure()
  let chain = s.find_model(0).unwrap().find_chain('A').unwrap()
  let matrix = @src.calc_sequence_property_matrix(chain)
  assert_eq(matrix.length(), 4)
  assert_eq(matrix[0].length(), 4)
  // Diagonal should be zero (self comparison)
  assert_true(matrix[0][0].abs() < 0.001)
  assert_true(matrix[1][1].abs() < 0.001)
  // Symmetric matrix
  assert_true((matrix[0][1] - matrix[1][0]).abs() < 0.001)
}
