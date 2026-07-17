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
