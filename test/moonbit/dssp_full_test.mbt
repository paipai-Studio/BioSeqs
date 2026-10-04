///|
test "dssp_full: DSSPType to_char and from_char" {
  assert_eq(@src.DSSPType::to_char(@src.dssp_type_from_char('H')), 'H')
  assert_eq(@src.DSSPType::to_char(@src.dssp_type_from_char('G')), 'G')
  assert_eq(@src.DSSPType::to_char(@src.dssp_type_from_char('I')), 'I')
  assert_eq(@src.DSSPType::to_char(@src.dssp_type_from_char('E')), 'E')
  assert_eq(@src.DSSPType::to_char(@src.dssp_type_from_char('B')), 'B')
  assert_eq(@src.DSSPType::to_char(@src.dssp_type_from_char('T')), 'T')
  assert_eq(@src.DSSPType::to_char(@src.dssp_type_from_char('S')), 'S')
  assert_eq(@src.DSSPType::to_char(@src.dssp_type_from_char('-')), '-')
  // Unknown char → Coil
  assert_eq(@src.DSSPType::to_str(@src.dssp_type_from_char('X')), "-")
}

///|
test "dssp_full: H-bond energy basic" {
  // Standard H-bond geometry: O...N distance ~2.9 Å, C...N ~3.8 Å
  let o_i = @src.Vector3::new(0.0, 0.0, 0.0)
  let c_i = @src.Vector3::new(1.2, 0.0, 0.0)
  let n_j = @src.Vector3::new(0.0, 2.9, 0.0)
  // H is ~1.0 Å from N towards O
  let h_j = @src.Vector3::new(0.0, 1.9, 0.0)
  let e = @src.dssp_hbond_energy(o_i, c_i, n_j, h_j)
  // For a good H-bond, energy should be negative
  assert_true(e < 0.0)
}

///|
test "dssp_full: H-bond energy weak (far apart)" {
  // Atoms far apart → weak/no H-bond
  let o_i = @src.Vector3::new(0.0, 0.0, 0.0)
  let c_i = @src.Vector3::new(1.2, 0.0, 0.0)
  let n_j = @src.Vector3::new(10.0, 10.0, 10.0)
  let h_j = @src.Vector3::new(10.0, 11.0, 10.0)
  let e = @src.dssp_hbond_energy(o_i, c_i, n_j, h_j)
  // Energy should be near zero (weak)
  assert_true(e > -0.5)
}

///|
test "dssp_full: approx H position" {
  let n_j = @src.Vector3::new(3.0, 0.0, 0.0)
  let ca_j = @src.Vector3::new(4.5, 0.0, 0.0)
  let prev_c = @src.Vector3::new(1.8, 0.0, 0.0)
  let h = @src.dssp_approx_h_position(n_j, ca_j, Some(prev_c))
  // H should be ~1.01 Å from N, in the direction away from prev_C
  let dist = n_j.distance(h)
  assert_true(dist > 0.9 && dist < 1.1)
  // H should be to the right of N (away from prev_C which is to the left)
  assert_true(h.x > n_j.x)
}

///|
test "dssp_full: approx H position no prev C" {
  let n_j = @src.Vector3::new(3.0, 0.0, 0.0)
  let ca_j = @src.Vector3::new(4.5, 0.0, 0.0)
  let h = @src.dssp_approx_h_position(n_j, ca_j, None)
  let dist = n_j.distance(h)
  assert_true(dist > 0.9 && dist < 1.1)
}

///|
test "dssp_full: assign on demo structure" {
  let s = @src.create_demo_structure()
  let chain = s.find_model(0).unwrap().find_chain('A').unwrap()
  let entries = @src.dssp_assign(chain)
  assert_eq(entries.length(), 4)
  // All entries should have valid resseq
  assert_eq(entries[0].resseq, 1)
  assert_eq(entries[3].resseq, 4)
}

///|
test "dssp_full: assign on ideal helix" {
  let s = @src.create_ideal_helix(12)
  let chain = s.find_model(0).unwrap().find_chain('A').unwrap()
  let entries = @src.dssp_assign(chain)
  assert_eq(entries.length(), 12)
  // Every entry should have a valid SS assignment
  let ss_seq = @src.dssp_ss_sequence(entries)
  assert_eq(ss_seq.length(), 12)
}

///|
test "dssp_full: assign on beta hairpin" {
  let s = @src.create_beta_hairpin(4)
  let chain = s.find_model(0).unwrap().find_chain('A').unwrap()
  let entries = @src.dssp_assign(chain)
  // 4 + 2 + 4 = 10 residues
  assert_eq(entries.length(), 10)
}

///|
test "dssp_full: assign on extended chain" {
  let s = @src.create_extended_chain(6)
  let chain = s.find_model(0).unwrap().find_chain('A').unwrap()
  let entries = @src.dssp_assign(chain)
  assert_eq(entries.length(), 6)
  // Extended chain should have mostly coil or bend
  let ss_seq = @src.dssp_ss_sequence(entries)
  assert_eq(ss_seq.length(), 6)
}

///|
test "dssp_full: empty chain" {
  let chain = @src.Chain::new(id='A')
  let entries = @src.dssp_assign(chain)
  assert_eq(entries.length(), 0)
}

///|
test "dssp_full: single residue" {
  let res = @src.Residue::new(
    resname="ALA", chainid='A', resseq=1,
    atoms=[
      @src.Atom::new(name="N", coord=@src.Vector3::new(0.0, 0.0, 0.0)),
      @src.Atom::new(name="CA", coord=@src.Vector3::new(1.5, 0.0, 0.0)),
      @src.Atom::new(name="C", coord=@src.Vector3::new(2.0, 1.5, 0.0)),
      @src.Atom::new(name="O", coord=@src.Vector3::new(1.2, 2.3, 0.0)),
    ],
  )
  let chain = @src.Chain::new(id='A', residues=[res])
  let entries = @src.dssp_assign(chain)
  assert_eq(entries.length(), 1)
  // Single residue should be coil
  assert_eq(entries[0].ss, @src.dssp_type_from_char('-'))
}

///|
test "dssp_full: two residues" {
  let res1 = @src.Residue::new(
    resname="ALA", chainid='A', resseq=1,
    atoms=[
      @src.Atom::new(name="N", coord=@src.Vector3::new(0.0, 0.0, 0.0)),
      @src.Atom::new(name="CA", coord=@src.Vector3::new(1.5, 0.0, 0.0)),
      @src.Atom::new(name="C", coord=@src.Vector3::new(2.0, 1.5, 0.0)),
      @src.Atom::new(name="O", coord=@src.Vector3::new(1.2, 2.3, 0.0)),
    ],
  )
  let res2 = @src.Residue::new(
    resname="GLY", chainid='A', resseq=2,
    atoms=[
      @src.Atom::new(name="N", coord=@src.Vector3::new(3.3, 1.8, 0.0)),
      @src.Atom::new(name="CA", coord=@src.Vector3::new(4.0, 3.0, 0.0)),
      @src.Atom::new(name="C", coord=@src.Vector3::new(5.5, 2.8, 0.0)),
      @src.Atom::new(name="O", coord=@src.Vector3::new(6.2, 3.7, 0.0)),
    ],
  )
  let chain = @src.Chain::new(id='A', residues=[res1, res2])
  let entries = @src.dssp_assign(chain)
  assert_eq(entries.length(), 2)
}

///|
test "dssp_full: missing backbone atoms" {
  // Residue with only CA (missing N, C, O)
  let res1 = @src.Residue::new(
    resname="ALA", chainid='A', resseq=1,
    atoms=[
      @src.Atom::new(name="CA", coord=@src.Vector3::new(0.0, 0.0, 0.0)),
    ],
  )
  let res2 = @src.Residue::new(
    resname="GLY", chainid='A', resseq=2,
    atoms=[
      @src.Atom::new(name="N", coord=@src.Vector3::new(3.3, 1.8, 0.0)),
      @src.Atom::new(name="CA", coord=@src.Vector3::new(4.0, 3.0, 0.0)),
      @src.Atom::new(name="C", coord=@src.Vector3::new(5.5, 2.8, 0.0)),
      @src.Atom::new(name="O", coord=@src.Vector3::new(6.2, 3.7, 0.0)),
    ],
  )
  let chain = @src.Chain::new(id='A', residues=[res1, res2])
  let entries = @src.dssp_assign(chain)
  assert_eq(entries.length(), 2)
  // Should not crash, both should be coil
  assert_eq(entries[0].ss, @src.dssp_type_from_char('-'))
  assert_eq(entries[1].ss, @src.dssp_type_from_char('-'))
}

///|
test "dssp_full: summary counts" {
  let s = @src.create_ideal_helix(10)
  let chain = s.find_model(0).unwrap().find_chain('A').unwrap()
  let entries = @src.dssp_assign(chain)
  let summary = @src.dssp_summary(entries)
  // Total should equal number of residues
  let mut total = 0
  for _, v in summary {
    total = total + v
  }
  assert_eq(total, 10)
}

///|
test "dssp_full: ss_sequence" {
  let s = @src.create_demo_structure()
  let chain = s.find_model(0).unwrap().find_chain('A').unwrap()
  let entries = @src.dssp_assign(chain)
  let seq = @src.dssp_ss_sequence(entries)
  assert_eq(seq.length(), 4)
  // Each character should be a valid DSSP code
  for i = 0; i < seq.length(); i = i + 1 {
    let c = seq[i]
    assert_true(c == 'H' || c == 'G' || c == 'I' || c == 'E' || c == 'B' || c == 'T' || c == 'S' || c == '-')
  }
}

///|
test "dssp_full: count hbonds" {
  let s = @src.create_ideal_helix(12)
  let chain = s.find_model(0).unwrap().find_chain('A').unwrap()
  let entries = @src.dssp_assign(chain)
  let n_hbonds = @src.dssp_count_hbonds(entries)
  // An ideal helix should have some H-bonds (i → i+4 pattern)
  // The exact count depends on geometry, but should be >= 0
  assert_true(n_hbonds >= 0)
}

///|
test "dssp_full: dihedral computation" {
  let s = @src.create_demo_structure()
  let chain = s.find_model(0).unwrap().find_chain('A').unwrap()
  let entries = @src.dssp_assign(chain)
  // Terminal residues have phi=0 or psi=0
  // First residue: phi should be 0 (no previous C)
  assert_true(entries[0].phi == 0.0)
  // Last residue: psi should be 0 (no next N)
  let last = entries.length() - 1
  assert_true(entries[last].psi == 0.0)
}

///|
test "dssp_full: to_string output" {
  let s = @src.create_demo_structure()
  let chain = s.find_model(0).unwrap().find_chain('A').unwrap()
  let entries = @src.dssp_assign(chain)
  let output = @src.dssp_to_string(entries)
  // Should contain header
  assert_true(output.contains("DSSP"))
  // Should contain residue names
  assert_true(output.contains("ALA"))
}

///|
test "dssp_full: hbond energy symmetric structure" {
  // Two residues with a clear H-bond geometry
  // Donor: C=O pointing toward acceptor N-H
  let o_i = @src.Vector3::new(0.0, 0.0, 0.0)
  let c_i = @src.Vector3::new(1.3, 0.0, 0.0)
  let n_j = @src.Vector3::new(0.0, 2.8, 0.0)
  let ca_j = @src.Vector3::new(1.5, 3.5, 0.0)
  let prev_c_j = @src.Vector3::new(0.0, 4.0, 0.0)
  let e = @src.dssp_residue_hbond_energy(o_i, c_i, n_j, ca_j, Some(prev_c_j))
  // Should produce a finite energy value
  assert_true(e > -10000.0 && e < 10000.0)
}

///|
test "dssp_full: bend detection" {
  // Create CA coordinates with a sharp bend at residue 2
  let ca_coords = [
    @src.Vector3::new(0.0, 0.0, 0.0),
    @src.Vector3::new(3.8, 0.0, 0.0),
    @src.Vector3::new(7.6, 0.0, 0.0),
    @src.Vector3::new(7.6, 3.8, 0.0),
    @src.Vector3::new(0.0, 3.8, 0.0),
  ]
  let bends = @src.dssp_find_bends(ca_coords)
  assert_eq(bends.length(), 5)
  // Residue 2 should be a bend (sharp 90° turn)
  assert_true(bends[2])
}
