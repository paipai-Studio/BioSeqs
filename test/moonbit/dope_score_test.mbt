///|
/// Tests for DOPE scoring.

///|
/// Build a test chain with specified residue types.
fn make_dope_chain(resnames : Array[String]) -> @src.Chain {
  let residues : Array[@src.Residue] = Array(capacity=0)
  for r = 0; r < resnames.length(); r = r + 1 {
    // Spread residues along x-axis with 3.8 Å spacing
    let x = r.to_double() * 3.8
    let ca = @src.Atom::new(
      name="CA",
      coord=@src.Vector3::new(x, 0.0, 0.0),
      resname=resnames[r],
      chainid='A',
      resseq=r + 1,
      element="C",
    )
    let res = @src.Residue::new(
      resname=resnames[r],
      chainid='A',
      resseq=r + 1,
      atoms=[ca],
    )
    residues.push(res)
  }
  @src.Chain::new(id='A', residues~)
}

///|
/// Build a compact chain (residues close together).
fn make_compact_chain(resnames : Array[String]) -> @src.Chain {
  let residues : Array[@src.Residue] = Array(capacity=0)
  for r = 0; r < resnames.length(); r = r + 1 {
    let angle = r.to_double() * 1.75
    let x = 3.0 * @math.cos(angle)
    let y = 3.0 * @math.sin(angle)
    let z = r.to_double() * 1.5
    let ca = @src.Atom::new(
      name="CA",
      coord=@src.Vector3::new(x, y, z),
      resname=resnames[r],
      chainid='A',
      resseq=r + 1,
      element="C",
    )
    let res = @src.Residue::new(
      resname=resnames[r],
      chainid='A',
      resseq=r + 1,
      atoms=[ca],
    )
    residues.push(res)
  }
  @src.Chain::new(id='A', residues~)
}

// ===== aa_index =====

test "dope_aa_index_standard" {
  assert_eq(@src.dope_aa_index("ALA"), 0)
  assert_eq(@src.dope_aa_index("VAL"), 19)
  assert_eq(@src.dope_aa_index("TRP"), 17)
  assert_eq(@src.dope_aa_index("CYS"), 4)
}

test "dope_aa_index_lowercase" {
  assert_eq(@src.dope_aa_index("ala"), 0)
  assert_eq(@src.dope_aa_index("cys"), 4)
}

test "dope_aa_index_nonstandard" {
  assert_eq(@src.dope_aa_index("XXX"), -1)
  assert_eq(@src.dope_aa_index("SEC"), -1)
}

// ===== pair_index =====

test "dope_pair_index_symmetric" {
  // pair_index(i, j) == pair_index(j, i)
  for i = 0; i < 20; i = i + 1 {
    for j = i; j < 20; j = j + 1 {
      assert_eq(@src.dope_pair_index(i, j), @src.dope_pair_index(j, i))
    }
  }
}

test "dope_pair_index_unique" {
  // All pairs should have unique indices
  let seen : Map[Int, Unit] = Map([])
  for i = 0; i < 20; i = i + 1 {
    for j = i; j < 20; j = j + 1 {
      let pidx = @src.dope_pair_index(i, j)
      assert_true(!seen.contains(pidx))
      seen[pidx] = ()
    }
  }
}

test "dope_pair_index_range" {
  // 210 unique pairs, indices 0..209
  for i = 0; i < 20; i = i + 1 {
    for j = i; j < 20; j = j + 1 {
      let pidx = @src.dope_pair_index(i, j)
      assert_true(pidx >= 0 && pidx < 210)
    }
  }
}

// ===== pair_compatibility =====

test "dope_pair_compat_hydrophobic" {
  // ILE-VAL: both hydrophobic → favorable
  let compat = @src.dope_pair_compatibility(9, 19) // ILE, VAL
  assert_true(compat > 1.0)
}

test "dope_pair_compat_same_charge" {
  // ARG-LYS: both positive → unfavorable
  let compat = @src.dope_pair_compatibility(1, 11) // ARG, LYS
  assert_true(compat < 0.0)
}

test "dope_pair_compat_opposite_charge" {
  // ARG-ASP: opposite charge → favorable
  let compat = @src.dope_pair_compatibility(1, 3) // ARG, ASP
  assert_true(compat > 1.0)
}

test "dope_pair_compat_cysteine" {
  // CYS-CYS: disulfide → very favorable
  let compat = @src.dope_pair_compatibility(4, 4) // CYS, CYS
  assert_true(compat > 2.0)
}

test "dope_pair_compat_aromatic" {
  // TRP-TYR: aromatic stacking → favorable
  let compat = @src.dope_pair_compatibility(17, 18) // TRP, TYR
  assert_true(compat > 1.0)
}

// ===== DopePotential::default =====

test "dope_potential_default" {
  let pot = @src.DopePotential::default()
  // Should have 210 pairs × 30 bins
  assert_eq(pot.energy.length(), 210)
  assert_eq(pot.energy[0].length(), 30)
  assert_true((pot.dist_cutoff - 15.0).abs() < 0.001)
  assert_true((pot.bin_size - 0.5).abs() < 0.001)
}

// ===== score_pair =====

test "dope_score_pair_steric" {
  let pot = @src.DopePotential::default()
  // Very short distance → steric clash (positive energy)
  let e = @src.DopePotential::score_pair(pot, "ALA", "VAL", 3.0)
  assert_true(e > 0.0)
}

test "dope_score_pair_long" {
  let pot = @src.DopePotential::default()
  // Long distance → neutral (0.0)
  let e = @src.DopePotential::score_pair(pot, "ALA", "VAL", 20.0)
  assert_true(e.abs() < 0.001)
}

test "dope_score_pair_nonstandard" {
  let pot = @src.DopePotential::default()
  let e = @src.DopePotential::score_pair(pot, "XXX", "YYY", 5.0)
  assert_true(e.abs() < 0.001)
}

test "dope_score_pair_hydrophobic_contact" {
  let pot = @src.DopePotential::default()
  // ILE-VAL at 5.0 Å: hydrophobic contact → negative (favorable)
  let e = @src.DopePotential::score_pair(pot, "ILE", "VAL", 5.0)
  assert_true(e < 0.0)
}

test "dope_score_pair_same_charge_close" {
  let pot = @src.DopePotential::default()
  // ARG-LYS at 4.0 Å: same charge close → positive (unfavorable)
  let e = @src.DopePotential::score_pair(pot, "ARG", "LYS", 4.0)
  assert_true(e > 0.0)
}

// ===== score (full chain) =====

test "dope_score_chain" {
  let pot = @src.DopePotential::default()
  let chain = make_compact_chain(["ALA", "VAL", "LEU", "ILE", "PHE", "TRP", "TYR", "MET"])
  match @src.DopePotential::score(pot, chain) {
    Some(result) => {
      assert_true(@src.DopeScore::get_num_pairs(result) > 0)
      // Per-residue should have 8 entries
      assert_eq(@src.DopeScore::get_per_residue(result).length(), 8)
    }
    None => assert_true(false)
  }
}

test "dope_score_compact_better" {
  // Compact chain (residues close) should score differently from extended
  let pot = @src.DopePotential::default()
  let resnames = ["ALA", "VAL", "LEU", "ILE", "PHE", "TRP", "TYR", "MET", "LEU", "VAL"]
  let extended = make_dope_chain(resnames)
  let compact = make_compact_chain(resnames)
  let s_ext = match @src.DopePotential::score(pot, extended) {
    Some(s) => @src.DopeScore::get_total(s)
    None => 0.0
  }
  let s_comp = match @src.DopePotential::score(pot, compact) {
    Some(s) => @src.DopeScore::get_total(s)
    None => 0.0
  }
  // They should be different
  assert_true((s_ext - s_comp).abs() > 0.01)
}

test "dope_score_normalized" {
  let pot = @src.DopePotential::default()
  let chain = make_dope_chain(["ALA", "VAL", "LEU", "ILE", "PHE", "TRP", "TYR", "MET"])
  match @src.DopePotential::score(pot, chain) {
    Some(result) => {
      let total = @src.DopeScore::get_total(result)
      let n = @src.DopeScore::get_num_pairs(result)
      let norm = @src.DopeScore::get_normalized(result)
      if n > 0 {
        assert_true((norm - total / n.to_double()).abs() < 0.001)
      }
    }
    None => assert_true(false)
  }
}

test "dope_score_too_short" {
  let pot = @src.DopePotential::default()
  let chain = make_dope_chain(["ALA"])
  match @src.DopePotential::score(pot, chain) {
    Some(_) => assert_true(false)
    None => assert_true(true)
  }
}

// ===== score_profile =====

test "dope_score_profile" {
  let pot = @src.DopePotential::default()
  let chain = make_dope_chain(["ALA", "VAL", "LEU", "ILE", "PHE", "TRP", "TYR", "MET"])
  match @src.DopePotential::score_profile(pot, chain) {
    Some(profile) => assert_eq(profile.length(), 8)
    None => assert_true(false)
  }
}

// ===== compare =====

test "dope_compare" {
  let pot = @src.DopePotential::default()
  let resnames = ["ALA", "VAL", "LEU", "ILE", "PHE", "TRP", "TYR", "MET", "LEU", "VAL"]
  let chain1 = make_dope_chain(resnames)
  let chain2 = make_compact_chain(resnames)
  match @src.DopePotential::compare(pot, chain1, chain2) {
    Some((s1, s2, diff)) => {
      assert_true((diff - (s1 - s2)).abs() < 0.001)
    }
    None => assert_true(false)
  }
}

// ===== rank_models =====

test "dope_rank_models" {
  let pot = @src.DopePotential::default()
  let resnames = ["ALA", "VAL", "LEU", "ILE", "PHE", "TRP", "TYR", "MET", "LEU", "VAL"]
  let chain1 = make_dope_chain(resnames)
  let chain2 = make_compact_chain(resnames)
  let chain3 = make_dope_chain(["GLY", "SER", "THR", "PRO", "GLY", "SER", "THR", "PRO", "GLY", "SER"])
  let ranked = @src.DopePotential::rank_models(pot, [chain1, chain2, chain3])
  assert_eq(ranked.length(), 3)
  // Should be sorted ascending (best first)
  assert_true(ranked[0].1 <= ranked[1].1)
  assert_true(ranked[1].1 <= ranked[2].1)
}

// ===== get_energy =====

test "dope_get_energy" {
  let pot = @src.DopePotential::default()
  let e = @src.DopePotential::get_energy(pot, 0, 0, 5) // ALA-ALA, bin 5
  // At bin 5, distance ~2.75 Å → steric clash
  assert_true(e > 0.0)
}

test "dope_get_energy_out_of_range" {
  let pot = @src.DopePotential::default()
  let e1 = @src.DopePotential::get_energy(pot, -1, 0, 5)
  assert_true(e1.abs() < 0.001)
  let e2 = @src.DopePotential::get_energy(pot, 0, 0, 100)
  assert_true(e2.abs() < 0.001)
}

// ===== worst/best residues =====

test "dope_worst_residues" {
  let pot = @src.DopePotential::default()
  let chain = make_dope_chain(["ALA", "VAL", "LEU", "ILE", "PHE", "TRP", "TYR", "MET"])
  match @src.DopePotential::score(pot, chain) {
    Some(result) => {
      let worst = @src.DopeScore::worst_residues(result, 3)
      assert_eq(worst.length(), 3)
    }
    None => assert_true(false)
  }
}

test "dope_best_residues" {
  let pot = @src.DopePotential::default()
  let chain = make_dope_chain(["ALA", "VAL", "LEU", "ILE", "PHE", "TRP", "TYR", "MET"])
  match @src.DopePotential::score(pot, chain) {
    Some(result) => {
      let best = @src.DopeScore::best_residues(result, 3)
      assert_eq(best.length(), 3)
    }
    None => assert_true(false)
  }
}

// ===== DopePotential::new with custom params =====

test "dope_potential_custom" {
  let pot = @src.DopePotential::new(dist_cutoff=10.0, bin_size=1.0, min_seq_sep=3)
  assert_true((pot.dist_cutoff - 10.0).abs() < 0.001)
  assert_true((pot.bin_size - 1.0).abs() < 0.001)
  assert_eq(pot.min_seq_sep, 3)
  // 10 bins
  assert_eq(pot.num_bins, 10)
}

test "dope_score_min_seq_sep" {
  // With min_seq_sep=3, more pairs should be counted
  let pot1 = @src.DopePotential::new(min_seq_sep=5)
  let pot2 = @src.DopePotential::new(min_seq_sep=3)
  let chain = make_dope_chain(["ALA", "VAL", "LEU", "ILE", "PHE", "TRP", "TYR", "MET"])
  let n1 = match @src.DopePotential::score(pot1, chain) {
    Some(s) => @src.DopeScore::get_num_pairs(s)
    None => 0
  }
  let n2 = match @src.DopePotential::score(pot2, chain) {
    Some(s) => @src.DopeScore::get_num_pairs(s)
    None => 0
  }
  // Lower min_seq_sep → more pairs
  assert_true(n2 >= n1)
}
