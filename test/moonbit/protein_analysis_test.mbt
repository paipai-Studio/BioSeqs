///|
/// Tests for the protein analysis module (protein_analysis.mbt).

// ============================================================================
// Hydropathy Tests
// ============================================================================

test "protein_hydropathy_basic" {
  let seq = "ACDEFGHIKLMNPQRSTVWY"
  let hydro = @src.protein_hydropathy(seq, window=7)
  assert_true(hydro.length() > 0)
  // All values should be finite
  let mut i = 0
  while i < hydro.length() {
    assert_true(hydro[i] > -10.0 && hydro[i] < 10.0)
    i = i + 1
  }
}

test "protein_hydropathy_short_sequence" {
  let seq = "ACG"
  let hydro = @src.protein_hydropathy(seq, window=7)
  assert_eq(hydro.length(), 0)
}

test "protein_hydropathy_custom_window" {
  let seq = "ACDEFGHIKLMNPQRSTVWY"
  let hydro5 = @src.protein_hydropathy(seq, window=5)
  let hydro7 = @src.protein_hydropathy(seq, window=7)
  assert_true(hydro5.length() > hydro7.length())
}

// ============================================================================
// Secondary Structure Prediction Tests
// ============================================================================

test "protein_secondary_structure_basic" {
  let seq = "ACDEFGHIKLMNPQRSTVWY"
  let ss = @src.protein_secondary_structure(seq)
  assert_eq(ss.length(), seq.length())
  // Each element should be H, E, T, or C
  let mut i = 0
  while i < ss.length() {
    let s = ss[i]
    assert_true(s == "H" || s == "E" || s == "T" || s == "C")
    i = i + 1
  }
}

test "protein_secondary_structure_all_chars" {
  // Test that all 20 amino acids are handled
  let seq = "ACDEFGHIKLMNPQRSTVWY"
  let ss = @src.protein_secondary_structure(seq)
  assert_eq(ss.length(), 20)
}

// ============================================================================
// Antigenicity Tests
// ============================================================================

test "protein_antigenicity_basic" {
  let seq = "ACDEFGHIKLMNPQRSTVWY"
  let antigen = @src.protein_antigenicity(seq)
  assert_eq(antigen.length(), seq.length())
  // Hopp-Woods scale ranges from about -0.07 to 0.37
  let mut i = 0
  while i < antigen.length() {
    assert_true(antigen[i] > -1.0 && antigen[i] < 1.0)
    i = i + 1
  }
}

test "protein_antigenicity_positive_residues" {
  // Known antigenic residues should have positive values
  let antigen = @src.protein_antigenicity("RKP")
  assert_true(antigen[0] > 0.0) // R = 0.37
  assert_true(antigen[1] > 0.0) // K = 0.16
  assert_true(antigen[2] > 0.0) // P = 0.18
}

// ============================================================================
// Transmembrane Prediction Tests
// ============================================================================

test "protein_transmembrane_basic" {
  let seq = "ACDEFGHIKLMNPQRSTVWY"
  let tm = @src.protein_transmembrane(seq, window=19)
  // Sequence is 20 chars, window is 19, should produce 2 values
  assert_eq(tm.length(), 2)
}

test "protein_transmembrane_short_sequence" {
  let seq = "ACG"
  let tm = @src.protein_transmembrane(seq, window=19)
  assert_eq(tm.length(), 0)
}

test "protein_transmembrane_propensity" {
  // Known transmembrane residues should have positive scores
  let v_tm = @src.protein_transmembrane_propensity('V')
  assert_true(v_tm > 0.0) // V = 1.19
  let i_tm = @src.protein_transmembrane_propensity('I')
  assert_true(i_tm > 0.0) // I = 1.0
}

// ============================================================================
// Amino Acid Composition Tests
// ============================================================================

test "protein_aa_composition_basic" {
  let seq = "AAL"
  let comp = @src.protein_aa_composition(seq)
  // A appears twice, L appears once
  let a_freq = comp['A']
  let l_freq = comp['L']
  assert_true(a_freq > 0.0)
  assert_true(l_freq > 0.0)
  // A should be twice as frequent as L
  assert_true(a_freq > l_freq)
}

test "protein_aa_composition_empty" {
  let seq = ""
  let comp = @src.protein_aa_composition(seq)
  // Empty sequence should have all zero frequencies
  assert_eq(comp.length(), 0)
}

// ============================================================================
// Dipeptide/Tripeptide Composition Tests
// ============================================================================

test "protein_dipeptide_composition_basic" {
  let seq = "ALA"
  let comp = @src.protein_dipeptide_composition(seq)
  assert_true(comp.length() > 0)
  let al_freq = comp["AL"]
  assert_true(al_freq > 0.0)
  let la_freq = comp["LA"]
  assert_true(la_freq > 0.0)
}

test "protein_dipeptide_composition_short" {
  let seq = "A"
  let comp = @src.protein_dipeptide_composition(seq)
  assert_eq(comp.length(), 0)
}

test "protein_tripeptide_composition_basic" {
  let seq = "ALA"
  let comp = @src.protein_tripeptide_composition(seq)
  assert_true(comp.length() > 0)
  let ala_freq = comp["ALA"]
  assert_true(ala_freq > 0.0)
}

test "protein_tripeptide_composition_short" {
  let seq = "AL"
  let comp = @src.protein_tripeptide_composition(seq)
  assert_eq(comp.length(), 0)
}

// ============================================================================
// Conservation Tests
// ============================================================================

test "protein_conservation_identical" {
  let alignment = ["AAAA", "AAAA", "AAAA"]
  let cons = @src.protein_conservation(alignment)
  assert_eq(cons.length(), 4)
  // All positions should have perfect conservation
  let mut i = 0
  while i < cons.length() {
    assert_true(cons[i] > 0.99)
    i = i + 1
  }
}

test "protein_conservation_different" {
  let alignment = ["AAAA", "CCCC", "GGGG"]
  let cons = @src.protein_conservation(alignment)
  assert_eq(cons.length(), 4)
  // All positions should have some conservation (different AAs at each position)
  let mut i = 0
  while i < cons.length() {
    assert_true(cons[i] < 0.7)
    i = i + 1
  }
}

test "protein_conservation_empty" {
  let alignment : Array[String] = []
  let cons = @src.protein_conservation(alignment)
  assert_eq(cons.length(), 0)
}

// ============================================================================
// Summary Report Tests
// ============================================================================

test "protein_summary_basic" {
  let seq = "ACDEFGHIKLMNPQRSTVWY"
  let summary = @src.protein_summary(seq)
  assert_true(summary.contains("Protein Sequence Analysis Report"))
  assert_true(summary.contains("Length: 20"))
  assert_true(summary.contains("Amino Acid Composition"))
  assert_true(summary.contains("Secondary Structure Prediction"))
}
