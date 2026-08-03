///|
/// Tests for Bio.CAPS cleaved amplified polymorphic sequences marker analysis.

// ---------------------------------------------------------------------------
// Data structure creation and accessors
// ---------------------------------------------------------------------------

test "caps_differential_cutsite_creation" {
  let dc = @src.CapsDifferentialCutsite::new(
    11, "EcoRI", [0], [1],
  )
  assert_eq(dc.start(), 11)
  assert_eq(dc.enzyme_name(), "EcoRI")
  assert_eq(dc.cuts_in().length(), 1)
  assert_eq(dc.cuts_in()[0], 0)
  assert_eq(dc.blocked_in().length(), 1)
  assert_eq(dc.blocked_in()[0], 1)
}

test "caps_differential_cutsite_to_string" {
  let dc = @src.CapsDifferentialCutsite::new(
    11, "EcoRI", [0], [1],
  )
  let s = dc.to_string()
  assert_true(s.contains("DifferentialCutsite"))
  assert_true(s.contains("pos=11"))
  assert_true(s.contains("enzyme=EcoRI"))
  assert_true(s.contains("cuts_in=[0]"))
  assert_true(s.contains("blocked_in=[1]"))
}

// ---------------------------------------------------------------------------
// CAPS map construction
// ---------------------------------------------------------------------------

test "caps_map_basic" {
  let seqs = @src.caps_sample_sequences()
  let enzymes = @src.caps_sample_enzymes()
  let m = @src.caps_map(seqs, enzymes)
  assert_eq(m.size(), 2)
  assert_eq(m.length(), 30)
  assert_true(m.enzymes().length() >= 1)
}

test "caps_map_named" {
  let seqs = @src.caps_sample_sequences()
  let enzymes = @src.caps_sample_enzymes()
  let m = @src.caps_map_named(seqs, ["strainA", "strainB"], enzymes)
  assert_eq(m.size(), 2)
  assert_eq(m.sequence_names()[0], "strainA")
  assert_eq(m.sequence_names()[1], "strainB")
}

test "caps_map_empty_sequences" {
  let enzymes = @src.caps_sample_enzymes()
  let m = @src.caps_map([], enzymes)
  assert_eq(m.size(), 0)
  assert_eq(m.length(), 0)
  assert_eq(m.dcut_count(), 0)
}

test "caps_map_single_sequence" {
  // Single sequence: no differential cutting possible
  let enzymes = @src.caps_sample_enzymes()
  let m = @src.caps_map(["AAAAAAAAGAGAATTCAAAAAAAAAAAAAA"], enzymes)
  assert_eq(m.size(), 1)
  assert_eq(m.dcut_count(), 0)
}

// ---------------------------------------------------------------------------
// Differential cutsite detection
// ---------------------------------------------------------------------------

test "caps_detects_ecori_differential" {
  let seqs = @src.caps_sample_sequences()
  let enzymes = @src.caps_sample_enzymes()
  let m = @src.caps_map(seqs, enzymes)
  // EcoRI (GAATTC) should cut seq0 but not seq1
  // seq0: ...GAATTC... at position 10, cut at 10+1=11
  // seq1: ...GAATTT... no cut
  let ecori_dcuts = m.get_dcuts_by_enzyme("EcoRI")
  assert_true(ecori_dcuts.length() >= 1)
  let dc = ecori_dcuts[0]
  assert_eq(dc.enzyme_name(), "EcoRI")
  assert_eq(dc.cuts_in()[0], 0)
  assert_eq(dc.blocked_in()[0], 1)
}

test "caps_no_differential_for_identical_sequences" {
  // Identical sequences: no differential cutting
  let enzymes = @src.caps_sample_enzymes()
  let m = @src.caps_map(
    ["AAAAAAAAGAGAATTCAAAAAAAAAAAAAA", "AAAAAAAAGAGAATTCAAAAAAAAAAAAAA"],
    enzymes,
  )
  // All sequences have the same cut sites, so no differential cutsites
  assert_eq(m.dcut_count(), 0)
}

test "caps_no_differential_without_enzyme_site" {
  // No enzyme recognition site in either sequence
  let enzymes = @src.caps_sample_enzymes()
  let m = @src.caps_map(
    ["ACGTACGTACGTACGTACGTACGTACGTAC", "ACGTACGTACGTACGTACGTACGTACGTAC"],
    enzymes,
  )
  assert_eq(m.dcut_count(), 0)
}

test "caps_multiple_differential_cutsites" {
  // Two sequences with multiple differential sites
  // seq0: GAATTC...AAGCTT (EcoRI + HindIII sites)
  // seq1: GAATTT...AAGCTT (only HindIII site)
  let enzymes = @src.caps_sample_enzymes()
  let m = @src.caps_map(
    [
      "AAAAAAAAGAGAATTCAAAAAGCTTAAAAAA",
      "AAAAAAAAGAGAATTTAAAAAGCTTAAAAAA",
    ],
    enzymes,
  )
  // EcoRI should produce a differential cutsite
  let ecori = m.get_dcuts_by_enzyme("EcoRI")
  assert_true(ecori.length() >= 1)
  // HindIII should NOT produce a differential cutsite (both sequences have it)
  let hindiii = m.get_dcuts_by_enzyme("HindIII")
  assert_eq(hindiii.length(), 0)
}

// ---------------------------------------------------------------------------
// Query methods
// ---------------------------------------------------------------------------

test "caps_get_dcuts_at_position" {
  let seqs = @src.caps_sample_sequences()
  let enzymes = @src.caps_sample_enzymes()
  let m = @src.caps_map(seqs, enzymes)
  // EcoRI cuts at position 11 (recognition at 10, cut at 10+1=11)
  let dcuts = m.get_dcuts_at_position(11)
  assert_true(dcuts.length() >= 1)
}

test "caps_has_dcuts_for_enzyme" {
  let seqs = @src.caps_sample_sequences()
  let enzymes = @src.caps_sample_enzymes()
  let m = @src.caps_map(seqs, enzymes)
  assert_true(m.has_dcuts_for_enzyme("EcoRI"))
  // HindIII has no site in the sample sequences
  assert_true(!m.has_dcuts_for_enzyme("HindIII"))
}

test "caps_dcut_count" {
  let seqs = @src.caps_sample_sequences()
  let enzymes = @src.caps_sample_enzymes()
  let m = @src.caps_map(seqs, enzymes)
  // At least one differential cutsite from EcoRI
  assert_true(m.dcut_count() >= 1)
}

// ---------------------------------------------------------------------------
// Formatting
// ---------------------------------------------------------------------------

test "caps_map_to_string" {
  let seqs = @src.caps_sample_sequences()
  let enzymes = @src.caps_sample_enzymes()
  let m = @src.caps_map(seqs, enzymes)
  let s = m.to_string()
  assert_true(s.contains("CAPSMap"))
  assert_true(s.contains("sequences=2"))
  assert_true(s.contains("length=30"))
}

test "caps_map_report" {
  let seqs = @src.caps_sample_sequences()
  let enzymes = @src.caps_sample_enzymes()
  let m = @src.caps_map(seqs, enzymes)
  let r = m.report()
  assert_true(r.contains("CAPSMap"))
  assert_true(r.contains("Differential cutsites:"))
  assert_true(r.contains("EcoRI"))
}

test "caps_map_report_no_dcuts" {
  let enzymes = @src.caps_sample_enzymes()
  let m = @src.caps_map(
    ["ACGTACGTACGTACGTACGTACGTACGTAC", "ACGTACGTACGTACGTACGTACGTACGTAC"],
    enzymes,
  )
  let r = m.report()
  assert_true(r.contains("No differential cutsites found."))
}

// ---------------------------------------------------------------------------
// Sample data
// ---------------------------------------------------------------------------

test "caps_sample_sequences_length" {
  let seqs = @src.caps_sample_sequences()
  assert_eq(seqs.length(), 2)
  // Both sequences should be the same length
  assert_eq(seqs[0].length(), seqs[1].length())
  assert_eq(seqs[0].length(), 30)
}

test "caps_sample_sequences_ecori_site" {
  let seqs = @src.caps_sample_sequences()
  // seq0 should contain GAATTC (EcoRI site)
  assert_true(seqs[0].contains("GAATTC"))
  // seq1 should contain GAATTT (no EcoRI site)
  assert_true(seqs[1].contains("GAATTT"))
  assert_true(!seqs[1].contains("GAATTC"))
}

test "caps_sample_enzymes" {
  let enzymes = @src.caps_sample_enzymes()
  // Should contain EcoRI and HindIII
  assert_true(enzymes.length() >= 2)
  let mut has_ecori = false
  let mut has_hindiii = false
  for e in enzymes {
    if e.name == "EcoRI" { has_ecori = true }
    if e.name == "HindIII" { has_hindiii = true }
  }
  assert_true(has_ecori)
  assert_true(has_hindiii)
}

// ---------------------------------------------------------------------------
// Edge cases
// ---------------------------------------------------------------------------

test "caps_three_sequences_differential" {
  // Three sequences: seq0 cut, seq1 not cut, seq2 cut
  let enzymes = @src.caps_sample_enzymes()
  let m = @src.caps_map(
    [
      "AAAAAAAAGAGAATTCAAAAAAAAAAAAAA", // has GAATTC
      "AAAAAAAAGAGAATTTAAAAAAAAAAAAAA", // has GAATTT (no cut)
      "AAAAAAAAGAGAATTCAAAAAAAAAAAAAA", // has GAATTC
    ],
    enzymes,
  )
  let ecori = m.get_dcuts_by_enzyme("EcoRI")
  assert_true(ecori.length() >= 1)
  let dc = ecori[0]
  // cuts_in should contain indices 0 and 2; blocked_in should contain 1
  assert_eq(dc.cuts_in().length(), 2)
  assert_eq(dc.blocked_in().length(), 1)
  assert_eq(dc.blocked_in()[0], 1)
}

test "caps_no_enzymes" {
  let seqs = @src.caps_sample_sequences()
  let m = @src.caps_map(seqs, [])
  assert_eq(m.dcut_count(), 0)
  assert_eq(m.enzymes().length(), 0)
}
