// Regression test for GitHub Issue #108.
// Seq.translate must resolve IUPAC-degenerate codons whose every expansion
// encodes the same amino acid (e.g. GCN → A, ATH → I), while still emitting
// 'X' when expansions disagree (e.g. TAN, NNN).

///|
test "issue108_GCN_to_A" {
  assert_eq(@src.Seq::new("GCN").translate().to_string(), "A")
}

///|
test "issue108_ATH_to_I" {
  assert_eq(@src.Seq::new("ATH").translate().to_string(), "I")
}

///|
test "issue108_TAN_still_X" {
  // TAN = TAA(*) / TAG(*) / TAC(Y) / TAT(Y) → not unanimous
  assert_eq(@src.Seq::new("TAN").translate().to_string(), "X")
}

///|
test "issue108_NNN_still_X" {
  assert_eq(@src.Seq::new("NNN").translate().to_string(), "X")
}

///|
test "issue108_consensus_L_YTR" {
  // YTR: CTA, CTG, TTA, TTG → all Leucine
  assert_eq(@src.Seq::new("YTR").translate().to_string(), "L")
}

///|
test "issue108_consensus_V_GTK" {
  // GTK: GTT, GTC, GTG → all Valine
  assert_eq(@src.Seq::new("GTK").translate().to_string(), "V")
}

///|
test "issue108_consensus_stop_early_phase" {
  // TAN contains two stops (TAA, TAG) and two Tyr (TAT, TAC) → not unanimous → X
  assert_eq(@src.Seq::new("TAN").translate().to_string(), "X")
}

///|
test "issue108_cds_standard_table_still_works" {
  // Sanity: exact standard codons untouched
  assert_eq(@src.Seq::new("ATGGCTGAC").translate().to_string(), "MAD")
}

///|
test "issue108_cds_with_to_stop" {
  assert_eq(
    @src.Seq::new("ATGGCTTAAGCT").translate(to_stop=true).to_string(),
    "MA",
  )
}

///|
test "issue108_lowercase_resolves" {
  assert_eq(@src.Seq::new("gcn").translate().to_string(), "A")
}

///|
test "issue108_non_iupac_char_yields_x" {
  // GCX is a single codon containing a non-IUPAC character: whole codon → X
  // (unchanged from pre-fix behavior).
  assert_eq(@src.Seq::new("GCX").translate().to_string(), "X")
}

///|
test "issue108_resolved_stop_with_to_stop" {
  // TAR → TAA/TAG both stop → resolved stop; to_stop stops before it.
  assert_eq(@src.Seq::new("GCNTAR").translate(to_stop=true).to_string(), "A")
}

///|
test "issue108_resolved_stop_default_symbol" {
  assert_eq(@src.Seq::new("GCNTAR").translate().to_string(), "A*")
}

///|
test "issue108_resolved_stop_in_cds_raises" {
  let raised = try {
    ignore(@src.Seq::new("ATGCCNTARAAA").translate(cds=true))
    false
  } catch {
    SeqError(_) => true
  }
  assert_true(raised)
}

///|
test "issue108_invalid_char_still_raises_in_cds" {
  // An invalid (non-nucleotide, non-IUPAC) char must still raise, not emit X.
  let raised = try {
    ignore(@src.Seq::new("ATG1CT").translate())
    false
  } catch {
    SeqError(_) => true
  }
  assert_true(raised)
}

// ============================================================
// Probe group: systematic consensus enumeration. Each probe pins
// one IUPAC combination to its biologically correct translation,
// guarding the expand-and-compare logic in resolve_degenerate_codon.
// ============================================================

///|
test "issue108_probe_agr_to_R" {
  // AGR = AGA|AGG → both Arginine.
  assert_eq(@src.Seq::new("AGR").translate().to_string(), "R")
}

///|
test "issue108_probe_agy_to_S" {
  // AGY = AGT|AGC → both Serine.
  assert_eq(@src.Seq::new("AGY").translate().to_string(), "S")
}

///|
test "issue108_probe_gkn_is_X_but_gtk_is_V" {
  // Counter-example: GKN's 2nd position is K(G|T), so it spans
  // GGN(Gly) + GTN(Val) → mixed → X.
  // The true Valine consensus is GTK = GTG|GTT (fixed 2nd base T).
  assert_eq(@src.Seq::new("GKN").translate().to_string(), "X")
  assert_eq(@src.Seq::new("GTK").translate().to_string(), "V")
}

///|
test "issue108_probe_mgr_to_R" {
  // MGR spans 6 codons: AGA/AGG/CGA/CGC/CGG/CGT → all Arginine.
  // Exercises the 2×3 expansion path with 6 table lookups.
  assert_eq(@src.Seq::new("MGR").translate().to_string(), "R")
}

///|
test "issue108_probe_ytg_to_L" {
  // YTG = CTG|TTG → both Leucine (differ from YTR test: G fixed here).
  assert_eq(@src.Seq::new("YTG").translate().to_string(), "L")
}

///|
test "issue108_probe_yat_is_X" {
  // Counter-example: YAT = CAT(H)|TAT(Y) → mixed → X.
  assert_eq(@src.Seq::new("YAT").translate().to_string(), "X")
}

///|
test "issue108_probe_ath_atw_atd" {
  // Base-set sensitivity on the 3rd IUPAC code (2nd base fixed T):
  //   ATH = ATA|ATC|ATT (H excludes G) → all I
  //   ATW = ATA|ATT (W = A|T) → still all I
  //   ATD = ATA|ATC|ATG|ATT (D = A|G|T) → ATG=M sneaks in → X
  assert_eq(@src.Seq::new("ATH").translate().to_string(), "I")
  assert_eq(@src.Seq::new("ATW").translate().to_string(), "I")
  assert_eq(@src.Seq::new("ATD").translate().to_string(), "X")
}

///|
test "issue108_probe_tar_resolved_stop" {
  // TAR = TAA|TAG → both stop → single '*' symbol.
  assert_eq(@src.Seq::new("TAR").translate().to_string(), "*")
}

///|
test "issue108_probe_rna_gar" {
  // RNA alphabet: GAR = GAA|GAG → Glutamate (U and R handled as masks).
  assert_eq(@src.Seq::new("GAR").translate().to_string(), "E")
}

///|
test "issue108_probe_mixed_case_codon" {
  // Case-insensitivity must hold across the degenerate path too.
  assert_eq(@src.Seq::new("GcN").translate().to_string(), "A")
}

///|
test "issue108_probe_cds_with_degenerate_internal" {
  // CDS mode with a resolvable degenerate codon: ATG|GCN|CCC|TAA.
  // The internal GCN must translate to A instead of aborting the CDS.
  assert_eq(
    @src.Seq::new("ATGGCNCCCTAA").translate(cds=true).to_string(),
    "MAP",
  )
}

///|
test "issue108_probe_to_stop_not_stopped_by_X" {
  // TAN is NOT resolvable (X), so it must not terminate a to_stop scan;
  // only a real/resolved stop codon does.
  assert_eq(
    @src.Seq::new("GCNTANAAA").translate(to_stop=true).to_string(),
    "AXK",
  )
}
