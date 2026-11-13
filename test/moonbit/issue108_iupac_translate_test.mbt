// Regression test for GitHub Issue #108.
// Seq.translate must resolve IUPAC-degenerate codons whose every expansion
// encodes the same amino acid (e.g. GCN → A, ATH → I), while still emitting
// 'X' when expansions disagree (e.g. TAN, NNN).

test "issue108_GCN_to_A" {
  assert_eq(@src.Seq::new("GCN").translate().to_string(), "A")
}

test "issue108_ATH_to_I" {
  assert_eq(@src.Seq::new("ATH").translate().to_string(), "I")
}

test "issue108_TAN_still_X" {
  // TAN = TAA(*) / TAG(*) / TAC(Y) / TAT(Y) → not unanimous
  assert_eq(@src.Seq::new("TAN").translate().to_string(), "X")
}

test "issue108_NNN_still_X" {
  assert_eq(@src.Seq::new("NNN").translate().to_string(), "X")
}

test "issue108_consensus_L_YTR" {
  // YTR: CTA, CTG, TTA, TTG → all Leucine
  assert_eq(@src.Seq::new("YTR").translate().to_string(), "L")
}

test "issue108_consensus_V_GTK" {
  // GTK: GTT, GTC, GTG → all Valine
  assert_eq(@src.Seq::new("GTK").translate().to_string(), "V")
}

test "issue108_consensus_stop_early_phase" {
  // TAN contains two stops (TAA, TAG) and two Tyr (TAT, TAC) → not unanimous → X
  assert_eq(@src.Seq::new("TAN").translate().to_string(), "X")
}

test "issue108_cds_standard_table_still_works" {
  // Sanity: exact standard codons untouched
  assert_eq(@src.Seq::new("ATGGCTGAC").translate().to_string(), "MAD")
}

test "issue108_cds_with_to_stop" {
  assert_eq(@src.Seq::new("ATGGCTTAAGCT").translate(to_stop=true).to_string(), "MA")
}
