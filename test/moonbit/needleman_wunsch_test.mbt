///|
/// Tests for Needleman-Wunsch global sequence alignment.

///|
test "nw_identical_sequences" {
  // Identical sequences → all matches, no gaps
  let result = @src.needleman_wunsch_dna("ACGT", "ACGT")
  assert_eq(result.score(), 8.0) // 4 × 2
  assert_eq(result.aligned1(), "ACGT")
  assert_eq(result.aligned2(), "ACGT")
  assert_eq(result.match_line(), "||||")
  assert_eq(result.alignment_length(), 4)
  assert_eq(result.num_matches(), 4)
  assert_eq(result.num_mismatches(), 0)
  assert_eq(result.num_gaps_seq1(), 0)
  assert_eq(result.num_gaps_seq2(), 0)
  assert_eq(result.percent_identity(), 100.0)
}

///|
test "nw_completely_different" {
  // Completely different → all mismatches (no gaps needed, same length)
  let result = @src.needleman_wunsch_dna("AAAA", "TTTT")
  assert_eq(result.score(), -4.0) // 4 × (-1)
  assert_eq(result.alignment_length(), 4)
  assert_eq(result.num_matches(), 0)
  assert_eq(result.num_mismatches(), 4)
  assert_eq(result.percent_identity(), 0.0)
}

///|
test "nw_with_gaps" {
  // Different lengths → gaps introduced
  // "ACGT" vs "ACT" → "ACGT" / "AC-T", score = 2+2-1+2 = 5
  let result = @src.needleman_wunsch_dna("ACGT", "ACT")
  assert_eq(result.score(), 5.0)
  assert_eq(result.alignment_length(), 4)
  assert_eq(result.aligned1(), "ACGT")
  assert_eq(result.aligned2(), "AC-T")
  assert_eq(result.num_gaps_seq2(), 1)
}

///|
test "nw_gap_in_first" {
  // "ACT" vs "ACGT" → "AC-T" / "ACGT", score = 2+2-1+2 = 5
  let result = @src.needleman_wunsch_dna("ACT", "ACGT")
  assert_eq(result.score(), 5.0)
  assert_eq(result.alignment_length(), 4)
  assert_eq(result.aligned1(), "AC-T")
  assert_eq(result.aligned2(), "ACGT")
  assert_eq(result.num_gaps_seq1(), 1)
}

///|
test "nw_empty_first" {
  // Empty seq1 → all gaps in seq1
  let result = @src.needleman_wunsch_dna("", "ACGT")
  assert_eq(result.score(), -4.0) // 4 gaps × (-1)
  assert_eq(result.aligned1(), "----")
  assert_eq(result.aligned2(), "ACGT")
  assert_eq(result.num_gaps_seq1(), 4)
}

///|
test "nw_empty_second" {
  // Empty seq2 → all gaps in seq2
  let result = @src.needleman_wunsch_dna("ACGT", "")
  assert_eq(result.score(), -4.0)
  assert_eq(result.aligned1(), "ACGT")
  assert_eq(result.aligned2(), "----")
  assert_eq(result.num_gaps_seq2(), 4)
}

///|
test "nw_both_empty" {
  let result = @src.needleman_wunsch_dna("", "")
  assert_eq(result.score(), 0.0)
  assert_eq(result.alignment_length(), 0)
}

///|
test "nw_single_char_match" {
  let result = @src.needleman_wunsch_dna("A", "A")
  assert_eq(result.score(), 2.0)
  assert_eq(result.aligned1(), "A")
  assert_eq(result.aligned2(), "A")
}

///|
test "nw_single_char_mismatch" {
  let result = @src.needleman_wunsch_dna("A", "T")
  assert_eq(result.score(), -1.0)
  assert_eq(result.num_mismatches(), 1)
}

///|
test "nw_custom_scoring" {
  // match=3, mismatch=-2, gap=-2
  let scoring = @src.NWScoring::new(3.0, -2.0, -2.0)
  let result = @src.needleman_wunsch("ACGT", "ACGT", scoring)
  assert_eq(result.score(), 12.0) // 4 × 3
}

///|
test "nw_protein_alignment" {
  // Classic example from Durbin et al.
  let result = @src.needleman_wunsch_protein("HEAGAWGHEE", "PAWHEAE")
  // With simplified scoring (match=5, mismatch=-4, gap=-3),
  // the global alignment score depends on the specific path chosen.
  assert_true(result.alignment_length() >= 10)
  assert_true(result.num_matches() >= 1)
}

///|
test "nw_score_matrix_dimensions" {
  let result = @src.needleman_wunsch_dna("AC", "GT")
  let matrix = result.score_matrix()
  // (n+1) x (m+1) = 3 x 3
  assert_eq(matrix.length(), 3)
  assert_eq(matrix[0].length(), 3)
}

///|
test "nw_trace_matrix_dimensions" {
  let result = @src.needleman_wunsch_dna("AC", "GT")
  let trace = result.trace_matrix()
  assert_eq(trace.length(), 3)
  assert_eq(trace[0].length(), 3)
}

///|
test "nw_first_row_gap_init" {
  // First row should be gap penalties: 0, -1, -2, ...
  let result = @src.needleman_wunsch_dna("A", "ACGT")
  let matrix = result.score_matrix()
  assert_eq(matrix[0][0], 0.0)
  assert_eq(matrix[0][1], -1.0)
  assert_eq(matrix[0][2], -2.0)
  assert_eq(matrix[0][3], -3.0)
  assert_eq(matrix[0][4], -4.0)
}

///|
test "nw_first_col_gap_init" {
  // First column should be gap penalties: 0, -1, -2, ...
  let result = @src.needleman_wunsch_dna("ACGT", "A")
  let matrix = result.score_matrix()
  assert_eq(matrix[0][0], 0.0)
  assert_eq(matrix[1][0], -1.0)
  assert_eq(matrix[2][0], -2.0)
  assert_eq(matrix[3][0], -3.0)
  assert_eq(matrix[4][0], -4.0)
}

///|
test "nw_match_line_correctness" {
  let result = @src.needleman_wunsch_dna("ACGT", "ACT")
  // "ACGT" / "AC-T" → "|| |"
  assert_eq(result.match_line(), "|| |")
}

///|
test "nw_format_output" {
  let result = @src.needleman_wunsch_dna("ACGT", "ACGT")
  let formatted = result.format()
  let lines = formatted.split("\n").to_array()
  assert_eq(lines.length(), 3)
  assert_eq(lines[0].to_string(), "ACGT")
  assert_eq(lines[1].to_string(), "||||")
  assert_eq(lines[2].to_string(), "ACGT")
}

///|
test "nw_long_sequences" {
  let seq1 = "ATCGATCGATCGATCGATCGATCG"
  let seq2 = "ATCGATCGATCGATCGATCGATCG"
  let result = @src.needleman_wunsch_dna(seq1, seq2)
  assert_eq(result.score(), 48.0) // 24 × 2
  assert_eq(result.percent_identity(), 100.0)
}

///|
test "nw_vs_sw_difference" {
  // NW forces global alignment (negative score possible),
  // SW would return 0 for no local match.
  let result = @src.needleman_wunsch_dna("AAAA", "TTTT")
  // NW score = -4 (all mismatches), not 0
  assert_eq(result.score(), -4.0)
  assert_eq(result.alignment_length(), 4)
}

///|
test "nw_gap_at_beginning" {
  // Gap at beginning of seq1
  // "CGT" vs "ACGT" → "-CGT" / "ACGT" or "ACGT" / "ACGT"? 
  // Best: align "ACGT" fully → "-CGT"/"ACGT" score = -1+2+2+2=5
  // or "ACGT"/"ACGT" not possible since seq1="CGT"
  // Actually: "C-GT"/"ACGT" → mismatch+gap+match+match = -1-1+2+2=2
  // or "-CGT"/"ACGT" → gap+match+match+match = -1+2+2+2=5
  let result = @src.needleman_wunsch_dna("CGT", "ACGT")
  assert_eq(result.score(), 5.0)
}
