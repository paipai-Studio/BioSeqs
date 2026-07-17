///|
/// Tests for Smith-Waterman local sequence alignment.

///|
test "sw_simple_match" {
  // Exact match: "ACGT" vs "ACGT" → score 8 (4×2)
  let result = @src.smith_waterman_dna("ACGT", "ACGT")
  assert_eq(result.score(), 8.0)
  assert_eq(result.aligned1(), "ACGT")
  assert_eq(result.aligned2(), "ACGT")
  assert_eq(result.match_line(), "||||")
  assert_eq(result.alignment_length(), 4)
  assert_eq(result.num_matches(), 4)
  assert_eq(result.percent_identity(), 100.0)
}

///|
test "sw_substring_match" {
  // seq2 is a substring of seq1
  let result = @src.smith_waterman_dna("GGTACGTAA", "ACGT")
  assert_eq(result.score(), 8.0)
  assert_eq(result.aligned1(), "ACGT")
  assert_eq(result.aligned2(), "ACGT")
  // seq1: "ACGT" starts at position 2, ends at 6 (0-based, end exclusive)
  let (s1, e1) = result.pos1()
  assert_eq(e1 - s1, 4) // alignment spans 4 characters in seq1
  // seq2: positions 0..4
  let (s2, e2) = result.pos2()
  assert_eq(s2, 0)
  assert_eq(e2, 4)
}

///|
test "sw_no_match" {
  // Completely different sequences → score 0, empty alignment
  let result = @src.smith_waterman_dna("AAAA", "TTTT")
  assert_eq(result.score(), 0.0)
  assert_eq(result.alignment_length(), 0)
}

///|
test "sw_with_gaps" {
  // "ACGT" vs "ACT" → aligns "ACGT" / "AC-T", score 5 (2+2-1+2)
  let result = @src.smith_waterman_dna("ACGT", "ACT")
  assert_eq(result.score(), 5.0)
  assert_eq(result.alignment_length(), 4)
  assert_eq(result.aligned1(), "ACGT")
  assert_eq(result.aligned2(), "AC-T")
}

///|
test "sw_partial_match" {
  // Two sequences sharing a common subsequence
  let result = @src.smith_waterman_dna("AACGGAATT", "AACGTT")
  // Best local alignment: "AACG" score 8, then mismatch/gaps reduce
  assert_true(result.score() > 0.0)
  assert_true(result.alignment_length() > 0)
  assert_true(result.num_matches() > 0)
}

///|
test "sw_protein_alignment" {
  let result = @src.smith_waterman_protein("HEAGAWGHEE", "PAWHEAE")
  // Known example from Durbin et al. textbook (with simplified scoring)
  // Should find "AWGHE" or "AWHE" region
  assert_true(result.score() > 0.0)
  assert_true(result.alignment_length() >= 3)
}

///|
test "sw_custom_scoring" {
  // Test with custom scoring: match=3, mismatch=-2, gap=-2
  let scoring = @src.Scoring::new(3.0, -2.0, -2.0)
  let result = @src.smith_waterman("ACGTACGT", "ACGT", scoring)
  assert_eq(result.score(), 12.0) // 4 matches × 3 = 12
  assert_eq(result.aligned1(), "ACGT")
  assert_eq(result.aligned2(), "ACGT")
}

///|
test "sw_empty_sequence" {
  let result = @src.smith_waterman_dna("", "ACGT")
  assert_eq(result.score(), 0.0)
  assert_eq(result.alignment_length(), 0)
}

///|
test "sw_single_char_match" {
  let result = @src.smith_waterman_dna("A", "A")
  assert_eq(result.score(), 2.0)
  assert_eq(result.aligned1(), "A")
  assert_eq(result.aligned2(), "A")
}

///|
test "sw_score_matrix_dimensions" {
  let result = @src.smith_waterman_dna("AC", "GT")
  let matrix = result.score_matrix()
  // (n+1) x (m+1) = 3 x 3
  assert_eq(matrix.length(), 3)
  assert_eq(matrix[0].length(), 3)
}

///|
test "sw_trace_matrix_dimensions" {
  let result = @src.smith_waterman_dna("AC", "GT")
  let trace = result.trace_matrix()
  assert_eq(trace.length(), 3)
  assert_eq(trace[0].length(), 3)
}

///|
test "sw_match_line_correctness" {
  // Test match line characters
  let result = @src.smith_waterman_dna("ACGT", "ACT")
  let line = result.match_line()
  // "ACGT" / "AC-T" → "|| |"
  assert_eq(line, "|| |")
}

///|
test "sw_percent_identity" {
  let result = @src.smith_waterman_dna("ACGT", "ACT")
  // 3 matches out of 4 positions = 75%
  assert_eq(result.percent_identity(), 75.0)
}

///|
test "sw_format_output" {
  let result = @src.smith_waterman_dna("ACGT", "ACGT")
  let formatted = result.format()
  // Should contain 3 lines separated by \n
  let lines = formatted.split("\n").to_array()
  assert_eq(lines.length(), 3)
  assert_eq(lines[0].to_string(), "ACGT")
  assert_eq(lines[1].to_string(), "||||")
  assert_eq(lines[2].to_string(), "ACGT")
}

///|
test "sw_long_sequence_alignment" {
  // Longer sequences with a shared motif
  let seq1 = "ATCGATCGATCGATCGATCGATCG"
  let seq2 = "GGGGATCGATCGGGGGATCGATCGGGGG"
  let result = @src.smith_waterman_dna(seq1, seq2)
  assert_true(result.score() > 0.0)
  assert_true(result.num_matches() >= 8)
}

///|
test "sw_repeated_motif" {
  // Find the best local alignment between sequences with repeats
  let result = @src.smith_waterman_dna("TACGTACGTACG", "ACGTACG")
  // Should align "ACGTACG" with "ACGTACG", score = 7×2 = 14
  assert_eq(result.score(), 14.0)
  assert_eq(result.num_matches(), 7)
}
