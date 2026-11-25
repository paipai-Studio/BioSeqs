///|
/// Tests for biotite_align_utils3.mbt — k-mer alphabets and gap utilities.
test "biotite3_kmer_alphabet_encode_decode" {
  let alphabet = @src.biotite3_kmer_alphabet(["A", "C"], 2)
  assert_eq(alphabet.kmer_symbols.length(), 4)
  assert_eq(@src.biotite3_kmer_alphabet_encode(alphabet, ["A", "C"]), 1)
  let decoded = @src.biotite3_kmer_alphabet_decode(alphabet, 1)
  assert_eq(decoded.length(), 2)
  assert_eq(decoded[0], "A")
  assert_eq(decoded[1], "C")
  assert_eq(@src.biotite3_kmer_alphabet_decode(alphabet, -1).length(), 0)
}

///|
test "biotite3_score_threshold_rule" {
  let rule = @src.biotite3_score_threshold_rule(
    ["A", "C"],
    [[0.0, 5.0], [5.0, 0.0]],
    3.0,
  )
  let similar = @src.biotite3_score_threshold_find_similar(rule, 0)
  assert_eq(similar.length(), 1)
  assert_eq(similar[0], 1)
}

///|
test "biotite3_mincode_select" {
  let positions = @src.biotite3_mincode_select(["A", "C", "G", "T"], 2, 2)
  assert_eq(positions.length(), 2)
  for p in positions {
    assert_true(p >= 0)
    assert_true(p <= 2)
  }
  assert_eq(@src.biotite3_mincode_select(["A"], 2, 1).length(), 0)
}

///|
test "biotite3_pairwise_identity_matrix" {
  let matrix = @src.biotite3_pairwise_identity_matrix(["ACGT", "ACGA", "TTTT"])
  assert_eq(matrix.length(), 3)
  assert_true(@moonbit.assert_num_close(matrix[0][0], 1.0, 1.0e-12))
  assert_true(@moonbit.assert_num_close(matrix[0][1], 0.75, 1.0e-12))
  assert_true(@moonbit.assert_num_close(matrix[0][2], 0.25, 1.0e-12))
}

///|
test "biotite3_terminal_gap_utilities" {
  let aligned = ["--AC--", "--AC--"]
  let (left, right) = @src.biotite3_find_terminal_gaps(aligned, "-")
  assert_eq(left, 2)
  assert_eq(right, 2)
  let trimmed = @src.biotite3_remove_terminal_gaps(aligned, "-")
  assert_eq(trimmed[0], "AC")
  assert_eq(trimmed[1], "AC")
  let ungapped = @src.biotite3_remove_all_gaps(["A-C", "ACG"], "-")
  assert_eq(ungapped[0], "AC")
  assert_eq(ungapped[1], "ACG")
}

///|
test "biotite3_alignment_score_and_symbols" {
  let score = @src.biotite3_alignment_score(
    ["AA", "AC"],
    [[2.0, -1.0], [-1.0, 2.0]],
    ["A", "C"],
  )
  assert_true(@moonbit.assert_num_close(score, 1.0, 1.0e-12))
  let codes = @src.biotite3_alignment_get_codes(["AC"], ["A", "C"])
  assert_eq(codes[0][0], 0)
  assert_eq(codes[0][1], 1)
  let symbols = @src.biotite3_alignment_get_symbols(["AC"])
  assert_eq(symbols[0][0], "A")
  assert_eq(symbols[0][1], "C")
}
