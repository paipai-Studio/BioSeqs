///|
test "pairwise_globalxx_identical" {
  let r = @src.pairwise_globalxx("ACGT", "ACGT", 2.0, -1.0, -5.0, -0.5)
  assert_eq(r.get_aligned_seq1(), "ACGT")
  assert_eq(r.get_aligned_seq2(), "ACGT")
  assert_true(r.get_score() > 0.0)
  assert_eq(r.get_mode(), "globalxx")
}

///|
test "pairwise_globalxx_with_gaps" {
  let r = @src.pairwise_globalxx("ACGT", "ACG", 2.0, -1.0, -5.0, -0.5)
  assert_eq(r.get_aligned_seq1().length(), r.get_aligned_seq2().length())
  assert_true(r.get_aligned_seq2().contains("-"))
}

///|
test "pairwise_globalxx_empty" {
  let r = @src.pairwise_globalxx("", "ACGT", 2.0, -1.0, -5.0, -0.5)
  assert_eq(r.get_aligned_seq1(), "")
  assert_eq(r.get_aligned_seq2(), "ACGT")
  assert_eq(r.get_score(), 0.0)
}

///|
test "pairwise_globalxx_both_empty" {
  let r = @src.pairwise_globalxx("", "", 2.0, -1.0, -5.0, -0.5)
  assert_eq(r.get_score(), 0.0)
}

///|
test "pairwise_localxx_identical" {
  let r = @src.pairwise_localxx("ACGTACGT", "ACGTACGT", 2.0, -1.0, -5.0, -0.5)
  assert_eq(r.get_aligned_seq1(), "ACGTACGT")
  assert_eq(r.get_aligned_seq2(), "ACGTACGT")
  assert_true(r.get_score() > 0.0)
}

///|
test "pairwise_localxx_substring" {
  let r = @src.pairwise_localxx("XXXACGTXXX", "ACGT", 2.0, -1.0, -5.0, -0.5)
  assert_eq(r.get_aligned_seq1(), "ACGT")
  assert_eq(r.get_aligned_seq2(), "ACGT")
}

///|
test "pairwise_localxx_no_match" {
  let r = @src.pairwise_localxx("AAAA", "TTTT", 2.0, -1.0, -5.0, -0.5)
  assert_eq(r.get_score(), 0.0)
}

///|
test "pairwise_localxx_empty" {
  let r = @src.pairwise_localxx("", "ACGT", 2.0, -1.0, -5.0, -0.5)
  assert_eq(r.get_score(), 0.0)
}

///|
test "pairwise_global_convenience" {
  let r = @src.pairwise_global("ACGT", "ACGT")
  assert_eq(r.get_aligned_seq1(), "ACGT")
  assert_eq(r.get_aligned_seq2(), "ACGT")
  assert_true(r.get_score() > 0.0)
}

///|
test "pairwise_local_convenience" {
  let r = @src.pairwise_local("ACGTACGT", "ACGT")
  assert_eq(r.get_aligned_seq2(), "ACGT")
  assert_true(r.get_score() > 0.0)
}

///|
test "simple_score_match" {
  let scorer = @src.simple_score(2.0, -1.0)
  assert_true((scorer("A", "A") - 2.0).abs() < 0.001)
  assert_true((scorer("A", "T") - (-1.0)).abs() < 0.001)
}

///|
test "simple_score_case_insensitive" {
  let scorer = @src.simple_score(2.0, -1.0)
  assert_true((scorer("A", "a") - 2.0).abs() < 0.001)
}

///|
test "identity_score" {
  let scorer = @src.identity_score()
  assert_true((scorer("A", "A") - 1.0).abs() < 0.001)
  assert_true((scorer("A", "T") - 0.0).abs() < 0.001)
}

///|
test "dna_matrix" {
  let m = @src.dna_matrix(2.0, -1.0)
  assert_true((m.get("AA").unwrap() - 2.0).abs() < 0.001)
  assert_true((m.get("AT").unwrap() - (-1.0)).abs() < 0.001)
  assert_eq(m.size(), 10) // 10 unique pairs
}

///|
test "matrix_score" {
  let m = @src.dna_matrix(2.0, -1.0)
  let scorer = @src.matrix_score(m, 0.0)
  assert_true((scorer("A", "A") - 2.0).abs() < 0.001)
  assert_true((scorer("A", "T") - (-1.0)).abs() < 0.001)
  assert_true((scorer("X", "X") - 0.0).abs() < 0.001) // default
}

///|
test "pairwise_globalms" {
  let m = @src.dna_matrix(2.0, -1.0)
  let r = @src.pairwise_globalms("ACGT", "ACGT", m, -5.0, -0.5)
  assert_eq(r.get_aligned_seq1(), "ACGT")
  assert_eq(r.get_aligned_seq2(), "ACGT")
  assert_true(r.get_score() > 0.0)
  assert_eq(r.get_mode(), "globalms")
}

///|
test "pairwise_localms" {
  let m = @src.dna_matrix(2.0, -1.0)
  let r = @src.pairwise_localms("XXXACGTXXX", "ACGT", m, -5.0, -0.5)
  assert_eq(r.get_aligned_seq1(), "ACGT")
  assert_eq(r.get_aligned_seq2(), "ACGT")
  assert_eq(r.get_mode(), "localms")
}

///|
test "pairwise_result_new" {
  let r = @src.pairwise_result_new("ACGT", "ACGT", 8.0, "test", 0, 4, 0, 4)
  assert_eq(r.get_aligned_seq1(), "ACGT")
  assert_eq(r.get_aligned_seq2(), "ACGT")
  assert_true((r.get_score() - 8.0).abs() < 0.001)
  assert_eq(r.get_mode(), "test")
}

///|
test "alignment_summary" {
  let r = @src.pairwise_global("ACGT", "ACGT")
  let s = @src.alignment_summary(r)
  assert_true(s.contains("Match"))
  assert_true(s.contains("Score"))
  assert_true(s.contains("Identity"))
}

///|
test "format_alignment" {
  let r = @src.pairwise_global("ACGT", "ACGT")
  let formatted = r.format_alignment(60)
  assert_true(formatted.contains("seq1"))
  assert_true(formatted.contains("seq2"))
  assert_true(formatted.contains("Score"))
}
