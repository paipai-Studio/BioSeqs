// Tests for Bio.Align.analysis module

test "AlnAnalysisResult creation" {
  let result = @src.AlnAnalysisResult::new(0.1, 0.05, 2.0, 0.1, 10, 5)
  assert_true(result.aln_get_dn() == 0.1)
  assert_true(result.aln_get_ds() == 0.05)
  assert_true(result.aln_get_ratio() == 2.0)
  assert_true(result.aln_get_ratio_sem() == 0.1)
  assert_true(result.aln_get_n_syn() == 10)
  assert_true(result.aln_get_n_nonsyn() == 5)
}

test "positive selection detection" {
  let result = @src.AlnAnalysisResult::new(0.1, 0.05, 2.0, 0.1, 10, 5)
  assert_true(result.aln_has_positive_selection())
  assert_false(result.aln_has_purifying_selection())
}

test "purifying selection detection" {
  let result = @src.AlnAnalysisResult::new(0.05, 0.1, 0.5, 0.1, 5, 10)
  assert_false(result.aln_has_positive_selection())
  assert_true(result.aln_has_purifying_selection())
}

test "aln_aa_from_codon" {
  assert_true(@src.aln_aa_from_codon("ATG") == "M")
  assert_true(@src.aln_aa_from_codon("AAA") == "K")
  assert_true(@src.aln_aa_from_codon("CCC") == "P")
  assert_true(@src.aln_aa_from_codon("TAA") == "*")
}

test "aln_is_synonymous" {
  assert_true(@src.aln_is_synonymous("AAA", 2, "G")) // AAG also codes for K
  assert_false(@src.aln_is_synonymous("AAA", 1, "T")) // ATA codes for I
  assert_true(@src.aln_is_synonymous("TTT", 2, "C")) // TTC also codes for F
}

test "aln_analyze_dn_ds identical sequences" {
  let seq1 = "ATGAAACCC"
  let seq2 = "ATGAAACCC"
  let result = @src.aln_analyze_dn_ds(seq1, seq2)
  assert_true(result.aln_get_dn() == 0.0)
  assert_true(result.aln_get_ds() == 0.0)
}

test "aln_analyze_dn_ds different sequences" {
  let seq1 = "ATGAAACCC"
  let seq2 = "ATGAAGCCC"
  let result = @src.aln_analyze_dn_ds(seq1, seq2)
  assert_true(result.aln_get_ds() > 0.0)
}

test "aln_jukes_cantor_distance identical" {
  let dist = @src.aln_jukes_cantor_distance("ACGT", "ACGT")
  assert_true(dist == 0.0)
}

test "aln_jukes_cantor_distance different" {
  let dist = @src.aln_jukes_cantor_distance("ACGT", "TCGT")
  assert_true(dist > 0.0)
}

test "aln_kimura_2p_distance identical" {
  let dist = @src.aln_kimura_2p_distance("ACGT", "ACGT")
  assert_true(dist == 0.0)
}

test "aln_kimura_2p_distance transition" {
  let dist = @src.aln_kimura_2p_distance("ACGT", "TCGT")
  assert_true(dist > 0.0)
}

test "aln_p_distance identical" {
  let dist = @src.aln_p_distance("ACGT", "ACGT")
  assert_true(dist == 0.0)
}

test "aln_p_distance" {
  let dist = @src.aln_p_distance("ACGT", "AATT")
  assert_true(dist > 0.4 && dist < 0.6)
}

test "aln_jukes_cantor_matrix" {
  let sequences = ["ACGT", "TCGT", "ACGA"]
  let matrix = @src.aln_jukes_cantor_matrix(sequences)
  assert_true(matrix.length() == 3)
  assert_true(matrix[0][0] == 0.0)
  assert_true(matrix[1][1] == 0.0)
  assert_true(matrix[0][1] > 0.0)
}

test "aln_calculate_dn_ds_matrix" {
  let sequences = ["ATGAAACCC", "ATGAAGCCC", "ATGAAACCC"]
  let matrix = @src.aln_calculate_dn_ds_matrix(sequences)
  assert_true(matrix.length() == 3)
  assert_true(matrix[0][0].aln_get_ratio() == 0.0)
}
