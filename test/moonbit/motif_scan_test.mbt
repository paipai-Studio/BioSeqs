///|
/// Tests for the FIMO-like motif scanning module (motif_scan.mbt).
///
/// Note: The PSSM struct is named `MotifScanPSSM` and the counts-to-PWM
/// builder is named `motif_scan_pwm_from_counts` to avoid collisions with the
/// pre-existing `PSSM` type alias and `pwm_from_counts` function in motifs.mbt.

///|
/// Helper: a simple 1-position PWM (A=0.9, C/G/T=0.05) for analytic p-value checks.
fn make_one_pos_pwm() -> @src.PositionWeightMatrix {
  @src.PositionWeightMatrix::new(
    motif_id="P1",
    motif_name="OnePos",
    matrix=[[0.9], [0.05], [0.05], [0.05]],
    background=[0.25, 0.25, 0.25, 0.25],
    strand="+",
  )
}

///|
test "pwm_construction_and_accessors" {
  let pwm = @src.PositionWeightMatrix::new(
    motif_id="M01",
    motif_name="MyMotif",
    matrix=[[0.9, 0.1], [0.05, 0.8], [0.03, 0.05], [0.02, 0.05]],
    background=[0.25, 0.25, 0.25, 0.25],
    strand="+-",
  )
  assert_eq(pwm.motif_id, "M01")
  assert_eq(pwm.motif_name, "MyMotif")
  assert_eq(pwm.matrix.length(), 4)
  assert_eq(pwm.matrix[0].length(), 2)
  assert_eq(pwm.background.length(), 4)
  assert_eq(pwm.strand, "+-")
}

///|
test "pwm_construction_default_background_and_strand" {
  let pwm = @src.PositionWeightMatrix::new(
    motif_id="M02",
    motif_name="DefaultMotif",
    matrix=[[0.25, 0.25], [0.25, 0.25], [0.25, 0.25], [0.25, 0.25]],
  )
  assert_eq(pwm.background[0], 0.25)
  assert_eq(pwm.background[3], 0.25)
  assert_eq(pwm.strand, "+-")
}

///|
test "pwm_from_counts_basic" {
  // 4 rows (A,C,G,T) x 3 columns
  let counts = [
    [9, 1, 1],
    [1, 8, 1],
    [1, 1, 7],
    [1, 2, 3],
  ]
  let pwm = @src.motif_scan_pwm_from_counts("MC1", "CountMotif", counts, 1.0)
  assert_eq(pwm.motif_id, "MC1")
  assert_eq(pwm.motif_name, "CountMotif")
  assert_eq(pwm.matrix.length(), 4)
  assert_eq(pwm.matrix[0].length(), 3)
  // Each column should sum to ~1.0 (with pseudocount).
  let col0 = pwm.matrix[0][0] + pwm.matrix[1][0] + pwm.matrix[2][0] + pwm.matrix[3][0]
  assert_true((col0 - 1.0).abs() < 0.001)
}

///|
test "pwm_from_counts_normalization" {
  // Counts: column 0 has total 10 + 4*pseudocount; each entry = (count+pseudo)/norm
  let counts = [
    [10, 0],
    [0, 10],
    [0, 0],
    [0, 0],
  ]
  let pwm = @src.motif_scan_pwm_from_counts("MC2", "Norm", counts, 1.0)
  // Column 0: total=10, norm=10+4=14. A = (10+1)/14 = 11/14
  let expected_a = 11.0 / 14.0
  assert_true((pwm.matrix[0][0] - expected_a).abs() < 0.001)
  // Column 1: total=10, norm=14. C = (10+1)/14 = 11/14
  let expected_c = 11.0 / 14.0
  assert_true((pwm.matrix[1][1] - expected_c).abs() < 0.001)
  // Zero-count entry: (0+1)/14
  let expected_zero = 1.0 / 14.0
  assert_true((pwm.matrix[1][0] - expected_zero).abs() < 0.001)
}

///|
test "pssm_construction_log_likelihood" {
  let pwm = make_one_pos_pwm()
  let pssm = @src.build_pssm(pwm)
  assert_eq(pssm.motif_id, "P1")
  assert_eq(pssm.scores.length(), 4)
  assert_eq(pssm.scores[0].length(), 1)
  // A: log2(0.9/0.25) = log2(3.6) ~ 1.848
  assert_true((pssm.scores[0][0] - 1.848).abs() < 0.01)
  // C: log2(0.05/0.25) = log2(0.2) ~ -2.322
  assert_true((pssm.scores[1][0] - (-2.322)).abs() < 0.01)
}

///|
test "pssm_construction_min_max_score" {
  let pwm = make_one_pos_pwm()
  let pssm = @src.build_pssm(pwm)
  // min_score = min of all bases at position 0 = -2.322
  assert_true((pssm.min_score - (-2.322)).abs() < 0.01)
  // max_score = max of all bases at position 0 = 1.848
  assert_true((pssm.max_score - 1.848).abs() < 0.01)
}

///|
test "pssm_carries_metadata" {
  let pwm = @src.PositionWeightMatrix::new(
    motif_id="M99",
    motif_name="Meta",
    matrix=[[0.9, 0.05, 0.05, 0.05, 0.9, 0.9],
           [0.05, 0.05, 0.05, 0.9, 0.05, 0.05],
           [0.05, 0.05, 0.9, 0.05, 0.05, 0.05],
           [0.05, 0.9, 0.05, 0.05, 0.05, 0.05]],
    background=[0.25, 0.25, 0.25, 0.25],
    strand="+-",
  )
  let pssm = @src.build_pssm(pwm)
  assert_eq(pssm.motif_id, "M99")
  assert_eq(pssm.motif_name, "Meta")
  assert_eq(pssm.background.length(), 4)
  assert_eq(pssm.strand, "+-")
}

///|
test "compute_score_single_position" {
  let pwm = make_one_pos_pwm()
  let pssm = @src.build_pssm(pwm)
  // Score of "A" = log2(0.9/0.25) ~ 1.848
  assert_true((@src.compute_score(pssm, "A") - 1.848).abs() < 0.01)
  // Score of "C" = log2(0.05/0.25) ~ -2.322
  assert_true((@src.compute_score(pssm, "C") - (-2.322)).abs() < 0.01)
}

///|
test "compute_score_perfect_match_equals_max" {
  let pwm = @src.motif_scan_sample_pwm()
  let pssm = @src.build_pssm(pwm)
  let score = @src.compute_score(pssm, "ATGCAA")
  assert_true((score - pssm.max_score).abs() < 0.001)
}

///|
test "compute_score_worst_match_equals_min" {
  let pwm = @src.motif_scan_sample_pwm()
  let pssm = @src.build_pssm(pwm)
  // CGTATG: each position picks a non-consensus (minimum-score) base.
  let score = @src.compute_score(pssm, "CGTATG")
  assert_true((score - pssm.min_score).abs() < 0.001)
}

///|
test "compute_score_lowercase_handled" {
  let pwm = @src.motif_scan_sample_pwm()
  let pssm = @src.build_pssm(pwm)
  let upper = @src.compute_score(pssm, "ATGCAA")
  let lower = @src.compute_score(pssm, "atgcaa")
  assert_true((upper - lower).abs() < 0.001)
}

///|
test "compute_pvalue_known_value_one_position" {
  // 1-position PWM: A scores ~1.848 (prob 0.25), C/G/T score ~-2.322 (prob 0.75).
  let pwm = make_one_pos_pwm()
  let pssm = @src.build_pssm(pwm)
  // P(score >= max_score) = P(A) = 0.25
  let pval_max = @src.compute_pvalue(pssm, pssm.max_score)
  assert_true((pval_max - 0.25).abs() < 0.02)
  // P(score >= 0) = P(A) = 0.25 (only A scores >= 0)
  let pval_zero = @src.compute_pvalue(pssm, 0.0)
  assert_true((pval_zero - 0.25).abs() < 0.02)
}

///|
test "compute_pvalue_min_score_is_one" {
  let pwm = make_one_pos_pwm()
  let pssm = @src.build_pssm(pwm)
  let pval = @src.compute_pvalue(pssm, pssm.min_score)
  assert_true((pval - 1.0).abs() < 0.001)
}

///|
test "compute_pvalue_above_max_is_zero" {
  let pwm = make_one_pos_pwm()
  let pssm = @src.build_pssm(pwm)
  let pval = @src.compute_pvalue(pssm, pssm.max_score + 1.0)
  assert_true((pval - 0.0).abs() < 0.001)
}

///|
test "compute_pvalue_monotonicity" {
  let pwm = @src.motif_scan_sample_pwm()
  let pssm = @src.build_pssm(pwm)
  // Higher score should yield lower (or equal) p-value.
  let pval_low = @src.compute_pvalue(pssm, pssm.min_score + 0.1)
  let pval_mid = @src.compute_pvalue(pssm, 0.0)
  let pval_high = @src.compute_pvalue(pssm, pssm.max_score)
  assert_true(pval_low >= pval_mid)
  assert_true(pval_mid >= pval_high)
}

///|
test "compute_pvalue_in_unit_range" {
  let pwm = @src.motif_scan_sample_pwm()
  let pssm = @src.build_pssm(pwm)
  for s in 0..<20 {
    let frac = s.to_double() / 20.0
    let score = pssm.min_score + frac * (pssm.max_score - pssm.min_score)
    let pval = @src.compute_pvalue(pssm, score)
    assert_true(pval >= 0.0)
    assert_true(pval <= 1.0)
  }
}

///|
test "compute_pvalue_perfect_match_small" {
  let pwm = @src.motif_scan_sample_pwm()
  let pssm = @src.build_pssm(pwm)
  let pval = @src.compute_pvalue(pssm, pssm.max_score)
  // Perfect match p-value ~ 0.25^6 = 0.000244
  assert_true(pval < 0.001)
  assert_true(pval > 0.0001)
}

///|
test "scan_sequence_forward_strand" {
  let pwm = @src.motif_scan_sample_pwm()
  let pssm = @src.build_pssm(pwm)
  // seq1 has ATGCAA at position 4 (forward).
  let matches = @src.scan_sequence(pssm, "GGGGATGCAAGGGG", "seq1", 0.001)
  let mut found_forward = false
  for i in 0..<matches.length() {
    if matches[i].strand == "+" && matches[i].start == 4 {
      assert_eq(matches[i].end, 10)
      assert_eq(matches[i].matched_sequence, "ATGCAA")
      assert_eq(matches[i].sequence_id, "seq1")
      found_forward = true
    }
  }
  assert_true(found_forward)
}

///|
test "scan_sequence_reverse_strand" {
  let pwm = @src.motif_scan_sample_pwm()
  let pssm = @src.build_pssm(pwm)
  // seq2 contains TTGCAT (revcomp of ATGCAA) at position 4 -> "-" match at 4.
  let matches = @src.scan_sequence(pssm, "GGGGTTGCATGGGG", "seq2", 0.001)
  let mut found_reverse = false
  for i in 0..<matches.length() {
    if matches[i].strand == "-" && matches[i].start == 4 {
      assert_eq(matches[i].end, 10)
      assert_eq(matches[i].matched_sequence, "ATGCAA")
      found_reverse = true
    }
  }
  assert_true(found_reverse)
}

///|
test "scan_sequence_both_strands" {
  let pwm = @src.motif_scan_sample_pwm()
  let pssm = @src.build_pssm(pwm)
  // Sequence with ATGCAA at pos 0 (forward) and TTGCAT at pos 10 (reverse).
  let matches = @src.scan_sequence(pssm, "ATGCAAGGGGTTGCAT", "both", 0.001)
  let mut found_plus = false
  let mut found_minus = false
  for i in 0..<matches.length() {
    if matches[i].strand == "+" && matches[i].start == 0 {
      found_plus = true
    }
    if matches[i].strand == "-" && matches[i].start == 10 {
      found_minus = true
    }
  }
  assert_true(found_plus)
  assert_true(found_minus)
}

///|
test "scan_sequence_forward_only_mode" {
  // A PWM with strand="+" should not produce reverse-strand matches.
  let pwm = @src.PositionWeightMatrix::new(
    motif_id="Mfwd",
    motif_name="FwdOnly",
    matrix=[[0.9, 0.05, 0.05, 0.05, 0.9, 0.9],
           [0.05, 0.05, 0.05, 0.9, 0.05, 0.05],
           [0.05, 0.05, 0.9, 0.05, 0.05, 0.05],
           [0.05, 0.9, 0.05, 0.05, 0.05, 0.05]],
    background=[0.25, 0.25, 0.25, 0.25],
    strand="+",
  )
  let pssm = @src.build_pssm(pwm)
  // seq2 only has a reverse-strand match; with strand="+" none should be found.
  let matches = @src.scan_sequence(pssm, "GGGGTTGCATGGGG", "seq2", 0.001)
  assert_eq(matches.length(), 0)
}

///|
test "scan_sequences_multiple" {
  let pwm = @src.motif_scan_sample_pwm()
  let pssm = @src.build_pssm(pwm)
  let seqs = @src.motif_scan_sample_sequences()
  let matches = @src.scan_sequences(pssm, seqs, 0.001)
  // seq1: 1 forward match; seq2: 1 reverse match; seq3: 2 forward matches.
  assert_eq(matches.length(), 4)
}

///|
test "motif_match_construction_and_accessors" {
  let m = @src.MotifMatch::new(
    motif_id="M1",
    motif_name="TF1",
    sequence_id="s1",
    start=5,
    end=11,
    strand="+",
    score=7.5,
    p_value=0.001,
    matched_sequence="ATGCAA",
  )
  assert_eq(m.motif_id, "M1")
  assert_eq(m.motif_name, "TF1")
  assert_eq(m.sequence_id, "s1")
  assert_eq(m.start, 5)
  assert_eq(m.end, 11)
  assert_eq(m.strand, "+")
  assert_true((m.score - 7.5).abs() < 0.001)
  assert_true((m.p_value - 0.001).abs() < 0.0001)
  assert_eq(m.matched_sequence, "ATGCAA")
}

///|
test "filter_matches_by_pvalue" {
  let m1 = @src.MotifMatch::new(
    motif_id="M",
    motif_name="N",
    sequence_id="s",
    start=0,
    end=4,
    strand="+",
    score=1.0,
    p_value=0.0001,
    matched_sequence="ACGT",
  )
  let m2 = @src.MotifMatch::new(
    motif_id="M",
    motif_name="N",
    sequence_id="s",
    start=1,
    end=5,
    strand="+",
    score=0.5,
    p_value=0.01,
    matched_sequence="CGTA",
  )
  let m3 = @src.MotifMatch::new(
    motif_id="M",
    motif_name="N",
    sequence_id="s",
    start=2,
    end=6,
    strand="-",
    score=2.0,
    p_value=0.0005,
    matched_sequence="TTAA",
  )
  let all = [m1, m2, m3]
  let filtered = @src.filter_matches(all, 0.001)
  assert_eq(filtered.length(), 2)
  // m1 (0.0001) and m3 (0.0005) pass; m2 (0.01) does not.
  assert_eq(filtered[0].start, 0)
  assert_eq(filtered[1].start, 2)
}

///|
test "motif_scan_summary_generation" {
  let m = @src.MotifMatch::new(
    motif_id="M1",
    motif_name="TF1",
    sequence_id="seq1",
    start=4,
    end=10,
    strand="+",
    score=11.09,
    p_value=0.000244,
    matched_sequence="ATGCAA",
  )
  let summary = @src.motif_scan_summary([m])
  assert_true(summary.contains("FIMO Motif Scan Summary"))
  assert_true(summary.contains("Total matches: 1"))
  assert_true(summary.contains("M1"))
  assert_true(summary.contains("seq1"))
  assert_true(summary.contains("ATGCAA"))
  assert_true(summary.contains("+"))
}

///|
test "motif_scan_summary_empty" {
  let summary = @src.motif_scan_summary([])
  assert_true(summary.contains("Total matches: 0"))
}

///|
test "sample_pwm_validation" {
  let pwm = @src.motif_scan_sample_pwm()
  assert_eq(pwm.motif_id, "M0001")
  assert_eq(pwm.motif_name, "SampleTF")
  assert_eq(pwm.matrix.length(), 4)
  assert_eq(pwm.matrix[0].length(), 6)
  // Consensus ATGCAA: pos0=A(0.9), pos1=T(row3=0.9), pos2=G(row2=0.9),
  // pos3=C(row1=0.9), pos4=A(0.9), pos5=A(0.9).
  assert_true((pwm.matrix[0][0] - 0.9).abs() < 0.001)
  assert_true((pwm.matrix[3][1] - 0.9).abs() < 0.001)
  assert_true((pwm.matrix[2][2] - 0.9).abs() < 0.001)
  assert_true((pwm.matrix[1][3] - 0.9).abs() < 0.001)
  assert_true((pwm.matrix[0][4] - 0.9).abs() < 0.001)
  assert_true((pwm.matrix[0][5] - 0.9).abs() < 0.001)
}

///|
test "sample_sequences_validation" {
  let seqs = @src.motif_scan_sample_sequences()
  assert_eq(seqs.length(), 3)
  assert_eq(seqs[0].0, "seq1")
  assert_eq(seqs[0].1, "GGGGATGCAAGGGG")
  // seq1 contains ATGCAA at position 4.
  assert_eq(seqs[0].1[4:10].to_owned(), "ATGCAA")
  // seq3 contains ATGCAA at positions 0 and 6.
  assert_eq(seqs[2].1[0:6].to_owned(), "ATGCAA")
  assert_eq(seqs[2].1[6:12].to_owned(), "ATGCAA")
}

///|
test "edge_case_empty_sequence" {
  let pwm = @src.motif_scan_sample_pwm()
  let pssm = @src.build_pssm(pwm)
  let matches = @src.scan_sequence(pssm, "", "empty", 0.001)
  assert_eq(matches.length(), 0)
}

///|
test "edge_case_sequence_shorter_than_motif" {
  let pwm = @src.motif_scan_sample_pwm()
  let pssm = @src.build_pssm(pwm)
  // Motif is 6 wide; sequence is 5.
  let matches = @src.scan_sequence(pssm, "ATGCA", "short", 0.001)
  assert_eq(matches.length(), 0)
}

///|
test "edge_case_all_same_base_sequence" {
  let pwm = @src.motif_scan_sample_pwm()
  let pssm = @src.build_pssm(pwm)
  // "AAAAAA" is exactly motif width; score should be computable and finite.
  let score = @src.compute_score(pssm, "AAAAAA")
  let pval = @src.compute_pvalue(pssm, score)
  assert_true(pval >= 0.0)
  assert_true(pval <= 1.0)
  // "AAAAAA" is not a perfect match (consensus is ATGCAA), so p-value is larger
  // than the perfect-match p-value.
  let perfect_pval = @src.compute_pvalue(pssm, pssm.max_score)
  assert_true(pval > perfect_pval)
}

///|
test "edge_case_perfect_match_scan" {
  let pwm = @src.motif_scan_sample_pwm()
  let pssm = @src.build_pssm(pwm)
  // Sequence exactly matching the consensus.
  let matches = @src.scan_sequence(pssm, "ATGCAA", "perfect", 0.001)
  assert_true(matches.length() >= 1)
  let mut found = false
  for i in 0..<matches.length() {
    if matches[i].start == 0 && matches[i].strand == "+" {
      assert_eq(matches[i].matched_sequence, "ATGCAA")
      assert_true(matches[i].p_value < 0.001)
      found = true
    }
  }
  assert_true(found)
}

///|
test "reverse_complement_basic" {
  assert_eq(@src.reverse_complement_sequence("ATGC"), "GCAT")
  assert_eq(@src.reverse_complement_sequence("AATTCCGG"), "CCGGAATT")
  assert_eq(@src.reverse_complement_sequence("AAAA"), "TTTT")
}

///|
test "reverse_complement_palindrome" {
  // ACGT is its own reverse complement.
  assert_eq(@src.reverse_complement_sequence("ACGT"), "ACGT")
  assert_eq(@src.reverse_complement_sequence("AATATT"), "AATATT")
}

///|
test "reverse_complement_empty_and_case" {
  assert_eq(@src.reverse_complement_sequence(""), "")
  // Lowercase input is handled (uppercased before complementing).
  assert_eq(@src.reverse_complement_sequence("atgc"), "GCAT")
}

///|
test "reverse_complement_non_acgt_passthrough" {
  // Non-ACGT characters pass through unchanged (uppercased).
  assert_eq(@src.reverse_complement_sequence("ANNT"), "ANNT")
}

///|
test "pwm_from_counts_empty" {
  let counts : Array[Array[Int]] = []
  let pwm = @src.motif_scan_pwm_from_counts("E", "Empty", counts, 1.0)
  assert_eq(pwm.matrix.length(), 0)
}

///|
test "compute_score_empty_pssm" {
  let pwm = @src.PositionWeightMatrix::new(
    motif_id="E",
    motif_name="E",
    matrix=[],
  )
  let pssm = @src.build_pssm(pwm)
  assert_true((@src.compute_score(pssm, "ACGT") - 0.0).abs() < 0.001)
  assert_true((@src.compute_pvalue(pssm, 0.0) - 1.0).abs() < 0.001)
}
