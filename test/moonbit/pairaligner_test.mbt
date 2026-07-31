///|
/// Test file for pairaligner module.

test "pairaligner_alignment_mode_global" {
  let mode = @src.pairaligner_global()
  let config = @src.PairwiseAlignerConfig::default_dna().set_mode(mode)
  assert_eq(config.mode, mode)
}

test "pairaligner_alignment_mode_local" {
  let mode = @src.pairaligner_local()
  let config = @src.PairwiseAlignerConfig::default_dna().set_mode(mode)
  assert_eq(config.mode, mode)
}

test "pairaligner_substitution_matrix_no_matrix" {
  let mat = @src.pairaligner_no_matrix()
  let config = @src.PairwiseAlignerConfig::default_dna().set_submatrix(mat)
  assert_eq(config.submatrix, mat)
}

test "pairaligner_substitution_matrix_blosum62" {
  let mat = @src.pairaligner_blosum62()
  let config = @src.PairwiseAlignerConfig::default_protein().set_submatrix(mat)
  assert_eq(config.submatrix, mat)
}

test "pairaligner_default_dna_config" {
  let config = @src.PairwiseAlignerConfig::default_dna()
  assert_eq(config.match_score, 1.0)
  assert_eq(config.mismatch_score, -1.0)
  assert_eq(config.gap_open, -1.0)
  assert_eq(config.gap_extend, -0.5)
  assert_eq(config.alphabet_type, "DNA")
  assert_eq(config.target_gap_open, -1.0)
  assert_eq(config.query_gap_open, -1.0)
}

test "pairaligner_default_protein_config" {
  let config = @src.PairwiseAlignerConfig::default_protein()
  assert_eq(config.gap_open, -10.0)
  assert_eq(config.gap_extend, -0.5)
  assert_eq(config.alphabet_type, "PROTEIN")
  assert_eq(config.submatrix, @src.pairaligner_blosum62())
}

test "pairaligner_config_setters" {
  let base = @src.PairwiseAlignerConfig::default_dna()
  let c1 = base.set_match_score(5.0)
  assert_eq(c1.match_score, 5.0)
  let c2 = c1.set_mismatch_score(-3.0)
  assert_eq(c2.mismatch_score, -3.0)
  let c3 = c2.set_gap_open(-5.0)
  assert_eq(c3.gap_open, -5.0)
  assert_eq(c3.target_gap_open, -5.0)
  assert_eq(c3.query_gap_open, -5.0)
  let c4 = c3.set_gap_extend(-1.0)
  assert_eq(c4.gap_extend, -1.0)
  assert_eq(c4.target_gap_extend, -1.0)
  let c5 = c4.set_alphabet_type("PROTEIN")
  assert_eq(c5.alphabet_type, "PROTEIN")
  assert_eq(base.match_score, 1.0)
  assert_eq(base.mismatch_score, -1.0)
}

test "pairaligner_align_global_dna" {
  let target = "ACGT"
  let query = "ACGT"
  let config = @src.PairwiseAlignerConfig::default_dna()
  let aln = @src.pairaligner_align(target, query, config=config)
  assert_eq(aln.aligned1(), "ACGT")
  assert_eq(aln.aligned2(), "ACGT")
  assert_eq(aln.identities(), 4)
  assert_eq(aln.alignment_length(), 4)
  let ml = aln.match_line
  assert_eq(ml.length(), 4)
  assert_eq(ml, "||||")
}

test "pairaligner_align_global_dna_with_mismatch" {
  let target = "ACGT"
  let query = "AGGT"
  let config = @src.PairwiseAlignerConfig::default_dna()
  let aln = @src.pairaligner_align(target, query, config=config)
  assert_eq(aln.alignment_length(), 4)
  let ml = aln.match_line
  assert_eq(ml[1:2], ".")
  assert_eq(aln.score(), 2.0)
}

test "pairaligner_align_local_dna" {
  let target = "XXXXACGTXXXX"
  let query = "ACGT"
  let config = @src.PairwiseAlignerConfig::default_dna()
    .set_mode(@src.pairaligner_local())
  let aln = @src.pairaligner_align(target, query, config=config)
  assert_eq(aln.aligned1(), "ACGT")
  assert_eq(aln.aligned2(), "ACGT")
  assert_eq(aln.identities(), 4)
  assert_eq(aln.score(), 4.0)
}

test "pairaligner_align_protein_blosum62" {
  let (target, query) = @src.pairaligner_sample_data()
  let config = @src.PairwiseAlignerConfig::default_protein()
  let aln = @src.pairaligner_align(target, query, config=config)
  assert_true(aln.alignment_length() > 0)
  assert_true(aln.score() > 0.0)
  assert_true(aln.identities() > 0)
}

test "pairaligner_aligned1_aligned2" {
  let aln = @src.pairaligner_align("ACGT", "ACGT")
  assert_eq(aln.aligned1(), aln.aligned_target)
  assert_eq(aln.aligned2(), aln.aligned_query)
}

test "pairaligner_alignment_length" {
  let aln = @src.pairaligner_align("AAAA", "AAAA")
  assert_eq(aln.alignment_length(), 4)
  assert_eq(aln.aligned1().length(), 4)
  assert_eq(aln.aligned2().length(), 4)
}

test "pairaligner_identities_and_identity_pct" {
  let aln = @src.pairaligner_align("ACGT", "ACGT")
  assert_eq(aln.identities(), 4)
  assert_eq(aln.identity_pct(), 100.0)
  let aln2 = @src.pairaligner_align("ACGT", "AGGT")
  assert_eq(aln2.identities(), 3)
  assert_eq(aln2.identity_pct(), 75.0)
}

test "pairaligner_gaps_count" {
  let config = @src.PairwiseAlignerConfig::default_dna()
    .set_gap_open(-2.0)
    .set_gap_extend(-1.0)
  let target = "AAACCC"
  let query = "AAA"
  let aln = @src.pairaligner_align(target, query, config=config)
  let gc = aln.gaps_count()
  assert_true(gc >= 3)
}

test "pairaligner_score_extraction" {
  let aln = @src.pairaligner_align("ACGT", "ACGT")
  assert_eq(aln.score(), 4.0)
  assert_eq(aln.score(), aln.score)
}

test "pairaligner_sample_data" {
  let (t, q) = @src.pairaligner_sample_data()
  assert_eq(t, "HEAGAWGHEE")
  assert_eq(q, "PAWHEAE")
  assert_true(t.length() > 0)
  assert_true(q.length() > 0)
}

test "pairaligner_affine_vs_linear" {
  let target = "AAAAAAAAAA"
  let query = "AAAA"
  let affine_config = @src.PairwiseAlignerConfig::default_dna()
    .set_mode(@src.pairaligner_global())
    .set_gap_open(-5.0)
    .set_gap_extend(-1.0)
  let affine_aln = @src.pairaligner_align(target, query, config=affine_config)
  let linear_config = @src.PairwiseAlignerConfig::default_dna()
    .set_mode(@src.pairaligner_global())
    .set_gap_open(-1.0)
    .set_gap_extend(-1.0)
  let linear_aln = @src.pairaligner_align(target, query, config=linear_config)
  assert_true(affine_aln.alignment_length() > 0)
  assert_true(linear_aln.alignment_length() > 0)
  assert_true(affine_aln.score() < 100.0)
  assert_true(linear_aln.score() < 100.0)
}

test "pairaligner_config_new_named" {
  let config = @src.PairwiseAlignerConfig::new(
    mode=@src.pairaligner_local(),
    match_score=3.0,
    mismatch_score=-2.0,
    gap_open=-3.0,
    gap_extend=-1.0,
    submatrix=@src.pairaligner_no_matrix(),
    alphabet_type="DNA",
  )
  assert_eq(config.match_score, 3.0)
  assert_eq(config.mismatch_score, -2.0)
  assert_eq(config.gap_open, -3.0)
  assert_eq(config.alphabet_type, "DNA")
}

test "pairaligner_local_mode_start_end_positions" {
  let target = "XXACGTYY"
  let query = "ACGT"
  let config = @src.PairwiseAlignerConfig::default_dna()
    .set_mode(@src.pairaligner_local())
  let aln = @src.pairaligner_align(target, query, config=config)
  assert_eq(aln.target_start, 2)
  assert_eq(aln.target_end, 6)
  assert_eq(aln.query_start, 0)
  assert_eq(aln.query_end, 4)
}
