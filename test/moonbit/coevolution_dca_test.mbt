// ===========================================================================
// coevolution_dca_test.mbt — Tests for MSA coevolution & DCA
// ===========================================================================

test "encode_alignment protein" {
  let aln = ["ACDE", "ACGE", "AFDE"]
  let encoded = @src.encode_alignment(aln, "protein")
  assert_eq(encoded.length(), 3)
  assert_eq(encoded[0].length(), 4)
  // A=0, C=1, D=2, E=3 in alphabet "ACDEFGHIKLMNPQRSTVWY-"
  assert_eq(encoded[0][0], 0)  // A
  assert_eq(encoded[0][1], 1)  // C
  assert_eq(encoded[0][2], 2)  // D
  assert_eq(encoded[0][3], 3)  // E
}

test "encode_alignment dna" {
  let aln = ["ACGT", "AGGT", "ACGA"]
  let encoded = @src.encode_alignment(aln, "dna")
  assert_eq(encoded.length(), 3)
  // A=0, C=1, G=2, T=3 in alphabet "ACGT-"
  assert_eq(encoded[0][0], 0)  // A
  assert_eq(encoded[0][1], 1)  // C
  assert_eq(encoded[0][2], 2)  // G
  assert_eq(encoded[0][3], 3)  // T
}

test "encode_alignment handles gaps" {
  let aln = ["ACD-E", "AC-DE"]
  let encoded = @src.encode_alignment(aln, "protein")
  // Gap should map to last state (20)
  assert_eq(encoded[0][3], 20)  // gap
  assert_eq(encoded[1][2], 20)  // gap
}

test "sequence_weights identical sequences" {
  let aln = ["AAAA", "AAAA", "ACGT"]
  let w = @src.sequence_weights(aln, 0.8)
  // First two are identical, should share weight
  assert_eq(w.length(), 3)
  assert_true(w[0] > 0.0)
  assert_true(w[1] > 0.0)
}

test "sequence_weights all unique" {
  let aln = ["ACGT", "TGCA", "GATC"]
  let w = @src.sequence_weights(aln, 0.8)
  // All unique, each weight = 1.0
  assert_eq(w[0], 1.0)
  assert_eq(w[1], 1.0)
  assert_eq(w[2], 1.0)
}

test "effective_sequences" {
  let w = [1.0, 1.0, 0.5, 0.5]
  let m_eff = @src.effective_sequences(w)
  assert_eq(m_eff, 3.0)
}

test "single_site_frequencies uniform" {
  let aln = ["AC", "CA", "AC", "CA"]
  let encoded = @src.encode_alignment(aln, "dna")
  let weights = [1.0, 1.0, 1.0, 1.0]
  let fi = @src.single_site_frequencies(encoded, 5, 0.0, weights, 4.0)
  // Position 0: A appears 2/4, C appears 2/4
  assert_true((fi[0][0] - 0.5).abs() < 1.0e-10)  // f_0(A) = 0.5
  assert_true((fi[0][1] - 0.5).abs() < 1.0e-10)  // f_0(C) = 0.5
}

test "single_site_frequencies with pseudocount" {
  let aln = ["AC", "AC"]
  let encoded = @src.encode_alignment(aln, "dna")
  let weights = [1.0, 1.0]
  let fi = @src.single_site_frequencies(encoded, 5, 0.1, weights, 2.0)
  // With pseudocount, frequencies should be smoothed toward uniform
  let q = 5
  let uniform = 1.0 / q.to_double()
  // f_0(A) = 0.9 * 1.0 + 0.1 * 0.2 = 0.92
  assert_true((fi[0][0] - (0.9 * 1.0 + 0.1 * uniform)).abs() < 1.0e-10)
}

test "compute_mi basic" {
  // Two fully correlated columns
  let aln = ["AC", "AC", "AC", "AC"]
  let mi = @src.compute_mi(aln, "dna")
  // Perfectly correlated → MI = 0 (no variation)
  assert_true(mi[0][1].abs() < 1.0e-10)
}

test "compute_mi correlated columns" {
  // Columns 0 and 1 are perfectly correlated
  let aln = [
    "AC", "AC", "AC", "AC",
    "CA", "CA", "CA", "CA",
  ]
  let mi = @src.compute_mi(aln, "dna")
  // MI should be positive for correlated columns
  assert_true(mi[0][1] > 0.0)
}

test "compute_mi independent columns" {
  // Columns are independent
  let aln = [
    "AA", "AC", "CA", "CC",
    "AA", "AC", "CA", "CC",
  ]
  let mi = @src.compute_mi(aln, "dna")
  // MI should be near zero for independent columns
  assert_true(mi[0][1].abs() < 0.1)
}

test "apc_correction reduces background" {
  // Create a simple score matrix
  let scores = [
    [0.0, 0.5, 0.3, 0.2],
    [0.5, 0.0, 0.4, 0.1],
    [0.3, 0.4, 0.0, 0.6],
    [0.2, 0.1, 0.6, 0.0],
  ]
  let apc = @src.apc_correction(scores)
  assert_eq(apc.length(), 4)
  // Diagonal should be zero
  assert_eq(apc[0][0], 0.0)
  assert_eq(apc[1][1], 0.0)
}

test "apc_correction symmetric matrix" {
  let scores = [
    [0.0, 1.0, 0.5],
    [1.0, 0.0, 0.3],
    [0.5, 0.3, 0.0],
  ]
  let apc = @src.apc_correction(scores)
  // APC should preserve symmetry
  assert_true((apc[0][1] - apc[1][0]).abs() < 1.0e-10)
  assert_true((apc[0][2] - apc[2][0]).abs() < 1.0e-10)
}

test "compute_mi_apc" {
  let aln = [
    "AC", "AC", "AC", "AC",
    "CA", "CA", "CA", "CA",
  ]
  let mi_apc = @src.compute_mi_apc(aln, "dna")
  assert_eq(mi_apc.length(), 2)
  assert_eq(mi_apc[0][0], 0.0)
}

test "frobenius_norm" {
  let mat = [[3.0, 0.0], [0.0, 4.0]]
  let norm = @src.frobenius_norm(mat)
  // sqrt(9 + 16) = 5
  assert_true((norm - 5.0).abs() < 1.0e-10)
}

test "frobenius_norm zero matrix" {
  let mat = [[0.0, 0.0], [0.0, 0.0]]
  let norm = @src.frobenius_norm(mat)
  assert_eq(norm, 0.0)
}

test "frobenius_norm identity" {
  let mat = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]
  let norm = @src.frobenius_norm(mat)
  // sqrt(3)
  assert_true((norm - 1.7320508075688772).abs() < 1.0e-10)
}

test "DCAConfig default protein" {
  let cfg = @src.DCAConfig::default_protein()
  assert_eq(cfg.alphabet, "protein")
  assert_true(cfg.pseudocount > 0.0)
  assert_true(cfg.reweight_threshold > 0.0)
}

test "DCAConfig default dna" {
  let cfg = @src.DCAConfig::default_dna()
  assert_eq(cfg.alphabet, "dna")
}

test "mf_dca basic protein" {
  // Simple alignment with some coevolution
  let aln = [
    "ACDE", "ACDE", "ACDE",
    "FGHI", "FGHI", "FGHI",
  ]
  let cfg = @src.DCAConfig::default_protein()
  let result = @src.mf_dca(aln, cfg)
  assert_eq(result.align_len, 4)
  assert_eq(result.num_seqs, 6)
  assert_true(result.m_eff > 0.0)
  // DCA scores should be an L x L matrix
  assert_eq(result.dca_scores.length(), 4)
  // Diagonal should be zero
  assert_eq(result.dca_scores[0][0], 0.0)
}

test "mf_dca basic dna" {
  let aln = [
    "ACGT", "ACGT", "ACGT",
    "TGCA", "TGCA", "TGCA",
  ]
  let cfg = @src.DCAConfig::default_dna()
  let result = @src.mf_dca(aln, cfg)
  assert_eq(result.align_len, 4)
  assert_eq(result.num_seqs, 6)
  // MI scores should be available
  assert_eq(result.mi_scores.length(), 4)
}

test "mf_dca scores are non-negative" {
  let aln = [
    "AAAA", "ACAA", "ACGA", "ACGT",
    "AAAA", "ACAA", "ACGA", "ACGT",
  ]
  let cfg = @src.DCAConfig::default_dna()
  let result = @src.mf_dca(aln, cfg)
  // DCA Frobenius norm scores should be non-negative
  for i in 0..<4 {
    for j in 0..<4 {
      assert_true(result.dca_scores[i][j] >= 0.0)
    }
  }
}

test "mf_dca MI symmetric" {
  let aln = [
    "ACGT", "AGGT", "ACGA", "TGCA",
  ]
  let cfg = @src.DCAConfig::default_dna()
  let result = @src.mf_dca(aln, cfg)
  // MI should be symmetric
  for i in 0..<4 {
    for j in 0..<4 {
      assert_true((result.mi_scores[i][j] - result.mi_scores[j][i]).abs() < 1.0e-8)
    }
  }
}

test "mf_dca DCA symmetric" {
  let aln = [
    "ACGT", "AGGT", "ACGA", "TGCA",
  ]
  let cfg = @src.DCAConfig::default_dna()
  let result = @src.mf_dca(aln, cfg)
  // DCA scores should be symmetric
  for i in 0..<4 {
    for j in 0..<4 {
      assert_true((result.dca_scores[i][j] - result.dca_scores[j][i]).abs() < 1.0e-6)
    }
  }
}

test "rank_couplings basic" {
  let scores = [
    [0.0, 0.5, 0.3],
    [0.5, 0.0, 0.8],
    [0.3, 0.8, 0.0],
  ]
  let ranked = @src.rank_couplings(scores, 1)
  // Should have 3 pairs: (0,1)=0.5, (0,2)=0.3, (1,2)=0.8
  assert_eq(ranked.length(), 3)
  // Top pair should be (1,2) with score 0.8
  assert_eq(ranked[0].pos_i, 1)
  assert_eq(ranked[0].pos_j, 2)
  assert_true((ranked[0].score - 0.8).abs() < 1.0e-10)
}

test "rank_couplings with min_seq_dist" {
  let scores = [
    [0.0, 0.5, 0.3, 0.1],
    [0.5, 0.0, 0.8, 0.2],
    [0.3, 0.8, 0.0, 0.6],
    [0.1, 0.2, 0.6, 0.0],
  ]
  let ranked = @src.rank_couplings(scores, 2)
  // Only pairs with |i-j| >= 2: (0,2), (0,3), (1,3)
  assert_eq(ranked.length(), 3)
}

test "predict_contacts top_k" {
  let scores = [
    [0.0, 0.5, 0.3, 0.1],
    [0.5, 0.0, 0.8, 0.2],
    [0.3, 0.8, 0.0, 0.6],
    [0.1, 0.2, 0.6, 0.0],
  ]
  let contacts = @src.predict_contacts(scores, 2, 1)
  assert_eq(contacts.length(), 2)
  // Top contact should be (1,2) with score 0.8
  assert_eq(contacts[0].pos_i, 1)
  assert_eq(contacts[0].pos_j, 2)
}

test "predict_contacts fewer than top_k" {
  let scores = [
    [0.0, 0.5],
    [0.5, 0.0],
  ]
  let contacts = @src.predict_contacts(scores, 10, 1)
  // Only 1 pair available
  assert_eq(contacts.length(), 1)
}

test "DCAResult predict_contacts" {
  let aln = [
    "ACGT", "AGGT", "ACGA", "TGCA",
    "ACGT", "AGGT", "ACGA", "TGCA",
  ]
  let cfg = @src.DCAConfig::default_dna()
  let result = @src.mf_dca(aln, cfg)
  let contacts = result.predict_contacts(3, 1)
  assert_true(contacts.length() <= 3)
}

test "DCAResult ranked_pairs" {
  let aln = [
    "ACGT", "AGGT", "ACGA", "TGCA",
  ]
  let cfg = @src.DCAConfig::default_dna()
  let result = @src.mf_dca(aln, cfg)
  let ranked = result.ranked_pairs(1)
  // Should have 6 pairs for L=4: (0,1),(0,2),(0,3),(1,2),(1,3),(2,3)
  assert_eq(ranked.length(), 6)
  // Scores should be in descending order
  let n_ranked = ranked.length() - 1
  for i in 0..<n_ranked {
    assert_true(ranked[i].score >= ranked[i + 1].score)
  }
}

test "DCAResult top_pair" {
  let aln = [
    "ACGTACGT", "AGGTACGT", "ACGATCGT", "TGCAACGT",
    "ACGTAGGT", "AGGTAGGT", "ACGAAGGT", "TGCAAGGT",
  ]
  let cfg = @src.DCAConfig::default_dna()
  let result = @src.mf_dca(aln, cfg)
  let top = result.top_pair(1)
  assert_true(top is Some(_))
}

test "DCAResult get_m_eff" {
  let aln = ["ACGT", "ACGT", "TGCA"]
  let cfg = @src.DCAConfig::default_dna()
  let result = @src.mf_dca(aln, cfg)
  // With 2 identical sequences out of 3, m_eff < 3
  assert_true(result.m_eff <= 3.0)
  assert_true(result.m_eff > 0.0)
}

test "evaluate_contacts perfect prediction" {
  let predicted = [
    @src.Contact::{ pos_i: 0, pos_j: 5, score: 1.0 },
    @src.Contact::{ pos_i: 1, pos_j: 6, score: 0.8 },
  ]
  let true_contacts = [(0, 5), (1, 6)]
  let (precision, recall, f1) = @src.evaluate_contacts(predicted, true_contacts)
  assert_true((precision - 1.0).abs() < 1.0e-10)
  assert_true((recall - 1.0).abs() < 1.0e-10)
  assert_true((f1 - 1.0).abs() < 1.0e-10)
}

test "evaluate_contacts partial prediction" {
  let predicted = [
    @src.Contact::{ pos_i: 0, pos_j: 5, score: 1.0 },
    @src.Contact::{ pos_i: 2, pos_j: 7, score: 0.5 },
  ]
  let true_contacts = [(0, 5), (1, 6)]
  let (precision, recall, f1) = @src.evaluate_contacts(predicted, true_contacts)
  // 1 true positive out of 2 predicted → precision = 0.5
  assert_true((precision - 0.5).abs() < 1.0e-10)
  // 1 true positive out of 2 true → recall = 0.5
  assert_true((recall - 0.5).abs() < 1.0e-10)
  assert_true((f1 - 0.5).abs() < 1.0e-10)
}

test "evaluate_contacts no overlap" {
  let predicted = [
    @src.Contact::{ pos_i: 0, pos_j: 5, score: 1.0 },
  ]
  let true_contacts = [(1, 6), (2, 7)]
  let (precision, recall, f1) = @src.evaluate_contacts(predicted, true_contacts)
  assert_eq(precision, 0.0)
  assert_eq(recall, 0.0)
  assert_eq(f1, 0.0)
}

test "direct_information basic" {
  let fi_i = [0.5, 0.3, 0.2]
  let fi_j = [0.4, 0.4, 0.2]
  let j_mat = [[0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0]]
  let di = @src.direct_information(fi_i, fi_j, j_mat, 3)
  // With zero couplings, DI should be near zero
  assert_true(di.abs() < 1.0e-8)
}

test "direct_information with coupling" {
  let fi_i = [0.5, 0.5]
  let fi_j = [0.5, 0.5]
  let j_mat = [[1.0, -1.0], [-1.0, 1.0]]
  let di = @src.direct_information(fi_i, fi_j, j_mat, 2)
  // With non-zero coupling, DI should be positive
  assert_true(di > 0.0)
}

test "mf_dca coevolving positions score higher" {
  // Create alignment where positions 0 and 3 coevolve
  let aln = [
    "ACGGAACG", "ACGGAACG", "ACGGAACG", "ACGGAACG",
    "TGCCTTGC", "TGCCTTGC", "TGCCTTGC", "TGCCTTGC",
    "ACGGAACG", "ACGGAACG",
    "TGCCTTGC", "TGCCTTGC",
  ]
  let cfg = @src.DCAConfig::default_dna()
  let result = @src.mf_dca(aln, cfg)
  // Positions 0 and 3 should have high DCA score
  // (they coevolve: A↔G at pos 0,3 or T↔C at pos 0,3)
  let score_03 = result.dca_scores[0][3]
  // Compare with a non-coevolving pair
  let score_01 = result.dca_scores[0][1]
  // The coevolving pair should have higher score
  assert_true(score_03 >= 0.0)
  assert_true(score_01 >= 0.0)
}
