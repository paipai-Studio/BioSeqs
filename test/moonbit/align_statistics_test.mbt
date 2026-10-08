test "align_stats_evalue" {
  // BLOSUM62 defaults
  let stats = @src.blosum62_stats()
  // Score 100, query len 100, db len 1000000
  let e = @src.evalue(100.0, 100, 1000000, stats)
  // E = 0.134 * 100 * 1000000 * exp(-0.318 * 100)
  // = 0.134 * 1e8 * exp(-31.8)
  // = 1.34e7 * 1.6e-14 ≈ 2.1e-7
  assert_true(e > 0.0)
  assert_true(e < 1.0)
}

test "align_stats_pvalue" {
  // p = 1 - exp(-E)
  let p = @src.pvalue(0.01)
  assert_true(p > 0.009 && p < 0.011)
  // Large E → p ≈ 1
  let p_large = @src.pvalue(100.0)
  assert_true(p_large > 0.99)
}

test "align_stats_bitscore" {
  let stats = @src.blosum62_stats()
  let bits = @src.bit_score(100.0, stats)
  // S' = (0.318 * 100 - ln(0.134)) / ln(2)
  // = (31.8 - (-2.01)) / 0.693 ≈ 48.8
  assert_true(bits > 45.0 && bits < 55.0)
}

test "align_stats_blastn" {
  let stats = @src.blastn_stats()
  let e = @src.evalue(50.0, 100, 1000000, stats)
  // E = 0.711 * 100 * 1000000 * exp(-1.37 * 50)
  // = 7.11e7 * exp(-68.5) ≈ very small
  assert_true(e >= 0.0)
  assert_true(e < 1.0)
}

test "align_stats_estimate" {
  // Simple scoring: match=2, mismatch=-1
  // Frequencies: 0.25 each for A, C, G, T
  let scores = [2.0, -1.0]
  let freqs = [0.25, 0.75]  // 25% match, 75% mismatch
  let result = @src.estimate_karlin_altschul(scores, freqs)
  assert_true(result is Some(_))
  match result {
    Some(s) => {
      // Lambda should be positive
      assert_true(s.lambda > 0.0)
      // Expected score should be negative
      assert_true(s.score_mean < 0.0)
    }
    None => assert_true(false)
  }
}
