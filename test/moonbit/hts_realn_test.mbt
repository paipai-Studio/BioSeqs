///|
test "probaln_glocal exact match" {
  // Reference and query are identical
  let ref = @src.seq_to_ints("ACGTACGT")
  let query = @src.seq_to_ints("ACGTACGT")
  let qual = [30, 30, 30, 30, 30, 30, 30, 30]
  let params = @src.ProbalnPar::make(0.001, 0.1, 10)
  let result = @src.probaln_glocal(ref, query, qual, params)
  // Score should be relatively good (low) for exact match
  assert_true(result.score < 50)
  // State should indicate matches (bit 0-1 = 0) for all positions
  for i = 0; i < 8; i = i + 1 {
    assert_eq(result.state[i] & 3, 0)
  }
}

///|
test "probaln_glocal one_mismatch" {
  let ref = @src.seq_to_ints("ACGTACGT")
  let query = @src.seq_to_ints("ACTTACGT")  // position 2: G->T
  let qual = [30, 30, 30, 30, 30, 30, 30, 30]
  let params = @src.ProbalnPar::make(0.001, 0.1, 10)
  let result = @src.probaln_glocal(ref, query, qual, params)
  // Should still align all positions as matches
  assert_eq(result.state.length(), 8)
}

///|
test "probaln_glocal insertion" {
  // Query has an extra base
  let ref = @src.seq_to_ints("ACGT")
  let query = @src.seq_to_ints("ACAGT")  // A inserted at position 2
  let qual = [30, 30, 30, 30, 30]
  let params = @src.ProbalnPar::make(0.001, 0.1, 10)
  let result = @src.probaln_glocal(ref, query, qual, params)
  assert_eq(result.state.length(), 5)
}

///|
test "probaln_glocal empty" {
  let ref : Array[Int] = []
  let query : Array[Int] = []
  let qual : Array[Int] = []
  let params = @src.ProbalnPar::make(0.001, 0.1, 10)
  let result = @src.probaln_glocal(ref, query, qual, params)
  assert_eq(result.score, 0)
  assert_eq(result.state.length(), 0)
}

///|
test "probaln_glocal short_match" {
  let ref = @src.seq_to_ints("ACTTC")
  let query = @src.seq_to_ints("ATTC")
  let qual = [30, 30, 30, 30]
  let params = @src.ProbalnPar::make(0.001, 0.1, 10)
  let result = @src.probaln_glocal(ref, query, qual, params)
  assert_eq(result.state.length(), 4)
  // Score should be reasonable
  assert_true(result.score >= 0)
}

///|
test "probaln_glocal quality_effect" {
  // Same sequences, different qualities
  let ref = @src.seq_to_ints("ACGTACGT")
  let query = @src.seq_to_ints("ACGTACGT")
  let qual_high = [40, 40, 40, 40, 40, 40, 40, 40]
  let qual_low = [10, 10, 10, 10, 10, 10, 10, 10]
  let params = @src.ProbalnPar::make(0.001, 0.1, 10)
  let result_high = @src.probaln_glocal(ref, query, qual_high, params)
  let result_low = @src.probaln_glocal(ref, query, qual_low, params)
  // Both should succeed
  assert_eq(result_high.state.length(), 8)
  assert_eq(result_low.state.length(), 8)
}

///|
test "base_to_int" {
  assert_eq(@src.base_to_int('A'), 0)
  assert_eq(@src.base_to_int('C'), 1)
  assert_eq(@src.base_to_int('G'), 2)
  assert_eq(@src.base_to_int('T'), 3)
  assert_eq(@src.base_to_int('N'), 4)
  assert_eq(@src.base_to_int('a'), 0)
  assert_eq(@src.base_to_int('c'), 1)
}

///|
test "seq_to_ints" {
  let result = @src.seq_to_ints("ACGTN")
  assert_eq(result.length(), 5)
  assert_eq(result[0], 0)
  assert_eq(result[1], 1)
  assert_eq(result[2], 2)
  assert_eq(result[3], 3)
  assert_eq(result[4], 4)
}

///|
test "errmod_init" {
  let em = @src.errmod_init(0.02)
  // fk[0] should be 1.0
  assert_true((em.fk[0] - 1.0).abs() < 0.000000001)
  // fk should be decreasing with depcorr > 0
  assert_true(em.fk[1] < 1.0)
}

///|
test "errmod_cal homozygous" {
  let em = @src.errmod_init(0.02)
  // All bases are A with quality 30
  // Encoding: quality<<5 | base = 30<<5 | 0 = 960
  let bases = [960, 960, 960, 960]
  let q = @src.errmod_cal(em, bases, 4)
  // Homozygous AA should have low likelihood (good)
  // Heterozygous should have higher likelihood (worse)
  assert_true(q[0] <= q[1])  // AA <= AC
  assert_true(q[0] <= q[5])  // AA <= CC
}

///|
test "errmod_cal mixed_bases" {
  let em = @src.errmod_init(0.02)
  // Mix of A and C bases
  // A with qual 30: 30<<5 | 0 = 960
  // C with qual 30: 30<<5 | 1 = 961
  let bases = [960, 960, 961, 961]
  let q = @src.errmod_cal(em, bases, 4)
  // Heterozygous AC should be better than homozygous AA or CC
  assert_true(q[1] <= q[0])  // AC <= AA
  assert_true(q[1] <= q[5])  // AC <= CC
}

///|
test "errmod_cal empty" {
  let em = @src.errmod_init(0.02)
  let q = @src.errmod_cal(em, [], 4)
  // All likelihoods should be 0
  for i = 0; i < 16; i = i + 1 {
    assert_true((q[i] - 0.0).abs() < 0.000000001)
  }
}

///|
test "sam_cap_mapq exact_match" {
  let ref = "ACGTACGTACGT"
  let query = "ACGTACGTACGT"
  let qual = [30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30]
  let cigar = [(0, 12)]  // 12M
  let result = @src.sam_cap_mapq(ref, query, qual, cigar, 0, 40)
  // No mismatches, should return high quality
  assert_true(result >= 0)
}

///|
test "sam_cap_mapq with_mismatch" {
  let ref = "ACGTACGTACGT"
  let query = "ACTTACGTACGT"  // mismatch at position 2
  let qual = [30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30]
  let cigar = [(0, 12)]
  let result = @src.sam_cap_mapq(ref, query, qual, cigar, 0, 40)
  // Should still be valid
  assert_true(result >= -1)
}

///|
test "sam_cap_mapq many_mismatches" {
  let ref = "ACGTACGTACGT"
  let query = "TTTTTTTTTTTT"  // many mismatches
  let qual = [30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30, 30]
  let cigar = [(0, 12)]
  let result = @src.sam_cap_mapq(ref, query, qual, cigar, 0, 40)
  // Should return -1 (too many mismatches)
  assert_eq(result, -1)
}
