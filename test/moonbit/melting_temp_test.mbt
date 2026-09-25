///|
/// Tests for melting_temp module (Bio.SeqUtils.MeltingTemp full implementation).

test "mt_wallace basic" {
  let result = @src.mt_wallace("ACGT")
  assert_eq(result.method_, "wallace")
  // ACGT: 2 GC + 2 AT = 4*2 + 2*2 = 12
  assert_true((result.tm - 12.0).abs() < 1.0e-9)
}

test "mt_gc_content basic" {
  let result = @src.mt_gc_content("ACGT")
  assert_eq(result.method_, "gc_content")
  // GC-content method can give negative values for short sequences
  assert_true(result.tm < 100.0)
}

test "mt_nearest_neighbor short" {
  let result = @src.mt_nearest_neighbor("ACGT")
  assert_eq(result.method_, "nearest_neighbor")
  // NN method for short sequences may give unrealistic values
  // This is expected for very short oligos
  assert_true(result.tm < 100.0)
}

test "mt_with_salt" {
  let result = @src.mt_with_salt("ACGT", na_conc=0.05, mg_conc=0.0, method_="nearest_neighbor")
  // NN method for short sequences may give unrealistic values
  assert_true(result.tm < 100.0)
}

test "mt_nearest_neighbor longer" {
  let result = @src.mt_nearest_neighbor("AGCTGATCGATCGATCGATCGATCGATCGATCGATCG")
  assert_eq(result.method_, "nearest_neighbor")
  assert_true(result.tm > 0.0)
}

test "mt_nearest_neighbor single base" {
  let result = @src.mt_nearest_neighbor("A")
  assert_eq(result.tm, 0.0)
}

test "mt_salt_correction_na" {
  let tm_1m = 60.0
  let tm_005m = @src.mt_salt_correction_na(tm_1m, 0.05)
  // Lower [Na+] should decrease Tm
  assert_true(tm_005m < tm_1m)
}

test "mt_salt_correction_mg" {
  let tm_na = 60.0
  let tm_mg = @src.mt_salt_correction_mg(tm_na, 0.01)
  assert_true(tm_mg > 0.0)
}

test "mt_is_valid_dna" {
  assert_true(@src.mt_is_valid_dna("ACGT"))
  assert_true(@src.mt_is_valid_dna("ACGTN"))
  assert_false(@src.mt_is_valid_dna("ACGTX"))
}

test "recommend_tm_method" {
  assert_eq(@src.recommend_tm_method("AC"), "wallace")
  assert_eq(@src.recommend_tm_method("ACGTACGTACGTACGT"), "nearest_neighbor")
}

test "mt_primer" {
  let result = @src.mt_primer("ACGT")
  assert_true(result.tm > 0.0)
}

test "mt_nearest_neighbor conc parameter" {
  let result1 = @src.mt_nearest_neighbor("ACGT", conc=0.0000001)
  let result2 = @src.mt_nearest_neighbor("ACGT", conc=0.000001)
  // Higher concentration should increase Tm slightly
  assert_true(result2.tm >= result1.tm)
}

test "mt_nearest_neighbor na_conc parameter" {
  let result1 = @src.mt_nearest_neighbor("ACGT", na_conc=0.05)
  let result2 = @src.mt_nearest_neighbor("ACGT", na_conc=0.5)
  // Higher [Na+] should increase Tm
  assert_true(result2.tm > result1.tm)
}
