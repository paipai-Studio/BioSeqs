///|
test "seq_complexity: shannon_entropy_bits_uniform" {
  let entropy = @src.shannon_entropy_bits("ATCGATCG")
  assert_true(entropy > 1.9 && entropy < 2.1)
}

///|
test "seq_complexity: shannon_entropy_bits_low" {
  let entropy = @src.shannon_entropy_bits("AAAAAAAA")
  assert_true(entropy < 0.01)
}

///|
test "seq_complexity: shannon_entropy_bits_empty" {
  let entropy = @src.shannon_entropy_bits("")
  assert_eq(entropy, 0.0)
}

///|
test "seq_complexity: linguistic_complexity" {
  let lc = @src.linguistic_complexity("ATCGATCG", k=2)
  assert_true(lc > 0.5)
  let lc2 = @src.linguistic_complexity("ATATATAT", k=2)
  assert_true(lc2 < lc)
}

///|
test "seq_complexity: linguistic_complexity_empty" {
  let lc = @src.linguistic_complexity("", k=2)
  assert_eq(lc, 0.0)
}

///|
test "seq_complexity: gc_content_percent" {
  assert_true(@src.gc_content_percent("GCGCGCGC") > 99.0)
  assert_true(@src.gc_content_percent("ATATATAT") < 1.0)
  assert_true(
    @src.gc_content_percent("ATCGATCG") > 49.0 &&
    @src.gc_content_percent("ATCGATCG") < 51.0,
  )
}

///|
test "seq_complexity: at_content" {
  assert_true(@src.at_content("ATATATAT") > 99.0)
  assert_true(@src.at_content("GCGCGCGC") < 1.0)
}

///|
test "seq_complexity: gc_skew" {
  let skew = @src.gc_skew("GGGGCC")
  assert_true(skew > 0.0)
  let skew2 = @src.gc_skew("CCCCGG")
  assert_true(skew2 < 0.0)
  let skew3 = @src.gc_skew("GCGC")
  assert_true(skew3 < 0.01 && skew3 > -0.01)
}

///|
test "seq_complexity: at_skew" {
  let skew = @src.at_skew("AAAATT")
  assert_true(skew > 0.0)
  let skew2 = @src.at_skew("TTTTAA")
  assert_true(skew2 < 0.0)
}

///|
test "seq_complexity: dust_score" {
  let dust_repeat = @src.dust_score("ATGATGATGATG")
  let dust_diverse = @src.dust_score("ATCGATCGATCG")
  assert_true(dust_repeat > 0.0)
  assert_true(dust_diverse > 0.0)
}

///|
test "seq_complexity: sequence_signature" {
  let sig = @src.sequence_signature("ATCG", k=2)
  assert_eq(sig.size(), 3)
  let freq = sig.get("AT")
  assert_true(freq is Some(_))
  assert_true(freq.unwrap() > 0.32 && freq.unwrap() < 0.34)
}

///|
test "seq_complexity: nucleotide_frequencies" {
  let (a, t, g, c) = @src.nucleotide_frequencies("ATCG")
  assert_true(a > 0.24 && a < 0.26)
  assert_true(t > 0.24 && t < 0.26)
  assert_true(g > 0.24 && g < 0.26)
  assert_true(c > 0.24 && c < 0.26)
}

///|
test "seq_complexity: nucleotide_frequencies_empty" {
  let (a, t, g, c) = @src.nucleotide_frequencies("")
  assert_eq(a, 0.0)
  assert_eq(t, 0.0)
  assert_eq(g, 0.0)
  assert_eq(c, 0.0)
}

///|
test "seq_complexity: chaos_game_representation" {
  let (x, y) = @src.chaos_game_representation("AAAA")
  assert_true(x < 0.1)
  assert_true(y < 0.1)
  let (x2, y2) = @src.chaos_game_representation("GGGG")
  assert_true(x2 > 0.9)
  assert_true(y2 > 0.9)
}

///|
test "seq_complexity: sequence_similarity_identical" {
  let sim = @src.sequence_similarity("ATCGATCG", "ATCGATCG", k=2)
  assert_true(sim > 0.99)
}

///|
test "seq_complexity: sequence_similarity_different" {
  let sim = @src.sequence_similarity("ATCGATCG", "GCTAGCTA", k=2)
  assert_true(sim < 1.0)
}

///|
test "seq_complexity: sequence_similarity_empty" {
  let sim = @src.sequence_similarity("", "", k=2)
  assert_eq(sim, 0.0)
}
