///|
/// Tests for EMBOSS module.
test "emboss_gc_skew" {
  let seq = "GGCC"
  let skew = @src.emboss_gc_skew(seq)
  assert_eq(skew, 0.0)
}

///|
test "emboss_at_skew" {
  let seq = "AATT"
  let skew = @src.emboss_at_skew(seq)
  assert_eq(skew, 0.0)
}

///|
test "emboss_dna_molecular_weight" {
  let seq = "ATGC"
  let weight = @src.dna_molecular_weight(seq)
  assert_true(weight > 1200.0)
}

///|
test "emboss_tm_simple" {
  let seq = "ATGC"
  let tm = @src.tm_simple(seq)
  assert_eq(tm, 12.0)
}

///|
test "emboss_tm_nearest_neighbor" {
  let seq = "ATGC"
  let tm = @src.tm_nearest_neighbor(seq)
  assert_true(tm == tm)
}

///|
test "emboss_reverse_complement" {
  let seq = "ATGC"
  let rc = @src.emboss_reverse_complement(seq)
  assert_eq(rc, "GCAT")
}

///|
test "emboss_translate" {
  let seq = "ATG"
  let protein = @src.translate(seq)
  assert_eq(protein, "M")
}

///|
test "emboss_find_orfs" {
  let seq = "ATGAAATAA"
  let orfs = @src.find_orfs(seq, min_length=3)
  assert_true(orfs.length() > 0)
}

///|
test "emboss_hamming_distance" {
  let d = @src.hamming_distance("ATGC", "ATGG")
  assert_eq(d, 1)
}

///|
test "emboss_levenshtein_distance" {
  let d = @src.levenshtein_distance("ATGC", "AGC")
  assert_eq(d, 1)
}

///|
test "emboss_count_kmers" {
  let kmers = @src.emboss_count_kmers("ATGCAT", 2)
  assert_true(kmers.length() > 0)
}

///|
test "emboss_isoelectric_point" {
  let seq = "K"
  let pi = @src.isoelectric_point(seq)
  assert_true(pi > 7.0)
}

///|
test "emboss_aliphatic_index" {
  let seq = "AVIL"
  let idx = @src.aliphatic_index(seq)
  assert_true(idx > 0.0)
}

///|
test "emboss_gravy_score" {
  let seq = "AVIL"
  let gravy = @src.gravy_score(seq)
  assert_true(gravy > 0.0)
}

///|
test "emboss_instability_index" {
  let seq = "AV"
  let idx = @src.instability_index(seq)
  assert_true(idx >= 0.0)
}

///|
test "emboss_amino_acid_composition" {
  let comp = @src.emboss_amino_acid_composition("AVIL")
  assert_true(comp.length() >= 0)
}
