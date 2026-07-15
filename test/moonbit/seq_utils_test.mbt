///|
/// SeqUtils module tests.

///|
test "molecular_weight_dna" {
  let weight = @src.sequtils_molecular_weight("ATCG", "dna")
  assert_true(weight > 0.0)
}

///|
test "molecular_weight_protein" {
  let weight = @src.sequtils_molecular_weight("AlaArg", "protein")
  assert_true(weight > 0.0)
}

///|
test "gc_content" {
  let gc = @src.gc_content("ATCG")
  assert_eq(gc, 0.5)
}

///|
test "gc_content_all_gc" {
  let gc = @src.gc_content("GCGCGC")
  assert_eq(gc, 1.0)
}

///|
test "melting_temperature_short" {
  let tm = @src.melting_temperature("ATCG")
  assert_true(tm > 0.0)
}

///|
test "melting_temperature_long" {
  let tm = @src.melting_temperature("ATCGATCGATCGATCGATCG")
  assert_true(tm > 0.0)
}

///|
test "seq3" {
  let seq = @src.seq3("AR")
  assert_eq(seq, "Ala Arg")
}

///|
test "seq1" {
  let seq = @src.seq1("AlaArg")
  assert_eq(seq, "AR")
}

///|
test "is_dna" {
  assert_true(@src.is_dna("ATCG"))
  assert_false(@src.is_dna("ATCGX"))
}

///|
test "is_rna" {
  assert_true(@src.is_rna("AUCG"))
  assert_false(@src.is_rna("ATCG"))
}

///|
test "is_protein" {
  assert_true(@src.is_protein("AlaArg"))
  assert_false(@src.is_protein("AlaArgXyz"))
}

///|
test "count_kmers" {
  let counts = @src.count_kmers("AAAA", 2)
  assert_true(counts.size() > 0)
}

///|
test "shannon_entropy" {
  let entropy = @src.sequtils_shannon_entropy("AAAA")
  assert_true(entropy >= 0.0)
}

///|
test "ambiguity_codes" {
  let expanded = @src.ambiguity_codes("N")
  assert_eq(expanded, "ATCG")
}

///|
test "codon_usage" {
  let usage = @src.codon_usage("ATGATG")
  assert_true(usage.size() > 0)
}

///|
test "reverse_complement" {
  let rc = @src.reverse_complement("ATCG")
  assert_eq(rc, "CGAT")
}