// Tests for Bio.FreqAnalysis module

test "fa_count_kmers single nucleotide" {
  let result = @src.fa_count_kmers("ACGT", 1)
  assert_true(result.fa_get_total_count() == 4)
  assert_true(result.fa_get_frequency("A") > 0.24 && result.fa_get_frequency("A") < 0.26)
}

test "fa_count_kmers dinucleotide" {
  let result = @src.fa_count_kmers("ACGTAC", 2)
  assert_true(result.fa_get_total_count() == 5)
}

test "fa_count_kmers empty" {
  let result = @src.fa_count_kmers("", 3)
  assert_true(result.fa_get_total_count() == 0)
}

test "fa_nucleotide_frequency" {
  let result = @src.fa_nucleotide_frequency("AATTCCGG")
  assert_true(result.fa_get_total_count() == 8)
  assert_true(result.fa_get_frequency("A") > 0.24 && result.fa_get_frequency("A") < 0.26)
}

test "fa_count_pattern" {
  let count = @src.fa_count_pattern("ATATAT", "AT")
  assert_true(count == 3)
}

test "fa_count_pattern no match" {
  let count = @src.fa_count_pattern("ACGT", "TT")
  assert_true(count == 0)
}

test "fa_codon_usage" {
  let result = @src.fa_codon_usage("ATGAAACCC")
  assert_true(result.fa_get_total_count() == 3)
}

test "fa_gc_content" {
  let result = @src.fa_gc_content("AAGGCC")
  assert_true(result > 0.66 && result < 0.68)
}

test "fa_gc_content all GC" {
  let result = @src.fa_gc_content("GCGC")
  assert_true(result == 1.0)
}

test "fa_at_content" {
  let result = @src.fa_at_content("AATTCC")
  assert_true(result > 0.66 && result < 0.68)
}

test "fa_find_motif" {
  let positions = @src.fa_find_motif("ATCGATCGATCG", "ATC")
  assert_true(positions.length() == 3)
  assert_true(positions[0] == 0)
  assert_true(positions[1] == 4)
  assert_true(positions[2] == 8)
}

test "fa_find_motif not found" {
  let positions = @src.fa_find_motif("ACGT", "TT")
  assert_true(positions.length() == 0)
}

test "fa_get_patterns" {
  let result = @src.fa_count_kmers("AAACCC", 1)
  let patterns = result.fa_get_patterns()
  assert_true(patterns.length() >= 2)
}
