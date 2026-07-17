///|
/// Population genetics module tests.

///|
test "allele_frequencies" {
  let pop = ["AA", "Aa", "aa"]
  let freqs = @src.allele_frequencies(pop)
  assert_true(freqs.size() > 0)
}

///|
test "genotype_frequencies" {
  let pop = ["AA", "Aa", "aa"]
  let freqs = @src.genotype_frequencies(pop)
  assert_true(freqs.size() > 0)
}

///|
test "hardy_weinberg_test" {
  let pop = ["AA", "Aa", "aa"]
  let chi_sq = @src.hardy_weinberg_test(pop)
  assert_true(chi_sq >= 0.0)
}

///|
test "heterozygosity" {
  let pop = ["AA", "Aa", "aa"]
  let het = @src.heterozygosity(pop)
  assert_true(het >= 0.0)
}

///|
test "observed_heterozygosity" {
  let pop = ["AA", "Aa", "aa"]
  let het = @src.observed_heterozygosity(pop)
  assert_true(het >= 0.0)
}

///|
test "inbreeding_coefficient" {
  let pop = ["AA", "Aa", "aa"]
  let f = @src.inbreeding_coefficient(pop)
  assert_true(f >= 0.0)
}

///|
test "allele_richness" {
  let pop = ["AA", "Aa", "aa"]
  let richness = @src.allele_richness(pop)
  assert_eq(richness, 2)
}

///|
test "watterson_theta" {
  let pop = ["AAAA", "AAAT"]
  let theta = @src.watterson_theta(pop)
  assert_true(theta >= 0.0)
}

///|
test "fst" {
  let pop1 = ["AA", "Aa", "aa"]
  let pop2 = ["AA", "AA", "Aa"]
  let f = @src.fst(pop1, pop2)
  assert_true(f >= 0.0)
}

///|
test "fst_multiple" {
  let pops = [["AA", "Aa"], ["AA", "AA"]]
  let f = @src.fst_multiple(pops)
  assert_true(f >= 0.0)
}

///|
test "pairwise_fst" {
  let pops = [["AA", "Aa"], ["AA", "AA"], ["Aa", "aa"]]
  let matrix = @src.pairwise_fst(pops)
  assert_eq(matrix.length(), 3)
}

///|
test "nucleotide_diversity" {
  let pop = ["AAAA", "AAAT", "AATT"]
  let pi = @src.nucleotide_diversity(pop)
  assert_true(pi >= 0.0)
}

///|
test "migrate_rate" {
  let m = @src.migrate_rate(0.1, 1000, 100)
  assert_true(m > 0.0)
}

///|
test "shannon_index" {
  let counts = [10, 10, 10]
  let h = @src.shannon_index(counts)
  assert_true(h > 0.0)
  assert_true(h <= 2.0)
}

///|
test "shannon_index_equal_abundance" {
  let counts = [100, 100, 100, 100]
  let h = @src.shannon_index(counts)
  assert_true((h - 2.0).abs() < 0.01)
}

///|
test "shannon_index_single_species" {
  let counts = [100]
  let h = @src.shannon_index(counts)
  assert_true(h == 0.0)
}

///|
test "simpson_index" {
  let counts = [10, 10, 10]
  let d = @src.simpson_index(counts)
  assert_true(d > 0.0)
  assert_true(d <= 1.0)
}

///|
test "simpson_index_equal_abundance" {
  let counts = [100, 100, 100]
  let d = @src.simpson_index(counts)
  assert_true(d > 0.5)
}

///|
test "simpson_reciprocal" {
  let counts = [10, 10, 10]
  let d = @src.simpson_reciprocal(counts)
  assert_true(d > 0.0)
}

///|
test "evenness" {
  let counts = [10, 10, 10]
  let e = @src.evenness(counts)
  assert_true(e >= 0.0)
  assert_true(e <= 1.0)
}

///|
test "evenness_max" {
  let counts = [100, 100, 100, 100]
  let e = @src.evenness(counts)
  assert_true((e - 1.0).abs() < 0.01)
}

///|
test "chao1" {
  let counts = [5, 4, 3, 2, 1, 1, 1]
  let est = @src.chao1(counts)
  assert_true(est >= counts.filter(fn(c) { c > 0 }).length().to_double())
}

///|
test "chao1_no_singletons" {
  let counts = [5, 4, 3, 2]
  let est = @src.chao1(counts)
  assert_true(est == 4.0)
}

///|
test "ace" {
  let counts = [10, 10, 5, 5, 3, 2, 1, 1, 1]
  let est = @src.ace(counts)
  assert_true(est >= 0.0)
}

///|
test "fishers_alpha" {
  let counts = [10, 10, 10, 5, 5, 5]
  let alpha = @src.fishers_alpha(counts)
  assert_true(alpha > 0.0)
}

///|
test "margalef_index" {
  let counts = [10, 10, 10, 10]
  let d = @src.margalef_index(counts)
  assert_true(d > 0.0)
}

///|
test "menhinick_index" {
  let counts = [10, 10, 10]
  let d = @src.menhinick_index(counts)
  assert_true(d > 0.0)
}

///|
test "diversity_empty_counts" {
  let empty = []
  assert_true(@src.shannon_index(empty) == 0.0)
  assert_true(@src.simpson_index(empty) == 0.0)
  assert_true(@src.chao1(empty) == 0.0)
}
