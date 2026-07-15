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