///|
/// Population genetics demo.
///
/// Run: moon run examples/popgen_demo/main.mbt
fn main {
  println("=== Population Genetics Demo ===\n")

  let population = ["AA", "Aa", "Aa", "aa", "AA", "Aa", "aa"]
  println("1. Population data:")
  println("   Genotypes: \{population}")

  println("\n2. Allele frequencies:")
  let allele_freqs = @bio.allele_frequencies(population)
  println("   Number of alleles: \{allele_freqs.size()}")

  println("\n3. Genotype frequencies:")
  let geno_freqs = @bio.genotype_frequencies(population)
  println("   Number of genotypes: \{geno_freqs.size()}")

  println("\n4. Hardy-Weinberg test:")
  let chi_sq = @bio.hardy_weinberg_test(population)
  println("   Chi-squared: \{chi_sq}")

  println("\n5. Heterozygosity:")
  let h_exp = @bio.heterozygosity(population)
  let h_obs = @bio.observed_heterozygosity(population)
  println("   Expected heterozygosity: \{h_exp}")
  println("   Observed heterozygosity: \{h_obs}")

  println("\n6. Inbreeding coefficient:")
  let f = @bio.inbreeding_coefficient(population)
  println("   F: \{f}")

  println("\n7. Allele richness:")
  let richness = @bio.allele_richness(population)
  println("   Richness: \{richness}")

  println("\n8. FST between populations:")
  let pop1 = ["AA", "Aa", "aa", "AA"]
  let pop2 = ["AA", "AA", "Aa", "Aa"]
  let fst_val = @bio.fst(pop1, pop2)
  println("   FST: \{fst_val}")

  println("\n9. Watterson's theta:")
  let sequences = ["ATCG", "ATCA", "ATGG", "ATGT"]
  let theta = @bio.watterson_theta(sequences)
  println("   Theta: \{theta}")

  println("\n10. Diversity indices:")
  let counts = [100, 100, 100, 100, 50, 50, 30, 20, 10, 10, 10]
  println("   Species counts: \{counts}")
  println("   Shannon index: \{@bio.shannon_index(counts)}")
  println("   Simpson index: \{@bio.simpson_index(counts)}")
  println("   Simpson reciprocal: \{@bio.simpson_reciprocal(counts)}")
  println("   Evenness: \{@bio.evenness(counts)}")
  println("   Chao1 estimate: \{@bio.chao1(counts)}")
  println("   ACE estimate: \{@bio.ace(counts)}")
  println("   Fisher's alpha: \{@bio.fishers_alpha(counts)}")
  println("   Margalef index: \{@bio.margalef_index(counts)}")
  println("   Menhinick index: \{@bio.menhinick_index(counts)}")

  println("\n=== Demo Complete ===")
}
