///|
/// Population genetics demo.

fn main {
  println("=== Population Genetics Demo ===\n")
  
  let population = ["AA", "Aa", "Aa", "aa", "AA", "Aa", "aa"]
  println("1. Population data:")
  println("   Genotypes: \{population}")
  
  println("\n2. Allele frequencies:")
  let allele_freqs = @src.allele_frequencies(population)
  println("   Alleles: \{allele_freqs.size()}")
  
  println("\n3. Genotype frequencies:")
  let geno_freqs = @src.genotype_frequencies(population)
  println("   Genotypes: \{geno_freqs.size()}")
  
  println("\n4. Hardy-Weinberg test:")
  let chi_sq = @src.hardy_weinberg_test(population)
  println("   Chi-squared: \{chi_sq}")
  
  println("\n5. Heterozygosity:")
  let h_exp = @src.heterozygosity(population)
  let h_obs = @src.observed_heterozygosity(population)
  println("   Expected heterozygosity: \{h_exp}")
  println("   Observed heterozygosity: \{h_obs}")
  
  println("\n6. Inbreeding coefficient:")
  let f = @src.inbreeding_coefficient(population)
  println("   F: \{f}")
  
  println("\n7. Allele richness:")
  let richness = @src.allele_richness(population)
  println("   Richness: \{richness}")
  
  println("\n8. FST between populations:")
  let pop1 = ["AA", "Aa", "aa", "AA"]
  let pop2 = ["AA", "AA", "Aa", "Aa"]
  let fst_val = @src.fst(pop1, pop2)
  println("   FST: \{fst_val}")
  
  println("\n9. Watterson's theta:")
  let sequences = ["ATCG", "ATCA", "ATGG", "ATGT"]
  let theta = @src.watterson_theta(sequences)
  println("   Theta: \{theta}")
  
  println("\n=== Demo Complete ===")
}