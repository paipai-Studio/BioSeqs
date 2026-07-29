///|
/// Microbiome Analysis Demo
///
/// This demo shows how to use the microbiome module for
/// analyzing microbial community data including alpha/beta diversity,
/// ordination, and differential abundance analysis.
fn main {
  println("=== Microbiome Analysis Demo ===")
  println("")

  println("1. Alpha Diversity Analysis:")
  let sample_counts = [10.0, 20.0, 15.0, 5.0, 25.0, 8.0, 12.0, 5.0]
  let alpha = @src.calc_alpha_diversity(sample_counts)
  println("   Sample counts: \{sample_counts}")
  println("   Observed OTUs: \{alpha.observed}")
  println("   Shannon index: \{alpha.shannon}")
  println("   Simpson index: \{alpha.simpson}")
  println("   Chao1 index: \{alpha.chao1}")
  println("   Pielou evenness: \{alpha.pielou_evenness}")
  println("   ACE index: \{alpha.ace}")
  println("   Fisher alpha: \{alpha.fisher_alpha}")
  println("")

  println("2. Beta Diversity Analysis:")
  let sample1 = [10.0, 20.0, 30.0, 40.0]
  let sample2 = [15.0, 25.0, 25.0, 35.0]
  let sample3 = [40.0, 30.0, 20.0, 10.0]
  println("   Sample 1: \{sample1}")
  println("   Sample 2: \{sample2}")
  println("   Sample 3: \{sample3}")
  println("   Bray-Curtis (1 vs 2): \{@src.bray_curtis(sample1, sample2)}")
  println("   Bray-Curtis (1 vs 3): \{@src.bray_curtis(sample1, sample3)}")
  println("   Jaccard distance (1 vs 2): \{@src.jaccard_distance(sample1, sample2)}")
  println("   Jensen-Shannon divergence (1 vs 2): \{@src.jensen_shannon_divergence(sample1, sample2)}")
  println("")

  println("3. Beta Diversity Matrix:")
  let otu_table_small = [
    [10.0, 20.0, 30.0, 40.0],
    [15.0, 25.0, 25.0, 35.0],
    [40.0, 30.0, 20.0, 10.0],
  ]
  let beta_matrix = @src.calc_beta_diversity_matrix(otu_table_small, "bray_curtis")
  println("   Distance matrix (Bray-Curtis):")
  let mut i = 0
  while i < beta_matrix.length() {
    println("     \{beta_matrix[i]}")
    i = i + 1
  }
  println("")

  println("4. PCoA (Principal Coordinates Analysis):")
  let otu_table = [
    [10.0, 20.0, 30.0, 40.0],
    [15.0, 25.0, 25.0, 35.0],
    [40.0, 30.0, 20.0, 10.0],
    [35.0, 32.0, 18.0, 15.0],
    [12.0, 22.0, 28.0, 38.0],
    [38.0, 29.0, 21.0, 12.0],
  ]
  let dist_for_pcoa = @src.calc_beta_diversity_matrix(otu_table, "bray_curtis")
  let pcoa_result = @src.pcoa(dist_for_pcoa)
  println("   Samples: \{otu_table.length()}")
  println("   OTUs per sample: \{otu_table[0].length()}")
  println("   Eigenvalues count: \{pcoa_result.eigenvalues.length()}")
  if pcoa_result.eigenvalues.length() > 1 {
    println("   Variance explained (PC1): \{pcoa_result.variance_explained[0]}")
    println("   Variance explained (PC2): \{pcoa_result.variance_explained[1]}")
  }
  println("")

  println("5. Differential Abundance Analysis:")
  let da_otu_table = [
    [10.0, 20.0, 30.0, 40.0],
    [12.0, 22.0, 28.0, 38.0],
    [11.0, 21.0, 29.0, 39.0],
    [13.0, 23.0, 27.0, 37.0],
    [40.0, 30.0, 20.0, 10.0],
    [38.0, 28.0, 22.0, 12.0],
    [42.0, 32.0, 18.0, 8.0],
    [39.0, 29.0, 21.0, 11.0],
  ]
  let taxa_names = ["OTU1", "OTU2", "OTU3", "OTU4"]
  let group1 = [0, 1, 2, 3]
  let group2 = [4, 5, 6, 7]
  let da_results = @src.differential_abundance(da_otu_table, taxa_names, group1, group2)
  println("   Group A samples: 4")
  println("   Group B samples: 4")
  println("   OTUs tested: \{da_results.length()}")
  println("   Results:")
  let mut k = 0
  while k < da_results.length() {
    let r = da_results[k]
    println("     \{r.taxa}: log2FC=\{r.log2_fold_change}, p=\{r.p_value}, adj_p=\{r.adjusted_p_value}")
    k = k + 1
  }
  println("")

  println("6. Example OTU table:")
  let example_otu = @src.create_example_otu_table()
  let example_taxa = @src.get_example_taxa_names()
  println("   Samples: \{example_otu.length()}")
  println("   Taxa: \{example_taxa.length()}")
  println("   Taxa names: \{example_taxa}")
  println("")

  println("=== Demo Complete ===")
}
