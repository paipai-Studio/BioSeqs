///|
/// Tests for microbiome module.

// ============================================================
// Alpha Diversity Tests
// ============================================================

test "calc_observed basic" {
  let counts : Array[Double] = [10.0, 20.0, 0.0, 5.0, 0.0, 15.0]
  let observed = @src.calc_observed(counts)
  assert_eq(observed, 4.0)
}

test "calc_observed all zeros" {
  let counts : Array[Double] = [0.0, 0.0, 0.0]
  let observed = @src.calc_observed(counts)
  assert_eq(observed, 0.0)
}

test "calc_shannon even" {
  let counts : Array[Double] = [10.0, 10.0, 10.0, 10.0]
  let shannon = @src.calc_shannon(counts)
  // For 4 equally abundant species: H = ln(4) ≈ 1.386
  assert_true(shannon > 1.3 && shannon < 1.5)
}

test "calc_shannon single taxon" {
  let counts : Array[Double] = [100.0]
  let shannon = @src.calc_shannon(counts)
  assert_eq(shannon, 0.0)
}

test "calc_simpson even" {
  let counts : Array[Double] = [10.0, 10.0, 10.0, 10.0]
  let simpson = @src.calc_simpson(counts)
  // For 4 equally abundant species: D = 1 - 4*(0.25^2) = 0.75
  assert_true(simpson > 0.7 && simpson < 0.8)
}

test "calc_simpson single taxon" {
  let counts : Array[Double] = [100.0]
  let simpson = @src.calc_simpson(counts)
  assert_eq(simpson, 0.0)
}

test "calc_inv_simpson even" {
  let counts : Array[Double] = [10.0, 10.0, 10.0, 10.0]
  let inv_simpson = @src.calc_inv_simpson(counts)
  // For 4 equally abundant species: 1/D = 4
  assert_true(inv_simpson > 3.5 && inv_simpson < 4.5)
}

test "calc_pielou_evenness even" {
  let counts : Array[Double] = [10.0, 10.0, 10.0, 10.0]
  let pielou = @src.calc_pielou_evenness(counts)
  // Even distribution should have Pielou's J close to 1
  assert_true(pielou > 0.95)
}

test "calc_pielou_evenness single taxon" {
  let counts : Array[Double] = [100.0]
  let pielou = @src.calc_pielou_evenness(counts)
  assert_eq(pielou, 0.0)
}

test "calc_chao1 basic" {
  let counts : Array[Double] = [10.0, 5.0, 3.0, 1.0, 1.0, 2.0, 2.0]
  let chao1 = @src.calc_chao1(counts)
  // Observed = 7, f1 = 2, f2 = 2
  // Chao1 = 7 + (2 * 1) / (2 * 3) = 7 + 2/6 ≈ 7.33
  assert_true(chao1 >= 7.0)
}

test "calc_fisher_alpha basic" {
  let counts : Array[Double] = [100.0, 50.0, 30.0, 20.0, 10.0]
  let alpha = @src.calc_fisher_alpha(counts)
  assert_true(alpha > 0.0)
}

test "calc_alpha_diversity all indices" {
  let counts : Array[Double] = [10.0, 20.0, 15.0, 5.0, 25.0]
  let div = @src.calc_alpha_diversity(counts)
  assert_eq(div.observed, 5.0)
  assert_true(div.shannon > 0.0)
  assert_true(div.simpson > 0.0)
  assert_true(div.inv_simpson > 0.0)
  assert_true(div.chao1 >= div.observed)
}

test "calc_alpha_diversity_table multiple samples" {
  let otu_table = @src.create_example_otu_table()
  let div_table = @src.calc_alpha_diversity_table(otu_table)
  assert_eq(div_table.length(), 6)
  assert_true(div_table[0].observed > 0.0)
  assert_true(div_table[5].shannon > 0.0)
}

// ============================================================
// Beta Diversity Tests
// ============================================================

test "bray_curtis identical" {
  let x : Array[Double] = [10.0, 20.0, 30.0]
  let y : Array[Double] = [10.0, 20.0, 30.0]
  let bc = @src.bray_curtis(x, y)
  assert_eq(bc, 0.0)
}

test "bray_curtis completely different" {
  let x : Array[Double] = [10.0, 0.0, 0.0]
  let y : Array[Double] = [0.0, 20.0, 30.0]
  let bc = @src.bray_curtis(x, y)
  // BC = (10 + 20 + 30) / (10 + 20 + 30) = 1.0
  assert_eq(bc, 1.0)
}

test "bray_curtis partial overlap" {
  let x : Array[Double] = [10.0, 20.0, 0.0]
  let y : Array[Double] = [10.0, 0.0, 20.0]
  let bc = @src.bray_curtis(x, y)
  assert_true(bc > 0.0 && bc < 1.0)
}

test "jaccard_distance identical" {
  let x : Array[Double] = [10.0, 20.0, 30.0]
  let y : Array[Double] = [5.0, 10.0, 15.0]
  let jaccard = @src.jaccard_distance(x, y)
  assert_eq(jaccard, 0.0)
}

test "jaccard_distance no overlap" {
  let x : Array[Double] = [10.0, 0.0, 0.0]
  let y : Array[Double] = [0.0, 20.0, 30.0]
  let jaccard = @src.jaccard_distance(x, y)
  assert_eq(jaccard, 1.0)
}

test "jaccard_distance partial" {
  let x : Array[Double] = [10.0, 20.0, 0.0]
  let y : Array[Double] = [10.0, 0.0, 20.0]
  let jaccard = @src.jaccard_distance(x, y)
  // Intersection = 1, Union = 3, J = 1 - 1/3 ≈ 0.667
  assert_true(jaccard > 0.6 && jaccard < 0.7)
}

test "jensen_shannon_divergence identical" {
  let x : Array[Double] = [10.0, 20.0, 30.0]
  let y : Array[Double] = [10.0, 20.0, 30.0]
  let jsd = @src.jensen_shannon_divergence(x, y)
  assert_eq(jsd, 0.0)
}

test "jensen_shannon_divergence symmetric" {
  let x : Array[Double] = [10.0, 20.0, 30.0]
  let y : Array[Double] = [30.0, 20.0, 10.0]
  let jsd_xy = @src.jensen_shannon_divergence(x, y)
  let jsd_yx = @src.jensen_shannon_divergence(y, x)
  assert_true((jsd_xy - jsd_yx).abs() < 0.0001)
}

test "calc_beta_diversity_matrix bray" {
  let otu_table = @src.create_example_otu_table()
  let dist_matrix = @src.calc_beta_diversity_matrix(otu_table, "bray")
  assert_eq(dist_matrix.length(), 6)
  assert_eq(dist_matrix[0].length(), 6)
  // Diagonal should be 0
  assert_eq(dist_matrix[0][0], 0.0)
  assert_eq(dist_matrix[3][3], 0.0)
  // Symmetric
  assert_true((dist_matrix[0][1] - dist_matrix[1][0]).abs() < 0.0001)
}

test "calc_beta_diversity_matrix jaccard" {
  let otu_table = @src.create_example_otu_table()
  let dist_matrix = @src.calc_beta_diversity_matrix(otu_table, "jaccard")
  assert_eq(dist_matrix.length(), 6)
  assert_true(dist_matrix[0][3] > 0.0)
}

test "calc_beta_diversity_matrix jsd" {
  let otu_table = @src.create_example_otu_table()
  let dist_matrix = @src.calc_beta_diversity_matrix(otu_table, "jsd")
  assert_eq(dist_matrix.length(), 6)
  assert_true(dist_matrix[0][1] >= 0.0)
}

// ============================================================
// PCoA Tests
// ============================================================

test "pcoa basic" {
  let otu_table = @src.create_example_otu_table()
  let dist_matrix = @src.calc_beta_diversity_matrix(otu_table, "bray")
  let pcoa_result = @src.pcoa(dist_matrix)
  
  assert_true(pcoa_result.eigenvalues.length() > 0)
  assert_true(pcoa_result.vectors.length() == 6)
  assert_true(pcoa_result.variance_explained.length() > 0)
  assert_true(pcoa_result.cumulative_variance.length() > 0)
}

test "pcoa eigenvalues positive" {
  let otu_table = @src.create_example_otu_table()
  let dist_matrix = @src.calc_beta_diversity_matrix(otu_table, "bray")
  let pcoa_result = @src.pcoa(dist_matrix)
  
  let mut i = 0
  while i < pcoa_result.eigenvalues.length() {
    assert_true(pcoa_result.eigenvalues[i] >= 0.0)
    i = i + 1
  }
}

test "pcoa has vectors" {
  let otu_table = @src.create_example_otu_table()
  let dist_matrix = @src.calc_beta_diversity_matrix(otu_table, "bray")
  let pcoa_result = @src.pcoa(dist_matrix)
  
  assert_true(pcoa_result.vectors.length() == 6)
  if pcoa_result.vectors.length() > 0 {
    assert_true(pcoa_result.vectors[0].length() > 0)
  }
}

test "pcoa single sample" {
  let dist_matrix : Array[Array[Double]] = [[0.0]]
  let pcoa_result = @src.pcoa(dist_matrix)
  assert_eq(pcoa_result.eigenvalues.length(), 0)
}

// ============================================================
// Differential Abundance Tests
// ============================================================

test "differential_abundance basic" {
  let otu_table = @src.create_example_otu_table()
  let taxa_names = @src.get_example_taxa_names()
  let group1 : Array[Int] = [0, 1, 2]
  let group2 : Array[Int] = [3, 4, 5]
  
  let results = @src.differential_abundance(otu_table, taxa_names, group1, group2)
  
  assert_eq(results.length(), 10)
  assert_eq(results[0].taxa, "Bacteroides")
  assert_true(results[0].p_value >= 0.0 && results[0].p_value <= 1.0)
  assert_true(results[0].adjusted_p_value >= results[0].p_value)
}

test "differential_abundance log2_fc direction" {
  let otu_table = @src.create_example_otu_table()
  let taxa_names = @src.get_example_taxa_names()
  let group1 : Array[Int] = [0, 1, 2]  // control
  let group2 : Array[Int] = [3, 4, 5]  // treatment
  
  let results = @src.differential_abundance(otu_table, taxa_names, group1, group2)
  
  // Taxon 0 (Bacteroides) is high in control, so log2FC should be negative
  assert_true(results[0].log2_fold_change < 0.0)
  // Taxon 1 (Prevotella) is high in treatment, so log2FC should be positive
  assert_true(results[1].log2_fold_change > 0.0)
}

test "differential_abundance empty input" {
  let otu_table : Array[Array[Double]] = []
  let taxa_names : Array[String] = []
  let group1 : Array[Int] = [0, 1]
  let group2 : Array[Int] = [2, 3]
  
  let results = @src.differential_abundance(otu_table, taxa_names, group1, group2)
  assert_eq(results.length(), 0)
}

// ============================================================
// Helper function tests
// ============================================================

test "create_example_otu_table dimensions" {
  let otu_table = @src.create_example_otu_table()
  assert_eq(otu_table.length(), 10)  // 10 taxa
  assert_eq(otu_table[0].length(), 6)  // 6 samples
}

test "get_example_taxa_names" {
  let names = @src.get_example_taxa_names()
  assert_eq(names.length(), 10)
  assert_eq(names[0], "Bacteroides")
  assert_eq(names[9], "Escherichia")
}
