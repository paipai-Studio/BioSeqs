///|
/// SingleR Demo - Cell type annotation by reference expression profiles.
/// Inspired by Bioconductor SingleR package.

///|
fn main {
  println("=== SingleR Cell Type Annotation Demo ===")
  println("Inspired by Bioconductor SingleR package")
  println("")
  
  println("1. Creating reference dataset with known cell type profiles...")
  let ref_data = @bio.single_r_create_reference_data()
  println("   Reference profiles: \{ref_data.n_profiles.to_string()}")
  println("   Cell types:")
  for ct in ref_data.cell_types {
    println("     - \{ct}")
  }
  println("   Genes per profile: \{ref_data.n_genes.to_string()}")
  println("")
  
  println("2. Creating synthetic single-cell data to annotate...")
  let test_data = @bio.single_r_create_test_data()
  println("   Cells: \{test_data.n_cells.to_string()}")
  println("   Genes: \{test_data.n_genes.to_string()}")
  println("")
  
  println("3. Running SingleR annotation (Spearman correlation)...")
  let params = @bio.SingleRParams::new()
  let result = @bio.single_r_annotate_cells(test_data, ref_data, params)
  println("   Annotation completed for \{result.cell_ids.length().to_string()} cells")
  println("")
  
  println("4. Annotation results summary:")
  println("   First 10 annotations:")
  let mut i = 0
  while i < 10 && i < result.cell_ids.length() {
    println("     Cell \{result.cell_ids[i]}: \{result.labels[i]} (score: \{result.scores[i].to_string()})")
    i = i + 1
  }
  println("")
  
  println("5. Summary statistics:")
  let summary = @bio.single_r_annotation_summary(result)
  let keys = summary.keys()
  for key in keys {
    println("   \{key}: \{summary.get(key).unwrap_or(0.0).to_string()}")
  }
  println("")
  
  println("6. Running SingleR with Pearson correlation...")
  let params_pearson = @bio.SingleRParams::with_method("pearson")
  let result_pearson = @bio.single_r_annotate_cells(test_data, ref_data, params_pearson)
  println("   Pearson annotation completed")
  println("   First 5 annotations (Pearson):")
  let mut j = 0
  while j < 5 && j < result_pearson.cell_ids.length() {
    println("     Cell \{result_pearson.cell_ids[j]}: \{result_pearson.labels[j]} (score: \{result_pearson.scores[j].to_string()})")
    j = j + 1
  }
  println("")
  
  println("7. Correlation analysis example:")
  let t_cell_expr = ref_data.profiles[0].expression
  let b_cell_expr = ref_data.profiles[1].expression
  let t_vs_b_spearman = @bio.single_r_spearman_correlation(t_cell_expr, b_cell_expr)
  let t_vs_self_spearman = @bio.single_r_spearman_correlation(t_cell_expr, t_cell_expr)
  let t_vs_b_pearson = @bio.single_r_pearson_correlation(t_cell_expr, b_cell_expr)
  println("   T cell vs B cell (Spearman): \{t_vs_b_spearman.to_string()}")
  println("   T cell vs T cell (Spearman): \{t_vs_self_spearman.to_string()}")
  println("   T cell vs B cell (Pearson): \{t_vs_b_pearson.to_string()}")
  println("")
  
  println("=== Demo completed successfully! ===")
  println("")
  println("Key takeaways:")
  println("- SingleR annotates cells by comparing expression profiles to reference datasets")
  println("- Spearman correlation is robust to outliers; Pearson is sensitive to magnitude")
  println("- Fine-tuning can resolve ambiguous assignments by comparing against related cell types")
  println("- The delta score (difference between top and 2nd best) indicates annotation confidence")
  println("- SingleR can be used with any reference dataset (bulk RNA-seq, sorted populations)")
}