///|
pub fn main() {
  println("=== monocle3 Demo ===")
  
  let (counts, gene_names, cell_names) = @src.create_example_monocle_data()
  
  println("\n1. Create CellDataSet:")
  let cds = @src.new_cell_data_set(counts, gene_names, cell_names)
  println("  Cells: " + cds.counts.length().to_string())
  println("  Genes: " + cds.gene_names.length().to_string())
  
  println("\n2. Preprocess CDS (normalization + PCA):")
  let preprocessed = @src.preprocess_cds(cds, 10)
  println("  Normalized: " + preprocessed.normalized.length().to_string() + " cells")
  println("  PCA dimensions: " + preprocessed.reduced_dimensions.length().to_string())
  
  println("\n3. Reduce dimension with UMAP:")
  let reduced = @src.reduce_dimension(preprocessed, "UMAP", 2)
  println("  Reduced dimensions: " + reduced.reduced_dimensions.length().to_string())
  
  println("\n4. Learn principal graph:")
  let with_graph = @src.learn_graph(reduced)
  println("  Cell partitions: " + with_graph.cell_partitions.length().to_string())
  println("  Principal graph nodes: " + with_graph.principal_graph.length().to_string())
  
  println("\n5. Order cells along trajectory:")
  let ordered = @src.order_cells(with_graph, [0])
  println("  Pseudotime assigned: " + ordered.pseudotime.length().to_string())
  
  let mut min_pt = ordered.pseudotime[0]
  let mut max_pt = ordered.pseudotime[0]
  let mut i = 0
  while i < ordered.pseudotime.length() {
    if ordered.pseudotime[i] < min_pt {
      min_pt = ordered.pseudotime[i]
    }
    if ordered.pseudotime[i] > max_pt {
      max_pt = ordered.pseudotime[i]
    }
    i = i + 1
  }
  println("  Pseudotime range: " + min_pt.to_string() + " - " + max_pt.to_string())
  
  println("\n6. Fit models to find pseudotime-dependent genes:")
  let test_genes = ["Gene_0", "Gene_1", "Gene_2", "Gene_20", "Gene_21", "Gene_40", "Gene_41"]
  let result = @src.fit_models(ordered, test_genes)
  println("  Tested genes: " + result.gene.length().to_string())
  
  println("\n7. Results summary:")
  i = 0
  while i < result.gene.length() {
    println("    " + result.gene[i] + ": slope=" + result.trend_statistic[i].to_string() + ", q-value=" + result.q_value[i].to_string())
    i = i + 1
  }
  
  println("\n8. Identify trajectory marker genes (q-value < 0.05):")
  let mut markers : Array[String] = Array::new()
  i = 0
  while i < result.gene.length() {
    if result.q_value[i] < 0.05 {
      markers.push(result.gene[i])
    }
    i = i + 1
  }
  
  println("  Trajectory markers: " + markers.length().to_string())
  if markers.length() > 0 {
    let mut j = 0
    while j < markers.length() {
      print(markers[j])
      if j < markers.length() - 1 {
        print(", ")
      }
      j = j + 1
    }
    println("")
  }
  
  println("\nmonocle3 demo completed!")
}