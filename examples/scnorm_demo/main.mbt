///|
/// SCnorm Single-cell RNA-seq Normalization Demo.
fn main {
  println("=== SCnorm Single-cell RNA-seq Normalization Demo ===\n")

  // Create simulated single-cell data
  println("1. Creating simulated scRNA-seq count matrix...")
  let genes = ["Gene_A", "Gene_B", "Gene_C", "Gene_D", "Gene_E"]
  let cell_ids = ["cell_1", "cell_2", "cell_3"]
  // Counts: genes (rows) x cells (columns)
  let counts : Array[Array[Double]] = [
    [10.0, 20.0, 30.0],  // Gene_A
    [5.0, 15.0, 25.0],   // Gene_B
    [100.0, 200.0, 300.0], // Gene_C (highly expressed)
    [50.0, 80.0, 120.0],  // Gene_D
    [1.0, 2.0, 3.0]       // Gene_E
  ]
  println("   Genes: \{genes.length()}, Cells: \{cell_ids.length()}")

  // Compute library sizes
  println("\n2. Computing library sizes per cell...")
  let lib_sizes = @src.sc_norm_library_sizes(counts)
  let mut s = 0
  while s < lib_sizes.length() {
    let ls = lib_sizes[s].to_string()
    println("   \{cell_ids[s]}: library_size=\{ls}")
    s = s + 1
  }

  // Normalize a single gene
  println("\n3. Normalizing Gene_C with SCnorm...")
  let params = @src.SCnormParams::new()
  let gene_id = "Gene_C"
  let gene_counts = counts[2]
  let norm_result = @src.sc_norm_normalize_gene(gene_id, gene_counts, lib_sizes, params)
  println("   Gene: \{norm_result.gene_id}")
  println("   Passed filter: \{norm_result.passed_filter}")
  println("   Normalized counts:")
  let mut n = 0
  while n < gene_counts.length() {
    let orig = gene_counts[n].to_string()
    let norm = norm_result.norm_counts[n].to_string()
    println("     \{cell_ids[n]}: original=\{orig}, normalized=\{norm}")
    n = n + 1
  }

  // Run full normalization pipeline
  println("\n4. Running full SCnorm normalization pipeline...")
  let full_result = @src.sc_norm_run(counts, genes, cell_ids, params)
  println("   Gene results: \{full_result.gene_results.length()}")
  println("   Filtered genes: \{full_result.gene_ids.length()}")

  // Print some gene results
  println("\n5. Normalization results for each gene:")
  let mut g = 0
  while g < full_result.gene_results.length() && g < 3 {
    let fit = full_result.gene_results[g]
    let slope_fit = fit.quant_fit.slope.to_string()
    let intercept_fit = fit.quant_fit.intercept.to_string()
    println("   \{fit.gene_id}: slope=\{slope_fit}, intercept=\{intercept_fit}")
    g = g + 1
  }

  println("\n=== SCnorm Demo Complete ===")
}
