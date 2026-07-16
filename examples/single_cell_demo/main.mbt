///|
fn main {
  println("=== Single Cell Analysis Demo ===")

  let counts = [
    [100.0, 200.0, 50.0, 150.0, 300.0],
    [120.0, 180.0, 60.0, 140.0, 280.0],
    [80.0, 220.0, 40.0, 160.0, 320.0],
    [90.0, 210.0, 45.0, 155.0, 310.0],
    [110.0, 190.0, 55.0, 145.0, 290.0],
  ]

  let gene_names = ["GeneA", "MT-CO1", "GeneB", "MT-ATP6", "GeneC"]

  println("\n1. Calculate QC Metrics:")
  let qc = @bio.calculate_qc_metrics(counts, gene_names)

  let mut i = 0
  while i < qc.n_umi.length() {
    println(
      "Cell " +
      (i + 1).to_string() +
      ": nUMI=" +
      qc.n_umi[i].to_string() +
      ", nGenes=" +
      qc.n_genes[i].to_string() +
      ", mito%=" +
      qc.mito_percent[i].to_string(),
    )
    i = i + 1
  }

  println("\n2. Filter Cells:")
  let filtered = @bio.filter_cells(counts, qc, 500, 1000, 3, 20.0)
  println("Filtered cells: " + filtered.length().to_string())
  
  let working_data = if filtered.length() > 0 { filtered } else { counts }
  
  println("\n3. Log Normalization:")
  let normalized = @bio.log_normalize(working_data, None)
  println(
    "Normalized matrix shape: " +
    normalized.length().to_string() +
    "x" +
    normalized[0].length().to_string(),
  )

  println("\n4. Find Variable Genes:")
  let top_genes = @bio.find_variable_genes(normalized, 3)
  println("Top variable gene indices: " + top_genes.length().to_string())

  let mut j = 0
  while j < top_genes.length() {
    println(top_genes[j].to_string() + " ")
    j = j + 1
  }
  println("")

  println("\n5. PCA Analysis:")
  let (pca_result, eigenvalues) = @bio.pca_power_iteration(normalized, 2)

  println("Eigenvalues:")
  j = 0
  while j < eigenvalues.length() {
    println("PC" + (j + 1).to_string() + ": " + eigenvalues[j].to_string())
    j = j + 1
  }

  println("\nPCA coordinates:")
  j = 0
  while j < pca_result.length() {
    println("Cell " + (j + 1).to_string() + ": [")
    let mut k = 0
    while k < pca_result[j].length() {
      println(pca_result[j][k].to_string())
      if k < pca_result[j].length() - 1 {
        println(", ")
      }
      k = k + 1
    }
    println("]")
    j = j + 1
  }

  println("\n6. SingleCellExperiment:")
  let assays = [normalized]
  let col_data = [["Cell1"], ["Cell2"], ["Cell3"], ["Cell4"], ["Cell5"]]
  let row_data = [["GeneA"], ["MT-CO1"], ["GeneB"], ["MT-ATP6"], ["GeneC"]]

  let sce = @bio.SingleCellExperiment::new(assays, col_data, row_data)
  println("Assays: " + sce.assays.length().to_string())
  println("Cells: " + sce.col_data.length().to_string())
  println("Genes: " + sce.row_data.length().to_string())

  println("\n=== Demo Complete ===")
}
