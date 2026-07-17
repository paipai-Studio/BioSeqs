///|
fn main {
  println("=== DESeq2 Differential Expression Analysis Demo ===")
  println("")

  let counts = [
    [100, 120, 200, 220],
    [50, 60, 55, 45],
    [10, 15, 80, 90],
    [200, 220, 190, 210],
    [30, 35, 150, 160],
    [80, 90, 75, 85],
    [5, 8, 100, 110],
    [150, 160, 140, 155],
  ]

  let row_names = [
    "GeneA", "GeneB", "GeneC", "GeneD", "GeneE", "GeneF", "GeneG", "GeneH",
  ]
  let col_names = ["Control1", "Control2", "Treat1", "Treat2"]
  let design = [[1.0, 0.0], [1.0, 0.0], [1.0, 1.0], [1.0, 1.0]]

  println("1. Creating DESeqDataSet...")
  println("   Genes: \{row_names.length()}")
  println("   Samples: \{col_names.length()}")
  println("")

  let dds = @bio.deseq_dataset(counts, row_names, col_names, design)

  println("2. Running DESeq2 analysis...")
  let results = @bio.results(dds)
  println("   Analysis completed!")
  println("")

  println("3. Results Summary:")
  println("   Gene\tbaseMean\tlog2FC")
  let mut i = 0
  while i < results.row_names.length() {
    println(
      "   \{results.row_names[i]}\t\{results.base_mean[i]}\t\{results.log2_fold_change[i]}",
    )
    i = i + 1
  }
  println("")

  println("4. Significant Genes (p_value < 0.05, |log2FC| > 1):")
  let sig_genes = @bio.significant_genes(results, 0.05, 1.0)
  if sig_genes.length() > 0 {
    let mut i = 0
    while i < sig_genes.length() {
      println("   - \\{sig_genes[i]}")
      i = i + 1
    }
  } else {
    println("   (No significant genes found)")
  }
  println("")

  println("=== Analysis Complete ===")
}
