///|
/// EnhancedVolcano module example - Volcano plots for differential expression.
///
/// Run: moon run examples/enhanced_volcano_demo
fn main {
  println("=== EnhancedVolcano Module Demo ===")
  println("")

  // Basic volcano plot creation
  println("=== Basic Volcano Plot ===")
  let ids = ["TP53", "BRCA1", "BRCA2", "MYC", "KRAS", "EGFR", "PTEN", "APC"]
  let lfcs = [2.5, -3.0, 1.8, 3.2, -1.5, 0.8, -2.2, 0.3]
  let pvals = [0.001, 0.005, 0.02, 0.0005, 0.01, 0.08, 0.002, 0.06]
  let result = @bio.volcano_plot(ids, lfcs, pvals, p_cutoff=0.05, fc_cutoff=1.0, title="Cancer DE Genes", x_label="log2FC", y_label="-log10(p)")

  println("Total genes: \{result.get_n_genes()}")
  println("Up-regulated: \{result.get_n_up()}")
  println("Down-regulated: \{result.get_n_down()}")
  println("Non-significant: \{result.get_n_nonsig()}")
  println("")

  // Classification
  println("=== Gene Classifications ===")
  let up_genes = result.get_up_genes()
  let mut i = 0
  while i < up_genes.length() {
    println("  UP: \{up_genes[i].label} (log2FC=@{up_genes[i].log2fc}, p=@{up_genes[i].p_val})")
    i = i + 1
  }

  let down_genes = result.get_down_genes()
  let mut j = 0
  while j < down_genes.length() {
    println("  DOWN: \{down_genes[j].label} (log2FC=@{down_genes[j].log2fc}, p=@{down_genes[j].p_val})")
    j = j + 1
  }
  println("")

  // Top genes
  println("=== Top 5 Significant Genes ===")
  let top = result.get_top_genes(5)
  let mut k = 0
  while k < top.length() {
    let g = top[k]
    println("  \{g.label}: log2FC=@{g.log2fc}, p=@{g.p_val}")
    k = k + 1
  }
  println("")

  // Significance thresholds
  println("=== Significance Thresholds ===")
  println("P-value cutoff: @{result.get_p_cutoff()}")
  println("Fold-change cutoff: |log2FC| >= @{result.get_fc_cutoff()}")
  println("")

  // -log10(p-value) calculation
  println("=== -log10(p-value) Calculations ===")
  println("p=0.05 -> -log10 = @{@bio.neg_log10_p(0.05)}")
  println("p=0.01 -> -log10 = @{@bio.neg_log10_p(0.01)}")
  println("p=0.001 -> -log10 = @{@bio.neg_log10_p(0.001)}")
  println("p=1e-10 -> -log10 = @{@bio.neg_log10_p(1e-10)}")
  println("")

  // Summary
  println("=== Summary ===")
  println(result.summary())

  // Sample data
  println("=== Sample Dataset ===")
  let sample = @bio.volcano_sample(50)
  println("Sample size: \{sample.get_n_genes()}")
  println("Significant: @{sample.get_n_up() + sample.get_n_down()}")
  println("  Up: \{sample.get_n_up()}, Down: \{sample.get_n_down()}")

  // Render ASCII
  println("")
  println("=== ASCII Volcano Plot ===")
  let ascii = sample.to_ascii(width=70, height=15)
  println(ascii)

  println("=== Done ===")
}