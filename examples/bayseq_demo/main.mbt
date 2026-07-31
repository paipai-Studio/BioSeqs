///|
/// baySeq module example - Bayesian differential expression analysis.
///
/// Run: moon run examples/bayseq_demo
fn main {
  println("=== Bioconductor baySeq Module Demo ===")
  println("Bayesian Analysis of RNA-seq Differential Expression")
  println("")

  // ---- Section 1: Load sample data ----
  println("=== 1. Load Sample RNA-seq Count Data ===")
  let (counts, groups, gene_names) = @bio.bayseq_sample_data()
  println("Input: " + counts.length().to_string() + " genes x " +
    counts[0].length().to_string() + " samples")
  println("Groups: 0 (control) x " + @bio.hts_count_group(groups, 0).to_string() +
    ", 1 (treatment) x " + @bio.hts_count_group(groups, 1).to_string())
  println("")

  // ---- Section 2: Estimate dispersion prior ----
  println("=== 2. Estimate Dispersion Prior (Gamma) ===")
  let (shape, rate) = @bio.bayseq_estimate_prior(counts, groups)
  println("Gamma shape: " + shape.to_string())
  println("Gamma rate: " + rate.to_string())
  println("")

  // ---- Section 3: Run baySeq test ----
  println("=== 3. Run baySeq Differential Expression Test ===")
  let result = @bio.bayseq_test(counts, groups, gene_names)
  println("Total genes analyzed: " + result.gene_results.length().to_string())
  println("DE genes (posterior >= " + result.alpha_level.to_string() + "): " +
    result.n_de_genes.to_string())
  println("")

  // ---- Section 4: Top DE genes ----
  println("=== 4. Top 10 DE Genes by |logFC| ===")
  let top = @bio.bayseq_get_top_de(result, 10)
  let mut ti = 0
  while ti < top.length() {
    let g = top[ti]
    let direction = if g.log_fold_change > 0.0 { "UP" } else { "DOWN" }
    println("  " + (ti + 1).to_string() + ". " + g.gene_name +
      " | logFC=" + g.log_fold_change.to_string() +
      " | Posterior=" + g.posterior_prob_de.to_string() +
      " | " + direction)
    ti = ti + 1
  }
  println("")

  // ---- Section 5: Significant DE genes ----
  println("=== 5. Significant DE Genes (posterior >= threshold) ===")
  let de_genes = @bio.bayseq_get_de_genes(result)
  if de_genes.length() > 0 {
    println("Total significant DE genes: " + de_genes.length().to_string())
    println("First 5 significant DE genes:")
    let mut di = 0
    while di < 5 && di < de_genes.length() {
      let g = de_genes[di]
      println("  " + g.gene_name +
        " | mean_ctrl=" + g.map_expression_a.to_string() +
        " | mean_treat=" + g.map_expression_b.to_string() +
        " | disp=" + g.dispersion.to_string())
      di = di + 1
    }
  } else {
    println("No significant DE genes found.")
  }
  println("")

  // ---- Section 6: Summary ----
  println("=== 6. Analysis Summary ===")
  println(@bio.bayseq_summary(result))
}
