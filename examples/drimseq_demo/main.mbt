///|
/// DRIMSeq module example - Differential transcript usage analysis.
///
/// Run: moon run examples/drimseq_demo
fn main {
  println("=== DRIMSeq Module Demo ===")
  println("Differential Transcript Usage Analysis")
  println("")

  // Create test data
  println("=== Sample Data ===")
  let counts = @bio.drimseq_sample_data()
  println("Total transcript counts: " + counts.length().to_string())
  let mut i = 0
  let genes : Array[String] = Array::new()
  while i < counts.length() {
    let g = counts[i].gene_id
    let mut found = false
    let mut j = 0
    while j < genes.length() {
      if genes[j] == g { found = true; break }
      j = j + 1
    }
    if !found { genes.push(g) }
    i = i + 1
  }
  println("Unique genes: " + genes.length().to_string())
  let mut g = 0
  while g < genes.length() {
    println("  " + genes[g])
    g = g + 1
  }
  println("")

  // Compute proportions
  println("=== Transcript Proportions ===")
  let gene_counts = @bio.drimseq_aggregate_by_gene(counts)
  let props = @bio.drimseq_compute_proportions(counts, gene_counts)
  println("Proportion records: " + props.length().to_string())
  let mut p = 0
  while p < props.length() && p < 6 {
    let rec = props[p]
    println("  " + rec.transcript_id + " (" + rec.gene_id + "): proportion=" + rec.proportion().to_string() + " (count=" + rec.count.to_string() + ")")
    p = p + 1
  }
  println("")

  // Filter counts
  println("=== Filtering ===")
  let filtered = @bio.drimseq_filter_counts(counts)
  println("Records after filtering: " + filtered.length().to_string())
  println("")

  // Run DRIMSeq analysis
  println("=== DRIMSeq Analysis ===")
  let config = @bio.DRIMSeqConfig::new()
  let result = @bio.drimseq_test_differential(counts, config=config)
  println(result.summary())
  println("")

  // Show results
  println("=== Gene Results ===")
  let top = result.get_top_genes(n=10)
  let mut t = 0
  while t < top.length() {
    let g = top[t]
    println("  " + g.gene_id + ": p=" + g.p_value.to_string() + ", adj_p=" + g.p_adj.to_string() + ", lr=" + g.lr_stat.to_string())
    t = t + 1
  }
  println("")

  // Significant results
  println("=== Significant Genes (alpha=" + config.alpha.to_string() + ") ===")
  let sig = result.get_significant()
  println("Number of significant genes: " + sig.length().to_string())
  let mut s = 0
  while s < sig.length() {
    let g = sig[s]
    println("  " + g.gene_id + ": delta=" + g.delta.to_string())
    s = s + 1
  }
  println("")

  // Config options
  println("=== Configuration Options ===")
  let cfg_low = @bio.DRIMSeqConfig::new().set_min_count(val=5.0)
  let result_low = @bio.drimseq_test_differential(counts, config=cfg_low)
  println("With min_count=5: tested=" + result_low.n_genes_tested.to_string() + ", filtered=" + result_low.n_genes_filtered.to_string())

  let cfg_alpha = @bio.DRIMSeqConfig::new().set_alpha(val=0.01)
  let result_alpha = @bio.drimseq_test_differential(counts, config=cfg_alpha)
  println("With alpha=0.01: significant=" + result_alpha.n_significant.to_string())
  println("")

  // Dirichlet-multinomial demonstration
  println("=== Dirichlet-Multinomial Log-Likelihood ===")
  let test_counts : Array[Double] = [50.0, 30.0, 20.0]
  let test_props : Array[Double] = [0.5, 0.3, 0.2]
  let ll = @bio.drimseq_dirichlet_multinomial(counts=test_counts, proportions=test_props)
  println("Log-likelihood: " + ll.to_string())
  println("")

  // Wald test demonstration
  println("=== Wald Test ===")
  let a : Array[Double] = [0.7, 0.3]
  let b : Array[Double] = [0.4, 0.6]
  let (stat, pval, delta) = @bio.drimseq_wald_test(a, b)
  println("  stat=" + stat.to_string() + ", pval=" + pval.to_string() + ", delta=" + delta.to_string())

  let c : Array[Double] = [0.5, 0.5]
  let d : Array[Double] = [0.5, 0.5]
  let (stat2, pval2, _) = @bio.drimseq_wald_test(c, d)
  println("  Same proportions: stat=" + stat2.to_string() + ", pval=" + pval2.to_string())
}
