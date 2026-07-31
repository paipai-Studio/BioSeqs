///|
/// GAGE module example - Generally Applicable Gene-set Enrichment analysis.
///
/// Run: moon run examples/gage_demo
fn main {
  println("=== GAGE Module Demo ===")
  println("Generally Applicable Gene-set Enrichment Analysis")
  println("")

  // ---- Section 1: Load sample data ----
  println("=== 1. Load Sample Data ===")
  let fcs = @bio.gage_sample_fold_changes()
  let sets = @bio.gage_sample_gene_sets()
  println("Gene fold changes: " + fcs.size().to_string())
  let mut si = 0
  while si < sets.length() {
    println("  " + sets[si].id + " (" + sets[si].name + "): " + sets[si].gene_ids.length().to_string() + " genes")
    si = si + 1
  }
  println("")

  // ---- Section 2: Run GAGE analysis ----
  println("=== 2. Run GAGE Analysis ===")
  let results = @bio.run_gage(fcs, sets, fdr_threshold=0.05)
  println(results.summary())
  println("")

  // ---- Section 3: Per-set results ----
  println("=== 3. Per-Set Results ===")
  let mut i = 0
  while i < results.results.length() {
    let r = results.results[i]
    let direction = if r.is_up { "UP" } else if r.is_down { "DOWN" } else { "NEUTRAL" }
    println("  " + r.gene_set_id + " (" + r.gene_set_name + "):")
    println("    n_genes=" + r.n_genes.to_string() + ", mean_fc=" + r.mean_fc.to_string() + ", stat=" + r.stat.to_string())
    println("    direction=" + direction + ", p_val=" + r.p_val.to_string() + ", p_adj=" + r.p_adj.to_string() + ", sig=" + r.significant.to_string())
    i = i + 1
  }
  println("")

  // ---- Section 4: Significant gene sets ----
  println("=== 4. Significant Gene Sets ===")
  let sig = results.get_significant()
  println("Total significant: " + sig.length().to_string())
  let mut j = 0
  while j < sig.length() {
    println("  " + sig[j].gene_set_id + ": p_adj=" + sig[j].p_adj.to_string())
    j = j + 1
  }
  println("")

  // ---- Section 5: Up-regulated gene sets ----
  println("=== 5. Up-regulated Gene Sets ===")
  let up = results.get_upregulated()
  println("Total up-regulated: " + up.length().to_string())
  let mut k = 0
  while k < up.length() {
    println("  " + up[k].gene_set_id + " (" + up[k].gene_set_name + "): mean_fc=" + up[k].mean_fc.to_string())
    k = k + 1
  }
  println("")

  // ---- Section 6: Down-regulated gene sets ----
  println("=== 6. Down-regulated Gene Sets ===")
  let down = results.get_downregulated()
  println("Total down-regulated: " + down.length().to_string())
  let mut m = 0
  while m < down.length() {
    println("  " + down[m].gene_set_id + " (" + down[m].gene_set_name + "): mean_fc=" + down[m].mean_fc.to_string())
    m = m + 1
  }
  println("")

  // ---- Section 7: Top ranked gene sets ----
  println("=== 7. Top Ranked Gene Sets ===")
  let top = results.get_top_sets(3)
  let mut n = 0
  while n < top.length() {
    println("  Rank " + (n + 1).to_string() + ": " + top[n].gene_set_id + " (p_adj=" + top[n].p_adj.to_string() + ")")
    n = n + 1
  }
  println("")

  // ---- Section 8: Paired test ----
  println("=== 8. Paired Test ===")
  let paired_results = @bio.run_gage(fcs, sets, fdr_threshold=0.05, paired=true)
  println(paired_results.summary())
  println("")

  // ---- Section 9: Strict FDR threshold ----
  println("=== 9. Strict FDR Threshold (0.001) ===")
  let strict_results = @bio.run_gage(fcs, sets, fdr_threshold=0.001)
  let strict_sig = strict_results.get_significant()
  println("Significant at FDR < 0.001: " + strict_sig.length().to_string())
  println("")

  // ---- Section 10: Permissive FDR threshold ----
  println("=== 10. Permissive FDR Threshold (1.0) ===")
  let perm_results = @bio.run_gage(fcs, sets, fdr_threshold=1.0)
  let perm_sig = perm_results.get_significant()
  println("Significant at FDR < 1.0: " + perm_sig.length().to_string())
  println("")

  println("=== Demo complete ===")
}
