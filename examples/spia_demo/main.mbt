///|
/// SPIA module example - Signaling Pathway Impact Analysis.
///
/// Run: moon run examples/spia_demo
fn main {
  println("=== SPIA Module Demo ===")
  println("Signaling Pathway Impact Analysis")
  println("")

  // ---- Section 1: Load sample pathways ----
  println("=== 1. Load Sample Pathways ===")
  let pathways = @bio.spia_sample_pathways()
  println("Total pathways: " + pathways.length().to_string())
  let mut pi = 0
  while pi < pathways.length() {
    let p = pathways[pi]
    println("  " + p.id + " (" + p.name + "): " + p.nodes.length().to_string() + " nodes")
    // Print nodes and their interactions
    let mut ni = 0
    while ni < p.nodes.length() {
      let node = p.nodes[ni]
      let mut downstream_str = "".to_string()
      let mut di = 0
      while di < node.downstream.length() {
        let arrow = if node.activation[di] > 0 { " ->" } else { " -|" }
        downstream_str = downstream_str + arrow + " " + node.downstream[di]
        if di < node.downstream.length() - 1 {
          downstream_str = downstream_str + ","
        }
        di = di + 1
      }
      if node.downstream.length() == 0 {
        downstream_str = "(terminal)"
      }
      println("    " + node.gene_id + " (" + node.gene_name + ") " + downstream_str)
      ni = ni + 1
    }
    pi = pi + 1
  }
  println("")

  // ---- Section 2: Load DE genes ----
  println("=== 2. Differentially Expressed Genes ===")
  let de = @bio.spia_sample_de_genes()
  println("Total DE genes: " + de.size().to_string())
  let keys = de.keys().collect()
  let mut ki = 0
  while ki < keys.length() {
    let gid = keys[ki]
    let lfc = de.get(gid).unwrap_or(0.0)
    let direction = if lfc > 0.0 { "UP" } else if lfc < 0.0 { "DOWN" } else { "NEUTRAL" }
    println("  " + gid + ": log2FC=" + lfc.to_string() + " (" + direction + ")")
    ki = ki + 1
  }
  println("")

  // ---- Section 3: Run SPIA ----
  println("=== 3. Run SPIA ===")
  let results = @bio.run_spia(pathways, de, 100, fdr_threshold=0.05)
  println(results.summary())
  println("")

  // ---- Section 4: Per-pathway results ----
  println("=== 4. Per-Pathway Results ===")
  let mut i = 0
  while i < results.results.length() {
    let r = results.results[i]
    let status = if r.activation_status > 0 { "Activated" } else if r.activation_status < 0 { "Inhibited" } else { "Neutral" }
    println("  " + r.pathway_id + " (" + r.pathway_name + "):")
    println("    n_genes=" + r.n_genes.to_string() + ", n_de=" + r.n_de.to_string())
    println("    p_or=" + r.p_or.to_string() + " (over-representation)")
    println("    pert_factor=" + r.pert_factor.to_string() + ", p_pert=" + r.p_pert.to_string())
    println("    p_combined=" + r.p_combined.to_string() + ", p_adj=" + r.p_adj.to_string())
    println("    status=" + status + ", significant=" + r.significant.to_string())
    i = i + 1
  }
  println("")

  // ---- Section 5: Activated pathways ----
  println("=== 5. Activated Pathways ===")
  let activated = results.get_activated()
  println("Total activated: " + activated.length().to_string())
  let mut j = 0
  while j < activated.length() {
    println("  " + activated[j].pathway_id + " (" + activated[j].pathway_name + "): pert_factor=" + activated[j].pert_factor.to_string())
    j = j + 1
  }
  println("")

  // ---- Section 6: Inhibited pathways ----
  println("=== 6. Inhibited Pathways ===")
  let inhibited = results.get_inhibited()
  println("Total inhibited: " + inhibited.length().to_string())
  let mut k = 0
  while k < inhibited.length() {
    println("  " + inhibited[k].pathway_id + " (" + inhibited[k].pathway_name + "): pert_factor=" + inhibited[k].pert_factor.to_string())
    k = k + 1
  }
  println("")

  // ---- Section 7: Top ranked pathways ----
  println("=== 7. Top Ranked Pathways ===")
  let top = results.get_top_pathways(2)
  let mut n = 0
  while n < top.length() {
    println("  Rank " + (n + 1).to_string() + ": " + top[n].pathway_id + " (" + top[n].pathway_name + "), p_combined=" + top[n].p_combined.to_string())
    n = n + 1
  }
  println("")

  // ---- Section 8: Significant pathways ----
  println("=== 8. Significant Pathways ===")
  let sig = results.get_significant()
  println("Total significant: " + sig.length().to_string())
  let mut m = 0
  while m < sig.length() {
    println("  " + sig[m].pathway_id + " (" + sig[m].pathway_name + "): p_adj=" + sig[m].p_adj.to_string())
    m = m + 1
  }
  println("")

  // ---- Section 9: Permissive threshold ----
  println("=== 9. Permissive Threshold (FDR < 1.0) ===")
  let perm_results = @bio.run_spia(pathways, de, 100, fdr_threshold=1.0)
  let perm_sig = perm_results.get_significant()
  println("Significant at FDR < 1.0: " + perm_sig.length().to_string())
  println("")

  // ---- Section 10: Strict threshold ----
  println("=== 10. Strict Threshold (FDR < 0.0) ===")
  let strict_results = @bio.run_spia(pathways, de, 100, fdr_threshold=0.0)
  let strict_sig = strict_results.get_significant()
  println("Significant at FDR < 0.0: " + strict_sig.length().to_string())
  println("")

  println("=== Demo complete ===")
}
