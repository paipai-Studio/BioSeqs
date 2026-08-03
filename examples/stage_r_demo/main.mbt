///|
/// stageR module demo - two-stage hypothesis testing with OFDR control.
///
/// Demonstrates the stage-wise testing paradigm that controls the Overall
/// False Discovery Rate (OFDR) at the gene level while allowing post-hoc
/// assessment of individual confirmation hypotheses (e.g. transcript-level tests).
///
/// Reference: Van den Berge et al. (2017) Genome Biology 18:151.
/// Run: moon run examples/stage_r_demo
fn main {
  println("=== stageR Module Demo ===")
  println("Two-stage hypothesis testing with OFDR control")
  println("Reference: Bioconductor stageR")
  println("")

  // ------------------------------------------------------------------
  // 1. Load sample data
  // ------------------------------------------------------------------
  println("=== 1. Sample Data ===")
  let (pScreen, pConfirmation) = @bio.stage_r_sample()
  println("Genes: " + pScreen.length().to_string())
  println("Confirmation hypotheses per gene: " + pConfirmation[0].length().to_string())
  println("Screening p-values: " + p_screen_to_string(pScreen))
  println("")

  // ------------------------------------------------------------------
  // 2. BH adjustment
  // ------------------------------------------------------------------
  println("=== 2. Benjamini-Hochberg FDR Adjustment ===")
  let bh_adj = @bio.stage_r_bh_adjust(pScreen)
  println("BH-adjusted screening p-values:")
  let mut i = 0
  while i < bh_adj.length() {
    println("  Gene " + i.to_string() + ": raw=" + pScreen[i].to_string() + ", adj=" + bh_adj[i].to_string())
    i = i + 1
  }
  println("")

  // ------------------------------------------------------------------
  // 3. Simes aggregation
  // ------------------------------------------------------------------
  println("=== 3. Simes Aggregation ===")
  let gene0_pvals = pConfirmation[0]
  println("Gene 0 confirmation p-values: " + p_conf_to_string(gene0_pvals))
  let simes = @bio.simes_aggregate(gene0_pvals)
  println("Simes aggregated p-value: " + simes.to_string())
  println("")

  // ------------------------------------------------------------------
  // 4. Stage-wise adjustment with Holm method
  // ------------------------------------------------------------------
  println("=== 4. Stage-wise Adjustment (Holm) ===")
  let config = @bio.StageRConfig::new(
    alpha=0.05,
    pScreenAdjusted=false,
    stageMethod=@bio.stage_r_method_holm(),
  )
  let result = @bio.stage_wise_adjustment(pScreen, pConfirmation, config)
  println("Total genes: " + result.nGenes.to_string())
  println("Significant genes (screening): " + result.nSignificantGenes.to_string())
  println("Adjusted alpha (R/G * alpha): " + result.alphaAdjusted.to_string())
  println("")

  // ------------------------------------------------------------------
  // 5. Significant genes and hypotheses
  // ------------------------------------------------------------------
  println("=== 5. Significant Results ===")
  let sig_genes = @bio.get_significant_genes(result)
  println("Significant gene indices: " + int_array_to_string(sig_genes))
  let sig_hyps = @bio.get_significant_hypotheses(result)
  println("Significant confirmation hypotheses: " + sig_hyps.length().to_string())
  let mut h = 0
  while h < sig_hyps.length() {
    let (g, j) = sig_hyps[h]
    println("  Gene " + g.to_string() + ", Hypothesis " + j.to_string())
    h = h + 1
  }
  println("")

  // ------------------------------------------------------------------
  // 6. Results matrix
  // ------------------------------------------------------------------
  println("=== 6. Results Matrix (gene x hypothesis+1) ===")
  let mat = @bio.get_results(result)
  println("Column 0 = gene-level, Columns 1..H = confirmation")
  let mut gi = 0
  while gi < mat.length() {
    let row_str = bool_row_to_string(mat[gi])
    println("  Gene " + gi.to_string() + ": " + row_str)
    gi = gi + 1
  }
  println("")

  // ------------------------------------------------------------------
  // 7. Compare methods
  // ------------------------------------------------------------------
  println("=== 7. Method Comparison ===")
  let methods = [
    ("none", @bio.stage_r_method_none()),
    ("holm", @bio.stage_r_method_holm()),
    ("dte", @bio.stage_r_method_dte()),
    ("dtu", @bio.stage_r_method_dtu()),
  ]
  let mut m = 0
  while m < methods.length() {
    let (name, method) = methods[m]
    let cfg = @bio.StageRConfig::new(alpha=0.05, pScreenAdjusted=false, stageMethod=method)
    let res = @bio.stage_wise_adjustment(pScreen, pConfirmation, cfg)
    let sigs = @bio.get_significant_hypotheses(res)
    println("  " + name + ": " + sigs.length().to_string() + " significant hypotheses")
    m = m + 1
  }
  println("")

  // ------------------------------------------------------------------
  // 8. Summary
  // ------------------------------------------------------------------
  println("=== 8. Summary ===")
  let summary = @bio.stage_r_summary(result)
  println(summary)
  println("=== Demo complete ===")
}

fn p_screen_to_string(arr : Array[Double]) -> String {
  let buf = StringBuilder::new()
  buf.write_string("[")
  let mut i = 0
  while i < arr.length() {
    if i > 0 {
      buf.write_string(", ")
    }
    buf.write_string(arr[i].to_string())
    i = i + 1
  }
  buf.write_string("]")
  buf.to_string()
}

fn p_conf_to_string(arr : Array[Double]) -> String {
  let buf = StringBuilder::new()
  buf.write_string("[")
  let mut i = 0
  while i < arr.length() {
    if i > 0 {
      buf.write_string(", ")
    }
    buf.write_string(arr[i].to_string())
    i = i + 1
  }
  buf.write_string("]")
  buf.to_string()
}

fn int_array_to_string(arr : Array[Int]) -> String {
  let buf = StringBuilder::new()
  buf.write_string("[")
  let mut i = 0
  while i < arr.length() {
    if i > 0 {
      buf.write_string(", ")
    }
    buf.write_string(arr[i].to_string())
    i = i + 1
  }
  buf.write_string("]")
  buf.to_string()
}

fn bool_row_to_string(row : Array[Bool]) -> String {
  let buf = StringBuilder::new()
  buf.write_string("[")
  let mut i = 0
  while i < row.length() {
    if i > 0 {
      buf.write_string(", ")
    }
    buf.write_string(if row[i] { "TRUE" } else { "FALSE" })
    i = i + 1
  }
  buf.write_string("]")
  buf.to_string()
}
