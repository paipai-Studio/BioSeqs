///|
/// decoupleR module example - Functional activity inference.
///
/// Run: moon run examples/decoupler_demo
fn main {
  println("=== Bioconductor decoupleR Module Demo ===")
  println("Functional Activity Inference from Prior Knowledge Networks")
  println("")

  // ---- Section 1: Load sample data ----
  println("=== 1. Load Sample Data ===")
  let (expr, samples, genes, pkn) = @src.decoupler_sample_data()
  println("Samples: " + samples.length().to_string())
  println("Genes: " + genes.length().to_string())
  println("PKN edges: " + pkn.edges.length().to_string())
  println("Regulators: " + pkn.regulators.length().to_string())
  println("Targets: " + pkn.targets.length().to_string())
  println("")

  // ---- Section 2: PKN edges ----
  println("=== 2. Prior Knowledge Network ===")
  let mut i = 0
  while i < pkn.edges.length() {
    let e = pkn.edges[i]
    let mode = if e.weight > 0.0 { "activation" } else { "inhibition" }
    println("  " + e.source + " -> " + e.target + " (" + mode + ", w=" + e.weight.to_string() + ")")
    i = i + 1
  }
  println("")

  // ---- Section 3: Run with WSum method ----
  println("=== 3. WSum Method ===")
  let wsum_result = @src.decoupler_run(expr, samples, genes, pkn, method=@src.decoupler_wsum_method())
  println("Method: " + wsum_result.method)
  println("Regulators: " + wsum_result.regulators.length().to_string())
  // Header row
  let mut hdr = "          "
  let mut si = 0
  while si < samples.length() {
    hdr = hdr + samples[si] + "   "
    si = si + 1
  }
  println(hdr)
  // Data rows
  let mut ri = 0
  while ri < wsum_result.regulators.length() {
    let mut row = wsum_result.regulators[ri] + "  "
    let mut sj = 0
    while sj < wsum_result.matrix[ri].length() {
      row = row + wsum_result.matrix[ri][sj].to_string() + "  "
      sj = sj + 1
    }
    println(row)
    ri = ri + 1
  }
  println("")

  // ---- Section 4: Run with ULM method ----
  println("=== 4. ULM Method ===")
  let ulm_result = @src.decoupler_run(expr, samples, genes, pkn, method=@src.decoupler_ulm_method())
  println("Method: " + ulm_result.method)
  let mut hdr2 = "          "
  let mut si2 = 0
  while si2 < samples.length() {
    hdr2 = hdr2 + samples[si2] + "   "
    si2 = si2 + 1
  }
  println(hdr2)
  let mut ri2 = 0
  while ri2 < ulm_result.regulators.length() {
    let mut row = ulm_result.regulators[ri2] + "  "
    let mut sj2 = 0
    while sj2 < ulm_result.matrix[ri2].length() {
      row = row + ulm_result.matrix[ri2][sj2].to_string() + "  "
      sj2 = sj2 + 1
    }
    println(row)
    ri2 = ri2 + 1
  }
  println("")

  // ---- Section 5: Run with MLM method ----
  println("=== 5. MLM Method ===")
  let mlm_result = @src.decoupler_run(expr, samples, genes, pkn, method=@src.decoupler_mlm_method())
  println("Method: " + mlm_result.method)
  let mut hdr3 = "          "
  let mut si3 = 0
  while si3 < samples.length() {
    hdr3 = hdr3 + samples[si3] + "   "
    si3 = si3 + 1
  }
  println(hdr3)
  let mut ri3 = 0
  while ri3 < mlm_result.regulators.length() {
    let mut row = mlm_result.regulators[ri3] + "  "
    let mut sj3 = 0
    while sj3 < mlm_result.matrix[ri3].length() {
      row = row + mlm_result.matrix[ri3][sj3].to_string() + "  "
      sj3 = sj3 + 1
    }
    println(row)
    ri3 = ri3 + 1
  }
  println("")

  // ---- Section 6: Top regulators ----
  println("=== 6. Top Regulators for S1 ===")
  let top = @src.decoupler_top_regulators(wsum_result, "S1", 3)
  let mut ti = 0
  while ti < top.length() {
    let s = top[ti]
    println("  " + (ti + 1).to_string() + ". " + s.regulator + " (score=" + s.score.to_string() + ")")
    ti = ti + 1
  }
  println("")

  // ---- Section 7: Compare methods ----
  println("=== 7. Compare All Methods ===")
  let methods = [
    (@src.decoupler_wsum_method(), "wsum"),
    (@src.decoupler_wmean_method(), "wmean"),
    (@src.decoupler_norm_method(), "norm"),
    (@src.decoupler_ulm_method(), "ulm"),
    (@src.decoupler_mlm_method(), "mlm"),
  ]
  let mut mi = 0
  while mi < methods.length() {
    let (m, name) = methods[mi]
    let r = @src.decoupler_run(expr, samples, genes, pkn, method=m)
    println("  " + name + ": " + r.scores.length().to_string() + " scores")
    mi = mi + 1
  }
  println("")

  println("=== Done ===")
}
