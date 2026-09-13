///|
fn main {
  println("=== affyPLM Demo ===")

  let ps = @src.affy_plm_sample_probeset()
  println("\n1. Probeset: " + ps.name)
  println("   Probes: " + ps.probes.length().to_string())
  println("   Samples: " + ps.probes[0].length().to_string())
  for i in 0..<ps.probes.length() {
    let mut row = ""
    for j in 0..<ps.probes[i].length() {
      row = row + ps.probes[i][j].to_string() + "  "
    }
    println("   probe " + i.to_string() + ": " + row)
  }

  println("\n2. Fitting PLM (Huber M-estimator, IRLS)")
  let fit = @src.fit_plm(ps)
  println("   Iterations: " + fit.n_iterations.to_string())
  println("   Sigma: " + fit.sigma.to_string())

  println("\n3. Chip effects (sample-level)")
  for j in 0..<fit.chip_effects.length() {
    println("   sample " + j.to_string() + ": " + fit.chip_effects[j].to_string())
  }

  println("\n4. Probe effects (probe-level)")
  for i in 0..<fit.probe_effects.length() {
    println("   probe " + i.to_string() + ": " + fit.probe_effects[i].to_string())
  }

  println("\n5. Residuals")
  for i in 0..<fit.residuals.length() {
    let mut row = ""
    for j in 0..<fit.residuals[i].length() {
      row = row + fit.residuals[i][j].to_string() + "  "
    }
    println("   probe " + i.to_string() + ": " + row)
  }

  println("\n6. QC statistics (NUSE / RLE)")
  let report = @src.affy_plm_qc([fit])
  println("   NUSE: " + report.nuse.to_string())
  println("   RLE:  " + report.rle.to_string())

  println("\n7. Tukey biweight estimator")
  let fit_tukey = @src.fit_plm(ps, estimator=@src.PlmEstimator::Tukey(4.685))
  println("   Iterations: " + fit_tukey.n_iterations.to_string())
  println("   Sigma: " + fit_tukey.sigma.to_string())
}