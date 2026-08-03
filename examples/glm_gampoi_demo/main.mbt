///|
/// glmGamPoi module demo - single-cell RNA-seq differential expression workflow.
///
/// Run: moon run examples/glm_gampoi_demo
fn main {
  println("=== glmGamPoi Module Demo ===")
  println("Single-cell RNA-seq differential expression analysis")
  println("Reference: Ahlmann-Eltze C & Huber W (2021) Bioinformatics 37(2):299-300")
  println("")

  // ------------------------------------------------------------------
  // 1. Create sample count data (5 genes x 8 samples)
  // ------------------------------------------------------------------
  println("=== 1. Create sample count data ===")
  let gene_names = ["Gene_A", "Gene_B", "Gene_C", "Gene_D", "Gene_E"]
  let sample_names = ["ctrl_1", "ctrl_2", "ctrl_3", "ctrl_4", "treat_1", "treat_2", "treat_3", "treat_4"]

  let counts : Array[Array[Double]] = [
    [10.0, 12.0, 11.0, 9.0, 50.0, 55.0, 48.0, 52.0],
    [5.0, 6.0, 5.0, 7.0, 15.0, 18.0, 14.0, 17.0],
    [20.0, 18.0, 22.0, 19.0, 21.0, 20.0, 23.0, 19.0],
    [50.0, 48.0, 52.0, 49.0, 20.0, 22.0, 18.0, 21.0],
    [2.0, 1.0, 3.0, 2.0, 1.0, 2.0, 1.0, 0.0],
  ]

  println("Genes: \{gene_names.length()}")
  println("Samples: \{sample_names.length()} (4 control + 4 treatment)")
  println("")

  // ------------------------------------------------------------------
  // 2. Create design matrix (intercept + treatment effect)
  // ------------------------------------------------------------------
  println("=== 2. Design matrix ===")
  let coef_names = ["Intercept", "Treatment"]
  let matrix : Array[Array[Double]] = [
    [1.0, 0.0],
    [1.0, 0.0],
    [1.0, 0.0],
    [1.0, 0.0],
    [1.0, 1.0],
    [1.0, 1.0],
    [1.0, 1.0],
    [1.0, 1.0],
  ]
  let design = @bio.GlmDesign::new(sample_names, coef_names, matrix)
  println("Samples (rows): \{design.n_samples()}")
  println("Coefficients (cols): \{design.n_coefs()}")
  println("")

  // ------------------------------------------------------------------
  // 3. Calculate size factors
  // ------------------------------------------------------------------
  println("=== 3. Size factors ===")
  let sf = @bio.glm_calculate_sf(counts)
  let mut si = 0
  while si < sf.length() {
    println("  \{sample_names[si]}: sf=\{sf[si]}")
    si = si + 1
  }
  println("")

  // ------------------------------------------------------------------
  // 4. Fit one gene and show coefficients
  // ------------------------------------------------------------------
  println("=== 4. Fit single gene (Gene_A) ===")
  let gene_a_counts = counts[0]
  let fit_a = @bio.glm_fit_one_gene(gene_a_counts, design, sf, "Gene_A")
  let mut ci = 0
  while ci < fit_a.coefficients.length() {
    println("  coef[\{ci}] (\{coef_names[ci]}): est=\{fit_a.coefficients[ci]} se=\{fit_a.std_errors[ci]}")
    ci = ci + 1
  }
  println("  dispersion: \{fit_a.dispersion}")
  println("  df_residual: \{fit_a.df_residual}")
  println("")

  // ------------------------------------------------------------------
  // 5. Test for differential expression (all genes)
  // ------------------------------------------------------------------
  println("=== 5. Differential expression test (Wald test) ===")
  let results = @bio.glm_test_de(counts, gene_names, design, sf, 1)
  let mut ri = 0
  while ri < results.length() {
    let r = results[ri]
    println("  \{r.to_string()}")
    ri = ri + 1
  }
  println("")

  // ------------------------------------------------------------------
  // 6. BH-FDR correction
  // ------------------------------------------------------------------
  println("=== 6. BH-FDR correction ===")
  let pvalues : Array[Double] = []
  let mut pi = 0
  while pi < results.length() {
    pvalues.push(results[pi].p_value())
    pi = pi + 1
  }
  let adj_pvalues = @bio.glm_bh_correct_array(pvalues)
  let mut ai = 0
  while ai < results.length() {
    println("  \{gene_names[ai]}: raw p=\{pvalues[ai]} adj p=\{adj_pvalues[ai]}")
    ai = ai + 1
  }
  println("")

  // ------------------------------------------------------------------
  // 7. Find significant genes
  // ------------------------------------------------------------------
  println("=== 7. Significant genes (adj_p < 0.05) ===")
  let sig_results = @bio.glm_significant(results, 0.05)
  println("Found \{sig_results.length()} significant gene(s)")
  let mut sri = 0
  while sri < sig_results.length() {
    let sr = sig_results[sri]
    println("  \{sr.to_string()}")
    sri = sri + 1
  }
  println("")

  // ------------------------------------------------------------------
  // 8. Linear algebra utility: solve 2x2 system
  // ------------------------------------------------------------------
  println("=== 8. Linear algebra: solve 2x2 system ===")
  println("  System: 3*x + 2*y = 7")
  println("          1*x + 4*y = 9")
  let a : Array[Array[Double]] = [[3.0, 2.0], [1.0, 4.0]]
  let b : Array[Double] = [7.0, 9.0]
  let x = @bio.glm_solve_linear_system(a, b, 2)
  println("  Solution: x=\{x[0]}, y=\{x[1]}")
  println("  Verify: 3*\{x[0]} + 2*\{x[1]} = \{3.0 * x[0] + 2.0 * x[1]} (expected 7)")
  println("  Verify: 1*\{x[0]} + 4*\{x[1]} = \{1.0 * x[0] + 4.0 * x[1]} (expected 9)")
  println("")

  // ------------------------------------------------------------------
  // 9. Statistical utilities
  // ------------------------------------------------------------------
  println("=== 9. Statistical utilities ===")

  println("  Standard normal CDF (glm_phi):")
  println("    phi(0.0) = \{@bio.glm_phi(0.0)} (expected ~0.5)")
  println("    phi(1.96) = \{@bio.glm_phi(1.96)} (expected ~0.975)")
  println("    phi(-1.96) = \{@bio.glm_phi(-1.96)} (expected ~0.025)")

  println("  Two-sided p-values (glm_normal_pvalue_two_sided):")
  println("    z=0.00 -> p=\{@bio.glm_normal_pvalue_two_sided(0.0)} (expected ~1.0)")
  println("    z=1.96 -> p=\{@bio.glm_normal_pvalue_two_sided(1.96)} (expected ~0.05)")
  println("    z=2.58 -> p=\{@bio.glm_normal_pvalue_two_sided(2.58)} (expected ~0.01)")

  println("  Manual z-test computation using @math.ln and @math.exp:")
  let manually_computed_z = @math.ln(50.0) - @math.ln(10.0)
  let manual_p = @bio.glm_normal_pvalue_two_sided(manually_computed_z)
  println("    ln(50) - ln(10) = \{manually_computed_z}")
  println("    p-value from manual z = \{manual_p}")
  println("")

  println("=== Demo complete ===")
}