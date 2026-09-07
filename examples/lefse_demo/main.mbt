///|
/// LEfSe (Linear Discriminant Analysis Effect Size) Demo
///
/// This demo shows how to use the LEfSe module for microbiome
/// biomarker discovery: creating abundance data with two classes,
/// running the full LEfSe pipeline, and printing significant biomarkers.
fn main {
  println("=== LEfSe Analysis Demo ===")
  println("")

  // -----------------------------------------------------------------------
  // 1. Create microbiome abundance data with two classes
  // -----------------------------------------------------------------------
  println("1. Creating microbiome abundance data:")
  let feature_names = [
    "Bacteroides", "Prevotella", "Faecalibacterium", "Roseburia",
    "Eubacterium", "Ruminococcus",
  ]
  let sample_names = [
    "H1", "H2", "H3", "H4", "H5",
    "D1", "D2", "D3", "D4", "D5",
  ]
  let classes = [
    "Healthy", "Healthy", "Healthy", "Healthy", "Healthy",
    "Disease", "Disease", "Disease", "Disease", "Disease",
  ]
  let subclasses : Array[String] = [
    "", "", "", "", "", "", "", "", "", "",
  ]
  // Bacteroides and Prevotella are elevated in Disease
  let features = [
    [10.0, 12.0, 11.0, 13.0, 10.0, 100.0, 102.0, 101.0, 103.0, 100.0],
    [20.0, 22.0, 21.0, 23.0, 20.0, 200.0, 202.0, 201.0, 203.0, 200.0],
    [50.0, 52.0, 51.0, 53.0, 50.0, 51.0, 50.0, 52.0, 49.0, 51.0],
    [30.0, 32.0, 31.0, 33.0, 30.0, 31.0, 30.0, 32.0, 29.0, 31.0],
    [40.0, 42.0, 41.0, 43.0, 40.0, 41.0, 40.0, 42.0, 39.0, 41.0],
    [25.0, 27.0, 26.0, 28.0, 25.0, 26.0, 25.0, 27.0, 24.0, 26.0],
  ]
  let input = @src.LEfSeInput::new(
    features, feature_names, sample_names, classes, subclasses,
  )
  println("   Features: \{input.feature_names.length()}")
  println("   Samples: \{input.sample_names.length()}")
  println("   Classes: Healthy (5), Disease (5)")
  println("")

  // -----------------------------------------------------------------------
  // 2. Run Kruskal-Wallis test directly
  // -----------------------------------------------------------------------
  println("2. Kruskal-Wallis test (Bacteroides):")
  let healthy_bact = [10.0, 12.0, 11.0, 13.0, 10.0]
  let disease_bact = [100.0, 102.0, 101.0, 103.0, 100.0]
  let (h, p_kw) = @src.lefse_kruskal_wallis([healthy_bact, disease_bact])
  println("   H statistic: \{h}")
  println("   p-value: \{p_kw}")
  println("")

  // -----------------------------------------------------------------------
  // 3. Run Wilcoxon rank-sum test directly
  // -----------------------------------------------------------------------
  println("3. Wilcoxon rank-sum test (Bacteroides):")
  let (u, p_wilc) = @src.lefse_wilcoxon_rank_sum(healthy_bact, disease_bact)
  println("   U statistic: \{u}")
  println("   p-value: \{p_wilc}")
  println("")

  // -----------------------------------------------------------------------
  // 4. Compute LDA score
  // -----------------------------------------------------------------------
  println("4. LDA effect size:")
  let means_bact = [11.2, 101.2] // Healthy mean, Disease mean
  let lda = @src.lefse_lda_score(means_bact, 1)
  println("   Class means: Healthy=\{means_bact[0]}, Disease=\{means_bact[1]}")
  println("   LDA score (best=Disease): \{lda}")
  println("")

  // -----------------------------------------------------------------------
  // 5. Run full LEfSe analysis
  // -----------------------------------------------------------------------
  println("5. Running full LEfSe pipeline:")
  let output = @src.lefse_run(input, alpha=0.05, lda_threshold=2.0)
  println("   Total features tested: \{output.results().length()}")
  println("   Significant biomarkers: \{output.n_significant()}")
  println("")

  // -----------------------------------------------------------------------
  // 6. Print significant biomarkers
  // -----------------------------------------------------------------------
  println("6. Significant biomarkers:")
  let mut found_sig = false
  for r in output.results() {
    if r.is_significant() {
      found_sig = true
      println("   \{r.feature_name()}: class=\{r.class_name()}, KW p=\{r.kw_pvalue()}, LDA=\{r.lda_score()}")
    }
  }
  if !found_sig {
    println("   (No significant biomarkers at default thresholds)")
  }
  println("")

  // -----------------------------------------------------------------------
  // 7. Print all per-feature results
  // -----------------------------------------------------------------------
  println("7. All per-feature results:")
  for r in output.results() {
    let sig_label = if r.is_significant() { " *" } else { "" }
    println("   \{r.feature_name()}: class=\{r.class_name()}, KW p=\{r.kw_pvalue()}, LDA=\{r.lda_score()}\{sig_label}")
  }
  println("")

  // -----------------------------------------------------------------------
  // 8. Run with sample data
  // -----------------------------------------------------------------------
  println("8. Using built-in sample data:")
  let sample_input = @src.lefse_sample_data()
  let sample_output = @src.lefse_run(sample_input)
  println("   Features: \{sample_input.feature_names.length()}")
  println("   Samples: \{sample_input.sample_names.length()}")
  println("   Significant biomarkers: \{sample_output.n_significant()}")
  for r in sample_output.results() {
    if r.is_significant() {
      println("   -> \{r.feature_name()} (\{r.class_name()}): LDA=\{r.lda_score()}, p=\{r.kw_pvalue()}")
    }
  }
  println("")

  // -----------------------------------------------------------------------
  // 9. Print result summary
  // -----------------------------------------------------------------------
  println("9. Result summary:")
  let summary = @src.lefse_result_summary(sample_output)
  println(summary)

  println("=== Demo Complete ===")
}
