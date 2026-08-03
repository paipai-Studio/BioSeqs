///|
/// ESTIMATE module demo - tumor microenvironment scoring.
///
/// Run: moon run examples/estimate_score_demo
fn main {
  println("=== ESTIMATE Module Demo ===")
  println("Estimation of STromal and Immune cells in MAlignant Tumours")
  println("Reference: Yoshihara K et al. (2013) Nature Communications 4:2612")
  println("")

  // ------------------------------------------------------------------
  // 1. Display built-in gene signatures
  // ------------------------------------------------------------------
  println("=== 1. Built-in gene signatures ===")
  let sig = @bio.est_default_signature()
  println("Stromal signature genes: \{sig.stromal_genes().length()}")
  println("Immune signature genes: \{sig.immune_genes().length()}")
  let stromal_preview = sig.stromal_genes()
  let immune_preview = sig.immune_genes()
  // Show first 5 of each
  let mut i = 0
  let mut st_preview = ""
  while i < 5 && i < stromal_preview.length() {
    if i > 0 { st_preview = st_preview + ", " }
    st_preview = st_preview + stromal_preview[i]
    i = i + 1
  }
  i = 0
  let mut im_preview = ""
  while i < 5 && i < immune_preview.length() {
    if i > 0 { im_preview = im_preview + ", " }
    im_preview = im_preview + immune_preview[i]
    i = i + 1
  }
  println("Stromal (first 5): \{st_preview}")
  println("Immune  (first 5): \{im_preview}")
  println("")

  // ------------------------------------------------------------------
  // 2. Load sample expression data
  // ------------------------------------------------------------------
  println("=== 2. Load sample expression data ===")
  let expr = @bio.est_sample_data()
  println("Genes: \{expr.gene_names().length()}")
  println("Samples: \{expr.sample_names().length()}")
  println("Sample names: \{expr.sample_names()}")
  println("")

  // ------------------------------------------------------------------
  // 3. Run ESTIMATE
  // ------------------------------------------------------------------
  println("=== 3. Run ESTIMATE scoring ===")
  let result = @bio.est_run(expr)
  println("Method: \{result.method_name()}")
  println("")

  // ------------------------------------------------------------------
  // 4. Display per-sample scores
  // ------------------------------------------------------------------
  println("=== 4. Per-sample ESTIMATE scores ===")
  for s in result.scores() {
    println("--- Sample: \{s.sample_id()} ---")
    println("  Stromal score: \{s.stromal_score()}")
    println("  Immune score:  \{s.immune_score()}")
    println("  ESTIMATE score: \{s.estimate_score()}")
    if s.tumor_purity().is_nan() {
      println("  Tumor purity:  NA (undefined)")
    } else {
      println("  Tumor purity:  \{s.tumor_purity()}")
    }
    println("")
  }

  // ------------------------------------------------------------------
  // 5. Compare infiltration profiles
  // ------------------------------------------------------------------
  println("=== 5. Infiltration profile comparison ===")
  let high_stromal = result.get_score("HighStromal")
  let high_immune = result.get_score("HighImmune")
  let high_both = result.get_score("HighBoth")
  let low_both = result.get_score("LowBoth")
  println("Sample        Stromal  Immune   ESTIMATE  Purity")
  println(
    "HighStromal   \{high_stromal.stromal_score()}    \{high_stromal.immune_score()}      \{high_stromal.estimate_score()}      \{purity_str(high_stromal.tumor_purity())}",
  )
  println(
    "HighImmune    \{high_immune.stromal_score()}    \{high_immune.immune_score()}      \{high_immune.estimate_score()}      \{purity_str(high_immune.tumor_purity())}",
  )
  println(
    "HighBoth      \{high_both.stromal_score()}    \{high_both.immune_score()}      \{high_both.estimate_score()}      \{purity_str(high_both.tumor_purity())}",
  )
  println(
    "LowBoth       \{low_both.stromal_score()}    \{low_both.immune_score()}      \{low_both.estimate_score()}      \{purity_str(low_both.tumor_purity())}",
  )
  println("")

  // ------------------------------------------------------------------
  // 6. Interpretation
  // ------------------------------------------------------------------
  println("=== 6. Interpretation ===")
  if high_stromal.stromal_score() > low_both.stromal_score() {
    println("-> HighStromal has higher stromal score than LowBoth (expected)")
  }
  if high_immune.immune_score() > low_both.immune_score() {
    println("-> HighImmune has higher immune score than LowBoth (expected)")
  }
  if high_both.estimate_score() > low_both.estimate_score() {
    println("-> HighBoth has higher ESTIMATE score than LowBoth (expected)")
  }
  if low_both.stromal_score() < 0.0 && low_both.immune_score() < 0.0 {
    println("-> LowBoth has negative scores (low infiltration, tumor-dominant)")
  }
  if !low_both.tumor_purity().is_nan() &&
    !high_both.tumor_purity().is_nan() &&
    low_both.tumor_purity() >= high_both.tumor_purity() {
    println("-> LowBoth has higher tumor purity (less infiltration)")
  }
  println("")

  // ------------------------------------------------------------------
  // 7. Custom signature example
  // ------------------------------------------------------------------
  println("=== 7. Custom gene signature example ===")
  let custom_sig = @bio.EstSignature::new(
    ["G1", "G2", "G3"],
    ["G3", "G4", "G5"],
  )
  let custom_expr = @bio.EstExpression::new(
    ["G1", "G2", "G3", "G4", "G5"],
    ["S1", "S2"],
    [
      [10.0, 1.0], // G1 - stromal, high in S1
      [8.0, 1.0], // G2 - stromal, high in S1
      [5.0, 5.0], // G3 - shared
      [1.0, 9.0], // G4 - immune, high in S2
      [1.0, 8.0], // G5 - immune, high in S2
    ],
  )
  let custom_result = @bio.est_run(custom_expr, signature=custom_sig)
  println("Custom signature result:")
  for s in custom_result.scores() {
    println(
      "  \{s.sample_id()}: stromal=\{s.stromal_score()}, immune=\{s.immune_score()}, estimate=\{s.estimate_score()}",
    )
  }
  println("")

  // ------------------------------------------------------------------
  // 8. Full result table
  // ------------------------------------------------------------------
  println("=== 8. Full result table ===")
  println(result.to_string())

  println("=== Demo Complete ===")
}

///|
/// Format tumor purity for display.
fn purity_str(p : Double) -> String {
  if p.is_nan() {
    "NA"
  } else {
    p.to_string()
  }
}
