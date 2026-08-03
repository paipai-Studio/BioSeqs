///|
/// ANCOMBC module demo - Microbiome differential abundance with bias correction.
///
/// Run: moon run examples/ancombc_demo
fn main {
  println("=== ANCOMBC Module Demo ===")
  println("Analysis of Composition of Microbiomes with Bias Correction")
  println("Reference: Lin & Peddada (2020) Nat Commun 11:5144")
  println("")

  // ------------------------------------------------------------------
  // 1. Load sample data
  // ------------------------------------------------------------------
  println("=== 1. Load Sample Data ===")
  let data = @bio.ancombc_sample_data()
  println("Samples: " + data.n_samples().to_string())
  println("Features: " + data.n_features().to_string())
  let groups = data.groups()
  println("Groups: " + groups.length().to_string())
  for g in groups {
    println("  - " + g)
  }
  println("")
  let lib_sizes = data.library_sizes()
  println("Library sizes:")
  for i in 0..<lib_sizes.length() {
    println("  " + data.sample_names()[i] + ": " + lib_sizes[i].to_string())
  }
  println("")

  // ------------------------------------------------------------------
  // 2. Sampling fractions
  // ------------------------------------------------------------------
  println("=== 2. Sampling Fractions ===")
  let fracs = @bio.ancombc_sampling_fractions(data)
  println("Sampling fractions:")
  for i in 0..<fracs.length() {
    println("  " + data.sample_names()[i] + ": " + fracs[i].to_string())
  }
  let mut total = 0.0
  for f in fracs {
    total = total + f
  }
  println("Sum (should be ~1.0): " + total.to_string())
  println("")

  // ------------------------------------------------------------------
  // 3. Reference feature selection
  // ------------------------------------------------------------------
  println("=== 3. Reference Feature ===")
  let ref_idx = @bio.ancombc_reference_feature(data)
  println("Reference feature index: " + ref_idx.to_string())
  println("Reference feature name: " + data.feature_names()[ref_idx])
  println("")

  // ------------------------------------------------------------------
  // 4. Differential abundance test
  // ------------------------------------------------------------------
  println("=== 4. Differential Abundance Test ===")
  let result = @bio.ancombc_test(data, "control", "treatment", 0.05)
  println("Comparison: " + result.group1() + " vs " + result.group2())
  println("Reference feature: " + result.reference_feature())
  println("Significant features: " + result.n_significant().to_string())
  println("")
  println("Per-feature results:")
  for r in result.results() {
    println(
      "  " +
        r.feature() +
        ": logFC=" +
        r.log_fold_change().to_string() +
        ", W=" +
        r.w_stat().to_string() +
        ", p=" +
        r.p_value().to_string() +
        ", q=" +
        r.q_value().to_string() +
        ", sig=" +
        r.significant().to_string(),
    )
  }
  println("")

  // ------------------------------------------------------------------
  // 5. Full summary
  // ------------------------------------------------------------------
  println("=== 5. Full Summary ===")
  let summary = @bio.ancombc_result_summary(result)
  println(summary)

  println("=== Demo complete ===")
}
