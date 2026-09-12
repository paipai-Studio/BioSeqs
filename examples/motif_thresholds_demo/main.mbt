///|
fn fmt4(value : Double) -> String {
  let scaled = (value * 10000.0).round() / 10000.0
  scaled.to_string()
}

///|
fn main raise {
  println("=== Bio.motifs.thresholds (ScoreDistribution) Demo ===")
  // Eight-position motif: nucleotide counts (rows A, C, G, T) from 10 sites.
  let counts = [
    [7, 1, 7, 2, 6, 1, 5, 3],
    [1, 6, 1, 1, 2, 7, 2, 4],
    [1, 2, 1, 6, 1, 1, 2, 2],
    [1, 1, 1, 1, 1, 1, 1, 1],
  ]

  println("\n1. Build the discrete score distribution (precision=1000/column)")
  let dist = @src.MotifScoreDistribution::from_counts(counts)
  println(
    "   grid: n=\{dist.n_points()} bins, min score=\{fmt4(dist.min_score())}, step=\{fmt4(dist.step())}"
  )
  println("   information content IC=\{fmt4(dist.information_content())} bits")
  let mo = dist.motif_density()
  let bg = dist.background_density()
  let mut mo_sum = 0.0
  let mut bg_sum = 0.0
  for i in 0..<mo.length() {
    mo_sum += mo[i]
    bg_sum += bg[i]
  }
  println("   density mass: motif=\{mo_sum} background=\{bg_sum}")

  println("\n2. Thresholds for target false positive rates (background hits)")
  for fpr in [0.001, 0.01, 0.05, 0.1] {
    println(
      "   FPR <= \{fpr}  ->  log-odds threshold \{fmt4(dist.threshold_fpr(fpr))}"
    )
  }

  println("\n3. Thresholds for target false negative rates (missed motifs)")
  for fnr in [0.001, 0.01, 0.05, 0.1] {
    println(
      "   FNR <= \{fnr}  ->  log-odds threshold \{fmt4(dist.threshold_fnr(fnr))}"
    )
  }

  println("\n4. Balanced error-rate thresholds and the patser criterion")
  for proportion in [0.5, 1.0, 2.0] {
    let (threshold, actual_fpr) = dist.threshold_balanced_with_rate(
      rate_proportion=proportion
    )
    println(
      "   FNR = FPR * \{proportion}: threshold=\{fmt4(threshold)} (actual FPR=\{fmt4(actual_fpr)})"
    )
  }
  println("   patser (FPR = 2^-IC): \{fmt4(dist.threshold_patser())}")

  println("\n5. Effect of an AT-rich background (A/T 30%, C/G 20%)")
  let biased = @src.MotifScoreDistribution::from_counts(
    counts,
    background=[0.30, 0.20, 0.20, 0.30]
  )
  println("   IC=\{fmt4(biased.information_content())} bits")
  println(
    "   balanced threshold=\{fmt4(biased.threshold_balanced())} (uniform: \{fmt4(dist.threshold_balanced())})"
  )
  println(
    "   patser threshold=\{fmt4(biased.threshold_patser())} (uniform: \{fmt4(dist.threshold_patser())})"
  )

  println("\n6. Pseudocounts rescue motifs with unobserved letters")
  let sparse = [
    [3, 1, 1, 2],
    [1, 0, 2, 1],
    [0, 2, 0, 1],
    [0, 1, 1, 0],
  ]
  let rejected = try {
    let _ = @src.MotifScoreDistribution::from_counts(sparse)
    false
  } catch {
    _ => true
  }
  println("   zero-count motif without pseudocount rejected: \{rejected}")
  let smoothed = @src.MotifScoreDistribution::from_counts(
    sparse,
    pseudocount=0.5,
    precision=500
  )
  println(
    "   with pseudocount 0.5: FPR 1% threshold=\{fmt4(smoothed.threshold_fpr(0.01))}"
  )

  println("\n7. Invalid matrices are rejected")
  let ragged = [
    [1.0, 2.0],
    [1.0],
    [1.0, 2.0],
    [1.0, 2.0],
  ]
  let bad_shape = try {
    let _ = @src.MotifScoreDistribution::from_pssm(ragged)
    false
  } catch {
    _ => true
  }
  println("   ragged PSSM rejected: \{bad_shape}")

  println("\n=== Demo Complete ===")
}
