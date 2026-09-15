///|
/// ashr demo - Adaptive shrinkage of noisy effect sizes.
/// Inspired by Bioconductor {ashr}.

fn main {
  println("=== ashr Demo: Adaptive Shrinkage (Empirical Bayes) ===")
  println("")

  let (estimates, ses) = @src.ashr_example_data()

  // Step 1: the raw estimates
  println("--- Step 1: Raw Effect Estimates ---")
  let mut i = 0
  while i < estimates.length() {
    println(
      "  effect " + i.to_string() + ": beta_hat = " +
      estimates[i].to_string() + "  (se = " + ses[i].to_string() + ")",
    )
    i = i + 1
  }
  println("")

  // Step 2: fit the unimodal prior by EM
  println("--- Step 2: Fit Unimodal Prior (EM) ---")
  let fit = @src.ashr_shrink(estimates, ses)
  println("  iterations: " + fit.iterations.to_string())
  println("  log-likelihood: " + fit.log_likelihood.to_string())
  let mode = @src.ashr_mode_component(fit)
  println(
    "  dominant prior component: sd = " + fit.prior_grid[mode].to_string() +
    " (weight " + fit.prior_weights[mode].to_string() + ")",
  )
  println("")

  // Step 3: shrinkage results
  println("--- Step 3: Shrunk Estimates vs Raw ---")
  i = 0
  while i < fit.n {
    println(
      "  effect " + i.to_string() + ": " + estimates[i].to_string() +
      " -> " + fit.posterior_mean[i].to_string() +
      "   [post.sd = " + fit.posterior_sd[i].to_string() +
      ", lfsr = " + fit.lfsr[i].to_string() + "]",
    )
    i = i + 1
  }
  println("")

  // Step 4: significance by local false sign rate
  println("--- Step 4: Significant Hits (lfsr < 0.05) ---")
  println("  " + @src.ashr_count_significant(fit, 0.05).to_string() + " of " + fit.n.to_string())
  println("")

  // Step 5: shrinkage factors
  println("--- Step 5: Shrinkage Factors ---")
  let factors = @src.ashr_shrinkage_factor(fit, estimates)
  i = 0
  while i < factors.length() {
    println(
      "  effect " + i.to_string() + ": factor = " + factors[i].to_string(),
    )
    i = i + 1
  }
  println("")

  // Step 6: z-score convenience wrapper
  println("--- Step 6: Z-score Shrinkage ---")
  let zfit = @src.ashr_shrink_z([3.0, -2.5, 0.4, -0.2])
  println("  significant z (lfsr < 0.05): " + @src.ashr_count_significant(zfit, 0.05).to_string())
  println(@src.ashr_summary(zfit))

  println("")
  println("=== ashr Demo Complete ===")
}
