///|
/// IHW (Independent Hypothesis Weighting) Demo.
/// Demonstrates basic IHW correction, custom configuration,
/// multiple covariates, Storey-style IHW, and comparison with
/// the Bonferroni baseline.

/// Helper: format an array of doubles into a short string
/// with up to `n` entries followed by "...".
fn format_array(arr : Array[Double], n : Int) -> String {
  let mut s = "["
  let mut i = 0
  let limit = if n < arr.length() { n } else { arr.length() }
  while i < limit {
    if i > 0 {
      s = s + ", "
    }
    s = s + arr[i].to_string()
    i = i + 1
  }
  if arr.length() > n {
    s = s + ", ..."
  }
  s = s + "]"
  s
}

/// Helper: count significant calls at level alpha in an
/// array of adjusted p-values.
fn count_significant(adj : Array[Double], alpha : Double) -> Int {
  let mut c = 0
  let mut i = 0
  while i < adj.length() {
    if adj[i] <= alpha {
      c = c + 1
    }
    i = i + 1
  }
  c
}

fn main {
  println("=== IHW (Independent Hypothesis Weighting) Demo ===\n")

  // ------------------------------------------------------------------
  // 1. Basic IHW correction with covariates
  // ------------------------------------------------------------------
  println("1. Basic IHW correction with covariates")
  println("---------------------------------------")
  let p_values : Array[Double] = [
    0.0001, 0.0020, 0.0100, 0.0300, 0.0400,
    0.0450, 0.0800, 0.1200, 0.3000, 0.5000,
    0.7000, 0.8500, 0.9000, 0.9500, 0.9900,
  ]
  let covariates : Array[Double] = [
    0.95, 0.90, 0.85, 0.70, 0.65,
    0.55, 0.40, 0.30, 0.25, 0.20,
    0.15, 0.10, 0.08, 0.05, 0.02,
  ]
  println("   Input p-values    : \{format_array(p_values, 5)}")
  println("   Input covariates  : \{format_array(covariates, 5)}")

  let result1 = @src.ihw(p_values, covariates, 0.05)
  println("   IHW attempts      : \{result1.n_attempts()}")
  println("   IHW weights       : \{format_array(result1.weights(), 5)}")
  println("   Adjusted p-values  : \{format_array(result1.adjusted_p_values(), 5)}")
  let n_sig1 = count_significant(result1.adjusted_p_values(), 0.05)
  println("   Significant (alpha=0.05): \{n_sig1}")

  // ------------------------------------------------------------------
  // 2. IHW with custom configuration (global vs local)
  // ------------------------------------------------------------------
  println("\n2. IHW with custom configuration (global and local)")
  println("---------------------------------------------------")

  let config_global = @src.IHWConfig::new(alpha=0.05, n_attempts=5, scale_type="global")
  let result_global = @src.ihw_with_config(p_values, covariates, config_global)
  println("   [global]  attempts     : \{result_global.n_attempts()}")
  println("   [global]  weights      : \{format_array(result_global.weights(), 5)}")
  println("   [global]  adj.p-values : \{format_array(result_global.adjusted_p_values(), 5)}")
  println("   [global]  significant  : \{count_significant(result_global.adjusted_p_values(), 0.05)}")

  let config_local = @src.IHWConfig::new(alpha=0.05, n_attempts=10, scale_type="local")
  let result_local = @src.ihw_with_config(p_values, covariates, config_local)
  println("   [local]   attempts     : \{result_local.n_attempts()}")
  println("   [local]   weights      : \{format_array(result_local.weights(), 5)}")
  println("   [local]   adj.p-values : \{format_array(result_local.adjusted_p_values(), 5)}")
  println("   [local]   significant  : \{count_significant(result_local.adjusted_p_values(), 0.05)}")

  // ------------------------------------------------------------------
  // 3. Multiple covariate IHW
  // ------------------------------------------------------------------
  println("\n3. Multiple covariate IHW")
  println("-------------------------")
  let cov1 : Array[Double] = [
    0.90, 0.85, 0.80, 0.70, 0.60,
    0.55, 0.45, 0.35, 0.25, 0.15,
    0.10, 0.08, 0.05, 0.03, 0.02,
  ]
  let cov2 : Array[Double] = [
    0.80, 0.82, 0.78, 0.65, 0.58,
    0.50, 0.42, 0.33, 0.22, 0.18,
    0.12, 0.09, 0.07, 0.04, 0.03,
  ]
  let cov3 : Array[Double] = [
    0.70, 0.75, 0.72, 0.68, 0.62,
    0.50, 0.40, 0.30, 0.20, 0.12,
    0.08, 0.06, 0.04, 0.03, 0.02,
  ]
  let multi_covs : Array[Array[Double]] = [cov1, cov2, cov3]
  let result_multi = @src.ihw_multiple(p_values, multi_covs, 0.05)
  println("   Covariates provided : 3 (cov1, cov2, cov3)")
  println("   Combined weights    : \{format_array(result_multi.weights(), 5)}")
  println("   Adjusted p-values   : \{format_array(result_multi.adjusted_p_values(), 5)}")
  println("   Significant         : \{count_significant(result_multi.adjusted_p_values(), 0.05)}")

  // ------------------------------------------------------------------
  // 4. Storey-style IHW
  // ------------------------------------------------------------------
  println("\n4. Storey-style IHW")
  println("-------------------")
  let storey_adj = @src.storey_ihw(p_values, covariates)
  println("   Storey adjusted p-values: \{format_array(storey_adj, 5)}")
  println("   Significant (alpha=0.05): \{count_significant(storey_adj, 0.05)}")

  // ------------------------------------------------------------------
  // 5. Compare IHW weights with Bonferroni baseline
  // ------------------------------------------------------------------
  println("\n5. Compare IHW with Bonferroni baseline")
  println("--------------------------------------")
  let bonf_adj = @src.bonferroni_ihw(p_values)
  println("   Bonferroni adj.p-values : \{format_array(bonf_adj, 5)}")
  println("   Bonferroni significant  : \{count_significant(bonf_adj, 0.05)}")
  println("   IHW adj.p-values        : \{format_array(result1.adjusted_p_values(), 5)}")
  println("   IHW significant          : \{count_significant(result1.adjusted_p_values(), 0.05)}")
  println("   Storey adj.p-values     : \{format_array(storey_adj, 5)}")
  println("   Storey significant      : \{count_significant(storey_adj, 0.05)}")

  println("\n   IHW vs Bonferroni weights:")
  let mut i = 0
  while i < result1.weights().length() && i < 5 {
    let w = result1.weights()[i]
    let equal = 1.0
    let marker = if w > equal { "higher" } else if w < equal { "lower" } else { "equal" }
    println(
      "      hypothesis \{i+1}: IHW weight = \{w}, Bonferroni weight = \{equal} (IHW \{marker})",
    )
    i = i + 1
  }

  println("\n=== IHW Demo Complete ===")
}
