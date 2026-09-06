///|
/// QValue Module Demo.
/// Demonstrates Storey's q-value FDR correction method,
/// pi0 estimation, bootstrap inference, and comparison
/// with the Benjamini-Hochberg procedure.

/// Format a Double to a String with a given number of decimal places.
fn fmt(value : Double, decimals : Int) -> String {
  let s = value.to_string()
  let parts = s.split(".").to_array()
  let int_part = if parts.length() > 0 { parts[0].to_owned() } else { "0" }
  let frac_part = if parts.length() > 1 { parts[1].to_owned() } else { "" }

  let mut padded_frac = frac_part
  if padded_frac.length() < decimals {
    padded_frac = padded_frac + "0".repeat(decimals - padded_frac.length())
  } else if padded_frac.length() > decimals {
    padded_frac = padded_frac[0:decimals].to_owned()
  }

  int_part + "." + padded_frac
}

/// Simple Benjamini-Hochberg FDR correction for comparison.
fn benjamini_hochberg(p_values : Array[Double]) -> Array[Double] {
  let n = p_values.length()
  if n == 0 {
    return []
  }

  let indexed : Array[Array[Double]] = Array::new()
  let mut i = 0
  while i < n {
    indexed.push([p_values[i], i.to_double()])
    i = i + 1
  }

  indexed.sort_by(fn(a : Array[Double], b : Array[Double]) -> Int {
    if a[0] < b[0] { -1 } else if a[0] > b[0] { 1 } else { 0 }
  })

  let adj_sorted : Array[Double] = Array::make(n, 1.0)
  let mut min_adj = 1.0
  let mut k = n - 1
  while k >= 0 {
    let p = indexed[k][0]
    let rank = (k + 1).to_double()
    let adj_val = p * n.to_double() / rank
    if adj_val < min_adj {
      min_adj = adj_val
    }
    adj_sorted[k] = if min_adj > 1.0 { 1.0 } else { min_adj }
    k = k - 1
  }

  let result : Array[Double] = Array::make(n, 0.0)
  let mut j = 0
  while j < n {
    let orig_idx = indexed[j][1].to_int()
    result[orig_idx] = adj_sorted[j]
    j = j + 1
  }

  result
}

fn main raise {
  println("=== QValue Module Demo ===\n")

  let p_values = [
    0.0001, 0.0005, 0.001, 0.002, 0.005,
    0.01, 0.02, 0.03, 0.04, 0.05,
    0.06, 0.07, 0.08, 0.09, 0.10,
    0.15, 0.20, 0.25, 0.30, 0.35,
    0.40, 0.45, 0.50, 0.55, 0.60,
    0.65, 0.70, 0.75, 0.80, 0.85,
    0.90, 0.95, 0.98,
  ]

  println("Input p-values (\{p_values.length()} tests):")
  let mut di = 0
  while di < p_values.length() {
    println("  p[\{di}] = \{p_values[di]}")
    di = di + 1
  }

  // ---------------------------------------------------------------------------
  // 1. Basic q-value calculation
  // ---------------------------------------------------------------------------
  println("\n1. Basic q-value calculation (Storey's method):")
  let result = @src.qvalue(p_values)
  let q_vals = result.q_values
  println("  Computed q-values:")
  let mut i = 0
  while i < p_values.length() {
    println("    p = \{fmt(p_values[i], 4)}  ->  q = \{fmt(q_vals[i], 4)}")
    i = i + 1
  }
  println("  Estimated pi0 = \{result.pi0}")
  println("  Selected lambda = \{result.lambda}")
  println("  Bootstrap gamma = \{result.gamma}")

  // ---------------------------------------------------------------------------
  // 2. pi0 estimation
  // ---------------------------------------------------------------------------
  println("\n2. pi0 estimation at different lambda values:")
  let lambdas = [0.05, 0.10, 0.20, 0.30, 0.40, 0.50, 0.60, 0.70, 0.80, 0.90]
  let mut li = 0
  while li < lambdas.length() {
    let pi0_est = @src.qvalue_estimate_pi0(p_values, lambdas[li])
    println("    lambda = \{fmt(lambdas[li], 2)}  ->  pi0 = \{fmt(pi0_est, 4)}")
    li = li + 1
  }
  println("  Final pi0 (auto-selected lambda = \{result.lambda}): \{result.pi0}")

  // ---------------------------------------------------------------------------
  // 3. Bootstrap pi0 estimation
  // ---------------------------------------------------------------------------
  println("\n3. Bootstrap pi0 estimation:")
  let boot_pi0_100 = @src.bootstrap_pi0(p_values, 100)
  println("  Bootstrap pi0 (100 iterations): \{fmt(boot_pi0_100, 4)}")
  let boot_pi0_500 = @src.bootstrap_pi0(p_values, 500)
  println("  Bootstrap pi0 (500 iterations): \{fmt(boot_pi0_500, 4)}")
  println("  Analytical pi0: \{fmt(result.pi0, 4)}")

  // ---------------------------------------------------------------------------
  // 4. Significance testing
  // ---------------------------------------------------------------------------
  println("\n4. Significance testing:")
  println("  Significance at q-value <= 0.05:")
  let sig_05 = result.significance
  let mut sig_count_05 = 0
  let mut si = 0
  while si < sig_05.length() {
    if sig_05[si] {
      sig_count_05 = sig_count_05 + 1
    }
    si = si + 1
  }
  println("    Significant: \{sig_count_05} / \{p_values.length()}")
  println("    Not significant: \{p_values.length() - sig_count_05} / \{p_values.length()}")

  println("\n  Significance at q-value <= 0.01:")
  let sig_01 = @src.qvalue_significance(result, 0.01)
  let mut sig_count_01 = 0
  let mut sj = 0
  while sj < sig_01.length() {
    if sig_01[sj] {
      sig_count_01 = sig_count_01 + 1
    }
    sj = sj + 1
  }
  println("    Significant: \{sig_count_01} / \{p_values.length()}")

  println("\n  Significance at q-value <= 0.10:")
  let sig_10 = @src.qvalue_significance(result, 0.10)
  let mut sig_count_10 = 0
  let mut sk = 0
  while sk < sig_10.length() {
    if sig_10[sk] {
      sig_count_10 = sig_count_10 + 1
    }
    sk = sk + 1
  }
  println("    Significant: \{sig_count_10} / \{p_values.length()}")

  // Show significance array details
  println("\n  Per-test significance (alpha=0.05):")
  let mut sti = 0
  while sti < p_values.length() {
    println("    Test \{sti}: p=\{fmt(p_values[sti], 4)}, q=\{fmt(q_vals[sti], 4)}, significant=\{sig_05[sti]}")
    sti = sti + 1
  }

  // ---------------------------------------------------------------------------
  // 5. Comparison with Benjamini-Hochberg
  // ---------------------------------------------------------------------------
  println("\n5. Comparison with Benjamini-Hochberg:")
  let bh_vals = benjamini_hochberg(p_values)
  println("  p-value       q-value (Storey)  BH-adjusted")
  println("  --------------------------------------------")
  let mut ci = 0
  while ci < p_values.length() {
    let p_str = fmt(p_values[ci], 4)
    let q_str = fmt(q_vals[ci], 4)
    let bh_str = fmt(bh_vals[ci], 4)
    println("  \{p_str}      \{q_str}          \{bh_str}")
    ci = ci + 1
  }

  let storey_sig_count = sig_count_05
  let mut bh_sig_count = 0
  let mut bi = 0
  while bi < bh_vals.length() {
    if bh_vals[bi] <= 0.05 {
      bh_sig_count = bh_sig_count + 1
    }
    bi = bi + 1
  }
  println("\n  Summary at FDR = 0.05:")
  println("    Storey's q-value: \{storey_sig_count} significant tests")
  println("    Benjamini-Hochberg: \{bh_sig_count} significant tests")

  // Show storey_fdr convenience function
  println("\n  storey_fdr() convenience function:")
  let storey_fdr_vals = @src.storey_fdr(p_values)
  let mut sf_ok = true
  let mut sfi = 0
  while sfi < q_vals.length() {
    let diff = storey_fdr_vals[sfi] - q_vals[sfi]
    let neg_diff = 0.0 - diff
    if diff > 0.0000000001 || neg_diff > 0.0000000001 {
      sf_ok = false
    }
    sfi = sfi + 1
  }
  println("    Results match qvalue(): \{sf_ok}")

  // Show qvalue_with_lambda
  println("\n  qvalue_with_lambda() with lambda=0.5:")
  let result_custom = @src.qvalue_with_lambda(p_values, 0.5)
  println("    pi0 = \{result_custom.pi0}")
  println("    lambda = \{result_custom.lambda}")
  let custom_q = result_custom.q_values
  let mut cj = 0
  while cj < p_values.length() {
    println("      p = \{fmt(p_values[cj], 4)}  ->  q = \{fmt(custom_q[cj], 4)}")
    cj = cj + 1
  }

  // ---------------------------------------------------------------------------
  // 6. Bioconductor pi0est(): df=3 smoothing spline vs closed-form bootstrap
  // ---------------------------------------------------------------------------
  println("\n6. pi0est() on a 50/50 signal/null mixture (Bioconductor parity):")
  let n_mix = 200
  let mix_p = Array::make(n_mix, 0.0)
  let mut mi = 0
  while mi < 100 {
    mix_p[mi] = (mi.to_double() + 0.5) / 100.0 * 0.002
    mi = mi + 1
  }
  mi = 0
  while mi < 100 {
    mix_p[100 + mi] = 0.01 + (mi.to_double() + 0.5) / 100.0 * 0.98
    mi = mi + 1
  }
  let pi0_smooth = @src.qvalue_pi0est(mix_p)
  let pi0_boot = @src.qvalue_pi0est(mix_p, method="bootstrap")
  println("  smoother  pi0 = \{fmt(pi0_smooth.pi0(), 4)} (cubic spline, df = 3)")
  println("  bootstrap pi0 = \{fmt(pi0_boot.pi0(), 4)} (closed-form MSE, deterministic)")
  println("  pi0(lambda) path at selected lambdas:")
  let path = pi0_smooth.pi0_lambda()
  let lams_out = pi0_smooth.lambdas()
  let mut pi4 = 0
  while pi4 < 4 {
    println("    lambda = \{fmt(lams_out[pi4], 2)} -> pi0 = \{fmt(path[pi4], 4)}")
    pi4 = pi4 + 1
  }
  println("    lambda = \{fmt(lams_out[18], 2)} -> pi0 = \{fmt(path[18], 4)} (spline endpoint)")

  // ---------------------------------------------------------------------------
  // 7. Local FDR (lfdr): probit transform + Gaussian KDE
  // ---------------------------------------------------------------------------
  println("\n7. Local FDR (lfdr) on the mixture (probit + KDE, monotone):")
  let lfdr_vals = @src.qvalue_lfdr(mix_p)
  let show_idx = [0, 20, 50, 99, 100, 150, 199]
  let mut shi = 0
  while shi < show_idx.length() {
    let idx = show_idx[shi]
    println("    p = \{fmt(mix_p[idx], 4)} -> lfdr = \{fmt(lfdr_vals[idx], 4)}")
    shi = shi + 1
  }

  // ---------------------------------------------------------------------------
  // 8. Empirical p-values (empPvals): pooled and test-specific
  // ---------------------------------------------------------------------------
  println("\n8. Empirical p-values (empPvals) from test statistics:")
  let ep = @src.empirical_pvalues([2.0, 5.0], [1.0, 2.0, 3.0, 4.0])
  println("  pooled:   stat = [2, 5], stat0 = [1, 2, 3, 4]")
  println("            p = [\{fmt(ep[0], 3)}, \{fmt(ep[1], 3)}] (ties: observed first; floor 1/m0 = 0.25)")
  let epr = @src.empirical_pvalues_rows([2.0, 5.0], [
    [1.0, 2.0, 3.0, 4.0], [1.0, 2.0, 3.0, 4.0]
  ])
  println("  per-test: each statistic vs its own 4 nulls")
  println("            p = [\{fmt(epr[0], 3)}, \{fmt(epr[1], 3)}] (3/4 = 0.75; 0/4 floored to 1/4 = 0.25)")

  // ---------------------------------------------------------------------------
  // 9. Full qvalue() analysis: pfdr correction and fdr.level significance
  // ---------------------------------------------------------------------------
  println("\n9. Full qvalue() analysis (Bioconductor parity):")
  let qa = @src.qvalue_analyze(mix_p, fdr_level=Some(0.05))
  let qaf = @src.qvalue_analyze(mix_p, pfdr=true, fdr_level=Some(0.05))
  let mut n_sig = 0
  let mut n_sig_pfdr = 0
  let mut qi = 0
  while qi < n_mix {
    if qa.significant()[qi] {
      n_sig = n_sig + 1
    }
    if qaf.significant()[qi] {
      n_sig_pfdr = n_sig_pfdr + 1
    }
    qi = qi + 1
  }
  println("  pi0 = \{fmt(qa.pi0(), 4)}, method = \{qa.pi0_method()}")
  println("  calls at FDR <= 0.05: \{n_sig} (plain q-value), \{n_sig_pfdr} (pFDR corrected)")
  println("  q at smallest p: plain = \{fmt(qa.q_values()[0], 4)}, pfdr = \{fmt(qaf.q_values()[0], 4)}")
  println(
    "  q at largest  p: plain = \{fmt(qa.q_values()[n_mix - 1], 4)}, pfdr = \{fmt(qaf.q_values()[n_mix - 1], 4)}",
  )

  println("\n=== QValue Module Demo Complete ===")
}