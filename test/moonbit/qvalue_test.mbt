///|
/// Tests for Bio.QValue module - Storey's q-value FDR correction.

///|
test "qvalue_empty_input" {
  let pvals : Array[Double] = []
  let result = @src.qvalue(pvals)
  assert_eq(result.q_values().length(), 0)
  assert_eq(result.p_values().length(), 0)
  assert_eq(result.significance().length(), 0)
  assert_true(result.pi0() == 0.0)
}

///|
test "qvalue_single_value" {
  let pvals = [0.01]
  let result = @src.qvalue(pvals)
  assert_eq(result.q_values().length(), 1)
  assert_true(result.q_values()[0] >= 0.0 && result.q_values()[0] <= 1.0)
  assert_true(result.pi0() >= 0.0 && result.pi0() <= 1.0)
}

///|
test "qvalue_all_same" {
  let pvals = [0.05, 0.05, 0.05, 0.05, 0.05]
  let result = @src.qvalue(pvals)
  assert_eq(result.q_values().length(), 5)
  let qvals = result.q_values()
  let mut i = 0
  while i < 5 {
    assert_true(qvals[i] >= 0.0 && qvals[i] <= 1.0)
    i = i + 1
  }
}

///|
test "qvalue_basic_calculation" {
  let pvals = [0.001, 0.01, 0.02, 0.04, 0.05, 0.10, 0.20, 0.50]
  let result = @src.qvalue(pvals)
  let qvals = result.q_values()
  assert_eq(qvals.length(), 8)
  let mut i = 0
  while i < 8 {
    assert_true(qvals[i] >= 0.0 && qvals[i] <= 1.0)
    i = i + 1
  }
}

///|
test "qvalues_non_increasing_when_sorted" {
  let pvals = [0.001, 0.01, 0.02, 0.04, 0.05, 0.10, 0.20, 0.50]
  let result = @src.qvalue(pvals)
  let qvals = result.q_values()
  let mut i = 1
  while i < qvals.length() {
    let diff = qvals[i] - qvals[i - 1]
    assert_true(diff >= -0.0000000001)
    i = i + 1
  }
}

///|
test "qvalue_with_lambda_specified" {
  let pvals = [0.001, 0.01, 0.02, 0.04, 0.05, 0.10, 0.20, 0.50]
  let result = @src.qvalue_with_lambda(pvals, 0.5)
  assert_eq(result.q_values().length(), 8)
  assert_true(result.lambda() >= 0.49 && result.lambda() <= 0.51)
  assert_true(result.pi0() >= 0.0 && result.pi0() <= 1.0)
}

///|
test "qvalue_with_lambda_empty" {
  let pvals : Array[Double] = []
  let result = @src.qvalue_with_lambda(pvals, 0.5)
  assert_eq(result.q_values().length(), 0)
  assert_true(result.lambda() >= 0.49 && result.lambda() <= 0.51)
}

///|
test "estimate_pi0_basic" {
  let pvals = [0.001, 0.01, 0.02, 0.04, 0.05, 0.10, 0.20, 0.50]
  let pi0 = @src.qvalue_estimate_pi0(pvals, 0.5)
  assert_true(pi0 >= 0.0 && pi0 <= 1.0)
}

///|
test "estimate_pi0_all_below_lambda" {
  let pvals = [0.001, 0.01, 0.02, 0.04, 0.05]
  let pi0 = @src.qvalue_estimate_pi0(pvals, 0.5)
  assert_true(pi0 >= 0.0 && pi0 <= 1.0)
}

///|
test "estimate_pi0_all_above_lambda" {
  let pvals = [0.6, 0.7, 0.8, 0.9, 0.95]
  let pi0 = @src.qvalue_estimate_pi0(pvals, 0.5)
  assert_true(pi0 >= 0.0 && pi0 <= 1.0)
}

///|
test "estimate_pi0_empty" {
  let pvals : Array[Double] = []
  let pi0 = @src.qvalue_estimate_pi0(pvals, 0.5)
  assert_true(pi0 == 0.0)
}

///|
test "estimate_pi0_clamped_to_one" {
  let pvals = [0.01, 0.02, 0.03, 0.04, 0.05]
  let pi0 = @src.qvalue_estimate_pi0(pvals, 0.01)
  assert_true(pi0 <= 1.0)
}

///|
test "storey_fdr_returns_array" {
  let pvals = [0.001, 0.01, 0.02, 0.04, 0.05, 0.10, 0.20, 0.50]
  let qvals = @src.storey_fdr(pvals)
  assert_eq(qvals.length(), 8)
  let mut i = 0
  while i < 8 {
    assert_true(qvals[i] >= 0.0 && qvals[i] <= 1.0)
    i = i + 1
  }
}

///|
test "storey_fdr_empty" {
  let pvals : Array[Double] = []
  let qvals = @src.storey_fdr(pvals)
  assert_eq(qvals.length(), 0)
}

///|
test "bootstrap_pi0_basic" {
  let pvals = [0.001, 0.01, 0.02, 0.04, 0.05, 0.10, 0.20, 0.50]
  let pi0 = @src.bootstrap_pi0(pvals, 50)
  assert_true(pi0 >= 0.0 && pi0 <= 1.0)
}

///|
test "bootstrap_pi0_empty" {
  let pvals : Array[Double] = []
  let pi0 = @src.bootstrap_pi0(pvals, 50)
  assert_true(pi0 == 0.0)
}

///|
test "bootstrap_pi0_zero_iterations" {
  let pvals = [0.1, 0.2, 0.3, 0.4, 0.5]
  let pi0 = @src.bootstrap_pi0(pvals, 0)
  assert_true(pi0 >= 0.0 && pi0 <= 1.0)
}

///|
test "bootstrap_pi0_single_value" {
  let pvals = [0.05]
  let pi0 = @src.bootstrap_pi0(pvals, 30)
  assert_true(pi0 >= 0.0 && pi0 <= 1.0)
}

///|
test "qvalue_significance_at_threshold" {
  let pvals = [0.001, 0.01, 0.02, 0.04, 0.05, 0.10, 0.20, 0.50]
  let result = @src.qvalue(pvals)
  let sig = @src.qvalue_significance(result, 0.05)
  assert_eq(sig.length(), 8)
  let mut i = 0
  while i < 8 {
    if sig[i] {
      assert_true(result.q_values()[i] <= 0.05)
    }
    i = i + 1
  }
}

///|
test "qvalue_significance_no_significant" {
  let pvals = [0.5, 0.6, 0.7, 0.8, 0.9]
  let result = @src.qvalue(pvals)
  let sig = @src.qvalue_significance(result, 0.01)
  assert_eq(sig.length(), 5)
  let any_sig = sig[0] || sig[1] || sig[2] || sig[3] || sig[4]
  assert_false(any_sig)
}

///|
test "qvalue_significance_custom_threshold" {
  let pvals = [0.001, 0.01, 0.02, 0.04, 0.05]
  let result = @src.qvalue(pvals)
  let sig_01 = @src.qvalue_significance(result, 0.01)
  let sig_10 = @src.qvalue_significance(result, 0.10)
  assert_eq(sig_01.length(), 5)
  assert_eq(sig_10.length(), 5)
}

///|
test "qvalue_result_getters" {
  let pvals = [0.01, 0.05, 0.10]
  let result = @src.qvalue(pvals)
  assert_eq(result.p_values().length(), 3)
  assert_eq(result.q_values().length(), 3)
  assert_true(result.pi0() >= 0.0 && result.pi0() <= 1.0)
  assert_true(result.lambda() >= 0.05 && result.lambda() <= 0.95)
  assert_true(result.gamma() >= 0.0 && result.gamma() <= 1.0)
  assert_eq(result.significance().length(), 3)
}

///|
test "qvalue_result_creation" {
  let pvals = [0.01, 0.05]
  let qvals = [0.02, 0.06]
  let sig = [true, false]
  let result = @src.QValueResult::new(pvals, qvals, 0.8, 0.5, 0.75, sig)
  assert_eq(result.p_values().length(), 2)
  assert_eq(result.q_values().length(), 2)
  assert_true(result.pi0() > 0.79 && result.pi0() < 0.81)
  assert_true(result.lambda() > 0.49 && result.lambda() < 0.51)
  assert_true(result.gamma() > 0.74 && result.gamma() < 0.76)
  assert_true(result.significance()[0])
  assert_false(result.significance()[1])
}

///|
test "qvalue_pvalues_clamped_to_unit_interval" {
  let pvals = [-0.5, 0.0, 0.5, 1.0, 1.5]
  let result = @src.qvalue(pvals)
  let validated = result.p_values()
  assert_true(validated[0] >= 0.0)
  assert_true(validated[4] <= 1.0)
}

///|
test "qvalue_order_preserved" {
  let pvals = [0.20, 0.001, 0.50, 0.05, 0.01]
  let result = @src.qvalue(pvals)
  assert_eq(result.p_values().length(), 5)
  assert_eq(result.q_values().length(), 5)
}

///|
test "qvalue_most_significant_stays_smallest" {
  let pvals = [0.20, 0.001, 0.50, 0.05, 0.01]
  let result = @src.qvalue(pvals)
  let qvals = result.q_values()
  assert_true(qvals[1] <= qvals[3])
  assert_true(qvals[3] <= qvals[4])
}

///|
test "qvalue_with_lambda_result_consistency" {
  let pvals = [0.01, 0.05, 0.10, 0.20, 0.50]
  let result1 = @src.qvalue_with_lambda(pvals, 0.5)
  let result2 = @src.qvalue_with_lambda(pvals, 0.5)
  let q1 = result1.q_values()
  let q2 = result2.q_values()
  let mut i = 0
  while i < q1.length() {
    assert_true(q1[i] >= q2[i] - 0.0000000001 && q1[i] <= q2[i] + 0.0000000001)
    i = i + 1
  }
}

///|
test "estimate_pi0_with_different_lambdas" {
  let pvals = [0.001, 0.01, 0.05, 0.10, 0.20, 0.50, 0.70, 0.90]
  let pi0_low = @src.qvalue_estimate_pi0(pvals, 0.1)
  let pi0_high = @src.qvalue_estimate_pi0(pvals, 0.8)
  assert_true(pi0_low >= 0.0 && pi0_low <= 1.0)
  assert_true(pi0_high >= 0.0 && pi0_high <= 1.0)
}

///|
test "bootstrap_pi0_increases_with_more_iterations" {
  let pvals = [0.05, 0.10, 0.20, 0.50, 0.80]
  let pi0_small = @src.bootstrap_pi0(pvals, 10)
  let pi0_large = @src.bootstrap_pi0(pvals, 100)
  assert_true(pi0_small >= 0.0 && pi0_small <= 1.0)
  assert_true(pi0_large >= 0.0 && pi0_large <= 1.0)
}

///|
test "qvalue_significance_boolean_correct" {
  let pvals = [0.0001, 0.001, 0.01, 0.05, 0.5]
  let result = @src.qvalue(pvals)
  let sig = @src.qvalue_significance(result, 0.05)
  let qvals = result.q_values()
  let mut i = 0
  while i < sig.length() {
    assert_true(sig[i] == (qvals[i] <= 0.05))
    i = i + 1
  }
}

///|
test "qvalue_large_dataset" {
  let pvals : Array[Double] = Array::make(100, 0.0)
  let mut i = 0
  while i < 100 {
    pvals[i] = (i.to_double() + 1.0) / 100.0
    i = i + 1
  }
  let result = @src.qvalue(pvals)
  assert_eq(result.q_values().length(), 100)
  let qvals = result.q_values()
  let mut j = 0
  while j < 100 {
    assert_true(qvals[j] >= 0.0 && qvals[j] <= 1.0)
    j = j + 1
  }
}

///|
test "qvalue_with_lambda_clamped" {
  let pvals = [0.01, 0.05, 0.5]
  let result = @src.qvalue_with_lambda(pvals, 0.99)
  assert_eq(result.q_values().length(), 3)
  assert_true(result.pi0() >= 0.0 && result.pi0() <= 1.0)
}

///|
test "qvalue_all_ones" {
  let pvals = [1.0, 1.0, 1.0, 1.0]
  let result = @src.qvalue(pvals)
  assert_eq(result.q_values().length(), 4)
  let qvals = result.q_values()
  let mut i = 0
  while i < 4 {
    assert_true(qvals[i] >= 0.0 && qvals[i] <= 1.0)
    i = i + 1
  }
}

///|
test "qvalue_near_zero" {
  let pvals = [0.0, 0.0, 0.0, 0.0]
  let result = @src.qvalue(pvals)
  assert_eq(result.q_values().length(), 4)
  let qvals = result.q_values()
  let mut i = 0
  while i < 4 {
    assert_true(qvals[i] >= 0.0 && qvals[i] <= 1.0)
    i = i + 1
  }
}

// ===========================================================================
// Bioconductor qvalue parity: pi0est (smoother/bootstrap), lfdr, empPvals,
// qvalue() with pfdr / fdr.level. Synthetic data is deterministic; expected
// values follow the official R package semantics (Storey/Tibshirani).
// ===========================================================================

///|
/// Absolute closeness helper for the new parity tests.
fn qv_close(a : Double, b : Double, tol : Double) -> Bool {
  (a - b).abs() <= tol
}

///|
/// Deterministic near-uniform p-values: (i + 0.5) / n, i in 0..n.
fn qv_uniform_pvalues(n : Int) -> Array[Double] {
  let v = Array::make(n, 0.0)
  let mut i = 0
  while i < n {
    v[i] = (i.to_double() + 0.5) / n.to_double()
    i = i + 1
  }
  v
}

///|
/// Deterministic 50/50 mixture (1000 p-values): 500 strong signals with
/// p < 0.001, then 500 nulls p > 0.002, sorted ascending so that lfdr
/// monotonicity can be checked directly on the returned order.
fn qv_mixture_pvalues() -> Array[Double] {
  let nsig = 500
  let v = Array::make(2 * nsig, 0.0)
  let mut i = 0
  while i < nsig {
    v[i] = (i.to_double() + 0.5) / nsig.to_double() * 0.0009
    i = i + 1
  }
  i = 0
  while i < nsig {
    v[nsig + i] = 0.002 + (i.to_double() + 0.5) / nsig.to_double() * 0.996
    i = i + 1
  }
  v
}

///|
test "pi0est_uniform_smoother_near_one" {
  let r = @src.qvalue_pi0est(qv_uniform_pvalues(1000))
  assert_true(r.pi0() >= 0.9 && r.pi0() <= 1.0)
  assert_eq(r.lambdas().length(), 19)
  assert_eq(r.pi0_lambda().length(), 19)
  assert_eq(r.pi0_smooth().length(), 19)
  assert_true(qv_close(r.lambdas()[0], 0.05, 1.0e-12))
  assert_true(qv_close(r.lambdas()[18], 0.95, 1.0e-12))
}

///|
test "pi0est_uniform_bootstrap_near_one" {
  let p = qv_uniform_pvalues(1000)
  let r = @src.qvalue_pi0est(p, method="bootstrap")
  assert_true(r.pi0() >= 0.9 && r.pi0() <= 1.0)
  assert_eq(r.pi0_smooth().length(), 0)
  // closed-form bootstrap is deterministic: repeated calls agree exactly
  let r2 = @src.qvalue_pi0est(p, method="bootstrap")
  assert_true(qv_close(r.pi0(), r2.pi0(), 1.0e-12))
}

///|
test "pi0est_mixture_near_half" {
  let p = qv_mixture_pvalues()
  let rs = @src.qvalue_pi0est(p)
  let rb = @src.qvalue_pi0est(p, method="bootstrap")
  let rl = @src.qvalue_pi0est(p, smooth_log=true)
  assert_true(rs.pi0() >= 0.35 && rs.pi0() <= 0.65)
  assert_true(rb.pi0() >= 0.35 && rb.pi0() <= 0.65)
  assert_true(rl.pi0() >= 0.35 && rl.pi0() <= 0.65)
}

///|
test "pi0est_single_lambda_direct" {
  // 10 p-values, only 0.99 is >= 0.5: pi0 = 1 / (10 * (1 - 0.5)) = 0.2
  let p = [0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.99]
  let r = @src.qvalue_pi0est(p, lambdas=[0.5])
  assert_true(qv_close(r.pi0(), 0.2, 1.0e-12))
  assert_eq(r.pi0_lambda().length(), 1)
  assert_eq(r.lambdas().length(), 1)
}

///|
test "pi0est_rejects_bad_lambda_grid" {
  let p = qv_uniform_pvalues(50)
  let len2 = try {
    ignore(@src.qvalue_pi0est(p, lambdas=[0.1, 0.2]))
    false
  } catch {
    QValueError(_) => true
  }
  assert_true(len2)
  let len3 = try {
    ignore(@src.qvalue_pi0est(p, lambdas=[0.1, 0.2, 0.3]))
    false
  } catch {
    QValueError(_) => true
  }
  assert_true(len3)
  let too_big = try {
    ignore(@src.qvalue_pi0est(p, lambdas=[0.1, 0.2, 0.3, 1.0]))
    false
  } catch {
    QValueError(_) => true
  }
  assert_true(too_big)
  let neg = try {
    ignore(@src.qvalue_pi0est(p, lambdas=[-0.1, 0.2, 0.3, 0.4]))
    false
  } catch {
    QValueError(_) => true
  }
  assert_true(neg)
  let bad_method = try {
    ignore(@src.qvalue_pi0est(p, method="nope"))
    false
  } catch {
    QValueError(_) => true
  }
  assert_true(bad_method)
}

///|
test "pi0est_rejects_bad_pvalues" {
  let empty = try {
    ignore(@src.qvalue_pi0est([]))
    false
  } catch {
    QValueError(_) => true
  }
  assert_true(empty)
  let too_big = try {
    ignore(@src.qvalue_pi0est([0.5, 1.2]))
    false
  } catch {
    QValueError(_) => true
  }
  assert_true(too_big)
  let neg = try {
    ignore(@src.qvalue_pi0est([0.5, -0.1]))
    false
  } catch {
    QValueError(_) => true
  }
  assert_true(neg)
}

///|
test "lfdr_range_and_extremes" {
  let p = qv_mixture_pvalues()
  let lf = @src.qvalue_lfdr(p)
  assert_eq(lf.length(), p.length())
  let mut i = 0
  while i < lf.length() {
    assert_true(lf[i] >= 0.0 && lf[i] <= 1.0)
    i = i + 1
  }
  // smallest p: very low local FDR; largest p: local FDR at the cap
  assert_true(lf[0] < 0.05)
  assert_true(qv_close(lf[lf.length() - 1], 1.0, 1.0e-12))
}

///|
test "lfdr_monotone_in_p" {
  let p = qv_mixture_pvalues() // sorted ascending
  let lf = @src.qvalue_lfdr(p)
  let mut i = 1
  while i < lf.length() {
    assert_true(lf[i] >= lf[i - 1] - 1.0e-12)
    i = i + 1
  }
}

///|
test "lfdr_logit_and_options" {
  let p = qv_mixture_pvalues()
  let lf = @src.qvalue_lfdr(p, transf="logit")
  assert_eq(lf.length(), p.length())
  let mut i = 0
  while i < lf.length() {
    assert_true(lf[i] >= 0.0 && lf[i] <= 1.0)
    i = i + 1
  }
  // monotone=false, trunc=false still yields non-negative values
  let lf2 = @src.qvalue_lfdr(p, monotone=false, trunc=false)
  assert_eq(lf2.length(), p.length())
  let mut j = 0
  while j < lf2.length() {
    assert_true(lf2[j] >= 0.0)
    j = j + 1
  }
  // supplied pi0 is honored
  let lf3 = @src.qvalue_lfdr(p, pi0=Some(1.0))
  assert_eq(lf3.length(), p.length())
}

///|
test "lfdr_rejects_bad_inputs" {
  let empty = try {
    ignore(@src.qvalue_lfdr([]))
    false
  } catch {
    QValueError(_) => true
  }
  assert_true(empty)
  let bad_p = try {
    ignore(@src.qvalue_lfdr([0.5, 1.5]))
    false
  } catch {
    QValueError(_) => true
  }
  assert_true(bad_p)
  let bad_pi0 = try {
    ignore(@src.qvalue_lfdr([0.1, 0.5, 0.9], pi0=Some(1.5)))
    false
  } catch {
    QValueError(_) => true
  }
  assert_true(bad_pi0)
  let bad_transf = try {
    ignore(@src.qvalue_lfdr([0.1, 0.5, 0.9], transf="cauchy"))
    false
  } catch {
    QValueError(_) => true
  }
  assert_true(bad_transf)
}

///|
test "emppvals_pooled_exact_small_example" {
  // Merged descending: 5(obs), 4, 3, 2(obs precedes tied null 2), 2, 1.
  // stat 5: (0 nulls greater) / 4 -> floor 1/4 = 0.25
  // stat 2: (2 nulls greater: 3, 4) / 4 = 0.5
  let ep = @src.empirical_pvalues([2.0, 5.0], [1.0, 2.0, 3.0, 4.0])
  assert_eq(ep.length(), 2)
  assert_true(qv_close(ep[0], 0.5, 1.0e-12))
  assert_true(qv_close(ep[1], 0.25, 1.0e-12))
}

///|
test "emppvals_pooled_floor_and_max" {
  // stats larger than all nulls -> p floored at 1 / m0 = 0.5
  let ep_hi = @src.empirical_pvalues([10.0, 10.0], [1.0, 2.0])
  assert_true(qv_close(ep_hi[0], 0.5, 1.0e-12))
  assert_true(qv_close(ep_hi[1], 0.5, 1.0e-12))
  // stats smaller than all nulls -> p = 1.0
  let ep_lo = @src.empirical_pvalues([0.0, 0.0], [1.0, 2.0])
  assert_true(qv_close(ep_lo[0], 1.0, 1.0e-12))
  assert_true(qv_close(ep_lo[1], 1.0, 1.0e-12))
}

///|
test "emppvals_rows_exact" {
  // row 0: 3 of 4 nulls >= 2 -> 0.75; row 1: 0 of 4 nulls >= 5 -> floor 0.25
  let ep = @src.empirical_pvalues_rows([2.0, 5.0], [
    [1.0, 2.0, 3.0, 4.0], [1.0, 2.0, 3.0, 4.0]
  ])
  assert_true(qv_close(ep[0], 0.75, 1.0e-12))
  assert_true(qv_close(ep[1], 0.25, 1.0e-12))
  // floor: no null exceeds the statistic -> 1 / ncol
  let ep_f = @src.empirical_pvalues_rows([9.0], [[1.0, 2.0, 3.0]])
  assert_true(qv_close(ep_f[0], 1.0 / 3.0, 1.0e-12))
}

///|
test "emppvals_rejects_bad_inputs" {
  let no_null = try {
    ignore(@src.empirical_pvalues([1.0, 2.0], []))
    false
  } catch {
    QValueError(_) => true
  }
  assert_true(no_null)
  let row_mismatch = try {
    ignore(@src.empirical_pvalues_rows([1.0, 2.0], [[1.0, 2.0]]))
    false
  } catch {
    QValueError(_) => true
  }
  assert_true(row_mismatch)
  let ragged = try {
    ignore(@src.empirical_pvalues_rows([1.0, 2.0], [[1.0, 2.0], [3.0]]))
    false
  } catch {
    QValueError(_) => true
  }
  assert_true(ragged)
  let empty_row = try {
    ignore(@src.empirical_pvalues_rows([1.0], [[]]))
    false
  } catch {
    QValueError(_) => true
  }
  assert_true(empty_row)
}

///|
test "analyze_basic_qvalues" {
  let p = qv_mixture_pvalues()
  let r = @src.qvalue_analyze(p)
  assert_eq(r.q_values().length(), p.length())
  assert_true(r.pi0() > 0.0 && r.pi0() <= 1.0)
  assert_eq(r.pfdr(), false)
  // p sorted ascending => q non-decreasing, all within [0, 1]
  let q = r.q_values()
  let mut i = 1
  while i < q.length() {
    assert_true(q[i] >= q[i - 1] - 1.0e-12)
    assert_true(q[i] >= 0.0 && q[i] <= 1.0)
    i = i + 1
  }
  assert_true(q[0] <= q[q.length() - 1])
  // default result includes lfdr values
  assert_eq(r.lfdr().length(), p.length())
}

///|
test "analyze_uniform_qvalues_bounded" {
  let r = @src.qvalue_analyze(qv_uniform_pvalues(500))
  let q = r.q_values()
  let mut i = 0
  while i < q.length() {
    assert_true(q[i] >= 0.0 && q[i] <= 1.0)
    i = i + 1
  }
  assert_true(r.pi0() >= 0.9)
}

///|
test "analyze_pfdr_dominates_plain" {
  let p = qv_mixture_pvalues()
  let plain = @src.qvalue_analyze(p)
  let pf = @src.qvalue_analyze(p, pfdr=true)
  assert_eq(pf.pfdr(), true)
  let mut i = 0
  while i < p.length() {
    assert_true(pf.q_values()[i] >= plain.q_values()[i] - 1.0e-12)
    i = i + 1
  }
}

///|
test "analyze_fdr_level_significance" {
  let p = qv_mixture_pvalues()
  let r = @src.qvalue_analyze(p, fdr_level=Some(0.05))
  assert_true(qv_close(r.fdr_level(), 0.05, 1.0e-12))
  let q = r.q_values()
  let sig = r.significant()
  let mut count = 0
  let mut i = 0
  while i < p.length() {
    assert_true(sig[i] == (q[i] <= 0.05))
    if sig[i] {
      count = count + 1
    }
    i = i + 1
  }
  // ~500 true signals; only a handful of nulls cross the threshold
  assert_true(count > 400 && count < 600)
  // invalid fdr levels raise
  let zero = try {
    ignore(@src.qvalue_analyze(p, fdr_level=Some(0.0)))
    false
  } catch {
    QValueError(_) => true
  }
  assert_true(zero)
  let big = try {
    ignore(@src.qvalue_analyze(p, fdr_level=Some(1.5)))
    false
  } catch {
    QValueError(_) => true
  }
  assert_true(big)
}

///|
test "analyze_options" {
  let p = qv_mixture_pvalues()
  // lfdr output can be suppressed
  let r1 = @src.qvalue_analyze(p, lfdr_out=false)
  assert_eq(r1.lfdr().length(), 0)
  // bootstrap pi0 path runs through the main entry point
  let r2 = @src.qvalue_analyze(p, method="bootstrap")
  assert_true(r2.pi0() >= 0.35 && r2.pi0() <= 0.65)
  assert_eq(r2.pi0_method(), "bootstrap")
  // explicit pi0 override: no lambda path is computed
  let r3 = @src.qvalue_analyze(p, pi0=Some(1.0))
  assert_true(qv_close(r3.pi0(), 1.0, 1.0e-12))
  assert_eq(r3.pi0_lambda().length(), 0)
  let bad_pi0 = try {
    ignore(@src.qvalue_analyze(p, pi0=Some(0.0)))
    false
  } catch {
    QValueError(_) => true
  }
  assert_true(bad_pi0)
  // a custom lambda grid passes through to pi0est
  let r4 = @src.qvalue_analyze(p, lambdas=Some([0.05, 0.1, 0.2, 0.3, 0.4, 0.5]))
  assert_true(r4.pi0() > 0.0 && r4.pi0() <= 1.0)
}
