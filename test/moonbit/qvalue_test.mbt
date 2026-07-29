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
  let result = @src.QValueResult::new(
    pvals,
    qvals,
    0.8,
    0.5,
    0.75,
    sig,
  )
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