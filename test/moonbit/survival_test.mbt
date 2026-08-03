///|
/// Survival analysis module tests.
/// Covers SurvRecord/SurvData data structures, Kaplan-Meier estimation,
/// log-rank test, Cox proportional hazards model, and statistical helpers.

// ===========================================================================
// SurvRecord creation and accessors
// ===========================================================================

///|
test "survival_record_creation" {
  let rec = @src.SurvRecord::new(5.0, true, "control", [55.0, 0.0])
  assert_true((rec.time() - 5.0).abs() < 0.001)
  assert_true(rec.event())
  assert_eq(rec.group(), "control")
  assert_eq(rec.covariates().length(), 2)
  assert_true((rec.covariates()[0] - 55.0).abs() < 0.001)
  assert_true((rec.covariates()[1] - 0.0).abs() < 0.001)
}

///|
test "survival_record_censored" {
  let rec = @src.SurvRecord::new(10.0, false, "treatment", [1.5])
  assert_true((rec.time() - 10.0).abs() < 0.001)
  assert_true(!rec.event())
  assert_eq(rec.group(), "treatment")
  assert_eq(rec.covariates().length(), 1)
}

///|
test "survival_record_empty_covariates" {
  let rec = @src.SurvRecord::new(3.0, true, "A", [])
  assert_eq(rec.covariates().length(), 0)
  assert_true((rec.time() - 3.0).abs() < 0.001)
}

// ===========================================================================
// SurvData creation (empty, from_records, add_record)
// ===========================================================================

///|
test "survival_data_empty" {
  let data = @src.SurvData::new()
  assert_eq(data.n_samples(), 0)
  assert_eq(data.n_groups(), 0)
  assert_eq(data.group_names().length(), 0)
  assert_eq(data.records().length(), 0)
}

///|
test "survival_data_from_records" {
  let records = [
    @src.SurvRecord::new(5.0, true, "A", [1.0]),
    @src.SurvRecord::new(10.0, false, "B", [2.0]),
    @src.SurvRecord::new(15.0, true, "A", [3.0]),
  ]
  let data = @src.SurvData::from_records(records)
  assert_eq(data.n_samples(), 3)
  assert_eq(data.n_groups(), 2)
  assert_eq(data.records().length(), 3)
  let names = data.group_names()
  assert_eq(names.length(), 2)
  assert_eq(names[0], "A")
  assert_eq(names[1], "B")
}

///|
test "survival_data_add_record" {
  let data = @src.SurvData::new()
  assert_eq(data.n_samples(), 0)
  assert_eq(data.n_groups(), 0)
  data.add_record(@src.SurvRecord::new(5.0, true, "A", []))
  assert_eq(data.n_samples(), 1)
  assert_eq(data.n_groups(), 1)
  data.add_record(@src.SurvRecord::new(10.0, false, "B", []))
  assert_eq(data.n_samples(), 2)
  assert_eq(data.n_groups(), 2)
  data.add_record(@src.SurvRecord::new(15.0, true, "A", []))
  assert_eq(data.n_samples(), 3)
  assert_eq(data.n_groups(), 2)
  let names = data.group_names()
  assert_eq(names[0], "A")
  assert_eq(names[1], "B")
}

// ===========================================================================
// SurvData properties (n_samples, n_groups, group_names, filter_by_group)
// ===========================================================================

///|
test "survival_data_properties" {
  let data = @src.survival_sample_data()
  assert_eq(data.n_samples(), 15)
  assert_eq(data.n_groups(), 2)
  let names = data.group_names()
  assert_eq(names.length(), 2)
  assert_eq(names[0], "control")
  assert_eq(names[1], "treatment")
}

///|
test "survival_data_filter_by_group" {
  let data = @src.survival_sample_data()
  let control = data.filter_by_group("control")
  let treatment = data.filter_by_group("treatment")
  assert_eq(control.length(), 7)
  assert_eq(treatment.length(), 8)
  for rec in control {
    assert_eq(rec.group(), "control")
  }
  for rec in treatment {
    assert_eq(rec.group(), "treatment")
  }
}

///|
test "survival_data_filter_missing_group" {
  let data = @src.survival_sample_data()
  let missing = data.filter_by_group("nonexistent")
  assert_eq(missing.length(), 0)
}

// ===========================================================================
// Kaplan-Meier basic: single group, all events (no censoring)
// ===========================================================================

///|
test "survival_km_basic" {
  let records = [
    @src.SurvRecord::new(1.0, true, "A", []),
    @src.SurvRecord::new(2.0, true, "A", []),
    @src.SurvRecord::new(3.0, true, "A", []),
    @src.SurvRecord::new(4.0, true, "A", []),
    @src.SurvRecord::new(5.0, true, "A", []),
  ]
  let data = @src.SurvData::from_records(records)
  let result = @src.survival_kaplan_meier(data, "A")
  // Survival should decrease to 0 with all events
  let last = result.survival[result.survival.length() - 1]
  assert_true((last - 0.0).abs() < 0.001)
  // Monotonically non-increasing
  let n = result.survival.length()
  for i = 1; i < n; i = i + 1 {
    assert_true(result.survival[i] <= result.survival[i - 1] + 1.0e-10)
  }
}

// ===========================================================================
// Kaplan-Meier with censoring
// ===========================================================================

///|
test "survival_km_with_censoring" {
  let records = [
    @src.SurvRecord::new(1.0, true, "A", []),
    @src.SurvRecord::new(2.0, false, "A", []),
    @src.SurvRecord::new(3.0, true, "A", []),
    @src.SurvRecord::new(4.0, false, "A", []),
  ]
  let data = @src.SurvData::from_records(records)
  let result = @src.survival_kaplan_meier(data, "A")
  // Survival should stay above 0 due to censoring
  let last = result.survival[result.survival.length() - 1]
  assert_true(last > 0.0)
  // At t=1: S = 1*(1-1/4) = 0.75
  assert_true((result.survival[1] - 0.75).abs() < 0.001)
  // At t=3: S = 0.75*(1-1/2) = 0.375
  assert_true((result.survival[2] - 0.375).abs() < 0.001)
}

// ===========================================================================
// Kaplan-Meier result fields (times, survival, n_at_risk, n_events)
// ===========================================================================

///|
test "survival_km_result_fields" {
  let records = [
    @src.SurvRecord::new(1.0, true, "A", []),
    @src.SurvRecord::new(2.0, true, "A", []),
    @src.SurvRecord::new(3.0, true, "A", []),
  ]
  let data = @src.SurvData::from_records(records)
  let result = @src.survival_kaplan_meier(data, "A")
  // t=0 (initial), 1, 2, 3
  assert_eq(result.times.length(), 4)
  assert_eq(result.survival.length(), 4)
  assert_eq(result.std_err.length(), 4)
  assert_eq(result.n_at_risk.length(), 4)
  assert_eq(result.n_events.length(), 4)
  assert_eq(result.group, "A")
  // Initial point
  assert_true((result.times[0] - 0.0).abs() < 0.001)
  assert_eq(result.n_at_risk[0], 3)
  assert_eq(result.n_events[0], 0)
  // t=1
  assert_eq(result.n_at_risk[1], 3)
  assert_eq(result.n_events[1], 1)
  // t=2
  assert_eq(result.n_at_risk[2], 2)
  assert_eq(result.n_events[2], 1)
  // t=3
  assert_eq(result.n_at_risk[3], 1)
  assert_eq(result.n_events[3], 1)
}

// ===========================================================================
// Kaplan-Meier survival starts at 1.0
// ===========================================================================

///|
test "survival_km_starts_at_one" {
  let data = @src.survival_sample_data()
  let result = @src.survival_kaplan_meier(data, "control")
  assert_true((result.survival[0] - 1.0).abs() < 0.001)
  assert_true((result.times[0] - 0.0).abs() < 0.001)
}

// ===========================================================================
// Kaplan-Meier Greenwood standard error
// ===========================================================================

///|
test "survival_km_std_err" {
  let data = @src.survival_sample_data()
  let result = @src.survival_kaplan_meier(data, "control")
  // std_err at t=0 should be 0
  assert_true((result.std_err[0] - 0.0).abs() < 0.001)
  // All std_err values should be non-negative
  for se in result.std_err {
    assert_true(se >= 0.0)
  }
}

// ===========================================================================
// Median survival
// ===========================================================================

///|
test "survival_median_survival" {
  let records = [
    @src.SurvRecord::new(1.0, true, "A", []),
    @src.SurvRecord::new(2.0, true, "A", []),
    @src.SurvRecord::new(3.0, true, "A", []),
    @src.SurvRecord::new(4.0, true, "A", []),
    @src.SurvRecord::new(5.0, true, "A", []),
  ]
  let data = @src.SurvData::from_records(records)
  let result = @src.survival_kaplan_meier(data, "A")
  // S=0.6 at t=2, S=0.4 at t=3 -> median = 3.0
  let med = @src.survival_median_survival(result)
  assert_true(med is Some(_))
  assert_true((med.unwrap() - 3.0).abs() < 0.001)
}

///|
test "survival_median_survival_none" {
  let records = [
    @src.SurvRecord::new(5.0, false, "A", []),
    @src.SurvRecord::new(10.0, false, "A", []),
  ]
  let data = @src.SurvData::from_records(records)
  let result = @src.survival_kaplan_meier(data, "A")
  let med = @src.survival_median_survival(result)
  assert_true(med is None)
}

// ===========================================================================
// Log-rank test: two groups with different survival
// ===========================================================================

///|
test "survival_log_rank_different" {
  let records = [
    @src.SurvRecord::new(5.0, true, "control", []),
    @src.SurvRecord::new(6.0, true, "control", []),
    @src.SurvRecord::new(7.0, true, "control", []),
    @src.SurvRecord::new(8.0, true, "control", []),
    @src.SurvRecord::new(10.0, true, "control", []),
    @src.SurvRecord::new(20.0, true, "treatment", []),
    @src.SurvRecord::new(25.0, true, "treatment", []),
    @src.SurvRecord::new(30.0, true, "treatment", []),
    @src.SurvRecord::new(35.0, true, "treatment", []),
    @src.SurvRecord::new(40.0, true, "treatment", []),
  ]
  let data = @src.SurvData::from_records(records)
  let result = @src.survival_log_rank(data, "control", "treatment")
  assert_true(result.chi_square > 0.0)
  assert_true(result.p_value < 0.05)
}

// ===========================================================================
// Log-rank test: two identical groups
// ===========================================================================

///|
test "survival_log_rank_identical" {
  let records = [
    @src.SurvRecord::new(5.0, true, "A", []),
    @src.SurvRecord::new(6.0, true, "A", []),
    @src.SurvRecord::new(7.0, true, "A", []),
    @src.SurvRecord::new(8.0, true, "A", []),
    @src.SurvRecord::new(10.0, true, "A", []),
    @src.SurvRecord::new(5.0, true, "B", []),
    @src.SurvRecord::new(6.0, true, "B", []),
    @src.SurvRecord::new(7.0, true, "B", []),
    @src.SurvRecord::new(8.0, true, "B", []),
    @src.SurvRecord::new(10.0, true, "B", []),
  ]
  let data = @src.SurvData::from_records(records)
  let result = @src.survival_log_rank(data, "A", "B")
  assert_true(result.chi_square < 0.001)
  assert_true(result.p_value > 0.95)
}

// ===========================================================================
// Log-rank result fields
// ===========================================================================

///|
test "survival_log_rank_fields" {
  let data = @src.survival_sample_data()
  let result = @src.survival_log_rank(data, "control", "treatment")
  assert_eq(result.df, 1)
  assert_eq(result.group1, "control")
  assert_eq(result.group2, "treatment")
  assert_eq(result.n1, 7)
  assert_eq(result.n2, 8)
  assert_true(result.chi_square >= 0.0)
  assert_true(result.p_value >= 0.0 && result.p_value <= 1.0)
}

// ===========================================================================
// Cox model: fit with sample data
// ===========================================================================

///|
test "survival_cox_fit_basic" {
  let data = @src.survival_sample_data()
  let result = @src.survival_cox_fit(data, ["age", "treatment"])
  assert_true(result.converged)
  assert_eq(result.n_samples, 15)
  assert_eq(result.coefficients.length(), 2)
  assert_eq(result.std_errors.length(), 2)
  assert_eq(result.hazard_ratios.length(), 2)
  assert_eq(result.z_values.length(), 2)
  assert_eq(result.p_values.length(), 2)
  assert_eq(result.covariate_names.length(), 2)
  assert_eq(result.covariate_names[0], "age")
  assert_eq(result.covariate_names[1], "treatment")
}

// ===========================================================================
// Cox model: coefficients, hazard ratios, p_values computed
// ===========================================================================

///|
test "survival_cox_coefficients_computed" {
  let data = @src.survival_sample_data()
  let result = @src.survival_cox_fit(data, ["age", "treatment"])
  for coef in result.coefficients {
    assert_true(!coef.is_nan())
  }
  for hr in result.hazard_ratios {
    assert_true(hr > 0.0)
  }
  for p in result.p_values {
    assert_true(p >= 0.0 && p <= 1.0)
  }
  for se in result.std_errors {
    assert_true(se >= 0.0)
  }
}

// ===========================================================================
// Cox model: hazard ratios = exp(coefficients)
// ===========================================================================

///|
test "survival_cox_hazard_ratios" {
  let data = @src.survival_sample_data()
  let result = @src.survival_cox_fit(data, ["age", "treatment"])
  let n = result.coefficients.length()
  for i = 0; i < n; i = i + 1 {
    let expected = @src.survival_exp(result.coefficients[i])
    assert_true((result.hazard_ratios[i] - expected).abs() < 0.001)
  }
}

// ===========================================================================
// Normal CDF
// ===========================================================================

///|
test "survival_normal_cdf" {
  // phi(0) = 0.5
  assert_true((@src.survival_normal_cdf(0.0) - 0.5).abs() < 0.001)
  // phi(1.96) ≈ 0.975
  assert_true((@src.survival_normal_cdf(1.96) - 0.975).abs() < 0.01)
  // phi(-1.96) ≈ 0.025
  assert_true((@src.survival_normal_cdf(-1.96) - 0.025).abs() < 0.01)
  // Monotonicity
  assert_true(@src.survival_normal_cdf(1.0) > @src.survival_normal_cdf(0.0))
  assert_true(@src.survival_normal_cdf(0.0) > @src.survival_normal_cdf(-1.0))
}

// ===========================================================================
// Normal p-value two-sided
// ===========================================================================

///|
test "survival_normal_p_value_two_sided" {
  // z=1.96 -> p ≈ 0.05
  assert_true((@src.survival_normal_p_value_two_sided(1.96) - 0.05).abs() < 0.01)
  // z=0 -> p = 1.0
  assert_true((@src.survival_normal_p_value_two_sided(0.0) - 1.0).abs() < 0.001)
  // p-value should be in [0, 1]
  assert_true(@src.survival_normal_p_value_two_sided(3.0) >= 0.0)
  assert_true(@src.survival_normal_p_value_two_sided(3.0) <= 1.0)
}

// ===========================================================================
// Chi-square p-value
// ===========================================================================

///|
test "survival_chi_square_p_value" {
  // df=1, chi_sq=3.84 ≈ 0.05
  assert_true((@src.survival_chi_square_p_value(3.84, 1) - 0.05).abs() < 0.01)
  // chi_sq=0 -> p=1.0
  assert_true((@src.survival_chi_square_p_value(0.0, 1) - 1.0).abs() < 0.001)
  // p-value in [0, 1] for general df
  assert_true(@src.survival_chi_square_p_value(10.0, 5) >= 0.0)
  assert_true(@src.survival_chi_square_p_value(10.0, 5) <= 1.0)
  // Larger chi_sq gives smaller p-value
  assert_true(@src.survival_chi_square_p_value(10.0, 1) < @src.survival_chi_square_p_value(1.0, 1))
}

// ===========================================================================
// Math helpers
// ===========================================================================

///|
test "survival_exp" {
  assert_true((@src.survival_exp(0.0) - 1.0).abs() < 0.001)
  assert_true((@src.survival_exp(1.0) - 2.71828).abs() < 0.01)
  assert_true((@src.survival_exp(-1.0) - 0.36788).abs() < 0.01)
}

///|
test "survival_ln" {
  assert_true((@src.survival_ln(1.0) - 0.0).abs() < 0.001)
  assert_true((@src.survival_ln(2.71828) - 1.0).abs() < 0.01)
  assert_true((@src.survival_ln(2.0) - 0.69315).abs() < 0.01)
}

///|
test "survival_abs" {
  assert_true((@src.survival_abs(-5.0) - 5.0).abs() < 0.001)
  assert_true((@src.survival_abs(5.0) - 5.0).abs() < 0.001)
  assert_true((@src.survival_abs(0.0) - 0.0).abs() < 0.001)
}

// ===========================================================================
// Summary functions
// ===========================================================================

///|
test "survival_km_summary" {
  let data = @src.survival_sample_data()
  let result = @src.survival_kaplan_meier(data, "control")
  let summary = @src.survival_km_summary(result)
  assert_true(summary.length() > 0)
  assert_true(summary.contains("Kaplan-Meier"))
  assert_true(summary.contains("control"))
}

///|
test "survival_logrank_summary" {
  let data = @src.survival_sample_data()
  let result = @src.survival_log_rank(data, "control", "treatment")
  let summary = @src.survival_logrank_summary(result)
  assert_true(summary.length() > 0)
  assert_true(summary.contains("Log-Rank"))
  assert_true(summary.contains("p-value"))
}

///|
test "survival_cox_summary" {
  let data = @src.survival_sample_data()
  let result = @src.survival_cox_fit(data, ["age", "treatment"])
  let summary = @src.survival_cox_summary(result)
  assert_true(summary.length() > 0)
  assert_true(summary.contains("Cox"))
  assert_true(summary.contains("age"))
  assert_true(summary.contains("treatment"))
}

///|
test "survival_km_to_string" {
  let data = @src.survival_sample_data()
  let result = @src.survival_kaplan_meier(data, "control")
  let s = @src.survival_km_to_string(result)
  assert_true(s.length() > 0)
  assert_true(s.contains("time"))
  assert_true(s.contains("survival"))
  assert_true(s.contains("control"))
}

// ===========================================================================
// Sample data validation
// ===========================================================================

///|
test "survival_sample_data" {
  let data = @src.survival_sample_data()
  assert_eq(data.n_samples(), 15)
  assert_eq(data.n_groups(), 2)
  let names = data.group_names()
  assert_eq(names.length(), 2)
  assert_eq(names[0], "control")
  assert_eq(names[1], "treatment")
  assert_eq(data.records().length(), 15)
  let control = data.filter_by_group("control")
  let treatment = data.filter_by_group("treatment")
  assert_eq(control.length(), 7)
  assert_eq(treatment.length(), 8)
}

// ===========================================================================
// Edge cases: all censored, single sample, empty data
// ===========================================================================

///|
test "survival_km_all_censored" {
  let records = [
    @src.SurvRecord::new(5.0, false, "A", []),
    @src.SurvRecord::new(10.0, false, "A", []),
    @src.SurvRecord::new(15.0, false, "A", []),
  ]
  let data = @src.SurvData::from_records(records)
  let result = @src.survival_kaplan_meier(data, "A")
  // Only the initial point at t=0, no events recorded
  assert_eq(result.times.length(), 1)
  assert_true((result.survival[0] - 1.0).abs() < 0.001)
  assert_eq(result.n_events[0], 0)
  // Median survival should be None
  let med = @src.survival_median_survival(result)
  assert_true(med is None)
}

///|
test "survival_km_single_sample" {
  let records = [@src.SurvRecord::new(5.0, true, "A", [])]
  let data = @src.SurvData::from_records(records)
  let result = @src.survival_kaplan_meier(data, "A")
  // t=0 (initial) and t=5 (event)
  assert_eq(result.times.length(), 2)
  assert_true((result.survival[0] - 1.0).abs() < 0.001)
  assert_true((result.survival[1] - 0.0).abs() < 0.001)
  assert_eq(result.n_at_risk[0], 1)
  assert_eq(result.n_at_risk[1], 1)
  assert_eq(result.n_events[1], 1)
}

///|
test "survival_km_empty_data" {
  let data = @src.SurvData::new()
  let result = @src.survival_kaplan_meier(data, "A")
  // Only the initial point with n_at_risk=0
  assert_eq(result.times.length(), 1)
  assert_true((result.survival[0] - 1.0).abs() < 0.001)
  assert_eq(result.n_at_risk[0], 0)
  assert_eq(result.n_events[0], 0)
}

///|
test "survival_cox_fit_empty" {
  let data = @src.SurvData::new()
  let result = @src.survival_cox_fit(data, ["x"])
  assert_eq(result.coefficients.length(), 0)
  assert_eq(result.n_samples, 0)
  assert_true(result.converged)
}

///|
test "survival_cox_fit_no_covariates" {
  let data = @src.survival_sample_data()
  let result = @src.survival_cox_fit(data, [])
  assert_eq(result.coefficients.length(), 0)
  assert_true(result.converged)
}
