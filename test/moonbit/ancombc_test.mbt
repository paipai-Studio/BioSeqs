///|
/// Tests for the ANCOMBC module (Bioconductor ANCOMBC equivalent).
///
/// Covers data containers (AncombcData, AncombcFeature, AncombcSample),
/// sampling fraction estimation, bias-corrected log-abundance,
/// reference feature selection, differential abundance testing,
/// BH-FDR correction, summary, and sample data validation.

// ===========================================================================
// 1. AncombcFeature
// ===========================================================================

///|
test "ancombc_feature_creation" {
  let f = @src.AncombcFeature::new("Bacteroides", [10, 20, 30])
  assert_eq(f.name(), "Bacteroides")
  assert_eq(f.counts().length(), 3)
  assert_eq(f.counts()[0], 10)
  assert_eq(f.counts()[2], 30)
}

///|
test "ancombc_feature_total" {
  let f = @src.AncombcFeature::new("Genus", [10, 20, 30])
  assert_eq(f.total(), 60)
}

///|
test "ancombc_feature_total_empty" {
  let f = @src.AncombcFeature::new("Empty", [])
  assert_eq(f.total(), 0)
}

// ===========================================================================
// 2. AncombcSample
// ===========================================================================

///|
test "ancombc_sample_creation" {
  let s = @src.AncombcSample::new("S1", "control")
  assert_eq(s.name(), "S1")
  assert_eq(s.group(), "control")
}

// ===========================================================================
// 3. AncombcData
// ===========================================================================

///|
test "ancombc_data_dimensions" {
  let data = @src.ancombc_sample_data()
  assert_eq(data.n_samples(), 8)
  assert_eq(data.n_features(), 6)
}

///|
test "ancombc_data_sample_names" {
  let data = @src.ancombc_sample_data()
  let names = data.sample_names()
  assert_eq(names[0], "S1")
  assert_eq(names[7], "S8")
}

///|
test "ancombc_data_sample_groups" {
  let data = @src.ancombc_sample_data()
  let groups = data.sample_groups()
  assert_eq(groups[0], "control")
  assert_eq(groups[4], "treatment")
}

///|
test "ancombc_data_feature_names" {
  let data = @src.ancombc_sample_data()
  let names = data.feature_names()
  assert_eq(names[0], "Bacteroides")
  assert_eq(names[5], "Ruminococcus")
}

///|
test "ancombc_data_get_count" {
  let data = @src.ancombc_sample_data()
  assert_eq(data.get_count(0, 0), 120)
  assert_eq(data.get_count(0, 4), 250)
}

///|
test "ancombc_data_groups" {
  let data = @src.ancombc_sample_data()
  let groups = data.groups()
  assert_eq(groups.length(), 2)
  assert_true(groups.contains("control"))
  assert_true(groups.contains("treatment"))
}

///|
test "ancombc_data_library_sizes" {
  let data = @src.ancombc_sample_data()
  let sizes = data.library_sizes()
  assert_eq(sizes.length(), 8)
  // S1: 120+80+50+30+40+25 = 345
  assert_eq(sizes[0], 345)
}

// ===========================================================================
// 4. Sampling fractions
// ===========================================================================

///|
test "ancombc_sampling_fractions_length" {
  let data = @src.ancombc_sample_data()
  let fracs = @src.ancombc_sampling_fractions(data)
  assert_eq(fracs.length(), 8)
}

///|
test "ancombc_sampling_fractions_sum_to_one" {
  let data = @src.ancombc_sample_data()
  let fracs = @src.ancombc_sampling_fractions(data)
  let mut sum = 0.0
  for f in fracs {
    sum = sum + f
  }
  assert_true(sum > 0.999)
  assert_true(sum < 1.001)
}

///|
test "ancombc_sampling_fractions_positive" {
  let data = @src.ancombc_sample_data()
  let fracs = @src.ancombc_sampling_fractions(data)
  for f in fracs {
    assert_true(f > 0.0)
  }
}

// ===========================================================================
// 5. Bias-corrected log-abundance
// ===========================================================================

///|
test "ancombc_bias_corrected_log_dimensions" {
  let data = @src.ancombc_sample_data()
  let log_data = @src.ancombc_bias_corrected_log(data)
  assert_eq(log_data.length(), 6) // n_features
  assert_eq(log_data[0].length(), 8) // n_samples
}

///|
test "ancombc_bias_corrected_log_finite" {
  let data = @src.ancombc_sample_data()
  let log_data = @src.ancombc_bias_corrected_log(data)
  for row in log_data {
    for v in row {
      assert_true(v == v) // not NaN
      assert_true(v > -100.0)
      assert_true(v < 100.0)
    }
  }
}

// ===========================================================================
// 6. Reference feature selection
// ===========================================================================

///|
test "ancombc_reference_feature_valid_index" {
  let data = @src.ancombc_sample_data()
  let ref_idx = @src.ancombc_reference_feature(data)
  assert_true(ref_idx >= 0)
  assert_true(ref_idx < data.n_features())
}

///|
test "ancombc_reference_feature_is_stable" {
  let data = @src.ancombc_sample_data()
  let ref_idx = @src.ancombc_reference_feature(data)
  // The reference should be one of the non-differential features (2-5)
  // since they have lower variance.
  let log_data = @src.ancombc_bias_corrected_log(data)
  let ref_var = {
    let mut s = 0.0
    let m = {
      let mut sum = 0.0
      for v in log_data[ref_idx] {
        sum = sum + v
      }
      sum / log_data[ref_idx].length().to_double()
    }
    for v in log_data[ref_idx] {
      let d = v - m
      s = s + d * d
    }
    s / (log_data[ref_idx].length() - 1).to_double()
  }
  // The reference should have the lowest variance.
  for i in 0..<data.n_features() {
    let var_i = {
      let mut s = 0.0
      let m = {
        let mut sum = 0.0
        for v in log_data[i] {
          sum = sum + v
        }
        sum / log_data[i].length().to_double()
      }
      for v in log_data[i] {
        let d = v - m
        s = s + d * d
      }
      s / (log_data[i].length() - 1).to_double()
    }
    assert_true(ref_var <= var_i + 0.000000001)
  }
}

// ===========================================================================
// 7. Differential abundance test
// ===========================================================================

///|
test "ancombc_test_basic" {
  let data = @src.ancombc_sample_data()
  let result = @src.ancombc_test(data, "control", "treatment", 0.05)
  assert_eq(result.results().length(), 6)
  assert_eq(result.group1(), "control")
  assert_eq(result.group2(), "treatment")
  assert_true(result.alpha() == 0.05)
}

///|
test "ancombc_test_reference_feature" {
  let data = @src.ancombc_sample_data()
  let result = @src.ancombc_test(data, "control", "treatment", 0.05)
  let ref_name = result.reference_feature()
  assert_true(ref_name.length() > 0)
}

///|
test "ancombc_test_results_have_fields" {
  let data = @src.ancombc_sample_data()
  let result = @src.ancombc_test(data, "control", "treatment", 0.05)
  for r in result.results() {
    assert_true(r.feature().length() > 0)
    // p-values should be in [0, 1]
    assert_true(r.p_value() >= 0.0)
    assert_true(r.p_value() <= 1.0)
    // q-values should be in [0, 1]
    assert_true(r.q_value() >= 0.0)
    assert_true(r.q_value() <= 1.0)
    // W stats should be non-negative
    assert_true(r.w_stat() >= 0)
  }
}

///|
test "ancombc_test_differential_features" {
  let data = @src.ancombc_sample_data()
  let result = @src.ancombc_test(data, "control", "treatment", 0.05)
  // Features 0 and 1 (Bacteroides, Prevotella) have larger count differences.
  let lfc0 = result.results()[0].log_fold_change().abs()
  let lfc1 = result.results()[1].log_fold_change().abs()
  let lfc2 = result.results()[2].log_fold_change().abs()
  let lfc3 = result.results()[3].log_fold_change().abs()
  // Differential features should have larger absolute log fold changes.
  assert_true(lfc0 > lfc2)
  assert_true(lfc0 > lfc3)
  assert_true(lfc1 > lfc2)
  assert_true(lfc1 > lfc3)
}

///|
test "ancombc_test_n_significant" {
  let data = @src.ancombc_sample_data()
  let result = @src.ancombc_test(data, "control", "treatment", 0.05)
  let n_sig = result.n_significant()
  assert_true(n_sig >= 0)
  assert_true(n_sig <= 6)
}

///|
test "ancombc_test_alpha_level" {
  let data = @src.ancombc_sample_data()
  let result = @src.ancombc_test(data, "control", "treatment", 0.1)
  assert_true(result.alpha() == 0.1)
}

// ===========================================================================
// 8. AncombcResult accessors
// ===========================================================================

///|
test "ancombc_result_accessors" {
  let r = @src.AncombcResult::new("test_feature", 1.5, 3, 0.01, 0.05, true)
  assert_eq(r.feature(), "test_feature")
  assert_true(r.log_fold_change() == 1.5)
  assert_eq(r.w_stat(), 3)
  assert_true(r.p_value() == 0.01)
  assert_true(r.q_value() == 0.05)
  assert_true(r.significant())
}

///|
test "ancombc_result_not_significant" {
  let r = @src.AncombcResult::new("feat", 0.1, 0, 0.5, 0.6, false)
  assert_false(r.significant())
}

// ===========================================================================
// 9. Summary
// ===========================================================================

///|
test "ancombc_result_summary_contains_info" {
  let data = @src.ancombc_sample_data()
  let result = @src.ancombc_test(data, "control", "treatment", 0.05)
  let summary = @src.ancombc_result_summary(result)
  assert_true(summary.contains("ANCOMBC"))
  assert_true(summary.contains("control"))
  assert_true(summary.contains("treatment"))
  assert_true(summary.contains("Bacteroides"))
}

// ===========================================================================
// 10. Sample data validation
// ===========================================================================

///|
test "ancombc_sample_data_groups_balanced" {
  let data = @src.ancombc_sample_data()
  let groups = data.sample_groups()
  let mut n_control = 0
  let mut n_treatment = 0
  for g in groups {
    if g == "control" {
      n_control = n_control + 1
    } else if g == "treatment" {
      n_treatment = n_treatment + 1
    }
  }
  assert_eq(n_control, 4)
  assert_eq(n_treatment, 4)
}

///|
test "ancombc_sample_data_counts_positive" {
  let data = @src.ancombc_sample_data()
  for i in 0..<data.n_features() {
    for j in 0..<data.n_samples() {
      assert_true(data.get_count(i, j) > 0)
    }
  }
}

///|
test "ancombc_sample_data_differential_pattern" {
  let data = @src.ancombc_sample_data()
  // Feature 0 (Bacteroides): control ~120, treatment ~250
  let ctrl_mean = (data.get_count(0, 0) + data.get_count(0, 1) + data.get_count(0, 2) + data.get_count(0, 3)).to_double() / 4.0
  let trt_mean = (data.get_count(0, 4) + data.get_count(0, 5) + data.get_count(0, 6) + data.get_count(0, 7)).to_double() / 4.0
  assert_true(trt_mean > ctrl_mean)
}

// ===========================================================================
// 11. Edge cases
// ===========================================================================

///|
test "ancombc_sampling_fractions_empty" {
  let samples : Array[String] = []
  let groups : Array[String] = []
  let features : Array[String] = []
  let counts : Array[Array[Int]] = []
  let data = @src.AncombcData::new(samples, groups, features, counts)
  let fracs = @src.ancombc_sampling_fractions(data)
  assert_eq(fracs.length(), 0)
}

///|
test "ancombc_test_qvalues_sorted" {
  let data = @src.ancombc_sample_data()
  let result = @src.ancombc_test(data, "control", "treatment", 0.05)
  // Collect p-values and q-values.
  let pvals : Array[Double] = []
  let qvals : Array[Double] = []
  for r in result.results() {
    pvals.push(r.p_value())
    qvals.push(r.q_value())
  }
  // q-values should be >= corresponding p-values (monotone after adjustment).
  for i in 0..<pvals.length() {
    assert_true(qvals[i] >= pvals[i] - 0.000000001)
  }
}
