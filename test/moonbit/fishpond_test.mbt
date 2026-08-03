///|
/// Tests for fishpond/Swish non-parametric differential expression.

// ---------------------------------------------------------------------------
// Data structure creation and accessors
// ---------------------------------------------------------------------------

test "fish_counts_creation" {
  let c = @src.fish_sample_counts()
  assert_eq(c.transcript_names().length(), 10)
  assert_eq(c.sample_names().length(), 16)
  assert_eq(c.condition().length(), 16)
  assert_eq(c.counts().length(), 10)
  assert_eq(c.counts()[0].length(), 16)
  assert_eq(c.condition()[0], "case")
  assert_eq(c.condition()[8], "control")
}

test "fish_result_creation" {
  let r = @src.FishResult::new(
    "TX1", 15.0, 1.5, 0.01, 0.05, "up",
  )
  assert_eq(r.transcript(), "TX1")
  assert_eq(r.statistic(), 15.0)
  assert_eq(r.log2_fold_change(), 1.5)
  assert_eq(r.p_value(), 0.01)
  assert_eq(r.fdr(), 0.05)
  assert_eq(r.direction(), "up")
}

// ---------------------------------------------------------------------------
// Mann-Whitney-Wilcoxon statistic
// ---------------------------------------------------------------------------

test "fish_mw_statistic_clear_separation" {
  // Cases all higher than controls
  let stat = @src.fish_mann_whitney_test(
    [10.0, 20.0, 30.0],
    [1.0, 2.0, 3.0],
  )
  // With complete separation, W = sum of ranks 4,5,6 = 15
  assert_eq(stat, 15.0)
}

test "fish_mw_statistic_identical_groups" {
  let stat = @src.fish_mann_whitney_test(
    [5.0, 5.0, 5.0],
    [5.0, 5.0, 5.0],
  )
  // All tied: average rank = 3.5 for each, sum for case = 3*3.5 = 10.5
  assert_eq(stat, 10.5)
}

test "fish_mw_statistic_overlap" {
  let stat = @src.fish_mann_whitney_test(
    [1.0, 3.0, 5.0],
    [2.0, 4.0, 6.0],
  )
  // Ranks: 1->1(case), 2->2(ctrl), 3->3(case), 4->4(ctrl), 5->5(case), 6->6(ctrl)
  // W_case = 1+3+5 = 9
  assert_eq(stat, 9.0)
}

test "fish_mw_statistic_empty_group" {
  let stat = @src.fish_mann_whitney_test([], [1.0, 2.0])
  assert_eq(stat, 0.0)
}

// ---------------------------------------------------------------------------
// Log2 fold change
// ---------------------------------------------------------------------------

test "fish_log2fc_positive" {
  let lfc = @src.fish_log2fc_test([8.0, 8.0], [2.0, 2.0])
  // log2(9) - log2(3) = log2(3) ≈ 1.585
  assert_true(lfc > 0.0)
}

test "fish_log2fc_negative" {
  let lfc = @src.fish_log2fc_test([2.0, 2.0], [8.0, 8.0])
  assert_true(lfc < 0.0)
}

test "fish_log2fc_zero_means" {
  let lfc = @src.fish_log2fc_test([0.0, 0.0], [0.0, 0.0])
  assert_eq(lfc, 0.0)
}

// ---------------------------------------------------------------------------
// Full Swish analysis
// ---------------------------------------------------------------------------

test "fish_swish_detects_upregulated" {
  let counts = @src.fish_sample_counts()
  let results = @src.fish_swish(counts, "case", 500)
  // TX1-3 should be up-regulated
  assert_eq(results[0].transcript(), "TX1")
  assert_eq(results[0].direction(), "up")
  assert_eq(results[1].transcript(), "TX2")
  assert_eq(results[1].direction(), "up")
  assert_eq(results[2].transcript(), "TX3")
  assert_eq(results[2].direction(), "up")
}

test "fish_swish_detects_downregulated" {
  let counts = @src.fish_sample_counts()
  let results = @src.fish_swish(counts, "case", 500)
  // TX4-6 should be down-regulated
  assert_eq(results[3].transcript(), "TX4")
  assert_eq(results[3].direction(), "down")
  assert_eq(results[4].transcript(), "TX5")
  assert_eq(results[4].direction(), "down")
  assert_eq(results[5].transcript(), "TX6")
  assert_eq(results[5].direction(), "down")
}

test "fish_swish_nonsignificant_stable" {
  let counts = @src.fish_sample_counts()
  let results = @src.fish_swish(counts, "case", 500)
  // TX7-10 should not be significant
  for i in 6..<10 {
    assert_eq(results[i].direction(), "ns")
  }
}

test "fish_swish_pvalues_in_range" {
  let counts = @src.fish_sample_counts()
  let results = @src.fish_swish(counts, "case", 500)
  for r in results {
    assert_true(r.p_value() >= 0.0)
    assert_true(r.p_value() <= 1.0)
  }
}

test "fish_swish_fdr_in_range" {
  let counts = @src.fish_sample_counts()
  let results = @src.fish_swish(counts, "case", 500)
  for r in results {
    assert_true(r.fdr() >= 0.0)
    assert_true(r.fdr() <= 1.0)
  }
}

test "fish_swish_empty_input" {
  let counts = @src.FishCounts::new([], [], [], [])
  let results = @src.fish_swish(counts, "case", 10)
  assert_eq(results.length(), 0)
}

test "fish_swish_too_few_samples" {
  let counts = @src.FishCounts::new(
    ["TX1"],
    ["s1"],
    ["case"],
    [[5.0]],
  )
  let results = @src.fish_swish(counts, "case", 10)
  assert_eq(results.length(), 0)
}

test "fish_swish_log2fc_sign_correct" {
  let counts = @src.fish_sample_counts()
  let results = @src.fish_swish(counts, "case", 500)
  // Up-regulated transcripts should have positive log2FC
  for i in 0..<3 {
    assert_true(results[i].log2_fold_change() > 0.0)
  }
  // Down-regulated transcripts should have negative log2FC
  for i in 3..<6 {
    assert_true(results[i].log2_fold_change() < 0.0)
  }
}

// ---------------------------------------------------------------------------
// Significant filtering
// ---------------------------------------------------------------------------

test "fish_significant_filter" {
  let counts = @src.fish_sample_counts()
  let results = @src.fish_swish(counts, "case", 500)
  let sig = @src.fish_significant(results, 0.05)
  // At least the up/down regulated transcripts should be significant
  assert_true(sig.length() >= 6)
}

test "fish_significant_strict_threshold" {
  let counts = @src.fish_sample_counts()
  let results = @src.fish_swish(counts, "case", 500)
  let sig = @src.fish_significant(results, 0.0)
  // No results should pass FDR < 0.0
  assert_eq(sig.length(), 0)
}

// ---------------------------------------------------------------------------
// Output formatting
// ---------------------------------------------------------------------------

test "fish_result_to_string" {
  let r = @src.FishResult::new(
    "TX1", 15.0, 1.5, 0.01, 0.05, "up",
  )
  let s = r.to_string()
  assert_true(s.contains("TX1"))
  assert_true(s.contains("up"))
}

test "fish_results_to_string" {
  let counts = @src.fish_sample_counts()
  let results = @src.fish_swish(counts, "case", 500)
  let s = @src.fish_results_to_string(results)
  assert_true(s.contains("Swish Results"))
  assert_true(s.contains("10 transcripts"))
}

// ---------------------------------------------------------------------------
// Sample data verification
// ---------------------------------------------------------------------------

test "fish_sample_data_up_pattern" {
  let c = @src.fish_sample_counts()
  // TX1 (transcript 0): case values (indices 0-7) should be higher than control (indices 8-15)
  let case_mean = (c.counts()[0][0] + c.counts()[0][1] + c.counts()[0][2] + c.counts()[0][3] + c.counts()[0][4] + c.counts()[0][5] + c.counts()[0][6] + c.counts()[0][7]) / 8.0
  let ctrl_mean = (c.counts()[0][8] + c.counts()[0][9] + c.counts()[0][10] + c.counts()[0][11] + c.counts()[0][12] + c.counts()[0][13] + c.counts()[0][14] + c.counts()[0][15]) / 8.0
  assert_true(case_mean > ctrl_mean)
}

test "fish_sample_data_down_pattern" {
  let c = @src.fish_sample_counts()
  // TX4 (transcript 3): case values should be lower than control
  let case_mean = (c.counts()[3][0] + c.counts()[3][1] + c.counts()[3][2] + c.counts()[3][3] + c.counts()[3][4] + c.counts()[3][5] + c.counts()[3][6] + c.counts()[3][7]) / 8.0
  let ctrl_mean = (c.counts()[3][8] + c.counts()[3][9] + c.counts()[3][10] + c.counts()[3][11] + c.counts()[3][12] + c.counts()[3][13] + c.counts()[3][14] + c.counts()[3][15]) / 8.0
  assert_true(case_mean < ctrl_mean)
}

test "fish_sample_data_ns_pattern" {
  let c = @src.fish_sample_counts()
  // TX7 (transcript 6): similar values in both groups
  let case_mean = (c.counts()[6][0] + c.counts()[6][1] + c.counts()[6][2] + c.counts()[6][3] + c.counts()[6][4] + c.counts()[6][5] + c.counts()[6][6] + c.counts()[6][7]) / 8.0
  let ctrl_mean = (c.counts()[6][8] + c.counts()[6][9] + c.counts()[6][10] + c.counts()[6][11] + c.counts()[6][12] + c.counts()[6][13] + c.counts()[6][14] + c.counts()[6][15]) / 8.0
  // Difference should be small
  let diff = if case_mean > ctrl_mean { case_mean - ctrl_mean } else { ctrl_mean - case_mean }
  assert_true(diff < 10.0)
}
