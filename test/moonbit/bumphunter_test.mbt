///|
/// Tests for bumphunter methylation/signal region detection.

// ---------------------------------------------------------------------------
// Data structure creation and accessors
// ---------------------------------------------------------------------------

test "bump_position_creation" {
  let p = @src.BumpPosition::new("chr1", 1000, [0.5, 0.6, 0.4, 0.3])
  assert_eq(p.chrom(), "chr1")
  assert_eq(p.pos(), 1000)
  assert_eq(p.values().length(), 4)
  assert_eq(p.values()[0], 0.5)
}

test "bump_region_creation" {
  let r = @src.BumpRegion::new("chr1", 1000, 2000, 3.5, 10.0, 0, 5)
  assert_eq(r.chrom(), "chr1")
  assert_eq(r.start(), 1000)
  assert_eq(r.end_(), 2000)
  assert_eq(r.value(), 3.5)
  assert_eq(r.area(), 10.0)
  assert_eq(r.index_start(), 0)
  assert_eq(r.index_end(), 5)
  assert_eq(r.length(), 1000)
}

test "bump_result_creation" {
  let r = @src.BumpResult::new(
    "chr1", 1000, 2000, 3.5, 10.0, 0.01, 0.05, 0, 5,
  )
  assert_eq(r.chrom(), "chr1")
  assert_eq(r.start(), 1000)
  assert_eq(r.end_(), 2000)
  assert_eq(r.value(), 3.5)
  assert_eq(r.area(), 10.0)
  assert_eq(r.p_value(), 0.01)
  assert_eq(r.fdr(), 0.05)
  assert_eq(r.index_start(), 0)
  assert_eq(r.index_end(), 5)
}

// ---------------------------------------------------------------------------
// t-statistic computation
// ---------------------------------------------------------------------------

test "bump_t_statistics_clear_difference" {
  // Position with clear difference between groups
  let positions = [
    @src.BumpPosition::new("chr1", 100, [1.0, 1.0, 5.0, 5.0]),
  ]
  let stats = @src.bump_compute_t_statistics_test(positions, 2)
  // Group1 mean=1, Group2 mean=5 => negative t-stat
  assert_true(stats[0] < 0.0)
}

test "bump_t_statistics_no_difference" {
  let positions = [
    @src.BumpPosition::new("chr1", 100, [3.0, 3.0, 3.0, 3.0]),
  ]
  let stats = @src.bump_compute_t_statistics_test(positions, 2)
  assert_eq(stats[0], 0.0)
}

// ---------------------------------------------------------------------------
// Smoothing
// ---------------------------------------------------------------------------

test "bump_smooth_basic" {
  let stats = [1.0, 2.0, 3.0, 4.0, 5.0]
  let smoothed = @src.bump_smooth_test(stats, 1)
  // With k=1, position 0 = mean(1,2) = 1.5
  assert_eq(smoothed[0], 1.5)
  // position 2 = mean(2,3,4) = 3.0
  assert_eq(smoothed[2], 3.0)
  // position 4 = mean(4,5) = 4.5
  assert_eq(smoothed[4], 4.5)
}

test "bump_smooth_preserves_constant" {
  let stats = [5.0, 5.0, 5.0, 5.0]
  let smoothed = @src.bump_smooth_test(stats, 1)
  for v in smoothed {
    assert_eq(v, 5.0)
  }
}

// ---------------------------------------------------------------------------
// Candidate bump finding
// ---------------------------------------------------------------------------

test "bump_find_candidates_basic" {
  let positions = [
    @src.BumpPosition::new("chr1", 100, [0.0]),
    @src.BumpPosition::new("chr1", 200, [3.0]),
    @src.BumpPosition::new("chr1", 300, [4.0]),
    @src.BumpPosition::new("chr1", 400, [0.0]),
  ]
  let smoothed = [0.0, 3.0, 4.0, 0.0]
  let bumps = @src.bump_find_candidates_test(positions, smoothed, 2.0)
  assert_eq(bumps.length(), 1)
  assert_eq(bumps[0].start(), 200)
  assert_eq(bumps[0].end_(), 300)
  assert_eq(bumps[0].value(), 4.0)
}

test "bump_find_candidates_none" {
  let positions = [
    @src.BumpPosition::new("chr1", 100, [0.0]),
    @src.BumpPosition::new("chr1", 200, [1.0]),
  ]
  let smoothed = [0.0, 1.0]
  let bumps = @src.bump_find_candidates_test(positions, smoothed, 2.0)
  assert_eq(bumps.length(), 0)
}

test "bump_find_candidates_negative_values" {
  let positions = [
    @src.BumpPosition::new("chr1", 100, [0.0]),
    @src.BumpPosition::new("chr1", 200, [-3.0]),
    @src.BumpPosition::new("chr1", 300, [-4.0]),
    @src.BumpPosition::new("chr1", 400, [0.0]),
  ]
  let smoothed = [0.0, -3.0, -4.0, 0.0]
  let bumps = @src.bump_find_candidates_test(positions, smoothed, 2.0)
  assert_eq(bumps.length(), 1)
  assert_true(bumps[0].value() < 0.0)
}

test "bump_find_candidates_chromosome_boundary" {
  let positions = [
    @src.BumpPosition::new("chr1", 100, [3.0]),
    @src.BumpPosition::new("chr2", 200, [4.0]),
  ]
  let smoothed = [3.0, 4.0]
  let bumps = @src.bump_find_candidates_test(positions, smoothed, 2.0)
  // Should be 2 separate bumps since they're on different chromosomes
  assert_eq(bumps.length(), 2)
  assert_eq(bumps[0].chrom(), "chr1")
  assert_eq(bumps[1].chrom(), "chr2")
}

// ---------------------------------------------------------------------------
// Full analysis pipeline
// ---------------------------------------------------------------------------

test "bump_hunt_detects_known_bump" {
  let data = @src.bump_sample_data()
  let results = @src.bump_hunt(data, 4, 2.0, 2, 20)
  // Should detect at least one bump in the known region
  assert_true(results.length() >= 1)
  // The bump should be on chr1
  assert_eq(results[0].chrom(), "chr1")
}

test "bump_hunt_empty_input" {
  let results = @src.bump_hunt([], 4, 2.0, 2, 10)
  assert_eq(results.length(), 0)
}

test "bump_hunt_zero_group_size" {
  let data = @src.bump_sample_data()
  let results = @src.bump_hunt(data, 0, 2.0, 2, 10)
  assert_eq(results.length(), 0)
}

test "bump_hunt_high_cutoff_finds_nothing" {
  let data = @src.bump_sample_data()
  // Very high cutoff should find no bumps
  let results = @src.bump_hunt(data, 4, 100.0, 2, 10)
  assert_eq(results.length(), 0)
}

test "bump_hunt_result_has_valid_pvalue" {
  let data = @src.bump_sample_data()
  let results = @src.bump_hunt(data, 4, 2.0, 2, 20)
  for r in results {
    assert_true(r.p_value() >= 0.0)
    assert_true(r.p_value() <= 1.0)
  }
}

test "bump_hunt_result_has_valid_fdr" {
  let data = @src.bump_sample_data()
  let results = @src.bump_hunt(data, 4, 2.0, 2, 20)
  for r in results {
    assert_true(r.fdr() >= 0.0)
    assert_true(r.fdr() <= 1.0)
  }
}

// ---------------------------------------------------------------------------
// Sample data
// ---------------------------------------------------------------------------

test "bump_sample_data_structure" {
  let data = @src.bump_sample_data()
  assert_eq(data.length(), 30)
  assert_eq(data[0].chrom(), "chr1")
  assert_eq(data[0].values().length(), 8)
  // Check positions are sequential
  assert_eq(data[0].pos(), 1000)
  assert_eq(data[1].pos(), 1100)
  assert_eq(data[29].pos(), 3900)
}

test "bump_sample_data_has_bump_signal" {
  let data = @src.bump_sample_data()
  // Cases (first 4 samples) should have higher values in positions 10-20
  let case_mean_at_5 = (data[5].values()[0] + data[5].values()[1] + data[5].values()[2] + data[5].values()[3]) / 4.0
  let case_mean_at_15 = (data[15].values()[0] + data[15].values()[1] + data[15].values()[2] + data[15].values()[3]) / 4.0
  // Position 15 (in bump region) should have higher case mean than position 5 (outside)
  assert_true(case_mean_at_15 > case_mean_at_5)
}

// ---------------------------------------------------------------------------
// Output formatting
// ---------------------------------------------------------------------------

test "bump_result_to_string" {
  let r = @src.BumpResult::new(
    "chr1", 1000, 2000, 3.5, 10.0, 0.01, 0.05, 0, 5,
  )
  let s = r.to_string()
  assert_true(s.contains("chr1"))
  assert_true(s.contains("1000"))
  assert_true(s.contains("2000"))
}

test "bump_results_to_string" {
  let results = [
    @src.BumpResult::new("chr1", 1000, 2000, 3.5, 10.0, 0.01, 0.05, 0, 5),
  ]
  let s = @src.bump_results_to_string(results)
  assert_true(s.contains("Bumphunter Results"))
  assert_true(s.contains("1 bumps"))
}

// ---------------------------------------------------------------------------
// Area computation
// ---------------------------------------------------------------------------

test "bump_area_is_sum_of_smoothed" {
  let positions = [
    @src.BumpPosition::new("chr1", 100, [0.0]),
    @src.BumpPosition::new("chr1", 200, [3.0]),
    @src.BumpPosition::new("chr1", 300, [5.0]),
    @src.BumpPosition::new("chr1", 400, [0.0]),
  ]
  let smoothed = [0.0, 3.0, 5.0, 0.0]
  let bumps = @src.bump_find_candidates_test(positions, smoothed, 2.0)
  assert_eq(bumps.length(), 1)
  // Area = 3.0 + 5.0 = 8.0
  assert_eq(bumps[0].area(), 8.0)
}

// ---------------------------------------------------------------------------
// Bump length
// ---------------------------------------------------------------------------

test "bump_length_calculation" {
  let r = @src.BumpRegion::new("chr1", 1000, 2500, 3.0, 9.0, 0, 3)
  assert_eq(r.length(), 1500)
}
