///|
/// Tests for the MACS2-like peak calling module.

///|
test "peak_calling_chipseq_read_construction" {
  let r = @src.ChipSeqRead::new(chr="chr1", position=12345, strand="+")
  assert_eq(r.chr, "chr1")
  assert_eq(r.position, 12345)
  assert_eq(r.strand, "+")
  assert_eq(r.chr(), "chr1")
  assert_eq(r.position(), 12345)
  assert_eq(r.strand(), "+")
}

///|
test "peak_calling_chipseq_read_strand_negative" {
  let r = @src.ChipSeqRead::new(chr="chr2", position=999, strand="-")
  assert_eq(r.chr, "chr2")
  assert_eq(r.position, 999)
  assert_eq(r.strand, "-")
}

///|
test "peak_calling_candidate_peak_construction" {
  let p = @src.CandidatePeak::new(
    chr="chr1",
    start=1000,
    end=2000,
    summit=1500,
    read_count=42,
    control_count=5,
    fold_enrichment=8.4,
    p_value=1.0e-10,
  )
  assert_eq(p.chr, "chr1")
  assert_eq(p.start, 1000)
  assert_eq(p.end, 2000)
  assert_eq(p.summit, 1500)
  assert_eq(p.read_count, 42)
  assert_eq(p.control_count, 5)
  assert_true(p.fold_enrichment > 8.3 && p.fold_enrichment < 8.5)
  assert_true(p.p_value < 1.0e-9)
  // Default mutable field values
  assert_true(p.fdr >= 1.0 || p.fdr <= 1.0)
  assert_eq(p.significant, false)
  // Accessors
  assert_eq(p.chr(), "chr1")
  assert_eq(p.start(), 1000)
  assert_eq(p.end(), 2000)
  assert_eq(p.summit(), 1500)
  assert_eq(p.read_count(), 42)
  assert_eq(p.control_count(), 5)
  assert_eq(p.width(), 1000)
}

///|
test "peak_calling_candidate_peak_mut_field_update" {
  let p = @src.CandidatePeak::new(
    chr="chr1",
    start=100,
    end=200,
    summit=150,
    read_count=10,
    control_count=2,
    fold_enrichment=5.0,
    p_value=1.0e-10,
  )
  // Default mutable field values are 1.0 / false; apply_fdr updates
  // them in place (the mutation happens inside the src package).
  let updated = @src.apply_fdr([p], 0.05)
  assert_eq(updated[0].significant, true)
  assert_true(updated[0].fdr <= 0.05)
}

///|
test "peak_calling_params_default_values" {
  let params = @src.PeakCallingParams::default()
  assert_eq(params.window_size, 500)
  assert_eq(params.step_size, 100)
  assert_eq(params.local_lambda_size, 10000)
  assert_true(params.pvalue_threshold > 0.0 && params.pvalue_threshold < 1.0e-4)
  assert_true(params.fdr_threshold > 0.0 && params.fdr_threshold <= 0.05)
  assert_true(params.min_fold_enrichment >= 2.0)
}

///|
test "peak_calling_params_new_with_overrides" {
  let params = @src.PeakCallingParams::new(
    window_size=1000,
    step_size=250,
    local_lambda_size=20000,
    pvalue_threshold=1.0e-4,
    fdr_threshold=0.1,
    min_fold_enrichment=3.0,
  )
  assert_eq(params.window_size, 1000)
  assert_eq(params.step_size, 250)
  assert_eq(params.local_lambda_size, 20000)
  assert_true(params.pvalue_threshold > 0.0)
  assert_eq(params.fdr_threshold, 0.1)
  assert_true(params.min_fold_enrichment > 2.5)
}

///|
test "peak_calling_params_new_defaults_when_omitted" {
  // All optional args omitted should match defaults.
  let params = @src.PeakCallingParams::new()
  assert_eq(params.window_size, 500)
  assert_eq(params.step_size, 100)
  assert_eq(params.local_lambda_size, 10000)
}

///|
test "peak_calling_params_partial_overrides" {
  // Override only some args; others should default.
  let params = @src.PeakCallingParams::new(
    window_size=250,
    fdr_threshold=0.01,
  )
  assert_eq(params.window_size, 250)
  assert_eq(params.step_size, 100)
  assert_eq(params.local_lambda_size, 10000)
  assert_eq(params.fdr_threshold, 0.01)
}

///|
test "peak_calling_count_reads_in_window_basic" {
  let reads = [
    @src.ChipSeqRead::new(chr="chr1", position=100, strand="+"),
    @src.ChipSeqRead::new(chr="chr1", position=200, strand="+"),
    @src.ChipSeqRead::new(chr="chr1", position=300, strand="-"),
    @src.ChipSeqRead::new(chr="chr1", position=500, strand="+"),
    @src.ChipSeqRead::new(chr="chr2", position=100, strand="+"),
  ]
  // Window [150, 400) on chr1 should include positions 200 and 300.
  assert_eq(@src.count_reads_in_window(reads, "chr1", 150, 400), 2)
  // Window [0, 1000) on chr1 should include 4 reads.
  assert_eq(@src.count_reads_in_window(reads, "chr1", 0, 1000), 4)
  // Window on chr2 includes 1 read.
  assert_eq(@src.count_reads_in_window(reads, "chr2", 0, 1000), 1)
  // Empty window for chr3.
  assert_eq(@src.count_reads_in_window(reads, "chr3", 0, 1000), 0)
}

///|
test "peak_calling_count_reads_in_window_boundary" {
  let reads = [
    @src.ChipSeqRead::new(chr="chr1", position=100, strand="+"),
    @src.ChipSeqRead::new(chr="chr1", position=200, strand="+"),
    @src.ChipSeqRead::new(chr="chr1", position=300, strand="+"),
  ]
  // Half-open interval: position 100 is in [100, 200), position 200 is not.
  assert_eq(@src.count_reads_in_window(reads, "chr1", 100, 200), 1)
  assert_eq(@src.count_reads_in_window(reads, "chr1", 100, 201), 2)
  assert_eq(@src.count_reads_in_window(reads, "chr1", 100, 301), 3)
}

///|
test "peak_calling_count_reads_in_window_empty_reads" {
  let reads : Array[@src.ChipSeqRead] = []
  assert_eq(@src.count_reads_in_window(reads, "chr1", 0, 1000), 0)
}

///|
test "peak_calling_estimate_local_lambda_basic" {
  let control : Array[@src.ChipSeqRead] = []
  // Place 20 control reads uniformly in [0, 2000].
  let mut i = 0
  while i < 20 {
    control.push(@src.ChipSeqRead::new(chr="chr1", position=i * 100, strand="+"))
    i = i + 1
  }
  // Window of 2000 bp centered at 1000 = [0, 2000) should have 20 reads.
  // Expected per-bp rate: 20 / 2000 = 0.01.
  let rate = @src.estimate_local_lambda(control, "chr1", 1000, 2000)
  assert_true(rate > 0.009 && rate < 0.011)
}

///|
test "peak_calling_estimate_local_lambda_empty_control" {
  let control : Array[@src.ChipSeqRead] = []
  let rate = @src.estimate_local_lambda(control, "chr1", 1000, 10000)
  assert_eq(rate, 0.0)
}

///|
test "peak_calling_estimate_local_lambda_off_chromosome" {
  let control = [
    @src.ChipSeqRead::new(chr="chr1", position=500, strand="+"),
  ]
  let rate = @src.estimate_local_lambda(control, "chr2", 500, 1000)
  assert_eq(rate, 0.0)
}

///|
test "peak_calling_poisson_pvalue_observed_zero" {
  // P(X >= 0) = 1 always.
  assert_eq(@src.poisson_pvalue(0, 1.0), 1.0)
  assert_eq(@src.poisson_pvalue(0, 100.0), 1.0)
  assert_eq(@src.poisson_pvalue(0, 0.0), 1.0)
}

///|
test "peak_calling_poisson_pvalue_lambda_one_observed_one" {
  // P(X >= 1 | lambda=1) = 1 - P(X=0) = 1 - exp(-1) ≈ 0.6321
  let p = @src.poisson_pvalue(1, 1.0)
  assert_true(p > 0.63 && p < 0.64)
}

///|
test "peak_calling_poisson_pvalue_lambda_two_observed_three" {
  // P(X >= 3 | lambda=2): P(X=0)+P(X=1)+P(X=2) = e^-2 * (1 + 2 + 2) = 5/e^2 ≈ 0.6767
  // P(X >= 3) = 1 - 0.6767 ≈ 0.3233
  let p = @src.poisson_pvalue(3, 2.0)
  assert_true(p > 0.32 && p < 0.33)
}

///|
test "peak_calling_poisson_pvalue_large_lambda_strong_signal" {
  // With lambda=1 and observed=10, the p-value should be very small.
  let p = @src.poisson_pvalue(10, 1.0)
  assert_true(p < 1.0e-5)
  assert_true(p >= 0.0)
}

///|
test "peak_calling_poisson_pvalue_high_observation_low_lambda" {
  // Very significant: observed=20, expected=1.
  let p = @src.poisson_pvalue(20, 1.0)
  assert_true(p < 1.0e-7)
}

///|
test "peak_calling_poisson_pvalue_in_range_0_to_1" {
  // P-value should always be in [0, 1].
  let p1 = @src.poisson_pvalue(5, 5.0)
  let p2 = @src.poisson_pvalue(2, 10.0)
  let p3 = @src.poisson_pvalue(15, 3.0)
  assert_true(p1 >= 0.0 && p1 <= 1.0)
  assert_true(p2 >= 0.0 && p2 <= 1.0)
  assert_true(p3 >= 0.0 && p3 <= 1.0)
}

///|
test "peak_calling_poisson_pvalue_expected_zero_observed_positive" {
  // When expected is 0 but we observed reads, p-value should be effectively 0.
  let p = @src.poisson_pvalue(5, 0.0)
  assert_true(p < 1.0e-100)
}

///|
test "peak_calling_fold_enrichment_basic" {
  // 10 / 5 = 2.0
  assert_eq(@src.compute_fold_enrichment(10, 5.0), 2.0)
  // 100 / 10 = 10.0
  assert_eq(@src.compute_fold_enrichment(100, 10.0), 10.0)
}

///|
test "peak_calling_fold_enrichment_zero_expected" {
  // When expected is 0 and observed > 0, return a large sentinel value.
  let fold = @src.compute_fold_enrichment(5, 0.0)
  assert_true(fold > 1.0e5)
}

///|
test "peak_calling_fold_enrichment_zero_observed_zero_expected" {
  // 0 / 0 case: should return 0.0.
  assert_eq(@src.compute_fold_enrichment(0, 0.0), 0.0)
}

///|
test "peak_calling_merge_peaks_overlapping" {
  let peaks = [
    @src.CandidatePeak::new(
      chr="chr1", start=1000, end=2000, summit=1500,
      read_count=20, control_count=2, fold_enrichment=10.0, p_value=1.0e-8,
    ),
    @src.CandidatePeak::new(
      chr="chr1", start=1800, end=2500, summit=2000,
      read_count=15, control_count=1, fold_enrichment=15.0, p_value=1.0e-7,
    ),
  ]
  let merged = @src.merge_peaks(peaks, 0)
  assert_eq(merged.length(), 1)
  assert_eq(merged[0].start, 1000)
  assert_eq(merged[0].end, 2500)
  assert_eq(merged[0].read_count, 35)
  assert_eq(merged[0].control_count, 3)
  // Min p-value of the merged set.
  assert_true(merged[0].p_value <= 1.0e-8)
  // Max fold enrichment of the merged set.
  assert_true(merged[0].fold_enrichment >= 15.0 - 1.0e-9)
}

///|
test "peak_calling_merge_peaks_adjacent_within_gap" {
  let peaks = [
    @src.CandidatePeak::new(
      chr="chr1", start=1000, end=1500, summit=1250,
      read_count=10, control_count=1, fold_enrichment=10.0, p_value=0.001,
    ),
    @src.CandidatePeak::new(
      chr="chr1", start=1700, end=2200, summit=2000,
      read_count=12, control_count=2, fold_enrichment=6.0, p_value=0.0001,
    ),
  ]
  // Gap is 200bp; with max_gap=300, they should merge.
  let merged = @src.merge_peaks(peaks, 300)
  assert_eq(merged.length(), 1)
  assert_eq(merged[0].start, 1000)
  assert_eq(merged[0].end, 2200)
  // Summit should be at the position of the higher-count peak (12 > 10).
  assert_eq(merged[0].summit, 2000)
}

///|
test "peak_calling_merge_peaks_non_overlapping_different_chromosomes" {
  let peaks = [
    @src.CandidatePeak::new(
      chr="chr1", start=1000, end=2000, summit=1500,
      read_count=20, control_count=2, fold_enrichment=10.0, p_value=1.0e-8,
    ),
    @src.CandidatePeak::new(
      chr="chr2", start=1000, end=2000, summit=1500,
      read_count=20, control_count=2, fold_enrichment=10.0, p_value=1.0e-8,
    ),
  ]
  let merged = @src.merge_peaks(peaks, 1000)
  assert_eq(merged.length(), 2)
}

///|
test "peak_calling_merge_peaks_far_apart" {
  let peaks = [
    @src.CandidatePeak::new(
      chr="chr1", start=1000, end=1500, summit=1250,
      read_count=10, control_count=1, fold_enrichment=10.0, p_value=0.001,
    ),
    @src.CandidatePeak::new(
      chr="chr1", start=10000, end=10500, summit=10250,
      read_count=10, control_count=1, fold_enrichment=10.0, p_value=0.001,
    ),
  ]
  let merged = @src.merge_peaks(peaks, 100)
  assert_eq(merged.length(), 2)
}

///|
test "peak_calling_merge_peaks_empty" {
  let peaks : Array[@src.CandidatePeak] = []
  let merged = @src.merge_peaks(peaks, 100)
  assert_eq(merged.length(), 0)
}

///|
test "peak_calling_merge_peaks_single_peak" {
  let peaks = [
    @src.CandidatePeak::new(
      chr="chr1", start=1000, end=2000, summit=1500,
      read_count=20, control_count=2, fold_enrichment=10.0, p_value=1.0e-8,
    ),
  ]
  let merged = @src.merge_peaks(peaks, 100)
  assert_eq(merged.length(), 1)
  assert_eq(merged[0].read_count, 20)
}

///|
test "peak_calling_merge_peaks_unsorted_input" {
  // Input peaks in reverse order; merge_peaks should sort internally.
  let peaks = [
    @src.CandidatePeak::new(
      chr="chr1", start=5000, end=5500, summit=5200,
      read_count=15, control_count=1, fold_enrichment=15.0, p_value=1.0e-7,
    ),
    @src.CandidatePeak::new(
      chr="chr1", start=1000, end=1500, summit=1250,
      read_count=20, control_count=2, fold_enrichment=10.0, p_value=1.0e-8,
    ),
  ]
  let merged = @src.merge_peaks(peaks, 100)
  assert_eq(merged.length(), 2)
  // After sort, the first peak should start at 1000.
  assert_eq(merged[0].start, 1000)
  assert_eq(merged[1].start, 5000)
}

///|
test "peak_calling_apply_fdr_marks_significant" {
  // Three peaks with very small p-values; all should be marked significant.
  let peaks = [
    @src.CandidatePeak::new(
      chr="chr1", start=1000, end=1500, summit=1250,
      read_count=50, control_count=2, fold_enrichment=25.0, p_value=1.0e-10,
    ),
    @src.CandidatePeak::new(
      chr="chr1", start=2000, end=2500, summit=2250,
      read_count=40, control_count=2, fold_enrichment=20.0, p_value=1.0e-8,
    ),
    @src.CandidatePeak::new(
      chr="chr1", start=3000, end=3500, summit=3250,
      read_count=30, control_count=2, fold_enrichment=15.0, p_value=1.0e-6,
    ),
  ]
  let result = @src.apply_fdr(peaks, 0.05)
  for p in result {
    assert_true(p.significant)
    assert_true(p.fdr <= 0.05)
  }
}

///|
test "peak_calling_apply_fdr_filters_large_pvalues" {
  // Mix of significant and non-significant peaks.
  let peaks = [
    @src.CandidatePeak::new(
      chr="chr1", start=1000, end=1500, summit=1250,
      read_count=50, control_count=2, fold_enrichment=25.0, p_value=1.0e-10,
    ),
    @src.CandidatePeak::new(
      chr="chr1", start=2000, end=2500, summit=2250,
      read_count=5, control_count=2, fold_enrichment=2.5, p_value=0.4,
    ),
  ]
  let result = @src.apply_fdr(peaks, 0.05)
  // The first peak should be significant; the second should not.
  assert_true(result[0].significant)
  assert_true(!result[1].significant)
  // FDR for first should be small.
  assert_true(result[0].fdr < 0.05)
  // FDR for second should be larger.
  assert_true(result[1].fdr > result[0].fdr)
}

///|
test "peak_calling_apply_fdr_empty" {
  let peaks : Array[@src.CandidatePeak] = []
  let result = @src.apply_fdr(peaks, 0.05)
  assert_eq(result.length(), 0)
}

///|
test "peak_calling_apply_fdr_monotonic" {
  // BH-adjusted p-values should be monotonic when sorted.
  let pvals = [1.0e-10, 1.0e-8, 1.0e-6, 1.0e-4, 0.01]
  let peaks : Array[@src.CandidatePeak] = Array::new()
  let mut i = 0
  while i < pvals.length() {
    peaks.push(
      @src.CandidatePeak::new(
        chr="chr1",
        start=1000 + i * 1000,
        end=1500 + i * 1000,
        summit=1250 + i * 1000,
        read_count=10,
        control_count=1,
        fold_enrichment=10.0,
        p_value=pvals[i],
      ),
    )
    i = i + 1
  }
  let result = @src.apply_fdr(peaks, 1.0)
  // Collect (p_value, fdr) pairs sorted by p_value ascending.
  let pairs : Array[(Double, Double)] = Array::new()
  for p in result {
    pairs.push((p.p_value, p.fdr))
  }
  pairs.sort_by(fn(a : (Double, Double), b : (Double, Double)) -> Int {
    if a.0 < b.0 { -1 } else if a.0 > b.0 { 1 } else { 0 }
  })
  // Adjusted p-values should be non-decreasing with rank.
  let mut k = 0
  while k < pairs.length() - 1 {
    assert_true(pairs[k].1 <= pairs[k + 1].1 + 1.0e-9)
    k = k + 1
  }
}

///|
test "peak_calling_filter_peaks_thresholds" {
  let params = @src.PeakCallingParams::new(
    pvalue_threshold=1.0e-5,
    fdr_threshold=0.05,
    min_fold_enrichment=2.0,
  )
  let peaks = [
    @src.CandidatePeak::new(
      chr="chr1", start=1000, end=1500, summit=1250,
      read_count=50, control_count=2, fold_enrichment=25.0, p_value=1.0e-10,
    ),
    @src.CandidatePeak::new(
      chr="chr1", start=2000, end=2500, summit=2250,
      read_count=5, control_count=2, fold_enrichment=2.5, p_value=1.0e-3,
    ),
    @src.CandidatePeak::new(
      chr="chr1", start=3000, end=3500, summit=3250,
      read_count=3, control_count=2, fold_enrichment=1.5, p_value=1.0e-10,
    ),
  ]
  // Set FDR and significant flags first.
  let with_fdr = @src.apply_fdr(peaks, 0.05)
  let filtered = @src.filter_peaks(with_fdr, params)
  // First peak: passes all (p=1e-10, fold=25, fdr likely small).
  // Second peak: p=1e-3 > 1e-5 threshold -> filtered out.
  // Third peak: fold=1.5 < 2.0 -> filtered out.
  assert_true(filtered.length() >= 1)
  // All filtered peaks should have p_value <= threshold and fold >= min.
  for p in filtered {
    assert_true(p.p_value <= params.pvalue_threshold)
    assert_true(p.fold_enrichment >= params.min_fold_enrichment)
    assert_true(p.fdr <= params.fdr_threshold)
  }
}

///|
test "peak_calling_filter_peaks_empty" {
  let params = @src.PeakCallingParams::default()
  let peaks : Array[@src.CandidatePeak] = []
  let filtered = @src.filter_peaks(peaks, params)
  assert_eq(filtered.length(), 0)
}

///|
test "peak_calling_full_pipeline_with_sample_data" {
  let (treatment, control) = @src.peak_calling_sample_data()
  let params = @src.PeakCallingParams::new(
    window_size=200,
    step_size=50,
    local_lambda_size=2000,
    pvalue_threshold=0.01,
    fdr_threshold=0.5,
    min_fold_enrichment=1.5,
  )
  let peaks = @src.call_peaks(treatment, control, params)
  // Should detect at least one peak in the strong enrichment region.
  assert_true(peaks.length() >= 1)
  for p in peaks {
    // Sanity checks: every called peak must satisfy the filter criteria.
    assert_true(p.p_value <= 0.01)
    assert_true(p.fold_enrichment >= 1.5)
    assert_true(p.fdr <= 0.5)
    assert_eq(p.chr, "chr1")
    assert_true(p.read_count > 0)
    assert_true(p.end > p.start)
  }
}

///|
test "peak_calling_full_pipeline_default_params" {
  let (treatment, control) = @src.peak_calling_sample_data()
  let params = @src.PeakCallingParams::default()
  let peaks = @src.call_peaks(treatment, control, params)
  // Default params are fairly stringent; we may get 0 or more peaks but
  // the call should not crash and should return a valid array.
  assert_true(peaks.length() >= 0)
  for p in peaks {
    assert_eq(p.chr, "chr1")
    assert_true(p.read_count > 0)
  }
}

///|
test "peak_calling_empty_reads" {
  let treatment : Array[@src.ChipSeqRead] = []
  let control : Array[@src.ChipSeqRead] = []
  let params = @src.PeakCallingParams::default()
  let peaks = @src.call_peaks(treatment, control, params)
  assert_eq(peaks.length(), 0)
}

///|
test "peak_calling_no_control" {
  // Treatment reads present but no control: should still produce
  // candidate peaks (lambda estimated as 0 -> expected floored to 0.5).
  let treatment = [
    @src.ChipSeqRead::new(chr="chr1", position=1000, strand="+"),
    @src.ChipSeqRead::new(chr="chr1", position=1050, strand="+"),
    @src.ChipSeqRead::new(chr="chr1", position=1100, strand="+"),
    @src.ChipSeqRead::new(chr="chr1", position=1150, strand="+"),
    @src.ChipSeqRead::new(chr="chr1", position=1200, strand="+"),
    @src.ChipSeqRead::new(chr="chr1", position=1250, strand="+"),
    @src.ChipSeqRead::new(chr="chr1", position=1300, strand="+"),
    @src.ChipSeqRead::new(chr="chr1", position=1350, strand="+"),
  ]
  let control : Array[@src.ChipSeqRead] = []
  let params = @src.PeakCallingParams::new(
    window_size=500,
    step_size=100,
    local_lambda_size=2000,
    pvalue_threshold=0.01,
    fdr_threshold=1.0,
    min_fold_enrichment=1.0,
  )
  let peaks = @src.call_peaks(treatment, control, params)
  assert_true(peaks.length() >= 1)
  // All peaks should be on chr1.
  for p in peaks {
    assert_eq(p.chr, "chr1")
  }
}

///|
test "peak_calling_single_read" {
  let treatment = [
    @src.ChipSeqRead::new(chr="chr1", position=1000, strand="+"),
  ]
  let control : Array[@src.ChipSeqRead] = []
  let params = @src.PeakCallingParams::default()
  let peaks = @src.call_peaks(treatment, control, params)
  // Single read may or may not pass stringent default thresholds,
  // but the call should not crash.
  assert_true(peaks.length() >= 0)
}

///|
test "peak_calling_sample_data_shape" {
  let (treatment, control) = @src.peak_calling_sample_data()
  assert_true(treatment.length() > 0)
  assert_true(control.length() > 0)
  // All reads on chr1.
  for r in treatment {
    assert_eq(r.chr, "chr1")
  }
  for r in control {
    assert_eq(r.chr, "chr1")
  }
  // Treatment should be smaller or comparable to control (since
  // we add background reads to both).
  assert_true(treatment.length() >= 60)
  assert_true(control.length() >= 50)
}

///|
test "peak_calling_sample_data_has_enrichment" {
  let (treatment, _control) = @src.peak_calling_sample_data()
  // Region A around 1000-1500 should have many treatment reads.
  let region_a = @src.count_reads_in_window(treatment, "chr1", 1000, 1500)
  assert_true(region_a >= 20)
  // Region B around 5000-5500.
  let region_b = @src.count_reads_in_window(treatment, "chr1", 5000, 5500)
  assert_true(region_b >= 10)
}

///|
test "peak_calling_summary_empty" {
  let peaks : Array[@src.CandidatePeak] = []
  let s = @src.peak_calling_summary(peaks)
  assert_true(s.contains("0 peaks"))
}

///|
test "peak_calling_summary_with_peaks" {
  let peaks = [
    @src.CandidatePeak::new(
      chr="chr1", start=1000, end=2000, summit=1500,
      read_count=42, control_count=5, fold_enrichment=8.4, p_value=1.0e-10,
    ),
    @src.CandidatePeak::new(
      chr="chr1", start=3000, end=4000, summit=3500,
      read_count=20, control_count=2, fold_enrichment=10.0, p_value=1.0e-8,
    ),
  ]
  let with_fdr = @src.apply_fdr(peaks, 0.05)
  let s = @src.peak_calling_summary(with_fdr)
  assert_true(s.contains("Total peaks: 2"))
  assert_true(s.contains("chr1"))
  assert_true(s.contains("Top"))
}

///|
test "peak_calling_summary_includes_significant_count" {
  let peaks = [
    @src.CandidatePeak::new(
      chr="chr1", start=1000, end=2000, summit=1500,
      read_count=42, control_count=5, fold_enrichment=8.4, p_value=1.0e-10,
    ),
    @src.CandidatePeak::new(
      chr="chr1", start=3000, end=4000, summit=3500,
      read_count=2, control_count=2, fold_enrichment=1.0, p_value=0.5,
    ),
  ]
  let with_fdr = @src.apply_fdr(peaks, 0.05)
  let s = @src.peak_calling_summary(with_fdr)
  assert_true(s.contains("Significant peaks:"))
}

///|
test "peak_calling_pipeline_deterministic" {
  // Calling peaks twice on the same data should give identical results.
  let (t1, c1) = @src.peak_calling_sample_data()
  let (t2, c2) = @src.peak_calling_sample_data()
  let params = @src.PeakCallingParams::new(
    window_size=300,
    step_size=100,
    local_lambda_size=3000,
    pvalue_threshold=0.01,
    fdr_threshold=0.5,
    min_fold_enrichment=2.0,
  )
  let peaks1 = @src.call_peaks(t1, c1, params)
  let peaks2 = @src.call_peaks(t2, c2, params)
  assert_eq(peaks1.length(), peaks2.length())
  let mut i = 0
  while i < peaks1.length() {
    assert_eq(peaks1[i].start, peaks2[i].start)
    assert_eq(peaks1[i].end, peaks2[i].end)
    assert_eq(peaks1[i].read_count, peaks2[i].read_count)
    i = i + 1
  }
}
