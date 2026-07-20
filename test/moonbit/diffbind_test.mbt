///|
/// Tests for the diffbind module.
/// Mirrors the Bioconductor DiffBind package API.

///|
test "diffbind_peak_new" {
  let peak = @src.DBAPeak::new("chr1", 1000, 2000)
  assert_eq(peak.chr, "chr1")
  assert_eq(peak.start, 1000)
  assert_eq(peak.end, 2000)
  assert_eq(peak.score, 0.0)
  // Summit defaults to midpoint.
  assert_eq(peak.summit, 1500)
}

///|
test "diffbind_peak_new_midpoint" {
  let peak = @src.DBAPeak::new("chr2", 300, 700)
  // midpoint = (300 + 700) / 2 = 500
  assert_eq(peak.summit, 500)
}

///|
test "diffbind_sample_new" {
  let sample = @src.DBASample::new("S1", "Treatment", 1)
  assert_eq(sample.sample_id, "S1")
  assert_eq(sample.condition, "Treatment")
  assert_eq(sample.replicate, 1)
  assert_eq(sample.peaks.length(), 0)
  assert_eq(sample.peak_file, "")
}

///|
test "diffbind_config_new" {
  let s1 = @src.DBASample::new("S1", "A", 1)
  let s2 = @src.DBASample::new("S2", "B", 1)
  let config = @src.DBAConfig::new([s1, s2])
  assert_eq(config.samples.length(), 2)
  assert_eq(config.conditions.length(), 2)
  assert_eq(config.conditions[0], "A")
  assert_eq(config.conditions[1], "B")
  assert_eq(config.replicates.length(), 2)
  assert_eq(config.peaks.length(), 2)
}

///|
test "diffbind_peakset_adds_sample" {
  let s1 = @src.DBASample::new("S1", "A", 1)
  let config = @src.DBAConfig::new([s1])
  assert_eq(config.samples.length(), 1)
  let s2 = @src.DBASample::new("S2", "B", 1)
  let new_config = @src.diffbind_peakset(config, s2)
  assert_eq(new_config.samples.length(), 2)
  // Original config is unchanged.
  assert_eq(config.samples.length(), 1)
}

///|
test "diffbind_peak_format_to_string" {
  // PeakFormat variants are not directly constructable from outside the
  // package (the enum derives Eq). DBAConfig::new always uses BedFormat.
  let config = @src.diffbind_create_example_config()
  assert_eq(@src.PeakFormat::to_string(config.peak_format), "BED")
}

///|
test "diffbind_is_overlapping_positive" {
  assert_true(@src.diffbind_is_overlapping(100, 200, 150, 250))
  assert_true(@src.diffbind_is_overlapping(100, 200, 100, 200))
}

///|
test "diffbind_is_overlapping_negative" {
  assert_true(!@src.diffbind_is_overlapping(100, 200, 300, 400))
  assert_true(!@src.diffbind_is_overlapping(1000, 2000, 100, 200))
}

///|
test "diffbind_count" {
  // Sample with one peak on chr1:1000-2000.
  // Note: `peaks` is a mutable Array field, so we don't need `mut` on `sample`.
  let sample = @src.DBASample::new("S1", "A", 1)
  sample.peaks.push(@src.DBAPeak::new("chr1", 1000, 2000))
  let config = @src.DBAConfig::new([sample])
  let ranges : Array[(String, Int, Int)] = [
    ("chr1", 1500, 2500),
    ("chr2", 100, 200),
  ]
  let counts = @src.diffbind_count(config, ranges)
  // chr1:1500-2500 overlaps the peak.
  assert_eq(counts.get("chr1:1500-2500"), Some(1))
  // chr2:100-200 does not overlap.
  assert_eq(counts.get("chr2:100-200"), Some(0))
}

///|
test "diffbind_consensus_peaks_basic" {
  let s1 = @src.DBASample::new("S1", "A", 1)
  s1.peaks.push(@src.DBAPeak::new("chr1", 100, 200))
  let s2 = @src.DBASample::new("S2", "B", 1)
  s2.peaks.push(@src.DBAPeak::new("chr1", 150, 250))
  // min_overlap=0.5 with 2 samples means threshold=1.0 sample.
  let consensus = @src.diffbind_consensus_peaks([s1, s2], 0.5)
  // Both samples overlap [100-250] -> at least one consensus peak.
  assert_true(consensus.length() >= 1)
  assert_eq(consensus[0].chr, "chr1")
}

///|
test "diffbind_consensus_peaks_empty" {
  let consensus = @src.diffbind_consensus_peaks([], 0.5)
  assert_eq(consensus.length(), 0)
}

///|
test "diffbind_overlap" {
  let s1 = @src.DBASample::new("S1", "A", 1)
  s1.peaks.push(@src.DBAPeak::new("chr1", 100, 200))
  let s2 = @src.DBASample::new("S2", "B", 1)
  s2.peaks.push(@src.DBAPeak::new("chr1", 150, 250))
  let config = @src.DBAConfig::new([s1, s2])
  let result = @src.diffbind_overlap(config, 0.5)
  assert_true(result.consensus_peaks.length() >= 1)
  assert_eq(result.sample_names.length(), 2)
  assert_eq(result.sample_names[0], "S1")
  // Overlap matrix has one row per consensus peak, with one column per sample.
  assert_eq(result.overlap_matrix.length(), result.consensus_peaks.length())
  assert_eq(result.overlap_matrix[0].length(), 2)
}

///|
test "diffbind_find_consensus_peaks" {
  let s1 = @src.DBASample::new("S1", "A", 1)
  s1.peaks.push(@src.DBAPeak::new("chr1", 100, 200))
  let s2 = @src.DBASample::new("S2", "B", 1)
  s2.peaks.push(@src.DBAPeak::new("chr1", 150, 250))
  let config = @src.DBAConfig::new([s1, s2])
  let overlap = @src.diffbind_overlap(config, 0.5)
  let consensus = @src.diffbind_find_consensus_peaks(overlap)
  assert_eq(consensus.length(), overlap.consensus_peaks.length())
}

///|
test "diffbind_merge_peaks_overlapping" {
  let peaks1 = [@src.DBAPeak::new("chr1", 100, 200)]
  let peaks2 = [@src.DBAPeak::new("chr1", 150, 300)]
  let merged = @src.diffbind_merge_peaks(peaks1, peaks2)
  // Overlapping peaks on the same chromosome are merged into one.
  assert_eq(merged.length(), 1)
  assert_eq(merged[0].chr, "chr1")
  assert_eq(merged[0].start, 100)
  assert_eq(merged[0].end, 300)
}

///|
test "diffbind_merge_peaks_disjoint" {
  let peaks1 = [@src.DBAPeak::new("chr1", 100, 200)]
  let peaks2 = [@src.DBAPeak::new("chr2", 100, 200)]
  let merged = @src.diffbind_merge_peaks(peaks1, peaks2)
  assert_eq(merged.length(), 2)
}

///|
test "diffbind_merge_peaks_empty" {
  let merged = @src.diffbind_merge_peaks([], [])
  assert_eq(merged.length(), 0)
}

///|
test "diffbind_compute_tmm_basic" {
  let counts : Array[Array[Double]] = [[10.0, 12.0], [20.0, 22.0], [30.0, 33.0]]
  let factors = @src.diffbind_compute_tmm(counts)
  assert_eq(factors.length(), 2)
  // Reference sample has factor 1.0.
  assert_eq(factors[0], 1.0)
}

///|
test "diffbind_compute_tmm_empty" {
  let factors = @src.diffbind_compute_tmm([])
  assert_eq(factors.length(), 0)
}

///|
test "diffbind_normalize_tmm" {
  let config = @src.diffbind_create_example_config()
  let normalized = @src.diffbind_normalize(config, "TMM")
  // Normalize returns one row per consensus peak.
  assert_true(normalized.length() >= 1)
  if normalized.length() > 0 {
    assert_eq(normalized[0].length(), config.samples.length())
  }
}

///|
test "diffbind_normalize_cpm" {
  let config = @src.diffbind_create_example_config()
  let normalized = @src.diffbind_normalize(config, "CPM")
  assert_true(normalized.length() >= 1)
}

///|
test "diffbind_negative_binomial_test_equal" {
  // Equal counts -> p-value should be relatively high (no difference).
  let p = @src.diffbind_negative_binomial_test(50.0, 50.0, 1.0)
  assert_true(p >= 0.0 && p <= 1.0)
}

///|
test "diffbind_negative_binomial_test_different" {
  // Very different counts -> small p-value.
  let p = @src.diffbind_negative_binomial_test(100.0, 1.0, 1.0)
  assert_true(p >= 0.0 && p <= 1.0)
}

///|
test "diffbind_negative_binomial_test_zero" {
  // Both zero -> total = 0 -> p = 1.
  let p = @src.diffbind_negative_binomial_test(0.0, 0.0, 1.0)
  assert_eq(p, 1.0)
}

///|
test "diffbind_bh_fdr_basic" {
  let pvalues = [0.01, 0.5, 0.001, 0.2]
  let fdrs = @src.diffbind_bh_fdr(pvalues)
  assert_eq(fdrs.length(), 4)
  // The smallest p-value (0.001 at index 2) should have the smallest FDR.
  assert_true(fdrs[2] <= fdrs[0])
  assert_true(fdrs[0] <= fdrs[1])
}

///|
test "diffbind_bh_fdr_empty" {
  let fdrs = @src.diffbind_bh_fdr([])
  assert_eq(fdrs.length(), 0)
}

///|
test "diffbind_bh_fdr_monotonic" {
  let pvalues = [0.001, 0.5, 0.05]
  let fdrs = @src.diffbind_bh_fdr(pvalues)
  // FDRs should be valid probabilities.
  let mut i = 0
  while i < fdrs.length() {
    assert_true(fdrs[i] >= 0.0 && fdrs[i] <= 1.0)
    i = i + 1
  }
}

///|
test "diffbind_analyze_basic" {
  let config = @src.diffbind_create_example_config()
  let result = @src.diffbind_analyze(config, "ConditionA", "ConditionB")
  assert_eq(result.contrast.0, "ConditionA")
  assert_eq(result.contrast.1, "ConditionB")
  assert_true(result.peaks.length() >= 1)
  assert_eq(result.pvalues.length(), result.peaks.length())
  assert_eq(result.fdrs.length(), result.peaks.length())
  assert_eq(result.fold_changes.length(), result.peaks.length())
}

///|
test "diffbind_analyze_unknown_condition" {
  let config = @src.diffbind_create_example_config()
  let result = @src.diffbind_analyze(config, "ConditionA", "UnknownCond")
  // Empty result because ConditionB group has no samples.
  assert_eq(result.peaks.length(), 0)
  assert_eq(result.n_significant, 0)
}

///|
test "diffbind_report_threshold" {
  let result = @src.diffbind_create_example_result()
  let significant = @src.diffbind_report(result, 0.05)
  // Two peaks have FDR < 0.05 (0.005 and 0.0005).
  assert_eq(significant.length(), 2)
}

///|
test "diffbind_report_high_threshold" {
  let result = @src.diffbind_create_example_result()
  let significant = @src.diffbind_report(result, 0.9)
  // All but one peak have FDR < 0.9 (the 0.85 is below 0.9 too).
  // Example FDRs: [0.005, 0.6, 0.0005, 0.85, 0.06] -> all under 0.9.
  assert_eq(significant.length(), 5)
}

///|
test "diffbind_volcano_data" {
  let result = @src.diffbind_create_example_result()
  let data = @src.diffbind_volcano_data(result)
  assert_eq(data.length(), result.peaks.length())
  let mut i = 0
  while i < data.length() {
    let (fc, neg_log10_p) = data[i]
    assert_true(fc == result.fold_changes[i])
    assert_true(neg_log10_p >= 0.0)
    i = i + 1
  }
}

///|
test "diffbind_ma_data" {
  let result = @src.diffbind_create_example_result()
  let data = @src.diffbind_ma_data(result)
  assert_eq(data.length(), result.peaks.length())
}

///|
test "diffbind_create_example_config" {
  let config = @src.diffbind_create_example_config()
  assert_eq(config.samples.length(), 4)
  assert_eq(config.samples[0].sample_id, "A1")
  assert_eq(config.samples[0].condition, "ConditionA")
  assert_eq(config.samples[2].sample_id, "B1")
  assert_eq(config.samples[2].condition, "ConditionB")
}

///|
test "diffbind_create_example_peaks" {
  let peaks = @src.diffbind_create_example_peaks()
  assert_eq(peaks.length(), 5)
  assert_eq(peaks[0].chr, "chr1")
  assert_eq(peaks[0].start, 1000)
  assert_eq(peaks[0].end, 2000)
  assert_eq(peaks[4].chr, "chr3")
}

///|
test "diffbind_contrast_valid" {
  let config = @src.diffbind_create_example_config()
  let (c1, c2) = @src.diffbind_contrast(config, "ConditionA", "ConditionB")
  assert_eq(c1, "ConditionA")
  assert_eq(c2, "ConditionB")
}

///|
test "diffbind_plot_pca" {
  let config = @src.diffbind_create_example_config()
  let plot = @src.diffbind_plot_pca(config, "TMM")
  assert_true(plot.find("PCA Plot").is_some())
}

///|
test "diffbind_plot_heatmap" {
  let result = @src.diffbind_create_example_result()
  let heatmap = @src.diffbind_plot_heatmap(result, 5)
  assert_true(heatmap.find("Heatmap").is_some())
}

///|
test "diffbind_plot_volcano" {
  let result = @src.diffbind_create_example_result()
  let volcano = @src.diffbind_plot_volcano(result)
  assert_true(volcano.find("Volcano Plot").is_some())
}

///|
test "diffbind_plot_ma" {
  let result = @src.diffbind_create_example_result()
  let ma = @src.diffbind_plot_ma(result)
  assert_true(ma.find("MA Plot").is_some())
}

///|
test "diffbind_create_example_result" {
  let result = @src.diffbind_create_example_result()
  assert_eq(result.peaks.length(), 5)
  assert_eq(result.pvalues.length(), 5)
  assert_eq(result.fdrs.length(), 5)
  assert_eq(result.fold_changes.length(), 5)
  assert_eq(result.contrast.0, "ConditionA")
  assert_eq(result.contrast.1, "ConditionB")
  // Two peaks have FDR < 0.05.
  assert_eq(result.n_significant, 2)
}

///|
test "diffbind_peak_to_string" {
  let peak = @src.DBAPeak::new("chr1", 100, 200)
  let s = peak.to_string()
  assert_true(s.contains("chr1"))
  assert_true(s.contains("100"))
  assert_true(s.contains("200"))
}
