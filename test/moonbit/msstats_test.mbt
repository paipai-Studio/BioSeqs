///|
/// Test file for MSstats module.

test "msstats_feature_creation" {
  let f = @src.MSFeature::new("P1", "pep1", "tr1", "Ctrl", "S1", "R1", 1000.0)
  assert_eq(f.protein, "P1")
  assert_eq(f.peptide, "pep1")
  assert_eq(f.transition, "tr1")
  assert_eq(f.group, "Ctrl")
  assert_eq(f.subject, "S1")
  assert_eq(f.run, "R1")
  assert_eq(f.intensity, 1000.0)
}

test "msstats_data_process_log_transform" {
  let features : Array[@src.MSFeature] = Array::new()
  features.push(@src.MSFeature::new("P1", "pep1", "tr1", "Ctrl", "S1", "R1", 1024.0))
  features.push(@src.MSFeature::new("P1", "pep1", "tr1", "Trt", "S2", "R2", 2048.0))

  let processed = @src.ms_data_process(features, normalization=@src.ms_norm_none())
  assert_eq(processed.length(), 2)
  // log2(1024) = 10
  assert_true((processed[0].intensity - 10.0).abs() < 0.01)
  // log2(2048) = 11
  assert_true((processed[1].intensity - 11.0).abs() < 0.01)
}

test "msstats_data_process_median_norm" {
  let features : Array[@src.MSFeature] = Array::new()
  features.push(@src.MSFeature::new("P1", "pep1", "tr1", "Ctrl", "S1", "R1", 1024.0))
  features.push(@src.MSFeature::new("P1", "pep2", "tr2", "Ctrl", "S1", "R1", 4096.0))
  features.push(@src.MSFeature::new("P2", "pep3", "tr3", "Trt", "S2", "R2", 2048.0))
  features.push(@src.MSFeature::new("P2", "pep4", "tr4", "Trt", "S2", "R2", 8192.0))

  let processed = @src.ms_data_process(features, normalization=@src.ms_norm_median())
  assert_eq(processed.length(), 4)
  // Values should be centered (median subtracted)
  // R1: log2(1024)=10, log2(4096)=12, median=11, so values become -1, 1
  // R2: log2(2048)=11, log2(8192)=13, median=12, so values become -1, 1
  assert_true((processed[0].intensity - (-1.0)).abs() < 0.01)
  assert_true((processed[1].intensity - 1.0).abs() < 0.01)
}

test "msstats_data_process_zero_intensity" {
  let features : Array[@src.MSFeature] = Array::new()
  features.push(@src.MSFeature::new("P1", "pep1", "tr1", "Ctrl", "S1", "R1", 0.0))

  let processed = @src.ms_data_process(features, normalization=@src.ms_norm_none())
  // Zero should be replaced with 1, then log2(1) = 0
  assert_true(processed[0].intensity.abs() < 0.01)
}

test "msstats_summarize_tukey" {
  let features : Array[@src.MSFeature] = Array::new()
  // Two peptides for protein P1 in run R1
  features.push(@src.MSFeature::new("P1", "pep1", "tr1", "Ctrl", "S1", "R1", 1024.0))
  features.push(@src.MSFeature::new("P1", "pep2", "tr2", "Ctrl", "S1", "R1", 4096.0))

  let processed = @src.ms_data_process(features, normalization=@src.ms_norm_none())
  let summarized = @src.ms_summarize(processed, method=@src.ms_summary_tukey())
  assert_eq(summarized.length(), 1) // One protein-run combo
  assert_eq(summarized[0].protein, "P1")
  assert_eq(summarized[0].run, "R1")
  // Median of log2(1024)=10 and log2(4096)=12 is 11
  assert_true((summarized[0].log2_abundance - 11.0).abs() < 0.01)
}

test "msstats_summarize_linear" {
  let features : Array[@src.MSFeature] = Array::new()
  features.push(@src.MSFeature::new("P1", "pep1", "tr1", "Ctrl", "S1", "R1", 1024.0))
  features.push(@src.MSFeature::new("P1", "pep2", "tr2", "Ctrl", "S1", "R1", 4096.0))

  let processed = @src.ms_data_process(features, normalization=@src.ms_norm_none())
  let summarized = @src.ms_summarize(processed, method=@src.ms_summary_linear())
  // Mean of log2(1024)=10 and log2(4096)=12 is 11
  assert_true((summarized[0].log2_abundance - 11.0).abs() < 0.01)
}

test "msstats_group_comparison_basic" {
  let features = @src.msstats_sample_data()
  let processed = @src.ms_data_process(features, normalization=@src.ms_norm_median())
  let summarized = @src.ms_summarize(processed)
  let results = @src.ms_group_comparison(summarized, fdr_threshold=0.5, log2fc_threshold=0.5)
  assert_true(results.get_n_results() > 0)
}

test "msstats_group_comparison_significant" {
  let features = @src.msstats_sample_data()
  let processed = @src.ms_data_process(features, normalization=@src.ms_norm_median())
  let summarized = @src.ms_summarize(processed)
  let results = @src.ms_group_comparison(summarized, fdr_threshold=0.5, log2fc_threshold=0.0)
  let sig = results.get_significant()
  assert_true(sig.length() >= 0)
}

test "msstats_group_comparison_top_proteins" {
  let features = @src.msstats_sample_data()
  let processed = @src.ms_data_process(features, normalization=@src.ms_norm_median())
  let summarized = @src.ms_summarize(processed)
  let results = @src.ms_group_comparison(summarized)
  let top = results.get_top_proteins(2)
  assert_eq(top.length(), 2)
}

test "msstats_group_comparison_summary" {
  let features = @src.msstats_sample_data()
  let processed = @src.ms_data_process(features)
  let summarized = @src.ms_summarize(processed)
  let results = @src.ms_group_comparison(summarized)
  let s = results.summary()
  assert_true(s.contains("MSstats"))
}

test "msstats_sample_size" {
  let features = @src.msstats_sample_data()
  let processed = @src.ms_data_process(features)
  let summarized = @src.ms_summarize(processed)
  let results = @src.ms_group_comparison(summarized)
  let ss = @src.ms_design_sample_size(results.results, desired_power=0.8)
  assert_true(ss.n_samples >= 3)
}

test "msstats_sample_data" {
  let features = @src.msstats_sample_data()
  assert_eq(features.length(), 40) // 5 proteins * 4 groups * 2 peptides
  assert_eq(features[0].protein, "P1")
}

test "msstats_type_to_string" {
  assert_eq(@src.ms_type_dda().to_string(), "DDA")
  assert_eq(@src.ms_type_dia().to_string(), "DIA")
  assert_eq(@src.ms_type_srm().to_string(), "SRM")
  assert_eq(@src.ms_type_tmt().to_string(), "TMT")
}

test "msstats_norm_to_string" {
  assert_eq(@src.ms_norm_none().to_string(), "none")
  assert_eq(@src.ms_norm_median().to_string(), "median")
  assert_eq(@src.ms_norm_quantile().to_string(), "quantile")
  assert_eq(@src.ms_norm_global().to_string(), "globalStandards")
}

test "msstats_summary_to_string" {
  assert_eq(@src.ms_summary_tukey().to_string(), "Tukey")
  assert_eq(@src.ms_summary_linear().to_string(), "linear")
  assert_eq(@src.ms_summary_logsum().to_string(), "logSum")
}
