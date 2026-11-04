///|
/// Tests for the minfi module.
/// Mirrors the Bioconductor minfi package API.

///|
test "minfi_compute_beta_basic" {
  // beta = M / (M + U + offset) = 800 / (800 + 200 + 100) = 800/1100.
  let beta = @src.minfi_compute_beta(800.0, 200.0, 100.0)
  assert_true(beta > 0.72 && beta < 0.73)
}

///|
test "minfi_compute_beta_zero_denominator" {
  // 0 + 0 + 0 = 0 -> denominator is 0 -> returns 0.
  let beta = @src.minfi_compute_beta(0.0, 0.0, 0.0)
  assert_eq(beta, 0.0)
}

///|
test "minfi_compute_beta_clamped" {
  // Negative inputs would make beta negative; function clamps to 0.
  let beta = @src.minfi_compute_beta(-100.0, 0.0, 0.0)
  assert_eq(beta, 0.0)
}

///|
test "minfi_compute_mvalue_basic" {
  // M = log2(800/200) = log2(4) = 2.
  let m = @src.minfi_compute_mvalue(800.0, 200.0)
  assert_true(m > 1.99 && m < 2.01)
}

///|
test "minfi_compute_mvalue_zero_inputs" {
  // Substitutes 1.0 for any input below 1.0 -> log2(1) = 0.
  let m = @src.minfi_compute_mvalue(0.0, 0.0)
  assert_eq(m, 0.0)
}

///|
test "minfi_mean_basic" {
  let arr = [1.0, 2.0, 3.0, 4.0, 5.0]
  assert_eq(@src.minfi_mean(arr), 3.0)
}

///|
test "minfi_mean_empty" {
  assert_eq(@src.minfi_mean([]), 0.0)
}

///|
test "minfi_median_odd" {
  // Odd count: median = middle element.
  let arr = [3.0, 1.0, 2.0, 5.0, 4.0]
  assert_eq(@src.minfi_median(arr), 3.0)
}

///|
test "minfi_median_even" {
  // Even count: median = average of two middle values.
  let arr = [1.0, 2.0, 3.0, 4.0]
  assert_eq(@src.minfi_median(arr), 2.5)
}

///|
test "minfi_variance_basic" {
  // Population variance of [1,2,3,4,5] mean=3 -> 2.5. Sample variance uses
  // n-1 -> 2.5*5/(5-1) = 2.5.
  let arr = [1.0, 2.0, 3.0, 4.0, 5.0]
  let v = @src.minfi_variance(arr)
  assert_true(v > 2.49 && v < 2.51)
}

///|
test "minfi_variance_single" {
  assert_eq(@src.minfi_variance([5.0]), 0.0)
}

///|
test "minfi_bh_fdr_basic" {
  let pvalues = [0.01, 0.5, 0.001, 0.2]
  let fdrs = @src.minfi_bh_fdr(pvalues)
  assert_eq(fdrs.length(), 4)
  // The smallest p-value (index 2) gets the smallest FDR.
  assert_true(fdrs[2] <= fdrs[0])
  // Adjusted values are valid probabilities.
  let mut i = 0
  while i < fdrs.length() {
    assert_true(fdrs[i] >= 0.0 && fdrs[i] <= 1.0)
    i = i + 1
  }
}

///|
test "minfi_bh_fdr_empty" {
  let fdrs = @src.minfi_bh_fdr([])
  assert_eq(fdrs.length(), 0)
}

///|
test "minfi_t_test_basic" {
  let x = [1.0, 2.0, 3.0, 4.0, 5.0]
  let y = [10.0, 11.0, 12.0, 13.0, 14.0]
  let (t, df) = @src.minfi_t_test(x, y)
  assert_eq(df, 8.0)
  // Mean difference is -9, t is negative.
  assert_true(t < 0.0)
}

///|
test "minfi_t_test_too_short" {
  // < 2 samples in either group -> (0, 0).
  let (t, df) = @src.minfi_t_test([1.0], [2.0, 3.0])
  assert_eq(t, 0.0)
  assert_eq(df, 0.0)
}

///|
test "minfi_t_test_equal_groups" {
  let x = [1.0, 2.0, 3.0]
  let y = [1.0, 2.0, 3.0]
  let (t, _df) = @src.minfi_t_test(x, y)
  // Same groups -> t = 0.
  assert_eq(t, 0.0)
}

///|
test "minfi_rgchannelset_new" {
  let red = [[100.0, 200.0], [300.0, 400.0]]
  let green = [[500.0, 600.0], [700.0, 800.0]]
  let rgset = @src.RGChannelSet::new(red, green, ["s1", "s2"], ["p1", "p2"])
  assert_eq(rgset.n_probes(), 2)
  assert_eq(rgset.n_samples(), 2)
  assert_eq(rgset.sample_names.length(), 2)
  assert_eq(rgset.probe_names.length(), 2)
}

///|
test "minfi_methylset_new" {
  let meth = [[100.0, 200.0], [300.0, 400.0]]
  let unmeth = [[10.0, 20.0], [30.0, 40.0]]
  let mset = @src.MethylSet::new(meth, unmeth, ["s1", "s2"], ["p1", "p2"])
  assert_eq(mset.n_probes(), 2)
  assert_eq(mset.n_samples(), 2)
}

///|
test "minfi_genomicratioset_new" {
  let beta = [[0.5, 0.6], [0.7, 0.8]]
  let mvalues = [[1.0, 2.0], [3.0, 4.0]]
  let grset = @src.GenomicRatioSet::new(beta, mvalues, ["s1", "s2"], [
    "p1", "p2",
  ])
  assert_eq(grset.n_probes(), 2)
  assert_eq(grset.n_samples(), 2)
}

///|
test "minfi_create_example_rgset" {
  let rgset = @src.minfi_create_example_rgset()
  assert_eq(rgset.n_probes(), 15)
  assert_eq(rgset.n_samples(), 4)
  assert_eq(rgset.sample_names.length(), 4)
  assert_eq(rgset.probe_names.length(), 15)
}

///|
test "minfi_create_example_grset" {
  let grset = @src.minfi_create_example_grset()
  assert_eq(grset.n_probes(), 15)
  assert_eq(grset.n_samples(), 4)
  // Beta matrix shape: probes x samples.
  assert_eq(grset.beta.length(), 15)
  assert_eq(grset.beta[0].length(), 4)
}

///|
test "minfi_create_example_dmp_result" {
  let dmp = @src.minfi_create_example_dmp_result()
  assert_eq(dmp.probe_names.length(), 10)
  assert_eq(dmp.chr.length(), 10)
  assert_eq(dmp.pos.length(), 10)
  assert_eq(dmp.pvalues.length(), 10)
  assert_eq(dmp.fdrs.length(), 10)
  assert_eq(dmp.beta_diffs.length(), 10)
  // 8 out of 10 have FDR < 0.05.
  assert_eq(dmp.n_dmp, 8)
}

///|
test "minfi_preprocess_noob" {
  let rgset = @src.minfi_create_example_rgset()
  let mset = @src.minfi_preprocess_noob(rgset)
  assert_eq(mset.n_probes(), rgset.n_probes())
  assert_eq(mset.n_samples(), rgset.n_samples())
  // Corrected values are floored at 1.0.
  let mut i = 0
  while i < mset.n_probes() {
    let mut j = 0
    while j < mset.n_samples() {
      assert_true(mset.methylated[i][j] >= 1.0)
      assert_true(mset.unmethylated[i][j] >= 1.0)
      j = j + 1
    }
    i = i + 1
  }
}

///|
test "minfi_preprocess_illumina" {
  let rgset = @src.minfi_create_example_rgset()
  let mset = @src.minfi_preprocess_illumina(rgset)
  assert_eq(mset.n_probes(), rgset.n_probes())
  assert_eq(mset.n_samples(), rgset.n_samples())
}

///|
test "minfi_preprocess_quantile" {
  let rgset = @src.minfi_create_example_rgset()
  let mset = @src.minfi_preprocess_quantile(rgset)
  assert_eq(mset.n_probes(), rgset.n_probes())
  assert_eq(mset.n_samples(), rgset.n_samples())
}

///|
test "minfi_preprocess_funnorm" {
  let rgset = @src.minfi_create_example_rgset()
  let grset = @src.minfi_preprocess_funnorm(rgset)
  assert_eq(grset.n_probes(), rgset.n_probes())
  assert_eq(grset.n_samples(), rgset.n_samples())
}

///|
test "minfi_ratio_convert" {
  let rgset = @src.minfi_create_example_rgset()
  let mset = @src.minfi_preprocess_noob(rgset)
  let grset = @src.minfi_ratio_convert(mset)
  assert_eq(grset.n_probes(), mset.n_probes())
  assert_eq(grset.n_samples(), mset.n_samples())
  // Beta values must be in [0, 1].
  let mut i = 0
  while i < grset.n_probes() {
    let mut j = 0
    while j < grset.n_samples() {
      assert_true(grset.beta[i][j] >= 0.0 && grset.beta[i][j] <= 1.0)
      j = j + 1
    }
    i = i + 1
  }
}

///|
test "minfi_get_beta" {
  let grset = @src.minfi_create_example_grset()
  let beta = @src.minfi_get_beta(grset)
  assert_eq(beta.length(), 15)
  assert_eq(beta[0].length(), 4)
}

///|
test "minfi_get_mvalues" {
  let grset = @src.minfi_create_example_grset()
  let mvalues = @src.minfi_get_mvalues(grset)
  assert_eq(mvalues.length(), 15)
  assert_eq(mvalues[0].length(), 4)
}

///|
test "minfi_detection_pvalue" {
  let rgset = @src.minfi_create_example_rgset()
  let pvalues = @src.minfi_detection_pvalue(rgset)
  assert_eq(pvalues.length(), rgset.n_probes())
  // Each p-value should be in [0, 1].
  let mut i = 0
  while i < pvalues.length() {
    assert_true(pvalues[i] >= 0.0 && pvalues[i] <= 1.0)
    i = i + 1
  }
}

///|
test "minfi_get_qc" {
  let rgset = @src.minfi_create_example_rgset()
  let qc = @src.minfi_get_qc(rgset)
  assert_eq(qc.sample_names.length(), rgset.n_samples())
  assert_eq(qc.median_meth.length(), rgset.n_samples())
  assert_eq(qc.median_unmeth.length(), rgset.n_samples())
  assert_eq(qc.detection_pvalues.length(), rgset.n_samples())
  assert_eq(qc.failed_probes.length(), rgset.n_samples())
  assert_eq(qc.bisulfite_conversion.length(), rgset.n_samples())
  assert_eq(qc.sex_predictions.length(), rgset.n_samples())
}

///|
test "minfi_estimate_sex" {
  let grset = @src.minfi_create_example_grset()
  let sexes = @src.minfi_estimate_sex(grset)
  assert_eq(sexes.length(), grset.n_samples())
  // The chrX probes for samples 0 and 1 have high beta (around 0.5).
  let mut i = 0
  while i < sexes.length() {
    assert_true(sexes[i] == "M" || sexes[i] == "F")
    i = i + 1
  }
}

///|
test "minfi_find_dmps" {
  let grset = @src.minfi_create_example_grset()
  // samples 0,1 -> group 0; samples 2,3 -> group 1.
  let groups = [0, 0, 1, 1]
  let result = @src.minfi_find_dmps(grset, groups)
  assert_eq(result.probe_names.length(), grset.n_probes())
  assert_eq(result.pvalues.length(), grset.n_probes())
  assert_eq(result.fdrs.length(), grset.n_probes())
  assert_eq(result.beta_diffs.length(), grset.n_probes())
}

///|
test "minfi_find_dmps_wrong_groups" {
  let grset = @src.minfi_create_example_grset()
  // Wrong number of groups (expected 4) -> empty result.
  let result = @src.minfi_find_dmps(grset, [0, 0, 1])
  assert_eq(result.probe_names.length(), 0)
}

///|
test "minfi_find_dmrs" {
  let dmp = @src.minfi_create_example_dmp_result()
  let dmr = @src.minfi_find_dmrs(dmp, 500)
  // At least one region should be formed from the cluster on chr1 around
  // 1000-1500 and the cluster around 5000-5200 and 9000-9100.
  assert_true(dmr.n_dmr > 0)
}

///|
test "minfi_filter_probes_low_threshold" {
  let grset = @src.minfi_create_example_grset()
  let filtered = @src.minfi_filter_probes(grset, 0.0)
  // No probes have stdev < 0.
  assert_eq(filtered.n_probes(), grset.n_probes())
}

///|
test "minfi_filter_probes_high_threshold" {
  let grset = @src.minfi_create_example_grset()
  let filtered = @src.minfi_filter_probes(grset, 100.0)
  // No probes have such large stdev; all should be removed.
  assert_eq(filtered.n_probes(), 0)
}

///|
test "minfi_drop_snp_probes" {
  // Build a GRSet with one rs-prefixed probe and one cg-prefixed probe.
  let beta = [[0.5, 0.6], [0.7, 0.8]]
  let mvalues = [[0.0, 1.0], [2.0, 3.0]]
  let grset = @src.GenomicRatioSet::new(beta, mvalues, ["s1", "s2"], [
    "cg00000001_chr1_1000", "rs1234_chr1_2000",
  ])
  let filtered = @src.minfi_drop_snp_probes(grset)
  assert_eq(filtered.n_probes(), 1)
  assert_eq(filtered.probe_names[0], "cg00000001_chr1_1000")
}

///|
test "minfi_map_to_genome" {
  let grset = @src.minfi_create_example_grset()
  let loc_map = @src.minfi_map_to_genome(grset)
  assert_eq(loc_map.length(), grset.n_probes())
  // Probe cg00000001_chr1_1000 maps to chr1 at position 1000.
  let loc = loc_map.get("cg00000001_chr1_1000")
  assert_true(loc is Some(_))
  let (chr, pos) = loc.unwrap()
  assert_eq(chr, "chr1")
  assert_eq(pos, 1000)
}

///|
test "minfi_map_to_genome_synthetic" {
  // Probe names without embedded location get a synthetic position on chr1.
  let beta = [[0.5]]
  let mvalues = [[0.0]]
  let grset = @src.GenomicRatioSet::new(beta, mvalues, ["s1"], ["cg_only"])
  let loc_map = @src.minfi_map_to_genome(grset)
  let loc = loc_map.get("cg_only")
  assert_true(loc is Some(_))
  let (chr, _pos) = loc.unwrap()
  assert_eq(chr, "chr1")
}
