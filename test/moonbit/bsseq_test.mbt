///|
/// Tests for the bsseq module.
/// Mirrors the Bioconductor bsseq package API.

///|
test "bsseq_bsrecord_new" {
  let rec = @src.BSRecord::new("chr1", 1000, "+", 8, 2)
  assert_eq(rec.chr, "chr1")
  assert_eq(rec.pos, 1000)
  assert_eq(rec.strand, "+")
  assert_eq(rec.meth, 8)
  assert_eq(rec.unmeth, 2)
  assert_eq(rec.context, "CpG")
}

///|
test "bsseq_bsrecord_new_with_context" {
  let rec = @src.BSRecord::new_with_context("chr2", 500, "-", 3, 7, "CHG")
  assert_eq(rec.chr, "chr2")
  assert_eq(rec.pos, 500)
  assert_eq(rec.strand, "-")
  assert_eq(rec.meth, 3)
  assert_eq(rec.unmeth, 7)
  assert_eq(rec.context, "CHG")
}

///|
test "bsseq_bsrecord_coverage" {
  let rec = @src.BSRecord::new("chr1", 100, "+", 6, 4)
  assert_eq(@src.BSRecord::coverage(rec), 10)
}

///|
test "bsseq_bsrecord_rate" {
  let rec = @src.BSRecord::new("chr1", 100, "+", 8, 2)
  let rate = @src.BSRecord::rate(rec)
  assert_true(rate > 0.79 && rate < 0.81)
}

///|
test "bsseq_bsrecord_to_string" {
  let rec = @src.BSRecord::new("chr1", 100, "+", 5, 5)
  let s = @src.BSRecord::to_string(rec)
  assert_true(s.contains("BSRecord"))
  assert_true(s.contains("chr1"))
  assert_true(s.contains("100"))
}

///|
test "bsseq_get_methylation_rate" {
  // meth=8, unmeth=2 -> rate = 0.8
  let rate = @src.bsseq_get_methylation_rate(8, 2)
  assert_true(rate > 0.79 && rate < 0.81)
}

///|
test "bsseq_get_methylation_rate_zero_total" {
  // meth=0, unmeth=0 -> rate = 0
  assert_eq(@src.bsseq_get_methylation_rate(0, 0), 0.0)
}

///|
test "bsseq_create_example_records" {
  let records = @src.bsseq_create_example_records()
  assert_eq(records.length(), 10)
  assert_eq(records[0].chr, "chr1")
  assert_eq(records[0].pos, 1000)
  assert_eq(records[0].meth, 8)
  assert_eq(records[0].unmeth, 2)
  assert_eq(records[0].context, "CpG")
}

///|
test "bsseq_create_example_bsdata" {
  let bsdata = @src.bsseq_create_example_bsdata()
  assert_eq(bsdata.records.length(), 10)
  assert_eq(bsdata.sample_names.length(), 1)
  assert_eq(bsdata.sample_names[0], "sample1")
  assert_eq(bsdata.coverage.length(), 10)
  assert_eq(bsdata.methylation.length(), 10)
}

///|
test "bsseq_bsdata_n_records_and_n_samples" {
  let bsdata = @src.bsseq_create_example_bsdata()
  assert_eq(@src.BSData::n_records(bsdata), 10)
  assert_eq(@src.BSData::n_samples(bsdata), 1)
}

///|
test "bsseq_get_coverage" {
  let bsdata = @src.bsseq_create_example_bsdata()
  let cov = @src.bsseq_get_coverage(bsdata)
  assert_eq(cov.length(), 10)
  // First record: meth=8, unmeth=2 -> cov=10
  assert_eq(cov[0], 10)
  // Second record: meth=9, unmeth=1 -> cov=10
  assert_eq(cov[1], 10)
}

///|
test "bsseq_get_methylation" {
  let bsdata = @src.bsseq_create_example_bsdata()
  let meth = @src.bsseq_get_methylation(bsdata)
  assert_eq(meth.length(), 10)
  // First record: meth=8, unmeth=2 -> rate=0.8
  assert_true(meth[0] > 0.79 && meth[0] < 0.81)
}

///|
test "bsseq_filter_by_coverage" {
  let bsdata = @src.bsseq_create_example_bsdata()
  // All records have coverage=10, so filtering 5..15 keeps all 10.
  let filtered = @src.bsseq_filter_by_coverage(bsdata, 5, 15)
  assert_eq(filtered.records.length(), 10)
}

///|
test "bsseq_filter_by_coverage_strict" {
  let bsdata = @src.bsseq_create_example_bsdata()
  // Coverage is always 10; filtering to >10 yields no records.
  let filtered = @src.bsseq_filter_by_coverage(bsdata, 11, 100)
  assert_eq(filtered.records.length(), 0)
}

///|
test "bsseq_filter_by_context" {
  // Build a BSData with mixed contexts.
  let mixed : Array[@src.BSRecord] = [
    @src.BSRecord::new_with_context("chr1", 1, "+", 5, 5, "CpG"),
    @src.BSRecord::new_with_context("chr1", 2, "+", 5, 5, "CHG"),
    @src.BSRecord::new_with_context("chr1", 3, "+", 5, 5, "CHH"),
    @src.BSRecord::new_with_context("chr1", 4, "+", 5, 5, "CpG"),
  ]
  let bsdata = @src.BSData::new(mixed, ["s1"])
  let cpg_only = @src.bsseq_filter_by_context(bsdata, "CpG")
  assert_eq(cpg_only.records.length(), 2)
  let chg_only = @src.bsseq_filter_by_context(bsdata, "CHG")
  assert_eq(chg_only.records.length(), 1)
  let chh_only = @src.bsseq_filter_by_context(bsdata, "CHH")
  assert_eq(chh_only.records.length(), 1)
}

///|
test "bsseq_smooth_methylation_basic" {
  // 5 positions, all coverage 10, meth values 0..4 -> rate = meth/10.
  let meth = [0, 1, 2, 3, 4]
  let cov = [10, 10, 10, 10, 10]
  let smoothed = @src.bsseq_smooth_methylation(meth, cov, 1)
  assert_eq(smoothed.length(), 5)
  // At position 0 (window=1): mean of rates at idx 0,1 = (0.0 + 0.1) / 2.
  assert_true(smoothed[0] >= 0.0 && smoothed[0] <= 0.1)
  // At position 2: mean of rates at idx 1,2,3 = (0.1 + 0.2 + 0.3) / 3 = 0.2.
  assert_true(smoothed[2] > 0.19 && smoothed[2] < 0.21)
}

///|
test "bsseq_smooth_methylation_empty" {
  let smoothed = @src.bsseq_smooth_methylation([], [], 2)
  assert_eq(smoothed.length(), 0)
}

///|
test "bsseq_bsmooth_single_chr" {
  let bsdata = @src.bsseq_create_example_bsdata()
  let result = @src.bsseq_bsmooth(bsdata, 2)
  assert_eq(result.pos.length(), 10)
  assert_eq(result.smoothed_meth.length(), 10)
  assert_eq(result.raw_meth.length(), 10)
  assert_eq(result.coverage.length(), 10)
  assert_eq(result.window_size, 2)
  assert_eq(result.chr, "chr1")
}

///|
test "bsseq_compute_tstat_basic" {
  // Two groups with mean methylation differences.
  let group1 = [0.9, 0.8, 0.85]
  let group2 = [0.1, 0.2, 0.15]
  let tstat = @src.bsseq_compute_tstat(group1, group2)
  assert_eq(tstat.t_stat.length(), 3)
  assert_eq(tstat.pvalue.length(), 3)
  assert_eq(tstat.fdr.length(), 3)
  assert_eq(tstat.mean_diff.length(), 3)
  // All differences should be positive and ~0.7
  assert_true(tstat.mean_diff[0] > 0.69 && tstat.mean_diff[0] < 0.81)
}

///|
test "bsseq_compute_tstat_empty" {
  let tstat = @src.bsseq_compute_tstat([], [])
  assert_eq(tstat.t_stat.length(), 0)
  assert_eq(tstat.pvalue.length(), 0)
}

///|
test "bsseq_compute_tstat_length_mismatch" {
  let tstat = @src.bsseq_compute_tstat([0.5, 0.6], [0.3])
  assert_eq(tstat.t_stat.length(), 0)
}

///|
test "bsseq_find_dmr_basic" {
  // Build a TStatResult with large |t| values to trigger DMR detection.
  let group1 = [0.9, 0.9, 0.1, 0.1, 0.9]
  let group2 = [0.1, 0.1, 0.9, 0.9, 0.1]
  let tstat = @src.bsseq_compute_tstat(group1, group2)
  let positions = [100, 200, 300, 400, 500]
  let dmr = @src.bsseq_find_dmr(tstat, positions, "chr1", 1.0, 200)
  assert_true(dmr.n_dmr >= 0)
  assert_eq(dmr.threshold, 1.0)
}

///|
test "bsseq_dmr_test_basic" {
  // group1: rates all 0.9, group2: rates all 0.1
  let meth1 = [9, 9, 9]
  let cov1 = [10, 10, 10]
  let meth2 = [1, 1, 1]
  let cov2 = [10, 10, 10]
  let (t, p) = @src.bsseq_dmr_test(meth1, cov1, meth2, cov2)
  // t should be large positive (group1 > group2)
  assert_true(t > 0.0)
  // p-value should be small for such a clear difference
  assert_true(p >= 0.0 && p <= 1.0)
}

///|
test "bsseq_dmr_test_empty" {
  let (t, p) = @src.bsseq_dmr_test([], [], [], [])
  assert_eq(t, 0.0)
  assert_eq(p, 1.0)
}

///|
test "bsseq_find_hyper_dmr" {
  // Construct a BsseqDMRResult with mixed directions.
  let regions : Array[@src.DMRegion] = [
    @src.DMRegion::new("chr1", 100, 200, 2, 0.5, 0.01, 0.02, "hyper"),
    @src.DMRegion::new("chr1", 300, 400, 2, -0.5, 0.01, 0.02, "hypo"),
    @src.DMRegion::new("chr1", 500, 600, 2, 0.3, 0.01, 0.02, "hyper"),
  ]
  let result = @src.BsseqDMRResult::new(regions, 0.05)
  let hyper = @src.bsseq_find_hyper_dmr(result)
  assert_eq(hyper.length(), 2)
}

///|
test "bsseq_find_hypo_dmr" {
  let regions : Array[@src.DMRegion] = [
    @src.DMRegion::new("chr1", 100, 200, 2, 0.5, 0.01, 0.02, "hyper"),
    @src.DMRegion::new("chr1", 300, 400, 2, -0.5, 0.01, 0.02, "hypo"),
    @src.DMRegion::new("chr1", 500, 600, 2, -0.3, 0.01, 0.02, "hypo"),
  ]
  let result = @src.BsseqDMRResult::new(regions, 0.05)
  let hypo = @src.bsseq_find_hypo_dmr(result)
  assert_eq(hypo.length(), 2)
}

///|
test "bsseq_compute_methylation_diff" {
  assert_eq(@src.bsseq_compute_methylation_diff(0.8, 0.3), 0.5)
  assert_true(@src.bsseq_compute_methylation_diff(0.1, 0.4) > -0.31 && @src.bsseq_compute_methylation_diff(0.1, 0.4) < -0.29)
  assert_eq(@src.bsseq_compute_methylation_diff(0.5, 0.5), 0.0)
}

///|
test "bsseq_merge_cpgs_basic" {
  let records : Array[@src.BSRecord] = [
    @src.BSRecord::new("chr1", 100, "+", 5, 5),
    @src.BSRecord::new("chr1", 105, "+", 3, 7),
    @src.BSRecord::new("chr1", 500, "+", 4, 6),
  ]
  let merged = @src.bsseq_merge_cpgs(records, 10)
  // 100 and 105 are within 10bp -> merged; 500 is separate.
  assert_eq(merged.length(), 2)
  // Merged record at pos 100: meth=8, unmeth=12.
  assert_eq(merged[0].meth, 8)
  assert_eq(merged[0].unmeth, 12)
  assert_eq(merged[0].pos, 100)
}

///|
test "bsseq_merge_cpgs_empty" {
  let merged = @src.bsseq_merge_cpgs([], 10)
  assert_eq(merged.length(), 0)
}

///|
test "bsseq_region_methylation" {
  let bsdata = @src.bsseq_create_example_bsdata()
  // Region covering positions 1000..1500 includes records at 1000,1100,1200,1300,1400,1500.
  let rate = @src.bsseq_region_methylation(bsdata, "chr1", 1000, 1500)
  // Coverage-weighted mean of (8/10, 9/10, 5/10, 3/10, 7/10, 6/10) -> total meth = 38, total cov = 60.
  assert_true(rate > 0.62 && rate < 0.64)
}

///|
test "bsseq_region_methylation_no_match" {
  let bsdata = @src.bsseq_create_example_bsdata()
  // Empty region -> 0.
  let rate = @src.bsseq_region_methylation(bsdata, "chrX", 1, 10)
  assert_eq(rate, 0.0)
}

///|
test "bsseq_normalize_none" {
  let bsdata = @src.bsseq_create_example_bsdata()
  let normalized = @src.bsseq_normalize(bsdata, "none")
  assert_eq(normalized.records.length(), bsdata.records.length())
  // "none" should preserve meth counts.
  assert_eq(normalized.records[0].meth, bsdata.records[0].meth)
  assert_eq(normalized.records[0].unmeth, bsdata.records[0].unmeth)
}

///|
test "bsseq_normalize_global" {
  let bsdata = @src.bsseq_create_example_bsdata()
  let normalized = @src.bsseq_normalize(bsdata, "global")
  assert_eq(normalized.records.length(), bsdata.records.length())
  // Coverage is preserved by global normalization.
  let mut i = 0
  while i < normalized.records.length() {
    let r = normalized.records[i]
    let orig = bsdata.records[i]
    assert_eq(r.meth + r.unmeth, orig.meth + orig.unmeth)
    i = i + 1
  }
}

///|
test "bsseq_combine_samples" {
  let cov = @src.bsseq_create_example_cov_matrix()
  assert_eq(cov.pos.length(), 10)
  assert_eq(cov.chr.length(), 10)
  // Two samples were combined.
  assert_eq(cov.cov.length(), 10)
  assert_eq(cov.cov[0].length(), 2)
  assert_eq(cov.meth[0].length(), 2)
}

///|
test "bsseq_summary" {
  let bsdata = @src.bsseq_create_example_bsdata()
  let s = @src.bsseq_summary(bsdata)
  assert_true(s.contains("BS-seq Data Summary"))
  assert_true(s.contains("Records: 10"))
  assert_true(s.contains("Samples: 1"))
}

///|
test "bsseq_bsdata_to_string" {
  let bsdata = @src.bsseq_create_example_bsdata()
  let s = @src.BSData::to_string(bsdata)
  assert_true(s.contains("BSData"))
  assert_true(s.contains("Records: 10"))
}

///|
test "bsseq_mean" {
  assert_eq(@src.bsseq_mean([1.0, 2.0, 3.0]), 2.0)
  assert_eq(@src.bsseq_mean([]), 0.0)
  assert_eq(@src.bsseq_mean([5.0]), 5.0)
}

///|
test "bsseq_median" {
  assert_eq(@src.bsseq_median([1.0, 2.0, 3.0]), 2.0)
  assert_eq(@src.bsseq_median([1.0, 2.0, 3.0, 4.0]), 2.5)
  assert_eq(@src.bsseq_median([]), 0.0)
}

///|
test "bsseq_variance" {
  // Variance of [1,2,3,4,5] is 2.5 (sample variance, n-1 denominator)
  let v = @src.bsseq_variance([1.0, 2.0, 3.0, 4.0, 5.0])
  assert_true(v > 2.49 && v < 2.51)
}

///|
test "bsseq_variance_singleton" {
  assert_eq(@src.bsseq_variance([5.0]), 0.0)
  assert_eq(@src.bsseq_variance([]), 0.0)
}

///|
test "bsseq_t_test_basic" {
  // Two clearly different groups.
  let x = [10.0, 12.0, 11.0, 13.0]
  let y = [1.0, 2.0, 3.0, 4.0]
  let (t, df) = @src.bsseq_t_test(x, y)
  // t should be large positive since x >> y.
  assert_true(t > 0.0)
  assert_true(df > 0.0)
}

///|
test "bsseq_t_test_short_input" {
  // Length < 2 returns (0, 0).
  let (t, df) = @src.bsseq_t_test([1.0], [2.0])
  assert_eq(t, 0.0)
  assert_eq(df, 0.0)
}

///|
test "bsseq_bh_fdr_basic" {
  let pvalues = [0.01, 0.04, 0.03, 0.20]
  let adj = @src.bsseq_bh_fdr(pvalues)
  assert_eq(adj.length(), 4)
  // Adjusted values should all be within [0,1].
  let mut i = 0
  while i < adj.length() {
    assert_true(adj[i] >= 0.0 && adj[i] <= 1.0)
    i = i + 1
  }
}

///|
test "bsseq_bh_fdr_empty" {
  let adj = @src.bsseq_bh_fdr([])
  assert_eq(adj.length(), 0)
}

///|
test "bsseq_parse_bed" {
  // Build a small BED-like methylation file.
  let content = "chrom\tstart\tend\tname\tmeth\tuneth\tstrand\tcontext\nchr1\t100\t101\tCpG_1\t8\t2\t+\tCpG\nchr1\t200\t201\tCpG_2\t4\t6\t+\tCpG\n"
  let bsdata = @src.bsseq_parse_bed(content)
  assert_eq(bsdata.records.length(), 2)
  assert_eq(bsdata.records[0].chr, "chr1")
  assert_eq(bsdata.records[0].pos, 100)
  assert_eq(bsdata.records[0].meth, 8)
  assert_eq(bsdata.records[0].unmeth, 2)
}

///|
test "bsseq_parse_bismark" {
  // Build Bismark methylation extractor output.
  let content = "Bismark methylation extractor\nchr1\t100\t+\t8\t2\tCpG\tCGT\nchr1\t200\t-\t4\t6\tCpG\tCGT\n"
  let bsdata = @src.bsseq_parse_bismark(content)
  assert_eq(bsdata.records.length(), 2)
  assert_eq(bsdata.records[0].chr, "chr1")
  assert_eq(bsdata.records[0].pos, 100)
  assert_eq(bsdata.records[0].meth, 8)
  assert_eq(bsdata.records[0].unmeth, 2)
}

///|
test "bsseq_dmregion_to_string" {
  let r = @src.DMRegion::new("chr1", 100, 200, 3, 0.5, 0.01, 0.02, "hyper")
  let s = @src.DMRegion::to_string(r)
  assert_true(s.contains("DMRegion"))
  assert_true(s.contains("chr1"))
  assert_true(s.contains("hyper"))
}

///|
test "bsseq_bsseq_dmr_result_to_string" {
  let regions : Array[@src.DMRegion] = [
    @src.DMRegion::new("chr1", 100, 200, 2, 0.5, 0.01, 0.02, "hyper"),
  ]
  let result = @src.BsseqDMRResult::new(regions, 0.05)
  let s = @src.BsseqDMRResult::to_string(result)
  assert_true(s.contains("BsseqDMRResult"))
  assert_true(s.contains("Total DMRs: 1"))
}

///|
test "bsseq_tstat_result_to_string" {
  let group1 = [0.9, 0.8, 0.85]
  let group2 = [0.1, 0.2, 0.15]
  let tstat = @src.bsseq_compute_tstat(group1, group2)
  let s = @src.TStatResult::to_string(tstat)
  assert_true(s.contains("TStatResult"))
  assert_true(s.contains("CpGs: 3"))
}

///|
test "bsseq_bsmooth_result_to_string" {
  let bsdata = @src.bsseq_create_example_bsdata()
  let result = @src.bsseq_bsmooth(bsdata, 2)
  let s = @src.BSmoothResult::to_string(result)
  assert_true(s.contains("BSmoothResult"))
  assert_true(s.contains("chr1"))
  assert_true(s.contains("Window size: 2"))
}

///|
test "bsseq_cov_matrix_to_string" {
  let cov = @src.bsseq_create_example_cov_matrix()
  let s = @src.CovMatrix::to_string(cov)
  assert_true(s.contains("CovMatrix"))
  assert_true(s.contains("Positions: 10"))
}

///|
test "bsseq_sort_doubles" {
  let sorted = @src.bsseq_sort_doubles([3.0, 1.0, 4.0, 1.5, 2.0])
  assert_eq(sorted.length(), 5)
  assert_eq(sorted[0], 1.0)
  assert_eq(sorted[1], 1.5)
  assert_eq(sorted[2], 2.0)
  assert_eq(sorted[3], 3.0)
  assert_eq(sorted[4], 4.0)
}

///|
test "bsseq_t_test_pval" {
  // t=0 -> p-value = 1.0
  assert_eq(@src.bsseq_t_test_pval(0.0, 10.0), 1.0)
  // df <= 0 -> p-value = 1.0
  assert_eq(@src.bsseq_t_test_pval(2.0, 0.0), 1.0)
}

///|
test "bsseq_bsrecord_to_string_format" {
  let rec = @src.BSRecord::new_with_context("chrX", 250, "-", 3, 4, "CHH")
  let s = @src.BSRecord::to_string(rec)
  assert_true(s.contains("chrX"))
  assert_true(s.contains("CHH"))
  assert_true(s.contains("250"))
}
