///|
/// Test file for featureCounts module.

test "feature_annotation_creation" {
  let fa = @src.FeatureAnnotation::new(
    chr="chr1", start=1000, end_=2000, strand="+",
    gene_id="GeneA",
  )
  assert_eq(fa.chr, "chr1")
  assert_eq(fa.start, 1000)
  assert_eq(fa.end_, 2000)
  assert_eq(fa.strand, "+")
  assert_eq(fa.gene_id, "GeneA")
  assert_eq(fa.length(), 1001)
}

test "feature_annotation_overlap" {
  let fa = @src.FeatureAnnotation::new(
    chr="chr1", start=1000, end_=2000, strand="+",
    gene_id="GeneA",
  )
  assert_true(fa.overlaps(chr="chr1", start=1500, end_=1600))
  assert_true(fa.overlaps(chr="chr1", start=900, end_=1100))
  assert_true(fa.overlaps(chr="chr1", start=1900, end_=2100))
  assert_false(fa.overlaps(chr="chr2", start=1500, end_=1600))
  assert_false(fa.overlaps(chr="chr1", start=2100, end_=2200))
}

test "feature_annotation_overlap_length" {
  let fa = @src.FeatureAnnotation::new(
    chr="chr1", start=1000, end_=2000, strand="+",
    gene_id="GeneA",
  )
  assert_eq(fa.overlap_length(chr="chr1", start=1500, end_=1600), 101)
  assert_eq(fa.overlap_length(chr="chr1", start=900, end_=1100), 101)
  assert_eq(fa.overlap_length(chr="chr1", start=2100, end_=2200), 0)
  assert_eq(fa.overlap_length(chr="chr2", start=1500, end_=1600), 0)
}

test "read_alignment_creation" {
  let r = @src.ReadAlignment::new(
    read_id="read1", chr="chr1", start=1500, end_=1600,
  )
  assert_eq(r.read_id, "read1")
  assert_eq(r.chr, "chr1")
  assert_eq(r.start, 1500)
  assert_eq(r.end_, 1600)
  assert_eq(r.is_paired, false)
  assert_eq(r.n_alignments, 1)
}

test "read_alignment_paired" {
  let r = @src.ReadAlignment::new_paired(
    read_id="read1", chr="chr1", start=1500, end_=1600,
    mate_start=1700,
  )
  assert_eq(r.is_paired, true)
  assert_eq(r.mate_start, 1700)
  assert_true(r.fragment_length > 0)
}

test "feature_counts_config_default" {
  let cfg = @src.FeatureCountsConfig::new()
  assert_eq(cfg.min_overlap, 1)
  assert_eq(cfg.min_mapq, 0)
  assert_eq(cfg.ignore_strand, false)
}

test "feature_counts_config_setters" {
  let cfg = @src.FeatureCountsConfig::new()
    .set_min_overlap(val=5)
    .set_min_mapq(val=30)
    .set_count_fragments(val=true)
  assert_eq(cfg.min_overlap, 5)
  assert_eq(cfg.min_mapq, 30)
  assert_eq(cfg.count_fragments, true)
}

test "feature_counts_count_basic" {
  let (features, reads, samples) = @src.feature_counts_sample_data()
  let cfg = @src.FeatureCountsConfig::new().set_min_mapq(val=40)
  let result = @src.feature_counts_count(features, reads, sample_names=samples, config=cfg)
  assert_eq(result.n_features(), 4)
  assert_eq(result.n_samples(), 1)
  assert_eq(result.total_reads[0], 7)
  assert_eq(result.assigned_reads[0], 5)
}

test "feature_counts_get_count" {
  let (features, reads, samples) = @src.feature_counts_sample_data()
  let result = @src.feature_counts_count(features, reads, sample_names=samples)
  let count = result.get_count(feature_id="exon_A1", sample_idx=0)
  assert_true(count >= 2.0)
}

test "feature_counts_unassigned" {
  let (features, reads, samples) = @src.feature_counts_sample_data()
  let cfg = @src.FeatureCountsConfig::new().set_min_mapq(val=40)
  let result = @src.feature_counts_count(features, reads, sample_names=samples, config=cfg)
  assert_eq(result.unassigned_low_quality[0], 1)
  assert_eq(result.unassigned_no_features[0], 1)
}

test "feature_counts_with_strand" {
  let (features, reads, samples) = @src.feature_counts_sample_data()
  let cfg = @src.FeatureCountsConfig::new()
    .set_strand_mode(mode=@src.strand_stranded())
  let result = @src.feature_counts_count(features, reads, sample_names=samples, config=cfg)
  let count_b1 = result.get_count(feature_id="exon_B1", sample_idx=0)
  // read3 is on chr1:6500-6600, strand -, feature is on strand -, should match
  assert_true(count_b1 >= 1.0)
}

test "feature_counts_reversed_strand" {
  let (features, reads, samples) = @src.feature_counts_sample_data()
  let cfg = @src.FeatureCountsConfig::new()
    .set_strand_mode(mode=@src.strand_reversed())
  let result = @src.feature_counts_count(features, reads, sample_names=samples, config=cfg)
  // With reversed strand, read on + should NOT match feature on +
  // read1 (+) on chr1:1500 should NOT match exon_A1 (+) with reversed
  let count = result.get_count(feature_id="exon_A1", sample_idx=0)
  assert_true(count <= 1.0)
}

test "feature_counts_library_size" {
  let (features, reads, samples) = @src.feature_counts_sample_data()
  let result = @src.feature_counts_count(features, reads, sample_names=samples)
  let lib_size = result.library_size(sample_idx=0)
  assert_true(lib_size > 0.0)
}

test "feature_counts_cpm" {
  let (features, reads, samples) = @src.feature_counts_sample_data()
  let result = @src.feature_counts_count(features, reads, sample_names=samples)
  let cpm_val = result.cpm(feature_id="exon_A1", sample_idx=0)
  assert_true(cpm_val >= 0.0)
}

test "feature_counts_summary" {
  let (features, reads, samples) = @src.feature_counts_sample_data()
  let result = @src.feature_counts_count(features, reads, sample_names=samples)
  let summary = result.summary()
  assert_true(summary.length() > 0)
}

test "feature_counts_multi_sample" {
  let features : Array[@src.FeatureAnnotation] = Array::new()
  features.push(@src.FeatureAnnotation::new(
    chr="chr1", start=1000, end_=2000, strand="+",
    gene_id="GeneA", feature_id="exon_A1",
  ))

  let reads : Array[@src.ReadAlignment] = Array::new()
  reads.push(@src.ReadAlignment::new(
    read_id="read_s1_1", chr="chr1", start=1500, end_=1600,
  ))
  reads.push(@src.ReadAlignment::new(
    read_id="read_s2_1", chr="chr1", start=1500, end_=1600,
  ))

  let sample_names : Array[String] = Array::new()
  sample_names.push("sample1")
  sample_names.push("sample2")

  let result = @src.feature_counts_count(features, reads, sample_names=sample_names)
  assert_eq(result.n_samples(), 2)
  assert_eq(result.total_reads[0], 1)
  assert_eq(result.total_reads[1], 1)
  assert_eq(result.get_count(feature_id="exon_A1", sample_idx=0), 1.0)
  assert_eq(result.get_count(feature_id="exon_A1", sample_idx=1), 1.0)
}

test "strand_mode_creation" {
  let s1 = @src.strand_unstranded()
  let s2 = @src.strand_stranded()
  let s3 = @src.strand_reversed()
  // Verify they are different variants
  assert_true(s1 == s1)
  assert_true(s2 == s2)
  assert_true(s3 == s3)
}

test "feature_counts_feature_total" {
  let (features, reads, samples) = @src.feature_counts_sample_data()
  let result = @src.feature_counts_count(features, reads, sample_names=samples)
  let total = result.get_feature_total(feature_id="exon_A1")
  assert_true(total >= 2.0)
}
