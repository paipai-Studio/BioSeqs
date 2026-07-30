///|
/// Test file for DRIMSeq module.

test "drimseq_transcript_count_creation" {
  let tc = @src.TranscriptCount::new(
    transcript_id="tx1", gene_id="GeneA",
    sample_id="s1", condition="ctrl",
    count=100.0, gene_count=150.0,
  )
  assert_eq(tc.transcript_id, "tx1")
  assert_eq(tc.gene_id, "GeneA")
  assert_eq(tc.count, 100.0)
  assert_true((tc.proportion() - 0.666666).abs() < 0.01)
}

test "drimseq_transcript_count_proportion_zero" {
  let tc = @src.TranscriptCount::new(
    transcript_id="tx1", gene_id="GeneA",
    sample_id="s1", condition="ctrl",
    count=100.0,
  )
  assert_eq(tc.proportion(), 0.0)
}

test "drimseq_filter_counts" {
  let counts : Array[@src.TranscriptCount] = Array::new()
  counts.push(@src.TranscriptCount::new(
    transcript_id="tx1", gene_id="G1",
    sample_id="s1", condition="c1",
    count=50.0, gene_count=100.0,
  ))
  counts.push(@src.TranscriptCount::new(
    transcript_id="tx2", gene_id="G2",
    sample_id="s1", condition="c1",
    count=5.0, gene_count=100.0,
  ))

  let filtered = @src.drimseq_filter_counts(counts)
  assert_eq(filtered.length(), 1)
  assert_eq(filtered[0].transcript_id, "tx1")
}

test "drimseq_aggregate_by_gene" {
  let counts : Array[@src.TranscriptCount] = Array::new()
  counts.push(@src.TranscriptCount::new(
    transcript_id="tx1", gene_id="G1",
    sample_id="s1", condition="c1", count=50.0,
  ))
  counts.push(@src.TranscriptCount::new(
    transcript_id="tx2", gene_id="G1",
    sample_id="s1", condition="c1", count=30.0,
  ))
  counts.push(@src.TranscriptCount::new(
    transcript_id="tx1", gene_id="G1",
    sample_id="s2", condition="c1", count=40.0,
  ))

  let aggregated = @src.drimseq_aggregate_by_gene(counts)
  assert_eq(aggregated.length(), 2) // G1-s1 and G1-s2
  assert_eq(aggregated[0].count, 80.0)
  assert_eq(aggregated[1].count, 40.0)
}

test "drimseq_compute_proportions" {
  let counts : Array[@src.TranscriptCount] = Array::new()
  counts.push(@src.TranscriptCount::new(
    transcript_id="tx1", gene_id="G1",
    sample_id="s1", condition="c1", count=50.0,
  ))
  counts.push(@src.TranscriptCount::new(
    transcript_id="tx2", gene_id="G1",
    sample_id="s1", condition="c1", count=50.0,
  ))

  let gene_counts : Array[@src.TranscriptCount] = Array::new()
  gene_counts.push(@src.TranscriptCount::new(
    transcript_id="__gene__", gene_id="G1",
    sample_id="s1", condition="c1", count=100.0,
  ))

  let props = @src.drimseq_compute_proportions(counts, gene_counts)
  assert_eq(props.length(), 2)
  assert_true((props[0].proportion() - 0.5).abs() < 0.01)
}

test "drimseq_wald_test_same" {
  let a : Array[Double] = Array::new()
  a.push(0.5); a.push(0.5)
  let b : Array[Double] = Array::new()
  b.push(0.5); b.push(0.5)

  let (stat, pval, _) = @src.drimseq_wald_test(a, b)
  assert_true(pval >= 0.99) // Should not be significant
}

test "drimseq_wald_test_different" {
  let a : Array[Double] = Array::new()
  a.push(0.8); a.push(0.2)
  let b : Array[Double] = Array::new()
  b.push(0.2); b.push(0.8)

  let (_, pval, delta) = @src.drimseq_wald_test(a, b)
  assert_true(pval < 0.05) // Should be significant
  assert_true(delta.length() == 2)
}

test "drimseq_result_summary" {
  let counts = @src.drimseq_sample_data()
  let result = @src.drimseq_test_differential(counts)
  let summary = result.summary()
  assert_true(summary.length() > 0)
}

test "drimseq_result_get_significant" {
  let counts = @src.drimseq_sample_data()
  let result = @src.drimseq_test_differential(counts)
  let sig = result.get_significant()
  assert_true(sig.length() >= 0)
}

test "drimseq_result_get_top_genes" {
  let counts = @src.drimseq_sample_data()
  let result = @src.drimseq_test_differential(counts)
  let top = result.get_top_genes(n=2)
  assert_true(top.length() <= 2)
}

test "drimseq_config_creation" {
  let cfg = @src.DRIMSeqConfig::new()
  assert_true(cfg.min_count > 0.0)
  assert_true(cfg.alpha > 0.0)
}

test "drimseq_config_setters" {
  let cfg = @src.DRIMSeqConfig::new()
    .set_min_count(val=5.0)
    .set_alpha(val=0.01)
  assert_eq(cfg.min_count, 5.0)
  assert_eq(cfg.alpha, 0.01)
}

test "drimseq_norm_creation" {
  let n1 = @src.drimseq_norm_none()
  let n2 = @src.drimseq_norm_tmm()
  let n3 = @src.drimseq_norm_sum()
  assert_true(n1 == n1)
  assert_true(n2 == n2)
  assert_true(n3 == n3)
}

test "drimseq_full_pipeline" {
  let counts = @src.drimseq_sample_data()
  let result = @src.drimseq_test_differential(counts)
  assert_true(result.n_genes_tested > 0)
  assert_true(result.n_genes_filtered > 0)
  assert_true(result.converged)
}

test "drimseq_empty_data" {
  let empty : Array[@src.TranscriptCount] = Array::new()
  let result = @src.drimseq_test_differential(empty)
  assert_eq(result.n_genes_tested, 0)
}
