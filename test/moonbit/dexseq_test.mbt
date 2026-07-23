///|
/// Tests for DEXSeq module.

test "ExonCount creation" {
  let ec = @src.ExonCount::new("gene1", "exon1", [10, 20, 30])
  assert_eq(ec.gene_id, "gene1")
  assert_eq(ec.exon_id, "exon1")
  assert_eq(ec.counts.length(), 3)
}

test "DEXSeqDataSet creation" {
  let exon_counts : Array[@src.ExonCount] = Array::new()
  exon_counts.push(@src.ExonCount::new("gene1", "exon1", [10, 12, 15]))
  exon_counts.push(@src.ExonCount::new("gene1", "exon2", [5, 6, 7]))
  let sample_groups = ["control", "control", "treatment"]
  let ds = @src.DEXSeqDataSet::new(exon_counts, sample_groups)
  assert_eq(ds.gene_ids.length(), 1)
  assert_eq(ds.gene_ids[0], "gene1")
}

test "dexseq_normalize_counts" {
  let ds = @src.create_example_dexseq_dataset()
  let normalized = @src.dexseq_normalize_counts(ds)
  assert_eq(normalized.exon_counts.length(), ds.exon_counts.length())
}

test "dexseq_test_for_exon_usage" {
  let ds = @src.create_example_dexseq_dataset()
  let results = @src.dexseq_test_for_exon_usage(ds)
  assert_true(results.length() > 0)
}

test "dexseq_filter_results" {
  let ds = @src.create_example_dexseq_dataset()
  let results = @src.dexseq_test_for_exon_usage(ds)
  let filtered = @src.dexseq_filter_results(results, 0.5, 0.1)
  assert_true(filtered.length() <= results.length())
}
