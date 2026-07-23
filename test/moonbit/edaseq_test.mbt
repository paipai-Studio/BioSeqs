///|
test "edaseq_params" {
  let params = @src.EDASeqParams::new()
  assert_true(params.loess_span > 0.0)
  assert_true(params.loess_span <= 1.0)
  assert_true(params.min_count > 0.0)
}

///|
test "edaseq_gene_anno" {
  let anno = @src.EDASeqGeneAnno::new("gene1", 0.45, 1000.0)
  assert_eq(anno.gene_id, "gene1")
  assert_true(anno.gc_content > 0.0)
  assert_true(anno.length > 0.0)
  assert_true(anno.eda_rpkm_scale() > 0.0)
}

///|
test "edaseq_gene_anno_min_length" {
  // Length below 200 should have effective length boosted to 200
  let anno = @src.EDASeqGeneAnno::new("gene1", 0.5, 100.0)
  assert_true(anno.effective_length >= 200.0)
}

///|
test "edaseq_dataset" {
  let gene_ids = ["gene1", "gene2"]
  let sample_ids = ["sample1", "sample2", "sample3"]
  let counts : Array[Array[Double]] = [
    [10.0, 20.0, 30.0],
    [50.0, 60.0, 70.0]
  ]
  let annotations = [
    @src.EDASeqGeneAnno::new("gene1", 0.45, 1000.0),
    @src.EDASeqGeneAnno::new("gene2", 0.55, 2000.0)
  ]
  let data = @src.EDASeqDataSet::new(gene_ids, sample_ids, counts, annotations)
  assert_eq(data.eda_n_genes(), 2)
  assert_eq(data.eda_n_samples(), 3)
}

///|
test "edaseq_gene_counts" {
  let gene_ids = ["gene1"]
  let sample_ids = ["s1", "s2"]
  let counts : Array[Array[Double]] = [[10.0, 20.0]]
  let annotations = [@src.EDASeqGeneAnno::new("gene1", 0.5, 500.0)]
  let data = @src.EDASeqDataSet::new(gene_ids, sample_ids, counts, annotations)
  let gene_counts = data.eda_gene_counts(0)
  assert_eq(gene_counts.length(), 2)
  assert_eq(gene_counts[0], 10.0)
}

///|
test "edaseq_sample_counts" {
  let gene_ids = ["gene1", "gene2"]
  let sample_ids = ["s1"]
  let counts : Array[Array[Double]] = [[10.0], [20.0]]
  let annotations = [
    @src.EDASeqGeneAnno::new("gene1", 0.5, 500.0),
    @src.EDASeqGeneAnno::new("gene2", 0.6, 1000.0)
  ]
  let data = @src.EDASeqDataSet::new(gene_ids, sample_ids, counts, annotations)
  let sample_counts = data.eda_sample_counts(0)
  assert_eq(sample_counts.length(), 2)
  assert_eq(sample_counts[0], 10.0)
  assert_eq(sample_counts[1], 20.0)
}

///|
test "edaseq_within_lane" {
  let gene_ids = ["gene1", "gene2", "gene3"]
  let sample_ids = ["s1", "s2"]
  let counts : Array[Array[Double]] = [
    [100.0, 200.0],
    [300.0, 400.0],
    [500.0, 600.0]
  ]
  let annotations = [
    @src.EDASeqGeneAnno::new("gene1", 0.4, 500.0),
    @src.EDASeqGeneAnno::new("gene2", 0.5, 1000.0),
    @src.EDASeqGeneAnno::new("gene3", 0.6, 2000.0)
  ]
  let data = @src.EDASeqDataSet::new(gene_ids, sample_ids, counts, annotations)
  let params = @src.EDASeqParams::new()
  let within_result = @src.eda_normalize_within_lane(data, params)
  assert_eq(within_result.norm_counts.length(), 3)
}

///|
test "edaseq_between_lane" {
  let gene_ids = ["gene1", "gene2"]
  let sample_ids = ["s1", "s2", "s3"]
  let counts : Array[Array[Double]] = [
    [100.0, 200.0, 300.0],
    [500.0, 600.0, 700.0]
  ]
  let annotations = [
    @src.EDASeqGeneAnno::new("gene1", 0.5, 500.0),
    @src.EDASeqGeneAnno::new("gene2", 0.55, 1000.0)
  ]
  let data = @src.EDASeqDataSet::new(gene_ids, sample_ids, counts, annotations)
  let params = @src.EDASeqParams::new()
  let between_result = @src.eda_normalize_between_lane(data, params)
  assert_true(between_result.norm_counts.length() > 0)
}

///|
test "edaseq_full_run" {
  let gene_ids = ["gene1", "gene2", "gene3"]
  let sample_ids = ["s1", "s2"]
  let counts : Array[Array[Double]] = [
    [100.0, 200.0],
    [300.0, 400.0],
    [500.0, 600.0]
  ]
  let annotations = [
    @src.EDASeqGeneAnno::new("gene1", 0.4, 500.0),
    @src.EDASeqGeneAnno::new("gene2", 0.5, 1000.0),
    @src.EDASeqGeneAnno::new("gene3", 0.6, 2000.0)
  ]
  let data = @src.EDASeqDataSet::new(gene_ids, sample_ids, counts, annotations)
  let params = @src.EDASeqParams::new()
  let result = @src.eda_run(data, params)
  assert_true(result.eda_n_genes() > 0)
  assert_true(result.eda_n_samples() > 0)
}

///|
test "edaseq_result_accessors" {
  let gene_ids = ["gene1", "gene2"]
  let sample_ids = ["s1", "s2"]
  let counts : Array[Array[Double]] = [
    [100.0, 200.0],
    [300.0, 400.0]
  ]
  let annotations = [
    @src.EDASeqGeneAnno::new("gene1", 0.5, 500.0),
    @src.EDASeqGeneAnno::new("gene2", 0.55, 1000.0)
  ]
  let data = @src.EDASeqDataSet::new(gene_ids, sample_ids, counts, annotations)
  let params = @src.EDASeqParams::new()
  let result = @src.eda_run(data, params)
  assert_true(result.eda_get_count(0, 0) > 0.0)
  assert_eq(result.eda_get_gene(0).length(), 2)
  assert_eq(result.eda_get_sample(0).length(), 2)
}
