///|
/// Tests for Singscore gene set scoring module.

test "singscore_create_example" {
  let (sample, spec) = @src.singscore_create_example()
  assert_true(sample.sample_id == "sample_1")
  assert_true(sample.gene_names.length() == 20)
  assert_true(spec.up_genes.length() > 0)
  assert_true(spec.down_genes.length() > 0)
}

test "singscore_score_basic" {
  let (sample, spec) = @src.singscore_create_example()
  let result = @src.singscore_score(sample, spec, "test_sample")
  assert_true(result.sample_id == "test_sample")
  assert_true(result.n_up_found > 0)
  assert_true(result.score >= -0.5 && result.score <= 0.5)
}

test "singscore_score_boundary" {
  // Score should always be in [-0.5, 0.5]
  let (sample, spec) = @src.singscore_create_example()
  let result = @src.singscore_score(sample, spec, "boundary_test")
  assert_true(result.score >= -0.5 - 1.0e-6)
  assert_true(result.score <= 0.5 + 1.0e-6)
}

test "singscore_empty_gene_set" {
  let genes = ["A", "B", "C"]
  let expression = [1.0, 2.0, 3.0]
  let sample = @src.SampleExpression::new("empty_test", genes, expression)
  let spec = @src.GeneSetSpec::new("empty_set", [], down_genes = [])
  let result = @src.singscore_score(sample, spec, "empty")
  assert_true(result.score == 0.0)
  assert_true(result.n_up_found == 0)
  assert_true(result.n_down_found == 0)
}

test "singscore_single_gene" {
  let genes = ["TP53"]
  let expression = [10.0]
  let sample = @src.SampleExpression::new("single_test", genes, expression)
  let spec = @src.GeneSetSpec::new("single_gene_set", ["TP53"])
  let result = @src.singscore_score(sample, spec, "single")
  assert_true(result.n_up_found == 1)
  assert_true(result.score >= -0.5 && result.score <= 0.5)
}

test "singscore_multiple_samples" {
  let (sample1, spec) = @src.singscore_create_example()
  let genes = sample1.gene_names
  let expression2 : Array[Double] = Array::new()
  for _g in genes {
    expression2.push(1.0)
  }
  let sample2 = @src.SampleExpression::new("sample_2", genes, expression2)
  let results = @src.singscore_score_multiple(sample2, [spec])
  assert_true(results.length() == 1)
  assert_true(results[0].score >= -0.5 && results[0].score <= 0.5)
}

test "singscore_missing_genes" {
  let genes = ["A", "B"]
  let expression = [1.0, 2.0]
  let sample = @src.SampleExpression::new("missing_test", genes, expression)
  let spec = @src.GeneSetSpec::new("missing_set", ["C", "D", "E"])
  let result = @src.singscore_score(sample, spec, "missing")
  assert_true(result.n_up_found == 0)
  assert_true(result.score == 0.0)
}

test "singscore_spec_creation" {
  let spec = @src.GeneSetSpec::new("test_set", ["G1", "G2"], down_genes = ["G3"])
  assert_true(spec.name == "test_set")
  assert_true(spec.up_genes.length() == 2)
  assert_true(spec.down_genes.length() == 1)
}
