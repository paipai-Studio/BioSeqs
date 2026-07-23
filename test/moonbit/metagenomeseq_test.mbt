///|
/// Tests for metagenomeSeq module.

test "MRexperiment creation" {
  let counts = [[10, 20], [30, 40]]
  let taxa_names = ["Bacteroides", "Firmicutes"]
  let sample_names = ["sample1", "sample2"]
  let obj = @src.MRexperiment::new(counts, taxa_names, sample_names)
  assert_eq(obj.taxa_names.length(), 2)
  assert_eq(obj.sample_names.length(), 2)
}

test "mg_normalize_counts" {
  let obj = @src.create_example_mrexperiment()
  let normalized = @src.mg_normalize_counts(obj)
  assert_eq(normalized.counts.length(), obj.counts.length())
}

test "mg_calculate_zero_inflation" {
  let obj = @src.create_example_mrexperiment()
  let zero_probs = @src.mg_calculate_zero_inflation(obj)
  assert_eq(zero_probs.length(), obj.taxa_names.length())
}

test "mg_test_zero_inflated" {
  let obj = @src.create_example_mrexperiment()
  let group = ["control", "control", "control", "treatment", "treatment", "treatment"]
  let results = @src.mg_test_zero_inflated(obj, group)
  assert_eq(results.length(), obj.taxa_names.length())
}

test "MRSampleData creation" {
  let sd = @src.MRSampleData::new("control", 100000)
  assert_eq(sd.group, "control")
  assert_eq(sd.depth, 100000)
}
