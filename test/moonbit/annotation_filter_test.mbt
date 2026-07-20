///|
test "af_create_example_annotations" {
  let annotations = @src.af_create_example_annotations()
  assert_eq(annotations.length(), 12)
}

///|
test "af_filter_chromosome" {
  let annotations = @src.af_create_example_annotations()
  let filtered = @src.af_filter_chromosome(annotations, "chr17")
  assert_eq(filtered.length(), 2)
}

///|
test "af_filter_biotype" {
  let annotations = @src.af_create_example_annotations()
  let filtered = @src.af_filter_biotype(annotations, "miRNA")
  assert_eq(filtered.length(), 2)
}

///|
test "af_filter_strand" {
  let annotations = @src.af_create_example_annotations()
  let filtered = @src.af_filter_strand(annotations, "-")
  assert_eq(filtered.length(), 3)
}

///|
test "af_filter_description_contains" {
  let annotations = @src.af_create_example_annotations()
  let filtered = @src.af_filter_description_contains(annotations, "cancer")
  assert_eq(filtered.length(), 2)
}

///|
test "af_filter_symbol_pattern" {
  let annotations = @src.af_create_example_annotations()
  let filtered = @src.af_filter_symbol_pattern(annotations, "BRCA")
  assert_eq(filtered.length(), 2)
}

///|
test "af_filter_biotype_in" {
  let annotations = @src.af_create_example_annotations()
  let filtered = @src.af_filter_biotype_in(annotations, ["miRNA", "lincRNA"])
  assert_eq(filtered.length(), 3)
}

///|
test "af_filter_region" {
  let annotations = @src.af_create_example_annotations()
  let filtered = @src.af_filter_region(annotations, "chr17", 7000000, 8000000)
  assert_eq(filtered.length(), 1)
}

///|
test "af_count_by_biotype" {
  let annotations = @src.af_create_example_annotations()
  let counts = @src.af_count_by_biotype(annotations)
  assert_true(counts.get("protein_coding").is_some())
}

///|
test "af_filter_summary" {
  let annotations = @src.af_create_example_annotations()
  let summary = @src.af_filter_summary(annotations)
  assert_true(summary.contains("Annotation"))
}
