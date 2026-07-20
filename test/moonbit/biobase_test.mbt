///|
test "biobase_expression_set_basic" {
  let assay = [[1.0, 2.0, 3.0], [4.0, 5.0, 6.0], [7.0, 8.0, 9.0]]
  let pheno : Map[String, Array[String]] = Map(
    [("sample_name", ["s1", "s2", "s3"])],
    capacity=1,
  )
  let feature : Map[String, Array[String]] = Map(
    [("feature_name", ["g1", "g2", "g3"])],
    capacity=1,
  )

  let eset = @src.ExpressionSet::new(assay, pheno, feature, "test", Map([]))

  assert_eq(@src.expression_set_nrow(eset), 3)
  assert_eq(@src.expression_set_ncol(eset), 3)
  assert_eq(@src.expression_set_get_annotation(eset), "test")
}

///|
test "biobase_expression_set_sample_names" {
  let assay = [[1.0, 2.0]]
  let pheno : Map[String, Array[String]] = Map(
    [("sample_name", ["sampleA", "sampleB"])],
    capacity=1,
  )
  let feature : Map[String, Array[String]] = Map([])

  let eset = @src.ExpressionSet::new(assay, pheno, feature, "", Map([]))
  let names = @src.expression_set_sample_names(eset)

  assert_eq(names[0], "sampleA")
  assert_eq(names[1], "sampleB")
}

///|
test "biobase_expression_set_subset" {
  let assay = [[1.0, 2.0, 3.0], [4.0, 5.0, 6.0], [7.0, 8.0, 9.0]]
  let pheno : Map[String, Array[String]] = Map(
    [("sample_name", ["s1", "s2", "s3"])],
    capacity=1,
  )
  let feature : Map[String, Array[String]] = Map(
    [("feature_name", ["g1", "g2", "g3"])],
    capacity=1,
  )

  let eset = @src.ExpressionSet::new(assay, pheno, feature, "", Map([]))
  let subset = @src.expression_set_subset(eset, [0, 2], [0, 2])

  assert_eq(@src.expression_set_nrow(subset), 2)
  assert_eq(@src.expression_set_ncol(subset), 2)
}

///|
test "biobase_expression_set_filter_by_row" {
  let assay = [[0.1, 0.2], [5.0, 6.0], [0.3, 0.4]]
  let pheno : Map[String, Array[String]] = Map([])
  let feature : Map[String, Array[String]] = Map([])

  let eset = @src.ExpressionSet::new(assay, pheno, feature, "", Map([]))
  let filtered = @src.expression_set_filter_by_row(eset, 1.0)

  assert_eq(@src.expression_set_nrow(filtered), 1)
}

///|
test "biobase_expression_set_summary" {
  let assay = [[1.0, 2.0]]
  let pheno : Map[String, Array[String]] = Map(
    [("sample_name", ["s1", "s2"])],
    capacity=1,
  )
  let feature : Map[String, Array[String]] = Map(
    [("feature_name", ["g1"])],
    capacity=1,
  )

  let eset = @src.ExpressionSet::new(assay, pheno, feature, "test", Map([]))
  let summary = @src.expression_set_summary(eset)

  assert_true(summary.contains("ExpressionSet Summary"))
  assert_true(summary.contains("Features: 1"))
  assert_true(summary.contains("Samples: 2"))
}
