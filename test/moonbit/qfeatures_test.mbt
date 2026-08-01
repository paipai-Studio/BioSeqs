///|
/// QFeatures module tests - quantitative proteomics data management.

///|
test "qf_assay_creation" {
  let row_names = ["r0", "r1", "r2"]
  let col_names = ["c0", "c1"]
  let data : Array[Array[Double]] = [
    [1.0, 2.0],
    [3.0, 4.0],
    [5.0, 6.0],
  ]
  let assay = @src.QfAssay::new("test", row_names, col_names, data)
  assert_eq(assay.name(), "test")
  assert_eq(assay.n_rows(), 3)
  assert_eq(assay.n_cols(), 2)
  assert_eq(assay.row_names().length(), 3)
  assert_eq(assay.col_names().length(), 2)
  assert_eq(assay.row_names()[0], "r0")
  assert_eq(assay.col_names()[1], "c1")
  assert_eq(assay.data().length(), 3)
  assert_eq(assay.data()[0].length(), 2)
}

///|
test "qf_assay_zeros" {
  let assay = @src.QfAssay::zeros("zero", 3, 4)
  assert_eq(assay.name(), "zero")
  assert_eq(assay.n_rows(), 3)
  assert_eq(assay.n_cols(), 4)
  assert_eq(assay.row_names().length(), 3)
  assert_eq(assay.col_names().length(), 4)
  assert_eq(assay.row_names()[0], "row0")
  assert_eq(assay.row_names()[2], "row2")
  assert_eq(assay.col_names()[0], "col0")
  assert_eq(assay.col_names()[3], "col3")
  for i in 0..<3 {
    for j in 0..<4 {
      assert_true((assay.get(i, j) - 0.0).abs() < 0.001)
    }
  }
}

///|
test "qf_assay_get_set" {
  let assay = @src.QfAssay::zeros("a", 2, 2)
  assert_true((assay.get(0, 0) - 0.0).abs() < 0.001)
  assay.set(0, 0, 1.5)
  assay.set(1, 1, 2.5)
  assert_true((assay.get(0, 0) - 1.5).abs() < 0.001)
  assert_true((assay.get(1, 1) - 2.5).abs() < 0.001)
  assert_true((assay.get(0, 1) - 0.0).abs() < 0.001)
  assert_true((assay.get(1, 0) - 0.0).abs() < 0.001)
}

///|
test "qf_assay_get_row_get_col" {
  let data : Array[Array[Double]] = [
    [1.0, 2.0, 3.0],
    [4.0, 5.0, 6.0],
  ]
  let assay = @src.QfAssay::new("a", ["r0", "r1"], ["c0", "c1", "c2"], data)
  let row0 = assay.get_row(0)
  let row1 = assay.get_row(1)
  assert_eq(row0.length(), 3)
  assert_eq(row1.length(), 3)
  assert_true((row0[0] - 1.0).abs() < 0.001)
  assert_true((row0[2] - 3.0).abs() < 0.001)
  assert_true((row1[0] - 4.0).abs() < 0.001)
  assert_true((row1[2] - 6.0).abs() < 0.001)
  let col0 = assay.get_col(0)
  let col2 = assay.get_col(2)
  assert_eq(col0.length(), 2)
  assert_eq(col2.length(), 2)
  assert_true((col0[0] - 1.0).abs() < 0.001)
  assert_true((col0[1] - 4.0).abs() < 0.001)
  assert_true((col2[0] - 3.0).abs() < 0.001)
  assert_true((col2[1] - 6.0).abs() < 0.001)
}

///|
test "qf_assay_row_col_stats" {
  let data : Array[Array[Double]] = [
    [1.0, 2.0, 3.0],
    [4.0, 5.0, 6.0],
  ]
  let assay = @src.QfAssay::new("a", ["r0", "r1"], ["c0", "c1", "c2"], data)
  let rmeans = assay.row_means()
  let cmeans = assay.col_means()
  let rsums = assay.row_sums()
  let csums = assay.col_sums()
  assert_eq(rmeans.length(), 2)
  assert_eq(cmeans.length(), 3)
  assert_true((rmeans[0] - 2.0).abs() < 0.001)
  assert_true((rmeans[1] - 5.0).abs() < 0.001)
  assert_true((cmeans[0] - 2.5).abs() < 0.001)
  assert_true((cmeans[1] - 3.5).abs() < 0.001)
  assert_true((cmeans[2] - 4.5).abs() < 0.001)
  assert_true((rsums[0] - 6.0).abs() < 0.001)
  assert_true((rsums[1] - 15.0).abs() < 0.001)
  assert_true((csums[0] - 5.0).abs() < 0.001)
  assert_true((csums[1] - 7.0).abs() < 0.001)
  assert_true((csums[2] - 9.0).abs() < 0.001)
}

///|
test "qf_assay_row_col_stats_with_nan" {
  let data : Array[Array[Double]] = [
    [1.0, Double::nan(), 3.0],
    [4.0, 5.0, Double::nan()],
  ]
  let assay = @src.QfAssay::new("a", ["r0", "r1"], ["c0", "c1", "c2"], data)
  let rmeans = assay.row_means()
  let rsums = assay.row_sums()
  // row 0: [1, NaN, 3] -> mean = (1+3)/2 = 2, sum = 4
  assert_true((rmeans[0] - 2.0).abs() < 0.001)
  assert_true((rsums[0] - 4.0).abs() < 0.001)
  // row 1: [4, 5, NaN] -> mean = (4+5)/2 = 4.5, sum = 9
  assert_true((rmeans[1] - 4.5).abs() < 0.001)
  assert_true((rsums[1] - 9.0).abs() < 0.001)
  let cmeans = assay.col_means()
  // col 1: [NaN, 5] -> mean = 5
  assert_true((cmeans[1] - 5.0).abs() < 0.001)
}

///|
test "qf_assay_to_string" {
  let data : Array[Array[Double]] = [[1.0, 2.0], [3.0, 4.0]]
  let assay = @src.QfAssay::new("myassay", ["r0", "r1"], ["c0", "c1"], data)
  let s = assay.to_string()
  assert_true(s.contains("myassay"))
  assert_true(s.contains("2 rows x 2 cols"))
}

///|
test "qf_feature_link_creation" {
  let link = @src.QfFeatureLink::new("psm", 0, "peptide", 0)
  assert_eq(link.from_assay(), "psm")
  assert_eq(link.from_row(), 0)
  assert_eq(link.to_assay(), "peptide")
  assert_eq(link.to_row(), 0)
}

///|
test "qf_feature_links_basic" {
  let links = @src.QfFeatureLinks::new()
  assert_eq(links.n_links(), 0)
  links.add_link(@src.QfFeatureLink::new("psm", 0, "peptide", 0))
  links.add_link(@src.QfFeatureLink::new("psm", 1, "peptide", 0))
  links.add_link(@src.QfFeatureLink::new("psm", 2, "peptide", 1))
  assert_eq(links.n_links(), 3)
  let from0 = links.get_links_from("psm", 0)
  assert_eq(from0.length(), 1)
  assert_eq(from0[0].to_assay(), "peptide")
  assert_eq(from0[0].to_row(), 0)
  let from2 = links.get_links_from("psm", 2)
  assert_eq(from2.length(), 1)
  assert_eq(from2[0].to_assay(), "peptide")
  assert_eq(from2[0].to_row(), 1)
  let from_other = links.get_links_from("psm", 99)
  assert_eq(from_other.length(), 0)
  let to0 = links.get_links_to("peptide", 0)
  assert_eq(to0.length(), 2)
  let to1 = links.get_links_to("peptide", 1)
  assert_eq(to1.length(), 1)
  let to_other = links.get_links_to("protein", 0)
  assert_eq(to_other.length(), 0)
}

///|
test "qf_features_empty" {
  let qf = @src.QFeatures::new()
  assert_eq(qf.n_assays(), 0)
  assert_eq(qf.assays().length(), 0)
  assert_eq(qf.assay_names().length(), 0)
  assert_eq(qf.links().n_links(), 0)
  assert_true(qf.get_assay("missing").is_empty())
}

///|
test "qf_features_add_assay_and_get" {
  let qf = @src.QFeatures::new()
  let a1 = @src.QfAssay::zeros("a1", 2, 3)
  let a2 = @src.QfAssay::zeros("a2", 4, 5)
  qf.add_assay(a1)
  qf.add_assay(a2)
  assert_eq(qf.n_assays(), 2)
  assert_eq(qf.assays().length(), 2)
  let names = qf.assay_names()
  assert_eq(names.length(), 2)
  assert_eq(names[0], "a1")
  assert_eq(names[1], "a2")
  match qf.get_assay("a1") {
    Some(a) => {
      assert_eq(a.n_rows(), 2)
      assert_eq(a.n_cols(), 3)
    }
    None => assert_true(false)
  }
  match qf.get_assay("a2") {
    Some(a) => assert_eq(a.n_rows(), 4)
    None => assert_true(false)
  }
  assert_true(qf.get_assay("nonexistent").is_empty())
}

///|
test "qf_features_add_link_and_counts" {
  let qf = @src.QFeatures::new()
  qf.add_assay(@src.QfAssay::zeros("psm", 6, 4))
  qf.add_assay(@src.QfAssay::zeros("peptide", 4, 4))
  qf.add_link("psm", 0, "peptide", 0)
  qf.add_link("psm", 1, "peptide", 0)
  qf.add_link("psm", 2, "peptide", 1)
  assert_eq(qf.links().n_links(), 3)
  assert_eq(qf.n_features("psm"), 6)
  assert_eq(qf.n_features("peptide"), 4)
  assert_eq(qf.n_samples("psm"), 4)
  assert_eq(qf.n_samples("peptide"), 4)
  assert_eq(qf.n_features("missing"), 0)
  assert_eq(qf.n_samples("missing"), 0)
}

///|
test "qf_aggregate_rows_sum" {
  let qf = @src.QFeatures::new()
  let data : Array[Array[Double]] = [
    [1.0, 2.0],
    [3.0, 4.0],
    [5.0, 6.0],
    [7.0, 8.0],
  ]
  qf.add_assay(
    @src.QfAssay::new("src", ["r0", "r1", "r2", "r3"], ["c0", "c1"], data),
  )
  let mapping : Array[(Int, String)] = [
    (0, "A"),
    (1, "A"),
    (2, "B"),
    (3, "B"),
  ]
  let qf2 = @src.qf_aggregate(qf, "src", "agg", mapping, "sum")
  assert_eq(qf2.n_assays(), 2)
  match qf2.get_assay("agg") {
    Some(a) => {
      assert_eq(a.n_rows(), 2)
      assert_eq(a.n_cols(), 2)
      // A: 1+3=4, 2+4=6
      assert_true((a.get(0, 0) - 4.0).abs() < 0.001)
      assert_true((a.get(0, 1) - 6.0).abs() < 0.001)
      // B: 5+7=12, 6+8=14
      assert_true((a.get(1, 0) - 12.0).abs() < 0.001)
      assert_true((a.get(1, 1) - 14.0).abs() < 0.001)
      assert_eq(a.row_names()[0], "A")
      assert_eq(a.row_names()[1], "B")
    }
    None => assert_true(false)
  }
  // Links should be added: 4 source rows -> 2 new rows
  assert_eq(qf2.links().n_links(), 4)
  let from0 = qf2.links().get_links_from("src", 0)
  assert_eq(from0.length(), 1)
  assert_eq(from0[0].to_assay(), "agg")
  assert_eq(from0[0].to_row(), 0)
}

///|
test "qf_aggregate_rows_mean" {
  let qf = @src.QFeatures::new()
  let data : Array[Array[Double]] = [
    [1.0, 2.0],
    [3.0, 4.0],
    [5.0, 6.0],
    [7.0, 8.0],
  ]
  qf.add_assay(
    @src.QfAssay::new("src", ["r0", "r1", "r2", "r3"], ["c0", "c1"], data),
  )
  let mapping : Array[(Int, String)] = [
    (0, "A"),
    (1, "A"),
    (2, "B"),
    (3, "B"),
  ]
  let qf2 = @src.qf_aggregate(qf, "src", "agg", mapping, "mean")
  match qf2.get_assay("agg") {
    Some(a) => {
      assert_eq(a.n_rows(), 2)
      // A: (1+3)/2=2, (2+4)/2=3
      assert_true((a.get(0, 0) - 2.0).abs() < 0.001)
      assert_true((a.get(0, 1) - 3.0).abs() < 0.001)
      // B: (5+7)/2=6, (6+8)/2=7
      assert_true((a.get(1, 0) - 6.0).abs() < 0.001)
      assert_true((a.get(1, 1) - 7.0).abs() < 0.001)
    }
    None => assert_true(false)
  }
}

///|
test "qf_aggregate_by_col_merge" {
  let qf = @src.QFeatures::new()
  let data : Array[Array[Double]] = [
    [1.0, 2.0, 3.0, 4.0],
    [5.0, 6.0, 7.0, 8.0],
  ]
  qf.add_assay(
    @src.QfAssay::new(
      "src",
      ["r0", "r1"],
      ["c0", "c1", "c2", "c3"],
      data,
    ),
  )
  let labels = ["rep1", "rep1", "rep2", "rep2"]
  // Sum aggregation.
  let qf_sum = @src.qf_aggregate_by_col(qf, "src", "agg_sum", labels, "sum")
  match qf_sum.get_assay("agg_sum") {
    Some(a) => {
      assert_eq(a.n_rows(), 2)
      assert_eq(a.n_cols(), 2)
      // row 0: [1+2, 3+4] = [3, 7]
      assert_true((a.get(0, 0) - 3.0).abs() < 0.001)
      assert_true((a.get(0, 1) - 7.0).abs() < 0.001)
      // row 1: [5+6, 7+8] = [11, 15]
      assert_true((a.get(1, 0) - 11.0).abs() < 0.001)
      assert_true((a.get(1, 1) - 15.0).abs() < 0.001)
      assert_eq(a.col_names()[0], "rep1")
      assert_eq(a.col_names()[1], "rep2")
    }
    None => assert_true(false)
  }
  // Mean aggregation.
  let qf_mean = @src.qf_aggregate_by_col(qf, "src", "agg_mean", labels, "mean")
  match qf_mean.get_assay("agg_mean") {
    Some(a) => {
      // row 0: [(1+2)/2, (3+4)/2] = [1.5, 3.5]
      assert_true((a.get(0, 0) - 1.5).abs() < 0.001)
      assert_true((a.get(0, 1) - 3.5).abs() < 0.001)
      // row 1: [(5+6)/2, (7+8)/2] = [5.5, 7.5]
      assert_true((a.get(1, 0) - 5.5).abs() < 0.001)
      assert_true((a.get(1, 1) - 7.5).abs() < 0.001)
    }
    None => assert_true(false)
  }
}

///|
test "qf_filter_features" {
  let qf = @src.QFeatures::new()
  let data : Array[Array[Double]] = [
    [1.0, 2.0],
    [3.0, 4.0],
    [5.0, 6.0],
    [7.0, 8.0],
  ]
  qf.add_assay(
    @src.QfAssay::new("a", ["r0", "r1", "r2", "r3"], ["c0", "c1"], data),
  )
  let qf2 = @src.qf_filter_features(qf, "a", [0, 2])
  match qf2.get_assay("a") {
    Some(a) => {
      assert_eq(a.n_rows(), 2)
      assert_eq(a.n_cols(), 2)
      assert_true((a.get(0, 0) - 1.0).abs() < 0.001)
      assert_true((a.get(0, 1) - 2.0).abs() < 0.001)
      assert_true((a.get(1, 0) - 5.0).abs() < 0.001)
      assert_true((a.get(1, 1) - 6.0).abs() < 0.001)
      assert_eq(a.row_names()[0], "r0")
      assert_eq(a.row_names()[1], "r2")
    }
    None => assert_true(false)
  }
}

///|
test "qf_filter_samples" {
  let qf = @src.QFeatures::new()
  let data : Array[Array[Double]] = [
    [1.0, 2.0, 3.0, 4.0],
    [5.0, 6.0, 7.0, 8.0],
  ]
  qf.add_assay(
    @src.QfAssay::new("a", ["r0", "r1"], ["c0", "c1", "c2", "c3"], data),
  )
  let qf2 = @src.qf_filter_samples(qf, "a", [1, 3])
  match qf2.get_assay("a") {
    Some(a) => {
      assert_eq(a.n_rows(), 2)
      assert_eq(a.n_cols(), 2)
      assert_true((a.get(0, 0) - 2.0).abs() < 0.001)
      assert_true((a.get(0, 1) - 4.0).abs() < 0.001)
      assert_true((a.get(1, 0) - 6.0).abs() < 0.001)
      assert_true((a.get(1, 1) - 8.0).abs() < 0.001)
      assert_eq(a.col_names()[0], "c1")
      assert_eq(a.col_names()[1], "c3")
    }
    None => assert_true(false)
  }
}

///|
test "qf_filter_na" {
  let qf = @src.QFeatures::new()
  let nan = Double::nan()
  let data : Array[Array[Double]] = [
    [1.0, 2.0, 3.0],
    [nan, nan, 6.0],
    [7.0, nan, 9.0],
  ]
  qf.add_assay(@src.QfAssay::new("a", ["r0", "r1", "r2"], ["c0", "c1", "c2"], data))
  // max_na_frac = 0.5 -> row 1 (2/3 NaN > 0.5) is dropped; row 2 (1/3 NaN) kept.
  let qf2 = @src.qf_filter_na(qf, "a", 0.5)
  match qf2.get_assay("a") {
    Some(a) => {
      assert_eq(a.n_rows(), 2)
      assert_eq(a.row_names()[0], "r0")
      assert_eq(a.row_names()[1], "r2")
    }
    None => assert_true(false)
  }
  // Stricter threshold 0.3 drops both row 1 (2/3) and row 2 (1/3 > 0.3).
  let qf3 = @src.qf_filter_na(qf, "a", 0.3)
  match qf3.get_assay("a") {
    Some(a) => assert_eq(a.n_rows(), 1)
    None => assert_true(false)
  }
}

///|
test "qf_filter_low_abundance" {
  let qf = @src.QFeatures::new()
  let data : Array[Array[Double]] = [
    [1.0, 2.0, 3.0],
    [0.5, 0.5, 0.5],
    [10.0, 20.0, 30.0],
  ]
  qf.add_assay(@src.QfAssay::new("a", ["r0", "r1", "r2"], ["c0", "c1", "c2"], data))
  // threshold = 1.0 -> row 1 (all 0.5 < 1.0) is dropped.
  let qf2 = @src.qf_filter_low_abundance(qf, "a", 1.0)
  match qf2.get_assay("a") {
    Some(a) => {
      assert_eq(a.n_rows(), 2)
      assert_eq(a.row_names()[0], "r0")
      assert_eq(a.row_names()[1], "r2")
    }
    None => assert_true(false)
  }
}

///|
test "qf_normalize_quantiles" {
  // Construct a matrix where columns have different distributions but
  // after quantile normalization they should share the same sorted values.
  let data : Array[Array[Double]] = [
    [5.0, 10.0],
    [2.0, 20.0],
    [10.0, 30.0],
  ]
  let assay = @src.QfAssay::new("a", ["r0", "r1", "r2"], ["c0", "c1"], data)
  let norm = @src.qf_normalize_quantiles(assay)
  // Compute sorted columns of normalized assay; they should be equal.
  let col0 = norm.get_col(0)
  let col1 = norm.get_col(1)
  col0.sort_by(fn(a : Double, b : Double) -> Int {
    if a < b { -1 } else if a > b { 1 } else { 0 }
  })
  col1.sort_by(fn(a : Double, b : Double) -> Int {
    if a < b { -1 } else if a > b { 1 } else { 0 }
  })
  for i in 0..<3 {
    assert_true((col0[i] - col1[i]).abs() < 0.001)
  }
  // Expected rank means: (2+10)/2=6, (5+20)/2=12.5, (10+30)/2=20.
  assert_true((col0[0] - 6.0).abs() < 0.001)
  assert_true((col0[1] - 12.5).abs() < 0.001)
  assert_true((col0[2] - 20.0).abs() < 0.001)
}

///|
test "qf_normalize_center_scale" {
  // Column 0: [1, 2, 3] -> mean=2, sample std=1
  // Column 1: [10, 20, 30] -> mean=20, sample std=10
  let data : Array[Array[Double]] = [
    [1.0, 10.0],
    [2.0, 20.0],
    [3.0, 30.0],
  ]
  let assay = @src.QfAssay::new("a", ["r0", "r1", "r2"], ["c0", "c1"], data)
  let norm = @src.qf_normalize_center_scale(assay, true, true)
  // Expected normalized values: col0 = [-1, 0, 1], col1 = [-1, 0, 1]
  for i in 0..<3 {
    let expected = i.to_double() - 1.0
    assert_true((norm.get(i, 0) - expected).abs() < 0.001)
    assert_true((norm.get(i, 1) - expected).abs() < 0.001)
  }
  // Verify column means ~0 and column std ~1.
  let cmeans = norm.col_means()
  assert_true(cmeans[0].abs() < 0.001)
  assert_true(cmeans[1].abs() < 0.001)
}

///|
test "qf_normalize_center_only" {
  let data : Array[Array[Double]] = [
    [1.0, 10.0],
    [3.0, 30.0],
  ]
  let assay = @src.QfAssay::new("a", ["r0", "r1"], ["c0", "c1"], data)
  let norm = @src.qf_normalize_center_scale(assay, true, false)
  // Col 0 mean=2 -> centered: [-1, 1]
  // Col 1 mean=20 -> centered: [-10, 10]
  assert_true((norm.get(0, 0) - -1.0).abs() < 0.001)
  assert_true((norm.get(1, 0) - 1.0).abs() < 0.001)
  assert_true((norm.get(0, 1) - -10.0).abs() < 0.001)
  assert_true((norm.get(1, 1) - 10.0).abs() < 0.001)
}

///|
test "qf_log_transform" {
  // log2(x+1): 1->1, 3->2, 7->3, 0->0
  let data : Array[Array[Double]] = [
    [1.0, 3.0],
    [7.0, 0.0],
  ]
  let assay = @src.QfAssay::new("a", ["r0", "r1"], ["c0", "c1"], data)
  let logged = @src.qf_log_transform(assay, 2.0, 1.0)
  assert_true((logged.get(0, 0) - 1.0).abs() < 0.001)
  assert_true((logged.get(0, 1) - 2.0).abs() < 0.001)
  assert_true((logged.get(1, 0) - 3.0).abs() < 0.001)
  assert_true((logged.get(1, 1) - 0.0).abs() < 0.001)
}

///|
test "qf_log_transform_preserves_nan" {
  let data : Array[Array[Double]] = [[Double::nan(), 1.0]]
  let assay = @src.QfAssay::new("a", ["r0"], ["c0", "c1"], data)
  let logged = @src.qf_log_transform(assay, 2.0, 1.0)
  assert_true(@src.qf_is_na(logged.get(0, 0)))
  assert_true((logged.get(0, 1) - 1.0).abs() < 0.001)
}

///|
test "qf_impute_knn" {
  // Row 0 has a NaN at column 1.
  // Distance(row0, row1) on shared columns [0, 2]: sqrt(((1-2)^2 + (3-4)^2)/2) = 1
  // Distance(row0, row2) on shared columns [0, 2]: sqrt(((1-3)^2 + (3-7)^2)/2) = sqrt(10)
  // So row 1 is the nearest neighbor; with k=1 the imputed value should be row1[1] = 5.0.
  let data : Array[Array[Double]] = [
    [1.0, Double::nan(), 3.0],
    [2.0, 5.0, 4.0],
    [3.0, 6.0, 7.0],
  ]
  let assay = @src.QfAssay::new("a", ["r0", "r1", "r2"], ["c0", "c1", "c2"], data)
  let imputed = @src.qf_impute_knn(assay, 1)
  assert_true(!@src.qf_is_na(imputed.get(0, 1)))
  assert_true((imputed.get(0, 1) - 5.0).abs() < 0.001)
  // Non-NaN values should be unchanged.
  assert_true((imputed.get(0, 0) - 1.0).abs() < 0.001)
  assert_true((imputed.get(1, 1) - 5.0).abs() < 0.001)
}

///|
test "qf_impute_mean" {
  // Col means (skip NaN): col0 = (1+4)/2 = 2.5, col1 = 5.0, col2 = (3+6)/2 = 4.5
  let data : Array[Array[Double]] = [
    [1.0, Double::nan(), 3.0],
    [4.0, 5.0, 6.0],
  ]
  let assay = @src.QfAssay::new("a", ["r0", "r1"], ["c0", "c1", "c2"], data)
  let imputed = @src.qf_impute_mean(assay)
  assert_true((imputed.get(0, 1) - 5.0).abs() < 0.001)
  // Non-NaN values unchanged.
  assert_true((imputed.get(0, 0) - 1.0).abs() < 0.001)
  assert_true((imputed.get(0, 2) - 3.0).abs() < 0.001)
  assert_true((imputed.get(1, 0) - 4.0).abs() < 0.001)
  assert_true((imputed.get(1, 1) - 5.0).abs() < 0.001)
  assert_true((imputed.get(1, 2) - 6.0).abs() < 0.001)
}

///|
test "qf_impute_zero" {
  let data : Array[Array[Double]] = [
    [1.0, Double::nan()],
    [Double::nan(), 2.0],
  ]
  let assay = @src.QfAssay::new("a", ["r0", "r1"], ["c0", "c1"], data)
  let imputed = @src.qf_impute_zero(assay)
  assert_true((imputed.get(0, 0) - 1.0).abs() < 0.001)
  assert_true((imputed.get(0, 1) - 0.0).abs() < 0.001)
  assert_true((imputed.get(1, 0) - 0.0).abs() < 0.001)
  assert_true((imputed.get(1, 1) - 2.0).abs() < 0.001)
}

///|
test "qf_count_na_and_is_na" {
  assert_true(@src.qf_is_na(Double::nan()))
  assert_true(!@src.qf_is_na(0.0))
  assert_true(!@src.qf_is_na(1.5))
  assert_true(!@src.qf_is_na(-1.0))
  let data : Array[Array[Double]] = [
    [1.0, Double::nan(), 3.0],
    [Double::nan(), Double::nan(), 6.0],
  ]
  let assay = @src.QfAssay::new("a", ["r0", "r1"], ["c0", "c1", "c2"], data)
  assert_eq(@src.qf_count_na(assay), 3)
  let by_row = @src.qf_count_na_by_row(assay)
  assert_eq(by_row.length(), 2)
  assert_eq(by_row[0], 1)
  assert_eq(by_row[1], 2)
}

///|
test "qf_summary" {
  let qf = @src.qf_sample_data()
  let s = @src.qf_summary(qf)
  assert_true(s.length() > 0)
  assert_true(s.contains("QFeatures Summary"))
  assert_true(s.contains("Number of assays: 3"))
  assert_true(s.contains("psm"))
  assert_true(s.contains("peptide"))
  assert_true(s.contains("protein"))
}

///|
test "qf_assay_summary" {
  let data : Array[Array[Double]] = [
    [1.0, 2.0],
    [3.0, 4.0],
  ]
  let assay = @src.QfAssay::new("myassay", ["r0", "r1"], ["c0", "c1"], data)
  let s = @src.qf_assay_summary(assay)
  assert_true(s.length() > 0)
  assert_true(s.contains("myassay"))
  assert_true(s.contains("2 rows x 2 cols"))
  assert_true(s.contains("NA count: 0 / 4"))
}

///|
test "qf_to_long_format" {
  let qf = @src.QFeatures::new()
  let data : Array[Array[Double]] = [
    [1.0, 2.0],
    [3.0, 4.0],
  ]
  qf.add_assay(@src.QfAssay::new("a", ["r0", "r1"], ["c0", "c1"], data))
  let long = @src.qf_to_long_format(qf, "a")
  assert_eq(long.length(), 4)
  assert_eq(long[0].0, "r0")
  assert_eq(long[0].1, "c0")
  assert_true((long[0].2 - 1.0).abs() < 0.001)
  assert_eq(long[3].0, "r1")
  assert_eq(long[3].1, "c1")
  assert_true((long[3].2 - 4.0).abs() < 0.001)
  // Non-existent assay returns empty.
  let empty = @src.qf_to_long_format(qf, "missing")
  assert_eq(empty.length(), 0)
}

///|
test "qf_sample_data_structure" {
  let qf = @src.qf_sample_data()
  assert_eq(qf.n_assays(), 3)
  let names = qf.assay_names()
  assert_eq(names[0], "psm")
  assert_eq(names[1], "peptide")
  assert_eq(names[2], "protein")
  // PSM: 6 rows x 4 cols
  match qf.get_assay("psm") {
    Some(a) => {
      assert_eq(a.n_rows(), 6)
      assert_eq(a.n_cols(), 4)
    }
    None => assert_true(false)
  }
  // Peptide: 4 rows x 4 cols
  match qf.get_assay("peptide") {
    Some(a) => {
      assert_eq(a.n_rows(), 4)
      assert_eq(a.n_cols(), 4)
    }
    None => assert_true(false)
  }
  // Protein: 3 rows x 4 cols
  match qf.get_assay("protein") {
    Some(a) => {
      assert_eq(a.n_rows(), 3)
      assert_eq(a.n_cols(), 4)
    }
    None => assert_true(false)
  }
  // Links connect psm->peptide (6 links) and peptide->protein (4 links).
  assert_eq(qf.links().n_links(), 10)
}

///|
test "qf_sample_data_values" {
  let qf = @src.qf_sample_data()
  // Verify peptide PepA = mean(PSM1, PSM2) per sample.
  // PSM1: [10, 12, 11, 13], PSM2: [20, 22, 21, 23]
  // PepA: [15, 17, 16, 18]
  let pep = match qf.get_assay("peptide") {
    Some(a) => a
    None => {
      assert_true(false)
      return
    }
  }
  assert_eq(pep.row_names()[0], "PepA")
  assert_true((pep.get(0, 0) - 15.0).abs() < 0.001)
  assert_true((pep.get(0, 1) - 17.0).abs() < 0.001)
  assert_true((pep.get(0, 2) - 16.0).abs() < 0.001)
  assert_true((pep.get(0, 3) - 18.0).abs() < 0.001)
  // Verify protein Prot1 = mean(PepA, PepB) per sample.
  // PepB = mean(PSM3, PSM4) = [35, 37, 36, 38]
  // Prot1 = [25, 27, 26, 28]
  let prot = match qf.get_assay("protein") {
    Some(a) => a
    None => {
      assert_true(false)
      return
    }
  }
  assert_eq(prot.row_names()[0], "Prot1")
  assert_true((prot.get(0, 0) - 25.0).abs() < 0.001)
  assert_true((prot.get(0, 1) - 27.0).abs() < 0.001)
  assert_true((prot.get(0, 2) - 26.0).abs() < 0.001)
  assert_true((prot.get(0, 3) - 28.0).abs() < 0.001)
}

///|
test "qf_edge_case_empty_qfeatures" {
  let qf = @src.QFeatures::new()
  assert_eq(qf.n_assays(), 0)
  assert_true(qf.get_assay("anything").is_empty())
  assert_eq(qf.n_features("anything"), 0)
  assert_eq(qf.n_samples("anything"), 0)
  // Aggregating from a missing assay should return the original qf unchanged.
  let mapping : Array[(Int, String)] = [(0, "g")]
  let qf2 = @src.qf_aggregate(qf, "missing", "new", mapping, "sum")
  assert_eq(qf2.n_assays(), 0)
  let qf3 = @src.qf_filter_features(qf, "missing", [0])
  assert_eq(qf3.n_assays(), 0)
}

///|
test "qf_edge_case_single_cell_assay" {
  let data : Array[Array[Double]] = [[42.0]]
  let assay = @src.QfAssay::new("single", ["r0"], ["c0"], data)
  assert_eq(assay.n_rows(), 1)
  assert_eq(assay.n_cols(), 1)
  assert_true((assay.get(0, 0) - 42.0).abs() < 0.001)
  assert_true((assay.row_means()[0] - 42.0).abs() < 0.001)
  assert_true((assay.col_means()[0] - 42.0).abs() < 0.001)
  assert_true((assay.row_sums()[0] - 42.0).abs() < 0.001)
  assert_true((assay.col_sums()[0] - 42.0).abs() < 0.001)
  let s = assay.to_string()
  assert_true(s.contains("1 rows x 1 cols"))
}

///|
test "qf_edge_case_all_nan_assay" {
  let nan = Double::nan()
  let data : Array[Array[Double]] = [
    [nan, nan],
    [nan, nan],
  ]
  let assay = @src.QfAssay::new("allnan", ["r0", "r1"], ["c0", "c1"], data)
  assert_eq(@src.qf_count_na(assay), 4)
  let by_row = @src.qf_count_na_by_row(assay)
  assert_eq(by_row[0], 2)
  assert_eq(by_row[1], 2)
  // Row means should be NaN when all values are NaN.
  let rmeans = assay.row_means()
  assert_true(@src.qf_is_na(rmeans[0]))
  assert_true(@src.qf_is_na(rmeans[1]))
  let cmeans = assay.col_means()
  assert_true(@src.qf_is_na(cmeans[0]))
  assert_true(@src.qf_is_na(cmeans[1]))
  // Row sums skip NaN, so they are 0.
  let rsums = assay.row_sums()
  assert_true((rsums[0] - 0.0).abs() < 0.001)
  assert_true((rsums[1] - 0.0).abs() < 0.001)
  // Impute zero should replace all NaN with 0.0.
  let imputed = @src.qf_impute_zero(assay)
  assert_eq(@src.qf_count_na(imputed), 0)
  assert_true((imputed.get(0, 0) - 0.0).abs() < 0.001)
  assert_true((imputed.get(1, 1) - 0.0).abs() < 0.001)
}
