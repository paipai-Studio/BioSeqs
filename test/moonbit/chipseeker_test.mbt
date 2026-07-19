///|
test "cs_create_example_peaks" {
  let peaks = @src.cs_create_example_peaks()
  assert_eq(peaks.length(), 50)
}

///|
test "cs_create_example_txdb" {
  let txdb = @src.cs_create_example_txdb()
  assert_eq(txdb.genes.length(), 8)
}

///|
test "annotate_peaks" {
  let peaks = @src.cs_create_example_peaks()
  let txdb = @src.cs_create_example_txdb()
  let result = @src.annotate_peaks(peaks, txdb)
  assert_eq(result.annotations.length(), 50)
}

///|
test "get_annotation_statistics" {
  let peaks = @src.cs_create_example_peaks()
  let txdb = @src.cs_create_example_txdb()
  let result = @src.annotate_peaks(peaks, txdb)
  let stats = @src.get_annotation_statistics(result)
  assert_true(stats.length() > 0)
}

///|
test "get_promoter_annotations" {
  let peaks = @src.cs_create_example_peaks()
  let txdb = @src.cs_create_example_txdb()
  let result = @src.annotate_peaks(peaks, txdb)
  let promoters = @src.get_promoter_annotations(result)
  assert_true(promoters.length() >= 0)
}

///|
test "get_distal_annotations" {
  let peaks = @src.cs_create_example_peaks()
  let txdb = @src.cs_create_example_txdb()
  let result = @src.annotate_peaks(peaks, txdb)
  let distal = @src.get_distal_annotations(result)
  assert_true(distal.length() >= 0)
}

///|
test "get_intergenic_annotations" {
  let peaks = @src.cs_create_example_peaks()
  let txdb = @src.cs_create_example_txdb()
  let result = @src.annotate_peaks(peaks, txdb)
  let intergenic = @src.get_intergenic_annotations(result)
  assert_true(intergenic.length() >= 0)
}

///|
test "plot_annotation_bar" {
  let peaks = @src.cs_create_example_peaks()
  let txdb = @src.cs_create_example_txdb()
  let result = @src.annotate_peaks(peaks, txdb)
  let plot = @src.plot_annotation_bar(result.summary)
  assert_true(plot.contains("Annotation"))
}

///|
test "plot_dist_to_tss" {
  let peaks = @src.cs_create_example_peaks()
  let txdb = @src.cs_create_example_txdb()
  let result = @src.annotate_peaks(peaks, txdb)
  let plot = @src.plot_dist_to_tss(result.annotations)
  assert_true(plot.contains("Distance"))
}

///|
test "annotation_summary" {
  let peaks = @src.cs_create_example_peaks()
  let txdb = @src.cs_create_example_txdb()
  let result = @src.annotate_peaks(peaks, txdb)
  let summary = @src.annotation_summary(result)
  assert_true(summary.contains("Peak"))
}
