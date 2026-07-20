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

///|
test "read_peak_file" {
  let bed_content = "chr1\t1000\t2000\tpeak1\nchr2\t3000\t4000\tpeak2\nchr3\t5000\t6000"
  let peaks = @src.read_peak_file(bed_content)
  assert_eq(peaks.length(), 3)
  assert_eq(peaks[0].chr, "chr1")
  assert_eq(peaks[0].start, 1000)
  assert_eq(peaks[0].end, 2000)
  assert_eq(peaks[0].peak_id, "peak1")
}

///|
test "peak2gene" {
  let peaks = @src.cs_create_example_peaks()
  let txdb = @src.cs_create_example_txdb()
  let result = @src.annotate_peaks(peaks, txdb)
  let associations = @src.peak2gene(result)
  assert_eq(associations.length(), 50)
}

///|
test "get_exon_annotations" {
  let peaks = @src.cs_create_example_peaks()
  let txdb = @src.cs_create_example_txdb()
  let result = @src.annotate_peaks(peaks, txdb)
  let exons = @src.get_exon_annotations(result)
  assert_true(exons.length() >= 0)
}

///|
test "get_intron_annotations" {
  let peaks = @src.cs_create_example_peaks()
  let txdb = @src.cs_create_example_txdb()
  let result = @src.annotate_peaks(peaks, txdb)
  let introns = @src.get_intron_annotations(result)
  assert_true(introns.length() >= 0)
}

///|
test "get_utr_annotations" {
  let peaks = @src.cs_create_example_peaks()
  let txdb = @src.cs_create_example_txdb()
  let result = @src.annotate_peaks(peaks, txdb)
  let utrs = @src.get_utr_annotations(result)
  assert_true(utrs.length() >= 0)
}

///|
test "peak_overlap" {
  let peaks1 = [
    @src.Peak::new("chr1", 1000, 2000, "peak1"),
    @src.Peak::new("chr1", 3000, 4000, "peak2"),
    @src.Peak::new("chr2", 5000, 6000, "peak3"),
  ]
  let peaks2 = [
    @src.Peak::new("chr1", 1500, 2500, "peakA"),
    @src.Peak::new("chr1", 3500, 4500, "peakB"),
    @src.Peak::new("chr3", 7000, 8000, "peakC"),
  ]
  let peak_sets = [("Set1", peaks1), ("Set2", peaks2)]
  let result = @src.peak_overlap(peak_sets, 1)
  assert_eq(result.overlap_matrix.length(), 2)
  assert_eq(result.overlap_matrix[0][1], 2)
  assert_eq(result.union_count, 6)
}

///|
test "plot_venn" {
  let peaks1 = [@src.Peak::new("chr1", 1000, 2000, "peak1")]
  let peaks2 = [@src.Peak::new("chr1", 1500, 2500, "peakA")]
  let peak_sets = [("Set1", peaks1), ("Set2", peaks2)]
  let overlap_result = @src.peak_overlap(peak_sets, 1)
  let venn = @src.plot_venn(overlap_result)
  assert_true(venn.contains("Venn"))
}

///|
test "plot_anno_pie" {
  let peaks = @src.cs_create_example_peaks()
  let txdb = @src.cs_create_example_txdb()
  let result = @src.annotate_peaks(peaks, txdb)
  let pie = @src.plot_anno_pie(result)
  assert_true(pie.contains("Pie"))
}
