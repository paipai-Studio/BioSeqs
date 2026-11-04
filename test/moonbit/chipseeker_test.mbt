///|
/// Tests for the ChIPseeker peak annotation module.

///|
test "cs_peak_construction_full" {
  let p = @src.Peak::new(
    chr="chr1",
    start=1000,
    end=2000,
    summit=1500,
    fold_enrichment=10.5,
    p_value=0.00001,
  )
  assert_eq(p.chr, "chr1")
  assert_eq(p.start, 1000)
  assert_eq(p.end, 2000)
  assert_eq(p.summit, 1500)
  assert_true(p.fold_enrichment == 10.5)
  assert_true(p.p_value == 0.00001)
}

///|
test "cs_peak_construction_optional" {
  // Only required args (chr, start) provided; rest default to zero.
  let p = @src.Peak::new(chr="chr2", start=500)
  assert_eq(p.chr, "chr2")
  assert_eq(p.start, 500)
  assert_eq(p.end, 0)
  assert_eq(p.summit, 0)
  assert_true(p.fold_enrichment == 0.0)
  assert_true(p.p_value == 0.0)
}

///|
test "cs_peak_accessors" {
  let p = @src.Peak::new(
    chr="chrX",
    start=100,
    end=200,
    summit=150,
    fold_enrichment=3.0,
    p_value=0.01,
  )
  assert_eq(p.get_chr(), "chrX")
  assert_eq(p.get_start(), 100)
  assert_eq(p.get_end(), 200)
  assert_eq(p.get_summit(), 150)
  assert_true(p.get_fold_enrichment() == 3.0)
  assert_true(p.get_p_value() == 0.01)
}

///|
test "cs_peak_set_end" {
  let p = @src.Peak::new(chr="chr1", start=100, end=200)
  assert_eq(p.end, 200)
  p.set_end(300)
  assert_eq(p.end, 300)
  assert_eq(p.get_end(), 300)
}

///|
test "cs_peak_position_with_summit" {
  let p = @src.Peak::new(chr="chr1", start=1000, end=2000, summit=1700)
  assert_eq(@src.cs_peak_position(p), 1700)
}

///|
test "cs_peak_position_without_summit" {
  // summit defaults to 0, so position should be midpoint.
  let p = @src.Peak::new(chr="chr1", start=1000, end=2001)
  assert_eq(@src.cs_peak_position(p), 1500)
}

///|
test "cs_gene_annotation_plus_strand" {
  let g = @src.GeneAnnotation::new(
    gene_id="G1",
    gene_name="Gene1",
    chr="chr1",
    gene_start=10000,
    gene_end=25000,
    strand="+",
  )
  assert_eq(g.gene_id, "G1")
  assert_eq(g.gene_name, "Gene1")
  assert_eq(g.chr, "chr1")
  assert_eq(g.gene_start, 10000)
  assert_eq(g.gene_end, 25000)
  assert_eq(g.strand, "+")
  // For + strand, TSS = gene_start.
  assert_eq(g.tss, 10000)
}

///|
test "cs_gene_annotation_minus_strand" {
  let g = @src.GeneAnnotation::new(
    gene_id="G2",
    gene_name="Gene2",
    chr="chr1",
    gene_start=50000,
    gene_end=70000,
    strand="-",
  )
  assert_eq(g.gene_id, "G2")
  assert_eq(g.strand, "-")
  // For - strand, TSS = gene_end.
  assert_eq(g.tss, 70000)
}

///|
test "cs_gene_annotation_accessors" {
  let g = @src.GeneAnnotation::new(
    gene_id="G3",
    gene_name="Gene3",
    chr="chr2",
    gene_start=1000,
    gene_end=2000,
    strand="+",
  )
  assert_eq(g.get_gene_id(), "G3")
  assert_eq(g.get_gene_name(), "Gene3")
  assert_eq(g.get_chr(), "chr2")
  assert_eq(g.get_gene_start(), 1000)
  assert_eq(g.get_gene_end(), 2000)
  assert_eq(g.get_strand(), "+")
  assert_eq(g.get_tss(), 1000)
}

///|
test "cs_tss_distance_positive" {
  // Peak downstream of TSS on + strand gene -> positive distance.
  let g = @src.GeneAnnotation::new(
    gene_id="G1",
    gene_name="Gene1",
    chr="chr1",
    gene_start=10000,
    gene_end=25000,
    strand="+",
  )
  let p = @src.Peak::new(chr="chr1", start=12000, end=12100, summit=12000)
  assert_eq(@src.compute_tss_distance(p, g), 2000)
}

///|
test "cs_tss_distance_negative" {
  // Peak upstream of TSS on + strand gene -> negative distance.
  let g = @src.GeneAnnotation::new(
    gene_id="G1",
    gene_name="Gene1",
    chr="chr1",
    gene_start=10000,
    gene_end=25000,
    strand="+",
  )
  let p = @src.Peak::new(chr="chr1", start=8000, end=8100, summit=8000)
  assert_eq(@src.compute_tss_distance(p, g), -2000)
}

///|
test "cs_tss_distance_zero" {
  // Peak exactly at TSS -> zero distance.
  let g = @src.GeneAnnotation::new(
    gene_id="G1",
    gene_name="Gene1",
    chr="chr1",
    gene_start=10000,
    gene_end=25000,
    strand="+",
  )
  let p = @src.Peak::new(chr="chr1", start=9900, end=10100, summit=10000)
  assert_eq(@src.compute_tss_distance(p, g), 0)
}

///|
test "cs_tss_distance_minus_strand" {
  // For - strand, TSS = gene_end. Distance is TSS - peak_pos.
  let g = @src.GeneAnnotation::new(
    gene_id="G2",
    gene_name="Gene2",
    chr="chr1",
    gene_start=50000,
    gene_end=70000,
    strand="-",
  )
  // Peak upstream (in genomic coords) of TSS => downstream in transcription => positive.
  let p1 = @src.Peak::new(chr="chr1", start=60000, end=60100, summit=60000)
  assert_eq(@src.compute_tss_distance(p1, g), 10000)
  // Peak downstream (in genomic coords) of TSS => upstream in transcription => negative.
  let p2 = @src.Peak::new(chr="chr1", start=75000, end=75100, summit=75000)
  assert_eq(@src.compute_tss_distance(p2, g), -5000)
}

///|
test "cs_feature_promoter" {
  let g = @src.GeneAnnotation::new(
    gene_id="G1",
    gene_name="Gene1",
    chr="chr1",
    gene_start=10000,
    gene_end=25000,
    strand="+",
  )
  // Peak 2000bp upstream of TSS, within default window 3000.
  let p = @src.Peak::new(chr="chr1", start=7900, end=8100, summit=8000)
  assert_eq(@src.assign_feature(p, g, 3000), "Promoter")
}

///|
test "cs_feature_5utr" {
  let g = @src.GeneAnnotation::new(
    gene_id="G1",
    gene_name="Gene1",
    chr="chr1",
    gene_start=10000,
    gene_end=25000,
    strand="+",
  )
  // 5'UTR = [TSS, TSS + 500] = [10000, 10500]. Use small window 100 so the
  // peak at 10300 is in 5'UTR but not in promoter.
  let p = @src.Peak::new(chr="chr1", start=10200, end=10400, summit=10300)
  assert_eq(@src.assign_feature(p, g, 100), "5'UTR")
}

///|
test "cs_feature_3utr" {
  let g = @src.GeneAnnotation::new(
    gene_id="G1",
    gene_name="Gene1",
    chr="chr1",
    gene_start=10000,
    gene_end=25000,
    strand="+",
  )
  // 3'UTR = [gene_end - 500, gene_end] = [24500, 25000].
  let p = @src.Peak::new(chr="chr1", start=24600, end=24700, summit=24700)
  assert_eq(@src.assign_feature(p, g, 3000), "3'UTR")
}

///|
test "cs_feature_exon" {
  let g = @src.GeneAnnotation::new(
    gene_id="G1",
    gene_name="Gene1",
    chr="chr1",
    gene_start=10000,
    gene_end=25000,
    strand="+",
  )
  // Exon 2 = [gene_end - 500 - 1000, gene_end - 500] = [23500, 24500].
  let p = @src.Peak::new(chr="chr1", start=23900, end=24100, summit=24000)
  assert_eq(@src.assign_feature(p, g, 3000), "Exon")
}

///|
test "cs_feature_intron" {
  let g = @src.GeneAnnotation::new(
    gene_id="G1",
    gene_name="Gene1",
    chr="chr1",
    gene_start=10000,
    gene_end=25000,
    strand="+",
  )
  // Intron = [TSS + 500 + 1000, gene_end - 500 - 1000] = [11500, 23500].
  // Peak at 14000 is well inside the intron and outside the promoter window.
  let p = @src.Peak::new(chr="chr1", start=13900, end=14100, summit=14000)
  assert_eq(@src.assign_feature(p, g, 3000), "Intron")
}

///|
test "cs_feature_downstream" {
  let g = @src.GeneAnnotation::new(
    gene_id="G1",
    gene_name="Gene1",
    chr="chr1",
    gene_start=10000,
    gene_end=25000,
    strand="+",
  )
  // Downstream = [gene_end, gene_end + 3000] = [25000, 28000].
  let p = @src.Peak::new(chr="chr1", start=25900, end=26100, summit=26000)
  assert_eq(@src.assign_feature(p, g, 3000), "Downstream")
}

///|
test "cs_feature_distal_intergenic" {
  let g = @src.GeneAnnotation::new(
    gene_id="G1",
    gene_name="Gene1",
    chr="chr1",
    gene_start=10000,
    gene_end=25000,
    strand="+",
  )
  // Peak at 50000 is far past the downstream window [25000, 28000].
  let p = @src.Peak::new(chr="chr1", start=49900, end=50100, summit=50000)
  assert_eq(@src.assign_feature(p, g, 3000), "Distal Intergenic")
}

///|
test "cs_feature_downstream_minus_strand" {
  let g = @src.GeneAnnotation::new(
    gene_id="G2",
    gene_name="Gene2",
    chr="chr1",
    gene_start=50000,
    gene_end=70000,
    strand="-",
  )
  // For - strand, downstream = [gene_start - 3000, gene_start] = [47000, 50000].
  let p = @src.Peak::new(chr="chr1", start=48900, end=49100, summit=49000)
  assert_eq(@src.assign_feature(p, g, 3000), "Downstream")
}

///|
test "cs_feature_5utr_minus_strand" {
  let g = @src.GeneAnnotation::new(
    gene_id="G2",
    gene_name="Gene2",
    chr="chr1",
    gene_start=50000,
    gene_end=70000,
    strand="-",
  )
  // For - strand, 5'UTR = [gene_end - 500, gene_end] = [69500, 70000].
  // Use small window 100 so peak at 69700 is in 5'UTR but not promoter.
  let p = @src.Peak::new(chr="chr1", start=69600, end=69800, summit=69700)
  assert_eq(@src.assign_feature(p, g, 100), "5'UTR")
}

///|
test "cs_find_nearest_gene" {
  let genes = @src.chipseeker_sample_genes()
  // Peak at chr1:8000 is nearest to GeneA (TSS=10000, dist=2000),
  // not GeneB (TSS=70000, dist=62000).
  let p = @src.Peak::new(chr="chr1", start=7900, end=8100, summit=8000)
  let g = @src.find_nearest_gene(p, genes)
  assert_eq(g.gene_id, "GENE1")
  assert_eq(g.gene_name, "GeneA")
}

///|
test "cs_find_nearest_gene_empty" {
  let p = @src.Peak::new(chr="chr1", start=1000, end=2000, summit=1500)
  let g = @src.find_nearest_gene(p, [])
  // Sentinel gene with empty gene_id.
  assert_eq(g.gene_id, "")
  assert_eq(g.gene_name, "")
}

///|
test "cs_find_nearest_gene_different_chr" {
  let genes = @src.chipseeker_sample_genes()
  // No gene on chr3 in the sample data.
  let p = @src.Peak::new(chr="chr3", start=10000, end=10200, summit=10100)
  let g = @src.find_nearest_gene(p, genes)
  assert_eq(g.gene_id, "")
}

///|
test "cs_annotate_single_peak" {
  let genes = @src.chipseeker_sample_genes()
  let p = @src.Peak::new(chr="chr1", start=7900, end=8100, summit=8000)
  let ann = @src.annotate_peak(p, genes, 3000)
  assert_eq(ann.chr, "chr1")
  assert_eq(ann.gene_id, "GENE1")
  assert_eq(ann.gene_name, "GeneA")
  assert_eq(ann.feature, "Promoter")
  assert_eq(ann.distance_to_tss, -2000)
  assert_eq(ann.tss, 10000)
  assert_true(ann.annotation.contains("Promoter"))
  assert_true(ann.annotation.contains("GeneA"))
}

///|
test "cs_annotate_multiple_peaks" {
  let peaks = @src.chipseeker_sample_peaks()
  let genes = @src.chipseeker_sample_genes()
  let anns = @src.annotate_peaks(peaks, genes, 3000)
  assert_eq(anns.length(), peaks.length())
  // Every annotation should have a non-empty feature string.
  for ann in anns {
    assert_true(ann.feature.length() > 0)
  }
}

///|
test "cs_annotation_summary" {
  let peaks = @src.chipseeker_sample_peaks()
  let genes = @src.chipseeker_sample_genes()
  let anns = @src.annotate_peaks(peaks, genes, 3000)
  let summary = @src.annotation_summary(anns)
  assert_eq(summary.total_peaks, peaks.length())
  // Counts should sum to total_peaks.
  let sum = summary.promoter +
    summary.utr5 +
    summary.utr3 +
    summary.exon +
    summary.intron +
    summary.downstream +
    summary.distal_intergenic
  assert_eq(sum, summary.total_peaks)
  // Sample data with default window has at least some promoter peaks.
  assert_true(summary.promoter > 0)
  // And at least one distal intergenic (the chr3 peak).
  assert_true(summary.distal_intergenic > 0)
}

///|
test "cs_filter_by_feature" {
  let peaks = @src.chipseeker_sample_peaks()
  let genes = @src.chipseeker_sample_genes()
  let anns = @src.annotate_peaks(peaks, genes, 3000)
  let promoters = @src.filter_by_feature(anns, "Promoter")
  for ann in promoters {
    assert_eq(ann.feature, "Promoter")
  }
  // Sample data has 4 promoter peaks (3 on chr1, 1 on chr2 at TSS).
  assert_true(promoters.length() >= 3)
}

///|
test "cs_filter_by_feature_intron" {
  let peaks = @src.chipseeker_sample_peaks()
  let genes = @src.chipseeker_sample_genes()
  let anns = @src.annotate_peaks(peaks, genes, 3000)
  let introns = @src.filter_by_feature(anns, "Intron")
  for ann in introns {
    assert_eq(ann.feature, "Intron")
  }
  assert_true(introns.length() >= 2)
}

///|
test "cs_filter_by_distance" {
  let peaks = @src.chipseeker_sample_peaks()
  let genes = @src.chipseeker_sample_genes()
  let anns = @src.annotate_peaks(peaks, genes, 3000)
  let close = @src.filter_by_distance(anns, 3000)
  for ann in close {
    assert_true(ann.distance_to_tss.abs() <= 3000)
  }
  // Promoter peaks (4 of them) all have |distance| <= 3000.
  assert_true(close.length() >= 4)
}

///|
test "cs_filter_by_distance_zero" {
  let peaks = @src.chipseeker_sample_peaks()
  let genes = @src.chipseeker_sample_genes()
  let anns = @src.annotate_peaks(peaks, genes, 3000)
  // Only the peak exactly at the TSS (chr2:20000) has |distance| == 0.
  let at_tss = @src.filter_by_distance(anns, 0)
  assert_true(at_tss.length() >= 1)
  for ann in at_tss {
    assert_eq(ann.distance_to_tss, 0)
  }
}

///|
test "cs_sample_peaks_validation" {
  let peaks = @src.chipseeker_sample_peaks()
  assert_true(peaks.length() >= 10)
  // All peaks should have a chromosome set.
  for p in peaks {
    assert_true(p.chr.length() > 0)
    assert_true(p.end >= p.start)
  }
}

///|
test "cs_sample_genes_validation" {
  let genes = @src.chipseeker_sample_genes()
  assert_eq(genes.length(), 3)
  // Check GeneA (+ strand).
  assert_eq(genes[0].gene_id, "GENE1")
  assert_eq(genes[0].strand, "+")
  assert_eq(genes[0].tss, genes[0].gene_start)
  // Check GeneB (- strand).
  assert_eq(genes[1].gene_id, "GENE2")
  assert_eq(genes[1].strand, "-")
  assert_eq(genes[1].tss, genes[1].gene_end)
  // Check GeneC (+ strand).
  assert_eq(genes[2].gene_id, "GENE3")
  assert_eq(genes[2].strand, "+")
  assert_eq(genes[2].tss, genes[2].gene_start)
}

///|
test "cs_promoter_window_custom" {
  let g = @src.GeneAnnotation::new(
    gene_id="G1",
    gene_name="Gene1",
    chr="chr1",
    gene_start=10000,
    gene_end=25000,
    strand="+",
  )
  // Peak 2000bp upstream of TSS.
  let p = @src.Peak::new(chr="chr1", start=7900, end=8100, summit=8000)
  // With window=1000, the peak is NOT in promoter.
  assert_eq(@src.assign_feature(p, g, 1000), "Distal Intergenic")
  // With window=3000, the peak IS in promoter.
  assert_eq(@src.assign_feature(p, g, 3000), "Promoter")
}

///|
test "cs_peak_at_tss" {
  let g = @src.GeneAnnotation::new(
    gene_id="G1",
    gene_name="Gene1",
    chr="chr1",
    gene_start=10000,
    gene_end=25000,
    strand="+",
  )
  let p = @src.Peak::new(chr="chr1", start=9900, end=10100, summit=10000)
  // Distance is 0, so always Promoter regardless of window (as long as window >= 0).
  assert_eq(@src.assign_feature(p, g, 0), "Promoter")
  assert_eq(@src.assign_feature(p, g, 3000), "Promoter")
  assert_eq(@src.compute_tss_distance(p, g), 0)
}

///|
test "cs_peak_far_away" {
  let genes = @src.chipseeker_sample_genes()
  // Peak on chr1 at 1,000,000 - far from all genes.
  let p = @src.Peak::new(chr="chr1", start=999900, end=1000100, summit=1000000)
  let g = @src.find_nearest_gene(p, genes)
  let dist = @src.compute_tss_distance(p, g)
  // Distance magnitude should be huge (peak ~930k bp from nearest TSS).
  assert_true(dist.abs() > 900000)
  assert_eq(@src.assign_feature(p, g, 3000), "Distal Intergenic")
}

///|
test "cs_peak_on_different_chr" {
  let genes = @src.chipseeker_sample_genes()
  let p = @src.Peak::new(chr="chrY", start=1000, end=2000, summit=1500)
  let g = @src.find_nearest_gene(p, genes)
  // No gene on chrY -> sentinel with empty gene_id.
  assert_eq(g.gene_id, "")
  assert_eq(@src.compute_tss_distance(p, g), 0)
  assert_eq(@src.assign_feature(p, g, 3000), "Distal Intergenic")
}

///|
test "cs_annotation_contains_details" {
  let genes = @src.chipseeker_sample_genes()
  let p = @src.Peak::new(chr="chr1", start=25900, end=26100, summit=26000)
  let ann = @src.annotate_peak(p, genes, 3000)
  assert_eq(ann.feature, "Downstream")
  // The annotation string should mention the gene name and "Downstream".
  assert_true(ann.annotation.contains("Downstream"))
  assert_true(ann.annotation.contains("GeneA"))
}

///|
test "cs_summary_string" {
  let peaks = @src.chipseeker_sample_peaks()
  let genes = @src.chipseeker_sample_genes()
  let anns = @src.annotate_peaks(peaks, genes, 3000)
  let summary = @src.annotation_summary(anns)
  let s = @src.chipseeker_summary_string(summary)
  assert_true(s.contains("ChIPseeker Annotation Summary"))
  assert_true(s.contains("Total peaks"))
  assert_true(s.contains("Promoter"))
  assert_true(s.contains("5'UTR"))
  assert_true(s.contains("3'UTR"))
  assert_true(s.contains("Exon"))
  assert_true(s.contains("Intron"))
  assert_true(s.contains("Downstream"))
  assert_true(s.contains("Distal Intergenic"))
}

///|
test "cs_peak_annotation_accessors" {
  let genes = @src.chipseeker_sample_genes()
  let p = @src.Peak::new(
    chr="chr1",
    start=7900,
    end=8100,
    summit=8000,
    fold_enrichment=5.0,
    p_value=0.001,
  )
  let ann = @src.annotate_peak(p, genes, 3000)
  assert_eq(ann.get_chr(), "chr1")
  assert_eq(ann.get_gene_id(), "GENE1")
  assert_eq(ann.get_gene_name(), "GeneA")
  assert_eq(ann.get_feature(), "Promoter")
  assert_eq(ann.get_distance_to_tss(), -2000)
  assert_eq(ann.get_tss(), 10000)
  assert_true(ann.get_annotation().contains("Promoter"))
  // The peak is preserved.
  let recovered = ann.get_peak()
  assert_eq(recovered.chr, "chr1")
  assert_eq(recovered.summit, 8000)
}

///|
test "cs_summary_accessors" {
  let summary = @src.AnnotationSummary::new(
    total_peaks=100,
    promoter=30,
    utr5=5,
    utr3=5,
    exon=20,
    intron=15,
    downstream=10,
    distal_intergenic=15,
  )
  assert_eq(summary.get_total_peaks(), 100)
  assert_eq(summary.get_promoter(), 30)
  assert_eq(summary.get_utr5(), 5)
  assert_eq(summary.get_utr3(), 5)
  assert_eq(summary.get_exon(), 20)
  assert_eq(summary.get_intron(), 15)
  assert_eq(summary.get_downstream(), 10)
  assert_eq(summary.get_distal_intergenic(), 15)
}

///|
test "cs_feature_priority_promoter_over_utr" {
  // When a peak is in BOTH the promoter window and the 5'UTR,
  // promoter should win (higher priority).
  let g = @src.GeneAnnotation::new(
    gene_id="G1",
    gene_name="Gene1",
    chr="chr1",
    gene_start=10000,
    gene_end=25000,
    strand="+",
  )
  // Peak at 10200 is in 5'UTR [10000, 10500] AND in promoter [7000, 13000].
  let p = @src.Peak::new(chr="chr1", start=10100, end=10300, summit=10200)
  assert_eq(@src.assign_feature(p, g, 3000), "Promoter")
}

///|
test "cs_annotate_empty_peaks" {
  let genes = @src.chipseeker_sample_genes()
  let anns = @src.annotate_peaks([], genes, 3000)
  assert_eq(anns.length(), 0)
  let summary = @src.annotation_summary(anns)
  assert_eq(summary.total_peaks, 0)
  assert_eq(summary.promoter, 0)
}

///|
test "cs_filter_by_feature_no_match" {
  let peaks = @src.chipseeker_sample_peaks()
  let genes = @src.chipseeker_sample_genes()
  let anns = @src.annotate_peaks(peaks, genes, 3000)
  // "Intron (non-exist)" - filtering for a feature that doesn't exist.
  let no_match = @src.filter_by_feature(anns, "Nonexistent")
  assert_eq(no_match.length(), 0)
}
