///|
/// ChIPseeker peak annotation demo.
///
/// Loads sample ChIP-seq peaks and gene annotations, annotates each peak
/// with its nearest gene and genomic feature, prints the per-peak
/// annotations, summarizes the feature distribution, and demonstrates
/// filtering by feature and by distance to TSS.

fn main {
  println("=== ChIPseeker Peak Annotation Demo ===")
  println("Inspired by the Bioconductor ChIPseeker package")
  println("")

  // 1. Load sample data.
  println("1. Loading sample data...")
  let peaks = @bio.chipseeker_sample_peaks()
  let genes = @bio.chipseeker_sample_genes()
  println("   Peaks loaded: \{peaks.length()}")
  println("   Genes loaded: \{genes.length()}")
  println("")

  // 2. Annotate peaks. Use the default promoter window of 3000bp.
  println("2. Annotating peaks (promoter window = 3000bp)...")
  let annotations = @bio.annotate_peaks(peaks, genes, 3000)
  println("   Annotated \{annotations.length()} peaks")
  println("")

  // 3. Display the first few per-peak annotations.
  println("3. Per-peak annotations (first 6):")
  println("   -----------------------------------------------")
  let mut i = 0
  while i < annotations.length() && i < 6 {
    let ann = annotations[i]
    let p = ann.get_peak()
    println("   Peak \{i + 1}: \{p.chr}:\{p.start}-\{p.end} (summit=\{p.summit})")
    println("     Gene      : \{ann.gene_name} (\{ann.gene_id})")
    println("     TSS       : \{ann.tss}")
    println("     Distance  : \{ann.distance_to_tss} (signed, + = downstream)")
    println("     Feature   : \{ann.feature}")
    println("     Annotation: \{ann.annotation}")
    println("")
    i = i + 1
  }

  // 4. Summary statistics.
  println("4. Annotation summary:")
  println("   --------------------")
  let summary = @bio.annotation_summary(annotations)
  println(@bio.chipseeker_summary_string(summary))

  // 5. Filter by feature: show promoter peaks.
  println("5. Filtering peaks by feature: 'Promoter'")
  let promoters = @bio.filter_by_feature(annotations, "Promoter")
  println("   Promoter peaks: \{promoters.length()}")
  for ann in promoters {
    let p = ann.get_peak()
    println("     - \{p.chr}:\{p.summit} -> \{ann.gene_name} (dist=\{ann.distance_to_tss})")
  }
  println("")

  // 6. Filter by feature: show exon peaks.
  println("6. Filtering peaks by feature: 'Exon'")
  let exons = @bio.filter_by_feature(annotations, "Exon")
  println("   Exon peaks: \{exons.length()}")
  for ann in exons {
    let p = ann.get_peak()
    println("     - \{p.chr}:\{p.summit} -> \{ann.gene_name} (dist=\{ann.distance_to_tss})")
  }
  println("")

  // 7. Filter by distance: peaks within 3000bp of a TSS.
  println("7. Filtering peaks by distance: |distance| <= 3000")
  let close = @bio.filter_by_distance(annotations, 3000)
  println("   Close peaks: \{close.length()}")
  for ann in close {
    let p = ann.get_peak()
    println("     - \{p.chr}:\{p.summit} -> \{ann.gene_name} (dist=\{ann.distance_to_tss}, \{ann.feature})")
  }
  println("")

  // 8. Demonstrate promoter window customization. With a smaller window
  //    (100bp), the 5'UTR of GeneA becomes accessible because the 5'UTR
  //    no longer sits entirely inside the promoter window.
  println("8. Promoter window customization (window=100bp):")
  println("   Annotating a peak inside the 5'UTR of GeneA (summit=10300)...")
  let peak_in_5utr = @bio.Peak::new(chr="chr1", start=10200, end=10400, summit=10300)
  let ann_small_window = @bio.annotate_peak(peak_in_5utr, genes, 100)
  println("   Feature with window=100 : \{ann_small_window.feature}")
  let ann_large_window = @bio.annotate_peak(peak_in_5utr, genes, 3000)
  println("   Feature with window=3000: \{ann_large_window.feature}")
  println("   (smaller window exposes the 5'UTR assignment)")
  println("")

  // 9. Edge case: peak on a chromosome with no genes.
  println("9. Edge case: peak on a chromosome with no genes (chrZ)")
  let orphan_peak = @bio.Peak::new(chr="chrZ", start=1000, end=2000, summit=1500)
  let orphan_ann = @bio.annotate_peak(orphan_peak, genes, 3000)
  println("   Gene id   : '\{orphan_ann.gene_id}' (empty -> no nearby gene)")
  println("   Feature   : \{orphan_ann.feature}")
  println("   Annotation: \{orphan_ann.annotation}")
  println("")

  println("=== ChIPseeker demo complete ===")
}
