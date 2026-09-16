///|
/// Demo for conumee module: copy number unified interface.

fn main {
  println("=== conumee Demo: Copy Number Analysis ===")
  println("")

  // Example 1: Create probe data
  println("--- Step 1: Probe Data ---")
  let probes = [
    @src.CnmProbe::new("chr1", 1000, 0.1),
    @src.CnmProbe::new("chr1", 2000, 0.2),
    @src.CnmProbe::new("chr1", 3000, 0.15),
    @src.CnmProbe::new("chr1", 4000, 1.5),
    @src.CnmProbe::new("chr1", 5000, 1.8),
    @src.CnmProbe::new("chr1", 6000, 1.6),
    @src.CnmProbe::new("chr1", 7000, -0.5),
    @src.CnmProbe::new("chr1", 8000, -0.8),
    @src.CnmProbe::new("chr2", 1000, 0.3),
    @src.CnmProbe::new("chr2", 2000, 0.4),
  ]
  let sample = @src.CnmSample::new("sample1", probes)
  println("Sample: " + sample.sample_id + ", probes=" + sample.probes.length().to_string())
  println("")

  // Example 2: CBS segmentation
  println("--- Step 2: CBS Segmentation ---")
  let segged = @src.cnm_segment(sample, alpha=0.05, nperm=1000)
  println("Segments found: " + segged.segments.length().to_string())
  for i in 0..<segged.segments.length() {
    let seg = segged.segments[i]
    println("  " + seg.chr + ":" + seg.start.to_string() + "-" + seg.end.to_string() +
            " CN=" + seg.copy_number.to_string() + " (n=" + seg.n_probes.to_string() + ")")
  }
  println("")

  // Example 3: Quality control
  println("--- Step 3: Quality Control ---")
  let qc = segged.qc
  println("Probes: " + qc.n_probes.to_string())
  println("Segments: " + qc.n_segments.to_string())
  println("Median logR: " + qc.median_logratio.to_string())
  println("MAD: " + qc.mad_logratio.to_string())
  println("Noise score: " + qc.noise_score.to_string())
  println("")

  // Example 4: Gene annotation
  println("--- Step 4: Gene Annotation ---")
  let genes = [
    ("chr1", 4000, 5000, "GENE1"),
    ("chr1", 7500, 8500, "GENE2"),
    ("chr2", 1000, 1500, "GENE3"),
  ]
  let annotated = @src.cnm_annotate_genes(segged, genes)
  for i in 0..<annotated.segments.length() {
    let seg = annotated.segments[i]
    if seg.genes.length() > 0 {
      println("  Segment " + seg.chr + ":" + seg.start.to_string() + " -> " +
              seg.genes.length().to_string() + " genes")
    }
  }
  println("")

  // Example 5: Filter segments
  println("--- Step 5: Filter CNVs ---")
  let gains = @src.cnm_filter_segments(annotated, gain_threshold=2.5, loss_threshold=1.5)
  println("CNV segments: " + gains.length().to_string())
  for i in 0..<gains.length() {
    let seg = gains[i]
    let label = if seg.copy_number >= 2.5 { "gain" } else { "loss" }
    println("  " + label + ": " + seg.chr + ":" + seg.start.to_string() + "-" + seg.end.to_string() +
            " CN=" + seg.copy_number.to_string())
  }
  println("")

  // Example 6: TSV export
  println("--- Step 6: TSV Export ---")
  let tsv = @src.cnm_segments_tsv(annotated.segments)
  println("TSV (first 200 chars):")
  let end = if tsv.length() > 200 { 200 } else { tsv.length() }
  println(tsv.substring(start=0, end=end))
  println("")

  // Example 7: ASCII plot
  println("--- Step 7: ASCII Plot ---")
  let plot = @src.cnm_ascii_plot(annotated.segments, "chr1")
  println(plot)
  println("")

  // Example 8: Compare two samples
  println("--- Step 8: Sample Comparison ---")
  let probes2 = [
    @src.CnmProbe::new("chr1", 1000, 0.2),
    @src.CnmProbe::new("chr1", 2000, 0.3),
    @src.CnmProbe::new("chr1", 3000, 0.25),
    @src.CnmProbe::new("chr1", 4000, 1.2),
    @src.CnmProbe::new("chr1", 5000, 1.5),
  ]
  let sample2 = @src.CnmSample::new("sample2", probes2)
  let seg2 = @src.cnm_segment(sample2)
  let cmp = @src.cnm_compare_samples(segged, seg2)
  println("Comparison: " + cmp.sample_a + " vs " + cmp.sample_b)
  println("  Shared segments: " + cmp.shared_segments.length().to_string())
  println("  Unique to " + cmp.sample_a + ": " + cmp.unique_a.length().to_string())
  println("  Unique to " + cmp.sample_b + ": " + cmp.unique_b.length().to_string())
  println("  Correlation: " + cmp.correlation.to_string())
  println("")

  println("=== conumee Demo Complete ===")
}
