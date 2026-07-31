///|
/// Gviz module example - Genome visualization track system.
///
/// Run: moon run examples/gviz_demo
fn main {
  println("=== Gviz Module Demo ===")
  println("Genome Visualization Track System")
  println("")

  // Create a plot
  println("=== Creating Genome Plot ===")
  let region = @bio.gviz_region("chr1", 1000, 5000)
  let plot = @bio.gviz_plot(region, title="BRCA1 Locus", width=80, height=25)
  println("Region: chr1:1000-5000")
  println("")

  // Add axis track
  println("=== Adding Axis Track ===")
  let axis = @bio.gviz_track("axis", @bio.track_type_genome_axis(), "chr1", 1000, 5000)
  axis.set_label("Position")
  plot.add_track(axis)
  println("Added GenomeAxisTrack")
  println("")

  // Add gene annotation track
  println("=== Adding Gene Annotation Track ===")
  let genes = @bio.gviz_track("genes", @bio.track_type_gene_region(), "chr1", 1000, 5000)
  genes.set_label("Gene Annotations")
  genes.set_color("blue")
  genes.add_feature(@bio.gviz_feature("g1", "chr1", 1200, 1800, @bio.track_strand_forward(), "exon", "BRCA1"))
  genes.add_feature(@bio.gviz_feature("g2", "chr1", 2000, 3500, @bio.track_strand_reverse(), "exon", "TP53"))
  genes.add_feature(@bio.gviz_feature("g3", "chr1", 4000, 4500, @bio.track_strand_forward(), "exon", "EGFR"))
  plot.add_track(genes)
  println("Added GeneRegionTrack with 3 features")
  println("")

  // Add data track
  println("=== Adding Coverage Data Track ===")
  let coverage = @bio.gviz_track("coverage", @bio.track_type_data(), "chr1", 1000, 5000)
  coverage.set_label("Read Coverage")
  coverage.set_color("red")
  let mut i = 0
  while i < 40 {
    let pos = 1000 + i * 100
    let val = 10.0 + @math.sin(i.to_double() * 0.5) * 5.0 + i.to_double() * 0.3
    coverage.add_data_point(pos, val)
    i = i + 1
  }
  plot.add_track(coverage)
  println("Added DataTrack with \{coverage.get_n_data_points()} data points")
  println("")

  // Add ideogram track
  println("=== Adding Ideogram Track ===")
  let ideo = @bio.gviz_track("ideogram", @bio.gviz_track_type_ideogram(), "chr1", 1000, 5000)
  ideo.set_label("Cytoband")
  plot.add_track(ideo)
  println("Added IdeogramTrack")
  println("")

  // Print summary
  println("=== Plot Summary ===")
  println(plot.summary())
  println("")

  // Render as ASCII
  println("=== ASCII Visualization ===")
  println(plot.to_ascii())
  println("")

  // Query features in a region
  println("=== Query Features in Region 1500-3000 ===")
  let features_in_region = genes.get_features_in_region(1500, 3000)
  println("Features found: \{features_in_region.length()}")
  let mut j = 0
  while j < features_in_region.length() {
    let f = features_in_region[j]
    println("  \{f.label}: \{f.chromosome}:\{f.start}-\{f.end_} (\{f.strand.to_string()})")
    j = j + 1
  }
  println("")

  // Use sample plot
  println("=== Sample Plot ===")
  let sample_plot = @bio.gviz_sample_plot()
  println(sample_plot.summary())
}
