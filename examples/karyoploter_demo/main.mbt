///|
/// karyoploteR module example - Karyotype visualization for genomic data.
///
/// Run: moon run examples/karyoploter_demo
fn main {
  println("=== karyoploteR Module Demo ===")
  println("")

  // Create a karyotype plot
  println("=== Creating Karyotype Plot ===")
  let plot = @bio.KaryotypePlot::new("hg38")
  println("Genome: \{plot.get_genome()}")
  println("Chromosomes: \{plot.get_chromosomes().length()}")
  println("")

  // Chromosome sizes
  println("=== Chromosome Sizes ===")
  let chroms = plot.get_chromosomes()
  let mut i = 0
  while i < 5 {
    let chr = chroms[i]
    let size = plot.get_chromosome_size(chr)
    println("  \{chr}: @{size} bp")
    i = i + 1
  }
  println("")

  // Add a data track
  println("=== Adding Data Tracks ===")
  let mut gc_track = @bio.KaryotypeTrack::new("gc", @bio.track_type_lines(), "chr1", 0, 248956422, "GC Content")
  gc_track = gc_track.set_color("#FF6600")

  let mut j = 0
  while j < 10 {
    let pos = j * 24895642 // evenly spaced positions on chr1
    let value = 0.35 + 0.2 * (j.to_double() / 10.0)
    gc_track = gc_track.add_point(@bio.TrackPoint::new("chr1", pos, value, ""))
    j = j + 1
  }
  let plot2 = plot.add_track(gc_track)
  println("Track 1 added: GC Content")

  // Add a second track
  let mut snp_track = @bio.KaryotypeTrack::new("snp", @bio.track_type_bars(), "chr1", 0, 248956422, "SNP Density")
  snp_track = snp_track.set_color("#003399")
  snp_track = snp_track.set_y_range(0.0, 1.0)

  let mut k = 0
  while k < 8 {
    let pos = k * 31119553
    let value = 0.2 + 0.5 * (k.to_double() / 8.0)
    snp_track = snp_track.add_point(@bio.TrackPoint::new("chr1", pos, value, ""))
    k = k + 1
  }
  let plot3 = plot2.add_track(snp_track)
  println("Track 2 added: SNP Density")

  // Add regions of interest
  println("")
  println("=== Adding Regions of Interest ===")
  let region1 = @bio.KaryotypeRegion::new("chr1", 11000000, 12000000, "BRCA2 region", "#FF0000")
  let region2 = @bio.KaryotypeRegion::new("chr1", 28000000, 30000000, "TP53 region", "#00FF00")
  let plot4 = plot3.add_region(region1).add_region(region2)
  println("Region 1: BRCA2 (chr1:11M-12M)")
  println("Region 2: TP53 (chr1:28M-30M)")
  println("")

  // Filter tracks by chromosome
  println("=== Chromosome Filtering ===")
  let chr1_tracks = plot4.tracks // just check what we have
  println("Total tracks: \{plot4.get_n_tracks()}")
  println("Total regions: \{plot4.get_n_regions()}")
  println("")

  // Generate ASCII visualization for chr1
  println("=== Chromosome 1 Visualization ===")
  let ascii = plot4.to_ascii("chr1", 70)
  println(ascii)

  // Summary
  println("=== Plot Summary ===")
  println(plot4.summary())

  // Sample
  println("=== Sample Karyotype ===")
  let sample = @bio.karyotype_sample()
  println(sample.summary())

  println("=== Done ===")
}