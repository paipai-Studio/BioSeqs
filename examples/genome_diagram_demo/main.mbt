///|
/// Genome Diagram Demo - Demonstrates genome diagram creation and SVG generation.

fn main {
  println("=== Genome Diagram Demo ===")

  let diagram = @bio.gd_create_diagram("Example Genome", 0, 5000)

  println("Created diagram: name=" + diagram.name() +
    ", length=" + @bio.gd_diagram_length(diagram).to_string())

  let diagram = @bio.gd_add_track(diagram, "Genes")
  let diagram = @bio.gd_add_track(diagram, "Exons")
  let diagram = @bio.gd_add_track(diagram, "Regulatory")
  println("Track count: " + @bio.gd_track_count(diagram).to_string())

  let diagram = @bio.gd_add_feature(
    diagram, 100, 600, "+", "BRCA1", @bio.gd_arrow_shape(),
  )
  let diagram = @bio.gd_add_feature(
    diagram, 800, 1400, "-", "TP53", @bio.gd_arrow_shape(),
  )
  let diagram = @bio.gd_add_feature(
    diagram, 1600, 1800, "+", "Promoter1", @bio.gd_diamond_shape(),
  )
  let diagram = @bio.gd_add_feature(
    diagram, 2000, 2800, "+", "EGFR", @bio.gd_arrow_shape(),
  )
  let diagram = @bio.gd_add_feature(
    diagram, 3000, 3200, "-", "Enhancer1", @bio.gd_terminators_shape(),
  )
  let diagram = @bio.gd_add_feature(
    diagram, 3500, 4200, "+", "MYC", @bio.gd_crossed_arrow_shape(),
  )
  let diagram = @bio.gd_add_feature(
    diagram, 4400, 4900, "-", "KRAS", @bio.gd_arrow_shape(),
  )
  println("Diagram feature count: " + @bio.gd_feature_count(diagram).to_string())

  let diagram = @bio.gd_add_track_feature(
    diagram, "Genes", 100, 600, "+", "BRCA1_gene",
  )
  let diagram = @bio.gd_add_track_feature(
    diagram, "Genes", 800, 1400, "-", "TP53_gene",
  )
  let diagram = @bio.gd_add_track_feature(
    diagram, "Genes", 2000, 2800, "+", "EGFR_gene",
  )
  let diagram = @bio.gd_add_track_feature(
    diagram, "Exons", 150, 300, "+", "Exon1",
  )
  let diagram = @bio.gd_add_track_feature(
    diagram, "Exons", 350, 500, "+", "Exon2",
  )
  let diagram = @bio.gd_add_track_feature(
    diagram, "Exons", 850, 1000, "-", "Exon3",
  )
  let diagram = @bio.gd_add_track_feature(
    diagram, "Exons", 1100, 1300, "-", "Exon4",
  )
  let diagram = @bio.gd_add_track_feature(
    diagram, "Regulatory", 1600, 1800, "+", "Promoter1",
  )
  let diagram = @bio.gd_add_track_feature(
    diagram, "Regulatory", 3000, 3200, "-", "Enhancer1",
  )
  println("Total feature count (including tracks): " +
    @bio.gd_feature_count(diagram).to_string())

  let overlapping = @bio.gd_find_overlapping_features(diagram, 500, 900)
  println("Features overlapping [500, 900]: " + overlapping.length().to_string())
  for feat in overlapping {
    println("  - " + feat.label() + " [" +
      feat.start().to_string() + ", " +
      feat.end().to_string() + "] " +
      feat.strand())
  }

  let style = @bio.DiagramStyle::new(
    scale=2.0,
    color_scheme="default",
  )
  let diagram = @bio.gd_set_style(diagram, style)

  let diagram = @bio.gd_label_features(diagram, 200)

  let svg = @bio.gd_to_svg_string(diagram, 900, 500)
  println("SVG output length: " + svg.length().to_string())
  println("SVG preview (first 200 chars):")
  if svg.length() > 200 {
    println(svg[0:200].to_owned())
  } else {
    println(svg)
  }

  println("")
  println("=== Diagram Statistics ===")
  println("Name: " + diagram.name())
  println("Range: [" + diagram.start().to_string() +
    ", " + diagram.end().to_string() + "]")
  println("Length: " + @bio.gd_diagram_length(diagram).to_string())
  println("Tracks: " + @bio.gd_track_count(diagram).to_string())
  println("Total features: " + @bio.gd_feature_count(diagram).to_string())

  let gene_features = @bio.gd_get_features(diagram, "Genes")
  println("Gene track features: " + gene_features.length().to_string())
  for feat in gene_features {
    println("  - " + feat.label() + " [" +
      feat.start().to_string() + ", " +
      feat.end().to_string() + "]")
  }

  println("")
  println("Demo completed successfully!")
}