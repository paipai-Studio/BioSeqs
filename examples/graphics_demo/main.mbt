///|
/// Graphics Demo - Demonstrates bioinformatics visualization.

fn main {
  println("=== Graphics Demo ===")
  
  println("")
  println("=== Sequence Logo ===")
  let seqs = [
    "ATGCCGTAA",
    "ATGCCGTAA",
    "ATGCCGTAA",
    "ATGCCGTAA",
    "ATGCCGTAA",
    "ATGCCGTAA",
    "ATGCCGTAA",
    "ATGCCGTAA",
    "ATGCCGTAA",
    "ATGCCGTAA"
  ]
  let logo = @bio.SeqLogo::new(seqs)
  let logo_output = @bio.seqlogo_generate_ascii(logo)
  println(logo_output)
  
  println("")
  println("=== Alignment Plot ===")
  let alignment_seqs = [
    ("Human", "ATGGCCGGTAAATGA"),
    ("Mouse", "ATGGCCGGTAAATGA"),
    ("Rat", "ATGGCCGGTAATAGA"),
    ("Chicken", "ATGGCCGGTAATGGA")
  ]
  let plot = @bio.AlignmentPlot::new(alignment_seqs)
  let alignment_output = @bio.alignment_plot_generate_ascii(plot)
  println(alignment_output)
  
  println("")
  println("=== Feature Plot ===")
  let features = [
    ("GeneA", 0, 100, "gene", "+"),
    ("Exon1", 0, 30, "exon", "+"),
    ("Intron1", 30, 70, "intron", "+"),
    ("Exon2", 70, 100, "exon", "+"),
    ("Promoter", -20, 0, "promoter", "+")
  ]
  let feat_plot = @bio.FeaturePlot::new(features, 150)
  let feature_output = @bio.feature_plot_generate_ascii(feat_plot)
  println(feature_output)
}