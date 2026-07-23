///|
/// Tests for Graphics module.

test "SeqLogo creation" {
  let seqs = ["ATGC"]
  let logo = @src.SeqLogo::new(seqs)
  assert_eq(logo.sequences.length(), 1)
}

test "LogoColumn creation" {
  let col = @src.LogoColumn::new(1)
  assert_eq(col.position, 1)
}

test "AlignmentPlot creation" {
  let seqs = [("Seq1", "ATGC")]
  let plot = @src.AlignmentPlot::new(seqs)
  assert_eq(plot.sequences.length(), 1)
}

test "FeaturePlot creation" {
  let features = [("Gene", 0, 10, "gene", "+")]
  let plot = @src.FeaturePlot::new(features, 100)
  assert_eq(plot.features.length(), 1)
}

test "seqlogo_calculate_columns" {
  let seqs = ["ATGC", "ATGC"]
  let columns = @src.seqlogo_calculate_columns(seqs)
  assert_eq(columns.length(), 4)
}

test "seqlogo_generate_ascii" {
  let seqs = ["ATGC", "ATGC", "ATGC"]
  let logo = @src.SeqLogo::new(seqs)
  let output = @src.seqlogo_generate_ascii(logo)
  assert_true(output.length() > 0)
}

test "alignment_plot_generate_ascii" {
  let seqs = [("Human", "ATGC"), ("Mouse", "ATGC")]
  let plot = @src.AlignmentPlot::new(seqs)
  let output = @src.alignment_plot_generate_ascii(plot)
  assert_true(output.length() > 0)
}

test "feature_plot_generate_ascii" {
  let features = [("Gene", 0, 10, "gene", "+")]
  let plot = @src.FeaturePlot::new(features, 50)
  let output = @src.feature_plot_generate_ascii(plot)
  assert_true(output.length() > 0)
}

test "create_example_seqlogo" {
  let logo = @src.create_example_seqlogo()
  assert_eq(logo.sequences.length(), 10)
}

test "create_example_alignment_plot" {
  let plot = @src.create_example_alignment_plot()
  assert_eq(plot.sequences.length(), 4)
}

test "create_example_feature_plot" {
  let plot = @src.create_example_feature_plot()
  assert_eq(plot.features.length(), 7)
}

test "get_default_colors" {
  let colors = @src.get_default_colors()
  assert_true(colors.contains("A"))
}

test "get_feature_colors" {
  let colors = @src.get_feature_colors()
  assert_true(colors.contains("exon"))
}