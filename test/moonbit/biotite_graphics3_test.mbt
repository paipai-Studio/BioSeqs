///|
/// Tests for biotite_graphics3.mbt — SVG graphics and color schemes.
test "biotite3_color_schemes" {
  let nucleotide = @src.biotite3_get_color_scheme("nucleotide")
  assert_eq(@src.biotite3_color_scheme_get(nucleotide, "A"), "#00FF00")
  assert_eq(@src.biotite3_color_scheme_get(nucleotide, "T"), "#DB3333")
  assert_eq(@src.biotite3_color_scheme_get(nucleotide, "X"), "#808080")
  let protein = @src.biotite3_get_color_scheme("protein")
  assert_eq(@src.biotite3_color_scheme_get(protein, "A"), "#DB3333")
  assert_eq(@src.biotite3_list_color_scheme_names().length(), 2)
}

///|
test "biotite3_sequence_logo_uses_alphabet_colors" {
  let logo = @src.biotite3_plot_sequence_logo(
    [[1.0, 0.0, 0.0, 0.0]],
    ["A", "C", "G", "T"],
    100,
  )
  assert_true(logo.has_prefix("<svg"))
  assert_true(logo.contains("#00FF00"))
  let protein_logo = @src.biotite3_plot_sequence_logo(
    [[1.0, 0.0, 0.0]],
    ["A", "C", "D"],
    100,
  )
  assert_true(protein_logo.contains("#DB3333"))
}

///|
test "biotite3_alignment_type_based_colors" {
  let scheme = @src.biotite3_get_color_scheme("nucleotide")
  let svg = @src.biotite3_plot_alignment_type_based(
    ["AC", "AG"],
    ["s1", "s2"],
    scheme,
  )
  assert_true(svg.contains("background:#0000FF"))
  assert_true(svg.contains("background:#FFFF00"))
}

///|
test "biotite3_plasmid_map_and_plotters" {
  let map = @src.biotite3_plot_plasmid_map(
    1000,
    [("gene", 0, 500, "#FF0000")],
    100,
  )
  assert_true(map.has_prefix("<svg"))
  assert_true(map.contains("gene"))
  let plotter = @src.biotite3_symbol_plotter(
    "rect",
    @src.biotite3_get_color_scheme("nucleotide"),
  )
  assert_true(
    @src.biotite3_symbol_plotter_plot(plotter, "A", 0, 0, 10).contains(
      "#00FF00",
    ),
  )
  let feature = @src.biotite3_coding_plotter()
  assert_true(
    @src.biotite3_feature_plotter_plot(feature, 0, 10, 0, 5).contains("CDS"),
  )
}

///|
test "biotite3_alignment_similarity_gradient" {
  let identical = @src.biotite3_plot_alignment_similarity_based(
    ["ACGT", "ACGT"],
    ["a", "b"],
  )
  assert_true(identical.contains("background:#006400"))
  let variable = @src.biotite3_plot_alignment_similarity_based(
    ["ACGT", "TGCA", "GGGG"],
    ["a", "b", "c"],
  )
  assert_true(variable.contains("background:#FFFFFF"))
}
