///|
/// Tests for the Biopython Bio.Graphics.KGML_vis-compatible renderer.

///|
fn load_sample_pathway() -> @src.KgmlPathway {
  match @src.parse_kgml(@src.kgml_sample_pathway()) {
    Some(p) => p
    None => abort("failed to parse sample KGML")
  }
}

///|
test "kgml canvas renders svg document" {
  let pathway = load_sample_pathway()
  let svg = @src.KgmlCanvas::new(pathway).to_svg(width=600, height=400)
  assert_true(svg.contains("<svg"))
  assert_true(svg.contains("</svg>"))
  assert_true(svg.contains("Glycolysis / Gluconeogenesis"))
}

///|
test "kgml canvas draws gene rectangles" {
  let pathway = load_sample_pathway()
  let svg = @src.KgmlCanvas::new(pathway).to_svg()
  // Gene entries use bgcolor #BFFFBF
  assert_true(svg.contains("#BFFFBF"))
  assert_true(svg.contains("<rect"))
}

///|
test "kgml canvas draws compound circles" {
  let pathway = load_sample_pathway()
  let svg = @src.KgmlCanvas::new(pathway).to_svg()
  assert_true(svg.contains("<circle"))
}

///|
test "kgml canvas draws map with dashed border" {
  let pathway = load_sample_pathway()
  let svg = @src.KgmlCanvas::new(pathway).to_svg()
  assert_true(svg.contains("stroke-dasharray=\"4,2\""))
  assert_true(svg.contains("TCA cycle"))
}

///|
test "kgml canvas draws relation lines" {
  let pathway = load_sample_pathway()
  let svg = @src.KgmlCanvas::new(pathway).to_svg()
  // ECrel color is #1f77b4
  assert_true(svg.contains("#1f77b4"))
  // maplink color is #7f7f7f
  assert_true(svg.contains("#7f7f7f"))
  assert_true(svg.contains("<line"))
}

///|
test "kgml canvas draws reaction arrows for complete pathway" {
  // Custom pathway where both substrate and product compounds have entries,
  // so the reaction arrow can be rendered.
  let xml =
    "<?xml version=\"1.0\"?>\n" +
    "<pathway name=\"path:ko00010\" org=\"ko\" number=\"00010\" title=\"Demo\">\n" +
    "<entry id=\"1\" name=\"cpd:C00001\" type=\"compound\">" +
    "<graphics name=\"C00001\" fgcolor=\"#000000\" bgcolor=\"#FFFFFF\" type=\"circle\" x=\"10\" y=\"10\" width=\"8\" height=\"8\"/>" +
    "</entry>\n" +
    "<entry id=\"2\" name=\"cpd:C00002\" type=\"compound\">" +
    "<graphics name=\"C00002\" fgcolor=\"#000000\" bgcolor=\"#FFFFFF\" type=\"circle\" x=\"100\" y=\"10\" width=\"8\" height=\"8\"/>" +
    "</entry>\n" +
    "<reaction name=\"rn:R00001\" type=\"irreversible\">" +
    "<substrate name=\"cpd:C00001\"/>" +
    "<product name=\"cpd:C00002\"/>" +
    "</reaction>\n" +
    "</pathway>"
  let svg = @src.kgml_to_svg(xml, width=200, height=100)
  match svg {
    Some(s) => assert_true(s.contains("<polygon"))
    None => assert_true(false)
  }
}

///|
test "kgml canvas options toggle labels" {
  let pathway = load_sample_pathway()
  let opts = @src.KgmlCanvasOptions::{
    ..@src.kgml_default_options(),
    show_labels: false,
  }
  let svg = @src.KgmlCanvas::new(pathway, options=opts).to_svg()
  // Entry graphics name "K00844" should not appear when labels are off.
  assert_true(!(svg.contains(">K00844<")))
}

///|
test "kgml canvas options toggle relations" {
  let pathway = load_sample_pathway()
  let opts = @src.KgmlCanvasOptions::{
    ..@src.kgml_default_options(),
    show_relations: false,
  }
  let svg = @src.KgmlCanvas::new(pathway, options=opts).to_svg()
  // Without relations the ECrel color should not appear.
  assert_true(!(svg.contains("#1f77b4")))
}

///|
test "kgml to svg convenience" {
  let svg = @src.kgml_to_svg(@src.kgml_sample_pathway(), width=400, height=300)
  assert_true(svg is Some(_))
}

///|
test "kgml entry summary counts by type" {
  let pathway = load_sample_pathway()
  let summary = @src.kgml_entry_summary(pathway)
  // 3 gene + 1 compound + 1 map
  let mut total = 0
  for kv in summary {
    total += kv.1
  }
  assert_eq(total, 5)
}

///|
test "kgml canvas default options relation colors" {
  let opts = @src.kgml_default_options()
  assert_eq(opts.relation_colors.get("ECrel"), Some("#1f77b4"))
  assert_eq(opts.relation_colors.get("PPrel"), Some("#2ca02c"))
  assert_eq(opts.relation_colors.get("maplink"), Some("#7f7f7f"))
}