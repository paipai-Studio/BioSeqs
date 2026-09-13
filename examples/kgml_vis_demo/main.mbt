///|
fn main {
  println("=== Bio.Graphics.KGML_vis Demo ===")

  let pathway = match @src.parse_kgml(@src.kgml_sample_pathway()) {
    Some(p) => p
    None => abort("failed to parse sample KGML")
  }

  println("\n1. Pathway metadata")
  println("  Name:  " + pathway.name())
  println("  Org:   " + pathway.org())
  println("  Title: " + pathway.title())

  println("\n2. Entry summary")
  let summary = @src.kgml_entry_summary(pathway)
  for kv in summary {
    println("  " + kv.0 + ": " + kv.1.to_string())
  }

  println("\n3. Render to SVG (default options)")
  let canvas = @src.KgmlCanvas::new(pathway)
  let svg = canvas.to_svg(width=600, height=400)
  println("  SVG length: " + svg.length().to_string() + " chars")
  println("  Contains <svg:    " + svg.contains("<svg").to_string())
  println("  Contains <rect:    " + svg.contains("<rect").to_string())
  println("  Contains <circle:  " + svg.contains("<circle").to_string())
  println("  Contains <line:    " + svg.contains("<line").to_string())
  println("  Contains <polygon: " + svg.contains("<polygon").to_string())
  println("  Contains <text:    " + svg.contains("<text").to_string())

  println("\n4. Render without labels")
  let opts = @src.KgmlCanvasOptions::{
    ..@src.kgml_default_options(),
    show_labels: false,
  }
  let svg_no_labels = @src.KgmlCanvas::new(pathway, options=opts).to_svg()
  println("  SVG length: " + svg_no_labels.length().to_string() + " chars")
  println("  Has 'K00844' text: " + svg_no_labels.contains(">K00844<").to_string())

  println("\n5. Render without relations or reactions")
  let opts_min = @src.KgmlCanvasOptions::{
    ..@src.kgml_default_options(),
    show_relations: false,
    show_reactions: false,
  }
  let svg_min = @src.KgmlCanvas::new(pathway, options=opts_min).to_svg()
  println("  SVG length: " + svg_min.length().to_string() + " chars")

  println("\n6. Convenience: kgml_to_svg from raw XML")
  match @src.kgml_to_svg(@src.kgml_sample_pathway(), width=400, height=300) {
    Some(s) => println("  Rendered " + s.length().to_string() + " chars")
    None => println("  Parse failed")
  }
}