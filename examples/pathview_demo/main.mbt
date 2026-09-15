///|
/// Demo for pathview module: KEGG pathway visualization.

fn main {
  println("=== pathview Demo: KEGG Pathway Visualization ===")
  println("")

  // Example 1: Build a simple pathway
  println("--- Step 1: Build Pathway ---")
  let pw0 = @src.pathview_new("hsa00010")
  let pw1 = @src.pathview_add_node(pw0, @src.PathviewNode::{ node_id: "n1", label: "HK1", node_type: "gene", x: 20.0, y: 50.0, width: 50.0, height: 22.0 })
  let pw2 = @src.pathview_add_node(pw1, @src.PathviewNode::{ node_id: "n2", label: "GPI", node_type: "gene", x: 100.0, y: 50.0, width: 50.0, height: 22.0 })
  let pw3 = @src.pathview_add_node(pw2, @src.PathviewNode::{ node_id: "n3", label: "PFKP", node_type: "gene", x: 180.0, y: 50.0, width: 50.0, height: 22.0 })
  let pw4 = @src.pathview_add_edge(pw3, @src.PathviewEdge::{ source: "n1", target: "n2", edge_type: "PPrel" })
  let pw5 = @src.pathview_add_edge(pw4, @src.PathviewEdge::{ source: "n2", target: "n3", edge_type: "PPrel" })
  let pw6 = @src.pathview_map_gene(pw5, "HK1", "n1")
  let pw7 = @src.pathview_map_gene(pw6, "GPI", "n2")
  let pw = @src.pathview_map_gene(pw7, "PFKP", "n3")
  println("Pathway: " + pw.pathway_id + ", nodes=" + pw.nodes.length().to_string() + ", edges=" + pw.edges.length().to_string())
  println("")

  // Example 2: Color mapping for logFC values
  println("--- Step 2: logFC Color Mapping ---")
  let up_color = @src.pathview_logfc_color(3.0)
  let down_color = @src.pathview_logfc_color(-3.0)
  let zero_color = @src.pathview_logfc_color(0.0)
  println("logFC=+3.0 -> " + @src.pathview_color_hex(up_color) + " (rgb=" + @src.pathview_color_rgb(up_color) + ")")
  println("logFC= 0.0 -> " + @src.pathview_color_hex(zero_color))
  println("logFC=-3.0 -> " + @src.pathview_color_hex(down_color) + " (rgb=" + @src.pathview_color_rgb(down_color) + ")")
  println("")

  // Example 3: Node statistics
  println("--- Step 3: Node Statistics ---")
  let vals = [("HK1", 2.5), ("GPI", -1.0), ("PFKP", 0.0)]
  let stats = @src.pathview_node_stats(pw, vals)
  println("Node stats: " + stats.length().to_string() + " nodes")
  for s in stats {
    println("  " + s.node_id + " " + "mean=" + s.mean_value.to_string() + " up=" + s.num_up.to_string() + " down=" + s.num_down.to_string() + " genes=" + s.num_genes.to_string())
  }
  println("")

  // Example 4: Render SVG
  println("--- Step 4: Render SVG ---")
  let svg = @src.pathview_render_svg(pw, vals)
  println("SVG length: " + svg.length().to_string() + " chars")
  println("Contains <svg>: " + svg.contains("<svg").to_string())
  println("Contains gene HK1: " + svg.contains("HK1").to_string())
  println("")

  // Example 5: Gene tally
  println("--- Step 5: Gene Tally ---")
  let (up_count, down_count) = @src.pathview_gene_tally(vals)
  println("Up-regulated: " + up_count.to_string() + ", Down-regulated: " + down_count.to_string())
  println("")

  // Example 6: Legend
  println("--- Step 6: Color Legend ---")
  let legend = @src.pathview_legend()
  println("Legend SVG length: " + legend.length().to_string() + " chars")
  println("")

  // Example 7: Find node
  println("--- Step 7: Find Node ---")
  match @src.pathview_find_node(pw, "n1") {
    Some(n) => println("Found node n1: " + n.label + " at (" + n.x.to_string() + ", " + n.y.to_string() + ")")
    None => println("Node n1 not found")
  }
  println("")

  println("=== pathview Demo Complete ===")
}
