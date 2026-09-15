///|
/// Tests for pathview module (KEGG pathway visualization).

fn pw_example() -> @src.PathviewPathway {
  let pw0 = @src.pathview_new("hsa00010")
  let pw1 = @src.pathview_add_node(pw0, @src.PathviewNode::{ node_id: "n1", label: "GENE1", node_type: "gene", x: 10.0, y: 10.0, width: 40.0, height: 20.0 })
  let pw2 = @src.pathview_add_node(pw1, @src.PathviewNode::{ node_id: "n2", label: "GENE2", node_type: "gene", x: 100.0, y: 10.0, width: 40.0, height: 20.0 })
  let pw3 = @src.pathview_add_edge(pw2, @src.PathviewEdge::{ source: "n1", target: "n2", edge_type: "PPrel" })
  let pw4 = @src.pathview_map_gene(pw3, "GENE1", "n1")
  @src.pathview_map_gene(pw4, "GENE2", "n2")
}

test "pathview_logfc_color_direction" {
  let up = @src.pathview_logfc_color(2.0)
  let down = @src.pathview_logfc_color(-2.0)
  let zero = @src.pathview_logfc_color(0.0)
  // up should be reddish (red > green)
  assert_true(up.red > up.green)
  // down should be greenish (green > red)
  assert_true(down.green > down.red)
  // zero should be gray (r == g == b)
  assert_eq(zero.red, zero.green)
  assert_eq(zero.green, zero.blue)
}

test "pathview_color_hex_format" {
  let c = @src.PathviewColor::{ red: 255, green: 0, blue: 0 }
  let hex = @src.pathview_color_hex(c)
  assert_eq(hex, "#FF0000")
}

test "pathview_color_rgb_format" {
  let c = @src.PathviewColor::{ red: 255, green: 0, blue: 0 }
  let rgb = @src.pathview_color_rgb(c)
  assert_eq(rgb, "rgb(255,0,0)")
}

test "pathview_node_stats" {
  let pw = pw_example()
  let vals = [("GENE1", 2.0), ("GENE2", -1.0)]
  let stats = @src.pathview_node_stats(pw, vals)
  assert_eq(stats.length(), 2)
  let mut found_up = false
  let mut found_down = false
  for s in stats {
    if s.mean_value > 0.0 { found_up = true }
    if s.mean_value < 0.0 { found_down = true }
  }
  assert_true(found_up)
  assert_true(found_down)
}

test "pathview_render_svg" {
  let pw = pw_example()
  let vals = [("GENE1", 2.5), ("GENE2", -2.5)]
  let svg = @src.pathview_render_svg(pw, vals)
  assert_true(svg.contains("<svg"))
  assert_true(svg.contains("</svg>"))
  assert_true(svg.contains("GENE1"))
}

test "pathview_gene_tally" {
  let vals = [("A", 2.0), ("B", -1.0), ("C", 3.0), ("D", 0.0)]
  let (up, down) = @src.pathview_gene_tally(vals)
  assert_eq(up, 2)
  assert_eq(down, 1)
}

test "pathview_legend" {
  let legend = @src.pathview_legend()
  assert_true(legend.contains("<svg"))
  assert_true(legend.contains("</svg>"))
}

test "pathview_find_node" {
  let pw = pw_example()
  match @src.pathview_find_node(pw, "n1") {
    Some(n) => assert_eq(n.label, "GENE1")
    None => assert_true(false)
  }
}
