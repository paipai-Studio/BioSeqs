///|
/// ggtree example: phylogenetic tree visualization layouts.
///
/// ggtree provides layout algorithms for phylogenetic tree visualization,
/// including rectangular (phylogram), radial, and unrooted layouts.
///
/// Run: moon run examples/ggtree_demo/main.mbt
fn main {
  // === Create a tree structure ===
  println("=== ggtree Phylogenetic Tree Visualization ===")
  let nodes = @bio.create_test_tree()
  let tree_data = @bio.create_tree_data_from_nodes(nodes)

  println("  Tree nodes: \{tree_data.nodes.length()}")

  // === Rectangular layout (phylogram) ===
  println("\n=== Rectangular Layout (Phylogram) ===")
  let rect_layout = @bio.layout_rectangular(tree_data)
  println(@bio.generate_tree_viz_data(rect_layout))

  // Display node positions
  println("  Node Positions:")
  for id in rect_layout.positions.keys() {
    let pos = rect_layout.positions[id]
    println("    \{id}: x=\{pos.0.to_string()}, y=\{pos.1.to_string()}")
  }

  // Leaf order
  println("\n  Leaf Order:")
  for i = 0; i < rect_layout.leaf_order.length(); i = i + 1 {
    println("    \{i}: \{rect_layout.leaf_order[i]}")
  }

  // Edge data
  println("\n  Edge Data:")
  let edges = @bio.generate_edge_data(rect_layout)
  for edge in edges {
    println("    \{edge}")
  }

  // Label data
  println("\n  Label Data:")
  let labels = @bio.generate_label_data(rect_layout)
  for label in labels {
    println("    \{label}")
  }

  // === Radial layout ===
  println("\n=== Radial Layout ===")
  let radial_layout = @bio.layout_radial(tree_data)
  println(@bio.generate_tree_viz_data(radial_layout))

  println("  Node Positions (radial):")
  for id in radial_layout.positions.keys() {
    let pos = radial_layout.positions[id]
    println("    \{id}: x=\{pos.0.to_string()}, y=\{pos.1.to_string()}, angle=\{pos.2.to_string()}")
  }

  // === Unrooted layout ===
  println("\n=== Unrooted Layout ===")
  let unrooted_layout = @bio.layout_unrooted(tree_data)
  println(@bio.generate_tree_viz_data(unrooted_layout))

  println("  Node Positions (unrooted):")
  for id in unrooted_layout.positions.keys() {
    let pos = unrooted_layout.positions[id]
    println("    \{id}: x=\{pos.0.to_string()}, y=\{pos.1.to_string()}")
  }

  // === Compare layouts ===
  println("\n=== Layout Comparison ===")
  println("  Rectangular: bounding box x=[\{rect_layout.min_x.to_string()}, \{rect_layout.max_x.to_string()}], y=[\{rect_layout.min_y.to_string()}, \{rect_layout.max_y.to_string()}]")
  println("  Radial:      bounding box x=[\{radial_layout.min_x.to_string()}, \{radial_layout.max_x.to_string()}], y=[\{radial_layout.min_y.to_string()}, \{radial_layout.max_y.to_string()}]")
  println("  Unrooted:    bounding box x=[\{unrooted_layout.min_x.to_string()}, \{unrooted_layout.max_x.to_string()}], y=[\{unrooted_layout.min_y.to_string()}, \{unrooted_layout.max_y.to_string()}]")

  // === Get specific node position ===
  println("\n=== Position Lookup ===")
  let a_pos = @bio.ggtree_get_position(rect_layout, "A")
  println("  Node 'A' (rectangular): x=\{a_pos.0.to_string()}, y=\{a_pos.1.to_string()}")

  let g_pos = @bio.ggtree_get_position(rect_layout, "G")
  println("  Node 'G' (rectangular): x=\{g_pos.0.to_string()}, y=\{g_pos.1.to_string()}")

  println("\nDone! Tree visualization layouts computed successfully.")
}
