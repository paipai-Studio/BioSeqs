///|
test "ggtree_create_test_tree" {
  let nodes = @src.create_test_tree()
  assert_eq(nodes.length() > 0, true)

  // Check root node
  let root = nodes[0]
  assert_eq(root.id, "root")
  assert_eq(root.label, "Root")
  assert_eq(root.is_leaf, false)

  // Check leaf nodes exist by counting them recursively
  let leaf_ids = @src.ggtree_collect_leaf_ids(root)
  assert_eq(leaf_ids.length() > 0, true)
}

///|
test "ggtree_tree_data_new" {
  let tree_data = @src.TreeData::new("(A:1,B:0.5):0.8;")
  assert_eq(tree_data.newick, "(A:1,B:0.5):0.8;")
}

///|
test "ggtree_count_leaves" {
  let nodes = @src.create_test_tree()
  let root = nodes[0]
  let n_leaves = @src.ggtree_count_leaves(root)
  assert_eq(n_leaves > 0, true)
}

///|
test "ggtree_y_positions" {
  let nodes = @src.create_test_tree()
  let root = nodes[0]
  let n_nodes = @src.ggtree_count_leaves(root) + @src.ggtree_collect_internal_ids(root).length()
  let positions : Map[String, (Double, Double, Double)] = Map([], capacity=n_nodes)
  let _ = @src.ggtree_compute_y_positions(root, 0.0, positions)

  // Root should have a y coordinate
  let root_pos = positions[root.id]
  assert_eq(root_pos.1 >= 0.0, true)

  // Leaf nodes should have sequential y positions
  let leaf_ids = @src.ggtree_collect_leaf_ids(root)
  for id in leaf_ids {
    let pos = positions[id]
    assert_eq(pos.1 >= 0.0, true)
  }
}

///|
test "ggtree_x_positions" {
  let nodes = @src.create_test_tree()
  let root = nodes[0]
  let n_nodes = @src.ggtree_count_leaves(root) + @src.ggtree_collect_internal_ids(root).length()
  let positions : Map[String, (Double, Double, Double)] = Map([], capacity=n_nodes)
  let _ = @src.ggtree_compute_y_positions(root, 0.0, positions)
  @src.ggtree_compute_x_positions(root, 0.0, positions)

  let root_pos = positions[root.id]
  assert_eq(root_pos.0 >= 0.0, true)
}

///|
test "ggtree_collect_leaf_ids" {
  let nodes = @src.create_test_tree()
  let leaf_order = @src.ggtree_collect_leaf_ids(nodes[0])
  assert_eq(leaf_order.length() > 0, true)
}

///|
test "ggtree_collect_internal_ids" {
  let nodes = @src.create_test_tree()
  let node_order = @src.ggtree_collect_internal_ids(nodes[0])
  assert_eq(node_order.length() > 0, true)
}

///|
test "ggtree_layout_rectangular" {
  let nodes = @src.create_test_tree()
  let tree_data = @src.create_tree_data_from_nodes(nodes)

  let layout = @src.layout_rectangular(tree_data)
  assert_eq(layout.layout_type, "rectangular")
  assert_eq(layout.nodes.length() > 0, true)
  assert_eq(layout.positions.length() > 0, true)
  assert_eq(layout.max_x >= layout.min_x, true)
  assert_eq(layout.max_y >= layout.min_y, true)
}

///|
test "ggtree_layout_rectangular_empty" {
  let tree_data = @src.TreeData::new("")
  let layout = @src.layout_rectangular(tree_data)
  assert_eq(layout.layout_type, "rectangular")
  assert_eq(layout.nodes.length(), 0)
}

///|
test "ggtree_layout_radial" {
  let nodes = @src.create_test_tree()
  let tree_data = @src.create_tree_data_from_nodes(nodes)

  let layout = @src.layout_radial(tree_data)
  assert_eq(layout.layout_type, "radial")
  assert_eq(layout.nodes.length() > 0, true)
}

///|
test "ggtree_layout_radial_empty" {
  let tree_data = @src.TreeData::new("")
  let layout = @src.layout_radial(tree_data)
  assert_eq(layout.layout_type, "radial")
  assert_eq(layout.nodes.length(), 0)
}

///|
test "ggtree_layout_unrooted" {
  let nodes = @src.create_test_tree()
  let tree_data = @src.create_tree_data_from_nodes(nodes)

  let layout = @src.layout_unrooted(tree_data)
  assert_eq(layout.layout_type, "unrooted")
  assert_eq(layout.nodes.length() > 0, true)
}

///|
test "ggtree_get_position" {
  let nodes = @src.create_test_tree()
  let tree_data = @src.create_tree_data_from_nodes(nodes)

  let layout = @src.layout_rectangular(tree_data)
  let pos = @src.ggtree_get_position(layout, "root")
  assert_eq(pos.0 >= 0.0, true)
  assert_eq(pos.1 >= 0.0, true)
}

///|
test "ggtree_generate_tree_viz_data" {
  let nodes = @src.create_test_tree()
  let tree_data = @src.create_tree_data_from_nodes(nodes)

  let layout = @src.layout_rectangular(tree_data)
  let viz_data = @src.generate_tree_viz_data(layout)
  assert_eq(viz_data.length() > 0, true)
}

///|
test "ggtree_generate_edge_data" {
  let nodes = @src.create_test_tree()
  let tree_data = @src.create_tree_data_from_nodes(nodes)

  let layout = @src.layout_rectangular(tree_data)
  let edges = @src.generate_edge_data(layout)
  assert_eq(edges.length() > 0, true)
}

///|
test "ggtree_generate_label_data" {
  let nodes = @src.create_test_tree()
  let tree_data = @src.create_tree_data_from_nodes(nodes)

  let layout = @src.layout_rectangular(tree_data)
  let labels = @src.generate_label_data(layout)
  assert_eq(labels.length() > 0, true)
}

///|
test "ggtree_layout_node" {
  let node = @src.LayoutNode::new("test_id", "TestLabel", 1.5, true)
  assert_eq(node.id, "test_id")
  assert_eq(node.label, "TestLabel")
  assert_eq(node.branch_length, 1.5)
  assert_eq(node.is_leaf, true)
  assert_eq(node.children.length(), 0)
}

///|
test "ggtree_tree_layout_new" {
  let nodes = @src.create_test_tree()
  let positions : Map[String, (Double, Double, Double)] = Map([], capacity=3)
  positions.set("root", (0.0, 0.0, 0.0))
  positions.set("A", (1.0, 1.0, 0.0))
  positions.set("B", (2.0, 2.0, 0.0))

  let layout = @src.TreeLayout::new(
    nodes,
    positions,
    ["B", "A"],
    ["root"],
    10.0,
    5.0,
    0.0,
    0.0,
    "rectangular",
  )
  assert_eq(layout.layout_type, "rectangular")
  assert_eq(layout.max_x, 10.0)
  assert_eq(layout.max_y, 5.0)
}
