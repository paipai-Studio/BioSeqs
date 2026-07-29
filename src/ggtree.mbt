///|
/// ggtree - Phylogenetic tree visualization layout.
/// Inspired by Bioconductor ggtree package.
/// Converts tree structures to visualization coordinates.

///|
/// Layout-specific tree node with immutable structure.
pub struct LayoutNode {
  pub id : String
  pub label : String
  pub branch_length : Double
  pub is_leaf : Bool
  pub children : Array[LayoutNode]
}

///|
pub fn LayoutNode::new(
  id : String,
  label : String,
  branch_length : Double,
  is_leaf : Bool,
) -> LayoutNode {
  LayoutNode::{
    id,
    label,
    branch_length,
    is_leaf,
    children: [],
  }
}

///|
pub struct TreeLayout {
  pub nodes : Array[LayoutNode]
  pub positions : Map[String, (Double, Double, Double)]
  pub leaf_order : Array[String]
  pub node_order : Array[String]
  pub max_x : Double
  pub max_y : Double
  pub min_x : Double
  pub min_y : Double
  pub layout_type : String
}

///|
pub fn TreeLayout::new(
  nodes : Array[LayoutNode],
  positions : Map[String, (Double, Double, Double)],
  leaf_order : Array[String],
  node_order : Array[String],
  max_x : Double,
  max_y : Double,
  min_x : Double,
  min_y : Double,
  layout_type : String,
) -> TreeLayout {
  TreeLayout::{
    nodes,
    positions,
    leaf_order,
    node_order,
    max_x,
    max_y,
    min_x,
    min_y,
    layout_type,
  }
}

///|
pub struct TreeData {
  pub newick : String
  pub root_id : String
  pub nodes : Array[LayoutNode]
}

///|
pub fn TreeData::new(newick : String) -> TreeData {
  let empty_nodes : Array[LayoutNode] = []
  TreeData::{ newick, root_id: "", nodes: empty_nodes }
}

///|
/// Create TreeData from existing nodes.
pub fn create_tree_data_from_nodes(nodes : Array[LayoutNode]) -> TreeData {
  let root_id = if nodes.length() > 0 { nodes[0].id } else { "" }
  TreeData::{ newick: "", root_id, nodes }
}

///|
/// Create a simple test tree structure manually.
pub fn create_test_tree() -> Array[LayoutNode] {
  let root = LayoutNode::new("root", "Root", 0.0, false)
  let node_a = LayoutNode::new("A", "A", 1.0, false)
  let node_b = LayoutNode::new("B", "B", 0.5, true)
  let node_c = LayoutNode::new("C", "C", 0.8, true)
  let node_d = LayoutNode::new("D", "D", 1.2, false)
  let node_e = LayoutNode::new("E", "E", 0.3, true)
  let node_f = LayoutNode::new("F", "F", 0.7, true)
  let node_g = LayoutNode::new("G", "G", 0.9, true)

  let children_d : Array[LayoutNode] = [node_e, node_f]
  let d_with_children = LayoutNode::{
    id: node_d.id,
    label: node_d.label,
    branch_length: node_d.branch_length,
    is_leaf: node_d.is_leaf,
    children: children_d,
  }

  let children_a : Array[LayoutNode] = [node_b, node_c]
  let a_with_children = LayoutNode::{
    id: node_a.id,
    label: node_a.label,
    branch_length: node_a.branch_length,
    is_leaf: node_a.is_leaf,
    children: children_a,
  }

  let children_root : Array[LayoutNode] = [a_with_children, d_with_children, node_g]
  let root_with_children = LayoutNode::{
    id: root.id,
    label: root.label,
    branch_length: root.branch_length,
    is_leaf: root.is_leaf,
    children: children_root,
  }

  [root_with_children]
}

///|
/// Compute the number of leaves under each node.
pub fn ggtree_count_leaves(node : LayoutNode) -> Int {
  if node.children.length() == 0 {
    1
  } else {
    let mut count = 0
    let mut i = 0
    while i < node.children.length() {
      count = count + ggtree_count_leaves(node.children[i])
      i = i + 1
    }
    count
  }
}

///|
/// Assign y positions and collect positions map.
/// Returns the next available y position.
pub fn ggtree_compute_y_positions(
  node : LayoutNode,
  y_pos : Double,
  positions : Map[String, (Double, Double, Double)],
) -> Double {
  if node.children.length() == 0 {
    // Leaf node
    positions.set(node.id, (0.0, y_pos, 0.0))
    y_pos + 1.0
  } else {
    // Internal node: process children first
    let mut current_y = y_pos
    let mut i = 0
    while i < node.children.length() {
      current_y = ggtree_compute_y_positions(node.children[i], current_y, positions)
      i = i + 1
    }
    // Center the internal node between its children
    let first_pos = positions[node.children[0].id]
    let last_pos = positions[node.children[node.children.length() - 1].id]
    let center_y = (first_pos.1 + last_pos.1) / 2.0
    positions.set(node.id, (0.0, center_y, 0.0))
    current_y
  }
}

///|
/// Assign x positions based on branch lengths.
pub fn ggtree_compute_x_positions(
  node : LayoutNode,
  x_offset : Double,
  positions : Map[String, (Double, Double, Double)],
) -> Unit {
  let y_pos = positions[node.id].1
  positions.set(node.id, (x_offset + node.branch_length, y_pos, 0.0))

  if node.children.length() > 0 {
    let node_x = x_offset + node.branch_length
    let mut i = 0
    while i < node.children.length() {
      ggtree_compute_x_positions(node.children[i], node_x, positions)
      i = i + 1
    }
  }
}

///|
/// Collect leaf node IDs in order.
pub fn ggtree_collect_leaf_ids(node : LayoutNode) -> Array[String] {
  let result : Array[String] = []
  ggtree_collect_leaf_ids_recursive(node, result)
  result
}

///|
fn ggtree_collect_leaf_ids_recursive(node : LayoutNode, result : Array[String]) -> Unit {
  if node.children.length() == 0 {
    result.push(node.id)
  } else {
    let mut i = 0
    while i < node.children.length() {
      ggtree_collect_leaf_ids_recursive(node.children[i], result)
      i = i + 1
    }
  }
}

///|
/// Collect internal node IDs in post-order.
pub fn ggtree_collect_internal_ids(node : LayoutNode) -> Array[String] {
  let result : Array[String] = []
  ggtree_collect_internal_ids_recursive(node, result)
  result
}

///|
fn ggtree_collect_internal_ids_recursive(node : LayoutNode, result : Array[String]) -> Unit {
  let mut i = 0
  while i < node.children.length() {
    ggtree_collect_internal_ids_recursive(node.children[i], result)
    i = i + 1
  }
  result.push(node.id)
}

///|
/// Compute bounding box from positions map.
fn ggtree_compute_bounds(
  nodes : Array[LayoutNode],
  positions : Map[String, (Double, Double, Double)],
) -> (Double, Double, Double, Double) {
  if nodes.length() == 0 {
    return (0.0, 0.0, 0.0, 0.0)
  }

  let mut min_x = positions[nodes[0].id].0
  let mut max_x = positions[nodes[0].id].0
  let mut min_y = positions[nodes[0].id].1
  let mut max_y = positions[nodes[0].id].1

  for node in nodes {
    if positions.contains(node.id) {
      let pos = positions[node.id]
      if pos.0 < min_x {
        min_x = pos.0
      }
      if pos.0 > max_x {
        max_x = pos.0
      }
      if pos.1 < min_y {
        min_y = pos.1
      }
      if pos.1 > max_y {
        max_y = pos.1
      }
    }
    // Recursively check children
    let mut i = 0
    while i < node.children.length() {
      let (child_min_x, child_max_x, child_min_y, child_max_y) = ggtree_compute_bounds([node.children[i]], positions)
      if child_min_x < min_x { min_x = child_min_x }
      if child_max_x > max_x { max_x = child_max_x }
      if child_min_y < min_y { min_y = child_min_y }
      if child_max_y > max_y { max_y = child_max_y }
      i = i + 1
    }
  }

  (min_x, max_x, min_y, max_y)
}

///|
/// Collect all node ids in the tree (for bounds computation).
fn ggtree_collect_all_ids(node : LayoutNode) -> Array[String] {
  let result : Array[String] = [node.id]
  if node.children.length() > 0 {
    let mut i = 0
    while i < node.children.length() {
      let child_ids = ggtree_collect_all_ids(node.children[i])
      for id in child_ids {
        result.push(id)
      }
      i = i + 1
    }
  }
  result
}

///|
/// Compute bounds from a flat list of positions.
fn ggtree_compute_bounds_flat(
  positions : Map[String, (Double, Double, Double)],
  all_ids : Array[String],
) -> (Double, Double, Double, Double) {
  if all_ids.length() == 0 {
    return (0.0, 0.0, 0.0, 0.0)
  }

  let first = positions[all_ids[0]]
  let mut min_x = first.0
  let mut max_x = first.0
  let mut min_y = first.1
  let mut max_y = first.1

  let mut i = 1
  while i < all_ids.length() {
    let pos = positions[all_ids[i]]
    if pos.0 < min_x { min_x = pos.0 }
    if pos.0 > max_x { max_x = pos.0 }
    if pos.1 < min_y { min_y = pos.1 }
    if pos.1 > max_y { max_y = pos.1 }
    i = i + 1
  }

  (min_x, max_x, min_y, max_y)
}

///|
/// Perform a rectangular (phylogram) layout on a tree.
pub fn layout_rectangular(tree_data : TreeData) -> TreeLayout {
  if tree_data.nodes.length() == 0 {
    return TreeLayout::{
      nodes: [],
      positions: Map([], capacity=0),
      leaf_order: [],
      node_order: [],
      max_x: 0.0,
      max_y: 0.0,
      min_x: 0.0,
      min_y: 0.0,
      layout_type: "rectangular",
    }
  }

  let root = tree_data.nodes[0]
  let n_nodes = ggtree_count_leaves(root) + (ggtree_collect_internal_ids(root).length())

  // Initialize positions map
  let positions : Map[String, (Double, Double, Double)] = Map([], capacity=n_nodes)

  // Step 1: Assign y positions
  let _ = ggtree_compute_y_positions(root, 0.0, positions)

  // Step 2: Assign x positions
  ggtree_compute_x_positions(root, 0.0, positions)

  // Step 3: Collect orders
  let leaf_order = ggtree_collect_leaf_ids(root)
  let node_order = ggtree_collect_internal_ids(root)

  // Step 4: Compute bounds
  let all_ids = ggtree_collect_all_ids(root)
  let (min_x, max_x, min_y, max_y) = ggtree_compute_bounds_flat(positions, all_ids)

  TreeLayout::{
    nodes: tree_data.nodes.copy(),
    positions,
    leaf_order,
    node_order,
    max_x,
    max_y,
    min_x,
    min_y,
    layout_type: "rectangular",
  }
}

///|
/// Perform a radial (circular) layout on a tree.
pub fn layout_radial(tree_data : TreeData) -> TreeLayout {
  if tree_data.nodes.length() == 0 {
    return TreeLayout::{
      nodes: [],
      positions: Map([], capacity=0),
      leaf_order: [],
      node_order: [],
      max_x: 0.0,
      max_y: 0.0,
      min_x: 0.0,
      min_y: 0.0,
      layout_type: "radial",
    }
  }

  // First compute rectangular layout
  let rect_layout = layout_rectangular(tree_data)

  // Convert to radial
  let radial_positions : Map[String, (Double, Double, Double)] = Map([], capacity=rect_layout.positions.length())
  let y_range = rect_layout.max_y - rect_layout.min_y
  let angle_scale = if y_range > 0.0 { 360.0 / y_range } else { 1.0 }

  for id in rect_layout.positions.keys() {
    let pos = rect_layout.positions[id]
    let x = pos.0
    let y = pos.1
    let angle = (y - rect_layout.min_y) * angle_scale - 90.0
    let angle_rad = angle * 3.141592653589793 / 180.0
    let new_x = x * @math.cos(angle_rad)
    let new_y = x * @math.sin(angle_rad)
    radial_positions.set(id, (new_x, new_y, angle))
  }

  let all_ids = ggtree_collect_all_ids(tree_data.nodes[0])
  let (min_x, max_x, min_y, max_y) = ggtree_compute_bounds_flat(radial_positions, all_ids)

  TreeLayout::{
    nodes: tree_data.nodes.copy(),
    positions: radial_positions,
    leaf_order: rect_layout.leaf_order,
    node_order: rect_layout.node_order,
    max_x,
    max_y,
    min_x,
    min_y,
    layout_type: "radial",
  }
}

///|
/// Perform an unrooted layout on a tree.
pub fn layout_unrooted(tree_data : TreeData) -> TreeLayout {
  if tree_data.nodes.length() == 0 {
    return TreeLayout::{
      nodes: [],
      positions: Map([], capacity=0),
      leaf_order: [],
      node_order: [],
      max_x: 0.0,
      max_y: 0.0,
      min_x: 0.0,
      min_y: 0.0,
      layout_type: "unrooted",
    }
  }

  // Compute rectangular layout
  let rect_layout = layout_rectangular(tree_data)

  // Shift to center at origin
  let center_x = (rect_layout.min_x + rect_layout.max_x) / 2.0
  let center_y = (rect_layout.min_y + rect_layout.max_y) / 2.0

  let shifted_positions : Map[String, (Double, Double, Double)] = Map([], capacity=rect_layout.positions.length())
  for id in rect_layout.positions.keys() {
    let pos = rect_layout.positions[id]
    shifted_positions.set(id, (pos.0 - center_x, pos.1 - center_y, pos.2))
  }

  let all_ids = ggtree_collect_all_ids(tree_data.nodes[0])
  let (min_x, max_x, min_y, max_y) = ggtree_compute_bounds_flat(shifted_positions, all_ids)

  TreeLayout::{
    nodes: tree_data.nodes.copy(),
    positions: shifted_positions,
    leaf_order: rect_layout.leaf_order,
    node_order: rect_layout.node_order,
    max_x,
    max_y,
    min_x,
    min_y,
    layout_type: "unrooted",
  }
}

///|
/// Get position for a node from layout.
pub fn ggtree_get_position(layout : TreeLayout, node_id : String) -> (Double, Double, Double) {
  if layout.positions.contains(node_id) {
    layout.positions[node_id]
  } else {
    (0.0, 0.0, 0.0)
  }
}

///|
/// Generate visualization data string for a tree layout.
pub fn generate_tree_viz_data(layout : TreeLayout) -> String {
  let result = "Tree Layout: \{layout.layout_type}\nBounding Box: x=[\{layout.min_x}, \{layout.max_x}], y=[\{layout.min_y}, \{layout.max_y}]\nNumber of leaf nodes: \{layout.leaf_order.length()}\nNumber of internal nodes: \{layout.node_order.length()}\n"
  result
}

///|
/// Generate edge data for tree visualization.
pub fn generate_edge_data(layout : TreeLayout) -> Array[String] {
  let edges : Array[String] = []
  if layout.nodes.length() > 0 {
    ggtree_generate_edges_recursive(layout.nodes[0], layout.positions, edges)
  }
  edges
}

///|
fn ggtree_generate_edges_recursive(
  node : LayoutNode,
  positions : Map[String, (Double, Double, Double)],
  edges : Array[String],
) -> Unit {
  if node.children.length() > 0 {
    let node_pos = if positions.contains(node.id) { positions[node.id] } else { (0.0, 0.0, 0.0) }
    let mut i = 0
    while i < node.children.length() {
      let child = node.children[i]
      let child_pos = if positions.contains(child.id) { positions[child.id] } else { (0.0, 0.0, 0.0) }
      let edge = "(\{node.id} -> \{child.id}: [\{node_pos.0}, \{node_pos.1}] -> [\{child_pos.0}, \{child_pos.1}])"
      edges.push(edge)
      ggtree_generate_edges_recursive(child, positions, edges)
      i = i + 1
    }
  }
}

///|
/// Generate label data for tree visualization.
pub fn generate_label_data(layout : TreeLayout) -> Array[String] {
  let labels : Array[String] = []
  if layout.nodes.length() > 0 {
    ggtree_generate_labels_recursive(layout.nodes[0], layout.positions, labels)
  }
  labels
}

///|
fn ggtree_generate_labels_recursive(
  node : LayoutNode,
  positions : Map[String, (Double, Double, Double)],
  labels : Array[String],
) -> Unit {
  if node.is_leaf {
    let pos = if positions.contains(node.id) { positions[node.id] } else { (0.0, 0.0, 0.0) }
    let label = "(\{node.id}, '\{node.label}', \{pos.0}, \{pos.1})"
    labels.push(label)
  }
  let mut i = 0
  while i < node.children.length() {
    ggtree_generate_labels_recursive(node.children[i], positions, labels)
    i = i + 1
  }
}
