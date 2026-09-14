///|
/// tidytree - A tidy approach to tree manipulation.
/// Inspired by Bioconductor tidytree package.
/// Provides tibble-like interface for phylogenetic tree operations.

///|
/// Tree tibble row: one row per node.
pub struct TreeTibbleRow {
  node : Int // node number
  parent : Int // parent node number (-1 for root)
  branch_length : Double
  label : String
  isTip : Bool
  x : Double // x coordinate (branch length cumulative)
  mut y : Double // y coordinate (for display)
  branch : Double // branch length
  angle : Double // angle for circular layout
  annotations : Map[String, String]
}

///|
/// Tree tibble: tidy representation of a tree.
pub struct TreeTibble {
  rows : Array[TreeTibbleRow]
  n_nodes : Int
  n_tips : Int
  root : Int
}

///|
/// Convert TreeioTreeData to TreeTibble (as_tibble).
pub fn tidy_as_tibble(tree : TreeioTreeData) -> TreeTibble {
  let rows : Array[TreeTibbleRow] = []
  let mut i = 0
  while i < tree.nodes.length() {
    let node = tree.nodes[i]
    // Compute x coordinate (cumulative branch length from root)
    let x = tidy_compute_x(tree, node.node_id)
    rows.push(TreeTibbleRow::{
      node: node.node_id,
      parent: node.parent,
      branch_length: node.branch_length,
      label: node.label,
      isTip: node.is_leaf,
      x,
      y: 0.0, // will be computed later
      branch: node.branch_length,
      angle: 0.0,
      annotations: node.annotations,
    })
    i = i + 1
  }
  // Compute y coordinates for tips
  let mut tip_count = 0
  let mut j = 0
  while j < rows.length() {
    if rows[j].isTip {
      rows[j].y = tip_count.to_double()
      tip_count = tip_count + 1
    }
    j = j + 1
  }
  // Compute y for internal nodes (average of children)
  let mut k = 0
  while k < rows.length() {
    if !rows[k].isTip {
      rows[k].y = tidy_compute_internal_y(tree, rows, rows[k].node)
    }
    k = k + 1
  }
  TreeTibble::{
    rows,
    n_nodes: tree.n_nodes,
    n_tips: tree.n_tips,
    root: tree.root,
  }
}

///|
/// Compute x coordinate (cumulative branch length from root).
fn tidy_compute_x(tree : TreeioTreeData, node_id : Int) -> Double {
  let mut x = 0.0
  let mut current = node_id
  while current != tree.root && current >= 0 {
    // Find node and add branch length
    let mut i = 0
    while i < tree.nodes.length() {
      if tree.nodes[i].node_id == current {
        x = x + tree.nodes[i].branch_length
        current = tree.nodes[i].parent
      }
      i = i + 1
    }
    if current < 0 {
      break
    }
  }
  x
}

///|
/// Compute y coordinate for internal node (average of children y).
fn tidy_compute_internal_y(
  tree : TreeioTreeData,
  rows : Array[TreeTibbleRow],
  node_id : Int,
) -> Double {
  // Find children
  let mut child_y_sum = 0.0
  let mut child_count = 0
  let mut i = 0
  while i < tree.nodes.length() {
    if tree.nodes[i].node_id == node_id {
      let mut ci = 0
      while ci < tree.nodes[i].children.length() {
        let child_id = tree.nodes[i].children[ci]
        // Find child y in rows
        let mut ri = 0
        while ri < rows.length() {
          if rows[ri].node == child_id {
            child_y_sum = child_y_sum + rows[ri].y
            child_count = child_count + 1
          }
          ri = ri + 1
        }
        ci = ci + 1
      }
    }
    i = i + 1
  }
  if child_count > 0 {
    child_y_sum / child_count.to_double()
  } else {
    0.0
  }
}

///|
/// Get all offspring (descendants) of a node.
pub fn tidy_offspring(tree : TreeioTreeData, node_id : Int) -> Array[Int] {
  let result : Array[Int] = []
  tidy_offspring_rec(tree, node_id, result)
  result
}

///|
fn tidy_offspring_rec(
  tree : TreeioTreeData,
  node_id : Int,
  result : Array[Int],
) -> Unit {
  let mut i = 0
  while i < tree.nodes.length() {
    if tree.nodes[i].node_id == node_id {
      let mut ci = 0
      while ci < tree.nodes[i].children.length() {
        let child_id = tree.nodes[i].children[ci]
        result.push(child_id)
        tidy_offspring_rec(tree, child_id, result)
        ci = ci + 1
      }
    }
    i = i + 1
  }
}

///|
/// Get direct children of a node.
pub fn tidy_child(tree : TreeioTreeData, node_id : Int) -> Array[Int] {
  let mut i = 0
  while i < tree.nodes.length() {
    if tree.nodes[i].node_id == node_id {
      return tree.nodes[i].children
    }
    i = i + 1
  }
  []
}

///|
/// Get parent of a node.
pub fn tidy_parent(tree : TreeioTreeData, node_id : Int) -> Int {
  let mut i = 0
  while i < tree.nodes.length() {
    if tree.nodes[i].node_id == node_id {
      return tree.nodes[i].parent
    }
    i = i + 1
  }
  -1
}

///|
/// Get sibling of a node (other children of parent).
pub fn tidy_sibling(tree : TreeioTreeData, node_id : Int) -> Array[Int] {
  let parent = tidy_parent(tree, node_id)
  if parent < 0 {
    return []
  }
  let siblings : Array[Int] = []
  let children = tidy_child(tree, parent)
  let mut i = 0
  while i < children.length() {
    if children[i] != node_id {
      siblings.push(children[i])
    }
    i = i + 1
  }
  siblings
}

///|
/// Get most recent common ancestor (MRCA) of two nodes.
pub fn tidy_mrca(tree : TreeioTreeData, node_a : Int, node_b : Int) -> Int {
  // Get ancestors of node_a (including itself)
  let ancestors_a : Array[Int] = []
  let mut current = node_a
  while current >= 0 {
    ancestors_a.push(current)
    current = tidy_parent(tree, current)
  }
  // Walk up from node_b and find first match
  current = node_b
  while current >= 0 {
    let mut i = 0
    while i < ancestors_a.length() {
      if ancestors_a[i] == current {
        return current
      }
      i = i + 1
    }
    current = tidy_parent(tree, current)
  }
  tree.root
}

///|
/// Check if a node is a leaf (tip).
pub fn tidy_is_leaf(tree : TreeioTreeData, node_id : Int) -> Bool {
  let mut i = 0
  while i < tree.nodes.length() {
    if tree.nodes[i].node_id == node_id {
      return tree.nodes[i].is_leaf
    }
    i = i + 1
  }
  false
}

///|
/// Get node depth (distance from root in number of edges).
pub fn tidy_node_depth(tree : TreeioTreeData, node_id : Int) -> Int {
  let mut depth = 0
  let mut current = node_id
  while current != tree.root && current >= 0 {
    current = tidy_parent(tree, current)
    depth = depth + 1
  }
  depth
}

///|
/// Get all ancestors of a node (path to root).
pub fn tidy_ancestors(tree : TreeioTreeData, node_id : Int) -> Array[Int] {
  let ancestors : Array[Int] = []
  let mut current = tidy_parent(tree, node_id)
  while current >= 0 {
    ancestors.push(current)
    current = tidy_parent(tree, current)
  }
  ancestors
}

///|
/// Get all tips (leaf nodes).
pub fn tidy_tips(tree : TreeioTreeData) -> Array[Int] {
  let tips : Array[Int] = []
  let mut i = 0
  while i < tree.nodes.length() {
    if tree.nodes[i].is_leaf {
      tips.push(tree.nodes[i].node_id)
    }
    i = i + 1
  }
  tips
}

///|
/// Get all internal nodes.
pub fn tidy_internal(tree : TreeioTreeData) -> Array[Int] {
  let internal : Array[Int] = []
  let mut i = 0
  while i < tree.nodes.length() {
    if !tree.nodes[i].is_leaf {
      internal.push(tree.nodes[i].node_id)
    }
    i = i + 1
  }
  internal
}

///|
/// Filter tibble by condition (isTip).
pub fn tidy_filter_tips(tibble : TreeTibble) -> TreeTibble {
  let rows : Array[TreeTibbleRow] = []
  let mut i = 0
  while i < tibble.rows.length() {
    if tibble.rows[i].isTip {
      rows.push(tibble.rows[i])
    }
    i = i + 1
  }
  TreeTibble::{ rows, n_nodes: rows.length(), n_tips: rows.length(), root: tibble.root }
}

///|
/// Filter tibble for internal nodes only.
pub fn tidy_filter_internal(tibble : TreeTibble) -> TreeTibble {
  let rows : Array[TreeTibbleRow] = []
  let mut i = 0
  while i < tibble.rows.length() {
    if !tibble.rows[i].isTip {
      rows.push(tibble.rows[i])
    }
    i = i + 1
  }
  TreeTibble::{ rows, n_nodes: rows.length(), n_tips: 0, root: tibble.root }
}

///|
/// Select specific nodes from tibble.
pub fn tidy_select(tibble : TreeTibble, node_ids : Array[Int]) -> TreeTibble {
  let rows : Array[TreeTibbleRow] = []
  let mut i = 0
  while i < tibble.rows.length() {
    let mut found = false
    let mut j = 0
    while j < node_ids.length() {
      if tibble.rows[i].node == node_ids[j] {
        found = true
      }
      j = j + 1
    }
    if found {
      rows.push(tibble.rows[i])
    }
    i = i + 1
  }
  let n_tips = {
    let mut count = 0
    let mut k = 0
    while k < rows.length() {
      if rows[k].isTip {
        count = count + 1
      }
      k = k + 1
    }
    count
  }
  TreeTibble::{ rows, n_nodes: rows.length(), n_tips, root: tibble.root }
}

///|
/// Get node by label.
pub fn tidy_get_node_by_label(tree : TreeioTreeData, label : String) -> Int {
  let mut i = 0
  while i < tree.nodes.length() {
    if tree.nodes[i].label == label {
      return tree.nodes[i].node_id
    }
    i = i + 1
  }
  -1
}

///|
/// Get label by node id.
pub fn tidy_get_label(tree : TreeioTreeData, node_id : Int) -> String {
  let mut i = 0
  while i < tree.nodes.length() {
    if tree.nodes[i].node_id == node_id {
      return tree.nodes[i].label
    }
    i = i + 1
  }
  ""
}

///|
/// Get branch length between node and its parent.
pub fn tidy_branch_length(tree : TreeioTreeData, node_id : Int) -> Double {
  let mut i = 0
  while i < tree.nodes.length() {
    if tree.nodes[i].node_id == node_id {
      return tree.nodes[i].branch_length
    }
    i = i + 1
  }
  0.0
}

///|
/// Calculate total branch length from node to root.
pub fn tidy_dist_to_root(tree : TreeioTreeData, node_id : Int) -> Double {
  let mut dist = 0.0
  let mut current = node_id
  while current != tree.root && current >= 0 {
    dist = dist + tidy_branch_length(tree, current)
    current = tidy_parent(tree, current)
  }
  dist
}

///|
/// Calculate patristic distance between two nodes.
pub fn tidy_patristic_distance(
  tree : TreeioTreeData,
  node_a : Int,
  node_b : Int,
) -> Double {
  let mrca = tidy_mrca(tree, node_a, node_b)
  let dist_a = tidy_dist_to_root(tree, node_a) - tidy_dist_to_root(tree, mrca)
  let dist_b = tidy_dist_to_root(tree, node_b) - tidy_dist_to_root(tree, mrca)
  dist_a + dist_b
}

///|
/// Reorder tree by post-order traversal.
pub fn tidy_post_order(tree : TreeioTreeData) -> Array[Int] {
  let result : Array[Int] = []
  tidy_post_order_rec(tree, tree.root, result)
  result
}

///|
fn tidy_post_order_rec(
  tree : TreeioTreeData,
  node_id : Int,
  result : Array[Int],
) -> Unit {
  let children = tidy_child(tree, node_id)
  let mut i = 0
  while i < children.length() {
    tidy_post_order_rec(tree, children[i], result)
    i = i + 1
  }
  result.push(node_id)
}

///|
/// Reorder tree by pre-order traversal.
pub fn tidy_pre_order(tree : TreeioTreeData) -> Array[Int] {
  let result : Array[Int] = []
  tidy_pre_order_rec(tree, tree.root, result)
  result
}

///|
fn tidy_pre_order_rec(
  tree : TreeioTreeData,
  node_id : Int,
  result : Array[Int],
) -> Unit {
  result.push(node_id)
  let children = tidy_child(tree, node_id)
  let mut i = 0
  while i < children.length() {
    tidy_pre_order_rec(tree, children[i], result)
    i = i + 1
  }
}

///|
/// Print tibble as a table (for debugging/display).
pub fn tidy_print(tibble : TreeTibble) -> String {
  let sb = StringBuilder::new()
  sb.write_string("node\tparent\tbranch_length\tlabel\tisTip\tx\ty\n")
  let mut i = 0
  while i < tibble.rows.length() {
    let r = tibble.rows[i]
    sb.write_string(r.node.to_string())
    sb.write_char('\t')
    sb.write_string(r.parent.to_string())
    sb.write_char('\t')
    sb.write_string(r.branch_length.to_string())
    sb.write_char('\t')
    sb.write_string(r.label)
    sb.write_char('\t')
    sb.write_string(if r.isTip { "true" } else { "false" })
    sb.write_char('\t')
    sb.write_string(r.x.to_string())
    sb.write_char('\t')
    sb.write_string(r.y.to_string())
    sb.write_char('\n')
    i = i + 1
  }
  sb.to_string()
}

///|
/// Create example tibble for testing.
pub fn tidy_example_tibble() -> TreeTibble {
  let tree = treeio_example_tree()
  tidy_as_tibble(tree)
}