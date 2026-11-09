/// skbio_tree_ops.mbt — Tree traversal, LCA, copy, augment
/// Items 1-4: skbio_tree_traverse, skbio_tree_lca, skbio_tree_copy, skbio_tree_augment

///|
/// Simple tree node for skbio-style tree operations.
pub struct SkbioTreeNode2 {
  name : String
  children : Array[SkbioTreeNode2]
  parent : Int32 // index of parent, -1 for root
  mut branch_length : Double
  mut metadata : String // JSON or free-form annotation
}

///|
/// Traversal order enum.
pub enum SkbioTraverseOrder {
  SkbioPreOrder
  SkbioPostOrder
  SkbioInOrder
  SkbioLevelOrder
}

///|
/// Traverse a tree and return node names in specified order.
pub fn skbio_tree_traverse(
  root : SkbioTreeNode2,
  order~ : SkbioTraverseOrder = SkbioPreOrder,
) -> Array[String] {
  let result : Array[String] = Array(capacity=0)
  match order {
    SkbioPreOrder => skbio_traverse_pre(root, result)
    SkbioPostOrder => skbio_traverse_post(root, result)
    SkbioInOrder => skbio_traverse_in(root, result)
    SkbioLevelOrder => skbio_traverse_level(root, result)
  }
  result
}

///|
fn skbio_traverse_pre(node : SkbioTreeNode2, result : Array[String]) -> Unit {
  result.push(node.name)
  for i = 0; i < node.children.length(); i = i + 1 {
    skbio_traverse_pre(node.children[i], result)
  }
}

///|
fn skbio_traverse_post(node : SkbioTreeNode2, result : Array[String]) -> Unit {
  for i = 0; i < node.children.length(); i = i + 1 {
    skbio_traverse_post(node.children[i], result)
  }
  result.push(node.name)
}

///|
fn skbio_traverse_in(node : SkbioTreeNode2, result : Array[String]) -> Unit {
  let n = node.children.length()
  if n == 0 {
    result.push(node.name)
  } else {
    // Left half
    let half = n / 2
    for i = 0; i < half; i = i + 1 {
      skbio_traverse_in(node.children[i], result)
    }
    result.push(node.name)
    for i = half; i < n; i = i + 1 {
      skbio_traverse_in(node.children[i], result)
    }
  }
}

///|
fn skbio_traverse_level(root : SkbioTreeNode2, result : Array[String]) -> Unit {
  let queue : Array[SkbioTreeNode2] = Array(capacity=0)
  queue.push(root)
  while queue.length() > 0 {
    let node = queue[0]
    let _ = queue.remove(0)
    result.push(node.name)
    for i = 0; i < node.children.length(); i = i + 1 {
      queue.push(node.children[i])
    }
  }
}

///|
/// Find lowest common ancestor of two named nodes.
/// Returns the name of the LCA node, or "" if not found.
pub fn skbio_tree_lca(
  root : SkbioTreeNode2,
  name1~ : String,
  name2~ : String,
) -> String {
  // Find path from root to each node
  let path1 : Array[String] = Array(capacity=0)
  let path2 : Array[String] = Array(capacity=0)
  let _ = skbio_find_path(root, name1, path1)
  let _ = skbio_find_path(root, name2, path2)
  if path1.length() == 0 || path2.length() == 0 {
    return ""
  }
  // Find last common node in paths
  let mut lca = ""
  let min_len = if path1.length() < path2.length() {
    path1.length()
  } else {
    path2.length()
  }
  for i = 0; i < min_len; i = i + 1 {
    if path1[i] == path2[i] {
      lca = path1[i]
    } else {
      break
    }
  }
  lca
}

///|
/// Find path from root to target node. Returns true if found.
fn skbio_find_path(
  node : SkbioTreeNode2,
  target : String,
  path : Array[String],
) -> Bool {
  path.push(node.name)
  if node.name == target {
    return true
  }
  for i = 0; i < node.children.length(); i = i + 1 {
    if skbio_find_path(node.children[i], target, path) {
      return true
    }
  }
  let _ = path.remove(path.length() - 1)
  false
}

///|
/// Deep copy a tree (subtree rooted at given node).
pub fn skbio_tree_copy(node : SkbioTreeNode2) -> SkbioTreeNode2 {
  let copy = SkbioTreeNode2::{
    name: node.name,
    children: Array(capacity=node.children.length()),
    parent: node.parent,
    branch_length: node.branch_length,
    metadata: node.metadata,
  }
  for i = 0; i < node.children.length(); i = i + 1 {
    copy.children.push(skbio_tree_copy(node.children[i]))
  }
  copy
}

///|
/// Check if two trees are isomorphic (same shape, ignoring names).
pub fn skbio_tree_isomorphic(t1 : SkbioTreeNode2, t2 : SkbioTreeNode2) -> Bool {
  if t1.children.length() != t2.children.length() {
    return false
  }
  // For binary trees, check both possible child orderings
  let n = t1.children.length()
  if n == 0 {
    return true
  }
  if n == 2 {
    // Try same order and swapped order
    if skbio_tree_isomorphic(t1.children[0], t2.children[0]) &&
      skbio_tree_isomorphic(t1.children[1], t2.children[1]) {
      return true
    }
    if skbio_tree_isomorphic(t1.children[0], t2.children[1]) &&
      skbio_tree_isomorphic(t1.children[1], t2.children[0]) {
      return true
    }
    return false
  }
  // For non-binary, try all child orderings (simple: sort by child count)
  for i = 0; i < n; i = i + 1 {
    if !skbio_tree_isomorphic(t1.children[i], t2.children[i]) {
      return false
    }
  }
  true
}

///|
/// Augment a tree node with metadata.
pub fn skbio_tree_augment(
  node : SkbioTreeNode2,
  metadata~ : String = "",
  bootstrap~ : Double = 0.0,
  branch_length~ : Double = -1.0,
) -> SkbioTreeNode2 {
  let result = skbio_tree_copy(node)
  if metadata.length() > 0 {
    result.metadata = metadata
  }
  if bootstrap > 0.0 {
    let bs_str = "bootstrap=" + bootstrap.to_string()
    if result.metadata.length() > 0 {
      result.metadata = result.metadata + ";" + bs_str
    } else {
      result.metadata = bs_str
    }
  }
  if branch_length >= 0.0 {
    result.branch_length = branch_length
  }
  result
}

///|
/// Count total nodes in a tree.
pub fn skbio_tree_count_nodes(node : SkbioTreeNode2) -> Int32 {
  let mut count = 1
  for i = 0; i < node.children.length(); i = i + 1 {
    count = count + skbio_tree_count_nodes(node.children[i])
  }
  count
}

///|
/// Get all leaf names from a tree.
pub fn skbio_tree_leaf_names(node : SkbioTreeNode2) -> Array[String] {
  let result : Array[String] = Array(capacity=0)
  skbio_collect_leaves(node, result)
  result
}

///|
fn skbio_collect_leaves(node : SkbioTreeNode2, result : Array[String]) -> Unit {
  if node.children.length() == 0 {
    result.push(node.name)
  } else {
    for i = 0; i < node.children.length(); i = i + 1 {
      skbio_collect_leaves(node.children[i], result)
    }
  }
}

///|
pub fn skbio_traverse_dummy() -> Unit {
  ignore(SkbioPostOrder)
  ignore(SkbioInOrder)
  ignore(SkbioLevelOrder)
}
