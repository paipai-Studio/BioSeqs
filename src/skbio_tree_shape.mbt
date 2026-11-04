/// Tree shape metrics: balance and imbalance indices.
///
/// scikit-bio equivalent: tree shape statistics
///
/// These metrics quantify the shape/balance of a phylogenetic tree
/// independent of branch lengths or taxon labels:
/// - Sackin index: sum of depths of all leaves
/// - Colless index: sum of |L(left) - L(right)| at each internal node
/// - Tree balance ratio: proportion of balanced internal nodes

///|
/// Tree node interface for shape metrics.
/// Uses a simple recursive structure to avoid coupling to specific tree types.
pub struct SkbioTreeNode {
  name : String
  children : Array[SkbioTreeNode]
  is_leaf : Bool
}

///|
pub fn SkbioTreeNode::leaf(name : String) -> SkbioTreeNode {
  { name, children: [], is_leaf: true, }
}

///|
pub fn SkbioTreeNode::internal(
  name : String,
  children : Array[SkbioTreeNode],
) -> SkbioTreeNode {
  { name, children, is_leaf: false, }
}

///|
pub fn SkbioTreeNode::name(self : SkbioTreeNode) -> String {
  self.name
}

///|
pub fn SkbioTreeNode::is_leaf(self : SkbioTreeNode) -> Bool {
  self.is_leaf
}

///|
pub fn SkbioTreeNode::children(self : SkbioTreeNode) -> Array[SkbioTreeNode] {
  self.children
}

///|
/// Count the number of leaves in a subtree.
fn skbio_shape_count_leaves(node : SkbioTreeNode) -> Int32 {
  if node.is_leaf {
    return 1
  }
  let mut count = 0
  for c in node.children {
    count = count + skbio_shape_count_leaves(c)
  }
  count
}

///|
/// Compute the Sackin index: sum of depths of all leaves.
/// Depth of root = 0, depth of a child = depth of parent + 1.
/// A balanced tree has lower Sackin index; a caterpillar (maximally imbalanced)
/// tree has the maximum Sackin index.
pub fn skbio_sackin_index(root : SkbioTreeNode) -> Int32 {
  let mut total = 0
  fn compute_depth(node : SkbioTreeNode, depth : Int32) -> Int32 {
    if node.is_leaf {
      return depth
    }
    let mut sum = 0
    for c in node.children {
      sum = sum + compute_depth(c, depth + 1)
    }
    sum
  }
  total = compute_depth(root, 0)
  total
}

///|
/// Compute the Colless index: sum of |L(left) - L(right)| at each internal node.
/// For binary trees only. A perfectly balanced tree has Colless = 0.
pub fn skbio_colless_index(root : SkbioTreeNode) -> Int32 {
  fn compute(node : SkbioTreeNode) -> Int32 {
    if node.is_leaf {
      return 0
    }
    let mut sum = 0
    // For binary nodes: |L(left) - L(right)|
    if node.children.length() == 2 {
      let l = skbio_shape_count_leaves(node.children[0])
      let r = skbio_shape_count_leaves(node.children[1])
      sum = sum + (l - r).abs()
    }
    // For multifurcating nodes: max - min leaf count among children
    if node.children.length() > 2 {
      let mut max_l = 0
      let mut min_l = 999999
      for c in node.children {
        let l = skbio_shape_count_leaves(c)
        if l > max_l {
          max_l = l
        }
        if l < min_l {
          min_l = l
        }
      }
      sum = sum + (max_l - min_l)
    }
    for c in node.children {
      sum = sum + compute(c)
    }
    sum
  }
  compute(root)
}

///|
/// Normalize the Colless index to range [0, 1].
/// C_norm = C / C_max, where C_max = (n-1)(n-2)/2 for n leaves.
pub fn skbio_colless_normalized(root : SkbioTreeNode) -> Double {
  let n = skbio_shape_count_leaves(root).to_double()
  let c = skbio_colless_index(root).to_double()
  if n <= 2.0 {
    return 0.0
  }
  let c_max = (n - 1.0) * (n - 2.0) / 2.0
  if c_max > 0.0 {
    c / c_max
  } else {
    0.0
  }
}

///|
/// Compute the number of cherries in a tree.
/// A cherry is a pair of leaves that share the same parent.
pub fn skbio_count_cherries(root : SkbioTreeNode) -> Int32 {
  fn compute(node : SkbioTreeNode) -> Int32 {
    if node.is_leaf {
      return 0
    }
    let mut sum = 0
    // Check if this node is a cherry (all children are leaves)
    if node.children.length() == 2 &&
      node.children[0].is_leaf &&
      node.children[1].is_leaf {
      sum = sum + 1
    }
    for c in node.children {
      sum = sum + compute(c)
    }
    sum
  }
  compute(root)
}

///|
/// Compute tree balance ratio: proportion of balanced internal nodes.
/// A balanced binary node has equal leaf counts in both subtrees.
pub fn skbio_balance_ratio(root : SkbioTreeNode) -> Double {
  fn compute(
    node : SkbioTreeNode,
    balanced : Int32,
    total : Int32,
  ) -> (Int32, Int32) {
    if node.is_leaf {
      return (balanced, total)
    }
    let mut b = balanced
    let mut t = total
    if node.children.length() == 2 {
      t = t + 1
      let l = skbio_shape_count_leaves(node.children[0])
      let r = skbio_shape_count_leaves(node.children[1])
      if l == r {
        b = b + 1
      }
    }
    for c in node.children {
      let (cb, ct) = compute(c, b, t)
      b = cb
      t = ct
    }
    (b, t)
  }
  let (balanced, total) = compute(root, 0, 0)
  if total == 0 {
    return 1.0
  }
  balanced.to_double() / total.to_double()
}

///|
/// Compute the total number of internal nodes.
pub fn skbio_count_internal(root : SkbioTreeNode) -> Int32 {
  fn compute(node : SkbioTreeNode) -> Int32 {
    if node.is_leaf {
      return 0
    }
    let mut sum = 1
    for c in node.children {
      sum = sum + compute(c)
    }
    sum
  }
  compute(root)
}

///|
/// Compute the maximum depth of the tree.
pub fn skbio_tree_depth(root : SkbioTreeNode) -> Int32 {
  fn compute(node : SkbioTreeNode, depth : Int32) -> Int32 {
    if node.is_leaf {
      return depth
    }
    let mut max_d = depth
    for c in node.children {
      let d = compute(c, depth + 1)
      if d > max_d {
        max_d = d
      }
    }
    max_d
  }
  compute(root, 0)
}

///|
/// Check if a tree is balanced (all leaves at the same depth).
pub fn skbio_is_balanced(root : SkbioTreeNode) -> Bool {
  fn compute(node : SkbioTreeNode, depth : Int32) -> (Int32, Bool) {
    if node.is_leaf {
      return (depth, true)
    }
    let mut leaf_depth = -1
    for c in node.children {
      let (d, ok) = compute(c, depth + 1)
      if !ok {
        return (0, false)
      }
      if leaf_depth == -1 {
        leaf_depth = d
      } else if d != leaf_depth {
        return (0, false)
      }
    }
    (leaf_depth, true)
  }
  let (_, ok) = compute(root, 0)
  ok
}
