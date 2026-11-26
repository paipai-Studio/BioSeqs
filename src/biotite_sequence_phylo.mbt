/// Biotite `sequence.phylo` — phylogenetic trees, distance-based
/// construction (UPGMA, neighbor joining), Newick I/O and tree
/// distance queries.
///
/// Node numbering follows the official Biotite semantics: leaves are
/// numbered `0..n-1` in the order they appear in the Newick notation,
/// intermediate nodes continue the numbering in the order they are
/// *closed* (post-order), and the root has the highest number.
/// `get_distance` accepts arbitrary node ids, not only leaves.

///|
/// A node of a phylogenetic tree.
///
/// Leaves carry an `index` (distance-matrix row/column, equal to the
/// node id for parsed trees); every node carries a `node_id` following
/// the Biotite numbering scheme.
pub struct BiotiteTreeNode {
  children : Array[BiotiteTreeNode]
  mut distance : Double
  name : String?
  index : Int32?
  mut node_id : Int32
  mut parent : BiotiteTreeNode?
}

///|
/// Create a tree node.
///
/// Either give `children` with `distances` (internal node) or `index`
/// (leaf). The parent pointers of all children are set to the new node.
pub fn biotite_tree_node_new(
  name~ : String? = None,
  index~ : Int32? = None,
  children~ : Array[BiotiteTreeNode] = [],
  distances~ : Array[Double] = [],
) -> BiotiteTreeNode {
  let node = BiotiteTreeNode::{
    children: [],
    distance: 0.0,
    name,
    index,
    node_id: -1,
    parent: None,
  }
  for i = 0; i < children.length(); i = i + 1 {
    let child = children[i]
    let d = if i < distances.length() { distances[i] } else { 0.0 }
    child.distance = d
    child.parent = Some(node)
    node.children.push(child)
  }
  node
}

///|
/// A phylogenetic tree with a root node.
pub struct BiotiteTree {
  root : BiotiteTreeNode
}

///|
/// Wrap a root node into a tree and assign node ids (Biotite scheme:
/// leaves keep their matrix index, intermediate nodes are numbered
/// from the leaf count upwards in post-order, so the root has the
/// highest id).
pub fn biotite_tree_new(root : BiotiteTreeNode) -> BiotiteTree {
  let tree = BiotiteTree::{ root, }
  biotite_tree_assign_ids(tree)
  tree
}

///|
/// Assign node ids in Biotite order (post-order, leaves first).
pub fn biotite_tree_assign_ids(tree : BiotiteTree) -> Unit {
  let n_leaves = tree.leaves().length()
  let counter : BiotitePosRef = BiotitePosRef::{ value: n_leaves, }
  biotite_assign_ids_rec(tree.root, counter)
}

///|
/// Post-order id assignment helper.
fn biotite_assign_ids_rec(
  node : BiotiteTreeNode,
  counter : BiotitePosRef,
) -> Unit {
  for child in node.children {
    biotite_assign_ids_rec(child, counter)
  }
  match node.index {
    Some(ix) =>
      if node.children.length() == 0 {
        node.node_id = ix
        return
      }
    None => ()
  }
  node.node_id = counter.value
  counter.value = counter.value + 1
}

///|
/// Whether the node is a leaf (no children).
pub fn BiotiteTreeNode::is_leaf(self : BiotiteTreeNode) -> Bool {
  self.children.length() == 0
}

///|
/// One-cell mutable cell for parser/count bookkeeping.
priv struct BiotitePosRef {
  mut value : Int32
}

///|
/// Ancestor chain from `node` up to the tree root (inclusive of both).
fn biotite_tree_node_ancestors(
  node : BiotiteTreeNode,
) -> Array[BiotiteTreeNode] {
  let out : Array[BiotiteTreeNode] = [node]
  let mut cur = node
  let mut done = false
  while !done {
    match cur.parent {
      Some(p) => {
        out.push(p)
        cur = p
      }
      None => done = true
    }
  }
  out
}

///|
/// Collect the leaves of the tree in depth-first (Newick appearance)
/// order, sorted by leaf index.
pub fn BiotiteTree::leaves(self : BiotiteTree) -> Array[BiotiteTreeNode] {
  let out : Array[BiotiteTreeNode] = []
  let stack : Array[BiotiteTreeNode] = [self.root]
  while stack.length() > 0 {
    let node = stack.pop().unwrap()
    if node.is_leaf() {
      if node.index is Some(_) {
        out.push(node)
      }
    } else {
      for child in node.children {
        stack.push(child)
      }
    }
  }
  out.sort_by(fn(a, b) {
    Int32::compare(a.index.unwrap_or(-1), b.index.unwrap_or(-1))
  })
  out
}

///|
/// Number of leaves.
pub fn BiotiteTree::n_leaves(self : BiotiteTree) -> Int32 {
  self.leaves().length()
}

///|
/// Find the node with the given node id (None when absent).
pub fn BiotiteTree::find_node(
  self : BiotiteTree,
  id : Int32,
) -> BiotiteTreeNode? {
  let stack : Array[BiotiteTreeNode] = [self.root]
  while stack.length() > 0 {
    let node = stack.pop().unwrap()
    if node.node_id == id {
      return Some(node)
    }
    for child in node.children {
      stack.push(child)
    }
  }
  None
}

///|
/// Bottom-up search for the two nodes and their path distances.
///
/// Returns (contains_a, dist_a, contains_b, dist_b, lca_done).
/// Once `lca_done` is true, dist_a + dist_b is the final node-node
/// distance via their most recent common ancestor and must be passed
/// upward unchanged. Otherwise dist_a / dist_b accumulate from the
/// nodes toward the current node (the pending LCA).
fn biotite_tree_walk(
  node : BiotiteTreeNode,
  ia : Int32,
  ib : Int32,
  topo : Bool,
) -> (Bool, Double, Bool, Double, Bool) {
  let mut fa = false
  let mut da = 0.0
  let mut fb = false
  let mut db = 0.0
  let mut lca_done = false
  if node.node_id == ia {
    fa = true
  }
  if node.node_id == ib {
    fb = true
  }
  for child in node.children {
    let (ca, dda, cb, ddb, done) = biotite_tree_walk(child, ia, ib, topo)
    if done {
      // LCA was below: final distance, pass through untouched.
      fa = true
      fb = true
      da = dda
      db = ddb
      lca_done = true
    } else {
      let step = if topo { 1.0 } else { child.distance }
      if ca {
        fa = true
        da = dda + step
      }
      if cb {
        fb = true
        db = ddb + step
      }
      // Both targets met within this child's subtree: the child is
      // the LCA; freeze the distance from here on.
      if ca && cb {
        lca_done = true
        da = dda
        db = ddb
      }
    }
  }
  (fa, da, fb, db, lca_done)
}

///|
/// Distance between two nodes via their most recent common ancestor.
///
/// `i` and `j` are node ids (leaves are `0..n-1`, intermediate nodes
/// continue the numbering, the root has the highest id).
/// `topological=true` counts each branch as length 1.
pub fn BiotiteTree::get_distance(
  self : BiotiteTree,
  i : Int32,
  j : Int32,
  topological~ : Bool = false,
) -> Double? {
  let (fa, da, fb, db, _done) = biotite_tree_walk(self.root, i, j, topological)
  if fa && fb {
    Some(da + db)
  } else {
    None
  }
}

///|
/// Pairwise distances between all leaves (n-by-n matrix, indexed by
/// leaf id).
pub fn BiotiteTree::get_leaf_distances(
  self : BiotiteTree,
  topological~ : Bool = false,
) -> Array[Array[Double]] {
  let leaves = self.leaves()
  let n = leaves.length()
  let out : Array[Array[Double]] = Array(capacity=n)
  for _i = 0; _i < n; _i = _i + 1 {
    out.push(Array::make(n, 0.0))
  }
  for i = 0; i < n; i = i + 1 {
    for j = i + 1; j < n; j = j + 1 {
      let id_i = leaves[i].node_id
      let id_j = leaves[j].node_id
      let d = match self.get_distance(id_i, id_j, topological~) {
        Some(v) => v
        None => 0.0
      }
      out[i][j] = d
      out[j][i] = d
    }
  }
  out
}

///|
/// Render the tree in Newick notation.
///
/// Mirrors the official `to_newick`: every node is suffixed with its
/// branch distance (the root gets `:0.0`); `include_distance=false`
/// omits the suffixes; `round_distance` fixes the number of decimals.
pub fn BiotiteTree::as_newick(
  self : BiotiteTree,
  include_distance~ : Bool = true,
  round_distance~ : Int32? = None,
) -> String {
  let sb = StringBuilder()
  biotite_newick_node(self.root, sb, include_distance, round_distance)
  let _ = sb.write(";")
  sb.to_string()
}

///|
fn biotite_newick_node(
  node : BiotiteTreeNode,
  sb : StringBuilder,
  include_distance : Bool,
  round_distance : Int32?,
) -> Unit {
  if node.is_leaf() {
    match node.name {
      Some(nm) => {
        let _ = sb.write(nm)
        ()
      }
      None =>
        match node.index {
          Some(ix) => {
            let _ = sb.write(ix.to_string())
            ()
          }
          None => {
            let _ = sb.write(node.node_id.to_string())
            ()
          }
        }
    }
  } else {
    let _ = sb.write("(")
    for i = 0; i < node.children.length(); i = i + 1 {
      if i > 0 {
        let _ = sb.write(",")
      }
      biotite_newick_node(
        node.children[i],
        sb,
        include_distance,
        round_distance,
      )
    }
    let _ = sb.write(")")
    match node.name {
      Some(nm) => {
        let _ = sb.write(nm)
        ()
      }
      None => ()
    }
  }
  if include_distance {
    let _ = sb.write(":")
    let _ = sb.write(biotite_format_dist(node.distance, round_distance))
  }
}

///|
/// Float formatting for Newick branch lengths: with a fixed number of
/// decimals when `round_distance` is given, otherwise the default
/// string form.
fn biotite_format_dist(d : Double, round_distance : Int32?) -> String {
  match round_distance {
    Some(digits) => {
      if digits <= 0 {
        return d.round().to_int().to_string()
      }
      let scale = @math.pow(10.0, digits.to_double())
      let scaled = (d * scale).round()
      let int_part = (scaled / scale).to_int()
      let frac = (scaled - int_part.to_double() * scale).to_int()
      let mut frac_str = frac.to_string()
      while frac_str.length() < digits {
        frac_str = "0" + frac_str
      }
      int_part.to_string() + "." + frac_str
    }
    None => d.to_string()
  }
}

///|
/// Parse a Newick string into a tree.
///
/// Leaves receive ids `0..n-1` in the order they appear; intermediate
/// nodes are numbered in the order they are closed (post-order), so
/// the root has the highest id — matching the official `from_newick`.
/// Leaf labels are kept as names; the leaf index equals the id.
pub fn biotite_tree_from_newick(newick : String) -> BiotiteTree {
  // Ignore any whitespace like the official parser.
  let joined = newick.trim().to_owned()
  let text = StringBuilder()
  for i = 0; i < joined.length(); i = i + 1 {
    let ch = joined[i]
    if ch != ' ' && ch != '\t' && ch != '\n' && ch != '\r' {
      let _ = text.write(joined[i:i + 1].to_owned())
    }
  }
  let mut clean = text.to_string()
  if clean.has_suffix(";") {
    clean = clean[0:clean.length() - 1].to_owned()
  }
  let pos = BiotitePosRef::{ value: 0, }
  let next_leaf = BiotitePosRef::{ value: 0, }
  let n_intermediate = BiotitePosRef::{ value: 0, }
  let root = biotite_parse_newick_node(clean, pos, next_leaf, n_intermediate)
  let n_leaves = next_leaf.value
  // Map the temporary negative intermediate ids to their final
  // numbers (first closed -> n_leaves, root -> highest).
  biotite_renumber_intermediate(root, n_leaves)
  BiotiteTree::{ root, }
}

///|
/// Convert temporary negative intermediate ids:
/// id = -1 - k  ->  n_leaves + k.
fn biotite_renumber_intermediate(
  node : BiotiteTreeNode,
  n_leaves : Int32,
) -> Unit {
  if node.node_id < 0 {
    node.node_id = n_leaves + (-node.node_id - 1)
  }
  for child in node.children {
    biotite_renumber_intermediate(child, n_leaves)
  }
}

///|
/// Recursive Newick node parser; `pos.value` points at the character
/// after the parsed node. Intermediate nodes get temporary negative
/// ids in closure order.
fn biotite_parse_newick_node(
  text : String,
  pos : BiotitePosRef,
  next_leaf : BiotitePosRef,
  n_intermediate : BiotitePosRef,
) -> BiotiteTreeNode {
  if pos.value < text.length() && text[pos.value] == '(' {
    pos.value = pos.value + 1
    let children : Array[BiotiteTreeNode] = []
    let distances : Array[Double] = []
    let mut parsing = true
    while parsing {
      let child = biotite_parse_newick_node(
        text, pos, next_leaf, n_intermediate,
      )
      children.push(child)
      distances.push(child.distance)
      if pos.value < text.length() && text[pos.value] == ',' {
        pos.value = pos.value + 1
      } else {
        parsing = false
      }
    }
    // consume ')'
    if pos.value < text.length() && text[pos.value] == ')' {
      pos.value = pos.value + 1
    }
    // optional internal name (discarded, like the official parser)
    let _name = biotite_parse_newick_name(text, pos)
    // optional branch length of this node
    let dist = biotite_parse_newick_dist(text, pos)
    let node = biotite_tree_node_new(children~, distances~)
    node.distance = dist
    node.node_id = -1 - n_intermediate.value
    n_intermediate.value = n_intermediate.value + 1
    node
  } else {
    // leaf
    let name = biotite_parse_newick_name(text, pos)
    let dist = biotite_parse_newick_dist(text, pos)
    let idx = next_leaf.value
    next_leaf.value = next_leaf.value + 1
    let node = biotite_tree_node_new(name~, index=Some(idx))
    node.distance = dist
    node.node_id = idx
    node
  }
}

///|
/// Read a label: all characters until ':' or ',' or ')' or end.
fn biotite_parse_newick_name(text : String, pos : BiotitePosRef) -> String? {
  let start = pos.value
  while pos.value < text.length() {
    let ch = text[pos.value]
    if ch == ':' || ch == ',' || ch == ')' {
      break
    }
    pos.value = pos.value + 1
  }
  let raw = text[start:pos.value].to_owned()
  let trimmed = raw.trim().to_owned()
  if trimmed.length() == 0 {
    None
  } else {
    Some(trimmed)
  }
}

///|
/// Read a branch length after ':' if present.
fn biotite_parse_newick_dist(text : String, pos : BiotitePosRef) -> Double {
  if pos.value < text.length() && text[pos.value] == ':' {
    pos.value = pos.value + 1
    let start = pos.value
    while pos.value < text.length() {
      let ch = text[pos.value]
      if ch == ',' || ch == ')' {
        break
      }
      pos.value = pos.value + 1
    }
    let raw = text[start:pos.value].to_owned().trim().to_owned()
    parse_double(raw).unwrap_or(0.0)
  } else {
    0.0
  }
}

///|
/// Convert a tree into a binary tree (official `as_binary`):
/// intermediate nodes with a single child are omitted (distances
/// summed), nodes with more than two children are split into a chain
/// of binary nodes. New intermediate nodes are numbered from the leaf
/// count upwards in post-order; the root has the highest id.
pub fn biotite_tree_as_binary(tree : BiotiteTree) -> BiotiteTree {
  let n_leaves = tree.leaves().length()
  let counter : BiotitePosRef = BiotitePosRef::{ value: n_leaves, }
  let (root, _extra) = biotite_as_binary_rec(tree.root, counter)
  // Rebuild parent pointers and node ids for the converted tree.
  let binary = BiotiteTree::{ root, }
  biotite_tree_assign_ids(binary)
  binary
}

///|
/// Recursive binary conversion; returns the converted subtree node and
/// the accumulated distance of omitted single-child nodes.
fn biotite_as_binary_rec(
  node : BiotiteTreeNode,
  counter : BiotitePosRef,
) -> (BiotiteTreeNode, Double) {
  if node.children.length() == 0 {
    return (node, 0.0)
  }
  if node.children.length() == 1 {
    let (child, extra) = biotite_as_binary_rec(node.children[0], counter)
    return (child, node.children[0].distance + extra)
  }
  // Convert every child first.
  let converted : Array[(BiotiteTreeNode, Double)] = []
  for child in node.children {
    converted.push(biotite_as_binary_rec(child, counter))
  }
  // Fold the children into a binary chain:
  // n1 = (c1, c2); n2 = (n1 with distance 0, c3); ...
  // Read the source distances BEFORE constructing the new node: the
  // converted children are the original nodes, and node_new would
  // overwrite their distances with the placeholders otherwise.
  let (first, extra0) = converted[0]
  let (second, extra1) = converted[1]
  let d0 = node.children[0].distance + extra0
  let d1 = node.children[1].distance + extra1
  let mut current = biotite_tree_node_new(children=[first, second], distances=[
    d0, d1,
  ])
  current.node_id = counter.value
  counter.value = counter.value + 1
  for k = 2; k < converted.length(); k = k + 1 {
    let (child, extra) = converted[k]
    let dk = node.children[k].distance + extra
    let parent = biotite_tree_node_new(children=[], distances=[])
    parent.node_id = counter.value
    counter.value = counter.value + 1
    // current joins with distance 0, the next child with its distance
    current.parent = Some(parent)
    parent.children.push(current)
    child.parent = Some(parent)
    child.distance = dk
    parent.children.push(child)
    current = parent
  }
  (current, 0.0)
}

///|
/// Build a tree with the UPGMA (average-linkage) method from an
/// n-by-n distance matrix.
pub fn biotite_upgma(distances : Array[Array[Double]]) -> BiotiteTree {
  let n = distances.length()
  if n == 0 {
    return biotite_tree_new(biotite_tree_node_new())
  }
  // Working matrix between clusters.
  let d : Array[Array[Double]] = Array(capacity=n)
  for i = 0; i < n; i = i + 1 {
    let row : Array[Double] = Array(capacity=n)
    for j = 0; j < n; j = j + 1 {
      row.push(distances[i][j])
    }
    d.push(row)
  }
  // Cluster records: member leaf indices, age, and the node built.
  let members : Array[Array[Int32]] = Array(capacity=n)
  let ages : Array[Double] = Array(capacity=n)
  let nodes : Array[BiotiteTreeNode?] = Array(capacity=n)
  let alive : Array[Bool] = Array(capacity=n)
  for i = 0; i < n; i = i + 1 {
    members.push([i])
    ages.push(0.0)
    nodes.push(Some(biotite_tree_node_new(index=Some(i))))
    alive.push(true)
  }
  let mut n_alive = n
  while n_alive > 1 {
    // Find the closest pair of live clusters.
    let mut best_i = -1
    let mut best_j = -1
    let mut best_d = 1.0e308
    for i = 0; i < n; i = i + 1 {
      if !alive[i] {
        continue
      }
      for j = i + 1; j < n; j = j + 1 {
        if !alive[j] {
          continue
        }
        if d[i][j] < best_d {
          best_d = d[i][j]
          best_i = i
          best_j = j
        }
      }
    }
    let i = best_i
    let j = best_j
    if i < 0 {
      break
    }
    let age = best_d / 2.0
    let child_a = nodes[i].unwrap()
    let child_b = nodes[j].unwrap()
    let new_node = biotite_tree_node_new(children=[child_a, child_b], distances=[
      age - ages[i],
      age - ages[j],
    ])
    // New cluster: average linkage over member pairs.
    for k = 0; k < n; k = k + 1 {
      if !alive[k] || k == i || k == j {
        continue
      }
      let mut sum = 0.0
      for a in members[i] {
        for b in members[k] {
          sum = sum + distances[a][b]
        }
      }
      for a in members[j] {
        for b in members[k] {
          sum = sum + distances[a][b]
        }
      }
      let cnt = (members[i].length() + members[j].length()).to_double() *
        members[k].length().to_double()
      let avg = if cnt > 0.0 { sum / cnt } else { 0.0 }
      d[i][k] = avg
      d[k][i] = avg
    }
    let merged : Array[Int32] = members[i].copy()
    for m in members[j] {
      merged.push(m)
    }
    members[i] = merged
    ages[i] = age
    nodes[i] = Some(new_node)
    alive[j] = false
    n_alive = n_alive - 1
  }
  let mut root_opt : BiotiteTreeNode? = None
  for i = 0; i < n; i = i + 1 {
    if alive[i] {
      root_opt = nodes[i]
    }
  }
  let root = match root_opt {
    Some(r) => r
    None => biotite_tree_node_new()
  }
  biotite_tree_new(root)
}

///|
/// Build a tree with the neighbor-joining method from an n-by-n
/// distance matrix.
pub fn biotite_neighbor_joining(
  distances : Array[Array[Double]],
) -> BiotiteTree {
  let n0 = distances.length()
  if n0 == 0 {
    return biotite_tree_new(biotite_tree_node_new())
  }
  if n0 == 1 {
    return biotite_tree_new(biotite_tree_node_new(index=Some(0)))
  }
  let mut d : Array[Array[Double]] = Array(capacity=n0)
  for i = 0; i < n0; i = i + 1 {
    let row : Array[Double] = Array(capacity=n0)
    for j = 0; j < n0; j = j + 1 {
      row.push(distances[i][j])
    }
    d.push(row)
  }
  let mut nodes : Array[BiotiteTreeNode?] = Array(capacity=n0)
  for i = 0; i < n0; i = i + 1 {
    nodes.push(Some(biotite_tree_node_new(index=Some(i))))
  }
  let mut m = n0
  while m > 2 {
    // Row sums r.
    let r : Array[Double] = Array(capacity=m)
    for i = 0; i < m; i = i + 1 {
      let mut s = 0.0
      for j = 0; j < m; j = j + 1 {
        if i != j {
          s = s + d[i][j]
        }
      }
      r.push(s)
    }
    // Minimal Q entry.
    let mut best_i = -1
    let mut best_j = -1
    let mut best_q = 1.0e308
    for i = 0; i < m; i = i + 1 {
      for j = i + 1; j < m; j = j + 1 {
        let q = (m - 2).to_double() * d[i][j] - r[i] - r[j]
        if q < best_q {
          best_q = q
          best_i = i
          best_j = j
        }
      }
    }
    let i = best_i
    let j = best_j
    if i < 0 {
      break
    }
    let dist_ij = d[i][j]
    let limb_i = 0.5 * dist_ij + (r[i] - r[j]) / (2.0 * (m - 2).to_double())
    let limb_j = dist_ij - limb_i
    let node_i = nodes[i].unwrap()
    let node_j = nodes[j].unwrap()
    let new_node = biotite_tree_node_new(children=[node_i, node_j], distances=[
      limb_i, limb_j,
    ])
    // Distances from the new cluster to all remaining ones.
    let new_row : Array[Double] = Array(capacity=m)
    for k = 0; k < m; k = k + 1 {
      if k == i || k == j {
        new_row.push(0.0)
      } else {
        new_row.push((d[i][k] + d[j][k] - dist_ij) / 2.0)
      }
    }
    // Compact the matrix: the new cluster becomes position 0, the
    // remaining original clusters follow in ascending order.
    let order : Array[Int32] = []
    for k = 0; k < m; k = k + 1 {
      if k != i && k != j {
        order.push(k)
      }
    }
    let next : Array[Array[Double]] = Array(capacity=m - 1)
    let row0 : Array[Double] = Array(capacity=m - 1)
    row0.push(0.0)
    for k in order {
      row0.push(new_row[k])
    }
    next.push(row0)
    for p = 0; p < order.length(); p = p + 1 {
      let a = order[p]
      let row : Array[Double] = Array(capacity=m - 1)
      row.push(new_row[a])
      for q = 0; q < order.length(); q = q + 1 {
        let b = order[q]
        row.push(if b == a { 0.0 } else { d[a][b] })
      }
      next.push(row)
    }
    d = next
    let new_nodes : Array[BiotiteTreeNode?] = Array(capacity=m - 1)
    new_nodes.push(Some(new_node))
    for k = 0; k < m; k = k + 1 {
      if k != i && k != j {
        new_nodes.push(nodes[k])
      }
    }
    nodes = new_nodes
    m = m - 1
  }
  // m == 2: connect the two remaining clusters via the root.
  let node_a = nodes[0].unwrap()
  let node_b = nodes[1].unwrap()
  let half = d[0][1] / 2.0
  let root = biotite_tree_node_new(children=[node_a, node_b], distances=[
    half, half,
  ])
  biotite_tree_new(root)
}

///|
pub fn _biotite_sequence_phylo_dummy() -> Unit {
  let t = biotite_tree_from_newick("((a:1,b:2):3,c:4);")
  let _ = t.as_newick()
  let _ = t.get_distance(0, 1)
  let _ = t.leaves()
  let _ = t.n_leaves()
  let _ = t.get_leaf_distances()
  let _ = biotite_tree_as_binary(t)
  let up = biotite_upgma([[0.0, 1.0, 2.0], [1.0, 0.0, 3.0], [2.0, 3.0, 0.0]])
  let _ = up.as_newick()
  let nj = biotite_neighbor_joining([
    [0.0, 1.0, 2.0],
    [1.0, 0.0, 3.0],
    [2.0, 3.0, 0.0],
  ])
  let _ = nj.as_newick()
  let node = biotite_tree_node_new(name=Some("x"))
  let _ = node.is_leaf()
  ignore(t.root)
  let _ = biotite_tree_node_ancestors(node)
  let _ = t.find_node(0)
}
