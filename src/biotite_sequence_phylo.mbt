/// Biotite `sequence.phylo` — phylogenetic trees, distance-based
/// construction (UPGMA, neighbor joining), Newick I/O and tree
/// distance queries.

///|
/// A node of a phylogenetic tree.
///
/// Leaves carry an `index` (matrix row/column); internal nodes carry an
/// optional `name` and a list of `children`, each child remembering its
/// branch `distance` to its parent.
pub struct BiotiteTreeNode {
  children : Array[BiotiteTreeNode]
  mut distance : Double
  name : String?
  index : Int32?
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
/// Wrap a root node into a tree.
pub fn biotite_tree_new(root : BiotiteTreeNode) -> BiotiteTree {
  BiotiteTree::{ root, }
}

///|
/// Whether the node is a leaf (no children).
pub fn BiotiteTreeNode::is_leaf(self : BiotiteTreeNode) -> Bool {
  self.children.length() == 0
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
/// Collect the leaves of the tree, sorted by leaf index.
///
/// A root without index and children (empty tree) contributes nothing.
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
/// Bottom-up search for the two leaves and their path distances.
///
/// Returns (contains_a, dist_a, contains_b, dist_b, lca_done).
/// Once `lca_done` is true, dist_a + dist_b is the final leaf-leaf
/// distance via their most recent common ancestor and must be passed
/// upward unchanged. Otherwise dist_a / dist_b accumulate from the
/// leaves toward the current node (the pending LCA).
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
  match node.index {
    Some(ix) => {
      if ix == ia {
        fa = true
      }
      if ix == ib {
        fb = true
      }
    }
    None => ()
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
/// Distance between two leaves via their most recent common ancestor.
///
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
/// Render the tree in Newick format.
pub fn BiotiteTree::as_newick(self : BiotiteTree) -> String {
  let sb = StringBuilder()
  biotite_newick_node(self.root, true, sb)
  let _ = sb.write(";")
  sb.to_string()
}

///|
fn biotite_newick_node(
  node : BiotiteTreeNode,
  is_root : Bool,
  sb : StringBuilder,
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
          None => ()
        }
    }
  } else {
    let _ = sb.write("(")
    for i = 0; i < node.children.length(); i = i + 1 {
      if i > 0 {
        let _ = sb.write(",")
      }
      biotite_newick_node(node.children[i], false, sb)
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
  if !is_root {
    let _ = sb.write(":")
    let _ = sb.write(biotite_format_dist(node.distance))
  }
}

///|
/// Compact float formatting for Newick branch lengths.
fn biotite_format_dist(d : Double) -> String {
  let s = d.to_string()
  if s.has_suffix(".0") {
    s[0:s.length() - 2].to_owned()
  } else {
    s
  }
}

///|
/// One-cell mutable cell for parser position bookkeeping.
priv struct BiotitePosRef {
  mut value : Int32
}

///|
/// Parse a Newick string into a tree.
///
/// Leaves receive indices in depth-first order (0-based), matching the
/// order they appear in the string.
pub fn biotite_tree_from_newick(newick : String) -> BiotiteTree {
  let cleaned = newick.trim().to_owned()
  let mut text = cleaned
  if text.has_suffix(";") {
    text = text[0:text.length() - 1].to_owned()
  }
  let pos = BiotitePosRef::{ value: 0, }
  let next_leaf = BiotitePosRef::{ value: 0, }
  let root = biotite_parse_newick_node(text, pos, next_leaf)
  BiotiteTree::{ root, }
}

///|
/// Recursive Newick node parser; `pos.value` points at the character
/// after the parsed node.
fn biotite_parse_newick_node(
  text : String,
  pos : BiotitePosRef,
  next_leaf : BiotitePosRef,
) -> BiotiteTreeNode {
  if pos.value < text.length() && text[pos.value] == '(' {
    pos.value = pos.value + 1
    let children : Array[BiotiteTreeNode] = []
    let distances : Array[Double] = []
    let mut parsing = true
    while parsing {
      let child = biotite_parse_newick_node(text, pos, next_leaf)
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
    // optional internal name
    let name = biotite_parse_newick_name(text, pos)
    // optional branch length of this node
    let dist = biotite_parse_newick_dist(text, pos)
    let node = biotite_tree_node_new(name~, children~, distances~)
    node.distance = dist
    node
  } else {
    // leaf
    let name = biotite_parse_newick_name(text, pos)
    let dist = biotite_parse_newick_dist(text, pos)
    let idx = next_leaf.value
    next_leaf.value = next_leaf.value + 1
    let node = biotite_tree_node_new(name~, index=Some(idx))
    node.distance = dist
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
}
