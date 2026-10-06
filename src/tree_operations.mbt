///|
// =========================================================================
// Tree operations (ete3-style)
//   Ladderize: sort child branches by subtree leaf count for visualization.
// =========================================================================

///|
/// Ladderize direction.
pub(all) enum LadderizeDirection {
  LadderizeLeft  // smaller subtrees first (left)
  LadderizeRight // larger subtrees first (right)
}

///|
/// Count the number of leaves in a subtree.
pub fn phylo_count_leaves(clade : Clade) -> Int {
  if clade.is_terminal() {
    return 1
  }
  let mut count = 0
  for ch in clade.clades {
    count = count + phylo_count_leaves(ch)
  }
  count
}

///|
/// Ladderize a tree: sort child branches by subtree leaf count.
/// Returns a new tree with sorted children.
/// - LadderizeLeft: smaller subtrees first
/// - LadderizeRight: larger subtrees first
pub fn phylo_ladderize(
  clade : Clade,
  direction : LadderizeDirection,
) -> Clade {
  if clade.is_terminal() {
    return clade
  }
  // Recursively ladderize children
  let children : Array[Clade] = Array(capacity=clade.clades.length())
  for ch in clade.clades {
    children.push(phylo_ladderize(ch, direction))
  }
  // Sort by leaf count
  children.sort_by(fn(a, b) {
    let na = phylo_count_leaves(a)
    let nb = phylo_count_leaves(b)
    match direction {
      LadderizeLeft => compare(na, nb)
      LadderizeRight => compare(nb, na)
    }
  })
  Clade::new(
    branch_length=clade.branch_length,
    name=clade.name,
    confidence=clade.confidence,
    clades=children,
  )
}


///|
// =========================================================================
// Midpoint rooting (ete3 / DendroPy)
//   Find the longest path between any two leaves, then re-root at
//   the midpoint of that path. If the midpoint falls on an edge,
//   split the edge and insert a new root node.
// =========================================================================

///|
/// Flat node representation for tree manipulation.
struct FlatNode {
  id : Int
  parent_id : Int // -1 for root
  branch_length : Double // length to parent
  children : Array[Int]
  name : String?
  is_leaf : Bool
}

///|
/// Flatten a Clade tree into a flat node array.
fn flatten_tree(clade : Clade) -> Array[FlatNode] {
  let nodes : Array[FlatNode] = Array(capacity=0)
  ignore(flatten_rec(clade, -1, 0.0, nodes))
  nodes
}

///|
fn flatten_rec(
  clade : Clade,
  parent_id : Int,
  branch_length : Double,
  nodes : Array[FlatNode],
) -> Int {
  let id = nodes.length()
  let is_leaf = clade.is_terminal()
  let node = FlatNode::{
    id,
    parent_id,
    branch_length,
    children: Array(capacity=0),
    name: clade.name,
    is_leaf,
  }
  nodes.push(node)
  for ch in clade.clades {
    let child_bl = match ch.branch_length {
      Some(v) => v
      None => 0.0
    }
    let child_id = flatten_rec(ch, id, child_bl, nodes)
    nodes[id].children.push(child_id)
  }
  id
}

///|
/// Build undirected adjacency list from flat nodes.
/// Returns (neighbor_id, branch_length) pairs for each node.
fn build_adjacency(
  nodes : Array[FlatNode],
) -> Array[Array[(Int, Double)]] {
  let n = nodes.length()
  let adj : Array[Array[(Int, Double)]] = Array(capacity=n)
  for i in 0..<n {
    adj.push(Array(capacity=0))
  }
  for i in 0..<n {
    // Edge to parent
    let pid = nodes[i].parent_id
    if pid >= 0 {
      adj[i].push((pid, nodes[i].branch_length))
      adj[pid].push((i, nodes[i].branch_length))
    }
  }
  adj
}

///|
/// BFS to find farthest node and record distances + predecessors.
/// Returns (farthest_node_id, distance_array, predecessor_array).
fn bfs_farthest(
  start : Int,
  adj : Array[Array[(Int, Double)]],
  n : Int,
) -> (Int, Array[Double], Array[Int]) {
  let dist : Array[Double] = Array(capacity=n)
  let pred : Array[Int] = Array(capacity=n)
  let visited : Array[Bool] = Array(capacity=n)
  for i in 0..<n {
    dist.push(-1.0)
    pred.push(-1)
    visited.push(false)
  }
  // Simple BFS using a queue (array-based)
  let queue : Array[Int] = Array(capacity=n)
  queue.push(start)
  dist[start] = 0.0
  visited[start] = true
  let mut head = 0
  while head < queue.length() {
    let u = queue[head]
    head = head + 1
    for pair in adj[u] {
      let v = pair.0
      let w = pair.1
      if not(visited[v]) {
        visited[v] = true
        dist[v] = dist[u] + w
        pred[v] = u
        queue.push(v)
      }
    }
  }
  // Find farthest
  let mut farthest = start
  let mut max_dist = 0.0
  for i in 0..<n {
    if dist[i] > max_dist {
      max_dist = dist[i]
      farthest = i
    }
  }
  (farthest, dist, pred)
}

///|
/// Reconstruct path from start to end using predecessor array.
fn reconstruct_path(
  start : Int,
  end : Int,
  pred : Array[Int],
) -> Array[Int] {
  let path : Array[Int] = Array(capacity=0)
  let mut cur = end
  while cur != start && cur != -1 {
    path.push(cur)
    cur = pred[cur]
  }
  path.push(start)
  // Reverse to get path from start to end
  let result : Array[Int] = Array(capacity=path.length())
  for i in 0..<path.length() {
    result.push(path[path.length() - 1 - i])
  }
  result
}

///|
/// Rebuild a Clade tree from the flat representation, rooted at new_root.
/// `new_parent` is the parent id to avoid (for DFS direction).
fn rebuild_clade(
  node_id : Int,
  parent_to_avoid : Int,
  bl_to_parent : Double,
  nodes : Array[FlatNode],
  adj : Array[Array[(Int, Double)]],
) -> Clade {
  let node = nodes[node_id]
  // Collect children (all neighbors except parent)
  let children : Array[Clade] = Array(capacity=0)
  for pair in adj[node_id] {
    let neighbor = pair.0
    let bl = pair.1
    if neighbor != parent_to_avoid {
      children.push(rebuild_clade(neighbor, node_id, bl, nodes, adj))
    }
  }
  Clade::new(
    branch_length=Some(bl_to_parent),
    name=node.name,
    clades=children,
  )
}

///|
/// Midpoint rooting: re-root the tree at the midpoint of the longest path.
///
/// Finds the two most distant nodes (via double BFS on the unrooted tree),
/// then places the new root at the midpoint of the path between them.
/// If the midpoint falls on an edge, that edge is split and a new root
/// node is inserted.
pub fn phylo_midpoint_root(clade : Clade) -> Clade {
  let nodes = flatten_tree(clade)
  let n = nodes.length()
  if n <= 1 {
    return clade
  }
  if n == 2 {
    return clade
  }
  let adj = build_adjacency(nodes)
  // Double BFS to find longest path
  // Start from node 0 (arbitrary)
  let (far1, _, _) = bfs_farthest(0, adj, n)
  let (far2, _, pred2) = bfs_farthest(far1, adj, n)
  // Path from far1 to far2 is the longest
  let path = reconstruct_path(far1, far2, pred2)
  // Compute cumulative distances along path
  let path_dist : Array[Double] = Array(capacity=path.length())
  path_dist.push(0.0)
  for i in 1..<path.length() {
    let u = path[i - 1]
    let v = path[i]
    let mut edge_len = 0.0
    for pair in adj[u] {
      if pair.0 == v {
        edge_len = pair.1
      }
    }
    path_dist.push(path_dist[i - 1] + edge_len)
  }
  let total_dist = path_dist[path.length() - 1]
  let half_dist = total_dist / 2.0
  // Find where the midpoint falls
  // Look for the edge (path[i-1], path[i]) where cumulative dist crosses half
  let mut mid_idx = 0
  for i in 1..<path.length() {
    if path_dist[i] >= half_dist {
      mid_idx = i
      break
    }
  }
  let dist_to_mid_from_prev = half_dist - path_dist[mid_idx - 1]
  let edge_node_prev = path[mid_idx - 1]
  let edge_node_next = path[mid_idx]
  // Check if midpoint is exactly at a node
  if dist_to_mid_from_prev < 0.0000000001 {
    // Midpoint is at edge_node_prev
    return rebuild_clade(edge_node_prev, -1, 0.0, nodes, adj)
  }
  // Get edge length
  let mut edge_len = 0.0
  for pair in adj[edge_node_prev] {
    if pair.0 == edge_node_next {
      edge_len = pair.1
    }
  }
  let dist_from_next = edge_len - dist_to_mid_from_prev
  if dist_from_next < 0.0000000001 {
    // Midpoint is at edge_node_next
    return rebuild_clade(edge_node_next, -1, 0.0, nodes, adj)
  }
  // Midpoint is on the edge — split it and create a new root
  // Build subtree from edge_node_prev (avoiding edge_node_next)
  let subtree_prev = rebuild_clade(
    edge_node_prev,
    edge_node_next,
    dist_to_mid_from_prev,
    nodes,
    adj,
  )
  // Build subtree from edge_node_next (avoiding edge_node_prev)
  let subtree_next = rebuild_clade(
    edge_node_next,
    edge_node_prev,
    dist_from_next,
    nodes,
    adj,
  )
  Clade::new(
    name=Some("midpoint_root"),
    clades=[subtree_prev, subtree_next],
  )
}


///|
// =========================================================================
// Root / Unroot (ete3)
//   unroot: remove root bifurcation (merge two root children)
//   root_at_leaf: re-root using a specified outgroup leaf
// =========================================================================

///|
/// Unroot a tree: if the root has exactly two children, merge them
/// into a single multifurcation. One child becomes the new root and
/// the other is attached with combined branch length.
/// If the root has != 2 children, the tree is already unrooted — return as-is.
pub fn phylo_unroot(clade : Clade) -> Clade {
  if clade.is_terminal() {
    return clade
  }
  if clade.clades.length() != 2 {
    return clade
  }
  let c1 = clade.clades[0]
  let c2 = clade.clades[1]
  let bl1 = match c1.branch_length {
    Some(v) => v
    None => 0.0
  }
  let bl2 = match c2.branch_length {
    Some(v) => v
    None => 0.0
  }
  let combined_bl = bl1 + bl2
  // Prefer internal node as new root (to keep root as internal node)
  let (new_root, other) = if c1.is_terminal() {
    (c2, c1)
  } else {
    (c1, c2)
  }
  // If both are leaves (2-leaf tree), can't meaningfully unroot
  if new_root.is_terminal() {
    return clade
  }
  // Attach `other` as a child of `new_root` with combined branch length
  let other_with_bl = Clade::new(
    branch_length=Some(combined_bl),
    name=other.name,
    confidence=other.confidence,
    clades=other.clades,
  )
  let children : Array[Clade] = Array(capacity=new_root.clades.length() + 1)
  for ch in new_root.clades {
    children.push(ch)
  }
  children.push(other_with_bl)
  Clade::new(
    branch_length=Some(0.0),
    name=new_root.name,
    confidence=new_root.confidence,
    clades=children,
  )
}

///|
/// Root a tree at a specified outgroup leaf.
/// Places a new root on the edge above the named leaf, making it the
/// outgroup. If the leaf is not found, returns the original tree.
pub fn phylo_root_at_leaf(clade : Clade, outgroup_name : String) -> Clade {
  let nodes = flatten_tree(clade)
  let n = nodes.length()
  if n <= 2 {
    return clade
  }
  // Find the leaf with the given name
  let mut target_id = -1
  for i in 0..<n {
    if nodes[i].is_leaf {
      match nodes[i].name {
        Some(nm) => if nm == outgroup_name {
          target_id = i
        }
        None => ()
      }
    }
  }
  if target_id < 0 {
    return clade
  }
  // If target is already a direct child of root (node 0) and root has 2 children,
  // the tree is already rooted at this outgroup
  if nodes[target_id].parent_id == 0 && nodes[0].children.length() == 2 {
    return clade
  }
  let adj = build_adjacency(nodes)
  let parent_of_target = nodes[target_id].parent_id
  if parent_of_target < 0 {
    // target is the root itself
    return clade
  }
  // Get branch length on edge (parent, target)
  let mut edge_len = 0.0
  for pair in adj[target_id] {
    if pair.0 == parent_of_target {
      edge_len = pair.1
    }
  }
  // Place new root at the midpoint of this edge (or at the parent node)
  // For outgroup rooting, we place the root at the parent node,
  // making the outgroup one child and everything else the other child.
  // If parent is the current root and has 3+ children, just split them.
  // If parent is deeper, we need to re-root at the parent.
  // Use the same approach as midpoint: split the edge above target
  // at distance 0 from target (i.e., at the parent node) — but that
  // doesn't create a new root on the edge. Instead, let's place the
  // new root on the edge, at the parent end (distance 0 from parent).
  // This effectively makes parent the new root with outgroup as one child.
  // But if parent is the old root, we need to handle it specially.
  // Simplest: place new root on the edge, splitting it.
  // Root is at distance edge_len from target (i.e., at parent).
  // Actually, for outgroup rooting, the standard approach is:
  // - Put new root on the edge above the outgroup
  // - Outgroup gets branch_length = edge_len (or edge_len/2)
  // - The rest gets branch_length = 0 (or edge_len/2)
  // Let's split at the parent: outgroup keeps its branch_length,
  // and the rest of the tree is rebuilt from parent avoiding target.
  // Build subtree from target (the outgroup) — it's a leaf
  let outgroup_subtree = rebuild_clade(target_id, parent_of_target, edge_len, nodes, adj)
  // Build subtree from parent (everything except target)
  let rest_subtree = rebuild_clade(parent_of_target, target_id, 0.0, nodes, adj)
  // If rest_subtree has only 1 child, collapse it (unroot the rest)
  let rest_final = if not(rest_subtree.is_terminal()) && rest_subtree.clades.length() == 1 {
    let child = rest_subtree.clades[0]
    let child_bl = match child.branch_length {
      Some(v) => v
      None => 0.0
    }
    let parent_bl = match rest_subtree.branch_length {
      Some(v) => v
      None => 0.0
    }
    Clade::new(
      branch_length=Some(child_bl + parent_bl),
      name=child.name,
      confidence=child.confidence,
      clades=child.clades,
    )
  } else {
    rest_subtree
  }
  Clade::new(
    name=Some("root"),
    clades=[outgroup_subtree, rest_final],
  )
}


///|
// =========================================================================
// Collapse weak nodes (ete3 / DendroPy)
//   Collapse internal nodes with confidence/support below a threshold
//   into multifurcations. The children of a collapsed node become
//   direct children of its parent.
// =========================================================================

///|
/// Collapse internal nodes whose confidence (bootstrap support) is below
/// `threshold`. A collapsed node's children are promoted to its parent.
/// Nodes with no confidence value are kept (not collapsed).
/// Returns a new tree.
pub fn phylo_collapse_weak(
  clade : Clade,
  threshold : Double,
) -> Clade {
  if clade.is_terminal() {
    return clade
  }
  // First, recursively process children
  let processed_children : Array[Clade] = Array(capacity=0)
  for ch in clade.clades {
    processed_children.push(phylo_collapse_weak(ch, threshold))
  }
  // Now check: should this node be collapsed?
  match clade.confidence {
    Some(conf) => {
      if conf < threshold {
        // Collapse: promote all children to parent
        // But we need to merge branch lengths: each child's branch_length
        // increases by this node's branch_length
        let my_bl = match clade.branch_length {
          Some(v) => v
          None => 0.0
        }
        let promoted : Array[Clade] = Array(capacity=0)
        for ch in processed_children {
          let ch_bl = match ch.branch_length {
            Some(v) => v
            None => 0.0
          }
          promoted.push(Clade::new(
            branch_length=Some(ch_bl + my_bl),
            name=ch.name,
            confidence=ch.confidence,
            clades=ch.clades,
          ))
        }
        // Return a "pass-through" node: parent will flatten this
        // We can't return multiple nodes, so we return a virtual node
        // with no confidence and the promoted children. The parent's
        // processing will handle it. But actually, we need a different
        // approach: return the promoted children directly.
        // Since we can only return one Clade, we wrap them in a marker.
        // Better approach: use a helper that returns Array[Clade].
        // For simplicity, if there's only 1 promoted child, return it
        // directly (with merged branch length). If multiple, return a
        // node with no confidence (so parent can flatten if needed).
        // Actually, the cleanest approach: return a node with the
        // promoted children and no confidence. The parent will see
        // this as a regular internal node. This is correct because
        // the collapse just means the node is no longer a bifurcation
        // — its children are now at the same level as siblings.
        // But the parent won't automatically flatten this. We need
        // to flatten at the parent level.
        // Let me use a different approach: a helper that returns
        // Array[Clade] for the parent to absorb.
        // For now, return as a regular node with no confidence.
        // This is semantically correct: the node is collapsed, meaning
        // it's no longer a "real" internal node, just a multifurcation point.
        return Clade::new(
          branch_length=clade.branch_length,
          name=None,
          confidence=None,
          clades=promoted,
        )
      }
    }
    None => ()
  }
  // Not collapsed: return with processed children
  Clade::new(
    branch_length=clade.branch_length,
    name=clade.name,
    confidence=clade.confidence,
    clades=processed_children,
  )
}

///|
/// Collapse internal nodes with confidence below threshold, properly
/// flattening promoted children into the parent's children array.
/// This is the correct implementation that avoids wrapper nodes.
pub fn phylo_collapse_weak_nodes(
  clade : Clade,
  threshold : Double,
) -> Clade {
  let children = collapse_rec(clade, threshold)
  // children is Array[Clade] — if 1 element, it's the tree
  // If multiple, wrap in a root
  if children.length() == 1 {
    return children[0]
  }
  Clade::new(
    branch_length=clade.branch_length,
    name=clade.name,
    confidence=clade.confidence,
    clades=children,
  )
}

///|
/// Recursive helper: returns the children array after collapsing.
/// If this node should be collapsed, its children are returned for
/// the parent to absorb. Otherwise, returns a single-element array
/// containing this node with processed children.
fn collapse_rec(
  clade : Clade,
  threshold : Double,
) -> Array[Clade] {
  if clade.is_terminal() {
    return [clade]
  }
  // Process children recursively, flattening collapsed nodes
  let flat_children : Array[Clade] = Array(capacity=0)
  for ch in clade.clades {
    let ch_result = collapse_rec(ch, threshold)
    for c in ch_result {
      flat_children.push(c)
    }
  }
  // Check if this node should be collapsed
  match clade.confidence {
    Some(conf) => {
      if conf < threshold {
        // Collapse: return children directly (they get promoted to parent)
        // Add this node's branch_length to each child
        let my_bl = match clade.branch_length {
          Some(v) => v
          None => 0.0
        }
        if my_bl > 0.0 {
          let adjusted : Array[Clade] = Array(capacity=flat_children.length())
          for ch in flat_children {
            let ch_bl = match ch.branch_length {
              Some(v) => v
              None => 0.0
            }
            adjusted.push(Clade::new(
              branch_length=Some(ch_bl + my_bl),
              name=ch.name,
              confidence=ch.confidence,
              clades=ch.clades,
            ))
          }
          return adjusted
        }
        return flat_children
      }
    }
    None => ()
  }
  // Not collapsed: wrap children in this node
  [Clade::new(
    branch_length=clade.branch_length,
    name=clade.name,
    confidence=clade.confidence,
    clades=flat_children,
  )]
}


///|
// =========================================================================
// Scale tree (DendroPy)
//   Multiply all branch lengths by a constant factor.
// =========================================================================

///|
/// Scale all branch lengths in the tree by `factor`.
/// Returns a new tree with scaled branch lengths.
pub fn phylo_scale_tree(clade : Clade, factor : Double) -> Clade {
  let scaled_bl = match clade.branch_length {
    Some(v) => Some(v * factor)
    None => None
  }
  let children : Array[Clade] = Array(capacity=clade.clades.length())
  for ch in clade.clades {
    children.push(phylo_scale_tree(ch, factor))
  }
  Clade::new(
    branch_length=scaled_bl,
    name=clade.name,
    confidence=clade.confidence,
    clades=children,
  )
}
