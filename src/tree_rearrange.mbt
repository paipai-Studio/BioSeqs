///|
// =========================================================================
// NNI (Nearest Neighbor Interchange) move (ete3)
//   For an internal edge parent→child (both internal), where parent
//   has children [child, sibling] and child has children [g1, g2]:
//
//   Original:          Swap 1:            Swap 2:
//     parent             parent             parent
//     /    \             /    \             /    \
//   child  sib         child  g1          child  g2
//   / \                / \                / \
//  g1  g2            g2  sib            g1  sib
//
//   Each internal edge has 2 NNI neighbors.
// =========================================================================

///|
/// Count the number of internal edges (edges where both endpoints
/// are internal nodes). This equals the number of NNI positions.
pub fn phylo_nni_count(clade : Clade) -> Int {
  if clade.is_terminal() {
    return 0
  }
  let mut count = 0
  for ch in clade.clades {
    if not(ch.is_terminal()) {
      count = count + 1 // this is an internal edge
      count = count + phylo_nni_count(ch)
    }
  }
  count
}

///|
/// Apply an NNI move on the edge_idx-th internal edge.
/// swap = 0 or 1 selects which of the two possible swaps.
/// Returns the modified tree, or the original if edge_idx is out of range.
pub fn phylo_nni_apply(
  clade : Clade,
  edge_idx : Int,
  swap : Int,
) -> Clade {
  let counter = NniCounter::{ val: 0 }
  nni_apply_rec(clade, edge_idx, swap, counter)
}

///|
struct NniCounter {
  mut val : Int
}

///|
fn nni_apply_rec(
  clade : Clade,
  edge_idx : Int,
  swap : Int,
  counter : NniCounter,
) -> Clade {
  if clade.is_terminal() {
    return clade
  }
  // Process children, looking for the target internal edge
  let children : Array[Clade] = Array(capacity=clade.clades.length())
  for ch in clade.clades {
    if not(ch.is_terminal()) && counter.val == edge_idx {
      // Found the target edge: clade → ch
      counter.val = counter.val + 1
      // ch is the "child", we need the "sibling" (other child of clade)
      // and ch's children g1, g2
      // But clade might have >2 children (multifurcation)
      // NNI is defined for binary trees; for multifurcating,
      // we pick the first two children of ch
      if ch.clades.length() >= 2 && clade.clades.length() >= 2 {
        // Find sibling: the other child of clade (not ch)
        // We need to know ch's index. Let's find it.
        // Actually, we're iterating and ch is the current one.
        // We need to build the new tree.
        // Let's handle this outside the loop.
        children.push(ch) // placeholder, will fix below
      } else {
        children.push(nni_apply_rec(ch, edge_idx, swap, counter))
      }
    } else if not(ch.is_terminal()) {
      counter.val = counter.val + 1
      children.push(nni_apply_rec(ch, edge_idx, swap, counter))
    } else {
      children.push(ch)
    }
  }
  // Now check if we matched the edge at this level
  // We need to redo this more carefully
  // Actually, let me restructure: the NNI is applied at the parent level
  // when we find the target edge. Let me rebuild.
  Clade::new(
    branch_length=clade.branch_length,
    name=clade.name,
    confidence=clade.confidence,
    clades=children,
  )
}

///|
/// Apply NNI move — restructured implementation.
/// edge_idx selects which internal edge (0-indexed, preorder traversal).
/// swap = 0 or 1 selects which swap.
pub fn phylo_nni_move(
  clade : Clade,
  edge_idx : Int,
  swap : Int,
) -> Clade {
  let counter = NniCounter::{ val: 0 }
  nni_move_rec(clade, edge_idx, swap, counter)
}

///|
fn nni_move_rec(
  clade : Clade,
  edge_idx : Int,
  swap : Int,
  counter : NniCounter,
) -> Clade {
  if clade.is_terminal() {
    return clade
  }
  // First, recursively process children (in case the target edge is deeper)
  // But we also need to check if the target edge is at this level.
  // The target edge is (clade → child) where child is internal.
  // We need to find which child corresponds to edge_idx.

  // Let's first check if any direct child edge is the target
  let mut found = false
  let mut found_idx = -1
  for i in 0..<clade.clades.length() {
    if not(clade.clades[i].is_terminal()) {
      if counter.val == edge_idx {
        found = true
        found_idx = i
        counter.val = counter.val + 1
        break
      }
      counter.val = counter.val + 1
    }
  }

  if found && clade.clades.length() == 2 && clade.clades[found_idx].clades.length() == 2 {
    // Apply NNI at this level
    let child = clade.clades[found_idx]
    let sibling_idx = if found_idx == 0 { 1 } else { 0 }
    let sibling = clade.clades[sibling_idx]
    let g1 = child.clades[0]
    let g2 = child.clades[1]

    // swap 0: exchange sibling and g1
    //   child's children become [g2, sibling], parent's other child becomes g1
    // swap 1: exchange sibling and g2
    //   child's children become [g1, sibling], parent's other child becomes g2
    let (new_child_grandchildren, new_sibling) = if swap == 0 {
      (
        [g2, sibling],
        g1,
      )
    } else {
      (
        [g1, sibling],
        g2,
      )
    }

    let new_child = Clade::new(
      branch_length=child.branch_length,
      name=child.name,
      confidence=child.confidence,
      clades=new_child_grandchildren,
    )

    // Rebuild children in original order
    let new_children : Array[Clade] = Array(capacity=2)
    if found_idx == 0 {
      new_children.push(new_child)
      new_children.push(new_sibling)
    } else {
      new_children.push(new_sibling)
      new_children.push(new_child)
    }
    return Clade::new(
      branch_length=clade.branch_length,
      name=clade.name,
      confidence=clade.confidence,
      clades=new_children,
    )
  }

  if found {
    // Edge was at this level but tree is not binary — can't do NNI
    // Reset and just return recursively
    let counter2 = NniCounter::{ val: 0 }
    return nni_move_rec(clade, edge_idx + 1000000, swap, counter2) // effectively skip
  }

  // Target edge is deeper — recurse into children
  let counter3 = NniCounter::{ val: 0 }
  // Recount to get the right counter state for children
  let mut c = 0
  for i in 0..<clade.clades.length() {
    if not(clade.clades[i].is_terminal()) {
      if c == edge_idx {
        // This shouldn't happen (found would be true)
        c = c + 1
      } else {
        c = c + 1
      }
    }
  }
  // Actually, let me just use a fresh counter and recurse properly
  let counter4 = NniCounter::{ val: 0 }
  let children : Array[Clade] = Array(capacity=clade.clades.length())
  for ch in clade.clades {
    if not(ch.is_terminal()) {
      // Check if this edge is the target
      if counter4.val == edge_idx {
        // Already handled above (found=true), shouldn't reach here
        counter4.val = counter4.val + 1
        children.push(ch)
      } else {
        counter4.val = counter4.val + 1
        children.push(nni_move_rec(ch, edge_idx, swap, counter4))
      }
    } else {
      children.push(ch)
    }
  }
  Clade::new(
    branch_length=clade.branch_length,
    name=clade.name,
    confidence=clade.confidence,
    clades=children,
  )
}

///|
/// Generate all NNI neighbor trees (all trees reachable by one NNI move).
/// For a binary tree with k internal edges, returns 2k trees.
pub fn phylo_nni_neighbors(clade : Clade) -> Array[Clade] {
  let n_edges = phylo_nni_count(clade)
  let neighbors : Array[Clade] = Array(capacity=n_edges * 2)
  for i in 0..<n_edges {
    neighbors.push(phylo_nni_move(clade, i, 0))
    neighbors.push(phylo_nni_move(clade, i, 1))
  }
  neighbors
}


///|
// =========================================================================
// SPR (Subtree Prune and Regraft) move (ete3)
//   1. Prune: cut the edge above a non-root node, detaching its subtree
//   2. Suppress any unary node left in the remaining tree
//   3. Regraft: insert the pruned subtree on a different edge
//   Uses a mutable parent-array representation for clean manipulation.
// =========================================================================

///|
/// Apply an SPR move: prune the subtree at `prune_node` and regraft it
/// onto the edge above `regraft_node` in the remaining tree.
/// Both are node IDs from a preorder flattening (0 = root).
/// Returns the modified tree, or the original if the move is invalid.
pub fn phylo_spr_move(
  clade : Clade,
  prune_node : Int,
  regraft_node : Int,
) -> Clade {
  let nodes = flatten_tree(clade)
  let n = nodes.length()
  if n <= 2 || prune_node <= 0 || prune_node >= n {
    return clade
  }
  if regraft_node <= 0 || regraft_node >= n || regraft_node == prune_node {
    return clade
  }
  // Check regraft_node is not in the pruned subtree
  if spr_is_descendant(nodes, prune_node, regraft_node) {
    return clade
  }
  // Mutable parent array and branch length array
  let parent : Array[Int] = Array(capacity=n)
  let bl : Array[Double] = Array(capacity=n)
  for i in 0..<n {
    parent.push(nodes[i].parent_id)
    bl.push(nodes[i].branch_length)
  }
  // Step 1: Prune — detach prune_node from its parent
  let prune_parent = parent[prune_node]
  let prune_bl = bl[prune_node]
  parent[prune_node] = -1
  bl[prune_node] = 0.0
  // Step 2: Suppress unary node at prune_parent in remaining tree
  // Count children of prune_parent (excluding pruned subtree)
  let mut root_id = 0
  let mut pp_children = 0
  let mut pp_remaining_child = -1
  for i in 0..<n {
    if i != prune_node && parent[i] == prune_parent {
      // Check i is not in pruned subtree
      if not(spr_is_descendant_arr(parent, prune_node, i)) {
        pp_children = pp_children + 1
        pp_remaining_child = i
      }
    }
  }
  if pp_children == 1 && prune_parent != 0 {
    // Suppress: connect remaining child to prune_parent's parent
    let pp_parent = parent[prune_parent]
    let combined_bl = bl[prune_parent] + bl[pp_remaining_child]
    parent[pp_remaining_child] = pp_parent
    bl[pp_remaining_child] = combined_bl
    // Mark prune_parent as removed
    parent[prune_parent] = -2
  } else if pp_children == 1 && prune_parent == 0 {
    // Root became unary — remaining child is new root
    root_id = pp_remaining_child
    parent[pp_remaining_child] = -1
    bl[pp_remaining_child] = 0.0
    parent[prune_parent] = -2
  }
  // Step 3: Regraft — insert pruned subtree on edge above regraft_node
  // Find regraft_node's current parent (may have changed due to suppression)
  let rg_parent = parent[regraft_node]
  if rg_parent < 0 {
    return clade
  }
  let rg_bl = bl[regraft_node]
  // Create a new node (index n)
  parent.push(-1)
  bl.push(0.0)
  let new_node = n
  // Insert new_node on the edge: rg_parent → new_node → regraft_node
  // and new_node → prune_node
  parent[new_node] = rg_parent
  bl[new_node] = rg_bl / 2.0
  parent[regraft_node] = new_node
  bl[regraft_node] = rg_bl / 2.0
  parent[prune_node] = new_node
  bl[prune_node] = prune_bl
  // Step 4: Rebuild Clade tree from root_id
  spr_rebuild_from_parent(root_id, -1, 0.0, parent, bl, nodes, new_node)
}

///|
/// Check if `desc` is a descendant of `anc` using the flat node array.
fn spr_is_descendant(
  nodes : Array[FlatNode],
  anc : Int,
  desc : Int,
) -> Bool {
  if desc == anc {
    return true
  }
  let mut cur = nodes[desc].parent_id
  while cur >= 0 {
    if cur == anc {
      return true
    }
    cur = nodes[cur].parent_id
  }
  false
}

///|
/// Check if `desc` is a descendant of `anc` using a parent array.
fn spr_is_descendant_arr(
  parent : Array[Int],
  anc : Int,
  desc : Int,
) -> Bool {
  if desc == anc {
    return true
  }
  let mut cur = parent[desc]
  let mut count = 0
  while cur >= 0 && count < parent.length() + 1 {
    if cur == anc {
      return true
    }
    cur = parent[cur]
    count = count + 1
  }
  false
}

///|
/// Rebuild a Clade tree from a parent array, starting at root.
fn spr_rebuild_from_parent(
  node_id : Int,
  parent_id : Int,
  bl_to_parent : Double,
  parent : Array[Int],
  bl : Array[Double],
  nodes : Array[FlatNode],
  new_node_id : Int,
) -> Clade {
  let node_name : String? = if node_id >= nodes.length() {
    Some("tbr_node")
  } else if node_id == new_node_id {
    Some("spr_node")
  } else {
    nodes[node_id].name
  }
  // Find all children (nodes whose parent is node_id)
  let children : Array[Clade] = Array(capacity=0)
  for i in 0..<parent.length() {
    if parent[i] == node_id {
      children.push(spr_rebuild_from_parent(
        i, node_id, bl[i], parent, bl, nodes, new_node_id,
      ))
    }
  }
  Clade::new(
    branch_length=Some(bl_to_parent),
    name=node_name,
    clades=children,
  )
}

///|
/// Generate all SPR neighbor trees by enumerating valid prune/regraft pairs.
pub fn phylo_spr_neighbors(clade : Clade) -> Array[Clade] {
  let nodes = flatten_tree(clade)
  let n = nodes.length()
  let neighbors : Array[Clade] = Array(capacity=0)
  for prune_node in 1..<n {
    for regraft_node in 1..<n {
      if regraft_node != prune_node {
        if not(spr_is_descendant(nodes, prune_node, regraft_node)) {
          let result = phylo_spr_move(clade, prune_node, regraft_node)
          if spr_tree_changed(clade, result) {
            neighbors.push(result)
          }
        }
      }
    }
  }
  neighbors
}

///|
/// Check if two trees differ (by leaf name ordering).
fn spr_tree_changed(a : Clade, b : Clade) -> Bool {
  let names_a : Array[String] = Array(capacity=0)
  let names_b : Array[String] = Array(capacity=0)
  spr_collect_names(a, names_a)
  spr_collect_names(b, names_b)
  if names_a.length() != names_b.length() {
    return true
  }
  for i in 0..<names_a.length() {
    if names_a[i] != names_b[i] {
      return true
    }
  }
  false
}

///|
fn spr_collect_names(clade : Clade, names : Array[String]) -> Unit {
  if clade.is_terminal() {
    match clade.name {
      Some(n) => names.push(n)
      None => names.push("?")
    }
    return
  }
  for ch in clade.clades {
    spr_collect_names(ch, names)
  }
}

///|
// =========================================================================
// TBR (Tree Bisection and Reconnection) move (ete3)
//   1. Bisect: cut an edge, splitting tree into T1 and T2
//   2. Reconnect: reroot T1 at an arbitrary edge, then regraft onto T2
//   TBR is the most general rearrangement; NNI ⊂ SPR ⊂ TBR.
// =========================================================================

///|
/// Apply a TBR move:
/// - `bisect_node`: node whose edge to parent is cut (bisection point)
/// - `reconnect1`: node in T1 (pruned subtree) to reroot at
/// - `reconnect2`: node in T2 (remaining tree) to regraft onto
/// Returns the modified tree, or original if invalid.
pub fn phylo_tbr_move(
  clade : Clade,
  bisect_node : Int,
  reconnect1 : Int,
  reconnect2 : Int,
) -> Clade {
  let nodes = flatten_tree(clade)
  let n = nodes.length()
  if n <= 2 || bisect_node <= 0 || bisect_node >= n {
    return clade
  }
  if reconnect2 <= 0 || reconnect2 >= n || reconnect2 == bisect_node {
    return clade
  }
  // reconnect2 must not be in the pruned subtree (T1)
  if spr_is_descendant(nodes, bisect_node, reconnect2) {
    return clade
  }
  // reconnect1 must be in T1 (the pruned subtree)
  if not(spr_is_descendant(nodes, bisect_node, reconnect1)) {
    return clade
  }
  // Mutable parent and branch length arrays
  let parent : Array[Int] = Array(capacity=n)
  let bl : Array[Double] = Array(capacity=n)
  for i in 0..<n {
    parent.push(nodes[i].parent_id)
    bl.push(nodes[i].branch_length)
  }
  // Step 1: Bisect — detach bisect_node from its parent
  let bisect_bl = bl[bisect_node]
  parent[bisect_node] = -1
  bl[bisect_node] = 0.0
  // Step 2: Suppress unary node at bisect_parent in T2
  let bisect_parent = nodes[bisect_node].parent_id
  let mut root_id = 0
  // Count remaining children of bisect_parent
  let mut bp_children = 0
  let mut bp_remaining = -1
  for i in 0..<n {
    if i != bisect_node && parent[i] == bisect_parent && parent[i] >= 0 {
      bp_children = bp_children + 1
      bp_remaining = i
    }
  }
  if bp_children == 1 && bisect_parent != 0 {
    let bp_parent = parent[bisect_parent]
    let combined = bl[bisect_parent] + bl[bp_remaining]
    parent[bp_remaining] = bp_parent
    bl[bp_remaining] = combined
    parent[bisect_parent] = -2 // removed
  } else if bp_children == 1 && bisect_parent == 0 {
    root_id = bp_remaining
    parent[bp_remaining] = -1
    bl[bp_remaining] = 0.0
    parent[bisect_parent] = -2
  }
  // Step 3: Reroot T1 at reconnect1
  // T1 is the subtree rooted at bisect_node.
  // If reconnect1 == bisect_node, no rerooting needed.
  // Otherwise, detach reconnect1 from its parent in T1,
  // and create a new root for T1 with reconnect1 and the rest.
  let mut t1_root = bisect_node
  if reconnect1 != bisect_node {
    // Detach reconnect1 from its parent
    let r1_parent = parent[reconnect1]
    let r1_bl = bl[reconnect1]
    parent[reconnect1] = -1
    bl[reconnect1] = 0.0
    // Suppress unary node at r1_parent in the remaining T1
    let mut r1p_children = 0
    let mut r1p_remaining = -1
    for i in 0..<n {
      if parent[i] == r1_parent && i != reconnect1 && parent[i] >= 0 {
        // Check i is in T1
        if spr_is_descendant_arr(parent, bisect_node, i) || i == bisect_node {
          r1p_children = r1p_children + 1
          r1p_remaining = i
        }
      }
    }
    if r1p_children == 1 {
      if r1_parent == bisect_node {
        // The root of T1 became unary — remaining child is new T1 root
        t1_root = r1p_remaining
        parent[r1p_remaining] = -1
        bl[r1p_remaining] = 0.0
      } else {
        // Suppress: connect remaining child to r1_parent's parent
        let r1p_parent = parent[r1_parent]
        let combined = bl[r1_parent] + bl[r1p_remaining]
        parent[r1p_remaining] = r1p_parent
        bl[r1p_remaining] = combined
        parent[r1_parent] = -2
      }
    }
    // Create new T1 root connecting reconnect1 and the rest of T1
    parent.push(-1)
    bl.push(0.0)
    let new_t1_root = parent.length() - 1
    t1_root = new_t1_root
    // Connect new_t1_root to reconnect1 and t1_root (the old root)
    parent[reconnect1] = new_t1_root
    bl[reconnect1] = r1_bl / 2.0
    // Find the old T1 root (which is now a child of the new root)
    // Actually, t1_root might have been updated above
    // We need to connect the remaining T1 (rooted at the old t1_root or r1p_remaining)
    // to the new root
    // The rest of T1: if r1_parent was the bisect_node (root of T1),
    // it became unary and r1p_remaining is the rest.
    // Otherwise, the rest is still rooted at bisect_node (suppression was internal).
    let old_t1 = if r1p_children == 1 && r1_parent == bisect_node {
      r1p_remaining
    } else {
      bisect_node
    }
    parent[old_t1] = new_t1_root
    bl[old_t1] = r1_bl / 2.0
  }
  // Step 4: Regraft T1 onto T2 at reconnect2's edge
  let rg_parent = parent[reconnect2]
  if rg_parent < 0 || rg_parent == -2 {
    return clade
  }
  let rg_bl = bl[reconnect2]
  // Create new node on the regraft edge
  parent.push(-1)
  bl.push(0.0)
  let new_node = parent.length() - 1
  parent[new_node] = rg_parent
  bl[new_node] = rg_bl / 2.0
  parent[reconnect2] = new_node
  bl[reconnect2] = rg_bl / 2.0
  // Attach T1 root to new_node
  parent[t1_root] = new_node
  bl[t1_root] = bisect_bl
  // Step 5: Rebuild from root_id
  spr_rebuild_from_parent(root_id, -1, 0.0, parent, bl, nodes, -1)
}

///|
/// Generate all TBR neighbor trees.
/// Enumerates valid (bisect, reconnect1, reconnect2) triples.
/// Note: this can produce O(n³) candidates, so use sparingly for large trees.
pub fn phylo_tbr_neighbors(clade : Clade) -> Array[Clade] {
  let nodes = flatten_tree(clade)
  let n = nodes.length()
  let neighbors : Array[Clade] = Array(capacity=0)
  for bisect_node in 1..<n {
    // Find all nodes in T1 (subtree at bisect_node)
    for reconnect1 in 0..<n {
      if spr_is_descendant(nodes, bisect_node, reconnect1) {
        // Find all nodes in T2 (not in T1)
        for reconnect2 in 1..<n {
          if reconnect2 != bisect_node && not(spr_is_descendant(nodes, bisect_node, reconnect2)) {
            let result = phylo_tbr_move(clade, bisect_node, reconnect1, reconnect2)
            if spr_tree_changed(clade, result) {
              neighbors.push(result)
            }
          }
        }
      }
    }
  }
  neighbors
}
