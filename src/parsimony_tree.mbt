///|
/// Maximum parsimony tree construction ported from Biopython
/// ``Bio.Phylo.TreeConstruction``:
///
/// - [ParsimonyScorer] combines the Fitch algorithm (no matrix) and the
///   Sankoff algorithm (with a [ParsimonyMatrix]); unrooted input trees are
///   rooted at their midpoint before scoring, matching Biopython.
/// - [NNITreeSearcher] performs Nearest Neighbor Interchange hill climbing.
/// - [ParsimonyTreeConstructor] runs the searcher from a user supplied
///   starting tree (default: a UPGMA tree built from identity distances).
/// - [tree_root_at_midpoint] / [Tree::root_at_midpoint] port
///   ``Tree.root_at_midpoint``.
///
/// All operations are non-destructive: inputs are copied internally.

///|
/// Error raised when parsimony tree construction input is invalid.
pub suberror ParsimonyTreeError {
  ParsimonyTreeError(String)
}

// ---------------------------------------------------------------------------
// Mutable working node (Clade is immutable-with-copy in this package)
// ---------------------------------------------------------------------------

priv struct NniNode {
  mut nid : Int
  name : String?
  mut branch_length : Double?
  confidence : Double?
  children : Array[NniNode]
}

fn nni_new_node(
  name : String?,
  branch_length : Double?,
  confidence : Double?,
  children : Array[NniNode],
) -> NniNode {
  NniNode::{ nid: -1, name, branch_length, confidence, children }
}

fn nni_from_clade(clade : Clade) -> NniNode {
  let children : Array[NniNode] = Array(capacity=0)
  for child in clade.clades {
    children.push(nni_from_clade(child))
  }
  nni_new_node(clade.name, clade.branch_length, clade.confidence, children)
}

fn nni_to_clade(node : NniNode) -> Clade {
  let children : Array[Clade] = Array(capacity=0)
  for child in node.children {
    children.push(nni_to_clade(child))
  }
  Clade::new(
    branch_length=node.branch_length,
    name=node.name,
    confidence=node.confidence,
    clades=children,
  )
}

fn nni_clone(node : NniNode) -> NniNode {
  let children : Array[NniNode] = Array(capacity=0)
  for child in node.children {
    children.push(nni_clone(child))
  }
  nni_new_node(node.name, node.branch_length, node.confidence, children)
}

fn nni_assign_ids(root : NniNode) -> Unit {
  let counter = Ref(0)
  fn visit(node : NniNode, counter : Ref[Int]) -> Unit {
    node.nid = counter.val
    counter.val = counter.val + 1
    for child in node.children {
      visit(child, counter)
    }
  }
  visit(root, counter)
}

fn nni_find(root : NniNode, nid : Int) -> NniNode? {
  if root.nid == nid {
    return Some(root)
  }
  for child in root.children {
    match nni_find(child, nid) {
      Some(found) => return Some(found)
      None => ()
    }
  }
  None
}

fn nni_collect_preorder(root : NniNode, out : Array[NniNode]) -> Unit {
  out.push(root)
  for child in root.children {
    nni_collect_preorder(child, out)
  }
}

///|
/// All terminals in preorder.
fn nni_terminals(root : NniNode) -> Array[NniNode] {
  let all : Array[NniNode] = Array(capacity=0)
  nni_collect_preorder(root, all)
  let tips : Array[NniNode] = Array(capacity=0)
  for node in all {
    if node.children.length() == 0 {
      tips.push(node)
    }
  }
  tips
}

///|
/// Internal (non-terminal) nodes in breadth-first ("level") order.
fn nni_nonterminals_level(root : NniNode) -> Array[NniNode] {
  let result : Array[NniNode] = Array(capacity=0)
  let queue : Array[NniNode] = [root]
  let mut head = 0
  while head < queue.length() {
    let node = queue[head]
    head = head + 1
    if node.children.length() > 0 {
      result.push(node)
      for child in node.children {
        queue.push(child)
      }
    }
  }
  result
}

///|
/// Map each node id (except the root) to its parent node.
fn nni_parent_map(root : NniNode) -> Map[Int, NniNode] {
  let parents : Map[Int, NniNode] = Map([])
  fn visit(node : NniNode, parents : Map[Int, NniNode]) -> Unit {
    for child in node.children {
      parents[child.nid] = node
      visit(child, parents)
    }
  }
  visit(root, parents)
  parents
}

fn nni_child_index(parent : NniNode, nid : Int) -> Int {
  let mut i = 0
  while i < parent.children.length() {
    if parent.children[i].nid == nid {
      return i
    }
    i = i + 1
  }
  -1
}

///|
/// Remove the child at ``idx`` (Array.remove returns the removed value).
fn nni_remove_at(children : Array[NniNode], idx : Int) -> Unit {
  let _ = children.remove(idx)
}

fn nni_bl_or_zero(node : NniNode) -> Double {
  match node.branch_length {
    Some(v) => v
    None => 0.0
  }
}

///|
/// Bifurcation check mirroring Biopython: the root may have 3 descendants,
/// every other internal node must have exactly 2.
fn nni_tree_bifurcating(root : NniNode) -> Bool {
  let root_n = root.children.length()
  if root_n != 0 && root_n != 2 && root_n != 3 {
    return false
  }
  fn visit(node : NniNode) -> Bool {
    let n = node.children.length()
    if n != 0 && n != 2 {
      return false
    }
    for child in node.children {
      if !visit(child) {
        return false
      }
    }
    true
  }
  for child in root.children {
    if !visit(child) {
      return false
    }
  }
  true
}

// ---------------------------------------------------------------------------
// NNI neighbor enumeration (ported from NNITreeSearcher._get_neighbors)
// ---------------------------------------------------------------------------

///|
/// Generate all NNI neighbor trees of a binary rooted tree. Two neighbors are
/// generated for each internal edge, in the same order as Biopython
/// (internal nodes visited breadth-first; root edges skipped).
pub fn nni_get_neighbors(tree : Tree) -> Array[Tree] {
  let work = nni_from_clade(tree.root)
  nni_assign_ids(work)
  nni_neighbors_of(work)
}

fn nni_neighbors_of(root : NniNode) -> Array[Tree] {
  let parents = nni_parent_map(root)
  let neighbors : Array[Tree] = Array(capacity=0)
  let root_child_ids : Array[Int] = Array(capacity=0)
  let internals = nni_nonterminals_level(root)
  for clade in internals {
    if clade.nid == root.nid {
      let left = clade.children[0]
      let right = clade.children[1]
      root_child_ids.push(left.nid)
      root_child_ids.push(right.nid)
      if left.children.length() > 0 && right.children.length() > 0 {
        let left_right = left.children[1]
        let right_left = right.children[0]
        let right_right = right.children[1]
        // neighbor 1 (left_left + right_right)
        nni_remove_at(left.children, 1)
        nni_remove_at(right.children, 1)
        left.children.push(right_right)
        right.children.push(left_right)
        neighbors.push(Tree::new(nni_to_clade(nni_clone(root)), rooted=true))
        // neighbor 2 (left_left + right_left)
        nni_remove_at(left.children, 1)
        nni_remove_at(right.children, 0)
        left.children.push(right_left)
        right.children.push(right_right)
        neighbors.push(Tree::new(nni_to_clade(nni_clone(root)), rooted=true))
        // change back (left_left + left_right)
        nni_remove_at(left.children, 1)
        nni_remove_at(right.children, 0)
        left.children.push(left_right)
        right.children.insert(0, right_left)
      }
    } else if root_child_ids.contains(clade.nid) {
      // skip root child
      continue
    } else {
      let left = clade.children[0]
      let right = clade.children[1]
      let parent = parents[clade.nid]
      if nni_child_index(parent, clade.nid) == 0 {
        let sister = parent.children[1]
        // neighbor 1 (parent + right)
        nni_remove_at(parent.children, 1)
        nni_remove_at(clade.children, 1)
        parent.children.push(right)
        clade.children.push(sister)
        neighbors.push(Tree::new(nni_to_clade(nni_clone(root)), rooted=true))
        // neighbor 2 (parent + left)
        nni_remove_at(parent.children, 1)
        nni_remove_at(clade.children, 0)
        parent.children.push(left)
        clade.children.push(right)
        neighbors.push(Tree::new(nni_to_clade(nni_clone(root)), rooted=true))
        // change back (parent + sister)
        nni_remove_at(parent.children, 1)
        nni_remove_at(clade.children, 0)
        parent.children.push(sister)
        clade.children.insert(0, left)
      } else {
        let sister = parent.children[0]
        // neighbor 1 (parent + right)
        nni_remove_at(parent.children, 0)
        nni_remove_at(clade.children, 1)
        parent.children.insert(0, right)
        clade.children.push(sister)
        neighbors.push(Tree::new(nni_to_clade(nni_clone(root)), rooted=true))
        // neighbor 2 (parent + left)
        nni_remove_at(parent.children, 0)
        nni_remove_at(clade.children, 0)
        parent.children.insert(0, left)
        clade.children.push(right)
        neighbors.push(Tree::new(nni_to_clade(nni_clone(root)), rooted=true))
        // change back (parent + sister)
        nni_remove_at(parent.children, 0)
        nni_remove_at(clade.children, 0)
        parent.children.insert(0, sister)
        clade.children.insert(0, left)
      }
    }
  }
  neighbors
}

// ---------------------------------------------------------------------------
// Rerooting (ported from Tree.root_with_outgroup / Tree.root_at_midpoint)
// ---------------------------------------------------------------------------

///|
/// Full node path from the root (included) to the node with ``nid``.
fn nni_root_path(root : NniNode, nid : Int, parents : Map[Int, NniNode]) -> Array[NniNode] {
  let chain : Array[NniNode] = Array(capacity=0)
  let mut current = match nni_find(root, nid) {
    Some(node) => node
    None => abort("nni: internal path error")
  }
  chain.push(current)
  while current.nid != root.nid {
    current = parents[current.nid]
    chain.push(current)
  }
  // reverse into a new array (tip -> root becomes root -> tip)
  let path : Array[NniNode] = Array(capacity=0)
  let mut i = chain.length() - 1
  while i >= 0 {
    path.push(chain[i])
    i = i - 1
  }
  path
}

///|
/// Port of ``Tree.root_with_outgroup``. When ``outgroup_branch_length`` is
/// ``None`` and the outgroup clade is internal, that clade becomes the new
/// (possibly trifurcating) root; otherwise a new bifurcating root is created.
fn nni_root_with_outgroup(
  root0 : NniNode,
  outgroup_nid : Int,
  outgroup_branch_length : Double?,
) -> NniNode {
  let root = nni_clone(root0)
  nni_assign_ids(root)
  let parents = nni_parent_map(root)
  let outgroup = match nni_find(root, outgroup_nid) {
    Some(node) => node
    None => abort("nni: outgroup not found in tree")
  }
  let path = nni_root_path(root, outgroup_nid, parents)
  if path.length() == 1 {
    // outgroup is the current root -- no change
    return root
  }
  let old_root_bl = root.branch_length
  let mut prev_blen = nni_bl_or_zero(outgroup)

  let pair : (NniNode, NniNode) =
    if outgroup.children.length() == 0 || outgroup_branch_length is Some(_) {
      // create a new root with a (possibly zero-length) branch to the outgroup
      let forced_bl = match outgroup_branch_length {
        Some(v) => v
        None => 0.0
      }
      outgroup.branch_length = Some(forced_bl)
      let root_children : Array[NniNode] = [outgroup]
      let nr = nni_new_node(None, old_root_bl, None, root_children)
      if path.length() == 2 {
        // outgroup is a direct child of the root
        (nr, nr)
      } else {
        let parent = path[path.length() - 2]
        let idx = nni_child_index(parent, outgroup.nid)
        nni_remove_at(parent.children, idx)
        let parent_old_bl = nni_bl_or_zero(parent)
        parent.branch_length = Some(prev_blen - forced_bl)
        prev_blen = parent_old_bl
        nr.children.insert(0, parent)
        (nr, parent)
      }
    } else {
      // use the internal outgroup node as the new (trifurcating) root
      outgroup.branch_length = old_root_bl
      (outgroup, outgroup)
    }
  let new_root = pair.0
  let mut new_parent = pair.1

  // reverse the branches along the lineage above the outgroup, stopping
  // below the original root (it is handled separately below)
  let upper_start = path.length() - 3
  let mut k = upper_start
  while k >= 1 {
    let parent = path[k]
    let idx = nni_child_index(parent, new_parent.nid)
    nni_remove_at(parent.children, idx)
    let parent_old_bl = nni_bl_or_zero(parent)
    parent.branch_length = Some(prev_blen)
    prev_blen = parent_old_bl
    new_parent.children.insert(0, parent)
    new_parent = parent
    k = k - 1
  }

  // handle the original root
  let old_root = root
  let out_is_root_child = path.length() == 2
  if out_is_root_child {
    let idx = nni_child_index(old_root, outgroup.nid)
    nni_remove_at(old_root.children, idx)
  } else {
    let idx = nni_child_index(old_root, new_parent.nid)
    nni_remove_at(old_root.children, idx)
  }
  if old_root.children.length() == 1 {
    // drop the old bifurcating root, preserving total branch length
    let ingroup = old_root.children[0]
    match ingroup.branch_length {
      Some(v) => ingroup.branch_length = Some(v + prev_blen)
      None => ingroup.branch_length = Some(prev_blen)
    }
    new_parent.children.insert(0, ingroup)
  } else {
    // keep the old trifurcating/polytomous root as an internal node
    old_root.branch_length = Some(prev_blen)
    new_parent.children.insert(0, old_root)
  }
  new_root
}

///|
/// Path of node ids from the root (excluded) to each node, with cumulative
/// branch length from the root.
fn nni_depths(
  root : NniNode,
) -> (Map[Int, Double], Map[Int, NniNode]) {
  let depths : Map[Int, Double] = Map([])
  let parents = nni_parent_map(root)
  depths[root.nid] = nni_bl_or_zero(root)
  let order : Array[NniNode] = Array(capacity=0)
  nni_collect_preorder(root, order)
  for node in order {
    if node.nid == root.nid {
      continue
    }
    let parent = parents[node.nid]
    depths[node.nid] = depths[parent.nid] + nni_bl_or_zero(node)
  }
  (depths, parents)
}

///|
/// Distance between two nodes along the tree, using root depths and LCA.
fn nni_pair_distance(
  root : NniNode,
  id_a : Int,
  id_b : Int,
  depths : Map[Int, Double],
  parents : Map[Int, NniNode],
) -> Double {
  let ancestors : Map[Int, Unit] = Map([])
  let mut cur = id_a
  while true {
    ancestors[cur] = ()
    if cur == root.nid {
      break
    }
    cur = parents[cur].nid
  }
  cur = id_b
  let mut lca = root.nid
  while true {
    if ancestors.contains(cur) {
      lca = cur
      break
    }
    cur = parents[cur].nid
  }
  depths[id_a] + depths[id_b] - 2.0 * depths[lca]
}

fn nni_root_at_midpoint(root0 : NniNode) -> NniNode {
  let root = nni_clone(root0)
  nni_assign_ids(root)
  let (depths, parents) = nni_depths(root)
  let tips = nni_terminals(root)

  // identify the diameter pair: preorder scan, strict improvement
  let mut max_distance = 0.0
  let mut tip1_id = tips[0].nid
  let mut tip2_id = tips[0].nid
  for t1 in tips {
    for t2 in tips {
      let d = nni_pair_distance(root, t1.nid, t2.nid, depths, parents)
      if d > max_distance {
        max_distance = d
        tip1_id = t1.nid
        tip2_id = t2.nid
      }
    }
  }

  // root with the first diameter tip as a zero-branch outgroup
  let mut rerooted = nni_root_with_outgroup(root, tip1_id, Some(0.0))
  nni_assign_ids(rerooted)
  let (_, new_parents) = nni_depths(rerooted)

  // nids were reassigned; locate the second diameter tip by terminal name
  let tip2_name = match nni_find(root, tip2_id) {
    Some(n) => n.name
    None => None
  }
  let new_tips = nni_terminals(rerooted)
  let mut new_tip2 = new_tips[0]
  for tip in new_tips {
    if tip.name == tip2_name {
      new_tip2 = tip
    }
  }

  // depth to travel from the ingroup side toward the outgroup tip
  let root_bl = nni_bl_or_zero(rerooted)
  let mut remainder = 0.5 * (max_distance - root_bl)

  // walk the root-to-tip2 path (root excluded), find the crossed edge
  let path = nni_root_path(rerooted, new_tip2.nid, new_parents)
  let mut outgroup_node = rerooted
  let mut outgroup_bl = 0.0
  let mut found = false
  let mut i = 1
  while i < path.length() {
    let node = path[i]
    remainder = remainder - nni_bl_or_zero(node)
    if remainder < 0.0 {
      outgroup_node = node
      outgroup_bl = -remainder
      found = true
      break
    }
    i = i + 1
  }
  if !found {
    abort("root_at_midpoint: failed to find the midpoint")
  }
  rerooted = nni_root_with_outgroup(rerooted, outgroup_node.nid, Some(outgroup_bl))
  rerooted
}

///|
/// Port of ``Tree.root_at_midpoint``: reroot at the midpoint of the longest
/// pairwise tip-to-tip distance. Returns a new rooted [Tree].
pub fn tree_root_at_midpoint(tree : Tree) -> Tree {
  let rooted = nni_root_at_midpoint(nni_from_clade(tree.root))
  Tree::new(
    nni_to_clade(rooted),
    rooted=true,
    name=tree.name,
    weight=tree.weight,
  )
}

///|
/// Return a copy rooted at the midpoint of its most distant taxa.
pub fn Tree::root_at_midpoint(self : Tree) -> Tree {
  tree_root_at_midpoint(self)
}

// ---------------------------------------------------------------------------
// ParsimonyScorer
// ---------------------------------------------------------------------------

///|
/// Parsimony scorer combining the Fitch algorithm (no matrix) and the
/// Sankoff algorithm (with a penalty matrix). Port of Biopython's
/// ``Bio.Phylo.TreeConstruction.ParsimonyScorer``.
pub struct ParsimonyScorer {
  matrix : ParsimonyMatrix?
}

///|
/// Create a parsimony scorer; without a matrix Fitch scoring is used.
pub fn ParsimonyScorer::new(matrix? : ParsimonyMatrix? = None) -> ParsimonyScorer {
  ParsimonyScorer::{ matrix }
}

///|
/// Return the integer parsimony score of the tree given the alignment
/// (taxon name -> sequence). Unrooted trees are rooted at their midpoint
/// in a private copy, exactly as Biopython does.
pub fn ParsimonyScorer::get_score(
  self : ParsimonyScorer,
  tree : Tree,
  alignment : Map[String, String],
) -> Int raise ParsimonyTreeError {
  let work_root = nni_from_clade(tree.root)
  nni_assign_ids(work_root)
  if !nni_tree_bifurcating(work_root) {
    raise ParsimonyTreeError("The tree provided should be bifurcating.")
  }
  let scored_tree =
    if tree.rooted {
      tree
    } else {
      Tree::new(
        nni_to_clade(nni_root_at_midpoint(work_root)),
        rooted=true,
        name=tree.name,
        weight=tree.weight,
      )
    }
  // validate taxon names
  let term_names : Array[String] = Array(capacity=0)
  for tip in scored_tree.root.get_terminals() {
    match tip.name {
      Some(n) => term_names.push(n)
      None =>
        raise ParsimonyTreeError(
          "Taxon names of the input tree should be the same with the alignment."
        )
    }
  }
  term_names.sort()
  let align_names : Array[String] = alignment.keys().collect()
  align_names.sort()
  if term_names.length() != align_names.length() {
    raise ParsimonyTreeError(
      "Taxon names of the input tree should be the same with the alignment."
    )
  }
  let mut ni = 0
  while ni < term_names.length() {
    if term_names[ni] != align_names[ni] {
      raise ParsimonyTreeError(
        "Taxon names of the input tree should be the same with the alignment."
      )
    }
    ni = ni + 1
  }
  let value = match self.matrix {
    Some(m) => sankoff_parsimony_score(scored_tree, alignment, m)
    None => fitch_parsimony_score(scored_tree, alignment)
  }
  value.to_int()
}

// ---------------------------------------------------------------------------
// NNITreeSearcher
// ---------------------------------------------------------------------------

///|
/// Tree searcher using the Nearest Neighbor Interchange (NNI) hill-climbing
/// algorithm. Port of ``Bio.Phylo.TreeConstruction.NNITreeSearcher``.
pub struct NNITreeSearcher {
  scorer : ParsimonyScorer
}

///|
/// Create an NNI searcher with the given parsimony scorer.
pub fn NNITreeSearcher::new(scorer : ParsimonyScorer) -> NNITreeSearcher {
  NNITreeSearcher::{ scorer }
}

///|
/// Search for the most parsimonious tree reachable by NNI moves from the
/// starting tree. The search stops when no neighbor has a strictly smaller
/// score; among improving neighbors, the last one achieving the running
/// minimum is accepted, matching Biopython.
pub fn NNITreeSearcher::search(
  self : NNITreeSearcher,
  starting_tree : Tree,
  alignment : Map[String, String],
) -> Tree raise ParsimonyTreeError {
  let mut best_tree =
    if starting_tree.rooted {
      starting_tree
    } else {
      starting_tree.root_at_midpoint()
    }
  while true {
    let mut best_score = self.scorer.get_score(best_tree, alignment)
    let temp = best_score
    for neighbor in nni_get_neighbors(best_tree) {
      let score = self.scorer.get_score(neighbor, alignment)
      if score < best_score {
        best_score = score
        best_tree = neighbor
      }
    }
    if best_score >= temp {
      break
    }
  }
  best_tree
}

// ---------------------------------------------------------------------------
// ParsimonyTreeConstructor
// ---------------------------------------------------------------------------

///|
/// Construct a maximum parsimony tree by running an [NNITreeSearcher] from a
/// starting tree. Port of
/// ``Bio.Phylo.TreeConstruction.ParsimonyTreeConstructor``.
pub struct ParsimonyTreeConstructor {
  searcher : NNITreeSearcher
  starting_tree : Tree?
}

///|
/// Create the constructor; when ``starting_tree`` is ``None`` a UPGMA tree
/// from identity distances is used.
pub fn ParsimonyTreeConstructor::new(
  searcher : NNITreeSearcher,
  starting_tree? : Tree? = None,
) -> ParsimonyTreeConstructor {
  ParsimonyTreeConstructor::{ searcher, starting_tree }
}

// ---------------------------------------------------------------------------
// Default starting tree: Biopython-style UPGMA from identity distances
// ---------------------------------------------------------------------------

///|
/// Byte-wise lexicographic string comparison (Python/C-locale ``sorted``).
/// MoonBit's native ``String`` order is length-first (shortlex), so the
/// default ``sort`` would order taxon names differently from Biopython's
/// alignment order and change UPGMA tie-breaking.
fn nni_str_cmp(a : String, b : String) -> Int {
  let n = if a.length() < b.length() { a.length() } else { b.length() }
  let mut cmp = 0
  let mut i = 0
  while i < n && cmp == 0 {
    let ca = a.unsafe_get(i).to_int()
    let cb = b.unsafe_get(i).to_int()
    if ca < cb {
      cmp = -1
    } else if ca > cb {
      cmp = 1
    }
    i = i + 1
  }
  if cmp != 0 {
    cmp
  } else if a.length() < b.length() {
    -1
  } else if a.length() > b.length() {
    1
  } else {
    0
  }
}

///|
/// UPGMA matching Biopython's ``DistanceTreeConstructor.upgma`` exactly:
/// equal (unweighted) averaging of the two merged rows, ties between minimum
/// pairs resolved in favour of the last lower-triangle pair encountered, and
/// the final root carries a branch length of 0.0. This differs from the
/// package's leaf-count-weighted [TreeConstructor] UPGMA, so it is kept
/// private to the parsimony constructor for parity with Biopython.
fn nni_biopython_upgma(
  matrix : Array[Array[Double]],
  names : Array[String],
) -> Tree {
  let n = names.length()
  let dm : Array[Array[Double]] = Array(capacity=0)
  for row in matrix {
    let copy : Array[Double] = Array(capacity=0)
    for v in row {
      copy.push(v)
    }
    dm.push(copy)
  }
  let clades : Array[Clade] = Array(capacity=0)
  for name in names {
    clades.push(Clade::new(name=Some(name)))
  }
  let heights : Array[Double] = Array(capacity=0)
  let mut k = 0
  while k < n {
    heights.push(0.0)
    k = k + 1
  }

  while dm.length() > 1 {
    let len = dm.length()
    let mut min_dist = dm[1][0]
    let mut min_i = 1
    let mut min_j = 0
    let mut i = 1
    while i < len {
      let mut j = 0
      while j < i {
        // ``>=`` reproduces Biopython: the last equal pair wins.
        if min_dist >= dm[i][j] {
          min_dist = dm[i][j]
          min_i = i
          min_j = j
        }
        j = j + 1
      }
      i = i + 1
    }
    let height = min_dist / 2.0
    let clade1 = Clade::new(
      branch_length=Some(height - heights[min_i]),
      name=clades[min_i].name,
      confidence=clades[min_i].confidence,
      clades=clades[min_i].clades,
    )
    let clade2 = Clade::new(
      branch_length=Some(height - heights[min_j]),
      name=clades[min_j].name,
      confidence=clades[min_j].confidence,
      clades=clades[min_j].clades,
    )
    let inner = Clade::new(clades=[clade1, clade2])

    // average the two rows (unweighted, as Biopython does)
    let mut k2 = 0
    while k2 < len {
      if k2 != min_i && k2 != min_j {
        let avg = (dm[min_i][k2] + dm[min_j][k2]) / 2.0
        dm[min_j][k2] = avg
        dm[k2][min_j] = avg
      }
      k2 = k2 + 1
    }

    clades[min_j] = inner
    let _ = clades.remove(min_i)
    heights[min_j] = height
    let _ = heights.remove(min_i)
    let _ = dm.remove(min_i)
    for row in dm {
      let _ = row.remove(min_i)
    }
  }

  let root = Clade::new(
    branch_length=Some(0.0),
    name=clades[0].name,
    confidence=clades[0].confidence,
    clades=clades[0].clades,
  )
  Tree::new(root, rooted=true)
}

///|
/// Build the parsimony tree for the given alignment (taxon name -> sequence).
pub fn ParsimonyTreeConstructor::build_tree(
  self : ParsimonyTreeConstructor,
  alignment : Map[String, String],
) -> Tree raise ParsimonyTreeError {
  let starting_tree = match self.starting_tree {
    Some(t) => t
    None => {
      // default: UPGMA tree from identity distances (may raise SeqError)
      try {
        let pairs : Array[(String, String)] = alignment.to_array()
        pairs.sort_by(fn(a, b) { nni_str_cmp(a.0, b.0) })
        let calculator = DistanceCalculator::new(model="identity")
        let matrix = calculator.get_distance(pairs)
        let names : Array[String] = pairs.map(fn(kv) { kv.0 })
        nni_biopython_upgma(matrix, names)
      } catch {
        SeqError(message) => raise ParsimonyTreeError(message)
      }
    }
  }
  self.searcher.search(starting_tree, alignment)
}
