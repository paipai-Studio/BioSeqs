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
