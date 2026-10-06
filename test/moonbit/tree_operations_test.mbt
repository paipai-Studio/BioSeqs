///|
/// Tests for tree operations (ladderize, count leaves).

///|
fn op_leaf(name : String) -> @src.Clade {
  @src.Clade::new(name=Some(name))
}

///|
fn op_internal(children : Array[@src.Clade]) -> @src.Clade {
  @src.Clade::new(name=Some("internal"), clades=children)
}

///|
/// Extract leaf names in order by in-order traversal.
fn op_leaf_names(clade : @src.Clade, names : Array[String]) -> Unit {
  if clade.is_terminal() {
    match clade.name {
      Some(n) => names.push(n)
      None => names.push("?")
    }
    return
  }
  for ch in clade.clades {
    op_leaf_names(ch, names)
  }
}

///|
test "count_leaves: single leaf" {
  assert_eq(@src.phylo_count_leaves(op_leaf("A")), 1)
}

///|
test "count_leaves: cherry" {
  let tree = op_internal([op_leaf("A"), op_leaf("B")])
  assert_eq(@src.phylo_count_leaves(tree), 2)
}

///|
test "count_leaves: balanced 8" {
  let tree = op_internal([
    op_internal([
      op_internal([op_leaf("A"), op_leaf("B")]),
      op_internal([op_leaf("C"), op_leaf("D")]),
    ]),
    op_internal([
      op_internal([op_leaf("E"), op_leaf("F")]),
      op_internal([op_leaf("G"), op_leaf("H")]),
    ]),
  ])
  assert_eq(@src.phylo_count_leaves(tree), 8)
}

///|
test "ladderize: leaf unchanged" {
  let leaf = op_leaf("A")
  let result = @src.phylo_ladderize(leaf, @src.LadderizeDirection::LadderizeLeft)
  assert_true(result.is_terminal())
}

///|
test "ladderize_left: small subtree first" {
  // (A, (B, (C, D))) — left child has 1 leaf, right has 3
  // Ladderize left: already small-first, should stay same order
  let tree = op_internal([
    op_leaf("A"),
    op_internal([op_leaf("B"), op_internal([op_leaf("C"), op_leaf("D")])]),
  ])
  let result = @src.phylo_ladderize(tree, @src.LadderizeDirection::LadderizeLeft)
  let names : Array[String] = Array(capacity=0)
  op_leaf_names(result, names)
  // A (1 leaf) should come before B,C,D (3 leaves)
  assert_eq(names[0], "A")
}

///|
test "ladderize_left: reorders large-first" {
  // ((B, (C, D)), A) — left child has 3 leaves, right has 1
  // Ladderize left: A (1 leaf) should come first
  let tree = op_internal([
    op_internal([op_leaf("B"), op_internal([op_leaf("C"), op_leaf("D")])]),
    op_leaf("A"),
  ])
  let result = @src.phylo_ladderize(tree, @src.LadderizeDirection::LadderizeLeft)
  let names : Array[String] = Array(capacity=0)
  op_leaf_names(result, names)
  assert_eq(names[0], "A")
}

///|
test "ladderize_right: large subtree first" {
  // (A, (B, (C, D))) — left has 1 leaf, right has 3
  // Ladderize right: (B,(C,D)) (3 leaves) should come first
  let tree = op_internal([
    op_leaf("A"),
    op_internal([op_leaf("B"), op_internal([op_leaf("C"), op_leaf("D")])]),
  ])
  let result = @src.phylo_ladderize(tree, @src.LadderizeDirection::LadderizeRight)
  let names : Array[String] = Array(capacity=0)
  op_leaf_names(result, names)
  // B,C,D should come before A
  assert_eq(names[3], "A")
}

///|
test "ladderize: preserves leaf count" {
  let tree = op_internal([
    op_internal([op_leaf("A"), op_leaf("B")]),
    op_internal([op_leaf("C"), op_internal([op_leaf("D"), op_leaf("E")])]),
  ])
  let result = @src.phylo_ladderize(tree, @src.LadderizeDirection::LadderizeLeft)
  assert_eq(@src.phylo_count_leaves(result), 5)
}

///|
test "ladderize: balanced tree unchanged" {
  // ((A,B),(C,D)) — both subtrees have 2 leaves, order preserved
  let tree = op_internal([
    op_internal([op_leaf("A"), op_leaf("B")]),
    op_internal([op_leaf("C"), op_leaf("D")]),
  ])
  let result = @src.phylo_ladderize(tree, @src.LadderizeDirection::LadderizeLeft)
  let names : Array[String] = Array(capacity=0)
  op_leaf_names(result, names)
  assert_eq(names[0], "A")
  assert_eq(names[1], "B")
  assert_eq(names[2], "C")
  assert_eq(names[3], "D")
}

///|
test "ladderize: recursive sorting" {
  // Root: ((C,D,E), (A,B)) — left has 3, right has 2
  // Ladderize left: (A,B) (2 leaves) first, then (C,D,E) (3 leaves)
  let tree = op_internal([
    op_internal([op_leaf("C"), op_leaf("D"), op_leaf("E")]),
    op_internal([op_leaf("A"), op_leaf("B")]),
  ])
  let result = @src.phylo_ladderize(tree, @src.LadderizeDirection::LadderizeLeft)
  let names : Array[String] = Array(capacity=0)
  op_leaf_names(result, names)
  // A,B should come first (2 leaves < 3 leaves)
  assert_eq(names[0], "A")
  assert_eq(names[1], "B")
  assert_eq(names[2], "C")
}

///|
test "ladderize: star tree stable" {
  // (A, B, C) — all leaves, all have 1 leaf, order preserved
  let tree = op_internal([op_leaf("A"), op_leaf("B"), op_leaf("C")])
  let result = @src.phylo_ladderize(tree, @src.LadderizeDirection::LadderizeLeft)
  let names : Array[String] = Array(capacity=0)
  op_leaf_names(result, names)
  assert_eq(names[0], "A")
  assert_eq(names[1], "B")
  assert_eq(names[2], "C")
}

///|
test "ladderize_left_right_inverse" {
  // Ladderize left then right should give reverse order at each node
  let tree = op_internal([
    op_leaf("A"),
    op_internal([op_leaf("B"), op_internal([op_leaf("C"), op_leaf("D")])]),
  ])
  let left = @src.phylo_ladderize(tree, @src.LadderizeDirection::LadderizeLeft)
  let right = @src.phylo_ladderize(tree, @src.LadderizeDirection::LadderizeRight)
  let names_left : Array[String] = Array(capacity=0)
  let names_right : Array[String] = Array(capacity=0)
  op_leaf_names(left, names_left)
  op_leaf_names(right, names_right)
  // A first in left, A last in right
  assert_eq(names_left[0], "A")
  assert_eq(names_right[3], "A")
}

///|
test "ladderize: preserves branch_length" {
  let tree = @src.Clade::new(
    name=Some("root"),
    clades=[
      @src.Clade::new(name=Some("A"), branch_length=Some(1.5)),
      @src.Clade::new(name=Some("B"), branch_length=Some(2.5)),
    ],
  )
  let result = @src.phylo_ladderize(tree, @src.LadderizeDirection::LadderizeLeft)
  // Both are leaves with 1 leaf each, order preserved
  assert_true(result.clades[0].branch_length.unwrap() - 1.5 < 0.0001)
  assert_true(result.clades[1].branch_length.unwrap() - 2.5 < 0.0001)
}
