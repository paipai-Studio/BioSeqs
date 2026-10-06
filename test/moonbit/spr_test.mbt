///|
/// Tests for SPR (Subtree Prune and Regraft) moves.

///|
fn spr_leaf(name : String) -> @src.Clade {
  @src.Clade::new(name=Some(name))
}

///|
fn spr_internal(children : Array[@src.Clade]) -> @src.Clade {
  @src.Clade::new(name=Some("internal"), clades=children)
}

///|
fn spr_count_leaves(clade : @src.Clade) -> Int {
  if clade.is_terminal() {
    return 1
  }
  let mut c = 0
  for ch in clade.clades {
    c = c + spr_count_leaves(ch)
  }
  c
}

///|
fn spr_leaf_names(clade : @src.Clade, names : Array[String]) -> Unit {
  if clade.is_terminal() {
    match clade.name {
      Some(n) => names.push(n)
      None => names.push("?")
    }
    return
  }
  for ch in clade.clades {
    spr_leaf_names(ch, names)
  }
}

///|
fn spr_is_binary(clade : @src.Clade) -> Bool {
  if clade.is_terminal() {
    return true
  }
  if clade.clades.length() != 2 {
    return false
  }
  spr_is_binary(clade.clades[0]) && spr_is_binary(clade.clades[1])
}

///|
fn spr_same_leaves(a : @src.Clade, b : @src.Clade) -> Bool {
  let names_a : Array[String] = Array(capacity=0)
  let names_b : Array[String] = Array(capacity=0)
  spr_leaf_names(a, names_a)
  spr_leaf_names(b, names_b)
  if names_a.length() != names_b.length() {
    return false
  }
  for n in names_a {
    let mut found = false
    for m in names_b {
      if n == m {
        found = true
      }
    }
    if not(found) {
      return false
    }
  }
  true
}

///|
/// 4-leaf balanced tree: ((A,B),(C,D))
fn spr_balanced4() -> @src.Clade {
  spr_internal([
    spr_internal([spr_leaf("A"), spr_leaf("B")]),
    spr_internal([spr_leaf("C"), spr_leaf("D")]),
  ])
}

///|
/// 5-leaf tree: ((A,B),(C,(D,E)))
fn spr_tree5() -> @src.Clade {
  spr_internal([
    spr_internal([spr_leaf("A"), spr_leaf("B")]),
    spr_internal([spr_leaf("C"), spr_internal([spr_leaf("D"), spr_leaf("E")])]),
  ])
}

///|
test "spr_move: preserves leaf count" {
  let tree = spr_balanced4()
  // Prune node 1 (first child of root), regraft on node 3
  let result = @src.phylo_spr_move(tree, 1, 3)
  assert_eq(spr_count_leaves(result), 4)
}

///|
test "spr_move: preserves leaf set" {
  let tree = spr_balanced4()
  let result = @src.phylo_spr_move(tree, 1, 3)
  assert_true(spr_same_leaves(tree, result))
}

///|
test "spr_move: preserves binary structure" {
  let tree = spr_tree5()
  let result = @src.phylo_spr_move(tree, 1, 3)
  assert_true(spr_is_binary(result))
}

///|
test "spr_move: invalid prune (root) returns original" {
  let tree = spr_balanced4()
  let result = @src.phylo_spr_move(tree, 0, 1)
  let names_orig : Array[String] = Array(capacity=0)
  let names_result : Array[String] = Array(capacity=0)
  spr_leaf_names(tree, names_orig)
  spr_leaf_names(result, names_result)
  assert_eq(names_orig.length(), names_result.length())
}

///|
test "spr_move: 5-leaf tree valid move" {
  let tree = spr_tree5()
  let result = @src.phylo_spr_move(tree, 2, 4)
  assert_eq(spr_count_leaves(result), 5)
  assert_true(spr_same_leaves(tree, result))
}

///|
test "spr_move: 6-leaf tree" {
  let tree = spr_internal([
    spr_internal([spr_leaf("A"), spr_leaf("B")]),
    spr_internal([
      spr_internal([spr_leaf("C"), spr_leaf("D")]),
      spr_internal([spr_leaf("E"), spr_leaf("F")]),
    ]),
  ])
  let result = @src.phylo_spr_move(tree, 1, 4)
  assert_eq(spr_count_leaves(result), 6)
  assert_true(spr_same_leaves(tree, result))
}

///|
test "spr_neighbors: 4-leaf tree produces neighbors" {
  let tree = spr_balanced4()
  let neighbors = @src.phylo_spr_neighbors(tree)
  // Should produce some neighbors
  assert_true(neighbors.length() > 0)
  for n in neighbors {
    assert_eq(spr_count_leaves(n), 4)
    assert_true(spr_same_leaves(tree, n))
  }
}

///|
test "spr_neighbors: 5-leaf tree all valid" {
  let tree = spr_tree5()
  let neighbors = @src.phylo_spr_neighbors(tree)
  assert_true(neighbors.length() > 0)
  for n in neighbors {
    assert_eq(spr_count_leaves(n), 5)
    assert_true(spr_same_leaves(tree, n))
  }
}

///|
test "spr_neighbors: cherry has no neighbors" {
  let tree = spr_internal([spr_leaf("A"), spr_leaf("B")])
  let neighbors = @src.phylo_spr_neighbors(tree)
  // 2-leaf tree has no valid SPR moves
  assert_eq(neighbors.length(), 0)
}

///|
test "spr_neighbors: NNI is a special case of SPR" {
  // NNI moves should be a subset of SPR neighbors
  let tree = spr_balanced4()
  let nni_neighbors = @src.phylo_nni_neighbors(tree)
  let spr_neighbors = @src.phylo_spr_neighbors(tree)
  // SPR should produce at least as many neighbors as NNI
  assert_true(spr_neighbors.length() >= nni_neighbors.length())
}

///|
test "spr_move: regraft on different subtree" {
  // ((A,B),(C,(D,E))) — prune (A,B) subtree, regraft on (D,E) edge
  let tree = spr_tree5()
  // Node layout (preorder):
  // 0: root, 1: (A,B), 2: A, 3: B, 4: (C,(D,E)), 5: C, 6: (D,E), 7: D, 8: E
  let result = @src.phylo_spr_move(tree, 1, 6)
  assert_eq(spr_count_leaves(result), 5)
  assert_true(spr_same_leaves(tree, result))
  assert_true(spr_is_binary(result))
}

///|
test "spr_neighbors: all neighbors are binary" {
  let tree = spr_tree5()
  let neighbors = @src.phylo_spr_neighbors(tree)
  for n in neighbors {
    assert_true(spr_is_binary(n))
  }
}
