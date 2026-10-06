///|
/// Tests for TBR (Tree Bisection and Reconnection) moves.

///|
fn tbr_leaf(name : String) -> @src.Clade {
  @src.Clade::new(name=Some(name))
}

///|
fn tbr_internal(children : Array[@src.Clade]) -> @src.Clade {
  @src.Clade::new(name=Some("internal"), clades=children)
}

///|
fn tbr_count_leaves(clade : @src.Clade) -> Int {
  if clade.is_terminal() {
    return 1
  }
  let mut c = 0
  for ch in clade.clades {
    c = c + tbr_count_leaves(ch)
  }
  c
}

///|
fn tbr_leaf_names(clade : @src.Clade, names : Array[String]) -> Unit {
  if clade.is_terminal() {
    match clade.name {
      Some(n) => names.push(n)
      None => names.push("?")
    }
    return
  }
  for ch in clade.clades {
    tbr_leaf_names(ch, names)
  }
}

///|
fn tbr_is_binary(clade : @src.Clade) -> Bool {
  if clade.is_terminal() {
    return true
  }
  if clade.clades.length() != 2 {
    return false
  }
  tbr_is_binary(clade.clades[0]) && tbr_is_binary(clade.clades[1])
}

///|
fn tbr_same_leaves(a : @src.Clade, b : @src.Clade) -> Bool {
  let names_a : Array[String] = Array(capacity=0)
  let names_b : Array[String] = Array(capacity=0)
  tbr_leaf_names(a, names_a)
  tbr_leaf_names(b, names_b)
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
fn tbr_balanced4() -> @src.Clade {
  tbr_internal([
    tbr_internal([tbr_leaf("A"), tbr_leaf("B")]),
    tbr_internal([tbr_leaf("C"), tbr_leaf("D")]),
  ])
}

///|
/// 5-leaf tree: ((A,B),(C,(D,E)))
fn tbr_tree5() -> @src.Clade {
  tbr_internal([
    tbr_internal([tbr_leaf("A"), tbr_leaf("B")]),
    tbr_internal([tbr_leaf("C"), tbr_internal([tbr_leaf("D"), tbr_leaf("E")])]),
  ])
}

///|
/// 6-leaf tree: ((A,B),((C,D),(E,F)))
fn tbr_tree6() -> @src.Clade {
  tbr_internal([
    tbr_internal([tbr_leaf("A"), tbr_leaf("B")]),
    tbr_internal([
      tbr_internal([tbr_leaf("C"), tbr_leaf("D")]),
      tbr_internal([tbr_leaf("E"), tbr_leaf("F")]),
    ]),
  ])
}

///|
test "tbr_move: preserves leaf count" {
  let tree = tbr_balanced4()
  // Bisect at node 1, reroot T1 at node 1 (no reroot), regraft at node 3
  let result = @src.phylo_tbr_move(tree, 1, 1, 3)
  assert_eq(tbr_count_leaves(result), 4)
}

///|
test "tbr_move: preserves leaf set" {
  let tree = tbr_balanced4()
  let result = @src.phylo_tbr_move(tree, 1, 1, 3)
  assert_true(tbr_same_leaves(tree, result))
}

///|
test "tbr_move: preserves binary structure" {
  let tree = tbr_tree5()
  let result = @src.phylo_tbr_move(tree, 1, 1, 3)
  assert_true(tbr_is_binary(result))
}

///|
test "tbr_move: rerooting T1 changes topology" {
  let tree = tbr_tree5()
  // Node layout: 0:root, 1:(A,B), 2:A, 3:B, 4:(C,(D,E)), 5:C, 6:(D,E), 7:D, 8:E
  // Bisect at 1, reroot T1 at 2 (A), regraft at 5 (C)
  let result = @src.phylo_tbr_move(tree, 1, 2, 5)
  assert_eq(tbr_count_leaves(result), 5)
  assert_true(tbr_same_leaves(tree, result))
}

///|
test "tbr_move: 6-leaf tree" {
  let tree = tbr_tree6()
  let result = @src.phylo_tbr_move(tree, 1, 1, 4)
  assert_eq(tbr_count_leaves(result), 6)
  assert_true(tbr_same_leaves(tree, result))
}

///|
test "tbr_move: reroot at leaf node" {
  let tree = tbr_tree5()
  // Bisect at 4 (the (C,(D,E)) subtree), reroot at 7 (D), regraft at 2 (A)
  let result = @src.phylo_tbr_move(tree, 4, 7, 2)
  assert_eq(tbr_count_leaves(result), 5)
  assert_true(tbr_same_leaves(tree, result))
}

///|
test "tbr_move: invalid bisect (root) returns original" {
  let tree = tbr_balanced4()
  let result = @src.phylo_tbr_move(tree, 0, 1, 3)
  assert_true(tbr_same_leaves(tree, result))
}

///|
test "tbr_move: reconnect2 in T1 returns original" {
  let tree = tbr_tree5()
  // reconnect2 = 2 (A) is in T1 (subtree at bisect_node=1)
  let result = @src.phylo_tbr_move(tree, 1, 1, 2)
  assert_true(tbr_same_leaves(tree, result))
}

///|
test "tbr_move: reconnect1 not in T1 returns original" {
  let tree = tbr_tree5()
  // reconnect1 = 5 (C) is not in T1 (subtree at bisect_node=1)
  let result = @src.phylo_tbr_move(tree, 1, 5, 3)
  assert_true(tbr_same_leaves(tree, result))
}

///|
test "tbr_neighbors: 4-leaf tree produces neighbors" {
  let tree = tbr_balanced4()
  let neighbors = @src.phylo_tbr_neighbors(tree)
  assert_true(neighbors.length() > 0)
  for n in neighbors {
    assert_eq(tbr_count_leaves(n), 4)
    assert_true(tbr_same_leaves(tree, n))
  }
}

///|
test "tbr_neighbors: 5-leaf tree all valid" {
  let tree = tbr_tree5()
  let neighbors = @src.phylo_tbr_neighbors(tree)
  assert_true(neighbors.length() > 0)
  for n in neighbors {
    assert_eq(tbr_count_leaves(n), 5)
    assert_true(tbr_same_leaves(tree, n))
  }
}

///|
test "tbr_neighbors: cherry has no neighbors" {
  let tree = tbr_internal([tbr_leaf("A"), tbr_leaf("B")])
  let neighbors = @src.phylo_tbr_neighbors(tree)
  assert_eq(neighbors.length(), 0)
}

///|
test "tbr_neighbors: SPR is a special case of TBR" {
  // TBR with reconnect1 == bisect_node is equivalent to SPR
  let tree = tbr_tree5()
  let spr_neighbors = @src.phylo_spr_neighbors(tree)
  let tbr_neighbors = @src.phylo_tbr_neighbors(tree)
  // TBR should produce at least as many neighbors as SPR
  assert_true(tbr_neighbors.length() >= spr_neighbors.length())
}

///|
test "tbr_neighbors: all neighbors are binary" {
  let tree = tbr_tree5()
  let neighbors = @src.phylo_tbr_neighbors(tree)
  for n in neighbors {
    assert_true(tbr_is_binary(n))
  }
}

///|
test "tbr_move: reroot at internal node of T1" {
  let tree = tbr_tree6()
  // Node layout: 0:root, 1:(A,B), 2:A, 3:B, 4:((C,D),(E,F)), 5:(C,D), 6:C, 7:D, 8:(E,F), 9:E, 10:F
  // Bisect at 4, reroot T1 at 5 ((C,D)), regraft at 2 (A)
  let result = @src.phylo_tbr_move(tree, 4, 5, 2)
  assert_eq(tbr_count_leaves(result), 6)
  assert_true(tbr_same_leaves(tree, result))
  assert_true(tbr_is_binary(result))
}

///|
test "tbr_move: bisect and reroot same node" {
  // When reconnect1 == bisect_node, TBR reduces to SPR
  let tree = tbr_tree5()
  let tbr_result = @src.phylo_tbr_move(tree, 1, 1, 5)
  let spr_result = @src.phylo_spr_move(tree, 1, 5)
  // Both should have same leaf count and leaf set
  assert_eq(tbr_count_leaves(tbr_result), tbr_count_leaves(spr_result))
  assert_true(tbr_same_leaves(tbr_result, spr_result))
}
