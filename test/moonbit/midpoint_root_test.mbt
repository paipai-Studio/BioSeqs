///|
/// Tests for midpoint rooting.

///|
fn mp_leaf(name : String, bl : Double) -> @src.Clade {
  @src.Clade::new(name=Some(name), branch_length=Some(bl))
}

///|
fn mp_leaf0(name : String) -> @src.Clade {
  @src.Clade::new(name=Some(name))
}

///|
fn mp_internal(children : Array[@src.Clade], bl : Double) -> @src.Clade {
  @src.Clade::new(name=Some("internal"), branch_length=Some(bl), clades=children)
}

///|
fn mp_internal0(children : Array[@src.Clade]) -> @src.Clade {
  @src.Clade::new(name=Some("internal"), clades=children)
}

///|
/// Collect all leaf names from a tree.
fn mp_leaf_names(clade : @src.Clade, names : Array[String]) -> Unit {
  if clade.is_terminal() {
    match clade.name {
      Some(n) => names.push(n)
      None => names.push("?")
    }
    return
  }
  for ch in clade.clades {
    mp_leaf_names(ch, names)
  }
}

///|
/// Count total leaves in a tree.
fn mp_count_leaves(clade : @src.Clade) -> Int {
  if clade.is_terminal() {
    return 1
  }
  let mut c = 0
  for ch in clade.clades {
    c = c + mp_count_leaves(ch)
  }
  c
}

///|
test "midpoint_root: single leaf" {
  let leaf = mp_leaf0("A")
  let result = @src.phylo_midpoint_root(leaf)
  assert_true(result.is_terminal())
}

///|
test "midpoint_root: two leaves" {
  // (A:1, B:1) — midpoint is at root, should stay the same
  let tree = mp_internal0([mp_leaf("A", 1.0), mp_leaf("B", 1.0)])
  let result = @src.phylo_midpoint_root(tree)
  assert_eq(mp_count_leaves(result), 2)
}

///|
test "midpoint_root: preserves leaf count" {
  // ((A:1,B:1):1,C:3) — longest path A-C = 1+1+3=5, midpoint at 2.5
  let tree = mp_internal0([
    mp_internal0([mp_leaf("A", 1.0), mp_leaf("B", 1.0)]),
    mp_leaf("C", 3.0),
  ])
  let result = @src.phylo_midpoint_root(tree)
  assert_eq(mp_count_leaves(result), 3)
}

///|
test "midpoint_root: balanced tree" {
  // ((A:1,B:1):1,(C:1,D:1):1) — all paths equal, midpoint at root
  let tree = mp_internal0([
    mp_internal0([mp_leaf("A", 1.0), mp_leaf("B", 1.0)]),
    mp_internal0([mp_leaf("C", 1.0), mp_leaf("D", 1.0)]),
  ])
  let result = @src.phylo_midpoint_root(tree)
  assert_eq(mp_count_leaves(result), 4)
  // All leaf names preserved
  let names : Array[String] = Array(capacity=0)
  mp_leaf_names(result, names)
  assert_eq(names.length(), 4)
}

///|
test "midpoint_root: caterpillar shifts root" {
  // (A:1,(B:1,(C:1,D:1):1):1)
  // Longest path: A to C (or D) = 1+1+1+1 = 4, midpoint at 2.0
  // Path A -> root -> (B,(C,D)) -> (C,D) -> C
  // At 2.0 from A: that's at the (B,(C,D)) node
  let tree = mp_internal0([
    mp_leaf("A", 1.0),
    mp_internal0([
      mp_leaf("B", 1.0),
      mp_internal0([mp_leaf("C", 1.0), mp_leaf("D", 1.0)]),
    ]),
  ])
  let result = @src.phylo_midpoint_root(tree)
  assert_eq(mp_count_leaves(result), 4)
  // All leaves preserved
  let names : Array[String] = Array(capacity=0)
  mp_leaf_names(result, names)
  assert_eq(names.length(), 4)
}

///|
test "midpoint_root: asymmetric branch lengths" {
  // (A:5, B:1) — longest path = 6, midpoint at 3.0 from A
  // Midpoint is on the edge from root to A (at 3.0 from A, 2.0 from root)
  let tree = mp_internal0([mp_leaf("A", 5.0), mp_leaf("B", 1.0)])
  let result = @src.phylo_midpoint_root(tree)
  assert_eq(mp_count_leaves(result), 2)
  let names : Array[String] = Array(capacity=0)
  mp_leaf_names(result, names)
  assert_eq(names.length(), 2)
}

///|
test "midpoint_root: 5-leaf tree" {
  // (((A:1,B:1):1,C:1):2,(D:1,E:1):1)
  let tree = mp_internal0([
    mp_internal0([
      mp_internal0([mp_leaf("A", 1.0), mp_leaf("B", 1.0)]),
      mp_leaf("C", 1.0),
    ]),
    mp_internal0([mp_leaf("D", 1.0), mp_leaf("E", 1.0)]),
  ])
  let result = @src.phylo_midpoint_root(tree)
  assert_eq(mp_count_leaves(result), 5)
  let names : Array[String] = Array(capacity=0)
  mp_leaf_names(result, names)
  assert_eq(names.length(), 5)
}

///|
test "midpoint_root: all leaves preserved" {
  // Larger tree — make sure no leaves are lost
  let tree = mp_internal0([
    mp_internal0([
      mp_internal0([mp_leaf("A", 2.0), mp_leaf("B", 3.0)]),
      mp_leaf("C", 1.0),
    ]),
    mp_internal0([
      mp_internal0([mp_leaf("D", 1.0), mp_leaf("E", 2.0)]),
      mp_internal0([mp_leaf("F", 1.0), mp_leaf("G", 1.0)]),
    ]),
  ])
  let result = @src.phylo_midpoint_root(tree)
  assert_eq(mp_count_leaves(result), 7)
  let names : Array[String] = Array(capacity=0)
  mp_leaf_names(result, names)
  assert_eq(names.length(), 7)
}

///|
test "midpoint_root: star tree" {
  // (A:1, B:2, C:3) — longest path A-C = 4, midpoint at 2.0 from A
  let tree = mp_internal0([
    mp_leaf("A", 1.0), mp_leaf("B", 2.0), mp_leaf("C", 3.0),
  ])
  let result = @src.phylo_midpoint_root(tree)
  assert_eq(mp_count_leaves(result), 3)
}

///|
test "midpoint_root: balanced with uniform lengths" {
  // Perfectly balanced 8-leaf tree, all branch lengths = 1
  // Longest path = 4 (depth 2 on each side), midpoint at root
  let tree = mp_internal0([
    mp_internal0([
      mp_internal0([mp_leaf("A", 1.0), mp_leaf("B", 1.0)]),
      mp_internal0([mp_leaf("C", 1.0), mp_leaf("D", 1.0)]),
    ]),
    mp_internal0([
      mp_internal0([mp_leaf("E", 1.0), mp_leaf("F", 1.0)]),
      mp_internal0([mp_leaf("G", 1.0), mp_leaf("H", 1.0)]),
    ]),
  ])
  let result = @src.phylo_midpoint_root(tree)
  assert_eq(mp_count_leaves(result), 8)
  let names : Array[String] = Array(capacity=0)
  mp_leaf_names(result, names)
  assert_eq(names.length(), 8)
}
