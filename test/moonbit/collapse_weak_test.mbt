///|
/// Tests for collapse weak nodes.

///|
fn cw_leaf(name : String) -> @src.Clade {
  @src.Clade::new(name=Some(name))
}

///|
fn cw_internal(children : Array[@src.Clade]) -> @src.Clade {
  @src.Clade::new(name=Some("internal"), clades=children)
}

///|
fn cw_internal_conf(children : Array[@src.Clade], conf : Double) -> @src.Clade {
  @src.Clade::new(name=Some("internal"), confidence=Some(conf), clades=children)
}

///|
fn cw_internal_bl(children : Array[@src.Clade], bl : Double) -> @src.Clade {
  @src.Clade::new(name=Some("internal"), branch_length=Some(bl), clades=children)
}

///|
fn cw_count_leaves(clade : @src.Clade) -> Int {
  if clade.is_terminal() {
    return 1
  }
  let mut c = 0
  for ch in clade.clades {
    c = c + cw_count_leaves(ch)
  }
  c
}

///|
/// Count internal nodes (non-leaf nodes).
fn cw_count_internal(clade : @src.Clade) -> Int {
  if clade.is_terminal() {
    return 0
  }
  let mut c = 1
  for ch in clade.clades {
    c = c + cw_count_internal(ch)
  }
  c
}

///|
fn cw_has_leaf(clade : @src.Clade, name : String) -> Bool {
  if clade.is_terminal() {
    match clade.name {
      Some(n) => return n == name
      None => return false
    }
  }
  for ch in clade.clades {
    if cw_has_leaf(ch, name) {
      return true
    }
  }
  false
}

///|
test "collapse: leaf unchanged" {
  let leaf = cw_leaf("A")
  let result = @src.phylo_collapse_weak_nodes(leaf, 50.0)
  assert_true(result.is_terminal())
}

///|
test "collapse: no collapse when all support high" {
  // All nodes have support >= 80, threshold = 50 → no collapse
  let tree = cw_internal([
    cw_internal_conf([cw_leaf("A"), cw_leaf("B")], 90.0),
    cw_internal_conf([cw_leaf("C"), cw_leaf("D")], 80.0),
  ])
  let result = @src.phylo_collapse_weak_nodes(tree, 50.0)
  assert_eq(cw_count_internal(result), 3)
  assert_eq(cw_count_leaves(result), 4)
}

///|
test "collapse: collapse low support node" {
  // ((A,B):30, (C,D):90) — threshold 50, first node collapses
  let tree = cw_internal([
    cw_internal_conf([cw_leaf("A"), cw_leaf("B")], 30.0),
    cw_internal_conf([cw_leaf("C"), cw_leaf("D")], 90.0),
  ])
  let result = @src.phylo_collapse_weak_nodes(tree, 50.0)
  assert_eq(cw_count_leaves(result), 4)
  // The low-support node is collapsed: A, B become direct children of root
  // Root should now have 3 children: A, B, (C,D)
  assert_eq(result.clades.length(), 3)
}

///|
test "collapse: collapse all low support" {
  // ((A,B):30, (C,D):40) — threshold 50, both collapse
  let tree = cw_internal([
    cw_internal_conf([cw_leaf("A"), cw_leaf("B")], 30.0),
    cw_internal_conf([cw_leaf("C"), cw_leaf("D")], 40.0),
  ])
  let result = @src.phylo_collapse_weak_nodes(tree, 50.0)
  assert_eq(cw_count_leaves(result), 4)
  // Both internal nodes collapsed: root has 4 children (star tree)
  assert_eq(result.clades.length(), 4)
}

///|
test "collapse: preserves all leaves" {
  let tree = cw_internal([
    cw_internal_conf([
      cw_internal_conf([cw_leaf("A"), cw_leaf("B")], 20.0),
      cw_leaf("C"),
    ], 90.0),
    cw_internal_conf([cw_leaf("D"), cw_leaf("E")], 30.0),
  ])
  let result = @src.phylo_collapse_weak_nodes(tree, 50.0)
  assert_eq(cw_count_leaves(result), 5)
  assert_true(cw_has_leaf(result, "A"))
  assert_true(cw_has_leaf(result, "B"))
  assert_true(cw_has_leaf(result, "C"))
  assert_true(cw_has_leaf(result, "D"))
  assert_true(cw_has_leaf(result, "E"))
}

///|
test "collapse: no confidence kept" {
  // Node without confidence should NOT be collapsed
  let tree = cw_internal([
    cw_internal([cw_leaf("A"), cw_leaf("B")]), // no confidence
    cw_internal_conf([cw_leaf("C"), cw_leaf("D")], 90.0),
  ])
  let result = @src.phylo_collapse_weak_nodes(tree, 50.0)
  assert_eq(cw_count_internal(result), 3)
}

///|
test "collapse: threshold zero collapses nothing" {
  let tree = cw_internal([
    cw_internal_conf([cw_leaf("A"), cw_leaf("B")], 10.0),
    cw_internal_conf([cw_leaf("C"), cw_leaf("D")], 20.0),
  ])
  let result = @src.phylo_collapse_weak_nodes(tree, 0.0)
  assert_eq(cw_count_internal(result), 3)
}

///|
test "collapse: threshold very high collapses all" {
  let tree = cw_internal([
    cw_internal_conf([cw_leaf("A"), cw_leaf("B")], 90.0),
    cw_internal_conf([cw_leaf("C"), cw_leaf("D")], 95.0),
  ])
  let result = @src.phylo_collapse_weak_nodes(tree, 100.0)
  assert_eq(cw_count_leaves(result), 4)
  // All internal nodes collapsed: star tree with 4 leaves
  assert_eq(result.clades.length(), 4)
}

///|
test "collapse: nested collapse" {
  // (((A,B):20, C):30, D) — threshold 50
  // Both internal nodes with support < 50 collapse
  // Result: (A, B, C, D) — star tree
  let tree = cw_internal([
    cw_internal_conf([
      cw_internal_conf([cw_leaf("A"), cw_leaf("B")], 20.0),
      cw_leaf("C"),
    ], 30.0),
    cw_leaf("D"),
  ])
  let result = @src.phylo_collapse_weak_nodes(tree, 50.0)
  assert_eq(cw_count_leaves(result), 4)
  assert_eq(result.clades.length(), 4)
}

///|
test "collapse: partial nested collapse" {
  // (((A,B):20, C):90, D) — threshold 50
  // Only the (A,B) node collapses; ((A,B),C) node stays (support 90)
  // Result: ((A,B,C):90, D) — root has 2 children, first has 3 children
  let tree = cw_internal([
    cw_internal_conf([
      cw_internal_conf([cw_leaf("A"), cw_leaf("B")], 20.0),
      cw_leaf("C"),
    ], 90.0),
    cw_leaf("D"),
  ])
  let result = @src.phylo_collapse_weak_nodes(tree, 50.0)
  assert_eq(cw_count_leaves(result), 4)
  // Root has 2 children: (A,B,C) and D
  assert_eq(result.clades.length(), 2)
  // First child should have 3 children (A, B, C)
  assert_true(not(result.clades[0].is_terminal()))
  assert_eq(result.clades[0].clades.length(), 3)
}

///|
test "collapse: branch lengths merged" {
  // ((A:1,B:1):2, C:3) with support 30, threshold 50
  // After collapse: A gets bl=1+2=3, B gets bl=1+2=3
  let tree = cw_internal([
    cw_internal_conf([
      @src.Clade::new(name=Some("A"), branch_length=Some(1.0)),
      @src.Clade::new(name=Some("B"), branch_length=Some(1.0)),
    ], 30.0),
  ])
  // Add branch_length to the internal node
  let tree_with_bl = @src.Clade::new(
    name=Some("root"),
    clades=[
      @src.Clade::new(
        name=Some("internal"),
        confidence=Some(30.0),
        branch_length=Some(2.0),
        clades=[
          @src.Clade::new(name=Some("A"), branch_length=Some(1.0)),
          @src.Clade::new(name=Some("B"), branch_length=Some(1.0)),
        ],
      ),
      @src.Clade::new(name=Some("C"), branch_length=Some(3.0)),
    ],
  )
  let result = @src.phylo_collapse_weak_nodes(tree_with_bl, 50.0)
  assert_eq(cw_count_leaves(result), 3)
  // Root should have 3 children: A (bl=3), B (bl=3), C (bl=3)
  assert_eq(result.clades.length(), 3)
  let mut found_merged = false
  for ch in result.clades {
    match (ch.name, ch.branch_length) {
      (Some(n), Some(bl)) => {
        if (n == "A" || n == "B") && (bl - 3.0).abs() < 0.0001 {
          found_merged = true
        }
      }
      _ => ()
    }
  }
  assert_true(found_merged)
}
