///|
/// Tests for root/unroot operations.

///|
fn ru_leaf(name : String, bl : Double) -> @src.Clade {
  @src.Clade::new(name=Some(name), branch_length=Some(bl))
}

///|
fn ru_leaf0(name : String) -> @src.Clade {
  @src.Clade::new(name=Some(name))
}

///|
fn ru_internal(children : Array[@src.Clade], bl : Double) -> @src.Clade {
  @src.Clade::new(name=Some("internal"), branch_length=Some(bl), clades=children)
}

///|
fn ru_internal0(children : Array[@src.Clade]) -> @src.Clade {
  @src.Clade::new(name=Some("internal"), clades=children)
}

///|
fn ru_count_leaves(clade : @src.Clade) -> Int {
  if clade.is_terminal() {
    return 1
  }
  let mut c = 0
  for ch in clade.clades {
    c = c + ru_count_leaves(ch)
  }
  c
}

///|
fn ru_leaf_names(clade : @src.Clade, names : Array[String]) -> Unit {
  if clade.is_terminal() {
    match clade.name {
      Some(n) => names.push(n)
      None => names.push("?")
    }
    return
  }
  for ch in clade.clades {
    ru_leaf_names(ch, names)
  }
}

///|
fn ru_has_leaf(clade : @src.Clade, name : String) -> Bool {
  if clade.is_terminal() {
    match clade.name {
      Some(n) => return n == name
      None => return false
    }
  }
  for ch in clade.clades {
    if ru_has_leaf(ch, name) {
      return true
    }
  }
  false
}

// =====================
// Unroot tests
// =====================

///|
test "unroot: leaf unchanged" {
  let leaf = ru_leaf0("A")
  let result = @src.phylo_unroot(leaf)
  assert_true(result.is_terminal())
}

///|
test "unroot: two-leaf tree unchanged" {
  // (A, B) — both children are leaves, can't unroot meaningfully
  let tree = ru_internal0([ru_leaf("A", 1.0), ru_leaf("B", 1.0)])
  let result = @src.phylo_unroot(tree)
  assert_eq(ru_count_leaves(result), 2)
}

///|
test "unroot: root with 3 children unchanged" {
  // Already unrooted (trifurcation at root)
  let tree = ru_internal0([ru_leaf("A", 1.0), ru_leaf("B", 1.0), ru_leaf("C", 1.0)])
  let result = @src.phylo_unroot(tree)
  assert_eq(result.clades.length(), 3)
}

///|
test "unroot: basic bifurcation" {
  // ((A,B),(C,D)) — root has 2 internal children
  // After unroot: one child becomes root with 3 children (A, B, other)
  let tree = ru_internal0([
    ru_internal0([ru_leaf("A", 1.0), ru_leaf("B", 1.0)]),
    ru_internal0([ru_leaf("C", 1.0), ru_leaf("D", 1.0)]),
  ])
  let result = @src.phylo_unroot(tree)
  assert_eq(ru_count_leaves(result), 4)
  // Root should now have 3 children (not 2)
  assert_eq(result.clades.length(), 3)
}

///|
test "unroot: preserves all leaves" {
  let tree = ru_internal0([
    ru_internal0([ru_leaf("A", 2.0), ru_leaf("B", 3.0)]),
    ru_internal0([ru_leaf("C", 1.0), ru_leaf("D", 1.0)]),
  ])
  let result = @src.phylo_unroot(tree)
  assert_eq(ru_count_leaves(result), 4)
  let names : Array[String] = Array(capacity=0)
  ru_leaf_names(result, names)
  assert_eq(names.length(), 4)
}

///|
test "unroot: combined branch length" {
  // ((A:1,B:1):2,(C:1,D:1):3)
  // After unroot: first child becomes root, second attached with bl=2+3=5
  let tree = ru_internal0([
    ru_internal([ru_leaf("A", 1.0), ru_leaf("B", 1.0)], 2.0),
    ru_internal([ru_leaf("C", 1.0), ru_leaf("D", 1.0)], 3.0),
  ])
  let result = @src.phylo_unroot(tree)
  assert_eq(ru_count_leaves(result), 4)
  // The merged child should have branch_length = 5.0
  let mut found_combined = false
  for ch in result.clades {
    match ch.branch_length {
      Some(v) => if (v - 5.0).abs() < 0.0001 {
        found_combined = true
      }
      None => ()
    }
  }
  assert_true(found_combined)
}

///|
test "unroot: caterpillar" {
  // (A,(B,(C,D))) — root has leaf A and internal (B,(C,D))
  // After unroot: (B,(C,D)) becomes root, A attached as sibling
  let tree = ru_internal0([
    ru_leaf("A", 1.0),
    ru_internal0([ru_leaf("B", 1.0), ru_internal0([ru_leaf("C", 1.0), ru_leaf("D", 1.0)])]),
  ])
  let result = @src.phylo_unroot(tree)
  assert_eq(ru_count_leaves(result), 4)
  // Root should have 3 children: B, (C,D), A
  assert_eq(result.clades.length(), 3)
}

// =====================
// Root at leaf tests
// =====================

///|
test "root_at_leaf: preserves leaf count" {
  let tree = ru_internal0([
    ru_internal0([ru_leaf("A", 1.0), ru_leaf("B", 1.0)]),
    ru_internal0([ru_leaf("C", 1.0), ru_leaf("D", 1.0)]),
  ])
  let result = @src.phylo_root_at_leaf(tree, "A")
  assert_eq(ru_count_leaves(result), 4)
}

///|
test "root_at_leaf: outgroup becomes direct child" {
  let tree = ru_internal0([
    ru_internal0([ru_leaf("A", 1.0), ru_leaf("B", 1.0)]),
    ru_internal0([ru_leaf("C", 1.0), ru_leaf("D", 1.0)]),
  ])
  let result = @src.phylo_root_at_leaf(tree, "C")
  // C should be a direct child of the new root
  assert_eq(result.clades.length(), 2)
  let mut c_is_direct = false
  for ch in result.clades {
    if ch.is_terminal() {
      match ch.name {
        Some(n) => if n == "C" {
          c_is_direct = true
        }
        None => ()
      }
    }
  }
  assert_true(c_is_direct)
}

///|
test "root_at_leaf: all leaves preserved" {
  let tree = ru_internal0([
    ru_internal0([
      ru_internal0([ru_leaf("A", 2.0), ru_leaf("B", 3.0)]),
      ru_leaf("C", 1.0),
    ]),
    ru_internal0([ru_leaf("D", 1.0), ru_leaf("E", 2.0)]),
  ])
  let result = @src.phylo_root_at_leaf(tree, "E")
  assert_eq(ru_count_leaves(result), 5)
  assert_true(ru_has_leaf(result, "A"))
  assert_true(ru_has_leaf(result, "B"))
  assert_true(ru_has_leaf(result, "C"))
  assert_true(ru_has_leaf(result, "D"))
  assert_true(ru_has_leaf(result, "E"))
}

///|
test "root_at_leaf: leaf not found returns original" {
  let tree = ru_internal0([ru_leaf("A", 1.0), ru_leaf("B", 1.0)])
  let result = @src.phylo_root_at_leaf(tree, "X")
  assert_eq(ru_count_leaves(result), 2)
}

///|
test "root_at_leaf: star tree" {
  // (A, B, C, D) — root at C
  let tree = ru_internal0([
    ru_leaf("A", 1.0), ru_leaf("B", 2.0), ru_leaf("C", 3.0), ru_leaf("D", 1.0),
  ])
  let result = @src.phylo_root_at_leaf(tree, "C")
  assert_eq(ru_count_leaves(result), 4)
  assert_eq(result.clades.length(), 2)
}

///|
test "root_at_leaf: caterpillar" {
  // (A,(B,(C,D))) — root at D
  let tree = ru_internal0([
    ru_leaf("A", 1.0),
    ru_internal0([ru_leaf("B", 1.0), ru_internal0([ru_leaf("C", 1.0), ru_leaf("D", 1.0)])]),
  ])
  let result = @src.phylo_root_at_leaf(tree, "D")
  assert_eq(ru_count_leaves(result), 4)
  assert_eq(result.clades.length(), 2)
  // D should be a direct child of root
  let mut d_is_direct = false
  for ch in result.clades {
    if ch.is_terminal() {
      match ch.name {
        Some(n) => if n == "D" {
          d_is_direct = true
        }
        None => ()
      }
    }
  }
  assert_true(d_is_direct)
}

///|
test "root_unroot_roundtrip: root then unroot" {
  // Start with unrooted (trifurcation), root at A, then unroot
  let tree = ru_internal0([
    ru_leaf("A", 1.0), ru_leaf("B", 1.0), ru_leaf("C", 1.0),
  ])
  let rooted = @src.phylo_root_at_leaf(tree, "A")
  let unrooted = @src.phylo_unroot(rooted)
  assert_eq(ru_count_leaves(unrooted), 3)
  // After unrooting, root should have 3 children again (trifurcation)
  // (or at least all 3 leaves preserved)
  assert_true(ru_has_leaf(unrooted, "A"))
  assert_true(ru_has_leaf(unrooted, "B"))
  assert_true(ru_has_leaf(unrooted, "C"))
}
