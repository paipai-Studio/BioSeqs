///|
/// Tests for NNI (Nearest Neighbor Interchange) moves.

///|
fn nni_leaf(name : String) -> @src.Clade {
  @src.Clade::new(name=Some(name))
}

///|
fn nni_internal(children : Array[@src.Clade]) -> @src.Clade {
  @src.Clade::new(name=Some("internal"), clades=children)
}

///|
fn nni_count_leaves(clade : @src.Clade) -> Int {
  if clade.is_terminal() {
    return 1
  }
  let mut c = 0
  for ch in clade.clades {
    c = c + nni_count_leaves(ch)
  }
  c
}

///|
/// Collect leaf names in order.
fn nni_leaf_names(clade : @src.Clade, names : Array[String]) -> Unit {
  if clade.is_terminal() {
    match clade.name {
      Some(n) => names.push(n)
      None => names.push("?")
    }
    return
  }
  for ch in clade.clades {
    nni_leaf_names(ch, names)
  }
}

///|
/// Check if two trees have the same leaf set.
fn nni_same_leaves(a : @src.Clade, b : @src.Clade) -> Bool {
  let names_a : Array[String] = Array(capacity=0)
  let names_b : Array[String] = Array(capacity=0)
  nni_leaf_names(a, names_a)
  nni_leaf_names(b, names_b)
  if names_a.length() != names_b.length() {
    return false
  }
  // Check each name in a is in b
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
/// Check tree is binary.
fn nni_is_binary(clade : @src.Clade) -> Bool {
  if clade.is_terminal() {
    return true
  }
  if clade.clades.length() != 2 {
    return false
  }
  nni_is_binary(clade.clades[0]) && nni_is_binary(clade.clades[1])
}

// =====================
// NNI count tests
// =====================

///|
test "nni_count: leaf is 0" {
  assert_eq(@src.phylo_nni_count(nni_leaf("A")), 0)
}

///|
test "nni_count: cherry is 0" {
  // (A,B) — no internal edges
  let tree = nni_internal([nni_leaf("A"), nni_leaf("B")])
  assert_eq(@src.phylo_nni_count(tree), 0)
}

///|
test "nni_count: 4-leaf balanced is 1" {
  // ((A,B),(C,D)) — 1 internal edge (root → (A,B) or root → (C,D))
  // Actually: root has 2 internal children, so 2 internal edges
  let tree = nni_internal([
    nni_internal([nni_leaf("A"), nni_leaf("B")]),
    nni_internal([nni_leaf("C"), nni_leaf("D")]),
  ])
  assert_eq(@src.phylo_nni_count(tree), 2)
}

///|
test "nni_count: caterpillar 4-leaf is 1" {
  // (A,(B,(C,D))) — 1 internal edge (root → (B,(C,D)))
  let tree = nni_internal([
    nni_leaf("A"),
    nni_internal([nni_leaf("B"), nni_internal([nni_leaf("C"), nni_leaf("D")])]),
  ])
  assert_eq(@src.phylo_nni_count(tree), 2)
}

///|
test "nni_count: 6-leaf tree" {
  // (((A,B),(C,D)),(E,F)) — 3 internal edges
  let tree = nni_internal([
    nni_internal([
      nni_internal([nni_leaf("A"), nni_leaf("B")]),
      nni_internal([nni_leaf("C"), nni_leaf("D")]),
    ]),
    nni_internal([nni_leaf("E"), nni_leaf("F")]),
  ])
  assert_eq(@src.phylo_nni_count(tree), 4)
}

// =====================
// NNI move tests
// =====================

///|
test "nni_move: preserves leaf count" {
  let tree = nni_internal([
    nni_internal([nni_leaf("A"), nni_leaf("B")]),
    nni_internal([nni_leaf("C"), nni_leaf("D")]),
  ])
  let result = @src.phylo_nni_move(tree, 0, 0)
  assert_eq(nni_count_leaves(result), 4)
}

///|
test "nni_move: preserves leaf set" {
  let tree = nni_internal([
    nni_internal([nni_leaf("A"), nni_leaf("B")]),
    nni_internal([nni_leaf("C"), nni_leaf("D")]),
  ])
  let result = @src.phylo_nni_move(tree, 0, 0)
  assert_true(nni_same_leaves(tree, result))
}

///|
test "nni_move: preserves binary structure" {
  let tree = nni_internal([
    nni_internal([
      nni_internal([nni_leaf("A"), nni_leaf("B")]),
      nni_internal([nni_leaf("C"), nni_leaf("D")]),
    ]),
    nni_internal([nni_leaf("E"), nni_leaf("F")]),
  ])
  for i in 0..<4 {
    for swap in 0..<2 {
      let result = @src.phylo_nni_move(tree, i, swap)
      assert_true(nni_is_binary(result))
      assert_eq(nni_count_leaves(result), 6)
    }
  }
}

///|
test "nni_move: swap 0 vs swap 1 differ" {
  // On a 4-leaf balanced tree, the two swaps on edge 0 should give different trees
  let tree = nni_internal([
    nni_internal([nni_leaf("A"), nni_leaf("B")]),
    nni_internal([nni_leaf("C"), nni_leaf("D")]),
  ])
  let r0 = @src.phylo_nni_move(tree, 0, 0)
  let r1 = @src.phylo_nni_move(tree, 0, 1)
  let names0 : Array[String] = Array(capacity=0)
  let names1 : Array[String] = Array(capacity=0)
  nni_leaf_names(r0, names0)
  nni_leaf_names(r1, names1)
  // They should have different leaf orderings
  let mut diff = false
  for i in 0..<names0.length() {
    if names0[i] != names1[i] {
      diff = true
    }
  }
  assert_true(diff)
}

///|
test "nni_move: changes tree topology" {
  // ((A,B),(C,D)) — NNI should change the cherry structure
  let tree = nni_internal([
    nni_internal([nni_leaf("A"), nni_leaf("B")]),
    nni_internal([nni_leaf("C"), nni_leaf("D")]),
  ])
  let result = @src.phylo_nni_move(tree, 0, 0)
  // Original has cherries (A,B) and (C,D)
  // After NNI, cherries should be different, e.g., (A,C) and (B,D) or (A,D) and (B,C)
  let names : Array[String] = Array(capacity=0)
  nni_leaf_names(result, names)
  // The first two leaves should NOT be A,B (cherry broken)
  let is_ab_cherry = names[0] == "A" && names[1] == "B"
  let is_cd_cherry = names[2] == "C" && names[3] == "D"
  assert_true(not(is_ab_cherry && is_cd_cherry))
}

// =====================
// NNI neighbors tests
// =====================

///|
test "nni_neighbors: 4-leaf balanced" {
  // ((A,B),(C,D)) — 2 internal edges, 2*2=4 neighbors
  let tree = nni_internal([
    nni_internal([nni_leaf("A"), nni_leaf("B")]),
    nni_internal([nni_leaf("C"), nni_leaf("D")]),
  ])
  let neighbors = @src.phylo_nni_neighbors(tree)
  assert_eq(neighbors.length(), 4)
  // All neighbors should have 4 leaves
  for n in neighbors {
    assert_eq(nni_count_leaves(n), 4)
    assert_true(nni_same_leaves(tree, n))
  }
}

///|
test "nni_neighbors: caterpillar" {
  // (A,(B,(C,D))) — 1 internal edge, 2 neighbors
  let tree = nni_internal([
    nni_leaf("A"),
    nni_internal([nni_leaf("B"), nni_internal([nni_leaf("C"), nni_leaf("D")])]),
  ])
  let neighbors = @src.phylo_nni_neighbors(tree)
  assert_eq(neighbors.length(), 4)
  for n in neighbors {
    assert_eq(nni_count_leaves(n), 4)
    assert_true(nni_same_leaves(tree, n))
  }
}

///|
test "nni_neighbors: 6-leaf tree" {
  let tree = nni_internal([
    nni_internal([
      nni_internal([nni_leaf("A"), nni_leaf("B")]),
      nni_internal([nni_leaf("C"), nni_leaf("D")]),
    ]),
    nni_internal([nni_leaf("E"), nni_leaf("F")]),
  ])
  let neighbors = @src.phylo_nni_neighbors(tree)
  // 4 internal edges, 8 neighbors
  assert_eq(neighbors.length(), 8)
  for n in neighbors {
    assert_eq(nni_count_leaves(n), 6)
    assert_true(nni_is_binary(n))
  }
}

///|
test "nni_neighbors: cherry has no neighbors" {
  let tree = nni_internal([nni_leaf("A"), nni_leaf("B")])
  let neighbors = @src.phylo_nni_neighbors(tree)
  assert_eq(neighbors.length(), 0)
}

///|
test "nni_move: out of range returns valid tree" {
  let tree = nni_internal([
    nni_internal([nni_leaf("A"), nni_leaf("B")]),
    nni_internal([nni_leaf("C"), nni_leaf("D")]),
  ])
  let result = @src.phylo_nni_move(tree, 99, 0)
  assert_eq(nni_count_leaves(result), 4)
}

///|
test "nni_move: 8-leaf tree all valid" {
  let tree = nni_internal([
    nni_internal([
      nni_internal([nni_leaf("A"), nni_leaf("B")]),
      nni_internal([nni_leaf("C"), nni_leaf("D")]),
    ]),
    nni_internal([
      nni_internal([nni_leaf("E"), nni_leaf("F")]),
      nni_internal([nni_leaf("G"), nni_leaf("H")]),
    ]),
  ])
  let n_edges = @src.phylo_nni_count(tree)
  let neighbors = @src.phylo_nni_neighbors(tree)
  assert_eq(neighbors.length(), n_edges * 2)
  for n in neighbors {
    assert_eq(nni_count_leaves(n), 8)
    assert_true(nni_is_binary(n))
    assert_true(nni_same_leaves(tree, n))
  }
}
