///|
/// Tests for cherry count.

///|
fn cherry_leaf(name : String) -> @src.Clade {
  @src.Clade::new(name=Some(name))
}

///|
fn cherry_internal(children : Array[@src.Clade]) -> @src.Clade {
  @src.Clade::new(name=Some("internal"), clades=children)
}

///|
fn cherry_balanced(n : Int) -> @src.Clade {
  if n <= 1 {
    return cherry_leaf("L")
  }
  let half = n / 2
  cherry_internal([cherry_balanced(half), cherry_balanced(n - half)])
}

///|
fn cherry_caterpillar(n : Int) -> @src.Clade {
  if n <= 1 {
    return cherry_leaf("L")
  }
  if n == 2 {
    return cherry_internal([cherry_leaf("L"), cherry_leaf("L")])
  }
  cherry_internal([cherry_leaf("L"), cherry_caterpillar(n - 1)])
}

///|
test "cherry_count: single leaf" {
  let leaf = cherry_leaf("A")
  assert_eq(@src.dendropy_cherry_count(leaf), 0)
}

///|
test "cherry_count: two leaves (one cherry)" {
  let tree = cherry_internal([cherry_leaf("A"), cherry_leaf("B")])
  assert_eq(@src.dendropy_cherry_count(tree), 1)
}

///|
test "cherry_count: balanced 4-leaf (two cherries)" {
  // ((A,B),(C,D)) — two cherries
  let tree = cherry_balanced(4)
  assert_eq(@src.dendropy_cherry_count(tree), 2)
}

///|
test "cherry_count: balanced 8-leaf (four cherries)" {
  // Perfectly balanced 8-leaf tree has 4 cherries
  let tree = cherry_balanced(8)
  assert_eq(@src.dendropy_cherry_count(tree), 4)
}

///|
test "cherry_count: caterpillar (one cherry)" {
  // Caterpillar always has exactly 1 cherry (the deepest pair)
  let tree = cherry_caterpillar(8)
  assert_eq(@src.dendropy_cherry_count(tree), 1)
}

///|
test "cherry_count: caterpillar 4-leaf" {
  let tree = cherry_caterpillar(4)
  assert_eq(@src.dendropy_cherry_count(tree), 1)
}

///|
test "cherry_count: star tree (no cherries)" {
  // Root with 3 leaves — no internal node has exactly 2 leaf children
  let tree = cherry_internal([cherry_leaf("A"), cherry_leaf("B"), cherry_leaf("C")])
  assert_eq(@src.dendropy_cherry_count(tree), 0)
}

///|
test "cherry_count: mixed tree" {
  // ((A,B),(C,(D,E))) — cherries: (A,B) and (D,E) = 2
  let tree = cherry_internal([
    cherry_internal([cherry_leaf("A"), cherry_leaf("B")]),
    cherry_internal([cherry_leaf("C"), cherry_internal([cherry_leaf("D"), cherry_leaf("E")])]),
  ])
  assert_eq(@src.dendropy_cherry_count(tree), 2)
}

///|
test "cherry_count: three cherries" {
  // (((A,B),(C,D)),(E,F)) — three cherries
  let tree = cherry_internal([
    cherry_internal([
      cherry_internal([cherry_leaf("A"), cherry_leaf("B")]),
      cherry_internal([cherry_leaf("C"), cherry_leaf("D")]),
    ]),
    cherry_internal([cherry_leaf("E"), cherry_leaf("F")]),
  ])
  assert_eq(@src.dendropy_cherry_count(tree), 3)
}

///|
test "cherry_count: no cherry in multifurcating" {
  // (A, (B, C, D)) — the inner node has 3 leaf children, not a cherry
  let tree = cherry_internal([
    cherry_leaf("A"),
    cherry_internal([cherry_leaf("B"), cherry_leaf("C"), cherry_leaf("D")]),
  ])
  assert_eq(@src.dendropy_cherry_count(tree), 0)
}

///|
test "cherry_count: balanced 16-leaf" {
  let tree = cherry_balanced(16)
  assert_eq(@src.dendropy_cherry_count(tree), 8)
}
