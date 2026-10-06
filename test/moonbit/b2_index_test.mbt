///|
/// Tests for B2 index (Shao & Sokal 1990).

///|
fn b2_leaf(name : String) -> @src.Clade {
  @src.Clade::new(name=Some(name))
}

///|
fn b2_internal(children : Array[@src.Clade]) -> @src.Clade {
  @src.Clade::new(name=Some("internal"), clades=children)
}

///|
fn b2_balanced(n : Int) -> @src.Clade {
  if n <= 1 {
    return b2_leaf("L")
  }
  let half = n / 2
  b2_internal([b2_balanced(half), b2_balanced(n - half)])
}

///|
fn b2_caterpillar(n : Int) -> @src.Clade {
  if n <= 1 {
    return b2_leaf("L")
  }
  if n == 2 {
    return b2_internal([b2_leaf("L"), b2_leaf("L")])
  }
  b2_internal([b2_leaf("L"), b2_caterpillar(n - 1)])
}

///|
test "b2: single leaf" {
  let leaf = b2_leaf("A")
  assert_true((@src.dendropy_b2(leaf) - 0.0).abs() < 0.0000000001)
}

///|
test "b2: two leaves (cherry)" {
  // Both leaves have p=0.5, B2 = -2*(0.5*log2(0.5)) = -2*(0.5*(-1)) = 1
  let tree = b2_internal([b2_leaf("A"), b2_leaf("B")])
  let b2 = @src.dendropy_b2(tree)
  assert_true((b2 - 1.0).abs() < 0.0000000001)
}

///|
test "b2: balanced 4-leaf" {
  // ((A,B),(C,D)): each leaf p=0.25, B2 = -4*(0.25*log2(0.25)) = -4*(0.25*(-2)) = 2
  let tree = b2_balanced(4)
  let b2 = @src.dendropy_b2(tree)
  assert_true((b2 - 2.0).abs() < 0.0000000001)
}

///|
test "b2: balanced 8-leaf" {
  // Perfectly balanced 8-leaf: each p=1/8, B2 = -8*(1/8*log2(1/8)) = 3
  let tree = b2_balanced(8)
  let b2 = @src.dendropy_b2(tree)
  assert_true((b2 - 3.0).abs() < 0.0000000001)
}

///|
test "b2: caterpillar has lower B2 than balanced" {
  // B2 is a balance index — higher = more balanced
  let balanced = b2_balanced(8)
  let caterpillar = b2_caterpillar(8)
  let b2_bal = @src.dendropy_b2(balanced)
  let b2_cat = @src.dendropy_b2(caterpillar)
  assert_true(b2_bal > b2_cat)
}

///|
test "b2: caterpillar 4-leaf manual" {
  // (A,(B,(C,D))): probabilities:
  // A: 1/2, B: 1/4, C: 1/8, D: 1/8
  // B2 = -(0.5*log2(0.5) + 0.25*log2(0.25) + 2*(0.125*log2(0.125)))
  //    = -(0.5*(-1) + 0.25*(-2) + 2*(0.125*(-3)))
  //    = -(-0.5 - 0.5 - 0.75) = 1.75
  let tree = b2_caterpillar(4)
  let b2 = @src.dendropy_b2(tree)
  assert_true((b2 - 1.75).abs() < 0.0000000001)
}

///|
test "b2: star tree 3 leaves" {
  // (A,B,C): each p=1/3, B2 = -3*(1/3*log2(1/3)) = log2(3)
  let tree = b2_internal([b2_leaf("A"), b2_leaf("B"), b2_leaf("C")])
  let b2 = @src.dendropy_b2(tree)
  let expected = @math.ln(3.0) / @math.ln(2.0)
  assert_true((b2 - expected).abs() < 0.0000000001)
}

///|
test "b2: natural log base" {
  // Same tree, base e: B2_e = B2_2 * ln(2)
  let tree = b2_caterpillar(4)
  let b2_2 = @src.dendropy_b2(tree)
  let b2_e = @src.dendropy_b2_base(tree, 2.718281828459045)
  assert_true((b2_e - b2_2 * @math.ln(2.0)).abs() < 0.0001)
}

///|
test "b2: balanced vs caterpillar 16-leaf" {
  let balanced = b2_balanced(16)
  let caterpillar = b2_caterpillar(16)
  let b2_bal = @src.dendropy_b2(balanced)
  let b2_cat = @src.dendropy_b2(caterpillar)
  // Balanced 16-leaf: B2 = 4
  assert_true((b2_bal - 4.0).abs() < 0.0000000001)
  assert_true(b2_bal > b2_cat)
}

///|
test "b2: multifurcating balanced" {
  // ((A,B,C),(D,E,F)): each leaf p=1/6, B2 = -6*(1/6*log2(1/6)) = log2(6)
  let tree = b2_internal([
    b2_internal([b2_leaf("A"), b2_leaf("B"), b2_leaf("C")]),
    b2_internal([b2_leaf("D"), b2_leaf("E"), b2_leaf("F")]),
  ])
  let b2 = @src.dendropy_b2(tree)
  let expected = @math.ln(6.0) / @math.ln(2.0)
  assert_true((b2 - expected).abs() < 0.0000000001)
}
