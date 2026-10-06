///|
/// Tests for s-shape statistic and normalized Colless/Sackin indices.

///|
/// Build a perfectly balanced binary tree with n leaves.
fn shape_make_balanced(n : Int) -> @src.Clade {
  if n <= 1 {
    return @src.Clade::new(name=Some("leaf"))
  }
  let half = n / 2
  let left = shape_make_balanced(half)
  let right = shape_make_balanced(n - half)
  @src.Clade::new(name=Some("internal"), clades=[left, right])
}

///|
/// Build a maximally imbalanced (caterpillar) tree with n leaves.
fn shape_make_caterpillar(n : Int) -> @src.Clade {
  if n <= 1 {
    return @src.Clade::new(name=Some("leaf"))
  }
  if n == 2 {
    return @src.Clade::new(
      name=Some("internal"),
      clades=[@src.Clade::new(name=Some("leaf")), @src.Clade::new(name=Some("leaf"))],
    )
  }
  @src.Clade::new(
    name=Some("internal"),
    clades=[
      @src.Clade::new(name=Some("leaf")),
      shape_make_caterpillar(n - 1),
    ],
  )
}

///|
test "shape_statistic: single leaf" {
  let leaf = @src.Clade::new(name=Some("A"))
  assert_eq(@src.dendropy_shape_statistic(leaf), 0.0)
}

///|
test "shape_statistic: two leaves (cherry)" {
  // One internal node with 2 leaves: log(2-1) = log(1) = 0
  let tree = shape_make_balanced(2)
  let s = @src.dendropy_shape_statistic(tree)
  assert_true((s - 0.0).abs() < 0.000000000001)
}

///|
test "shape_statistic: balanced 4-leaf tree" {
  // ((A,B),(C,D)) — 3 internal nodes
  // Root: n=4, log(3); two children: n=2 each, log(1)=0
  // Total = log(3) + 0 + 0 = log(3)
  let tree = shape_make_balanced(4)
  let s = @src.dendropy_shape_statistic(tree)
  let expected = @math.ln(3.0)
  assert_true((s - expected).abs() < 0.0000000001)
}

///|
test "shape_statistic: caterpillar 4-leaf tree" {
  // (A,(B,(C,D))) — 3 internal nodes
  // Root: n=4, log(3); next: n=3, log(2); next: n=2, log(1)=0
  // Total = log(3) + log(2) + 0 = log(6)
  let tree = shape_make_caterpillar(4)
  let s = @src.dendropy_shape_statistic(tree)
  let expected = @math.ln(3.0) + @math.ln(2.0)
  assert_true((s - expected).abs() < 0.0000000001)
}

///|
test "shape_statistic: balanced vs caterpillar" {
  // Balanced tree should have lower s-shape (more balanced)
  let balanced = shape_make_balanced(8)
  let caterpillar = shape_make_caterpillar(8)
  let s_bal = @src.dendropy_shape_statistic(balanced)
  let s_cat = @src.dendropy_shape_statistic(caterpillar)
  assert_true(s_cat > s_bal)
}

///|
test "shape_statistic_yule: balanced 8-leaf" {
  let tree = shape_make_balanced(8)
  let s = @src.dendropy_shape_statistic_yule(tree)
  // Should be a finite number (negative = more balanced than Yule)
  assert_true(s < 0.0) // balanced tree should be negative
}

///|
test "shape_statistic_yule: caterpillar 8-leaf" {
  let tree = shape_make_caterpillar(8)
  let s = @src.dendropy_shape_statistic_yule(tree)
  // Caterpillar is highly imbalanced, should be positive
  assert_true(s > 0.0)
}

///|
test "shape_statistic_pda: balanced 8-leaf" {
  let tree = shape_make_balanced(8)
  let s = @src.dendropy_shape_statistic_pda(tree)
  // Should be a finite number
  assert_true(s < 0.0) // balanced tree should be negative under PDA too
}

///|
test "shape_statistic_pda: caterpillar 8-leaf" {
  let tree = shape_make_caterpillar(8)
  let s = @src.dendropy_shape_statistic_pda(tree)
  assert_true(s > 0.0) // imbalanced tree should be positive
}

///|
test "colless_yule: balanced tree near zero" {
  // Perfectly balanced tree should have Colless = 0, so yule norm ≈ negative
  let tree = shape_make_balanced(8)
  let c = @src.dendropy_colless_yule(tree)
  // Colless raw = 0, so normalized = -(n*log(n) + n*(gamma-1-log2))/n
  assert_true(c < 0.0)
}

///|
test "colless_yule: caterpillar positive" {
  let tree = shape_make_caterpillar(8)
  let c = @src.dendropy_colless_yule(tree)
  // Caterpillar has max Colless, should be positive
  assert_true(c > 0.0)
}

///|
test "colless_pda: balanced zero" {
  let tree = shape_make_balanced(8)
  let c = @src.dendropy_colless_pda(tree)
  // Colless raw = 0 for balanced
  assert_true((c - 0.0).abs() < 0.0000000001)
}

///|
test "colless_pda: caterpillar positive" {
  let tree = shape_make_caterpillar(8)
  let c = @src.dendropy_colless_pda(tree)
  assert_true(c > 0.0)
}

///|
test "sackin_yule: balanced tree" {
  let tree = shape_make_balanced(8)
  let s = @src.dendropy_sackin_yule(tree)
  // Balanced tree should have negative normalized Sackin
  assert_true(s < 0.0)
}

///|
test "sackin_yule: caterpillar positive" {
  let tree = shape_make_caterpillar(8)
  let s = @src.dendropy_sackin_yule(tree)
  // Caterpillar should have positive normalized Sackin
  assert_true(s > 0.0)
}

///|
test "sackin_pda: balanced vs caterpillar" {
  let balanced = shape_make_balanced(8)
  let caterpillar = shape_make_caterpillar(8)
  let s_bal = @src.dendropy_sackin_pda(balanced)
  let s_cat = @src.dendropy_sackin_pda(caterpillar)
  assert_true(s_cat > s_bal)
}

///|
test "sackin_pda: positive values" {
  let tree = shape_make_balanced(4)
  let s = @src.dendropy_sackin_pda(tree)
  assert_true(s > 0.0)
}

///|
test "colless_yule: matches manual calculation" {
  // 4-leaf caterpillar: Colless = 2, n = 4
  // Ic' = (2 - 4*log(4) - 4*(gamma - 1 - log(2))) / 4
  let tree = shape_make_caterpillar(4)
  let c = @src.dendropy_colless_yule(tree)
  let gamma = 0.5772156649015329
  let expected = (3.0 - 4.0 * @math.ln(4.0) - 4.0 * (gamma - 1.0 - @math.ln(2.0))) / 4.0
  assert_true((c - expected).abs() < 0.0000000001)
}

///|
test "sackin_yule: matches manual calculation" {
  // 4-leaf balanced: Sackin = 8 (each leaf at depth 2), n = 4
  // Is' = (8 - 2*4*(1/2 + 1/3 + 1/4)) / 4
  let tree = shape_make_balanced(4)
  let s = @src.dendropy_sackin_yule(tree)
  let h = 1.0 / 2.0 + 1.0 / 3.0 + 1.0 / 4.0
  let expected = (8.0 - 2.0 * 4.0 * h) / 4.0
  assert_true((s - expected).abs() < 0.0000000001)
}
