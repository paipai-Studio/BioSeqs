///|
/// Tests for Colless-like index (Mir, Rotger & Rossello 2018).

///|
fn cl_make_leaf(name : String) -> @src.Clade {
  @src.Clade::new(name=Some(name))
}

///|
fn cl_make_internal(children : Array[@src.Clade]) -> @src.Clade {
  @src.Clade::new(name=Some("internal"), clades=children)
}

///|
/// Balanced binary tree with n leaves.
fn cl_balanced(n : Int) -> @src.Clade {
  if n <= 1 {
    return cl_make_leaf("L")
  }
  let half = n / 2
  cl_make_internal([cl_balanced(half), cl_balanced(n - half)])
}

///|
/// Caterpillar (maximally imbalanced) tree with n leaves.
fn cl_caterpillar(n : Int) -> @src.Clade {
  if n <= 1 {
    return cl_make_leaf("L")
  }
  if n == 2 {
    return cl_make_internal([cl_make_leaf("L"), cl_make_leaf("L")])
  }
  cl_make_internal([cl_make_leaf("L"), cl_caterpillar(n - 1)])
}

///|
test "colless_like: single leaf" {
  let leaf = cl_make_leaf("A")
  let cl = @src.dendropy_colless_like(leaf, @src.CollessLikeFSize::FExp, @src.CollessLikeDissim::DMdm)
  assert_true((cl - 0.0).abs() < 0.0000000001)
}

///|
test "colless_like: balanced tree is zero" {
  // Perfectly balanced tree: all balance values should be 0
  let tree = cl_balanced(8)
  let cl = @src.dendropy_colless_like(tree, @src.CollessLikeFSize::FExp, @src.CollessLikeDissim::DMdm)
  assert_true((cl - 0.0).abs() < 0.0000000001)
}

///|
test "colless_like: caterpillar positive" {
  let tree = cl_caterpillar(8)
  let cl = @src.dendropy_colless_like(tree, @src.CollessLikeFSize::FExp, @src.CollessLikeDissim::DMdm)
  assert_true(cl > 0.0)
}

///|
test "colless_like: caterpillar vs balanced" {
  let balanced = cl_balanced(8)
  let caterpillar = cl_caterpillar(8)
  let cl_bal = @src.dendropy_colless_like(balanced, @src.CollessLikeFSize::FExp, @src.CollessLikeDissim::DMdm)
  let cl_cat = @src.dendropy_colless_like(caterpillar, @src.CollessLikeFSize::FExp, @src.CollessLikeDissim::DMdm)
  assert_true(cl_cat > cl_bal)
}

///|
test "colless_like: exp vs ln weight" {
  // Different weight functions should give different values for imbalanced tree
  let tree = cl_caterpillar(6)
  let cl_exp = @src.dendropy_colless_like(tree, @src.CollessLikeFSize::FExp, @src.CollessLikeDissim::DMdm)
  let cl_ln = @src.dendropy_colless_like(tree, @src.CollessLikeFSize::FLn, @src.CollessLikeDissim::DMdm)
  assert_true(cl_exp > 0.0)
  assert_true(cl_ln > 0.0)
  // They should differ (exp grows faster)
  assert_true((cl_exp - cl_ln).abs() > 0.001)
}

///|
test "colless_like: mdm vs var vs sd" {
  let tree = cl_caterpillar(6)
  let cl_mdm = @src.dendropy_colless_like(tree, @src.CollessLikeFSize::FExp, @src.CollessLikeDissim::DMdm)
  let cl_var = @src.dendropy_colless_like(tree, @src.CollessLikeFSize::FExp, @src.CollessLikeDissim::DVar)
  let cl_sd = @src.dendropy_colless_like(tree, @src.CollessLikeFSize::FExp, @src.CollessLikeDissim::DSd)
  assert_true(cl_mdm > 0.0)
  assert_true(cl_var > 0.0)
  assert_true(cl_sd > 0.0)
  // sd = sqrt(var), but they're summed over nodes so not directly related
  // Just check all are positive for imbalanced tree
}

///|
test "colless_like: multifurcating tree" {
  // Star tree: root with 4 leaves (maximally multifurcating)
  let star = cl_make_internal([
    cl_make_leaf("A"), cl_make_leaf("B"), cl_make_leaf("C"), cl_make_leaf("D"),
  ])
  // All children have same f.size, so balance = 0
  let cl = @src.dendropy_colless_like(star, @src.CollessLikeFSize::FExp, @src.CollessLikeDissim::DMdm)
  assert_true((cl - 0.0).abs() < 0.0000000001)
}

///|
test "colless_like: asymmetric multifurcating" {
  // Root with 3 children: one leaf, one cherry, one leaf
  //   (A, (B, C), D)
  let tree = cl_make_internal([
    cl_make_leaf("A"),
    cl_make_internal([cl_make_leaf("B"), cl_make_leaf("C")]),
    cl_make_leaf("D"),
  ])
  // The cherry node has balance 0 (two identical children)
  // The root has 3 children with different f.sizes: leaf vs cherry vs leaf
  // Leaf f.size = f(0) = e^0 = 1
  // Cherry f.size = f(2) + 2*f(0) = e^2 + 2
  // So root balance = mdm(1, e^2+2, 1) > 0
  let cl = @src.dendropy_colless_like(tree, @src.CollessLikeFSize::FExp, @src.CollessLikeDissim::DMdm)
  assert_true(cl > 0.0)
}

///|
test "colless_like: two-leaf cherry zero" {
  // (A, B) — single cherry, both children have same f.size
  let tree = cl_make_internal([cl_make_leaf("A"), cl_make_leaf("B")])
  let cl = @src.dendropy_colless_like(tree, @src.CollessLikeFSize::FExp, @src.CollessLikeDissim::DMdm)
  assert_true((cl - 0.0).abs() < 0.0000000001)
}

///|
test "colless_like: ln weight balanced zero" {
  let tree = cl_balanced(8)
  let cl = @src.dendropy_colless_like(tree, @src.CollessLikeFSize::FLn, @src.CollessLikeDissim::DVar)
  assert_true((cl - 0.0).abs() < 0.0000000001)
}

///|
test "colless_like: 4-leaf caterpillar manual" {
  // (A, (B, (C, D))) with f=exp, D=mdm
  // Leaf f.size = e^0 = 1
  // Cherry (C,D): f.size = e^2 + 2*1 = e^2 + 2
  //   balance = mdm(1, 1) = 0
  // Node (B, (C,D)): f.size = e^2 + 1 + (e^2+2) = 2*e^2 + 3
  //   children f.sizes: [1, e^2+2]
  //   balance = mdm(1, e^2+2) = |1 - median| + |e^2+2 - median| / 2
  //   median = (1 + e^2+2)/2 = (e^2+3)/2
  //   = (|(e^2+1)/2| + |(e^2+1)/2|) / 2 = (e^2+1)/2
  // Root (A, (B,(C,D))): children f.sizes: [1, 2*e^2+3]
  //   balance = mdm(1, 2*e^2+3) = (2*e^2+2)/2 = e^2+1
  // Total = 0 + (e^2+1)/2 + (e^2+1)
  let tree = cl_caterpillar(4)
  let cl = @src.dendropy_colless_like(tree, @src.CollessLikeFSize::FExp, @src.CollessLikeDissim::DMdm)
  let e2 = @math.exp(2.0)
  let expected = 0.0 + (e2 + 1.0) / 2.0 + (e2 + 1.0)
  assert_true((cl - expected).abs() < 0.0001)
}
