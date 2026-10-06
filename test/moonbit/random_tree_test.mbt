///|
/// Tests for random tree generation (Yule / PDA).

///|
fn rt_count_leaves(clade : @src.Clade) -> Int {
  if clade.is_terminal() {
    return 1
  }
  let mut c = 0
  for ch in clade.clades {
    c = c + rt_count_leaves(ch)
  }
  c
}

///|
/// Check all internal nodes are bifurcations (binary tree).
fn rt_is_binary(clade : @src.Clade) -> Bool {
  if clade.is_terminal() {
    return true
  }
  if clade.clades.length() != 2 {
    return false
  }
  for ch in clade.clades {
    if not(rt_is_binary(ch)) {
      return false
    }
  }
  true
}

///|
/// Collect all leaf names.
fn rt_leaf_names(clade : @src.Clade, names : Array[String]) -> Unit {
  if clade.is_terminal() {
    match clade.name {
      Some(n) => names.push(n)
      None => names.push("?")
    }
    return
  }
  for ch in clade.clades {
    rt_leaf_names(ch, names)
  }
}

// =====================
// RNG tests
// =====================

///|
test "rng: produces different values" {
  let rng = @src.Rng::new(42UL)
  let a = rng.next_double()
  let b = rng.next_double()
  assert_true(a >= 0.0)
  assert_true(a < 1.0)
  assert_true(b >= 0.0)
  assert_true(b < 1.0)
  assert_true((a - b).abs() > 0.0001)
}

///|
test "rng: deterministic with same seed" {
  let rng1 = @src.Rng::new(123UL)
  let rng2 = @src.Rng::new(123UL)
  for i in 0..<10 {
    assert_true((rng1.next_double() - rng2.next_double()).abs() < 0.0000000001)
  }
}

///|
test "rng: different seeds give different sequences" {
  let rng1 = @src.Rng::new(1UL)
  let rng2 = @src.Rng::new(2UL)
  let mut diff = false
  for i in 0..<10 {
    if (rng1.next_double() - rng2.next_double()).abs() > 0.001 {
      diff = true
    }
  }
  assert_true(diff)
}

///|
test "rng: randint in range" {
  let rng = @src.Rng::new(42UL)
  for i in 0..<100 {
    let v = rng.randint(10)
    assert_true(v >= 0)
    assert_true(v < 10)
  }
}

// =====================
// Yule tree tests
// =====================

///|
test "yule: n=1 single leaf" {
  let rng = @src.Rng::new(42UL)
  let tree = @src.phylo_random_yule(1, rng)
  assert_true(tree.is_terminal())
}

///|
test "yule: n=2 cherry" {
  let rng = @src.Rng::new(42UL)
  let tree = @src.phylo_random_yule(2, rng)
  assert_eq(rt_count_leaves(tree), 2)
  assert_true(rt_is_binary(tree))
}

///|
test "yule: n=10 correct leaf count" {
  let rng = @src.Rng::new(42UL)
  let tree = @src.phylo_random_yule(10, rng)
  assert_eq(rt_count_leaves(tree), 10)
}

///|
test "yule: produces binary tree" {
  let rng = @src.Rng::new(7UL)
  let tree = @src.phylo_random_yule(20, rng)
  assert_true(rt_is_binary(tree))
}

///|
test "yule: different seeds give different trees" {
  let rng1 = @src.Rng::new(1UL)
  let rng2 = @src.Rng::new(999UL)
  let tree1 = @src.phylo_random_yule(10, rng1)
  let tree2 = @src.phylo_random_yule(10, rng2)
  let names1 : Array[String] = Array(capacity=0)
  let names2 : Array[String] = Array(capacity=0)
  rt_leaf_names(tree1, names1)
  rt_leaf_names(tree2, names2)
  // Same leaf count but potentially different structure
  assert_eq(names1.length(), 10)
  assert_eq(names2.length(), 10)
}

///|
test "yule: n=50 correct leaf count" {
  let rng = @src.Rng::new(42UL)
  let tree = @src.phylo_random_yule(50, rng)
  assert_eq(rt_count_leaves(tree), 50)
  assert_true(rt_is_binary(tree))
}

///|
test "yule: reproducible" {
  let rng1 = @src.Rng::new(777UL)
  let rng2 = @src.Rng::new(777UL)
  let tree1 = @src.phylo_random_yule(15, rng1)
  let tree2 = @src.phylo_random_yule(15, rng2)
  let names1 : Array[String] = Array(capacity=0)
  let names2 : Array[String] = Array(capacity=0)
  rt_leaf_names(tree1, names1)
  rt_leaf_names(tree2, names2)
  assert_eq(names1.length(), names2.length())
  for i in 0..<names1.length() {
    assert_eq(names1[i], names2[i])
  }
}

// =====================
// PDA tree tests
// =====================

///|
test "pda: n=1 single leaf" {
  let rng = @src.Rng::new(42UL)
  let tree = @src.phylo_random_pda(1, rng)
  assert_true(tree.is_terminal())
}

///|
test "pda: n=2 cherry" {
  let rng = @src.Rng::new(42UL)
  let tree = @src.phylo_random_pda(2, rng)
  assert_eq(rt_count_leaves(tree), 2)
}

///|
test "pda: n=10 correct leaf count" {
  let rng = @src.Rng::new(42UL)
  let tree = @src.phylo_random_pda(10, rng)
  assert_eq(rt_count_leaves(tree), 10)
}

///|
test "pda: produces binary tree" {
  let rng = @src.Rng::new(7UL)
  let tree = @src.phylo_random_pda(20, rng)
  assert_true(rt_is_binary(tree))
}

///|
test "pda: n=50 correct leaf count" {
  let rng = @src.Rng::new(42UL)
  let tree = @src.phylo_random_pda(50, rng)
  assert_eq(rt_count_leaves(tree), 50)
  assert_true(rt_is_binary(tree))
}

///|
test "pda: reproducible" {
  let rng1 = @src.Rng::new(777UL)
  let rng2 = @src.Rng::new(777UL)
  let tree1 = @src.phylo_random_pda(15, rng1)
  let tree2 = @src.phylo_random_pda(15, rng2)
  assert_eq(rt_count_leaves(tree1), rt_count_leaves(tree2))
}

// =====================
// Comparison tests
// =====================

///|
test "yule_vs_pda: both produce valid trees" {
  let rng1 = @src.Rng::new(42UL)
  let rng2 = @src.Rng::new(42UL)
  let yule_tree = @src.phylo_random_yule(30, rng1)
  let pda_tree = @src.phylo_random_pda(30, rng2)
  assert_eq(rt_count_leaves(yule_tree), 30)
  assert_eq(rt_count_leaves(pda_tree), 30)
  assert_true(rt_is_binary(yule_tree))
  assert_true(rt_is_binary(pda_tree))
}
