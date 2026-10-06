///|
/// Tests for scale tree.

///|
fn st_leaf(name : String, bl : Double) -> @src.Clade {
  @src.Clade::new(name=Some(name), branch_length=Some(bl))
}

///|
fn st_internal(children : Array[@src.Clade], bl : Double) -> @src.Clade {
  @src.Clade::new(name=Some("internal"), branch_length=Some(bl), clades=children)
}

///|
/// Collect all branch lengths from a tree.
fn st_collect_bls(clade : @src.Clade, bls : Array[Double]) -> Unit {
  match clade.branch_length {
    Some(v) => bls.push(v)
    None => ()
  }
  for ch in clade.clades {
    st_collect_bls(ch, bls)
  }
}

///|
test "scale: leaf" {
  let leaf = st_leaf("A", 2.0)
  let result = @src.phylo_scale_tree(leaf, 3.0)
  match result.branch_length {
    Some(v) => assert_true((v - 6.0).abs() < 0.0000000001)
    None => assert_true(false)
  }
}

///|
test "scale: basic tree" {
  // ((A:1,B:2):3,C:4) * 2 = ((A:2,B:4):6,C:8)
  let tree = st_internal([
    st_internal([st_leaf("A", 1.0), st_leaf("B", 2.0)], 3.0),
    st_leaf("C", 4.0),
  ], 0.0)
  let result = @src.phylo_scale_tree(tree, 2.0)
  let bls : Array[Double] = Array(capacity=0)
  st_collect_bls(result, bls)
  // Should have [6, 2, 4, 8] (root bl, A bl, B bl, C bl)
  assert_eq(bls.length(), 5)
  let mut has2 = false
  let mut has4 = false
  let mut has6 = false
  let mut has8 = false
  for bl in bls {
    if (bl - 2.0).abs() < 0.0001 { has2 = true }
    if (bl - 4.0).abs() < 0.0001 { has4 = true }
    if (bl - 6.0).abs() < 0.0001 { has6 = true }
    if (bl - 8.0).abs() < 0.0001 { has8 = true }
  }
  assert_true(has2)
  assert_true(has4)
  assert_true(has6)
  assert_true(has8)
}

///|
test "scale: factor 1 unchanged" {
  let tree = st_internal([st_leaf("A", 1.5), st_leaf("B", 2.5)], 3.5)
  let result = @src.phylo_scale_tree(tree, 1.0)
  let bls_orig : Array[Double] = Array(capacity=0)
  let bls_new : Array[Double] = Array(capacity=0)
  st_collect_bls(tree, bls_orig)
  st_collect_bls(result, bls_new)
  assert_eq(bls_orig.length(), bls_new.length())
  for i in 0..<bls_orig.length() {
    assert_true((bls_orig[i] - bls_new[i]).abs() < 0.0000000001)
  }
}

///|
test "scale: factor 0 zeros all" {
  let tree = st_internal([st_leaf("A", 1.0), st_leaf("B", 2.0)], 3.0)
  let result = @src.phylo_scale_tree(tree, 0.0)
  let bls : Array[Double] = Array(capacity=0)
  st_collect_bls(result, bls)
  for bl in bls {
    assert_true((bl - 0.0).abs() < 0.0000000001)
  }
}

///|
test "scale: fractional factor" {
  let tree = st_internal([st_leaf("A", 4.0), st_leaf("B", 6.0)], 2.0)
  let result = @src.phylo_scale_tree(tree, 0.5)
  let bls : Array[Double] = Array(capacity=0)
  st_collect_bls(result, bls)
  // Should have [1, 2, 3]
  let mut has1 = false
  let mut has2 = false
  let mut has3 = false
  for bl in bls {
    if (bl - 1.0).abs() < 0.0001 { has1 = true }
    if (bl - 2.0).abs() < 0.0001 { has2 = true }
    if (bl - 3.0).abs() < 0.0001 { has3 = true }
  }
  assert_true(has1)
  assert_true(has2)
  assert_true(has3)
}

///|
test "scale: preserves structure" {
  let tree = st_internal([
    st_internal([st_leaf("A", 1.0), st_leaf("B", 1.0)], 1.0),
    st_internal([st_leaf("C", 1.0), st_leaf("D", 1.0)], 1.0),
  ], 1.0)
  let result = @src.phylo_scale_tree(tree, 10.0)
  // Should still have 4 leaves and 3 internal nodes
  let lc = @src.phylo_count_leaves(result)
  assert_eq(lc, 4)
}

///|
test "scale: negative factor" {
  let tree = st_internal([st_leaf("A", 1.0), st_leaf("B", 2.0)], 3.0)
  let result = @src.phylo_scale_tree(tree, -1.0)
  let bls : Array[Double] = Array(capacity=0)
  st_collect_bls(result, bls)
  // All should be negative
  for bl in bls {
    assert_true(bl < 0.0)
  }
}
