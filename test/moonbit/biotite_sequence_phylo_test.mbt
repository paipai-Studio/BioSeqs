///|
/// Tests for biotite_sequence_phylo.mbt.
test "tree from newick and back" {
  let t = @src.biotite_tree_from_newick("((a:1.0,b:2.0):3.0,c:4.0);")
  assert_eq(t.n_leaves(), 3)
  let newick = t.as_newick()
  assert_true(newick.has_prefix("(("))
  assert_true(newick.has_suffix(";"))
  assert_true(newick.contains("a:1"))
  assert_true(newick.contains("c:4"))
}

///|
test "tree get_distance" {
  let t = @src.biotite_tree_from_newick("((a:1,b:2):3,c:4);")
  match t.get_distance(0, 1) {
    Some(d) => assert_true(@moonbit.assert_num_close(d, 3.0, 1.0e-12))
    None => fail("distance a-b should exist")
  }
  match t.get_distance(0, 2) {
    Some(d) => assert_true(@moonbit.assert_num_close(d, 8.0, 1.0e-12))
    None => fail("distance a-c should exist")
  }
  // topological distances count branches as 1
  match t.get_distance(0, 1, topological=true) {
    Some(d) => assert_true(@moonbit.assert_num_close(d, 2.0, 1.0e-12))
    None => fail("topological a-b should exist")
  }
  match t.get_distance(0, 2, topological=true) {
    Some(d) => assert_true(@moonbit.assert_num_close(d, 3.0, 1.0e-12))
    None => fail("topological a-c should exist")
  }
  // missing leaf index
  assert_eq(t.get_distance(0, 99), None)
}

///|
test "tree leaves sorted by index" {
  let t = @src.biotite_tree_from_newick("(c:4,(a:1,b:2):3);")
  let leaves = t.leaves()
  assert_eq(leaves.length(), 3)
  assert_eq(leaves[0].index, Some(0))
  assert_eq(leaves[1].index, Some(1))
  assert_eq(leaves[2].index, Some(2))
  match leaves[0].name {
    Some(n) => assert_eq(n, "c")
    None => fail("leaf should keep its name")
  }
}

///|
test "node constructors" {
  let leaf = @src.biotite_tree_node_new(index=Some(0), name=Some("leaf"))
  assert_true(leaf.is_leaf())
  let inner = @src.biotite_tree_node_new(name=Some("root"), children=[leaf], distances=[
    0.5,
  ])
  assert_false(inner.is_leaf())
  assert_eq(inner.children.length(), 1)
  match leaf.parent {
    Some(_) => ()
    None => fail("parent pointer should be set")
  }
  assert_true(@moonbit.assert_num_close(leaf.distance, 0.5, 1.0e-12))
}

///|
test "upgma symmetric tree" {
  let t = @src.biotite_upgma([[0.0, 1.0, 2.0], [1.0, 0.0, 2.0], [2.0, 2.0, 0.0]])
  assert_eq(t.n_leaves(), 3)
  // symmetric: leaves 0 and 1 are closest, so they share a parent
  match t.get_distance(0, 1) {
    Some(d) => assert_true(@moonbit.assert_num_close(d, 1.0, 1.0e-9))
    None => fail("upgma distance 0-1 should exist")
  }
  match t.get_distance(0, 2) {
    Some(d) => assert_true(d > 1.0)
    None => fail("upgma distance 0-2 should exist")
  }
  // ultrametric: all leaves equidistant from root
  match (t.get_distance(0, 2), t.get_distance(1, 2)) {
    (Some(a), Some(b)) => assert_true(@moonbit.assert_num_close(a, b, 1.0e-9))
    _ => fail("ultrametric distances should exist")
  }
}

///|
test "upgma hand-computed merge distances" {
  let t = @src.biotite_upgma([[0.0, 1.0, 2.0], [1.0, 0.0, 3.0], [2.0, 3.0, 0.0]])
  // merge {0,1} at age 0.5; root at age (2+3)/2/2 = 1.25
  match t.get_distance(0, 1) {
    Some(d) => assert_true(@moonbit.assert_num_close(d, 1.0, 1.0e-9))
    None => fail("merge distance should exist")
  }
  match t.get_distance(0, 2) {
    Some(d) => assert_true(@moonbit.assert_num_close(d, 2.5, 1.0e-9))
    None => fail("root distance should exist")
  }
}

///|
test "neighbor joining three taxa" {
  let t = @src.biotite_neighbor_joining([
    [0.0, 1.0, 2.0],
    [1.0, 0.0, 3.0],
    [2.0, 3.0, 0.0],
  ])
  assert_eq(t.n_leaves(), 3)
  match t.get_distance(0, 1) {
    Some(d) => assert_true(@moonbit.assert_num_close(d, 1.0, 1.0e-9))
    None => fail("nj distance 0-1 should exist")
  }
  match t.get_distance(0, 2) {
    Some(d) => assert_true(@moonbit.assert_num_close(d, 2.0, 1.0e-9))
    None => fail("nj distance 0-2 should exist")
  }
  match t.get_distance(1, 2) {
    Some(d) => assert_true(@moonbit.assert_num_close(d, 3.0, 1.0e-9))
    None => fail("nj distance 1-2 should exist")
  }
}

///|
test "neighbor joining star tree" {
  let t = @src.biotite_neighbor_joining([
    [0.0, 2.0, 2.0, 2.0],
    [2.0, 0.0, 2.0, 2.0],
    [2.0, 2.0, 0.0, 2.0],
    [2.0, 2.0, 2.0, 0.0],
  ])
  assert_eq(t.n_leaves(), 4)
  // all pairwise distances preserved up to tolerance
  match t.get_distance(0, 1) {
    Some(d) => assert_true(@moonbit.assert_num_close(d, 2.0, 1.0e-6))
    None => fail("star distance should exist")
  }
}

///|
test "edge cases empty and single" {
  let empty = @src.biotite_upgma([])
  assert_eq(empty.n_leaves(), 0)
  let single = @src.biotite_neighbor_joining([[0.0]])
  assert_eq(single.n_leaves(), 1)
  let nj_empty = @src.biotite_neighbor_joining([])
  assert_eq(nj_empty.n_leaves(), 0)
}
