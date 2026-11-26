///|
/// Tests for biotite_sequence_phylo.mbt — verified against the official
/// Biotite doctest scenarios (integer node numbering, root distance).
test "official from_newick node numbering" {
  // Official doctest: from_newick("((0:5.0,1:7.0):3.0,2:10.0);")
  // nodes [0,1,2,3,4] — root is node 4.
  let t = @src.biotite_tree_from_newick("((0:5.0,1:7.0):3.0,2:10.0);")
  assert_eq(t.n_leaves(), 3)
  match t.find_node(0) {
    Some(n) => assert_eq(n.node_id, 0)
    None => fail("leaf 0 must exist")
  }
  match t.find_node(4) {
    Some(n) =>
      // root: highest id, no parent
      match n.parent {
        Some(_) => fail("node 4 should be the root")
        None => ()
      }
    None => fail("node 4 must exist")
  }
  match t.find_node(3) {
    Some(n) =>
      // intermediate node closing first: parent of leaves 0 and 1
      assert_eq(n.children.length(), 2)
    None => fail("node 3 must exist")
  }
}

///|
test "official get_distance doctest" {
  let t = @src.biotite_tree_from_newick("((0:5.0,1:7.0):3.0,2:10.0);")
  match t.get_distance(0, 1) {
    Some(d) => assert_true(@moonbit.assert_num_close(d, 12.0, 1.0e-12))
    None => fail("distance 0-1 should exist")
  }
  match t.get_distance(0, 2) {
    Some(d) => assert_true(@moonbit.assert_num_close(d, 18.0, 1.0e-12))
    None => fail("distance 0-2 should exist")
  }
  match t.get_distance(0, 2, topological=true) {
    Some(d) => assert_true(@moonbit.assert_num_close(d, 3.0, 1.0e-12))
    None => fail("topological distance should exist")
  }
}

///|
test "official leaf distances matrix" {
  let t = @src.biotite_tree_from_newick("((0:5.0,1:7.0):3.0,2:10.0);")
  let m = t.get_leaf_distances()
  assert_eq(m.length(), 3)
  assert_true(@moonbit.assert_num_close(m[0][1], 12.0, 1.0e-12))
  assert_true(@moonbit.assert_num_close(m[0][2], 18.0, 1.0e-12))
  assert_true(@moonbit.assert_num_close(m[1][2], 20.0, 1.0e-12))
  assert_true(@moonbit.assert_num_close(m[2][0], 18.0, 1.0e-12))
}

///|
test "official to_newick root distance" {
  let t = @src.biotite_tree_from_newick("((0:5.0,1:7.0):3.0,2:10.0);")
  let newick = t.as_newick(round_distance=Some(1))
  // Official doctest: ((0:5.0,1:7.0):3.0,2:10.0):0.0;
  assert_eq(newick, "((0:5.0,1:7.0):3.0,2:10.0):0.0;")
  let no_dist = t.as_newick(include_distance=false)
  assert_eq(no_dist, "((0,1),2);")
}

///|
test "newick intermediate node numbering by closure order" {
  // ((0,1),(2,3)): leaves 0..3, first closed inner = 4, second = 5,
  // root = 5? The outer bracket closes last -> root = 5? No: two inner
  // nodes close first (4 and 5), the root closes last (6).
  let t = @src.biotite_tree_from_newick("((0,1),(2,3));")
  assert_eq(t.n_leaves(), 4)
  match t.find_node(4) {
    Some(n) => assert_eq(n.children.length(), 2)
    None => fail("node 4 must exist")
  }
  // root must have the highest id among all nodes
  let root = t.root
  let max_id = root.node_id
  let stack : Array[@src.BiotiteTreeNode] = [root]
  while stack.length() > 0 {
    let n = stack.pop().unwrap()
    assert_true(n.node_id <= max_id)
    for c in n.children {
      stack.push(c)
    }
  }
}

///|
test "newick roundtrip with names" {
  let t = @src.biotite_tree_from_newick("((a:1,b:2):3,c:4);")
  let newick = t.as_newick()
  assert_true(newick.has_prefix("((a:1"))
  assert_true(newick.has_suffix(";"))
  assert_true(newick.has_suffix("):0.0;") || newick.contains("c:4"))
  let t2 = @src.biotite_tree_from_newick(newick)
  assert_eq(t2.n_leaves(), 3)
  match (t.get_distance(0, 1), t2.get_distance(0, 1)) {
    (Some(d1), Some(d2)) =>
      assert_true(@moonbit.assert_num_close(d1, d2, 1.0e-12))
    _ => fail("distances should exist")
  }
}

///|
test "newick distance default zero" {
  let t = @src.biotite_tree_from_newick("((0,1),2);")
  // missing distances parse as 0.0
  match t.get_distance(0, 1) {
    Some(d) => assert_true(@moonbit.assert_num_close(d, 0.0, 1.0e-12))
    None => fail("distance should exist")
  }
}

///|
test "node constructors and parent pointers" {
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
  match t.get_distance(0, 1) {
    Some(d) => assert_true(@moonbit.assert_num_close(d, 1.0, 1.0e-9))
    None => fail("upgma distance 0-1 should exist")
  }
  match t.get_distance(0, 2) {
    Some(d) => assert_true(d > 1.0)
    None => fail("upgma distance 0-2 should exist")
  }
  match (t.get_distance(0, 2), t.get_distance(1, 2)) {
    (Some(a), Some(b)) => assert_true(@moonbit.assert_num_close(a, b, 1.0e-9))
    _ => fail("ultrametric distances should exist")
  }
}

///|
test "upgma hand-computed merge distances" {
  let t = @src.biotite_upgma([[0.0, 1.0, 2.0], [1.0, 0.0, 3.0], [2.0, 3.0, 0.0]])
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
  match t.get_distance(0, 1) {
    Some(d) => assert_true(@moonbit.assert_num_close(d, 2.0, 1.0e-6))
    None => fail("star distance should exist")
  }
}

///|
test "as binary splits multi-way nodes" {
  // Official doctest: (0:1.0,1:2.0,2:3.0) -> ((0:1.0,1:2.0):0.0,2:3.0)
  let t = @src.biotite_tree_from_newick("(0:1.0,1:2.0,2:3.0);")
  let binary = @src.biotite_tree_as_binary(t)
  let newick = binary.as_newick(round_distance=Some(1))
  assert_eq(newick, "((0:1.0,1:2.0):0.0,2:3.0):0.0;")
}

///|
test "as binary omits single-child nodes" {
  // Official doctest: ((0:1.0,1:2.0):5.0):3.0 -> (0:1.0,1:2.0)
  let t = @src.biotite_tree_from_newick("((0:1.0,1:2.0):5.0):3.0;")
  let binary = @src.biotite_tree_as_binary(t)
  let newick = binary.as_newick(round_distance=Some(1))
  assert_eq(newick, "(0:1.0,1:2.0):0.0;")
}

///|
test "as binary preserves leaf distances" {
  let t = @src.biotite_tree_from_newick("(0:1.0,1:2.0,2:3.0);")
  let binary = @src.biotite_tree_as_binary(t)
  match binary.get_distance(0, 2) {
    Some(d) => assert_true(@moonbit.assert_num_close(d, 4.0, 1.0e-9))
    None => fail("distance 0-2 should exist")
  }
  match binary.get_distance(1, 2) {
    Some(d) => assert_true(@moonbit.assert_num_close(d, 5.0, 1.0e-9))
    None => fail("distance 1-2 should exist")
  }
}
