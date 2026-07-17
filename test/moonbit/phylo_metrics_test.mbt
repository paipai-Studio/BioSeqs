///|
test "phylo_metrics: total_branch_length" {
  let tree = @src.create_demo_tree()
  let total = tree.total_branch_length()
  assert_true(total > 4.9 && total < 5.1)
}

///|
test "phylo_metrics: max_depth" {
  let tree = @src.create_demo_tree()
  let depth = tree.max_depth()
  assert_true(depth > 1.4 && depth < 1.6)
}

///|
test "phylo_metrics: get_leaf_names" {
  let tree = @src.create_demo_tree()
  let root = tree.get_root()
  let leaves = root.get_leaf_names()
  assert_eq(leaves.length(), 4)
  assert_true(leaves.contains("A"))
  assert_true(leaves.contains("B"))
  assert_true(leaves.contains("C"))
  assert_true(leaves.contains("D"))
}

///|
test "phylo_metrics: colless_index" {
  let tree = @src.create_demo_tree()
  let index = tree.colless_index()
  assert_eq(index, 0)
}

///|
test "phylo_metrics: patristic_distance" {
  let tree = @src.create_demo_tree()
  let dist = tree.patristic_distance("A", "B")
  assert_true(dist is Some(_))
  assert_true(dist.unwrap() > 1.9 && dist.unwrap() < 2.1)

  let dist2 = tree.patristic_distance("A", "C")
  assert_true(dist2 is Some(_))
  assert_true(dist2.unwrap() > 2.9 && dist2.unwrap() < 3.1)
}

///|
test "phylo_metrics: patristic_distance_matrix" {
  let tree = @src.create_demo_tree()
  let matrix = tree.patristic_distance_matrix()
  assert_eq(matrix.length(), 6)
}

///|
test "phylo_metrics: count_internal" {
  let tree = @src.create_demo_tree()
  let count = tree.count_internal()
  assert_eq(count, 3)
}

///|
test "phylo_metrics: get_bipartitions" {
  let tree = @src.create_demo_tree()
  let root = tree.get_root()
  let biparts = root.get_bipartitions()
  assert_true(biparts.length() > 0)
}

///|
test "phylo_metrics: robinson_foulds_distance_same_tree" {
  let tree1 = @src.create_demo_tree()
  let tree2 = @src.create_demo_tree()
  let rf = @src.robinson_foulds_distance(tree1, tree2)
  assert_eq(rf, 0)
}

///|
test "phylo_metrics: robinson_foulds_distance_different_tree" {
  let tree1 = @src.create_demo_tree()
  let tree2 = @src.create_demo_tree2()
  let rf = @src.robinson_foulds_distance(tree1, tree2)
  assert_true(rf > 0)
}
