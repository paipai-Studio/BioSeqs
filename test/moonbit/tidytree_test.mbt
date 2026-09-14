///|
/// Tests for tidytree module.

///|
test "tidytree as_tibble" {
  let tree = @src.treeio_example_tree()
  let tibble = @src.tidy_as_tibble(tree)
  assert_eq(tibble.n_nodes, tree.n_nodes)
  assert_eq(tibble.n_tips, tree.n_tips)
}

///|
test "tidytree offspring" {
  let tree = @src.treeio_example_tree()
  // Root's offspring should be all other nodes
  let offspring = @src.tidy_offspring(tree, tree.root)
  assert_true(offspring.length() >= 2)
}

///|
test "tidytree child" {
  let tree = @src.treeio_example_tree()
  let children = @src.tidy_child(tree, tree.root)
  assert_true(children.length() >= 2)
}

///|
test "tidytree parent" {
  let tree = @src.treeio_example_tree()
  // Root's parent should be -1
  assert_eq(@src.tidy_parent(tree, tree.root), -1)
}

///|
test "tidytree sibling" {
  let tree = @src.treeio_example_tree()
  // Find a tip and check it has siblings
  let tips = @src.tidy_tips(tree)
  if tips.length() >= 2 {
    let siblings = @src.tidy_sibling(tree, tips[0])
    assert_true(siblings.length() >= 1)
  }
}

///|
test "tidytree mrca" {
  let tree = @src.treeio_example_tree()
  let tips = @src.tidy_tips(tree)
  if tips.length() >= 2 {
    let mrca = @src.tidy_mrca(tree, tips[0], tips[1])
    assert_true(mrca >= 0)
  }
}

///|
test "tidytree is leaf" {
  let tree = @src.treeio_example_tree()
  let tips = @src.tidy_tips(tree)
  if tips.length() > 0 {
    assert_true(@src.tidy_is_leaf(tree, tips[0]))
  }
  // Root should not be a leaf
  assert_true(!@src.tidy_is_leaf(tree, tree.root))
}

///|
test "tidytree node depth" {
  let tree = @src.treeio_example_tree()
  // Root depth should be 0
  assert_eq(@src.tidy_node_depth(tree, tree.root), 0)
  // Tips should have depth > 0
  let tips = @src.tidy_tips(tree)
  if tips.length() > 0 {
    assert_true(@src.tidy_node_depth(tree, tips[0]) > 0)
  }
}

///|
test "tidytree ancestors" {
  let tree = @src.treeio_example_tree()
  let tips = @src.tidy_tips(tree)
  if tips.length() > 0 {
    let ancestors = @src.tidy_ancestors(tree, tips[0])
    assert_true(ancestors.length() >= 1)
    // Last ancestor should be root
    assert_eq(ancestors[ancestors.length() - 1], tree.root)
  }
}

///|
test "tidytree tips" {
  let tree = @src.treeio_example_tree()
  let tips = @src.tidy_tips(tree)
  assert_eq(tips.length(), tree.n_tips)
}

///|
test "tidytree internal" {
  let tree = @src.treeio_example_tree()
  let internal = @src.tidy_internal(tree)
  assert_eq(internal.length(), tree.n_nodes - tree.n_tips)
}

///|
test "tidytree filter tips" {
  let tree = @src.treeio_example_tree()
  let tibble = @src.tidy_as_tibble(tree)
  let tips_only = @src.tidy_filter_tips(tibble)
  assert_eq(tips_only.n_tips, tibble.n_tips)
}

///|
test "tidytree filter internal" {
  let tree = @src.treeio_example_tree()
  let tibble = @src.tidy_as_tibble(tree)
  let internal_only = @src.tidy_filter_internal(tibble)
  assert_eq(internal_only.n_tips, 0)
}

///|
test "tidytree get node by label" {
  let tree = @src.treeio_example_tree()
  let node_id = @src.tidy_get_node_by_label(tree, "A")
  assert_true(node_id >= 0)
}

///|
test "tidytree get label" {
  let tree = @src.treeio_example_tree()
  let node_id = @src.tidy_get_node_by_label(tree, "A")
  assert_eq(@src.tidy_get_label(tree, node_id), "A")
}

///|
test "tidytree branch length" {
  let tree = @src.treeio_example_tree()
  let node_id = @src.tidy_get_node_by_label(tree, "A")
  let blen = @src.tidy_branch_length(tree, node_id)
  assert_true((blen - 0.1).abs() < 0.001)
}

///|
test "tidytree dist to root" {
  let tree = @src.treeio_example_tree()
  // Root distance to root is 0
  assert_true(@src.tidy_dist_to_root(tree, tree.root) < 0.001)
}

///|
test "tidytree patristic distance" {
  let tree = @src.treeio_example_tree()
  let tips = @src.tidy_tips(tree)
  if tips.length() >= 2 {
    let dist = @src.tidy_patristic_distance(tree, tips[0], tips[1])
    assert_true(dist > 0.0)
  }
}

///|
test "tidytree post order" {
  let tree = @src.treeio_example_tree()
  let order = @src.tidy_post_order(tree)
  assert_eq(order.length(), tree.n_nodes)
  // Last node should be root
  assert_eq(order[order.length() - 1], tree.root)
}

///|
test "tidytree pre order" {
  let tree = @src.treeio_example_tree()
  let order = @src.tidy_pre_order(tree)
  assert_eq(order.length(), tree.n_nodes)
  // First node should be root
  assert_eq(order[0], tree.root)
}

///|
test "tidytree print" {
  let tibble = @src.tidy_example_tibble()
  let output = @src.tidy_print(tibble)
  assert_true(output.contains("node"))
  assert_true(output.contains("parent"))
}

///|
test "tidytree example tibble" {
  let tibble = @src.tidy_example_tibble()
  assert_true(tibble.n_nodes > 0)
  assert_true(tibble.n_tips > 0)
}