///|
test "tree_io_parse_newick" {
  let newick = "(A:0.1,B:0.2,(C:0.3,D:0.4):0.5);"
  let tree = @src.parse_newick(newick)

  assert_true(tree.children.length() > 0)
}

///|
test "tree_io_write_newick" {
  let newick = "(A,B);"
  let tree = @src.parse_newick(newick)
  let result = @src.write_newick(tree, false)

  assert_true(result.length() > 0)
}

///|
test "tree_io_get_leaves" {
  let newick = "(A,B,(C,D));"
  let tree = @src.parse_newick(newick)
  let leaves = @src.get_leaves(tree)

  assert_eq(leaves.length(), 4)
}

///|
test "tree_io_get_leaf_count" {
  let newick = "(A,B,(C,D));"
  let tree = @src.parse_newick(newick)

  assert_eq(@src.get_leaf_count(tree), 4)
}

///|
test "tree_io_get_node_count" {
  let newick = "(A,B,(C,D));"
  let tree = @src.parse_newick(newick)

  assert_eq(@src.get_node_count(tree), 5)
}

///|
test "tree_io_find_node" {
  let newick = "(A,B,(C,D));"
  let tree = @src.parse_newick(newick)

  match @src.find_node(tree, "C") {
    Some(node) => assert_eq(node.name, "C")
    None => assert_true(false)
  }
}

///|
test "tree_io_total_branch_length" {
  let newick = "(A:0.1,B:0.2);"
  let tree = @src.parse_newick(newick)
  let total = @src.total_branch_length(tree)

  assert_true(total > 0.0)
}

///|
test "tree_io_get_tree_depth" {
  let newick = "(A:0.1,(B:0.2,C:0.3):0.4);"
  let tree = @src.parse_newick(newick)
  let depth = @src.get_tree_depth(tree)

  assert_true(depth > 0.0)
}
