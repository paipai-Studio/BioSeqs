///|
test "phylo_newick_parse_simple" {
  let newick = "(A,B);"
  let tree = @src.parse_newick(newick)

  assert_eq(@src.get_leaf_count(tree), 2)
  let leaves = @src.get_leaves(tree)
  assert_eq(leaves.length(), 2)
}

///|
test "phylo_newick_parse_nested" {
  let newick = "(((A,B),(C,D)),(E,F,G));"
  let tree = @src.parse_newick(newick)

  assert_eq(@src.get_leaf_count(tree), 7)
}

///|
test "phylo_newick_parse_with_branch_lengths" {
  let newick = "(A:0.1,B:0.2,(C:0.3,D:0.4):0.5);"
  let tree = @src.parse_newick(newick)

  assert_eq(@src.get_leaf_count(tree), 4)

  let found = @src.find_node(tree, "A")
  assert_eq(found.is_some(), true)
}

///|
test "phylo_newick_parse_no_semicolon" {
  let newick = "(A,B)"
  let tree = @src.parse_newick(newick)

  assert_eq(@src.get_leaf_count(tree), 2)
}

///|
test "phylo_newick_roundtrip" {
  let original = "(A:0.1,B:0.2,(C:0.3,D:0.4):0.5);"
  let tree = @src.parse_newick(original)
  let result = @src.write_newick(tree, true)

  let tree2 = @src.parse_newick(result)
  assert_eq(@src.get_leaf_count(tree2), 4)
}

///|
test "phylo_get_leaves" {
  let newick = "((A,B),(C,D));"
  let tree = @src.parse_newick(newick)
  let leaves = @src.get_leaves(tree)

  assert_eq(leaves.length(), 4)
}

///|
test "phylo_get_node_count" {
  let newick = "((A,B),(C,D));"
  let tree = @src.parse_newick(newick)

  assert_eq(@src.get_node_count(tree), 5)
}

///|
test "phylo_get_tree_depth" {
  let newick = "(A:0.1,B:0.2);"
  let tree = @src.parse_newick(newick)

  let depth = @src.get_tree_depth(tree)
  assert_true(depth > 0.0)
}

///|
test "phylo_total_branch_length" {
  let newick = "(A:0.1,B:0.2);"
  let tree = @src.parse_newick(newick)

  let total = @src.total_branch_length(tree)
  assert_true(total > 0.0)
}

///|
test "phylo_find_node" {
  let newick = "((A,B),(C,D));"
  let tree = @src.parse_newick(newick)

  let found = @src.find_node(tree, "A")
  assert_eq(found.is_some(), true)

  let not_found = @src.find_node(tree, "Z")
  assert_eq(not_found.is_none(), true)
}

///|
test "phylo_prune_tree" {
  let newick = "((A,B),(C,D));"
  let tree = @src.parse_newick(newick)

  let pruned = @src.prune_tree(tree, ["A", "B"])
  assert_eq(pruned.is_some(), true)

  match pruned {
    Some(p) => assert_eq(@src.get_leaf_count(p), 2)
    None => assert_true(false)
  }
}

///|
test "phylo_parse_nhx" {
  let nhx = "(A[&&NHX:name=A],B[&&NHX:name=B]);"
  let tree = @src.parse_nhx(nhx)

  assert_eq(@src.get_leaf_count(tree), 2)
}
