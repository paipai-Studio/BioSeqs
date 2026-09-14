///|
/// Tests for treeio module.

///|
test "treeio parse newick simple" {
  let tree = @src.treeio_parse_newick("((A:0.1,B:0.2):0.3,C:0.4);")
  assert_true(tree.n_tips >= 2)
  assert_true(tree.n_nodes >= 3)
  assert_eq(tree.format, "newick")
}

///|
test "treeio parse newick with labels" {
  let tree = @src.treeio_parse_newick("((A:0.1,B:0.2):0.3,C:0.4);")
  assert_true(tree.tip_labels.contains("A"))
  assert_true(tree.tip_labels.contains("B"))
  assert_true(tree.tip_labels.contains("C"))
}

///|
test "treeio write newick" {
  let tree = @src.treeio_parse_newick("((A:0.1,B:0.2):0.3,C:0.4);")
  let newick = @src.treeio_write_newick(tree)
  assert_true(newick.contains("A"))
  assert_true(newick.contains("B"))
  assert_true(newick.contains("C"))
  assert_true(newick.contains(";"))
}

///|
test "treeio parse beast tree" {
  let tree = @src.treeio_parse_beast("((A[&height=0.1]:0.1,B[&height=0.1]:0.2)[&height=0.3]:0.3,C[&height=0.1]:0.4);")
  assert_eq(tree.format, "beast")
  assert_true(tree.n_tips >= 2)
}

///|
test "treeio parse jplace" {
  let jplace = @src.treeio_parse_jplace("{\"tree\":\"((A:0.1,B:0.2):0.3,C:0.4);\",\"placements\":[]}")
  assert_true(jplace.tree.contains("A"))
}

///|
test "treeio parse r8s" {
  let content = "tree1 = ((A:0.1,B:0.2):0.3,C:0.4);\ntree2 = ((A:0.2,B:0.3):0.4,C:0.5);"
  let trees = @src.treeio_parse_r8s(content)
  assert_eq(trees.length(), 2)
}

///|
test "treeio parse phylip" {
  let content = "3\nA 0.1\nB 0.2\nC 0.3\n((A:0.1,B:0.2):0.3,C:0.4);"
  let tree = @src.treeio_parse_phylip(content)
  assert_eq(tree.format, "phylip")
  assert_true(tree.n_tips >= 2)
}

///|
test "treeio tree stats" {
  let tree = @src.treeio_example_tree()
  let (n_tips, n_internal, n_edges) = @src.treeio_tree_stats(tree)
  assert_true(n_tips >= 2)
  assert_true(n_internal >= 1)
  assert_true(n_edges >= 2)
}

///|
test "treeio is binary" {
  let tree = @src.treeio_example_tree()
  assert_true(@src.treeio_is_binary(tree))
}

///|
test "treeio is rooted" {
  let tree = @src.treeio_example_tree()
  assert_true(@src.treeio_is_rooted(tree))
}

///|
test "treeio total length" {
  let tree = @src.treeio_example_tree()
  let total = @src.treeio_total_length(tree)
  // ((A:0.1,B:0.2):0.3,C:0.4) = 0.1+0.2+0.3+0.4 = 1.0
  assert_true((total - 1.0).abs() < 0.001)
}

///|
test "treeio max depth" {
  let tree = @src.treeio_example_tree()
  let depth = @src.treeio_max_depth(tree)
  assert_true(depth >= 2)
}

///|
test "treeio example tree" {
  let tree = @src.treeio_example_tree()
  assert_true(tree.n_tips >= 2)
  assert_true(tree.root >= 0)
}

///|
test "treeio example beast" {
  let tree = @src.treeio_example_beast()
  assert_eq(tree.format, "beast")
}