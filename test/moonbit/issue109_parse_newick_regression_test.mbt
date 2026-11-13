// Regression tests for GitHub Issue #109: parse_newick was flattening
// nested topology and dropping internal branch lengths.

test "issue109_nested_topology_two_children" {
  let tree = @src.parse_newick("(A:0.1,(B:0.2,C:0.3):0.4);")
  assert_eq(tree.children.length(), 2)
  // Root has two children: leaf A and internal node (B,C).
  let a = tree.children[0]
  let internal_bc = tree.children[1]
  assert_eq(a.name, "A")
  assert_true((a.branch_length - 0.1).abs() < 0.000000001)
  assert_eq(internal_bc.children.length(), 2)
  assert_true((internal_bc.branch_length - 0.4).abs() < 0.000000001)
  // Verify total branch length: 0.1 + 0.4 + 0.2 + 0.3 = 1.0
  let total = @src.total_branch_length(tree)
  assert_true((total - 1.0).abs() < 0.000000001)
  // Max depth: 0.1 (A) vs 0.4 + 0.3 = 0.7 (C) — so 0.7
  let depth = @src.get_tree_depth(tree)
  assert_true((depth - 0.7).abs() < 0.000000001)
}

test "issue109_roundtrip_write" {
  let tree = @src.parse_newick("(A:0.1,(B:0.2,C:0.3):0.4);")
  let written = @src.write_newick(tree, true)
  // Re-parse the serialized form and check invariants again.
  let reparsed = @src.parse_newick(written + ";")
  assert_eq(reparsed.children.length(), 2)
  let total = @src.total_branch_length(reparsed)
  assert_true((total - 1.0).abs() < 0.000000001)
}

test "issue109_old_tests_still_pass_parseable" {
  let tree = @src.parse_newick("(A:0.1,B:0.2,(C:0.3,D:0.4):0.5);")
  // Now root has 3 children: A, B, and internal (C,D)
  assert_eq(tree.children.length(), 3)
  assert_eq(@src.get_leaf_count(tree), 4)
  // Nodes: root + A + B + internal(C,D) + C + D = 6
  assert_eq(@src.get_node_count(tree), 6)
  // Leaves depth: A=0.1, B=0.2, C=0.5+0.3=0.8, D=0.5+0.4=0.9 -> max 0.9
  let depth = @src.get_tree_depth(tree)
  assert_true((depth - 0.9).abs() < 0.000000001)
  // Total: 0.1+0.2+0.5+0.3+0.4 = 1.5
  let total = @src.total_branch_length(tree)
  assert_true((total - 1.5).abs() < 0.000000001)
}

test "issue109_simple_two_leaves" {
  let tree = @src.parse_newick("(A,B);")
  assert_eq(tree.children.length(), 2)
  assert_eq(@src.get_leaf_count(tree), 2)
}

test "issue109_unparenthesized_single_root" {
  let tree = @src.parse_newick("F:0.5;")
  assert_eq(tree.children.length(), 0)
  assert_eq(tree.name, "F")
  assert_true((tree.branch_length - 0.5).abs() < 0.000000001)
}

test "issue109_root_label_and_length" {
  let tree = @src.parse_newick("(A,B)C:0.3;")
  assert_eq(tree.children.length(), 2)
  assert_eq(tree.name, "C")
  assert_true((tree.branch_length - 0.3).abs() < 0.000000001)
}
