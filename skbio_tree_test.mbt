///|
/// Tests for scikit-bio Tree module.

///|
test "tree_node_creation" {
  let node = TreeNode::new(name=Some("root"))
  assert_eq(node.name, Some("root"))
  assert_true(node.is_tip())
}

///|
test "tree_node_with_children" {
  let child1 = TreeNode::new(name=Some("child1"))
  let child2 = TreeNode::new(name=Some("child2"))
  let children = Array::new()
  children.push(child1)
  children.push(child2)
  let parent = TreeNode::new(name=Some("parent"), children=children)
  
  assert_eq(parent.name, Some("parent"))
  assert_eq(parent.num_children(), 2)
  assert_true(parent.has_children())
  assert_true(!parent.is_tip())
}

///|
test "newick_parse_simple" {
  let newick = "(A,B);"
  let tree = TreeNode::from_newick(newick)
  
  assert_eq(tree.children.length(), 2)
  assert_eq(tree.children[0].name, Some("A"))
  assert_eq(tree.children[1].name, Some("B"))
}

///|
test "newick_parse_with_lengths" {
  let newick = "(A:0.1,B:0.2)root:0.3;"
  let tree = TreeNode::from_newick(newick)
  
  assert_eq(tree.name, Some("root"))
  assert_eq(tree.length, Some(0.3))
  assert_eq(tree.children[0].length, Some(0.1))
  assert_eq(tree.children[1].length, Some(0.2))
}

///|
test "newick_parse_nested" {
  let newick = "((A,B)C,D)E;"
  let tree = TreeNode::from_newick(newick)
  
  assert_eq(tree.name, Some("E"))
  assert_eq(tree.children.length(), 2)
  assert_eq(tree.children[0].name, Some("C"))
  assert_eq(tree.children[0].children.length(), 2)
  assert_eq(tree.children[0].children[0].name, Some("A"))
  assert_eq(tree.children[0].children[1].name, Some("B"))
  assert_eq(tree.children[1].name, Some("D"))
}

///|
test "newick_to_newick" {
  let original = "(A:0.1,B:0.2)root;"
  let tree = TreeNode::from_newick(original)
  let output = tree.to_newick()
  
  assert_true(output.has_suffix(";"))
}

///|
test "preorder_traversal" {
  let newick = "((A,B)C,D)E;"
  let tree = TreeNode::from_newick(newick)
  
  let nodes = tree.preorder()
  assert_eq(nodes.length(), 5)
  assert_eq(nodes[0].name, Some("E"))
}

///|
test "postorder_traversal" {
  let newick = "((A,B)C,D)E;"
  let tree = TreeNode::from_newick(newick)
  
  let nodes = tree.postorder()
  assert_eq(nodes.length(), 5)
}

///|
test "levelorder_traversal" {
  let newick = "((A,B)C,D)E;"
  let tree = TreeNode::from_newick(newick)
  
  let nodes = tree.levelorder()
  assert_eq(nodes.length(), 5)
}

///|
test "tips" {
  let newick = "((A,B)C,D)E;"
  let tree = TreeNode::from_newick(newick)
  
  let tips = tree.tips()
  assert_eq(tips.length(), 3)
}

///|
test "count_tips" {
  let newick = "((A,B)C,D)E;"
  let tree = TreeNode::from_newick(newick)
  
  assert_eq(tree.count_tips(), 3)
}

///|
test "count_nodes" {
  let newick = "((A,B)C,D)E;"
  let tree = TreeNode::from_newick(newick)
  
  assert_eq(tree.count(), 5)
}

///|
test "find_node" {
  let newick = "((A,B)C,D)E;"
  let tree = TreeNode::from_newick(newick)
  
  let found = tree.find("C")
  assert_true(found is Some(_))
  assert_eq(found.unwrap().name, Some("C"))
}

///|
test "find_not_found" {
  let newick = "((A,B)C,D)E;"
  let tree = TreeNode::from_newick(newick)
  
  let found = tree.find("X")
  assert_true(found is None)
}

///|
test "total_length" {
  let newick = "(A:0.1,B:0.2)root:0.3;"
  let tree = TreeNode::from_newick(newick)
  
  let total = tree.total_length()
  assert_true((total - 0.6).abs() < 0.0001)
}

///|
test "lca" {
  let newick = "((A,B)C,D)E;"
  let tree = TreeNode::from_newick(newick)
  
  let names = Array::new()
  names.push("A")
  names.push("B")
  let ancestor = tree.lca(names)
  
  assert_true(ancestor is Some(_))
  assert_eq(ancestor.unwrap().name, Some("C"))
}

///|
test "distance" {
  let newick = "(A:0.1,B:0.2)C:0.3;"
  let tree = TreeNode::from_newick(newick)
  
  let dist = tree.distance("A", "B")
  assert_true(dist is Some(_))
  assert_true((dist.unwrap() - 0.3).abs() < 0.0001)
}
