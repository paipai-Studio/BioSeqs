///|
/// Phylo.Consensus module tests

test "ConsensusNode::new" {
  let node = @src.ConsensusNode::new("A", false)
  assert_eq(node.name, "A")
  assert_true(node.is_leaf())
}

test "ConsensusNode::add_child" {
  let parent = @src.ConsensusNode::new("internal", true)
  let child = @src.ConsensusNode::new("A", false)
  let parent_with_child = parent.add_child(child)
  assert_eq(parent_with_child.children.length(), 1)
}

test "ConsensusNode::get_leaves" {
  let root = @src.ConsensusNode::new("root", true)
  let child1 = @src.ConsensusNode::new("A", false)
  let child2 = @src.ConsensusNode::new("B", false)
  root.children.push(child1)
  root.children.push(child2)
  let leaves = root.get_leaves()
  assert_eq(leaves.length(), 2)
}

test "parse_newick" {
  let newick = "((A,B),(C,D));"
  match @src.parse_consensus_newick(newick) {
    Some(root) => assert_true(root.is_internal)
    None => assert_true(false)
  }
}

test "newick_to_tree" {
  let tree = @src.newick_to_tree("((A,B),(C,D));")
  match tree.root {
    Some(root) => assert_eq(root.children.length(), 2)
    None => assert_true(false)
  }
}

test "ConsensusTree::to_newick" {
  let tree = @src.newick_to_tree("((A,B),(C,D));")
  let newick = tree.to_newick()
  assert_true(newick.length() > 0)
}

test "build_majority_consensus" {
  let trees = @src.create_example_trees()
  let consensus = @src.build_majority_consensus(trees)
  match consensus.root {
    Some(root) => assert_true(root.children.length() > 0)
    None => assert_true(false)
  }
}

test "build_strict_consensus" {
  let trees = @src.create_example_trees()
  let consensus = @src.build_strict_consensus(trees)
  match consensus.root {
    Some(_) => ()
    None => ()
  }
}

test "get_all_splits" {
  let tree = @src.newick_to_tree("((A,B),(C,D));")
  match tree.root {
    Some(root) => {
      let leaves = root.get_leaves()
      let splits = @src.get_all_splits(root, leaves)
      assert_true(splits.length() > 0)
    }
    None => assert_true(false)
  }
}

test "Split::normalize" {
  let split = @src.Split::new(["B", "A"])
  let normalized = split.normalize()
  assert_eq(normalized.taxa[0], "A")
  assert_eq(normalized.taxa[1], "B")
}

test "Split::hash" {
  let split = @src.Split::new(["A", "B"])
  let h = split.hash()
  assert_eq(h, "A,B")
}

test "calculate_consensus_support" {
  let trees = @src.create_example_trees()
  let split = @src.Split::new(["A", "B"])
  let support = @src.calculate_consensus_support(trees, split)
  assert_true(support >= 0.0 && support <= 1.0)
}

test "get_all_splits_from_trees" {
  let trees = @src.create_example_trees()
  let splits = @src.get_all_splits_from_trees(trees)
  assert_true(splits.length() > 0)
}

test "create_example_trees" {
  let trees = @src.create_example_trees()
  assert_eq(trees.length(), 5)
}

test "create_simple_tree" {
  let tree = @src.create_simple_tree()
  match tree.root {
    Some(root) => assert_eq(root.children.length(), 2)
    None => assert_true(false)
  }
}

test "build_consensus_empty" {
  let trees : Array[@src.ConsensusTree] = Array::new()
  let consensus = @src.build_consensus(trees, 0.5)
  match consensus.root {
    Some(_) => assert_true(false)
    None => ()
  }
}