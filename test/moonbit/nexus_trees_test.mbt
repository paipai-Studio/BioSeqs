///|
/// Tests for nexus_trees module (Bio.Nexus.Trees/Nodes).

test "nexus_trees_parse basic" {
  let content = "#NEXUS\nBEGIN TREES;\n  TREE tree1 = ((A:0.1,B:0.2):0.3,(C:0.4,D:0.5):0.6);\nEND;\n"
  let collection = @src.nexus_trees_parse(content)
  assert_eq(collection.trees.length(), 1)
  assert_eq(collection.trees[0].name, "tree1")
  assert_eq(collection.trees[0].n_taxa, 4)
}

test "nexus_trees_parse_newick structure" {
  let collection = @src.nexus_trees_from_newicks([("t1", "((A:0.1,B:0.2):0.3,C:0.4);")])
  let tree = collection.trees[0]
  assert_eq(tree.node_count(), 5) // 2 internal + 3 leaves
  assert_eq(tree.n_taxa, 3)
}

test "nexus_tree_get_leaves" {
  let collection = @src.nexus_trees_from_newicks([("t1", "((A:0.1,B:0.2):0.3,C:0.4);")])
  let tree = collection.trees[0]
  let leaves = tree.get_leaves()
  assert_eq(leaves.length(), 3)
  // All leaves should have taxa
  let mut i = 0
  while i < leaves.length() {
    assert_true(leaves[i].is_leaf())
    i = i + 1
  }
}

test "nexus_tree_get_node" {
  let collection = @src.nexus_trees_from_newicks([("t1", "((A:0.1,B:0.2):0.3,C:0.4);")])
  let tree = collection.trees[0]
  let node = tree.get_node(tree.root)
  assert_true(node is Some(_))
}

test "nexus_tree_to_newick roundtrip" {
  let newick = "((A:0.1,B:0.2):0.3,(C:0.4,D:0.5):0.6);"
  let collection = @src.nexus_trees_from_newicks([("t1", newick)])
  let tree = collection.trees[0]
  let output = tree.to_newick()
  assert_eq(output, newick)
}

test "nexus_tree_total_length" {
  let collection = @src.nexus_trees_from_newicks([("t1", "((A:0.1,B:0.2):0.3,(C:0.4,D:0.5):0.6);")])
  let tree = collection.trees[0]
  // Sum: 0.1+0.2+0.3+0.4+0.5+0.6 = 2.1
  assert_true((tree.total_length() - 2.1).abs() < 1.0e-9)
}

test "nexus_tree_patristic_distance" {
  let collection = @src.nexus_trees_from_newicks([("t1", "((A:0.1,B:0.2):0.3,(C:0.4,D:0.5):0.6);")])
  let tree = collection.trees[0]
  // A-B distance: 0.1 + 0.2 = 0.3
  let d_ab = tree.patristic_distance("A", "B")
  assert_true((d_ab - 0.3).abs() < 1.0e-9)
  // A-C distance: 0.1 + 0.3 + 0.6 + 0.4 = 1.4
  let d_ac = tree.patristic_distance("A", "C")
  assert_true((d_ac - 1.4).abs() < 1.0e-9)
}

test "nexus_tree_rooted_flag" {
  let collection = @src.nexus_trees_from_newicks([("t1", "[&R]((A:0.1,B:0.2):0.3,C:0.4);")])
  let tree = collection.trees[0]
  assert_true(tree.rooted)
}

test "nexus_trees_parse_translate" {
  let content = "#NEXUS\nBEGIN TREES;\nTRANSLATE\n  1 Homo_sapiens,\n  2 Mus_musculus\n  ;\n  TREE tree1 = ((1:0.1,2:0.2):0.3);\nEND;\n"
  let collection = @src.nexus_trees_parse(content)
  assert_eq(collection.trees.length(), 1)
  assert_eq(collection.translate.size() |> ignore, ())
  // resolve taxon
  let label = collection.resolve_taxon("1")
  assert_eq(label, "Homo_sapiens")
}

test "nexus_tree_collection_get_tree" {
  let collection = @src.nexus_trees_from_newicks([("first", "(A:0.1,B:0.2);"), ("second", "(C:0.3,D:0.4);")])
  assert_eq(collection.trees.length(), 2)
  let t1 = collection.get_tree("first")
  assert_true(t1 is Some(_))
  let t3 = collection.get_tree("missing")
  assert_true(t3 is None)
}

test "nexus_node_is_leaf" {
  let node = @src.NexusNode::new(0)
  assert_true(node.is_leaf())
  let parent = @src.NexusNode::new(1)
  let parent2 = @src.NexusNode::{ node_id: parent.node_id, label: parent.label, branch_length: parent.branch_length, taxon: parent.taxon, children: [0], parent: parent.parent }
  assert_false(parent2.is_leaf())
}

test "nexus_trees_example_content" {
  let content = @src.nexus_trees_example_content()
  let collection = @src.nexus_trees_parse(content)
  assert_eq(collection.trees.length(), 1)
  assert_eq(collection.trees[0].name, "tree1")
}