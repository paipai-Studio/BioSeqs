///|
/// Tests for phylo_nexml module.

// ============================================================================
// NeXML data structure creation tests
// ============================================================================

///|
test "nexml_state_creation" {
  let state = @src.NeXMLState::new("s1", "A")
  assert_eq(state.id, "s1")
  assert_eq(state.symbol, "A")
}

///|
test "nexml_character_creation" {
  let states : Array[String] = ["A", "C", "G", "T"]
  let char = @src.NeXMLCharacter::new("c1", "pos1", states)
  assert_eq(char.id, "c1")
  assert_eq(char.label, "pos1")
  assert_eq(char.states.length(), 4)
}

///|
test "nexml_otu_creation" {
  let otu = @src.NeXMLOTU::new("otu1", "Homo_sapiens")
  assert_eq(otu.id, "otu1")
  assert_eq(otu.label, "Homo_sapiens")
}

///|
test "nexml_otus_creation" {
  let otus_arr : Array[@src.NeXMLOTU] = []
  let otus = @src.NeXMLOTUs::new("otus1", otus_arr)
  assert_eq(otus.id, "otus1")
  assert_eq(otus.otus.length(), 0)
}

///|
test "nexml_node_creation" {
  let node = @src.NeXMLNode::new("n1", "Node1", "otu1")
  assert_eq(node.id, "n1")
  assert_eq(node.label, "Node1")
  assert_eq(node.otu, "otu1")
}

///|
test "nexml_edge_creation" {
  let edge = @src.NeXMLEdge::new("e1", "n1", "n2")
  assert_eq(edge.id, "e1")
  assert_eq(edge.source, "n1")
  assert_eq(edge.target, "n2")
}

///|
test "nexml_tree_creation" {
  let nodes : Array[@src.NeXMLNode] = []
  let edges : Array[@src.NeXMLEdge] = []
  let tree = @src.NeXMLTree::new(
    "tree1",
    "Test Tree",
    "FloatTree",
    true,
    nodes,
    edges,
  )
  assert_eq(tree.id, "tree1")
  assert_eq(tree.name, "Test Tree")
  assert_eq(tree.tree_type, "FloatTree")
  assert_true(tree.rooted)
}

///|
test "nexml_trees_creation" {
  let trees_arr : Array[@src.NeXMLTree] = []
  let trees = @src.NeXMLTrees::new("trees1", "otus1", trees_arr)
  assert_eq(trees.id, "trees1")
  assert_eq(trees.otus_id, "otus1")
  assert_eq(trees.trees.length(), 0)
}

///|
test "nexml_document_creation" {
  let doc = @src.NeXMLDocument::new(version="0.9", generator="BioSeqs")
  assert_eq(doc.version, "0.9")
  assert_eq(doc.generator, "BioSeqs")
  assert_eq(doc.otus_blocks.length(), 0)
  assert_eq(doc.characters_blocks.length(), 0)
  assert_eq(doc.trees_blocks.length(), 0)
}

// ============================================================================
// NeXML serialization tests
// ============================================================================

///|
test "nexml_to_xml_simple" {
  let doc = @src.NeXMLDocument::new(version="0.9", generator="TestGen")
  let xml = @src.nexml_to_xml(doc)
  assert_true(xml.length() > 0)
  assert_true(xml.has_prefix("<?xml"))
}

// ============================================================================
// Utility function tests
// ============================================================================

///|
test "nexml_get_otu_count" {
  let doc = @src.NeXMLDocument::new()
  assert_eq(@src.nexml_get_otu_count(doc), 0)
}

///|
test "nexml_get_tree_count" {
  let doc = @src.NeXMLDocument::new()
  assert_eq(@src.nexml_get_tree_count(doc), 0)
}

///|
test "nexml_get_all_otu_labels" {
  let doc = @src.NeXMLDocument::new()
  let labels = @src.nexml_get_all_otu_labels(doc)
  assert_eq(labels.length(), 0)
}

///|
test "nexml_get_tree_names" {
  let doc = @src.NeXMLDocument::new()
  let names = @src.nexml_get_tree_names(doc)
  assert_eq(names.length(), 0)
}

// ============================================================================
// Tree to Newick conversion tests
// ============================================================================

///|
test "nexml_tree_to_newick_empty" {
  let nodes : Array[@src.NeXMLNode] = []
  let edges : Array[@src.NeXMLEdge] = []
  let tree = @src.NeXMLTree::new("t1", "Test", "FloatTree", true, nodes, edges)
  let newick = @src.nexml_tree_to_newick(tree)
  assert_eq(newick, ";")
}

///|
test "nexml_tree_to_newick_simple" {
  let n1 = @src.NeXMLNode::new("n1", "root", "")
  let n2 = @src.NeXMLNode::new("n2", "A", "")
  let n3 = @src.NeXMLNode::new("n3", "B", "")
  let nodes = [n1, n2, n3]
  let e1 = @src.NeXMLEdge::new("e1", "n1", "n2", length=Some(0.5))
  let e2 = @src.NeXMLEdge::new("e2", "n1", "n3", length=Some(0.3))
  let edges = [e1, e2]
  let tree = @src.NeXMLTree::new("t1", "Test", "FloatTree", true, nodes, edges)
  let newick = @src.nexml_tree_to_newick(tree)
  assert_true(newick.length() > 0)
  assert_true(newick.has_prefix("("))
  assert_true(newick.has_suffix(";"))
}
