///|
/// Tests for the CDAO (Comparative Data Analysis Ontology) format module.
/// Tests CDAO namespace constants, data structures, parsing, writing, and conversion.

// ============================================================================
// Namespace Constants
// ============================================================================

///|
test "cdao_namespace_constant" {
  let ns = @src.cdao_namespace()
  assert_eq(ns, "http://purl.obolibrary.org/obo/cdao.owl#")
}

///|
test "cdao_rdf_namespace_constant" {
  let ns = @src.cdao_rdf_namespace()
  assert_eq(ns, "http://www.w3.org/1999/02/22-rdf-syntax-ns#")
}

///|
test "cdao_rdfs_namespace_constant" {
  let ns = @src.cdao_rdfs_namespace()
  assert_eq(ns, "http://www.w3.org/2000/01/rdf-schema#")
}

// ============================================================================
// Data Structure Construction
// ============================================================================

///|
test "cdao_tree_construction" {
  let t = @src.Cdaotree::new("#tree1", rooted=true, root_node_id="#node1", name=Some("MyTree"))
  assert_eq(t.id, "#tree1")
  assert_eq(t.rooted, true)
  assert_eq(t.root_node_id, "#node1")
  assert_eq(t.name.unwrap(), "MyTree")
}

///|
test "cdao_tree_default_values" {
  let t = @src.Cdaotree::new("#tree2")
  assert_eq(t.id, "#tree2")
  assert_eq(t.rooted, false)
  assert_eq(t.root_node_id, "")
  assert_eq(t.name.is_some(), false)
}

///|
test "cdao_node_construction" {
  let n = @src.CdaoNode::new("#node1", children=["#node2", "#node3"], parent_id=Some("#node0"))
  assert_eq(n.id, "#node1")
  assert_eq(n.children.length(), 2)
  assert_eq(n.parent_id.unwrap(), "#node0")
}

///|
test "cdao_node_default_values" {
  let n = @src.CdaoNode::new("#node2")
  assert_eq(n.id, "#node2")
  assert_eq(n.children.length(), 0)
  assert_eq(n.parent_id.is_some(), false)
  assert_eq(n.tu_id.is_some(), false)
  assert_eq(n.branch_length.is_some(), false)
}

///|
test "cdao_tu_construction" {
  let tu = @src.CdaoTU::new("#tu1", label=Some("Homo sapiens"))
  assert_eq(tu.id, "#tu1")
  assert_eq(tu.label.unwrap(), "Homo sapiens")
}

///|
test "cdao_tu_default_values" {
  let tu = @src.CdaoTU::new("#tu2")
  assert_eq(tu.id, "#tu2")
  assert_eq(tu.label.is_some(), false)
}

///|
test "cdao_document_new_empty" {
  let doc = @src.CdaoDocument::new()
  assert_eq(doc.trees.length(), 0)
  assert_eq(doc.nodes.size(), 0)
  assert_eq(doc.tus.size(), 0)
}

// ============================================================================
// Sample Data Tests
// ============================================================================

///|
test "cdao_sample_document_structure" {
  let doc = @src.cdao_sample_document()
  assert_eq(doc.trees.length(), 1)
  assert_eq(doc.trees[0].id, "#tree1")
  assert_eq(doc.trees[0].rooted, true)
  assert_eq(doc.trees[0].root_node_id, "#node1")
  // 7 nodes (1 root + 2 internal + 4 leaves)
  assert_eq(doc.nodes.size(), 7)
  // 4 TUs
  assert_eq(doc.tus.size(), 4)
}

///|
test "cdao_sample_document_tu_labels" {
  let doc = @src.cdao_sample_document()
  assert_eq(doc.tus["#tu1"].label.unwrap(), "A")
  assert_eq(doc.tus["#tu2"].label.unwrap(), "B")
  assert_eq(doc.tus["#tu3"].label.unwrap(), "C")
  assert_eq(doc.tus["#tu4"].label.unwrap(), "D")
}

///|
test "cdao_sample_document_node_relationships" {
  let doc = @src.cdao_sample_document()
  // Root node has 2 children
  let root = doc.nodes["#node1"]
  assert_eq(root.children.length(), 2)
  assert_eq(root.children[0], "#node2")
  assert_eq(root.children[1], "#node3")
  // Internal node #node2 has 2 children and parent #node1
  let n2 = doc.nodes["#node2"]
  assert_eq(n2.children.length(), 2)
  assert_eq(n2.parent_id.unwrap(), "#node1")
  // Leaf node #node4 has TU
  let n4 = doc.nodes["#node4"]
  assert_eq(n4.children.length(), 0)
  assert_eq(n4.tu_id.unwrap(), "#tu1")
  assert_eq(n4.branch_length.unwrap(), 1.0)
}

///|
test "cdao_sample_xml_contains_key_elements" {
  let xml = @src.cdao_sample_xml()
  assert_true(xml.contains("rdf:RDF"))
  assert_true(xml.contains("cdao:Tree"))
  assert_true(xml.contains("cdao:Node"))
  assert_true(xml.contains("cdao:TU"))
  assert_true(xml.contains("cdao:has_Root"))
  assert_true(xml.contains("cdao:has_Child"))
  assert_true(xml.contains("cdao:belongs_to_TU"))
  assert_true(xml.contains("rdfs:label"))
}

// ============================================================================
// CDAO Parsing Tests
// ============================================================================

///|
test "cdao_parse_sample_xml" {
  let xml = @src.cdao_sample_xml()
  let doc = @src.cdao_parse(xml)
  assert_eq(doc.trees.length(), 1)
  assert_eq(doc.trees[0].id, "#tree1")
  assert_eq(doc.trees[0].root_node_id, "#node1")
  // Should have parsed nodes
  assert_true(doc.nodes.size() > 0)
}

///|
test "cdao_parse_extracts_tu_labels" {
  let xml = @src.cdao_sample_xml()
  let doc = @src.cdao_parse(xml)
  // Check that at least some TUs were parsed
  assert_true(doc.tus.size() > 0)
  // Check for label A
  let mut found_a = false
  let tu_keys = doc.tus.keys()
  let tu_arr : Array[String] = tu_keys.collect()
  for tu_id in tu_arr {
    let tu = doc.tus[tu_id]
    match tu.label {
      Some(l) => if l == "A" { found_a = true }
      None => ()
    }
  }
  assert_eq(found_a, true)
}

///|
test "cdao_parse_extracts_tree_root" {
  let xml = @src.cdao_sample_xml()
  let doc = @src.cdao_parse(xml)
  assert_eq(doc.trees.length(), 1)
  assert_eq(doc.trees[0].root_node_id, "#node1")
}

///|
test "cdao_parse_node_children" {
  let xml = @src.cdao_sample_xml()
  let doc = @src.cdao_parse(xml)
  // Root node should have children
  let root = doc.nodes.get("#node1")
  assert_eq(root.is_some(), true)
  let r = root.unwrap()
  assert_eq(r.children.length(), 2)
}

// ============================================================================
// CDAO to Tree Conversion Tests
// ============================================================================

///|
test "cdao_to_trees_from_sample_document" {
  let doc = @src.cdao_sample_document()
  let trees = @src.cdao_to_trees(doc)
  assert_eq(trees.length(), 1)
  let tree = trees[0]
  assert_eq(tree.rooted, true)
  // Root should have 2 children
  assert_eq(tree.root.clades.length(), 2)
}

///|
test "cdao_to_trees_preserves_leaf_labels" {
  let doc = @src.cdao_sample_document()
  let trees = @src.cdao_to_trees(doc)
  let tree = trees[0]
  // Check that leaf labels are preserved
  let left = tree.root.clades[0]
  let right = tree.root.clades[1]
  assert_eq(left.clades.length(), 2)
  assert_eq(right.clades.length(), 2)
  // Left subtree leaves: A and B
  assert_eq(left.clades[0].name.unwrap(), "A")
  assert_eq(left.clades[1].name.unwrap(), "B")
  // Right subtree leaves: C and D
  assert_eq(right.clades[0].name.unwrap(), "C")
  assert_eq(right.clades[1].name.unwrap(), "D")
}

///|
test "cdao_to_trees_preserves_branch_lengths" {
  let doc = @src.cdao_sample_document()
  let trees = @src.cdao_to_trees(doc)
  let tree = trees[0]
  // Leaf nodes should have branch_length = 1.0
  let left_leaf = tree.root.clades[0].clades[0]
  assert_eq(left_leaf.branch_length.unwrap(), 1.0)
}

// ============================================================================
// CDAO Writer Tests
// ============================================================================

///|
test "cdao_write_produces_valid_xml" {
  let tree = @src.parsimony_sample_tree()
  let xml = @src.cdao_write(tree)
  assert_true(xml.contains("<?xml"))
  assert_true(xml.contains("rdf:RDF"))
  assert_true(xml.contains("</rdf:RDF>"))
}

///|
test "cdao_write_contains_tree_element" {
  let tree = @src.parsimony_sample_tree()
  let xml = @src.cdao_write(tree)
  assert_true(xml.contains("cdao:Tree"))
  assert_true(xml.contains("cdao:has_Root"))
}

///|
test "cdao_write_contains_node_elements" {
  let tree = @src.parsimony_sample_tree()
  let xml = @src.cdao_write(tree)
  assert_true(xml.contains("cdao:Node"))
  assert_true(xml.contains("cdao:has_Child"))
}

///|
test "cdao_write_contains_tu_elements" {
  let tree = @src.parsimony_sample_tree()
  let xml = @src.cdao_write(tree)
  assert_true(xml.contains("cdao:TU"))
  assert_true(xml.contains("rdfs:label"))
  // The sample tree has leaves A, B, C, D
  assert_true(xml.contains("A"))
  assert_true(xml.contains("B"))
  assert_true(xml.contains("C"))
  assert_true(xml.contains("D"))
}

///|
test "cdao_write_roundtrip" {
  // Write a tree to CDAO, parse it back, and verify structure
  let tree = @src.parsimony_sample_tree()
  let xml = @src.cdao_write(tree)
  let doc = @src.cdao_parse(xml)
  let trees = @src.cdao_to_trees(doc)
  // Should have at least one tree
  assert_eq(trees.length(), 1)
}

// ============================================================================
// CDAO Document Summary Tests
// ============================================================================

///|
test "cdao_document_summary_format" {
  let doc = @src.cdao_sample_document()
  let summary = @src.cdao_document_summary(doc)
  assert_true(summary.contains("CDAO Document Summary"))
  assert_true(summary.contains("Trees: 1"))
  assert_true(summary.contains("Nodes: 7"))
  assert_true(summary.contains("TUs: 4"))
}

///|
test "cdao_document_summary_empty" {
  let doc = @src.CdaoDocument::new()
  let summary = @src.cdao_document_summary(doc)
  assert_true(summary.contains("Trees: 0"))
  assert_true(summary.contains("Nodes: 0"))
  assert_true(summary.contains("TUs: 0"))
}

// ============================================================================
// Edge Cases
// ============================================================================

///|
test "cdao_parse_empty_xml" {
  let doc = @src.cdao_parse("")
  assert_eq(doc.trees.length(), 0)
  assert_eq(doc.nodes.size(), 0)
}

///|
test "cdao_to_trees_empty_document" {
  let doc = @src.CdaoDocument::new()
  let trees = @src.cdao_to_trees(doc)
  assert_eq(trees.length(), 0)
}

///|
test "cdao_node_fields_via_constructor" {
  let n = @src.CdaoNode::new("#test", parent_id=Some("#parent"), tu_id=Some("#tu1"), branch_length=Some(2.5), label=Some("TestNode"))
  assert_eq(n.parent_id.unwrap(), "#parent")
  assert_eq(n.tu_id.unwrap(), "#tu1")
  assert_eq(n.branch_length.unwrap(), 2.5)
  assert_eq(n.label.unwrap(), "TestNode")
}

///|
test "cdao_node_children_via_constructor" {
  let n = @src.CdaoNode::new("#parent", children=["#child1", "#child2"])
  assert_eq(n.children.length(), 2)
  assert_eq(n.children[0], "#child1")
  assert_eq(n.children[1], "#child2")
}

///|
test "cdao_tu_label_via_constructor" {
  let tu = @src.CdaoTU::new("#tu1", label=Some("E. coli"))
  assert_eq(tu.label.unwrap(), "E. coli")
}
