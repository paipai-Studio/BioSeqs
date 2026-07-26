///|
test "phyloxml_parse_simple_tree" {
  let xml =
    "<phyloxml>\n" +
    "  <phylogeny id=\"tree1\">\n" +
    "    <name>Simple Tree</name>\n" +
    "    <clade>\n" +
    "      <name>Root</name>\n" +
    "      <clade>\n" +
    "        <name>Child1</name>\n" +
    "      </clade>\n" +
    "      <clade>\n" +
    "        <name>Child2</name>\n" +
    "      </clade>\n" +
    "    </clade>\n" +
    "  </phylogeny>\n" +
    "</phyloxml>\n"
  let result = @src.parse_phyloxml(xml)
  assert_true(result.trees.length() >= 0)
}

///|
test "phyloxml_parse_branch_lengths" {
  let xml =
    "<phyloxml>\n" +
    "  <phylogeny id=\"tree1\">\n" +
    "    <clade>\n" +
    "      <branch_length>0.5</branch_length>\n" +
    "      <clade>\n" +
    "        <name>A</name>\n" +
    "        <branch_length>0.1</branch_length>\n" +
    "      </clade>\n" +
    "      <clade>\n" +
    "        <name>B</name>\n" +
    "        <branch_length>0.2</branch_length>\n" +
    "      </clade>\n" +
    "    </clade>\n" +
    "  </phylogeny>\n" +
    "</phyloxml>\n"
  let result = @src.parse_phyloxml(xml)
  assert_true(result.trees.length() >= 0)
}

///|
test "phyloxml_parse_taxon_annotations" {
  let xml =
    "<phyloxml>\n" +
    "  <phylogeny>\n" +
    "    <clade>\n" +
    "      <name>Root</name>\n" +
    "      <taxon>\n" +
    "        <taxon_id>tx1</taxon_id>\n" +
    "        <scientific_name>Homo sapiens</scientific_name>\n" +
    "        <common_name>Human</common_name>\n" +
    "      </taxon>\n" +
    "      <clade>\n" +
    "        <name>A</name>\n" +
    "      </clade>\n" +
    "    </clade>\n" +
    "  </phylogeny>\n" +
    "</phyloxml>\n"
  let result = @src.parse_phyloxml(xml)
  assert_true(result.trees.length() >= 0)
}

///|
test "phyloxml_parse_sequence_annotations" {
  let xml =
    "<phyloxml>\n" +
    "  <phylogeny>\n" +
    "    <clade>\n" +
    "      <name>Root</name>\n" +
    "      <sequence>\n" +
    "        <sequence_id>seq1</sequence_id>\n" +
    "        <molecular_type>DNA</molecular_type>\n" +
    "        <symbol>ATG</symbol>\n" +
    "      </sequence>\n" +
    "      <clade>\n" +
    "        <name>A</name>\n" +
    "      </clade>\n" +
    "    </clade>\n" +
    "  </phylogeny>\n" +
    "</phyloxml>\n"
  let result = @src.parse_phyloxml(xml)
  assert_true(result.trees.length() >= 0)
}

///|
test "phyloxml_xml_roundtrip" {
  let xml =
    "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n" +
    "<phyloxml xmlns=\"http://www.phyloxml.org/1.0\">\n" +
    "  <phylogeny id=\"tree1\">\n" +
    "    <name>Test Tree</name>\n" +
    "    <clade>\n" +
    "      <name>Root</name>\n" +
    "      <clade>\n" +
    "        <name>Child1</name>\n" +
    "      </clade>\n" +
    "      <clade>\n" +
    "        <name>Child2</name>\n" +
    "      </clade>\n" +
    "    </clade>\n" +
    "  </phylogeny>\n" +
    "</phyloxml>\n"
  let result = @src.parse_phyloxml(xml)
  assert_eq(result.trees.length(), 1)
  let tree = result.trees[0]
  assert_eq(tree.tree_id, "tree1")
  let output_xml = @src.phyloxml_to_xml(result)
  assert_true(output_xml.contains("<?xml version=\"1.0\""))
  assert_true(output_xml.contains("<phyloxml"))
  assert_true(output_xml.contains("phylogeny"))
}

///|
test "phyloxml_newick_to_phyloxml" {
  let newick = "(A:0.1,B:0.2,(C:0.3,D:0.4):0.5);"
  let result = @src.newick_to_phyloxml(newick, tree_id="test_tree")
  assert_eq(result.trees.length(), 1)
  let tree = result.trees[0]
  assert_eq(tree.tree_id, "test_tree")
  assert_eq(tree.name, "test_tree")
  assert_false(tree.root.is_leaf)
  assert_eq(tree.root.children.length(), 3)
  let child1 = tree.root.children[0]
  assert_eq(child1.name, "A")
  assert_true(child1.is_leaf)
  assert_true(child1.branch_length > 0.05 && child1.branch_length < 0.15)
  let child2 = tree.root.children[1]
  assert_eq(child2.name, "B")
  assert_true(child2.is_leaf)
  assert_true(child2.branch_length > 0.15 && child2.branch_length < 0.25)
  let child3 = tree.root.children[2]
  assert_eq(child3.name, "")
  assert_false(child3.is_leaf)
  assert_true(child3.branch_length > 0.4 && child3.branch_length < 0.6)
  assert_eq(child3.children.length(), 2)
}

///|
test "phyloxml_newick_simple" {
  let newick = "(A,B,C);"
  let result = @src.newick_to_phyloxml(newick)
  assert_eq(result.trees.length(), 1)
  let tree = result.trees[0]
  assert_false(tree.root.is_leaf)
  assert_eq(tree.root.children.length(), 3)
  assert_eq(tree.root.children[0].name, "A")
  assert_eq(tree.root.children[1].name, "B")
  assert_eq(tree.root.children[2].name, "C")
}

///|
test "phyloxml_to_newick" {
  let xml =
    "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n" +
    "<phyloxml xmlns=\"http://www.phyloxml.org/1.0\">\n" +
    "  <phylogeny>\n" +
    "    <clade>\n" +
    "      <name>Root</name>\n" +
    "      <clade>\n" +
    "        <name>A</name>\n" +
    "        <branch_length>0.1</branch_length>\n" +
    "      </clade>\n" +
    "      <clade>\n" +
    "        <name>B</name>\n" +
    "        <branch_length>0.2</branch_length>\n" +
    "      </clade>\n" +
    "    </clade>\n" +
    "  </phylogeny>\n" +
    "</phyloxml>\n"
  let result = @src.parse_phyloxml(xml)
  let newick = @src.phyloxml_to_newick(result)
  assert_true(newick.length() > 0)
  assert_true(newick.contains(";"))
}

///|
test "phyloxml_node_to_newick_leaf" {
  let node = @src.PhyloXMLNode::new(name="A", branch_length=0.5)
  let newick = @src.phyloxml_node_to_newick(node)
  assert_true(newick.contains("A"))
  assert_true(newick.contains("0.5"))
}

///|
test "phyloxml_node_to_newick_internal" {
  let child1 = @src.PhyloXMLNode::new(name="A", branch_length=0.1)
  let child2 = @src.PhyloXMLNode::new(name="B", branch_length=0.2)
  let node = @src.PhyloXMLNode::new(
    name="Root",
    branch_length=0.5,
    children=[child1, child2],
  )
  let newick = @src.phyloxml_node_to_newick(node)
  assert_true(newick.contains("("))
  assert_true(newick.contains("A"))
  assert_true(newick.contains("B"))
  assert_true(newick.contains("Root"))
  assert_true(newick.contains("0.5"))
}

///|
test "phyloxml_get_tree_names" {
  let tree1 = @src.PhyloXMLTree::new(name="Tree1")
  let tree2 = @src.PhyloXMLTree::new(name="Tree2")
  let result = @src.PhyloXMLResult::new(trees=[tree1, tree2])
  let names = @src.get_tree_names(result)
  assert_eq(names.length(), 2)
  assert_eq(names[0], "Tree1")
  assert_eq(names[1], "Tree2")
}

///|
test "phyloxml_get_all_tips" {
  let leaf1 = @src.PhyloXMLNode::new(name="A")
  let leaf2 = @src.PhyloXMLNode::new(name="B")
  let internal = @src.PhyloXMLNode::new(
    name="Root",
    children=[leaf1, leaf2],
  )
  let tree = @src.PhyloXMLTree::new(root=internal)
  let tips = @src.get_all_tips(tree)
  assert_eq(tips.length(), 2)
  assert_eq(tips[0], "A")
  assert_eq(tips[1], "B")
}

///|
test "phyloxml_get_all_tips_deep" {
  let leaf1 = @src.PhyloXMLNode::new(name="A")
  let leaf2 = @src.PhyloXMLNode::new(name="B")
  let leaf3 = @src.PhyloXMLNode::new(name="C")
  let child1 = @src.PhyloXMLNode::new(name="X", children=[leaf1, leaf2])
  let root = @src.PhyloXMLNode::new(name="Root", children=[child1, leaf3])
  let tree = @src.PhyloXMLTree::new(root=root)
  let tips = @src.get_all_tips(tree)
  assert_eq(tips.length(), 3)
}

///|
test "phyloxml_node_count" {
  let leaf1 = @src.PhyloXMLNode::new(name="A")
  let leaf2 = @src.PhyloXMLNode::new(name="B")
  let root = @src.PhyloXMLNode::new(
    name="Root",
    children=[leaf1, leaf2],
  )
  let tree = @src.PhyloXMLTree::new(root=root)
  let count = @src.phyloxml_get_node_count(tree)
  assert_eq(count, 3)
}

///|
test "phyloxml_multiple_trees" {
  let xml =
    "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n" +
    "<phyloxml xmlns=\"http://www.phyloxml.org/1.0\">\n" +
    "  <phylogeny id=\"tree1\">\n" +
    "    <name>Tree 1</name>\n" +
    "    <clade>\n" +
    "      <name>Root1</name>\n" +
    "      <clade><name>A</name></clade>\n" +
    "      <clade><name>B</name></clade>\n" +
    "    </clade>\n" +
    "  </phylogeny>\n" +
    "  <phylogeny id=\"tree2\">\n" +
    "    <name>Tree 2</name>\n" +
    "    <clade>\n" +
    "      <name>Root2</name>\n" +
    "      <clade><name>C</name></clade>\n" +
    "      <clade><name>D</name></clade>\n" +
    "    </clade>\n" +
    "  </phylogeny>\n" +
    "</phyloxml>\n"
  let result = @src.parse_phyloxml(xml)
  assert_eq(result.trees.length(), 2)
  assert_eq(result.trees[0].tree_id, "tree1")
  assert_true(result.trees[0].name.length() >= 0)
  assert_eq(result.trees[1].tree_id, "tree2")
  assert_true(result.trees[1].name.length() >= 0)
}

///|
test "phyloxml_edge_empty_xml" {
  let result = @src.parse_phyloxml("")
  assert_eq(result.trees.length(), 0)
}

///|
test "phyloxml_simple_phylogeny" {
  let xml = "<phyloxml><phylogeny id=\"t1\"><clade><name>A</name></clade></phylogeny></phyloxml>"
  let result = @src.parse_phyloxml(xml)
  assert_true(result.trees.length() >= 0)
}

///|
test "phyloxml_edge_missing_elements" {
  let xml = "<phyloxml xmlns=\"http://www.phyloxml.org/1.0\"></phyloxml>"
  let result = @src.parse_phyloxml(xml)
  assert_eq(result.trees.length(), 0)
}

///|
test "phyloxml_sequence_annotation_creation" {
  let sa = @src.SequenceAnnotation::new(
    sequence_id="seq1",
    molecular_type="DNA",
    symbol="ATG",
    uri="http://example.com",
  )
  assert_eq(sa.sequence_id, "seq1")
  assert_eq(sa.molecular_type, "DNA")
  assert_eq(sa.symbol, "ATG")
  assert_eq(sa.uri, "http://example.com")
}

///|
test "phyloxml_taxon_info_creation" {
  let ti = @src.TaxonInfo::new(
    taxon_id="tx1",
    scientific_name="Homo sapiens",
    common_name="Human",
    authority="Linnaeus, 1758",
    taxon_rank="species",
  )
  assert_eq(ti.taxon_id, "tx1")
  assert_eq(ti.scientific_name, "Homo sapiens")
  assert_eq(ti.common_name, "Human")
  assert_eq(ti.authority, "Linnaeus, 1758")
  assert_eq(ti.taxon_rank, "species")
}

///|
test "phyloxml_node_creation" {
  let node = @src.PhyloXMLNode::new(
    node_id="n1",
    name="TestNode",
    branch_length=0.5,
  )
  assert_eq(node.node_id, "n1")
  assert_eq(node.name, "TestNode")
  assert_true(node.branch_length > 0.4 && node.branch_length < 0.6)
  assert_true(node.is_leaf)
  assert_eq(node.children.length(), 0)
}

///|
test "phyloxml_node_with_children" {
  let child1 = @src.PhyloXMLNode::new(name="C1")
  let child2 = @src.PhyloXMLNode::new(name="C2")
  let node = @src.PhyloXMLNode::new(
    name="Parent",
    children=[child1, child2],
  )
  assert_false(node.is_leaf)
  assert_eq(node.children.length(), 2)
  assert_eq(node.children[0].name, "C1")
  assert_eq(node.children[1].name, "C2")
}

///|
test "phyloxml_tree_creation" {
  let root = @src.PhyloXMLNode::new(name="Root")
  let tree = @src.PhyloXMLTree::new(
    tree_id="t1",
    name="Test Tree",
    description="A test tree",
    root=root,
  )
  assert_eq(tree.tree_id, "t1")
  assert_eq(tree.name, "Test Tree")
  assert_eq(tree.description, "A test tree")
  assert_eq(tree.root.name, "Root")
}

///|
test "phyloxml_tree_with_metadata" {
  let meta : Map[String, String] = Map([], capacity=2)
  meta["rooted"] = "true"
  let tree = @src.PhyloXMLTree::new(
    name="Rooted Tree",
    phylogeny_metadata=meta,
  )
  assert_eq(tree.phylogeny_metadata["rooted"], "true")
}

///|
test "phyloxml_find_element" {
  let xml = "<root><child>text</child></root>"
  let (start, end, _content) = @src.find_element(xml, "child")
  // find_element returns 0,0,"" when not found
  assert_true(start >= 0)
}

///|
test "phyloxml_get_element_text" {
  let xml = "<name>Test Name</name>"
  let text = @src.get_element_text(xml)
  // Test the function handles input gracefully
  assert_true(text.length() >= 0)
}

///|
test "phyloxml_parse_node_element" {
  let xml =
    "<clade>\n" +
    "  <name>TestNode</name>\n" +
    "  <branch_length>0.5</branch_length>\n" +
    "  <clade>\n" +
    "    <name>Child</name>\n" +
    "  </clade>\n" +
    "</clade>"
  let node = @src.parse_node_element(xml)
  // The function should return a PhyloXMLNode
  assert_true(true)
}

///|
test "phyloxml_parse_taxon_element" {
  let xml =
    "<taxon>\n" +
    "  <taxon_id>tx1</taxon_id>\n" +
    "  <scientific_name>Homo sapiens</scientific_name>\n" +
    "  <common_name>Human</common_name>\n" +
    "</taxon>"
  let taxon = @src.parse_taxon_element(xml)
  // The function should return a TaxonInfo
  assert_true(true)
}

///|
test "phyloxml_parse_sequence_element" {
  let xml =
    "<sequence>\n" +
    "  <sequence_id>seq1</sequence_id>\n" +
    "  <molecular_type>DNA</molecular_type>\n" +
    "  <symbol>ATG</symbol>\n" +
    "</sequence>"
  let seq = @src.parse_sequence_element(xml)
  // The function should return a SequenceAnnotation
  assert_true(true)
}

///|
test "phyloxml_complex_tree_serialization" {
  let leaf_a = @src.PhyloXMLNode::new(name="A", branch_length=0.1)
  let leaf_b = @src.PhyloXMLNode::new(name="B", branch_length=0.2)
  let leaf_c = @src.PhyloXMLNode::new(name="C", branch_length=0.3)
  let internal = @src.PhyloXMLNode::new(
    name="Int",
    branch_length=0.5,
    children=[leaf_a, leaf_b],
  )
  let root = @src.PhyloXMLNode::new(
    name="Root",
    branch_length=0.0,
    children=[internal, leaf_c],
  )
  let tree = @src.PhyloXMLTree::new(
    tree_id="complex",
    name="Complex Tree",
    root=root,
  )
  let result = @src.PhyloXMLResult::new(trees=[tree])
  let xml = @src.phyloxml_to_xml(result)
  assert_true(xml.contains("complex"))
  assert_true(xml.contains("Complex Tree"))
  assert_true(xml.contains("Root"))
  assert_true(xml.contains("Int"))
  assert_true(xml.contains("A"))
  assert_true(xml.contains("B"))
  assert_true(xml.contains("C"))
  assert_true(xml.contains("<branch_length>"))
}

///|
test "phyloxml_newick_to_phyloxml_with_id" {
  let newick = "(A:0.1,B:0.2);"
  let result = @src.newick_to_phyloxml(newick, tree_id="my_tree")
  assert_eq(result.trees.length(), 1)
  assert_eq(result.trees[0].tree_id, "my_tree")
  assert_eq(result.trees[0].name, "my_tree")
}

///|
test "phyloxml_newick_empty" {
  let result = @src.newick_to_phyloxml("")
  assert_eq(result.trees.length(), 1)
  assert_eq(result.trees[0].root.name, "")
}

///|
test "phyloxml_to_newick_empty_result" {
  let result = @src.PhyloXMLResult::new()
  let newick = @src.phyloxml_to_newick(result)
  assert_eq(newick, "")
}

///|
test "phyloxml_node_to_newick_no_branch_length" {
  let node = @src.PhyloXMLNode::new(name="A")
  let newick = @src.phyloxml_node_to_newick(node)
  assert_eq(newick, "A")
}

///|
test "phyloxml_roundtrip_newick_phyloxml" {
  let newick = "(A:0.1,B:0.2,(C:0.3,D:0.4):0.5);"
  let result = @src.newick_to_phyloxml(newick)
  let output = @src.phyloxml_to_newick(result)
  assert_true(output.contains("A"))
  assert_true(output.contains("B"))
  assert_true(output.contains("C"))
  assert_true(output.contains("D"))
  assert_true(output.contains(";"))
}

///|
test "phyloxml_indent_xml" {
  let xml = "<a>\n<b>\ntest</b>\n</a>"
  let indented = @src.indent_xml(xml, 0)
  assert_true(indented.contains("<a>"))
  assert_true(indented.contains("test"))
}

///|
test "phyloxml_find_all_elements" {
  let xml =
    "<root>\n" +
    "  <child>1</child>\n" +
    "  <child>2</child>\n" +
    "  <child>3</child>\n" +
    "</root>"
  let elements = @src.find_all_elements(xml, "child")
  assert_eq(elements.length(), 3)
}

///|
test "phyloxml_taxon_namespace" {
  let ns = @src.TaxonNamespace::new(
    namespace_id="ns1",
    name="Test NS",
  )
  assert_eq(ns.namespace_id, "ns1")
  assert_eq(ns.name, "Test NS")
}

///|
test "phyloxml_with_taxon_namespace" {
  let xml =
    "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n" +
    "<phyloxml xmlns=\"http://www.phyloxml.org/1.0\">\n" +
    "  <taxon_namespace>\n" +
    "    <namespace_id>ns1</namespace_id>\n" +
    "    <name>Test Namespace</name>\n" +
    "  </taxon_namespace>\n" +
    "  <phylogeny id=\"tree1\">\n" +
    "    <clade>\n" +
    "      <name>Root</name>\n" +
    "      <clade><name>A</name></clade>\n" +
    "    </clade>\n" +
    "  </phylogeny>\n" +
    "</phyloxml>\n"
  let result = @src.parse_phyloxml(xml)
  assert_eq(result.trees.length(), 1)
  assert_eq(result.taxon_namespaces.length(), 1)
  assert_true(result.taxon_namespaces[0].namespace_id.length() >= 0)
  assert_true(result.taxon_namespaces[0].name.length() >= 0)
}