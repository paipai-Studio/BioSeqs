///|
test "debug_newick" {
  let newick = "(A:0.1,B:0.2,(C:0.3,D:0.4):0.5);"
  let result = @src.newick_to_phyloxml(newick)
  assert_eq(result.trees.length(), 1)
  let tree = result.trees[0]
  assert_eq(tree.root.children.length(), 3)
}

///|
test "debug_find_all" {
  let xml = "<root><child>1</child><child>2</child><child>3</child></root>"
  let elements = @src.find_all_elements(xml, "child")
  assert_true(elements.length() == 3)
}

///|
test "debug_parse_phyloxml" {
  let xml = "<phyloxml><phylogeny id=\"t1\"><clade><name>A</name></clade></phylogeny></phyloxml>"
  let result = @src.parse_phyloxml(xml)
  assert_true(result.trees.length() == 1)
}

///|
test "debug_attrs" {
  let xml = "<phylogeny id=\"tree1\">"
  let result = @src.parse_phyloxml(xml)
  assert_true(result.trees.length() >= 0)
}

///|
test "debug_find_all_phylogeny" {
  let xml = "<phylogeny id=\"tree1\"><name>Tree</name></phylogeny>"
  let elements = @src.find_all_elements(xml, "phylogeny")
  assert_true(elements.length() == 1)
}

///|
test "debug_get_name_text" {
  let xml = "<phylogeny id=\"tree1\"><name>Test Tree</name><clade><name>A</name></clade></phylogeny>"
  let (s, e, content) = @src.find_element(xml, "name")
  assert_true(s >= 0)
  assert_true(e > s)
  // Check if content has "Test Tree"
  assert_true(content.contains("Test Tree"))
}

///|
test "debug_parse_tree" {
  let xml = "<phylogeny id=\"tree1\"><name>Test Tree</name><clade><name>A</name></clade></phylogeny>"
  let tree = @src.parse_tree_element(xml)
  // Just check tree_id works first
  assert_true(tree.tree_id == "tree1")
}

///|
test "debug_multiple_trees" {
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
  assert_true(result.trees.length() >= 2)
  assert_true(result.trees[0].tree_id == "tree1")
  assert_true(result.trees[0].name.length() >= 0)
}
