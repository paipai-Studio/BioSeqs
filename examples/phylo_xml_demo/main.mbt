///|
/// PhyloXML Demo
/// Demonstrates Bio.Phylo.PhyloXML functionality:
///   - Parse PhyloXML text
///   - Display tree names, node structures, branch lengths
///   - Show taxon annotations and sequence annotations
///   - Demonstrate XML serialization (round-trip)
///   - Convert Newick to PhyloXML and back
///   - Get all tips, node counts
///
/// Run: moon run examples/phylo_xml_demo/

fn main {
  println("=== PhyloXML Demo ===")

  // =========================================================================
  // 1. Create sample PhyloXML text
  // =========================================================================
  println("")
  println("--- 1. Sample PhyloXML Input ---")

  let phyloxml_text =
    "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n" +
    "<phyloxml xmlns=\"http://www.phyloxml.org/1.0\">\n" +
    "  <phylogeny id=\"tree1\" rooted=\"true\">\n" +
    "    <name>Species Tree</name>\n" +
    "    <description>A phylogenetic tree of species</description>\n" +
    "    <clade>\n" +
    "      <name>root</name>\n" +
    "      <branch_length>0.0</branch_length>\n" +
    "      <taxon>\n" +
    "        <taxon_id>root_taxon</taxon_id>\n" +
    "        <scientific_name>Life</scientific_name>\n" +
    "        <common_name>Root</common_name>\n" +
    "      </taxon>\n" +
    "      <clade>\n" +
    "        <name>mammals</name>\n" +
    "        <branch_length>0.12</branch_length>\n" +
    "        <clade>\n" +
    "          <name>Human</name>\n" +
    "          <branch_length>0.05</branch_length>\n" +
    "          <taxon>\n" +
    "            <taxon_id>9606</taxon_id>\n" +
    "            <scientific_name>Homo sapiens</scientific_name>\n" +
    "            <common_name>Human</common_name>\n" +
    "            <taxon_rank>species</taxon_rank>\n" +
    "          </taxon>\n" +
    "          <sequence>\n" +
    "            <sequence_id>human_seq</sequence_id>\n" +
    "            <molecular_type>DNA</molecular_type>\n" +
    "            <symbol>Hs</symbol>\n" +
    "          </sequence>\n" +
    "          <annotation>\n" +
    "            <name>gc_content</name>\n" +
    "            <value>0.41</value>\n" +
    "          </annotation>\n" +
    "        </clade>\n" +
    "        <clade>\n" +
    "          <name>Mouse</name>\n" +
    "          <branch_length>0.08</branch_length>\n" +
    "          <taxon>\n" +
    "            <taxon_id>10090</taxon_id>\n" +
    "            <scientific_name>Mus musculus</scientific_name>\n" +
    "            <common_name>Mouse</common_name>\n" +
    "          </taxon>\n" +
    "        </clade>\n" +
    "      </clade>\n" +
    "      <clade>\n" +
    "        <name>yeast</name>\n" +
    "        <branch_length>0.30</branch_length>\n" +
    "        <taxon>\n" +
    "          <taxon_id>4932</taxon_id>\n" +
    "          <scientific_name>Saccharomyces cerevisiae</scientific_name>\n" +
    "          <common_name>Baker's yeast</common_name>\n" +
    "        </taxon>\n" +
    "      </clade>\n" +
    "    </clade>\n" +
    "  </phylogeny>\n" +
    "  <phylogeny id=\"tree2\">\n" +
    "    <name>Gene Tree</name>\n" +
    "    <clade>\n" +
    "      <name>GeneA</name>\n" +
    "      <branch_length>0.01</branch_length>\n" +
    "    </clade>\n" +
    "  </phylogeny>\n" +
    "</phyloxml>\n"

  println(phyloxml_text)

  // =========================================================================
  // 2. Parse using parse_phyloxml()
  // =========================================================================
  println("")
  println("--- 2. Parse PhyloXML ---")

  let result = @bio.parse_phyloxml(phyloxml_text)
  println("Parse complete!")
  println("Number of trees: \{result.trees.length()}")

  // =========================================================================
  // 3. Display tree names
  // =========================================================================
  println("")
  println("--- 3. Tree Names ---")

  let tree_names = @bio.get_tree_names(result)
  for i = 0; i < tree_names.length(); i = i + 1 {
    println("  Tree \{i + 1}: \{tree_names[i]}")
  }

  // =========================================================================
  // 4. Display tree structure and node details
  // =========================================================================
  println("")
  println("--- 4. Tree Structure Details ---")

  for tree_idx = 0; tree_idx < result.trees.length(); tree_idx = tree_idx + 1 {
    let tree = result.trees[tree_idx]
    println("")
    println("Tree \{tree_idx + 1}: '\{tree.name}'")
    println("  ID: '\{tree.tree_id}'")
    println("  Description: '\{tree.description}'")

    let rooted = tree.phylogeny_metadata.get("rooted").unwrap_or("false")
    println("  Rooted: \{rooted}")

    let node_count = @bio.phyloxml_get_node_count(tree)
    println("  Total nodes: \{node_count}")

    let tips = @bio.get_all_tips(tree)
    println("  Tips (\{tips.length()}):")
    for tip in tips {
      println("    - \{tip}")
    }

    println("  Root node name: '\{tree.root.name}'")
    println("  Root branch length: \{tree.root.branch_length}")
    println("  Root is leaf: \{tree.root.is_leaf}")
    println("  Root children count: \{tree.root.children.length()}")

    display_node(tree.root, "  ")
  }

  // =========================================================================
  // 5. Show taxon and sequence annotations
  // =========================================================================
  println("")
  println("--- 5. Taxon & Sequence Annotations ---")

  let tree1 = result.trees[0]
  println("Examining '\{tree1.name}':")
  display_annotations(tree1.root, "")

  // =========================================================================
  // 6. XML serialization (round-trip)
  // =========================================================================
  println("")
  println("--- 6. XML Serialization (Round-Trip) ---")

  let serialized = @bio.phyloxml_to_xml(result)
  println("Serialized XML:")
  println(serialized)

  let reparsed = @bio.parse_phyloxml(serialized)
  println("Re-parsed tree count: \{reparsed.trees.length()}")
  println("Re-parsed tree names: \{@bio.get_tree_names(reparsed).join(", ")}")

  // =========================================================================
  // 7. Convert Newick to PhyloXML and back
  // =========================================================================
  println("")
  println("--- 7. Newick <-> PhyloXML Conversion ---")

  let newick = "((Human:0.05,Mouse:0.08)mammals:0.12,Yeast:0.30)root:0.0;"
  println("Input Newick: \{newick}")

  let newick_result = @bio.newick_to_phyloxml(newick, tree_id="newick_tree")
  println("")
  println("Converted to PhyloXML:")
  println("  Trees: \{newick_result.trees.length()}")
  println("  Tree name: '\{newick_result.trees[0].name}'")
  println("  Tree ID: '\{newick_result.trees[0].tree_id}'")

  let newick_xml = @bio.phyloxml_to_xml(newick_result)
  println("")
  println("Generated PhyloXML:")
  println(newick_xml)

  let back_to_newick = @bio.phyloxml_to_newick(newick_result)
  println("")
  println("Converted back to Newick: \{back_to_newick}")

  // =========================================================================
  // 8. Get all tips and node counts
  // =========================================================================
  println("")
  println("--- 8. Tree Utility Functions ---")

  for tree_idx = 0; tree_idx < result.trees.length(); tree_idx = tree_idx + 1 {
    let tree = result.trees[tree_idx]
    let tips = @bio.get_all_tips(tree)
    let node_count = @bio.phyloxml_get_node_count(tree)
    println("Tree '\{tree.name}':")
    println("  Tips: \{tips.length()}")
    println("  Nodes: \{node_count}")
  }

  // =========================================================================
  // 9. Construct tree programmatically with annotations
  // =========================================================================
  println("")
  println("--- 9. Programmatically Constructed Tree ---")

  let human_node = @bio.PhyloXMLNode::new(
    name="Human",
    branch_length=0.05,
    taxon=Some(@bio.TaxonInfo::new(
      taxon_id="9606",
      scientific_name="Homo sapiens",
      common_name="Human",
    )),
    sequence_annotation=Some(@bio.SequenceAnnotation::new(
      sequence_id="human_seq",
      molecular_type="DNA",
      symbol="Hs",
    )),
  )

  let chimp_node = @bio.PhyloXMLNode::new(
    name="Chimpanzee",
    branch_length=0.06,
    taxon=Some(@bio.TaxonInfo::new(
      taxon_id="9598",
      scientific_name="Pan troglodytes",
      common_name="Chimpanzee",
    )),
  )

  let mammals_node = @bio.PhyloXMLNode::new(
    name="Mammalia",
    branch_length=0.10,
    children=[human_node, chimp_node],
    annotations=Map([("group", "mammals"), ("description", "Primate clade")], capacity=2),
  )

  let tree_root = @bio.PhyloXMLNode::new(
    name="root",
    branch_length=0.0,
    children=[mammals_node],
  )

  let constructed_tree = @bio.PhyloXMLTree::new(
    tree_id="constructed",
    name="Constructed Tree",
    description="Tree built programmatically",
    root=tree_root,
  )

  let constructed_result = @bio.PhyloXMLResult::new(
    trees=[constructed_tree],
  )

  println("Constructed tree XML:")
  let constructed_xml = @bio.phyloxml_to_xml(constructed_result)
  println(constructed_xml)

  // =========================================================================
  // 10. Edge cases
  // =========================================================================
  println("")
  println("--- 10. Edge Cases ---")

  let empty = @bio.parse_phyloxml("")
  println("Empty text: \{empty.trees.length()} trees")

  let simple = @bio.parse_phyloxml(
    "<?xml version=\"1.0\"?>\n<phyloxml>\n  <phylogeny><name>Empty</name></phylogeny>\n</phyloxml>",
  )
  println("Tree with no clades: '\{simple.trees[0].name}'")
  println("  Root name: '\{simple.trees[0].root.name}'")
  println("  Root is leaf: \{simple.trees[0].root.is_leaf}")

  println("")
  println("=== Demo complete ===")
}

///|
/// Recursively display node details including children.
fn display_node(node : @bio.PhyloXMLNode, indent : String) -> Unit {
  println("\{indent}Node: '\{node.name}'")
  println("\{indent}  branch_length: \{node.branch_length}")
  println("\{indent}  is_leaf: \{node.is_leaf}")
  println("\{indent}  children: \{node.children.length()}")

  match node.taxon {
    Some(t) => {
      println("\{indent}  taxon:")
      println("\{indent}    ID: '\{t.taxon_id}'")
      println("\{indent}    Scientific: '\{t.scientific_name}'")
      println("\{indent}    Common: '\{t.common_name}'")
      println("\{indent}    Rank: '\{t.taxon_rank}'")
    }
    None => ()
  }

  match node.sequence_annotation {
    Some(s) => {
      println("\{indent}  sequence:")
      println("\{indent}    ID: '\{s.sequence_id}'")
      println("\{indent}    Type: '\{s.molecular_type}'")
      println("\{indent}    Symbol: '\{s.symbol}'")
    }
    None => ()
  }

  let anns = node.annotations
  if anns.length() > 0 {
    println("\{indent}  annotations (\{anns.length()}):")
    for key in anns.keys() {
      println("\{indent}    \{key} = '\{anns[key]}'")
    }
  }

  for child in node.children {
    display_node(child, indent + "  ")
  }
}

///|
/// Recursively display annotations on nodes.
fn display_annotations(node : @bio.PhyloXMLNode, indent : String) -> Unit {
  if node.name.length() > 0 {
    println("\{indent}'\{node.name}':")
  }
  match node.taxon {
    Some(t) => {
      if t.scientific_name.length() > 0 {
        println("\{indent}  Taxon: \{t.scientific_name} (\{t.common_name})")
      }
    }
    None => ()
  }
  match node.sequence_annotation {
    Some(s) => {
      println("\{indent}  Sequence: \{s.sequence_id} [\{s.molecular_type}]")
    }
    None => ()
  }
  let anns = node.annotations
  for key in anns.keys() {
    println("\{indent}  Annotation: \{key} = '\{anns[key]}'")
  }
  for child in node.children {
    display_annotations(child, indent + "  ")
  }
}