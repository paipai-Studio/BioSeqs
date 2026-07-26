///|
/// Phylo NeXML Demo - NeXML format creation and serialization.
///
/// Run: moon run examples/phylo_nexml_demo/main.mbt
fn main {
  println("=== Phylo NeXML Demo ===\n")

  // 1. Create OTUs
  println("1. Creating OTUs block:")
  let otu1 = @bio.NeXMLOTU::new("otu1", "Homo_sapiens")
  let otu2 = @bio.NeXMLOTU::new("otu2", "Pan_troglodytes")
  let otu3 = @bio.NeXMLOTU::new("otu3", "Gorilla_gorilla")
  let otus = [otu1, otu2, otu3]
  let otus_block = @bio.NeXMLOTUs::new("otus1", otus)
  println("   OTUs block ID: \{otus_block.id}")
  println("   Number of OTUs: \{otus_block.otus.length()}")

  // 2. Create a tree
  println("\n2. Creating a phylogenetic tree:")
  let node1 = @bio.NeXMLNode::new("n1", "root", "")
  let node2 = @bio.NeXMLNode::new("n2", "Human", "otu1")
  let node3 = @bio.NeXMLNode::new("n3", "Chimp", "otu2")
  let node4 = @bio.NeXMLNode::new("n4", "Gorilla", "otu3")
  let nodes = [node1, node2, node3, node4]
  let edge1 = @bio.NeXMLEdge::new("e1", "n1", "n2", length=Some(0.1))
  let edge2 = @bio.NeXMLEdge::new("e2", "n1", "n3", length=Some(0.15))
  let edge3 = @bio.NeXMLEdge::new("e3", "n1", "n4", length=Some(0.2))
  let edges = [edge1, edge2, edge3]
  let tree = @bio.NeXMLTree::new(
    "tree1",
    "Primates Tree",
    "FloatTree",
    true,
    nodes,
    edges,
  )
  println("   Tree ID: \{tree.id}")
  println("   Tree name: \{tree.name}")
  println("   Number of nodes: \{tree.nodes.length()}")
  println("   Number of edges: \{tree.edges.length()}")

  // 3. Create trees block
  let trees = [tree]
  let trees_block = @bio.NeXMLTrees::new("trees1", "otus1", trees)
  println("   Trees block ID: \{trees_block.id}")

  // 4. Create full NeXML document
  println("\n3. Creating NeXML document:")
  let doc = @bio.NeXMLDocument::new(
    version="0.9",
    generator="BioSeqs",
  )
  println("   Version: \{doc.version}")
  println("   Generator: \{doc.generator}")
  println("   OTU blocks: \{doc.otus_blocks.length()}")
  println("   Tree blocks: \{doc.trees_blocks.length()}")

  // 5. Convert tree to Newick format
  println("\n4. Converting tree to Newick format:")
  let newick = @bio.nexml_tree_to_newick(tree)
  println("   Newick: \{newick}")

  // 6. Serialize to XML
  println("\n5. Serializing to NeXML format:")
  let xml = @bio.nexml_to_xml(doc)
  println("   XML length: \{xml.length()} characters")
  println("   First 200 chars:")
  let preview = if xml.length() > 200 { xml[0:200].to_owned() + "..." } else { xml }
  println("   \{preview}")

  // 7. Utility functions
  println("\n6. Utility functions:")
  println("   Total OTU count: \{@bio.nexml_get_otu_count(doc)}")
  println("   Total tree count: \{@bio.nexml_get_tree_count(doc)}")
  let labels = @bio.nexml_get_all_otu_labels(doc)
  println("   All OTU labels: \{labels.join(", ")}")

  println("\n=== Demo Complete ===")
}
