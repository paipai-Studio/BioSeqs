///|
/// CDAO (Comparative Data Analysis Ontology) format demo.
///
/// Demonstrates the CDAO RDF/XML format for phylogenetic data:
///   - CDAO document structure (trees, nodes, TUs)
///   - Parsing CDAO RDF/XML
///   - Writing trees to CDAO format
///   - Converting between CDAO and BioSeqs Tree
///
/// Inspired by Biopython's Bio.Phylo.CDAO and Bio.Phylo.CDAOIO.
/// Run: moon run examples/phylo_cdao_demo
fn main {
  println("=== CDAO Module Demo ===")
  println("Comparative Data Analysis Ontology format")
  println("Reference: Biopython Bio.Phylo.CDAO / CDAOIO")
  println("")

  // ------------------------------------------------------------------
  // 1. Namespace constants
  // ------------------------------------------------------------------
  println("=== 1. CDAO Namespace Constants ===")
  println("CDAO namespace: " + @bio.cdao_namespace())
  println("RDF namespace:  " + @bio.cdao_rdf_namespace())
  println("RDFS namespace: " + @bio.cdao_rdfs_namespace())
  println("")

  // ------------------------------------------------------------------
  // 2. Sample CDAO document
  // ------------------------------------------------------------------
  println("=== 2. Sample CDAO Document ===")
  let doc = @bio.cdao_sample_document()
  println(@bio.cdao_document_summary(doc))
  println("")

  // Show node relationships
  println("Node relationships:")
  let node_keys = doc.nodes.keys()
  let node_arr : Array[String] = node_keys.collect()
  let mut i = 0
  while i < node_arr.length() {
    let node_id = node_arr[i]
    let node = doc.nodes[node_id]
    let parent_str = match node.parent_id {
      Some(p) => p
      None => "(root)"
    }
    let tu_str = match node.tu_id {
      Some(t) => " TU=" + t
      None => ""
    }
    let bl_str = match node.branch_length {
      Some(bl) => " bl=" + bl.to_string()
      None => ""
    }
    let child_count = node.children.length()
    println("  " + node_id + " parent=" + parent_str + " children=" + child_count.to_string() + tu_str + bl_str)
    i = i + 1
  }
  println("")

  // ------------------------------------------------------------------
  // 3. Parse CDAO XML
  // ------------------------------------------------------------------
  println("=== 3. Parse CDAO RDF/XML ===")
  let xml = @bio.cdao_sample_xml()
  println("Parsing sample CDAO XML...")
  let parsed_doc = @bio.cdao_parse(xml)
  println("Parsed document summary:")
  println(@bio.cdao_document_summary(parsed_doc))
  println("")

  // ------------------------------------------------------------------
  // 4. Convert CDAO to BioSeqs Tree
  // ------------------------------------------------------------------
  println("=== 4. Convert CDAO to Tree ===")
  let trees = @bio.cdao_to_trees(doc)
  println("Converted " + trees.length().to_string() + " tree(s) from CDAO document")
  if trees.length() > 0 {
    let tree = trees[0]
    println("Tree rooted: " + tree.rooted.to_string())
    println("Root children: " + tree.root.clades.length().to_string())
    // Print leaf names
    let terminals = tree.root.get_terminals()
    println("Leaf taxa (" + terminals.length().to_string() + "):")
    let mut t = 0
    while t < terminals.length() {
      let name = match terminals[t].name {
        Some(n) => n
        None => "(unnamed)"
      }
      println("  " + name)
      t = t + 1
    }
  }
  println("")

  // ------------------------------------------------------------------
  // 5. Write tree to CDAO format
  // ------------------------------------------------------------------
  println("=== 5. Write Tree to CDAO Format ===")
  let sample_tree = @bio.parsimony_sample_tree()
  println("Tree topology: ((A, B), (C, D))")
  println("")
  println("Generating CDAO RDF/XML...")
  let cdao_xml = @bio.cdao_write(sample_tree)
  println("Output XML (" + cdao_xml.length().to_string() + " chars):")
  println("---")
  // Print first 500 characters
  let display_len = if cdao_xml.length() > 500 { 500 } else { cdao_xml.length() }
  println(cdao_xml[0:display_len].to_owned())
  if cdao_xml.length() > 500 {
    println("... (truncated)")
  }
  println("---")
  println("")

  // ------------------------------------------------------------------
  // 6. Roundtrip: Write then parse
  // ------------------------------------------------------------------
  println("=== 6. Roundtrip: Write -> Parse ===")
  let roundtrip_doc = @bio.cdao_parse(cdao_xml)
  let roundtrip_trees = @bio.cdao_to_trees(roundtrip_doc)
  println("Roundtrip results:")
  println("  Trees: " + roundtrip_trees.length().to_string())
  if roundtrip_trees.length() > 0 {
    let rt_tree = roundtrip_trees[0]
    let rt_terminals = rt_tree.root.get_terminals()
    println("  Leaf taxa: " + rt_terminals.length().to_string())
    let mut rt = 0
    while rt < rt_terminals.length() {
      let name = match rt_terminals[rt].name {
        Some(n) => n
        None => "(unnamed)"
      }
      println("    " + name)
      rt = rt + 1
    }
  }
  println("")

  // ------------------------------------------------------------------
  // 7. Build custom CDAO document
  // ------------------------------------------------------------------
  println("=== 7. Build Custom CDAO Document ===")
  let custom_doc = @bio.CdaoDocument::new()
  custom_doc.trees.push(@bio.Cdaotree::new("#mytree", rooted=true, root_node_id="#n1"))
  // Simple tree: ((X, Y), Z)
  custom_doc.nodes["#n1"] = @bio.CdaoNode::new("#n1", children=["#n2", "#n3"])
  custom_doc.nodes["#n2"] = @bio.CdaoNode::new("#n2", children=["#n4", "#n5"], parent_id=Some("#n1"))
  custom_doc.nodes["#n3"] = @bio.CdaoNode::new("#n3", parent_id=Some("#n1"), tu_id=Some("#tu3"), branch_length=Some(0.5))
  custom_doc.nodes["#n4"] = @bio.CdaoNode::new("#n4", parent_id=Some("#n2"), tu_id=Some("#tu1"), branch_length=Some(1.0))
  custom_doc.nodes["#n5"] = @bio.CdaoNode::new("#n5", parent_id=Some("#n2"), tu_id=Some("#tu2"), branch_length=Some(1.0))
  custom_doc.tus["#tu1"] = @bio.CdaoTU::new("#tu1", label=Some("X"))
  custom_doc.tus["#tu2"] = @bio.CdaoTU::new("#tu2", label=Some("Y"))
  custom_doc.tus["#tu3"] = @bio.CdaoTU::new("#tu3", label=Some("Z"))
  println("Custom document summary:")
  println(@bio.cdao_document_summary(custom_doc))
  let custom_trees = @bio.cdao_to_trees(custom_doc)
  println("Converted to " + custom_trees.length().to_string() + " tree(s)")
  if custom_trees.length() > 0 {
    let ct = custom_trees[0]
    let ct_terminals = ct.root.get_terminals()
    println("Leaf taxa:")
    let mut ci = 0
    while ci < ct_terminals.length() {
      let name = match ct_terminals[ci].name {
        Some(n) => n
        None => "(unnamed)"
      }
      println("  " + name)
      ci = ci + 1
    }
  }
  println("")

  println("=== CDAO Demo Complete ===")
}
