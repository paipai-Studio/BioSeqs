///|
/// Bio.Phylo.NeXML - NeXML format parsing and serialization.
/// NeXML is an XML-based format for phylogenetic data, including
/// trees, character matrices, and OTU definitions.

// ============================================================================
// Data Structures
// ============================================================================

///|
/// NeXML state (for discrete characters).
pub struct NeXMLState {
  id : String
  symbol : String
} derive(Eq)

///|
pub extend NeXMLState with Eq::{equal, not_equal}

///|
/// NeXML character definition.
pub struct NeXMLCharacter {
  id : String
  label : String
  states : Array[String]
} derive(Eq)

///|
pub extend NeXMLCharacter with Eq::{equal, not_equal}

///|
/// NeXML characters matrix.
pub struct NeXMLCharacters {
  id : String
  otus_id : String
  matrix_type : String
  characters : Array[NeXMLCharacter]
  matrix : Map[String, Array[String]]
} derive(Eq)

///|
pub extend NeXMLCharacters with Eq::{equal, not_equal}

///|
/// NeXML OTU (Operational Taxonomic Unit).
pub struct NeXMLOTU {
  id : String
  label : String
} derive(Eq)

///|
pub extend NeXMLOTU with Eq::{equal, not_equal}

///|
/// NeXML OTUs block.
pub struct NeXMLOTUs {
  id : String
  otus : Array[NeXMLOTU]
} derive(Eq)

///|
pub extend NeXMLOTUs with Eq::{equal, not_equal}

///|
/// NeXML tree node.
pub struct NeXMLNode {
  id : String
  label : String
  otu : String
  length : Double?
} derive(Eq)

///|
pub extend NeXMLNode with Eq::{equal, not_equal}

///|
/// NeXML tree edge.
pub struct NeXMLEdge {
  id : String
  source : String
  target : String
  length : Double?
} derive(Eq)

///|
pub extend NeXMLEdge with Eq::{equal, not_equal}

///|
/// NeXML tree (edges + nodes).
pub struct NeXMLTree {
  id : String
  name : String
  tree_type : String
  rooted : Bool
  nodes : Array[NeXMLNode]
  edges : Array[NeXMLEdge]
} derive(Eq)

///|
pub extend NeXMLTree with Eq::{equal, not_equal}

///|
/// NeXML trees block.
pub struct NeXMLTrees {
  id : String
  otus_id : String
  trees : Array[NeXMLTree]
} derive(Eq)

///|
pub extend NeXMLTrees with Eq::{equal, not_equal}

///|
/// NeXML document root.
pub struct NeXMLDocument {
  version : String
  generator : String
  otus_blocks : Array[NeXMLOTUs]
  characters_blocks : Array[NeXMLCharacters]
  trees_blocks : Array[NeXMLTrees]
} derive(Eq)

///|
pub extend NeXMLDocument with Eq::{equal, not_equal}

// ============================================================================
// Constructors
// ============================================================================

///|
pub fn NeXMLState::new(id : String, symbol : String) -> NeXMLState {
  NeXMLState::{ id, symbol, }
}

///|
pub fn NeXMLCharacter::new(
  id : String,
  label : String,
  states : Array[String],
) -> NeXMLCharacter {
  NeXMLCharacter::{ id, label, states, }
}

///|
pub fn NeXMLCharacters::new(
  id : String,
  otus_id : String,
  matrix_type : String,
  characters : Array[NeXMLCharacter],
  matrix : Map[String, Array[String]],
) -> NeXMLCharacters {
  NeXMLCharacters::{ id, otus_id, matrix_type, characters, matrix, }
}

///|
pub fn NeXMLOTU::new(id : String, label : String) -> NeXMLOTU {
  NeXMLOTU::{ id, label, }
}

///|
pub fn NeXMLOTUs::new(id : String, otus : Array[NeXMLOTU]) -> NeXMLOTUs {
  NeXMLOTUs::{ id, otus, }
}

///|
pub fn NeXMLNode::new(
  id : String,
  label : String,
  otu : String,
  length~ : Double? = None,
) -> NeXMLNode {
  NeXMLNode::{ id, label, otu, length, }
}

///|
pub fn NeXMLEdge::new(
  id : String,
  source : String,
  target : String,
  length~ : Double? = None,
) -> NeXMLEdge {
  NeXMLEdge::{ id, source, target, length, }
}

///|
pub fn NeXMLTree::new(
  id : String,
  name : String,
  tree_type : String,
  rooted : Bool,
  nodes : Array[NeXMLNode],
  edges : Array[NeXMLEdge],
) -> NeXMLTree {
  NeXMLTree::{ id, name, tree_type, rooted, nodes, edges, }
}

///|
pub fn NeXMLTrees::new(
  id : String,
  otus_id : String,
  trees : Array[NeXMLTree],
) -> NeXMLTrees {
  NeXMLTrees::{ id, otus_id, trees, }
}

///|
pub fn NeXMLDocument::new(
  version~ : String = "0.9",
  generator~ : String = "",
) -> NeXMLDocument {
  NeXMLDocument::{
    version,
    generator,
    otus_blocks: [],
    characters_blocks: [],
    trees_blocks: [],
  }
}

// ============================================================================
// Helper functions
// ============================================================================

// ============================================================================
// NeXML Parsing
// ============================================================================

// ============================================================================
// NeXML Serialization
// ============================================================================

///|
/// Serialize a NeXML document to XML string.
pub fn nexml_to_xml(doc : NeXMLDocument) -> String {
  let mut xml = "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n"
  xml = xml + "<nex:nexml\n"
  xml = xml + "  version=\"" + doc.version + "\"\n"
  xml = xml + "  generator=\"" + doc.generator + "\"\n"
  xml = xml + "  xmlns:nex=\"http://www.nexml.org/2009\"\n"
  xml = xml + "  xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\"\n"
  xml = xml + "  xmlns:xml=\"http://www.w3.org/XML/1998/namespace\">\n"

  // Write OTUs blocks
  for i = 0; i < doc.otus_blocks.length(); i = i + 1 {
    let otus_block = doc.otus_blocks[i]
    xml = xml + "  <otus id=\"" + otus_block.id + "\">\n"
    for j = 0; j < otus_block.otus.length(); j = j + 1 {
      let otu = otus_block.otus[j]
      xml = xml +
        "    <otu id=\"" +
        otu.id +
        "\" label=\"" +
        otu.label +
        "\"/>\n"
    }
    xml = xml + "  </otus>\n"
  }

  // Write characters blocks
  for i = 0; i < doc.characters_blocks.length(); i = i + 1 {
    let chars_block = doc.characters_blocks[i]
    xml = xml +
      "  <characters id=\"" +
      chars_block.id +
      "\" otus=\"" +
      chars_block.otus_id +
      "\">\n"
    for j = 0; j < chars_block.characters.length(); j = j + 1 {
      let ch = chars_block.characters[j]
      xml = xml +
        "    <char id=\"" +
        ch.id +
        "\" label=\"" +
        ch.label +
        "\"/>\n"
    }
    xml = xml + "    <matrix>\n"
    for otu in chars_block.matrix.keys() {
      let row_opt = chars_block.matrix.get(otu)
      if row_opt is Some(_) {
        let row = row_opt.unwrap()
        xml = xml + "      <row otu=\"" + otu + "\">\n"
        for k = 0; k < row.length(); k = k + 1 {
          xml = xml + "        <cell state=\"" + row[k] + "\"/>\n"
        }
        xml = xml + "      </row>\n"
      }
    }
    xml = xml + "    </matrix>\n"
    xml = xml + "  </characters>\n"
  }

  // Write trees blocks
  for i = 0; i < doc.trees_blocks.length(); i = i + 1 {
    let trees_block = doc.trees_blocks[i]
    xml = xml +
      "  <trees id=\"" +
      trees_block.id +
      "\" otus=\"" +
      trees_block.otus_id +
      "\">\n"
    for j = 0; j < trees_block.trees.length(); j = j + 1 {
      let tree = trees_block.trees[j]
      let rooted_str = if tree.rooted { "true" } else { "false" }
      xml = xml +
        "    <tree id=\"" +
        tree.id +
        "\" label=\"" +
        tree.name +
        "\" rooted=\"" +
        rooted_str +
        "\">\n"
      for k = 0; k < tree.nodes.length(); k = k + 1 {
        let node = tree.nodes[k]
        xml = xml + "      <node id=\"" + node.id + "\""
        if node.label.length() > 0 {
          xml = xml + " label=\"" + node.label + "\""
        }
        if node.otu.length() > 0 {
          xml = xml + " otu=\"" + node.otu + "\""
        }
        match node.length {
          Some(len) => xml = xml + " length=\"" + len.to_string() + "\""
          None => ()
        }
        xml = xml + "/>\n"
      }
      for k = 0; k < tree.edges.length(); k = k + 1 {
        let edge = tree.edges[k]
        xml = xml +
          "      <edge id=\"" +
          edge.id +
          "\" source=\"" +
          edge.source +
          "\" target=\"" +
          edge.target +
          "\""
        match edge.length {
          Some(len) => xml = xml + " length=\"" + len.to_string() + "\""
          None => ()
        }
        xml = xml + "/>\n"
      }
      xml = xml + "    </tree>\n"
    }
    xml = xml + "  </trees>\n"
  }

  xml = xml + "</nex:nexml>\n"

  xml
}

// ============================================================================
// Conversion to Newick
// ============================================================================

///|
/// Convert a NeXML tree to Newick format string.
pub fn nexml_tree_to_newick(tree : NeXMLTree) -> String {
  if tree.nodes.length() == 0 {
    return ";"
  }

  // Build adjacency list
  let children_map : Map[String, Array[(String, Double?)]] = Map(
    [],
    capacity=16,
  )
  let parent_map : Map[String, String] = Map([], capacity=16)

  for i = 0; i < tree.edges.length(); i = i + 1 {
    let edge = tree.edges[i]
    let existing_opt = children_map.get(edge.source)
    let existing = match existing_opt {
      Some(arr) => arr
      None => []
    }
    let new_children = existing.copy()
    new_children.push((edge.target, edge.length))
    children_map.set(edge.source, new_children)
    parent_map.set(edge.target, edge.source)
  }

  // Find root (node with no parent)
  let mut root_id = ""
  for i = 0; i < tree.nodes.length(); i = i + 1 {
    let node = tree.nodes[i]
    if parent_map.get(node.id) is None {
      root_id = node.id
      break
    }
  }

  if root_id.length() == 0 && tree.nodes.length() > 0 {
    root_id = tree.nodes[0].id
  }

  // Build Newick recursively
  let (newick_str, _) = build_newick(tree, root_id, children_map)
  newick_str + ";"
}

///|
fn build_newick(
  tree : NeXMLTree,
  node_id : String,
  children_map : Map[String, Array[(String, Double?)]],
) -> (String, Bool) {
  let children_opt = children_map.get(node_id)

  // Find node label
  let mut node_label = ""
  let mut _node_length : Double? = None
  for i = 0; i < tree.nodes.length(); i = i + 1 {
    if tree.nodes[i].id == node_id {
      node_label = tree.nodes[i].label
      _node_length = tree.nodes[i].length
      break
    }
  }

  match children_opt {
    Some(children) =>
      if children.length() == 0 {
        // Leaf node
        (node_label, true)
      } else {
        // Internal node
        let mut result = "("
        for i = 0; i < children.length(); i = i + 1 {
          let (child_id, edge_length) = children[i]
          let (child_str, _) = build_newick(tree, child_id, children_map)
          if i > 0 {
            result = result + ","
          }
          match edge_length {
            Some(len) => result = result + child_str + ":" + len.to_string()
            None => result = result + child_str
          }
        }
        result = result + ")" + node_label
        (result, false)
      }
    None =>
      // Leaf node (no children in map)
      (node_label, true)
  }
}

// ============================================================================
// Utility functions
// ============================================================================

///|
/// Get all OTU labels from a NeXML document.
pub fn nexml_get_all_otu_labels(doc : NeXMLDocument) -> Array[String] {
  let labels : Array[String] = Array(capacity=0)
  for i = 0; i < doc.otus_blocks.length(); i = i + 1 {
    let otus_block = doc.otus_blocks[i]
    for j = 0; j < otus_block.otus.length(); j = j + 1 {
      labels.push(otus_block.otus[j].label)
    }
  }
  labels
}

///|
/// Get all tree names from a NeXML document.
pub fn nexml_get_tree_names(doc : NeXMLDocument) -> Array[String] {
  let names : Array[String] = Array(capacity=0)
  for i = 0; i < doc.trees_blocks.length(); i = i + 1 {
    let trees_block = doc.trees_blocks[i]
    for j = 0; j < trees_block.trees.length(); j = j + 1 {
      names.push(trees_block.trees[j].name)
    }
  }
  names
}

///|
/// Get the number of trees in a NeXML document.
pub fn nexml_get_tree_count(doc : NeXMLDocument) -> Int32 {
  let mut count = 0
  for i = 0; i < doc.trees_blocks.length(); i = i + 1 {
    count = count + doc.trees_blocks[i].trees.length()
  }
  count
}

///|
/// Get the number of OTUs in a NeXML document.
pub fn nexml_get_otu_count(doc : NeXMLDocument) -> Int32 {
  let mut count = 0
  for i = 0; i < doc.otus_blocks.length(); i = i + 1 {
    count = count + doc.otus_blocks[i].otus.length()
  }
  count
}
