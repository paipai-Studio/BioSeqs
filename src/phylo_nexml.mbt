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
/// NeXML character definition.
pub struct NeXMLCharacter {
  id : String
  label : String
  states : Array[String]
} derive(Eq)

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
/// NeXML OTU (Operational Taxonomic Unit).
pub struct NeXMLOTU {
  id : String
  label : String
} derive(Eq)

///|
/// NeXML OTUs block.
pub struct NeXMLOTUs {
  id : String
  otus : Array[NeXMLOTU]
} derive(Eq)

///|
/// NeXML tree node.
pub struct NeXMLNode {
  id : String
  label : String
  otu : String
  length : Double?
} derive(Eq)

///|
/// NeXML tree edge.
pub struct NeXMLEdge {
  id : String
  source : String
  target : String
  length : Double?
} derive(Eq)

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
/// NeXML trees block.
pub struct NeXMLTrees {
  id : String
  otus_id : String
  trees : Array[NeXMLTree]
} derive(Eq)

///|
/// NeXML document root.
pub struct NeXMLDocument {
  version : String
  generator : String
  otus_blocks : Array[NeXMLOTUs]
  characters_blocks : Array[NeXMLCharacters]
  trees_blocks : Array[NeXMLTrees]
} derive(Eq)

// ============================================================================
// Constructors
// ============================================================================

///|
pub fn NeXMLState::new(id : String, symbol : String) -> NeXMLState {
  NeXMLState::{ id, symbol }
}

///|
pub fn NeXMLCharacter::new(
  id : String,
  label : String,
  states : Array[String],
) -> NeXMLCharacter {
  NeXMLCharacter::{ id, label, states }
}

///|
pub fn NeXMLCharacters::new(
  id : String,
  otus_id : String,
  matrix_type : String,
  characters : Array[NeXMLCharacter],
  matrix : Map[String, Array[String]],
) -> NeXMLCharacters {
  NeXMLCharacters::{ id, otus_id, matrix_type, characters, matrix }
}

///|
pub fn NeXMLOTU::new(id : String, label : String) -> NeXMLOTU {
  NeXMLOTU::{ id, label }
}

///|
pub fn NeXMLOTUs::new(id : String, otus : Array[NeXMLOTU]) -> NeXMLOTUs {
  NeXMLOTUs::{ id, otus }
}

///|
pub fn NeXMLNode::new(
  id : String,
  label : String,
  otu : String,
  length? : Double? = None,
) -> NeXMLNode {
  NeXMLNode::{ id, label, otu, length }
}

///|
pub fn NeXMLEdge::new(
  id : String,
  source : String,
  target : String,
  length? : Double? = None,
) -> NeXMLEdge {
  NeXMLEdge::{ id, source, target, length }
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
  NeXMLTree::{ id, name, tree_type, rooted, nodes, edges }
}

///|
pub fn NeXMLTrees::new(
  id : String,
  otus_id : String,
  trees : Array[NeXMLTree],
) -> NeXMLTrees {
  NeXMLTrees::{ id, otus_id, trees }
}

///|
pub fn NeXMLDocument::new(
  version? : String = "0.9",
  generator? : String = "",
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

///|
fn nexml_parse_float(s : String) -> Double? {
  if s == "." || s.length() == 0 {
    return None
  }
  let mut result = 0.0
  let mut sign = 1.0
  let mut i = 0
  let zero = '0'.to_int()
  
  if s.unsafe_get(0).to_int() == '-'.to_int() {
    sign = -1.0
    i = 1
  }
  
  let mut has_dot = false
  let mut decimal_place = 1.0
  
  while i < s.length() {
    let c = s.unsafe_get(i).to_int()
    if c == '.'.to_int() {
      has_dot = true
      i = i + 1
      continue
    }
    if c >= zero && c <= zero + 9 {
      let digit = (c - zero).to_double()
      result = result * 10.0 + digit
      if has_dot {
        decimal_place = decimal_place * 10.0
      }
    }
    i = i + 1
  }
  
  Some(result * sign / decimal_place)
}

// ============================================================================
// NeXML Parsing
// ============================================================================

///|
/// Parse a NeXML document from XML string.
pub fn parse_nexml(xml_content : String) -> NeXMLDocument {
  let mut doc = NeXMLDocument::new()
  
  // Check if there's a nexml element
  let (nexml_start, _, _) = find_element(xml_content, "nexml")
  if nexml_start == 0 && xml_content.length() > 0 {
    // Try to get version from root
    let nexml_str = xml_content[0:].to_owned()
    let version_attr = get_attribute(nexml_str, "version")
    let generator_attr = get_attribute(nexml_str, "generator")
    if version_attr.length() > 0 {
      doc = NeXMLDocument::{ ..doc, version: version_attr }
    }
    if generator_attr.length() > 0 {
      doc = NeXMLDocument::{ ..doc, generator: generator_attr }
    }
  }
  
  // Parse OTUs blocks
  doc = parse_otus_blocks(xml_content, doc)
  
  // Parse characters blocks
  doc = parse_characters_blocks(xml_content, doc)
  
  // Parse trees blocks
  doc = parse_trees_blocks(xml_content, doc)
  
  doc
}

///|
fn parse_otus_blocks(xml : String, doc : NeXMLDocument) -> NeXMLDocument {
  let otus_elements = find_all_elements(xml, "otus")
  let new_otus_blocks = doc.otus_blocks.copy()
  
  for i = 0; i < otus_elements.length(); i = i + 1 {
    let (start, end) = otus_elements[i]
    let otus_xml = xml[start:end].to_owned()
    let otus_id = get_attribute(otus_xml, "id")
    
    let otu_list : Array[NeXMLOTU] = Array::new()
    let otu_elements = find_all_elements(otus_xml, "otu")
    
    for j = 0; j < otu_elements.length(); j = j + 1 {
      let (otu_start, otu_end) = otu_elements[j]
      let otu_xml = otus_xml[otu_start:otu_end].to_owned()
      let id = get_attribute(otu_xml, "id")
      let label = get_attribute(otu_xml, "label")
      if label.length() == 0 {
        otu_list.push(NeXMLOTU::new(id, id))
      } else {
        otu_list.push(NeXMLOTU::new(id, label))
      }
    }
    
    new_otus_blocks.push(NeXMLOTUs::new(otus_id, otu_list))
  }
  
  NeXMLDocument::{ ..doc, otus_blocks: new_otus_blocks }
}

///|
fn parse_characters_blocks(xml : String, doc : NeXMLDocument) -> NeXMLDocument {
  let chars_elements = find_all_elements(xml, "characters")
  let new_chars_blocks = doc.characters_blocks.copy()
  
  for i = 0; i < chars_elements.length(); i = i + 1 {
    let (start, end) = chars_elements[i]
    let chars_xml = xml[start:end].to_owned()
    let chars_id = get_attribute(chars_xml, "id")
    let otus_id = get_attribute(chars_xml, "otus")
    let matrix_type = get_attribute(chars_xml, "xsi:type")
    if matrix_type.length() == 0 {
      let xtype = get_attribute(chars_xml, "type")
      let final_type = if xtype.length() > 0 { xtype } else { "StandardCells" }
      
      let char_list : Array[NeXMLCharacter] = Array::new()
      let char_elements = find_all_elements(chars_xml, "char")
      
      for j = 0; j < char_elements.length(); j = j + 1 {
        let (c_start, c_end) = char_elements[j]
        let char_xml = chars_xml[c_start:c_end].to_owned()
        let cid = get_attribute(char_xml, "id")
        let clabel = get_attribute(char_xml, "label")
        char_list.push(NeXMLCharacter::new(cid, clabel, []))
      }
      
      let matrix : Map[String, Array[String]] = Map([], capacity=16)
      let row_elements = find_all_elements(chars_xml, "row")
      
      for j = 0; j < row_elements.length(); j = j + 1 {
        let (r_start, r_end) = row_elements[j]
        let row_xml = chars_xml[r_start:r_end].to_owned()
        let otu = get_attribute(row_xml, "otu")
        let cell_elements = find_all_elements(row_xml, "cell")
        let row_data : Array[String] = Array::new()
        
        for k = 0; k < cell_elements.length(); k = k + 1 {
          let (c_start2, c_end2) = cell_elements[k]
          let cell_xml = row_xml[c_start2:c_end2].to_owned()
          let state = get_attribute(cell_xml, "state")
          row_data.push(state)
        }
        
        matrix.set(otu, row_data)
      }
      
      new_chars_blocks.push(
        NeXMLCharacters::new(chars_id, otus_id, final_type, char_list, matrix)
      )
    } else {
      let char_list : Array[NeXMLCharacter] = Array::new()
      let char_elements = find_all_elements(chars_xml, "char")
      
      for j = 0; j < char_elements.length(); j = j + 1 {
        let (c_start, c_end) = char_elements[j]
        let char_xml = chars_xml[c_start:c_end].to_owned()
        let cid = get_attribute(char_xml, "id")
        let clabel = get_attribute(char_xml, "label")
        char_list.push(NeXMLCharacter::new(cid, clabel, []))
      }
      
      let matrix : Map[String, Array[String]] = Map([], capacity=16)
      let row_elements = find_all_elements(chars_xml, "row")
      
      for j = 0; j < row_elements.length(); j = j + 1 {
        let (r_start, r_end) = row_elements[j]
        let row_xml = chars_xml[r_start:r_end].to_owned()
        let otu = get_attribute(row_xml, "otu")
        let cell_elements = find_all_elements(row_xml, "cell")
        let row_data : Array[String] = Array::new()
        
        for k = 0; k < cell_elements.length(); k = k + 1 {
          let (c_start2, c_end2) = cell_elements[k]
          let cell_xml = row_xml[c_start2:c_end2].to_owned()
          let state = get_attribute(cell_xml, "state")
          row_data.push(state)
        }
        
        matrix.set(otu, row_data)
      }
      
      new_chars_blocks.push(
        NeXMLCharacters::new(chars_id, otus_id, matrix_type, char_list, matrix)
      )
    }
  }
  
  NeXMLDocument::{ ..doc, characters_blocks: new_chars_blocks }
}

///|
fn parse_trees_blocks(xml : String, doc : NeXMLDocument) -> NeXMLDocument {
  let trees_elements = find_all_elements(xml, "trees")
  let new_trees_blocks = doc.trees_blocks.copy()
  
  for i = 0; i < trees_elements.length(); i = i + 1 {
    let (start, end) = trees_elements[i]
    let trees_xml = xml[start:end].to_owned()
    let trees_id = get_attribute(trees_xml, "id")
    let otus_id = get_attribute(trees_xml, "otus")
    
    let tree_list : Array[NeXMLTree] = Array::new()
    let tree_elements = find_all_elements(trees_xml, "tree")
    
    for j = 0; j < tree_elements.length(); j = j + 1 {
      let (t_start, t_end) = tree_elements[j]
      let tree_xml = trees_xml[t_start:t_end].to_owned()
      let tid = get_attribute(tree_xml, "id")
      let tname = get_attribute(tree_xml, "label")
      let rooted_str = get_attribute(tree_xml, "rooted")
      let rooted = rooted_str == "true"
      let ttype = get_attribute(tree_xml, "xsi:type")
      let final_type = if ttype.length() > 0 { ttype } else { "FloatTree" }
      
      let node_list : Array[NeXMLNode] = Array::new()
      let node_elements = find_all_elements(tree_xml, "node")
      
      for k = 0; k < node_elements.length(); k = k + 1 {
        let (n_start, n_end) = node_elements[k]
        let node_xml = tree_xml[n_start:n_end].to_owned()
        let nid = get_attribute(node_xml, "id")
        let nlabel = get_attribute(node_xml, "label")
        let notu = get_attribute(node_xml, "otu")
        
        let length_str = get_attribute(node_xml, "length")
        let nlength : Double? = if length_str.length() > 0 {
          nexml_parse_float(length_str)
        } else {
          None
        }
        
        node_list.push(NeXMLNode::new(nid, nlabel, notu, length=nlength))
      }
      
      let edge_list : Array[NeXMLEdge] = Array::new()
      let edge_elements = find_all_elements(tree_xml, "edge")
      
      for k = 0; k < edge_elements.length(); k = k + 1 {
        let (e_start, e_end) = edge_elements[k]
        let edge_xml = tree_xml[e_start:e_end].to_owned()
        let eid = get_attribute(edge_xml, "id")
        let esource = get_attribute(edge_xml, "source")
        let etarget = get_attribute(edge_xml, "target")
        
        let length_str = get_attribute(edge_xml, "length")
        let elength : Double? = if length_str.length() > 0 {
          nexml_parse_float(length_str)
        } else {
          None
        }
        
        edge_list.push(NeXMLEdge::new(eid, esource, etarget, length=elength))
      }
      
      tree_list.push(
        NeXMLTree::new(tid, tname, final_type, rooted, node_list, edge_list)
      )
    }
    
    new_trees_blocks.push(NeXMLTrees::new(trees_id, otus_id, tree_list))
  }
  
  NeXMLDocument::{ ..doc, trees_blocks: new_trees_blocks }
}

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
      xml = xml + "    <otu id=\"" + otu.id + "\" label=\"" + otu.label + "\"/>\n"
    }
    xml = xml + "  </otus>\n"
  }
  
  // Write characters blocks
  for i = 0; i < doc.characters_blocks.length(); i = i + 1 {
    let chars_block = doc.characters_blocks[i]
    xml = xml + "  <characters id=\"" + chars_block.id +
          "\" otus=\"" + chars_block.otus_id + "\">\n"
    for j = 0; j < chars_block.characters.length(); j = j + 1 {
      let ch = chars_block.characters[j]
      xml = xml + "    <char id=\"" + ch.id + "\" label=\"" + ch.label + "\"/>\n"
    }
    xml = xml + "    <matrix>\n"
    for otu in chars_block.matrix.keys() {
      let row_opt = chars_block.matrix.get(otu)
      if row_opt.is_some() {
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
    xml = xml + "  <trees id=\"" + trees_block.id +
          "\" otus=\"" + trees_block.otus_id + "\">\n"
    for j = 0; j < trees_block.trees.length(); j = j + 1 {
      let tree = trees_block.trees[j]
      let rooted_str = if tree.rooted { "true" } else { "false" }
      xml = xml + "    <tree id=\"" + tree.id +
            "\" label=\"" + tree.name +
            "\" rooted=\"" + rooted_str + "\">\n"
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
        xml = xml + "      <edge id=\"" + edge.id +
              "\" source=\"" + edge.source +
              "\" target=\"" + edge.target + "\""
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
  let children_map : Map[String, Array[(String, Double?)]] = Map([], capacity=16)
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
    if parent_map.get(node.id).is_none() {
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
  let mut node_length : Double? = None
  for i = 0; i < tree.nodes.length(); i = i + 1 {
    if tree.nodes[i].id == node_id {
      node_label = tree.nodes[i].label
      node_length = tree.nodes[i].length
      break
    }
  }
  
  match children_opt {
    Some(children) => {
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
    }
    None => {
      // Leaf node (no children in map)
      (node_label, true)
    }
  }
}

// ============================================================================
// Utility functions
// ============================================================================

///|
/// Get all OTU labels from a NeXML document.
pub fn nexml_get_all_otu_labels(doc : NeXMLDocument) -> Array[String] {
  let labels : Array[String] = Array::new()
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
  let names : Array[String] = Array::new()
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
pub fn nexml_get_tree_count(doc : NeXMLDocument) -> Int {
  let mut count = 0
  for i = 0; i < doc.trees_blocks.length(); i = i + 1 {
    count = count + doc.trees_blocks[i].trees.length()
  }
  count
}

///|
/// Get the number of OTUs in a NeXML document.
pub fn nexml_get_otu_count(doc : NeXMLDocument) -> Int {
  let mut count = 0
  for i = 0; i < doc.otus_blocks.length(); i = i + 1 {
    count = count + doc.otus_blocks[i].otus.length()
  }
  count
}
