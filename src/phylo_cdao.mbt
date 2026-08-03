///|
/// Bio.Phylo.CDAO - CDAO (Comparative Data Analysis Ontology) format support.
/// Inspired by Biopython's Bio.Phylo.CDAO and Bio.Phylo.CDAOIO.
///
/// CDAO is an RDF-based format for representing phylogenetic data using the
/// Comparative Data Analysis Ontology. It serializes trees as RDF triples
/// with CDAO ontology terms.
///
/// Key CDAO ontology terms:
///   - cdao:Tree          - A phylogenetic tree
///   - cdao:Node          - A node in the tree
///   - cdao:Edge          - An edge (branch) in the tree
///   - cdao:has_Root      - Links a tree to its root node
///   - cdao:has_Ancestor  - Links a node to its parent node
///   - cdao:belongs_to_TU - Links a node to its taxonomic unit
///   - cdao:TU            - A taxonomic unit (OTU/leaf label)
///
/// This module provides a simplified CDAO parser/serializer that does not
/// require an RDF triple store (librdf). It parses the RDF/XML structure
/// directly to extract tree topology and labels.
///
/// Reference: http://evolutionaryontology.org/cdao

// ============================================================================
// CDAO Namespace Constants
// ============================================================================

///|
/// CDAO namespace URI.
pub fn cdao_namespace() -> String {
  "http://purl.obolibrary.org/obo/cdao.owl#"
}

///|
/// RDF namespace URI.
pub fn cdao_rdf_namespace() -> String {
  "http://www.w3.org/1999/02/22-rdf-syntax-ns#"
}

///|
/// RDFS namespace URI.
pub fn cdao_rdfs_namespace() -> String {
  "http://www.w3.org/2000/01/rdf-schema#"
}

// ============================================================================
// CDAO Data Structures
// ============================================================================

///|
/// Represents a CDAO tree with metadata.
pub struct Cdaotree {
  id : String
  rooted : Bool
  root_node_id : String
  name : String?
}

///|
/// Construct a Cdaotree.
pub fn Cdaotree::new(
  id : String,
  rooted? : Bool = false,
  root_node_id? : String = "",
  name? : String? = None
) -> Cdaotree {
  Cdaotree::{ id, rooted, root_node_id, name }
}

///|
/// Represents a CDAO node.
pub struct CdaoNode {
  children : Array[String]
  mut parent_id : String?
  mut tu_id : String?
  mut branch_length : Double?
  mut label : String?
  id : String
}

///|
/// Construct a CdaoNode.
pub fn CdaoNode::new(
  id : String,
  children? : Array[String] = [],
  parent_id? : String? = None,
  tu_id? : String? = None,
  branch_length? : Double? = None,
  label? : String? = None
) -> CdaoNode {
  CdaoNode::{ id, children, parent_id, tu_id, branch_length, label }
}

///|
/// Represents a CDAO Taxonomic Unit (TU).
pub struct CdaoTU {
  id : String
  mut label : String?
}

///|
/// Construct a CdaoTU.
pub fn CdaoTU::new(id : String, label? : String? = None) -> CdaoTU {
  CdaoTU::{ id, label }
}

///|
/// Represents a complete CDAO document (parsed RDF graph).
pub struct CdaoDocument {
  trees : Array[Cdaotree]
  nodes : Map[String, CdaoNode]
  tus : Map[String, CdaoTU]
}

///|
/// Construct an empty CdaoDocument.
pub fn CdaoDocument::new() -> CdaoDocument {
  CdaoDocument::{
    trees: [],
    nodes: Map([], capacity=16),
    tus: Map([], capacity=16),
  }
}

// ============================================================================
// CDAO Parser (RDF/XML to CdaoDocument)
// ============================================================================

///|
/// Parse a CDAO RDF/XML string into a CdaoDocument.
pub fn cdao_parse(xml : String) -> CdaoDocument {
  let doc = CdaoDocument::new()
  let triples = cdao_parse_rdf_xml(xml)
  cdao_build_document(doc, triples)
  doc
}

///|
/// An RDF triple (subject, predicate, object).
struct CdaoTriple {
  subject : String
  predicate : String
  object : String
  is_literal : Bool
}

///|
/// Parse RDF/XML into a list of triples.
fn cdao_parse_rdf_xml(xml : String) -> Array[CdaoTriple] {
  let triples : Array[CdaoTriple] = Array::new()
  let normalized = cdao_normalize_xml(xml)
  cdao_extract_elements(normalized, triples)
  triples
}

///|
/// Normalize XML by removing excessive whitespace.
fn cdao_normalize_xml(xml : String) -> String {
  let mut result = ""
  let mut in_tag = false
  let mut i = 0
  while i < xml.length() {
    let c = xml[i:i + 1].to_owned()
    if c == "<" {
      in_tag = true
      result = result + c
    } else if c == ">" {
      in_tag = false
      result = result + c
    } else if in_tag && (c == " " || c == "\t" || c == "\n" || c == "\r") {
      if result.length() > 0 {
        let last = result[result.length() - 1:result.length()].to_owned()
        if last != " " && last != "<" {
          result = result + " "
        }
      }
    } else if !in_tag && (c == "\t" || c == "\r" || c == "\n") {
      result = result + " "
    } else {
      result = result + c
    }
    i = i + 1
  }
  result
}

///|
/// Extract RDF elements and their properties from normalized XML.
fn cdao_extract_elements(xml : String, triples : Array[CdaoTriple]) -> Unit {
  let mut pos = 0
  while pos < xml.length() {
    let tag_start = cdao_find_substring(xml, "<", pos)
    if tag_start < 0 {
      break
    }
    let tag_end = cdao_find_substring(xml, ">", tag_start)
    if tag_end < 0 {
      break
    }
    let tag_content = xml[tag_start + 1:tag_end].to_owned()
    if tag_content[0:1].to_owned() == "/" || tag_content[0:4].to_owned() == "!--" || tag_content.contains("rdf:RDF") {
      pos = tag_end + 1
      continue
    }
    let is_self_closing = tag_content[tag_content.length() - 1:tag_content.length()].to_owned() == "/"
    let tag_name = cdao_get_tag_name(tag_content)
    let about = cdao_get_attribute(tag_content, "rdf:about")
    if tag_name.contains("cdao:Tree") || tag_name.contains("cdao:Node") || tag_name.contains("cdao:TU") || tag_name.contains("rdf:Description") {
      let subject = match about {
        Some(a) => a
        None => ""
      }
      if subject.length() > 0 {
        triples.push(CdaoTriple::{
          subject,
          predicate: "rdf:type",
          object: tag_name,
          is_literal: false,
        })
        if is_self_closing {
          pos = tag_end + 1
        } else {
          let close_tag = "</" + tag_name + ">"
          let close_pos = cdao_find_substring(xml, close_tag, tag_end + 1)
          if close_pos >= 0 {
            let inner = xml[tag_end + 1:close_pos].to_owned()
            cdao_parse_children(subject, inner, triples)
            pos = close_pos + close_tag.length()
          } else {
            pos = tag_end + 1
          }
        }
      } else {
        pos = tag_end + 1
      }
    } else {
      pos = tag_end + 1
    }
  }
}

///|
/// Parse child elements of an RDF resource element.
fn cdao_parse_children(
  subject : String,
  inner : String,
  triples : Array[CdaoTriple]
) -> Unit {
  let mut pos = 0
  while pos < inner.length() {
    let tag_start = cdao_find_substring(inner, "<", pos)
    if tag_start < 0 {
      break
    }
    let tag_end = cdao_find_substring(inner, ">", tag_start)
    if tag_end < 0 {
      break
    }
    let tag_content = inner[tag_start + 1:tag_end].to_owned()
    if tag_content[0:1].to_owned() == "/" || tag_content[0:4].to_owned() == "!--" {
      pos = tag_end + 1
      continue
    }
    let is_self_closing = tag_content[tag_content.length() - 1:tag_content.length()].to_owned() == "/"
    let tag_name = cdao_get_tag_name(tag_content)
    let resource = cdao_get_attribute(tag_content, "rdf:resource")
    let close_tag = "</" + tag_name + ">"
    let close_pos = cdao_find_substring(inner, close_tag, tag_end + 1)
    let text_content = if is_self_closing {
      ""
    } else if close_pos >= 0 {
      inner[tag_end + 1:close_pos].trim().to_string()
    } else {
      ""
    }
    if resource is Some(res) {
      triples.push(CdaoTriple::{
        subject,
        predicate: tag_name,
        object: res,
        is_literal: false,
      })
    } else if text_content.length() > 0 {
      triples.push(CdaoTriple::{
        subject,
        predicate: tag_name,
        object: text_content,
        is_literal: true,
      })
    }
    if is_self_closing {
      pos = tag_end + 1
    } else if close_pos >= 0 {
      pos = close_pos + close_tag.length()
    } else {
      pos = tag_end + 1
    }
  }
}

// ============================================================================
// CDAO Document Builder (triples to CdaoDocument)
// ============================================================================

///|
/// Build a CdaoDocument from parsed RDF triples.
fn cdao_build_document(doc : CdaoDocument, triples : Array[CdaoTriple]) -> Unit {
  let mut i = 0
  while i < triples.length() {
    let t = triples[i]
    if t.predicate == "rdf:type" {
      if t.object.contains("cdao:Tree") {
        // Tree type
      } else if t.object.contains("cdao:Node") {
        doc.nodes[t.subject] = CdaoNode::new(t.subject)
      } else if t.object.contains("cdao:TU") {
        doc.tus[t.subject] = CdaoTU::new(t.subject)
      }
    }
    i = i + 1
  }
  // Handle resources that have CDAO properties but no explicit type
  let mut j = 0
  while j < triples.length() {
    let t = triples[j]
    if !doc.nodes.contains(t.subject) && !doc.tus.contains(t.subject) {
      if t.predicate.contains("cdao:") {
        if t.predicate.contains("has_Child") || t.predicate.contains("has_Ancestor") || t.predicate.contains("belongs_to_TU") {
          doc.nodes[t.subject] = CdaoNode::new(t.subject)
        }
      }
    }
    j = j + 1
  }
  // Populate properties
  let mut k = 0
  while k < triples.length() {
    let t = triples[k]
    if t.predicate == "rdf:type" {
      // Skip
    } else if t.predicate.contains("has_Root") {
      let mut found = false
      let mut idx = 0
      while idx < doc.trees.length() {
        if doc.trees[idx].id == t.subject {
          found = true
          break
        }
        idx = idx + 1
      }
      if !found {
        doc.trees.push(Cdaotree::new(t.subject, rooted=true, root_node_id=t.object))
      }
    } else if t.predicate.contains("has_Child") || t.predicate.contains("has_Descendant") {
      let node = doc.nodes.get(t.subject)
      match node {
        Some(n) => n.children.push(t.object)
        None => ()
      }
    } else if t.predicate.contains("has_Ancestor") || t.predicate.contains("has_Parent") {
      let node = doc.nodes.get(t.subject)
      match node {
        Some(n) => n.parent_id = Some(t.object)
        None => ()
      }
    } else if t.predicate.contains("belongs_to_TU") {
      let node = doc.nodes.get(t.subject)
      match node {
        Some(n) => n.tu_id = Some(t.object)
        None => ()
      }
    } else if t.predicate == "rdfs:label" || t.predicate.contains("has_Label") {
      let tu = doc.tus.get(t.subject)
      match tu {
        Some(tu_val) => tu_val.label = Some(t.object)
        None => {
          let node = doc.nodes.get(t.subject)
          match node {
            Some(n) => n.label = Some(t.object)
            None => ()
          }
        }
      }
    } else if t.predicate.contains("has_branch_length") {
      let node = doc.nodes.get(t.subject)
      match node {
        Some(n) => n.branch_length = Some(cdao_parse_double(t.object))
        None => ()
      }
    }
    k = k + 1
  }
}

// ============================================================================
// CDAO to BioSeqs Tree conversion
// ============================================================================

///|
/// Convert a CdaoDocument to an array of BioSeqs Tree objects.
pub fn cdao_to_trees(doc : CdaoDocument) -> Array[Tree] {
  let result : Array[Tree] = Array::new()
  let mut i = 0
  while i < doc.trees.length() {
    let cdao_tree = doc.trees[i]
    let root_clade = cdao_build_clade(doc, cdao_tree.root_node_id)
    let tree = Tree::new(
      root_clade,
      rooted=cdao_tree.rooted,
      name=cdao_tree.name,
    )
    result.push(tree)
    i = i + 1
  }
  result
}

///|
/// Recursively build a Clade from CDAO nodes.
fn cdao_build_clade(doc : CdaoDocument, node_id : String) -> Clade {
  let node = doc.nodes.get(node_id)
  match node {
    Some(n) => {
      let label = match n.tu_id {
        Some(tu_id) =>
          match doc.tus.get(tu_id) {
            Some(tu) =>
              match tu.label {
                Some(l) => Some(l)
                None => Some(tu_id)
              }
            None => Some(tu_id)
          }
        None =>
          match n.label {
            Some(l) => Some(l)
            None => None
          }
      }
      let children : Array[Clade] = Array::new()
      let mut ci = 0
      while ci < n.children.length() {
        children.push(cdao_build_clade(doc, n.children[ci]))
        ci = ci + 1
      }
      Clade::new(
        branch_length=n.branch_length,
        name=label,
        clades=children,
      )
    }
    None => Clade::new(name=Some(node_id))
  }
}

// ============================================================================
// CDAO Writer (Tree to RDF/XML)
// ============================================================================

///|
/// Mutable counter holder for CDAO writer.
struct CdaoWriteState {
  mut node_counter : Int
  mut tu_counter : Int
  tu_map : Map[String, String]
}

///|
/// Serialize a Tree to CDAO RDF/XML format.
pub fn cdao_write(tree : Tree) -> String {
  let state = CdaoWriteState::{
    node_counter: 1,
    tu_counter: 0,
    tu_map: Map([], capacity=16),
  }
  let mut output = ""
  output = output + "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n"
  output = output + "<rdf:RDF\n"
  output = output + "  xmlns:rdf=\"" + cdao_rdf_namespace() + "\"\n"
  output = output + "  xmlns:rdfs=\"" + cdao_rdfs_namespace() + "\"\n"
  output = output + "  xmlns:cdao=\"" + cdao_namespace() + "\">\n"
  output = output + "  <cdao:Tree rdf:about=\"#tree1\">\n"
  output = output + "    <cdao:has_Root rdf:resource=\"#node1\"/>\n"
  output = output + "  </cdao:Tree>\n"
  let (node_output, _) = cdao_write_clade(tree.root, "#node1", state)
  output = output + node_output
  // TU elements
  let tu_keys = state.tu_map.keys()
  let tu_arr : Array[String] = tu_keys.collect()
  let mut ti = 0
  while ti < tu_arr.length() {
    let tu_label = tu_arr[ti]
    let tu_id = state.tu_map[tu_label]
    output = output + "  <cdao:TU rdf:about=\"" + tu_id + "\">\n"
    output = output + "    <rdfs:label>" + cdao_escape_xml(tu_label) + "</rdfs:label>\n"
    output = output + "  </cdao:TU>\n"
    ti = ti + 1
  }
  output = output + "</rdf:RDF>\n"
  output
}

///|
/// Recursively write a clade and its children as CDAO nodes.
/// Returns (output_string, updated_state).
fn cdao_write_clade(
  clade : Clade,
  node_uri : String,
  state : CdaoWriteState
) -> (String, CdaoWriteState) {
  let mut output = ""
  output = output + "  <cdao:Node rdf:about=\"" + node_uri + "\">\n"
  match clade.branch_length {
    Some(bl) => {
      output = output + "    <cdao:has_branch_length>" + bl.to_string() + "</cdao:has_branch_length>\n"
    }
    None => ()
  }
  if clade.is_terminal() {
    let label = match clade.name {
      Some(n) => n
      None => "unknown"
    }
    let tu_id_opt = state.tu_map.get(label)
    let tu_uri = match tu_id_opt {
      Some(id) => id
      None => {
        state.tu_counter = state.tu_counter + 1
        let new_id = "#tu" + state.tu_counter.to_string()
        state.tu_map[label] = new_id
        new_id
      }
    }
    output = output + "    <cdao:belongs_to_TU rdf:resource=\"" + tu_uri + "\"/>\n"
  }
  // Collect child URIs
  let child_uris : Array[String] = Array::new()
  let mut ci = 0
  while ci < clade.clades.length() {
    state.node_counter = state.node_counter + 1
    let child_uri = "#node" + state.node_counter.to_string()
    child_uris.push(child_uri)
    output = output + "    <cdao:has_Child rdf:resource=\"" + child_uri + "\"/>\n"
    ci = ci + 1
  }
  output = output + "  </cdao:Node>\n"
  // Recursively write children
  let mut di = 0
  while di < clade.clades.length() {
    let (child_output, _) = cdao_write_clade(clade.clades[di], child_uris[di], state)
    output = output + child_output
    di = di + 1
  }
  (output, state)
}

// ============================================================================
// CDAO Utility Functions
// ============================================================================

///|
/// Create a sample CDAO document for testing.
/// Tree topology: ((A, B), (C, D))
pub fn cdao_sample_document() -> CdaoDocument {
  let doc = CdaoDocument::new()
  doc.trees.push(Cdaotree::new("#tree1", rooted=true, root_node_id="#node1"))
  doc.nodes["#node1"] = CdaoNode::new("#node1", children=["#node2", "#node3"])
  doc.nodes["#node2"] = CdaoNode::new("#node2", children=["#node4", "#node5"], parent_id=Some("#node1"))
  doc.nodes["#node3"] = CdaoNode::new("#node3", children=["#node6", "#node7"], parent_id=Some("#node1"))
  doc.nodes["#node4"] = CdaoNode::new("#node4", parent_id=Some("#node2"), tu_id=Some("#tu1"), branch_length=Some(1.0))
  doc.nodes["#node5"] = CdaoNode::new("#node5", parent_id=Some("#node2"), tu_id=Some("#tu2"), branch_length=Some(1.0))
  doc.nodes["#node6"] = CdaoNode::new("#node6", parent_id=Some("#node3"), tu_id=Some("#tu3"), branch_length=Some(1.0))
  doc.nodes["#node7"] = CdaoNode::new("#node7", parent_id=Some("#node3"), tu_id=Some("#tu4"), branch_length=Some(1.0))
  doc.tus["#tu1"] = CdaoTU::new("#tu1", label=Some("A"))
  doc.tus["#tu2"] = CdaoTU::new("#tu2", label=Some("B"))
  doc.tus["#tu3"] = CdaoTU::new("#tu3", label=Some("C"))
  doc.tus["#tu4"] = CdaoTU::new("#tu4", label=Some("D"))
  doc
}

///|
/// Create a sample CDAO RDF/XML string for testing.
pub fn cdao_sample_xml() -> String {
  "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n" +
  "<rdf:RDF xmlns:rdf=\"http://www.w3.org/1999/02/22-rdf-syntax-ns#\"\n" +
  "         xmlns:rdfs=\"http://www.w3.org/2000/01/rdf-schema#\"\n" +
  "         xmlns:cdao=\"http://purl.obolibrary.org/obo/cdao.owl#\">\n" +
  "  <cdao:Tree rdf:about=\"#tree1\">\n" +
  "    <cdao:has_Root rdf:resource=\"#node1\"/>\n" +
  "  </cdao:Tree>\n" +
  "  <cdao:Node rdf:about=\"#node1\">\n" +
  "    <cdao:has_Child rdf:resource=\"#node2\"/>\n" +
  "    <cdao:has_Child rdf:resource=\"#node3\"/>\n" +
  "  </cdao:Node>\n" +
  "  <cdao:Node rdf:about=\"#node2\">\n" +
  "    <cdao:has_Ancestor rdf:resource=\"#node1\"/>\n" +
  "    <cdao:has_Child rdf:resource=\"#node4\"/>\n" +
  "    <cdao:has_Child rdf:resource=\"#node5\"/>\n" +
  "  </cdao:Node>\n" +
  "  <cdao:Node rdf:about=\"#node3\">\n" +
  "    <cdao:has_Ancestor rdf:resource=\"#node1\"/>\n" +
  "    <cdao:has_Child rdf:resource=\"#node6\"/>\n" +
  "    <cdao:has_Child rdf:resource=\"#node7\"/>\n" +
  "  </cdao:Node>\n" +
  "  <cdao:Node rdf:about=\"#node4\">\n" +
  "    <cdao:has_Ancestor rdf:resource=\"#node2\"/>\n" +
  "    <cdao:belongs_to_TU rdf:resource=\"#tu1\"/>\n" +
  "  </cdao:Node>\n" +
  "  <cdao:Node rdf:about=\"#node5\">\n" +
  "    <cdao:has_Ancestor rdf:resource=\"#node2\"/>\n" +
  "    <cdao:belongs_to_TU rdf:resource=\"#tu2\"/>\n" +
  "  </cdao:Node>\n" +
  "  <cdao:Node rdf:about=\"#node6\">\n" +
  "    <cdao:has_Ancestor rdf:resource=\"#node3\"/>\n" +
  "    <cdao:belongs_to_TU rdf:resource=\"#tu3\"/>\n" +
  "  </cdao:Node>\n" +
  "  <cdao:Node rdf:about=\"#node7\">\n" +
  "    <cdao:has_Ancestor rdf:resource=\"#node3\"/>\n" +
  "    <cdao:belongs_to_TU rdf:resource=\"#tu4\"/>\n" +
  "  </cdao:Node>\n" +
  "  <cdao:TU rdf:about=\"#tu1\">\n" +
  "    <rdfs:label>A</rdfs:label>\n" +
  "  </cdao:TU>\n" +
  "  <cdao:TU rdf:about=\"#tu2\">\n" +
  "    <rdfs:label>B</rdfs:label>\n" +
  "  </cdao:TU>\n" +
  "  <cdao:TU rdf:about=\"#tu3\">\n" +
  "    <rdfs:label>C</rdfs:label>\n" +
  "  </cdao:TU>\n" +
  "  <cdao:TU rdf:about=\"#tu4\">\n" +
  "    <rdfs:label>D</rdfs:label>\n" +
  "  </cdao:TU>\n" +
  "</rdf:RDF>\n"
}

///|
/// Get a summary of a CdaoDocument.
pub fn cdao_document_summary(doc : CdaoDocument) -> String {
  let mut s = "=== CDAO Document Summary ===\n"
  s = s + "Trees: " + doc.trees.length().to_string() + "\n"
  s = s + "Nodes: " + doc.nodes.size().to_string() + "\n"
  s = s + "TUs: " + doc.tus.size().to_string() + "\n"
  let mut i = 0
  while i < doc.trees.length() {
    let tree = doc.trees[i]
    s = s + "\nTree " + tree.id + " (rooted=" + tree.rooted.to_string() + ")\n"
    s = s + "  Root: " + tree.root_node_id + "\n"
    i = i + 1
  }
  s
}

// ============================================================================
// Internal XML Helpers
// ============================================================================

///|
/// Find the position of a substring, starting from a given index.
fn cdao_find_substring(s : String, sub : String, start : Int) -> Int {
  if start >= s.length() {
    return -1
  }
  let sub_len = sub.length()
  if sub_len == 0 {
    return start
  }
  let mut i = start
  while i <= s.length() - sub_len {
    let mut found = true
    let mut j = 0
    while j < sub_len {
      if s[i + j:i + j + 1].to_owned() != sub[j:j + 1].to_owned() {
        found = false
        break
      }
      j = j + 1
    }
    if found {
      return i
    }
    i = i + 1
  }
  -1
}

///|
/// Get the tag name from a tag content string.
fn cdao_get_tag_name(tag_content : String) -> String {
  let mut name = ""
  let mut i = 0
  while i < tag_content.length() {
    let c = tag_content[i:i + 1].to_owned()
    if c == " " || c == "/" || c == "\t" || c == "\n" {
      break
    }
    name = name + c
    i = i + 1
  }
  name
}

///|
/// Get an attribute value from a tag content string.
fn cdao_get_attribute(tag_content : String, attr_name : String) -> String? {
  let search = attr_name + "=\""
  let pos = cdao_find_substring(tag_content, search, 0)
  if pos < 0 {
    return None
  }
  let value_start = pos + search.length()
  let value_end = cdao_find_substring(tag_content, "\"", value_start)
  if value_end < 0 {
    return None
  }
  Some(tag_content[value_start:value_end].to_owned())
}

///|
/// Parse a double from a string.
fn cdao_parse_double(s : String) -> Double {
  let trimmed = s.trim().to_string()
  let neg = trimmed[0:1].to_owned() == "-"
  let digits = if neg { trimmed[1:trimmed.length()].to_string() } else { trimmed }
  let result = parse_double(digits)
  match result {
    Some(d) => if neg { -d } else { d }
    None => 0.0
  }
}

///|
/// Escape special XML characters.
fn cdao_escape_xml(s : String) -> String {
  let mut result = ""
  let mut i = 0
  while i < s.length() {
    let c = s[i:i + 1].to_owned()
    match c {
      "<" => result = result + "&lt;"
      ">" => result = result + "&gt;"
      "&" => result = result + "&amp;"
      "\"" => result = result + "&quot;"
      "'" => result = result + "&apos;"
      _ => result = result + c
    }
    i = i + 1
  }
  result
}
