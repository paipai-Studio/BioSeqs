///|
/// Bio.Phylo.PhyloXML - XML-based phylogenetic tree format.
/// Provides parsing, generation, and conversion of PhyloXML format
/// with support for rich taxon and sequence annotations.

// ============================================================================
// Data Structures
// ============================================================================

///|
/// Sequence annotation associated with a PhyloXML node.
pub struct SequenceAnnotation {
  sequence_id : String
  molecular_type : String
  symbol : String
  uri : String
} derive(Eq)

///|
/// Taxonomic information associated with a PhyloXML node.
pub struct TaxonInfo {
  taxon_id : String
  scientific_name : String
  common_name : String
  authority : String
  taxon_rank : String
} derive(Eq)

///|
/// Taxon namespace definition for PhyloXML.
pub struct TaxonNamespace {
  namespace_id : String
  name : String
} derive(Eq)

///|
/// A node (clade) in a PhyloXML tree.
pub struct PhyloXMLNode {
  node_id : String
  name : String
  branch_length : Double
  children : Array[PhyloXMLNode]
  is_leaf : Bool
  annotations : Map[String, String]
  sequence_annotation : SequenceAnnotation?
  taxon : TaxonInfo?
} derive(Eq)

///|
/// A phylogenetic tree in PhyloXML format.
pub struct PhyloXMLTree {
  tree_id : String
  name : String
  description : String
  root : PhyloXMLNode
  phylogeny_metadata : Map[String, String]
  taxon_namespaces : Array[TaxonNamespace]
} derive(Eq)

///|
/// Top-level result of PhyloXML parsing.
pub struct PhyloXMLResult {
  trees : Array[PhyloXMLTree]
  metadata : Map[String, String]
  taxon_namespaces : Array[TaxonNamespace]
} derive(Eq)

// ============================================================================
// Constructors
// ============================================================================

///|
pub fn SequenceAnnotation::new(
  sequence_id? : String = "",
  molecular_type? : String = "",
  symbol? : String = "",
  uri? : String = "",
) -> SequenceAnnotation {
  SequenceAnnotation::{ sequence_id, molecular_type, symbol, uri }
}

///|
pub fn TaxonInfo::new(
  taxon_id? : String = "",
  scientific_name? : String = "",
  common_name? : String = "",
  authority? : String = "",
  taxon_rank? : String = "",
) -> TaxonInfo {
  TaxonInfo::{ taxon_id, scientific_name, common_name, authority, taxon_rank }
}

///|
pub fn TaxonNamespace::new(
  namespace_id? : String = "",
  name? : String = "",
) -> TaxonNamespace {
  TaxonNamespace::{ namespace_id, name }
}

///|
pub fn PhyloXMLNode::new(
  node_id? : String = "",
  name? : String = "",
  branch_length? : Double = 0.0,
  children? : Array[PhyloXMLNode] = [],
  annotations? : Map[String, String] = Map([], capacity=8),
  sequence_annotation? : SequenceAnnotation? = None,
  taxon? : TaxonInfo? = None,
) -> PhyloXMLNode {
  let is_leaf = children.length() == 0
  PhyloXMLNode::{
    node_id,
    name,
    branch_length,
    children,
    is_leaf,
    annotations,
    sequence_annotation,
    taxon,
  }
}

///|
pub fn PhyloXMLTree::new(
  tree_id? : String = "",
  name? : String = "",
  description? : String = "",
  root? : PhyloXMLNode = PhyloXMLNode::new(),
  phylogeny_metadata? : Map[String, String] = Map([], capacity=8),
  taxon_namespaces? : Array[TaxonNamespace] = [],
) -> PhyloXMLTree {
  PhyloXMLTree::{
    tree_id,
    name,
    description,
    root,
    phylogeny_metadata,
    taxon_namespaces,
  }
}

///|
pub fn PhyloXMLResult::new(
  trees? : Array[PhyloXMLTree] = [],
  metadata? : Map[String, String] = Map([], capacity=8),
  taxon_namespaces? : Array[TaxonNamespace] = [],
) -> PhyloXMLResult {
  PhyloXMLResult::{ trees, metadata, taxon_namespaces }
}

// ============================================================================
// XML Low-Level Helpers
// ============================================================================

///|
/// Strip namespace prefix from a tag name (e.g., "px:name" -> "name").
fn xml_strip_n(tag : String) -> String {
  let colon = tag.find(":")
  if colon.is_some() {
    tag[colon.unwrap() + 1:].to_owned()
  } else {
    tag
  }
}

///|
/// Parse attributes from an opening tag's attribute string.
/// Input is the portion after "<tagname" and before ">".
fn xml_parse_attrs(attr_str : String) -> Map[String, String] {
  let result : Map[String, String] = Map([], capacity=8)
  let trimmed = trim(attr_str)
  if trimmed.length() == 0 {
    return result
  }
  let n = trimmed.length()
  let mut i = 0
  let space = ' '.to_int()
  let eq = '='.to_int()
  let gt = '>'.to_int()
  let dquote = '"'.to_int()
  let squote = '\''.to_int()
  while i < n {
    while i < n && trimmed.unsafe_get(i).to_int() == space {
      i = i + 1
    }
    if i >= n {
      break
    }
    let key_start = i
    while i < n && trimmed.unsafe_get(i).to_int() != eq &&
      trimmed.unsafe_get(i).to_int() != space {
      i = i + 1
    }
    let key = trimmed[key_start:i].to_owned()
    while i < n && trimmed.unsafe_get(i).to_int() == space {
      i = i + 1
    }
    if i >= n || trimmed.unsafe_get(i).to_int() != eq {
      result[key] = ""
      continue
    }
    i = i + 1
    while i < n && trimmed.unsafe_get(i).to_int() == space {
      i = i + 1
    }
    if i >= n {
      result[key] = ""
      break
    }
    let ch = trimmed.unsafe_get(i).to_int()
    if ch == dquote || ch == squote {
      i = i + 1
      let val_start = i
      while i < n && trimmed.unsafe_get(i).to_int() != ch {
        i = i + 1
      }
      let val = trimmed[val_start:i].to_owned()
      result[key] = val
      if i < n {
        i = i + 1
      }
    } else {
      let val_start = i
      while i < n && trimmed.unsafe_get(i).to_int() != space &&
        trimmed.unsafe_get(i).to_int() != gt {
        i = i + 1
      }
      result[key] = trimmed[val_start:i].to_owned()
    }
  }
  result
}

///|
/// Find the next occurrence of an opening tag for the given tag name.
/// Returns (position_after_tag_start, attributes, is_self_closing) or None.
fn xml_find_open_tag(
  xml : String,
  tag_name : String,
  start : Int,
) -> (Int, Map[String, String], Bool)? {
  let search_tag = tag_name
  let n = xml.length()
  let mut i = start
  while i < n {
    if char_at(xml, i) == "<" {
      let rest = xml[i + 1:].to_owned()
      let rest_len = rest.length()
      let mut j = 0
      while j < rest_len && char_at(rest, j) != ">" && char_at(rest, j) != "/" &&
        char_at(rest, j) != " " && char_at(rest, j) != "\t" {
        j = j + 1
      }
      let found_tag = xml_strip_n(rest[0:j].to_owned())
      if found_tag == search_tag {
        let after_tag = rest[j:].to_owned()
        let is_self_closing = after_tag.has_prefix("/>")
        let attr_str = if is_self_closing {
          after_tag[0:after_tag.length() - 2].to_owned()
        } else {
          let gt_pos = after_tag.find(">")
          if gt_pos.is_some() {
            after_tag[0:gt_pos.unwrap()].to_owned()
          } else {
            after_tag
          }
        }
        let attrs = xml_parse_attrs(attr_str)
        let tag_end = if is_self_closing {
          i + 1 + j + 1 + attr_str.length() + 1
        } else {
          i + 1 + j + 1 + attr_str.length() + 1
        }
        return Some((i, attrs, is_self_closing))
      }
      i = i + 1
    } else {
      i = i + 1
    }
  }
  None
}

///|
/// Find the matching closing tag for a given opening tag, handling nesting.
/// Returns the position of the matching "<" or None.
fn xml_find_close_tag(
  xml : String,
  tag_name : String,
  open_start : Int,
) -> Int? {
  let n = xml.length()
  let search_tag = tag_name
  let mut depth = 1
  let mut i = open_start + 1
  while i < n {
    if char_at(xml, i) == "<" {
      let rest = xml[i + 1:].to_owned()
      if rest.has_prefix("/") {
        let after_slash = rest[1:].to_owned()
        let mut j = 0
        let rlen = after_slash.length()
        while j < rlen && char_at(after_slash, j) != ">" &&
          char_at(after_slash, j) != " " && char_at(after_slash, j) != "\t" {
          j = j + 1
        }
        let close_tag = xml_strip_n(after_slash[0:j].to_owned())
        if close_tag == search_tag {
          depth = depth - 1
          if depth == 0 {
            return Some(i)
          }
        }
        i = i + 1
      } else if rest.has_prefix("!") {
        let gt = rest.find(">")
        if gt.is_some() {
          i = i + 1 + gt.unwrap() + 1
        } else {
          return None
        }
      } else {
        let mut j = 0
        let rlen = rest.length()
        while j < rlen && char_at(rest, j) != ">" &&
          char_at(rest, j) != "/" && char_at(rest, j) != " " &&
          char_at(rest, j) != "\t" {
          j = j + 1
        }
        let found_tag = xml_strip_n(rest[0:j].to_owned())
        let is_self_closing = j < rlen && char_at(rest, j) == "/"
        if found_tag == search_tag && !is_self_closing {
          depth = depth + 1
        }
        i = i + 1
      }
    } else {
      i = i + 1
    }
  }
  None
}

///|
/// Get the text content between an element's start and end tags.
fn xml_get_text_content(
  xml : String,
  open_start : Int,
  close_start : Int,
  tag_name : String,
) -> String {
  let n = xml.length()
  let search_tag = tag_name
  let tag_end = find_tag_end(xml, open_start, search_tag)
  if tag_end >= close_start {
    return ""
  }
  let content = xml[tag_end:close_start].to_owned()
  let trimmed = trim(content)
  if trimmed.length() == 0 {
    return ""
  }
  if has_nested_elements(trimmed) {
    return ""
  }
  trimmed
}

///|
/// Find the position after the closing ">" of the tag that starts at open_start.
fn find_tag_end(xml : String, open_start : Int, tag_name : String) -> Int {
  let search_tag = tag_name
  let n = xml.length()
  let mut i = open_start + 1
  while i < n {
    if char_at(xml, i) == "<" {
      let rest = xml[i + 1:].to_owned()
      let mut j = 0
      let rlen = rest.length()
      while j < rlen && char_at(rest, j) != ">" &&
        char_at(rest, j) != "/" && char_at(rest, j) != " " {
        j = j + 1
      }
      let found_tag = xml_strip_n(rest[0:j].to_owned())
      if found_tag == search_tag {
        let after = rest[j:].to_owned()
        if after.has_prefix("/>") {
          return i + j + 2 + 1
        }
        let gt = after.find(">")
        if gt.is_some() {
          return i + j + 1 + gt.unwrap() + 1
        }
      }
      i = i + 1
    } else if char_at(xml, i) == ">" {
      return i + 1
    } else {
      i = i + 1
    }
  }
  n
}

///|
/// Check if a string contains nested XML elements (child tags).
fn has_nested_elements(s : String) -> Bool {
  s.find("<").is_some()
}

///|
/// Extract a complete XML element by tag name.
/// Returns (attributes_map, children_text) or (empty_map, "").
fn xml_extract_element(xml : String, tag_name : String) -> (Map[String, String], String) {
  let open_opt = xml_find_open_tag(xml, tag_name, 0)
  if open_opt.is_none() {
    return (Map([], capacity=4), "")
  }
  let (open_start, attrs, is_self_closing) = open_opt.unwrap()
  if is_self_closing {
    return (attrs, "")
  }
  let close_opt = xml_find_close_tag(xml, tag_name, open_start)
  if close_opt.is_none() {
    return (attrs, "")
  }
  let close_start = close_opt.unwrap()
  let tag_end = find_tag_end(xml, open_start, tag_name)
  let children_text = xml[tag_end:close_start].to_owned()
  (attrs, children_text)
}

///|
/// Get text content between an open tag and close tag positions.
fn xml_get_text_between(xml : String, open_tag : Int, close_tag : Int) -> String {
  if open_tag >= close_tag {
    return ""
  }
  let content = xml[open_tag:close_tag].to_owned()
  trim(content)
}

// ============================================================================
// XML Parser
// ============================================================================

///|
/// Find first XML element by tag name.
/// Returns (start_pos, end_pos, content_string) or (0, 0, "").
pub fn find_element(
  xml : String,
  tag : String,
) -> (Int, Int, String) {
  let open_opt = xml_find_open_tag(xml, tag, 0)
  if open_opt.is_none() {
    return (0, 0, "")
  }
  let (open_start, _, is_self_closing) = open_opt.unwrap()
  if is_self_closing {
    let tag_end = find_tag_end(xml, open_start, tag)
    return (open_start, tag_end, "")
  }
  let close_opt = xml_find_close_tag(xml, tag, open_start)
  if close_opt.is_none() {
    return (open_start, xml.length(), "")
  }
  let close_start = close_opt.unwrap()
  let tag_end = find_tag_end(xml, open_start, tag)
  let content = xml[tag_end:close_start].to_owned()
  let end_pos = find_tag_end(xml, close_start, tag)
  return (open_start, end_pos, content)
}

///|
/// Find all XML elements by tag name (non-overlapping, top-level only within the search range).
/// Returns array of (start_pos, end_pos) tuples.
pub fn find_all_elements(
  xml : String,
  tag : String,
) -> Array[(Int, Int)] {
  let result : Array[(Int, Int)] = Array::new()
  let mut search_start = 0
  let n = xml.length()
  while search_start < n {
    let open_opt = xml_find_open_tag(xml, tag, search_start)
    if open_opt.is_none() {
      break
    }
    let (open_start, _, is_self_closing) = open_opt.unwrap()
    if is_self_closing {
      let tag_end = find_tag_end(xml, open_start, tag)
      result.push((open_start, tag_end))
      search_start = tag_end
    } else {
      let close_opt = xml_find_close_tag(xml, tag, open_start)
      if close_opt.is_none() {
        break
      }
      let close_start = close_opt.unwrap()
      let close_end = find_tag_end(xml, close_start, tag)
      result.push((open_start, close_end))
      search_start = close_end
      if search_start >= n {
        break
      }
    }
  }
  result
}

///|
/// Get an attribute value from an element at the given position.
pub fn get_attribute(
  xml : String,
  attr_name : String,
) -> String {
  let trimmed = trim(xml)
  if !trimmed.has_prefix("<") {
    return ""
  }
  let gt_pos = trimmed.find(">")
  if gt_pos.is_none() {
    return ""
  }
  let tag_part = trimmed[1:gt_pos.unwrap()].to_owned()
  let attrs = xml_parse_attrs(tag_part)
  attrs.get(attr_name).unwrap_or("")
}

///|
/// Get the text content of an element (leaf text only, no nested elements).
/// Handles both XML elements (with tags) and plain text content.
pub fn get_element_text(xml : String) -> String {
  let trimmed = trim(xml)
  if trimmed.length() == 0 {
    return ""
  }
  // If the content doesn't have XML tags, it's plain text
  if !trimmed.has_prefix("<") {
    return trimmed
  }
  let open_gt = trimmed.find(">")
  if open_gt.is_none() {
    return ""
  }
  let start = open_gt.unwrap() + 1
  let n = trimmed.length()
  let mut last_lt = -1
  let mut i = n - 1
  while i >= start {
    if char_at(trimmed, i) == "<" {
      last_lt = i
      break
    }
    i = i - 1
  }
  if last_lt < 0 {
    return ""
  }
  let content = trimmed[start:last_lt].to_owned()
  trim(content)
}

///|
/// Parse a PhyloXML string into a PhyloXMLResult.
pub fn parse_phyloxml(xml_text : String) -> PhyloXMLResult {
  let trimmed = trim(xml_text)
  if trimmed.length() == 0 {
    return PhyloXMLResult::new()
  }

  let metadata : Map[String, String] = Map([], capacity=8)
  let namespaces : Array[TaxonNamespace] = Array::new()
  let trees : Array[PhyloXMLTree] = Array::new()

  let tree_blocks = find_all_elements(trimmed, "phylogeny")
  for block in tree_blocks {
    let (start, end) = block
    let tree_xml = trimmed[start:end].to_owned()
    let tree = parse_tree_element(tree_xml)
    trees.push(tree)
  }

  let ns_blocks = find_all_elements(trimmed, "taxon_namespace")
  for ns_block in ns_blocks {
    let (s, e) = ns_block
    let ns_xml = trimmed[s:e].to_owned()
    let ns = parse_taxon_namespace_element(ns_xml)
    namespaces.push(ns)
  }

  let id_blocks = find_all_elements(trimmed, "phylogeny_id")
  for id_block in id_blocks {
    let (s, e) = id_block
    let id_xml = trimmed[s:e].to_owned()
    let text = get_element_text(id_xml)
    if text.length() > 0 {
      metadata["phylogeny_id"] = text
    }
  }

  PhyloXMLResult::{ trees, metadata, taxon_namespaces: namespaces }
}

///|
/// Parse a PhyloXML file (content) into a PhyloXMLResult.
pub fn parse_phyloxml_file(file_path : String) -> PhyloXMLResult {
  parse_phyloxml(file_path)
}

///|
/// Parse a single phylogeny tree element.
pub fn parse_tree_element(xml : String) -> PhyloXMLTree {
  let trimmed = trim(xml)
  if trimmed.length() == 0 {
    return PhyloXMLTree::new()
  }

  let attrs_start = trimmed.find("<phylogeny")
  let attrs_end = if attrs_start.is_some() {
    let gt_opt = trimmed[attrs_start.unwrap():].find(">")
    if gt_opt.is_some() { Some(attrs_start.unwrap() + gt_opt.unwrap() + 1) } else { None }
  } else {
    None
  }

  let attrs_str = if attrs_start.is_some() && attrs_end.is_some() {
    trimmed[attrs_start.unwrap():attrs_end.unwrap()].to_owned()
  } else {
    ""
  }
  let attrs = xml_parse_attrs(attrs_str)

  let tree_id = attrs.get("id").unwrap_or("")
  let rooted = attrs.get("rooted").unwrap_or("") == "true"

  let name_text = parse_child_text(trimmed, "name")
  let desc_text = parse_child_text(trimmed, "description")
  let root_node = parse_child_node(trimmed, "clade")

  let meta : Map[String, String] = Map([], capacity=4)
  if rooted {
    meta["rooted"] = "true"
  }

  PhyloXMLTree::{
    tree_id,
    name: name_text,
    description: desc_text,
    root: root_node,
    phylogeny_metadata: meta,
    taxon_namespaces: [],
  }
}

///|
/// Parse a clade/node element recursively.
pub fn parse_node_element(xml : String) -> PhyloXMLNode {
  parse_child_node(xml, "clade")
}

///|
/// Parse taxon annotation from an XML element.
pub fn parse_taxon_element(xml : String) -> TaxonInfo {
  parse_taxon_info_from(if has_element(xml, "taxon") {
    extract_child(xml, "taxon")
  } else {
    xml
  })
}

///|
/// Parse sequence annotation from an XML element.
pub fn parse_sequence_element(xml : String) -> SequenceAnnotation {
  parse_seq_annotation_from(if has_element(xml, "sequence") {
    extract_child(xml, "sequence")
  } else {
    xml
  })
}

// ============================================================================
// Private Parsing Helpers
// ============================================================================

///|
fn has_element(xml : String, tag : String) -> Bool {
  let (s, e, _) = find_element(xml, tag)
  s != 0 || e != 0
}

///|
fn extract_child(xml : String, tag : String) -> String {
  let (s, e, _) = find_element(xml, tag)
  if s == 0 && e == 0 {
    return ""
  }
  xml[s:e].to_owned()
}

///|
fn parse_child_text(xml : String, tag : String) -> String {
  let (_, _, content) = find_element(xml, tag)
  if content.length() == 0 {
    get_element_text(extract_child(xml, tag))
  } else {
    get_element_text(content)
  }
}

///|
fn parse_child_node(xml : String, tag : String) -> PhyloXMLNode {
  let (s, e, content) = find_element(xml, tag)
  if s == 0 && e == 0 {
    return PhyloXMLNode::new()
  }
  let node_xml = if content.length() > 0 { content } else { xml[s:e].to_owned() }
  parse_node_from_content(node_xml)
}

///|
fn parse_node_from_content(xml : String) -> PhyloXMLNode {
  let name_text = parse_child_text(xml, "name")
  let bl_text = parse_child_text(xml, "branch_length")
  let bl_val = match parse_double(bl_text) {
    Some(v) => v
    None => 0.0
  }

  let children : Array[PhyloXMLNode] = Array::new()
  let mut search_start = 0
  let n = xml.length()
  while search_start < n {
    let open_opt = xml_find_open_tag(xml, "clade", search_start)
    if open_opt.is_none() {
      break
    }
    let (open_start, _, is_self_closing) = open_opt.unwrap()
    if is_self_closing {
      let tag_end = find_tag_end(xml, open_start, "clade")
      let child_xml = xml[open_start:tag_end].to_owned()
      children.push(parse_node_from_content(child_xml))
      search_start = tag_end
    } else {
      let close_opt = xml_find_close_tag(xml, "clade", open_start)
      if close_opt.is_none() {
        break
      }
      let close_start = close_opt.unwrap()
      let open_tag_end = find_tag_end(xml, open_start, "clade")
      let close_tag_end = find_tag_end(xml, close_start, "clade")
      let inner = xml[open_tag_end:close_start].to_owned()
      if inner.length() > 0 {
        let child = parse_node_from_content(inner)
        children.push(child)
      }
      search_start = close_tag_end
    }
  }

  let taxon = if has_element(xml, "taxon") {
    Some(parse_taxon_info_from(extract_child(xml, "taxon")))
  } else {
    None
  }

  let seq_ann = if has_element(xml, "sequence") {
    Some(parse_seq_annotation_from(extract_child(xml, "sequence")))
  } else {
    None
  }

  let anns : Map[String, String] = Map([], capacity=4)
  let ann_blocks = find_all_elements(xml, "annotation")
  for ann_block in ann_blocks {
    let (as_, ae) = ann_block
    let ann_xml = xml[as_:ae].to_owned()
    let ann_name = parse_child_text(ann_xml, "name")
    let ann_value = parse_child_text(ann_xml, "value")
    if ann_name.length() > 0 {
      anns[ann_name] = ann_value
    }
  }

  PhyloXMLNode::{
    node_id: "",
    name: name_text,
    branch_length: bl_val,
    children,
    is_leaf: children.length() == 0,
    annotations: anns,
    sequence_annotation: seq_ann,
    taxon,
  }
}
///|
fn parse_taxon_info_from(xml : String) -> TaxonInfo {
  TaxonInfo::{
    taxon_id: parse_child_text(xml, "taxon_id"),
    scientific_name: parse_child_text(xml, "scientific_name"),
    common_name: parse_child_text(xml, "common_name"),
    authority: parse_child_text(xml, "authority"),
    taxon_rank: parse_child_text(xml, "taxon_rank"),
  }
}

///|
fn parse_seq_annotation_from(xml : String) -> SequenceAnnotation {
  SequenceAnnotation::{
    sequence_id: parse_child_text(xml, "sequence_id"),
    molecular_type: parse_child_text(xml, "molecular_type"),
    symbol: parse_child_text(xml, "symbol"),
    uri: parse_child_text(xml, "uri"),
  }
}

///|
fn parse_taxon_namespace_element(xml : String) -> TaxonNamespace {
  TaxonNamespace::{
    namespace_id: parse_child_text(xml, "namespace_id"),
    name: parse_child_text(xml, "name"),
  }
}

// ============================================================================
// XML Writer
// ============================================================================

///|
/// Serialize a PhyloXMLResult to an XML string.
pub fn phyloxml_to_xml(result : PhyloXMLResult) -> String {
  let buf = StringBuilder::new()
  buf.write_string("<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n")
  buf.write_string("<phyloxml xmlns=\"http://www.phyloxml.org/1.0\">\n")

  for ns in result.taxon_namespaces {
    buf.write_string("  <taxon_namespace>\n")
    buf.write_string("    <namespace_id>")
    buf.write_string(ns.namespace_id)
    buf.write_string("</namespace_id>\n")
    buf.write_string("    <name>")
    buf.write_string(ns.name)
    buf.write_string("</name>\n")
    buf.write_string("  </taxon_namespace>\n")
  }

  for tree in result.trees {
    let tree_xml = tree_to_xml(tree)
    buf.write_string(tree_xml)
  }

  buf.write_string("</phyloxml>\n")
  buf.to_string()
}

///|
/// Serialize a single PhyloXMLTree to XML.
pub fn tree_to_xml(tree : PhyloXMLTree) -> String {
  let buf = StringBuilder::new()
  let indent = "  "
  buf.write_string(indent)
  buf.write_string("<phylogeny")
  if tree.tree_id.length() > 0 {
    buf.write_string(" id=\"")
    buf.write_string(tree.tree_id)
    buf.write_string("\"")
  }
  let is_rooted = tree.phylogeny_metadata.get("rooted").unwrap_or("")
  if is_rooted == "true" {
    buf.write_string(" rooted=\"true\"")
  }
  buf.write_string(">\n")

  if tree.name.length() > 0 {
    buf.write_string(indent)
    buf.write_string("  <name>")
    buf.write_string(tree.name)
    buf.write_string("</name>\n")
  }
  if tree.description.length() > 0 {
    buf.write_string(indent)
    buf.write_string("  <description>")
    buf.write_string(tree.description)
    buf.write_string("</description>\n")
  }

  let node_xml = node_to_xml(tree.root, 2)
  buf.write_string(node_xml)

  buf.write_string(indent)
  buf.write_string("</phylogeny>\n")
  buf.to_string()
}

///|
/// Serialize a PhyloXMLNode to XML recursively with indentation.
pub fn node_to_xml(node : PhyloXMLNode, indent : Int) -> String {
  let buf = StringBuilder::new()
  let prefix = make_indent(indent)
  let prefix2 = make_indent(indent + 1)

  buf.write_string(prefix)
  buf.write_string("<clade>\n")

  if node.name.length() > 0 {
    buf.write_string(prefix2)
    buf.write_string("<name>")
    buf.write_string(node.name)
    buf.write_string("</name>\n")
  }

  if node.branch_length != 0.0 {
    buf.write_string(prefix2)
    buf.write_string("<branch_length>")
    buf.write_string(node.branch_length.to_string())
    buf.write_string("</branch_length>\n")
  }

  match node.taxon {
    Some(t) => {
      buf.write_string(prefix2)
      buf.write_string("<taxon>\n")
      if t.taxon_id.length() > 0 {
        buf.write_string(prefix2)
        buf.write_string("  <taxon_id>")
        buf.write_string(t.taxon_id)
        buf.write_string("</taxon_id>\n")
      }
      if t.scientific_name.length() > 0 {
        buf.write_string(prefix2)
        buf.write_string("  <scientific_name>")
        buf.write_string(t.scientific_name)
        buf.write_string("</scientific_name>\n")
      }
      if t.common_name.length() > 0 {
        buf.write_string(prefix2)
        buf.write_string("  <common_name>")
        buf.write_string(t.common_name)
        buf.write_string("</common_name>\n")
      }
      if t.authority.length() > 0 {
        buf.write_string(prefix2)
        buf.write_string("  <authority>")
        buf.write_string(t.authority)
        buf.write_string("</authority>\n")
      }
      if t.taxon_rank.length() > 0 {
        buf.write_string(prefix2)
        buf.write_string("  <taxon_rank>")
        buf.write_string(t.taxon_rank)
        buf.write_string("</taxon_rank>\n")
      }
      buf.write_string(prefix2)
      buf.write_string("</taxon>\n")
    }
    None => ()
  }

  match node.sequence_annotation {
    Some(s) => {
      buf.write_string(prefix2)
      buf.write_string("<sequence>\n")
      if s.sequence_id.length() > 0 {
        buf.write_string(prefix2)
        buf.write_string("  <sequence_id>")
        buf.write_string(s.sequence_id)
        buf.write_string("</sequence_id>\n")
      }
      if s.molecular_type.length() > 0 {
        buf.write_string(prefix2)
        buf.write_string("  <molecular_type>")
        buf.write_string(s.molecular_type)
        buf.write_string("</molecular_type>\n")
      }
      if s.symbol.length() > 0 {
        buf.write_string(prefix2)
        buf.write_string("  <symbol>")
        buf.write_string(s.symbol)
        buf.write_string("</symbol>\n")
      }
      if s.uri.length() > 0 {
        buf.write_string(prefix2)
        buf.write_string("  <uri>")
        buf.write_string(s.uri)
        buf.write_string("</uri>\n")
      }
      buf.write_string(prefix2)
      buf.write_string("</sequence>\n")
    }
    None => ()
  }

  for key in node.annotations.keys() {
    let val = node.annotations[key]
    buf.write_string(prefix2)
    buf.write_string("<annotation>\n")
    buf.write_string(prefix2)
    buf.write_string("  <name>")
    buf.write_string(key)
    buf.write_string("</name>\n")
    buf.write_string(prefix2)
    buf.write_string("  <value>")
    buf.write_string(val)
    buf.write_string("</value>\n")
    buf.write_string(prefix2)
    buf.write_string("</annotation>\n")
  }

  for child in node.children {
    let child_xml = node_to_xml(child, indent + 1)
    buf.write_string(child_xml)
  }

  buf.write_string(prefix)
  buf.write_string("</clade>\n")
  buf.to_string()
}

///|
/// Add proper indentation to XML string.
pub fn indent_xml(xml : String, indent : Int) -> String {
  let prefix = make_indent(indent)
  let lines = xml.split("\n").to_array()
  let buf = StringBuilder::new()
  let mut i = 0
  while i < lines.length() {
    let line = lines[i].to_owned()
    let trimmed = trim(line)
    if trimmed.length() > 0 {
      buf.write_string(prefix)
      buf.write_string(trimmed)
    }
    if i < lines.length() - 1 {
      buf.write_string("\n")
    }
    i = i + 1
  }
  buf.to_string()
}

///|
fn make_indent(level : Int) -> String {
  let mut s = ""
  let mut i = 0
  while i < level {
    s = s + "  "
    i = i + 1
  }
  s
}

// ============================================================================
// Newick <-> PhyloXML Conversion
// ============================================================================

///|
/// Convert Newick format string to PhyloXMLResult.
pub fn newick_to_phyloxml(
  newick_text : String,
  tree_id? : String = "",
) -> PhyloXMLResult {
  let node = parse_newick_to_node(newick_text)
  let tree = PhyloXMLTree::{
    tree_id,
    name: if tree_id.length() > 0 { tree_id } else { "Newick Tree" },
    description: "Converted from Newick format",
    root: node,
    phylogeny_metadata: Map([("rooted", "true")], capacity=2),
    taxon_namespaces: [],
  }
  PhyloXMLResult::{ trees: [tree], metadata: Map([], capacity=4), taxon_namespaces: [] }
}

///|
fn parse_newick_to_node(content : String) -> PhyloXMLNode {
  let trimmed = trim(content).replace(old=";", new="")
  let n = trimmed.length()
  if n == 0 {
    return PhyloXMLNode::new()
  }

  let root = parse_newick_subtree(trimmed, 0)
  root
}

///|
fn parse_newick_subtree(s : String, start : Int) -> PhyloXMLNode {
  let trimmed = trim(s).replace(old=";", new="")
  let n = trimmed.length()
  if n == 0 {
    return PhyloXMLNode::new()
  }

  if char_at(trimmed, 0) == "(" {
    let children : Array[PhyloXMLNode] = Array::new()
    let mut depth = 1
    let mut child_start = 1
    let mut i = 1
    while i < n {
      let c = char_at(trimmed, i)
      if c == "(" {
        depth = depth + 1
      } else if c == ")" {
        depth = depth - 1
        if depth == 0 {
          let child_str = trimmed[child_start:i].to_owned()
          let child = parse_newick_subtree(child_str, 0)
          children.push(child)
          let after = trimmed[i + 1:].to_owned()
          let colon_pos = after.find(":")
          let name_str = if colon_pos.is_some() {
            after[0:colon_pos.unwrap()].to_owned()
          } else {
            let comma_pos = after.find(",")
            let close_pos = after.find(")")
            let semi_pos = after.find(";")
            let end = find_min_pos([comma_pos, close_pos, semi_pos])
            if end.is_some() { after[0:end.unwrap()].to_owned() } else { after }
          }
          let bl_str = if colon_pos.is_some() {
            let after_colon = after[colon_pos.unwrap() + 1:].to_owned()
            let end_pos = find_first_of(after_colon, [",", ")", ";"])
            if end_pos.is_some() { after_colon[0:end_pos.unwrap()].to_owned() } else { after_colon }
          } else {
            ""
          }
          let bl = match parse_double(trim(bl_str)) {
            Some(v) => v
            None => 0.0
          }
          return PhyloXMLNode::{
            node_id: "",
            name: trim(name_str),
            branch_length: bl,
            children,
            is_leaf: false,
            annotations: Map([], capacity=4),
            sequence_annotation: None,
            taxon: None,
          }
        }
      } else if c == "," && depth == 1 {
        let child_str = trimmed[child_start:i].to_owned()
        let child = parse_newick_subtree(child_str, 0)
        children.push(child)
        child_start = i + 1
      }
      i = i + 1
    }
    PhyloXMLNode::{
      node_id: "",
      name: "",
      branch_length: 0.0,
      children,
      is_leaf: false,
      annotations: Map([], capacity=4),
      sequence_annotation: None,
      taxon: None,
    }
  } else {
    let colon_pos = trimmed.find(":")
    let name = if colon_pos.is_some() {
      trim(trimmed[0:colon_pos.unwrap()].to_owned())
    } else {
      trim(trimmed)
    }
    let bl_str = if colon_pos.is_some() {
      trim(trimmed[colon_pos.unwrap() + 1:].to_owned())
    } else {
      ""
    }
    let bl = match parse_double(bl_str) {
      Some(v) => v
      None => 0.0
    }
    PhyloXMLNode::{
      node_id: "",
      name,
      branch_length: bl,
      children: [],
      is_leaf: true,
      annotations: Map([], capacity=4),
      sequence_annotation: None,
      taxon: None,
    }
  }
}

///|
fn find_min_pos(positions : Array[Int?]) -> Int? {
  let mut min_val : Int? = None
  for p in positions {
    if p.is_some() {
      let v = p.unwrap()
      if min_val.is_none() || v < min_val.unwrap() {
        min_val = Some(v)
      }
    }
  }
  min_val
}

///|
fn find_first_of(s : String, chars : Array[String]) -> Int? {
  let mut result : Int? = None
  for ch in chars {
    let pos = s.find(ch)
    if pos.is_some() {
      if result.is_none() || pos.unwrap() < result.unwrap() {
        result = pos
      }
    }
  }
  result
}

///|
/// Convert PhyloXMLResult to Newick format string.
pub fn phyloxml_to_newick(result : PhyloXMLResult) -> String {
  let buf = StringBuilder::new()
  let mut i = 0
  for tree in result.trees {
    if i > 0 {
      buf.write_string("\n")
    }
    let newick = phyloxml_node_to_newick(tree.root)
    buf.write_string(newick)
    buf.write_string(";")
    i = i + 1
  }
  buf.to_string()
}

///|
/// Recursively convert a PhyloXMLNode to Newick string.
pub fn phyloxml_node_to_newick(node : PhyloXMLNode) -> String {
  if node.children.length() == 0 {
    let mut result = node.name
    if node.branch_length != 0.0 {
      result = result + ":" + node.branch_length.to_string()
    }
    result
  } else {
    let buf = StringBuilder::new()
    buf.write_string("(")
    let mut i = 0
    for child in node.children {
      if i > 0 {
        buf.write_string(",")
      }
      buf.write_string(phyloxml_node_to_newick(child))
      i = i + 1
    }
    buf.write_string(")")
    buf.write_string(node.name)
    if node.branch_length != 0.0 {
      buf.write_string(":")
      buf.write_string(node.branch_length.to_string())
    }
    buf.to_string()
  }
}

// ============================================================================
// Tree Utility Functions
// ============================================================================

///|
/// Get all tree names from a PhyloXMLResult.
pub fn get_tree_names(result : PhyloXMLResult) -> Array[String] {
  let names : Array[String] = Array::new()
  for tree in result.trees {
    names.push(tree.name)
  }
  names
}

///|
/// Get all leaf (tip) names from a PhyloXMLTree.
pub fn get_all_tips(tree : PhyloXMLTree) -> Array[String] {
  let tips : Array[String] = Array::new()
  collect_tips(tree.root, tips)
  tips
}

///|
fn collect_tips(node : PhyloXMLNode, tips : Array[String]) -> Unit {
  if node.is_leaf {
    tips.push(node.name)
  } else {
    for child in node.children {
      collect_tips(child, tips)
    }
  }
}

///|
/// Count the total number of nodes in a PhyloXMLTree.
pub fn phyloxml_get_node_count(tree : PhyloXMLTree) -> Int {
  let mut count = 0
  let stack : Array[PhyloXMLNode] = Array::new()
  stack.push(tree.root)
  while stack.length() > 0 {
    let current = match stack.pop() {
      Some(v) => v
      None => tree.root
    }
    count = count + 1
    for child in current.children {
      stack.push(child)
    }
  }
  count
}