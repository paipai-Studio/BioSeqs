///|
/// KGML - KEGG Markup Language pathway map format parser.
///
/// References:
///   - Biopython Bio.KEGG.KGML.KGML_parser
///   - https://www.kegg.jp/kegg/xml/
///
/// KGML is an XML format for representing KEGG pathway maps. A pathway
/// contains entries (genes, compounds, maps), relations (interactions),
/// and reactions (biochemical transformations).

// ============================================================================
// Data structures
// ============================================================================

///|
/// Graphics descriptor for a pathway entry.
pub struct KgmlGraphics {
  x : Double
  y : Double
  width : Double
  height : Double
  gtype : String
  name : String
  fgcolor : String
  bgcolor : String
}

///|
/// Create graphics descriptor.
pub fn KgmlGraphics::new(
  x : Double,
  y : Double,
  width : Double,
  height : Double,
  gtype : String,
  name : String,
  fgcolor : String,
  bgcolor : String,
) -> KgmlGraphics {
  { x, y, width, height, gtype, name, fgcolor, bgcolor }
}

///|
/// Get x coordinate.
pub fn KgmlGraphics::x(self : KgmlGraphics) -> Double {
  self.x
}

///|
/// Get y coordinate.
pub fn KgmlGraphics::y(self : KgmlGraphics) -> Double {
  self.y
}

///|
/// Get width.
pub fn KgmlGraphics::width(self : KgmlGraphics) -> Double {
  self.width
}

///|
/// Get height.
pub fn KgmlGraphics::height(self : KgmlGraphics) -> Double {
  self.height
}

///|
/// Get graphics type.
pub fn KgmlGraphics::gtype(self : KgmlGraphics) -> String {
  self.gtype
}

///|
/// Get display name.
pub fn KgmlGraphics::name(self : KgmlGraphics) -> String {
  self.name
}

///|
/// Get foreground color.
pub fn KgmlGraphics::fgcolor(self : KgmlGraphics) -> String {
  self.fgcolor
}

///|
/// Get background color.
pub fn KgmlGraphics::bgcolor(self : KgmlGraphics) -> String {
  self.bgcolor
}

///|
/// A pathway entry (node): gene, compound, enzyme, map, etc.
pub struct KgmlEntry {
  id : Int
  name : String
  etype : String
  mut reaction : String
  mut link : String
  mut graphics : KgmlGraphics?
  components : Array[Int]
}

///|
/// Create a pathway entry.
pub fn KgmlEntry::new(
  id : Int,
  name : String,
  etype : String,
) -> KgmlEntry {
  {
    id,
    name,
    etype,
    reaction: "",
    link: "",
    graphics: None,
    components: [],
  }
}

///|
/// Get entry ID.
pub fn KgmlEntry::id(self : KgmlEntry) -> Int {
  self.id
}

///|
/// Get entry name (KEGG identifiers).
pub fn KgmlEntry::name(self : KgmlEntry) -> String {
  self.name
}

///|
/// Get entry type.
pub fn KgmlEntry::etype(self : KgmlEntry) -> String {
  self.etype
}

///|
/// Get reaction identifiers.
pub fn KgmlEntry::reaction(self : KgmlEntry) -> String {
  self.reaction
}

///|
/// Get URL link.
pub fn KgmlEntry::link(self : KgmlEntry) -> String {
  self.link
}

///|
/// Get graphics descriptor (if any).
pub fn KgmlEntry::graphics(self : KgmlEntry) -> KgmlGraphics? {
  self.graphics
}

///|
/// Get component entry IDs.
pub fn KgmlEntry::components(self : KgmlEntry) -> Array[Int] {
  self.components
}

///|
/// A subtype of a relation (e.g. activation, inhibition).
pub struct KgmlSubType {
  name : String
  value : String
}

///|
/// Create a relation subtype.
pub fn KgmlSubType::new(name : String, value : String) -> KgmlSubType {
  { name, value }
}

///|
/// Get subtype name.
pub fn KgmlSubType::name(self : KgmlSubType) -> String {
  self.name
}

///|
/// Get subtype value.
pub fn KgmlSubType::value(self : KgmlSubType) -> String {
  self.value
}

///|
/// A relation between two pathway entries.
pub struct KgmlRelation {
  entry1 : Int
  entry2 : Int
  rtype : String
  subtypes : Array[KgmlSubType]
}

///|
/// Create a relation.
pub fn KgmlRelation::new(
  entry1 : Int,
  entry2 : Int,
  rtype : String,
) -> KgmlRelation {
  { entry1, entry2, rtype, subtypes: [] }
}

///|
/// Get source entry ID.
pub fn KgmlRelation::entry1(self : KgmlRelation) -> Int {
  self.entry1
}

///|
/// Get target entry ID.
pub fn KgmlRelation::entry2(self : KgmlRelation) -> Int {
  self.entry2
}

///|
/// Get relation type.
pub fn KgmlRelation::rtype(self : KgmlRelation) -> String {
  self.rtype
}

///|
/// Get subtypes.
pub fn KgmlRelation::subtypes(self : KgmlRelation) -> Array[KgmlSubType] {
  self.subtypes
}

///|
/// A reaction in the pathway.
pub struct KgmlReaction {
  name : String
  rtype : String
  substrates : Array[String]
  products : Array[String]
}

///|
/// Create a reaction.
pub fn KgmlReaction::new(name : String, rtype : String) -> KgmlReaction {
  { name, rtype, substrates: [], products: [] }
}

///|
/// Get reaction name.
pub fn KgmlReaction::name(self : KgmlReaction) -> String {
  self.name
}

///|
/// Get reaction type.
pub fn KgmlReaction::rtype(self : KgmlReaction) -> String {
  self.rtype
}

///|
/// Get substrate names.
pub fn KgmlReaction::substrates(self : KgmlReaction) -> Array[String] {
  self.substrates
}

///|
/// Get product names.
pub fn KgmlReaction::products(self : KgmlReaction) -> Array[String] {
  self.products
}

///|
/// A KEGG pathway map (KGML root element).
pub struct KgmlPathway {
  name : String
  org : String
  number : String
  title : String
  mut image : String
  mut link : String
  mut entries : Array[KgmlEntry]
  mut relations : Array[KgmlRelation]
  mut reactions : Array[KgmlReaction]
}

///|
/// Create a pathway.
pub fn KgmlPathway::new(
  name : String,
  org : String,
  number : String,
  title : String,
) -> KgmlPathway {
  {
    name,
    org,
    number,
    title,
    image: "",
    link: "",
    entries: [],
    relations: [],
    reactions: [],
  }
}

///|
/// Get pathway name.
pub fn KgmlPathway::name(self : KgmlPathway) -> String {
  self.name
}

///|
/// Get organism code.
pub fn KgmlPathway::org(self : KgmlPathway) -> String {
  self.org
}

///|
/// Get pathway number.
pub fn KgmlPathway::number(self : KgmlPathway) -> String {
  self.number
}

///|
/// Get pathway title.
pub fn KgmlPathway::title(self : KgmlPathway) -> String {
  self.title
}

///|
/// Get image URL.
pub fn KgmlPathway::image(self : KgmlPathway) -> String {
  self.image
}

///|
/// Get link URL.
pub fn KgmlPathway::link(self : KgmlPathway) -> String {
  self.link
}

///|
/// Get all entries.
pub fn KgmlPathway::entries(self : KgmlPathway) -> Array[KgmlEntry] {
  self.entries
}

///|
/// Get all relations.
pub fn KgmlPathway::relations(self : KgmlPathway) -> Array[KgmlRelation] {
  self.relations
}

///|
/// Get all reactions.
pub fn KgmlPathway::reactions(self : KgmlPathway) -> Array[KgmlReaction] {
  self.reactions
}

///|
/// Get entry by ID.
pub fn KgmlPathway::get_entry(self : KgmlPathway, id : Int) -> KgmlEntry? {
  for e in self.entries {
    if e.id == id { return Some(e) }
  }
  None
}

///|
/// Get entries by type (e.g. "gene", "compound", "map").
pub fn KgmlPathway::get_entries_by_type(
  self : KgmlPathway,
  etype : String,
) -> Array[KgmlEntry] {
  let result : Array[KgmlEntry] = []
  for e in self.entries {
    if e.etype == etype { result.push(e) }
  }
  result
}

///|
/// Get all relations involving a specific entry ID.
pub fn KgmlPathway::get_relations_for_entry(
  self : KgmlPathway,
  id : Int,
) -> Array[KgmlRelation] {
  let result : Array[KgmlRelation] = []
  for r in self.relations {
    if r.entry1 == id || r.entry2 == id { result.push(r) }
  }
  result
}

// ============================================================================
// String parsing utilities
// ============================================================================

///|
/// Find substring starting from a given position. Returns byte offset or None.
fn kgml_find_from(s : String, sub : String, from : Int) -> Int? {
  let n = s.length()
  let sub_len = sub.length()
  if sub_len == 0 || from >= n { return None }
  let mut i = from
  while i + sub_len <= n {
    let mut match_found = true
    let mut j = 0
    while j < sub_len {
      if s[i + j] != sub[j] {
        match_found = false
        break
      }
      j = j + 1
    }
    if match_found { return Some(i) }
    i = i + 1
  }
  None
}

///|
/// Extract an XML attribute value from a tag string.
/// Returns "" if the attribute is not found.
fn kgml_attr(tag : String, attr_name : String) -> String {
  let pattern = attr_name + "=\""
  match tag.find(pattern) {
    Some(start) => {
      let val_start = start + pattern.length()
      match kgml_find_from(tag, "\"", val_start) {
        Some(e) => tag[val_start:e].to_owned()
        None => ""
      }
    }
    None => {
      let pattern2 = attr_name + "='"
      match tag.find(pattern2) {
        Some(start) => {
          let val_start = start + pattern2.length()
          match kgml_find_from(tag, "'", val_start) {
            Some(e) => tag[val_start:e].to_owned()
            None => ""
          }
        }
        None => ""
      }
    }
  }
}

///|
/// Parse an integer from a string, returning default on failure.
fn kgml_parse_int(s : String, default : Int) -> Int {
  if s.length() == 0 { return default }
  let trimmed = s.trim().to_owned()
  let result = parse_int(trimmed)
  // parse_int returns 0 on failure; distinguish empty from real 0
  if trimmed.length() == 0 { return default }
  result
}

///|
/// Parse a double from a string, returning default on failure.
fn kgml_parse_double(s : String, default : Double) -> Double {
  if s.length() == 0 { return default }
  let trimmed = s.trim().to_owned()
  match parse_double(trimmed) {
    Some(d) => d
    None => default
  }
}

///|
/// Extract an integer XML attribute, returns default if not found.
fn kgml_attr_int(tag : String, attr_name : String, default : Int) -> Int {
  kgml_parse_int(kgml_attr(tag, attr_name), default)
}

///|
/// Extract a double XML attribute, returns default if not found.
fn kgml_attr_double(tag : String, attr_name : String, default : Double) -> Double {
  kgml_parse_double(kgml_attr(tag, attr_name), default)
}

///|
/// Check if a string starts with a prefix (ignoring leading whitespace).
fn kgml_starts_with(tag : String, prefix : String) -> Bool {
  let trimmed = tag.trim_start()
  trimmed.has_prefix(prefix)
}

///|
/// Check if a tag is self-closing (ends with '/').
fn kgml_is_self_closing(tag : String) -> Bool {
  tag.trim_end().has_suffix("/")
}

// ============================================================================
// KGML parser
// ============================================================================

///|
/// Parse a KGML XML string into a KgmlPathway.
/// Returns None if the input is not valid KGML.
pub fn parse_kgml(content : String) -> KgmlPathway? {
  let pathway_start = content.find("<pathway")
  match pathway_start {
    None => None
    Some(ps) => {
      match kgml_find_from(content, ">", ps) {
        None => None
        Some(te) => {
          let pathway_tag = content[ps:te].to_owned()
          let pw_name = kgml_attr(pathway_tag, "name")
          let pw_org = kgml_attr(pathway_tag, "org")
          let pw_number = kgml_attr(pathway_tag, "number")
          let pw_title = kgml_attr(pathway_tag, "title")
          let pw_image = kgml_attr(pathway_tag, "image")
          let pw_link = kgml_attr(pathway_tag, "link")
          let pathway = KgmlPathway::new(pw_name, pw_org, pw_number, pw_title)
          pathway.image = pw_image
          pathway.link = pw_link
          let entries : Array[KgmlEntry] = []
          let relations : Array[KgmlRelation] = []
          let reactions : Array[KgmlReaction] = []
          kgml_parse_children(content, te + 1, entries, relations, reactions)
          pathway.entries = entries
          pathway.relations = relations
          pathway.reactions = reactions
          Some(pathway)
        }
      }
    }
  }
}

///|
/// Parse child elements of <pathway>.
fn kgml_parse_children(
  content : String,
  start : Int,
  entries : Array[KgmlEntry],
  relations : Array[KgmlRelation],
  reactions : Array[KgmlReaction],
) -> Unit {
  let n = content.length()
  let mut pos = start
  while pos < n {
    match kgml_find_from(content, "<", pos) {
      None => break
      Some(lt_pos) => {
        if lt_pos + 1 < n {
          let next_ch = content[lt_pos + 1 : lt_pos + 2].to_owned()
          if next_ch == "/" || next_ch == "?" || next_ch == "!" {
            match kgml_find_from(content, ">", lt_pos) {
              None => break
              Some(g) => {
                pos = g + 1
                continue
              }
            }
          }
        }
        match kgml_find_from(content, ">", lt_pos) {
          None => break
          Some(g) => {
            let tag_content = content[lt_pos:g].to_owned()
            let self_closing = kgml_is_self_closing(tag_content)
            // Extract element name
            let space_idx = tag_content.find(" ")
            let elem_name : String = match space_idx {
              Some(sp) => tag_content[1:sp].trim().to_owned()
              None => {
                if self_closing {
                  tag_content[1:tag_content.length() - 1].trim().to_owned()
                } else {
                  tag_content[1:].trim().to_owned()
                }
              }
            }
            if elem_name == "entry" {
              let entry = kgml_parse_entry(content, lt_pos, g, self_closing)
              entries.push(entry)
            } else if elem_name == "relation" {
              let rel = kgml_parse_relation(content, lt_pos, g, self_closing)
              relations.push(rel)
            } else if elem_name == "reaction" {
              let rxn = kgml_parse_reaction(content, lt_pos, g, self_closing)
              reactions.push(rxn)
            }
            if self_closing {
              pos = g + 1
            } else {
              let close_tag = "</" + elem_name + ">"
              match kgml_find_from(content, close_tag, g + 1) {
                None => pos = g + 1
                Some(cp) => pos = cp + close_tag.length()
              }
            }
          }
        }
      }
    }
  }
}

///|
/// Parse an <entry> element.
fn kgml_parse_entry(
  content : String,
  lt_pos : Int,
  gt_pos : Int,
  self_closing : Bool,
) -> KgmlEntry {
  let tag = content[lt_pos:gt_pos].to_owned()
  let id = kgml_attr_int(tag, "id", 0)
  let name = kgml_attr(tag, "name")
  let etype = kgml_attr(tag, "type")
  let reaction = kgml_attr(tag, "reaction")
  let link = kgml_attr(tag, "link")
  let entry = KgmlEntry::new(id, name, etype)
  entry.reaction = reaction
  entry.link = link
  if !self_closing {
    let mut pos = gt_pos + 1
    let end_pos = match kgml_find_from(content, "</entry>", pos) {
      Some(cp) => cp
      None => content.length()
    }
    while pos < end_pos {
      match kgml_find_from(content, "<", pos) {
        None => break
        Some(lp) => {
          match kgml_find_from(content, ">", lp) {
            None => break
            Some(gp) => {
              let child_tag = content[lp:gp].to_owned()
              if kgml_starts_with(child_tag, "<graphics") {
                let gfx = kgml_parse_graphics(child_tag)
                entry.graphics = Some(gfx)
              } else if kgml_starts_with(child_tag, "<component") {
                let cid = kgml_attr_int(child_tag, "id", 0)
                entry.components.push(cid)
              }
              pos = gp + 1
            }
          }
        }
      }
    }
  }
  entry
}

///|
/// Parse a <graphics> element.
fn kgml_parse_graphics(tag : String) -> KgmlGraphics {
  let x = kgml_attr_double(tag, "x", 0.0)
  let y = kgml_attr_double(tag, "y", 0.0)
  let w = kgml_attr_double(tag, "width", 0.0)
  let h = kgml_attr_double(tag, "height", 0.0)
  let gtype = kgml_attr(tag, "type")
  let name = kgml_attr(tag, "name")
  let fg = kgml_attr(tag, "fgcolor")
  let bg = kgml_attr(tag, "bgcolor")
  let fg_final = if fg.length() == 0 { "#000000" } else { fg }
  let bg_final = if bg.length() == 0 { "#FFFFFF" } else { bg }
  KgmlGraphics::new(x, y, w, h, gtype, name, fg_final, bg_final)
}

///|
/// Parse a <relation> element and its <subtype> children.
fn kgml_parse_relation(
  content : String,
  lt_pos : Int,
  gt_pos : Int,
  self_closing : Bool,
) -> KgmlRelation {
  let tag = content[lt_pos:gt_pos].to_owned()
  let e1 = kgml_attr_int(tag, "entry1", 0)
  let e2 = kgml_attr_int(tag, "entry2", 0)
  let rtype = kgml_attr(tag, "type")
  let rel = KgmlRelation::new(e1, e2, rtype)
  if !self_closing {
    let mut pos = gt_pos + 1
    let end_pos = match kgml_find_from(content, "</relation>", pos) {
      Some(cp) => cp
      None => content.length()
    }
    while pos < end_pos {
      match kgml_find_from(content, "<", pos) {
        None => break
        Some(lp) => {
          match kgml_find_from(content, ">", lp) {
            None => break
            Some(gp) => {
              let child_tag = content[lp:gp].to_owned()
              if kgml_starts_with(child_tag, "<subtype") {
                let sn = kgml_attr(child_tag, "name")
                let sv = kgml_attr(child_tag, "value")
                rel.subtypes.push(KgmlSubType::new(sn, sv))
              }
              pos = gp + 1
            }
          }
        }
      }
    }
  }
  rel
}

///|
/// Parse a <reaction> element and its substrate/product children.
fn kgml_parse_reaction(
  content : String,
  lt_pos : Int,
  gt_pos : Int,
  self_closing : Bool,
) -> KgmlReaction {
  let tag = content[lt_pos:gt_pos].to_owned()
  let name = kgml_attr(tag, "name")
  let rtype = kgml_attr(tag, "type")
  let rxn = KgmlReaction::new(name, rtype)
  if !self_closing {
    let mut pos = gt_pos + 1
    let end_pos = match kgml_find_from(content, "</reaction>", pos) {
      Some(cp) => cp
      None => content.length()
    }
    while pos < end_pos {
      match kgml_find_from(content, "<", pos) {
        None => break
        Some(lp) => {
          match kgml_find_from(content, ">", lp) {
            None => break
            Some(gp) => {
              let child_tag = content[lp:gp].to_owned()
              if kgml_starts_with(child_tag, "<substrate") {
                let sname = kgml_attr(child_tag, "name")
                rxn.substrates.push(sname)
              } else if kgml_starts_with(child_tag, "<product") {
                let pname = kgml_attr(child_tag, "name")
                rxn.products.push(pname)
              }
              pos = gp + 1
            }
          }
        }
      }
    }
  }
  rxn
}

// ============================================================================
// Sample KGML data for testing
// ============================================================================

///|
/// Create a sample KGML pathway for testing.
/// Represents a simplified glycolysis pathway with 5 entries, 3 relations,
/// and 2 reactions.
pub fn kgml_sample_pathway() -> String {
  let mut s = "<?xml version=\"1.0\"?>\n"
  s = s + "<pathway name=\"path:ko00010\" org=\"ko\" number=\"00010\""
  s = s + " title=\"Glycolysis / Gluconeogenesis\""
  s = s + " image=\"https://www.kegg.jp/kegg/pathway/ko/ko00010.png\""
  s = s + " link=\"https://www.kegg.jp/kegg-bin/show_pathway?ko00010\">\n"
  // Entry 1: gene K00844
  s = s + "<entry id=\"1\" name=\"ko:K00844\" type=\"gene\" reaction=\"rn:R01786\""
  s = s + " link=\"https://www.kegg.jp/dbget-bin/get_brite?ko:K00844\">"
  s = s + "<graphics name=\"K00844\" fgcolor=\"#000000\" bgcolor=\"#BFFFBF\""
  s = s + " type=\"rect\" x=\"100\" y=\"200\" width=\"46\" height=\"17\"/>"
  s = s + "</entry>\n"
  // Entry 2: gene K01810
  s = s + "<entry id=\"2\" name=\"ko:K01810\" type=\"gene\" reaction=\"rn:R01070\""
  s = s + " link=\"https://www.kegg.jp/dbget-bin/get_brite?ko:K01810\">"
  s = s + "<graphics name=\"K01810\" fgcolor=\"#000000\" bgcolor=\"#BFFFBF\""
  s = s + " type=\"rect\" x=\"200\" y=\"200\" width=\"46\" height=\"17\"/>"
  s = s + "</entry>\n"
  // Entry 3: compound C00118
  s = s + "<entry id=\"3\" name=\"cpd:C00118\" type=\"compound\""
  s = s + " link=\"https://www.kegg.jp/dbget-bin/www_bget?cpd:C00118\">"
  s = s + "<graphics name=\"C00118\" fgcolor=\"#000000\" bgcolor=\"#FFFFFF\""
  s = s + " type=\"circle\" x=\"150\" y=\"250\" width=\"8\" height=\"8\"/>"
  s = s + "</entry>\n"
  // Entry 4: gene K01623
  s = s + "<entry id=\"4\" name=\"ko:K01623\" type=\"gene\" reaction=\"rn:R00200\""
  s = s + " link=\"https://www.kegg.jp/dbget-bin/get_brite?ko:K01623\">"
  s = s + "<graphics name=\"K01623\" fgcolor=\"#000000\" bgcolor=\"#BFFFBF\""
  s = s + " type=\"rect\" x=\"300\" y=\"200\" width=\"46\" height=\"17\"/>"
  s = s + "</entry>\n"
  // Entry 5: map
  s = s + "<entry id=\"5\" name=\"path:ko00020\" type=\"map\""
  s = s + " link=\"https://www.kegg.jp/kegg-bin/show_pathway?ko00020\">"
  s = s + "<graphics name=\"TCA cycle\" fgcolor=\"#000000\" bgcolor=\"#FFFFFF\""
  s = s + " type=\"roundrect\" x=\"400\" y=\"200\" width=\"46\" height=\"17\"/>"
  s = s + "</entry>\n"
  // Relations
  s = s + "<relation entry1=\"1\" entry2=\"2\" type=\"ECrel\">"
  s = s + "<subtype name=\"compound\" value=\"C00118\"/>"
  s = s + "</relation>\n"
  s = s + "<relation entry1=\"2\" entry2=\"4\" type=\"ECrel\">"
  s = s + "<subtype name=\"compound\" value=\"C00118\"/>"
  s = s + "</relation>\n"
  s = s + "<relation entry1=\"4\" entry2=\"5\" type=\"maplink\">"
  s = s + "<subtype name=\"map\" value=\"00020\"/>"
  s = s + "</relation>\n"
  // Reactions
  s = s + "<reaction name=\"rn:R01786\" type=\"irreversible\">"
  s = s + "<substrate name=\"cpd:C00118\"/>"
  s = s + "<product name=\"cpd:C00197\"/>"
  s = s + "</reaction>\n"
  s = s + "<reaction name=\"rn:R01070\" type=\"reversible\">"
  s = s + "<substrate name=\"cpd:C00197\"/>"
  s = s + "<substrate name=\"cpd:C00031\"/>"
  s = s + "<product name=\"cpd:C00118\"/>"
  s = s + "<product name=\"cpd:C00022\"/>"
  s = s + "</reaction>\n"
  s = s + "</pathway>"
  s
}

// ============================================================================
// Serialization (to KGML string)
// ============================================================================

///|
/// Serialize a KgmlPathway back to KGML XML string.
pub fn kgml_to_string(pathway : KgmlPathway) -> String {
  let mut result = "<?xml version=\"1.0\"?>\n"
  result = result +
    "<pathway name=\"" +
    pathway.name +
    "\" org=\"" +
    pathway.org +
    "\" number=\"" +
    pathway.number +
    "\" title=\"" +
    pathway.title +
    "\""
  if pathway.image.length() > 0 {
    result = result + " image=\"" + pathway.image + "\""
  }
  if pathway.link.length() > 0 {
    result = result + " link=\"" + pathway.link + "\""
  }
  result = result + ">\n"
  // Entries
  for e in pathway.entries {
    result = result + "  <entry id=\"" + e.id.to_string() + "\""
    result = result + " name=\"" + e.name + "\""
    result = result + " type=\"" + e.etype + "\""
    if e.reaction.length() > 0 {
      result = result + " reaction=\"" + e.reaction + "\""
    }
    if e.link.length() > 0 {
      result = result + " link=\"" + e.link + "\""
    }
    match e.graphics {
      Some(g) => {
        result = result + ">\n"
        result = result +
          "    <graphics name=\"" +
          g.name +
          "\" fgcolor=\"" +
          g.fgcolor +
          "\" bgcolor=\"" +
          g.bgcolor +
          "\""
        result = result + " type=\"" + g.gtype + "\""
        result = result +
          " x=\"" +
          g.x.to_string() +
          "\" y=\"" +
          g.y.to_string() +
          "\""
        result = result +
          " width=\"" +
          g.width.to_string() +
          "\" height=\"" +
          g.height.to_string() +
          "\"/>\n"
        result = result + "  </entry>\n"
      }
      None => result = result + "/>\n"
    }
  }
  // Relations
  for r in pathway.relations {
    result = result +
      "  <relation entry1=\"" +
      r.entry1.to_string() +
      "\" entry2=\"" +
      r.entry2.to_string() +
      "\""
    result = result + " type=\"" + r.rtype + "\""
    if r.subtypes.length() == 0 {
      result = result + "/>\n"
    } else {
      result = result + ">\n"
      for st in r.subtypes {
        result = result +
          "    <subtype name=\"" +
          st.name +
          "\" value=\"" +
          st.value +
          "\"/>\n"
      }
      result = result + "  </relation>\n"
    }
  }
  // Reactions
  for rxn in pathway.reactions {
    result = result + "  <reaction name=\"" + rxn.name + "\""
    result = result + " type=\"" + rxn.rtype + "\""
    if rxn.substrates.length() == 0 && rxn.products.length() == 0 {
      result = result + "/>\n"
    } else {
      result = result + ">\n"
      for s in rxn.substrates {
        result = result + "    <substrate name=\"" + s + "\"/>\n"
      }
      for p in rxn.products {
        result = result + "    <product name=\"" + p + "\"/>\n"
      }
      result = result + "  </reaction>\n"
    }
  }
  result = result + "</pathway>\n"
  result
}
