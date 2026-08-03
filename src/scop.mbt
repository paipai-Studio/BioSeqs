///|
/// SCOP - Structural Classification of Proteins database parser.
///
/// References:
///   - Biopython Bio.SCOP package
///   - Murzin AG, et al. (1995) SCOP: a structural classification
///     of proteins database. J Mol Biol 247:536-540.
///   - https://scop.berkeley.edu/
///
/// The SCOP hierarchy classifies protein structures into a tree:
///   root -> class -> fold -> superfamily -> family -> protein -> species -> domain
///
/// Three file formats are parsed:
///   - dir.cla.scop.txt: classification (sid -> hierarchy sunids)
///   - dir.des.scop.txt: description (sunid -> name/description)
///   - dir.hie.scop.txt: hierarchy (sunid -> parent/children)

// ============================================================================
// Node type constants and helpers
// ============================================================================

///|
/// SCOP node type codes (2-letter abbreviations).
/// - "ro": root
/// - "cl": class
/// - "cf": fold
/// - "sf": superfamily
/// - "fa": family
/// - "dm": protein
/// - "sp": species
/// - "px": domain (leaf)
pub fn scop_node_type_name(code : String) -> String {
  match code {
    "ro" => "root"
    "cl" => "class"
    "cf" => "fold"
    "sf" => "superfamily"
    "fa" => "family"
    "dm" => "protein"
    "sp" => "species"
    "px" => "domain"
    _ => "unknown"
  }
}

///|
/// Convert a long node type name to its 2-letter code.
pub fn scop_node_type_code(name : String) -> String {
  match name {
    "root" => "ro"
    "class" => "cl"
    "fold" => "cf"
    "superfamily" => "sf"
    "family" => "fa"
    "protein" => "dm"
    "species" => "sp"
    "domain" => "px"
    _ => "??"
  }
}

///|
/// Canonical ordering of node types from root to leaf.
/// Used for short-circuiting ascendent/descendent searches.
pub fn scop_node_type_order(code : String) -> Int {
  match code {
    "ro" => 0
    "cl" => 1
    "cf" => 2
    "sf" => 3
    "fa" => 4
    "dm" => 5
    "sp" => 6
    "px" => 7
    _ => 8
  }
}

// ============================================================================
// Residues (PDB residue range specification)
// ============================================================================

///|
/// A single chain fragment: (chainID, startResID, endResID).
/// Empty strings indicate "whole chain" or "unspecified".
pub struct ScopFragment {
  chain : String
  start : String
  end_ : String
}

///|
/// Create a chain fragment.
pub fn ScopFragment::new(
  chain : String,
  start : String,
  end_ : String,
) -> ScopFragment {
  { chain, start, end_ }
}

///|
/// Get chain identifier.
pub fn ScopFragment::chain(self : ScopFragment) -> String {
  self.chain
}

///|
/// Get start residue ID.
pub fn ScopFragment::start(self : ScopFragment) -> String {
  self.start
}

///|
/// Get end residue ID.
pub fn ScopFragment::end_(self : ScopFragment) -> String {
  self.end_
}

///|
/// A collection of PDB residues defining a SCOP domain.
/// Format examples:
///   "1bba A:10-20,B:" -> pdbid=1bba, fragments=[(A,10,20),(B,"","")]
///   "-"               -> whole sequence (empty fragments)
pub struct ScopResidues {
  pdbid : String
  fragments : Array[ScopFragment]
}

///|
/// Create a residues descriptor.
pub fn ScopResidues::new(
  pdbid : String,
  fragments : Array[ScopFragment],
) -> ScopResidues {
  { pdbid, fragments }
}

///|
/// Get PDB identifier.
pub fn ScopResidues::pdbid(self : ScopResidues) -> String {
  self.pdbid
}

///|
/// Get chain fragments.
pub fn ScopResidues::fragments(self : ScopResidues) -> Array[ScopFragment] {
  self.fragments
}

///|
/// Parse a residue specification string into a ScopResidues.
///
/// Format: "[pdbid ]chain:start-end,chain:start-end,..."
/// Special cases:
///   "-" means whole sequence (no fragments)
///   "(-)" means whole sequence with one unnamed chain
pub fn scop_parse_residues(s : String) -> ScopResidues {
  let trimmed = s.trim().to_owned()
  if trimmed == "-" || trimmed.length() == 0 {
    return ScopResidues::new("", [])
  }
  if trimmed == "(-)" {
    return ScopResidues::new("", [ScopFragment::new("", "", "")])
  }
  // Optional pdbid prefix followed by a space, then fragments
  let space_idx = trimmed.find(" ")
  let mut pdbid = ""
  let mut frag_str = ""
  match space_idx {
    Some(sp) => {
      pdbid = trimmed[0:sp].to_owned()
      frag_str = trimmed[sp + 1:].to_owned()
    }
    None => {
      // No pdbid, entire string is fragments
      pdbid = ""
      frag_str = trimmed
    }
  }
  let fragments : Array[ScopFragment] = []
  // Split by comma
  let parts = split_by_char(frag_str, ','.to_int())
  for part in parts {
    let p = part.trim().to_owned()
    if p.length() == 0 {
      continue
    }
    // Format: chain:start-end  OR  chain:
    let colon_idx = p.find(":")
    match colon_idx {
      None => {
        // No colon: treat the whole thing as a chain with no range
        fragments.push(ScopFragment::new(p, "", ""))
      }
      Some(ci) => {
        let chain = p[0:ci].to_owned()
        let range = p[ci + 1:].to_owned()
        if range.length() == 0 {
          // "chain:" -> whole chain
          fragments.push(ScopFragment::new(chain, "", ""))
        } else {
          // "chain:start-end"
          let dash_idx = range.find("-")
          match dash_idx {
            None => {
              fragments.push(ScopFragment::new(chain, range, ""))
            }
            Some(di) => {
              let start = range[0:di].to_owned()
              let end_ = range[di + 1:].to_owned()
              fragments.push(ScopFragment::new(chain, start, end_))
            }
          }
        }
      }
    }
  }
  ScopResidues::new(pdbid, fragments)
}

///|
/// Serialize a ScopResidues back to its string representation.
pub fn scop_residues_to_string(res : ScopResidues) -> String {
  if res.fragments.length() == 0 {
    return "-"
  }
  let mut s = ""
  if res.pdbid.length() > 0 {
    s = res.pdbid + " "
  }
  let mut first = true
  for frag in res.fragments {
    if !first {
      s = s + ","
    }
    first = false
    s = s + frag.chain + ":"
    if frag.start.length() > 0 || frag.end_.length() > 0 {
      s = s + frag.start + "-" + frag.end_
    }
  }
  s
}

// ============================================================================
// File record types
// ============================================================================

///|
/// A record from dir.cla.scop.txt (classification file).
/// Format: sid\tpdbid\tresidues\tsccs\tsunid\thierarchy
pub struct ClaRecord {
  sid : String
  pdbid : String
  residues : ScopResidues
  sccs : String
  sunid : Int
  hierarchy : Map[String, Int] // node type code -> sunid
}

///|
/// Create a CLA record.
pub fn ClaRecord::new(
  sid : String,
  pdbid : String,
  residues : ScopResidues,
  sccs : String,
  sunid : Int,
  hierarchy : Map[String, Int],
) -> ClaRecord {
  { sid, pdbid, residues, sccs, sunid, hierarchy }
}

///|
/// Get sid.
pub fn ClaRecord::sid(self : ClaRecord) -> String {
  self.sid
}

///|
/// Get PDB id.
pub fn ClaRecord::pdbid(self : ClaRecord) -> String {
  self.pdbid
}

///|
/// Get residues.
pub fn ClaRecord::residues(self : ClaRecord) -> ScopResidues {
  self.residues
}

///|
/// Get SCCS string.
pub fn ClaRecord::sccs(self : ClaRecord) -> String {
  self.sccs
}

///|
/// Get sunid.
pub fn ClaRecord::sunid(self : ClaRecord) -> Int {
  self.sunid
}

///|
/// Get hierarchy map.
pub fn ClaRecord::hierarchy(self : ClaRecord) -> Map[String, Int] {
  self.hierarchy
}

///|
/// A record from dir.des.scop.txt (description file).
/// Format: sunid\tnodetype\tsccs\tname\tdescription
pub struct DesRecord {
  sunid : Int
  nodetype : String
  sccs : String
  name : String
  description : String
}

///|
/// Create a DES record.
pub fn DesRecord::new(
  sunid : Int,
  nodetype : String,
  sccs : String,
  name : String,
  description : String,
) -> DesRecord {
  { sunid, nodetype, sccs, name, description }
}

///|
/// Get sunid.
pub fn DesRecord::sunid(self : DesRecord) -> Int {
  self.sunid
}

///|
/// Get node type.
pub fn DesRecord::nodetype(self : DesRecord) -> String {
  self.nodetype
}

///|
/// Get SCCS.
pub fn DesRecord::sccs(self : DesRecord) -> String {
  self.sccs
}

///|
/// Get name (sid for domains, "-" otherwise).
pub fn DesRecord::name(self : DesRecord) -> String {
  self.name
}

///|
/// Get description.
pub fn DesRecord::description(self : DesRecord) -> String {
  self.description
}

///|
/// A record from dir.hie.scop.txt (hierarchy file).
/// Format: sunid\tparent\tchildren
pub struct HieRecord {
  sunid : Int
  parent : Int // -1 for root
  children : Array[Int]
}

///|
/// Create a HIE record.
pub fn HieRecord::new(sunid : Int, parent : Int, children : Array[Int]) -> HieRecord {
  { sunid, parent, children }
}

///|
/// Get sunid.
pub fn HieRecord::sunid(self : HieRecord) -> Int {
  self.sunid
}

///|
/// Get parent sunid (-1 for root).
pub fn HieRecord::parent(self : HieRecord) -> Int {
  self.parent
}

///|
/// Get children sunids.
pub fn HieRecord::children(self : HieRecord) -> Array[Int] {
  self.children
}

// ============================================================================
// File parsers
// ============================================================================

///|
/// Split a string by tab character.
fn scop_split_tabs(s : String) -> Array[String] {
  split_by_char(s, '\t'.to_int())
}

///|
/// Parse a single line from dir.cla.scop.txt into a ClaRecord.
/// Returns None for comment lines (starting with '#') or malformed lines.
pub fn scop_parse_cla_line(line : String) -> ClaRecord? {
  if line.length() == 0 {
    return None
  }
  let first_ch = line[0:1].to_owned()
  if first_ch == "#" {
    return None
  }
  let fields = scop_split_tabs(line)
  if fields.length() < 6 {
    return None
  }
  let sid = fields[0]
  let pdbid = fields[1]
  let residues = scop_parse_residues(fields[2])
  let sccs = fields[3]
  let sunid = parse_int(fields[4].trim().to_owned())
  // hierarchy: "cl=48724,cf=48725,sf=49265,fa=49266,dm=49267,sp=49268,px=21953"
  let hierarchy : Map[String, Int] = Map([], capacity=8)
  let pairs = split_by_char(fields[5], ','.to_int())
  for pair in pairs {
    let eq_idx = pair.find("=")
    match eq_idx {
      Some(ei) => {
        let key = pair[0:ei].to_owned()
        let val_str = pair[ei + 1:].to_owned()
        let val = parse_int(val_str)
        hierarchy.set(key, val)
      }
      None => continue
    }
  }
  Some(ClaRecord::new(sid, pdbid, residues, sccs, sunid, hierarchy))
}

///|
/// Parse all CLA records from a multi-line string.
pub fn scop_parse_cla(content : String) -> Array[ClaRecord] {
  let records : Array[ClaRecord] = []
  let lines = split_by_char(content, '\n'.to_int())
  for line in lines {
    match scop_parse_cla_line(line) {
      Some(r) => records.push(r)
      None => continue
    }
  }
  records
}

///|
/// Parse a single line from dir.des.scop.txt into a DesRecord.
/// Returns None for comment lines or malformed lines.
pub fn scop_parse_des_line(line : String) -> DesRecord? {
  if line.length() == 0 {
    return None
  }
  let first_ch = line[0:1].to_owned()
  if first_ch == "#" {
    return None
  }
  let fields = scop_split_tabs(line)
  if fields.length() < 5 {
    return None
  }
  let sunid = parse_int(fields[0].trim().to_owned())
  let nodetype = fields[1].trim().to_owned()
  let sccs = fields[2].trim().to_owned()
  let name = fields[3].trim().to_owned()
  let description = fields[4].trim().to_owned()
  Some(DesRecord::new(sunid, nodetype, sccs, name, description))
}

///|
/// Parse all DES records from a multi-line string.
pub fn scop_parse_des(content : String) -> Array[DesRecord] {
  let records : Array[DesRecord] = []
  let lines = split_by_char(content, '\n'.to_int())
  for line in lines {
    match scop_parse_des_line(line) {
      Some(r) => records.push(r)
      None => continue
    }
  }
  records
}

///|
/// Parse a single line from dir.hie.scop.txt into a HieRecord.
/// Returns None for comment lines or malformed lines.
pub fn scop_parse_hie_line(line : String) -> HieRecord? {
  if line.length() == 0 {
    return None
  }
  let first_ch = line[0:1].to_owned()
  if first_ch == "#" {
    return None
  }
  let fields = scop_split_tabs(line)
  if fields.length() < 3 {
    return None
  }
  let sunid = parse_int(fields[0].trim().to_owned())
  let parent_str = fields[1].trim().to_owned()
  let parent = if parent_str == "-" { -1 } else { parse_int(parent_str) }
  let children : Array[Int] = []
  let child_str = fields[2].trim().to_owned()
  if child_str != "-" {
    let parts = split_by_char(child_str, ','.to_int())
    for p in parts {
      let p_trimmed = p.trim().to_owned()
      if p_trimmed.length() > 0 {
        children.push(parse_int(p_trimmed))
      }
    }
  }
  Some(HieRecord::new(sunid, parent, children))
}

///|
/// Parse all HIE records from a multi-line string.
pub fn scop_parse_hie(content : String) -> Array[HieRecord] {
  let records : Array[HieRecord] = []
  let lines = split_by_char(content, '\n'.to_int())
  for line in lines {
    match scop_parse_hie_line(line) {
      Some(r) => records.push(r)
      None => continue
    }
  }
  records
}

// ============================================================================
// Hierarchy node and Scop tree
// ============================================================================

///|
/// A node in the SCOP hierarchy tree.
/// For domain nodes (type "px"), `sid` and `residues` are populated.
pub struct ScopNode {
  sunid : Int
  mut parent : Int // -1 for root, else parent sunid
  mut children : Array[Int] // child sunids
  sccs : String
  nodetype : String
  description : String
  mut sid : String // domain identifier (for type "px")
  mut residues : ScopResidues? // residue spec (for type "px")
  mut pdbid : String // PDB id (for type "px")
}

///|
/// Create a SCOP node.
pub fn ScopNode::new(
  sunid : Int,
  nodetype : String,
  sccs : String,
  description : String,
) -> ScopNode {
  {
    sunid,
    parent: -1,
    children: [],
    sccs,
    nodetype,
    description,
    sid: "",
    residues: None,
    pdbid: "",
  }
}

///|
/// Get sunid.
pub fn ScopNode::sunid(self : ScopNode) -> Int {
  self.sunid
}

///|
/// Get parent sunid (-1 for root).
pub fn ScopNode::parent(self : ScopNode) -> Int {
  self.parent
}

///|
/// Get children sunids.
pub fn ScopNode::children(self : ScopNode) -> Array[Int] {
  self.children
}

///|
/// Get SCCS string.
pub fn ScopNode::sccs(self : ScopNode) -> String {
  self.sccs
}

///|
/// Get node type code.
pub fn ScopNode::nodetype(self : ScopNode) -> String {
  self.nodetype
}

///|
/// Get description.
pub fn ScopNode::description(self : ScopNode) -> String {
  self.description
}

///|
/// Get domain sid (for domain nodes).
pub fn ScopNode::sid(self : ScopNode) -> String {
  self.sid
}

///|
/// Get residues (for domain nodes).
pub fn ScopNode::residues(self : ScopNode) -> ScopResidues? {
  self.residues
}

///|
/// Get PDB id (for domain nodes).
pub fn ScopNode::pdbid(self : ScopNode) -> String {
  self.pdbid
}

///|
/// Check if this node is a domain (type "px").
pub fn ScopNode::is_domain(self : ScopNode) -> Bool {
  self.nodetype == "px"
}

///|
/// The SCOP hierarchy database.
pub struct Scop {
  root : ScopNode
  sunid_dict : Map[Int, ScopNode]
  sid_dict : Map[String, ScopNode]
  domains : Array[ScopNode]
}

///|
/// Build a Scop database from CLA, DES, and HIE record lists.
///
/// Construction algorithm (in-memory):
///   1. Create root node (sunid=0, type="ro").
///   2. Parse DES records: create a node for each, store in sunid_dict.
///      For domain nodes (type "px"), set sid from record.name.
///   3. Parse HIE records: wire parent/children relationships.
///   4. Parse CLA records: set residues, pdbid, and populate sid_dict.
pub fn Scop::new(
  cla_records : Array[ClaRecord],
  des_records : Array[DesRecord],
  hie_records : Array[HieRecord],
) -> Scop {
  let sunid_dict : Map[Int, ScopNode] = Map([], capacity=des_records.length() + 1)
  let sid_dict : Map[String, ScopNode] = Map([], capacity=cla_records.length())
  let domains : Array[ScopNode] = []
  // 1. Create root
  let root = ScopNode::new(0, "ro", "", "SCOP Root")
  sunid_dict.set(0, root)
  // 2. Parse DES records
  for rec in des_records {
    let node = ScopNode::new(rec.sunid, rec.nodetype, rec.sccs, rec.description)
    if rec.nodetype == "px" {
      node.sid = rec.name
    }
    sunid_dict.set(rec.sunid, node)
  }
  // 3. Parse HIE records (must set back into map since structs are value types)
  for rec in hie_records {
    match sunid_dict.get(rec.sunid) {
      None => continue
      Some(node) => {
        node.parent = rec.parent
        node.children = rec.children
        sunid_dict.set(rec.sunid, node)
      }
    }
  }
  // 4. Parse CLA records (must set back into map)
  for rec in cla_records {
    match sunid_dict.get(rec.sunid) {
      None => continue
      Some(node) => {
        node.residues = Some(rec.residues)
        node.pdbid = rec.pdbid
        if node.sid.length() == 0 {
          node.sid = rec.sid
        }
        sunid_dict.set(rec.sunid, node)
        if node.sid.length() > 0 {
          sid_dict.set(node.sid, node)
        }
      }
    }
  }
  // Collect domains from the fully-built map
  for rec in des_records {
    if rec.nodetype == "px" {
      match sunid_dict.get(rec.sunid) {
        Some(node) => domains.push(node)
        None => continue
      }
    }
  }
  // Get the updated root from the map
  let final_root = match sunid_dict.get(0) {
    Some(r) => r
    None => root
  }
  { root: final_root, sunid_dict, sid_dict, domains }
}

///|
/// Create an empty Scop with just a root node.
pub fn scop_empty() -> Scop {
  let root = ScopNode::new(0, "ro", "", "SCOP Root")
  let sunid_dict : Map[Int, ScopNode] = Map([], capacity=1)
  sunid_dict.set(0, root)
  let sid_dict : Map[String, ScopNode] = Map([], capacity=0)
  let domains : Array[ScopNode] = []
  { root, sunid_dict, sid_dict, domains }
}

///|
/// Get the root node.
pub fn Scop::root(self : Scop) -> ScopNode {
  self.root
}

///|
/// Get a node by sunid. Returns None if not found.
pub fn Scop::get_node_by_sunid(self : Scop, sunid : Int) -> ScopNode? {
  self.sunid_dict.get(sunid)
}

///|
/// Get a domain by sid. Returns None if not found.
pub fn Scop::get_domain_by_sid(self : Scop, sid : String) -> ScopNode? {
  self.sid_dict.get(sid)
}

///|
/// Get all domain nodes.
pub fn Scop::get_domains(self : Scop) -> Array[ScopNode] {
  self.domains
}

///|
/// Get the parent of a node. Returns None for root or unknown.
pub fn Scop::get_parent(self : Scop, node : ScopNode) -> ScopNode? {
  if node.parent < 0 {
    return None
  }
  self.sunid_dict.get(node.parent)
}

///|
/// Get the children of a node as ScopNode array.
pub fn Scop::get_children(self : Scop, node : ScopNode) -> Array[ScopNode] {
  let result : Array[ScopNode] = []
  for cid in node.children {
    match self.sunid_dict.get(cid) {
      Some(child) => result.push(child)
      None => continue
    }
  }
  result
}

///|
/// Get the ascendent of a node of the given type.
/// Walks up the parent chain until a node of the specified type is found.
/// `node_type` can be a 2-letter code (e.g. "fa") or long name (e.g. "family").
/// Returns None if not found.
pub fn Scop::get_ascendent(
  self : Scop,
  node : ScopNode,
  node_type : String,
) -> ScopNode? {
  // Normalize to 2-letter code
  let code = if node_type.length() == 2 {
    node_type
  } else {
    scop_node_type_code(node_type)
  }
  let target_order = scop_node_type_order(code)
  let mut current = node.parent
  while current >= 0 {
    match self.sunid_dict.get(current) {
      None => break
      Some(n) => {
        if n.nodetype == code {
          return Some(n)
        }
        // Short-circuit: if we've gone above the target level, stop
        if scop_node_type_order(n.nodetype) < target_order {
          break
        }
        current = n.parent
      }
    }
  }
  None
}

///|
/// Get all descendents of a node that are of the given type.
/// Walks the subtree breadth-first.
/// `node_type` can be a 2-letter code or long name.
pub fn Scop::get_descendents(
  self : Scop,
  node : ScopNode,
  node_type : String,
) -> Array[ScopNode] {
  let code = if node_type.length() == 2 {
    node_type
  } else {
    scop_node_type_code(node_type)
  }
  let target_order = scop_node_type_order(code)
  let result : Array[ScopNode] = []
  // BFS queue of sunids
  let queue : Array[Int] = []
  for cid in node.children {
    queue.push(cid)
  }
  while queue.length() > 0 {
    let sunid = queue.remove(0)
    match self.sunid_dict.get(sunid) {
      None => continue
      Some(n) => {
        if n.nodetype == code {
          result.push(n)
        }
        // Only descend further if this node is above the target level
        if scop_node_type_order(n.nodetype) < target_order {
          for cid in n.children {
            queue.push(cid)
          }
        }
      }
    }
  }
  result
}

// ============================================================================
// SCCS comparison
// ============================================================================

///|
/// Compare two SCOP concise classification strings (SCCS).
/// Format: "letter.fold.superfamily.family" (e.g. "a.1.1.2").
/// The first component is compared lexically; the rest numerically.
/// Returns: -1 if sccs1 < sccs2, 0 if equal, 1 if sccs1 > sccs2.
pub fn scop_cmp_sccs(sccs1 : String, sccs2 : String) -> Int {
  let parts1 = split_by_char(sccs1, '.'.to_int())
  let parts2 = split_by_char(sccs2, '.'.to_int())
  let n1 = parts1.length()
  let n2 = parts2.length()
  let n = if n1 < n2 { n1 } else { n2 }
  for i in 0..<n {
    let p1 = parts1[i].trim().to_owned()
    let p2 = parts2[i].trim().to_owned()
    if i == 0 {
      // First component: lexical comparison
      if p1 < p2 {
        return -1
      }
      if p1 > p2 {
        return 1
      }
    } else {
      // Numeric comparison
      let v1 = parse_int(p1)
      let v2 = parse_int(p2)
      if v1 < v2 {
        return -1
      }
      if v1 > v2 {
        return 1
      }
    }
  }
  // All compared components are equal; shorter string sorts first
  if n1 < n2 {
    return -1
  }
  if n1 > n2 {
    return 1
  }
  0
}

// ============================================================================
// Serialization (to file format strings)
// ============================================================================

///|
/// Serialize a HieRecord to a string.
pub fn scop_hie_to_string(rec : HieRecord) -> String {
  let parent_str = if rec.parent < 0 { "-" } else { rec.parent.to_string() }
  let children_str = if rec.children.length() == 0 {
    "-"
  } else {
    let mut s = ""
    let mut first = true
    for c in rec.children {
      if !first {
        s = s + ","
      }
      first = false
      s = s + c.to_string()
    }
    s
  }
  rec.sunid.to_string() + "\t" + parent_str + "\t" + children_str
}

///|
/// Serialize a DesRecord to a string.
pub fn scop_des_to_string(rec : DesRecord) -> String {
  let name_str = if rec.name.length() == 0 { "-" } else { rec.name }
  rec.sunid.to_string() +
  "\t" +
  rec.nodetype +
  "\t" +
  rec.sccs +
  "\t" +
  name_str +
  "\t" +
  rec.description
}

///|
/// Serialize a ClaRecord to a string.
pub fn scop_cla_to_string(rec : ClaRecord) -> String {
  let mut hier_str = ""
  let mut first = true
  // Output in canonical order: cl, cf, sf, fa, dm, sp, px
  let order = ["cl", "cf", "sf", "fa", "dm", "sp", "px"]
  for key in order {
    match rec.hierarchy.get(key) {
      Some(v) => {
        if !first {
          hier_str = hier_str + ","
        }
        first = false
        hier_str = hier_str + key + "=" + v.to_string()
      }
      None => continue
    }
  }
  rec.sid +
  "\t" +
  rec.pdbid +
  "\t" +
  scop_residues_to_string(rec.residues) +
  "\t" +
  rec.sccs +
  "\t" +
  rec.sunid.to_string() +
  "\t" +
  hier_str
}

///|
/// Serialize a Scop's hierarchy to a full HIE file content string.
pub fn Scop::write_hie(self : Scop) -> String {
  // Collect all sunids and sort
  let sunids : Array[Int] = []
  for k in self.sunid_dict.keys() {
    sunids.push(k)
  }
  // Simple bubble sort (ascending)
  let n = sunids.length()
  let mut i = 0
  while i < n {
    let mut j = 0
    while j < n - 1 - i {
      if sunids[j] > sunids[j + 1] {
        let tmp = sunids[j]
        sunids[j] = sunids[j + 1]
        sunids[j + 1] = tmp
      }
      j = j + 1
    }
    i = i + 1
  }
  let mut result = ""
  for sid in sunids {
    match self.sunid_dict.get(sid) {
      None => continue
      Some(node) => {
        let rec = HieRecord::new(node.sunid, node.parent, node.children)
        if result.length() > 0 {
          result = result + "\n"
        }
        result = result + scop_hie_to_string(rec)
      }
    }
  }
  result
}

///|
/// Serialize a Scop's descriptions to a full DES file content string.
pub fn Scop::write_des(self : Scop) -> String {
  let sunids : Array[Int] = []
  for k in self.sunid_dict.keys() {
    sunids.push(k)
  }
  let n = sunids.length()
  let mut i = 0
  while i < n {
    let mut j = 0
    while j < n - 1 - i {
      if sunids[j] > sunids[j + 1] {
        let tmp = sunids[j]
        sunids[j] = sunids[j + 1]
        sunids[j + 1] = tmp
      }
      j = j + 1
    }
    i = i + 1
  }
  let mut result = ""
  for sid in sunids {
    match self.sunid_dict.get(sid) {
      None => continue
      Some(node) => {
        let name = if node.nodetype == "px" { node.sid } else { "-" }
        let rec = DesRecord::new(
          node.sunid, node.nodetype, node.sccs, name, node.description,
        )
        if result.length() > 0 {
          result = result + "\n"
        }
        result = result + scop_des_to_string(rec)
      }
    }
  }
  result
}

///|
/// Serialize a Scop's classifications to a full CLA file content string.
pub fn Scop::write_cla(self : Scop) -> String {
  let mut result = ""
  for node in self.domains {
    match node.residues {
      None => continue
      Some(res) => {
        let hierarchy : Map[String, Int] = Map([], capacity=8)
        // Walk up the parent chain to build the hierarchy
        let mut current = node.parent
        while current >= 0 {
          match self.sunid_dict.get(current) {
            None => break
            Some(n) => {
              if n.nodetype != "ro" {
                hierarchy.set(n.nodetype, n.sunid)
              }
              current = n.parent
            }
          }
        }
        let rec = ClaRecord::new(
          node.sid, node.pdbid, res, node.sccs, node.sunid, hierarchy,
        )
        if result.length() > 0 {
          result = result + "\n"
        }
        result = result + scop_cla_to_string(rec)
      }
    }
  }
  result
}

// ============================================================================
// Sample data for testing
// ============================================================================

///|
/// Create a sample DES file content for testing.
/// A small SCOP subset with 1 class, 1 fold, 1 superfamily, 1 family,
/// 1 protein, 1 species, and 2 domains.
pub fn scop_sample_des() -> String {
  let mut s = ""
  s = s + "0\tro\t-\t-\tSCOP Root\n"
  s = s + "100\tcl\ta\t-\tAll alpha proteins\n"
  s = s + "200\tcf\ta.1\t-\tGlobin-like\n"
  s = s + "300\tsf\ta.1.1\t-\tGlobin-like\n"
  s = s + "400\tfa\ta.1.1.1\t-\tGlobins\n"
  s = s + "500\tdm\ta.1.1.1.1\t-\tHemoglobin, alpha-chain\n"
  s = s + "600\tsp\ta.1.1.1.1.1\t-\tHuman (Homo sapiens)\n"
  s = s + "1000\tpx\ta.1.1.1.1.1.1\td1hba_\t1hba Hemoglobin alpha chain\n"
  s = s + "1001\tpx\ta.1.1.1.1.1.2\td1hbb_\t1hba Hemoglobin beta chain\n"
  s
}

///|
/// Create a sample HIE file content for testing.
pub fn scop_sample_hie() -> String {
  let mut s = ""
  s = s + "0\t-\t100\n"
  s = s + "100\t0\t200\n"
  s = s + "200\t100\t300\n"
  s = s + "300\t200\t400\n"
  s = s + "400\t300\t500\n"
  s = s + "500\t400\t600\n"
  s = s + "600\t500\t1000,1001\n"
  s = s + "1000\t600\t-\n"
  s = s + "1001\t600\t-\n"
  s
}

///|
/// Create a sample CLA file content for testing.
pub fn scop_sample_cla() -> String {
  let mut s = ""
  s = s + "d1hba_\t1hba\tA:1-141\ta.1.1.1.1.1.1\t1000\tcl=100,cf=200,sf=300,fa=400,dm=500,sp=600,px=1000\n"
  s = s + "d1hbb_\t1hba\tB:1-146\ta.1.1.1.1.1.2\t1001\tcl=100,cf=200,sf=300,fa=400,dm=500,sp=600,px=1001\n"
  s
}
