// Bio.SeqIO.GfaIO - Parser and writer for the GFA (Graph Fragment Assembly)
// format, mirroring BioPython's Bio.SeqIO.GfaIO.
//
// GFA is a tab-separated text format for genome assembly graphs. Each line
// starts with a record-type character:
//   H  Header        e.g.  H\tVN:Z:1.0
//   S  Segment       e.g.  S\tseq1\tACGTACGT\tLN:i:8
//   L  Link          e.g.  L\tseq1\t+\tseq2\t-\t8M
//   C  Containment   e.g.  C\tseq1\t+\tseq2\t+\t5\t8M
//   P  Path          e.g.  P\tpath1\tseq1+,seq2-\t8M,4M
//   W  Walk (GFA2)   e.g.  W\twalk1\tseq1\t0\t1\t0\t8
//
// Optional tags have the form  tag:type:value  where type is one of
// Z (string), i (int), f (float), J (json), H (hex), B (byte array).

// ============================================================================
// GfaTag - an optional tag on a GFA record
// ============================================================================

///|
/// An optional tag attached to a GFA record, of the form `name:type:value`
/// (e.g. `LN:i:100`).
pub struct GfaTag {
  name : String
  tag_type : String
  value : String
}

///|
/// Construct a GfaTag from its name, type, and value.
pub fn GfaTag::new(name : String, tag_type : String, value : String) -> GfaTag {
  GfaTag::{ name, tag_type, value }
}

///|
/// Tag name (e.g. "LN").
pub fn GfaTag::name(self : GfaTag) -> String {
  self.name
}

///|
/// Tag type code (e.g. "i", "Z", "f").
pub fn GfaTag::tag_type(self : GfaTag) -> String {
  self.tag_type
}

///|
/// Tag value as a raw string.
pub fn GfaTag::value(self : GfaTag) -> String {
  self.value
}

///|
/// Render the tag as `name:type:value`.
pub fn GfaTag::to_string(self : GfaTag) -> String {
  self.name + ":" + self.tag_type + ":" + self.value
}

// ============================================================================
// GfaHeader - an H record
// ============================================================================

///|
/// A GFA header record (H line). Stores the format version (from the VN tag)
/// and any additional tags.
pub struct GfaHeader {
  version : String
  tags : Array[GfaTag]
}

///|
/// Construct a GfaHeader with the given version and tags.
pub fn GfaHeader::new(version : String, tags : Array[GfaTag]) -> GfaHeader {
  GfaHeader::{ version, tags }
}

///|
/// Format version string (e.g. "1.0"), or empty if unspecified.
pub fn GfaHeader::version(self : GfaHeader) -> String {
  self.version
}

///|
/// All tags on this header (including VN if present).
pub fn GfaHeader::tags(self : GfaHeader) -> Array[GfaTag] {
  self.tags.copy()
}

// ============================================================================
// GfaSegment - an S record
// ============================================================================

///|
/// A GFA segment (S line): a named sequence with optional tags.
/// `sequence` is "*" when the sequence is not stored in the graph.
pub struct GfaSegment {
  name : String
  sequence : String
  tags : Array[GfaTag]
}

///|
/// Construct a GfaSegment with the given name, sequence, and tags.
pub fn GfaSegment::new(
  name : String,
  sequence : String,
  tags : Array[GfaTag],
) -> GfaSegment {
  GfaSegment::{ name, sequence, tags }
}

///|
/// Segment name.
pub fn GfaSegment::name(self : GfaSegment) -> String {
  self.name
}

///|
/// Segment sequence string ("*" if absent).
pub fn GfaSegment::sequence(self : GfaSegment) -> String {
  self.sequence
}

///|
/// Tags attached to this segment.
pub fn GfaSegment::tags(self : GfaSegment) -> Array[GfaTag] {
  self.tags.copy()
}

///|
/// Human-readable representation of this segment.
pub fn GfaSegment::to_string(self : GfaSegment) -> String {
  let buf = StringBuilder::new()
  buf.write_string("GfaSegment(name=")
  buf.write_string(self.name)
  buf.write_string(", sequence=")
  buf.write_string(self.sequence)
  buf.write_string(", tags=[")
  for i in 0..<self.tags.length() {
    if i > 0 {
      buf.write_string(", ")
    }
    buf.write_string(self.tags[i].to_string())
  }
  buf.write_string("])")
  buf.to_string()
}

// ============================================================================
// GfaLink - an L record
// ============================================================================

///|
/// A GFA link (L line) connecting two segments with orientations and an
/// overlap CIGAR string.
pub struct GfaLink {
  from_segment : String
  from_orient : String
  to_segment : String
  to_orient : String
  overlap : String
  tags : Array[GfaTag]
}

///|
/// Construct a GfaLink.
pub fn GfaLink::new(
  from_segment : String,
  from_orient : String,
  to_segment : String,
  to_orient : String,
  overlap : String,
  tags : Array[GfaTag],
) -> GfaLink {
  GfaLink::{
    from_segment,
    from_orient,
    to_segment,
    to_orient,
    overlap,
    tags,
  }
}

///|
/// Name of the source segment.
pub fn GfaLink::from_segment(self : GfaLink) -> String {
  self.from_segment
}

///|
/// Orientation of the source segment ("+" or "-").
pub fn GfaLink::from_orient(self : GfaLink) -> String {
  self.from_orient
}

///|
/// Name of the destination segment.
pub fn GfaLink::to_segment(self : GfaLink) -> String {
  self.to_segment
}

///|
/// Orientation of the destination segment ("+" or "-").
pub fn GfaLink::to_orient(self : GfaLink) -> String {
  self.to_orient
}

///|
/// Overlap CIGAR string (e.g. "8M").
pub fn GfaLink::overlap(self : GfaLink) -> String {
  self.overlap
}

///|
/// Tags attached to this link.
pub fn GfaLink::tags(self : GfaLink) -> Array[GfaTag] {
  self.tags.copy()
}

// ============================================================================
// GfaContainment - a C record
// ============================================================================

///|
/// A GFA containment (C line): one segment is fully contained within another
/// at a given position with an overlap CIGAR.
pub struct GfaContainment {
  container : String
  container_orient : String
  contained : String
  contained_orient : String
  position : Int
  overlap : String
  tags : Array[GfaTag]
}

///|
/// Construct a GfaContainment.
pub fn GfaContainment::new(
  container : String,
  container_orient : String,
  contained : String,
  contained_orient : String,
  position : Int,
  overlap : String,
  tags : Array[GfaTag],
) -> GfaContainment {
  GfaContainment::{
    container,
    container_orient,
    contained,
    contained_orient,
    position,
    overlap,
    tags,
  }
}

///|
/// Name of the containing segment.
pub fn GfaContainment::container(self : GfaContainment) -> String {
  self.container
}

///|
/// Orientation of the containing segment.
pub fn GfaContainment::container_orient(self : GfaContainment) -> String {
  self.container_orient
}

///|
/// Name of the contained segment.
pub fn GfaContainment::contained(self : GfaContainment) -> String {
  self.contained
}

///|
/// Orientation of the contained segment.
pub fn GfaContainment::contained_orient(self : GfaContainment) -> String {
  self.contained_orient
}

///|
/// 0-based start position of the contained segment within the container.
pub fn GfaContainment::position(self : GfaContainment) -> Int {
  self.position
}

///|
/// Overlap CIGAR string.
pub fn GfaContainment::overlap(self : GfaContainment) -> String {
  self.overlap
}

///|
/// Tags attached to this containment.
pub fn GfaContainment::tags(self : GfaContainment) -> Array[GfaTag] {
  self.tags.copy()
}

// ============================================================================
// GfaPath - a P record
// ============================================================================

///|
/// A GFA path (P line): an ordered list of oriented segments with per-step
/// overlap CIGARs.
pub struct GfaPath {
  name : String
  segment_list : Array[String]
  overlaps : Array[String]
  tags : Array[GfaTag]
}

///|
/// Construct a GfaPath.
pub fn GfaPath::new(
  name : String,
  segment_list : Array[String],
  overlaps : Array[String],
  tags : Array[GfaTag],
) -> GfaPath {
  GfaPath::{ name, segment_list, overlaps, tags }
}

///|
/// Path name.
pub fn GfaPath::name(self : GfaPath) -> String {
  self.name
}

///|
/// Ordered segment list with orientations (e.g. ["seq1+", "seq2-"]).
pub fn GfaPath::segment_list(self : GfaPath) -> Array[String] {
  self.segment_list.copy()
}

///|
/// Per-step overlap CIGARs (e.g. ["8M", "4M"]).
pub fn GfaPath::overlaps(self : GfaPath) -> Array[String] {
  self.overlaps.copy()
}

///|
/// Tags attached to this path.
pub fn GfaPath::tags(self : GfaPath) -> Array[GfaTag] {
  self.tags.copy()
}

// ============================================================================
// GfaGraph - a complete GFA document
// ============================================================================

///|
/// A complete GFA graph: headers, segments, links, containments, and paths.
pub struct GfaGraph {
  headers : Array[GfaHeader]
  segments : Array[GfaSegment]
  links : Array[GfaLink]
  containments : Array[GfaContainment]
  paths : Array[GfaPath]
}

///|
/// Construct an empty GfaGraph.
pub fn GfaGraph::new() -> GfaGraph {
  GfaGraph::{
    headers: Array::new(),
    segments: Array::new(),
    links: Array::new(),
    containments: Array::new(),
    paths: Array::new(),
  }
}

///|
/// All header records.
pub fn GfaGraph::headers(self : GfaGraph) -> Array[GfaHeader] {
  self.headers.copy()
}

///|
/// All segment records.
pub fn GfaGraph::segments(self : GfaGraph) -> Array[GfaSegment] {
  self.segments.copy()
}

///|
/// All link records.
pub fn GfaGraph::links(self : GfaGraph) -> Array[GfaLink] {
  self.links.copy()
}

///|
/// All containment records.
pub fn GfaGraph::containments(self : GfaGraph) -> Array[GfaContainment] {
  self.containments.copy()
}

///|
/// All path records.
pub fn GfaGraph::paths(self : GfaGraph) -> Array[GfaPath] {
  self.paths.copy()
}

///|
/// Append a header to the graph.
pub fn GfaGraph::add_header(self : GfaGraph, header : GfaHeader) -> Unit {
  self.headers.push(header)
}

///|
/// Append a segment to the graph.
pub fn GfaGraph::add_segment(self : GfaGraph, segment : GfaSegment) -> Unit {
  self.segments.push(segment)
}

///|
/// Append a link to the graph.
pub fn GfaGraph::add_link(self : GfaGraph, link : GfaLink) -> Unit {
  self.links.push(link)
}

///|
/// Append a containment to the graph.
pub fn GfaGraph::add_containment(
  self : GfaGraph,
  cont : GfaContainment,
) -> Unit {
  self.containments.push(cont)
}

///|
/// Append a path to the graph.
pub fn GfaGraph::add_path(self : GfaGraph, path : GfaPath) -> Unit {
  self.paths.push(path)
}

// ============================================================================
// Parsing helpers
// ============================================================================

///|
/// Split `s` on `delim` and return an array of owned strings. (MoonBit's
/// `String::split` returns an iterator of StringView; this helper materializes
/// the pieces into `Array[String]`.)
fn gfa_split_strings(s : String, delim : String) -> Array[String] {
  let views = s.split(delim).to_array()
  let out : Array[String] = Array::new()
  for i in 0..<views.length() {
    out.push(views[i].to_owned())
  }
  out
}

///|
/// Parse a single tag field of the form `name:type:value` (e.g. "LN:i:100")
/// into a GfaTag. Returns None when the field does not contain two ':'.
pub fn gfa_parse_tag(field : String) -> GfaTag? {
  let first = field.find(":")
  match first {
    None => None
    Some(p1) =>
      if p1 == 0 || p1 + 1 >= field.length() {
        None
      } else {
        let name = field[0:p1].to_owned()
        let rest = field[p1 + 1:].to_owned()
        let second = rest.find(":")
        match second {
          None => None
          Some(p2) => {
            let tag_type = rest[0:p2].to_owned()
            let value = rest[p2 + 1:].to_owned()
            Some(GfaTag::{ name, tag_type, value })
          }
        }
      }
  }
}

///|
/// Parse all tag fields in `fields` starting from index `start`. Non-tag
/// fields (those lacking two ':') are silently skipped.
fn gfa_parse_tags(fields : Array[String], start : Int) -> Array[GfaTag] {
  let tags = Array::new()
  for i in start..<fields.length() {
    match gfa_parse_tag(fields[i]) {
      Some(t) => tags.push(t)
      None => ()
    }
  }
  tags
}

///|
/// Parse an H line's tab-split fields (with fields[0] == "H") into a
/// GfaHeader. The version is taken from the VN tag if present.
pub fn gfa_parse_header(fields : Array[String]) -> GfaHeader {
  let tags = gfa_parse_tags(fields, 1)
  let mut version = ""
  for i in 0..<tags.length() {
    if tags[i].name == "VN" {
      version = tags[i].value
    }
  }
  GfaHeader::{ version, tags }
}

///|
/// Parse an S line's tab-split fields (with fields[0] == "S") into a
/// GfaSegment. fields[1] is the name, fields[2] is the sequence (defaults to
/// "*" when absent), and remaining fields are tags.
pub fn gfa_parse_segment(fields : Array[String]) -> GfaSegment {
  let n = fields.length()
  let name = if n > 1 { fields[1] } else { "" }
  let sequence = if n > 2 { fields[2] } else { "*" }
  let tags = if n > 3 { gfa_parse_tags(fields, 3) } else { Array::new() }
  GfaSegment::{ name, sequence, tags }
}

///|
/// Parse an L line's tab-split fields (with fields[0] == "L") into a GfaLink.
pub fn gfa_parse_link(fields : Array[String]) -> GfaLink {
  let n = fields.length()
  let from_segment = if n > 1 { fields[1] } else { "" }
  let from_orient = if n > 2 { fields[2] } else { "+" }
  let to_segment = if n > 3 { fields[3] } else { "" }
  let to_orient = if n > 4 { fields[4] } else { "+" }
  let overlap = if n > 5 { fields[5] } else { "*" }
  let tags = if n > 6 { gfa_parse_tags(fields, 6) } else { Array::new() }
  GfaLink::{
    from_segment,
    from_orient,
    to_segment,
    to_orient,
    overlap,
    tags,
  }
}

///|
/// Parse a C line's tab-split fields (with fields[0] == "C") into a
/// GfaContainment. fields[5] is the 0-based position.
pub fn gfa_parse_containment(fields : Array[String]) -> GfaContainment {
  let n = fields.length()
  let container = if n > 1 { fields[1] } else { "" }
  let container_orient = if n > 2 { fields[2] } else { "+" }
  let contained = if n > 3 { fields[3] } else { "" }
  let contained_orient = if n > 4 { fields[4] } else { "+" }
  let position = if n > 5 { parse_int(fields[5]) } else { 0 }
  let overlap = if n > 6 { fields[6] } else { "*" }
  let tags = if n > 7 { gfa_parse_tags(fields, 7) } else { Array::new() }
  GfaContainment::{
    container,
    container_orient,
    contained,
    contained_orient,
    position,
    overlap,
    tags,
  }
}

///|
/// Parse a P line's tab-split fields (with fields[0] == "P") into a GfaPath.
/// fields[2] is a comma-separated list of oriented segments and fields[3] is
/// a comma-separated list of overlap CIGARs.
pub fn gfa_parse_path(fields : Array[String]) -> GfaPath {
  let n = fields.length()
  let name = if n > 1 { fields[1] } else { "" }
  let segment_list = if n > 2 && fields[2].length() > 0 {
    gfa_split_strings(fields[2], ",")
  } else {
    Array::new()
  }
  let overlaps = if n > 3 && fields[3].length() > 0 {
    gfa_split_strings(fields[3], ",")
  } else {
    Array::new()
  }
  let tags = if n > 4 { gfa_parse_tags(fields, 4) } else { Array::new() }
  GfaPath::{ name, segment_list, overlaps, tags }
}

///|
/// Parse a single GFA line and dispatch on its record-type character. The
/// line is parsed for validation but the resulting record is not accumulated
/// anywhere; use `gfa_parse` to build a complete graph.
pub fn gfa_parse_line(line : String) -> Unit {
  let n = line.length()
  let cr : UInt16 = '\r'.to_int().to_uint16()
  let mut end = n
  if n > 0 && line.unsafe_get(n - 1) == cr {
    end = n - 1
  }
  if end == 0 {
    return
  }
  let body = line[0:end].to_owned()
  let fields = gfa_split_strings(body, "\t")
  if fields.length() == 0 {
    return
  }
  let rtype = fields[0]
  match rtype {
    "H" => {
      let _ = gfa_parse_header(fields)
    }
    "S" => {
      let _ = gfa_parse_segment(fields)
    }
    "L" => {
      let _ = gfa_parse_link(fields)
    }
    "C" => {
      let _ = gfa_parse_containment(fields)
    }
    "P" => {
      let _ = gfa_parse_path(fields)
    }
    _ => ()
  }
}

///|
/// Parse a full GFA document into a GfaGraph. Blank lines are skipped and
/// unknown record types are ignored.
pub fn gfa_parse(content : String) -> GfaGraph {
  let graph = GfaGraph::new()
  let lines = gfa_split_strings(content, "\n")
  let cr : UInt16 = '\r'.to_int().to_uint16()
  for line in lines {
    let n = line.length()
    let mut end = n
    if n > 0 && line.unsafe_get(n - 1) == cr {
      end = n - 1
    }
    if end == 0 {
      continue
    }
    let body = line[0:end].to_owned()
    let fields = gfa_split_strings(body, "\t")
    if fields.length() == 0 {
      continue
    }
    let rtype = fields[0]
    match rtype {
      "H" => graph.headers.push(gfa_parse_header(fields))
      "S" => graph.segments.push(gfa_parse_segment(fields))
      "L" => graph.links.push(gfa_parse_link(fields))
      "C" => graph.containments.push(gfa_parse_containment(fields))
      "P" => graph.paths.push(gfa_parse_path(fields))
      _ => ()
    }
  }
  graph
}

// ============================================================================
// Writing functions
// ============================================================================

///|
/// Render a GfaTag as `name:type:value`.
pub fn gfa_tag_to_string(tag : GfaTag) -> String {
  tag.to_string()
}

///|
/// Render a GfaHeader as an H line (without trailing newline).
pub fn gfa_header_to_line(header : GfaHeader) -> String {
  let buf = StringBuilder::new()
  buf.write_string("H")
  for i in 0..<header.tags.length() {
    buf.write_string("\t")
    buf.write_string(header.tags[i].to_string())
  }
  buf.to_string()
}

///|
/// Render a GfaSegment as an S line (without trailing newline).
pub fn gfa_segment_to_line(segment : GfaSegment) -> String {
  let buf = StringBuilder::new()
  buf.write_string("S\t")
  buf.write_string(segment.name)
  buf.write_string("\t")
  buf.write_string(segment.sequence)
  for i in 0..<segment.tags.length() {
    buf.write_string("\t")
    buf.write_string(segment.tags[i].to_string())
  }
  buf.to_string()
}

///|
/// Render a GfaLink as an L line (without trailing newline).
pub fn gfa_link_to_line(link : GfaLink) -> String {
  let buf = StringBuilder::new()
  buf.write_string("L\t")
  buf.write_string(link.from_segment)
  buf.write_string("\t")
  buf.write_string(link.from_orient)
  buf.write_string("\t")
  buf.write_string(link.to_segment)
  buf.write_string("\t")
  buf.write_string(link.to_orient)
  buf.write_string("\t")
  buf.write_string(link.overlap)
  for i in 0..<link.tags.length() {
    buf.write_string("\t")
    buf.write_string(link.tags[i].to_string())
  }
  buf.to_string()
}

///|
/// Render a GfaContainment as a C line (without trailing newline).
pub fn gfa_containment_to_line(cont : GfaContainment) -> String {
  let buf = StringBuilder::new()
  buf.write_string("C\t")
  buf.write_string(cont.container)
  buf.write_string("\t")
  buf.write_string(cont.container_orient)
  buf.write_string("\t")
  buf.write_string(cont.contained)
  buf.write_string("\t")
  buf.write_string(cont.contained_orient)
  buf.write_string("\t")
  buf.write_string(cont.position.to_string())
  buf.write_string("\t")
  buf.write_string(cont.overlap)
  for i in 0..<cont.tags.length() {
    buf.write_string("\t")
    buf.write_string(cont.tags[i].to_string())
  }
  buf.to_string()
}

///|
/// Render a GfaPath as a P line (without trailing newline).
pub fn gfa_path_to_line(path : GfaPath) -> String {
  let buf = StringBuilder::new()
  buf.write_string("P\t")
  buf.write_string(path.name)
  buf.write_string("\t")
  for i in 0..<path.segment_list.length() {
    if i > 0 {
      buf.write_string(",")
    }
    buf.write_string(path.segment_list[i])
  }
  buf.write_string("\t")
  for i in 0..<path.overlaps.length() {
    if i > 0 {
      buf.write_string(",")
    }
    buf.write_string(path.overlaps[i])
  }
  for i in 0..<path.tags.length() {
    buf.write_string("\t")
    buf.write_string(path.tags[i].to_string())
  }
  buf.to_string()
}

///|
/// Serialize a complete GfaGraph to GFA text. Headers, segments, links,
/// containments, and paths are emitted in that order, one record per line.
pub fn gfa_to_string(graph : GfaGraph) -> String {
  let buf = StringBuilder::new()
  for i in 0..<graph.headers.length() {
    buf.write_string(gfa_header_to_line(graph.headers[i]))
    buf.write_string("\n")
  }
  for i in 0..<graph.segments.length() {
    buf.write_string(gfa_segment_to_line(graph.segments[i]))
    buf.write_string("\n")
  }
  for i in 0..<graph.links.length() {
    buf.write_string(gfa_link_to_line(graph.links[i]))
    buf.write_string("\n")
  }
  for i in 0..<graph.containments.length() {
    buf.write_string(gfa_containment_to_line(graph.containments[i]))
    buf.write_string("\n")
  }
  for i in 0..<graph.paths.length() {
    buf.write_string(gfa_path_to_line(graph.paths[i]))
    buf.write_string("\n")
  }
  buf.to_string()
}

// ============================================================================
// Utility functions
// ============================================================================

///|
/// Number of segments in the graph.
pub fn gfa_graph_n_segments(graph : GfaGraph) -> Int {
  graph.segments.length()
}

///|
/// Number of links in the graph.
pub fn gfa_graph_n_links(graph : GfaGraph) -> Int {
  graph.links.length()
}

///|
/// Number of paths in the graph.
pub fn gfa_graph_n_paths(graph : GfaGraph) -> Int {
  graph.paths.length()
}

///|
/// Look up a segment by name, returning `Some(seg)` if found.
pub fn gfa_get_segment(graph : GfaGraph, name : String) -> GfaSegment? {
  let mut found : GfaSegment? = None
  for i in 0..<graph.segments.length() {
    if graph.segments[i].name == name {
      found = Some(graph.segments[i])
      break
    }
  }
  found
}

///|
/// Convert all segments with real sequence data (not "*") to SeqRecord
/// objects. The segment name becomes the record id, name, and description.
pub fn gfa_get_segments_as_records(graph : GfaGraph) -> Array[SeqRecord] {
  let records = Array::new()
  for i in 0..<graph.segments.length() {
    let seg = graph.segments[i]
    let seq_str = if seg.sequence == "*" { "" } else { seg.sequence }
    records.push(
      SeqRecord::new(
        Seq::new(seq_str),
        id=seg.name,
        name=seg.name,
        description=seg.name,
      ),
    )
  }
  records
}

///|
/// Return the reverse of a link: the destination becomes the source (with its
/// orientation preserved) and vice versa. Tags are copied unchanged.
pub fn gfa_reverse_link(link : GfaLink) -> GfaLink {
  GfaLink::{
    from_segment: link.to_segment,
    from_orient: link.to_orient,
    to_segment: link.from_segment,
    to_orient: link.from_orient,
    overlap: link.overlap,
    tags: link.tags.copy(),
  }
}

///|
/// Build a short human-readable summary of the graph: counts of each record
/// type and total sequence length.
pub fn gfa_graph_summary(graph : GfaGraph) -> String {
  let mut total_len = 0
  for i in 0..<graph.segments.length() {
    let seq = graph.segments[i].sequence
    if seq != "*" {
      total_len = total_len + seq.length()
    }
  }
  let buf = StringBuilder::new()
  buf.write_string("GfaGraph(headers=")
  buf.write_string(graph.headers.length().to_string())
  buf.write_string(", segments=")
  buf.write_string(graph.segments.length().to_string())
  buf.write_string(", links=")
  buf.write_string(graph.links.length().to_string())
  buf.write_string(", containments=")
  buf.write_string(graph.containments.length().to_string())
  buf.write_string(", paths=")
  buf.write_string(graph.paths.length().to_string())
  buf.write_string(", total_sequence_length=")
  buf.write_string(total_len.to_string())
  buf.write_string(")")
  buf.to_string()
}

// ============================================================================
// Sample data
// ============================================================================

///|
/// Return a sample GfaGraph with 3 segments, 2 links, and 1 path, suitable
/// for tests and demos.
pub fn gfa_sample_graph() -> GfaGraph {
  let graph = GfaGraph::new()
  // Header
  let header_tags : Array[GfaTag] = [
    GfaTag::new("VN", "Z", "1.0"),
    GfaTag::new("OR", "Z", "sample"),
  ]
  graph.add_header(GfaHeader::new("1.0", header_tags))
  // Segments
  graph.add_segment(
    GfaSegment::new("seq1", "ACGTACGT", [GfaTag::new("LN", "i", "8")]),
  )
  graph.add_segment(
    GfaSegment::new("seq2", "TTTTGGGG", [GfaTag::new("LN", "i", "8")]),
  )
  graph.add_segment(
    GfaSegment::new("seq3", "CCCCAAAA", [GfaTag::new("LN", "i", "8")]),
  )
  // Links
  graph.add_link(GfaLink::new("seq1", "+", "seq2", "-", "4M", []))
  graph.add_link(GfaLink::new("seq2", "+", "seq3", "+", "4M", []))
  // Path
  graph.add_path(
    GfaPath::new("path1", ["seq1+", "seq2-", "seq3+"], ["4M", "4M"], []),
  )
  graph
}
