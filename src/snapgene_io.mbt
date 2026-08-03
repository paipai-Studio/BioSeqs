// Bio.SeqIO.SnapGeneIO - Parser and writer for the SnapGene .dna binary format,
// mirroring BioPython's Bio.SeqIO.SnapGeneIO.
//
// SnapGene .dna is a binary format used by SnapGene software for storing
// plasmid and sequence data. The file consists of a sequence of packets.
//
// Packet layout (all multi-byte integers are BIG-ENDIAN):
//   Byte 0    : packet type (1 byte)
//   Bytes 1-4 : data length (uint32 big-endian)
//   Bytes 5.. : data (length bytes)
//
// Packet types:
//   0x00 : Comments  (JSON-like: type, comments, exportVersion, ...)
//   0x01 : Features  (XML-like: feature definitions with segments/qualifiers)
//   0x02 : Sequence  (raw DNA sequence as ASCII bytes)
//   0x03 : Primers   (binary: primer definitions with sequences/templates)
//   0x04 : History   (binary: history tree nodes)
//   0x05 : Notes     (XML-like: notes/annotations)
//   0x06 : Alignable sequence
//   0x07 : Alignable sequence
//   0x08-0x0B : additional/unknown packet types
//
// The first packet in a .dna file is always a Comments packet (type 0x00),
// which contains the sequence type (DNA/RNA/protein) and other metadata.

// ============================================================================
// Data structures
// ============================================================================

///|
/// The molecular type stored in the SnapGene comments packet.
pub enum SnapgeneSeqType {
  SnapgeneDna
  SnapgeneRna
  SnapgeneProtein
  SnapgeneUnknown
} derive(Eq)

///|
/// Convert a sequence type to its string name.
pub fn SnapgeneSeqType::to_string(self : SnapgeneSeqType) -> String {
  match self {
    SnapgeneDna => "DNA"
    SnapgeneRna => "RNA"
    SnapgeneProtein => "protein"
    SnapgeneUnknown => "unknown"
  }
}

///|
/// Parse a sequence type from a string.
pub fn SnapgeneSeqType::from_string(s : String) -> SnapgeneSeqType {
  match s {
    "DNA" => SnapgeneDna
    "RNA" => SnapgeneRna
    "protein" => SnapgeneProtein
    _ => SnapgeneUnknown
  }
}

///|
/// A SnapGene feature (annotation on the sequence).
pub struct SnapgeneFeature {
  name : String
  type_ : String // e.g. "CDS", "promoter", "terminator", "rep_origin"
  direction : Int // 1 = forward, -1 = reverse, 0 = none
  segments : Array[(Int, Int)] // list of (start, end) 0-based inclusive segments
  qualifiers : Map[String, String]
}

///|
/// Construct a SnapgeneFeature.
pub fn SnapgeneFeature::new(
  name~ : String,
  type_~ : String,
  direction~ : Int,
  segments~ : Array[(Int, Int)],
  qualifiers~ : Map[String, String],
) -> SnapgeneFeature {
  SnapgeneFeature::{ name, type_, direction, segments, qualifiers }
}

///|
/// A SnapGene primer.
pub struct SnapgenePrimer {
  name : String
  sequence : String
  description : String
  bind_to : Int // 1 = binds top strand, 2 = binds bottom strand
  binds_at_position : Int // 0-based binding position
}

///|
/// Construct a SnapgenePrimer.
pub fn SnapgenePrimer::new(
  name~ : String,
  sequence~ : String,
  description~ : String,
  bind_to~ : Int,
  binds_at_position~ : Int,
) -> SnapgenePrimer {
  SnapgenePrimer::{
    name,
    sequence,
    description,
    bind_to,
    binds_at_position,
  }
}

///|
/// A parsed SnapGene .dna file.
pub struct SnapgeneFile {
  pub mut seq_type : SnapgeneSeqType
  pub mut comments : String // raw comments text
  pub mut export_version : String
  pub mut sequence : String
  features : Array[SnapgeneFeature]
  primers : Array[SnapgenePrimer]
  pub mut notes : String // raw notes text
  pub mut is_circular : Bool
  extra_packets : Array[(Int, Array[Int])] // (type, raw_bytes) for unknown types
}

///|
/// Construct an empty SnapgeneFile.
pub fn SnapgeneFile::new() -> SnapgeneFile {
  SnapgeneFile::{
    seq_type: SnapgeneDna,
    comments: "",
    export_version: "",
    sequence: "",
    features: [],
    primers: [],
    notes: "",
    is_circular: false,
    extra_packets: [],
  }
}

///|
/// Set the sequence type.
pub fn SnapgeneFile::set_seq_type(self : SnapgeneFile, t : SnapgeneSeqType) -> Unit {
  self.seq_type = t
}

///|
/// Set the comments text.
pub fn SnapgeneFile::set_comments(self : SnapgeneFile, s : String) -> Unit {
  self.comments = s
}

///|
/// Set the export version.
pub fn SnapgeneFile::set_export_version(self : SnapgeneFile, s : String) -> Unit {
  self.export_version = s
}

///|
/// Set the DNA/RNA sequence.
pub fn SnapgeneFile::set_sequence(self : SnapgeneFile, s : String) -> Unit {
  self.sequence = s
}

///|
/// Set the notes text.
pub fn SnapgeneFile::set_notes(self : SnapgeneFile, s : String) -> Unit {
  self.notes = s
}

///|
/// Set whether the sequence is circular.
pub fn SnapgeneFile::set_circular(self : SnapgeneFile, b : Bool) -> Unit {
  self.is_circular = b
}

///|
/// Add a feature to the file.
pub fn SnapgeneFile::add_feature(self : SnapgeneFile, f : SnapgeneFeature) -> Unit {
  self.features.push(f)
}

///|
/// Add a primer to the file.
pub fn SnapgeneFile::add_primer(self : SnapgeneFile, p : SnapgenePrimer) -> Unit {
  self.primers.push(p)
}

// ============================================================================
// Binary helpers (hex string <-> byte array)
// ============================================================================

///|
/// Convert an integer 0-15 to a lowercase hex digit character.
fn sg_hex_digit(n : Int) -> Char {
  if n < 10 {
    Char::from_int(48 + n) // '0'..'9'
  } else {
    Char::from_int(87 + n) // 'a'..'f'
  }
}

///|
/// Convert a hex digit character's ASCII code to its integer value 0-15.
fn sg_hex_value(c : Int) -> Int {
  if c >= 48 && c <= 57 {
    c - 48
  } else if c >= 65 && c <= 70 {
    c - 55
  } else if c >= 97 && c <= 102 {
    c - 87
  } else {
    0
  }
}

///|
/// Read a uint32 big-endian from a byte array at the given offset.
fn sg_read_u32_be(data : Array[Int], offset : Int) -> Int {
  (data[offset] & 0xFF) * 16777216 +
  (data[offset + 1] & 0xFF) * 65536 +
  (data[offset + 2] & 0xFF) * 256 +
  (data[offset + 3] & 0xFF)
}

///|
/// Write a uint32 big-endian into a byte array at the given offset.
fn sg_write_u32_be(data : Array[Int], offset : Int, value : Int) -> Unit {
  data[offset] = (value >> 24) & 0xFF
  data[offset + 1] = (value >> 16) & 0xFF
  data[offset + 2] = (value >> 8) & 0xFF
  data[offset + 3] = value & 0xFF
}

///|
/// Convert a hex string to a byte array.
fn sg_hex_to_bytes(hex : String) -> Array[Int] {
  let n = hex.length()
  let bytes = Array::make(n / 2, 0)
  for i in 0..<bytes.length() {
    let high = sg_hex_value(hex.unsafe_get(i * 2).to_int())
    let low = sg_hex_value(hex.unsafe_get(i * 2 + 1).to_int())
    bytes[i] = high * 16 + low
  }
  bytes
}

///|
/// Convert a byte array to a lowercase hex string.
fn sg_bytes_to_hex(data : Array[Int]) -> String {
  let sb = StringBuilder::new()
  for i in 0..<data.length() {
    let b = data[i] & 0xFF
    sb.write_char(sg_hex_digit((b >> 4) & 0xF))
    sb.write_char(sg_hex_digit(b & 0xF))
  }
  sb.to_string()
}

///|
/// Convert ASCII byte array slice to a string.
fn sg_ascii_to_string(data : Array[Int], offset : Int, length : Int) -> String {
  let sb = StringBuilder::new()
  for i in 0..<length {
    let b = data[offset + i] & 0xFF
    sb.write_char(Char::from_int(b))
  }
  sb.to_string()
}

///|
/// Convert a string to an ASCII byte array.
fn sg_string_to_ascii(s : String) -> Array[Int] {
  let n = s.length()
  let out = Array::make(n, 0)
  for i in 0..<n {
    out[i] = s.unsafe_get(i).to_int() & 0xFF
  }
  out
}

// ============================================================================
// Comments packet parsing (JSON-like)
// ============================================================================

///|
/// Parse the comments packet JSON-like text to extract metadata.
/// Returns (seq_type_string, comments_text, export_version, is_circular).
fn sg_parse_comments(text : String) -> (SnapgeneSeqType, String, String, Bool) {
  let mut seq_type = SnapgeneDna
  let mut comments = ""
  let mut export_version = ""
  let mut is_circular = false
  // Simple JSON-like field extraction (not a full JSON parser)
  let type_val = sg_json_extract_string(text, "type")
  if type_val != "" {
    seq_type = SnapgeneSeqType::from_string(type_val)
  }
  comments = sg_json_extract_string(text, "comments")
  export_version = sg_json_extract_string(text, "exportVersion")
  // Check for "circular" boolean field
  if text.contains("\"circular\":true") || text.contains("\"circular\": true") {
    is_circular = true
  }
  (seq_type, comments, export_version, is_circular)
}

///|
/// Extract a string value for a key from simple JSON-like text.
/// Returns "" if not found.
fn sg_json_extract_string(text : String, key : String) -> String {
  let pattern = "\"" + key + "\":"
  let mut pos = text.find(pattern)
  match pos {
    None => return ""
    Some(_) => ()
  }
  pos = text.find(pattern)
  let start = match pos {
    Some(p) => p + pattern.length()
    None => return ""
  }
  // Skip whitespace
  let chars = text.to_array()
  let mut i = start
  while i < chars.length() && (chars[i] == ' ' || chars[i] == '\t') {
    i = i + 1
  }
  if i >= chars.length() || chars[i] != '"' {
    return ""
  }
  i = i + 1
  let sb = StringBuilder::new()
  while i < chars.length() && chars[i] != '"' {
    if chars[i] == '\\' && i + 1 < chars.length() {
      // Escape sequence
      i = i + 1
      match chars[i] {
        'n' => sb.write_char('\n')
        't' => sb.write_char('\t')
        'r' => sb.write_char('\r')
        '"' => sb.write_char('"')
        '\\' => sb.write_char('\\')
        '/' => sb.write_char('/')
        c => sb.write_char(c)
      }
    } else {
      sb.write_char(chars[i])
    }
    i = i + 1
  }
  sb.to_string()
}

///|
/// Build a comments packet JSON-like text from metadata.
fn sg_build_comments(
  seq_type : SnapgeneSeqType,
  comments : String,
  export_version : String,
  is_circular : Bool,
) -> String {
  let type_str = seq_type.to_string()
  let circ_str = if is_circular { "true" } else { "false" }
  let ev = if export_version == "" { "1.0" } else { export_version }
  "{\"comments\":\"" + sg_json_escape(comments) + "\",\"type\":\"" + type_str + "\",\"exportVersion\":\"" + ev + "\",\"circular\":" + circ_str + "}"
}

///|
/// Escape a string for JSON.
fn sg_json_escape(s : String) -> String {
  let sb = StringBuilder::new()
  for c in s.to_array() {
    match c {
      '"' => {
        sb.write_string("\\\"")
      }
      '\\' => {
        sb.write_string("\\\\")
      }
      '\n' => {
        sb.write_string("\\n")
      }
      '\t' => {
        sb.write_string("\\t")
      }
      '\r' => {
        sb.write_string("\\r")
      }
      c => sb.write_char(c)
    }
  }
  sb.to_string()
}

// ============================================================================
// Features packet parsing (XML-like)
// ============================================================================

///|
/// Parse the features packet XML-like text into an array of features.
fn sg_parse_features(text : String) -> Array[SnapgeneFeature] {
  let features = Array::new()
  // Features are stored in <Feature> elements
  let segments = text.split("<Feature").to_array()
  for seg in segments {
    let seg_str = seg.to_string()
    if seg_str.length() > 0 && seg_str[0:1].to_string() != "<" {
      // Try to parse this feature segment
      let feature = sg_parse_one_feature(seg_str)
      match feature {
        Some(f) => features.push(f)
        None => continue
      }
    }
  }
  features
}

///|
/// Parse a single <Feature> element text.
fn sg_parse_one_feature(text : String) -> SnapgeneFeature? {
  // Extract type attribute
  let type_ = sg_xml_attr(text, "type")
  if type_ == "" {
    return None
  }
  // Extract name attribute
  let name = sg_xml_attr(text, "name")
  // Extract directionality attribute
  let direction = match sg_xml_attr(text, "directionality") {
    "1" => 1
    "-1" => -1
    _ => 0
  }
  // Extract segments
  let segs = Array::new()
  let seg_parts = text.split("<Segment").to_array()
  for sp in seg_parts {
    let sp_str = sp.to_string()
    if sp_str.contains("range=\"") {
      let range = sg_xml_attr(sp_str, "range")
      if range != "" {
        let parsed = sg_parse_range(range)
        match parsed {
          Some(p) => segs.push(p)
          None => continue
        }
      }
    }
  }
  // Extract qualifiers from <Q> elements
  let qualifiers : Map[String, String] = Map([], capacity=0)
  let q_parts = text.split("<Q").to_array()
  for qp in q_parts {
    let qp_str = qp.to_string()
    if qp_str.contains("name=\"") && qp_str.contains("content=\"") {
      let qname = sg_xml_attr(qp_str, "name")
      let qcontent = sg_xml_attr(qp_str, "content")
      if qname != "" {
        qualifiers[qname] = qcontent
      }
    }
  }
  Some(
    SnapgeneFeature::new(
      name=name,
      type_=type_,
      direction=direction,
      segments=segs,
      qualifiers=qualifiers,
    ),
  )
}

///|
/// Extract an XML attribute value: name="value"
fn sg_xml_attr(text : String, attr : String) -> String {
  let pattern = attr + "=\""
  let pos = text.find(pattern)
  match pos {
    None => return ""
    Some(p) => {
      let start = p + pattern.length()
      let chars = text.to_array()
      let sb = StringBuilder::new()
      let mut i = start
      while i < chars.length() && chars[i] != '"' {
        sb.write_char(chars[i])
        i = i + 1
      }
      sb.to_string()
    }
  }
}

///|
/// Parse a range string "start..end" or "(start,end)" into (start, end).
fn sg_parse_range(s : String) -> (Int, Int)? {
  // SnapGene ranges can be "1234..5678" format
  let parts = s.split("..").to_array()
  if parts.length() == 2 {
    let start = parts[0].trim().to_string()
    let end = parts[1].trim().to_string()
    // Validate that both parts are numeric (non-empty digits)
    if start.length() == 0 || end.length() == 0 {
      return None
    }
    if !sg_all_digits(start) || !sg_all_digits(end) {
      return None
    }
    let sv = parse_int(start)
    let ev = parse_int(end)
    return Some((sv - 1, ev - 1)) // 1-based to 0-based
  }
  None
}

///|
/// Check whether a string consists only of digit characters.
fn sg_all_digits(s : String) -> Bool {
  if s.length() == 0 {
    return false
  }
  for c in s.to_array() {
    if c < '0' || c > '9' {
      return false
    }
  }
  true
}

///|
/// Build the features packet XML-like text from an array of features.
fn sg_build_features(features : Array[SnapgeneFeature]) -> String {
  let sb = StringBuilder::new()
  sb.write_string("<Features>")
  for f in features {
    let dir_str = if f.direction == 1 { "1" } else if f.direction == -1 { "-1" } else { "0" }
    sb.write_string("<Feature name=\"" + sg_xml_escape(f.name) + "\" type=\"" + sg_xml_escape(f.type_) + "\" directionality=\"" + dir_str + "\">")
    for seg in f.segments {
      let (start, end) = seg
      sb.write_string("<Segment range=\"" + (start + 1).to_string() + ".." + (end + 1).to_string() + "\" />")
    }
    let keys = f.qualifiers.keys()
    for k in keys {
      sb.write_string("<Q name=\"" + sg_xml_escape(k) + "\" content=\"" + sg_xml_escape(f.qualifiers[k]) + "\" />")
    }
    sb.write_string("</Feature>")
  }
  sb.write_string("</Features>")
  sb.to_string()
}

///|
/// Escape a string for XML attribute values.
fn sg_xml_escape(s : String) -> String {
  let sb = StringBuilder::new()
  for c in s.to_array() {
    match c {
      '"' => sb.write_string("&quot;")
      '<' => sb.write_string("&lt;")
      '>' => sb.write_string("&gt;")
      '&' => sb.write_string("&amp;")
      c => sb.write_char(c)
    }
  }
  sb.to_string()
}

// ============================================================================
// Primers packet parsing (binary)
// ============================================================================

///|
/// Parse the primers packet binary data into an array of primers.
fn sg_parse_primers(data : Array[Int]) -> Array[SnapgenePrimer] {
  let primers = Array::new()
  let n = data.length()
  let mut offset = 0
  while offset + 4 <= n {
    // Primer ID (uint32, often sequential)
    let _id = sg_read_u32_be(data, offset)
    offset = offset + 4
    // Name length (uint32)
    if offset + 4 > n {
      break
    }
    let name_len = sg_read_u32_be(data, offset)
    offset = offset + 4
    if offset + name_len > n {
      break
    }
    let name = sg_ascii_to_string(data, offset, name_len)
    offset = offset + name_len
    // Description length (uint32)
    if offset + 4 > n {
      break
    }
    let desc_len = sg_read_u32_be(data, offset)
    offset = offset + 4
    if offset + desc_len > n {
      break
    }
    let description = sg_ascii_to_string(data, offset, desc_len)
    offset = offset + desc_len
    // Sequence length (uint32)
    if offset + 4 > n {
      break
    }
    let seq_len = sg_read_u32_be(data, offset)
    offset = offset + 4
    if offset + seq_len > n {
      break
    }
    let sequence = sg_ascii_to_string(data, offset, seq_len)
    offset = offset + seq_len
    // Binding info: bind_to (uint32) + binds_at_position (uint32)
    if offset + 8 > n {
      break
    }
    let bind_to = sg_read_u32_be(data, offset)
    offset = offset + 4
    let binds_at = sg_read_u32_be(data, offset)
    offset = offset + 4
    primers.push(
      SnapgenePrimer::new(
        name=name,
        sequence=sequence,
        description=description,
        bind_to=bind_to,
        binds_at_position=binds_at,
      ),
    )
  }
  primers
}

///|
/// Build the primers packet binary data from an array of primers.
fn sg_build_primers(primers : Array[SnapgenePrimer]) -> Array[Int] {
  let out = Array::new()
  for i in 0..<primers.length() {
    let p = primers[i]
    // Primer ID
    let id_bytes = Array::make(4, 0)
    sg_write_u32_be(id_bytes, 0, i + 1)
    for b in id_bytes {
      out.push(b)
    }
    // Name
    let name_bytes = sg_string_to_ascii(p.name)
    let nl = Array::make(4, 0)
    sg_write_u32_be(nl, 0, name_bytes.length())
    for b in nl {
      out.push(b)
    }
    for b in name_bytes {
      out.push(b)
    }
    // Description
    let desc_bytes = sg_string_to_ascii(p.description)
    let dl = Array::make(4, 0)
    sg_write_u32_be(dl, 0, desc_bytes.length())
    for b in dl {
      out.push(b)
    }
    for b in desc_bytes {
      out.push(b)
    }
    // Sequence
    let seq_bytes = sg_string_to_ascii(p.sequence)
    let sl = Array::make(4, 0)
    sg_write_u32_be(sl, 0, seq_bytes.length())
    for b in sl {
      out.push(b)
    }
    for b in seq_bytes {
      out.push(b)
    }
    // Binding info
    let bt = Array::make(4, 0)
    sg_write_u32_be(bt, 0, p.bind_to)
    for b in bt {
      out.push(b)
    }
    let bp = Array::make(4, 0)
    sg_write_u32_be(bp, 0, p.binds_at_position)
    for b in bp {
      out.push(b)
    }
  }
  out
}

// ============================================================================
// Main parsing and serialization
// ============================================================================

///|
/// Parse a SnapGene .dna file from a hex string.
pub fn snapgene_parse(hex : String) -> SnapgeneFile {
  let bytes = sg_hex_to_bytes(hex)
  snapgene_from_bytes(bytes)
}

///|
/// Parse a SnapGene .dna file from a byte array.
pub fn snapgene_from_bytes(data : Array[Int]) -> SnapgeneFile {
  let file = SnapgeneFile::new()
  let n = data.length()
  let mut offset = 0
  let mut first_packet = true
  while offset + 5 <= n {
    let pkt_type = data[offset] & 0xFF
    let pkt_len = sg_read_u32_be(data, offset + 1)
    offset = offset + 5
    if offset + pkt_len > n {
      break
    }
    let pkt_data_offset = offset
    match pkt_type {
      0 => {
        // Comments packet
        let text = sg_ascii_to_string(data, pkt_data_offset, pkt_len)
        let (st, comments, ev, circ) = sg_parse_comments(text)
        file.seq_type = st
        file.comments = comments
        file.export_version = ev
        file.is_circular = circ
      }
      1 => {
        // Features packet
        let text = sg_ascii_to_string(data, pkt_data_offset, pkt_len)
        let feats = sg_parse_features(text)
        for f in feats {
          file.features.push(f)
        }
      }
      2 => {
        // Sequence packet (raw DNA as ASCII)
        let seq = sg_ascii_to_string(data, pkt_data_offset, pkt_len)
        file.sequence = seq
      }
      3 => {
        // Primers packet
        let primer_bytes = Array::new()
        for i in 0..<pkt_len {
          primer_bytes.push(data[pkt_data_offset + i])
        }
        let prims = sg_parse_primers(primer_bytes)
        for p in prims {
          file.primers.push(p)
        }
      }
      5 => {
        // Notes packet
        let text = sg_ascii_to_string(data, pkt_data_offset, pkt_len)
        file.notes = text
      }
      _ => {
        // Unknown packet type — store raw bytes
        let raw = Array::new()
        for i in 0..<pkt_len {
          raw.push(data[pkt_data_offset + i])
        }
        file.extra_packets.push((pkt_type, raw))
      }
    }
    offset = offset + pkt_len
    first_packet = false
  }
  file
}

///|
/// Serialize a SnapgeneFile to a hex string.
pub fn snapgene_write(file : SnapgeneFile) -> String {
  let bytes = snapgene_to_bytes(file)
  sg_bytes_to_hex(bytes)
}

///|
/// Serialize a SnapgeneFile to a byte array.
pub fn snapgene_to_bytes(file : SnapgeneFile) -> Array[Int] {
  let out = Array::new()
  // Packet 0: Comments
  let comments_text = sg_build_comments(
    file.seq_type,
    file.comments,
    file.export_version,
    file.is_circular,
  )
  let comments_bytes = sg_string_to_ascii(comments_text)
  out.push(0)
  let cl = Array::make(4, 0)
  sg_write_u32_be(cl, 0, comments_bytes.length())
  for b in cl {
    out.push(b)
  }
  for b in comments_bytes {
    out.push(b)
  }
  // Packet 2: Sequence
  let seq_bytes = sg_string_to_ascii(file.sequence)
  out.push(2)
  let sl = Array::make(4, 0)
  sg_write_u32_be(sl, 0, seq_bytes.length())
  for b in sl {
    out.push(b)
  }
  for b in seq_bytes {
    out.push(b)
  }
  // Packet 1: Features (if any)
  if file.features.length() > 0 {
    let features_text = sg_build_features(file.features)
    let features_bytes = sg_string_to_ascii(features_text)
    out.push(1)
    let fl = Array::make(4, 0)
    sg_write_u32_be(fl, 0, features_bytes.length())
    for b in fl {
      out.push(b)
    }
    for b in features_bytes {
      out.push(b)
    }
  }
  // Packet 3: Primers (if any)
  if file.primers.length() > 0 {
    let primer_bytes = sg_build_primers(file.primers)
    out.push(3)
    let pl = Array::make(4, 0)
    sg_write_u32_be(pl, 0, primer_bytes.length())
    for b in pl {
      out.push(b)
    }
    for b in primer_bytes {
      out.push(b)
    }
  }
  // Packet 5: Notes (if any)
  if file.notes != "" {
    let notes_bytes = sg_string_to_ascii(file.notes)
    out.push(5)
    let nl = Array::make(4, 0)
    sg_write_u32_be(nl, 0, notes_bytes.length())
    for b in nl {
      out.push(b)
    }
    for b in notes_bytes {
      out.push(b)
    }
  }
  // Extra packets
  for pair in file.extra_packets {
    let (pkt_type, raw) = pair
    out.push(pkt_type)
    let el = Array::make(4, 0)
    sg_write_u32_be(el, 0, raw.length())
    for b in el {
      out.push(b)
    }
    for b in raw {
      out.push(b)
    }
  }
  out
}

// ============================================================================
// SeqRecord conversion
// ============================================================================

///|
/// Convert a SnapgeneFile to a SeqRecord.
pub fn snapgene_to_seqrecord(file : SnapgeneFile) -> SeqRecord {
  let id = if file.comments != "" { file.comments } else { "snapgene" }
  let rec = SeqRecord::new(
    Seq::new(file.sequence),
    id=id,
    name=id,
    description=file.comments,
  )
  rec
}

///|
/// Convert a SeqRecord to a SnapgeneFile (with default type DNA).
pub fn snapgene_from_seqrecord(rec : SeqRecord) -> SnapgeneFile {
  let file = SnapgeneFile::new()
  file.seq_type = SnapgeneDna
  file.comments = rec.description
  file.sequence = rec.seq.data
  file
}

// ============================================================================
// Utility functions
// ============================================================================

///|
/// Count the number of features of a given type.
pub fn snapgene_count_features_by_type(
  file : SnapgeneFile,
  type_name : String,
) -> Int {
  let mut count = 0
  for f in file.features {
    if f.type_ == type_name {
      count = count + 1
    }
  }
  count
}

///|
/// Get the total sequence length.
pub fn snapgene_seq_length(file : SnapgeneFile) -> Int {
  file.sequence.length()
}

///|
/// Generate a summary string for a SnapgeneFile.
pub fn snapgene_summary(file : SnapgeneFile) -> String {
  let sb = StringBuilder::new()
  sb.write_string("SnapGene File Summary:\n")
  sb.write_string("  Type: " + file.seq_type.to_string() + "\n")
  sb.write_string("  Circular: " + file.is_circular.to_string() + "\n")
  sb.write_string("  Export Version: " + file.export_version + "\n")
  sb.write_string("  Sequence Length: " + snapgene_seq_length(file).to_string() + " bp\n")
  sb.write_string("  Features: " + file.features.length().to_string() + "\n")
  sb.write_string("  Primers: " + file.primers.length().to_string() + "\n")
  if file.features.length() > 0 {
    sb.write_string("  Feature Types:\n")
    let type_counts : Map[String, Int] = Map([], capacity=0)
    for f in file.features {
      let t = f.type_
      match type_counts.get(t) {
        Some(c) => type_counts[t] = c + 1
        None => type_counts[t] = 1
      }
    }
    let keys = type_counts.keys()
    for k in keys {
      sb.write_string("    " + k + ": " + type_counts[k].to_string() + "\n")
    }
  }
  sb.to_string()
}

///|
/// Create a sample SnapgeneFile for testing/demonstration.
pub fn snapgene_sample() -> SnapgeneFile {
  let file = SnapgeneFile::new()
  file.seq_type = SnapgeneDna
  file.comments = "Sample plasmid"
  file.export_version = "1.0"
  file.is_circular = true
  file.sequence = "ATGGCTAGCAATGGGTCGCGATAGCGATAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGC"
  file.features.push(
    SnapgeneFeature::new(
      name="ORF1",
      type_="CDS",
      direction=1,
      segments=[(0, 8)],
      qualifiers=Map([], capacity=0),
    ),
  )
  file.features.push(
    SnapgeneFeature::new(
      name="Promoter1",
      type_="promoter",
      direction=1,
      segments=[(40, 50)],
      qualifiers=Map([], capacity=0),
    ),
  )
  file.primers.push(
    SnapgenePrimer::new(
      name="Fwd_Primer",
      sequence="ATGGCTAGCAAT",
      description="Forward primer for ORF1",
      bind_to=1,
      binds_at_position=0,
    ),
  )
  file.primers.push(
    SnapgenePrimer::new(
      name="Rev_Primer",
      sequence="GCTAGCTAGCT",
      description="Reverse primer for ORF1",
      bind_to=2,
      binds_at_position=58,
    ),
  )
  file
}
