///|
/// SAM header tag (key:value pair on a header line).
/// Ported from htslib header.c sam_hrec_tag_t.
pub struct SamHrecTag {
  key : String
  value : String
} derive(Eq)
pub extend SamHrecTag with Eq::{equal, not_equal}

///|
/// A parsed SAM header line (e.g. "@SQ\tSN:chr1\tLN:100").
/// Ported from htslib header.c sam_hrec_type_t.
pub struct SamHrecLine {
  type_str : String
  tags : Array[SamHrecTag]
} derive(Eq)
pub extend SamHrecLine with Eq::{equal, not_equal}

///|
/// Parsed SAM header records.
/// Ported from htslib header.c sam_hrecs_t.
pub struct SamHrecs {
  lines : Array[SamHrecLine]
  ref_names : Array[String]
  ref_lengths : Array[Int]
  rg_ids : Array[String]
  pg_ids : Array[String]
}

///|
/// Parse a single SAM header line into a SamHrecLine.
fn parse_header_line(line : String) -> SamHrecLine? {
  if line.length() < 3 {
    return None
  }
  if line[0] != '@' {
    return None
  }
  let type_str = line[1:3].to_owned()
  let t0 = type_str[0]
  let t1 = type_str[1]
  if !(t0 >= 'A' && t0 <= 'Z' && t1 >= 'A' && t1 <= 'Z') {
    return None
  }
  let tags : Array[SamHrecTag] = []
  let mut pos = 3
  let n = line.length()
  while pos < n {
    if line[pos] == '\t' {
      pos = pos + 1
    }
    if pos >= n {
      break
    }
    let mut end = pos
    while end < n && line[end] != '\t' {
      end = end + 1
    }
    let field = line[pos:end].to_owned()
    let mut colon = -1
    for i = 0; i < field.length(); i = i + 1 {
      if field[i] == ':' {
        colon = i
        break
      }
    }
    if colon > 0 && colon < field.length() - 1 {
      tags.push(SamHrecTag::{
        key: field[:colon].to_owned(),
        value: field[colon + 1:].to_owned(),
      })
    } else if field.length() > 0 {
      tags.push(SamHrecTag::{ key: field, value: "" })
    }
    pos = end
  }
  Some(SamHrecLine::{ type_str, tags })
}

///|
/// Rebuild indexes from lines.
fn rebuild_indexes(lines : Array[SamHrecLine]) -> SamHrecs {
  let ref_names : Array[String] = []
  let ref_lengths : Array[Int] = []
  let rg_ids : Array[String] = []
  let pg_ids : Array[String] = []
  for i = 0; i < lines.length(); i = i + 1 {
    let line = lines[i]
    if line.type_str == "SQ" {
      let mut name = ""
      let mut lenv = 0
      for j = 0; j < line.tags.length(); j = j + 1 {
        if line.tags[j].key == "SN" {
          name = line.tags[j].value
        } else if line.tags[j].key == "LN" {
          lenv = parse_int(line.tags[j].value)
        }
      }
      if name.length() > 0 {
        ref_names.push(name)
        ref_lengths.push(lenv)
      }
    } else if line.type_str == "RG" {
      for j = 0; j < line.tags.length(); j = j + 1 {
        if line.tags[j].key == "ID" {
          rg_ids.push(line.tags[j].value)
          break
        }
      }
    } else if line.type_str == "PG" {
      for j = 0; j < line.tags.length(); j = j + 1 {
        if line.tags[j].key == "ID" {
          pg_ids.push(line.tags[j].value)
          break
        }
      }
    }
  }
  SamHrecs::{ lines, ref_names, ref_lengths, rg_ids, pg_ids }
}

///|
/// Parse SAM header text into structured records.
/// Ported from htslib sam_hrecs_parse_lines().
pub fn sam_hrecs_parse(header_text : String) -> SamHrecs {
  let lines : Array[SamHrecLine] = []
  let mut pos = 0
  let n = header_text.length()
  while pos < n {
    let mut end = pos
    while end < n && header_text[end] != '\n' && header_text[end] != '\r' {
      end = end + 1
    }
    if end > pos {
      let line = header_text[pos:end].to_owned()
      match parse_header_line(line) {
        Some(hline) => lines.push(hline)
        None => ()
      }
    }
    pos = end
    if pos < n && header_text[pos] == '\r' {
      pos = pos + 1
    }
    if pos < n && header_text[pos] == '\n' {
      pos = pos + 1
    }
  }
  rebuild_indexes(lines)
}

///|
/// Rebuild a single header line to text.
/// Ported from htslib build_header_line().
pub fn SamHrecLine::to_text(self : SamHrecLine) -> String {
  let sb = StringBuilder()
  sb.write_char('@')
  sb.write_string(self.type_str)
  for i = 0; i < self.tags.length(); i = i + 1 {
    sb.write_char('\t')
    sb.write_string(self.tags[i].key)
    sb.write_char(':')
    sb.write_string(self.tags[i].value)
  }
  sb.to_string()
}

///|
/// Rebuild full header text from parsed records.
/// Ported from htslib sam_hrecs_rebuild_text().
pub fn SamHrecs::to_text(self : SamHrecs) -> String {
  let sb = StringBuilder()
  for i = 0; i < self.lines.length(); i = i + 1 {
    sb.write_string(self.lines[i].to_text())
    sb.write_char('\n')
  }
  sb.to_string()
}

///|
/// Number of reference sequences (@SQ lines).
/// Ported from htslib sam_hdr_nref().
pub fn SamHrecs::nref(self : SamHrecs) -> Int {
  self.ref_names.length()
}

///|
/// Reference name to tid (0-based index).
/// Ported from htslib sam_hdr_name2tid().
pub fn SamHrecs::name2tid(self : SamHrecs, name : String) -> Int {
  for i = 0; i < self.ref_names.length(); i = i + 1 {
    if self.ref_names[i] == name {
      return i
    }
  }
  -1
}

///|
/// tid to reference name.
/// Ported from htslib sam_hdr_tid2name().
pub fn SamHrecs::tid2name(self : SamHrecs, tid : Int) -> String? {
  if tid < 0 || tid >= self.ref_names.length() {
    return None
  }
  Some(self.ref_names[tid])
}

///|
/// tid to reference length.
/// Ported from htslib sam_hdr_tid2len().
pub fn SamHrecs::tid2len(self : SamHrecs, tid : Int) -> Int {
  if tid < 0 || tid >= self.ref_lengths.length() {
    return 0
  }
  self.ref_lengths[tid]
}

///|
/// Count lines of a given type (e.g. "SQ", "RG", "PG").
/// Ported from htslib sam_hdr_count_lines().
pub fn SamHrecs::count_lines(self : SamHrecs, type_str : String) -> Int {
  let mut count = 0
  for i = 0; i < self.lines.length(); i = i + 1 {
    if self.lines[i].type_str == type_str {
      count = count + 1
    }
  }
  count
}

///|
/// Find a line by type and ID key/value.
/// Ported from htslib sam_hdr_find_line_id().
pub fn SamHrecs::find_line_id(
  self : SamHrecs,
  type_str : String,
  id_key : String,
  id_value : String
) -> SamHrecLine? {
  for i = 0; i < self.lines.length(); i = i + 1 {
    let line = self.lines[i]
    if line.type_str != type_str {
      continue
    }
    if id_key.length() == 0 {
      return Some(line)
    }
    for j = 0; j < line.tags.length(); j = j + 1 {
      if line.tags[j].key == id_key && line.tags[j].value == id_value {
        return Some(line)
      }
    }
  }
  None
}

///|
/// Find a line by type and position (0-based among lines of that type).
/// Ported from htslib sam_hdr_find_line_pos().
pub fn SamHrecs::find_line_pos(
  self : SamHrecs,
  type_str : String,
  position : Int
) -> SamHrecLine? {
  let mut idx = 0
  for i = 0; i < self.lines.length(); i = i + 1 {
    if self.lines[i].type_str == type_str {
      if idx == position {
        return Some(self.lines[i])
      }
      idx = idx + 1
    }
  }
  None
}

///|
/// Find a tag value from a line identified by type and ID.
/// Ported from htslib sam_hdr_find_tag_id().
pub fn SamHrecs::find_tag_id(
  self : SamHrecs,
  type_str : String,
  id_key : String,
  id_value : String,
  tag_key : String
) -> String? {
  match self.find_line_id(type_str, id_key, id_value) {
    Some(line) => {
      for i = 0; i < line.tags.length(); i = i + 1 {
        if line.tags[i].key == tag_key {
          return Some(line.tags[i].value)
        }
      }
      None
    }
    None => None
  }
}

///|
/// Find a tag value from a line identified by type and position.
/// Ported from htslib sam_hdr_find_tag_pos().
pub fn SamHrecs::find_tag_pos(
  self : SamHrecs,
  type_str : String,
  position : Int,
  tag_key : String
) -> String? {
  match self.find_line_pos(type_str, position) {
    Some(line) => {
      for i = 0; i < line.tags.length(); i = i + 1 {
        if line.tags[i].key == tag_key {
          return Some(line.tags[i].value)
        }
      }
      None
    }
    None => None
  }
}

///|
/// Add a header line from text.
/// Ported from htslib sam_hdr_add_line().
pub fn SamHrecs::add_line(self : SamHrecs, line_text : String) -> SamHrecs {
  match parse_header_line(line_text) {
    Some(hline) => {
      let new_lines = []
      for i = 0; i < self.lines.length(); i = i + 1 {
        new_lines.push(self.lines[i])
      }
      new_lines.push(hline)
      rebuild_indexes(new_lines)
    }
    None => self
  }
}

///|
/// Add multiple header lines from text (newline-separated).
/// Ported from htslib sam_hdr_add_lines().
pub fn SamHrecs::add_lines(self : SamHrecs, lines_text : String) -> SamHrecs {
  let mut result = self
  let mut pos = 0
  let n = lines_text.length()
  while pos < n {
    let mut end = pos
    while end < n && lines_text[end] != '\n' && lines_text[end] != '\r' {
      end = end + 1
    }
    if end > pos {
      let line = lines_text[pos:end].to_owned()
      result = result.add_line(line)
    }
    pos = end
    if pos < n && lines_text[pos] == '\r' {
      pos = pos + 1
    }
    if pos < n && lines_text[pos] == '\n' {
      pos = pos + 1
    }
  }
  result
}

///|
/// Remove a line by type and ID key/value.
/// Ported from htslib sam_hdr_remove_line_id().
pub fn SamHrecs::remove_line_id(
  self : SamHrecs,
  type_str : String,
  id_key : String,
  id_value : String
) -> SamHrecs {
  let new_lines = []
  let mut removed = false
  for i = 0; i < self.lines.length(); i = i + 1 {
    let line = self.lines[i]
    if !removed && line.type_str == type_str {
      if id_key.length() == 0 {
        removed = true
        continue
      }
      let mut is_match = false
      for j = 0; j < line.tags.length(); j = j + 1 {
        if line.tags[j].key == id_key && line.tags[j].value == id_value {
          is_match = true
          break
        }
      }
      if is_match {
        removed = true
        continue
      }
    }
    new_lines.push(line)
  }
  if !removed {
    return self
  }
  rebuild_indexes(new_lines)
}

///|
/// Remove a line by type and position.
/// Ported from htslib sam_hdr_remove_line_pos().
pub fn SamHrecs::remove_line_pos(
  self : SamHrecs,
  type_str : String,
  position : Int
) -> SamHrecs {
  let new_lines = []
  let mut idx = 0
  let mut removed = false
  for i = 0; i < self.lines.length(); i = i + 1 {
    let line = self.lines[i]
    if !removed && line.type_str == type_str {
      if idx == position {
        removed = true
        continue
      }
      idx = idx + 1
    }
    new_lines.push(line)
  }
  if !removed {
    return self
  }
  rebuild_indexes(new_lines)
}

///|
/// Update a tag value in a line identified by type and ID.
/// If the tag doesn't exist, it's added. If new_value is empty, the tag is removed.
/// Ported from htslib sam_hdr_update_line().
pub fn SamHrecs::update_line(
  self : SamHrecs,
  type_str : String,
  id_key : String,
  id_value : String,
  tag_key : String,
  new_value : String
) -> SamHrecs {
  let new_lines = []
  let mut found = false
  for i = 0; i < self.lines.length(); i = i + 1 {
    let line = self.lines[i]
    if !found && line.type_str == type_str {
      let mut matches = id_key.length() == 0
      if id_key.length() > 0 {
        for j = 0; j < line.tags.length(); j = j + 1 {
          if line.tags[j].key == id_key && line.tags[j].value == id_value {
            matches = true
            break
          }
        }
      }
      if matches {
        found = true
        let new_tags = []
        let mut tag_exists = false
        for j = 0; j < line.tags.length(); j = j + 1 {
          if line.tags[j].key == tag_key {
            tag_exists = true
            if new_value.length() > 0 {
              new_tags.push(SamHrecTag::{ key: tag_key, value: new_value })
            }
          } else {
            new_tags.push(line.tags[j])
          }
        }
        if !tag_exists && new_value.length() > 0 {
          new_tags.push(SamHrecTag::{ key: tag_key, value: new_value })
        }
        new_lines.push(SamHrecLine::{
          type_str: line.type_str,
          tags: new_tags,
        })
        continue
      }
    }
    new_lines.push(line)
  }
  if !found {
    return self
  }
  rebuild_indexes(new_lines)
}

///|
/// Remove a tag from a line identified by type and ID.
/// Ported from htslib sam_hdr_remove_tag_id().
pub fn SamHrecs::remove_tag_id(
  self : SamHrecs,
  type_str : String,
  id_key : String,
  id_value : String,
  tag_key : String
) -> SamHrecs {
  self.update_line(type_str, id_key, id_value, tag_key, "")
}

///|
/// Remove all lines of a given type except the one matching ID.
/// Ported from htslib sam_hdr_remove_except().
pub fn SamHrecs::remove_except(
  self : SamHrecs,
  type_str : String,
  id_key : String,
  id_value : String
) -> SamHrecs {
  let new_lines = []
  for i = 0; i < self.lines.length(); i = i + 1 {
    let line = self.lines[i]
    if line.type_str == type_str {
      let mut matches = false
      for j = 0; j < line.tags.length(); j = j + 1 {
        if line.tags[j].key == id_key && line.tags[j].value == id_value {
          matches = true
          break
        }
      }
      if matches {
        new_lines.push(line)
      }
    } else {
      new_lines.push(line)
    }
  }
  rebuild_indexes(new_lines)
}

///|
/// Line index: find the position of a line with given type and key/value.
/// Returns 0-based index among lines of that type, or -1 if not found.
/// Ported from htslib sam_hdr_line_index().
pub fn SamHrecs::line_index(
  self : SamHrecs,
  type_str : String,
  id_key : String,
  id_value : String
) -> Int {
  let mut idx = 0
  for i = 0; i < self.lines.length(); i = i + 1 {
    let line = self.lines[i]
    if line.type_str == type_str {
      if id_key.length() == 0 {
        return idx
      }
      for j = 0; j < line.tags.length(); j = j + 1 {
        if line.tags[j].key == id_key && line.tags[j].value == id_value {
          return idx
        }
      }
      idx = idx + 1
    }
  }
  -1
}

///|
/// Get all line texts of a given type.
pub fn SamHrecs::lines_of_type(
  self : SamHrecs,
  type_str : String
) -> Array[String] {
  let result = []
  for i = 0; i < self.lines.length(); i = i + 1 {
    if self.lines[i].type_str == type_str {
      result.push(self.lines[i].to_text())
    }
  }
  result
}

///|
/// Get all reference names.
pub fn SamHrecs::ref_seq_names(self : SamHrecs) -> Array[String] {
  let result = []
  for i = 0; i < self.ref_names.length(); i = i + 1 {
    result.push(self.ref_names[i])
  }
  result
}

///|
/// Get all read group IDs.
pub fn SamHrecs::read_group_ids(self : SamHrecs) -> Array[String] {
  let result = []
  for i = 0; i < self.rg_ids.length(); i = i + 1 {
    result.push(self.rg_ids[i])
  }
  result
}

///|
/// Get all program IDs.
pub fn SamHrecs::program_ids(self : SamHrecs) -> Array[String] {
  let result = []
  for i = 0; i < self.pg_ids.length(); i = i + 1 {
    result.push(self.pg_ids[i])
  }
  result
}

///|
/// Check if a read group ID exists.
pub fn SamHrecs::has_read_group(self : SamHrecs, rg_id : String) -> Bool {
  for i = 0; i < self.rg_ids.length(); i = i + 1 {
    if self.rg_ids[i] == rg_id {
      return true
    }
  }
  false
}

///|
/// Check if a program ID exists.
pub fn SamHrecs::has_program(self : SamHrecs, pg_id : String) -> Bool {
  for i = 0; i < self.pg_ids.length(); i = i + 1 {
    if self.pg_ids[i] == pg_id {
      return true
    }
  }
  false
}

///|
/// Get the sort order from @HD line (SO tag).
pub fn SamHrecs::sort_order(self : SamHrecs) -> String {
  match self.find_tag_id("HD", "", "", "SO") {
    Some(v) => v
    None => "unknown"
  }
}

///|
/// Get the group order from @HD line (GO tag).
pub fn SamHrecs::group_order(self : SamHrecs) -> String {
  match self.find_tag_id("HD", "", "", "GO") {
    Some(v) => v
    None => "none"
  }
}

///|
/// Get the sequencing technology from @RG line (PL tag).
pub fn SamHrecs::platform(self : SamHrecs, rg_id : String) -> String {
  match self.find_tag_id("RG", "ID", rg_id, "PL") {
    Some(v) => v
    None => "unknown"
  }
}
