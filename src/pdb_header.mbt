///|
/// Bio.PDB.ParsePDBHeader - PDB metadata record parsing.
/// Parses PDB file header records including HEADER, TITLE, COMPOUND,
/// SOURCE, KEYWDS, REMARK, AUTH, JRNL, and DBREF.

///|
/// PDB header information containing metadata about the structure.
pub struct PDBHeader {
  pdb_id : String
  title : String
  compound : CompoundInfo
  source : SourceInfo
  keywords : Array[String]
  resolution : Double
  rfactor : Double
  date_deposited : String
  authors : Array[String]
  journal : String
  remark_lines : Array[String]
  dbrefs : Array[DBRefEntry]
  chain_ids : Array[String]
}

///|
/// Construct a PDBHeader with default values.
pub fn PDBHeader::new(
  pdb_id? : String = "",
  title? : String = "",
  compound? : CompoundInfo = CompoundInfo::new(),
  source? : SourceInfo = SourceInfo::new(),
  keywords? : Array[String] = [],
  resolution? : Double = 0.0,
  rfactor? : Double = 0.0,
  date_deposited? : String = "",
  authors? : Array[String] = [],
  journal? : String = "",
  remark_lines? : Array[String] = [],
  dbrefs? : Array[DBRefEntry] = [],
  chain_ids? : Array[String] = [],
) -> PDBHeader {
  PDBHeader::{
    pdb_id,
    title,
    compound,
    source,
    keywords,
    resolution,
    rfactor,
    date_deposited,
    authors,
    journal,
    remark_lines,
    dbrefs,
    chain_ids,
  }
}

///|
/// Compound information from COMPOUND record.
pub struct CompoundInfo {
  molecule_id : String
  chain_ids : Array[String]
  compound_name : String
  formula : String
  ka : String
  kb : String
  ec_number : String
  engineered : Bool
  mutation : Bool
}

///|
/// Construct a CompoundInfo with default values.
pub fn CompoundInfo::new(
  molecule_id? : String = "",
  chain_ids? : Array[String] = [],
  compound_name? : String = "",
  formula? : String = "",
  ka? : String = "",
  kb? : String = "",
  ec_number? : String = "",
  engineered? : Bool = false,
  mutation? : Bool = false,
) -> CompoundInfo {
  CompoundInfo::{
    molecule_id,
    chain_ids,
    compound_name,
    formula,
    ka,
    kb,
    ec_number,
    engineered,
    mutation,
  }
}

///|
/// Source organism information from SOURCE record.
pub struct SourceInfo {
  source : String
  organism_scientific : String
  organism_common : String
  strain : String
  cell_line : String
  tissue : String
  tax_id : String
}

///|
/// Construct a SourceInfo with default values.
pub fn SourceInfo::new(
  source? : String = "",
  organism_scientific? : String = "",
  organism_common? : String = "",
  strain? : String = "",
  cell_line? : String = "",
  tissue? : String = "",
  tax_id? : String = "",
) -> SourceInfo {
  SourceInfo::{
    source,
    organism_scientific,
    organism_common,
    strain,
    cell_line,
    tissue,
    tax_id,
  }
}

///|
/// Database cross-reference entry from DBREF record.
pub struct DBRefEntry {
  db_name : String
  db_accession : String
  db_id : String
}

///|
/// Construct a DBRefEntry with default values.
pub fn DBRefEntry::new(
  db_name? : String = "",
  db_accession? : String = "",
  db_id? : String = "",
) -> DBRefEntry {
  DBRefEntry::{ db_name, db_accession, db_id }
}

///|
/// Convert StringView to String with whitespace trimming.
fn trim_sv(sv : StringView) -> String {
  sv.trim().to_owned()
}

///|
/// Convert a string to uppercase.
fn to_upper(s : String) -> String {
  let result = StringBuilder::new()
  let mut i = 0
  while i < s.length() {
    let code = s.unsafe_get(i).to_int()
    if code >= 97 && code <= 122 {
      let ch = match (code - 32).to_char() {
        Some(c) => c
        None => 'A'
      }
      result.write_char(ch)
    } else {
      let ch = match code.to_char() {
        Some(c) => c
        None => ' '
      }
      result.write_char(ch)
    }
    i = i + 1
  }
  result.to_string()
}

///|
/// Split PDB text into lines, removing empty lines.
fn split_pdb_lines(text : String) -> Array[String] {
  let parts = text.split("\n").to_array()
  let result = Array::new()
  for part in parts {
    let line = part.to_string()
    if line.length() > 0 {
      result.push(line)
    }
  }
  result
}

///|
/// Find the start and end indices of a record section in the lines array.
/// Returns (-1, -1) if the record is not found.
fn extract_record_section(lines : Array[String], record_name : String) -> (Int, Int) {
  let mut start = -1
  let mut end = -1
  let mut i = 0
  while i < lines.length() {
    let line = lines[i]
    if line.length() >= 6 && trim_sv(line[0:6]) == record_name {
      if start == -1 {
        start = i
      }
      end = i
    } else if start != -1 {
      break
    }
    i = i + 1
  }
  (start, end)
}

///|
/// Parse semicolon-separated fields into an array of trimmed strings.
fn parse_semicolon_fields(text : String) -> Array[String] {
  let parts = text.split(";").to_array()
  let result = Array::new()
  for part in parts {
    let trimmed = trim_sv(part)
    if trimmed.length() > 0 {
      result.push(trimmed)
    }
  }
  result
}

///|
/// Extract a fixed-width field from a PDB line using 1-indexed column positions.
/// Moonbit slices are exclusive on the end: line[a:b] gives chars at indices a..b-1.
/// So to extract columns start_col..end_col (inclusive), we use
///   line[start_col-1 : end_col]
/// which gives indices start_col-1 .. end_col-1 (inclusive), length = end_col - start_col + 1.
fn trim_pdb_field(line : String, start_col : Int, end_col : Int) -> String {
  let start_idx = start_col - 1
  let end_idx = end_col
  if line.length() < start_idx {
    return ""
  }
  let actual_end = if line.length() < end_idx { line.length() } else { end_idx }
  trim_sv(line[start_idx:actual_end])
}

///|
/// Clean whitespace from PDB text.
fn clean_pdb_text(text : String) -> String {
  text.trim().to_owned()
}

///|
/// Parse a double value from a string, returning 0.0 on failure.
fn parse_double_safe(s : String) -> Double {
  let trimmed = s.trim().to_owned()
  if trimmed.length() == 0 {
    return 0.0
  }
  let mut result = 0.0
  let mut decimal = 0.0
  let mut decimal_places = 0
  let mut is_negative = false
  let zero = '0'.to_int()
  let mut i = 0

  if trimmed.unsafe_get(0).to_int() == 45 {
    is_negative = true
    i = 1
  }

  while i < trimmed.length() {
    let c = trimmed.unsafe_get(i).to_int()
    if c == 46 {
      i = i + 1
      while i < trimmed.length() {
        let dc = trimmed.unsafe_get(i).to_int()
        if dc >= zero && dc <= zero + 9 {
          decimal = decimal * 10.0 + (dc - zero).to_double()
          decimal_places = decimal_places + 1
        }
        i = i + 1
      }
      break
    }
    if c >= zero && c <= zero + 9 {
      result = result * 10.0 + (c - zero).to_double()
    }
    i = i + 1
  }

  while decimal_places > 0 {
    decimal = decimal / 10.0
    decimal_places = decimal_places - 1
  }

  let final_result = result + decimal
  if is_negative {
    -final_result
  } else {
    final_result
  }
}

///|
/// Parse a single COMPOUND or SOURCE field line into key-value pairs.
/// Returns an array of (key, value) tuples.
fn parse_record_fields(text : String) -> Array[(String, String)] {
  let fields = parse_semicolon_fields(text)
  let result : Array[(String, String)] = []
  for field in fields {
    let cleaned = strip_field_prefix(field)
    if cleaned.length() == 0 {
      continue
    }
    let colon_idx = find_colon(cleaned)
    if colon_idx >= 0 {
      let key = trim_sv(cleaned[0:colon_idx])
      let value = trim_sv(cleaned[colon_idx + 1:])
      if key.length() > 0 {
        result.push((key, value))
      }
    } else {
      result.push(("", cleaned))
    }
  }
  result
}

///|
/// Strip leading field prefix (numbers like "1 ", "2 ", etc.).
fn strip_field_prefix(field : String) -> String {
  let mut i = 0
  while i < field.length() {
    let c = field.unsafe_get(i).to_int()
    if c >= '0'.to_int() && c <= '9'.to_int() {
      i = i + 1
    } else if c == ' '.to_int() {
      i = i + 1
    } else {
      break
    }
  }
  if i > 0 {
    field[i:].to_string()
  } else {
    field
  }
}

///|
/// Find the first colon position in a string.
fn find_colon(s : String) -> Int {
  let mut i = 0
  while i < s.length() {
    if s.unsafe_get(i).to_int() == 58 {
      return i
    }
    i = i + 1
  }
  -1
}

///|
/// Join text content from a multi-line record section.
/// Extracts text from cols 11-80 (index 10-80) of each line.
fn join_record_text(lines : Array[String], start : Int, end : Int) -> String {
  let result = StringBuilder::new()
  let mut i = start
  while i <= end && i < lines.length() {
    let line = lines[i]
    if line.length() > 10 {
      let text = line[10:].to_string()
      result.write_string(text)
    }
    i = i + 1
  }
  result.to_string()
}

///|
/// Parse the chain IDs string (e.g., "A, B, C") into an array of strings.
fn parse_chain_ids_str(s : String) -> Array[String] {
  let parts = s.split(",").to_array()
  let result = Array::new()
  for part in parts {
    let trimmed = trim_sv(part)
    if trimmed.length() > 0 {
      result.push(trimmed)
    }
  }
  result
}

// === Public API Functions ===

///|
/// Parse full PDB file header from text content.
pub fn parse_pdb_header(pdb_text : String) -> PDBHeader {
  let lines = split_pdb_lines(pdb_text)
  parse_pdb_header_lines(lines)
}

///|
/// Parse PDB header from pre-split lines.
pub fn parse_pdb_header_lines(lines : Array[String]) -> PDBHeader {
  let mut pdb_id = ""
  let mut date_deposited = ""
  let mut title_start = ""

  let mut i = 0
  while i < lines.length() {
    let line = lines[i]
    if line.length() >= 6 && line[0:6] == "HEADER" {
      let (id, date, ts) = parse_header_record(line)
      pdb_id = id
      date_deposited = date
      title_start = ts
      break
    }
    i = i + 1
  }

  let title = parse_title_record(lines, 0)
  let compound = parse_compound_record(lines, 0)
  let source = parse_source_record(lines, 0)
  let keywords = parse_keywords_record(lines, 0)
  let remark_lines = parse_remark_records(lines, 0)
  let resolution = parse_resolution_from_lines(lines)
  let rfactor = parse_rfactor_from_lines(lines)
  let authors = parse_authors_record(lines, 0)
  let journal = parse_journal_record(lines, 0)
  let dbrefs = parse_dbref_record(lines, 0)
  let chain_ids = extract_chain_ids_from_compound(compound)

  let final_chain_ids = if chain_ids.length() > 0 {
    chain_ids
  } else {
    extract_chain_ids_from_dbref(dbrefs)
  }

  PDBHeader::new(
    pdb_id=pdb_id,
    title=title,
    compound=compound,
    source=source,
    keywords=keywords,
    resolution=resolution,
    rfactor=rfactor,
    date_deposited=date_deposited,
    authors=authors,
    journal=journal,
    remark_lines=remark_lines,
    dbrefs=dbrefs,
    chain_ids=final_chain_ids,
  )
}

///|
/// Parse HEADER record line.
/// Returns (pdb_id, date, title_start) where title_start is the classification field.
pub fn parse_header_record(line : String) -> (String, String, String) {
  let pdb_id = trim_pdb_field(line, 50, 54)
  let date = trim_pdb_field(line, 19, 28)
  let title_start = trim_pdb_field(line, 13, 17)
  (pdb_id, date, title_start)
}

///|
/// Parse multi-line TITLE record.
pub fn parse_title_record(lines : Array[String], _start_idx : Int) -> String {
  let (start, end) = extract_record_section(lines, "TITLE")
  if start < 0 {
    return ""
  }
  let text = join_record_text(lines, start, end)
  clean_pdb_text(text)
}

///|
/// Parse multi-line COMPOUND record.
pub fn parse_compound_record(lines : Array[String], _start_idx : Int) -> CompoundInfo {
  let (start, end) = extract_record_section(lines, "COMPND")
  if start < 0 {
    return CompoundInfo::new()
  }
  let text = join_record_text(lines, start, end)
  let fields = parse_record_fields(text)

  let mut molecule_id = ""
  let mut chain_ids : Array[String] = []
  let mut compound_name = ""
  let mut formula = ""
  let mut ka = ""
  let mut kb = ""
  let mut ec_number = ""
  let mut engineered = false
  let mut mutation = false

  let mut fi = 0
  while fi < fields.length() {
    let (key, value) = fields[fi]
    let upper_key = to_upper(key)
    if upper_key == "MOL_ID" {
      molecule_id = value
    } else if upper_key == "CHAIN" {
      chain_ids = parse_chain_ids_str(value)
    } else if upper_key == "COMPOUND" || upper_key == "COMPOUND_NAME" {
      compound_name = value
    } else if upper_key == "FORMULA" {
      formula = value
    } else if upper_key == "KA" {
      ka = value
    } else if upper_key == "KB" {
      kb = value
    } else if upper_key == "EC" || upper_key == "EC_NUMBER" {
      ec_number = value
    } else if upper_key == "ENGINEERED" {
      engineered = to_upper(value) == "YES" || value == "Y"
    } else if upper_key == "MUTATION" {
      mutation = to_upper(value) == "YES" || value == "Y"
    }
    fi = fi + 1
  }

  CompoundInfo::new(
    molecule_id=molecule_id,
    chain_ids=chain_ids,
    compound_name=compound_name,
    formula=formula,
    ka=ka,
    kb=kb,
    ec_number=ec_number,
    engineered=engineered,
    mutation=mutation,
  )
}

///|
/// Parse multi-line SOURCE record.
pub fn parse_source_record(lines : Array[String], _start_idx : Int) -> SourceInfo {
  let (start, end) = extract_record_section(lines, "SOURCE")
  if start < 0 {
    return SourceInfo::new()
  }
  let text = join_record_text(lines, start, end)
  let fields = parse_record_fields(text)

  let mut source = ""
  let mut organism_scientific = ""
  let mut organism_common = ""
  let mut strain = ""
  let mut cell_line = ""
  let mut tissue = ""
  let mut tax_id = ""

  let mut fi = 0
  while fi < fields.length() {
    let (key, value) = fields[fi]
    let upper_key = to_upper(key)
    if upper_key == "ORGANISM" || upper_key == "ORGANISM_SCIENTIFIC" {
      organism_scientific = value
    } else if upper_key == "COMMON_NAME" || upper_key == "ORGANISM_COMMON" {
      organism_common = value
    } else if upper_key == "STRAIN" {
      strain = value
    } else if upper_key == "CELL_LINE" {
      cell_line = value
    } else if upper_key == "TISSUE" {
      tissue = value
    } else if upper_key == "TAX_ID" {
      tax_id = value
    } else if upper_key == "SOURCE" {
      if source.length() == 0 {
        source = value
      }
    }
    fi = fi + 1
  }

  SourceInfo::new(
    source=source,
    organism_scientific=organism_scientific,
    organism_common=organism_common,
    strain=strain,
    cell_line=cell_line,
    tissue=tissue,
    tax_id=tax_id,
  )
}

///|
/// Parse multi-line KEYWDS record.
pub fn parse_keywords_record(lines : Array[String], _start_idx : Int) -> Array[String] {
  let (start, end) = extract_record_section(lines, "KEYWDS")
  if start < 0 {
    return []
  }
  let text = join_record_text(lines, start, end)
  let parts = text.split(";").to_array()
  let result = Array::new()
  for part in parts {
    let trimmed = trim_sv(part)
    if trimmed.length() > 0 {
      result.push(trimmed)
    }
  }
  result
}

///|
/// Parse REMARK 2 RESOLUTION record line.
pub fn parse_resolution_record(line : String) -> Double {
  let idx = find_keyword(line, "RESOLUTION.")
  if idx < 0 {
    return 0.0
  }
  let after = line[idx + 11:].to_string()
  let trimmed = after.trim().to_owned()
  if trimmed.length() == 0 {
    return 0.0
  }
  parse_first_number(trimmed)
}

///|
/// Parse REMARK 3 R FACTOR record line.
pub fn parse_rfactor_record(line : String) -> Double {
  let idx = find_keyword(line, "R FACTOR.")
  if idx < 0 {
    return 0.0
  }
  let after = line[idx + 9:].to_string()
  let trimmed = after.trim().to_owned()
  if trimmed.length() == 0 {
    return 0.0
  }
  parse_first_number(trimmed)
}

///|
/// Parse multi-line AUTH record.
/// Returns array of author names (semicolon-separated).
pub fn parse_authors_record(lines : Array[String], _start_idx : Int) -> Array[String] {
  let (start, end) = extract_record_section(lines, "AUTH")
  if start < 0 {
    return []
  }
  let text = join_record_text(lines, start, end)
  let parts = text.split(";").to_array()
  let result = Array::new()
  for part in parts {
    let trimmed = trim_sv(part)
    if trimmed.length() > 0 {
      result.push(trimmed)
    }
  }
  result
}

///|
/// Parse multi-line JRNL record.
pub fn parse_journal_record(lines : Array[String], _start_idx : Int) -> String {
  let (start, end) = extract_record_section(lines, "JRNL")
  if start < 0 {
    return ""
  }
  let text = join_record_text(lines, start, end)
  clean_pdb_text(text)
}

///|
/// Parse DBREF records.
pub fn parse_dbref_record(lines : Array[String], _start_idx : Int) -> Array[DBRefEntry] {
  let result : Array[DBRefEntry] = []
  let (start, end) = extract_record_section(lines, "DBREF")
  if start < 0 {
    return result
  }
  let mut i = start
  while i <= end && i < lines.length() {
    let line = lines[i]
    if line.length() >= 56 {
      let db_name = trim_pdb_field(line, 13, 18)
      let db_accession = trim_pdb_field(line, 24, 33)
      let db_id = trim_pdb_field(line, 49, 56)
      result.push(DBRefEntry::new(
        db_name=db_name,
        db_accession=db_accession,
        db_id=db_id,
      ))
    }
    i = i + 1
  }
  result
}

///|
/// Parse all REMARK records into an array of strings.
pub fn parse_remark_records(lines : Array[String], _start_idx : Int) -> Array[String] {
  let result : Array[String] = []
  let (start, end) = extract_record_section(lines, "REMARK")
  if start < 0 {
    return result
  }
  let mut i = start
  while i <= end && i < lines.length() {
    let line = lines[i]
    if line.length() > 11 {
      let remark_text = trim_sv(line[11:])
      result.push(remark_text)
    }
    i = i + 1
  }
  result
}

///|
/// Extract chain IDs from compound information.
pub fn extract_chain_ids(header : PDBHeader) -> Array[String] {
  header.compound.chain_ids.copy()
}

// === Internal helper functions ===

///|
/// Parse resolution value from all REMARK 2 lines.
fn parse_resolution_from_lines(lines : Array[String]) -> Double {
  let mut resolution = 0.0
  let mut i = 0
  while i < lines.length() {
    let line = lines[i]
    if line.length() >= 6 && line[0:6] == "REMARK" {
      if line.length() >= 10 && trim_sv(line[6:10]) == "2" {
        resolution = parse_resolution_record(line)
      }
    }
    i = i + 1
  }
  resolution
}

///|
/// Parse rfactor value from all REMARK 3 lines.
fn parse_rfactor_from_lines(lines : Array[String]) -> Double {
  let mut rfactor = 0.0
  let mut i = 0
  while i < lines.length() {
    let line = lines[i]
    if line.length() >= 6 && line[0:6] == "REMARK" {
      if line.length() >= 10 && trim_sv(line[6:10]) == "3" {
        rfactor = parse_rfactor_record(line)
      }
    }
    i = i + 1
  }
  rfactor
}

///|
/// Find the position of a keyword in a line (case-insensitive).
fn find_keyword(line : String, keyword : String) -> Int {
  let line_upper = to_upper(line)
  let keyword_upper = to_upper(keyword)
  let kw_len = keyword_upper.length()
  if line_upper.length() < kw_len {
    return -1
  }
  let mut i = 0
  while i <= line_upper.length() - kw_len {
    if line_upper[i:i + kw_len] == keyword_upper {
      return i
    }
    i = i + 1
  }
  -1
}

///|
/// Parse the first numeric value from a string.
fn parse_first_number(s : String) -> Double {
  let n = s.length()
  if n == 0 {
    return 0.0
  }

  let mut i = 0
  while i < n {
    let c = s.unsafe_get(i).to_int()
    if (c >= '0'.to_int() && c <= '9'.to_int()) || c == 45 {
      break
    }
    i = i + 1
  }

  if i >= n {
    return 0.0
  }

  let start = i
  let mut has_dot = false
  while i < n {
    let c = s.unsafe_get(i).to_int()
    if c >= '0'.to_int() && c <= '9'.to_int() {
      i = i + 1
    } else if c == 46 && !has_dot {
      has_dot = true
      i = i + 1
    } else {
      break
    }
  }

  if i == start {
    return 0.0
  }

  let num_str = s[start:i].to_string()
  parse_double_safe(num_str)
}

///|
/// Extract chain IDs from compound info.
fn extract_chain_ids_from_compound(compound : CompoundInfo) -> Array[String] {
  compound.chain_ids.copy()
}

///|
/// Extract chain IDs from DBREF entries as fallback.
fn extract_chain_ids_from_dbref(dbrefs : Array[DBRefEntry]) -> Array[String] {
  let result : Array[String] = []
  result
}