///|
/// GEO SOFT format parser - Gene Expression Omnibus Simple Omnibus Format.
///
/// References:
///   - Biopython Bio.Geo module
///   - Barrett T et al. (2013) NCBI GEO: archive for functional genomics
///     data sets--update. Nucleic Acids Res 41(Database issue):D991-5.
///   - Edgar R, Domrachev M, Lash AE (2002) Gene Expression Omnibus: NCBI
///     gene expression and hybridization array data repository. Nucleic
///     Acids Res 30(1):207-10.
///
/// SOFT (Simple Omnibus Format in Text) is used by NCBI GEO to distribute
/// platform, series, and sample records. Each record begins with a caret
/// entity marker (^), followed by attribute lines (!), column definitions
/// (#), and data rows.
///
/// All public identifiers are prefixed with `geo_` to avoid conflicts.

// ============================================================================
// Data structures
// ============================================================================

///|
/// Type of GEO SOFT entity.
pub enum GeoEntityType {
  Platform
  Sample
  Series
  Dataset
  Unknown
}

///|
/// Create a Platform entity type.
pub fn geo_platform() -> GeoEntityType {
  Platform
}

///|
/// Create a Sample entity type.
pub fn geo_sample() -> GeoEntityType {
  Sample
}

///|
/// Create a Series entity type.
pub fn geo_series() -> GeoEntityType {
  Series
}

///|
/// Create a Dataset entity type.
pub fn geo_dataset() -> GeoEntityType {
  Dataset
}

///|
/// Create an Unknown entity type.
pub fn geo_unknown() -> GeoEntityType {
  Unknown
}

///|
/// Convert entity type to string.
pub fn GeoEntityType::to_string(self : GeoEntityType) -> String {
  match self {
    Platform => "Platform"
    Sample => "Sample"
    Series => "Series"
    Dataset => "Dataset"
    Unknown => "Unknown"
  }
}

///|
/// A key-value attribute (lines starting with !).
pub struct GeoAttribute {
  name : String
  value : String
}

///|
/// Create a GeoAttribute.
pub fn GeoAttribute::new(name : String, value : String) -> GeoAttribute {
  { name, value }
}

///|
/// Get attribute name.
pub fn GeoAttribute::name(self : GeoAttribute) -> String {
  self.name
}

///|
/// Get attribute value.
pub fn GeoAttribute::value(self : GeoAttribute) -> String {
  self.value
}

///|
/// A column definition (lines starting with #).
pub struct GeoColumnDef {
  name : String
  description : String
}

///|
/// Create a GeoColumnDef.
pub fn GeoColumnDef::new(name : String, description : String) -> GeoColumnDef {
  { name, description }
}

///|
/// Get column name.
pub fn GeoColumnDef::name(self : GeoColumnDef) -> String {
  self.name
}

///|
/// Get column description.
pub fn GeoColumnDef::description(self : GeoColumnDef) -> String {
  self.description
}

///|
/// A data table within a GEO record.
pub struct GeoDataTable {
  columns : Array[GeoColumnDef]
  rows : Array[Array[String]]
}

///|
/// Create a GeoDataTable.
pub fn GeoDataTable::new(
  columns : Array[GeoColumnDef],
  rows : Array[Array[String]],
) -> GeoDataTable {
  { columns, rows }
}

///|
/// Get column definitions.
pub fn GeoDataTable::columns(self : GeoDataTable) -> Array[GeoColumnDef] {
  self.columns
}

///|
/// Get data rows.
pub fn GeoDataTable::rows(self : GeoDataTable) -> Array[Array[String]] {
  self.rows
}

///|
/// Number of columns.
pub fn GeoDataTable::n_columns(self : GeoDataTable) -> Int {
  self.columns.length()
}

///|
/// Number of rows.
pub fn GeoDataTable::n_rows(self : GeoDataTable) -> Int {
  self.rows.length()
}

///|
/// Get a column index by name. Returns -1 if not found.
pub fn GeoDataTable::column_index(self : GeoDataTable, name : String) -> Int {
  for i in 0..<self.columns.length() {
    if self.columns[i].name == name {
      return i
    }
  }
  -1
}

///|
/// Get values from a specific column as an array of strings.
pub fn GeoDataTable::get_column(
  self : GeoDataTable,
  name : String,
) -> Array[String] {
  let idx = self.column_index(name)
  let result : Array[String] = []
  if idx < 0 {
    return result
  }
  for row in self.rows {
    if idx < row.length() {
      result.push(row[idx])
    }
  }
  result
}

///|
/// A GEO SOFT record: platform, sample, or series.
pub struct GeoRecord {
  entity_type : GeoEntityType
  accession : String
  attributes : Array[GeoAttribute]
  data_table : GeoDataTable
}

///|
/// Create a GeoRecord.
pub fn GeoRecord::new(
  entity_type : GeoEntityType,
  accession : String,
  attributes : Array[GeoAttribute],
  data_table : GeoDataTable,
) -> GeoRecord {
  { entity_type, accession, attributes, data_table }
}

///|
/// Get entity type.
pub fn GeoRecord::entity_type(self : GeoRecord) -> GeoEntityType {
  self.entity_type
}

///|
/// Get accession ID.
pub fn GeoRecord::accession(self : GeoRecord) -> String {
  self.accession
}

///|
/// Get attributes.
pub fn GeoRecord::attributes(self : GeoRecord) -> Array[GeoAttribute] {
  self.attributes
}

///|
/// Get data table.
pub fn GeoRecord::data_table(self : GeoRecord) -> GeoDataTable {
  self.data_table
}

///|
/// Get an attribute value by name. Returns None if not found.
pub fn GeoRecord::get_attribute(self : GeoRecord, name : String) -> String? {
  for attr in self.attributes {
    if attr.name == name {
      return Some(attr.value)
    }
  }
  None
}

///|
/// Get the title attribute (common to most records).
/// Searches for "title", "Platform_title", "Sample_title", or "Series_title".
pub fn GeoRecord::title(self : GeoRecord) -> String {
  match self.get_attribute("title") {
    Some(v) => v
    None =>
      match self.get_attribute("Platform_title") {
        Some(v) => v
        None =>
          match self.get_attribute("Sample_title") {
            Some(v) => v
            None =>
              match self.get_attribute("Series_title") {
                Some(v) => v
                None => ""
              }
          }
      }
  }
}

///|
/// Top-level container for parsed SOFT data.
pub struct GeoFile {
  records : Array[GeoRecord]
}

///|
/// Create a GeoFile.
pub fn GeoFile::new(records : Array[GeoRecord]) -> GeoFile {
  { records }
}

///|
/// Get all records.
pub fn GeoFile::records(self : GeoFile) -> Array[GeoRecord] {
  self.records
}

///|
/// Number of records.
pub fn GeoFile::n_records(self : GeoFile) -> Int {
  self.records.length()
}

///|
/// Get records by entity type.
pub fn GeoFile::get_by_type(
  self : GeoFile,
  etype : GeoEntityType,
) -> Array[GeoRecord] {
  let result : Array[GeoRecord] = []
  for r in self.records {
    match (r.entity_type, etype) {
      (Platform, Platform) => result.push(r)
      (Sample, Sample) => result.push(r)
      (Series, Series) => result.push(r)
      (Dataset, Dataset) => result.push(r)
      _ => continue
    }
  }
  result
}

///|
/// Get a record by accession ID.
pub fn GeoFile::get_by_accession(
  self : GeoFile,
  accession : String,
) -> GeoRecord? {
  for r in self.records {
    if r.accession == accession {
      return Some(r)
    }
  }
  None
}

// ============================================================================
// Parsing
// ============================================================================

///|
/// Tokenize a line by tab character.
fn geo_split_tab(line : String) -> Array[String] {
  split_by_char(line, '\t'.to_int())
}

///|
/// Determine entity type from the caret marker keyword.
fn geo_entity_type_from_keyword(keyword : String) -> GeoEntityType {
  if keyword == "PLATFORM" {
    Platform
  } else if keyword == "SAMPLE" {
    Sample
  } else if keyword == "SERIES" {
    Series
  } else if keyword == "DATASET" {
    Dataset
  } else {
    Unknown
  }
}

///|
/// Parse an attribute line: !Name = value or !Name_value = value
/// Returns (name, value) or None if not a valid attribute.
fn geo_parse_attribute_line(line : String) -> GeoAttribute? {
  // line starts with "!"
  let rest = line[1:line.length()].to_string()
  // Find "=" position
  let eq_pos = geo_find_char(rest, '='.to_int())
  if eq_pos < 0 {
    return None
  }
  let name = trim(rest[0:eq_pos].to_string())
  let value = trim(rest[eq_pos + 1:rest.length()].to_string())
  Some(GeoAttribute::new(name, value))
}

///|
/// Find the first occurrence of a character in a string. Returns -1 if not found.
fn geo_find_char(s : String, ch : Int) -> Int {
  let n = s.length()
  let mut i = 0
  while i < n {
    if s.unsafe_get(i).to_int() == ch {
      return i
    }
    i = i + 1
  }
  -1
}

///|
/// Parse a column definition line: #NAME = description
fn geo_parse_column_line(line : String) -> GeoColumnDef? {
  // line starts with "#"
  let rest = line[1:line.length()].to_string()
  let eq_pos = geo_find_char(rest, '='.to_int())
  if eq_pos < 0 {
    return None
  }
  let name = trim(rest[0:eq_pos].to_string())
  let desc = trim(rest[eq_pos + 1:rest.length()].to_string())
  Some(GeoColumnDef::new(name, desc))
}

///|
/// Parse a complete SOFT file from text input.
pub fn geo_parse(text : String) -> GeoFile {
  let lines = split_by_char(text, '\n'.to_int())
  let records : Array[GeoRecord] = []
  let mut idx = 0
  while idx < lines.length() {
    let raw_line = lines[idx].to_string()
    let line = trim(raw_line)
    if line.length() == 0 {
      idx = idx + 1
      continue
    }
    // Check for entity marker (^PLATFORM, ^SAMPLE, etc.)
    if line.unsafe_get(0).to_int() == '^'.to_int() {
      let rest = line[1:line.length()].to_string()
      let eq_pos = geo_find_char(rest, '='.to_int())
      if eq_pos >= 0 {
        let keyword = trim(rest[0:eq_pos].to_string())
        let accession = trim(rest[eq_pos + 1:rest.length()].to_string())
        let etype = geo_entity_type_from_keyword(keyword)
        let (record, next) = geo_parse_record(
          lines,
          idx + 1,
          etype,
          accession,
        )
        records.push(record)
        idx = next
      } else {
        idx = idx + 1
      }
    } else {
      idx = idx + 1
    }
  }
  GeoFile::new(records)
}

///|
/// Parse a single GEO record starting after the entity marker line.
fn geo_parse_record(
  lines : Array[String],
  start : Int,
  etype : GeoEntityType,
  accession : String,
) -> (GeoRecord, Int) {
  let attributes : Array[GeoAttribute] = []
  let columns : Array[GeoColumnDef] = []
  let rows : Array[Array[String]] = []
  let mut in_table = false
  let mut idx = start
  while idx < lines.length() {
    let raw_line = lines[idx].to_string()
    let line = trim(raw_line)
    if line.length() == 0 {
      idx = idx + 1
      continue
    }
    // Next entity marker means this record ends
    if line.unsafe_get(0).to_int() == '^'.to_int() {
      break
    }
    if line.unsafe_get(0).to_int() == '!'.to_int() {
      // Attribute line
      match geo_parse_attribute_line(line) {
        Some(attr) => attributes.push(attr)
        None => ()
      }
      in_table = false
    } else if line.unsafe_get(0).to_int() == '#'.to_int() {
      // Column definition
      match geo_parse_column_line(line) {
        Some(col) => columns.push(col)
        None => ()
      }
      in_table = true
    } else {
      // Data row (if we're in a table section with columns defined)
      if in_table && columns.length() > 0 {
        let fields = geo_split_tab(line)
        rows.push(fields)
      }
    }
    idx = idx + 1
  }
  (
    GeoRecord::new(
      etype,
      accession,
      attributes,
      GeoDataTable::new(columns, rows),
    ),
    idx,
  )
}

// ============================================================================
// Sample data
// ============================================================================

///|
/// Generate sample SOFT format text with a platform and two samples.
pub fn geo_sample_text() -> String {
  let mut s = "^PLATFORM = GPL1\n"
  s = s + "!Platform_title = Test Affymetrix Platform\n"
  s = s + "!Platform_organism = Homo sapiens\n"
  s = s + "!Platform_technology = in_situ_oligonucleotide\n"
  s = s + "#ID = Probe Set ID\n"
  s = s + "#VALUE = Expression Value\n"
  s = s + "1007_s_at\t8.5\n"
  s = s + "117_at\t7.2\n"
  s = s + "1294_at\t9.1\n"
  s = s + "^SAMPLE = GSM1\n"
  s = s + "!Sample_title = Control Sample 1\n"
  s = s + "!Sample_source = brain tissue\n"
  s = s + "!Sample_characteristics = normal\n"
  s = s + "#ID = Probe Set ID\n"
  s = s + "#VALUE = Expression Value\n"
  s = s + "1007_s_at\t6.2\n"
  s = s + "117_at\t5.8\n"
  s = s + "1294_at\t7.4\n"
  s = s + "^SAMPLE = GSM2\n"
  s = s + "!Sample_title = Treatment Sample 1\n"
  s = s + "!Sample_source = brain tissue\n"
  s = s + "!Sample_characteristics = treated\n"
  s = s + "#ID = Probe Set ID\n"
  s = s + "#VALUE = Expression Value\n"
  s = s + "1007_s_at\t9.1\n"
  s = s + "117_at\t8.7\n"
  s = s + "1294_at\t10.2\n"
  s
}

// ============================================================================
// Utilities
// ============================================================================

///|
/// Format a GeoAttribute as a readable string.
pub fn GeoAttribute::to_string(self : GeoAttribute) -> String {
  self.name + " = " + self.value
}

///|
/// Format a GeoDataTable as a readable string.
pub fn GeoDataTable::to_string(self : GeoDataTable) -> String {
  let mut s = "DataTable(" + self.n_columns().to_string() + " cols, " +
    self.n_rows().to_string() + " rows)\n"
  // Column headers
  s = s + "  "
  for i in 0..<self.columns.length() {
    if i > 0 {
      s = s + "\t"
    }
    s = s + self.columns[i].name
  }
  s = s + "\n"
  // First few rows
  let max_rows = if self.n_rows() < 5 { self.n_rows() } else { 5 }
  for i in 0..<max_rows {
    s = s + "  "
    for j in 0..<self.rows[i].length() {
      if j > 0 {
        s = s + "\t"
      }
      s = s + self.rows[i][j]
    }
    s = s + "\n"
  }
  if self.n_rows() > max_rows {
    s = s + "  ... (" + (self.n_rows() - max_rows).to_string() + " more rows)\n"
  }
  s
}

///|
/// Format a GeoRecord as a readable string.
pub fn GeoRecord::to_string(self : GeoRecord) -> String {
  let mut s = "GeoRecord(type=" + self.entity_type.to_string() + ", accession=" +
    self.accession + ")\n"
  s = s + "Attributes (" + self.attributes.length().to_string() + "):\n"
  for attr in self.attributes {
    s = s + "  " + attr.to_string() + "\n"
  }
  s = s + self.data_table.to_string()
  s
}

///|
/// Format a GeoFile as a readable summary.
pub fn GeoFile::to_string(self : GeoFile) -> String {
  let mut s = "GEO SOFT File (" + self.n_records().to_string() + " records)\n"
  for r in self.records {
    s = s + r.to_string() + "\n"
  }
  s
}
