///|
/// Bio.UniProt.GOA - Gene Association Format (GAF) parser.
///
/// References:
///   - Biopython Bio.UniProt.GOA module
///   - Gene Ontology Consortium GAF specification
///     (http://geneontology.org/docs/go-annotation-file-gaf-format-2.2/)
///   - Ashburner M et al. (2000) Gene Ontology: tool for the unification of
///     biology. Nat Genet 25(1):25-9.
///   - The Gene Ontology Consortium (2019) The Gene Ontology Resource: 20
///     years and still GOing strong. Nucleic Acids Res 47(D1):D330-D338.
///
/// GAF (Gene Association Format) is the tab-separated text format used by
/// the Gene Ontology consortium to distribute functional annotations of
/// gene products. Header lines begin with "!" (e.g. "!gaf-version: 2.2").
/// GAF 2.2 data rows have 17 tab-separated columns.
///
/// All public identifiers are prefixed with `goa_` to avoid conflicts.

// ============================================================================
// Data structures
// ============================================================================

///|
/// The three GO aspects (namespaces), plus an Unknown fallback.
pub enum GafAspect {
  MolecularFunction
  BiologicalProcess
  CellularComponent
  Unknown
}

///|
/// Parse a single-character aspect code into a GafAspect.
/// "F" -> MolecularFunction, "P" -> BiologicalProcess,
/// "C" -> CellularComponent, anything else -> Unknown.
pub fn GafAspect::from_string(s : String) -> GafAspect {
  if s == "F" {
    MolecularFunction
  } else if s == "P" {
    BiologicalProcess
  } else if s == "C" {
    CellularComponent
  } else {
    Unknown
  }
}

///|
/// Convert a GafAspect back to its single-letter code.
pub fn GafAspect::to_string(self : GafAspect) -> String {
  match self {
    MolecularFunction => "F"
    BiologicalProcess => "P"
    CellularComponent => "C"
    Unknown => "?"
  }
}

///|
/// Human-readable description of a GafAspect.
pub fn GafAspect::description(self : GafAspect) -> String {
  match self {
    MolecularFunction => "Molecular Function"
    BiologicalProcess => "Biological Process"
    CellularComponent => "Cellular Component"
    Unknown => "Unknown"
  }
}

///|
/// A single GAF annotation record with all 17 columns of GAF 2.2.
pub struct GafRecord {
  db : String
  db_object_id : String
  db_object_symbol : String
  qualifier : String
  go_id : String
  db_reference : String
  evidence_code : String
  with_from : String
  aspect : String
  db_object_name : String
  db_object_synonyms : String
  db_object_type : String
  taxon : String
  date : String
  assigned_by : String
  annotation_extension : String
  gene_product_form_id : String
}

///|
/// Construct a GafRecord with all 17 fields.
pub fn GafRecord::new(
  db~ : String,
  db_object_id~ : String,
  db_object_symbol~ : String,
  qualifier~ : String,
  go_id~ : String,
  db_reference~ : String,
  evidence_code~ : String,
  with_from~ : String,
  aspect~ : String,
  db_object_name~ : String,
  db_object_synonyms~ : String,
  db_object_type~ : String,
  taxon~ : String,
  date~ : String,
  assigned_by~ : String,
  annotation_extension~ : String,
  gene_product_form_id~ : String,
) -> GafRecord {
  GafRecord::{
    db,
    db_object_id,
    db_object_symbol,
    qualifier,
    go_id,
    db_reference,
    evidence_code,
    with_from,
    aspect,
    db_object_name,
    db_object_synonyms,
    db_object_type,
    taxon,
    date,
    assigned_by,
    annotation_extension,
    gene_product_form_id,
  }
}

///|
/// Get the DB (source database) column.
pub fn GafRecord::db(self : GafRecord) -> String {
  self.db
}

///|
/// Get the DB Object ID column.
pub fn GafRecord::db_object_id(self : GafRecord) -> String {
  self.db_object_id
}

///|
/// Get the DB Object Symbol (gene/protein name) column.
pub fn GafRecord::db_object_symbol(self : GafRecord) -> String {
  self.db_object_symbol
}

///|
/// Get the Qualifier column (e.g. "NOT", "enables").
pub fn GafRecord::qualifier(self : GafRecord) -> String {
  self.qualifier
}

///|
/// Get the GO ID column (e.g. "GO:0003674").
pub fn GafRecord::go_id(self : GafRecord) -> String {
  self.go_id
}

///|
/// Get the DB Reference column (pipe-separated, e.g. "PMID:12345").
pub fn GafRecord::db_reference(self : GafRecord) -> String {
  self.db_reference
}

///|
/// Get the Evidence Code column (e.g. "IDA", "IEA").
pub fn GafRecord::evidence_code(self : GafRecord) -> String {
  self.evidence_code
}

///|
/// Get the With/From column (pipe-separated supporting evidence).
pub fn GafRecord::with_from(self : GafRecord) -> String {
  self.with_from
}

///|
/// Get the Aspect column ("F", "P", or "C").
pub fn GafRecord::aspect(self : GafRecord) -> String {
  self.aspect
}

///|
/// Get the DB Object Name column (full name).
pub fn GafRecord::db_object_name(self : GafRecord) -> String {
  self.db_object_name
}

///|
/// Get the DB Object Synonyms column (pipe-separated).
pub fn GafRecord::db_object_synonyms(self : GafRecord) -> String {
  self.db_object_synonyms
}

///|
/// Get the DB Object Type column (e.g. "protein", "gene").
pub fn GafRecord::db_object_type(self : GafRecord) -> String {
  self.db_object_type
}

///|
/// Get the Taxon column (e.g. "taxon:9606|taxon:9606").
pub fn GafRecord::taxon(self : GafRecord) -> String {
  self.taxon
}

///|
/// Get the Date column (YYYYMMDD format).
pub fn GafRecord::date(self : GafRecord) -> String {
  self.date
}

///|
/// Get the Assigned By column (e.g. "UniProt").
pub fn GafRecord::assigned_by(self : GafRecord) -> String {
  self.assigned_by
}

///|
/// Get the Annotation Extension column (pipe-separated).
pub fn GafRecord::annotation_extension(self : GafRecord) -> String {
  self.annotation_extension
}

///|
/// Get the Gene Product Form ID column.
pub fn GafRecord::gene_product_form_id(self : GafRecord) -> String {
  self.gene_product_form_id
}

///|
/// Return the GO ID without the "GO:" prefix (e.g. "0003674").
pub fn GafRecord::go_id_short(self : GafRecord) -> String {
  let id = self.go_id
  if starts_with(id, 0, "GO:") {
    id[3:id.length()].to_owned()
  } else {
    id
  }
}

///|
/// Split the qualifier field on "|" into an array of strings.
pub fn GafRecord::qualifiers(self : GafRecord) -> Array[String] {
  goa_split_pipe(self.qualifier)
}

///|
/// Split the synonyms field on "|" into an array of strings.
pub fn GafRecord::synonyms(self : GafRecord) -> Array[String] {
  goa_split_pipe(self.db_object_synonyms)
}

///|
/// Split the DB Reference field on "|" into an array of strings.
pub fn GafRecord::references(self : GafRecord) -> Array[String] {
  goa_split_pipe(self.db_reference)
}

///|
/// Extract numeric NCBI taxon IDs from the taxon column.
/// For "taxon:9606|taxon:9606" returns ["9606", "9606"].
pub fn GafRecord::taxon_ids(self : GafRecord) -> Array[String] {
  let result : Array[String] = Array::new()
  let parts = split_by_char(self.taxon, '|'.to_int())
  let mut i = 0
  while i < parts.length() {
    let p = trim(parts[i])
    let id = if starts_with(p, 0, "taxon:") {
      p[6:p.length()].to_owned()
    } else {
      p
    }
    if id.length() > 0 {
      result.push(id)
    }
    i = i + 1
  }
  result
}

///|
/// Format a GafRecord as a single-line readable representation.
pub fn GafRecord::to_string(self : GafRecord) -> String {
  self.db +
  ":" +
  self.db_object_id +
  " | " +
  self.db_object_symbol +
  " | " +
  self.go_id +
  " (" +
  self.evidence_code +
  ") | " +
  self.aspect +
  " | " +
  self.taxon
}

///|
/// Container holding all GAF records and header metadata.
pub struct GoaDatabase {
  records : Array[GafRecord]
  version : String
  created_by : String
}

///|
/// Create an empty GoaDatabase.
pub fn GoaDatabase::new() -> GoaDatabase {
  GoaDatabase::{ records: Array::new(), version: "", created_by: "" }
}

///|
/// Get all GAF records.
pub fn GoaDatabase::records(self : GoaDatabase) -> Array[GafRecord] {
  self.records
}

///|
/// Get the GAF version string (e.g. "2.2").
pub fn GoaDatabase::version(self : GoaDatabase) -> String {
  self.version
}

///|
/// Get the name of the tool/source that created the GAF file.
pub fn GoaDatabase::created_by(self : GoaDatabase) -> String {
  self.created_by
}

///|
/// Number of records in the database.
pub fn GoaDatabase::n_records(self : GoaDatabase) -> Int {
  self.records.length()
}

// ============================================================================
// Parsing
// ============================================================================

///|
/// Find the first occurrence of a character in a string. Returns -1 if not found.
fn goa_find_char(s : String, ch : Int) -> Int {
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
/// Split a pipe-separated string into trimmed, non-empty parts.
fn goa_split_pipe(s : String) -> Array[String] {
  let result : Array[String] = Array::new()
  if s.length() == 0 {
    return result
  }
  let parts = split_by_char(s, '|'.to_int())
  let mut i = 0
  while i < parts.length() {
    let p = trim(parts[i])
    if p.length() > 0 {
      result.push(p)
    }
    i = i + 1
  }
  result
}

///|
/// Parse a GAF header line of the form "!key: value".
/// Returns (key, value) or None if the line is not a valid header.
pub fn goa_parse_header_line(line : String) -> (String, String)? {
  let trimmed = trim(line)
  if !starts_with(trimmed, 0, "!") {
    return None
  }
  let rest = trimmed[1:trimmed.length()].to_owned()
  let colon_pos = goa_find_char(rest, ':'.to_int())
  if colon_pos < 0 {
    return None
  }
  let key = trim(rest[0:colon_pos].to_owned())
  let value = trim(rest[colon_pos + 1:rest.length()].to_owned())
  Some((key, value))
}

///|
/// Parse a single GAF data line into a GafRecord.
/// Returns None for blank lines or header lines (starting with "!").
/// Lines with fewer than 17 columns are right-padded with empty strings.
pub fn goa_parse_line(line : String) -> GafRecord? {
  let trimmed = trim(line)
  if trimmed.length() == 0 {
    return None
  }
  if starts_with(trimmed, 0, "!") {
    return None
  }
  let fields = split_by_char(trimmed, '\t'.to_int())
  while fields.length() < 17 {
    fields.push("")
  }
  Some(
    GafRecord::new(
      db=fields[0],
      db_object_id=fields[1],
      db_object_symbol=fields[2],
      qualifier=fields[3],
      go_id=fields[4],
      db_reference=fields[5],
      evidence_code=fields[6],
      with_from=fields[7],
      aspect=fields[8],
      db_object_name=fields[9],
      db_object_synonyms=fields[10],
      db_object_type=fields[11],
      taxon=fields[12],
      date=fields[13],
      assigned_by=fields[14],
      annotation_extension=fields[15],
      gene_product_form_id=fields[16],
    ),
  )
}

///|
/// Parse full GAF text content into a GoaDatabase.
/// Header lines "!gaf-version:" and "!generated-by:" populate the database
/// metadata; all data lines are parsed into records.
pub fn goa_parse(content : String) -> GoaDatabase {
  let records : Array[GafRecord] = Array::new()
  let mut version = ""
  let mut created_by = ""
  let lines = split_by_char(content, '\n'.to_int())
  let mut i = 0
  while i < lines.length() {
    let line = trim(lines[i])
    if line.length() == 0 {
      i = i + 1
      continue
    }
    if starts_with(line, 0, "!") {
      match goa_parse_header_line(line) {
        Some((key, value)) => {
          if key == "gaf-version" {
            version = value
          } else if key == "generated-by" {
            created_by = value
          }
        }
        None => ()
      }
      i = i + 1
      continue
    }
    match goa_parse_line(line) {
      Some(record) => records.push(record)
      None => ()
    }
    i = i + 1
  }
  GoaDatabase::{ records, version, created_by }
}

// ============================================================================
// Filtering
// ============================================================================

///|
/// Return all records matching the given GO ID (e.g. "GO:0003674").
pub fn goa_filter_by_go_id(db : GoaDatabase, go_id : String) -> Array[GafRecord] {
  let result : Array[GafRecord] = Array::new()
  let recs = db.records
  let mut i = 0
  while i < recs.length() {
    if recs[i].go_id == go_id {
      result.push(recs[i])
    }
    i = i + 1
  }
  result
}

///|
/// Return all records whose aspect matches the given GafAspect.
pub fn goa_filter_by_aspect(
  db : GoaDatabase,
  aspect : GafAspect,
) -> Array[GafRecord] {
  let target = aspect.to_string()
  let result : Array[GafRecord] = Array::new()
  let recs = db.records
  let mut i = 0
  while i < recs.length() {
    if recs[i].aspect == target {
      result.push(recs[i])
    }
    i = i + 1
  }
  result
}

///|
/// Return all records with the given evidence code (e.g. "IDA").
pub fn goa_filter_by_evidence(
  db : GoaDatabase,
  evidence : String,
) -> Array[GafRecord] {
  let result : Array[GafRecord] = Array::new()
  let recs = db.records
  let mut i = 0
  while i < recs.length() {
    if recs[i].evidence_code == evidence {
      result.push(recs[i])
    }
    i = i + 1
  }
  result
}

///|
/// Return all records annotated to the given numeric taxon ID
/// (e.g. "9606" matches a record whose taxon column is "taxon:9606|...").
pub fn goa_filter_by_taxon(
  db : GoaDatabase,
  taxon_id : String,
) -> Array[GafRecord] {
  let result : Array[GafRecord] = Array::new()
  let recs = db.records
  let mut i = 0
  while i < recs.length() {
    let ids = recs[i].taxon_ids()
    let mut j = 0
    let mut found = false
    while j < ids.length() {
      if ids[j] == taxon_id {
        found = true
        break
      }
      j = j + 1
    }
    if found {
      result.push(recs[i])
    }
    i = i + 1
  }
  result
}

///|
/// Return all records with the given DB Object ID (e.g. "Q12345").
pub fn goa_filter_by_db_object_id(
  db : GoaDatabase,
  object_id : String,
) -> Array[GafRecord] {
  let result : Array[GafRecord] = Array::new()
  let recs = db.records
  let mut i = 0
  while i < recs.length() {
    if recs[i].db_object_id == object_id {
      result.push(recs[i])
    }
    i = i + 1
  }
  result
}

// ============================================================================
// Unique value extraction
// ============================================================================

///|
/// Collect all unique GO IDs found across the database, preserving first-seen order.
pub fn goa_unique_go_ids(db : GoaDatabase) -> Array[String] {
  let result : Array[String] = Array::new()
  let seen : Map[String, Bool] = Map([], capacity=0)
  let recs = db.records
  let mut i = 0
  while i < recs.length() {
    let gid = recs[i].go_id
    if !seen.contains(gid) {
      seen[gid] = true
      result.push(gid)
    }
    i = i + 1
  }
  result
}

///|
/// Collect all unique evidence codes found across the database.
pub fn goa_unique_evidence_codes(db : GoaDatabase) -> Array[String] {
  let result : Array[String] = Array::new()
  let seen : Map[String, Bool] = Map([], capacity=0)
  let recs = db.records
  let mut i = 0
  while i < recs.length() {
    let code = recs[i].evidence_code
    if !seen.contains(code) {
      seen[code] = true
      result.push(code)
    }
    i = i + 1
  }
  result
}

///|
/// Collect all unique numeric taxon IDs found across the database.
pub fn goa_unique_taxon_ids(db : GoaDatabase) -> Array[String] {
  let result : Array[String] = Array::new()
  let seen : Map[String, Bool] = Map([], capacity=0)
  let recs = db.records
  let mut i = 0
  while i < recs.length() {
    let ids = recs[i].taxon_ids()
    let mut j = 0
    while j < ids.length() {
      let id = ids[j]
      if !seen.contains(id) {
        seen[id] = true
        result.push(id)
      }
      j = j + 1
    }
    i = i + 1
  }
  result
}

// ============================================================================
// Counting / aggregation
// ============================================================================

///|
/// Count records grouped by aspect code. Keys are the raw aspect strings
/// ("F", "P", "C", ...).
pub fn goa_count_by_aspect(db : GoaDatabase) -> Map[String, Int] {
  let counts : Map[String, Int] = Map([], capacity=0)
  let recs = db.records
  let mut i = 0
  while i < recs.length() {
    let aspect = recs[i].aspect
    if counts.contains(aspect) {
      counts[aspect] = counts.get(aspect).unwrap() + 1
    } else {
      counts[aspect] = 1
    }
    i = i + 1
  }
  counts
}

///|
/// Count records grouped by evidence code.
pub fn goa_count_by_evidence(db : GoaDatabase) -> Map[String, Int] {
  let counts : Map[String, Int] = Map([], capacity=0)
  let recs = db.records
  let mut i = 0
  while i < recs.length() {
    let code = recs[i].evidence_code
    if counts.contains(code) {
      counts[code] = counts.get(code).unwrap() + 1
    } else {
      counts[code] = 1
    }
    i = i + 1
  }
  counts
}

///|
/// Build a human-readable summary string of the database contents.
pub fn goa_to_summary(db : GoaDatabase) -> String {
  let mut s = "GOA Database Summary\n"
  s = s + "===================\n"
  s = s + "Total records: " + db.n_records().to_string() + "\n"
  s = s + "GAF version: " + db.version + "\n"
  s = s + "Created by: " + db.created_by + "\n"
  s = s + "\nCounts by aspect:\n"
  let aspect_counts = goa_count_by_aspect(db)
  let aspect_entries = aspect_counts.iter().collect()
  let mut i = 0
  while i < aspect_entries.length() {
    let (aspect, count) = aspect_entries[i]
    let desc = GafAspect::description(GafAspect::from_string(aspect))
    s = s + "  " + aspect + " (" + desc + "): " + count.to_string() + "\n"
    i = i + 1
  }
  s = s + "\nCounts by evidence code:\n"
  let ev_counts = goa_count_by_evidence(db)
  let ev_entries = ev_counts.iter().collect()
  let mut j = 0
  while j < ev_entries.length() {
    let (code, count) = ev_entries[j]
    s = s + "  " + code + ": " + count.to_string() + "\n"
    j = j + 1
  }
  s = s + "\nUnique GO IDs: " + goa_unique_go_ids(db).length().to_string() + "\n"
  s = s + "Unique evidence codes: " +
    goa_unique_evidence_codes(db).length().to_string() + "\n"
  s = s + "Unique taxon IDs: " + goa_unique_taxon_ids(db).length().to_string() + "\n"
  s
}

// ============================================================================
// Sample data
// ============================================================================

///|
/// Return a small sample GoaDatabase with ~10 records covering different
/// aspects (F/P/C), evidence codes (IDA/IEA/EXP/ISS/IPI), and taxa
/// (9606 human, 4932 yeast, 7227 fly).
pub fn goa_sample_database() -> GoaDatabase {
  let records : Array[GafRecord] = Array::new()
  // Record 1: PROT1, molecular function, IDA, human
  records.push(
    GafRecord::new(
      db="UniProtKB",
      db_object_id="Q12345",
      db_object_symbol="PROT1",
      qualifier="enables",
      go_id="GO:0003674",
      db_reference="PMID:12345",
      evidence_code="IDA",
      with_from="GO:0005515",
      aspect="F",
      db_object_name="Protein 1",
      db_object_synonyms="P1|PROT-1",
      db_object_type="protein",
      taxon="taxon:9606|taxon:9606",
      date="20210115",
      assigned_by="UniProt",
      annotation_extension="",
      gene_product_form_id="",
    ),
  )
  // Record 2: PROT1, biological process, IDA, human
  records.push(
    GafRecord::new(
      db="UniProtKB",
      db_object_id="Q12345",
      db_object_symbol="PROT1",
      qualifier="involved_in",
      go_id="GO:0008150",
      db_reference="PMID:12345",
      evidence_code="IDA",
      with_from="",
      aspect="P",
      db_object_name="Protein 1",
      db_object_synonyms="P1|PROT-1",
      db_object_type="protein",
      taxon="taxon:9606|taxon:9606",
      date="20210115",
      assigned_by="UniProt",
      annotation_extension="",
      gene_product_form_id="",
    ),
  )
  // Record 3: PROT1, cellular component, IDA, human
  records.push(
    GafRecord::new(
      db="UniProtKB",
      db_object_id="Q12345",
      db_object_symbol="PROT1",
      qualifier="located_in",
      go_id="GO:0005737",
      db_reference="PMID:12345",
      evidence_code="IDA",
      with_from="",
      aspect="C",
      db_object_name="Protein 1",
      db_object_synonyms="P1|PROT-1",
      db_object_type="protein",
      taxon="taxon:9606|taxon:9606",
      date="20210115",
      assigned_by="UniProt",
      annotation_extension="",
      gene_product_form_id="",
    ),
  )
  // Record 4: TP53, molecular function, EXP, human
  records.push(
    GafRecord::new(
      db="UniProtKB",
      db_object_id="P04637",
      db_object_symbol="TP53",
      qualifier="enables",
      go_id="GO:0003700",
      db_reference="PMID:67890",
      evidence_code="EXP",
      with_from="GO:0005634",
      aspect="F",
      db_object_name="Cellular tumor antigen p53",
      db_object_synonyms="p53|LFS1",
      db_object_type="protein",
      taxon="taxon:9606|taxon:9606",
      date="20210220",
      assigned_by="UniProt",
      annotation_extension="",
      gene_product_form_id="",
    ),
  )
  // Record 5: TP53, biological process, IEA, human
  records.push(
    GafRecord::new(
      db="UniProtKB",
      db_object_id="P04637",
      db_object_symbol="TP53",
      qualifier="involved_in",
      go_id="GO:0006915",
      db_reference="PMID:67890",
      evidence_code="IEA",
      with_from="",
      aspect="P",
      db_object_name="Cellular tumor antigen p53",
      db_object_synonyms="p53|LFS1",
      db_object_type="protein",
      taxon="taxon:9606|taxon:9606",
      date="20210301",
      assigned_by="InterPro",
      annotation_extension="",
      gene_product_form_id="",
    ),
  )
  // Record 6: YFG1, molecular function, ISS, yeast
  records.push(
    GafRecord::new(
      db="UniProtKB",
      db_object_id="P00549",
      db_object_symbol="YFG1",
      qualifier="enables",
      go_id="GO:0005524",
      db_reference="PMID:11111",
      evidence_code="ISS",
      with_from="UniProtKB:Q12345",
      aspect="F",
      db_object_name="Your favorite gene 1",
      db_object_synonyms="YFG-1",
      db_object_type="protein",
      taxon="taxon:4932|taxon:4932",
      date="20210410",
      assigned_by="SGD",
      annotation_extension="",
      gene_product_form_id="",
    ),
  )
  // Record 7: YFG1, cellular component, IDA, yeast
  records.push(
    GafRecord::new(
      db="UniProtKB",
      db_object_id="P00549",
      db_object_symbol="YFG1",
      qualifier="located_in",
      go_id="GO:0005634",
      db_reference="PMID:11111",
      evidence_code="IDA",
      with_from="",
      aspect="C",
      db_object_name="Your favorite gene 1",
      db_object_synonyms="YFG-1",
      db_object_type="protein",
      taxon="taxon:4932|taxon:4932",
      date="20210410",
      assigned_by="SGD",
      annotation_extension="",
      gene_product_form_id="",
    ),
  )
  // Record 8: HSP70, biological process, IEA, fly
  records.push(
    GafRecord::new(
      db="UniProtKB",
      db_object_id="P02825",
      db_object_symbol="Hsp70",
      qualifier="involved_in",
      go_id="GO:0006457",
      db_reference="PMID:22222",
      evidence_code="IEA",
      with_from="",
      aspect="P",
      db_object_name="Heat shock protein 70",
      db_object_synonyms="Hsp-70|Hsp70",
      db_object_type="protein",
      taxon="taxon:7227|taxon:7227",
      date="20210505",
      assigned_by="FlyBase",
      annotation_extension="",
      gene_product_form_id="",
    ),
  )
  // Record 9: HSP70, cellular component, IDA, fly
  records.push(
    GafRecord::new(
      db="UniProtKB",
      db_object_id="P02825",
      db_object_symbol="Hsp70",
      qualifier="located_in",
      go_id="GO:0005829",
      db_reference="PMID:22222",
      evidence_code="IDA",
      with_from="",
      aspect="C",
      db_object_name="Heat shock protein 70",
      db_object_synonyms="Hsp-70|Hsp70",
      db_object_type="protein",
      taxon="taxon:7227|taxon:7227",
      date="20210505",
      assigned_by="FlyBase",
      annotation_extension="",
      gene_product_form_id="",
    ),
  )
  // Record 10: RIB1, molecular function, IPI, human
  records.push(
    GafRecord::new(
      db="UniProtKB",
      db_object_id="P05388",
      db_object_symbol="RPL5",
      qualifier="enables",
      go_id="GO:0003735",
      db_reference="PMID:33333",
      evidence_code="IPI",
      with_from="GO:0005515",
      aspect="F",
      db_object_name="60S ribosomal protein L5",
      db_object_synonyms="RPL-5|RIB1",
      db_object_type="protein",
      taxon="taxon:9606|taxon:9606",
      date="20210601",
      assigned_by="UniProt",
      annotation_extension="",
      gene_product_form_id="",
    ),
  )
  GoaDatabase::{ records, version: "2.2", created_by: "UniProt" }
}

// ============================================================================
// Writing
// ============================================================================

///|
/// Serialize a single GafRecord as a tab-separated GAF line (no trailing newline).
pub fn goa_record_to_gaf_line(record : GafRecord) -> String {
  record.db +
  "\t" +
  record.db_object_id +
  "\t" +
  record.db_object_symbol +
  "\t" +
  record.qualifier +
  "\t" +
  record.go_id +
  "\t" +
  record.db_reference +
  "\t" +
  record.evidence_code +
  "\t" +
  record.with_from +
  "\t" +
  record.aspect +
  "\t" +
  record.db_object_name +
  "\t" +
  record.db_object_synonyms +
  "\t" +
  record.db_object_type +
  "\t" +
  record.taxon +
  "\t" +
  record.date +
  "\t" +
  record.assigned_by +
  "\t" +
  record.annotation_extension +
  "\t" +
  record.gene_product_form_id
}

///|
/// Serialize an entire GoaDatabase as GAF text, including header lines.
pub fn goa_to_gaf(db : GoaDatabase) -> String {
  let mut s = "!gaf-version: " + db.version + "\n"
  if db.created_by.length() > 0 {
    s = s + "!generated-by: " + db.created_by + "\n"
  }
  let recs = db.records
  let mut i = 0
  while i < recs.length() {
    s = s + goa_record_to_gaf_line(recs[i]) + "\n"
    i = i + 1
  }
  s
}
