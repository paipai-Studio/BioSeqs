///|
/// Bio.SearchIO.InterproscanIO - InterProScan output parser.
///
/// References:
///   - Biopython Bio.SearchIO.InterproscanIO module
///   - InterProScan documentation: https://interproscan-docs.readthedocs.io/
///   - Jones P et al. (2014) InterProScan 5: genome-scale protein function
///     classification. Bioinformatics 30(9):1236-40.
///
/// InterProScan scans protein sequences against InterPro's protein signature
/// databases (Pfam, SMART, PROSITE, etc.) and reports matches in TSV or XML
/// format. This module parses the 14-column TSV output format.
///
/// TSV columns:
///   1. protein ID
///   2. MD5 hash of sequence
///   3. sequence length
///   4. analysis (database name, e.g. "Pfam")
///   5. signature accession
///   6. signature description
///   7. start position
///   8. stop position
///   9. score
///  10. status (e.g. "T" for true)
///  11. date
///  12. IPR accession
///  13. IPR description
///  14. GO terms (pipe-separated)
///
/// Lines starting with "#" are comments. Lines where the analysis column is
/// "No hits" indicate no matches and are skipped. Missing IPR and GO term
/// fields are represented as "-" in the file and normalized to empty values.
///
/// All public identifiers are prefixed with `interproscan_` to avoid conflicts.
/// Internal helpers use the `ips_` prefix.

// ============================================================================
// Data structure
// ============================================================================

///|
/// A single InterProScan TSV record with all 14 columns.
pub struct InterproScanRecord {
  protein_id : String
  md5 : String
  seq_length : Int
  analysis : String
  signature_acc : String
  signature_desc : String
  start_pos : Int
  stop_pos : Int
  score : Double
  status : String
  date : String
  ipr_acc : String
  ipr_desc : String
  go_terms : Array[String]
}

///|
/// Construct an InterproScanRecord with all fields.
pub fn InterproScanRecord::new(
  protein_id~ : String,
  md5~ : String,
  seq_length~ : Int,
  analysis~ : String,
  signature_acc~ : String,
  signature_desc~ : String,
  start_pos~ : Int,
  stop_pos~ : Int,
  score~ : Double,
  status~ : String,
  date~ : String,
  ipr_acc~ : String,
  ipr_desc~ : String,
  go_terms~ : Array[String],
) -> InterproScanRecord {
  InterproScanRecord::{
    protein_id,
    md5,
    seq_length,
    analysis,
    signature_acc,
    signature_desc,
    start_pos,
    stop_pos,
    score,
    status,
    date,
    ipr_acc,
    ipr_desc,
    go_terms,
  }
}

// ============================================================================
// Accessors
// ============================================================================

///|
/// Get the protein ID column.
pub fn InterproScanRecord::protein_id(self : InterproScanRecord) -> String {
  self.protein_id
}

///|
/// Get the MD5 hash column.
pub fn InterproScanRecord::md5(self : InterproScanRecord) -> String {
  self.md5
}

///|
/// Get the sequence length column.
pub fn InterproScanRecord::seq_length(self : InterproScanRecord) -> Int {
  self.seq_length
}

///|
/// Get the analysis (database name) column.
pub fn InterproScanRecord::analysis(self : InterproScanRecord) -> String {
  self.analysis
}

///|
/// Get the signature accession column.
pub fn InterproScanRecord::signature_acc(self : InterproScanRecord) -> String {
  self.signature_acc
}

///|
/// Get the signature description column.
pub fn InterproScanRecord::signature_desc(self : InterproScanRecord) -> String {
  self.signature_desc
}

///|
/// Get the start position column.
pub fn InterproScanRecord::start_pos(self : InterproScanRecord) -> Int {
  self.start_pos
}

///|
/// Get the stop position column.
pub fn InterproScanRecord::stop_pos(self : InterproScanRecord) -> Int {
  self.stop_pos
}

///|
/// Get the score column.
pub fn InterproScanRecord::score(self : InterproScanRecord) -> Double {
  self.score
}

///|
/// Get the status column.
pub fn InterproScanRecord::status(self : InterproScanRecord) -> String {
  self.status
}

///|
/// Get the date column.
pub fn InterproScanRecord::date(self : InterproScanRecord) -> String {
  self.date
}

///|
/// Get the IPR accession column.
pub fn InterproScanRecord::ipr_acc(self : InterproScanRecord) -> String {
  self.ipr_acc
}

///|
/// Get the IPR description column.
pub fn InterproScanRecord::ipr_desc(self : InterproScanRecord) -> String {
  self.ipr_desc
}

///|
/// Get the GO terms array.
pub fn InterproScanRecord::go_terms(self : InterproScanRecord) -> Array[String] {
  self.go_terms
}

///|
/// Human-readable string representation of the record.
pub fn InterproScanRecord::to_string(self : InterproScanRecord) -> String {
  let go_str = ips_join(self.go_terms, "|")
  "InterproScanRecord(protein_id=" +
  self.protein_id +
  ", analysis=" +
  self.analysis +
  ", signature=" +
  self.signature_acc +
  ", start=" +
  self.start_pos.to_string() +
  ", stop=" +
  self.stop_pos.to_string() +
  ", score=" +
  self.score.to_string() +
  ", ipr=" +
  self.ipr_acc +
  ", go_terms=" +
  go_str +
  ")"
}

// ============================================================================
// Internal helpers (ips_ prefix)
// ============================================================================

///|
/// Parse a score string into Double. Returns 0.0 for "-", "", or "NaN".
fn ips_parse_score(s : String) -> Double {
  if s == "-" || s == "" || s == "NaN" {
    return 0.0
  }
  match parse_double(s) {
    Some(d) => d
    None => 0.0
  }
}

///|
/// Normalize a field value: "-" becomes "" (for missing IPR/GO fields).
fn ips_normalize(s : String) -> String {
  if s == "-" {
    ""
  } else {
    s
  }
}

///|
/// Split a pipe-separated GO terms string into an array.
/// Returns an empty array for "-", "", or whitespace-only values.
/// Empty segments (e.g. from trailing "|") are filtered out.
fn ips_split_go_terms(s : String) -> Array[String] {
  let normalized = ips_normalize(s)
  if normalized.length() == 0 {
    return Array::new()
  }
  let parts = split_by_char(normalized, '|'.to_int())
  let result : Array[String] = Array::new()
  for p in parts {
    if p.length() > 0 {
      result.push(p)
    }
  }
  result
}

///|
/// Join an array of strings with a separator.
fn ips_join(arr : Array[String], sep : String) -> String {
  if arr.length() == 0 {
    return ""
  }
  let mut s = arr[0]
  for i in 1..<arr.length() {
    s = s + sep + arr[i]
  }
  s
}

///|
/// Parse a single InterProScan TSV line into a record.
/// Returns None for blank lines, comment lines (starting with "#"),
/// and "No hits" lines. Lines with fewer than 14 columns are padded.
fn ips_parse_line(line : String) -> InterproScanRecord? {
  let trimmed = trim(line)
  if trimmed.length() == 0 {
    return None
  }
  if starts_with(trimmed, 0, "#") {
    return None
  }
  let fields = split_by_tab(trimmed)
  // Lines where the analysis column (index 3) is "No hits" indicate no matches.
  if fields.length() >= 4 && fields[3] == "No hits" {
    return None
  }
  // Pad to 14 fields if necessary.
  while fields.length() < 14 {
    fields.push("")
  }
  let protein_id = fields[0]
  let md5 = fields[1]
  let seq_length = parse_int(fields[2])
  let analysis = fields[3]
  let signature_acc = fields[4]
  let signature_desc = fields[5]
  let start_pos = parse_int(fields[6])
  let stop_pos = parse_int(fields[7])
  let score = ips_parse_score(fields[8])
  let status = fields[9]
  let date = fields[10]
  let ipr_acc = ips_normalize(fields[11])
  let ipr_desc = ips_normalize(fields[12])
  let go_terms = ips_split_go_terms(fields[13])
  Some(
    InterproScanRecord::{
      protein_id,
      md5,
      seq_length,
      analysis,
      signature_acc,
      signature_desc,
      start_pos,
      stop_pos,
      score,
      status,
      date,
      ipr_acc,
      ipr_desc,
      go_terms,
    },
  )
}

// ============================================================================
// Parsing
// ============================================================================

///|
/// Parse InterProScan TSV output into an array of records.
/// This is the main entry point and defaults to TSV format.
pub fn parse_interproscan(content : String) -> Array[InterproScanRecord] {
  parse_interproscan_tsv(content)
}

///|
/// Explicit TSV parser for InterProScan output.
/// Skips comment lines (starting with "#"), blank lines, and "No hits" lines.
pub fn parse_interproscan_tsv(content : String) -> Array[InterproScanRecord] {
  let records : Array[InterproScanRecord] = Array::new()
  let lines = split_by_char(content, '\n'.to_int())
  for line in lines {
    match ips_parse_line(line) {
      Some(rec) => records.push(rec)
      None => ()
    }
  }
  records
}

// ============================================================================
// Sample data
// ============================================================================

///|
/// Return sample InterProScan TSV data for testing and demonstration.
/// Contains 8 records spanning 3 proteins, 5 analyses, and 8 unique signatures.
pub fn interproscan_sample_data() -> String {
  let t = "\t"
  "# InterProScan output (sample)\n" +
  "# version 5.0\n" +
  "sp|P12345|PROT_HUMAN" + t + "abc123def456" + t + "500" + t + "Pfam" + t + "PF00001" + t + "p53 DNA-binding domain" + t + "1" + t + "100" + t + "150.5" + t + "T" + t + "01-Jan-2024" + t + "IPR000001" + t + "p53, DNA-binding domain" + t + "GO:0003674|GO:0008150\n" +
  "sp|P12345|PROT_HUMAN" + t + "abc123def456" + t + "500" + t + "SMART" + t + "SM00100" + t + "p53 domain" + t + "1" + t + "100" + t + "120.0" + t + "T" + t + "01-Jan-2024" + t + "IPR000001" + t + "p53, DNA-binding domain" + t + "GO:0003674\n" +
  "sp|P12345|PROT_HUMAN" + t + "abc123def456" + t + "500" + t + "PROSITEPATTERNS" + t + "PS00001" + t + "p53 pattern" + t + "50" + t + "60" + t + "0.0" + t + "T" + t + "01-Jan-2024" + t + "-" + t + "-" + t + "-\n" +
  "sp|Q9Y6K9|XYZ_HUMAN" + t + "def789ghi012" + t + "300" + t + "Profile" + t + "PS50122" + t + "Helicase C-terminal" + t + "10" + t + "250" + t + "200.3" + t + "T" + t + "02-Jan-2024" + t + "IPR000002" + t + "Helicase, C-terminal" + t + "GO:0005524|GO:0006280|GO:0003676\n" +
  "sp|Q9Y6K9|XYZ_HUMAN" + t + "def789ghi012" + t + "300" + t + "Pfam" + t + "PF00002" + t + "Kinase domain" + t + "5" + t + "150" + t + "180.0" + t + "T" + t + "02-Jan-2024" + t + "IPR000003" + t + "Kinase domain" + t + "GO:0004672\n" +
  "sp|O15143|ABC_HUMAN" + t + "ghi345jkl678" + t + "600" + t + "Pfam" + t + "PF00003" + t + "ABC transporter" + t + "100" + t + "300" + t + "110.5" + t + "T" + t + "03-Jan-2024" + t + "IPR000004" + t + "ABC transporter-like" + t + "GO:0005524|GO:0016887\n" +
  "sp|O15143|ABC_HUMAN" + t + "ghi345jkl678" + t + "600" + t + "HAMMER" + t + "HM00001" + t + "Hammerhead ribozyme" + t + "1" + t + "50" + t + "90.0" + t + "T" + t + "03-Jan-2024" + t + "IPR000005" + t + "Hammerhead ribozyme" + t + "-\n" +
  "sp|O15143|ABC_HUMAN" + t + "ghi345jkl678" + t + "600" + t + "SMART" + t + "SM00150" + t + "DnaJ domain" + t + "400" + t + "500" + t + "85.0" + t + "T" + t + "03-Jan-2024" + t + "IPR000006" + t + "DnaJ domain" + t + "GO:0044183\n"
}

// ============================================================================
// Filtering
// ============================================================================

///|
/// Filter records by analysis (database) name, e.g. "Pfam", "SMART".
pub fn interproscan_filter_by_analysis(
  records : Array[InterproScanRecord],
  analysis : String,
) -> Array[InterproScanRecord] {
  let result : Array[InterproScanRecord] = Array::new()
  for r in records {
    if r.analysis == analysis {
      result.push(r)
    }
  }
  result
}

///|
/// Filter records by protein ID.
pub fn interproscan_filter_by_protein(
  records : Array[InterproScanRecord],
  protein_id : String,
) -> Array[InterproScanRecord] {
  let result : Array[InterproScanRecord] = Array::new()
  for r in records {
    if r.protein_id == protein_id {
      result.push(r)
    }
  }
  result
}

// ============================================================================
// Unique value extraction
// ============================================================================

///|
/// Collect all unique protein IDs, preserving first-seen order.
pub fn interproscan_get_unique_proteins(
  records : Array[InterproScanRecord],
) -> Array[String] {
  let result : Array[String] = Array::new()
  let seen : Map[String, Bool] = Map([], capacity=0)
  for r in records {
    if !seen.contains(r.protein_id) {
      seen[r.protein_id] = true
      result.push(r.protein_id)
    }
  }
  result
}

///|
/// Collect all unique signature accessions, preserving first-seen order.
pub fn interproscan_get_unique_signatures(
  records : Array[InterproScanRecord],
) -> Array[String] {
  let result : Array[String] = Array::new()
  let seen : Map[String, Bool] = Map([], capacity=0)
  for r in records {
    if !seen.contains(r.signature_acc) {
      seen[r.signature_acc] = true
      result.push(r.signature_acc)
    }
  }
  result
}

///|
/// Collect all unique GO terms across all records, preserving first-seen order.
pub fn interproscan_get_go_terms(
  records : Array[InterproScanRecord],
) -> Array[String] {
  let result : Array[String] = Array::new()
  let seen : Map[String, Bool] = Map([], capacity=0)
  for r in records {
    for g in r.go_terms {
      if !seen.contains(g) {
        seen[g] = true
        result.push(g)
      }
    }
  }
  result
}

// ============================================================================
// Summary
// ============================================================================

///|
/// Build a human-readable summary string of the records.
/// Includes total records, unique proteins, signatures, and analyses.
pub fn interproscan_summary(records : Array[InterproScanRecord]) -> String {
  let total = records.length()
  let proteins = interproscan_get_unique_proteins(records)
  let signatures = interproscan_get_unique_signatures(records)
  let analyses : Array[String] = Array::new()
  let seen : Map[String, Bool] = Map([], capacity=0)
  for r in records {
    if !seen.contains(r.analysis) {
      seen[r.analysis] = true
      analyses.push(r.analysis)
    }
  }
  let mut s = "InterProScan Summary\n"
  s = s + "====================\n"
  s = s + "Total records: " + total.to_string() + "\n"
  s = s + "Unique proteins: " + proteins.length().to_string() + "\n"
  s = s + "Unique signatures: " + signatures.length().to_string() + "\n"
  s = s + "Unique analyses: " + analyses.length().to_string() + "\n"
  s = s + "Analyses: " + ips_join(analyses, ", ") + "\n"
  s
}

// ============================================================================
// Grouping
// ============================================================================

///|
/// Group records by protein ID, preserving first-seen order of proteins.
/// Returns an array of (protein_id, records) tuples.
pub fn interproscan_group_by_protein(
  records : Array[InterproScanRecord],
) -> Array[(String, Array[InterproScanRecord])] {
  let groups : Array[(String, Array[InterproScanRecord])] = Array::new()
  let indices : Map[String, Int] = Map([], capacity=0)
  for r in records {
    let pid = r.protein_id
    if indices.contains(pid) {
      let idx = indices.get(pid).unwrap()
      let (_, arr) = groups[idx]
      arr.push(r)
    } else {
      let idx = groups.length()
      indices[pid] = idx
      let arr : Array[InterproScanRecord] = Array::new()
      arr.push(r)
      groups.push((pid, arr))
    }
  }
  groups
}
