// Bio.SeqIO.ImgtIO - Parser and writer for the IMGT sequence format,
// mirroring BioPython's Bio.SeqIO.ImgtIO.
//
// IMGT is a FASTA-like format used by the international ImMunoGeneTics
// information system. The header line starts with '>' followed by specific
// IMGT fields separated by '|'. For IMGT/HLA there are 13 fields:
//
//   1.  Accession number
//   2.  Sequence ID / allele name
//   3.  Official name
//   4.  Region
//   5.  Date of last update
//   6.  Date of creation
//   7.  Origin species
//   8.  Gene designation
//   9.  IMGT group
//   10. Specific locus
//   11. Sequence length (with 'bp' suffix)
//   12. Nucleotide sequence type (e.g. " genomic-DNA")
//   13. Additional comments
//
// Example header:
//   >HLA00001|A*01:01:01:01|A*01:01:01:01|ORF|2015-03-31|1996-01-19|human|HLA-A|I|MHC|2703 bp| genomic-DNA|
//
// Sequence lines follow the header until the next '>' or end of input.

///|
/// Default line width for wrapping IMGT sequence output (matches FASTA).
let imgt_line_width : Int = 60

// ============================================================================
// ImgtHeader - the 13 pipe-separated fields of an IMGT header
// ============================================================================

///|
/// The 13 pipe-separated fields of an IMGT/HLA header line.
pub struct ImgtHeader {
  accession : String
  seq_id : String
  official_name : String
  region : String
  date_update : String
  date_creation : String
  species : String
  gene : String
  imgt_group : String
  locus : String
  seq_length : String
  seq_type : String
  comments : String
}

///|
/// Construct an ImgtHeader from its 13 fields.
pub fn ImgtHeader::new(
  accession : String,
  seq_id : String,
  official_name : String,
  region : String,
  date_update : String,
  date_creation : String,
  species : String,
  gene : String,
  imgt_group : String,
  locus : String,
  seq_length : String,
  seq_type : String,
  comments : String,
) -> ImgtHeader {
  ImgtHeader::{
    accession,
    seq_id,
    official_name,
    region,
    date_update,
    date_creation,
    species,
    gene,
    imgt_group,
    locus,
    seq_length,
    seq_type,
    comments,
  }
}

///|
/// Field 1: accession number.
pub fn ImgtHeader::accession(self : ImgtHeader) -> String {
  self.accession
}

///|
/// Field 2: sequence ID / allele name.
pub fn ImgtHeader::seq_id(self : ImgtHeader) -> String {
  self.seq_id
}

///|
/// Field 3: official name.
pub fn ImgtHeader::official_name(self : ImgtHeader) -> String {
  self.official_name
}

///|
/// Field 4: region.
pub fn ImgtHeader::region(self : ImgtHeader) -> String {
  self.region
}

///|
/// Field 5: date of last update.
pub fn ImgtHeader::date_update(self : ImgtHeader) -> String {
  self.date_update
}

///|
/// Field 6: date of creation.
pub fn ImgtHeader::date_creation(self : ImgtHeader) -> String {
  self.date_creation
}

///|
/// Field 7: origin species.
pub fn ImgtHeader::species(self : ImgtHeader) -> String {
  self.species
}

///|
/// Field 8: gene designation.
pub fn ImgtHeader::gene(self : ImgtHeader) -> String {
  self.gene
}

///|
/// Field 9: IMGT group.
pub fn ImgtHeader::imgt_group(self : ImgtHeader) -> String {
  self.imgt_group
}

///|
/// Field 10: specific locus.
pub fn ImgtHeader::locus(self : ImgtHeader) -> String {
  self.locus
}

///|
/// Field 11: sequence length (with 'bp' suffix).
pub fn ImgtHeader::seq_length(self : ImgtHeader) -> String {
  self.seq_length
}

///|
/// Field 12: nucleotide sequence type.
pub fn ImgtHeader::seq_type(self : ImgtHeader) -> String {
  self.seq_type
}

///|
/// Field 13: additional comments.
pub fn ImgtHeader::comments(self : ImgtHeader) -> String {
  self.comments
}

///|
/// Render the header as pipe-separated fields (without the leading '>').
pub fn ImgtHeader::to_string(self : ImgtHeader) -> String {
  self.accession +
  "|" +
  self.seq_id +
  "|" +
  self.official_name +
  "|" +
  self.region +
  "|" +
  self.date_update +
  "|" +
  self.date_creation +
  "|" +
  self.species +
  "|" +
  self.gene +
  "|" +
  self.imgt_group +
  "|" +
  self.locus +
  "|" +
  self.seq_length +
  "|" +
  self.seq_type +
  "|" +
  self.comments
}

// ============================================================================
// ImgtRecord - a header plus its sequence
// ============================================================================

///|
/// A single IMGT record: a 13-field header and the associated nucleotide or
/// protein sequence.
pub struct ImgtRecord {
  header : ImgtHeader
  sequence : String
}

///|
/// Construct an ImgtRecord from a header and a sequence string.
pub fn ImgtRecord::new(header : ImgtHeader, sequence : String) -> ImgtRecord {
  ImgtRecord::{ header, sequence }
}

///|
/// The record's IMGT header.
pub fn ImgtRecord::header(self : ImgtRecord) -> ImgtHeader {
  self.header
}

///|
/// The record's sequence string.
pub fn ImgtRecord::sequence(self : ImgtRecord) -> String {
  self.sequence
}

///|
/// The record id: the header's seq_id field.
pub fn ImgtRecord::id(self : ImgtRecord) -> String {
  self.header.seq_id
}

///|
/// The record description: the full header rendered as pipe-separated fields
/// (without the leading '>').
pub fn ImgtRecord::description(self : ImgtRecord) -> String {
  self.header.to_string()
}

///|
/// Human-readable representation of this record.
pub fn ImgtRecord::to_string(self : ImgtRecord) -> String {
  let buf = StringBuilder::new()
  buf.write_string("ImgtRecord(id=")
  buf.write_string(self.header.seq_id)
  buf.write_string(", accession=")
  buf.write_string(self.header.accession)
  buf.write_string(", gene=")
  buf.write_string(self.header.gene)
  buf.write_string(", species=")
  buf.write_string(self.header.species)
  buf.write_string(", length=")
  buf.write_string(self.sequence.length().to_string())
  buf.write_string(")")
  buf.to_string()
}

// ============================================================================
// Parsing
// ============================================================================

///|
/// Split `s` on `delim` and return an array of owned strings. (MoonBit's
/// `String::split` returns an iterator of StringView; this helper materializes
/// the pieces into `Array[String]`.)
fn imgt_split_strings(s : String, delim : String) -> Array[String] {
  let views = s.split(delim).to_array()
  let out : Array[String] = Array::new()
  for i in 0..<views.length() {
    out.push(views[i].to_owned())
  }
  out
}

///|
/// Split an IMGT header line by '|', stripping a leading '>' if present.
/// Returns the array of field strings (in order, possibly fewer or more than
/// 13).
pub fn imgt_parse_line_to_fields(line : String) -> Array[String] {
  let n = line.length()
  let body = if n > 0 && line.unsafe_get(0) == '>'.to_int().to_uint16() {
    line[1:].to_owned()
  } else {
    line
  }
  imgt_split_strings(body, "|")
}

///|
/// Parse an IMGT header line (starting with '>') into an ImgtHeader.
/// Returns None when the line does not start with '>'. Fewer than 13 fields
/// are padded with empty strings; extra fields are ignored.
pub fn imgt_parse_header(line : String) -> ImgtHeader? {
  let n = line.length()
  if n == 0 || line.unsafe_get(0) != '>'.to_int().to_uint16() {
    return None
  }
  let fields = imgt_parse_line_to_fields(line)
  let count = fields.length()
  let acc = imgt_field_or_empty(fields, count, 0)
  let sid = imgt_field_or_empty(fields, count, 1)
  let on = imgt_field_or_empty(fields, count, 2)
  let reg = imgt_field_or_empty(fields, count, 3)
  let du = imgt_field_or_empty(fields, count, 4)
  let dc = imgt_field_or_empty(fields, count, 5)
  let sp = imgt_field_or_empty(fields, count, 6)
  let gn = imgt_field_or_empty(fields, count, 7)
  let ig = imgt_field_or_empty(fields, count, 8)
  let lo = imgt_field_or_empty(fields, count, 9)
  let sl = imgt_field_or_empty(fields, count, 10)
  let st = imgt_field_or_empty(fields, count, 11)
  let cm = imgt_field_or_empty(fields, count, 12)
  Some(
    ImgtHeader::{
      accession: acc,
      seq_id: sid,
      official_name: on,
      region: reg,
      date_update: du,
      date_creation: dc,
      species: sp,
      gene: gn,
      imgt_group: ig,
      locus: lo,
      seq_length: sl,
      seq_type: st,
      comments: cm,
    },
  )
}

///|
/// Return `fields[i]` if `i < count`, else "".
fn imgt_field_or_empty(fields : Array[String], count : Int, i : Int) -> String {
  if i < count {
    fields[i]
  } else {
    ""
  }
}

///|
/// Parse full IMGT file content into an array of ImgtRecord objects. Header
/// lines start with '>'; sequence lines are concatenated until the next
/// header.
pub fn imgt_parse(content : String) -> Array[ImgtRecord] {
  let records = Array::new()
  let lines = content.split("\n").to_array()
  let mut current_header : ImgtHeader? = None
  let mut current_seq = StringBuilder::new()
  let gt : UInt16 = '>'.to_int().to_uint16()
  let cr : UInt16 = '\r'.to_int().to_uint16()
  let sp : UInt16 = ' '.to_int().to_uint16()
  let tab : UInt16 = '\t'.to_int().to_uint16()
  for line in lines {
    let n = line.length()
    let mut trimmed_len = n
    if n > 0 && line.unsafe_get(n - 1) == cr {
      trimmed_len = n - 1
    }
    if trimmed_len == 0 {
      continue
    }
    if line.unsafe_get(0) == gt {
      match current_header {
        Some(h) =>
          records.push(ImgtRecord::{ header: h, sequence: current_seq.to_string() })
        None => ()
      }
      match imgt_parse_header(line[0:trimmed_len].to_owned()) {
        Some(h) => current_header = Some(h)
        None => current_header = None
      }
      current_seq = StringBuilder::new()
    } else {
      for i = 0; i < trimmed_len; i = i + 1 {
        let c = line.unsafe_get(i)
        if c != sp && c != tab {
          current_seq.write_char(c.unsafe_to_char())
        }
      }
    }
  }
  match current_header {
    Some(h) =>
      records.push(ImgtRecord::{ header: h, sequence: current_seq.to_string() })
    None => ()
  }
  records
}

// ============================================================================
// Writing
// ============================================================================

///|
/// Render an ImgtHeader as an IMGT header line: ">field1|field2|...|field13"
/// (without a trailing newline).
pub fn imgt_header_to_string(header : ImgtHeader) -> String {
  ">" + header.to_string()
}

///|
/// Render an ImgtRecord as full IMGT text: the header line followed by the
/// sequence wrapped at `imgt_line_width` characters per line. A trailing
/// newline is included.
pub fn imgt_record_to_string(record : ImgtRecord) -> String {
  let buf = StringBuilder::new()
  buf.write_string(imgt_header_to_string(record.header))
  buf.write_string("\n")
  let seq = record.sequence
  let n = seq.length()
  let mut i = 0
  while i < n {
    let end = if i + imgt_line_width < n { i + imgt_line_width } else { n }
    buf.write_string(seq[i:end].to_owned())
    buf.write_string("\n")
    i = end
  }
  buf.to_string()
}

///|
/// Serialize an array of ImgtRecord objects to IMGT text.
pub fn imgt_write(records : Array[ImgtRecord]) -> String {
  let buf = StringBuilder::new()
  for i in 0..<records.length() {
    buf.write_string(imgt_record_to_string(records[i]))
  }
  buf.to_string()
}

// ============================================================================
// Utility functions
// ============================================================================

///|
/// Number of IMGT records.
pub fn imgt_n_records(records : Array[ImgtRecord]) -> Int {
  records.length()
}

///|
/// Find the first record with the given accession number, or None.
pub fn imgt_get_by_accession(
  records : Array[ImgtRecord],
  accession : String,
) -> ImgtRecord? {
  let mut found : ImgtRecord? = None
  for i in 0..<records.length() {
    if records[i].header.accession == accession {
      found = Some(records[i])
      break
    }
  }
  found
}

///|
/// Return all records whose origin species matches exactly.
pub fn imgt_get_by_species(
  records : Array[ImgtRecord],
  species : String,
) -> Array[ImgtRecord] {
  let out = Array::new()
  for i in 0..<records.length() {
    if records[i].header.species == species {
      out.push(records[i])
    }
  }
  out
}

///|
/// Return all records whose gene designation matches exactly.
pub fn imgt_get_by_gene(
  records : Array[ImgtRecord],
  gene : String,
) -> Array[ImgtRecord] {
  let out = Array::new()
  for i in 0..<records.length() {
    if records[i].header.gene == gene {
      out.push(records[i])
    }
  }
  out
}

///|
/// Return all records whose specific locus matches exactly.
pub fn imgt_filter_by_locus(
  records : Array[ImgtRecord],
  locus : String,
) -> Array[ImgtRecord] {
  let out = Array::new()
  for i in 0..<records.length() {
    if records[i].header.locus == locus {
      out.push(records[i])
    }
  }
  out
}

///|
/// Return the set of distinct species values across all records, in order of
/// first appearance.
pub fn imgt_unique_species(records : Array[ImgtRecord]) -> Array[String] {
  let out = Array::new()
  for i in 0..<records.length() {
    let s = records[i].header.species
    if !string_array_contains(out, s) {
      out.push(s)
    }
  }
  out
}

///|
/// Return the set of distinct gene designations across all records, in order
/// of first appearance.
pub fn imgt_unique_genes(records : Array[ImgtRecord]) -> Array[String] {
  let out = Array::new()
  for i in 0..<records.length() {
    let g = records[i].header.gene
    if !string_array_contains(out, g) {
      out.push(g)
    }
  }
  out
}

///|
/// Build a short human-readable summary of the records: counts and distinct
/// species/genes.
pub fn imgt_summary(records : Array[ImgtRecord]) -> String {
  let species = imgt_unique_species(records)
  let genes = imgt_unique_genes(records)
  let buf = StringBuilder::new()
  buf.write_string("IMGT records: ")
  buf.write_string(records.length().to_string())
  buf.write_string(", species: ")
  buf.write_string(species.length().to_string())
  buf.write_string(", genes: ")
  buf.write_string(genes.length().to_string())
  buf.write_string("\nspecies: [")
  for i in 0..<species.length() {
    if i > 0 {
      buf.write_string(", ")
    }
    buf.write_string(species[i])
  }
  buf.write_string("]\ngenes: [")
  for i in 0..<genes.length() {
    if i > 0 {
      buf.write_string(", ")
    }
    buf.write_string(genes[i])
  }
  buf.write_string("]")
  buf.to_string()
}

// ============================================================================
// SeqRecord integration
// ============================================================================

///|
/// Convert ImgtRecord objects to SeqRecord objects. The id and name become the
/// header's seq_id; the description becomes the full pipe-separated header.
pub fn imgt_to_seq_records(records : Array[ImgtRecord]) -> Array[SeqRecord] {
  let out = Array::new()
  for i in 0..<records.length() {
    let r = records[i]
    out.push(
      SeqRecord::new(
        Seq::new(r.sequence),
        id=r.header.seq_id,
        name=r.header.seq_id,
        description=r.header.to_string(),
      ),
    )
  }
  out
}

///|
/// Build ImgtRecord objects from SeqRecord objects using a minimal IMGT
/// header: accession and seq_id are set to the record id, the seq_length is
/// derived from the sequence length, and all other fields are empty.
pub fn imgt_from_seq_records(records : Array[SeqRecord]) -> Array[ImgtRecord] {
  let out = Array::new()
  for i in 0..<records.length() {
    let rec = records[i]
    let seq_str = rec.seq.to_string()
    let seq_len_str = seq_str.length().to_string() + " bp"
    let header = ImgtHeader::{
      accession: rec.id,
      seq_id: rec.id,
      official_name: "",
      region: "",
      date_update: "",
      date_creation: "",
      species: "",
      gene: "",
      imgt_group: "",
      locus: "",
      seq_length: seq_len_str,
      seq_type: "",
      comments: rec.description,
    }
    out.push(ImgtRecord::{ header, sequence: seq_str })
  }
  out
}

// ============================================================================
// Sample data
// ============================================================================

///|
/// Return 3 sample IMGT records: 2 HLA-A alleles (human) and 1 HLA-B allele
/// (human), suitable for tests and demos.
pub fn imgt_sample_records() -> Array[ImgtRecord] {
  let out = Array::new()
  // HLA-A*01:01:01:01 (human)
  out.push(
    ImgtRecord::{
      header: ImgtHeader::{
        accession: "HLA00001",
        seq_id: "A*01:01:01:01",
        official_name: "A*01:01:01:01",
        region: "ORF",
        date_update: "2015-03-31",
        date_creation: "1996-01-19",
        species: "human",
        gene: "HLA-A",
        imgt_group: "I",
        locus: "MHC",
        seq_length: "2703 bp",
        seq_type: " genomic-DNA",
        comments: "",
      },
      sequence: "ATGGCCGTCATGGCGCCCCGAACCCTCCTCCTC",
    },
  )
  // HLA-A*02:01:01:01 (human)
  out.push(
    ImgtRecord::{
      header: ImgtHeader::{
        accession: "HLA00002",
        seq_id: "A*02:01:01:01",
        official_name: "A*02:01:01:01",
        region: "ORF",
        date_update: "2016-05-12",
        date_creation: "1996-01-19",
        species: "human",
        gene: "HLA-A",
        imgt_group: "I",
        locus: "MHC",
        seq_length: "2709 bp",
        seq_type: " genomic-DNA",
        comments: "",
      },
      sequence: "ATGGCCGTCATGGCGCCCCGAACCCTCCTCCTG",
    },
  )
  // HLA-B*07:02:01 (human)
  out.push(
    ImgtRecord::{
      header: ImgtHeader::{
        accession: "HLA00003",
        seq_id: "B*07:02:01",
        official_name: "B*07:02:01",
        region: "ORF",
        date_update: "2017-07-18",
        date_creation: "1996-01-19",
        species: "human",
        gene: "HLA-B",
        imgt_group: "I",
        locus: "MHC",
        seq_length: "2711 bp",
        seq_type: " genomic-DNA",
        comments: "",
      },
      sequence: "ATGGCGCCCCGAACCCTCCTCCTGCTGCTCTC",
    },
  )
  out
}

// ============================================================================
// Private helpers
// ============================================================================

///|
/// Return true if `arr` contains the string `value`.
fn string_array_contains(arr : Array[String], value : String) -> Bool {
  let mut found = false
  for i in 0..<arr.length() {
    if arr[i] == value {
      found = true
      break
    }
  }
  found
}
