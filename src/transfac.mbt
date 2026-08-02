// Bio.motifs.transfac - TRANSFAC motif format parsing and writing.
// Inspired by Biopython's Bio.motifs.transfac module.
//
// Parses TRANSFAC database files containing transcription factor binding
// site matrices. TRANSFAC format uses two-letter field codes and includes
// position frequency matrices (PFM) for each motif.
//
// TRANSFAC format example:
//   VV  VERSION January 2024
//   XX
//   //
//   AC  T00001
//   XX
//   ID  AP-1
//   XX
//   DE  activator protein 1
//   XX
//   BF  Jun
//   XX
//   P0  A C G T
//   01  4 0 0 0
//   02  0 0 4 0
//   03  0 4 0 0
//   04  0 0 0 4
//   XX
//   //
//
// Reference:
//   Wingender E et al. (2000) TRANSFAC: an integrated system for gene
//   expression regulation. Nucleic Acids Research 28:316-319.

// ===========================================================================
// Data structures
// ===========================================================================

///|
/// A TRANSFAC reference entry.
pub struct TransfacReference {
  pub mut number : String // RN field (reference number)
  pub mut authors : String // RA field (authors)
  pub mut title : String // RT field (title)
  pub mut journal : String // RL field (journal/data)
  pub mut pubmed_id : String // RX field (PubMed ID)
}

///|
/// Create a new TransfacReference.
pub fn TransfacReference::new() -> TransfacReference {
  TransfacReference::{
    number: "",
    authors: "",
    title: "",
    journal: "",
    pubmed_id: "",
  }
}

///|
/// A TRANSFAC motif with position frequency matrix and metadata.
pub struct TransfacMotifFull {
  pub accession : String // AC field (accession number)
  pub mut identifier : String // ID field (identifier)
  pub mut description : String // DE field (short description)
  pub mut binding_factors : Array[String] // BF field (binding factors)
  pub mut comments : String // CC field (comments)
  pub mut copyright : String // CO field (copyright notice)
  pub mut name : String // NA field (name of binding factor)
  pub mut species : String // OS field (species/taxon)
  pub mut taxonomic_class : String // OC field (taxonomic classification)
  pub mut pfm : Array[Array[Int]] // position frequency matrix (width x 4 for A,C,G,T)
  pub mut alphabet : Array[String] // alphabet (default: A,C,G,T)
  pub mut width : Int // motif width
  pub mut references : Array[TransfacReference] // references list
  pub mut binding_sites : Array[String] // BS field (binding sites)
  pub mut external_refs : Array[String] // DR field (external database references)
  pub mut creation_date : String // DT field (date)
  pub mut motif_type : String // TY field (type)
}

///|
/// Create a new TransfacMotifFull.
pub fn TransfacMotifFull::new(accession~ : String = "") -> TransfacMotifFull {
  TransfacMotifFull::{
    accession,
    identifier: "",
    description: "",
    binding_factors: Array::new(),
    comments: "",
    copyright: "",
    name: "",
    species: "",
    taxonomic_class: "",
    pfm: Array::new(),
    alphabet: ["A", "C", "G", "T"],
    width: 0,
    references: Array::new(),
    binding_sites: Array::new(),
    external_refs: Array::new(),
    creation_date: "",
    motif_type: "",
  }
}

///|
/// Get the consensus sequence of the PFM.
pub fn TransfacMotifFull::consensus(self : TransfacMotifFull) -> String {
  let result = StringBuilder::new()
  for col = 0; col < self.width; col = col + 1 {
    let mut max_count = -1
    let mut max_base = "N"
    for row = 0; row < self.alphabet.length(); row = row + 1 {
      if col < self.pfm[row].length() {
        let count = self.pfm[row][col]
        if count > max_count {
          max_count = count
          max_base = self.alphabet[row]
        }
      }
    }
    result.write_string(max_base)
  }
  result.to_string()
}

///|
/// Get the total count at a position (sum of all bases).
pub fn TransfacMotifFull::total_at(self : TransfacMotifFull, position : Int) -> Int {
  if position < 0 || position >= self.width {
    return 0
  }
  let mut total = 0
  for row = 0; row < self.pfm.length(); row = row + 1 {
    if position < self.pfm[row].length() {
      total = total + self.pfm[row][position]
    }
  }
  total
}

///|
/// Get the frequency of a specific base at a position.
pub fn TransfacMotifFull::frequency(
  self : TransfacMotifFull,
  base : String,
  position : Int,
) -> Double {
  let row_idx = self.alphabet_index(base)
  if row_idx < 0 || position < 0 || position >= self.width {
    return 0.0
  }
  if row_idx >= self.pfm.length() || position >= self.pfm[row_idx].length() {
    return 0.0
  }
  let total = self.total_at(position)
  if total == 0 {
    return 0.0
  }
  self.pfm[row_idx][position].to_double() / total.to_double()
}

///|
/// Find the index of a base in the alphabet.
pub fn TransfacMotifFull::alphabet_index(self : TransfacMotifFull, base : String) -> Int {
  let mut i = 0
  while i < self.alphabet.length() {
    if self.alphabet[i] == base {
      return i
    }
    i = i + 1
  }
  -1
}

///|
/// A TRANSFAC record (file-level container for multiple motifs).
pub struct TransfacRecordFull {
  pub mut version : String // VV field (version)
  pub motifs : Array[TransfacMotifFull] // list of motifs
}

///|
/// Create a new empty TransfacRecordFull.
pub fn TransfacRecordFull::new() -> TransfacRecordFull {
  TransfacRecordFull::{ version: "", motifs: Array::new() }
}

///|
/// Get the number of motifs.
pub fn TransfacRecordFull::num_motifs(self : TransfacRecordFull) -> Int {
  self.motifs.length()
}

///|
/// Find a motif by accession number.
pub fn TransfacRecordFull::find(self : TransfacRecordFull, accession : String) -> TransfacMotifFull? {
  for m in self.motifs {
    if m.accession == accession {
      return Some(m)
    }
  }
  None
}

// ===========================================================================
// Parsing
// ===========================================================================

///|
/// Parse TRANSFAC format text into a TransfacRecordFull.
pub fn transfac_parse_full(text : String) -> TransfacRecordFull {
  let record = TransfacRecordFull::new()
  let lines = text.split("\n").to_array()
  // Use local mutable variables for building the current motif
  let mut has_motif = false
  let mut m_accession = ""
  let mut m_identifier = ""
  let mut m_description = ""
  let mut m_name = ""
  let mut m_species = ""
  let mut m_taxonomic_class = ""
  let mut m_comments = ""
  let mut m_copyright = ""
  let mut m_creation_date = ""
  let mut m_motif_type = ""
  let mut m_binding_factors : Array[String] = Array::new()
  let mut m_binding_sites : Array[String] = Array::new()
  let mut m_external_refs : Array[String] = Array::new()
  let mut m_alphabet : Array[String] = ["A", "C", "G", "T"]
  let mut m_pfm : Array[Array[Int]] = Array::new()
  let mut m_references : Array[TransfacReference] = Array::new()
  let mut in_matrix = false
  // Reference building variables
  let mut has_ref = false
  let mut ref_number = ""
  let mut ref_authors = ""
  let mut ref_title = ""
  let mut ref_journal = ""
  let mut ref_pubmed = ""
  for line in lines {
    let trimmed = line.trim().to_string()
    if trimmed.length() == 0 {
      continue
    }
    // Get the two-letter field code
    let field = if trimmed.length() >= 2 {
      trimmed[0:2].to_string()
    } else {
      ""
    }
    let value = if trimmed.length() > 2 {
      trimmed[2:].trim().to_string()
    } else {
      ""
    }
    // Version
    if field == "VV" {
      record.version = value
      continue
    }
    // Empty line (XX)
    if field == "XX" {
      continue
    }
    // Helper to finalize current reference
    if field != "RA" && field != "RT" && field != "RL" && field != "RX" && has_ref {
      let r = TransfacReference::new()
      r.number = ref_number
      r.authors = ref_authors
      r.title = ref_title
      r.journal = ref_journal
      r.pubmed_id = ref_pubmed
      m_references.push(r)
      has_ref = false
      ref_number = ""
      ref_authors = ""
      ref_title = ""
      ref_journal = ""
      ref_pubmed = ""
    }
    // End of motif (//)
    if field == "//" {
      if has_motif {
        let m = TransfacMotifFull::new(accession=m_accession)
        m.identifier = m_identifier
        m.description = m_description
        m.name = m_name
        m.species = m_species
        m.taxonomic_class = m_taxonomic_class
        m.comments = m_comments
        m.copyright = m_copyright
        m.creation_date = m_creation_date
        m.motif_type = m_motif_type
        m.binding_factors = m_binding_factors
        m.binding_sites = m_binding_sites
        m.external_refs = m_external_refs
        m.alphabet = m_alphabet
        m.pfm = m_pfm
        m.width = transfac_matrix_width(m_pfm)
        m.references = m_references
        record.motifs.push(m)
      }
      // Reset all motif variables
      has_motif = false
      m_accession = ""
      m_identifier = ""
      m_description = ""
      m_name = ""
      m_species = ""
      m_taxonomic_class = ""
      m_comments = ""
      m_copyright = ""
      m_creation_date = ""
      m_motif_type = ""
      m_binding_factors = Array::new()
      m_binding_sites = Array::new()
      m_external_refs = Array::new()
      m_alphabet = ["A", "C", "G", "T"]
      m_pfm = Array::new()
      m_references = Array::new()
      in_matrix = false
      continue
    }
    // Accession (AC) - start of a new motif
    if field == "AC" {
      if has_motif {
        let m = TransfacMotifFull::new(accession=m_accession)
        m.identifier = m_identifier
        m.description = m_description
        m.name = m_name
        m.species = m_species
        m.taxonomic_class = m_taxonomic_class
        m.comments = m_comments
        m.copyright = m_copyright
        m.creation_date = m_creation_date
        m.motif_type = m_motif_type
        m.binding_factors = m_binding_factors
        m.binding_sites = m_binding_sites
        m.external_refs = m_external_refs
        m.alphabet = m_alphabet
        m.pfm = m_pfm
        m.width = transfac_matrix_width(m_pfm)
        m.references = m_references
        record.motifs.push(m)
        m_binding_factors = Array::new()
        m_binding_sites = Array::new()
        m_external_refs = Array::new()
        m_references = Array::new()
        m_pfm = Array::new()
      }
      m_accession = value
      has_motif = true
      in_matrix = false
      continue
    }
    if !has_motif {
      continue
    }
    // Identifier (ID)
    if field == "ID" {
      m_identifier = value
      continue
    }
    // Description (DE)
    if field == "DE" {
      if m_description.length() > 0 {
        m_description = m_description + " " + value
      } else {
        m_description = value
      }
      continue
    }
    // Binding factors (BF)
    if field == "BF" {
      m_binding_factors.push(value)
      continue
    }
    // Comments (CC)
    if field == "CC" {
      if m_comments.length() > 0 {
        m_comments = m_comments + " " + value
      } else {
        m_comments = value
      }
      continue
    }
    // Copyright (CO)
    if field == "CO" {
      m_copyright = value
      continue
    }
    // Name (NA)
    if field == "NA" {
      m_name = value
      continue
    }
    // Species (OS)
    if field == "OS" {
      m_species = value
      continue
    }
    // Taxonomic classification (OC)
    if field == "OC" {
      m_taxonomic_class = value
      continue
    }
    // Date (DT)
    if field == "DT" {
      m_creation_date = value
      continue
    }
    // Type (TY)
    if field == "TY" {
      m_motif_type = value
      continue
    }
    // External references (DR)
    if field == "DR" {
      m_external_refs.push(value)
      continue
    }
    // Binding sites (BS)
    if field == "BS" {
      m_binding_sites.push(value)
      continue
    }
    // Reference fields
    if field == "RN" {
      has_ref = true
      ref_number = value
      continue
    }
    if field == "RA" && has_ref {
      ref_authors = value
      continue
    }
    if field == "RT" && has_ref {
      ref_title = value
      continue
    }
    if field == "RL" && has_ref {
      ref_journal = value
      continue
    }
    if field == "RX" && has_ref {
      ref_pubmed = value
      continue
    }
    // Position matrix header (P0)
    if field == "P0" {
      in_matrix = true
      // Parse alphabet from header (e.g., "A C G T")
      let header_parts = value.split(" ").to_array()
      let alphabet : Array[String] = Array::new()
      for p in header_parts {
        let s = p.trim().to_string()
        if s.length() > 0 {
          alphabet.push(s)
        }
      }
      if alphabet.length() > 0 {
        m_alphabet = alphabet
        // Initialize PFM rows (one row per alphabet letter)
        // Use loop push to create independent arrays (Array::make shares references)
        m_pfm = Array::new()
        let mut k = 0
        while k < alphabet.length() {
          m_pfm.push(Array::new())
          k = k + 1
        }
      }
      continue
    }
    // Matrix data rows (start with digits)
    if in_matrix && field.length() >= 2 && field[0:1].to_string() >= "0" && field[0:1].to_string() <= "9" {
      // Parse the row number and counts
      let parts = trimmed.split(" ").to_array()
      let counts : Array[Int] = Array::new()
      let mut first = true
      for p in parts {
        let s = p.trim().to_string()
        if first {
          // Skip the row number (01, 02, etc.)
          first = false
          continue
        }
        if s.length() > 0 {
          counts.push(parse_int(s))
        }
      }
      // Add counts to each alphabet row
      let mut i = 0
      while i < counts.length() && i < m_pfm.length() {
        m_pfm[i].push(counts[i])
        i = i + 1
      }
      continue
    }
  }
  // Save last motif if not terminated by //
  if has_motif {
    let m = TransfacMotifFull::new(accession=m_accession)
    m.identifier = m_identifier
    m.description = m_description
    m.name = m_name
    m.species = m_species
    m.taxonomic_class = m_taxonomic_class
    m.comments = m_comments
    m.copyright = m_copyright
    m.creation_date = m_creation_date
    m.motif_type = m_motif_type
    m.binding_factors = m_binding_factors
    m.binding_sites = m_binding_sites
    m.external_refs = m_external_refs
    m.alphabet = m_alphabet
    m.pfm = m_pfm
    m.width = transfac_matrix_width(m_pfm)
    m.references = m_references
    record.motifs.push(m)
  }
  record
}

///|
/// Calculate the width of a TRANSFAC PFM matrix.
fn transfac_matrix_width(pfm : Array[Array[Int]]) -> Int {
  if pfm.length() == 0 {
    return 0
  }
  pfm[0].length()
}

// ===========================================================================
// Serialization
// ===========================================================================

///|
/// Serialize a TransfacMotifFull to TRANSFAC format.
pub fn transfac_motif_full_to_string(motif : TransfacMotifFull) -> String {
  let sb = StringBuilder::new()
  sb.write_string("AC  " + motif.accession + "\nXX\n")
  if motif.identifier.length() > 0 {
    sb.write_string("ID  " + motif.identifier + "\nXX\n")
  }
  if motif.name.length() > 0 {
    sb.write_string("NA  " + motif.name + "\nXX\n")
  }
  if motif.description.length() > 0 {
    sb.write_string("DE  " + motif.description + "\nXX\n")
  }
  if motif.species.length() > 0 {
    sb.write_string("OS  " + motif.species + "\nXX\n")
  }
  if motif.taxonomic_class.length() > 0 {
    sb.write_string("OC  " + motif.taxonomic_class + "\nXX\n")
  }
  for bf in motif.binding_factors {
    sb.write_string("BF  " + bf + "\nXX\n")
  }
  if motif.pfm.length() > 0 {
    // Header
    sb.write_string("P0  ")
    let mut first = true
    for a in motif.alphabet {
      if !first {
        sb.write_char(' ')
      }
      sb.write_string(a)
      first = false
    }
    sb.write_char('\n')
    // Matrix rows
    let width = if motif.pfm.length() > 0 { motif.pfm[0].length() } else { 0 }
    for col = 0; col < width; col = col + 1 {
      let row_num = (col + 1).to_string()
      let padded = if row_num.length() < 2 { "0" + row_num } else { row_num }
      sb.write_string(padded + " ")
      for row = 0; row < motif.pfm.length(); row = row + 1 {
        if col < motif.pfm[row].length() {
          if row > 0 {
            sb.write_char(' ')
          }
          sb.write_string(motif.pfm[row][col].to_string())
        }
      }
      sb.write_char('\n')
    }
    sb.write_string("XX\n")
  }
  if motif.comments.length() > 0 {
    sb.write_string("CC  " + motif.comments + "\nXX\n")
  }
  sb.write_string("//\n")
  sb.to_string()
}

///|
/// Serialize a TransfacRecordFull to TRANSFAC format.
pub fn transfac_record_full_to_string(record : TransfacRecordFull) -> String {
  let sb = StringBuilder::new()
  if record.version.length() > 0 {
    sb.write_string("VV  " + record.version + "\nXX\n//\n")
  }
  for motif in record.motifs {
    sb.write_string(transfac_motif_full_to_string(motif))
  }
  sb.to_string()
}

// ===========================================================================
// Sample data
// ===========================================================================

///|
/// Create a sample TRANSFAC format text for testing.
pub fn transfac_sample_text_full() -> String {
  "VV  January 2024\nXX\n//\nAC  T00001\nXX\nID  AP-1\nXX\nDE  activator protein 1\nXX\nNA  AP-1\nXX\nOS  Homo sapiens\nXX\nBF  Jun\nXX\nBF  Fos\nXX\nP0  A C G T\n01  4 0 0 0\n02  0 0 4 0\n03  0 4 0 0\n04  0 0 0 4\n05  3 1 0 0\nXX\nCC  AP-1 transcription factor binding site\nXX\n//\nAC  T00002\nXX\nID  NF-kappaB\nXX\nDE  nuclear factor kappa B\nXX\nNA  NF-kappaB\nXX\nOS  Mus musculus\nXX\nP0  A C G T\n01  0 0 8 0\n02  8 0 0 0\n03  0 0 0 8\n04  0 8 0 0\nXX\n//\n"
}

///|
/// Create a sample TransfacMotifFull for testing.
pub fn transfac_sample_motif_full() -> TransfacMotifFull {
  let m = TransfacMotifFull::new(accession="T00001")
  m.identifier = "AP-1"
  m.description = "activator protein 1"
  m.name = "AP-1"
  m.species = "Homo sapiens"
  m.binding_factors.push("Jun")
  m.binding_factors.push("Fos")
  m.pfm = [
    [4, 0, 0, 0, 3],
    [0, 0, 4, 0, 1],
    [0, 4, 0, 0, 0],
    [0, 0, 0, 4, 0],
  ]
  m.width = 5
  m.comments = "AP-1 transcription factor binding site"
  m
}
