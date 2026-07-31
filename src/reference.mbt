///|
/// Bio.Reference module - Literature reference objects for SeqRecord.
/// Inspired by Biopython Bio.Reference.
///
/// Provides data structures for representing publications
/// associated with sequence records (GenBank, EMBL, etc.).
///
/// References:
///   - https://biopython.org/wiki/Reference

// ============================================================================
// Trait & Data Structures
// ============================================================================

///|
/// Abstract reference trait.
pub trait AbstractReference {
  /// Get the title of the reference.
  fn title(self : Self) -> String
  /// Get the authors string.
  fn authors(self : Self) -> String
  /// Get the journal string.
  fn journal(self : Self) -> String
  /// Get the year as a string.
  fn year(self : Self) -> String
  /// Get the reference as a formatted citation.
  fn citation(self : Self) -> String
  /// Get the PubMed ID.
  fn pubmed_id(self : Self) -> String
  /// Get the DOI.
  fn doi(self : Self) -> String
  /// Get the reference type (journal, book, patent, etc.).
  fn reference_type(self : Self) -> String
}

///|
/// Concrete literature reference implementation.
/// Corresponds to Biopython's Reference class.
pub struct BioReference {
  mut title_val : String
  mut authors_val : String
  mut journal_val : String
  mut year_val : String
  mut pubmed_id_val : String
  mut doi_val : String
  mut type_val : String
  locations : Array[(Int, Int)]
  mut comment : String
  mut data_pgr : String
}

///|
/// Create a new BioReference.
pub fn BioReference::new(
  title? : String = "",
  authors? : String = "",
  journal? : String = "",
  year? : String = "",
) -> BioReference {
  BioReference::{
    title_val: title,
    authors_val: authors,
    journal_val: journal,
    year_val: year,
    pubmed_id_val: "",
    doi_val: "",
    type_val: "journal article",
    locations: Array::new(),
    comment: "",
    data_pgr: "",
  }
}

///|
/// Create a reference with PubMed ID.
pub fn BioReference::with_pubmed(
  title : String,
  authors : String,
  journal : String,
  year : String,
  pubmed_id : String,
) -> BioReference {
  let r = BioReference::new(title=title, authors=authors, journal=journal, year=year)
  r.pubmed_id_val = pubmed_id
  r
}

///|
/// Create a reference with DOI.
pub fn BioReference::with_doi(
  title : String,
  authors : String,
  journal : String,
  year : String,
  doi : String,
) -> BioReference {
  let r = BioReference::new(title=title, authors=authors, journal=journal, year=year)
  r.doi_val = doi
  r
}

pub impl AbstractReference for BioReference with fn title(self : BioReference) -> String { self.title_val }
pub impl AbstractReference for BioReference with fn authors(self : BioReference) -> String { self.authors_val }
pub impl AbstractReference for BioReference with fn journal(self : BioReference) -> String { self.journal_val }
pub impl AbstractReference for BioReference with fn year(self : BioReference) -> String { self.year_val }
pub impl AbstractReference for BioReference with fn pubmed_id(self : BioReference) -> String { self.pubmed_id_val }
pub impl AbstractReference for BioReference with fn doi(self : BioReference) -> String { self.doi_val }
pub impl AbstractReference for BioReference with fn reference_type(self : BioReference) -> String { self.type_val }
pub impl AbstractReference for BioReference with fn citation(self : BioReference) -> String {
  let mut s = self.authors_val
  if s.length() > 0 { s = s + ". " }
  s = s + self.title_val
  if s.length() > 0 { s = s + ". " }
  s = s + self.journal_val
  if self.year_val.length() > 0 {
    s = s + " (" + self.year_val + ")"
  }
  if self.pubmed_id_val.length() > 0 {
    s = s + " PMID: " + self.pubmed_id_val
  }
  if self.doi_val.length() > 0 {
    s = s + " DOI: " + self.doi_val
  }
  s
}

// ============================================================================
// BioReference Methods
// ============================================================================

///|
/// Set the title.
pub fn BioReference::set_title(self : BioReference, title : String) -> Unit {
  self.title_val = title
}

///|
/// Set the authors.
pub fn BioReference::set_authors(self : BioReference, authors : String) -> Unit {
  self.authors_val = authors
}

///|
/// Set the journal.
pub fn BioReference::set_journal(self : BioReference, journal : String) -> Unit {
  self.journal_val = journal
}

///|
/// Set the year.
pub fn BioReference::set_year(self : BioReference, year : String) -> Unit {
  self.year_val = year
}

///|
/// Set the PubMed ID.
pub fn BioReference::set_pubmed_id(self : BioReference, id : String) -> Unit {
  self.pubmed_id_val = id
}

///|
/// Set the DOI.
pub fn BioReference::set_doi(self : BioReference, doi : String) -> Unit {
  self.doi_val = doi
}

///|
/// Set the reference type.
pub fn BioReference::set_type(self : BioReference, type_ : String) -> Unit {
  self.type_val = type_
}

///|
/// Add a location (pair of start/end positions).
pub fn BioReference::add_location(self : BioReference, start : Int, end_ : Int) -> Unit {
  self.locations.push((start, end_))
}

///|
/// Get number of locations.
pub fn BioReference::get_n_locations(self : BioReference) -> Int {
  self.locations.length()
}

///|
/// Get location by index.
pub fn BioReference::get_location(self : BioReference, idx : Int) -> (Int, Int) {
  if idx >= 0 && idx < self.locations.length() {
    self.locations[idx]
  } else {
    (0, 0)
  }
}

///|
/// Set comment.
pub fn BioReference::set_comment(self : BioReference, comment : String) -> Unit {
  self.comment = comment
}

///|
/// Get comment.
pub fn BioReference::get_comment(self : BioReference) -> String {
  self.comment
}

///|
/// Set PGR (Protein Group Reference) string.
pub fn BioReference::set_pgr(self : BioReference, pgr : String) -> Unit {
  self.data_pgr = pgr
}

///|
/// Get PGR string.
pub fn BioReference::get_pgr(self : BioReference) -> String {
  self.data_pgr
}

// ============================================================================
// ReferenceList
// ============================================================================

///|
/// Container for a list of references.
pub struct ReferenceList {
  refs : Array[BioReference]
}

///|
/// Create a new reference list.
pub fn ReferenceList::new() -> ReferenceList {
  ReferenceList::{ refs: Array::new() }
}

///|
/// Add a reference.
pub fn ReferenceList::add(self : ReferenceList, ref : BioReference) -> Unit {
  self.refs.push(ref)
}

///|
/// Get number of references.
pub fn ReferenceList::count(self : ReferenceList) -> Int {
  self.refs.length()
}

///|
/// Get reference by index.
pub fn ReferenceList::get(self : ReferenceList, idx : Int) -> BioReference {
  if idx >= 0 && idx < self.refs.length() {
    self.refs[idx]
  } else {
    BioReference::new()
  }
}

///|
/// Get all references.
pub fn ReferenceList::get_all(self : ReferenceList) -> Array[BioReference] {
  self.refs.copy()
}

///|
/// Get references by author name.
pub fn ReferenceList::by_author(self : ReferenceList, author : String) -> Array[BioReference] {
  let result : Array[BioReference] = Array::new()
  let mut i = 0
  while i < self.refs.length() {
    if self.refs[i].authors_val.contains(author) {
      result.push(self.refs[i])
    }
    i = i + 1
  }
  result
}

///|
/// Get references by year.
pub fn ReferenceList::by_year(self : ReferenceList, year : String) -> Array[BioReference] {
  let result : Array[BioReference] = Array::new()
  let mut i = 0
  while i < self.refs.length() {
    if self.refs[i].year_val == year {
      result.push(self.refs[i])
    }
    i = i + 1
  }
  result
}

///|
/// Get references with PubMed IDs.
pub fn ReferenceList::with_pubmed(self : ReferenceList) -> Array[BioReference] {
  let result : Array[BioReference] = Array::new()
  let mut i = 0
  while i < self.refs.length() {
    if self.refs[i].pubmed_id_val.length() > 0 {
      result.push(self.refs[i])
    }
    i = i + 1
  }
  result
}

///|
/// Format all references as a bibliography.
pub fn ReferenceList::bibliography(self : ReferenceList) -> String {
  let mut s = ""
  let mut i = 0
  while i < self.refs.length() {
    s = s + (i + 1).to_string() + ". " + self.refs[i].citation() + "\n"
    i = i + 1
  }
  s
}

// ============================================================================
// Sample Data
// ============================================================================

///|
/// Create sample references for testing.
pub fn reference_sample_data() -> ReferenceList {
  let list = ReferenceList::new()
  let ref1 = BioReference::with_pubmed(
    "The Human Genome: A Complete Sequence",
    "Nurk S, Koren S, Rhie A",
    "Science",
    "2022",
    "36189102",
  )
  ref1.add_location(1, 250000000)
  ref1.set_comment("Telomere-to-telomere assembly")
  list.add(ref1)

  let ref2 = BioReference::with_doi(
    "CRISPR-Cas9 genome editing",
    "Doudna JA, Charpentier E",
    "Science",
    "2014",
    "10.1126/science.1258096",
  )
  ref2.add_location(100, 5000)
  list.add(ref2)

  let ref3 = BioReference::new(
    title="AlphaFold: Protein Structure Prediction",
    authors="Jumper J, Evans R, Pritzel A",
    journal="Nature",
    year="2021",
  )
  ref3.set_type("journal article")
  list.add(ref3)

  list
}
