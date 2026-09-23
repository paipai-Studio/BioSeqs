///|
/// Parser for the ExPASy Cellosaurus flat-file format.
///
/// Mirrors the record model and parsing behavior of Biopython's
/// `Bio.ExPASy.cellosaurus` module.

///|
/// Error raised for malformed or ambiguous Cellosaurus input.
pub suberror CellosaurusError {
  CellosaurusError(String)
}

///|
/// A Cellosaurus `DR` cross-reference.
pub struct CellosaurusCrossReference {
  database : String
  accession : String
} derive(Eq)
pub extend CellosaurusCrossReference with Eq::{equal, not_equal}

///|
/// Create a Cellosaurus cross-reference.
pub fn CellosaurusCrossReference::new(
  database~ : String,
  accession~ : String,
) -> CellosaurusCrossReference {
  CellosaurusCrossReference::{ database, accession, }
}

///|
/// A typed Cellosaurus cell-line record.
///
/// Fields correspond to the official two-letter line codes:
/// ID, AC, AS, SY, DR, RX, WW, CC, ST, DI, OX, HI, OI, SX, AG, CA, and DT.
pub struct CellosaurusRecord {
  identifier : String
  accession : String
  secondary_accessions : String
  synonyms : String
  cross_references : Array[CellosaurusCrossReference]
  reference_identifiers : Array[String]
  web_pages : Array[String]
  comments : Array[String]
  str_profiles : Array[String]
  diseases : Array[String]
  species : Array[String]
  hierarchy : Array[String]
  same_individual : Array[String]
  sex : String
  age : String
  category : String
  date : String
} derive(Eq)
pub extend CellosaurusRecord with Eq::{equal, not_equal}

///|
/// Create a Cellosaurus record.
pub fn CellosaurusRecord::new(
  identifier~ : String,
  accession? : String = "",
  secondary_accessions? : String = "",
  synonyms? : String = "",
  cross_references? : Array[CellosaurusCrossReference] = [],
  reference_identifiers? : Array[String] = [],
  web_pages? : Array[String] = [],
  comments? : Array[String] = [],
  str_profiles? : Array[String] = [],
  diseases? : Array[String] = [],
  species? : Array[String] = [],
  hierarchy? : Array[String] = [],
  same_individual? : Array[String] = [],
  sex? : String = "",
  age? : String = "",
  category? : String = "",
  date? : String = "",
) -> CellosaurusRecord {
  CellosaurusRecord::{
    identifier,
    accession,
    secondary_accessions,
    synonyms,
    cross_references,
    reference_identifiers,
    web_pages,
    comments,
    str_profiles,
    diseases,
    species,
    hierarchy,
    same_individual,
    sex,
    age,
    category,
    date,
  }
}

///|
/// Return secondary accessions split on semicolons.
pub fn CellosaurusRecord::secondary_accession_list(
  self : CellosaurusRecord,
) -> Array[String] {
  cellosaurus_split_semicolon_list(self.secondary_accessions)
}

///|
/// Return cell-line synonyms split on semicolons.
pub fn CellosaurusRecord::synonym_list(
  self : CellosaurusRecord,
) -> Array[String] {
  cellosaurus_split_semicolon_list(self.synonyms)
}

///|
/// Return accessions referenced by the requested external database.
pub fn CellosaurusRecord::cross_references_for(
  self : CellosaurusRecord,
  database : String,
) -> Array[String] {
  let result : Array[String] = []
  for reference in self.cross_references {
    if reference.database == database {
      result.push(reference.accession)
    }
  }
  result
}

///|
/// Test whether any OX line contains the requested species text.
pub fn CellosaurusRecord::has_species(
  self : CellosaurusRecord,
  species_name : String,
) -> Bool {
  for origin in self.species {
    if origin.contains(species_name) {
      return true
    }
  }
  false
}

///|
/// Return a compact representation similar to Biopython's `Record.__repr__`.
pub fn CellosaurusRecord::repr(self : CellosaurusRecord) -> String {
  if self.identifier.length() == 0 {
    "CellosaurusRecord ( )"
  } else if self.accession.length() == 0 {
    "CellosaurusRecord (" + self.identifier + ")"
  } else {
    "CellosaurusRecord (" + self.identifier + ", " + self.accession + ")"
  }
}

///|
/// Serialize a record to canonical Cellosaurus flat-file text.
pub fn CellosaurusRecord::to_string(self : CellosaurusRecord) -> String {
  let output = StringBuilder::new()
  cellosaurus_write_field(output, "ID", self.identifier)
  cellosaurus_write_field(output, "AC", self.accession)
  cellosaurus_write_field(output, "AS", self.secondary_accessions)
  cellosaurus_write_field(output, "SY", self.synonyms)
  for reference in self.cross_references {
    cellosaurus_write_field(
      output,
      "DR",
      reference.database + "; " + reference.accession,
    )
  }
  cellosaurus_write_fields(output, "RX", self.reference_identifiers)
  cellosaurus_write_fields(output, "WW", self.web_pages)
  cellosaurus_write_fields(output, "CC", self.comments)
  cellosaurus_write_fields(output, "ST", self.str_profiles)
  cellosaurus_write_fields(output, "DI", self.diseases)
  cellosaurus_write_fields(output, "OX", self.species)
  cellosaurus_write_fields(output, "HI", self.hierarchy)
  cellosaurus_write_fields(output, "OI", self.same_individual)
  cellosaurus_write_field(output, "SX", self.sex)
  cellosaurus_write_field(output, "AG", self.age)
  cellosaurus_write_field(output, "CA", self.category)
  cellosaurus_write_field(output, "DT", self.date)
  output.write_string("//\n")
  output.to_string()
}

///|
/// Parse all Cellosaurus records from flat-file text.
///
/// Header lines before the first ID record are ignored. A started record must
/// end with `//`; malformed DR lines and unterminated records raise
/// `CellosaurusError`.
pub fn cellosaurus_parse(
  text : String,
) -> Array[CellosaurusRecord] raise CellosaurusError {
  let records : Array[CellosaurusRecord] = []
  let mut current : CellosaurusBuilder? = None
  let lines = text.split("\n").to_array()
  for line_index, line_view in lines {
    let mut line = line_view.to_owned()
    if line.length() > 0 &&
      line.unsafe_get(line.length() - 1) == '\r'.to_int().to_uint16() {
      line = line[0:line.length() - 1].to_owned()
    }
    if line.length() < 2 {
      continue
    }
    let code = line[0:2].to_owned()
    let value = if line.length() > 5 { line[5:].trim().to_owned() } else { "" }
    if code == "ID" {
      match current {
        Some(_) =>
          raise CellosaurusError(
            "Cellosaurus record before line " +
            (line_index + 1).to_string() +
            " is missing its // terminator",
          )
        None => ()
      }
      let builder = CellosaurusBuilder::new()
      builder.identifier = value
      current = Some(builder)
    } else if code == "//" {
      match current {
        Some(builder) => {
          records.push(builder.finish())
          current = None
        }
        None => ()
      }
    } else {
      match current {
        Some(builder) =>
          cellosaurus_apply_field(
            builder,
            code,
            value,
            line_number=line_index + 1,
          )
        None => ()
      }
    }
  }
  match current {
    Some(_) =>
      raise CellosaurusError(
        "Unexpected end of Cellosaurus input: record is missing its // terminator",
      )
    None => ()
  }
  records
}

///|
/// Read zero or one Cellosaurus record.
///
/// Returns `None` for empty input and raises `CellosaurusError` if more than
/// one record is present.
pub fn cellosaurus_read(
  text : String,
) -> CellosaurusRecord? raise CellosaurusError {
  let records = cellosaurus_parse(text)
  if records.length() == 0 {
    None
  } else if records.length() == 1 {
    Some(records[0])
  } else {
    raise CellosaurusError("More than one cell line record found")
  }
}

///|
/// Return a small valid Cellosaurus record for examples.
pub fn cellosaurus_sample_text() -> String {
  "ID   #15310-LN\n" +
  "AC   CVCL_E548\n" +
  "SY   15310-LN; TER461; TER-461; Extract 519\n" +
  "DR   dbMHC; 48439\n" +
  "DR   ECACC; 94050311\n" +
  "DR   Wikidata; Q54398957\n" +
  "WW   http://pathology.ucla.edu/workfiles/360cx.pdf\n" +
  "CC   Population: Caucasian; French Canadian.\n" +
  "CC   Cell type: B-cell; CL=CL_0000236.\n" +
  "OX   NCBI_TaxID=9606; ! Homo sapiens (Human)\n" +
  "SX   Female\n" +
  "AG   Age unspecified\n" +
  "CA   Transformed cell line\n" +
  "DT   Created: 22-10-12; Last updated: 30-01-24; Version: 18\n" +
  "//\n"
}

///|
struct CellosaurusBuilder {
  mut identifier : String
  mut accession : String
  mut secondary_accessions : String
  mut synonyms : String
  cross_references : Array[CellosaurusCrossReference]
  reference_identifiers : Array[String]
  web_pages : Array[String]
  comments : Array[String]
  str_profiles : Array[String]
  diseases : Array[String]
  species : Array[String]
  hierarchy : Array[String]
  same_individual : Array[String]
  mut sex : String
  mut age : String
  mut category : String
  mut date : String
}

///|
fn CellosaurusBuilder::new() -> CellosaurusBuilder {
  CellosaurusBuilder::{
    identifier: "",
    accession: "",
    secondary_accessions: "",
    synonyms: "",
    cross_references: [],
    reference_identifiers: [],
    web_pages: [],
    comments: [],
    str_profiles: [],
    diseases: [],
    species: [],
    hierarchy: [],
    same_individual: [],
    sex: "",
    age: "",
    category: "",
    date: "",
  }
}

///|
fn CellosaurusBuilder::finish(self : CellosaurusBuilder) -> CellosaurusRecord {
  CellosaurusRecord::new(
    identifier=self.identifier,
    accession=self.accession,
    secondary_accessions=self.secondary_accessions,
    synonyms=self.synonyms,
    cross_references=self.cross_references,
    reference_identifiers=self.reference_identifiers,
    web_pages=self.web_pages,
    comments=self.comments,
    str_profiles=self.str_profiles,
    diseases=self.diseases,
    species=self.species,
    hierarchy=self.hierarchy,
    same_individual=self.same_individual,
    sex=self.sex,
    age=self.age,
    category=self.category,
    date=self.date,
  )
}

///|
fn cellosaurus_apply_field(
  builder : CellosaurusBuilder,
  code : String,
  value : String,
  line_number~ : Int,
) -> Unit raise CellosaurusError {
  match code {
    "AC" => builder.accession = builder.accession + value
    "AS" => builder.secondary_accessions = builder.secondary_accessions + value
    "SY" => builder.synonyms = builder.synonyms + value
    "DR" =>
      match value.find(";") {
        Some(separator) => {
          let database = value[0:separator].trim().to_owned()
          let accession = value[separator + 1:].trim().to_owned()
          if database.length() == 0 || accession.length() == 0 {
            raise CellosaurusError(
              "Invalid DR field at line " + line_number.to_string(),
            )
          }
          builder.cross_references.push(
            CellosaurusCrossReference::new(database~, accession~),
          )
        }
        None =>
          raise CellosaurusError(
            "Invalid DR field at line " +
            line_number.to_string() +
            ": expected 'database; accession'",
          )
      }
    "RX" => builder.reference_identifiers.push(value)
    "WW" => builder.web_pages.push(value)
    "CC" => builder.comments.push(value)
    "ST" => builder.str_profiles.push(value)
    "DI" => builder.diseases.push(value)
    "OX" => builder.species.push(value)
    "HI" => builder.hierarchy.push(value)
    "OI" => builder.same_individual.push(value)
    "SX" => builder.sex = builder.sex + value
    "AG" => builder.age = builder.age + value
    "CA" => builder.category = builder.category + value
    "DT" => builder.date = builder.date + value
    _ => ()
  }
}

///|
fn cellosaurus_split_semicolon_list(text : String) -> Array[String] {
  let result : Array[String] = []
  for part in text.split(";") {
    let value = part.trim().to_owned()
    if value.length() > 0 {
      result.push(value)
    }
  }
  result
}

///|
fn cellosaurus_write_field(
  output : StringBuilder,
  code : String,
  value : String,
) -> Unit {
  if value.length() > 0 {
    output.write_string(code)
    output.write_string("   ")
    output.write_string(value)
    output.write_string("\n")
  }
}

///|
fn cellosaurus_write_fields(
  output : StringBuilder,
  code : String,
  values : Array[String],
) -> Unit {
  for value in values {
    cellosaurus_write_field(output, code, value)
  }
}
