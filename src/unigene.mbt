///|
/// Parser and writer for NCBI UniGene flat files.
///
/// The data model follows Biopython's `Bio.UniGene` module while retaining
/// unknown key-value fields found inside SEQUENCE, PROTSIM, STS, and TXMAP
/// lines.

///|
/// Error raised for malformed or ambiguous UniGene input.
pub suberror UniGeneError {
  UniGeneError(String)
}

///|
/// A key-value field not covered by the standard UniGene schema.
pub struct UniGeneExtraField {
  key : String
  value : String
} derive(Eq)
pub extend UniGeneExtraField with Eq::{equal, not_equal}

///|
pub fn UniGeneExtraField::new(
  key~ : String,
  value~ : String,
) -> UniGeneExtraField {
  UniGeneExtraField::{ key, value, }
}

///|
/// One `SEQUENCE` line from a UniGene record.
pub struct UniGeneSequence {
  accession : String
  nucleotide_id : String
  protein_id : String
  clone_id : String
  read_end : String
  library_id : String
  mgc : String
  sequence_type : String
  trace : String
  is_image : Bool
  image_id : String
  extra_fields : Array[UniGeneExtraField]
} derive(Eq)
pub extend UniGeneSequence with Eq::{equal, not_equal}

///|
pub fn UniGeneSequence::new(
  accession? : String = "",
  nucleotide_id? : String = "",
  protein_id? : String = "",
  clone_id? : String = "",
  read_end? : String = "",
  library_id? : String = "",
  mgc? : String = "",
  sequence_type? : String = "",
  trace? : String = "",
  extra_fields? : Array[UniGeneExtraField] = [],
) -> UniGeneSequence {
  let image_id = unigene_image_id(clone_id)
  UniGeneSequence::{
    accession,
    nucleotide_id,
    protein_id,
    clone_id,
    read_end,
    library_id,
    mgc,
    sequence_type,
    trace,
    is_image: image_id.length() > 0,
    image_id,
    extra_fields,
  }
}

///|
/// Return a field by its uppercase UniGene key.
pub fn UniGeneSequence::field(self : UniGeneSequence, key : String) -> String? {
  match key {
    "ACC" => unigene_nonempty(self.accession)
    "NID" => unigene_nonempty(self.nucleotide_id)
    "PID" => unigene_nonempty(self.protein_id)
    "CLONE" => unigene_nonempty(self.clone_id)
    "END" => unigene_nonempty(self.read_end)
    "LID" => unigene_nonempty(self.library_id)
    "MGC" => unigene_nonempty(self.mgc)
    "SEQTYPE" => unigene_nonempty(self.sequence_type)
    "TRACE" => unigene_nonempty(self.trace)
    _ => unigene_extra_value(self.extra_fields, key)
  }
}

///|
/// Serialize this sequence entry without the fixed-width `SEQUENCE` tag.
pub fn UniGeneSequence::to_value(self : UniGeneSequence) -> String {
  let fields : Array[String] = []
  unigene_push_field(fields, "ACC", self.accession)
  unigene_push_field(fields, "NID", self.nucleotide_id)
  unigene_push_field(fields, "PID", self.protein_id)
  unigene_push_field(fields, "CLONE", self.clone_id)
  unigene_push_field(fields, "END", self.read_end)
  unigene_push_field(fields, "LID", self.library_id)
  unigene_push_field(fields, "MGC", self.mgc)
  unigene_push_field(fields, "SEQTYPE", self.sequence_type)
  unigene_push_field(fields, "TRACE", self.trace)
  for field in self.extra_fields {
    unigene_push_field(fields, field.key, field.value)
  }
  fields.join("; ")
}

///|
/// One `PROTSIM` protein-similarity line.
pub struct UniGeneProteinSimilarity {
  organism : String
  protein_gi : String
  protein_id : String
  percent : String
  alignment_length : String
  extra_fields : Array[UniGeneExtraField]
} derive(Eq)
pub extend UniGeneProteinSimilarity with Eq::{equal, not_equal}

///|
pub fn UniGeneProteinSimilarity::new(
  organism? : String = "",
  protein_gi? : String = "",
  protein_id? : String = "",
  percent? : String = "",
  alignment_length? : String = "",
  extra_fields? : Array[UniGeneExtraField] = [],
) -> UniGeneProteinSimilarity {
  UniGeneProteinSimilarity::{
    organism,
    protein_gi,
    protein_id,
    percent,
    alignment_length,
    extra_fields,
  }
}

///|
/// Return a field by its uppercase UniGene key.
pub fn UniGeneProteinSimilarity::field(
  self : UniGeneProteinSimilarity,
  key : String,
) -> String? {
  match key {
    "ORG" => unigene_nonempty(self.organism)
    "PROTGI" => unigene_nonempty(self.protein_gi)
    "PROTID" => unigene_nonempty(self.protein_id)
    "PCT" => unigene_nonempty(self.percent)
    "ALN" => unigene_nonempty(self.alignment_length)
    _ => unigene_extra_value(self.extra_fields, key)
  }
}

///|
pub fn UniGeneProteinSimilarity::to_value(
  self : UniGeneProteinSimilarity,
) -> String {
  let fields : Array[String] = []
  unigene_push_field(fields, "ORG", self.organism)
  unigene_push_field(fields, "PROTGI", self.protein_gi)
  unigene_push_field(fields, "PROTID", self.protein_id)
  unigene_push_field(fields, "PCT", self.percent)
  unigene_push_field(fields, "ALN", self.alignment_length)
  for field in self.extra_fields {
    unigene_push_field(fields, field.key, field.value)
  }
  fields.join("; ")
}

///|
/// One `STS` sequence-tagged-site line.
pub struct UniGeneSTS {
  accession : String
  unists : String
  extra_fields : Array[UniGeneExtraField]
} derive(Eq)
pub extend UniGeneSTS with Eq::{equal, not_equal}

///|
pub fn UniGeneSTS::new(
  accession? : String = "",
  unists? : String = "",
  extra_fields? : Array[UniGeneExtraField] = [],
) -> UniGeneSTS {
  UniGeneSTS::{ accession, unists, extra_fields, }
}

///|
pub fn UniGeneSTS::field(self : UniGeneSTS, key : String) -> String? {
  match key {
    "ACC" => unigene_nonempty(self.accession)
    "UNISTS" => unigene_nonempty(self.unists)
    _ => unigene_extra_value(self.extra_fields, key)
  }
}

///|
pub fn UniGeneSTS::to_value(self : UniGeneSTS) -> String {
  let fields : Array[String] = []
  unigene_push_field(fields, "ACC", self.accession)
  unigene_push_field(fields, "UNISTS", self.unists)
  for field in self.extra_fields {
    unigene_push_field(fields, field.key, field.value)
  }
  fields.join(" ")
}

///|
/// One `TXMAP` transcript-map line.
pub struct UniGeneTranscriptMap {
  marker : String
  radiation_hybrid_panel : String
  extra_fields : Array[UniGeneExtraField]
} derive(Eq)
pub extend UniGeneTranscriptMap with Eq::{equal, not_equal}

///|
pub fn UniGeneTranscriptMap::new(
  marker? : String = "",
  radiation_hybrid_panel? : String = "",
  extra_fields? : Array[UniGeneExtraField] = [],
) -> UniGeneTranscriptMap {
  UniGeneTranscriptMap::{ marker, radiation_hybrid_panel, extra_fields, }
}

///|
pub fn UniGeneTranscriptMap::field(
  self : UniGeneTranscriptMap,
  key : String,
) -> String? {
  match key {
    "MARKER" => unigene_nonempty(self.marker)
    "RHPANEL" => unigene_nonempty(self.radiation_hybrid_panel)
    _ => unigene_extra_value(self.extra_fields, key)
  }
}

///|
pub fn UniGeneTranscriptMap::to_value(self : UniGeneTranscriptMap) -> String {
  let fields : Array[String] = []
  unigene_push_field(fields, "MARKER", self.marker)
  unigene_push_field(fields, "RHPANEL", self.radiation_hybrid_panel)
  for field in self.extra_fields {
    unigene_push_field(fields, field.key, field.value)
  }
  fields.join("; ")
}

///|
/// A complete UniGene cluster record.
pub struct UniGeneRecord {
  identifier : String
  species : String
  title : String
  symbol : String
  cytoband : String
  expression : Array[String]
  restricted_expression : Array[String]
  genomic_terminus : String
  gene_id : String
  locuslink : String
  homology : Bool?
  chromosome : String
  protein_similarities : Array[UniGeneProteinSimilarity]
  sequences : Array[UniGeneSequence]
  sts : Array[UniGeneSTS]
  transcript_maps : Array[UniGeneTranscriptMap]
  declared_sequence_count : Int
} derive(Eq)
pub extend UniGeneRecord with Eq::{equal, not_equal}

///|
pub fn UniGeneRecord::new(
  identifier~ : String,
  species? : String = "",
  title? : String = "",
  symbol? : String = "",
  cytoband? : String = "",
  expression? : Array[String] = [],
  restricted_expression? : Array[String] = [],
  genomic_terminus? : String = "",
  gene_id? : String = "",
  locuslink? : String = "",
  homology? : Bool? = None,
  chromosome? : String = "",
  protein_similarities? : Array[UniGeneProteinSimilarity] = [],
  sequences? : Array[UniGeneSequence] = [],
  sts? : Array[UniGeneSTS] = [],
  transcript_maps? : Array[UniGeneTranscriptMap] = [],
) -> UniGeneRecord {
  UniGeneRecord::{
    identifier,
    species: if species.length() > 0 {
      species
    } else {
      unigene_species(identifier)
    },
    title,
    symbol,
    cytoband,
    expression,
    restricted_expression,
    genomic_terminus,
    gene_id,
    locuslink,
    homology,
    chromosome,
    protein_similarities,
    sequences,
    sts,
    transcript_maps,
    declared_sequence_count: sequences.length(),
  }
}

///|
/// Return all sequences with the requested `SEQTYPE`.
pub fn UniGeneRecord::sequences_of_type(
  self : UniGeneRecord,
  sequence_type : String,
) -> Array[UniGeneSequence] {
  let result : Array[UniGeneSequence] = []
  for sequence in self.sequences {
    if sequence.sequence_type == sequence_type {
      result.push(sequence)
    }
  }
  result
}

///|
/// Return protein similarities for an NCBI taxonomy identifier.
pub fn UniGeneRecord::protein_similarities_for(
  self : UniGeneRecord,
  organism : String,
) -> Array[UniGeneProteinSimilarity] {
  let result : Array[UniGeneProteinSimilarity] = []
  for similarity in self.protein_similarities {
    if similarity.organism == organism {
      result.push(similarity)
    }
  }
  result
}

///|
/// Return sequences whose CLONE value identifies an IMAGE clone.
pub fn UniGeneRecord::image_sequences(
  self : UniGeneRecord,
) -> Array[UniGeneSequence] {
  let result : Array[UniGeneSequence] = []
  for sequence in self.sequences {
    if sequence.is_image {
      result.push(sequence)
    }
  }
  result
}

///|
/// Find a sequence by GenBank/EMBL/DDBJ accession.
pub fn UniGeneRecord::find_sequence(
  self : UniGeneRecord,
  accession : String,
) -> UniGeneSequence? {
  for sequence in self.sequences {
    if sequence.accession == accession {
      return Some(sequence)
    }
  }
  None
}

///|
/// Return a compact record description.
pub fn UniGeneRecord::summary(self : UniGeneRecord) -> String {
  "UniGeneRecord(" +
  self.identifier +
  ", symbol=" +
  self.symbol +
  ", sequences=" +
  self.sequences.length().to_string() +
  ", protein_similarities=" +
  self.protein_similarities.length().to_string() +
  ")"
}

///|
/// Serialize this record to canonical fixed-width UniGene text.
pub fn UniGeneRecord::to_string(self : UniGeneRecord) -> String {
  let output = StringBuilder::new()
  unigene_write_line(output, "ID", self.identifier)
  unigene_write_line(output, "TITLE", self.title)
  unigene_write_line(output, "GENE", self.symbol)
  unigene_write_line(output, "CYTOBAND", self.cytoband)
  unigene_write_line(output, "GENE_ID", self.gene_id)
  unigene_write_line(output, "LOCUSLINK", self.locuslink)
  match self.homology {
    Some(value) =>
      unigene_write_line(output, "HOMOL", if value { "YES" } else { "NO" })
    None => ()
  }
  if self.expression.length() > 0 {
    unigene_write_line(output, "EXPRESS", self.expression.join("| "))
  }
  if self.restricted_expression.length() > 0 {
    unigene_write_line(
      output,
      "RESTR_EXPR",
      self.restricted_expression.join("| "),
    )
  }
  unigene_write_line(output, "GNM_TERMINUS", self.genomic_terminus)
  unigene_write_line(output, "CHROMOSOME", self.chromosome)
  for similarity in self.protein_similarities {
    unigene_write_line(output, "PROTSIM", similarity.to_value())
  }
  for transcript_map in self.transcript_maps {
    unigene_write_line(output, "TXMAP", transcript_map.to_value())
  }
  unigene_write_line(output, "SCOUNT", self.sequences.length().to_string())
  for sequence in self.sequences {
    unigene_write_line(output, "SEQUENCE", sequence.to_value())
  }
  for site in self.sts {
    unigene_write_line(output, "STS", site.to_value())
  }
  output.write_string("//\n")
  output.to_string()
}

///|
/// Parse all records in UniGene flat-file text.
pub fn unigene_parse(text : String) -> Array[UniGeneRecord] raise UniGeneError {
  let records : Array[UniGeneRecord] = []
  let mut current : UniGeneBuilder? = None
  for line_index, line_view in text.split("\n") {
    let mut line = line_view.to_owned()
    if line.length() > 0 &&
      line.unsafe_get(line.length() - 1) == '\r'.to_int().to_uint16() {
      line = line[0:line.length() - 1].to_owned()
    }
    if line.trim().length() == 0 {
      continue
    }
    if line.trim().to_owned() == "//" {
      match current {
        Some(builder) => {
          records.push(builder.finish())
          current = None
        }
        None =>
          raise UniGeneError(
            "Unexpected UniGene terminator at line " +
            (line_index + 1).to_string(),
          )
      }
      continue
    }
    let (tag, value) = unigene_split_line(line, line_number=line_index + 1)
    if tag == "ID" {
      match current {
        Some(_) =>
          raise UniGeneError(
            "UniGene record before line " +
            (line_index + 1).to_string() +
            " is missing its // terminator",
          )
        None => ()
      }
      if value.length() == 0 {
        raise UniGeneError(
          "Empty UniGene ID at line " + (line_index + 1).to_string(),
        )
      }
      let builder = UniGeneBuilder::new()
      builder.identifier = value
      current = Some(builder)
    } else {
      match current {
        Some(builder) =>
          unigene_apply_field(builder, tag, value, line_number=line_index + 1)
        None =>
          raise UniGeneError(
            "UniGene field before ID at line " + (line_index + 1).to_string(),
          )
      }
    }
  }
  match current {
    Some(_) =>
      raise UniGeneError(
        "Unexpected end of UniGene input: record is missing its // terminator",
      )
    None => ()
  }
  records
}

///|
/// Read exactly one UniGene record.
pub fn unigene_read(text : String) -> UniGeneRecord raise UniGeneError {
  let records = unigene_parse(text)
  if records.length() == 0 {
    raise UniGeneError("No UniGene record found")
  } else if records.length() > 1 {
    raise UniGeneError("More than one UniGene record found")
  }
  records[0]
}

///|
/// Return a compact valid UniGene record for examples and tests.
pub fn unigene_sample_text() -> String {
  "ID          Hs.12345\n" +
  "TITLE       Example kinase family member\n" +
  "GENE        EXK1\n" +
  "CYTOBAND    7q31.2\n" +
  "GENE_ID     12345\n" +
  "LOCUSLINK   12345\n" +
  "HOMOL       YES\n" +
  "EXPRESS     brain| heart| adult\n" +
  "RESTR_EXPR  brain| adult\n" +
  "GNM_TERMINUS T\n" +
  "CHROMOSOME  7\n" +
  "PROTSIM     ORG=9606; PROTGI=123456; PROTID=NP_000001.1; PCT=99.50; ALN=401\n" +
  "TXMAP       MARKER=D7S1234; RHPANEL=GB4\n" +
  "SCOUNT      3\n" +
  "SEQUENCE    ACC=NM_000001.2; NID=g123455; PID=g123456; SEQTYPE=mRNA\n" +
  "SEQUENCE    ACC=AA000001.1; NID=g223456; CLONE=IMAGE:123456; END=5'; LID=100; SEQTYPE=EST; TRACE=900001\n" +
  "SEQUENCE    ACC=AA000002.1; NID=g223457; CLONE=clone-2; END=3'; LID=100; MGC=600; SEQTYPE=EST\n" +
  "STS         ACC=G12345 UNISTS=76543\n" +
  "//\n"
}

///|
priv struct UniGeneBuilder {
  mut identifier : String
  mut title : String
  mut symbol : String
  mut cytoband : String
  mut expression : Array[String]
  mut restricted_expression : Array[String]
  mut genomic_terminus : String
  mut gene_id : String
  mut locuslink : String
  mut homology : Bool?
  mut chromosome : String
  protein_similarities : Array[UniGeneProteinSimilarity]
  sequences : Array[UniGeneSequence]
  sts : Array[UniGeneSTS]
  transcript_maps : Array[UniGeneTranscriptMap]
  mut declared_sequence_count : Int
  mut has_sequence_count : Bool
}

///|
fn UniGeneBuilder::new() -> UniGeneBuilder {
  UniGeneBuilder::{
    identifier: "",
    title: "",
    symbol: "",
    cytoband: "",
    expression: [],
    restricted_expression: [],
    genomic_terminus: "",
    gene_id: "",
    locuslink: "",
    homology: None,
    chromosome: "",
    protein_similarities: [],
    sequences: [],
    sts: [],
    transcript_maps: [],
    declared_sequence_count: 0,
    has_sequence_count: false,
  }
}

///|
fn UniGeneBuilder::finish(
  self : UniGeneBuilder,
) -> UniGeneRecord raise UniGeneError {
  if !self.has_sequence_count {
    raise UniGeneError(
      "UniGene record " + self.identifier + " is missing its SCOUNT field",
    )
  }
  if self.declared_sequence_count != self.sequences.length() {
    raise UniGeneError(
      "UniGene record " +
      self.identifier +
      " declares " +
      self.declared_sequence_count.to_string() +
      " sequences but contains " +
      self.sequences.length().to_string(),
    )
  }
  UniGeneRecord::{
    identifier: self.identifier,
    species: unigene_species(self.identifier),
    title: self.title,
    symbol: self.symbol,
    cytoband: self.cytoband,
    expression: self.expression,
    restricted_expression: self.restricted_expression,
    genomic_terminus: self.genomic_terminus,
    gene_id: self.gene_id,
    locuslink: self.locuslink,
    homology: self.homology,
    chromosome: self.chromosome,
    protein_similarities: self.protein_similarities,
    sequences: self.sequences,
    sts: self.sts,
    transcript_maps: self.transcript_maps,
    declared_sequence_count: self.declared_sequence_count,
  }
}

///|
fn unigene_apply_field(
  builder : UniGeneBuilder,
  tag : String,
  value : String,
  line_number~ : Int,
) -> Unit raise UniGeneError {
  match tag {
    "TITLE" => builder.title = value
    "GENE" => builder.symbol = value
    "CYTOBAND" => builder.cytoband = value
    "EXPRESS" => builder.expression = unigene_split_list(value, "|")
    "RESTR_EXPR" =>
      builder.restricted_expression = unigene_split_list(value, "|")
    "GNM_TERMINUS" => builder.genomic_terminus = value
    "GENE_ID" => builder.gene_id = value
    "LOCUSLINK" => builder.locuslink = value
    "HOMOL" =>
      if value == "YES" {
        builder.homology = Some(true)
      } else if value == "NO" {
        builder.homology = Some(false)
      } else {
        raise UniGeneError(
          "Invalid HOMOL value at line " +
          line_number.to_string() +
          ": expected YES or NO",
        )
      }
    "CHROMOSOME" => builder.chromosome = value
    "PROTSIM" =>
      builder.protein_similarities.push(
        unigene_parse_protein_similarity(value, line_number~),
      )
    "TXMAP" =>
      builder.transcript_maps.push(
        unigene_parse_transcript_map(value, line_number~),
      )
    "SCOUNT" => {
      if builder.has_sequence_count {
        raise UniGeneError(
          "Duplicate SCOUNT field at line " + line_number.to_string(),
        )
      }
      builder.declared_sequence_count = unigene_parse_nonnegative_int(
        value,
        context="SCOUNT at line " + line_number.to_string(),
      )
      builder.has_sequence_count = true
    }
    "SEQUENCE" =>
      builder.sequences.push(unigene_parse_sequence(value, line_number~))
    "STS" => builder.sts.push(unigene_parse_sts(value, line_number~))
    _ =>
      raise UniGeneError(
        "Unknown UniGene tag " + tag + " at line " + line_number.to_string(),
      )
  }
}

///|
fn unigene_parse_sequence(
  value : String,
  line_number~ : Int,
) -> UniGeneSequence raise UniGeneError {
  let fields = unigene_parse_fields(
    value,
    separator=";",
    context="SEQUENCE at line " + line_number.to_string(),
  )
  let extras : Array[UniGeneExtraField] = []
  for field in fields {
    if !unigene_is_known_sequence_key(field.key) {
      extras.push(field)
    }
  }
  UniGeneSequence::new(
    accession=unigene_field_or_empty(fields, "ACC"),
    nucleotide_id=unigene_field_or_empty(fields, "NID"),
    protein_id=unigene_field_or_empty(fields, "PID"),
    clone_id=unigene_field_or_empty(fields, "CLONE"),
    read_end=unigene_field_or_empty(fields, "END"),
    library_id=unigene_field_or_empty(fields, "LID"),
    mgc=unigene_field_or_empty(fields, "MGC"),
    sequence_type=unigene_field_or_empty(fields, "SEQTYPE"),
    trace=unigene_field_or_empty(fields, "TRACE"),
    extra_fields=extras,
  )
}

///|
fn unigene_parse_protein_similarity(
  value : String,
  line_number~ : Int,
) -> UniGeneProteinSimilarity raise UniGeneError {
  let fields = unigene_parse_fields(
    value,
    separator=";",
    context="PROTSIM at line " + line_number.to_string(),
  )
  let extras : Array[UniGeneExtraField] = []
  for field in fields {
    if !unigene_is_known_protein_key(field.key) {
      extras.push(field)
    }
  }
  UniGeneProteinSimilarity::new(
    organism=unigene_field_or_empty(fields, "ORG"),
    protein_gi=unigene_field_or_empty(fields, "PROTGI"),
    protein_id=unigene_field_or_empty(fields, "PROTID"),
    percent=unigene_field_or_empty(fields, "PCT"),
    alignment_length=unigene_field_or_empty(fields, "ALN"),
    extra_fields=extras,
  )
}

///|
fn unigene_parse_sts(
  value : String,
  line_number~ : Int,
) -> UniGeneSTS raise UniGeneError {
  let fields = unigene_parse_fields(
    value,
    separator=" ",
    context="STS at line " + line_number.to_string(),
  )
  let extras : Array[UniGeneExtraField] = []
  for field in fields {
    if field.key != "ACC" && field.key != "UNISTS" {
      extras.push(field)
    }
  }
  UniGeneSTS::new(
    accession=unigene_field_or_empty(fields, "ACC"),
    unists=unigene_field_or_empty(fields, "UNISTS"),
    extra_fields=extras,
  )
}

///|
fn unigene_parse_transcript_map(
  value : String,
  line_number~ : Int,
) -> UniGeneTranscriptMap raise UniGeneError {
  let fields = unigene_parse_fields(
    value,
    separator=";",
    context="TXMAP at line " + line_number.to_string(),
  )
  let extras : Array[UniGeneExtraField] = []
  for field in fields {
    if field.key != "MARKER" && field.key != "RHPANEL" {
      extras.push(field)
    }
  }
  UniGeneTranscriptMap::new(
    marker=unigene_field_or_empty(fields, "MARKER"),
    radiation_hybrid_panel=unigene_field_or_empty(fields, "RHPANEL"),
    extra_fields=extras,
  )
}

///|
fn unigene_split_line(
  line : String,
  line_number~ : Int,
) -> (String, String) raise UniGeneError {
  if line.length() < 12 {
    raise UniGeneError(
      "Malformed fixed-width UniGene line at " + line_number.to_string(),
    )
  }
  let tag = line[0:12].trim().to_owned()
  let value = line[12:].trim().to_owned()
  if tag.length() == 0 {
    raise UniGeneError("Missing UniGene tag at line " + line_number.to_string())
  }
  (tag, value)
}

///|
fn unigene_parse_fields(
  value : String,
  separator~ : String,
  context~ : String,
) -> Array[UniGeneExtraField] raise UniGeneError {
  let result : Array[UniGeneExtraField] = []
  for part_view in value.split(separator) {
    let part = part_view.trim().to_owned()
    if part.length() == 0 {
      continue
    }
    match part.find("=") {
      Some(position) => {
        let key = part[0:position].trim().to_owned()
        let field_value = part[position + 1:].trim().to_owned()
        if key.length() == 0 || field_value.length() == 0 {
          raise UniGeneError("Malformed key-value field in " + context)
        }
        result.push(UniGeneExtraField::new(key~, value=field_value))
      }
      None =>
        raise UniGeneError(
          "Malformed key-value field '" + part + "' in " + context,
        )
    }
  }
  result
}

///|
fn unigene_split_list(value : String, separator : String) -> Array[String] {
  let result : Array[String] = []
  for part in value.split(separator) {
    let item = part.trim().to_owned()
    if item.length() > 0 {
      result.push(item)
    }
  }
  result
}

///|
fn unigene_species(identifier : String) -> String {
  match identifier.find(".") {
    Some(position) => identifier[0:position].to_owned()
    None => identifier
  }
}

///|
fn unigene_image_id(clone_id : String) -> String {
  if clone_id.has_prefix("IMAGE:") && clone_id.length() > 6 {
    clone_id[6:].trim().to_owned()
  } else if clone_id.has_prefix("IMAGE ") && clone_id.length() > 6 {
    clone_id[6:].trim().to_owned()
  } else {
    ""
  }
}

///|
fn unigene_field_or_empty(
  fields : Array[UniGeneExtraField],
  key : String,
) -> String {
  match unigene_extra_value(fields, key) {
    Some(value) => value
    None => ""
  }
}

///|
fn unigene_extra_value(
  fields : Array[UniGeneExtraField],
  key : String,
) -> String? {
  for field in fields {
    if field.key == key {
      return Some(field.value)
    }
  }
  None
}

///|
fn unigene_nonempty(value : String) -> String? {
  if value.length() > 0 {
    Some(value)
  } else {
    None
  }
}

///|
fn unigene_is_known_sequence_key(key : String) -> Bool {
  key == "ACC" ||
  key == "NID" ||
  key == "PID" ||
  key == "CLONE" ||
  key == "END" ||
  key == "LID" ||
  key == "MGC" ||
  key == "SEQTYPE" ||
  key == "TRACE"
}

///|
fn unigene_is_known_protein_key(key : String) -> Bool {
  key == "ORG" ||
  key == "PROTGI" ||
  key == "PROTID" ||
  key == "PCT" ||
  key == "ALN"
}

///|
fn unigene_parse_nonnegative_int(
  value : String,
  context~ : String,
) -> Int raise UniGeneError {
  if value.length() == 0 {
    raise UniGeneError("Missing integer value for " + context)
  }
  for index = 0; index < value.length(); index = index + 1 {
    let character = value.unsafe_get(index).to_int()
    if character < '0'.to_int() || character > '9'.to_int() {
      raise UniGeneError("Invalid non-negative integer for " + context)
    }
  }
  parse_int(value)
}

///|
fn unigene_push_field(
  fields : Array[String],
  key : String,
  value : String,
) -> Unit {
  if value.length() > 0 {
    fields.push(key + "=" + value)
  }
}

///|
fn unigene_write_line(
  output : StringBuilder,
  tag : String,
  value : String,
) -> Unit {
  if value.length() == 0 {
    return
  }
  output.write_string(tag)
  for index = tag.length(); index < 12; index = index + 1 {
    output.write_char(' ')
  }
  output.write_string(value)
  output.write_string("\n")
}
