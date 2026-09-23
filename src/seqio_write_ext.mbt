///|
/// SeqIO write-side extension adapters (design 2.2 (6)).
///
/// Each `write_*_records` function maps an `Array[SeqRecord]` onto the
/// existing record-level writers for its format; where no record-level
/// writer exists (pir, swiss, phd) a minimal serializer is implemented
/// here. Parsing logic is never duplicated and every format round-trips
/// through its unified parse counterpart.
///
/// Formats (spec 5.3):
///   - embl:       SeqRecord -> InsdcRecord -> write_embl_record
///   - pir:        minimal PIR/NBRF flat writer (">P1;<id>" headers)
///   - fasta-2line: delegates to write_fasta_2line (wired in seqio.mbt)
///   - seqxml:     seqxml_from_seq_records + seqxml_to_xml
///   - swiss:      minimal Swiss-Prot flat writer
///   - ig:         ig_from_seq_records + ig_records_to_string
///   - imgt:       SeqRecord -> ImgtRecord -> imgt_write
///   - phd:        SeqRecord -> PhdFile-equivalent blocks, requires
///                 letter_annotations["phred_quality"]
///
/// Empty record arrays produce an empty string, consistent with the
/// existing write_fasta behaviour (spec 5.3.3.3).

///|
/// Read a record annotation with a fallback default.
fn writeext_annotation(
  rec : SeqRecord,
  key : String,
  fallback : String,
) -> String {
  match rec.annotations.get(key) {
    Some(v) => v
    None => fallback
  }
}

///|
/// Split a comma-joined annotation value back into its items.
fn writeext_split_annotation(value : String) -> Array[String] {
  let items : Array[String] = Array(capacity=0)
  for part in value.split(", ") {
    let trimmed = part.trim().to_string()
    if trimmed.length() > 0 {
      items.push(trimmed)
    }
  }
  items
}

///|
/// Write records in Swiss-Prot flat file format. A minimal entry is
/// emitted per record (ID/AC/DE/OS/KW/SQ plus the sequence block), which
/// round-trips through parse_swiss_records.
pub fn write_swiss_records(
  records : Array[SeqRecord],
) -> String raise SeqIOError {
  let buf = StringBuilder::new()
  for rec in records {
    let seq = rec.seq.to_string()
    buf.write_string("ID   ")
    buf.write_string(rec.id)
    buf.write_string("              Standard;       ")
    buf.write_string(seq.length().to_string())
    buf.write_string(" AA.\n")
    let accessions = writeext_split_annotation(
      writeext_annotation(rec, "accessions", ""),
    )
    if accessions.length() > 0 {
      buf.write_string("AC   ")
      buf.write_string(accessions.join("; "))
      buf.write_string(";\n")
    }
    if rec.description.length() > 0 {
      buf.write_string("DE   RecName: Full=")
      buf.write_string(rec.description)
      buf.write_string(";\n")
    }
    let organism = writeext_annotation(rec, "organism", "")
    if organism.length() > 0 {
      buf.write_string("OS   ")
      buf.write_string(organism)
      buf.write_string(".\n")
    }
    let keywords = writeext_split_annotation(
      writeext_annotation(rec, "keywords", ""),
    )
    if keywords.length() > 0 {
      buf.write_string("KW   ")
      buf.write_string(keywords.join("; "))
      buf.write_string(".\n")
    }
    buf.write_string("SQ   SEQUENCE   ")
    buf.write_string(seq.length().to_string())
    buf.write_string(" AA;\n")
    // Sequence in 10-character groups, 60 characters per line.
    let mut i = 0
    while i < seq.length() {
      let groups : Array[String] = Array(capacity=0)
      let mut j = i
      while j < seq.length() && j < i + 60 {
        let end = if j + 10 < seq.length() && j + 10 <= i + 60 {
          j + 10
        } else if j + 10 <= i + 60 {
          seq.length()
        } else {
          i + 60
        }
        groups.push(seq[j:end].to_owned())
        j = end
      }
      buf.write_string("     ")
      buf.write_string(groups.join(" "))
      buf.write_char('\n')
      i = i + 60
    }
    buf.write_string("//\n")
  }
  buf.to_string()
}

///|
/// Write records in PIR/NBRF format: ">P1;<id>" header, description line,
/// and the sequence wrapped at 60 columns with a trailing '*'.
pub fn write_pir_records(records : Array[SeqRecord]) -> String raise SeqIOError {
  let buf = StringBuilder::new()
  for rec in records {
    let seq = rec.seq.to_string()
    buf.write_string(">P1;")
    buf.write_string(rec.id)
    buf.write_char('\n')
    buf.write_string(rec.description)
    buf.write_char('\n')
    let mut i = 0
    while i < seq.length() {
      let end = if i + 60 < seq.length() { i + 60 } else { seq.length() }
      buf.write_string(seq[i:end].to_owned())
      if end == seq.length() {
        buf.write_string("*")
      }
      buf.write_char('\n')
      i = end
    }
  }
  buf.to_string()
}

///|
/// Render a SeqFeature location as an INSDC location string.
fn writeext_location_string(feature : SeqFeature) -> String {
  match feature.location {
    Some(loc) => {
      let start = loc.start_int()
      let end = loc.end_int()
      if start == end {
        start.to_string()
      } else {
        start.to_string() + ".." + end.to_string()
      }
    }
    None => ""
  }
}

///|
/// Write records in EMBL format by adapting each SeqRecord onto an
/// InsdcRecord and reusing the existing record-level writer. Header
/// metadata is sourced from annotations; features are preserved through
/// the InsdcFeature feature table.
pub fn write_embl_records(
  records : Array[SeqRecord],
) -> String raise SeqIOError {
  let buf = StringBuilder::new()
  for rec in records {
    let seq = rec.seq.to_string()
    let features : Array[InsdcFeature] = Array(capacity=0)
    for feature in rec.features {
      let location = writeext_location_string(feature)
      if location.length() == 0 {
        continue
      }
      features.push(InsdcFeature::{
        feature_type: feature.feature_type,
        location_str: location,
        qualifiers: feature.qualifiers,
      })
    }
    let insdc = InsdcRecord::{
      locus_name: rec.name,
      sequence_length: seq.length(),
      molecule_type: writeext_annotation(rec, "molecule_type", "DNA"),
      topology: writeext_annotation(rec, "topology", "linear"),
      division: writeext_annotation(rec, "data_file_division", "UNC"),
      date: writeext_annotation(rec, "date", "01-JAN-1980"),
      definition: rec.description,
      accession: writeext_annotation(rec, "accession", rec.id),
      version: writeext_annotation(rec, "version", rec.id),
      organism: writeext_annotation(rec, "organism", ""),
      features,
      sequence: seq,
    }
    buf.write_string(write_embl_record(insdc))
  }
  buf.to_string()
}

///|
/// Write records in IMGT format by rebuilding the 13-field header from
/// annotations (the same keys parse_imgt_records emits) and reusing the
/// existing imgt_write serializer.
pub fn write_imgt_records(
  records : Array[SeqRecord],
) -> String raise SeqIOError {
  let imgt_records : Array[ImgtRecord] = Array(capacity=0)
  for rec in records {
    let header = ImgtHeader::new(
      writeext_annotation(rec, "accession", rec.id),
      writeext_annotation(rec, "seq_id", rec.id),
      writeext_annotation(rec, "official_name", rec.name),
      writeext_annotation(rec, "region", ""),
      writeext_annotation(rec, "date_update", ""),
      writeext_annotation(rec, "date_creation", ""),
      writeext_annotation(rec, "species", ""),
      writeext_annotation(rec, "gene", ""),
      writeext_annotation(rec, "imgt_group", ""),
      writeext_annotation(rec, "locus", ""),
      rec.seq.length().to_string(),
      writeext_annotation(rec, "seq_type", ""),
      writeext_annotation(rec, "comments", ""),
    )
    imgt_records.push(ImgtRecord::new(header, rec.seq.to_string()))
  }
  imgt_write(imgt_records)
}

///|
/// Write records in IntelliGenetics format, reusing the existing
/// SeqRecord -> IgRecord conversion and file serializer.
pub fn write_ig_records(records : Array[SeqRecord]) -> String raise SeqIOError {
  ig_records_to_string(ig_from_seq_records(records), None)
}

///|
/// Write records in SeqXML format, reusing the existing SeqRecord ->
/// SeqXmlDocument conversion and serializer.
pub fn write_seqxml_records(
  records : Array[SeqRecord],
) -> String raise SeqIOError {
  seqxml_to_xml(seqxml_from_seq_records(records))
}

///|
/// Write records in PHD format. Every record must carry per-base Phred
/// qualities in letter_annotations["phred_quality"]; a missing (or
/// length-mismatched) quality track raises before any output text is
/// produced for that record (spec 5.3.1.4 / 5.3.3.2).
pub fn write_phd_records(records : Array[SeqRecord]) -> String raise SeqIOError {
  let buf = StringBuilder::new()
  for rec in records {
    let seq = rec.seq.to_string()
    let qualities = match rec.phred_quality() {
      Some(q) => q
      None =>
        raise SeqIOError("record '" + rec.id + "' lacks phred quality scores")
    }
    if qualities.length() != seq.length() {
      raise SeqIOError(
        "record '" +
        rec.id +
        "' phred quality length does not match sequence length",
      )
    }
    buf.write_string("BEGIN_SEQUENCE ")
    buf.write_string(rec.id)
    buf.write_string("\n\nBEGIN_COMMENT\n")
    buf.write_string("CHROMAT_FILE: \n")
    buf.write_string("PHRED_VERSION: \n")
    buf.write_string("PHRED_PARAMETER_FILE: \n")
    buf.write_string("TRACE_INDEX_MIN: 0\n")
    buf.write_string("TRACE_INDEX_MAX: 0\n")
    buf.write_string("TRIM: 0 0 0.0\n")
    buf.write_string("CHEMISTRY: \n")
    buf.write_string("DYE: \n")
    buf.write_string("END_COMMENT\n\nBEGIN_DNA\n")
    for i in 0..<seq.length() {
      let base = seq.unsafe_get(i).unsafe_to_char().to_string()
      buf.write_string(base)
      buf.write_char(' ')
      buf.write_string(qualities[i].to_string())
      buf.write_string(" 0\n")
    }
    buf.write_string("END_DNA\n\nEND_SEQUENCE\n")
  }
  buf.to_string()
}
