///|
/// SeqIO text/container unification adapters.
///
/// Each `parse_*_records` function delegates to the existing standalone
/// parser for its format and normalizes the parser-specific record types
/// into `SeqRecord`, so that the unified `seqio_parse` entry point returns
/// the same data model for every format. Parsing logic is never duplicated:
/// the adapters only map fields and wrap error messages with the format
/// name.
///
/// Annotation normalization rules (design 1.2.3):
///   - swiss:        description/keywords/organism/... go to annotations,
///                   features become SeqFeature entries
///   - uniprot-xml:  bridged through the existing uniprot_to_seqrecord
///   - seqxml:       bridged through the existing seqxml_to_seq_records
///   - imgt:         the 13 header fields go to annotations
///   - ig:           bridged through the existing ig_to_seq_records
///   - insdc:        GenBank/EMBL records keep their feature table through
///                   SeqRecord.features
///   - phd:          per-read base qualities go to
///                   letter_annotations["phred_quality"]
///   - ace:          one record per contig consensus, read counts etc. in
///                   annotations
///   - gfa:          bridged through gfa_get_segments_as_records
///   - pdb-seqres / pdb-atom: parsers already produce SeqRecord
///
/// Empty content (or content without records) yields an empty array, not
/// an error (spec 5.1.3.4).

///|
/// Set an annotation key only when the value is non-empty (keeps the
/// annotation map clean and deterministic).
fn unify_set_annotation(
  annotations : Map[String, String],
  key : String,
  value : String,
) -> Unit {
  if value.length() > 0 {
    annotations[key] = value
  }
}

///|
/// Normalize a Swiss-Prot record into a SeqRecord.
fn unify_swiss_to_record(entry : SwissProtRecord) -> SeqRecord raise SeqIOError {
  let annotations : Map[String, String] = Map([], capacity=16)
  unify_set_annotation(
    annotations,
    "accessions",
    entry.accession_numbers.join(", "),
  )
  unify_set_annotation(annotations, "gene_names", entry.gene_names.join(", "))
  unify_set_annotation(annotations, "organism", entry.organism)
  unify_set_annotation(annotations, "organism_id", entry.organism_id)
  unify_set_annotation(annotations, "taxonomy", entry.taxonomy.join(", "))
  unify_set_annotation(annotations, "keywords", entry.keywords.join(", "))
  unify_set_annotation(annotations, "seq_length", entry.seq_length.to_string())
  unify_set_annotation(annotations, "seq_mass", entry.seq_mass.to_string())
  unify_set_annotation(annotations, "ec_numbers", entry.ec_numbers.join(", "))
  unify_set_annotation(annotations, "domains", entry.domains.join(", "))
  unify_set_annotation(annotations, "function", entry.function)
  unify_set_annotation(
    annotations,
    "subcellular_location",
    entry.subcellular_location.join(", "),
  )
  unify_set_annotation(
    annotations,
    "tissue_specificity",
    entry.tissue_specificity,
  )
  unify_set_annotation(annotations, "date_created", entry.date_created)
  unify_set_annotation(annotations, "date_modified", entry.date_modified)
  unify_set_annotation(annotations, "version", entry.version.to_string())
  let features : Array[SeqFeature] = Array(capacity=0)
  for f in entry.features {
    let quals : Map[String, Array[String]] = Map(
      [],
      capacity=f.qualifiers.length(),
    )
    for k, v in f.qualifiers {
      quals[k] = [v]
    }
    let location : FeatureLocation? = if f.start > 0 || f.end > 0 {
      match (try? FeatureLocation::new_exact(f.start, f.end)) {
        Ok(loc) => Some(loc)
        Err(_) => None
      }
    } else {
      None
    }
    features.push(
      SeqFeature::new(
        location~,
        feature_type=f.feature_type,
        id=f.description,
        qualifiers=quals,
      ),
    )
  }
  let record = SeqRecord::new(
    Seq::new(entry.sequence),
    id=entry.entry_name,
    name=entry.entry_name,
    description=entry.protein_name,
  )
  for k, v in annotations {
    record.annotations[k] = v
  }
  for feature in features {
    record.features.push(feature)
  }
  record
}

///|
/// Parse Swiss-Prot flat file content into SeqRecords (one record per
/// "//"-delimited entry). Raises "Malformed swiss entry" when no entry can
/// be parsed.
pub fn parse_swiss_records(
  content : String,
) -> Array[SeqRecord] raise SeqIOError {
  if content.trim() == "" {
    return []
  }
  let entries = parse_swissprot_all(content)
  if entries.length() == 0 {
    raise SeqIOError("Malformed swiss entry")
  }
  let records : Array[SeqRecord] = Array(capacity=0)
  for entry in entries {
    records.push(unify_swiss_to_record(entry))
  }
  records
}

///|
/// Parse UniProt XML content into SeqRecords, bridging the existing
/// `parse_uniprot_xml` + `uniprot_to_seqrecord` pair.
pub fn parse_uniprot_xml_records(
  content : String,
) -> Array[SeqRecord] raise SeqIOError {
  if content.trim() == "" {
    return []
  }
  let entries = parse_uniprot_xml(content)
  let records : Array[SeqRecord] = Array(capacity=0)
  for entry in entries {
    let (id, description, sequence, annotations) = uniprot_to_seqrecord(entry)
    let record = SeqRecord::new(Seq::new(sequence), id~, name=id, description~)
    for k, v in annotations {
      record.annotations[k] = v
    }
    records.push(record)
  }
  records
}

///|
/// Parse SeqXML content into SeqRecords, bridging the existing
/// `seqxml_parse` + `seqxml_to_seq_records` pair.
pub fn parse_seqxml_records(
  content : String,
) -> Array[SeqRecord] raise SeqIOError {
  if content.trim() == "" {
    return []
  }
  seqxml_to_seq_records(seqxml_parse(content))
}

///|
/// Parse IMGT content into SeqRecords: one record per IMGT entry with the
/// 13 header fields mirrored into annotations (spec design 1.2.3).
pub fn parse_imgt_records(
  content : String,
) -> Array[SeqRecord] raise SeqIOError {
  if content.trim() == "" {
    return []
  }
  let entries = imgt_parse(content)
  let records : Array[SeqRecord] = Array(capacity=0)
  for entry in entries {
    let header = entry.header()
    let record = SeqRecord::new(
      Seq::new(entry.sequence()),
      id=entry.id(),
      name=entry.id(),
      description=entry.description(),
    )
    unify_set_annotation(record.annotations, "accession", header.accession)
    unify_set_annotation(record.annotations, "seq_id", header.seq_id)
    unify_set_annotation(
      record.annotations,
      "official_name",
      header.official_name,
    )
    unify_set_annotation(record.annotations, "region", header.region)
    unify_set_annotation(record.annotations, "date_update", header.date_update)
    unify_set_annotation(
      record.annotations,
      "date_creation",
      header.date_creation,
    )
    unify_set_annotation(record.annotations, "species", header.species)
    unify_set_annotation(record.annotations, "gene", header.gene)
    unify_set_annotation(record.annotations, "imgt_group", header.imgt_group)
    unify_set_annotation(record.annotations, "locus", header.locus)
    unify_set_annotation(record.annotations, "seq_length", header.seq_length)
    unify_set_annotation(record.annotations, "seq_type", header.seq_type)
    unify_set_annotation(record.annotations, "comments", header.comments)
    records.push(record)
  }
  records
}

///|
/// Parse IntelliGenetics content into SeqRecords, bridging the existing
/// `ig_parse` + `ig_to_seq_records` pair.
pub fn parse_ig_records(content : String) -> Array[SeqRecord] raise SeqIOError {
  if content.trim() == "" {
    return []
  }
  ig_to_seq_records(ig_parse(content))
}

///|
/// Split multi-record INSDC content (GenBank or EMBL) on the "//" record
/// terminator lines. Every chunk keeps its trailing terminator so the
/// underlying single-record parsers behave unchanged.
fn unify_split_insdc(content : String) -> Array[String] {
  let raw_lines = content.split("\n").to_array()
  let chunks : Array[String] = Array(capacity=0)
  let mut buf = StringBuilder()
  for raw in raw_lines {
    let len = raw.length()
    let line = if len > 0 &&
      raw.unsafe_get(len - 1) == '\r'.to_int().to_uint16() {
      raw[0:len - 1].to_owned()
    } else {
      raw.to_owned()
    }
    let trimmed = line.trim()
    buf.write_string(line)
    buf.write_char('\n')
    if trimmed == "//" {
      chunks.push(buf.to_string())
      buf = StringBuilder()
    }
  }
  let tail = buf.to_string()
  if tail.trim() != "" {
    chunks.push(tail)
  }
  chunks
}

///|
/// Normalize an INSDC (GenBank/EMBL) record into a SeqRecord, keeping the
/// feature table accessible through SeqRecord.features.
fn unify_insdc_to_record(record : InsdcRecord) -> SeqRecord {
  let id = if record.accession.length() > 0 {
    record.accession
  } else {
    record.locus_name
  }
  let seq_record = SeqRecord::new(
    Seq::new(record.sequence),
    id~,
    name=record.locus_name,
    description=record.definition,
  )
  unify_set_annotation(
    seq_record.annotations,
    "sequence_length",
    record.sequence_length.to_string(),
  )
  unify_set_annotation(
    seq_record.annotations,
    "molecule_type",
    record.molecule_type,
  )
  unify_set_annotation(seq_record.annotations, "topology", record.topology)
  unify_set_annotation(seq_record.annotations, "division", record.division)
  unify_set_annotation(seq_record.annotations, "date", record.date)
  unify_set_annotation(seq_record.annotations, "accession", record.accession)
  unify_set_annotation(seq_record.annotations, "version", record.version)
  unify_set_annotation(seq_record.annotations, "organism", record.organism)
  for feature in record.features {
    let quals : Map[String, Array[String]] = Map(
      [],
      capacity=feature.qualifiers.length(),
    )
    for k, v in feature.qualifiers {
      quals[k] = v.copy()
    }
    seq_record.features.push(
      SeqFeature::new(
        feature_type=feature.feature_type,
        id=feature.location_str,
        qualifiers=quals,
      ),
    )
  }
  seq_record
}

///|
/// Parse INSDC content (GenBank or EMBL records, mixed per record) into
/// SeqRecords. Each record must start with a LOCUS (GenBank) or ID (EMBL)
/// header line; the feature table is preserved through SeqRecord.features.
pub fn parse_insdc_records(
  content : String,
) -> Array[SeqRecord] raise SeqIOError {
  if content.trim() == "" {
    return []
  }
  let chunks = unify_split_insdc(content)
  let records : Array[SeqRecord] = Array(capacity=0)
  for chunk in chunks {
    let mut first = ""
    for line in chunk.split("\n") {
      let trimmed = line.trim().to_owned()
      if trimmed != "" {
        first = trimmed
        break
      }
    }
    if first == "" {
      continue
    }
    let record = if first.unsafe_get(0) == 'L' && first.has_prefix("LOCUS") {
      unify_insdc_to_record(parse_genbank_record(chunk))
    } else if first.has_prefix("ID") {
      unify_insdc_to_record(parse_embl_record(chunk))
    } else {
      raise SeqIOError("Unrecognized INSDC record: " + first)
    }
    records.push(record)
  }
  records
}

///|
/// Parse PHD content into SeqRecords: one record per read with per-base
/// qualities stored in letter_annotations["phred_quality"].
pub fn parse_phd_records(content : String) -> Array[SeqRecord] raise SeqIOError {
  if content.trim() == "" {
    return []
  }
  let file = phd_parse(content)
  let records : Array[SeqRecord] = Array(capacity=0)
  for read in file.reads {
    let bases = StringBuilder()
    let qualities : Array[Int] = Array(capacity=0)
    for base in read.bases {
      bases.write_string(base.base)
      qualities.push(base.quality)
    }
    let record = SeqRecord::new(
      Seq::new(bases.to_string()),
      id=read.name,
      name=read.name,
      description=read.name,
    )
    record.letter_annotations["phred_quality"] = qualities
    records.push(record)
  }
  records
}

///|
/// Parse ACE content into SeqRecords: one record per contig consensus
/// sequence, with the contig name as id and read counts etc. in
/// annotations. Contig base qualities are attached as phred_quality when
/// their length matches the consensus.
pub fn parse_ace_records(content : String) -> Array[SeqRecord] raise SeqIOError {
  if content.trim() == "" {
    return []
  }
  let data = ace_parse(content)
  let records : Array[SeqRecord] = Array(capacity=0)
  for contig in data.contigs {
    let sequence = ace_consensus_sequence(contig)
    let record = SeqRecord::new(
      Seq::new(sequence),
      id=contig.contig_name,
      name=contig.contig_name,
      description=contig.contig_name,
    )
    unify_set_annotation(
      record.annotations,
      "reads_count",
      contig.reads.length().to_string(),
    )
    unify_set_annotation(
      record.annotations,
      "length",
      contig.length.to_string(),
    )
    unify_set_annotation(
      record.annotations,
      "reverse_complemented",
      if contig.reverse_complemented {
        "true"
      } else {
        "false"
      },
    )
    if contig.base_qualities.length() == sequence.length() &&
      sequence.length() > 0 {
      record.letter_annotations["phred_quality"] = contig.base_qualities.copy()
    }
    records.push(record)
  }
  records
}

///|
/// Parse GFA content into SeqRecords, bridging the existing
/// `gfa_parse` + `gfa_get_segments_as_records` pair.
pub fn parse_gfa_records(content : String) -> Array[SeqRecord] raise SeqIOError {
  if content.trim() == "" {
    return []
  }
  gfa_get_segments_as_records(gfa_parse(content))
}

///|
/// Parse PDB SEQRES content into SeqRecords (bridged; the underlying
/// parser already produces SeqRecord).
pub fn parse_pdb_seqres_records(
  content : String,
) -> Array[SeqRecord] raise SeqIOError {
  if content.trim() == "" {
    return []
  }
  pdb_seqres_parser(content)
}

///|
/// Parse PDB ATOM content into SeqRecords (bridged; the underlying parser
/// already produces SeqRecord).
pub fn parse_pdb_atom_records(
  content : String,
) -> Array[SeqRecord] raise SeqIOError {
  if content.trim() == "" {
    return []
  }
  pdb_atom_parser(content)
}
