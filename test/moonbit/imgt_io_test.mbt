///|
/// Tests for the IMGT I/O module (Bio.SeqIO.ImgtIO equivalent).
///
/// These tests exercise the public API exposed in `src/imgt_io.mbt`:
/// the ImgtHeader (13 pipe-separated fields) and ImgtRecord structs,
/// the parser (`imgt_parse`, `imgt_parse_header`, `imgt_parse_line_to_fields`),
/// the serializers (`imgt_header_to_string`, `imgt_record_to_string`,
/// `imgt_write`), the query helpers (`imgt_get_by_accession`,
/// `imgt_get_by_species`, `imgt_get_by_gene`, `imgt_filter_by_locus`,
/// `imgt_unique_species`, `imgt_unique_genes`, `imgt_summary`,
/// `imgt_n_records`), the SeqRecord integration (`imgt_to_seq_records` /
/// `imgt_from_seq_records`), and the sample data (`imgt_sample_records`).

// ============================================================================
// 1. ImgtHeader - creation, all 13 accessors, to_string
// ============================================================================

///|
test "imgt_header_new_and_all_accessors" {
  let h = @src.ImgtHeader::new(
    "HLA00001",
    "A*01:01:01:01",
    "A*01:01:01:01",
    "ORF",
    "2015-03-31",
    "1996-01-19",
    "human",
    "HLA-A",
    "I",
    "MHC",
    "2703 bp",
    " genomic-DNA",
    "no comments",
  )
  assert_eq(h.accession(), "HLA00001")
  assert_eq(h.seq_id(), "A*01:01:01:01")
  assert_eq(h.official_name(), "A*01:01:01:01")
  assert_eq(h.region(), "ORF")
  assert_eq(h.date_update(), "2015-03-31")
  assert_eq(h.date_creation(), "1996-01-19")
  assert_eq(h.species(), "human")
  assert_eq(h.gene(), "HLA-A")
  assert_eq(h.imgt_group(), "I")
  assert_eq(h.locus(), "MHC")
  assert_eq(h.seq_length(), "2703 bp")
  assert_eq(h.seq_type(), " genomic-DNA")
  assert_eq(h.comments(), "no comments")
}

///|
test "imgt_header_to_string" {
  let h = @src.ImgtHeader::new(
    "HLA00001",
    "A*01:01:01:01",
    "A*01:01:01:01",
    "ORF",
    "2015-03-31",
    "1996-01-19",
    "human",
    "HLA-A",
    "I",
    "MHC",
    "2703 bp",
    " genomic-DNA",
    "",
  )
  assert_eq(
    h.to_string(),
    "HLA00001|A*01:01:01:01|A*01:01:01:01|ORF|2015-03-31|1996-01-19|human|HLA-A|I|MHC|2703 bp| genomic-DNA|",
  )
}

///|
test "imgt_header_to_string_all_fields_populated" {
  let h = @src.ImgtHeader::new(
    "ACC1",
    "SID1",
    "ON1",
    "R1",
    "D1",
    "D2",
    "SP1",
    "G1",
    "IG1",
    "L1",
    "100 bp",
    "DNA",
    "C1",
  )
  let s = h.to_string()
  // 13 fields separated by 12 pipes.
  assert_eq(s, "ACC1|SID1|ON1|R1|D1|D2|SP1|G1|IG1|L1|100 bp|DNA|C1")
}

///|
test "imgt_header_to_string_empty_fields" {
  let h = @src.ImgtHeader::new("", "", "", "", "", "", "", "", "", "", "", "", "")
  assert_eq(h.to_string(), "||||||||||||")
}

// ============================================================================
// 2. ImgtRecord - creation, header, sequence, id, description, to_string
// ============================================================================

///|
test "imgt_record_new_and_accessors" {
  let header = @src.ImgtHeader::new(
    "HLA00001",
    "A*01:01:01:01",
    "A*01:01:01:01",
    "ORF",
    "2015-03-31",
    "1996-01-19",
    "human",
    "HLA-A",
    "I",
    "MHC",
    "2703 bp",
    " genomic-DNA",
    "",
  )
  let rec = @src.ImgtRecord::new(header, "ATGGCCGTCATGGCGCCCCGAACCCTCCTCCTC")
  assert_eq(rec.header().accession(), "HLA00001")
  assert_eq(rec.sequence(), "ATGGCCGTCATGGCGCCCCGAACCCTCCTCCTC")
}

///|
test "imgt_record_id" {
  let header = @src.ImgtHeader::new(
    "HLA00001",
    "A*01:01:01:01",
    "",
    "",
    "",
    "",
    "",
    "",
    "",
    "",
    "",
    "",
    "",
  )
  let rec = @src.ImgtRecord::new(header, "ACGT")
  assert_eq(rec.id(), "A*01:01:01:01")
}

///|
test "imgt_record_description" {
  let header = @src.ImgtHeader::new(
    "HLA00001",
    "A*01:01:01:01",
    "A*01:01:01:01",
    "ORF",
    "2015-03-31",
    "1996-01-19",
    "human",
    "HLA-A",
    "I",
    "MHC",
    "2703 bp",
    " genomic-DNA",
    "",
  )
  let rec = @src.ImgtRecord::new(header, "ACGT")
  // description() returns the full pipe-separated header (without '>').
  assert_eq(rec.description(), header.to_string())
}

///|
test "imgt_record_to_string" {
  let header = @src.ImgtHeader::new(
    "HLA00001",
    "A*01:01:01:01",
    "",
    "",
    "",
    "",
    "human",
    "HLA-A",
    "",
    "",
    "",
    "",
    "",
  )
  let rec = @src.ImgtRecord::new(header, "ATGGCC")
  let s = rec.to_string()
  assert_true(s.contains("ImgtRecord(id=A*01:01:01:01"))
  assert_true(s.contains("accession=HLA00001"))
  assert_true(s.contains("gene=HLA-A"))
  assert_true(s.contains("species=human"))
  assert_true(s.contains("length=6"))
}

// ============================================================================
// 3. imgt_parse_header - valid, invalid, short
// ============================================================================

///|
test "imgt_parse_header_valid" {
  let line = ">HLA00001|A*01:01:01:01|A*01:01:01:01|ORF|2015-03-31|1996-01-19|human|HLA-A|I|MHC|2703 bp| genomic-DNA|"
  match @src.imgt_parse_header(line) {
    Some(h) => {
      assert_eq(h.accession(), "HLA00001")
      assert_eq(h.seq_id(), "A*01:01:01:01")
      assert_eq(h.official_name(), "A*01:01:01:01")
      assert_eq(h.region(), "ORF")
      assert_eq(h.date_update(), "2015-03-31")
      assert_eq(h.date_creation(), "1996-01-19")
      assert_eq(h.species(), "human")
      assert_eq(h.gene(), "HLA-A")
      assert_eq(h.imgt_group(), "I")
      assert_eq(h.locus(), "MHC")
      assert_eq(h.seq_length(), "2703 bp")
      assert_eq(h.seq_type(), " genomic-DNA")
      assert_eq(h.comments(), "")
    }
    None => assert_true(false)
  }
}

///|
test "imgt_parse_header_invalid_no_gt" {
  // A line that does not start with '>' returns None.
  match @src.imgt_parse_header("HLA00001|A*01:01:01:01|ORF") {
    Some(_) => assert_true(false)
    None => assert_true(true)
  }
}

///|
test "imgt_parse_header_empty" {
  match @src.imgt_parse_header("") {
    Some(_) => assert_true(false)
    None => assert_true(true)
  }
}

///|
test "imgt_parse_header_short_header_pads_empty" {
  // A header with only 2 fields: the remaining 11 fields should be empty.
  match @src.imgt_parse_header(">HLA00001|A*01:01:01:01") {
    Some(h) => {
      assert_eq(h.accession(), "HLA00001")
      assert_eq(h.seq_id(), "A*01:01:01:01")
      assert_eq(h.official_name(), "")
      assert_eq(h.region(), "")
      assert_eq(h.species(), "")
      assert_eq(h.gene(), "")
      assert_eq(h.locus(), "")
      assert_eq(h.comments(), "")
    }
    None => assert_true(false)
  }
}

///|
test "imgt_parse_header_bare_gt" {
  // A header with just '>' and no fields.
  match @src.imgt_parse_header(">") {
    Some(h) => {
      assert_eq(h.accession(), "")
      assert_eq(h.seq_id(), "")
    }
    None => assert_true(false)
  }
}

// ============================================================================
// 4. imgt_parse_line_to_fields - split by |
// ============================================================================

///|
test "imgt_parse_line_to_fields_basic" {
  let fields = @src.imgt_parse_line_to_fields(">HLA00001|A*01:01:01:01|ORF")
  assert_eq(fields.length(), 3)
  assert_eq(fields[0], "HLA00001")
  assert_eq(fields[1], "A*01:01:01:01")
  assert_eq(fields[2], "ORF")
}

///|
test "imgt_parse_line_to_fields_without_gt" {
  // A leading '>' is stripped if present, but the function also works without it.
  let fields = @src.imgt_parse_line_to_fields("HLA00001|A*01:01:01:01")
  assert_eq(fields.length(), 2)
  assert_eq(fields[0], "HLA00001")
  assert_eq(fields[1], "A*01:01:01:01")
}

///|
test "imgt_parse_line_to_fields_empty" {
  let fields = @src.imgt_parse_line_to_fields("")
  assert_eq(fields.length(), 1)
  assert_eq(fields[0], "")
}

///|
test "imgt_parse_line_to_fields_trailing_pipe" {
  // A trailing pipe produces an extra empty field.
  let fields = @src.imgt_parse_line_to_fields(">a|b|")
  assert_eq(fields.length(), 3)
  assert_eq(fields[0], "a")
  assert_eq(fields[1], "b")
  assert_eq(fields[2], "")
}

// ============================================================================
// 5. imgt_parse - full content, single record, empty content
// ============================================================================

///|
test "imgt_parse_two_records" {
  let content = ">HLA00001|A*01:01:01:01|A*01:01:01:01|ORF|2015-03-31|1996-01-19|human|HLA-A|I|MHC|2703 bp| genomic-DNA|\nATGGCCGTCATGGCGCCCCGAACCCTCCTCCTC\n>HLA00003|B*07:02:01|B*07:02:01|ORF|2017-07-18|1996-01-19|human|HLA-B|I|MHC|2711 bp| genomic-DNA|\nATGGCGCCCCGAACCCTCCTCCTGCTGCTCTC\n"
  let records = @src.imgt_parse(content)
  assert_eq(records.length(), 2)
  assert_eq(records[0].id(), "A*01:01:01:01")
  assert_eq(records[0].sequence(), "ATGGCCGTCATGGCGCCCCGAACCCTCCTCCTC")
  assert_eq(records[0].header().gene(), "HLA-A")
  assert_eq(records[1].id(), "B*07:02:01")
  assert_eq(records[1].sequence(), "ATGGCGCCCCGAACCCTCCTCCTGCTGCTCTC")
  assert_eq(records[1].header().gene(), "HLA-B")
}

///|
test "imgt_parse_single_record" {
  let content = ">HLA00001|A*01:01:01:01|||||human|HLA-A||MHC|| genomic-DNA|\nACGTACGT\n"
  let records = @src.imgt_parse(content)
  assert_eq(records.length(), 1)
  assert_eq(records[0].id(), "A*01:01:01:01")
  assert_eq(records[0].sequence(), "ACGTACGT")
  assert_eq(records[0].header().species(), "human")
}

///|
test "imgt_parse_empty_content" {
  let records = @src.imgt_parse("")
  assert_eq(records.length(), 0)
}

///|
test "imgt_parse_blank_lines_skipped" {
  let content = "\n\n>HLA00001|A*01:01:01:01|||||human|HLA-A||MHC|||\nACGT\n\n"
  let records = @src.imgt_parse(content)
  assert_eq(records.length(), 1)
  assert_eq(records[0].sequence(), "ACGT")
}

///|
test "imgt_parse_multiline_sequence" {
  // Sequence lines are concatenated until the next header.
  let content = ">HLA00001|sid1||||||||||||\nACGT\nTTTT\nGGGG\n"
  let records = @src.imgt_parse(content)
  assert_eq(records.length(), 1)
  assert_eq(records[0].sequence(), "ACGTTTTTGGGG")
}

///|
test "imgt_parse_strips_whitespace_in_sequence" {
  // Spaces and tabs within sequence lines are removed.
  let content = ">HLA00001|sid1||||||||||||\nAC GT\tAC GT\n"
  let records = @src.imgt_parse(content)
  assert_eq(records.length(), 1)
  assert_eq(records[0].sequence(), "ACGTACGT")
}

///|
test "imgt_parse_no_trailing_newline" {
  let content = ">HLA00001|sid1||||||||||||\nACGT"
  let records = @src.imgt_parse(content)
  assert_eq(records.length(), 1)
  assert_eq(records[0].sequence(), "ACGT")
}

// ============================================================================
// 6. imgt_header_to_string - round-trip with parse
// ============================================================================

///|
test "imgt_header_to_string_basic" {
  let h = @src.ImgtHeader::new(
    "HLA00001",
    "A*01:01:01:01",
    "",
    "",
    "",
    "",
    "",
    "",
    "",
    "",
    "",
    "",
    "",
  )
  assert_eq(@src.imgt_header_to_string(h), ">HLA00001|A*01:01:01:01|||||||||||")
}

///|
test "imgt_header_to_string_round_trip_with_parse" {
  let h = @src.ImgtHeader::new(
    "HLA00001",
    "A*01:01:01:01",
    "A*01:01:01:01",
    "ORF",
    "2015-03-31",
    "1996-01-19",
    "human",
    "HLA-A",
    "I",
    "MHC",
    "2703 bp",
    " genomic-DNA",
    "",
  )
  let s = @src.imgt_header_to_string(h)
  match @src.imgt_parse_header(s) {
    Some(h2) => {
      assert_eq(h2.accession(), h.accession())
      assert_eq(h2.seq_id(), h.seq_id())
      assert_eq(h2.official_name(), h.official_name())
      assert_eq(h2.region(), h.region())
      assert_eq(h2.date_update(), h.date_update())
      assert_eq(h2.date_creation(), h.date_creation())
      assert_eq(h2.species(), h.species())
      assert_eq(h2.gene(), h.gene())
      assert_eq(h2.imgt_group(), h.imgt_group())
      assert_eq(h2.locus(), h.locus())
      assert_eq(h2.seq_length(), h.seq_length())
      assert_eq(h2.seq_type(), h.seq_type())
      assert_eq(h2.comments(), h.comments())
    }
    None => assert_true(false)
  }
}

// ============================================================================
// 7. imgt_record_to_string - includes header and wrapped sequence
// ============================================================================

///|
test "imgt_record_to_string_includes_header_and_sequence" {
  let header = @src.ImgtHeader::new(
    "HLA00001",
    "A*01:01:01:01",
    "",
    "",
    "",
    "",
    "",
    "",
    "",
    "",
    "",
    "",
    "",
  )
  let rec = @src.ImgtRecord::new(header, "ATGGCCGTC")
  let s = @src.imgt_record_to_string(rec)
  // The header line starts with '>'.
  assert_true(s.contains(">HLA00001|A*01:01:01:01|"))
  // The sequence is on the next line.
  assert_true(s.contains("ATGGCCGTC"))
}

///|
test "imgt_record_to_string_wraps_long_sequence" {
  // A 120-character sequence should be wrapped into two 60-character lines.
  let seq = "ATGGCCGTCATGGCGCCCCGAACCCTCCTCCTGCTGCTCTCTGGGGCCCCGGGGCCCCGGGGCCCCGGATGGCCGTCATGGCGCCCCGAACCCTCCTCCTGCTGCTCTCTGGGGCCCCGGGGCCCCGGGGCCCCGG"
  let header = @src.ImgtHeader::new("A", "B", "", "", "", "", "", "", "", "", "", "", "")
  let rec = @src.ImgtRecord::new(header, seq)
  let s = @src.imgt_record_to_string(rec)
  // The first line is the header (starts with '>'); the next two lines are
  // each 60 characters of sequence, then a trailing newline.
  let lines = s.split("\n").to_array()
  assert_true(lines[0].starts_with(">"))
  assert_eq(lines[1].length(), 60)
  assert_eq(lines[2].length(), 60)
}

///|
test "imgt_record_to_string_empty_sequence" {
  let header = @src.ImgtHeader::new("A", "B", "", "", "", "", "", "", "", "", "", "", "")
  let rec = @src.ImgtRecord::new(header, "")
  let s = @src.imgt_record_to_string(rec)
  // Header line followed by a newline; no sequence line content.
  assert_true(s.contains(">A|B|"))
}

// ============================================================================
// 8. imgt_write - round-trip with parse
// ============================================================================

///|
test "imgt_write_round_trip_with_parse" {
  let records = @src.imgt_sample_records()
  let text = @src.imgt_write(records)
  let reparsed = @src.imgt_parse(text)
  assert_eq(reparsed.length(), records.length())
  for i in 0..<records.length() {
    assert_eq(reparsed[i].id(), records[i].id())
    assert_eq(reparsed[i].sequence(), records[i].sequence())
    assert_eq(reparsed[i].header().accession(), records[i].header().accession())
    assert_eq(reparsed[i].header().gene(), records[i].header().gene())
    assert_eq(reparsed[i].header().species(), records[i].header().species())
  }
}

///|
test "imgt_write_empty_array" {
  assert_eq(@src.imgt_write([]), "")
}

// ============================================================================
// 9. imgt_n_records
// ============================================================================

///|
test "imgt_n_records_sample" {
  let records = @src.imgt_sample_records()
  assert_eq(@src.imgt_n_records(records), 3)
}

///|
test "imgt_n_records_empty" {
  assert_eq(@src.imgt_n_records([]), 0)
}

///|
test "imgt_n_records_single" {
  let h = @src.ImgtHeader::new("A", "B", "", "", "", "", "", "", "", "", "", "", "")
  let rec = @src.ImgtRecord::new(h, "ACGT")
  assert_eq(@src.imgt_n_records([rec]), 1)
}

// ============================================================================
// 10. imgt_get_by_accession - found and not found
// ============================================================================

///|
test "imgt_get_by_accession_found" {
  let records = @src.imgt_sample_records()
  match @src.imgt_get_by_accession(records, "HLA00001") {
    Some(rec) => {
      assert_eq(rec.id(), "A*01:01:01:01")
      assert_eq(rec.header().gene(), "HLA-A")
    }
    None => assert_true(false)
  }
}

///|
test "imgt_get_by_accession_found_third" {
  let records = @src.imgt_sample_records()
  match @src.imgt_get_by_accession(records, "HLA00003") {
    Some(rec) => assert_eq(rec.id(), "B*07:02:01")
    None => assert_true(false)
  }
}

///|
test "imgt_get_by_accession_not_found" {
  let records = @src.imgt_sample_records()
  match @src.imgt_get_by_accession(records, "MISSING") {
    Some(_) => assert_true(false)
    None => assert_true(true)
  }
}

///|
test "imgt_get_by_accession_empty_array" {
  match @src.imgt_get_by_accession([], "anything") {
    Some(_) => assert_true(false)
    None => assert_true(true)
  }
}

// ============================================================================
// 11. imgt_get_by_species - filter results
// ============================================================================

///|
test "imgt_get_by_species_human" {
  let records = @src.imgt_sample_records()
  let result = @src.imgt_get_by_species(records, "human")
  assert_eq(result.length(), 3)
}

///|
test "imgt_get_by_species_no_match" {
  let records = @src.imgt_sample_records()
  let result = @src.imgt_get_by_species(records, "mouse")
  assert_eq(result.length(), 0)
}

///|
test "imgt_get_by_species_empty_array" {
  let result = @src.imgt_get_by_species([], "human")
  assert_eq(result.length(), 0)
}

// ============================================================================
// 12. imgt_get_by_gene - filter results
// ============================================================================

///|
test "imgt_get_by_gene_hla_a" {
  let records = @src.imgt_sample_records()
  let result = @src.imgt_get_by_gene(records, "HLA-A")
  assert_eq(result.length(), 2)
  assert_eq(result[0].id(), "A*01:01:01:01")
  assert_eq(result[1].id(), "A*02:01:01:01")
}

///|
test "imgt_get_by_gene_hla_b" {
  let records = @src.imgt_sample_records()
  let result = @src.imgt_get_by_gene(records, "HLA-B")
  assert_eq(result.length(), 1)
  assert_eq(result[0].id(), "B*07:02:01")
}

///|
test "imgt_get_by_gene_no_match" {
  let records = @src.imgt_sample_records()
  let result = @src.imgt_get_by_gene(records, "HLA-C")
  assert_eq(result.length(), 0)
}

// ============================================================================
// 13. imgt_filter_by_locus
// ============================================================================

///|
test "imgt_filter_by_locus_mhc" {
  let records = @src.imgt_sample_records()
  let result = @src.imgt_filter_by_locus(records, "MHC")
  assert_eq(result.length(), 3)
}

///|
test "imgt_filter_by_locus_no_match" {
  let records = @src.imgt_sample_records()
  let result = @src.imgt_filter_by_locus(records, "NON-MHC")
  assert_eq(result.length(), 0)
}

// ============================================================================
// 14. imgt_unique_species and imgt_unique_genes
// ============================================================================

///|
test "imgt_unique_species_sample" {
  let records = @src.imgt_sample_records()
  let species = @src.imgt_unique_species(records)
  assert_eq(species.length(), 1)
  assert_eq(species[0], "human")
}

///|
test "imgt_unique_species_empty" {
  let species = @src.imgt_unique_species([])
  assert_eq(species.length(), 0)
}

///|
test "imgt_unique_genes_sample" {
  let records = @src.imgt_sample_records()
  let genes = @src.imgt_unique_genes(records)
  assert_eq(genes.length(), 2)
  assert_eq(genes[0], "HLA-A")
  assert_eq(genes[1], "HLA-B")
}

///|
test "imgt_unique_genes_empty" {
  let genes = @src.imgt_unique_genes([])
  assert_eq(genes.length(), 0)
}

///|
test "imgt_unique_species_preserves_first_appearance_order" {
  // Create records with multiple species to verify order.
  let h1 = @src.ImgtHeader::new("A", "s1", "", "", "", "", "zebra", "", "", "", "", "", "")
  let h2 = @src.ImgtHeader::new("B", "s2", "", "", "", "", "mouse", "", "", "", "", "", "")
  let h3 = @src.ImgtHeader::new("C", "s3", "", "", "", "", "zebra", "", "", "", "", "", "")
  let recs = [
    @src.ImgtRecord::new(h1, "ACGT"),
    @src.ImgtRecord::new(h2, "TTTT"),
    @src.ImgtRecord::new(h3, "GGGG"),
  ]
  let species = @src.imgt_unique_species(recs)
  assert_eq(species.length(), 2)
  assert_eq(species[0], "zebra")
  assert_eq(species[1], "mouse")
}

// ============================================================================
// 15. imgt_summary - non-empty string
// ============================================================================

///|
test "imgt_summary_sample" {
  let records = @src.imgt_sample_records()
  let summary = @src.imgt_summary(records)
  assert_true(summary.length() > 0)
  assert_true(summary.contains("IMGT records: 3"))
  assert_true(summary.contains("species: 1"))
  assert_true(summary.contains("genes: 2"))
  assert_true(summary.contains("human"))
  assert_true(summary.contains("HLA-A"))
  assert_true(summary.contains("HLA-B"))
}

///|
test "imgt_summary_empty" {
  let summary = @src.imgt_summary([])
  assert_true(summary.contains("IMGT records: 0"))
  assert_true(summary.contains("species: 0"))
  assert_true(summary.contains("genes: 0"))
}

// ============================================================================
// 16. SeqRecord conversion - to and from
// ============================================================================

///|
test "imgt_to_seq_records" {
  let records = @src.imgt_sample_records()
  let seq_records = @src.imgt_to_seq_records(records)
  assert_eq(seq_records.length(), 3)
  // id and name come from the header's seq_id.
  assert_eq(seq_records[0].id, "A*01:01:01:01")
  assert_eq(seq_records[0].name, "A*01:01:01:01")
  // description is the full pipe-separated header.
  assert_eq(seq_records[0].description, records[0].header().to_string())
  assert_eq(seq_records[0].seq.to_string(), records[0].sequence())
  assert_eq(seq_records[1].id, "A*02:01:01:01")
  assert_eq(seq_records[2].id, "B*07:02:01")
}

///|
test "imgt_from_seq_records" {
  let seq_records = [
    @src.SeqRecord::new(@src.Seq::new("ACGTACGT"), id="rec1", description="desc1"),
    @src.SeqRecord::new(@src.Seq::new("TTTTGGGG"), id="rec2", description="desc2"),
  ]
  let imgt_records = @src.imgt_from_seq_records(seq_records)
  assert_eq(imgt_records.length(), 2)
  // accession and seq_id are both set to the record id.
  assert_eq(imgt_records[0].header().accession(), "rec1")
  assert_eq(imgt_records[0].header().seq_id(), "rec1")
  // seq_length is derived from the sequence length.
  assert_eq(imgt_records[0].header().seq_length(), "8 bp")
  // comments holds the original description.
  assert_eq(imgt_records[0].header().comments(), "desc1")
  assert_eq(imgt_records[0].sequence(), "ACGTACGT")
  assert_eq(imgt_records[1].header().seq_length(), "8 bp")
  assert_eq(imgt_records[1].sequence(), "TTTTGGGG")
}

///|
test "imgt_seq_records_round_trip" {
  let records = @src.imgt_sample_records()
  let seq_records = @src.imgt_to_seq_records(records)
  let back = @src.imgt_from_seq_records(seq_records)
  assert_eq(back.length(), records.length())
  // The sequences survive the round trip.
  for i in 0..<records.length() {
    assert_eq(back[i].sequence(), records[i].sequence())
  }
}

// ============================================================================
// 17. Sample records validation - 3 records, 2 HLA-A + 1 HLA-B
// ============================================================================

///|
test "imgt_sample_records_count" {
  let records = @src.imgt_sample_records()
  assert_eq(records.length(), 3)
}

///|
test "imgt_sample_records_hla_a_count" {
  let records = @src.imgt_sample_records()
  let hla_a = @src.imgt_get_by_gene(records, "HLA-A")
  assert_eq(hla_a.length(), 2)
}

///|
test "imgt_sample_records_hla_b_count" {
  let records = @src.imgt_sample_records()
  let hla_b = @src.imgt_get_by_gene(records, "HLA-B")
  assert_eq(hla_b.length(), 1)
}

///|
test "imgt_sample_records_first_record" {
  let records = @src.imgt_sample_records()
  let r = records[0]
  assert_eq(r.header().accession(), "HLA00001")
  assert_eq(r.id(), "A*01:01:01:01")
  assert_eq(r.header().official_name(), "A*01:01:01:01")
  assert_eq(r.header().region(), "ORF")
  assert_eq(r.header().date_update(), "2015-03-31")
  assert_eq(r.header().date_creation(), "1996-01-19")
  assert_eq(r.header().species(), "human")
  assert_eq(r.header().gene(), "HLA-A")
  assert_eq(r.header().imgt_group(), "I")
  assert_eq(r.header().locus(), "MHC")
  assert_eq(r.header().seq_length(), "2703 bp")
  assert_eq(r.header().seq_type(), " genomic-DNA")
  assert_eq(r.sequence(), "ATGGCCGTCATGGCGCCCCGAACCCTCCTCCTC")
}

///|
test "imgt_sample_records_second_record" {
  let records = @src.imgt_sample_records()
  let r = records[1]
  assert_eq(r.header().accession(), "HLA00002")
  assert_eq(r.id(), "A*02:01:01:01")
  assert_eq(r.header().gene(), "HLA-A")
  assert_eq(r.sequence(), "ATGGCCGTCATGGCGCCCCGAACCCTCCTCCTG")
}

///|
test "imgt_sample_records_third_record" {
  let records = @src.imgt_sample_records()
  let r = records[2]
  assert_eq(r.header().accession(), "HLA00003")
  assert_eq(r.id(), "B*07:02:01")
  assert_eq(r.header().gene(), "HLA-B")
  assert_eq(r.sequence(), "ATGGCGCCCCGAACCCTCCTCCTGCTGCTCTC")
}

///|
test "imgt_sample_records_all_human" {
  let records = @src.imgt_sample_records()
  for i in 0..<records.length() {
    assert_eq(records[i].header().species(), "human")
  }
}

///|
test "imgt_sample_records_all_mhc_locus" {
  let records = @src.imgt_sample_records()
  for i in 0..<records.length() {
    assert_eq(records[i].header().locus(), "MHC")
  }
}

// ============================================================================
// 18. Edge cases - empty array, single record, header with missing fields
// ============================================================================

///|
test "edge_case_empty_array_n_records" {
  assert_eq(@src.imgt_n_records([]), 0)
}

///|
test "edge_case_empty_array_write" {
  assert_eq(@src.imgt_write([]), "")
}

///|
test "edge_case_single_record" {
  let h = @src.ImgtHeader::new("A", "B", "", "", "", "", "", "", "", "", "", "", "")
  let rec = @src.ImgtRecord::new(h, "ACGT")
  let text = @src.imgt_write([rec])
  let reparsed = @src.imgt_parse(text)
  assert_eq(reparsed.length(), 1)
  assert_eq(reparsed[0].id(), "B")
  assert_eq(reparsed[0].sequence(), "ACGT")
}

///|
test "edge_case_header_with_missing_fields" {
  // Only 3 fields provided; the rest default to empty.
  match @src.imgt_parse_header(">ACC|SID|ORF") {
    Some(h) => {
      assert_eq(h.accession(), "ACC")
      assert_eq(h.seq_id(), "SID")
      assert_eq(h.official_name(), "ORF")
      assert_eq(h.region(), "")
      assert_eq(h.date_update(), "")
      assert_eq(h.species(), "")
      assert_eq(h.gene(), "")
      assert_eq(h.locus(), "")
      assert_eq(h.comments(), "")
    }
    None => assert_true(false)
  }
}

///|
test "edge_case_parse_header_only_no_sequence" {
  // A header line with no sequence following it still produces a record with
  // an empty sequence.
  let records = @src.imgt_parse(">HLA00001|sid1||||||||||||\n")
  assert_eq(records.length(), 1)
  assert_eq(records[0].id(), "sid1")
  assert_eq(records[0].sequence(), "")
}

///|
test "edge_case_to_seq_records_empty" {
  let seq_records = @src.imgt_to_seq_records([])
  assert_eq(seq_records.length(), 0)
}

///|
test "edge_case_from_seq_records_empty" {
  let imgt_records = @src.imgt_from_seq_records([])
  assert_eq(imgt_records.length(), 0)
}
