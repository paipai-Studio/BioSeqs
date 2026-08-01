///|
/// Tests for Bio.UniProt.GOA (Gene Association Format / GAF) parser.

// ---------------------------------------------------------------------------
// GafAspect enum
// ---------------------------------------------------------------------------

test "goa_aspect_from_string_F" {
  let a = @src.GafAspect::from_string("F")
  assert_eq(@src.GafAspect::to_string(a), "F")
  assert_eq(@src.GafAspect::description(a), "Molecular Function")
}

test "goa_aspect_from_string_P" {
  let a = @src.GafAspect::from_string("P")
  assert_eq(@src.GafAspect::to_string(a), "P")
  assert_eq(@src.GafAspect::description(a), "Biological Process")
}

test "goa_aspect_from_string_C" {
  let a = @src.GafAspect::from_string("C")
  assert_eq(@src.GafAspect::to_string(a), "C")
  assert_eq(@src.GafAspect::description(a), "Cellular Component")
}

test "goa_aspect_from_string_unknown" {
  // Unknown inputs (e.g. "X", "", "foo") map to Unknown.
  let a1 = @src.GafAspect::from_string("X")
  let a2 = @src.GafAspect::from_string("")
  let a3 = @src.GafAspect::from_string("foo")
  assert_eq(@src.GafAspect::to_string(a1), "?")
  assert_eq(@src.GafAspect::to_string(a2), "?")
  assert_eq(@src.GafAspect::to_string(a3), "?")
  assert_eq(@src.GafAspect::description(a1), "Unknown")
  assert_eq(@src.GafAspect::description(a2), "Unknown")
  assert_eq(@src.GafAspect::description(a3), "Unknown")
}

test "goa_aspect_to_string_roundtrip" {
  // from_string -> to_string should round-trip for the three valid codes.
  assert_eq(@src.GafAspect::to_string(@src.GafAspect::from_string("F")), "F")
  assert_eq(@src.GafAspect::to_string(@src.GafAspect::from_string("P")), "P")
  assert_eq(@src.GafAspect::to_string(@src.GafAspect::from_string("C")), "C")
}

// ---------------------------------------------------------------------------
// GafRecord creation and accessors
// ---------------------------------------------------------------------------

test "goa_record_creation" {
  let r = @src.GafRecord::new(
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
  )
  assert_eq(r.db(), "UniProtKB")
  assert_eq(r.db_object_id(), "Q12345")
  assert_eq(r.db_object_symbol(), "PROT1")
  assert_eq(r.qualifier(), "enables")
  assert_eq(r.go_id(), "GO:0003674")
  assert_eq(r.db_reference(), "PMID:12345")
  assert_eq(r.evidence_code(), "IDA")
  assert_eq(r.with_from(), "GO:0005515")
  assert_eq(r.aspect(), "F")
  assert_eq(r.db_object_name(), "Protein 1")
  assert_eq(r.db_object_synonyms(), "P1|PROT-1")
  assert_eq(r.db_object_type(), "protein")
  assert_eq(r.taxon(), "taxon:9606|taxon:9606")
  assert_eq(r.date(), "20210115")
  assert_eq(r.assigned_by(), "UniProt")
  assert_eq(r.annotation_extension(), "")
  assert_eq(r.gene_product_form_id(), "")
}

test "goa_record_creation_all_fields_populated" {
  // Record with every field populated, including the optional extension/form id.
  let r = @src.GafRecord::new(
    db="UniProtKB",
    db_object_id="P12345",
    db_object_symbol="GENE_X",
    qualifier="NOT|enables",
    go_id="GO:0005515",
    db_reference="PMID:111|PMID:222",
    evidence_code="IPI",
    with_from="UniProtKB:Q99999",
    aspect="F",
    db_object_name="Gene X protein",
    db_object_synonyms="GX|XGENE",
    db_object_type="protein",
    taxon="taxon:9606|taxon:9606",
    date="20230101",
    assigned_by="UniProt",
    annotation_extension="part_of(CL:0000001)",
    gene_product_form_id="UniProtKB:P12345-1",
  )
  assert_eq(r.annotation_extension(), "part_of(CL:0000001)")
  assert_eq(r.gene_product_form_id(), "UniProtKB:P12345-1")
  assert_eq(r.qualifier(), "NOT|enables")
  assert_eq(r.db_reference(), "PMID:111|PMID:222")
  assert_eq(r.db_object_synonyms(), "GX|XGENE")
}

// ---------------------------------------------------------------------------
// GafRecord utility methods
// ---------------------------------------------------------------------------

test "goa_record_go_id_short" {
  let r = @src.GafRecord::new(
    db="UniProtKB",
    db_object_id="Q12345",
    db_object_symbol="PROT1",
    qualifier="enables",
    go_id="GO:0003674",
    db_reference="PMID:12345",
    evidence_code="IDA",
    with_from="",
    aspect="F",
    db_object_name="",
    db_object_synonyms="",
    db_object_type="protein",
    taxon="taxon:9606",
    date="20210115",
    assigned_by="UniProt",
    annotation_extension="",
    gene_product_form_id="",
  )
  assert_eq(r.go_id_short(), "0003674")
}

test "goa_record_go_id_short_no_prefix" {
  // A GO ID without "GO:" prefix should be returned unchanged.
  let r = @src.GafRecord::new(
    db="UniProtKB",
    db_object_id="Q12345",
    db_object_symbol="PROT1",
    qualifier="enables",
    go_id="0003674",
    db_reference="",
    evidence_code="IDA",
    with_from="",
    aspect="F",
    db_object_name="",
    db_object_synonyms="",
    db_object_type="protein",
    taxon="taxon:9606",
    date="20210115",
    assigned_by="UniProt",
    annotation_extension="",
    gene_product_form_id="",
  )
  assert_eq(r.go_id_short(), "0003674")
}

test "goa_record_qualifiers_split" {
  let r = @src.GafRecord::new(
    db="UniProtKB",
    db_object_id="Q12345",
    db_object_symbol="PROT1",
    qualifier="NOT|enables",
    db_reference="",
    go_id="GO:0003674",
    evidence_code="IDA",
    with_from="",
    aspect="F",
    db_object_name="",
    db_object_synonyms="",
    db_object_type="protein",
    taxon="taxon:9606",
    date="20210115",
    assigned_by="UniProt",
    annotation_extension="",
    gene_product_form_id="",
  )
  let qs = r.qualifiers()
  assert_eq(qs.length(), 2)
  assert_eq(qs[0], "NOT")
  assert_eq(qs[1], "enables")
}

test "goa_record_qualifiers_single" {
  let r = @src.GafRecord::new(
    db="UniProtKB",
    db_object_id="Q12345",
    db_object_symbol="PROT1",
    qualifier="enables",
    db_reference="",
    go_id="GO:0003674",
    evidence_code="IDA",
    with_from="",
    aspect="F",
    db_object_name="",
    db_object_synonyms="",
    db_object_type="protein",
    taxon="taxon:9606",
    date="20210115",
    assigned_by="UniProt",
    annotation_extension="",
    gene_product_form_id="",
  )
  let qs = r.qualifiers()
  assert_eq(qs.length(), 1)
  assert_eq(qs[0], "enables")
}

test "goa_record_qualifiers_empty" {
  let r = @src.GafRecord::new(
    db="UniProtKB",
    db_object_id="Q12345",
    db_object_symbol="PROT1",
    qualifier="",
    db_reference="",
    go_id="GO:0003674",
    evidence_code="IDA",
    with_from="",
    aspect="F",
    db_object_name="",
    db_object_synonyms="",
    db_object_type="protein",
    taxon="taxon:9606",
    date="20210115",
    assigned_by="UniProt",
    annotation_extension="",
    gene_product_form_id="",
  )
  assert_eq(r.qualifiers().length(), 0)
}

test "goa_record_synonyms_split" {
  let r = @src.GafRecord::new(
    db="UniProtKB",
    db_object_id="P04637",
    db_object_symbol="TP53",
    qualifier="enables",
    go_id="GO:0003700",
    db_reference="",
    evidence_code="EXP",
    with_from="",
    aspect="F",
    db_object_name="Cellular tumor antigen p53",
    db_object_synonyms="p53|LFS1",
    db_object_type="protein",
    taxon="taxon:9606",
    date="20210220",
    assigned_by="UniProt",
    annotation_extension="",
    gene_product_form_id="",
  )
  let syns = r.synonyms()
  assert_eq(syns.length(), 2)
  assert_eq(syns[0], "p53")
  assert_eq(syns[1], "LFS1")
}

test "goa_record_references_split" {
  let r = @src.GafRecord::new(
    db="UniProtKB",
    db_object_id="Q12345",
    db_object_symbol="PROT1",
    qualifier="enables",
    go_id="GO:0003674",
    db_reference="PMID:111|PMID:222|GO_REF:000001",
    evidence_code="IDA",
    with_from="",
    aspect="F",
    db_object_name="",
    db_object_synonyms="",
    db_object_type="protein",
    taxon="taxon:9606",
    date="20210115",
    assigned_by="UniProt",
    annotation_extension="",
    gene_product_form_id="",
  )
  let refs = r.references()
  assert_eq(refs.length(), 3)
  assert_eq(refs[0], "PMID:111")
  assert_eq(refs[1], "PMID:222")
  assert_eq(refs[2], "GO_REF:000001")
}

test "goa_record_taxon_ids_with_prefix" {
  // "taxon:9606|taxon:9606" yields ["9606", "9606"] (duplicates preserved).
  let r = @src.GafRecord::new(
    db="UniProtKB",
    db_object_id="Q12345",
    db_object_symbol="PROT1",
    qualifier="enables",
    go_id="GO:0003674",
    db_reference="",
    evidence_code="IDA",
    with_from="",
    aspect="F",
    db_object_name="",
    db_object_synonyms="",
    db_object_type="protein",
    taxon="taxon:9606|taxon:9606",
    date="20210115",
    assigned_by="UniProt",
    annotation_extension="",
    gene_product_form_id="",
  )
  let ids = r.taxon_ids()
  assert_eq(ids.length(), 2)
  assert_eq(ids[0], "9606")
  assert_eq(ids[1], "9606")
}

test "goa_record_taxon_ids_bare_numeric" {
  // Numeric taxon IDs without the "taxon:" prefix should be returned as-is.
  let r = @src.GafRecord::new(
    db="UniProtKB",
    db_object_id="Q12345",
    db_object_symbol="PROT1",
    qualifier="enables",
    go_id="GO:0003674",
    db_reference="",
    evidence_code="IDA",
    with_from="",
    aspect="F",
    db_object_name="",
    db_object_synonyms="",
    db_object_type="protein",
    taxon="9606",
    date="20210115",
    assigned_by="UniProt",
    annotation_extension="",
    gene_product_form_id="",
  )
  let ids = r.taxon_ids()
  assert_eq(ids.length(), 1)
  assert_eq(ids[0], "9606")
}

test "goa_record_taxon_ids_multiple_distinct" {
  // Two distinct taxon IDs (interactor case).
  let r = @src.GafRecord::new(
    db="UniProtKB",
    db_object_id="Q12345",
    db_object_symbol="PROT1",
    qualifier="enables",
    go_id="GO:0003674",
    db_reference="",
    evidence_code="IDA",
    with_from="",
    aspect="F",
    db_object_name="",
    db_object_synonyms="",
    db_object_type="protein",
    taxon="taxon:9606|taxon:4932",
    date="20210115",
    assigned_by="UniProt",
    annotation_extension="",
    gene_product_form_id="",
  )
  let ids = r.taxon_ids()
  assert_eq(ids.length(), 2)
  assert_eq(ids[0], "9606")
  assert_eq(ids[1], "4932")
}

test "goa_record_to_string" {
  let r = @src.GafRecord::new(
    db="UniProtKB",
    db_object_id="Q12345",
    db_object_symbol="PROT1",
    qualifier="enables",
    go_id="GO:0003674",
    db_reference="",
    evidence_code="IDA",
    with_from="",
    aspect="F",
    db_object_name="",
    db_object_synonyms="",
    db_object_type="protein",
    taxon="taxon:9606|taxon:9606",
    date="20210115",
    assigned_by="UniProt",
    annotation_extension="",
    gene_product_form_id="",
  )
  let s = r.to_string()
  assert_true(s.contains("UniProtKB"))
  assert_true(s.contains("Q12345"))
  assert_true(s.contains("PROT1"))
  assert_true(s.contains("GO:0003674"))
  assert_true(s.contains("IDA"))
  assert_true(s.contains("F"))
}

// ---------------------------------------------------------------------------
// Header line parsing
// ---------------------------------------------------------------------------

test "goa_parse_header_line_valid" {
  let h = @src.goa_parse_header_line("!gaf-version: 2.2")
  assert_true(h is Some(_))
  match h {
    Some((k, v)) => {
      assert_eq(k, "gaf-version")
      assert_eq(v, "2.2")
    }
    None => assert_true(false)
  }
}

test "goa_parse_header_line_generated_by" {
  let h = @src.goa_parse_header_line("!generated-by: UniProt")
  assert_true(h is Some(_))
  match h {
    Some((k, v)) => {
      assert_eq(k, "generated-by")
      assert_eq(v, "UniProt")
    }
    None => assert_true(false)
  }
}

test "goa_parse_header_line_with_extra_spaces" {
  let h = @src.goa_parse_header_line("!  gaf-version  :   2.2  ")
  assert_true(h is Some(_))
  match h {
    Some((k, v)) => {
      assert_eq(k, "gaf-version")
      assert_eq(v, "2.2")
    }
    None => assert_true(false)
  }
}

test "goa_parse_header_line_not_header" {
  // A line that does not start with "!" is not a header.
  let h = @src.goa_parse_header_line("UniProtKB\tQ12345")
  assert_true(h is None)
}

test "goa_parse_header_line_no_colon" {
  // A header line with no colon is invalid.
  let h = @src.goa_parse_header_line("!this-has-no-colon")
  assert_true(h is None)
}

test "goa_parse_header_line_empty" {
  let h = @src.goa_parse_header_line("")
  assert_true(h is None)
}

test "goa_parse_header_line_empty_value" {
  // A header line with a colon but empty value is valid (value is "").
  let h = @src.goa_parse_header_line("!gaf-version:")
  assert_true(h is Some(_))
  match h {
    Some((k, v)) => {
      assert_eq(k, "gaf-version")
      assert_eq(v, "")
    }
    None => assert_true(false)
  }
}

// ---------------------------------------------------------------------------
// Data line parsing
// ---------------------------------------------------------------------------

test "goa_parse_line_full_17_columns" {
  let line = "UniProtKB\tQ12345\tPROT1\tenables\tGO:0003674\tPMID:12345\tIDA\tGO:0005515\tF\tProtein 1\tP1|PROT-1\tprotein\ttaxon:9606|taxon:9606\t20210115\tUniProt\t\t"
  let r = @src.goa_parse_line(line)
  assert_true(r is Some(_))
  match r {
    Some(rec) => {
      assert_eq(rec.db(), "UniProtKB")
      assert_eq(rec.db_object_id(), "Q12345")
      assert_eq(rec.db_object_symbol(), "PROT1")
      assert_eq(rec.qualifier(), "enables")
      assert_eq(rec.go_id(), "GO:0003674")
      assert_eq(rec.db_reference(), "PMID:12345")
      assert_eq(rec.evidence_code(), "IDA")
      assert_eq(rec.with_from(), "GO:0005515")
      assert_eq(rec.aspect(), "F")
      assert_eq(rec.db_object_name(), "Protein 1")
      assert_eq(rec.db_object_synonyms(), "P1|PROT-1")
      assert_eq(rec.db_object_type(), "protein")
      assert_eq(rec.taxon(), "taxon:9606|taxon:9606")
      assert_eq(rec.date(), "20210115")
      assert_eq(rec.assigned_by(), "UniProt")
      assert_eq(rec.annotation_extension(), "")
      assert_eq(rec.gene_product_form_id(), "")
    }
    None => assert_true(false)
  }
}

test "goa_parse_line_short_padded" {
  // A line with fewer than 17 columns should be right-padded with empty strings.
  let line = "UniProtKB\tQ12345\tPROT1\tenables\tGO:0003674"
  let r = @src.goa_parse_line(line)
  assert_true(r is Some(_))
  match r {
    Some(rec) => {
      assert_eq(rec.db(), "UniProtKB")
      assert_eq(rec.db_object_id(), "Q12345")
      assert_eq(rec.db_object_symbol(), "PROT1")
      assert_eq(rec.qualifier(), "enables")
      assert_eq(rec.go_id(), "GO:0003674")
      assert_eq(rec.db_reference(), "")
      assert_eq(rec.evidence_code(), "")
      assert_eq(rec.aspect(), "")
      assert_eq(rec.taxon(), "")
      assert_eq(rec.assigned_by(), "")
    }
    None => assert_true(false)
  }
}

test "goa_parse_line_empty" {
  let r = @src.goa_parse_line("")
  assert_true(r is None)
}

test "goa_parse_line_whitespace_only" {
  let r = @src.goa_parse_line("     ")
  assert_true(r is None)
}

test "goa_parse_line_header_returns_none" {
  // Header lines (starting with "!") should not be parsed as data records.
  let r = @src.goa_parse_line("!gaf-version: 2.2")
  assert_true(r is None)
}

test "goa_parse_line_comment_returns_none" {
  // Any line starting with "!" is treated as a header/comment, not data.
  let r = @src.goa_parse_line("!some comment line")
  assert_true(r is None)
}

// ---------------------------------------------------------------------------
// Full content parsing
// ---------------------------------------------------------------------------

test "goa_parse_with_header_and_data" {
  let content = "!gaf-version: 2.2\n" +
    "!generated-by: UniProt\n" +
    "UniProtKB\tQ12345\tPROT1\tenables\tGO:0003674\tPMID:12345\tIDA\tGO:0005515\tF\tProtein 1\tP1|PROT-1\tprotein\ttaxon:9606|taxon:9606\t20210115\tUniProt\t\t\n" +
    "UniProtKB\tP04637\tTP53\tenables\tGO:0003700\tPMID:67890\tEXP\tGO:0005634\tF\tCellular tumor antigen p53\tp53|LFS1\tprotein\ttaxon:9606|taxon:9606\t20210220\tUniProt\t\t\n"
  let db = @src.goa_parse(content)
  assert_eq(db.n_records(), 2)
  assert_eq(db.version(), "2.2")
  assert_eq(db.created_by(), "UniProt")
  let recs = db.records()
  assert_eq(recs[0].db_object_symbol(), "PROT1")
  assert_eq(recs[1].db_object_symbol(), "TP53")
}

test "goa_parse_empty_input" {
  let db = @src.goa_parse("")
  assert_eq(db.n_records(), 0)
  assert_eq(db.version(), "")
  assert_eq(db.created_by(), "")
}

test "goa_parse_only_headers" {
  let content = "!gaf-version: 2.2\n" + "!generated-by: SGD\n"
  let db = @src.goa_parse(content)
  assert_eq(db.n_records(), 0)
  assert_eq(db.version(), "2.2")
  assert_eq(db.created_by(), "SGD")
}

test "goa_parse_skips_blank_lines" {
  let content = "!gaf-version: 2.2\n" +
    "\n" +
    "UniProtKB\tQ12345\tPROT1\tenables\tGO:0003674\tPMID:12345\tIDA\t\tF\tProtein 1\t\tprotein\ttaxon:9606\t20210115\tUniProt\t\t\n" +
    "\n" +
    "UniProtKB\tP04637\tTP53\tenables\tGO:0003700\tPMID:67890\tEXP\t\tF\tCellular tumor antigen p53\t\tprotein\ttaxon:9606\t20210220\tUniProt\t\t\n"
  let db = @src.goa_parse(content)
  assert_eq(db.n_records(), 2)
}

test "goa_parse_ignores_unknown_header_keys" {
  // Header keys other than gaf-version / generated-by should not overwrite metadata.
  let content = "!gaf-version: 2.2\n" +
    "!generated-by: UniProt\n" +
    "!some-other-key: whatever\n" +
    "UniProtKB\tQ12345\tPROT1\tenables\tGO:0003674\tPMID:12345\tIDA\t\tF\tProtein 1\t\tprotein\ttaxon:9606\t20210115\tUniProt\t\t\n"
  let db = @src.goa_parse(content)
  assert_eq(db.n_records(), 1)
  assert_eq(db.version(), "2.2")
  assert_eq(db.created_by(), "UniProt")
}

// ---------------------------------------------------------------------------
// GoaDatabase properties
// ---------------------------------------------------------------------------

test "goa_database_new_empty" {
  let db = @src.GoaDatabase::new()
  assert_eq(db.n_records(), 0)
  assert_eq(db.records().length(), 0)
  assert_eq(db.version(), "")
  assert_eq(db.created_by(), "")
}

test "goa_database_records_accessor" {
  let db = @src.goa_sample_database()
  let recs = db.records()
  assert_eq(recs.length(), db.n_records())
}

// ---------------------------------------------------------------------------
// Filter functions
// ---------------------------------------------------------------------------

test "goa_filter_by_go_id" {
  let db = @src.goa_sample_database()
  // GO:0003674 (molecular_function) appears once in the sample database.
  let hits = @src.goa_filter_by_go_id(db, "GO:0003674")
  assert_eq(hits.length(), 1)
  assert_eq(hits[0].db_object_symbol(), "PROT1")
  assert_eq(hits[0].go_id(), "GO:0003674")
}

test "goa_filter_by_go_id_no_match" {
  let db = @src.goa_sample_database()
  let hits = @src.goa_filter_by_go_id(db, "GO:9999999")
  assert_eq(hits.length(), 0)
}

test "goa_filter_by_aspect_molecular_function" {
  let db = @src.goa_sample_database()
  // F records: PROT1, TP53, YFG1, RPL5 (4 total).
  let hits = @src.goa_filter_by_aspect(db, @src.GafAspect::from_string("F"))
  assert_eq(hits.length(), 4)
  for h in hits {
    assert_eq(h.aspect(), "F")
  }
}

test "goa_filter_by_aspect_biological_process" {
  let db = @src.goa_sample_database()
  // P records: PROT1, TP53, Hsp70 (3 total).
  let hits = @src.goa_filter_by_aspect(db, @src.GafAspect::from_string("P"))
  assert_eq(hits.length(), 3)
  for h in hits {
    assert_eq(h.aspect(), "P")
  }
}

test "goa_filter_by_aspect_cellular_component" {
  let db = @src.goa_sample_database()
  // C records: PROT1, YFG1, Hsp70 (3 total).
  let hits = @src.goa_filter_by_aspect(db, @src.GafAspect::from_string("C"))
  assert_eq(hits.length(), 3)
  for h in hits {
    assert_eq(h.aspect(), "C")
  }
}

test "goa_filter_by_aspect_unknown" {
  let db = @src.goa_sample_database()
  // No Unknown aspect records in the sample database.
  let hits = @src.goa_filter_by_aspect(db, @src.GafAspect::from_string("X"))
  assert_eq(hits.length(), 0)
}

test "goa_filter_by_evidence_ida" {
  let db = @src.goa_sample_database()
  // IDA records: PROT1 (x3), YFG1 (C), Hsp70 (C) = 5 total.
  let hits = @src.goa_filter_by_evidence(db, "IDA")
  assert_eq(hits.length(), 5)
  for h in hits {
    assert_eq(h.evidence_code(), "IDA")
  }
}

test "goa_filter_by_evidence_iea" {
  let db = @src.goa_sample_database()
  // IEA records: TP53 (P), Hsp70 (P) = 2 total.
  let hits = @src.goa_filter_by_evidence(db, "IEA")
  assert_eq(hits.length(), 2)
}

test "goa_filter_by_evidence_no_match" {
  let db = @src.goa_sample_database()
  let hits = @src.goa_filter_by_evidence(db, "NONEXISTENT")
  assert_eq(hits.length(), 0)
}

test "goa_filter_by_taxon_human" {
  let db = @src.goa_sample_database()
  // taxon:9606 records: PROT1 (x3), TP53 (x2), RPL5 = 6 total.
  let hits = @src.goa_filter_by_taxon(db, "9606")
  assert_eq(hits.length(), 6)
}

test "goa_filter_by_taxon_yeast" {
  let db = @src.goa_sample_database()
  // taxon:4932 records: YFG1 (x2) = 2 total.
  let hits = @src.goa_filter_by_taxon(db, "4932")
  assert_eq(hits.length(), 2)
}

test "goa_filter_by_taxon_fly" {
  let db = @src.goa_sample_database()
  // taxon:7227 records: Hsp70 (x2) = 2 total.
  let hits = @src.goa_filter_by_taxon(db, "7227")
  assert_eq(hits.length(), 2)
}

test "goa_filter_by_taxon_no_match" {
  let db = @src.goa_sample_database()
  let hits = @src.goa_filter_by_taxon(db, "0000000")
  assert_eq(hits.length(), 0)
}

test "goa_filter_by_db_object_id" {
  let db = @src.goa_sample_database()
  // Q12345 (PROT1) has 3 records (F, P, C).
  let hits = @src.goa_filter_by_db_object_id(db, "Q12345")
  assert_eq(hits.length(), 3)
  for h in hits {
    assert_eq(h.db_object_id(), "Q12345")
    assert_eq(h.db_object_symbol(), "PROT1")
  }
}

test "goa_filter_by_db_object_id_multiple" {
  let db = @src.goa_sample_database()
  // P04637 (TP53) has 2 records (F, P).
  let hits = @src.goa_filter_by_db_object_id(db, "P04637")
  assert_eq(hits.length(), 2)
  for h in hits {
    assert_eq(h.db_object_symbol(), "TP53")
  }
}

test "goa_filter_by_db_object_id_no_match" {
  let db = @src.goa_sample_database()
  let hits = @src.goa_filter_by_db_object_id(db, "ZZZZZZ")
  assert_eq(hits.length(), 0)
}

// ---------------------------------------------------------------------------
// Unique value extraction
// ---------------------------------------------------------------------------

test "goa_unique_go_ids" {
  let db = @src.goa_sample_database()
  // Every record in the sample database has a distinct GO ID.
  let ids = @src.goa_unique_go_ids(db)
  assert_eq(ids.length(), 10)
  // First-seen order: first record's GO ID is GO:0003674.
  assert_eq(ids[0], "GO:0003674")
}

test "goa_unique_go_ids_dedup" {
  // Build a database with duplicate GO IDs to test deduplication.
  let rec = @src.GafRecord::new(
    db="UniProtKB",
    db_object_id="Q12345",
    db_object_symbol="PROT1",
    qualifier="enables",
    go_id="GO:0003674",
    db_reference="",
    evidence_code="IDA",
    with_from="",
    aspect="F",
    db_object_name="",
    db_object_synonyms="",
    db_object_type="protein",
    taxon="taxon:9606",
    date="20210115",
    assigned_by="UniProt",
    annotation_extension="",
    gene_product_form_id="",
  )
  let content = @src.goa_record_to_gaf_line(rec) + "\n" +
    @src.goa_record_to_gaf_line(rec) + "\n"
  let db = @src.goa_parse(content)
  let ids = @src.goa_unique_go_ids(db)
  assert_eq(ids.length(), 1)
  assert_eq(ids[0], "GO:0003674")
}

test "goa_unique_evidence_codes" {
  let db = @src.goa_sample_database()
  // Sample DB uses IDA, EXP, IEA, ISS, IPI = 5 distinct evidence codes.
  let codes = @src.goa_unique_evidence_codes(db)
  assert_eq(codes.length(), 5)
}

test "goa_unique_taxon_ids" {
  let db = @src.goa_sample_database()
  // Sample DB uses 9606, 4932, 7227 = 3 distinct taxon IDs.
  let ids = @src.goa_unique_taxon_ids(db)
  assert_eq(ids.length(), 3)
  assert_true(ids.contains("9606"))
  assert_true(ids.contains("4932"))
  assert_true(ids.contains("7227"))
}

// ---------------------------------------------------------------------------
// Count functions
// ---------------------------------------------------------------------------

test "goa_count_by_aspect" {
  let db = @src.goa_sample_database()
  let counts = @src.goa_count_by_aspect(db)
  // F: 4, P: 3, C: 3.
  assert_eq(counts.get("F"), Some(4))
  assert_eq(counts.get("P"), Some(3))
  assert_eq(counts.get("C"), Some(3))
  // No Unknown aspect records.
  assert_eq(counts.get("?"), None)
}

test "goa_count_by_aspect_empty_db" {
  let db = @src.GoaDatabase::new()
  let counts = @src.goa_count_by_aspect(db)
  assert_eq(counts.get("F"), None)
  assert_eq(counts.get("P"), None)
  assert_eq(counts.get("C"), None)
}

test "goa_count_by_evidence" {
  let db = @src.goa_sample_database()
  let counts = @src.goa_count_by_evidence(db)
  // IDA: 5, EXP: 1, IEA: 2, ISS: 1, IPI: 1.
  assert_eq(counts.get("IDA"), Some(5))
  assert_eq(counts.get("EXP"), Some(1))
  assert_eq(counts.get("IEA"), Some(2))
  assert_eq(counts.get("ISS"), Some(1))
  assert_eq(counts.get("IPI"), Some(1))
  // Non-existent evidence code returns None.
  assert_eq(counts.get("TAS"), None)
}

test "goa_count_by_evidence_empty_db" {
  let db = @src.GoaDatabase::new()
  let counts = @src.goa_count_by_evidence(db)
  assert_eq(counts.get("IDA"), None)
}

// ---------------------------------------------------------------------------
// Summary generation
// ---------------------------------------------------------------------------

test "goa_to_summary_has_basic_fields" {
  let db = @src.goa_sample_database()
  let s = @src.goa_to_summary(db)
  assert_true(s.contains("GOA Database Summary"))
  assert_true(s.contains("Total records: 10"))
  assert_true(s.contains("GAF version: 2.2"))
  assert_true(s.contains("Created by: UniProt"))
}

test "goa_to_summary_has_aspect_section" {
  let db = @src.goa_sample_database()
  let s = @src.goa_to_summary(db)
  assert_true(s.contains("Counts by aspect"))
  assert_true(s.contains("Molecular Function"))
  assert_true(s.contains("Biological Process"))
  assert_true(s.contains("Cellular Component"))
}

test "goa_to_summary_has_evidence_section" {
  let db = @src.goa_sample_database()
  let s = @src.goa_to_summary(db)
  assert_true(s.contains("Counts by evidence code"))
  assert_true(s.contains("IDA"))
  assert_true(s.contains("EXP"))
  assert_true(s.contains("IEA"))
  assert_true(s.contains("ISS"))
  assert_true(s.contains("IPI"))
}

test "goa_to_summary_has_unique_counts" {
  let db = @src.goa_sample_database()
  let s = @src.goa_to_summary(db)
  assert_true(s.contains("Unique GO IDs: 10"))
  assert_true(s.contains("Unique evidence codes: 5"))
  assert_true(s.contains("Unique taxon IDs: 3"))
}

test "goa_to_summary_empty_db" {
  let db = @src.GoaDatabase::new()
  let s = @src.goa_to_summary(db)
  assert_true(s.contains("Total records: 0"))
  assert_true(s.contains("Unique GO IDs: 0"))
}

// ---------------------------------------------------------------------------
// Writing functions
// ---------------------------------------------------------------------------

test "goa_record_to_gaf_line_roundtrip" {
  let r = @src.GafRecord::new(
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
  )
  let line = @src.goa_record_to_gaf_line(r)
  // Should be 17 tab-separated fields (16 tabs).
  let fields = line.split("\t").to_array()
  assert_eq(fields.length(), 17)
  assert_eq(fields[0].to_string(), "UniProtKB")
  assert_eq(fields[1].to_string(), "Q12345")
  assert_eq(fields[2].to_string(), "PROT1")
  assert_eq(fields[3].to_string(), "enables")
  assert_eq(fields[4].to_string(), "GO:0003674")
  assert_eq(fields[5].to_string(), "PMID:12345")
  assert_eq(fields[6].to_string(), "IDA")
  assert_eq(fields[7].to_string(), "GO:0005515")
  assert_eq(fields[8].to_string(), "F")
  assert_eq(fields[9].to_string(), "Protein 1")
  assert_eq(fields[10].to_string(), "P1|PROT-1")
  assert_eq(fields[11].to_string(), "protein")
  assert_eq(fields[12].to_string(), "taxon:9606|taxon:9606")
  assert_eq(fields[13].to_string(), "20210115")
  assert_eq(fields[14].to_string(), "UniProt")
  assert_eq(fields[15].to_string(), "")
  assert_eq(fields[16].to_string(), "")
}

test "goa_record_to_gaf_line_parse_roundtrip" {
  // Serializing and re-parsing should produce an equivalent record.
  let r = @src.GafRecord::new(
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
  )
  let line = @src.goa_record_to_gaf_line(r)
  let parsed = @src.goa_parse_line(line)
  assert_true(parsed is Some(_))
  match parsed {
    Some(rec) => {
      assert_eq(rec.db(), r.db())
      assert_eq(rec.db_object_id(), r.db_object_id())
      assert_eq(rec.db_object_symbol(), r.db_object_symbol())
      assert_eq(rec.qualifier(), r.qualifier())
      assert_eq(rec.go_id(), r.go_id())
      assert_eq(rec.evidence_code(), r.evidence_code())
      assert_eq(rec.aspect(), r.aspect())
      assert_eq(rec.taxon(), r.taxon())
      assert_eq(rec.assigned_by(), r.assigned_by())
    }
    None => assert_true(false)
  }
}

test "goa_to_gaf_includes_header" {
  let db = @src.goa_sample_database()
  let text = @src.goa_to_gaf(db)
  assert_true(text.contains("!gaf-version: 2.2"))
  assert_true(text.contains("!generated-by: UniProt"))
}

test "goa_to_gaf_includes_all_records" {
  let db = @src.goa_sample_database()
  let text = @src.goa_to_gaf(db)
  // The GAF text should contain every sample record's GO ID.
  assert_true(text.contains("GO:0003674"))
  assert_true(text.contains("GO:0008150"))
  assert_true(text.contains("GO:0005737"))
  assert_true(text.contains("GO:0003700"))
  assert_true(text.contains("GO:0006915"))
  assert_true(text.contains("GO:0005524"))
  assert_true(text.contains("GO:0005634"))
  assert_true(text.contains("GO:0006457"))
  assert_true(text.contains("GO:0005829"))
  assert_true(text.contains("GO:0003735"))
}

test "goa_to_gaf_roundtrip" {
  // Serializing the sample database and re-parsing should preserve its
  // record count and metadata.
  let db = @src.goa_sample_database()
  let text = @src.goa_to_gaf(db)
  let db2 = @src.goa_parse(text)
  assert_eq(db2.n_records(), db.n_records())
  assert_eq(db2.version(), db.version())
  assert_eq(db2.created_by(), db.created_by())
}

test "goa_to_gaf_empty_db" {
  let db = @src.GoaDatabase::new()
  let text = @src.goa_to_gaf(db)
  // Empty database still emits a gaf-version header line.
  assert_true(text.contains("!gaf-version:"))
  // No "generated-by" header when created_by is empty.
  assert_false(text.contains("!generated-by:"))
}

// ---------------------------------------------------------------------------
// Sample database validation
// ---------------------------------------------------------------------------

test "goa_sample_database_has_10_records" {
  let db = @src.goa_sample_database()
  assert_eq(db.n_records(), 10)
}

test "goa_sample_database_metadata" {
  let db = @src.goa_sample_database()
  assert_eq(db.version(), "2.2")
  assert_eq(db.created_by(), "UniProt")
}

test "goa_sample_database_has_3_aspects" {
  let db = @src.goa_sample_database()
  let counts = @src.goa_count_by_aspect(db)
  assert_eq(counts.get("F"), Some(4))
  assert_eq(counts.get("P"), Some(3))
  assert_eq(counts.get("C"), Some(3))
}

test "goa_sample_database_has_3_taxa" {
  let db = @src.goa_sample_database()
  let taxa = @src.goa_unique_taxon_ids(db)
  assert_eq(taxa.length(), 3)
  assert_true(taxa.contains("9606"))
  assert_true(taxa.contains("4932"))
  assert_true(taxa.contains("7227"))
}

test "goa_sample_database_has_5_evidence_codes" {
  let db = @src.goa_sample_database()
  let codes = @src.goa_unique_evidence_codes(db)
  assert_eq(codes.length(), 5)
  assert_true(codes.contains("IDA"))
  assert_true(codes.contains("EXP"))
  assert_true(codes.contains("IEA"))
  assert_true(codes.contains("ISS"))
  assert_true(codes.contains("IPI"))
}

test "goa_sample_database_first_record" {
  let db = @src.goa_sample_database()
  let rec = db.records()[0]
  assert_eq(rec.db(), "UniProtKB")
  assert_eq(rec.db_object_id(), "Q12345")
  assert_eq(rec.db_object_symbol(), "PROT1")
  assert_eq(rec.go_id(), "GO:0003674")
  assert_eq(rec.evidence_code(), "IDA")
  assert_eq(rec.aspect(), "F")
  assert_eq(rec.taxon(), "taxon:9606|taxon:9606")
}

test "goa_sample_database_last_record" {
  let db = @src.goa_sample_database()
  let rec = db.records()[9]
  assert_eq(rec.db_object_id(), "P05388")
  assert_eq(rec.db_object_symbol(), "RPL5")
  assert_eq(rec.go_id(), "GO:0003735")
  assert_eq(rec.evidence_code(), "IPI")
  assert_eq(rec.aspect(), "F")
}

test "goa_sample_database_all_uniprot_db" {
  let db = @src.goa_sample_database()
  for r in db.records() {
    assert_eq(r.db(), "UniProtKB")
  }
}

test "goa_sample_database_unique_go_ids_count" {
  let db = @src.goa_sample_database()
  // Each of the 10 records has a distinct GO ID.
  let ids = @src.goa_unique_go_ids(db)
  assert_eq(ids.length(), 10)
}

// ---------------------------------------------------------------------------
// Edge cases
// ---------------------------------------------------------------------------

test "goa_parse_line_single_column" {
  // A single-column line should be padded to 17 empty fields and parsed.
  let r = @src.goa_parse_line("UniProtKB")
  assert_true(r is Some(_))
  match r {
    Some(rec) => {
      assert_eq(rec.db(), "UniProtKB")
      assert_eq(rec.db_object_id(), "")
      assert_eq(rec.go_id(), "")
      assert_eq(rec.aspect(), "")
    }
    None => assert_true(false)
  }
}

test "goa_parse_line_extra_columns_kept_in_field" {
  // A line with more than 17 columns: the parser splits on all tabs and uses
  // only the first 17 fields; extra columns are silently ignored.
  let line = "UniProtKB\tQ12345\tPROT1\tenables\tGO:0003674\tPMID:12345\tIDA\tGO:0005515\tF\tProtein 1\tP1\tprotein\ttaxon:9606\t20210115\tUniProt\t\textra1\textra2"
  let r = @src.goa_parse_line(line)
  assert_true(r is Some(_))
  match r {
    Some(rec) => {
      assert_eq(rec.db(), "UniProtKB")
      assert_eq(rec.db_object_id(), "Q12345")
      assert_eq(rec.go_id(), "GO:0003674")
    }
    None => assert_true(false)
  }
}

test "goa_parse_line_with_trailing_newline" {
  // Trailing whitespace/newlines are trimmed before parsing.
  let line = "UniProtKB\tQ12345\tPROT1\tenables\tGO:0003674\tPMID:12345\tIDA\tGO:0005515\tF\tProtein 1\tP1\tprotein\ttaxon:9606\t20210115\tUniProt\t\t\n"
  let r = @src.goa_parse_line(line)
  assert_true(r is Some(_))
  match r {
    Some(rec) => {
      assert_eq(rec.db_object_symbol(), "PROT1")
      assert_eq(rec.go_id(), "GO:0003674")
    }
    None => assert_true(false)
  }
}

test "goa_parse_only_comment_lines" {
  // A line starting with "!" is treated as a header/comment regardless of
  // whether it has a colon, and never produces a data record.
  let content = "!gaf-version: 2.2\n" +
    "!generated-by: UniProt\n" +
    "!this is just a comment with no colon\n"
  let db = @src.goa_parse(content)
  assert_eq(db.n_records(), 0)
  assert_eq(db.version(), "2.2")
  assert_eq(db.created_by(), "UniProt")
}

test "goa_parse_mixed_valid_and_invalid_lines" {
  let content = "!gaf-version: 2.2\n" +
    "!generated-by: UniProt\n" +
    "\n" +
    "!comment line\n" +
    "\n" +
    "UniProtKB\tQ12345\tPROT1\tenables\tGO:0003674\tPMID:12345\tIDA\t\tF\tProtein 1\t\tprotein\ttaxon:9606\t20210115\tUniProt\t\t\n" +
    "\n" +
    "   \n" +
    "UniProtKB\tP04637\tTP53\tenables\tGO:0003700\tPMID:67890\tEXP\t\tF\tCellular tumor antigen p53\t\tprotein\ttaxon:9606\t20210220\tUniProt\t\t\n"
  let db = @src.goa_parse(content)
  assert_eq(db.n_records(), 2)
  assert_eq(db.version(), "2.2")
  assert_eq(db.created_by(), "UniProt")
}

test "goa_filter_on_empty_database" {
  let db = @src.GoaDatabase::new()
  assert_eq(@src.goa_filter_by_go_id(db, "GO:0003674").length(), 0)
  assert_eq(@src.goa_filter_by_aspect(db, @src.GafAspect::from_string("F")).length(), 0)
  assert_eq(@src.goa_filter_by_evidence(db, "IDA").length(), 0)
  assert_eq(@src.goa_filter_by_taxon(db, "9606").length(), 0)
  assert_eq(@src.goa_filter_by_db_object_id(db, "Q12345").length(), 0)
}

test "goa_unique_on_empty_database" {
  let db = @src.GoaDatabase::new()
  assert_eq(@src.goa_unique_go_ids(db).length(), 0)
  assert_eq(@src.goa_unique_evidence_codes(db).length(), 0)
  assert_eq(@src.goa_unique_taxon_ids(db).length(), 0)
}
