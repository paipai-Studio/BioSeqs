///|
/// Tests for Bio.ExPASy.cellosaurus-compatible parsing.

///|
test "cellosaurus_parse_official_sample_fields" {
  let records = @src.cellosaurus_parse(@src.cellosaurus_sample_text())
  assert_eq(records.length(), 1)
  let record = records[0]
  assert_eq(record.identifier, "#15310-LN")
  assert_eq(record.accession, "CVCL_E548")
  assert_eq(record.sex, "Female")
  assert_eq(record.age, "Age unspecified")
  assert_eq(record.category, "Transformed cell line")
  assert_eq(record.comments.length(), 2)
  assert_eq(record.species.length(), 1)
}

///|
test "cellosaurus_parse_cross_references" {
  let record = @src.cellosaurus_parse(@src.cellosaurus_sample_text())[0]
  assert_eq(record.cross_references.length(), 3)
  assert_eq(record.cross_references[0].database, "dbMHC")
  assert_eq(record.cross_references[0].accession, "48439")
  assert_eq(record.cross_references[2].database, "Wikidata")
  assert_eq(record.cross_references[2].accession, "Q54398957")
}

///|
test "cellosaurus_list_helpers" {
  let record = @src.CellosaurusRecord::new(
    identifier="XP3OS",
    secondary_accessions="CVCL_F511; CVCL_TEST",
    synonyms="Xeroderma Pigmentosum 3 OSaka; GM04314; GM4314",
  )
  assert_eq(record.secondary_accession_list(), ["CVCL_F511", "CVCL_TEST"])
  assert_eq(record.synonym_list(), [
    "Xeroderma Pigmentosum 3 OSaka", "GM04314", "GM4314",
  ])
}

///|
test "cellosaurus_cross_reference_filter" {
  let text = "ID   Example\n" +
    "AC   CVCL_0001\n" +
    "DR   JCRB; JCRB0303\n" +
    "DR   JCRB; KURB1002\n" +
    "DR   Wikidata; Q1\n" +
    "//\n"
  let record = @src.cellosaurus_parse(text)[0]
  assert_eq(record.cross_references_for("JCRB"), ["JCRB0303", "KURB1002"])
  assert_eq(record.cross_references_for("ATCC"), [])
}

///|
test "cellosaurus_species_and_repr" {
  let record = @src.cellosaurus_parse(@src.cellosaurus_sample_text())[0]
  assert_true(record.has_species("Homo sapiens"))
  assert_false(record.has_species("Mus musculus"))
  assert_eq(record.repr(), "CellosaurusRecord (#15310-LN, CVCL_E548)")
}

///|
test "cellosaurus_parse_multiple_records" {
  let text = "ID   XP3OS\n" +
    "AC   CVCL_3245\n" +
    "AS   CVCL_F511\n" +
    "RX   PubMed=832273;\n" +
    "ST   Amelogenin: X\n" +
    "DI   ORDO; Orphanet_910; Xeroderma pigmentosum\n" +
    "OX   NCBI_TaxID=9606; ! Homo sapiens (Human)\n" +
    "//\n" +
    "ID   1-5c-4\n" +
    "AC   CVCL_2260\n" +
    "HI   CVCL_0030 ! HeLa\n" +
    "OI   CVCL_0002 ! Same donor\n" +
    "CA   Cancer cell line\n" +
    "//\n"
  let records = @src.cellosaurus_parse(text)
  assert_eq(records.length(), 2)
  assert_eq(records[0].secondary_accessions, "CVCL_F511")
  assert_eq(records[0].reference_identifiers, ["PubMed=832273;"])
  assert_eq(records[0].str_profiles, ["Amelogenin: X"])
  assert_eq(records[0].diseases.length(), 1)
  assert_eq(records[1].identifier, "1-5c-4")
  assert_eq(records[1].hierarchy, ["CVCL_0030 ! HeLa"])
  assert_eq(records[1].same_individual, ["CVCL_0002 ! Same donor"])
}

///|
test "cellosaurus_read_single" {
  match @src.cellosaurus_read(@src.cellosaurus_sample_text()) {
    Some(record) => assert_eq(record.accession, "CVCL_E548")
    None => fail("expected one Cellosaurus record")
  }
}

///|
test "cellosaurus_read_empty" {
  assert_true(@src.cellosaurus_read("") is None)
  assert_true(@src.cellosaurus_read("//\n") is None)
}

///|
test "cellosaurus_read_multiple_error" {
  let text = "ID   First\nAC   CVCL_0001\n//\n" +
    "ID   Second\nAC   CVCL_0002\n//\n"
  let raised = try {
    ignore(@src.cellosaurus_read(text))
    false
  } catch {
    CellosaurusError(_) => true
  }
  assert_true(raised)
}

///|
test "cellosaurus_unterminated_record_error" {
  let raised = try {
    ignore(@src.cellosaurus_parse("ID   Incomplete\nAC   CVCL_0001\n"))
    false
  } catch {
    CellosaurusError(_) => true
  }
  assert_true(raised)
}

///|
test "cellosaurus_new_id_before_terminator_error" {
  let text = "ID   First\nAC   CVCL_0001\nID   Second\nAC   CVCL_0002\n//\n"
  let raised = try {
    ignore(@src.cellosaurus_parse(text))
    false
  } catch {
    CellosaurusError(_) => true
  }
  assert_true(raised)
}

///|
test "cellosaurus_invalid_cross_reference_error" {
  let raised = try {
    ignore(@src.cellosaurus_parse("ID   Bad\nDR   missing delimiter\n//\n"))
    false
  } catch {
    CellosaurusError(_) => true
  }
  assert_true(raised)
}

///|
test "cellosaurus_header_unknown_and_crlf" {
  let text = "Cellosaurus release header\r\n" +
    "__________\r\n" +
    "//\r\n" +
    "ID   CRLF example\r\n" +
    "AC   CVCL_TEST\r\n" +
    "ZZ   ignored extension\r\n" +
    "OX   NCBI_TaxID=10090; ! Mus musculus (Mouse)\r\n" +
    "//\r\n"
  let records = @src.cellosaurus_parse(text)
  assert_eq(records.length(), 1)
  assert_eq(records[0].identifier, "CRLF example")
  assert_eq(records[0].accession, "CVCL_TEST")
  assert_true(records[0].has_species("Mus musculus"))
}

///|
test "cellosaurus_repeated_single_value_fields_concatenate" {
  let text = "ID   Split\n" +
    "AC   CVCL_\n" +
    "AC   1234\n" +
    "SY   First; \n" +
    "SY   Second\n" +
    "//\n"
  let record = @src.cellosaurus_parse(text)[0]
  assert_eq(record.accession, "CVCL_1234")
  assert_eq(record.synonyms, "First;Second")
}

///|
test "cellosaurus_serialization_round_trip" {
  let original = @src.cellosaurus_parse(@src.cellosaurus_sample_text())[0]
  let serialized = original.to_string()
  let reparsed = @src.cellosaurus_parse(serialized)
  assert_eq(reparsed.length(), 1)
  assert_true(reparsed[0] == original)
}

///|
test "cellosaurus_empty_constructor_repr" {
  let record = @src.CellosaurusRecord::new(identifier="")
  assert_eq(record.repr(), "CellosaurusRecord ( )")
  assert_eq(record.to_string(), "//\n")
}
