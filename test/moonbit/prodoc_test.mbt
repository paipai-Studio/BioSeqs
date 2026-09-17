///|
/// Tests for the Bio.ExPASy.Prodoc port. Expected record values mirror
/// Biopython's Prodoc.parse/read output on the same sample documents
/// (cross-checked against Biopython 1.88).

///|
test "prodoc_parse_two_records" {
  let records = @src.prodoc_parse(@src.prodoc_sample_content()) catch {
    _ => abort("sample prosite.doc should parse")
  }
  assert_eq(records.length(), 2)

  let r0 = records[0]
  assert_eq(r0.accession, "PDOC00001")
  assert_eq(r0.prosite_refs.length(), 2)
  assert_eq(r0.prosite_refs[0].0, "PS00001")
  assert_eq(r0.prosite_refs[0].1, "ASN_GLYCOSYLATION")
  assert_eq(r0.prosite_refs[1].0, "PS00007")
  assert_eq(r0.prosite_refs[1].1, "PROTEIN_KINASE_ATP")
  assert_eq(r0.text, "N-glycosylation site documentation.\nSecond paragraph line.\n")
  assert_eq(r0.references.length(), 2)
  assert_eq(r0.references[0].number, "1")
  assert_eq(r0.references[0].authors, "Hunter, T., Plowman, G.D.")
  assert_eq(
    r0.references[0].citation,
    "In \"Methods in Enzymology\", Academic Press, New York, (1991).",
  )
  assert_eq(r0.references[1].number, "E1")
  assert_eq(r0.references[1].authors, "")
  assert_eq(r0.references[1].citation, "https://www.ebi.ac.uk/interpro/entry/PDOC00001/")

  let r1 = records[1]
  assert_eq(r1.accession, "PDOC00002")
  assert_eq(r1.prosite_refs.length(), 1)
  assert_eq(r1.prosite_refs[0].0, "PS00002")
  assert_eq(r1.prosite_refs[0].1, "GLYCOSAMINOGLYCAN")
  assert_eq(r1.text, "Glycosaminoglycan attachment site.\n")
  assert_eq(r1.references.length(), 1)
  assert_eq(r1.references[0].number, "1")
  assert_eq(r1.references[0].authors, "Bourdon, M.A., Krusius, T., Ruoslahti, E.")
  assert_eq(r1.references[0].citation, "Identification and synthesis of a recognition signal.")
}

///|
test "prodoc_parse_skips_copyright" {
  let records = @src.prodoc_parse(@src.prodoc_sample_content()) catch {
    _ => abort("sample prosite.doc should parse")
  }
  // The copyright box lines are not part of the free text.
  assert_true(!@src.contains(records[0].text, "Copyright"))
}

///|
test "prodoc_read_single_record" {
  let content =
    "{PDOC00001}\n{PS00001; ASN_GLYCOSYLATION}\n{BEGIN}\nSome text.\n{END}\n"
  let r = @src.prodoc_read(content) catch {
    _ => abort("single record should parse")
  }
  match r {
    Some(record) => {
      assert_eq(record.accession, "PDOC00001")
      assert_eq(record.text, "Some text.\n")
    }
    None => abort("expected one record")
  }
}

///|
test "prodoc_read_empty_document" {
  let r = @src.prodoc_read("") catch {
    _ => abort("empty document should not raise")
  }
  match r {
    None => ()
    Some(_) => abort("expected None for empty document")
  }
}

///|
test "prodoc_read_multiple_records_raises" {
  let r = try {
    ignore(@src.prodoc_read(@src.prodoc_sample_content()))
    false
  } catch {
    ProdocError(msg) => msg == "More than one Prodoc record found"
  }
  assert_true(r)
}

///|
test "prodoc_bad_accession_raises" {
  let r = try {
    ignore(@src.prodoc_parse("BOGUS\n{BEGIN}\n{END}\n"))
    false
  } catch {
    ProdocError(_) => true
  }
  assert_true(r)
}

///|
test "prodoc_read_trailing_separator_is_extra_content" {
  // A single record followed by its "//" terminator still counts as extra
  // content for read(), matching Biopython's single-readline() check.
  let content = "{PDOC00001}\n{BEGIN}\nSome text.\n{END}\n//\n"
  let r = try {
    ignore(@src.prodoc_read(content))
    false
  } catch {
    ProdocError(msg) => msg == "More than one Prodoc record found"
  }
  assert_true(r)
}