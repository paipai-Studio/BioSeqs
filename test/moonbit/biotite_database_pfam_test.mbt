///|
/// Tests for biotite_database_pfam.mbt.
test "biotite_pfam_urls" {
  let entry = @src.biotite_pfam_entry_url("PF00069")
  assert_true(entry.contains("PF00069"))
  assert_true(entry.contains("interpro"))
  let alignment = @src.biotite_pfam_alignment_url("PF00069")
  assert_true(alignment.contains("PF00069"))
  assert_true(alignment.contains("stockholm"))
}

///|
test "biotite_pfam_parse_fields" {
  let json = "{\"metadata\":{\"accession\":\"PF00069\",\"description\":\"Protein kinase domain\"}}"
  match @src.biotite_pfam_parse_accession(json) {
    Some(acc) => assert_eq(acc, "PF00069")
    None => fail("expected accession")
  }
  match @src.biotite_pfam_parse_description(json) {
    Some(desc) => assert_eq(desc, "Protein kinase domain")
    None => fail("expected description")
  }
}

///|
test "biotite_pfam_missing_fields" {
  assert_true(@src.biotite_pfam_parse_description("{}") is None)
  assert_true(@src.biotite_pfam_parse_accession("{}") is None)
}
