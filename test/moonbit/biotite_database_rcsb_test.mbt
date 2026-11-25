///|
/// Tests for biotite_database_rcsb.mbt.
test "biotite_rcsb_download_urls" {
  assert_eq(
    @src.biotite_rcsb_fetch_pdb_url("1ake"),
    "https://files.rcsb.org/download/1ake.pdb",
  )
  assert_true(@src.biotite_rcsb_fetch_cif_url("1ake").has_suffix("1ake.cif"))
  assert_true(@src.biotite_rcsb_entry_api_url("1ake").has_suffix("/1ake"))
  assert_true(@src.biotite_rcsb_assembly_api_url("1ake").has_suffix("/1ake"))
}

///|
test "biotite_rcsb_parse_assembly_ids" {
  let json = "{\"assembly_id\":\"1\"},{\"assembly_id\":\"2\"}"
  let ids = @src.biotite_rcsb_parse_assembly_ids(json)
  assert_eq(ids.length(), 2)
  assert_eq(ids[0], "1")
  assert_eq(ids[1], "2")
  assert_eq(@src.biotite_rcsb_parse_assembly_ids("{}").length(), 0)
}

///|
test "biotite_rcsb_parse_title" {
  let json = "{\"struct\":{\"title\":\"A test structure\"}}"
  match @src.biotite_rcsb_parse_title(json) {
    Some(title) => assert_eq(title, "A test structure")
    None => fail("expected title")
  }
  assert_true(@src.biotite_rcsb_parse_title("{}") is None)
}
