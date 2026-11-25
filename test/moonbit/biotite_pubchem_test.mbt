///|
/// Tests for biotite_pubchem.mbt — query URL builders and JSON parsing.
test "biotite_pubchem_query_urls" {
  let by_name = @src.biotite_pubchem_query_url(
    @src.BiotitePubChemByName("aspirin"),
    ["MolecularFormula"],
  )
  assert_true(by_name.contains("name/aspirin"))
  assert_true(by_name.contains("properties=MolecularFormula"))
  let by_cid = @src.biotite_pubchem_query_url(@src.BiotitePubChemByCid(2244), [])
  assert_true(by_cid.contains("cid/2244"))
  assert_true(by_cid.contains("/JSON"))
}

///|
test "biotite_pubchem_properties_and_sdf_urls" {
  let props = @src.biotite_pubchem_properties_url(2244, [
    "MolecularFormula", "MolecularWeight",
  ])
  assert_true(props.contains("/cid/2244/property/"))
  assert_true(props.contains("MolecularWeight"))
  let sdf = @src.biotite_pubchem_sdf_url(2244)
  assert_true(sdf.contains("2244"))
  assert_true(sdf.contains("SDF"))
}

///|
test "biotite_pubchem_json_parsing" {
  let json = "{\"MolecularFormula\":\"C9H8O4\",\"CID\":2244}"
  match @src.biotite_pubchem_parse_json_string(json, "MolecularFormula") {
    Some(v) => assert_eq(v, "C9H8O4")
    None => fail("expected string property")
  }
  match @src.biotite_pubchem_parse_json_number(json, "CID") {
    Some(v) => assert_true(@moonbit.assert_num_close(v, 2244.0, 1.0e-9))
    None => fail("expected numeric property")
  }
}

///|
test "biotite_pubchem_json_missing_key" {
  let json = "{\"CID\":2244}"
  assert_true(@src.biotite_pubchem_parse_json_string(json, "Missing") is None)
  assert_true(@src.biotite_pubchem_parse_json_number(json, "Missing") is None)
}
