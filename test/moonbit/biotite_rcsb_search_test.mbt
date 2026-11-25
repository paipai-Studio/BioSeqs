///|
/// Tests for biotite_rcsb_search.mbt — RCSB query construction.
test "biotite_rcsb_query_json_and_url" {
  let query = @src.BiotiteRcsbSearchQuery::new()
  let with_method = query.add_condition(
    @src.BiotiteRcsbExperimentalMethod("X-RAY DIFFRACTION"),
  )
  let json = with_method.build_json()
  assert_true(json.contains("X-RAY DIFFRACTION"))
  assert_true(json.contains("\"return_type\":\"entry\""))
  let url = with_method.build_url()
  assert_true(url.has_prefix("https://search.rcsb.org/rcsbsearch/v2/query"))
}

///|
test "biotite_rcsb_resolution_condition" {
  let query = @src.BiotiteRcsbSearchQuery::new().add_condition(
    @src.BiotiteRcsbResolution(1.0, 2.5),
  )
  let json = query.build_json()
  assert_true(json.contains("resolution_combined"))
  assert_true(json.contains("\"from\":1"))
  assert_true(json.contains("\"to\":2.5"))
}

///|
test "biotite_rcsb_sort_and_return_type" {
  let query = @src.BiotiteRcsbSearchQuery::new()
    .add_condition(@src.BiotiteRcsbOrganism("Escherichia coli"))
    .set_return_type("polymer_entity")
    .set_sort("score", sort_dir="desc")
  let json = query.build_json()
  assert_true(json.contains("polymer_entity"))
  assert_true(json.contains("\"sort_by\":\"score\""))
  assert_true(json.contains("\"direction\":\"desc\""))
}

///|
test "biotite_rcsb_convenience_builders" {
  assert_true(@src.biotite_rcsb_search_by_method("NMR").contains("NMR"))
  assert_true(@src.biotite_rcsb_search_by_resolution(1.5, 3.0).contains("1.5"))
}

///|
test "biotite_rcsb_empty_query" {
  let json = @src.BiotiteRcsbSearchQuery::new().build_json()
  assert_true(json.contains("\"nodes\":[]"))
}
