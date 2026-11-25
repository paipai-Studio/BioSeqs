///|
/// Tests for biotite_database_dob.mbt.
test "biotite_dob_entry_url" {
  let client = @src.biotite_dob_client("https://example.org")
  assert_eq(
    @src.biotite_dob_entry_url(client, "abc"),
    "https://example.org/entry/abc",
  )
}

///|
test "biotite_dob_search_url_encodes_query" {
  let client = @src.biotite_dob_client("https://example.org")
  let url = @src.biotite_dob_search_url(client, "a b")
  assert_true(url.contains("q=a%20b"))
}

///|
test "biotite_dob_client_custom_base" {
  let client = @src.biotite_dob_client("http://localhost:8080")
  assert_true(
    @src.biotite_dob_entry_url(client, "x").has_prefix("http://localhost:8080"),
  )
}
