///|
/// Tests for biotite_database_fetch.mbt.
test "biotite_database_fetch_url_dispatch" {
  assert_eq(
    @src.biotite_database_fetch_url("rcsb", "1ake"),
    "https://files.rcsb.org/download/1ake.pdb",
  )
  assert_true(
    @src.biotite_database_fetch_url("pdb", "1ake").has_suffix("1ake.pdb"),
  )
  assert_true(
    @src.biotite_database_fetch_url("pfam", "PF00069").contains("PF00069"),
  )
  assert_eq(@src.biotite_database_fetch_url("unknown", "x"), "")
}

///|
test "biotite_database_fetch_url_with_base" {
  let url = @src.biotite_database_fetch_url_with_base(
    "dob", "abc", "https://example.org",
  )
  assert_eq(url, "https://example.org/entry/abc")
  assert_eq(
    @src.biotite_database_fetch_url_with_base("rcsb", "1ake", "ignored"),
    "https://files.rcsb.org/download/1ake.pdb",
  )
}

///|
test "biotite_database_cache_path" {
  let path = @src.biotite_database_cache_path(
    "/tmp/cache", "https://files.rcsb.org/download/1ake.pdb",
  )
  assert_true(path.has_prefix("/tmp/cache/"))
  assert_true(path.has_suffix(".dat"))
  assert_true(!path.contains("://"))
}
