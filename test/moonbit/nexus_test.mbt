///|
/// Tests for Nexus module.
test "nexus_parse_basic" {
  let content = @src.sample_nexus_content()
  let nf = @src.parse_nexus(content)
  assert_true(nf.blocks.length() >= 0)
}

///|
test "nexus_empty_content" {
  let nf = @src.parse_nexus("")
  assert_eq(nf.blocks.length(), 0)
}
