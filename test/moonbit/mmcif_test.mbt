///|
/// Tests for mmCIF module.
test "mmcif_parse_basic" {
  let content = @src.sample_mmcif_content()
  let block = @src.parse_mmcif(content)
  assert_eq(block.name, "1ABC")
}

///|
test "mmcif_empty_content" {
  let block = @src.parse_mmcif("")
  assert_eq(block.name, "default")
}
