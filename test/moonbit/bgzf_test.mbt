///|
/// Test BGZF decompression functionality.

///|
/// Test bgzf_decompress with empty data.
test "bgzf_decompress_empty" {
  let result = @src.bgzf_decompress("")
  assert_eq(result, "")
}

///|
/// Test is_bgzf function.
test "bgzf_is_bgzf" {
  let bgzf_magic = "\x1f\x8b\x08\x04"
  assert_true(@src.is_bgzf(bgzf_magic))

  let not_bgzf = "BAM\x01"
  assert_true(!@src.is_bgzf(not_bgzf))
}
