///|
test "checksum_gcg_basic" {
  let seq = "ACDEFGHIKLMNPQRSTVWY"
  let result = @src.checksum_gcg(seq)
  assert_eq(result.method, "gcg")
  assert_true(result.checksum.length() == 4)
}

///|
test "checksum_gcg_simple" {
  let seq = "ABC"
  let result = @src.checksum_gcg(seq)
  assert_eq(result.method, "gcg")
  assert_true(result.checksum.length() == 4)
}

///|
test "checksum_seguid_basic" {
  let seq = "ACGT"
  let result = @src.checksum_seguid(seq)
  assert_eq(result.method, "seguid")
  assert_true(result.checksum.length() > 0)
}

///|
test "checksum_seguid_dna" {
  let seq = "AGCTTTTCATTCTGACTGCAACGGGCAATATGTCTCTGTGTGGATTAAAAAAAGAGTGTCTGATAGCAGC"
  let result = @src.checksum_seguid(seq)
  assert_eq(result.method, "seguid")
  assert_true(result.checksum.length() > 0)
}

///|
test "create_example_sequence" {
  let seq = @src.create_example_sequence()
  assert_eq(seq, "ACDEFGHIKLMNPQRSTVWY")
}

///|
test "verify_checksum_gcg" {
  let seq = "ACDEFGHIKLMNPQRSTVWY"
  let result = @src.checksum_gcg(seq)
  let verified = @src.verify_checksum(seq, result.checksum, "gcg")
  assert_true(verified)
}

///|
test "verify_checksum_seguid" {
  let seq = "ACGT"
  let result = @src.checksum_seguid(seq)
  let verified = @src.verify_checksum(seq, result.checksum, "seguid")
  assert_true(verified)
}