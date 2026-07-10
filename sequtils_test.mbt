///|
test "crc32_basic" {
  let seq = "AGCTTTTCATTCTGACTGCAACGGGCAATATGTCTCTGTGTGGATTAAAAAAAGAGTGTCTGATAGCAGC"
  let result = crc32(seq)
  let expected = -47322090
  assert_eq(result & 0xFFFFFFFF, expected & 0xFFFFFFFF)
}

///|
test "crc32_simple" {
  let abc_expected = -1551695032
  let agct_expected = -1066667478
  assert_eq(crc32("ABC") & 0xFFFFFFFF, abc_expected & 0xFFFFFFFF)
  assert_eq(crc32("AGCT") & 0xFFFFFFFF, agct_expected & 0xFFFFFFFF)
}

///|
test "gcg_basic" {
  let seq = "AGCTTTTCATTCTGACTGCAACGGGCAATATGTCTCTGTGTGGATTAAAAAAAGAGTGTCTGATAGCAGC"
  let result = gcg(seq)
  assert_eq(result, 5493)
}

///|
test "seguid_basic" {
  let seq = "AGCTTTTCATTCTGACTGCAACGGGCAATATGTCTCTGTGTGGATTAAAAAAAGAGTGTCTGATAGCAGC"
  let result = seguid(seq)
  assert_eq(result, "kspO3Jbmvk9xX3jVv8QOSapllIo")
}

///|
test "sha1_basic" {
  let hash = sha1_string("abc")
  assert_eq(hash.length(), 20)
  assert_eq(hash[0], 169)
  assert_eq(hash[1], 153)
  assert_eq(hash[2], 62)
}

test "sha1_h0_after_round" {
  let hash = sha1_string("abc")
  let expected_h0_bytes = [169, 153, 62, 54]
  assert_eq(hash[0], expected_h0_bytes[0])
  assert_eq(hash[1], expected_h0_bytes[1])
  assert_eq(hash[2], expected_h0_bytes[2])
  assert_eq(hash[3], expected_h0_bytes[3])
}

test "sha1_raw_bytes" {
  let bytes = [97, 98, 99]
  let hash = sha1(bytes)
  assert_eq(hash[0], 169)
  assert_eq(hash[1], 153)
  assert_eq(hash[2], 62)
}

test "base64_test" {
  let bytes = [169, 153, 62, 54, 71, 6, 129, 106, 186, 62, 37, 113, 120, 80, 194, 108, 156, 208, 216, 157]
  let result = bytes_to_base64(bytes)
  assert_eq(result, "q71HzGOFBjaYuzVZUDKe7APTKFU")
}

test "string_to_bytes_test" {
  let bytes = string_to_bytes("abc")
  assert_eq(bytes[0], 97)
  assert_eq(bytes[1], 98)
  assert_eq(bytes[2], 99)
}

///|
test "gc_simple" {
  let result = gc("AGCT")
  assert_true(result > 49.9 && result < 50.1)
}

///|
test "gc_full" {
  let seq = "AGCTTTTCATTCTGACTGCAACGGGCAATATGTCTCTGTGTGGATTAAAAAAAGAGTGTCTGATAGCAGC"
  let result = gc(seq)
  assert_true(result > 41.0 && result < 42.0)
}

///|
test "protein_molecular_weight" {
  let protein_seq = "MAEGEITTFTALTEKFNLPPGNYKKPKLLYCSNGGHFLRILPDGTVDGTILTRLRLY"
  let protein = ProteinAnalysis::new(protein_seq)
  let result = protein.molecular_weight()
  assert_true(result > 6380.0 && result < 6390.0)
}

///|
test "protein_aromaticity" {
  let protein_seq = "MAEGEITTFTALTEKFNLPPGNYKKPKLLYCSNGGHFLRILPDGTVDGTILTRLRLY"
  let protein = ProteinAnalysis::new(protein_seq)
  let result = protein.aromaticity()
  assert_true(result > 0.1 && result < 0.11)
}

///|
test "protein_gravy" {
  let protein_seq = "MAEGEITTFTALTEKFNLPPGNYKKPKLLYCSNGGHFLRILPDGTVDGTILTRLRLY"
  let protein = ProteinAnalysis::new(protein_seq)
  let result = protein.gravy()
  assert_true(result > -0.2 && result < -0.1)
}

///|
test "tm_wallace_simple" {
  let result = tm_wallace("AAAAATTTTT")
  assert_eq(result, 20.0)
}

///|
test "tm_wallace_gc" {
  let result = tm_wallace("GGCCGGCC")
  assert_eq(result, 32.0)
}

///|
test "tm_gc_simple" {
  let result = tm_gc("AGCTAGCT")
  assert_true(result > 5.0 && result < 6.0)
}

///|
test "tm_gc_long" {
  let seq = "AGCTTTTCATTCTGACTGCAACGGGCAATATGTCTCTGTGTGGATTAAAAAAAGAGTGTCTGATAGCAGC"
  let result = tm_gc(seq)
  assert_true(result > 68.0 && result < 69.0)
}