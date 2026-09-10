// Whitebox tests for CRC32/CRC64 checksums (spec 5.6).
// Reference values computed with Biopython 1.88 (crc32 = binascii.crc32,
// crc64 = Bio.SeqUtils.CheckSum.crc64) on the upper-cased sequence, since
// our implementations upper-case before hashing (spec 5.6.1.3).

///|
test "checksum_crc32_matches_biopython" {
  let cases : Array[(String, String)] = [
    ("", "0"),
    ("A", "3554254475"),
    ("ACGT", "2735644658"),
    ("ATGC", "968142784"),
    ("ACGTACGTACGT", "20049947"),
    ("ACDEFGHIKLMNPQRSTVWY", "149026042"),
    (
      "MKWVTFISLLFLFSSAYSRGVFRRDTHKSEIAHRFKDLGEEHFKGLVLIAFSQYLQQCPFDEHVK",
      "4288159397",
    ),
  ]
  for case in cases {
    let (seq, expected) = case
    let result = @src.checksum_crc32(seq)
    assert_eq(result.checksum, expected)
    assert_eq(result.method, "crc32")
  }
}

///|
test "checksum_crc64_matches_biopython" {
  let cases : Array[(String, String)] = [
    ("", "CRC-0000000000000000"),
    ("A", "CRC-6DB0000000000000"),
    ("ACGT", "CRC-71A87EBDB0000000"),
    ("ATGC", "CRC-6EB861ADB0000000"),
    ("ACGTACGTACGT", "CRC-C4FBB762C4A87EBD"),
    ("ACDEFGHIKLMNPQRSTVWY", "CRC-F0D170F3DFA2290A"),
    (
      "MKWVTFISLLFLFSSAYSRGVFRRDTHKSEIAHRFKDLGEEHFKGLVLIAFSQYLQQCPFDEHVK",
      "CRC-475F23DE3366B7A2",
    ),
  ]
  for case in cases {
    let (seq, expected) = case
    let result = @src.checksum_crc64(seq)
    assert_eq(result.checksum, expected)
    assert_eq(result.method, "crc64")
  }
}

///|
test "checksum_crc_case_insensitive" {
  // Two sequences differing only in case produce identical checksums
  // (spec 5.6.1.3)
  assert_eq(
    @src.checksum_crc32("acgt").checksum,
    @src.checksum_crc32("ACGT").checksum,
  )
  assert_eq(
    @src.checksum_crc64("acgt").checksum,
    @src.checksum_crc64("ACGT").checksum,
  )
  assert_eq(
    @src.checksum_crc64("acgtACGTacgt").checksum,
    @src.checksum_crc64("ACGTACGTACGT").checksum,
  )
}

///|
test "checksum_crc_empty_is_deterministic" {
  assert_eq(@src.checksum_crc32("").checksum, @src.checksum_crc32("").checksum)
  assert_eq(@src.checksum_crc64("").checksum, @src.checksum_crc64("").checksum)
  assert_eq(@src.checksum_crc32("").checksum, "0")
  assert_eq(@src.checksum_crc64("").checksum, "CRC-0000000000000000")
}

///|
test "checksum_by_name_dispatch" {
  let seq = "ACGT"
  assert_eq(@src.checksum_by_name(seq, "crc32").checksum, @src.checksum_crc32(seq).checksum)
  assert_eq(@src.checksum_by_name(seq, "crc64").checksum, @src.checksum_crc64(seq).checksum)
  assert_eq(@src.checksum_by_name(seq, "CRC32").checksum, @src.checksum_crc32(seq).checksum)
  assert_eq(@src.checksum_by_name(seq, "CRC64").checksum, @src.checksum_crc64(seq).checksum)
  // Existing algorithms keep working through the dispatcher
  assert_eq(@src.checksum_by_name(seq, "gcg").checksum, @src.checksum_gcg(seq).checksum)
  assert_eq(
    @src.checksum_by_name(seq, "seguid").checksum,
    @src.checksum_seguid(seq).checksum,
  )
  let raised = try {
    ignore(@src.checksum_by_name(seq, "md5"))
    false
  } catch {
    SeqIOError(msg) => msg == "Unknown checksum algorithm: md5"
  }
  assert_true(raised)
}