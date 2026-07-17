///|
/// Test BAM format parsing.
/// Uses base64-encoded raw BAM data for testing.

///|
/// Base64-encoded minimal BAM file for testing.
/// Generated from:
/// @HD	VN:1.6	SO:coordinate
/// @SQ	SN:chr1	LN:1000
/// read1	0	chr1	100	60	8M	*	0	0	AAAAAAAA	!!!!!!!!
let test_bam_base64 = "QkFNAS0AAABASEQJVk46MS42CVNPOmNvb3JkaW5hdGUKQFNRCVNOOmNocjEJTE46MTAwMAoBAAAABQAAAGNocjEA6AMAADkAAAAAAAAAZAAAAAAAAAA8AAAGAAEACAAAAP////8AAAAAAAAAAHJlYWQxAIAAAAAAAAAAAAAAAAAAAAA="

///|
test "bam_parse_header" {
  let bam = @src.parse_bam_base64(test_bam_base64)

  assert_eq(bam.num_references(), 1)
  let ref_names = bam.get_reference_names()
  assert_eq(ref_names.length(), 1)
  assert_eq(ref_names[0], "chr1")
  assert_eq(bam.get_reference_length("chr1"), Some(1000))
}

///|
test "bam_parse_record" {
  let bam = @src.parse_bam_base64(test_bam_base64)

  assert_eq(bam.len(), 1)

  let record = bam.records[0]
  assert_eq(record.qname, "read1")
  assert_eq(record.flag, 0)
  assert_eq(record.ref_id, 0)
  assert_eq(record.pos, 100)
  assert_eq(record.mapq, 60)
  assert_eq(record.l_seq, 8)
}

///|
test "bam_cigar_parse" {
  let bam = @src.parse_bam_base64(test_bam_base64)

  let record = bam.records[0]
  assert_eq(record.cigar.length(), 1)
  assert_eq(record.cigar[0].length, 8)
}

///|
test "bam_seq_decode" {
  let bam = @src.parse_bam_base64(test_bam_base64)

  let record = bam.records[0]
  assert_eq(record.seq, "AAAAAAAA")
}

///|
test "bam_qual_decode" {
  let bam = @src.parse_bam_base64(test_bam_base64)

  let record = bam.records[0]
  assert_eq(record.qual, "!!!!!!!!")
}

///|
test "bam_flag_checks" {
  let bam = @src.parse_bam_base64(test_bam_base64)

  let record = bam.records[0]
  assert_true(record.is_paired() == false)
  assert_true(record.is_proper_pair() == false)
  assert_true(record.is_unmapped() == false)
  assert_true(record.is_reverse() == false)
}

///|
test "bam_get_reference_name" {
  let bam = @src.parse_bam_base64(test_bam_base64)

  let record = bam.records[0]
  let ref_name = record.get_reference_name(bam.header)
  assert_eq(ref_name, "chr1")
}

///|
test "bam_alignment_stats" {
  let bam = @src.parse_bam_base64(test_bam_base64)

  let record = bam.records[0]
  assert_eq(record.get_end(), 107)
  assert_eq(record.get_alignment_length(), 8)
}

///|
test "bam_to_sam_string" {
  let bam = @src.parse_bam_base64(test_bam_base64)

  let record = bam.records[0]
  let sam_str = record.to_sam_string(bam.header)
  assert_true(sam_str.find("read1").is_some())
  assert_true(sam_str.find("chr1").is_some())
  assert_true(sam_str.find("8M").is_some())
}

///|
test "bam_to_sam" {
  let bam = @src.parse_bam_base64(test_bam_base64)

  let sam_str = bam.to_sam()
  assert_true(sam_str.find("@HD").is_some())
  assert_true(sam_str.find("@SQ").is_some())
  assert_true(sam_str.find("read1").is_some())
}

///|
test "bam_base64_decode" {
  let encoded = "YWJj"
  let decoded = @src.base64_decode(encoded)
  assert_eq(decoded.length(), 3)
  assert_eq(decoded.unsafe_get(0).to_int(), 97)
  assert_eq(decoded.unsafe_get(1).to_int(), 98)
  assert_eq(decoded.unsafe_get(2).to_int(), 99)
}

///|
test "bam_debug_magic" {
  let raw_data = @src.base64_decode(test_bam_base64)
  assert_eq(raw_data.length() >= 4, true)

  let magic_ok = @src.is_bam_magic(raw_data)
  assert_true(magic_ok)
}
