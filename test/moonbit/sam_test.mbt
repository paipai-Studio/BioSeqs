///|
/// Test data: a simple SAM file for testing.
/// Each sequence is exactly 80 bases (20 repeats of "ACGT").
let test_sam_content = "@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:1000\n@SQ\tSN:chr2\tLN:2000\n@RG\tID:test\tSM:sample1\n@CO\tThis is a comment\nread1\t0\tchr1\t100\t60\t80M\t*\t0\t0\tACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACG\t!\"#$%&'()*+,-./0123456789:;<=>?@ABCDEFGHIJKLMNOPQRSTUVWXYZ[\\]^_`abcdefghijklmnopqrstuvwxyz{|}~!\tNM:i:0\tMD:Z:80\tAS:i:80\nread2\t16\tchr1\t200\t40\t50M20D30M\t*\t0\t0\tACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACG\t!\"#$%&'()*+,-./0123456789:;<=>?@ABCDEFGHIJKLMNOPQRSTUVWXYZ[\\]^_`abcd\tNM:i:2\tMD:Z:48A20D28C2\nread3\t73\tchr2\t50\t50\t10S60M10S\t=\t100\t150\tCGTACGTACGACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGATGCATGCATG\t\"#$%&'()*+,-./0123456789:;<=>?@ABCDEFGHIJKLMNOPQRSTUVWXYZ[\\]^_`abcdefghijklmnopqrstuvwxyz{|}~!\"#$%&'()*+,-./012\tNM:i:1\nread4\t4\t*\t0\t0\t*\t*\t0\t0\tACGTACGTACGTACGTACG\t!\"#$%&'()*+,-./0123456789:;\n"

///|
test "sam_parse_basic" {
  let sam = @src.parse_sam(test_sam_content)

  assert_eq(sam.num_headers(), 5)
  assert_eq(sam.len(), 4)

  let seq_names = sam.get_sequence_names()
  assert_eq(seq_names.length(), 2)
  assert_true(seq_names[0] == "chr1" || seq_names[1] == "chr1")
  assert_true(seq_names[0] == "chr2" || seq_names[1] == "chr2")

  assert_eq(sam.get_sequence_length("chr1"), Some(1000))
  assert_eq(sam.get_sequence_length("chr2"), Some(2000))
  assert_eq(sam.get_sequence_length("chr3"), None)
}

///|
test "sam_parse_record" {
  let sam = @src.parse_sam(test_sam_content)
  let records = sam.records

  let r1 = records[0]
  assert_eq(r1.qname, "read1")
  assert_eq(r1.flag, 0)
  assert_eq(r1.rname, "chr1")
  assert_eq(r1.pos, 100)
  assert_eq(r1.mapq, 60)

  let r4 = records[3]
  assert_eq(r4.qname, "read4")
  assert_eq(r4.flag, 4)
  assert_eq(r4.rname, "*")
  assert_eq(r4.pos, 0)
  assert_eq(r4.mapq, 0)
}

///|
test "sam_cigar_parse" {
  let cigar1 = @src.parse_cigar("80M")
  assert_eq(cigar1.length(), 1)
  assert_eq(cigar1[0].length, 80)

  let cigar2 = @src.parse_cigar("50M20D30M")
  assert_eq(cigar2.length(), 3)
  assert_eq(cigar2[0].length, 50)
  assert_eq(cigar2[1].length, 20)
  assert_eq(cigar2[2].length, 30)

  let cigar3 = @src.parse_cigar("10S60M10S")
  assert_eq(cigar3.length(), 3)
  assert_eq(cigar3[0].length, 10)
  assert_eq(cigar3[1].length, 60)
  assert_eq(cigar3[2].length, 10)

  let cigar4 = @src.parse_cigar("*")
  assert_eq(cigar4.length(), 0)
}

///|
test "sam_flag_checks" {
  let sam = @src.parse_sam(test_sam_content)
  let records = sam.records

  let r1 = records[0]
  assert_true(r1.is_paired() == false)
  assert_true(r1.is_proper_pair() == false)
  assert_true(r1.is_unmapped() == false)
  assert_true(r1.is_reverse() == false)

  let r2 = records[1]
  assert_true(r2.is_paired() == false)
  assert_true(r2.is_reverse() == true)

  let r3 = records[2]
  assert_true(r3.is_paired() == true)
  assert_true(r3.is_read1() == true)

  let r4 = records[3]
  assert_true(r4.is_unmapped() == true)
}

///|
test "sam_alignment_stats" {
  let sam = @src.parse_sam(test_sam_content)
  let records = sam.records

  let r1 = records[0]
  assert_eq(r1.get_end(), 179)
  assert_eq(r1.get_alignment_length(), 80)

  let r2 = records[1]
  assert_eq(r2.get_end(), 299)
  assert_eq(r2.get_alignment_length(), 80)

  let r3 = records[2]
  assert_eq(r3.get_end(), 109)
  assert_eq(r3.get_alignment_length(), 60)
}

///|
test "sam_tags" {
  let sam = @src.parse_sam(test_sam_content)
  let records = sam.records

  let r1 = records[0]
  assert_eq(r1.tags.get("NM"), Some("0"))
  assert_eq(r1.tags.get("MD"), Some("80"))
  assert_eq(r1.tags.get("AS"), Some("80"))

  let r2 = records[1]
  assert_eq(r2.tags.get("NM"), Some("2"))
  assert_eq(r2.tags.get("MD"), Some("48A20D28C2"))
}
