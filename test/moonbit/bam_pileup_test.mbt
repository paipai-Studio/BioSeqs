///|
/// Tests for BAM pileup, cross-validated against pysam `pileup(stepper="all")`.

///|
let pileup_bam_base64 = "QkFNAS0AAABASEQJVk46MS42CVNPOmNvb3JkaW5hdGUKQFNRCVNOOmNocjEJTE46MTAwMAoBAAAABQAAAGNocjEA6AMAADcAAAAAAAAACAAAAAY8SRICAAAABgAAAP//////////AAAAAHJlYWQ1ACQAAABAAAAAEkgS////////MwAAAAAAAAAKAAAABjxJEgEAAAAGAAAA//////////8AAAAAcmVhZDEAYAAAABJIEv///////zAAAAAAAAAACgAAAAY8SRIBAAAABAAAAP//////////AAAAAHJlYWQ2AEAAAAASSP////87AAAAAAAAAAsAAAAGPEkSAwAAAAYAAAD//////////wAAAAByZWFkMwAwAAAAEQAAACAAAAASSBL///////86AAAAAAAAAAwAAAAGPEkSAwAAAAUAAAD//////////wAAAAByZWFkMgAgAAAAEgAAADAAAAASSBD//////zgAAAAAAAAADQAAAAY8SRIDAAAABAAAAP//////////AAAAAHJlYWQ0ACAAAAATAAAAIAAAABJI/////w=="

///|
fn pl_is(pr : @src.PileupRead, qp : Int?, is_del : Bool, is_refskip : Bool, indel : Int) -> Bool {
  pr.query_position == qp &&
  pr.is_del == is_del &&
  pr.is_refskip == is_refskip &&
  pr.indel == indel
}

///|
test "bam_pileup_columns" {
  let bam = @src.parse_bam_base64(pileup_bam_base64)
  let cols = bam.pileup("chr1", 0, 25)
  assert_eq(cols.length(), 10)
  let mut i = 0
  while i < 10 {
    assert_eq(cols[i].reference_pos, 8 + i)
    i = i + 1
  }
}

///|
test "bam_pileup_simple_match" {
  let bam = @src.parse_bam_base64(pileup_bam_base64)
  let cols = bam.pileup("chr1", 0, 25)
  // pos=10: read5(q=4, soft-clipped 2 + 2 matched), read1(q=0), read6(q=0)
  let c = cols[2] // pos=10
  assert_eq(c.nsegments, 3)
  assert_true(pl_is(c.pileups[0], Some(4), false, false, 0)) // read5
  assert_true(pl_is(c.pileups[1], Some(0), false, false, 0)) // read1
  assert_true(pl_is(c.pileups[2], Some(0), false, false, 0)) // read6
}

///|
test "bam_pileup_insertion_deletion" {
  let bam = @src.parse_bam_base64(pileup_bam_base64)
  let cols = bam.pileup("chr1", 0, 25)
  // pos=13: read3 has insertion (indel=+1), read2 has deletion (indel=-1)
  let c = cols[5] // pos=13
  assert_eq(c.nsegments, 5)
  assert_true(pl_is(c.pileups[2], Some(2), false, false, 1))  // read3 (3M1I2M)
  assert_true(pl_is(c.pileups[3], Some(1), false, false, -1)) // read2 (2M1D3M)
  assert_true(pl_is(c.pileups[4], Some(0), false, false, 0))  // read4
  // pos=14: read2 is a deletion; read3 skips the inserted base (q jumps 2->4)
  let d = cols[6] // pos=14
  assert_true(pl_is(d.pileups[0], Some(4), false, false, 0))  // read1
  assert_true(pl_is(d.pileups[1], Some(4), false, false, 0))  // read3 (insertion consumed q=3)
  assert_true(pl_is(d.pileups[2], None, true, false, 0))      // read2 deletion
  assert_true(pl_is(d.pileups[3], Some(1), false, false, 0))  // read4
}

///|
test "bam_pileup_refskip" {
  let bam = @src.parse_bam_base64(pileup_bam_base64)
  let cols = bam.pileup("chr1", 0, 25)
  // pos=15: read4 is a refskip (N) -> is_del=1 AND is_refskip=1 (htslib quirk)
  let d = cols[7] // pos=15
  assert_eq(d.nsegments, 4)
  assert_true(pl_is(d.pileups[0], Some(5), false, false, 0)) // read1
  assert_true(pl_is(d.pileups[1], Some(5), false, false, 0)) // read3
  assert_true(pl_is(d.pileups[2], Some(2), false, false, 0)) // read2 (2M1D3M, after D)
  assert_true(pl_is(d.pileups[3], None, true, true, 0))      // read4 refskip
  // pos=16 and 17: only read2 and read4 remain
  let e = cols[8] // pos=16
  assert_eq(e.nsegments, 2)
  assert_true(pl_is(e.pileups[0], Some(3), false, false, 0)) // read2
  assert_true(pl_is(e.pileups[1], Some(2), false, false, 0)) // read4
}
///|
test "bam_query_accessors" {
  let bam = @src.parse_bam_base64(pileup_bam_base64)
  // read5 (index 0): 2S4M, seq ACGTAC
  let r5 = bam.records[0]
  assert_eq(r5.query_sequence(), "ACGTAC")
  assert_eq(r5.infer_query_length(), 6) // 2S + 4M
  assert_eq(r5.query_length(), 6)
  // read3 (index 3): 3M1I2M
  let r3 = bam.records[3]
  assert_eq(r3.query_sequence(), "ACGTAC")
  assert_eq(r3.infer_query_length(), 6) // 3M + 1I + 2M
  // read2 (index 4): 2M1D3M, D excluded from query length
  let r2 = bam.records[4]
  assert_eq(r2.query_sequence(), "ACGTA")
  assert_eq(r2.infer_query_length(), 5) // 2M + 3M
  // read4 (index 5): 2M1N2M, N excluded
  let r4 = bam.records[5]
  assert_eq(r4.query_sequence(), "ACGT")
  assert_eq(r4.infer_query_length(), 4) // 2M + 2M
}
