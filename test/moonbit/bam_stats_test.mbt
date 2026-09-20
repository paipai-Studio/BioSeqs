///|
/// Tests for samtools flagstat + idxstats, cross-validated against the
/// authoritative semantics in samtools `bam_stat.c` / `bam_index.c`.

///|
let flagstat_bam_base64 = "QkFNAS0AAABASEQJVk46MS42CVNPOmNvb3JkaW5hdGUKQFNRCVNOOmNocjEJTE46MTAwMAoBAAAABQAAAGNocjEA6AMAADAAAAAAAAAACgAAAAY8SRIBAEMABAAAAAAAAAD/////AAAAAHBhaXIxAEAAAAASSP////8wAAAAAAAAAAoAAAAGPEkSAQCDAAQAAAAAAAAA/////wAAAABwYWlyMgBAAAAAEkj/////MAAAAAAAAAAUAAAABjxJEgEASQAEAAAA//////////8AAAAAc2luZzEAQAAAABJI/////zAAAAAAAAAAFAAAAAY8SRIBAIkABAAAAP//////////AAAAAHNpbmcyAEAAAAASSP////80AAAAAAAAAB4AAAAKPEkSAQAAAQQAAAD//////////wAAAABzZWNvbmRhcnkAQAAAABJI/////y4AAAAAAAAAKAAAAAQ8SRIBAAAEBAAAAP//////////AAAAAGR1cABAAAAAEkj/////MwAAAP//////////CTxIEgEABAAEAAAA//////////8AAAAAdW5tYXBwZWQAQAAAABJI/////zAAAAAAAAAAMgAAAAY8SRIBAAAABAAAAP//////////AAAAAHBsYWluAEAAAAASSP////8="

///|
test "bam_flagstat" {
  let bam = @src.parse_bam_base64(flagstat_bam_base64)
  let s = bam.flagstat()
  assert_eq(s.total, 8)
  assert_eq(s.primary, 7)
  assert_eq(s.secondary, 1)
  assert_eq(s.supplementary, 0)
  assert_eq(s.duplicates, 1)
  assert_eq(s.primary_duplicates, 1)
  assert_eq(s.mapped, 7)
  assert_eq(s.primary_mapped, 6)
  assert_eq(s.paired, 4)
  assert_eq(s.read1, 2)
  assert_eq(s.read2, 2)
  assert_eq(s.properly_paired, 2)
  assert_eq(s.both_mapped, 2)
  assert_eq(s.singletons, 2)
  assert_eq(s.mate_diff_chr, 0)
  assert_eq(s.mate_diff_chr_mapq5, 0)
}

///|
test "bam_idxstats" {
  let bam = @src.parse_bam_base64(flagstat_bam_base64)
  let r = bam.idxstats()
  assert_eq(r.length(), 2)
  assert_eq(r[0].ref_name, "chr1")
  assert_eq(r[0].ref_length, 1000)
  assert_eq(r[0].mapped, 7)
  assert_eq(r[0].unmapped, 0)
  assert_eq(r[1].ref_name, "*")
  assert_eq(r[1].ref_length, 0)
  assert_eq(r[1].mapped, 0)
  assert_eq(r[1].unmapped, 1)
}