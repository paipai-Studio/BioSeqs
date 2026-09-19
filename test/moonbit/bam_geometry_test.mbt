///|
/// Tests for BAM record geometry methods (reference_end / get_reference_positions
/// / get_overlap / get_aligned_pairs), cross-validated against pysam
/// `AlignedSegment`.

///|
let bam_geom_base64 = "QkFNAS0AAABASEQJVk46MS42CVNPOmNvb3JkaW5hdGUKQFNRCVNOOmNocjEJTE46MTAwMAoBAAAABQAAAGNocjEA6AMAAC0AAAAAAAAACgAAAAM8SRIBAAAABAAAAP//////////AAAAAHIxAEAAAAASSP////87AAAAAAAAAAoAAAADPEkSAwAAAAgAAAD//////////wAAAAByMgBAAAAAIQAAACAAAAASSBJI//////////83AAAAAAAAAAoAAAADPEkSAwAAAAUAAAD//////////wAAAAByMwAwAAAAEgAAACAAAAASSED//////zwAAAAAAAAACgAAAAM8SRIEAAAABgAAAP//////////AAAAAHI0ACAAAAARAAAAEAAAACQAAAASSBL///////8+AAAAAAAAAAoAAAADPEkSBAAAAAcAAAD//////////wAAAAByNQAkAAAAMAAAABMAAAAgAAAAEkgSgP////////9AAAAAAAAAAAoAAAADPEkSBQAAAAYAAAD//////////wAAAAByNgAUAAAAIAAAABIAAAAQAAAAJAAAABJIEv///////zsAAAAAAAAACgAAAAM8SRIDAAAACAAAAP//////////AAAAAHI3AFcAAAAYAAAAIAAAABJIEkj//////////w=="

///|
fn opt_ints(a : Array[Int?], b : Array[Int?]) -> Bool {
  if a.length() != b.length() {
    return false
  }
  for i = 0; i < a.length(); i = i + 1 {
    if a[i] != b[i] {
      return false
    }
  }
  true
}

///|
test "bam_reference_end" {
  let bam = @src.parse_bam_base64(bam_geom_base64)
  assert_eq(bam.records[0].reference_end(), 14) // 4M
  assert_eq(bam.records[1].reference_end(), 16) // 4M2I2M
  assert_eq(bam.records[2].reference_end(), 16) // 3M1D2M
  assert_eq(bam.records[3].reference_end(), 13) // 2M1I1M2S
  assert_eq(bam.records[4].reference_end(), 16) // 2S3M1N2M
  assert_eq(bam.records[5].reference_end(), 14) // 1S2M1D1M2S
  assert_eq(bam.records[6].reference_end(), 18) // 5=1X2M
  assert_eq(bam.records[0].reference_start(), 10)
}

///|
test "bam_get_reference_positions" {
  let bam = @src.parse_bam_base64(bam_geom_base64)
  // 4M2I2M: insertions skipped by default, None placeholders with full_length
  assert_true(opt_ints(bam.records[1].get_reference_positions(false), [Some(10), Some(11), Some(12), Some(13), Some(14), Some(15)]))
  assert_true(opt_ints(bam.records[1].get_reference_positions(true), [Some(10), Some(11), Some(12), Some(13), None, None, Some(14), Some(15)]))
  // 3M1D2M: deletion position excluded
  assert_true(opt_ints(bam.records[2].get_reference_positions(false), [Some(10), Some(11), Some(12), Some(14), Some(15)]))
  // 2S3M1N2M: soft clips None, skip excluded
  assert_true(opt_ints(bam.records[4].get_reference_positions(true), [None, None, Some(10), Some(11), Some(12), Some(14), Some(15)]))
  // 5=1X2M: = and X emit positions
  assert_true(opt_ints(bam.records[6].get_reference_positions(false), [Some(10), Some(11), Some(12), Some(13), Some(14), Some(15), Some(16), Some(17)]))
}

///|
test "bam_get_overlap" {
  let bam = @src.parse_bam_base64(bam_geom_base64)
  assert_eq(bam.records[0].get_overlap(8, 20), 4) // 4M within [8,20)
  assert_eq(bam.records[0].get_overlap(12, 14), 2)
  assert_eq(bam.records[0].get_overlap(20, 30), 0)
  assert_eq(bam.records[1].get_overlap(8, 20), 6) // 4M2I2M
  assert_eq(bam.records[2].get_overlap(8, 20), 5) // deletion not counted
  assert_eq(bam.records[4].get_overlap(8, 20), 5) // skip not counted
}

///|
test "bam_get_aligned_pairs" {
  let bam = @src.parse_bam_base64(bam_geom_base64)
  // 4M2I2M: insertions as (query, None); matches_only filters them out
  let ap = bam.records[1].get_aligned_pairs(false)
  assert_eq(ap.length(), 8)
  assert_eq(ap[0], (Some(0), Some(10)))
  assert_eq(ap[4], (Some(4), None))
  assert_eq(ap[6], (Some(6), Some(14)))
  let apm = bam.records[1].get_aligned_pairs(true)
  assert_eq(apm.length(), 6)
  assert_eq(apm[4], (Some(6), Some(14)))
  // 3M1D2M: deletion as (None, ref)
  let ap2 = bam.records[2].get_aligned_pairs(false)
  assert_eq(ap2[3], (None, Some(13)))
  assert_eq(ap2[4], (Some(3), Some(14)))
  // 1S2M1D1M2S: soft clips as (query, None)
  let ap6 = bam.records[5].get_aligned_pairs(false)
  assert_eq(ap6[0], (Some(0), None))
  assert_eq(ap6[3], (None, Some(12)))
}