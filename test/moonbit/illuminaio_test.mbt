///|
/// Tests for the illuminaio IDAT raw file parser. The fixture is a synthetic
/// IDAT built with struct.pack following illuminaio's `readIDAT_nonenc`, and
/// the expected values come from an independent Python transcription.

///|
test "idat_parse_header" {
  let raw = @src.base64_decode(@src.idat_sample_base64())
  let r = @src.parse_idat(raw)
  assert_eq(r.version, 3)
  assert_eq(r.n_fields, 7)
  assert_eq(r.n_snps_read, 3)
}

///|
test "idat_parse_quant_fields" {
  let raw = @src.base64_decode(@src.idat_sample_base64())
  let r = @src.parse_idat(raw)
  assert_eq(r.illumina_id.length(), 3)
  assert_eq(r.illumina_id[0], 101)
  assert_eq(r.illumina_id[1], 202)
  assert_eq(r.illumina_id[2], 303)
  assert_eq(r.mean[0], 1000)
  assert_eq(r.mean[1], 2000)
  assert_eq(r.mean[2], 3000)
  assert_eq(r.sd[0], 10)
  assert_eq(r.sd[1], 20)
  assert_eq(r.sd[2], 30)
  assert_eq(r.n_beads[0], 5)
  assert_eq(r.n_beads[1], 6)
  assert_eq(r.n_beads[2], 7)
}

///|
test "idat_parse_metadata" {
  let raw = @src.base64_decode(@src.idat_sample_base64())
  let r = @src.parse_idat(raw)
  assert_eq(r.barcode, "1234567890")
  assert_eq(r.chip_type, "HumanMethylation450")
}

///|
test "idat_bad_magic" {
  let r = try {
    ignore(@src.parse_idat("XXXX along string here"))
    false
  } catch {
    IdatError(_) => true
  }
  assert_true(r)
}