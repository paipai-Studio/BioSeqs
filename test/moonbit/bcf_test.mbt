///|
/// Tests for the BCF (binary VCF) parser. The fixture is a real BCF produced
/// by htslib (via pysam), embedded as base64; decoded values were
/// cross-checked against pysam.

///|
fn bcf_fixture() -> String {
  @src.base64_decode(
    "QkNGAgJLAgAAIyNmaWxlZm9ybWF0PVZDRnY0LjIKIyNGSUxURVI9PElEPVBBU1MsRGVzY3JpcHRpb249IkFsbCBmaWx0ZXJzIHBhc3NlZCIsSURYPTA+CiMjY29udGlnPTxJRD1jaHIxLGxlbmd0aD0xMDAwLElEWD0wPgojI2NvbnRpZz08SUQ9Y2hyMixsZW5ndGg9NTAwLElEWD0xPgojI0ZJTFRFUj08SUQ9cTEwLERlc2NyaXB0aW9uPSJRdWFsaXR5ID4gMTAiLElEWD0xPgojI0lORk89PElEPURQLE51bWJlcj0xLFR5cGU9SW50ZWdlcixEZXNjcmlwdGlvbj0iRGVwdGgiLElEWD0yPgojI0lORk89PElEPUFGLE51bWJlcj1BLFR5cGU9RmxvYXQsRGVzY3JpcHRpb249IkFsbGVsZSBmcmVxIixJRFg9Mz4KIyNJTkZPPTxJRD1GTEFHWCxOdW1iZXI9MCxUeXBlPUZsYWcsRGVzY3JpcHRpb249IkEgZmxhZyIsSURYPTQ+CiMjRk9STUFUPTxJRD1HVCxOdW1iZXI9MSxUeXBlPVN0cmluZyxEZXNjcmlwdGlvbj0iR2Vub3R5cGUiLElEWD01PgojI0ZPUk1BVD08SUQ9QUQsTnVtYmVyPVIsVHlwZT1JbnRlZ2VyLERlc2NyaXB0aW9uPSJBbGxlbGUgZGVwdGhzIixJRFg9Nj4KI0NIUk9NCVBPUwlJRAlSRUYJQUxUCVFVQUwJRklMVEVSCUlORk8JRk9STUFUCVMxCgAwAAAACgAAAAAAAABlAAAAAQAAAAAA9EEDAAIAAQAAAjdyczEXQRdDEQARAhEMEQMVAACAPhEEABEFIQIEEQYhCgIwAAAACwAAAAEAAADHAAAAAQAAAAEAgH8CAAMAAQAAAgcXRxdUF0ERARECEQURAyXNzMw9zcxMPhEFIQQGEQYxAAMC",
  )
}

///|
test "parse_bcf_header" {
  let h = @src.parse_bcf_header(bcf_fixture()) catch {
    _ => abort("BCF header should parse")
  }
  assert_eq(h.contig_names.length(), 2)
  assert_eq(h.contig_names[0], "chr1")
  assert_eq(h.contig_names[1], "chr2")
  assert_eq(h.contig_lengths[0].unwrap_or(-1), 1000)
  assert_eq(h.id_names.length(), 7)
  assert_eq(h.id_names[0], "PASS")
  assert_eq(h.id_names[1], "q10")
  assert_eq(h.id_names[2], "DP")
  assert_eq(h.id_names[3], "AF")
  assert_eq(h.id_names[4], "FLAGX")
  assert_eq(h.id_names[5], "GT")
  assert_eq(h.id_names[6], "AD")
  assert_eq(h.filter_ids, [0, 1])
  assert_eq(h.info_ids, [2, 3, 4])
  assert_eq(h.format_ids, [5, 6])
  assert_eq(h.sample_names.length(), 1)
  assert_eq(h.sample_names[0], "S1")
}

///|
test "parse_bcf_records" {
  let f = @src.parse_bcf(bcf_fixture()) catch {
    _ => abort("BCF should parse")
  }
  assert_eq(f.records.length(), 2)

  let r = f.records[0]
  assert_eq(r.rid, 0)
  assert_eq(r.chrom, "chr1")
  assert_eq(r.pos, 101)
  assert_eq(r.rlen, 1)
  assert_eq(r.id, "rs1")
  assert_eq(r.ref_seq, "A")
  assert_eq(r.alt.length(), 1)
  assert_eq(r.alt[0], "C")
  match r.qual {
    Some(q) => assert_true((q - 30.5).abs() < 1.0e-4)
    None => abort("expected qual")
  }
  assert_eq(r.filter_ids, [0])
  // INFO: DP (12), AF (0.25), FLAGX (flag)
  assert_eq(r.info.length(), 3)
  assert_eq(r.info[0].key_id, 2)
  match r.info[0].values[0] {
    BcfInt(v) => assert_eq(v, 12)
    _ => abort("expected DP int")
  }
  assert_eq(r.info[1].key_id, 3)
  match r.info[1].values[0] {
    BcfFloat(x) => assert_true((x - 0.25).abs() < 1.0e-6)
    _ => abort("expected AF float")
  }
  assert_eq(r.info[2].key_id, 4)
  match r.info[2].values[0] {
    BcfFlag => ()
    _ => abort("expected FLAGX flag")
  }
  // FORMAT: GT [0/1] -> [2, 4]; AD [10, 2]
  assert_eq(r.formats.length(), 2)
  assert_eq(r.formats[0].key_id, 5)
  assert_eq(r.formats[0].values_per_sample, 2)
  match r.formats[0].values[0] {
    BcfInt(v) => {
      assert_eq(v, 2)
      assert_eq(@src.bcf_gt_allele(v), Some(0))
    }
    _ => abort("expected GT int")
    _ => abort("expected GT int")
  }
  assert_eq(r.formats[1].key_id, 6)
  assert_eq(r.formats[1].values_per_sample, 2)
  match (r.formats[1].values[0], r.formats[1].values[1]) {
    (BcfInt(a), BcfInt(b)) => {
      assert_eq(a, 10)
      assert_eq(b, 2)
    }
    _ => abort("expected AD ints")
  }

  // record 2
  let r2 = f.records[1]
  assert_eq(r2.rid, 1)
  assert_eq(r2.chrom, "chr2")
  assert_eq(r2.pos, 199)
  assert_eq(r2.rlen, 1)
  assert_eq(r2.id, "")
  assert_eq(r2.ref_seq, "G")
  assert_eq(r2.alt.length(), 2)
  assert_eq(r2.alt[0], "T")
  assert_eq(r2.alt[1], "A")
  match r2.qual {
    None => ()
    Some(_) => abort("expected missing qual")
  }
  assert_eq(r2.filter_ids, [1])
  assert_eq(r2.info.length(), 2)
  assert_eq(r2.info[0].key_id, 2)
  match r2.info[0].values[0] {
    BcfInt(v) => assert_eq(v, 5)
    _ => abort("expected DP int")
  }
  assert_eq(r2.info[1].key_id, 3)
  match (r2.info[1].values[0], r2.info[1].values[1]) {
    (BcfFloat(x), BcfFloat(y)) => {
      assert_true((x - 0.1).abs() < 1.0e-6)
      assert_true((y - 0.2).abs() < 1.0e-6)
    }
    _ => abort("expected AF floats")
  }
  match r2.formats[0].values[0] {
    BcfInt(v) => {
      assert_eq(v, 4)
      assert_eq(@src.bcf_gt_allele(v), Some(1))
    }
    _ => abort("expected GT int")
  }
  match (r2.formats[1].values[0], r2.formats[1].values[1], r2.formats[1].values[2]) {
    (BcfInt(a), BcfInt(b), BcfInt(c)) => {
      assert_eq(a, 0)
      assert_eq(b, 3)
      assert_eq(c, 2)
    }
    _ => abort("expected AD ints")
  }
}

///|
test "parse_bcf_bad_magic" {
  let r = try {
    ignore(@src.parse_bcf_header("XXXXX"))
    false
  } catch {
    BcfError(msg) => msg == "Invalid BCF2 magic string"
  }
  assert_true(r)
}
///|
fn bcf_val_eq(a : @src.BcfValue, b : @src.BcfValue) -> Bool {
  match (a, b) {
    (@src.BcfFlag, @src.BcfFlag) => true
    (@src.BcfInt(x), @src.BcfInt(y)) => x == y
    (@src.BcfFloat(x), @src.BcfFloat(y)) => assert_close(x, y, 1.0e-6)
    (@src.BcfStr(x), @src.BcfStr(y)) => x == y
    _ => false
  }
}

///|
test "bcf_roundtrip" {
  let f = @src.parse_bcf(bcf_fixture()) catch {
    _ => abort("BCF should parse")
  }
  let bytes = f.to_bcf_bytes()
  let g = @src.parse_bcf(bytes) catch {
    _ => abort("round-tripped BCF should parse")
  }

  assert_eq(g.header.text, f.header.text)
  assert_eq(g.header.contig_names, f.header.contig_names)
  assert_eq(g.header.id_names, f.header.id_names)
  assert_eq(g.header.filter_ids, f.header.filter_ids)
  assert_eq(g.header.info_ids, f.header.info_ids)
  assert_eq(g.header.format_ids, f.header.format_ids)
  assert_eq(g.header.sample_names, f.header.sample_names)
  assert_eq(g.records.length(), f.records.length())

  for i = 0; i < f.records.length(); i = i + 1 {
    let a = f.records[i]
    let b = g.records[i]
    assert_eq(a.rid, b.rid)
    assert_eq(a.pos, b.pos)
    assert_eq(a.rlen, b.rlen)
    assert_eq(a.id, b.id)
    assert_eq(a.ref_seq, b.ref_seq)
    assert_eq(a.alt, b.alt)
    assert_eq(a.filter_ids, b.filter_ids)
    assert_eq(a.n_sample, b.n_sample)
    assert_eq(a.n_fmt, b.n_fmt)
    match (a.qual, b.qual) {
      (None, None) => ()
      (Some(qa), Some(qb)) => assert_true(assert_close(qa, qb, 1.0e-5))
      _ => assert_true(false)
    }
    assert_eq(a.info.length(), b.info.length())
    for j = 0; j < a.info.length(); j = j + 1 {
      assert_eq(a.info[j].key_id, b.info[j].key_id)
      assert_eq(a.info[j].values.length(), b.info[j].values.length())
      for k = 0; k < a.info[j].values.length(); k = k + 1 {
        assert_true(bcf_val_eq(a.info[j].values[k], b.info[j].values[k]))
      }
    }
    assert_eq(a.formats.length(), b.formats.length())
    for j = 0; j < a.formats.length(); j = j + 1 {
      assert_eq(a.formats[j].key_id, b.formats[j].key_id)
      assert_eq(a.formats[j].values_per_sample, b.formats[j].values_per_sample)
      assert_eq(a.formats[j].values.length(), b.formats[j].values.length())
      for k = 0; k < a.formats[j].values.length(); k = k + 1 {
        assert_true(bcf_val_eq(a.formats[j].values[k], b.formats[j].values[k]))
      }
    }
  }
}
