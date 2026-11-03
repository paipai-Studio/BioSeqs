// ===========================================================================
// airr_test.mbt — AIRR TSV 读写测试
// ===========================================================================

///|
fn airr_test_record(id : String, cdr3 : String, aa : String) -> @src.AirrRecord {
  {
    sequence_id: id,
    locus: "TRB",
    v_call: "TRBV01",
    d_call: "TRBD1",
    j_call: "TRBJ02",
    cdr3,
    cdr3_aa: aa,
    productive: true,
    v_identity: 0.95,
    j_identity: 1.0,
    junction_length: cdr3.length(),
  }
}

///|
test "airr write header and rows" {
  let records : Array[@src.AirrRecord] = [
    airr_test_record("r1", "TGTACGTGG", "CTW"),
    airr_test_record("r2", "TGTTTATGG", "CLW"),
  ]
  let text = @src.airr_write(records)
  let lines = text.split("\n").to_array()
  assert_eq(
    lines[0].to_owned(),
    "sequence_id\tlocus\tv_call\td_call\tj_call\tcdr3\tcdr3_aa\tproductive\tv_identity\tj_identity\tjunction_length",
  )
  assert_eq(
    lines[1].to_owned(),
    "r1\tTRB\tTRBV01\tTRBD1\tTRBJ02\tTGTACGTGG\tCTW\tT\t0.95\t1.0\t9",
  )
  assert_eq(
    lines[2].to_owned(),
    "r2\tTRB\tTRBV01\tTRBD1\tTRBJ02\tTGTTTATGG\tCLW\tT\t0.95\t1.0\t9",
  )
}

///|
test "airr write parse roundtrip" {
  let records : Array[@src.AirrRecord] = [
    airr_test_record("r1", "TGTACGTGG", "CTW"),
    airr_test_record("r2", "TGTTTATGG", "CLW"),
  ]
  let text = @src.airr_write(records)
  let parsed = @src.airr_parse(text)
  assert_eq(parsed.length(), 2)
  assert_eq(parsed[0].sequence_id, "r1")
  assert_eq(parsed[0].locus, "TRB")
  assert_eq(parsed[0].v_call, "TRBV01")
  assert_eq(parsed[0].d_call, "TRBD1")
  assert_eq(parsed[0].j_call, "TRBJ02")
  assert_eq(parsed[0].cdr3, "TGTACGTGG")
  assert_eq(parsed[0].cdr3_aa, "CTW")
  assert_true(parsed[0].productive)
  assert_true(parsed[0].v_identity > 0.94 && parsed[0].v_identity < 0.96)
  assert_true(parsed[0].j_identity > 0.99)
  assert_eq(parsed[0].junction_length, 9)
  assert_eq(parsed[1].sequence_id, "r2")
}

///|
test "airr parse tolerant to missing columns" {
  let text = "sequence_id\tv_call\tcdr3\nr1\tTRBV01\tTGTACGTGG\n"
  let parsed = @src.airr_parse(text)
  assert_eq(parsed.length(), 1)
  assert_eq(parsed[0].v_call, "TRBV01")
  assert_eq(parsed[0].cdr3, "TGTACGTGG")
  assert_eq(parsed[0].j_call, "")
  assert_false(parsed[0].productive)
  assert_eq(parsed[0].junction_length, 9)
}

///|
test "airr parse column order arbitrary" {
  let text = "cdr3_aa\tsequence_id\tproductive\nCTW\tr1\tT\n"
  let parsed = @src.airr_parse(text)
  assert_eq(parsed[0].sequence_id, "r1")
  assert_eq(parsed[0].cdr3_aa, "CTW")
  assert_true(parsed[0].productive)
}

///|
test "airr parse missing sequence_id raises" {
  let text = "v_call\tcdr3\nTRBV01\tTGTACGTGG\n"
  let mut raised = false
  try {
    let parsed = @src.airr_parse(text)
    ignore(parsed)
  } catch {
    _ => raised = true
  }
  assert_true(raised)
}

///|
test "airr parse empty and blank lines" {
  let text = "sequence_id\tcdr3\n\nr1\tTGTACGTGG\n\n"
  let parsed = @src.airr_parse(text)
  assert_eq(parsed.length(), 1)
  assert_eq(parsed[0].sequence_id, "r1")
}

///|
test "airr from vdj chain conversion" {
  let db = @src.vdj_default_trb_db()
  ignore(db)
  let chain : @src.VdjChain = {
    sequence_id: "c1",
    v_call: "TRBV05",
    d_call: "TRBD2",
    j_call: "TRBJ01",
    v_identity: 0.9,
    j_identity: 0.88,
    cdr3_nt: "TGTACGTGG",
    cdr3_aa: "CTW",
    in_frame: true,
    stop_codon: false,
    productive: true,
  }
  let rec = @src.AirrRecord::from_vdj_chain(chain, locus="TRB")
  assert_eq(rec.sequence_id, "c1")
  assert_eq(rec.locus, "TRB")
  assert_eq(rec.v_call, "TRBV05")
  assert_eq(rec.cdr3, "TGTACGTGG")
  assert_eq(rec.cdr3_aa, "CTW")
  assert_true(rec.productive)
  assert_eq(rec.junction_length, 9)
}
