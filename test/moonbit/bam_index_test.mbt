///|
/// Tests for the BAI (BAM index) parser and region query. The fixture is an
/// independently generated index (Python struct.pack), embedded as base64, so
/// the parser is checked against a source of truth it did not produce.

///|
fn bai_fixture() -> String {
  @src.base64_decode(
    "QkFJAQIAAAACAAAAAAAAAAEAAABkAAAAAAAAAMgAAAAAAAAAAwAAAAIAAAAsAQAAAAAAAJABAAAAAAAA9AEAAAAAAABYAgAAAAAAAAMAAABkAAAAAAAAAGQAAAAAAAAALAEAAAAAAAABAAAAAQAAAAEAAADoAwAAAAAAAEwEAAAAAAAAAQAAAOgDAAAAAAAABwAAAAAAAAA=",
  )
}

///|
test "bai_reg2bins_basic" {
  let bins = @src.bai_reg2bins(0, 100)
  assert_eq(bins.length(), 6)
  assert_eq(bins[0], 0)
  assert_eq(bins[1], 1)
  assert_eq(bins[2], 9)
  assert_eq(bins[3], 73)
  assert_eq(bins[4], 585)
  assert_eq(bins[5], 4681)
}

///|
test "bai_reg2bins_multi_level" {
  let bins = @src.bai_reg2bins(500000, 520000)
  assert_eq(bins.length(), 7)
  assert_eq(bins[0], 0)
  assert_eq(bins[1], 1)
  assert_eq(bins[2], 9)
  assert_eq(bins[3], 73)
  assert_eq(bins[4], 588)
  assert_eq(bins[5], 4711)
  assert_eq(bins[6], 4712)
}

///|
test "bai_reg2bins_empty" {
  assert_eq(@src.bai_reg2bins(10, 10).length(), 0)
  assert_eq(@src.bai_reg2bins(5, 5).length(), 0)
  assert_eq(@src.bai_reg2bins(7, 3).length(), 0)
}

///|
test "parse_bai_structure" {
  let index = @src.parse_bai(bai_fixture()) catch {
    _ => abort("BAI fixture should parse")
  }
  assert_eq(index.n_ref, 2)
  assert_eq(index.references.length(), 2)

  let r0 = index.references[0]
  assert_eq(r0.linear.length(), 3)
  assert_eq(r0.linear[0].to_int(), 100)
  assert_eq(r0.linear[2].to_int(), 300)
  let bin0 = r0.bins.get(0).unwrap_or([])
  assert_eq(bin0.length(), 1)
  assert_eq(bin0[0].chunk_beg.to_int(), 100)
  assert_eq(bin0[0].chunk_end.to_int(), 200)
  let bin3 = r0.bins.get(3).unwrap_or([])
  assert_eq(bin3.length(), 2)
  assert_eq(bin3[1].chunk_beg.to_int(), 500)
  assert_eq(bin3[1].chunk_end.to_int(), 600)

  let r1 = index.references[1]
  assert_eq(r1.linear.length(), 1)
  assert_eq(r1.linear[0].to_int(), 1000)
  let bin1 = r1.bins.get(1).unwrap_or([])
  assert_eq(bin1[0].chunk_beg.to_int(), 1000)
  assert_eq(bin1[0].chunk_end.to_int(), 1100)

  match index.n_no_coor {
    Some(n) => assert_eq(n.to_int(), 7)
    None => abort("expected n_no_coor")
  }
}

///|
test "bai_query_returns_overlapping_chunks" {
  let index = @src.parse_bai(bai_fixture()) catch {
    _ => abort("BAI fixture should parse")
  }
  let q = @src.bai_query(index, 0, 0, 100)
  assert_eq(q.min_off.to_int(), 100)
  assert_eq(q.chunks.length(), 1)
  assert_eq(q.chunks[0].chunk_beg.to_int(), 100)
  assert_eq(q.chunks[0].chunk_end.to_int(), 200)

  let q2 = @src.bai_query(index, 1, 0, 100)
  assert_eq(q2.min_off.to_int(), 1000)
  assert_eq(q2.chunks.length(), 1)
  assert_eq(q2.chunks[0].chunk_beg.to_int(), 1000)
}

///|
test "bai_query_invalid" {
  let index = @src.parse_bai(bai_fixture()) catch {
    _ => abort("BAI fixture should parse")
  }
  let r1 = try {
    ignore(@src.bai_query(index, 9, 0, 100))
    false
  } catch {
    BaiIndexError(_) => true
  }
  assert_true(r1)
  let r2 = try {
    ignore(@src.bai_query(index, 0, 100, 50))
    false
  } catch {
    BaiIndexError(_) => true
  }
  assert_true(r2)
}

///|
test "parse_bai_bad_magic" {
  let r = try {
    ignore(@src.parse_bai("XXXX"))
    false
  } catch {
    BaiIndexError(msg) => msg == "Not a BAI index: bad magic bytes"
  }
  assert_true(r)
}