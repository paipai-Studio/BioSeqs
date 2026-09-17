///|
/// Tests for the Tabix (.tbi) index parser and region query. The fixture is
/// an independently generated index (Python struct.pack), embedded as base64.

///|
fn tbi_fixture() -> String {
  @src.base64_decode(
    "VEJJAQIAAAACAAAAAQAAAAIAAAAAAAAAIwAAAAAAAAAKAAAAY2hyMQBjaHIyAAEAAAAAAAAAAQAAADIAAAAAAAAAPAAAAAAAAAABAAAAMgAAAAAAAAABAAAAAAAAAAEAAABGAAAAAAAAAFAAAAAAAAAAAQAAAEYAAAAAAAAA",
  )
}

///|
test "parse_tabix_structure" {
  let index = @src.parse_tabix(tbi_fixture()) catch {
    _ => abort("TBI fixture should parse")
  }
  assert_eq(index.n_ref, 2)
  assert_eq(index.format, 2)
  assert_eq(index.col_seq, 1)
  assert_eq(index.col_beg, 2)
  assert_eq(index.col_end, 0)
  assert_eq(index.meta, 35)
  assert_eq(index.skip, 0)
  assert_eq(index.names.length(), 2)
  assert_eq(index.names[0], "chr1")
  assert_eq(index.names[1], "chr2")

  let r0 = index.references[0]
  assert_eq(r0.linear.length(), 1)
  assert_eq(r0.linear[0].to_int(), 50)
  let bin0 = r0.bins.get(0).unwrap_or([])
  assert_eq(bin0.length(), 1)
  assert_eq(bin0[0].chunk_beg.to_int(), 50)
  assert_eq(bin0[0].chunk_end.to_int(), 60)

  let r1 = index.references[1]
  assert_eq(r1.linear[0].to_int(), 70)
  let bin0b = r1.bins.get(0).unwrap_or([])
  assert_eq(bin0b[0].chunk_beg.to_int(), 70)
  assert_eq(bin0b[0].chunk_end.to_int(), 80)
}

///|
test "tabix_query_returns_overlapping_chunks" {
  let index = @src.parse_tabix(tbi_fixture()) catch {
    _ => abort("TBI fixture should parse")
  }
  let q = @src.tabix_query(index, 0, 0, 100)
  assert_eq(q.min_off.to_int(), 50)
  assert_eq(q.chunks.length(), 1)
  assert_eq(q.chunks[0].chunk_beg.to_int(), 50)
  assert_eq(q.chunks[0].chunk_end.to_int(), 60)
}

///|
test "tabix_query_invalid" {
  let index = @src.parse_tabix(tbi_fixture()) catch {
    _ => abort("TBI fixture should parse")
  }
  let r = try {
    ignore(@src.tabix_query(index, -1, 0, 100))
    false
  } catch {
    TabixError(_) => true
  }
  assert_true(r)
}

///|
test "parse_tabix_bad_magic" {
  let r = try {
    ignore(@src.parse_tabix("XXXX"))
    false
  } catch {
    TabixError(msg) => msg == "Not a Tabix index: bad magic bytes"
  }
  assert_true(r)
}