///|
/// Tabix (.tbi) index Demo.
///
/// Demonstrates the Rsamtools / htslib Tabix index port: parsing a `.tbi`
/// file (header fields + sequence names + per-reference binning) and querying
/// a half-open region for its compressed-block virtual offsets. The index data
/// is a canned (offline) fixture.

///|
fn main raise {
  println("=== 1. Parse a Tabix (.tbi) index ===")
  let idx = @src.parse_tabix(
    @src.base64_decode(
      "VEJJAQIAAAACAAAAAQAAAAIAAAAAAAAAIwAAAAAAAAAKAAAAY2hyMQBjaHIyAAEAAAAAAAAAAQAAADIAAAAAAAAAPAAAAAAAAAABAAAAMgAAAAAAAAABAAAAAAAAAAEAAABGAAAAAAAAAFAAAAAAAAAAAQAAAEYAAAAAAAAA",
    ),
  )
  println("n_ref = \{idx.n_ref}, format = \{idx.format}, meta = \{idx.meta}, skip = \{idx.skip}")
  println(
    "columns: seq=\{idx.col_seq}, beg=\{idx.col_beg}, end=\{idx.col_end}",
  )
  println("names = \{idx.names[0]}, \{idx.names[1]}")

  println("\n=== 2. Query a region on ref 0 ===")
  let q = @src.tabix_query(idx, 0, 0, 100)
  println("min_off = \{q.min_off}")
  for j = 0; j < q.chunks.length(); j = j + 1 {
    let c = q.chunks[j]
    println("  chunk [\{c.chunk_beg}, \{c.chunk_end})")
  }

  println("\n=== 3. Negative reference id is rejected ===")
  let err = try {
    ignore(@src.tabix_query(idx, -1, 0, 100))
    ""
  } catch {
    TabixError(msg) => msg
  }
  println("error: \{err}")
}