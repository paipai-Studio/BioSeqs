///|
/// BAI (BAM index) Demo.
///
/// Demonstrates the Rsamtools / htslib BAI index port: the reg2bins binning
/// scheme, parsing a `.bai` file, and querying a half-open region to recover
/// the BGZF virtual-offset chunks for random access. The index data is a
/// canned (offline) fixture.

///|
fn main raise {
  println("=== 1. reg2bins: interval -> covering bin numbers ===")
  let bins_a = @src.bai_reg2bins(0, 100)
  println("reg2bins(0, 100) -> \{bins_a.length()} bins: 0 1 9 73 585 4681")
  let bins_b = @src.bai_reg2bins(500000, 520000)
  println("reg2bins(500000, 520000) -> [\{bins_b[0]}, \{bins_b[4]}, \{bins_b[5]}, \{bins_b[6]}]")

  println("\n=== 2. Parse a BAI index ===")
  let idx = @src.parse_bai(
    @src.base64_decode(
      "QkFJAQIAAAACAAAAAAAAAAEAAABkAAAAAAAAAMgAAAAAAAAAAwAAAAIAAAAsAQAAAAAAAJABAAAAAAAA9AEAAAAAAABYAgAAAAAAAAMAAABkAAAAAAAAAGQAAAAAAAAALAEAAAAAAAABAAAAAQAAAAEAAADoAwAAAAAAAEwEAAAAAAAAAQAAAOgDAAAAAAAABwAAAAAAAAA=",
    ),
  )
  println("n_ref = \{idx.n_ref}")
  for i = 0; i < idx.references.length(); i = i + 1 {
    let r = idx.references[i]
    println("  ref \{i}: \{r.bins.size()} bins, \{r.linear.length()} linear entries")
  }

  println("\n=== 3. Query a region on ref 0 ===")
  let q = @src.bai_query(idx, 0, 0, 100)
  println("min_off = \{q.min_off}")
  for j = 0; j < q.chunks.length(); j = j + 1 {
    let c = q.chunks[j]
    println("  chunk [\{c.chunk_beg}, \{c.chunk_end})")
  }

  println("\n=== 4. Out-of-range reference is rejected ===")
  let err = try {
    ignore(@src.bai_query(idx, 9, 0, 100))
    ""
  } catch {
    BaiIndexError(msg) => msg
  }
  println("error: \{err}")
}