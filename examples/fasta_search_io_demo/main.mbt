///|
/// fasta_search_io module demo - FASTA -m8/-m9 search output parsing.
///
/// Demonstrates parsing of `fasta36` `-m 8` (compact tabular) and
/// `-m 9` (tabular with comment header) output, including metadata
/// extraction and query/hit/HSP aggregation.
///
/// Reference: Biopython Bio.SearchIO.FastaIO
/// Run: moon run examples/fasta_search_io_demo
fn main {
  println("=== fasta_search_io Module Demo ===")
  println("FASTA -m8/-m9 search output parsing")
  println("Reference: Biopython Bio.SearchIO.FastaIO")
  println("")

  // ------------------------------------------------------------------
  // 1. Parse -m 8 (compact tabular) sample
  // ------------------------------------------------------------------
  println("=== 1. Parse -m 8 Sample ===")
  let m8 = @bio.fasta_m8_sample()
  let m8_results = @bio.parse_fasta_m8(m8)
  println("Number of queries: " + m8_results.length().to_string())
  let mut qi = 0
  while qi < m8_results.length() {
    let qr = m8_results[qi]
    println("Query: " + qr.id + "  hits: " + qr.hits.length().to_string())
    let mut hi = 0
    while hi < qr.hits.length() {
      let h = qr.hits[hi]
      println(
        "  Hit: " + h.id +
        "  bitscore=" + h.bitscore.to_string() +
        "  evalue=" + h.evalue.to_string() +
        "  hsps=" + h.hsps.length().to_string(),
      )
      let mut si = 0
      while si < h.hsps.length() {
        let frag = h.hsps[si].fragments[0]
        println(
          "    HSP: query=" + frag.query_start.to_string() + "-" + frag.query_end.to_string() +
          "  subject=" + frag.hit_start.to_string() + "-" + frag.hit_end.to_string() +
          "  aln_len=" + frag.aln_span.to_string(),
        )
        si = si + 1
      }
      hi = hi + 1
    }
    qi = qi + 1
  }
  println("")

  // ------------------------------------------------------------------
  // 2. Parse -m 9 (tabular with comments) sample
  // ------------------------------------------------------------------
  println("=== 2. Parse -m 9 Sample ===")
  let m9 = @bio.fasta_m9_sample()
  let m9_result = @bio.parse_fasta_m9(m9)
  let meta = m9_result.meta
  println("Metadata:")
  println("  Program: " + meta.program)
  println("  Version: " + meta.version)
  println("  Database: " + meta.database)
  println("  Comment lines: " + meta.comments.length().to_string())
  println("Queries parsed: " + m9_result.queries.length().to_string())
  println("")

  // ------------------------------------------------------------------
  // 3. Custom -m 8 input
  // ------------------------------------------------------------------
  println("=== 3. Custom -m 8 Input ===")
  let custom =
    "q1\ts1\t95.0\t100\t5\t0\t1\t100\t1\t100\t1e-50\t200.3\n" +
    "q1\ts2\t80.0\t80\t16\t2\t1\t80\t1\t80\t1e-30\t150.0\n"
  let custom_results = @bio.parse_fasta_m8(custom)
  println("Queries: " + custom_results.length().to_string())
  println(
    "First query hits: " + custom_results[0].hits.length().to_string(),
  )
  println(
    "First hit id: " + custom_results[0].hits[0].id +
    "  bitscore=" + custom_results[0].hits[0].bitscore.to_string(),
  )
  println("")

  println("=== Demo Complete ===")
}
