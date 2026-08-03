///|
/// hmmer_io module demo - HMMER3 domtblout and text output parsing.
///
/// Demonstrates parsing of HMMER3 `--domtblout` (domain table) and
/// human-readable text output from `hmmscan`/`hmmsearch`, including
/// query/hit/HSP/domain extraction and bit-score/E-value reporting.
///
/// Reference: Biopython Bio.SearchIO.HmmerIO
/// Run: moon run examples/hmmer_io_demo
fn main {
  println("=== hmmer_io Module Demo ===")
  println("HMMER3 domtblout and text output parsing")
  println("Reference: Biopython Bio.SearchIO.HmmerIO")
  println("")

  // ------------------------------------------------------------------
  // 1. Parse domtblout (domain table) sample
  // ------------------------------------------------------------------
  println("=== 1. Parse domtblout Sample ===")
  let domtab = @bio.hmmer3_domtab_sample()
  let domtab_results = @bio.parse_hmmer3_domtab(domtab)
  println("Number of query results: " + domtab_results.length().to_string())
  let mut qi = 0
  while qi < domtab_results.length() {
    let qr = domtab_results[qi]
    println("Query: " + qr.id + " (program=" + qr.program + ")")
    println("  Hits: " + qr.hits.length().to_string())
    println("  Total HSPs: " + @bio.hmmer_io_count_hsps(qr).to_string())
    let mut hi = 0
    while hi < qr.hits.length() {
      let h = qr.hits[hi]
      println(
        "    Hit: " + h.id + "  bitscore=" + h.bitscore.to_string() +
        "  evalue=" + h.evalue.to_string() +
        "  hsps=" + h.hsps.length().to_string(),
      )
      hi = hi + 1
    }
    qi = qi + 1
  }
  println("")

  // ------------------------------------------------------------------
  // 2. Parse human-readable text output sample
  // ------------------------------------------------------------------
  println("=== 2. Parse Text Output Sample ===")
  let text = @bio.hmmer3_text_sample()
  let text_results = @bio.parse_hmmer3_text(text)
  println("Number of query results: " + text_results.length().to_string())
  let mut ti = 0
  while ti < text_results.length() {
    let qr = text_results[ti]
    println(
      "Query: " + qr.id + "  length=" + qr.seq_len.to_string() +
      "  hits=" + qr.hits.length().to_string(),
    )
    let mut hi = 0
    while hi < qr.hits.length() {
      let h = qr.hits[hi]
      println("  >> " + h.id + "  " + h.description)
      let mut si = 0
      while si < h.hsps.length() {
        let hsp = h.hsps[si]
        let frag = hsp.fragments[0]
        println(
          "     domain " + (si + 1).to_string() +
          ": bitscore=" + hsp.bitscore.to_string() +
          "  evalue=" + hsp.evalue.to_string() +
          "  ali=" + frag.query_start.to_string() + "-" + frag.query_end.to_string(),
        )
        si = si + 1
      }
      hi = hi + 1
    }
    ti = ti + 1
  }
  println("")

  // ------------------------------------------------------------------
  // 3. Custom domtblout input
  // ------------------------------------------------------------------
  println("=== 3. Custom domtblout Input ===")
  let custom =
    "# header\n" +
    "target1 acc1 queryA qacc1 1e-50 200.0 1.0 1 2 1e-50 1e-53 100.0 1.0 1 80 1 80 1 80 1.0 1 0 0 1 1 1 1 desc1\n" +
    "target1 acc1 queryA qacc1 1e-50 200.0 1.0 2 2 1e-40 1e-43 100.0 1.0 81 150 81 150 81 150 1.0 1 0 0 1 1 1 1 desc1\n"
  let custom_results = @bio.parse_hmmer3_domtab(custom)
  println("Queries: " + custom_results.length().to_string())
  println(
    "First query id: " + custom_results[0].id +
    "  hits: " + custom_results[0].hits.length().to_string(),
  )
  println(
    "Multi-domain hit HSP count: " + custom_results[0].hits[0].hsps.length().to_string(),
  )
  println("")

  println("=== Demo Complete ===")
}
