///|
/// SearchIO demo - demonstrates parsing HMMER3 tabular, BLAT PSL formats
/// and converting BLAST results to SearchIO model.
///
/// Run with: moon run examples/search_io_demo/main.mbt
fn main {
  println("=== SearchIO Demo ===\n")

  println("1. Parsing HMMER3 Tabular Format:")
  let hmmer3_data = "# hmmscan 3.3.2 (March 2023)\n" +
    "# query name           target name        accession   E-value  score  bias\n" +
    "PF00001.28            sp|Q9Y2W8|A1BG_HUMAN  Q9Y2W8     1e-50    200    5\n" +
    "PF00001.28            sp|P12345|ABC_HUMAN   P12345     1e-30    150    2\n" +
    "PF00001.28            sp|X99999|XYZ_HUMAN   X99999     0.5      50     1\n"

  let hmmer3_results = @bio.parse_hmmer3_tab(hmmer3_data)
  println("   Results found: " + hmmer3_results.length().to_string())
  if hmmer3_results.length() > 0 {
    let qr = hmmer3_results[0]
    println("   Query: " + qr.id)
    println("   Hits: " + qr.hits.length().to_string())

    let mut i = 0
    while i < qr.hits.length() {
      let hit = qr.hits[i]
      println("   Hit " + (i + 1).to_string() + ": " + hit.id)
      println("     E-value: " + hit.evalue.to_string())
      println("     Bitscore: " + hit.bitscore.to_string())
      i = i + 1
    }
  }

  println("\n2. Parsing BLAT PSL Format:")
  let psl_data = "100\t0\t0\t0\t0\t0\t0\t0\t+\tquery1\t1000\t0\t100\ttarget1\t2000\t0\t100\t1\t100,\t0,\t0,\n" +
    "80\t0\t0\t0\t0\t0\t0\t0\t-\tquery1\t1000\t100\t180\ttarget2\t1500\t50\t130\t1\t80,\t100,\t50,\n"

  let psl_results = @bio.parse_blat_psl(psl_data)
  println("   Results found: " + psl_results.length().to_string())
  if psl_results.length() > 0 {
    let qr = psl_results[0]
    println("   Query: " + qr.id)
    println("   Query length: " + qr.seq_len.to_string())
    println("   Hits: " + qr.hits.length().to_string())

    let mut i = 0
    while i < qr.hits.length() {
      let hit = qr.hits[i]
      println("   Hit " + (i + 1).to_string() + ": " + hit.id)
      i = i + 1
    }
  }

  println("\n3. Counting HSPs:")
  if hmmer3_results.length() > 0 {
    let count = @bio.count_hsps(hmmer3_results[0])
    println("   Total HSPs: " + count.to_string())
  }

  println("\n4. Getting Top Hits (by E-value):")
  if hmmer3_results.length() > 0 {
    let tops = @bio.top_hits(hmmer3_results[0], 2, 0.01)
    println("   Top 2 hits (E-value < 0.01):")

    let mut i = 0
    while i < tops.length() {
      let hit = tops[i]
      println(
        "   " +
        (i + 1).to_string() +
        ". " +
        hit.id +
        " (E-value: " +
        hit.evalue.to_string() +
        ")",
      )
      i = i + 1
    }
  }

  println("\n=== Demo Complete ===")
}
