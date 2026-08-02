///|
/// InterProScan Demo - Demonstrates parsing InterProScan TSV output.
///
/// Run: moon run examples/interproscan_demo
fn main {
  println("=== InterProScan Demo ===")
  println("InterProScan TSV output parser (14-column format)")
  println("Reference: InterProScan 5 genome-scale protein function classification")
  println("")

  // ------------------------------------------------------------------
  // 1. Parse sample InterProScan TSV data
  // ------------------------------------------------------------------
  println("=== 1. Parse sample data ===")
  let data = @bio.interproscan_sample_data()
  let records = @bio.parse_interproscan(data)
  println("Parsed \{records.length()} records from sample data")
  println("")

  // ------------------------------------------------------------------
  // 2. Display first few records with key fields
  // ------------------------------------------------------------------
  println("=== 2. First few records ===")
  let show_count = if records.length() < 3 { records.length() } else { 3 }
  for i in 0..<show_count {
    let r = records[i]
    println("  Record \{i + 1}:")
    println("    protein_id: \{r.protein_id()}")
    println("    analysis: \{r.analysis()}")
    println("    signature: \{r.signature_acc()} - \{r.signature_desc()}")
    println("    position: \{r.start_pos()}-\{r.stop_pos()}")
    println("    score: \{r.score()}")
    println("    IPR: \{r.ipr_acc()} (\{r.ipr_desc()})")
    println("    GO terms: \{r.go_terms().length()}")
  }
  println("")

  // ------------------------------------------------------------------
  // 3. Filter by analysis (Pfam)
  // ------------------------------------------------------------------
  println("=== 3. Filter by analysis: Pfam ===")
  let pfam = @bio.interproscan_filter_by_analysis(records, "Pfam")
  println("Found \{pfam.length()} Pfam record(s)")
  for r in pfam {
    println("  \{r.protein_id()} -> \{r.signature_acc()} (score=\{r.score()})")
  }
  println("")

  // ------------------------------------------------------------------
  // 4. Filter by protein ID
  // ------------------------------------------------------------------
  println("=== 4. Filter by protein: sp|P12345|PROT_HUMAN ===")
  let prot = @bio.interproscan_filter_by_protein(records, "sp|P12345|PROT_HUMAN")
  println("Found \{prot.length()} record(s)")
  for r in prot {
    println("  \{r.analysis()} | \{r.signature_acc()} | pos=\{r.start_pos()}-\{r.stop_pos()}")
  }
  println("")

  // ------------------------------------------------------------------
  // 5. Show all GO terms
  // ------------------------------------------------------------------
  println("=== 5. All GO terms ===")
  let go = @bio.interproscan_get_go_terms(records)
  println("Unique GO terms count: \{go.length()}")
  for g in go {
    println("  \{g}")
  }
  println("")

  // ------------------------------------------------------------------
  // 6. Show unique proteins and signatures
  // ------------------------------------------------------------------
  println("=== 6. Unique proteins and signatures ===")
  let proteins = @bio.interproscan_get_unique_proteins(records)
  let signatures = @bio.interproscan_get_unique_signatures(records)
  println("Unique proteins (\{proteins.length()}):")
  for p in proteins {
    println("  \{p}")
  }
  println("Unique signatures (\{signatures.length()}):")
  for s in signatures {
    println("  \{s}")
  }
  println("")

  // ------------------------------------------------------------------
  // 7. Group records by protein
  // ------------------------------------------------------------------
  println("=== 7. Group by protein ===")
  let groups = @bio.interproscan_group_by_protein(records)
  for g in groups {
    let (pid, recs) = g
    println("  \{pid}: \{recs.length()} record(s)")
    for r in recs {
      println("    - \{r.analysis()} / \{r.signature_acc()}")
    }
  }
  println("")

  // ------------------------------------------------------------------
  // 8. Generate summary
  // ------------------------------------------------------------------
  println("=== 8. Summary ===")
  println(@bio.interproscan_summary(records))

  println("=== Demo complete ===")
}
