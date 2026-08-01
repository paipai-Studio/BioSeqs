///|
/// GOA (Gene Association Format) demo - Gene Ontology Annotation parser.
///
/// Run: moon run examples/goa_demo
fn main {
  println("=== GOA (Gene Association Format) Demo ===")
  println("Gene Ontology Annotation file parser (GAF 2.2)")
  println("Reference: Gene Ontology Consortium GAF specification")
  println("")

  // ------------------------------------------------------------------
  // 1. Load sample GAF database
  // ------------------------------------------------------------------
  println("=== 1. Load sample GAF database ===")
  let db = @bio.goa_sample_database()
  println("Loaded sample GOA database")
  println("")

  // ------------------------------------------------------------------
  // 2. Print database summary (n_records, version, created_by)
  // ------------------------------------------------------------------
  println("=== 2. Database summary ===")
  println("n_records: \{db.n_records()}")
  println("version: \{db.version()}")
  println("created_by: \{db.created_by()}")
  println("")

  // ------------------------------------------------------------------
  // 3. Show first few records with key fields
  // ------------------------------------------------------------------
  println("=== 3. First few records ===")
  let records = db.records()
  let show_count = if records.length() < 3 { records.length() } else { 3 }
  for i in 0..<show_count {
    let r = records[i]
    println("  Record \{i + 1}:")
    println("    db_object_symbol: \{r.db_object_symbol()}")
    println("    go_id: \{r.go_id()}")
    println("    evidence_code: \{r.evidence_code()}")
    println("    aspect: \{r.aspect()}")
  }
  println("")

  // ------------------------------------------------------------------
  // 4. Filter by aspect (MolecularFunction) and show results
  // ------------------------------------------------------------------
  println("=== 4. Filter by aspect: Molecular Function ===")
  let mf = @bio.goa_filter_by_aspect(db, @bio.GafAspect::from_string("F"))
  println("Found \{mf.length()} Molecular Function record(s)")
  for r in mf {
    println("  \{r.db_object_symbol()} -> \{r.go_id()} (\{r.evidence_code()})")
  }
  println("")

  // ------------------------------------------------------------------
  // 5. Filter by evidence code (IEA) and show results
  // ------------------------------------------------------------------
  println("=== 5. Filter by evidence code: IEA ===")
  let iea = @bio.goa_filter_by_evidence(db, "IEA")
  println("Found \{iea.length()} IEA record(s)")
  for r in iea {
    println("  \{r.db_object_symbol()} -> \{r.go_id()} [aspect=\{r.aspect()}]")
  }
  println("")

  // ------------------------------------------------------------------
  // 6. Filter by taxon (9606 - human) and show count
  // ------------------------------------------------------------------
  println("=== 6. Filter by taxon: 9606 (human) ===")
  let human = @bio.goa_filter_by_taxon(db, "9606")
  println("Found \{human.length()} human record(s)")
  println("")

  // ------------------------------------------------------------------
  // 7. Show unique GO IDs count
  // ------------------------------------------------------------------
  println("=== 7. Unique GO IDs ===")
  let go_ids = @bio.goa_unique_go_ids(db)
  println("Unique GO IDs count: \{go_ids.length()}")
  for gid in go_ids {
    println("  \{gid}")
  }
  println("")

  // ------------------------------------------------------------------
  // 8. Show count by aspect using the summary function
  // ------------------------------------------------------------------
  println("=== 8. Count by aspect ===")
  let counts = @bio.goa_count_by_aspect(db)
  let entries = counts.iter().collect()
  for i in 0..<entries.length() {
    let (aspect, count) = entries[i]
    let desc = @bio.GafAspect::description(@bio.GafAspect::from_string(aspect))
    println("  \{aspect} (\{desc}): \{count}")
  }
  println("")
  println("--- Full summary ---")
  println(@bio.goa_to_summary(db))

  // ------------------------------------------------------------------
  // 9. Demonstrate writing: convert a record to GAF line format
  // ------------------------------------------------------------------
  println("=== 9. Write record to GAF line ===")
  let first = records[0]
  let line = @bio.goa_record_to_gaf_line(first)
  println("GAF line: \{line}")
  println("")

  // ------------------------------------------------------------------
  // 10. Full round-trip: parse a GAF string, then write back
  // ------------------------------------------------------------------
  println("=== 10. Full round-trip: parse GAF string then write back ===")
  let gaf_text = @bio.goa_to_gaf(db)
  let reparsed = @bio.goa_parse(gaf_text)
  println("Original n_records: \{db.n_records()}")
  println("Reparsed n_records: \{reparsed.n_records()}")
  let rewritten = @bio.goa_to_gaf(reparsed)
  println("Round-trip preserved: \{gaf_text == rewritten}")
  println("")

  println("=== Demo complete ===")
}
