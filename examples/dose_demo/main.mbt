///|
/// Bioconductor DOSE demo — Disease Ontology enrichment analysis.
///
/// Ported from the Bioconductor DOSE package (Yu et al., 2020),
/// this demo shows Fisher exact enrichment over a synthetic disease
/// term database, BH–FDR correction, and fold enrichment / odds ratio.

fn main {
  println("=== Bioconductor DOSE Demo ===")
  println("Disease Ontology enrichment analysis via Fisher's exact test")
  println("")

  // ---------------------------------------------------------------------
  // 1. Build a synthetic disease-term database.
  //    Each term contains a mix of query and background genes so we
  //    get clean enrichment signals.
  // ---------------------------------------------------------------------
  println("1. Build synthetic disease-term database")
  let db : Map[String, Array[String]] = Map([])

  // TERM_A: 28 of 180 query genes are in this term → strong enrichment
  let ta : Array[String] = Array::new()
  let mut i = 1
  while i <= 600 { ta.push("TA" + i.to_string()); i = i + 1 }
  let qa = ["Q001","Q002","Q003","Q004","Q005","Q006","Q007","Q008",
            "Q009","Q010","Q011","Q012","Q013","Q014","Q015","Q016",
            "Q017","Q018","Q019","Q020","Q021","Q022","Q023","Q024",
            "Q025","Q026","Q027","Q028"]
  i = 0
  while i < qa.length() { ta[i] = qa[i]; i = i + 1 }
  db["TERM_A (Cancer Signaling)"] = ta

  // TERM_B: 15 of 180 → weak enrichment
  let tb : Array[String] = Array::new()
  i = 1
  while i <= 900 { tb.push("TB" + i.to_string()); i = i + 1 }
  let qb = ["Q001","Q002","Q003","Q004","Q005","Q006","Q007","Q008",
            "Q009","Q010","Q011","Q012","Q013","Q014","Q015"]
  i = 0
  while i < qb.length() { tb[i] = qb[i]; i = i + 1 }
  db["TERM_B (Cell Cycle)"] = tb

  // TERM_C: 12 of 180 → moderate enrichment
  let tc : Array[String] = Array::new()
  i = 1
  while i <= 400 { tc.push("TC" + i.to_string()); i = i + 1 }
  let qc = ["Q001","Q002","Q003","Q004","Q005","Q006","Q007","Q008",
            "Q009","Q010","Q011","Q012"]
  i = 0
  while i < qc.length() { tc[i] = qc[i]; i = i + 1 }
  db["TERM_C (Apoptosis)"] = tc

  // TERM_D: 8 of 180 in a huge 1500-gene term → depletion
  let td : Array[String] = Array::new()
  i = 1
  while i <= 1500 { td.push("TD" + i.to_string()); i = i + 1 }
  let qd = ["Q001","Q002","Q003","Q004","Q005","Q006","Q007","Q008"]
  i = 0
  while i < qd.length() { td[i] = qd[i]; i = i + 1 }
  db["TERM_D (Ubiquitous Housekeeping)"] = td

  // TERM_E: 15 of 180 in a small 250-gene term → strong enrichment
  let te : Array[String] = Array::new()
  i = 1
  while i <= 250 { te.push("TE" + i.to_string()); i = i + 1 }
  let qe = ["Q010","Q011","Q012","Q013","Q014","Q015","Q016","Q017",
            "Q018","Q019","Q020","Q021","Q022","Q023","Q024"]
  i = 0
  while i < qe.length() { te[i] = qe[i]; i = i + 1 }
  db["TERM_E (DNA Repair)"] = te

  println("   Number of disease terms: " + db.keys().collect().length().to_string())

  // ---------------------------------------------------------------------
  // 2. Build query gene list and background.
  // ---------------------------------------------------------------------
  println("\n2. Query and background")
  let query : Array[String] = Array::new()
  i = 1
  while i <= 180 {
    // Pad to 3 digits so "Q1" → "Q001", matching the term-gene names.
    let s = i.to_string()
    let suffix = if i < 10 {
      "00" + s
    } else if i < 100 {
      "0" + s
    } else {
      s
    }
    query.push("Q" + suffix)
    i = i + 1
  }
  let bg : Array[String] = Array::new()
  let keys = db.keys().collect()
  i = 0
  while i < keys.length() {
    let genes = db[keys[i]]
    let mut j = 0
    while j < genes.length() { bg.push(genes[j]); j = j + 1 }
    i = i + 1
  }
  // Add query genes to background too.
  i = 0
  while i < query.length() { bg.push(query[i]); i = i + 1 }
  // Pad to 12000 with unique markers.
  i = 1
  while bg.length() < 12000 { bg.push("BG" + i.to_string()); i = i + 1 }
  println("   Query size:     " + query.length().to_string())
  println("   Background size: " + bg.length().to_string())

  // ---------------------------------------------------------------------
  // 3. Run Fisher-exact enrichment with BH–FDR correction.
  // ---------------------------------------------------------------------
  println("\n3. Enrichment analysis (Fisher exact + BH–FDR)")
  let result = @src.enrich_disease_fisher(query, db, bg)
  println("   Terms returned:  " + result.terms.length().to_string())
  println("   Summary:         " + @src.dose_summarize(result))

  // ---------------------------------------------------------------------
  // 4. Print top 3 results.
  // ---------------------------------------------------------------------
  println("\n4. Top 3 enriched terms")
  let top = @src.dose_top_terms(result, 3)
  i = 0
  while i < top.terms.length() {
    let t = top.terms[i]
    println("   #" + (i + 1).to_string() + "  " + t.term_id)
    println("      overlap=" + t.overlap_count.to_string() +
           "/" + t.gene_count.to_string() +
           "  FE=" + t.fold_enrichment.to_string() +
           "  OR=" + t.odds_ratio.to_string())
    println("      raw p=" + t.p_value.to_string() +
           "  adj p=" + t.adjusted_p_value.to_string())
    i = i + 1
  }

  // ---------------------------------------------------------------------
  // 5. Show Fisher-exact helpers directly.
  // ---------------------------------------------------------------------
  println("\n5. Direct Fisher-exact / hypergeometric calls")
  // 2×2 table: a=28, b=152, c=572, d=11248  (TERM_A enrichment)
  let r = @src.fisher_exact_greater(28, 152, 572, 11248)
  println("   Fisher greater (TERM_A): p=" + r.p_value.to_string() +
          "  OR=" + r.odds_ratio.to_string())

  // Hypergeometric: k=12, query=180, term=400, bg=12000  (TERM_C)
  let p_hyp = @src.dose_hypergeometric_test(12, 180, 400, 12000)
  println("   Hypergeometric (TERM_C): p=" + p_hyp.to_string())

  println("\n=== DOSE demo complete ===")
}
