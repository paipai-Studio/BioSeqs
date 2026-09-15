///|
/// Demo for mirna_targets module: miRNA target prediction.

fn main {
  println("=== mirna_targets Demo: miRNA Target Prediction ===")
  println("")

  // Example 1: Extract miRNA seed (DNA alphabet)
  println("--- Step 1: Extract miRNA Seed ---")
  // Note: Use DNA alphabet (T instead of U) for sequences
  let mirna = "TGGAATGTAAAGAAGTATGTAT"  // hsa-miR-1 (DNA)
  let seed = @src.mir_extract_seed(mirna)
  let seed6 = @src.mir_extract_seed6(mirna)
  println("miRNA: " + mirna)
  println("Seed (pos 2-8): " + seed)
  println("Seed6 (pos 2-7): " + seed6)
  println("")

  // Example 2: Find seed matches in a 3'UTR
  println("--- Step 2: Find Seed Matches ---")
  // miRNA seed = "GGAATGT", rc = "ACATTCC"
  // UTR with "ACATTCC" preceded by A -> 8mer
  let utr = "AATACATTCCGGGUUU"
  let sites = @src.mir_find_seed_matches(mirna, utr, mir_id="hsa-miR-1", target_id="TEST_GENE")
  println("3'UTR: " + utr)
  println("Sites found: " + sites.length().to_string())
  for s in sites {
    println("  pos=" + s.utr_pos.to_string() + " len=" + s.site_len.to_string() +
            " type=" + @src.mir_seed_type_label(s.seed_type) +
            " score=" + s.score.to_string())
  }
  println("")

  // Example 3: Predict targets for multiple miRNAs
  println("--- Step 3: Predict Targets ---")
  let mirnas = [
    ("hsa-let-7a", "TGAGGTAGTAGGTTGTATAGTT"),
    ("hsa-miR-1", "TGGAATGTAAAGAAGTATGTAT"),
  ]
  let utrs = [
    ("GENE1", "AACTACCTCCTCAACATACAACCTACTACCTCATT"),
    ("GENE2", "AATACATTCCGGGUUU"),
  ]
  let results = @src.mir_predict_targets(mirnas, utrs, min_score=0.1)
  println("Results: " + results.length().to_string())
  for r in results {
    println("  " + r.mir_id + " -> " + r.target_id + " (score=" + r.total_score.to_string() + ")")
  }
  println("")

  // Example 4: Export sites
  println("--- Step 4: Export Sites ---")
  let tsv = @src.mir_export_sites(results)
  println(tsv)
  println("")

  // Example 5: Count by seed type
  println("--- Step 5: Count by Seed Type ---")
  let counts = @src.mir_count_by_seed_type(results)
  println("8mer: " + counts.get("8mer").unwrap_or(0).to_string())
  println("7mer-m8: " + counts.get("7mer-m8").unwrap_or(0).to_string())
  println("7mer-A1: " + counts.get("7mer-A1").unwrap_or(0).to_string())
  println("6mer: " + counts.get("6mer").unwrap_or(0).to_string())
}
