///|
/// msigdbr demo - MSigDB gene set collections, GMT I/O and ID mapping.
/// Inspired by Bioconductor {msigdbr}.

fn main {
  println("=== msigdbr Demo: MSigDB Gene Set Collections ===")
  println("")

  // Step 1: built-in catalog overview
  println("--- Step 1: Built-in Catalog ---")
  let catalog = @src.msig_catalog()
  for col in catalog {
    println(
      "  [" + col.category.code() + "] " + col.subcategory +
      " - " + col.size().to_string() + " sets",
    )
  }
  println("  total sets: " + @src.msig_catalog_size(catalog).to_string())
  println("")

  // Step 2: inspect a hallmark set
  println("--- Step 2: Inspect a Hallmark Set ---")
  match @src.msig_lookup(catalog, "HALLMARK_P53_PATHWAY") {
    Some(gs) => {
      println("  " + gs.name + " (" + gs.description + ")")
      let buf = StringBuilder()
      let mut j = 0
      while j < gs.gene_ids.length() {
        if j > 0 {
          buf.write_string(", ")
        }
        buf.write_string(gs.gene_ids[j])
        j = j + 1
      }
      println("  genes (" + gs.size().to_string() + "): " + buf.to_string())
    }
    None => println("  not found")
  }
  println("")

  // Step 3: GMT round-trip
  println("--- Step 3: GMT Round-trip ---")
  let gmt = @src.msig_to_gmt(@src.msig_hallmark())
  let first_line = gmt.split("\n").next()
  match first_line {
    Some(line) => {
      let owned = line.to_owned()
      let end = if owned.length() > 80 { 80 } else { owned.length() }
      println("  GMT line 1 (truncated): " + owned.substring(start=0, end=end) + "...")
    }
    None => println("  empty GMT")
  }
  let parsed = @src.msig_from_gmt(gmt, @src.MsigHallmark)
  println("  parsed back: " + parsed.size().to_string() + " sets")
  println("")

  // Step 4: symbol -> Entrez mapping
  println("--- Step 4: Gene ID Mapping ---")
  let table = @src.msig_id_table([("TP53", "7157"), ("BRCA1", "672"), ("BAX", "581")])
  let mapped = @src.msig_symbol_to_entrez(@src.msig_ontology_go_bp(), table)
  match @src.msig_lookup([mapped], "GO_APOPTOTIC_PROCESS") {
    Some(gs) => {
      let buf = StringBuilder()
      let mut j = 0
      while j < gs.gene_ids.length() {
        if j > 0 {
          buf.write_string(", ")
        }
        buf.write_string(gs.gene_ids[j])
        j = j + 1
      }
      println("  GO_APOPTOTIC_PROCESS (entrez): " + buf.to_string())
    }
    None => println("  not found")
  }
  println("")

  // Step 5: gene frequency (hub genes)
  println("--- Step 5: Hub Genes Across Sets ---")
  let freq = @src.msig_gene_frequency(@src.msig_hallmark())
  println("  distinct genes in H collection: " + freq.size().to_string())
  match freq.get("TP53") {
    Some(n) => println("  TP53 appears in " + n.to_string() + " sets")
    None => println("  TP53 not present")
  }
  println("")

  // Step 6: summary
  println("--- Step 6: Collection Summary ---")
  println(@src.msig_summary(@src.msig_hallmark()))

  println("")
  println("=== msigdbr Demo Complete ===")
}
