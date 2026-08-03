///|
/// SMART domain database format parser demo.
///
/// Demonstrates the SMART protein domain analysis module for:
///   - Parsing SMART domain annotations
///   - Querying domains by name and evalue
///   - Domain classification (SMART, Pfam, signal peptide, transmembrane)
///   - Summary generation
///
/// Inspired by Biopython's Bio.Smart module.
/// Run: moon run examples/smart_demo
fn main {
  println("=== SMART Domain Parser Demo ===")
  println("Protein domain annotation analysis")
  println("Reference: Biopython Bio.Smart")
  println("")

  // ------------------------------------------------------------------
  // 1. Parse SMART sample data
  // ------------------------------------------------------------------
  println("=== 1. Parsing SMART Sample Data ===")
  let content = @bio.smart_sample()
  println("Input content (first 100 chars):")
  println(content[0:100].to_owned())
  println("")

  let result = @bio.smart_parse(content)
  println("Parsed results:")
  println("  Sequence ID: " + result.sequence_id)
  println("  Sequence Length: " + result.sequence_length.to_string())
  println("  Total Domains: " + @bio.smart_total_domains(result).to_string())
  println("  SMART Domains: " + result.domains.length().to_string())
  println("  Pfam Domains: " + result.pfam_domains.length().to_string())
  println("  Signal Peptides: " + result.signal_peptides.length().to_string())
  println("  Transmembrane: " + result.transmembrane.length().to_string())
  println("")

  // ------------------------------------------------------------------
  // 2. List all domains
  // ------------------------------------------------------------------
  println("=== 2. Domain Details ===")
  let mut i = 0
  while i < result.domains.length() {
    let d = result.domains[i]
    println("  Domain " + (i + 1).to_string() + ": " + d.domain_name)
    println("    Type: " + d.domain_type)
    println("    Location: " + d.start.to_string() + " - " + d.end.to_string())
    println("    E-value: " + d.evalue.to_string())
    println("    Description: " + d.description)
    println("")
    i = i + 1
  }

  // ------------------------------------------------------------------
  // 3. Query domains by name
  // ------------------------------------------------------------------
  println("=== 3. Query: Find 'ABC' domains ===")
  let abc_domains = @bio.smart_find_domains(result, "ABC")
  println("Found " + abc_domains.length().to_string() + " domains matching 'ABC':")
  let mut q = 0
  while q < abc_domains.length() {
    println("  - " + abc_domains[q].domain_name + " (evalue: " + abc_domains[q].evalue.to_string() + ")")
    q = q + 1
  }
  println("")

  // ------------------------------------------------------------------
  // 4. Filter by evalue
  // ------------------------------------------------------------------
  println("=== 4. Filter by E-value (<= 1e-10) ===")
  let significant = @bio.smart_filter_evalue(result, 1.0e-10)
  println("Found " + significant.length().to_string() + " significant domains:")
  let mut s = 0
  while s < significant.length() {
    println("  - " + significant[s].domain_name + " (evalue: " + significant[s].evalue.to_string() + ")")
    s = s + 1
  }
  println("")

  // ------------------------------------------------------------------
  // 5. Best domain hit
  // ------------------------------------------------------------------
  println("=== 5. Best Domain Hit ===")
  match @bio.smart_best_domain(result) {
    Some(best) => {
      println("Best domain: " + best.domain_name)
      println("  E-value: " + best.evalue.to_string())
      println("  Location: " + best.start.to_string() + " - " + best.end.to_string())
    }
    None => println("  No domains found")
  }
  println("")

  // ------------------------------------------------------------------
  // 6. Summary report
  // ------------------------------------------------------------------
  println("=== 6. Summary Report ===")
  let summary = @bio.smart_summary(result)
  println(summary)

  println("=== SMART Demo Complete ===")
}
