///|
/// Prosite Demo - Demonstrates protein motif database search.

fn main {
  println("=== Prosite Demo ===")
  
  let patterns = @bio.prosite_create_example_patterns()
  println("Loaded \{patterns.length()} Prosite patterns")
  
  let pattern = @bio.prosite_get_pattern("PS00001")
  match pattern {
    Some(p) => println("Pattern PS00001: \{p.name} - \{p.pattern}")
    None => println("Pattern PS00001 not found")
  }
  
  println("")
  println("Searching for PKC phosphorylation site ([ST]-X-[RK]):")
  let matches = @bio.prosite_search("[ST]-X-[RK]", "AASTKKST")
  println("  Found \{matches.length()} matches")
  
  println("")
  println("Scanning sequence with all patterns:")
  let seq = "AASTKKSTGXXGKKST"
  let scan_results = @bio.prosite_scan(seq, patterns)
  println("  Sequence: \{seq}")
  println("  Found \{scan_results.length()} motif matches")
  
  println("")
  println("Prosite entry:")
  let entry = @bio.prosite_create_example_entry()
  println("  Accession: \{entry.accession}")
  println("  Name: \{entry.name}")
  println("  Pattern: \{entry.pattern}")
}