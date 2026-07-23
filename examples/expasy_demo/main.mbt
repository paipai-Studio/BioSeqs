///|
/// ExPASy Demo - Demonstrates protein analysis tools interface.

fn main {
  println("=== ExPASy Demo ===")
  
  let entry = @bio.create_example_expasy_entry()
  println("Accession: \{entry.accession}")
  println("Name: \{entry.name}")
  println("Keywords: \{entry.keywords.length()}")
  
  println("")
  println("Protein analysis:")
  let results = @bio.expasy_analyze_protein("AVG")
  println("  Molecular weight: \{results.get("molecular_weight").unwrap_or(0.0)}")
  println("  GRAVY: \{results.get("gravy").unwrap_or(0.0)}")
  println("  Instability index: \{results.get("instability_index").unwrap_or(0.0)}")
  println("  Isoelectric point: \{results.get("isoelectric_point").unwrap_or(0.0)}")
  
  println("")
  println("Enzyme entry:")
  let enzyme = @bio.create_example_enzyme_entry()
  println("  EC number: \{enzyme.ec_number}")
  println("  Name: \{enzyme.name}")
  
  println("")
  println("EC parsing:")
  let (class, subclass, subsubclass, serial) = @bio.enzyme_parse_ec("1.1.1.1")
  println("  Class: \{class}, Subclass: \{subclass}, Subsubclass: \{subsubclass}, Serial: \{serial}")
}