///|
/// SeqXML demo - XML-based sequence exchange format.
///
/// Run: moon run examples/seqxml_demo
fn main {
  println("=== SeqXML Demo ===")
  println("XML format for sequence data exchange")
  println("Reference: BioPython Bio.SeqIO.SeqXmlIO")
  println("")

  // ------------------------------------------------------------------
  // 1. Create a sample SeqXML document
  // ------------------------------------------------------------------
  println("=== 1. Create sample SeqXML document ===")
  let doc = @bio.seqxml_sample_document()
  println("Created sample SeqXML document")
  println("")

  // ------------------------------------------------------------------
  // 2. Print document info (n_entries, species_name, source_db)
  // ------------------------------------------------------------------
  println("=== 2. Document info ===")
  println("n_entries: \{doc.n_entries()}")
  println("species_name: \{doc.species_name()}")
  println("source_db: \{doc.source_db()}")
  println("")

  // ------------------------------------------------------------------
  // 3. Show all entries (id, type, description, sequence)
  // ------------------------------------------------------------------
  println("=== 3. All entries ===")
  let entries = doc.entries()
  for e in entries {
    let type_desc = @bio.SeqXmlType::description(e.seq_type())
    println("  id: \{e.id()}")
    println("    type: \{type_desc}")
    println("    description: \{e.description()}")
    println("    sequence: \{e.sequence()}")
  }
  println("")

  // ------------------------------------------------------------------
  // 4. Parse a SeqXML string and show results
  // ------------------------------------------------------------------
  println("=== 4. Parse SeqXML string ===")
  let xml_str = "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n<seqXML>\n  <species name=\"Mus musculus\" ncbiTaxID=\"10090\"/>\n  <sourceDB name=\"Ensembl\"/>\n  <entry id=\"parsed1\" desc=\"Parsed DNA\">\n    <dna>GGGATCCC</dna>\n  </entry>\n  <entry id=\"parsed2\" desc=\"Parsed protein\">\n    <protein>MKLV</protein>\n  </entry>\n</seqXML>\n"
  let parsed = @bio.seqxml_parse(xml_str)
  println("Parsed \{parsed.n_entries()} entries")
  println("species_name: \{parsed.species_name()}")
  println("source_db: \{parsed.source_db()}")
  for e in parsed.entries() {
    let type_desc = @bio.SeqXmlType::description(e.seq_type())
    println("  \{e.id()} [\{type_desc}]: \{e.sequence()}")
  }
  println("")

  // ------------------------------------------------------------------
  // 5. Write document back to XML and print
  // ------------------------------------------------------------------
  println("=== 5. Write document back to XML ===")
  let xml_out = @bio.seqxml_to_xml(doc)
  println(xml_out)

  // ------------------------------------------------------------------
  // 6. Convert to SeqRecord array and show details
  // ------------------------------------------------------------------
  println("=== 6. Convert to SeqRecord array ===")
  let seq_records = @bio.seqxml_to_seq_records(doc)
  println("Converted to \{seq_records.length()} SeqRecord(s)")
  for rec in seq_records {
    println("  \{rec.repr()}")
  }
  println("")

  // ------------------------------------------------------------------
  // 7. Create from SeqRecord array and verify round-trip
  // ------------------------------------------------------------------
  println("=== 7. Create from SeqRecord array (round-trip) ===")
  let rebuilt = @bio.seqxml_from_seq_records(seq_records)
  println("Rebuilt document n_entries: \{rebuilt.n_entries()}")
  for e in rebuilt.entries() {
    let type_desc = @bio.SeqXmlType::description(e.seq_type())
    println("  \{e.id()} [\{type_desc}]: \{e.sequence()}")
  }
  println("")

  println("=== Demo complete ===")
}
