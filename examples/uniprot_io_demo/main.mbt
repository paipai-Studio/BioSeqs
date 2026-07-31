///|
/// UniProt XML Format IO Demo.
///
/// Demonstrates parsing UniProt XML entries and extracting annotations.
///
/// Run: moon run examples/uniprot_io_demo/main.mbt
fn main {
  println("=== UniProt XML IO Demo ===\n")

  // Example UniProt XML content
  let uniprot_xml = "<uniprot>\n" +
    "    <entry>\n" +
    "      <accession>P12345</accession>\n" +
    "      <protein><recommendedName><fullName>Hemoglobin subunit alpha</fullName></recommendedName></protein>\n" +
    "      <gene><name>HBA1</name></gene>\n" +
    "      <organism><ncbi_taxonomy_id>9606</ncbi_taxonomy_id><scientificName>Homo sapiens</scientificName></organism>\n" +
    "      <sequence length=\"141\" mass=\"15126.0\">VLSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSFPTTKTYFPHFDLSHGSAQVKGHGKKVADALTNAVAHVDDMPNALSALSDLHAHKLRVDPVNFKLLSHCLLVTLAAHLPAEFTPAVHASLDKFLASVSTVLTSKYR</sequence>\n" +
    "      <dbReference type=\"PDB\" id=\"1A3N\">\n" +
    "        <property type=\"id\" value=\"1A3N\"/>\n" +
    "      </dbReference>\n" +
    "      <dbReference type=\"Pfam\" id=\"PF00042\">\n" +
    "        <property type=\"id\" value=\"PF00042\"/>\n" +
    "      </dbReference>\n" +
    "      <feature>\n" +
    "        <location><begin position=\"1\"/><end position=\"141\"/></location>\n" +
    "        <evidence>\n" +
    "          <source>test</source>\n" +
    "          <type>test</type>\n" +
    "        </evidence>\n" +
    "      </feature>\n" +
    "      <comment type=\"subcellular location\">\n" +
    "        <subcellularLocation>\n" +
    "          <location>Cytoplasm</location>\n" +
    "        </subcellularLocation>\n" +
    "      </comment>\n" +
    "    </entry>\n" +
    "  </uniprot>"

  // Parse UniProt XML
  println("1. Parse UniProt XML:")
  let entries = @bio.parse_uniprot_xml(uniprot_xml)
  println("   Found \{entries.length()} entries")

  if entries.length() > 0 {
    let entry = entries[0]
    println("\n2. Entry details:")
    println("   Accession: \{entry.accession}")
    println("   Protein: \{entry.protein_name}")
    println("   Organism: \{entry.organism}")
    println("   Sequence length: \{entry.seq_length}")
    let seq_str = entry.sequence
    let prefix = if seq_str.length() > 20 { seq_str.substring(start=0, end=20) + "..." } else { seq_str }
    println("   Sequence (first 20 chars): \{prefix}")

    println("\n3. Gene names:")
    for g in entry.gene_names {
      println("   - \{g}")
    }

    println("\n4. Database references:")
    for db in entry.db_references {
      println("   - \{db.database}: \{db.id}")
    }

    println("\n5. Features:")
    for f in entry.features {
      println("   - \{f.feature_type}: \{f.description}")
    }

    println("\n6. Convert to SeqRecord-compatible format:")
    let (header, id, desc, kwargs) = @bio.uniprot_to_seqrecord(entry)
    println("   Header: \{header}")
    println("   ID: \{id}")
    println("   Description: \{desc}")
  }

  println("\n=== Demo Complete ===")
}
