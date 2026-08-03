///|
/// INSDC (GenBank/EMBL) IO Demo.
///
/// Demonstrates Bio.SeqIO.InsdcIO module: parsing and writing GenBank
/// and EMBL records, feature tables, and qualifier extraction.
/// Inspired by BioPython's Bio.SeqIO.InsdcIO module.

fn main {
  println("=== INSDC IO Demo ===")
  println("Bio.SeqIO.InsdcIO: GenBank/EMBL parsing and writing\n")

  // === 1. Construct a GenBank format text ===
  println("=== 1. Construct GenBank Text ===")
  let genbank_text = "LOCUS       TESTSEQ     60 bp    DNA     linear   UNK 01-JAN-2024\n" +
    "DEFINITION  Test sequence for InsdcIO demo.\n" +
    "ACCESSION   TEST0001\n" +
    "VERSION     TEST0001.1\n" +
    "ORGANISM  Escherichia coli\n" +
    "FEATURES         Location/Qualifiers\n" +
    "     gene            1..60\n" +
    "                     /gene=\"rpoB\"\n" +
    "                     /locus_tag=\"ECK_1234\"\n" +
    "     CDS             1..60\n" +
    "                     /gene=\"rpoB\"\n" +
    "                     /product=\"DNA-directed RNA polymerase subunit beta\"\n" +
    "                     /translation=\"ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT\"\n" +
    "ORIGIN\n" +
    "        1 acgtacgtac gtacgtacgt acgtacgtac gtacgtacgt acgtacgtac gtacgtacgt\n" +
    "//\n"
  println("  Constructed GenBank record with LOCUS/DEFINITION/ACCESSION/")
  println("  VERSION/FEATURES(gene+CDS)/ORIGIN sections.")

  // === 2. Parse the GenBank record ===
  println("\n=== 2. Parse GenBank Record ===")
  let record = @bio.parse_genbank_record(genbank_text)
  println("  parse_genbank_record: OK")

  // === 3. Print parsed metadata ===
  println("\n=== 3. Parsed Metadata ===")
  println("  locus_name:       \{record.locus_name}")
  println("  accession:        \{record.accession}")
  println("  version:          \{record.version}")
  println("  organism:         \{record.organism}")
  println("  sequence length:  \{record.sequence.length()}")
  println("  molecule_type:    \{record.molecule_type}")
  println("  topology:         \{record.topology}")
  println("  division:         \{record.division}")
  println("  definition:       \{record.definition}")
  println("  features count:   \{record.features.length()}")

  // === 4. Iterate features, print type/location/qualifiers ===
  println("\n=== 4. Features ===")
  let mut fi = 0
  while fi < record.features.length() {
    let f = record.features[fi]
    let loc = @bio.parse_insdc_location(f.location_str)
    println("  Feature [\{fi}] type=\{f.feature_type}, location=\{loc}")
    let qkeys = f.qualifiers.keys().collect()
    qkeys.sort()
    let mut qi = 0
    while qi < qkeys.length() {
      let k = qkeys[qi]
      let vals = f.qualifiers.get_or_default(k, [])
      println("    /\{k} = \{vals}")
      qi = qi + 1
    }
    fi = fi + 1
  }

  // === 5. Extract CDS /gene, /product, /translation ===
  println("\n=== 5. Extract CDS Qualifiers ===")
  fi = 0
  while fi < record.features.length() {
    let f = record.features[fi]
    if f.feature_type == "CDS" {
      let gene = @bio.extract_gene_name(f).unwrap_or("(none)")
      let product = @bio.extract_feature_qualifier(f, "product").unwrap_or("(none)")
      let translation = @bio.extract_translation(f).unwrap_or("(none)")
      println("  /gene:        \{gene}")
      println("  /product:     \{product}")
      println("  /translation: \{translation}")
    }
    fi = fi + 1
  }

  // === 6. Write GenBank record back ===
  println("\n=== 6. Write GenBank Record ===")
  let genbank_rewritten = @bio.write_genbank_record(record)
  println(genbank_rewritten)

  // === 7. Construct and parse an EMBL record ===
  println("=== 7. Parse EMBL Record ===")
  let embl_text = "ID   TESTSEQ; SV 1; linear; DNA; STD; UNK; 60 BP.\n" +
    "XX\n" +
    "AC   TEST0001;\n" +
    "XX\n" +
    "DE   Test sequence for InsdcIO demo.\n" +
    "XX\n" +
    "OS   Escherichia coli\n" +
    "XX\n" +
    "FH   Key             Location/Qualifiers\n" +
    "FH\n" +
    "FT   gene            1..60\n" +
    "FT                   /gene=\"rpoB\"\n" +
    "FT                   /locus_tag=\"ECK_1234\"\n" +
    "FT   CDS             1..60\n" +
    "FT                   /gene=\"rpoB\"\n" +
    "FT                   /product=\"DNA-directed RNA polymerase subunit beta\"\n" +
    "FT                   /translation=\"ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT\"\n" +
    "XX\n" +
    "SQ   Sequence 60 BP; 15 A; 15 C; 15 G; 15 T; 0 other;\n" +
    "     acgtacgtac gtacgtacgt acgtacgtac gtacgtacgt acgtacgtac gtacgtacgt\n" +
    "//\n"
  let embl_record = @bio.parse_embl_record(embl_text)
  println("  parse_embl_record: OK")
  println("  locus_name:       \{embl_record.locus_name}")
  println("  accession:        \{embl_record.accession}")
  println("  version:          \{embl_record.version}")
  println("  organism:         \{embl_record.organism}")
  println("  sequence length:  \{embl_record.sequence.length()}")
  println("  molecule_type:    \{embl_record.molecule_type}")
  println("  topology:         \{embl_record.topology}")
  println("  features count:   \{embl_record.features.length()}")

  // === 8. Write EMBL record back ===
  println("\n=== 8. Write EMBL Record ===")
  let embl_rewritten = @bio.write_embl_record(embl_record)
  println(embl_rewritten)

  println("=== INSDC IO Demo Complete ===")
}
