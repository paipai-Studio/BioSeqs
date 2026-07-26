///|
/// PDB Header Demo
/// Demonstrates Bio.PDB.ParsePDBHeader functionality:
///   - Parse PDB header text
///   - Display PDB ID, title, compound info, source info
///   - Display resolution, R-factor, keywords, chain IDs
///   - Display authors, DBREF entries
///   - Show individual record parsing
///   - Demonstrate edge cases
///
/// Run: moon run examples/pdb_header_demo/

fn main {
  println("=== PDB Header Demo ===")

  // =========================================================================
  // 1. Create sample PDB header text
  // =========================================================================
  println("")
  println("--- 1. Sample PDB Header ---")

  let header_line =
    "HEADER      TRANS 05-JAN-21                      1ABC                           "
  let pdb_text =
    header_line + "\n" +
    "TITLE     CRYSTAL STRUCTURE OF A TRANSFERASE ENZYME COMPLEXED WITH\n" +
    "TITLE     ITS SUBSTRATE AND COFACTOR\n" +
    "COMPND    1 ABC; A; ENGINEERED; MUTATION;\n" +
    "COMPND    2 DEF; B;\n" +
    "SOURCE    1 ORGANISM Escherichia coli; COMMON_NAME Bacterium;\n" +
    "SOURCE    2 ORGANISM Homo sapiens; TAX_ID 9606;\n" +
    "KEYWDS    TRANSFERASE, METABOLISM, ENZYME;\n" +
    "KEYWDS    CRYSTAL STRUCTURE;\n" +
    "REMARK   1 RESOLUTION. 2.10\n" +
    "REMARK   3 R FACTOR. 0.195\n" +
    "AUTHOR    JOHN DOE, JANE SMITH, BOB JONES;\n" +
    "JRNL      Article title of the publication;\n" +
    "DBREF     1ABC A  UNP   P00000   1-300;\n" +
    "DBREF     1ABC A  PDB   1XYZ    1-300;\n"

  println("Input PDB text:")
  println(pdb_text)

  // =========================================================================
  // 2. Parse using parse_pdb_header()
  // =========================================================================
  println("")
  println("--- 2. Parse Full Header ---")

  let header = @bio.parse_pdb_header(pdb_text)
  println("PDB ID: \{header.pdb_id}")
  println("Date Deposited: \{header.date_deposited}")
  println("Title: \{header.title}")

  println("")
  println("--- 3. Compound Info ---")
  let compound = header.compound
  println("Molecule ID: \{compound.molecule_id}")
  println("Chain IDs (\{compound.chain_ids.length()}):")
  for chain in compound.chain_ids {
    println("  - \{chain}")
  }
  println("Compound Name: \{compound.compound_name}")
  println("Formula: \{compound.formula}")
  println("KA: \{compound.ka}")
  println("KB: \{compound.kb}")
  println("EC Number: \{compound.ec_number}")
  println("Engineered: \{compound.engineered}")
  println("Mutation: \{compound.mutation}")

  println("")
  println("--- 4. Source Info ---")
  let source = header.source
  println("Source: \{source.source}")
  println("Organism (scientific): \{source.organism_scientific}")
  println("Organism (common): \{source.organism_common}")
  println("Strain: \{source.strain}")
  println("Cell Line: \{source.cell_line}")
  println("Tissue: \{source.tissue}")
  println("Tax ID: \{source.tax_id}")

  println("")
  println("--- 5. Quality Metrics ---")
  println("Resolution: \{header.resolution} A")
  println("R-factor: \{header.rfactor}")

  println("")
  println("--- 6. Keywords ---")
  let keywords = header.keywords
  println("Keywords (\{keywords.length()}):")
  for kw in keywords {
    println("  - \{kw}")
  }

  println("")
  println("--- 7. Chain IDs (via extract_chain_ids) ---")
  let chain_ids = @bio.extract_chain_ids(header)
  println("Chain IDs (\{chain_ids.length()}):")
  for chain in chain_ids {
    println("  - \{chain}")
  }

  println("")
  println("--- 8. Authors ---")
  let authors = header.authors
  println("Authors (\{authors.length()}):")
  for author in authors {
    println("  - \{author}")
  }

  println("")
  println("--- 9. Journal ---")
  println("Journal: \{header.journal}")

  println("")
  println("--- 10. DBREF Entries ---")
  let dbrefs = header.dbrefs
  println("DBREF entries (\{dbrefs.length()}):")
  for dbref in dbrefs {
    println("  DB: \{dbref.db_name} | Accession: \{dbref.db_accession} | ID: \{dbref.db_id}")
  }

  println("")
  println("--- 11. Remark Lines ---")
  let remarks = header.remark_lines
  println("Remark lines (\{remarks.length()}):")
  for remark in remarks {
    println("  \{remark}")
  }

  // =========================================================================
  // 12. Individual record parsing
  // =========================================================================
  println("")
  println("--- 12. Individual Record Parsing ---")

  let header_line =
    "HEADER      TRANS 05-JAN-21                      1ABC                           "
  let (pdb_id, date, classification) = @bio.parse_header_record(header_line)
  println("parse_header_record:")
  println("  PDB ID: \{pdb_id}")
  println("  Date: \{date}")
  println("  Classification: \{classification}")

  let remark2 = "REMARK   2 RESOLUTION. 1.85"
  let resolution = @bio.parse_resolution_record(remark2)
  println("")
  println("parse_resolution_record: \{resolution} A")

  let remark3 = "REMARK   3 R FACTOR. 0.172"
  let rfactor = @bio.parse_rfactor_record(remark3)
  println("parse_rfactor_record: \{rfactor}")

  // =========================================================================
  // 13. Edge cases
  // =========================================================================
  println("")
  println("--- 13. Edge Cases ---")

  println("Empty text:")
  let empty_header = @bio.parse_pdb_header("")
  println("  PDB ID: '\{empty_header.pdb_id}'")
  println("  Resolution: \{empty_header.resolution}")
  println("  Keywords count: \{empty_header.keywords.length()}")

  println("")
  println("Minimal PDB (only HEADER record):")
  let minimal_text =
    "HEADER      LIGAN 01-JAN-20                      2XYZ                           \n"
  let minimal = @bio.parse_pdb_header(minimal_text)
  println("  PDB ID: '\{minimal.pdb_id}'")
  println("  Title: '\{minimal.title}'")
  println("  Compound name: '\{minimal.compound.compound_name}'")
  println("  Resolution: \{minimal.resolution}")
  println("  R-factor: \{minimal.rfactor}")
  println("  Chain IDs: \{minimal.chain_ids.length()}")

  println("")
  println("Malformed REMARK (no keyword):")
  let malformed = "REMARK   2 THIS LINE HAS NO RESOLUTION DATA\n"
  let mal_res = @bio.parse_resolution_record(malformed)
  println("  Resolution (expected 0.0): \{mal_res}")

  let malformed_r = "REMARK   3 NO RFACTOR HERE\n"
  let mal_rf = @bio.parse_rfactor_record(malformed_r)
  println("  R-factor (expected 0.0): \{mal_rf}")

  println("")
  println("=== Demo complete ===")
}