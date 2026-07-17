///|
/// SwissProt example: parse and analyze UniProt/Swiss-Prot records.
///
/// Run: moon run examples/swissprot_demo/main.mbt
fn main {
  println("=== SwissProt / UniProt Record Parsing Demo ===")
  println("")

  // === Example 1: Sample record ===
  println("--- 1. Sample SwissProt Record ---")

  let record = @bio.sample_swissprot_record()

  println("  Entry name: \{record.entry_name}")
  match record.primary_accession() {
    Some(acc) => println("  Primary accession: \{acc}")
    None => println("  Primary accession: (none)")
  }
  println("  Protein name: \{record.protein_name}")
  println("  Organism: \{record.organism} (TaxID: \{record.organism_id})")
  println("  Sequence length: \{record.seq_length} AA")
  println("  Molecular weight: \{record.seq_mass} Da")
  println("")

  // === Example 2: Gene names ===
  println("--- 2. Gene Names ---")
  let mut genes_str = "  Genes: "
  for i = 0; i < record.gene_names.length(); i = i + 1 {
    if i > 0 {
      genes_str = genes_str + ", "
    }
    genes_str = genes_str + record.gene_names[i]
  }
  println(genes_str)
  match record.primary_gene() {
    Some(g) => println("  Primary gene: \{g}")
    None => println("  Primary gene: (none)")
  }
  println("")

  // === Example 3: Taxonomy ===
  println("--- 3. Taxonomy ---")
  let mut tax_str = "  Lineage: "
  for i = 0; i < record.taxonomy.length(); i = i + 1 {
    if i > 0 {
      tax_str = tax_str + " > "
    }
    tax_str = tax_str + record.taxonomy[i]
  }
  println(tax_str)
  println("")

  // === Example 4: Keywords ===
  println("--- 4. Keywords ---")
  let mut kw_str = "  "
  for i = 0; i < record.keywords.length(); i = i + 1 {
    if i > 0 {
      kw_str = kw_str + ", "
    }
    kw_str = kw_str + record.keywords[i]
  }
  println(kw_str)
  println("")

  // === Example 5: EC numbers ===
  println("--- 5. EC Numbers ---")
  for i = 0; i < record.ec_numbers.length(); i = i + 1 {
    println("  EC \{i + 1}: \{record.ec_numbers[i]}")
  }
  println("")

  // === Example 6: Features ===
  println("--- 6. Features ---")
  println("  Total features: \{record.features.length()}")

  let domains = record.get_features("DOMAIN")
  println("  Domains: \{domains.length()}")
  for d in domains {
    println("    - \{d.description}: \{d.start}-\{d.end}")
  }

  let act_sites = record.get_features("ACT_SITE")
  println("  Active sites: \{act_sites.length()}")
  for s in act_sites {
    println("    - \{s.description}: residue \{s.start}")
  }

  let binding = record.get_features("BINDING")
  println("  Binding sites: \{binding.length()}")
  for b in binding {
    println("    - \{b.description}: \{b.start}-\{b.end}")
  }
  println("")

  // === Example 7: Function and localization ===
  println("--- 7. Function & Localization ---")
  println("  Function: \{record.function}")
  let mut loc_str = "  Subcellular location: "
  for i = 0; i < record.subcellular_location.length(); i = i + 1 {
    if i > 0 {
      loc_str = loc_str + ", "
    }
    loc_str = loc_str + record.subcellular_location[i]
  }
  println(loc_str)
  println("  Tissue specificity: \{record.tissue_specificity}")
  println("")

  // === Example 8: References ===
  println("--- 8. References ---")
  println("  Total references: \{record.references.length()}")
  for r in record.references {
    println("  Reference \{r.number}:")
    println("    Title: \{r.title}")
    let mut auth_str = "    Authors: "
    for i = 0; i < r.authors.length(); i = i + 1 {
      if i > 0 {
        auth_str = auth_str + ", "
      }
      auth_str = auth_str + r.authors[i]
    }
    println(auth_str)
    println("    Journal: \{r.journal}")
    println("    PubMed: \{r.pubmed_id}")
    println("    DOI: \{r.doi}")
  }
  println("")

  // === Example 9: Parse from text ===
  println("--- 9. Parse from UniProt Flat File ---")

  let sample_text = "ID   HBA_HUMAN              Standard;       142 AA.\nAC   P69905; P01922; Q6B0B5;\nDE   RecName: Full=Hemoglobin subunit alpha; EC=1.2.1.-;\nDE   AltName: Full=Alpha-globin;\nDE   AltName: Full=Hemoglobin alpha chain;\nGN   Name=HBA1; Synonyms=HBA, HBAA;\nOS   Homo sapiens (Human).\nOX   NCBI_TaxID=9606;\nOC   Eukaryota; Metazoa; Chordata; Craniata; Vertebrata; Euteleostomi;\nOC   Mammalia; Eutheria; Euarchontoglires; Primates; Haplorrhini;\nOC   Catarrhini; Hominoidea; Hominidae; Homo.\nKW   Acetylation; Alternative splicing; Blood group;\nKW   Cadmium; Developmental protein; Diacetylation;\nKW   Glycation; Hemoglobin; Iron; Metal-binding; Methylation;\nKW   Oxygen transport; Phosphoprotein; Transport.\nFT   HELIX           16      36\nFT   HELIX           42      57\nFT   HELIX           61      80\nFT   HELIX           85      96\nFT   HELIX          106     126\nFT   BINDING         58      63       Heme.\nFT   BINDING         87      87       Heme iron.\nFT   ACT_SITE        87      87       Proximal His.\nRN   [1]\nRX   PUBMED=13576275; DOI=10.1016/0006-291x(59)90510-7;\nRA   Schroeder W.A., Shelton J.R., Shelton J.B., Robberson B.,\nRA   Apell G.;\nRT   \"Some chemical and physical properties of human hemoglobin A2 and\nRT   of its subunits.\";\nRL   Biochim. Biophys. Acta 36:453-463(1959).\nCC   -!- FUNCTION: Involved in oxygen transport from the lung to the\nCC   various peripheral tissues.\nCC   -!- SUBCELLULAR LOCATION: Cytoplasm.\nCC   -!- TISSUE SPECIFICITY: Red blood cells.\nDT   01-NOV-1986, integrated into UniProtKB/Swiss-Prot.\nDT   23-JAN-2007, sequence version 3.\nDT   18-MAY-2024, entry version 226.\nSQ   SEQUENCE  142 AA;  15258 MW;  8F05A60E7D966B73 CRC64;\n     MVLSPADKTN VKAAWGKVGA HAGEYGAEAL ERMFLSFPTT KTYFPHFDL\n     SHGSAQVKGH GKKVADALTN AVAHVDDMPN ALSALSDLHA HKLRVDPVNF\n     KLLSHCLLVT LAAHLPAEFT PAVHASLDKFLASVSTVLTS KYR\n//\n"

  let parsed = @bio.parse_swissprot(sample_text)
  match parsed {
    Some(rec) => {
      println("  Parsed successfully!")
      println("  Entry: \{rec.entry_name}")
      println("  Accessions: \{rec.accession_numbers.length()} total")
      println("  Protein: \{rec.protein_name}")
      println("  Organism: \{rec.organism}")
      println("  Sequence length: \{rec.seq_length} AA")
      println("  Keywords: \{rec.keywords.length()}")
      println("  Features: \{rec.features.length()}")
      println("  References: \{rec.references.length()}")
      println("  Function: \{rec.function}")
    }
    None => println("  Failed to parse")
  }
  println("")

  println("=== Demo Complete ===")
}
