///|
/// IMGT module demo - ImMunoGeneTics HLA sequence format parser and writer.
///
/// Run: moon run examples/imgt_demo
fn main {
  println("=== IMGT Module Demo ===")
  println("ImMunoGeneTics information system HLA sequence format")
  println("Reference: Robinson J et al. (2015) Nucleic Acids Res 43:D423-31")
  println("")

  // ------------------------------------------------------------------
  // 1. Load sample IMGT records
  // ------------------------------------------------------------------
  println("=== 1. Load sample IMGT records ===")
  let records = @bio.imgt_sample_records()
  println("Loaded \{records.length()} sample IMGT record(s)")
  println("")

  // ------------------------------------------------------------------
  // 2. Print summary (n_records, unique species, unique genes)
  // ------------------------------------------------------------------
  println("=== 2. Summary ===")
  let n_records = @bio.imgt_n_records(records)
  let species = @bio.imgt_unique_species(records)
  let genes = @bio.imgt_unique_genes(records)
  println("n_records: \{n_records}")
  println("unique species: \{species.length()}")
  for i in 0..<species.length() {
    println("  - \{species[i]}")
  }
  println("unique genes: \{genes.length()}")
  for i in 0..<genes.length() {
    println("  - \{genes[i]}")
  }
  println("")
  println("Built-in summary:")
  println(@bio.imgt_summary(records))
  println("")

  // ------------------------------------------------------------------
  // 3. Show all records (accession, seq_id, gene, species, sequence length)
  // ------------------------------------------------------------------
  println("=== 3. All records ===")
  for i in 0..<records.length() {
    let rec = records[i]
    let h = rec.header()
    println("Record \{i}:")
    println("  accession: \{h.accession()}")
    println("  seq_id: \{h.seq_id()}")
    println("  gene: \{h.gene()}")
    println("  species: \{h.species()}")
    println("  locus: \{h.locus()}")
    println("  seq_length (header): \{h.seq_length()}")
    println("  actual sequence length: \{rec.sequence().length()}")
    println("  to_string: \{rec.to_string()}")
  }
  println("")

  // ------------------------------------------------------------------
  // 4. Filter by gene (HLA-A) and show results
  // ------------------------------------------------------------------
  println("=== 4. Filter by gene (HLA-A) ===")
  let hla_a = @bio.imgt_get_by_gene(records, "HLA-A")
  println("Found \{hla_a.length()} record(s) with gene=HLA-A")
  for i in 0..<hla_a.length() {
    let rec = hla_a[i]
    let h = rec.header()
    println("  \{i}: accession=\{h.accession()}, seq_id=\{h.seq_id()}, length=\{rec.sequence().length()}")
  }
  println("")
  // Also demonstrate filtering by species.
  println("Filter by species (human):")
  let human = @bio.imgt_get_by_species(records, "human")
  println("Found \{human.length()} record(s) with species=human")
  println("")

  // ------------------------------------------------------------------
  // 5. Parse an IMGT string and verify
  // ------------------------------------------------------------------
  println("=== 5. Parse an IMGT string ===")
  let imgt_string = ">HLA99999|C*01:02:01|C*01:02:01|ORF|2020-01-01|2000-01-01|human|HLA-C|I|MHC|2730 bp| genomic-DNA|test record\nATGGCGCCCCGAACCCTCCTCCTG\n"
  println("Input IMGT string:")
  println(imgt_string)
  let parsed = @bio.imgt_parse(imgt_string)
  println("Parsed \{parsed.length()} record(s)")
  if parsed.length() > 0 {
    let p = parsed[0]
    let h = p.header()
    println("  accession: \{h.accession()}")
    println("  seq_id: \{h.seq_id()}")
    println("  gene: \{h.gene()}")
    println("  species: \{h.species()}")
    println("  sequence: \{p.sequence()}")
    println("  sequence length: \{p.sequence().length()}")
  }
  println("")

  // ------------------------------------------------------------------
  // 6. Write records to IMGT format and show output
  // ------------------------------------------------------------------
  println("=== 6. Write records to IMGT format ===")
  let output = @bio.imgt_write(records)
  println("Serialized \{records.length()} record(s) to IMGT text:")
  // Print each line of the output.
  let lines = output.split("\n").to_array()
  for i in 0..<lines.length() {
    if lines[i].length() > 0 {
      println("  |\{lines[i]}")
    }
  }
  println("")
  // Show a single record serialization as well.
  println("Single record serialization (record 0):")
  println(@bio.imgt_record_to_string(records[0]))
  println("")

  // ------------------------------------------------------------------
  // 7. Convert to SeqRecord array and show
  // ------------------------------------------------------------------
  println("=== 7. Convert to SeqRecord array ===")
  let seq_records = @bio.imgt_to_seq_records(records)
  println("Converted \{seq_records.length()} ImgtRecord(s) to SeqRecord(s)")
  for i in 0..<seq_records.length() {
    let rec = seq_records[i]
    println("  Record \{i}: id=\{rec.id}, name=\{rec.name}, length=\{rec.seq.to_string().length()}")
    println("    description (truncated): \{rec.description}")
  }
  println("")

  // ------------------------------------------------------------------
  // 8. Search by accession number
  // ------------------------------------------------------------------
  println("=== 8. Search by accession number ===")
  // Search for an accession that exists.
  let acc_to_find = "HLA00002"
  match @bio.imgt_get_by_accession(records, acc_to_find) {
    Some(rec) => {
      let h = rec.header()
      println("Found accession '\{acc_to_find}':")
      println("  seq_id: \{h.seq_id()}")
      println("  gene: \{h.gene()}")
      println("  species: \{h.species()}")
      println("  sequence length: \{rec.sequence().length()}")
    }
    None => println("Accession '\{acc_to_find}' not found")
  }
  println("")
  // Search for an accession that does NOT exist.
  let missing_acc = "HLA99999"
  match @bio.imgt_get_by_accession(records, missing_acc) {
    Some(_) => println("Unexpected: found accession '\{missing_acc}'")
    None => println("Accession '\{missing_acc}' not found (expected)")
  }
  println("")

  println("=== Demo complete ===")
}
