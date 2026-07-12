///|
/// Sequence I/O example: parse & write FASTA, FASTQ, GenBank.
///
/// Run: moon run examples/seqio_demo/main.mbt
fn main {
  // === Parse FASTA ===
  let fasta_content = ">seq1 Description for seq1\nATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG\n>seq2 Another sequence\nACGTACGTACGTACGT\n>seq3\nATATATATATAT"

  let records = @bio.parse_fasta(fasta_content)
  println("=== FASTA Parsing ===")
  for i = 0; i < records.length(); i = i + 1 {
    let r = records[i]
    println("  [\{r.id}] len=\{r.length()} desc=\"\{r.description}\"")
    println("    seq: \{r.seq.to_string()}")
  }

  // === Write FASTA ===
  let fasta_out = @bio.write_fasta(records)
  println("\n=== FASTA Output ===")
  println(fasta_out)

  // === Parse FASTQ ===
  let fastq_content = "@seq1 Description for seq1\nATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG\n+\n" +
    "I".repeat(39) +
    "\n@seq2 Another read\nACGTACGT\n+\n" +
    "I".repeat(8)

  try {
    let fq_records = @bio.parse_fastq(fastq_content)
    println("=== FASTQ Parsing ===")
    for i = 0; i < fq_records.length(); i = i + 1 {
      let r = fq_records[i]
      let qual = r.phred_quality()
      match qual {
        Some(q) =>
          println("  [\{r.id}] len=\{r.length()} qual_len=\{q.length()}")
        None => println("  [\{r.id}] len=\{r.length()} no_quality")
      }
    }

    // === Write FASTQ ===
    let fastq_out = @bio.write_fastq(fq_records)
    println("\n=== FASTQ Output ===")
    println(fastq_out)
  } catch {
    _ => println("FASTQ error occurred")
  }

  // === SeqIO unified interface ===
  try {
    let unified = @bio.seqio_parse(fasta_content, "fasta")
    println("=== SeqIO parse (fasta) ===")
    println("  Records: \{unified.length()}")
    for i = 0; i < unified.length(); i = i + 1 {
      println("  [\{unified[i].id}] \{unified[i].seq.to_string()}")
    }
  } catch {
    _ => println("SeqIO error occurred")
  }

  // === SeqIO to_dict ===
  let dict = @bio.seqio_to_dict(records)
  println("\n=== SeqIO to_dict ===")
  let keys = dict.keys()
  println("  keys: \{keys}")

  // === Parse GenBank ===
  let genbank_content = "LOCUS       NM_001 48 bp DNA linear PRI 01-JAN-2024\nDEFINITION  Test gene for comparison.\nACCESSION   NM_001\nVERSION     NM_001.1\nORIGIN\nATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGC\n//"

  let gb_records = @bio.parse_genbank(genbank_content)
  println("\n=== GenBank Parsing ===")
  for i = 0; i < gb_records.length(); i = i + 1 {
    let r = gb_records[i]
    println("  [\{r.id}] len=\{r.length()} desc=\"\{r.description}\"")
    println("    seq: \{r.seq.to_string()}")
  }

  // === FASTA indexing (pyfaidx style) ===
  let fasta = @bio.Fasta::from_content(fasta_content)
  println("\n=== FASTA Index ===")
  let names = fasta.get_names()
  println("  Indexed sequences: \{names}")

  // Fetch subsequence: seq1[3..15]
  try {
    let sub = fasta.fetch("seq1", 3, 15)
    match sub {
      Some(s) => println("  seq1[3:15]: \{s.to_string()}")
      None => println("  seq1[3:15]: not found")
    }
  } catch {
    _ => println("  fetch error occurred")
  }

  // Get full sequence
  let full = fasta.get_seq("seq2")
  match full {
    Some(s) => println("  seq2 full: \{s.to_string()}")
    None => println("  seq2 full: not found")
  }
}
