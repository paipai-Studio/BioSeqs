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

  // === New: EMBL Format Parsing ===
  println("\n=== EMBL Format Parsing (NEW) ===")
  let embl_content =
    "ID   HSBGLOD; SV 1; linear; genomic DNA; STD; HUM; 20 BP.\n" +
    "XX\n" +
    "AC   M16496;\n" +
    "XX\n" +
    "DE   Human beta-globin gene region.\n" +
    "XX\n" +
    "SQ   Sequence 20 BP; 5 A; 5 C; 5 G; 5 T; 0 other;\n" +
    "     atgcatgcat gcatgcatgc\n" +
    "//\n"
  try {
    let embl_records = @bio.seqio_parse(embl_content, "embl")
    for i = 0; i < embl_records.length(); i = i + 1 {
      let r = embl_records[i]
      println("  [\{r.id}] name=\{r.name} len=\{r.length()}")
      println("    desc: \"\{r.description}\"")
      println("    seq:  \{r.seq.to_string()}")
    }
  } catch {
    _ => println("  EMBL parse error")
  }

  // === New: PIR/NBRF Format Parsing ===
  println("\n=== PIR/NBRF Format Parsing (NEW) ===")
  let pir_content =
    ">P1;ALB_HUMAN\n" +
    "Serum albumin precursor.\n" +
    "MKWVTFISLLFLFSSAYS*\n" +
    ">P1;HBA_HUMAN\n" +
    "Hemoglobin subunit alpha.\n" +
    "MVLSPADKTNVKAAWGKVGA*"
  try {
    let pir_records = @bio.seqio_parse(pir_content, "pir")
    println("  Parsed PIR records: \{pir_records.length()}")
    for i = 0; i < pir_records.length(); i = i + 1 {
      let r = pir_records[i]
      println("  [\{r.id}] len=\{r.length()}")
      println("    desc: \"\{r.description}\"")
    }
  } catch {
    _ => println("  PIR parse error")
  }

  // === New: Tab Format Read/Write ===
  println("\n=== Tab-delimited Format (NEW) ===")
  let tab_content = "gene1\tATGCATGC\ngene2\tGGGGCCCC\ngene3\tTTTTAAAA"
  try {
    let tab_records = @bio.seqio_parse(tab_content, "tab")
    println("  Parsed Tab records: \{tab_records.length()}")
    for i = 0; i < tab_records.length(); i = i + 1 {
      println("    \{tab_records[i].id}\t<\{tab_records[i].seq.to_string()}>")
    }
    // Write Tab format back
    let tab_back = @bio.seqio_write(tab_records, "tab")
    println("\n  Tab-format output:")
    println("\{tab_back}")
  } catch {
    _ => println("  Tab parse error")
  }

  // === New: GenBank Format Writing ===
  println("\n=== GenBank Format Writing (NEW) ===")
  try {
    let write_test = [
      @bio.SeqRecord::new(
        @bio.Seq::new("ATGCATGCATGCATGCATGC"),
        id="GB001",
        name="DEMO",
        description="Demo record for GenBank writer",
      ),
    ]
    let gb_text = @bio.seqio_write(write_test, "genbank")
    println("\{gb_text}")
  } catch {
    _ => println("  GenBank write error")
  }

  // === New: seqio_convert Format Conversion ===
  println("\n=== seqio_convert: FASTA -> Tab (NEW) ===")
  try {
    let converted = @bio.seqio_convert(fasta_content, "fasta", "tab")
    println("\{converted}")
  } catch {
    _ => println("  Format conversion error")
  }

  println("\n=== seqio_convert: Tab -> FASTA (NEW) ===")
  try {
    let conv2 = @bio.seqio_convert("x1\tAAAA\nx2\tCCCC", "tab", "fasta")
    println("\{conv2}")
  } catch {
    _ => println("  Conversion error 2")
  }

  println("\n=== seqio_convert: GenBank Round-trip (NEW) ===")
  try {
    // FASTA -> GenBank(write) -> GenBank(parse) -> FASTA
    let recs1 = @bio.seqio_parse(fasta_content, "fasta")
    let as_gb = @bio.seqio_write(recs1, "genbank")
    let roundtrip = @bio.seqio_parse(as_gb, "genbank")
    println("  Original: \{recs1.length()} records")
    println("  After round-trip: \{roundtrip.length()} records")
    if roundtrip.length() > 0 {
      println("  First record seq: \{roundtrip[0].seq.to_string()}")
    }
  } catch {
    _ => println("  Round-trip error occurred")
  }

  // === New: seqio_read (single record) ===
  println("\n=== seqio_read: single record (NEW) ===")
  try {
    let only_one = @bio.seqio_read("solo\tACGTACGT", "tab")
    println("  id: \{only_one.id} seq: \{only_one.seq.to_string()}")
  } catch {
    _ => println("  seqio_read error")
  }
}
