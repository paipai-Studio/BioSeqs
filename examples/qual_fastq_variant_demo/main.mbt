///|
/// QUAL / FASTQ variant example: QUAL format, FASTQ variant encodings,
/// quality score conversion, and PairedFastaQualIterator.
///
/// Run: moon run examples/qual_fastq_variant_demo/main.mbt
fn main {
  // === QUAL Format Parsing ===
  println("=== QUAL Format Parsing ===")
  let qual_content = ">seq1\n10 20 30 40\n>seq2\n50 60 70 80\n"
  let qual_records = @bio.parse_qual(qual_content)
  println("  Parsed \{qual_records.length()} QUAL records")
  for i = 0; i < qual_records.length(); i = i + 1 {
    let r = qual_records[i]
    let quals = r.letter_annotations["phred_quality"]
    println("  [\{r.id}] len=\{r.length()} quals=" + array_to_string(quals))
  }

  // === QUAL Writing (round-trip) ===
  println("\n=== QUAL Write (round-trip) ===")
  try {
    let qual_out = @bio.write_qual(qual_records)
    println("  Output:")
    println(qual_out)
  } catch {
    _ => println("  QUAL write error")
  }

  // === FASTQ Illumina (offset 64) ===
  println("=== FASTQ Illumina (offset 64) ===")
  // '@' = 64 (ASCII), so '@' - 64 = 0, 'A'=65-64=1, 'B'=66-64=2, 'C'=67-64=3
  let illumina_fastq = "@seq1\nACGT\n+\n@ABC\n"
  try {
    let illumina_records = @bio.parse_fastq_illumina(illumina_fastq)
    println("  Parsed \{illumina_records.length()} Illumina FASTQ records")
    for i = 0; i < illumina_records.length(); i = i + 1 {
      let r = illumina_records[i]
      let quals = r.letter_annotations["phred_quality"]
      println("  [\{r.id}] seq=\{r.seq.to_string()} quals=" + array_to_string(quals))
    }

    // Write back to Illumina FASTQ
    let illumina_out = @bio.write_fastq_illumina(illumina_records)
    println("  Written back:")
    println(illumina_out)
  } catch {
    _ => println("  Illumina FASTQ error occurred")
  }

  // === FASTQ Solexa (Solexa scores, offset 64) ===
  println("=== FASTQ Solexa (Solexa scores, offset 64) ===")
  // ';' = 59 (ASCII), so ';' - 64 = -5 (minimum Solexa score)
  let solexa_fastq = "@seq1\nACGT\n+\n;>?@\n"
  try {
    let solexa_records = @bio.parse_fastq_solexa(solexa_fastq)
    println("  Parsed \{solexa_records.length()} Solexa FASTQ records")
    for i = 0; i < solexa_records.length(); i = i + 1 {
      let r = solexa_records[i]
      let quals = r.letter_annotations["phred_quality"]
      println("  [\{r.id}] seq=\{r.seq.to_string()} quals=" + array_to_string(quals))
    }

    // Write back to Solexa FASTQ
    let solexa_out = @bio.write_fastq_solexa(solexa_records)
    println("  Written back:")
    println(solexa_out)
  } catch {
    _ => println("  Solexa FASTQ error occurred")
  }

  // === Quality Score Conversion ===
  println("=== Quality Score Conversion ===")
  let test_phred = [40.0, 30.0, 20.0, 10.0, 5.0, 3.0, 2.0, 1.0, 0.0]
  println("  Phred -> Solexa -> Phred:")
  for q in test_phred {
    let solexa = @bio.solexa_quality_from_phred(q)
    let back = @bio.phred_quality_from_solexa(solexa)
    let q_str = double_to_str(q)
    let s_str = double_to_str(solexa)
    let b_str = double_to_str(back)
    println("  " + q_str + " -> " + s_str + " -> " + b_str)
  }

  // === PairedFastaQualIterator ===
  println("\n=== PairedFastaQualIterator ===")
  let fasta_content = ">seq1\nACGT\n>seq2\nTGCA\n"
  let qual_for_pair = ">seq1\n10 20 30 40\n>seq2\n50 60 70 80\n"
  let paired_records = @bio.paired_fasta_qual_iterator(fasta_content, qual_for_pair)
  println("  Paired \{paired_records.length()} records")
  for i = 0; i < paired_records.length(); i = i + 1 {
    let r = paired_records[i]
    let quals = r.letter_annotations["phred_quality"]
    println("  [\{r.id}] seq=\{r.seq.to_string()} quals=" + array_to_string(quals))
  }

  // === SeqIO Unified Interface (QUAL, FASTQ variants) ===
  println("\n=== SeqIO Unified Interface ===")
  try {
    let qual_via_seqio = @bio.seqio_parse(qual_content, "qual")
    println("  seqio_parse(\"qual\"): \{qual_via_seqio.length()} records")
  } catch {
    _ => println("  seqio_parse(\"qual\") error")
  }

  try {
    let illumina_via_seqio = @bio.seqio_parse(illumina_fastq, "fastq-illumina")
    println("  seqio_parse(\"fastq-illumina\"): \{illumina_via_seqio.length()} records")
  } catch {
    _ => println("  seqio_parse(\"fastq-illumina\") error")
  }

  try {
    let solexa_via_seqio = @bio.seqio_parse(solexa_fastq, "fastq-solexa")
    println("  seqio_parse(\"fastq-solexa\"): \{solexa_via_seqio.length()} records")
  } catch {
    _ => println("  seqio_parse(\"fastq-solexa\") error")
  }

  // === SeqIO convert: FASTQ -> QUAL ===
  println("\n=== SeqIO convert: FASTQ -> QUAL ===")
  let fastq_sanger = "@seq1\nACGT\n+\n!\"#$\n"
  try {
    let fastq_to_qual = @bio.seqio_convert(fastq_sanger, "fastq", "qual")
    println("  Converted FASTQ to QUAL:")
    println(fastq_to_qual)
  } catch {
    _ => println("  FASTQ -> QUAL conversion error")
  }

  // === SeqIO convert: QUAL -> FASTQ ===
  println("\n=== SeqIO convert: QUAL -> FASTQ ===")
  try {
    let qual_to_fastq = @bio.seqio_convert(qual_content, "qual", "fastq")
    println("  Converted QUAL to FASTQ:")
    println(qual_to_fastq)
  } catch {
    _ => println("  QUAL -> FASTQ conversion error")
  }
}

///|
/// Helper: convert an Array[Int] to a display string like "[10, 20, 30]"
fn array_to_string(arr : Array[Int]) -> String {
  let buf = StringBuilder()
  buf.write_char('[')
  for i = 0; i < arr.length(); i = i + 1 {
    if i > 0 {
      buf.write_string(", ")
    }
    buf.write_string(arr[i].to_string())
  }
  buf.write_char(']')
  buf.to_string()
}

///|
/// Helper: format a Double to 4 decimal places
fn double_to_str(v : Double) -> String {
  let s = v.to_string()
  // Find decimal point and truncate to 4 places
  match s.find(".") {
    Some(pos) => {
      let end = (pos + 5).min(s.length())
      if end < s.length() {
        s[0:end].to_owned()
      } else {
        s
      }
    }
    None => s
  }
}