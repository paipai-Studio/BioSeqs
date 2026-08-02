///|
/// Exonerate output parsing demo - demonstrates parsing vulgar and cigar
/// output formats from the Exonerate spliced alignment tool.
fn main {
  println("=== Exonerate Parsing Demo ===\n")

  println("1. Parsing vulgar format:")
  let vulgar_data = @bio.exonerate_sample_data()
  let vulgar_records = @bio.parse_exonerate_vulgar(vulgar_data)
  println("   Parsed \{vulgar_records.length()} vulgar records")
  let mut i = 0
  while i < vulgar_records.length() {
    let rec = vulgar_records[i]
    println("   [\{i + 1}] \{rec.query_name()} \{rec.query_start()}-\{rec.query_end()}\{rec.query_strand()} -> \{rec.target_name()} \{rec.target_start()}-\{rec.target_end()}\{rec.target_strand()} score=\{rec.score().to_int()}")
    let mut j = 0
    while j < rec.alignment_blocks().length() {
      let b = rec.alignment_blocks()[j]
      println("       block \{j}: code=\{b.code()} q=\{b.query_length()} t=\{b.target_length()} score=\{b.block_score().to_int()}")
      j = j + 1
    }
    i = i + 1
  }

  println("\n2. Parsing cigar format:")
  let cigar_data = @bio.exonerate_sample_cigar_data()
  let cigar_records = @bio.parse_exonerate_cigar(cigar_data)
  println("   Parsed \{cigar_records.length()} cigar records")
  i = 0
  while i < cigar_records.length() {
    let rec = cigar_records[i]
    println("   [\{i + 1}] \{rec.query_name()} -> \{rec.target_name()} score=\{rec.score().to_int()} blocks=\{rec.alignment_blocks().length()}")
    i = i + 1
  }

  println("\n3. Auto-detecting format (mixed input):")
  let mixed = "vulgar: mixed_q 0 100 + mixed_t 0 95 + 1850 M 80 80 80 I 0 0 2000 M 20 15 20\n" +
    "cigar: mixed_q2 0 50 + mixed_t2 0 50 + 850 M 50\n"
  let mixed_records = @bio.parse_exonerate(mixed)
  println("   Parsed \{mixed_records.length()} records from mixed input")
  i = 0
  while i < mixed_records.length() {
    println("   [\{i + 1}] \{mixed_records[i].query_name()}")
    i = i + 1
  }

  println("\n4. Filtering by score (>= 1000):")
  let filtered = @bio.exonerate_filter_by_score(vulgar_records, 1000.0)
  println("   \{filtered.length()} of \{vulgar_records.length()} records passed")
  i = 0
  while i < filtered.length() {
    println("   [\{i + 1}] \{filtered[i].query_name()} score=\{filtered[i].score().to_int()}")
    i = i + 1
  }

  println("\n5. Per-record alignment statistics:")
  i = 0
  while i < vulgar_records.length() {
    let rec = vulgar_records[i]
    let aln_len = @bio.exonerate_total_alignment_length(rec)
    let intron_len = @bio.exonerate_total_intron_length(rec)
    let n_introns = @bio.exonerate_count_introns(rec)
    let n_matches = @bio.exonerate_count_matches(rec)
    println("   \{rec.query_name()}: matches=\{n_matches} aln_len=\{aln_len} introns=\{n_introns} intron_len=\{intron_len}")
    i = i + 1
  }

  println("\n6. Reconstructing vulgar and cigar strings:")
  let first = vulgar_records[0]
  println("   vulgar: \{@bio.exonerate_to_vulgar_string(first)}")
  println("   cigar:  \{@bio.exonerate_to_cigar_string(first)}")

  println("\n7. Summary:")
  let summary = @bio.exonerate_summary(vulgar_records)
  println(summary)

  println("=== Demo Complete ===")
}
