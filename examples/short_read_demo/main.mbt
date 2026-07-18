///|
fn main {
  println("=== ShortRead Demo ===")
  
  let sequences = [
    "ACGTACGTAGCTAGCTAGATCGGAAGAGC",
    "ACGTACGTAGCTAGCTAGATCGGAAGAGC",
    "GGGGCCCCAAAATTTTGGGGCCCCAAAA",
    "ACGTNNNNAGCTAGCTAGATCGGAAGAGC",
    "TTTTCCCCGGGGAAAATTTTCCCCGGGG",
  ]
  let qualities = [
    "IIIIIIIIIIIIIIIIIIIIIIIIIIII",
    "IIIIIIIIIIIIIIIIIIIIIIIIIIII",
    "IIIIIIIIIIIIIIIIIIIIIIIIIIII",
    "IIIIIIII########IIIIIIIIIIII",
    "IIIIIIIIIIIIIIIIIIIIIIIIIIII",
  ]
  
  println("\n1. Quality Assessment Summary:")
  let stats = @src.qa_summary(sequences, qualities)
  println(stats.to_string())
  
  println("\n2. Per-base quality (first 10 cycles):")
  let mut i = 0
  while i < 10 && i < stats.quality_by_cycle.length() {
    println("  Cycle " + (i + 1).to_string() + ": " + stats.quality_by_cycle[i].to_string())
    i = i + 1
  }
  
  println("\n3. Trim tail (quality threshold 20):")
  let seq = "ACGTACGT##"
  let qual = "IIIIIIII##"
  let (trimmed_seq, trimmed_qual) = @src.trim_tail(seq, qual, 20)
  println("  Original: " + seq + " (" + seq.length().to_string() + "bp)")
  println("  Trimmed:  " + trimmed_seq + " (" + trimmed_seq.length().to_string() + "bp)")
  
  println("\n4. Trim leading (quality threshold 20):")
  let seq2 = "##ACGTACGT"
  let qual2 = "##IIIIIIII"
  let (trimmed_seq2, trimmed_qual2) = @src.trim_leading(seq2, qual2, 20)
  println("  Original: " + seq2 + " (" + seq2.length().to_string() + "bp)")
  println("  Trimmed:  " + trimmed_seq2 + " (" + trimmed_seq2.length().to_string() + "bp)")
  
  println("\n5. Trim adapter (AGATCGGAAGAGC):")
  let adapter = "AGATCGGAAGAGC"
  let (trimmed_adapter, _) = @src.trim_adapter(sequences[0], qualities[0], adapter, 0)
  println("  Original: " + sequences[0] + " (" + sequences[0].length().to_string() + "bp)")
  println("  Trimmed:  " + trimmed_adapter + " (" + trimmed_adapter.length().to_string() + "bp)")
  
  println("\n6. Filter reads (min quality 30, min length 20, max N 2):")
  let (filtered_seqs, filtered_quals) = @src.filter_reads(sequences, qualities, 30.0, 20, 2)
  println("  Input reads:  " + sequences.length().to_string())
  println("  Filtered reads: " + filtered_seqs.length().to_string())
  
  println("\n7. FastQC Report:")
  let report = @src.generate_fastqc_report(sequences, qualities)
  println(report.to_string())
  
  println("\nShortRead demo completed!")
}
