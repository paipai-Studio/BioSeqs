///|
/// NGS quality trimming demo.
///
/// Demonstrates parsing FASTQ data, quality trimming reads, adapter
/// removal, length and GC filtering, batch processing, and
/// trimming statistics computation.
///
/// Run: moon run examples/seq_quality_trim_demo/main.mbt
fn main {
  println("=== NGS Quality Trimming Demo ===")
  println("")

  // Sample FASTQ data with various quality profiles
  let fastq_data =
    "@read1 high_quality\nACGTACGTACGTACGT\n+\nIIIIIIIIIIIIIIII\n@read2 adapter_contaminated\nACGTACGTAGATCGGAA\n+\nIIIIIIIIIIIIIIIII\n@read3 low_quality_ends\nACGTACGTACGT\n+\nIIIIIIII####\n@read4 poly_A_tail\nACGTACGGGGAAAA\n+\nIIIIIIIIIIIIII\n@read5 short_read\nACGT\n+\nIIII\n@read6 extreme_gc\nGGGGGGGGGGGGCCCC\n+\nIIIIIIIIIIIIIIII\n"

  println("1. Parsing FASTQ data")
  let reads = @bio.sqt_fastq_parse(fastq_data)
  println("   Reads parsed: \{reads.length()}")
  let mut i = 0
  while i < reads.length() {
    let r = reads[i]
    println("     \{r.id}: seq_len=\{r.sequence.length()}, qual_len=\{r.quality.length()}")
    i = i + 1
  }

  println("")
  println("2. Quality trimming (min_quality=20)")
  println("   Trimming low-quality bases from both ends...")
  let seq = "ACGTACGTACGT"
  let qual = "IIIIIIII####"
  let (trimmed_seq, _) = @bio.sqt_quality_trim(seq, qual, min_qual=20)
  println("   Original seq:   \{seq} (len=\{seq.length()})")
  println("   Original qual:  \{qual}")
  println("   Trimmed seq:    \{trimmed_seq} (len=\{trimmed_seq.length()})")

  println("")
  println("3. Adapter removal")
  let seq_adapter = "ACGTACGTAGATCGGAA"
  let qual_adapter = "IIIIIIIIIIIIIIIII"
  let adapter = "AGATCGGAA"
  println("   Adapter: 5'-AGATCGGAA-3' (Illumina TruSeq adapter)")
  println("   Before:   \{seq_adapter}")
  let (clean_seq, _) = @bio.sqt_trim_adapter(
    seq_adapter,
    qual_adapter,
    adapter,
    allowed_mismatches=0,
  )
  println("   After:    \{clean_seq} (adapter removed)")

  println("")
  println("4. Poly-A / Poly-T trimming")
  let seq_polya = "ACGTACGGGGAAAA"
  let qual_polya = "IIIIIIIIIIIIII"
  println("   Before:  \{seq_polya}")
  let (poly_clean_seq, _) = @bio.sqt_trim_poly_a(seq_polya, qual_polya, min_tail=3)
  println("   After:   \{poly_clean_seq} (poly-A tail removed)")

  let seq_polyt = "TTTTACGTACGT"
  println("   Before:  \{seq_polyt}")
  let (polyt_clean_seq, _) = @bio.sqt_trim_poly_a(
    seq_polyt,
    "IIIIIIIIIIII",
    min_tail=3,
  )
  println("   After:   \{polyt_clean_seq} (poly-T tail removed)")

  println("")
  println("5. Filtering by length (min_length=8)")
  let short_seq = "ACGT"
  let short_qual = "IIII"
  let long_seq = "ACGTACGTACGT"
  let long_qual = "IIIIIIIIIIII"

  let short_result = @bio.sqt_trim_by_length(short_seq, short_qual, min_length=8)
  match short_result {
    Some(_) => println("   Short read (len=4): KEPT")
    None => println("   Short read (len=4): DISCARDED (below min_length)")
  }

  let long_result = @bio.sqt_trim_by_length(long_seq, long_qual, min_length=8)
  match long_result {
    Some(_) => println("   Long read (len=12): KEPT")
    None => println("   Long read (len=12): DISCARDED")
  }

  println("")
  println("6. Filtering by GC content (gc_min=0.3, gc_max=0.7)")
  let gc_high_seq = "GGGGGGGGGGGGCCCC"
  let gc_high_qual = "IIIIIIIIIIIIIIII"
  let gc_at_seq = "AAAAAAAAAATTTTTT"
  let gc_at_qual = "IIIIIIIIIIIIIIII"

  let gc_high_result = @bio.sqt_trim_by_gc(
    gc_high_seq,
    gc_high_qual,
    gc_min=0.3,
    gc_max=0.7,
  )
  match gc_high_result {
    Some(_) => println("   GC-high read: KEPT")
    None => println("   GC-high read (83% GC): DISCARDED (above gc_max)")
  }

  let gc_at_result = @bio.sqt_trim_by_gc(
    gc_at_seq,
    gc_at_qual,
    gc_min=0.3,
    gc_max=0.7,
  )
  match gc_at_result {
    Some(_) => println("   AT-rich read: KEPT")
    None => println("   AT-rich read (0% GC): DISCARDED (below gc_min)")
  }

  println("")
  println("7. Batch trimming pipeline")
  println("   Running full trimming pipeline on all reads...")
  let qual_config = @bio.QualityConfig::new(20, 4, 4, 20)
  let adapter_config = @bio.AdapterConfig::new(["AGATCGGAA"], 4, 0)

  let trim_results = @bio.sqt_trim_reads(reads, qual_config, Some(adapter_config))
  println("   Trim results:")
  i = 0
  while i < trim_results.length() {
    let tr = trim_results[i]
    let discard_note = if tr.trimmed_seq.length() == 0 { " [DISCARDED]" } else { "" }
    println(
      "     \{tr.read_id}: \{tr.original_seq} -> \{tr.trimmed_seq} (type=\{tr.trim_type})\{discard_note}",
    )
    i = i + 1
  }

  println("")
  println("8. Computing trimming statistics")
  let stats = @bio.sqt_compute_stats(trim_results, reads)
  println("   Total reads:    \{stats.total_reads}")
  println("   Kept reads:     \{stats.kept_reads}")
  println("   Trimmed reads:  \{stats.trimmed_reads}")
  println("   Discarded reads: \{stats.discarded_reads}")
  println("   Total bp before: \{stats.total_bp_before}")
  println("   Total bp after:  \{stats.total_bp_after}")
  println("   GC before:       \{stats.gc_content_before}")
  println("   GC after:        \{stats.gc_content_after}")

  println("")
  println("9. Serializing trimmed reads back to FASTQ")
  let kept_reads : Array[@bio.FastqRead] = []
  i = 0
  while i < trim_results.length() {
    let tr = trim_results[i]
    if tr.trimmed_seq.length() > 0 {
      kept_reads.push(@bio.FastqRead::new(
        tr.read_id,
        tr.trimmed_seq,
        tr.trimmed_qual,
      ))
    }
    i = i + 1
  }
  let fastq_out = @bio.sqt_fastq_serialize(kept_reads)
  println("   Output FASTQ (\{kept_reads.length()} reads):")
  let out_lines = fastq_out.split("\n").to_array()
  let mut l = 0
  while l < out_lines.length() && l < 12 {
    println("     \{out_lines[l]}")
    l = l + 1
  }
  if out_lines.length() > 12 {
    println("     ... (\{out_lines.length()} total lines)")
  }

  println("")
  println("=== Demo Complete ===")
}