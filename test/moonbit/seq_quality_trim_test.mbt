///|
/// SeqQualityTrim module tests - quality trimming, adapter removal,
/// poly-A trimming, length/GC filtering, batch trimming, FASTQ I/O,
/// and statistics computation.

///|
test "sqt_quality_trim_3_prime_basic" {
  let seq = "ACGTACGT"
  let qual = "IIIIII##"
  let (t_seq, t_qual) = @src.sqt_quality_trim(seq, qual, min_qual=20)
  assert_eq(t_seq, "ACGTAC")
  assert_eq(t_seq.length(), 6)
  assert_eq(t_qual.length(), 6)
}

///|
test "sqt_quality_trim_5_prime_low_quality" {
  let seq = "ACGTACGT"
  let qual = "##IIIIII"
  let (t_seq, _) = @src.sqt_quality_trim(seq, qual, min_qual=20)
  assert_eq(t_seq, "GTACGT")
  assert_eq(t_seq.length(), 6)
}

///|
test "sqt_quality_trim_both_ends_low" {
  let seq = "ACGTACGTACGT"
  let qual = "##IIIIIIII##"
  let (t_seq, _) = @src.sqt_quality_trim(seq, qual, min_qual=20)
  assert_eq(t_seq, "GTACGTAC")
  assert_eq(t_seq.length(), 8)
}

///|
test "sqt_quality_trim_all_high_quality" {
  let seq = "ACGTACGT"
  let qual = "IIIIIIII"
  let (t_seq, _) = @src.sqt_quality_trim(seq, qual, min_qual=20)
  assert_eq(t_seq, "ACGTACGT")
  assert_eq(t_seq.length(), 8)
}

///|
test "sqt_quality_trim_all_low_quality" {
  let seq = "ACGTACGT"
  let qual = "########"
  let (t_seq, _) = @src.sqt_quality_trim(seq, qual, min_qual=20)
  assert_eq(t_seq.length(), 0)
}

///|
test "sqt_quality_trim_custom_threshold" {
  let seq = "ACGTACGT"
  let qual = "IIIIIIII"
  let (t_seq, _) = @src.sqt_quality_trim(seq, qual, min_qual=10)
  assert_eq(t_seq, "ACGTACGT")
}

///|
test "sqt_quality_trim_high_threshold" {
  let seq = "ACGTACGT"
  let qual = "IIIIIIII"
  let (t_seq, _) = @src.sqt_quality_trim(seq, qual, min_qual=41)
  assert_eq(t_seq.length(), 0)
}

///|
test "sqt_quality_trim_empty_sequence" {
  let seq = ""
  let qual = ""
  let (t_seq, t_qual) = @src.sqt_quality_trim(seq, qual, min_qual=20)
  assert_eq(t_seq.length(), 0)
  assert_eq(t_qual.length(), 0)
}

///|
test "sqt_quality_trim_quality_string_preserved" {
  let seq = "ACGTACGTACGT"
  let qual = "IIIIIIIIII##"
  let (t_seq, t_qual) = @src.sqt_quality_trim(seq, qual, min_qual=20)
  assert_eq(t_seq.length(), t_qual.length())
  assert_eq(t_qual, "IIIIIIIIII")
}

///|
test "sqt_quality_trim_single_base" {
  let seq = "A"
  let qual = "I"
  let (t_seq, _) = @src.sqt_quality_trim(seq, qual, min_qual=20)
  assert_eq(t_seq, "A")
}

///|
test "sqt_quality_trim_mixed_quality_pattern" {
  let seq = "ACGTACGTACGT"
  let qual = "II##IIII##II"
  let (t_seq, t_qual) = @src.sqt_quality_trim(seq, qual, min_qual=20)
  assert_eq(t_seq.length(), t_qual.length())
  assert_true(t_seq.length() > 0)
}

///|
test "sqt_quality_trim_no_trim_needed" {
  let seq = "GGCCGGCC"
  let qual = "IIIIIIII"
  let (t_seq, _) = @src.sqt_quality_trim(seq, qual, min_qual=30)
  assert_eq(t_seq, "GGCCGGCC")
}

///|
test "sqt_trim_adapter_perfect_match" {
  let seq = "ACGTACGTAGATCGGAA"
  let qual = "IIIIIIIIIIIIIIII"
  let adapter = "AGATCGGAA"
  let (t_seq, _) = @src.sqt_trim_adapter(
    seq, qual, adapter, allowed_mismatches=0,
  )
  assert_eq(t_seq, "ACGTACGT")
  assert_eq(t_seq.length(), 8)
}

///|
test "sqt_trim_adapter_with_mismatches" {
  let seq = "ACGTACGTAGATCGGAA"
  let qual = "IIIIIIIIIIIIIIII"
  let adapter = "AGATCGGAA"
  let (t_seq, _) = @src.sqt_trim_adapter(
    seq, qual, adapter, allowed_mismatches=2,
  )
  assert_eq(t_seq, "ACGTACGT")
}

///|
test "sqt_trim_adapter_no_match" {
  let seq = "ACGTACGTACGTACGT"
  let qual = "IIIIIIIIIIIIIIII"
  let adapter = "GGGGGGGG"
  let (t_seq, _) = @src.sqt_trim_adapter(
    seq, qual, adapter, allowed_mismatches=0,
  )
  assert_eq(t_seq, "ACGTACGTACGTACGT")
  assert_eq(t_seq.length(), 16)
}

///|
test "sqt_trim_adapter_empty_adapter" {
  let seq = "ACGTACGT"
  let qual = "IIIIIIII"
  let (t_seq, _) = @src.sqt_trim_adapter(
    seq, qual, "", allowed_mismatches=0,
  )
  assert_eq(t_seq, "ACGTACGT")
}

///|
test "sqt_trim_adapter_short_sequence" {
  let seq = "AC"
  let qual = "II"
  let adapter = "AGATCGGAA"
  let (t_seq, _) = @src.sqt_trim_adapter(
    seq, qual, adapter, allowed_mismatches=0,
  )
  assert_eq(t_seq, "AC")
}

///|
test "sqt_trim_adapter_at_beginning" {
  let seq = "AGATCGGAAACGTACGT"
  let qual = "IIIIIIIIIIIIIIII"
  let adapter = "AGATCGGAA"
  let (t_seq, _) = @src.sqt_trim_adapter(
    seq, qual, adapter, allowed_mismatches=0,
  )
  assert_eq(t_seq.length(), 0)
}

///|
test "sqt_trim_adapter_partial_overlap" {
  let seq = "ACGTACGTAGATCGGAA"
  let qual = "IIIIIIIIIIIIIIII"
  let adapter = "AGATCGGAA"
  let (t_seq, _) = @src.sqt_trim_adapter(
    seq, qual, adapter, allowed_mismatches=1,
  )
  assert_eq(t_seq, "ACGTACGT")
}

///|
test "sqt_trim_adapter_at_3_prime_end" {
  let seq = "ACGTACGTACGTAGATCGGAA"
  let qual = "IIIIIIIIIIIIIIIIIIIII"
  let adapter = "AGATCGGAA"
  let (t_seq, _) = @src.sqt_trim_adapter(
    seq, qual, adapter, allowed_mismatches=0,
  )
  assert_eq(t_seq, "ACGTACGTACGT")
}

///|
test "sqt_trim_poly_a_3_prime" {
  let seq = "ACGTACGGGAAAA"
  let qual = "IIIIIIIIIIIIII"
  let (t_seq, _) = @src.sqt_trim_poly_a(seq, qual, min_tail=3)
  assert_eq(t_seq, "ACGTACGGG")
}

///|
test "sqt_trim_poly_t_3_prime" {
  let seq = "ACGTACGGGTTTT"
  let qual = "IIIIIIIIIIIIII"
  let (t_seq, _) = @src.sqt_trim_poly_a(seq, qual, min_tail=3)
  assert_eq(t_seq, "ACGTACGGG")
}

///|
test "sqt_trim_poly_a_5_prime" {
  let seq = "AAAACGTACGT"
  let qual = "IIIIIIIIIII"
  let (t_seq, _) = @src.sqt_trim_poly_a(seq, qual, min_tail=3)
  assert_eq(t_seq, "CGTACGT")
}

///|
test "sqt_trim_poly_a_no_tail" {
  let seq = "ACGTACGT"
  let qual = "IIIIIIII"
  let (t_seq, _) = @src.sqt_trim_poly_a(seq, qual, min_tail=3)
  assert_eq(t_seq, "ACGTACGT")
}

///|
test "sqt_trim_poly_a_all_poly" {
  let seq = "AAAAAAAA"
  let qual = "IIIIIIII"
  let (t_seq, _) = @src.sqt_trim_poly_a(seq, qual, min_tail=3)
  assert_eq(t_seq.length(), 0)
}

///|
test "sqt_trim_poly_a_combined_at" {
  let seq = "AAATTTACGTAA"
  let qual = "IIIIIIIIIIII"
  let (t_seq, _) = @src.sqt_trim_poly_a(seq, qual, min_tail=3)
  assert_eq(t_seq, "CG")
}

///|
test "sqt_trim_poly_a_min_tail_threshold" {
  let seq = "AAACGTAAA"
  let qual = "IIIIIIIII"
  let (t_seq, _) = @src.sqt_trim_poly_a(seq, qual, min_tail=3)
  assert_eq(t_seq, "CG")
}

///|
test "sqt_trim_by_length_keep" {
  let seq = "ACGTACGTACGTACGT"
  let qual = "IIIIIIIIIIIIIIII"
  let result = @src.sqt_trim_by_length(seq, qual, min_length=5)
  assert_true(result is Some)
}

///|
test "sqt_trim_by_length_discard" {
  let seq = "ACGT"
  let qual = "IIII"
  let result = @src.sqt_trim_by_length(seq, qual, min_length=10)
  assert_true(result is None)
}

///|
test "sqt_trim_by_length_exact_boundary" {
  let seq = "ACGTACGT"
  let qual = "IIIIIIII"
  let result = @src.sqt_trim_by_length(seq, qual, min_length=8)
  assert_true(result is Some)
}

///|
test "sqt_trim_by_length_zero_length" {
  let seq = ""
  let qual = ""
  let result = @src.sqt_trim_by_length(seq, qual, min_length=0)
  assert_true(result is Some)
}

///|
test "sqt_trim_by_gc_within_range" {
  let seq = "ACGTACGT"
  let qual = "IIIIIIII"
  let result = @src.sqt_trim_by_gc(seq, qual, gc_min=0.3, gc_max=0.7)
  assert_true(result is Some)
}

///|
test "sqt_trim_by_gc_out_of_range" {
  let seq = "GGGGGGGG"
  let qual = "IIIIIIII"
  let result = @src.sqt_trim_by_gc(seq, qual, gc_min=0.3, gc_max=0.7)
  assert_true(result is None)
}

///|
test "sqt_trim_by_gc_extreme_at" {
  let seq = "AAAAAAAA"
  let qual = "IIIIIIII"
  let result = @src.sqt_trim_by_gc(seq, qual, gc_min=0.3, gc_max=0.7)
  assert_true(result is None)
}

///|
test "sqt_trim_by_gc_mid_range" {
  let seq = "ACGTACGTACGT"
  let qual = "IIIIIIIIIIII"
  let result = @src.sqt_trim_by_gc(seq, qual, gc_min=0.4, gc_max=0.6)
  assert_true(result is Some)
}

///|
test "sqt_trim_by_gc_boundary_values" {
  let seq = "ACGT"
  let qual = "IIII"
  let result = @src.sqt_trim_by_gc(seq, qual, gc_min=0.5, gc_max=0.5)
  assert_true(result is Some)
}

///|
test "sqt_trim_reads_batch_with_quality_config" {
  let reads = [
    @src.FastqRead::new("read1", "ACGTACGT", "IIIIIIII"),
    @src.FastqRead::new("read2", "ACGTACGTAGATCGGAA", "IIIIIIIIIIIIIIIII"),
  ]
  let qual_config = @src.QualityConfig::new(20, 10, 4, 20)
  let results = @src.sqt_trim_reads(reads, qual_config, None)
  assert_eq(results.length(), 2)
  assert_eq(results[0].trimmed_seq, "ACGTACGT")
  assert_eq(results[0].trim_type, "none")
}

///|
test "sqt_trim_reads_batch_with_adapter_config" {
  let reads = [
    @src.FastqRead::new("read1", "ACGTACGTAGATCGGAA", "IIIIIIIIIIIIIIIII"),
  ]
  let qual_config = @src.QualityConfig::new(20, 4, 4, 20)
  let adapter_config = @src.AdapterConfig::new(["AGATCGGAA"], 4, 0)
  let results = @src.sqt_trim_reads(
    reads, qual_config, Some(adapter_config),
  )
  assert_eq(results.length(), 1)
  assert_eq(results[0].trimmed_seq, "ACGTACGT")
  assert_true(results[0].trim_type == "adapter:AGATCGGAA")
}

///|
test "sqt_trim_reads_quality_trim_only" {
  let reads = [
    @src.FastqRead::new("r1", "ACGTACGT", "IIII####"),
  ]
  let qual_config = @src.QualityConfig::new(20, 4, 4, 20)
  let results = @src.sqt_trim_reads(reads, qual_config, None)
  assert_eq(results.length(), 1)
  assert_eq(results[0].trimmed_seq.length(), 4)
}

///|
test "sqt_trim_reads_poly_a_trim" {
  let reads = [
    @src.FastqRead::new("r1", "ACGTACGGGAAAA", "IIIIIIIIIIIIII"),
  ]
  let qual_config = @src.QualityConfig::new(20, 4, 4, 20)
  let results = @src.sqt_trim_reads(reads, qual_config, None)
  assert_eq(results.length(), 1)
  assert_eq(results[0].trimmed_seq, "ACGTACGGG")
}

///|
test "sqt_trim_reads_adapter_at_end" {
  let reads = [
    @src.FastqRead::new("r1", "ACGTACGTACGTAGATCGGAA", "IIIIIIIIIIIIIIIIIIIII"),
  ]
  let qual_config = @src.QualityConfig::new(20, 4, 4, 20)
  let adapter_config = @src.AdapterConfig::new(["AGATCGGAA"], 4, 0)
  let results = @src.sqt_trim_reads(
    reads, qual_config, Some(adapter_config),
  )
  assert_eq(results.length(), 1)
  assert_eq(results[0].trimmed_seq, "ACGTACGTACGT")
}

///|
test "sqt_trim_reads_empty_batch" {
  let reads : Array[@src.FastqRead] = []
  let qual_config = @src.QualityConfig::new(20, 10, 4, 20)
  let results = @src.sqt_trim_reads(reads, qual_config, None)
  assert_eq(results.length(), 0)
}

///|
test "sqt_trim_reads_preserves_original_data" {
  let reads = [
    @src.FastqRead::new("read_orig", "ACGTACGT", "IIIIIIII"),
  ]
  let qual_config = @src.QualityConfig::new(20, 4, 4, 20)
  let results = @src.sqt_trim_reads(reads, qual_config, None)
  assert_eq(results.length(), 1)
  assert_eq(results[0].original_seq, "ACGTACGT")
  assert_eq(results[0].original_qual, "IIIIIIII")
  assert_eq(results[0].read_id, "read_orig")
}

///|
test "sqt_compute_stats_basic" {
  let reads = [
    @src.FastqRead::new("r1", "ACGTACGT", "IIIIIIII"),
    @src.FastqRead::new("r2", "TTTTTTTT", "IIIIIIII"),
  ]
  let results = [
    @src.TrimResult::new("r1", "ACGTACGT", "ACGTACGT", "IIIIIIII", "IIIIIIII", "none"),
    @src.TrimResult::new("r2", "TTTTTTTT", "", "IIIIIIII", "", "poly_A"),
  ]
  let stats = @src.sqt_compute_stats(results, reads)
  assert_eq(stats.total_reads, 2)
  assert_eq(stats.kept_reads, 1)
  assert_eq(stats.discarded_reads, 1)
  assert_eq(stats.total_bp_before, 16)
  assert_eq(stats.total_bp_after, 8)
}

///|
test "sqt_compute_stats_gc_content" {
  let reads = [
    @src.FastqRead::new("gc_high", "GGGGCCCC", "IIIIIIII"),
    @src.FastqRead::new("at_high", "AAAATTTT", "IIIIIIII"),
  ]
  let results = [
    @src.TrimResult::new("gc_high", "GGGGCCCC", "GGGGCCCC", "IIIIIIII", "IIIIIIII", "none"),
    @src.TrimResult::new("at_high", "AAAATTTT", "AAAATTTT", "IIIIIIII", "IIIIIIII", "none"),
  ]
  let stats = @src.sqt_compute_stats(results, reads)
  assert_true(stats.gc_content_before > 0.4)
  assert_true(stats.gc_content_before < 0.6)
}

///|
test "sqt_compute_stats_empty" {
  let reads : Array[@src.FastqRead] = []
  let results : Array[@src.TrimResult] = []
  let stats = @src.sqt_compute_stats(results, reads)
  assert_eq(stats.total_reads, 0)
  assert_eq(stats.kept_reads, 0)
  assert_eq(stats.total_bp_before, 0)
  assert_eq(stats.gc_content_before, 0.0)
}

///|
test "sqt_compute_stats_all_discarded" {
  let reads = [
    @src.FastqRead::new("r1", "ACGT", "IIII"),
  ]
  let results = [
    @src.TrimResult::new("r1", "ACGT", "", "IIII", "", "quality"),
  ]
  let stats = @src.sqt_compute_stats(results, reads)
  assert_eq(stats.total_reads, 1)
  assert_eq(stats.kept_reads, 0)
  assert_eq(stats.discarded_reads, 1)
  assert_eq(stats.total_bp_after, 0)
  assert_eq(stats.gc_content_after, 0.0)
}

///|
test "sqt_fastq_parse_single_record" {
  let content = "@read1 desc\nACGTACGT\n+\nIIIIIIII\n"
  let reads = @src.sqt_fastq_parse(content)
  assert_eq(reads.length(), 1)
  assert_eq(reads[0].id, "read1")
  assert_eq(reads[0].sequence, "ACGTACGT")
  assert_eq(reads[0].quality, "IIIIIIII")
}

///|
test "sqt_fastq_parse_multiple_records" {
  let content = "@r1\nACGT\n+\nIIII\n@r2\nTGCA\n+\nIIII\n"
  let reads = @src.sqt_fastq_parse(content)
  assert_eq(reads.length(), 2)
  assert_eq(reads[0].id, "r1")
  assert_eq(reads[0].sequence, "ACGT")
  assert_eq(reads[1].id, "r2")
  assert_eq(reads[1].sequence, "TGCA")
}

///|
test "sqt_fastq_parse_with_description" {
  let content = "@read1 some description here\nACGT\n+\nIIII\n"
  let reads = @src.sqt_fastq_parse(content)
  assert_eq(reads.length(), 1)
  assert_eq(reads[0].id, "read1")
}

///|
test "sqt_fastq_parse_empty_content" {
  let content = ""
  let reads = @src.sqt_fastq_parse(content)
  assert_eq(reads.length(), 0)
}

///|
test "sqt_fastq_parse_malformed" {
  let content = "not a fastq file\n"
  let reads = @src.sqt_fastq_parse(content)
  assert_eq(reads.length(), 0)
}

///|
test "sqt_fastq_parse_crlf" {
  let content = "@r1 desc\r\nACGT\r\n+\r\nIIII\r\n"
  let reads = @src.sqt_fastq_parse(content)
  assert_eq(reads.length(), 1)
  assert_eq(reads[0].id, "r1")
  assert_eq(reads[0].sequence, "ACGT")
}

///|
test "sqt_fastq_parse_multiple_with_variable_length" {
  let content = "@r1\nACGTACGT\n+\nIIIIIIII\n@r2\nACGT\n+\nIIII\n"
  let reads = @src.sqt_fastq_parse(content)
  assert_eq(reads.length(), 2)
  assert_eq(reads[0].sequence.length(), 8)
  assert_eq(reads[1].sequence.length(), 4)
}

///|
test "sqt_fastq_serialize_single" {
  let reads = [
    @src.FastqRead::new("test1", "ACGT", "IIII"),
  ]
  let output = @src.sqt_fastq_serialize(reads)
  assert_eq(output, "@test1\nACGT\n+\nIIII\n")
}

///|
test "sqt_fastq_serialize_multiple" {
  let reads = [
    @src.FastqRead::new("r1", "AAAA", "IIII"),
    @src.FastqRead::new("r2", "CCCC", "IIII"),
  ]
  let output = @src.sqt_fastq_serialize(reads)
  assert_eq(output, "@r1\nAAAA\n+\nIIII\n@r2\nCCCC\n+\nIIII\n")
}

///|
test "sqt_fastq_roundtrip" {
  let original = "@sample1 desc\nACGTACGT\n+\nIIIIIIII\n@sample2\nTGCA\n+\nIIII\n"
  let reads = @src.sqt_fastq_parse(original)
  assert_eq(reads.length(), 2)
  let serialized = @src.sqt_fastq_serialize(reads)
  let reparsed = @src.sqt_fastq_parse(serialized)
  assert_eq(reparsed.length(), 2)
  assert_eq(reparsed[0].id, "sample1")
  assert_eq(reparsed[0].sequence, "ACGTACGT")
  assert_eq(reparsed[1].id, "sample2")
  assert_eq(reparsed[1].sequence, "TGCA")
}