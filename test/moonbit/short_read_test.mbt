///|
test "short_read_qa_summary_basic" {
  let sequences = ["ACGTACGT", "ACGTACGT", "ACGTACGT"]
  let qualities = ["IIIIIIII", "IIIIIIII", "IIIIIIII"]
  let stats = @src.qa_summary(sequences, qualities)
  assert_eq(stats.n_reads, 3)
  assert_eq(stats.read_length, 8)
  assert_true(stats.gc_content > 0.4 && stats.gc_content < 0.6)
  assert_true(stats.mean_quality > 30.0)
}

///|
test "short_read_qa_summary_empty" {
  let sequences : Array[String] = []
  let qualities : Array[String] = []
  let stats = @src.qa_summary(sequences, qualities)
  assert_eq(stats.n_reads, 0)
  assert_eq(stats.read_length, 0)
}

///|
test "short_read_qa_summary_gc_content" {
  let sequences = ["GGGGCCCC", "GGGGCCCC"]
  let qualities = ["IIIIIIII", "IIIIIIII"]
  let stats = @src.qa_summary(sequences, qualities)
  assert_true(stats.gc_content > 0.99)
}

///|
test "short_read_qa_summary_dup_rate" {
  let sequences = ["ACGTACGT", "ACGTACGT", "ACGTACGT"]
  let qualities = ["IIIIIIII", "IIIIIIII", "IIIIIIII"]
  let stats = @src.qa_summary(sequences, qualities)
  assert_true(stats.dup_rate > 0.5)
}

///|
test "short_read_trim_tail_basic" {
  let seq = "ACGTACGT"
  let qual = "IIIIII##"
  let (trimmed_seq, trimmed_qual) = @src.trim_tail(seq, qual, 20)
  assert_eq(trimmed_seq.length(), 6)
  assert_eq(trimmed_qual.length(), 6)
}

///|
test "short_read_trim_tail_all_low_quality" {
  let seq = "ACGTACGT"
  let qual = "########"
  let (trimmed_seq, trimmed_qual) = @src.trim_tail(seq, qual, 20)
  assert_eq(trimmed_seq.length(), 0)
  assert_eq(trimmed_qual.length(), 0)
}

///|
test "short_read_trim_leading_basic" {
  let seq = "ACGTACGT"
  let qual = "##IIIIII"
  let (trimmed_seq, trimmed_qual) = @src.trim_leading(seq, qual, 20)
  assert_eq(trimmed_seq.length(), 6)
  assert_eq(trimmed_qual.length(), 6)
}

///|
test "short_read_trim_adapter_basic" {
  let seq = "ACGTACGTAGATCGGAA"
  let qual = "IIIIIIIIIIIIIIII"
  let adapter = "AGATCGGAA"
  let (trimmed_seq, trimmed_qual) = @src.trim_adapter(seq, qual, adapter, 0)
  assert_eq(trimmed_seq.length(), 8)
}

///|
test "short_read_trim_adapter_no_match" {
  let seq = "ACGTACGTACGTACGT"
  let qual = "IIIIIIIIIIIIIIII"
  let adapter = "GGGGGGGG"
  let (trimmed_seq, trimmed_qual) = @src.trim_adapter(seq, qual, adapter, 0)
  assert_eq(trimmed_seq.length(), 16)
}

///|
test "short_read_filter_reads_basic" {
  let sequences = ["ACGTACGT", "ACGTNNNN", "AAAAAAAA"]
  let qualities = ["IIIIIIII", "IIIIIIII", "########"]
  let (filtered_seqs, filtered_quals) = @src.filter_reads(sequences, qualities, 30.0, 5, 2)
  assert_eq(filtered_seqs.length(), 1)
  assert_eq(filtered_quals.length(), 1)
}

///|
test "short_read_filter_reads_length" {
  let sequences = ["ACGT", "ACGTACGTACGT"]
  let qualities = ["IIII", "IIIIIIIIIIII"]
  let (filtered_seqs, filtered_quals) = @src.filter_reads(sequences, qualities, 20.0, 8, 10)
  assert_eq(filtered_seqs.length(), 1)
}

///|
test "short_read_fastqc_report" {
  let sequences = ["ACGTACGT", "ACGTACGT", "ACGTACGT"]
  let qualities = ["IIIIIIII", "IIIIIIII", "IIIIIIII"]
  let report = @src.generate_fastqc_report(sequences, qualities)
  assert_eq(report.status, "PASS")
  assert_eq(report.per_base_quality.length(), 8)
  assert_eq(report.per_base_n_content.length(), 8)
}

///|
test "short_read_fastqc_report_fail" {
  let sequences = ["ACGTACGT"]
  let qualities = ["########"]
  let report = @src.generate_fastqc_report(sequences, qualities)
  assert_eq(report.status, "FAIL")
}

///|
test "short_read_quality_by_cycle" {
  let sequences = ["ACGTACGT", "ACGTACGT"]
  let qualities = ["IIIIIIII", "IIIIIIII"]
  let stats = @src.qa_summary(sequences, qualities)
  assert_eq(stats.quality_by_cycle.length(), 8)
  assert_true(stats.quality_by_cycle[0] > 30.0)
}

///|
test "short_read_n_count_by_cycle" {
  let sequences = ["ACGTNCGT", "ACGTNCGT"]
  let qualities = ["IIIIIIII", "IIIIIIII"]
  let stats = @src.qa_summary(sequences, qualities)
  assert_eq(stats.n_count_by_cycle[4], 2)
}
