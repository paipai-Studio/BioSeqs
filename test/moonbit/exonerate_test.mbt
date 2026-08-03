///|
/// Tests for Bio.SearchIO.ExonerateIO module (vulgar and cigar parsing).

///|
/// AlignmentBlock construction and accessors.
test "exonerate_alignment_block_new" {
  let block = @src.AlignmentBlock::new("M", 100, 100, 50.0)
  assert_eq(block.code(), "M")
  assert_eq(block.query_length(), 100)
  assert_eq(block.target_length(), 100)
  assert_eq(block.block_score(), 50.0)
}

///|
/// AlignmentBlock with intron code.
test "exonerate_alignment_block_intron" {
  let block = @src.AlignmentBlock::new("I", 0, 0, 2000.0)
  assert_eq(block.code(), "I")
  assert_eq(block.query_length(), 0)
  assert_eq(block.target_length(), 0)
  assert_eq(block.block_score(), 2000.0)
}

///|
/// AlignmentBlock with splice site code.
test "exonerate_alignment_block_splice" {
  let block = @src.AlignmentBlock::new("5", 0, 0, 0.0)
  assert_eq(block.code(), "5")
  assert_eq(block.query_length(), 0)
  assert_eq(block.target_length(), 0)
}

///|
/// ExonerateRecord construction and accessors.
test "exonerate_record_new" {
  let blocks : Array[@src.AlignmentBlock] = [
    @src.AlignmentBlock::new("M", 50, 50, 50.0),
    @src.AlignmentBlock::new("I", 0, 0, 1000.0),
  ]
  let rec = @src.ExonerateRecord::new(
    "query1", 0, 100, "+",
    "target1", 0, 95, "+",
    1850.0, blocks,
  )
  assert_eq(rec.query_name(), "query1")
  assert_eq(rec.query_start(), 0)
  assert_eq(rec.query_end(), 100)
  assert_eq(rec.query_strand(), "+")
  assert_eq(rec.target_name(), "target1")
  assert_eq(rec.target_start(), 0)
  assert_eq(rec.target_end(), 95)
  assert_eq(rec.target_strand(), "+")
  assert_eq(rec.score(), 1850.0)
  assert_eq(rec.alignment_blocks().length(), 2)
}

///|
/// Parse a single vulgar line.
test "exonerate_parse_vulgar_single" {
  let data = "vulgar: query1 0 100 + target1 0 95 + 1850 M 80 80 80 I 0 0 2000 M 20 15 20\n"
  let records = @src.parse_exonerate_vulgar(data)
  assert_eq(records.length(), 1)
  let rec = records[0]
  assert_eq(rec.query_name(), "query1")
  assert_eq(rec.query_start(), 0)
  assert_eq(rec.query_end(), 100)
  assert_eq(rec.query_strand(), "+")
  assert_eq(rec.target_name(), "target1")
  assert_eq(rec.target_start(), 0)
  assert_eq(rec.target_end(), 95)
  assert_eq(rec.target_strand(), "+")
  assert_eq(rec.score(), 1850.0)
  assert_eq(rec.alignment_blocks().length(), 3)
  assert_eq(rec.alignment_blocks()[0].code(), "M")
  assert_eq(rec.alignment_blocks()[0].query_length(), 80)
  assert_eq(rec.alignment_blocks()[1].code(), "I")
  assert_eq(rec.alignment_blocks()[1].block_score(), 2000.0)
  assert_eq(rec.alignment_blocks()[2].code(), "M")
  assert_eq(rec.alignment_blocks()[2].query_length(), 20)
}

///|
/// Parse multiple vulgar lines.
test "exonerate_parse_vulgar_multiple" {
  let records = @src.parse_exonerate_vulgar(@src.exonerate_sample_data())
  assert_eq(records.length(), 3)
  assert_eq(records[0].query_name(), "query1")
  assert_eq(records[1].query_name(), "query2")
  assert_eq(records[2].query_name(), "query3")
}

///|
/// Parse vulgar with splice site blocks (5 and 3 codes with score field).
test "exonerate_parse_vulgar_splice" {
  let data = "vulgar: q1 0 200 + t1 0 200 + 3200 M 100 100 100 5 0 0 0 I 0 0 1500 3 0 0 0 M 100 100 100\n"
  let records = @src.parse_exonerate_vulgar(data)
  assert_eq(records.length(), 1)
  let rec = records[0]
  assert_eq(rec.alignment_blocks().length(), 5)
  assert_eq(rec.alignment_blocks()[0].code(), "M")
  assert_eq(rec.alignment_blocks()[1].code(), "5")
  assert_eq(rec.alignment_blocks()[2].code(), "I")
  assert_eq(rec.alignment_blocks()[2].block_score(), 1500.0)
  assert_eq(rec.alignment_blocks()[3].code(), "3")
  assert_eq(rec.alignment_blocks()[4].code(), "M")
}

///|
/// Parse a single cigar line.
test "exonerate_parse_cigar_single" {
  let data = "cigar: query1 0 100 + target1 0 95 + 1850 M 80 I 2000 M 20\n"
  let records = @src.parse_exonerate_cigar(data)
  assert_eq(records.length(), 1)
  let rec = records[0]
  assert_eq(rec.query_name(), "query1")
  assert_eq(rec.target_name(), "target1")
  assert_eq(rec.score(), 1850.0)
  assert_eq(rec.alignment_blocks().length(), 3)
  assert_eq(rec.alignment_blocks()[0].code(), "M")
  assert_eq(rec.alignment_blocks()[0].query_length(), 80)
  assert_eq(rec.alignment_blocks()[0].target_length(), 80)
  assert_eq(rec.alignment_blocks()[1].code(), "I")
  assert_eq(rec.alignment_blocks()[1].block_score(), 2000.0)
  assert_eq(rec.alignment_blocks()[1].query_length(), 0)
  assert_eq(rec.alignment_blocks()[2].code(), "M")
  assert_eq(rec.alignment_blocks()[2].query_length(), 20)
}

///|
/// Parse multiple cigar lines.
test "exonerate_parse_cigar_multiple" {
  let records = @src.parse_exonerate_cigar(@src.exonerate_sample_cigar_data())
  assert_eq(records.length(), 3)
  assert_eq(records[0].query_name(), "query1")
  assert_eq(records[2].query_name(), "query3")
}

///|
/// Auto-detect vulgar format.
test "exonerate_parse_auto_vulgar" {
  let data = "vulgar: query1 0 100 + target1 0 95 + 1850 M 80 80 80 I 0 0 2000 M 20 15 20\n"
  let records = @src.parse_exonerate(data)
  assert_eq(records.length(), 1)
  assert_eq(records[0].query_name(), "query1")
  assert_eq(records[0].alignment_blocks().length(), 3)
}

///|
/// Auto-detect cigar format.
test "exonerate_parse_auto_cigar" {
  let data = "cigar: query1 0 100 + target1 0 95 + 1850 M 80 I 2000 M 20\n"
  let records = @src.parse_exonerate(data)
  assert_eq(records.length(), 1)
  assert_eq(records[0].query_name(), "query1")
  assert_eq(records[0].alignment_blocks().length(), 3)
}

///|
/// Auto-detect mixed vulgar and cigar lines.
test "exonerate_parse_auto_mixed" {
  let data = "vulgar: q1 0 100 + t1 0 95 + 1850 M 80 80 80 I 0 0 2000 M 20 15 20\n" +
    "cigar: q2 0 50 + t2 0 50 + 850 M 50\n"
  let records = @src.parse_exonerate(data)
  assert_eq(records.length(), 2)
  assert_eq(records[0].query_name(), "q1")
  assert_eq(records[1].query_name(), "q2")
}

///|
/// Empty input returns no records.
test "exonerate_parse_empty" {
  let records = @src.parse_exonerate("")
  assert_eq(records.length(), 0)
  let records_v = @src.parse_exonerate_vulgar("")
  assert_eq(records_v.length(), 0)
  let records_c = @src.parse_exonerate_cigar("")
  assert_eq(records_c.length(), 0)
}

///|
/// Non-exonerate lines are ignored.
test "exonerate_parse_ignores_other_lines" {
  let data = "header line\n" +
    "vulgar: q1 0 50 + t1 0 50 + 850 M 50 50 50\n" +
    "some comment\n" +
    "\n"
  let records = @src.parse_exonerate(data)
  assert_eq(records.length(), 1)
  assert_eq(records[0].query_name(), "q1")
}

///|
/// Filter records by score.
test "exonerate_filter_by_score" {
  let records = @src.parse_exonerate_vulgar(@src.exonerate_sample_data())
  assert_eq(records.length(), 3)
  let filtered = @src.exonerate_filter_by_score(records, 1000.0)
  assert_eq(filtered.length(), 2)
  assert_eq(filtered[0].query_name(), "query1")
  assert_eq(filtered[1].query_name(), "query2")
}

///|
/// Filter records by a high threshold returns none.
test "exonerate_filter_by_score_high" {
  let records = @src.parse_exonerate_vulgar(@src.exonerate_sample_data())
  let filtered = @src.exonerate_filter_by_score(records, 5000.0)
  assert_eq(filtered.length(), 0)
}

///|
/// Total alignment length sums M blocks' query lengths.
test "exonerate_total_alignment_length" {
  let data = "vulgar: q1 0 100 + t1 0 95 + 1850 M 80 80 80 I 0 0 2000 M 20 15 20\n"
  let records = @src.parse_exonerate_vulgar(data)
  let total = @src.exonerate_total_alignment_length(records[0])
  assert_eq(total, 100)
}

///|
/// Total intron length sums I blocks' scores.
test "exonerate_total_intron_length" {
  let data = "vulgar: q1 0 200 + t1 0 200 + 3200 M 100 100 100 5 0 0 0 I 0 0 1500 3 0 0 0 I 0 0 500 M 100 100 100\n"
  let records = @src.parse_exonerate_vulgar(data)
  let total = @src.exonerate_total_intron_length(records[0])
  assert_eq(total, 2000)
}

///|
/// Count intron blocks.
test "exonerate_count_introns" {
  let data = "vulgar: q1 0 200 + t1 0 200 + 3200 M 100 100 100 5 0 0 0 I 0 0 1500 3 0 0 0 I 0 0 500 M 100 100 100\n"
  let records = @src.parse_exonerate_vulgar(data)
  let count = @src.exonerate_count_introns(records[0])
  assert_eq(count, 2)
}

///|
/// Count match blocks.
test "exonerate_count_matches" {
  let data = "vulgar: q1 0 100 + t1 0 95 + 1850 M 80 80 80 I 0 0 2000 M 20 15 20\n"
  let records = @src.parse_exonerate_vulgar(data)
  let count = @src.exonerate_count_matches(records[0])
  assert_eq(count, 2)
}

///|
/// Record with no introns.
test "exonerate_no_introns" {
  let data = "vulgar: q3 0 50 - t3 500 550 + 850 M 50 50 50\n"
  let records = @src.parse_exonerate_vulgar(data)
  assert_eq(records.length(), 1)
  assert_eq(@src.exonerate_count_introns(records[0]), 0)
  assert_eq(@src.exonerate_total_intron_length(records[0]), 0)
  assert_eq(@src.exonerate_count_matches(records[0]), 1)
  assert_eq(@src.exonerate_total_alignment_length(records[0]), 50)
}

///|
/// Reconstruct vulgar string from a record.
test "exonerate_to_vulgar_string" {
  let data = "vulgar: q1 0 100 + t1 0 95 + 1850 M 80 80 80 I 0 0 2000 M 20 15 20\n"
  let records = @src.parse_exonerate_vulgar(data)
  let s = @src.exonerate_to_vulgar_string(records[0])
  assert_eq(s, "M 80 80 80 I 0 0 2000 M 20 15 20")
}

///|
/// Reconstruct cigar string from a record.
test "exonerate_to_cigar_string" {
  let data = "cigar: q1 0 100 + t1 0 95 + 1850 M 80 I 2000 M 20\n"
  let records = @src.parse_exonerate_cigar(data)
  let s = @src.exonerate_to_cigar_string(records[0])
  assert_eq(s, "M 80 I 2000 M 20")
}

///|
/// Reconstruct cigar string from a vulgar-parsed record (intron length used).
test "exonerate_to_cigar_string_from_vulgar" {
  let data = "vulgar: q1 0 100 + t1 0 95 + 1850 M 80 80 80 I 0 0 2000 M 20 15 20\n"
  let records = @src.parse_exonerate_vulgar(data)
  let s = @src.exonerate_to_cigar_string(records[0])
  assert_eq(s, "M 80 I 2000 M 20")
}

///|
/// Summary string contains expected content.
test "exonerate_summary" {
  let records = @src.parse_exonerate_vulgar(@src.exonerate_sample_data())
  let summary = @src.exonerate_summary(records)
  assert_true(summary.contains("Exonerate records: 3"))
  assert_true(summary.contains("query1"))
  assert_true(summary.contains("query2"))
  assert_true(summary.contains("query3"))
  assert_true(summary.contains("score=1850"))
}

///|
/// Sample data is non-empty and parseable.
test "exonerate_sample_data" {
  let data = @src.exonerate_sample_data()
  assert_true(data.length() > 0)
  let records = @src.parse_exonerate_vulgar(data)
  assert_eq(records.length(), 3)
}

///|
/// Sample cigar data is non-empty and parseable.
test "exonerate_sample_cigar_data" {
  let data = @src.exonerate_sample_cigar_data()
  assert_true(data.length() > 0)
  let records = @src.parse_exonerate_cigar(data)
  assert_eq(records.length(), 3)
}

///|
/// Negative strand is parsed correctly.
test "exonerate_negative_strand" {
  let data = "vulgar: q3 0 50 - t3 500 550 + 850 M 50 50 50\n"
  let records = @src.parse_exonerate_vulgar(data)
  assert_eq(records[0].query_strand(), "-")
  assert_eq(records[0].target_strand(), "+")
}

///|
/// Cigar with only a single match block (edge case).
test "exonerate_cigar_single_block" {
  let data = "cigar: q1 0 50 + t1 0 50 + 850 M 50\n"
  let records = @src.parse_exonerate_cigar(data)
  assert_eq(records.length(), 1)
  assert_eq(records[0].alignment_blocks().length(), 1)
  assert_eq(@src.exonerate_total_alignment_length(records[0]), 50)
  assert_eq(@src.exonerate_count_matches(records[0]), 1)
}
