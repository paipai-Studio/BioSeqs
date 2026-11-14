// Regression tests for GitHub Issue #107: parse_fastq silently discarded
// incomplete trailing records. All FASTQ batch parsers (Sanger, Illumina,
// Solexa) must now raise SeqIOError when non-blank lines remain after the
// last complete 4-line record, while tolerating trailing blank lines.

///|
test "issue107_parse_fastq_truncated_raises" {
  let raised = try {
    ignore(@src.parse_fastq("@s1\nACGT\n+\nIIII\nEXTRA_LINE"))
    false
  } catch {
    SeqIOError(_) => true
  }
  assert_true(raised)
}

///|
test "issue107_error_message_reports_remaining_line_count" {
  let msg = try {
    ignore(@src.parse_fastq("@s1\nACGT\n+\nIIII\nEXTRA_LINE"))
    ""
  } catch {
    SeqIOError(m) => m
  }
  // Exactly 1 leftover line should be reported (not the total line count).
  assert_true(msg.contains("got 1 remaining lines"))
}

///|
test "issue107_parse_fastq_ok_with_trailing_blank_lines" {
  // Trailing blank lines are the normal result of a final newline; no error.
  let records = @src.parse_fastq("@s1\nACGT\n+\nIIII\n\n\n")
  assert_eq(records.length(), 1)
}

///|
test "issue107_parse_fastq_illumina_truncated_raises" {
  let raised = try {
    ignore(@src.parse_fastq_illumina("@s1\nACGT\n+\n@@@@\nEXTRA"))
    false
  } catch {
    SeqIOError(_) => true
  }
  assert_true(raised)
}

///|
test "issue107_parse_fastq_solexa_truncated_raises" {
  let raised = try {
    ignore(@src.parse_fastq_solexa("@s1\nACGT\n+\n@@@@\nEXTRA"))
    false
  } catch {
    SeqIOError(_) => true
  }
  assert_true(raised)
}

///|
test "issue107_complete_records_still_parse" {
  let records = @src.parse_fastq(
    "@s1 d\nACGT\n+s1 d\nIIII\n@s2\nTTTT\n+\nJJJJ\n",
  )
  assert_eq(records.length(), 2)
  assert_eq(records[0].seq.to_string(), "ACGT")
  assert_eq(records[1].seq.to_string(), "TTTT")
}

// ============================================================
// Probe group: truncation shapes — every prefix length below a full
// record must raise, so silent data loss is impossible.
// ============================================================

///|
test "issue107_probe_single_line_raises" {
  // Header alone: 1 of 4 lines.
  let raised = try {
    ignore(@src.parse_fastq("@s1"))
    false
  } catch {
    SeqIOError(_) => true
  }
  assert_true(raised)
}

///|
test "issue107_probe_two_lines_raises" {
  // Header + sequence: 2 of 4 lines (missing separator and quality).
  let raised = try {
    ignore(@src.parse_fastq("@s1\nACGT"))
    false
  } catch {
    SeqIOError(_) => true
  }
  assert_true(raised)
}

///|
test "issue107_probe_three_lines_raises" {
  // Header + sequence + separator: 3 of 4 lines (missing quality).
  let raised = try {
    ignore(@src.parse_fastq("@s1\nACGT\n+"))
    false
  } catch {
    SeqIOError(_) => true
  }
  assert_true(raised)
}

///|
test "issue107_probe_exactly_four_lines_ok" {
  // Exactly one complete record without a trailing newline must parse.
  let records = @src.parse_fastq("@s1\nACGT\n+\nIIII")
  assert_eq(records.length(), 1)
  assert_eq(records[0].seq.to_string(), "ACGT")
}

// ============================================================
// Probe group: line-ending and whitespace tolerance at the tail.
// ============================================================

///|
test "issue107_probe_crlf_tail_not_flagged" {
  // A Windows-style trailing newline splits into "\r" + ""; neither may
  // be mistaken for an incomplete record.
  let records = @src.parse_fastq("@s1\nACGT\n+\nIIII\r\n")
  assert_eq(records.length(), 1)
}

///|
test "issue107_probe_whitespace_only_tail_raises" {
  // Documented strict behavior: a whitespace-only leftover line is NOT
  // treated as blank — it is a truncated record remnant and raises.
  // (No legal FASTQ record component can consist solely of spaces.)
  let raised = try {
    ignore(@src.parse_fastq("@s1\nACGT\n+\nIIII\n   "))
    false
  } catch {
    SeqIOError(_) => true
  }
  assert_true(raised)
}

///|
test "issue107_probe_multi_record_then_truncation" {
  // Two complete records followed by a 2-line orphan must raise —
  // data before the truncation must not mask the error.
  let raised = try {
    ignore(@src.parse_fastq("@a\nAAAA\n+\n@@@@\n@b\nTTTT\n+\n@@@@\n@c\nGG"))
    false
  } catch {
    SeqIOError(_) => true
  }
  assert_true(raised)
}

///|
test "issue107_probe_illumina_trailing_blank_ok" {
  // The parse_fastq_offset path (Illumina) must tolerate trailing blank
  // lines just like the Sanger parser.
  let records = @src.parse_fastq_illumina("@s1\nACGT\n+\n@@@@\n\n")
  assert_eq(records.length(), 1)
}
