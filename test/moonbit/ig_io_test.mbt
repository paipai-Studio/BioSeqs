///|
/// Tests for the IgIO module (Bio.SeqIO.IgIO equivalent).
///
/// These tests exercise the public API exposed in `src/ig_io.mbt`:
/// the IgRecord struct and its accessors, the parser (`ig_parse`),
/// the serializers (`ig_record_to_string`, `ig_record_to_string_width`,
/// `ig_records_to_string`, `ig_to_string`), the SeqRecord interop
/// (`ig_to_seq_record`, `ig_to_seq_records`, `ig_from_seq_record`,
/// `ig_from_seq_records`), and the sample data helpers.

// ============================================================================
// 1. IgRecord - creation and accessors
// ============================================================================

///|
test "ig_record_new_and_accessors" {
  let r = @src.IgRecord::new(
    "A_U455",
    "HIV-1 group M subtype A",
    "ATGGCAGCTGTGATGAAGCAGAGACGGGTAAGAGCTC",
  )
  assert_eq(r.title(), "A_U455")
  assert_eq(r.comment(), "HIV-1 group M subtype A")
  assert_eq(r.sequence(), "ATGGCAGCTGTGATGAAGCAGAGACGGGTAAGAGCTC")
}

///|
test "ig_record_id_name_description_match_title" {
  let r = @src.IgRecord::new("B_HXB2R", "comment", "ATGC")
  assert_eq(r.id(), "B_HXB2R")
  assert_eq(r.name(), "B_HXB2R")
  assert_eq(r.description(), "B_HXB2R")
}

///|
test "ig_record_length" {
  let r = @src.IgRecord::new("X", "", "ATGCATGC")
  assert_eq(r.length(), 8)
}

///|
test "ig_record_length_empty_sequence" {
  let r = @src.IgRecord::new("X", "", "")
  assert_eq(r.length(), 0)
}

///|
test "ig_record_to_string_contains_id_and_length" {
  let r = @src.IgRecord::new("A_U455", "comment", "ATGCATGC")
  let s = r.to_string()
  assert_true(s.contains("IgRecord(id=A_U455"))
  assert_true(s.contains("length=8"))
}

///|
test "ig_record_to_string_shows_comment_chars_when_present" {
  let r = @src.IgRecord::new("X", "five", "ATGC")
  let s = r.to_string()
  assert_true(s.contains("comment=4 chars"))
}

///|
test "ig_record_to_string_no_comment_when_empty" {
  let r = @src.IgRecord::new("X", "", "ATGC")
  let s = r.to_string()
  assert_false(s.contains("comment="))
}

// ============================================================================
// 2. Parsing - basic
// ============================================================================

///|
test "ig_parse_single_record" {
  let text = ";\nA_U455\nATGGCAGCTGTGATGAAGCAGAGACGGGTAAGAGCTC\n1\n"
  let records = @src.ig_parse(text)
  assert_eq(records.length(), 1)
  assert_eq(records[0].title(), "A_U455")
  assert_eq(records[0].sequence(), "ATGGCAGCTGTGATGAAGCAGAGACGGGTAAGAGCTC")
}

///|
test "ig_parse_strips_trailing_terminator" {
  let text = ";\nA_U455\nATGCATGCATGC\n1\n"
  let records = @src.ig_parse(text)
  assert_eq(records[0].sequence(), "ATGCATGCATGC")
}

///|
test "ig_parse_multiple_records" {
  let text = ";\nA_U455\nATGGCAGCTGTG\n1\n;\nB_HXB2R\nATGGGTGCTGTG\n1\n"
  let records = @src.ig_parse(text)
  assert_eq(records.length(), 2)
  assert_eq(records[0].title(), "A_U455")
  assert_eq(records[1].title(), "B_HXB2R")
}

///|
test "ig_parse_collects_comment_lines" {
  let text = ";HIV-1 group M\n;subtype A\nA_U455\nATGCATGC\n1\n"
  let records = @src.ig_parse(text)
  assert_eq(records.length(), 1)
  let comment = records[0].comment()
  assert_true(comment.contains("HIV-1 group M"))
  assert_true(comment.contains("subtype A"))
}

///|
test "ig_parse_comment_with_leading_space" {
  let text = "; HIV-1 group M\nA_U455\nATGCATGC\n1\n"
  let records = @src.ig_parse(text)
  assert_eq(records[0].comment(), "HIV-1 group M")
}

///|
test "ig_parse_empty_comment_line" {
  let text = ";\nA_U455\nATGCATGC\n1\n"
  let records = @src.ig_parse(text)
  assert_eq(records.length(), 1)
  assert_eq(records[0].comment(), "")
}

///|
test "ig_parse_removes_internal_whitespace_from_sequence" {
  let text = ";\nA_U455\nATGC ATGC ATGC\n1\n"
  let records = @src.ig_parse(text)
  assert_eq(records[0].sequence(), "ATGCATGCATGC")
}

///|
test "ig_parse_multiline_sequence" {
  let text = ";\nA_U455\nATGCATGC\nATGCATGC\n1\n"
  let records = @src.ig_parse(text)
  assert_eq(records[0].sequence(), "ATGCATGCATGCATGC")
}

// ============================================================================
// 3. Parsing - file headers
// ============================================================================

///|
test "ig_parse_skips_file_header_lines" {
  let text = ";;\n;; This is a file header\n;\nA_U455\nATGCATGC\n1\n"
  let records = @src.ig_parse(text)
  assert_eq(records.length(), 1)
  assert_eq(records[0].title(), "A_U455")
}

///|
test "ig_parse_file_header_only_no_records" {
  let text = ";;\n;; header only\n"
  let records = @src.ig_parse(text)
  assert_eq(records.length(), 0)
}

///|
test "ig_parse_empty_string" {
  let records = @src.ig_parse("")
  assert_eq(records.length(), 0)
}

// ============================================================================
// 4. Parsing - edge cases
// ============================================================================

///|
test "ig_parse_no_terminator" {
  // Sequence without trailing "1".
  let text = ";\nA_U455\nATGCATGC\n"
  let records = @src.ig_parse(text)
  assert_eq(records.length(), 1)
  assert_eq(records[0].sequence(), "ATGCATGC")
}

///|
test "ig_parse_crlf_line_endings" {
  let text = ";\r\nA_U455\r\nATGCATGC\r\n1\r\n"
  let records = @src.ig_parse(text)
  assert_eq(records.length(), 1)
  assert_eq(records[0].title(), "A_U455")
  assert_eq(records[0].sequence(), "ATGCATGC")
}

///|
test "ig_parse_handles_blank_lines_between_records" {
  let text = ";\nA_U455\nATGCATGC\n1\n\n;\nB_HXB2R\nATGCATGC\n1\n"
  let records = @src.ig_parse(text)
  assert_eq(records.length(), 2)
}

///|
test "ig_parse_sequence_with_tabs" {
  let text = ";\nA_U455\nATGC\tATGC\tATGC\n1\n"
  let records = @src.ig_parse(text)
  assert_eq(records[0].sequence(), "ATGCATGCATGC")
}

// ============================================================================
// 5. Serialization - ig_record_to_string
// ============================================================================

///|
test "ig_record_to_string_basic" {
  let r = @src.IgRecord::new("A_U455", "", "ATGCATGC")
  let s = @src.ig_record_to_string(r)
  assert_true(s.contains("A_U455"))
  assert_true(s.contains("ATGCATGC"))
  assert_true(s.contains("1"))
}

///|
test "ig_record_to_string_includes_comment_lines" {
  let r = @src.IgRecord::new("A_U455", "HIV-1 group M", "ATGCATGC")
  let s = @src.ig_record_to_string(r)
  // Comment lines start with ";"
  assert_true(s.contains("; HIV-1 group M"))
}

///|
test "ig_record_to_string_empty_comment_emits_semicolon_line" {
  let r = @src.IgRecord::new("A_U455", "", "ATGCATGC")
  let s = @src.ig_record_to_string(r)
  // When comment is empty, a minimal ";" line is emitted.
  assert_true(s.contains(";\n"))
}

///|
test "ig_record_to_string_wraps_long_sequence" {
  let long_seq = "ATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGC"
  let r = @src.IgRecord::new("X", "", long_seq)
  let s = @src.ig_record_to_string(r)
  // The default width is 60; the sequence should be split into multiple lines.
  let lines = s.split("\n").to_array()
  let mut max_len = 0
  for i in 0..<lines.length() {
    let line = lines[i].to_owned()
    // Skip non-sequence lines (comment, title, terminator).
    if line.length() > 0 &&
      line.unsafe_get(0) != ';'.to_int().to_uint16() &&
      line != "X" &&
      line != "1" {
      if line.length() > max_len {
        max_len = line.length()
      }
    }
  }
  assert_true(max_len <= 60)
}

///|
test "ig_record_to_string_width_zero_disables_wrapping" {
  let long_seq = "ATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGC"
  let r = @src.IgRecord::new("X", "", long_seq)
  let s = @src.ig_record_to_string_width(r, 0)
  assert_true(s.contains(long_seq))
}

///|
test "ig_record_to_string_multiline_comment" {
  let r = @src.IgRecord::new("X", "line1\nline2\nline3", "ATGC")
  let s = @src.ig_record_to_string(r)
  assert_true(s.contains("; line1"))
  assert_true(s.contains("; line2"))
  assert_true(s.contains("; line3"))
}

// ============================================================================
// 6. Serialization - ig_records_to_string
// ============================================================================

///|
test "ig_records_to_string_multiple_records" {
  let r1 = @src.IgRecord::new("A", "", "ATGC")
  let r2 = @src.IgRecord::new("B", "", "ATGC")
  let records = [r1, r2]
  let s = @src.ig_records_to_string(records, None)
  assert_true(s.contains("A\n"))
  assert_true(s.contains("B\n"))
}

///|
test "ig_records_to_string_with_file_header" {
  let r = @src.IgRecord::new("A", "", "ATGC")
  let header = [";;", ";; Sample file"]
  let s = @src.ig_records_to_string([r], Some(header))
  assert_true(s.contains(";;\n"))
  assert_true(s.contains(";; Sample file\n"))
}

///|
test "ig_records_to_string_header_lines_prefixed_with_semicolons" {
  let r = @src.IgRecord::new("A", "", "ATGC")
  let header = ["Custom header"]
  let s = @src.ig_records_to_string([r], Some(header))
  assert_true(s.contains(";;Custom header\n"))
}

///|
test "ig_records_to_string_empty_records" {
  let s = @src.ig_records_to_string([], None)
  assert_eq(s, "")
}

// ============================================================================
// 7. Round-trip: parse -> serialize -> parse
// ============================================================================

///|
test "ig_round_trip_single_record" {
  let original = @src.IgRecord::new("A_U455", "test comment", "ATGCATGCATGC")
  let text = @src.ig_record_to_string(original)
  let parsed = @src.ig_parse(text)
  assert_eq(parsed.length(), 1)
  assert_eq(parsed[0].title(), "A_U455")
  assert_eq(parsed[0].sequence(), "ATGCATGCATGC")
  assert_eq(parsed[0].comment(), "test comment")
}

///|
test "ig_round_trip_multiple_records" {
  let r1 = @src.IgRecord::new("A", "comment A", "ATGCATGC")
  let r2 = @src.IgRecord::new("B", "comment B", "ATGCATGCATGC")
  let records = [r1, r2]
  let text = @src.ig_records_to_string(records, None)
  let parsed = @src.ig_parse(text)
  assert_eq(parsed.length(), 2)
  assert_eq(parsed[0].title(), "A")
  assert_eq(parsed[1].title(), "B")
  assert_eq(parsed[0].sequence(), "ATGCATGC")
  assert_eq(parsed[1].sequence(), "ATGCATGCATGC")
}

///|
test "ig_round_trip_with_file_header" {
  let r = @src.IgRecord::new("X", "comment", "ATGCATGC")
  let header = [";;", ";; Test file"]
  let text = @src.ig_records_to_string([r], Some(header))
  let parsed = @src.ig_parse(text)
  assert_eq(parsed.length(), 1)
  assert_eq(parsed[0].title(), "X")
}

///|
test "ig_round_trip_long_sequence" {
  let long_seq = "ATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGC"
  let r = @src.IgRecord::new("LONG", "", long_seq)
  let text = @src.ig_record_to_string(r)
  let parsed = @src.ig_parse(text)
  assert_eq(parsed[0].sequence(), long_seq)
}

// ============================================================================
// 8. SeqRecord interop
// ============================================================================

///|
test "ig_to_seq_record_basic" {
  let r = @src.IgRecord::new("A_U455", "comment", "ATGCATGC")
  let sr = @src.ig_to_seq_record(r)
  assert_eq(sr.id, "A_U455")
  assert_eq(sr.name, "A_U455")
  assert_eq(sr.seq.to_string(), "ATGCATGC")
}

///|
test "ig_to_seq_records_multiple" {
  let r1 = @src.IgRecord::new("A", "", "ATGC")
  let r2 = @src.IgRecord::new("B", "", "ATGCAT")
  let records = [r1, r2]
  let srs = @src.ig_to_seq_records(records)
  assert_eq(srs.length(), 2)
  assert_eq(srs[0].id, "A")
  assert_eq(srs[1].id, "B")
}

///|
test "ig_from_seq_record_basic" {
  let sr = @src.SeqRecord::new(
    @src.Seq::new("ATGCATGC"),
    id="A_U455",
    name="A_U455",
    description="A_U455",
  )
  let r = @src.ig_from_seq_record(sr)
  assert_eq(r.title(), "A_U455")
  assert_eq(r.sequence(), "ATGCATGC")
}

///|
test "ig_from_seq_record_with_description_as_comment" {
  let sr = @src.SeqRecord::new(
    @src.Seq::new("ATGCATGC"),
    id="A_U455",
    name="A_U455",
    description="Some description",
  )
  let r = @src.ig_from_seq_record(sr)
  assert_eq(r.comment(), "Some description")
}

///|
test "ig_from_seq_records_multiple" {
  let sr1 = @src.SeqRecord::new(
    @src.Seq::new("ATGC"),
    id="A",
    name="A",
    description="A",
  )
  let sr2 = @src.SeqRecord::new(
    @src.Seq::new("ATGCAT"),
    id="B",
    name="B",
    description="B",
  )
  let records = @src.ig_from_seq_records([sr1, sr2])
  assert_eq(records.length(), 2)
  assert_eq(records[0].title(), "A")
  assert_eq(records[1].title(), "B")
}

///|
test "ig_seq_record_round_trip" {
  let original = @src.IgRecord::new("X", "", "ATGCATGC")
  let sr = @src.ig_to_seq_record(original)
  let back = @src.ig_from_seq_record(sr)
  assert_eq(back.title(), original.title())
  assert_eq(back.sequence(), original.sequence())
}

// ============================================================================
// 9. Sample data
// ============================================================================

///|
test "ig_sample_text_has_records" {
  let text = @src.ig_sample_text()
  let records = @src.ig_parse(text)
  assert_true(records.length() >= 2)
}

///|
test "ig_sample_records_nonempty" {
  let records = @src.ig_sample_records()
  assert_true(records.length() >= 2)
  assert_eq(records[0].title(), "A_U455")
  assert_eq(records[1].title(), "B_HXB2R")
}

///|
test "ig_sample_records_have_sequences" {
  let records = @src.ig_sample_records()
  for i in 0..<records.length() {
    assert_true(records[i].sequence().length() > 0)
  }
}

///|
test "ig_sample_file_header_has_two_lines" {
  let header = @src.ig_sample_file_header()
  assert_eq(header.length(), 2)
  // Both lines should start with ";;".
  for i in 0..<header.length() {
    assert_true(header[i].length() >= 2)
    assert_eq(header[i].unsafe_get(0), ';'.to_int().to_uint16())
    assert_eq(header[i].unsafe_get(1), ';'.to_int().to_uint16())
  }
}

///|
test "ig_sample_record_basic" {
  let r = @src.ig_sample_record()
  assert_eq(r.title(), "A_U455")
  assert_true(r.sequence().length() > 0)
  assert_true(r.comment().length() > 0)
}

///|
test "ig_sample_record_2_basic" {
  let r = @src.ig_sample_record_2()
  assert_eq(r.title(), "B_HXB2R")
  assert_true(r.sequence().length() > 0)
  assert_true(r.comment().length() > 0)
}

///|
test "ig_sample_records_round_trip" {
  let records = @src.ig_sample_records()
  let header = @src.ig_sample_file_header()
  let text = @src.ig_records_to_string(records, Some(header))
  let parsed = @src.ig_parse(text)
  assert_eq(parsed.length(), records.length())
  for i in 0..<parsed.length() {
    assert_eq(parsed[i].title(), records[i].title())
    assert_eq(parsed[i].sequence(), records[i].sequence())
  }
}

// ============================================================================
// 10. Integration with sample text
// ============================================================================

///|
test "ig_sample_text_starts_with_file_header" {
  let text = @src.ig_sample_text()
  assert_true(text.length() >= 2)
  assert_eq(text.unsafe_get(0), ';'.to_int().to_uint16())
  assert_eq(text.unsafe_get(1), ';'.to_int().to_uint16())
}

///|
test "ig_parse_sample_text_finds_expected_titles" {
  let records = @src.ig_parse(@src.ig_sample_text())
  let titles = records.map(fn(r) { r.title() })
  assert_true(titles.contains("A_U455"))
  assert_true(titles.contains("B_HXB2R"))
}
