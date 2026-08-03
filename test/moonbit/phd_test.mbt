///|
/// Tests for Bio.Sequencing.Phd Phred quality file parser.

// ---------------------------------------------------------------------------
// Data structure creation and accessors
// ---------------------------------------------------------------------------

test "phd_base_creation" {
  let b = @src.PhdBase::new("A", 35, 1)
  assert_eq(b.base(), "A")
  assert_eq(b.quality(), 35)
  assert_eq(b.peak_position(), 1)
}

test "phd_comment_creation" {
  let c = @src.PhdComment::new(
    "read.ab1", "0.020425.c", "/etc/phred.dat",
    0, 1000, 0, 10, 0.05, "term", "big",
  )
  assert_eq(c.chromat_file(), "read.ab1")
  assert_eq(c.phred_version(), "0.020425.c")
  assert_eq(c.phred_parameter_file(), "/etc/phred.dat")
  assert_eq(c.trace_index_min(), 0)
  assert_eq(c.trace_index_max(), 1000)
  assert_eq(c.trim_start(), 0)
  assert_eq(c.trim_end(), 10)
  assert_eq(c.trim_prob(), 0.05)
  assert_eq(c.chemistry(), "term")
  assert_eq(c.dye(), "big")
}

test "phd_read_creation_and_sequence" {
  let bases = [
    @src.PhdBase::new("A", 35, 1),
    @src.PhdBase::new("C", 40, 2),
    @src.PhdBase::new("G", 45, 3),
  ]
  let comment = @src.PhdComment::new(
    "read.ab1", "0.020425.c", "/etc/phred.dat",
    0, 1000, 0, 3, 0.05, "term", "big",
  )
  let read = @src.PhdRead::new("test_read", bases, comment)
  assert_eq(read.name(), "test_read")
  assert_eq(read.length(), 3)
  assert_eq(read.sequence(), "ACG")
  assert_eq(read.quality(), [35, 40, 45])
  assert_eq(read.peak_positions(), [1, 2, 3])
}

test "phd_read_empty_sequence" {
  let comment = @src.PhdComment::new(
    "", "", "", 0, 0, 0, 0, 0.0, "", "",
  )
  let read = @src.PhdRead::new("empty", [], comment)
  assert_eq(read.length(), 0)
  assert_eq(read.sequence(), "")
  assert_eq(read.quality().length(), 0)
}

test "phd_file_creation" {
  let comment = @src.PhdComment::new(
    "read.ab1", "0.020425.c", "", 0, 0, 0, 0, 0.0, "", "",
  )
  let read1 = @src.PhdRead::new("r1", [@src.PhdBase::new("A", 30, 1)], comment)
  let read2 = @src.PhdRead::new("r2", [@src.PhdBase::new("T", 40, 1)], comment)
  let file = @src.PhdFile::new([read1, read2])
  assert_eq(file.num_reads(), 2)
  assert_true(file.get(0).unwrap().name() == "r1")
  assert_true(file.get(1).unwrap().name() == "r2")
  assert_true(file.get(2) is None)
}

// ---------------------------------------------------------------------------
// Sample data
// ---------------------------------------------------------------------------

test "phd_sample_text_has_begin_sequence" {
  let text = @src.phd_sample_text()
  assert_true(text.contains("BEGIN_SEQUENCE"))
  assert_true(text.contains("END_SEQUENCE"))
}

test "phd_sample_text_has_comment_block" {
  let text = @src.phd_sample_text()
  assert_true(text.contains("BEGIN_COMMENT"))
  assert_true(text.contains("END_COMMENT"))
  assert_true(text.contains("CHROMAT_FILE"))
  assert_true(text.contains("PHRED_VERSION"))
}

test "phd_sample_text_has_dna_block" {
  let text = @src.phd_sample_text()
  assert_true(text.contains("BEGIN_DNA"))
  assert_true(text.contains("END_DNA"))
}

// ---------------------------------------------------------------------------
// Parsing
// ---------------------------------------------------------------------------

test "phd_parse_single_read" {
  let text = @src.phd_sample_text()
  let file = @src.phd_parse(text)
  assert_eq(file.num_reads(), 1)
  let read = file.get(0).unwrap()
  assert_eq(read.name(), "sample_read")
  assert_eq(read.length(), 10)
  assert_eq(read.sequence(), "ACGTACGTAC")
}

test "phd_parse_comment_fields" {
  let text = @src.phd_sample_text()
  let file = @src.phd_parse(text)
  let read = file.get(0).unwrap()
  let c = read.comment()
  assert_eq(c.chromat_file(), "sample_read.ab1")
  assert_eq(c.phred_version(), "0.020425.c")
  assert_eq(c.trace_index_min(), 0)
  assert_eq(c.trace_index_max(), 1000)
  assert_eq(c.trim_start(), 0)
  assert_eq(c.trim_end(), 10)
  assert_eq(c.chemistry(), "term")
  assert_eq(c.dye(), "big")
}

test "phd_parse_quality_scores" {
  let text = @src.phd_sample_text()
  let file = @src.phd_parse(text)
  let read = file.get(0).unwrap()
  let quals = read.quality()
  assert_eq(quals.length(), 10)
  assert_eq(quals[0], 35)
  assert_eq(quals[1], 40)
  assert_eq(quals[9], 44)
}

test "phd_parse_peak_positions" {
  let text = @src.phd_sample_text()
  let file = @src.phd_parse(text)
  let read = file.get(0).unwrap()
  let peaks = read.peak_positions()
  assert_eq(peaks.length(), 10)
  assert_eq(peaks[0], 1)
  assert_eq(peaks[9], 10)
}

test "phd_parse_multiple_reads" {
  let text = @src.phd_sample_text() + "\n" + @src.phd_sample_text()
  let file = @src.phd_parse(text)
  assert_eq(file.num_reads(), 2)
}

test "phd_parse_empty_input" {
  let file = @src.phd_parse("")
  assert_eq(file.num_reads(), 0)
}

test "phd_parse_minimal_phd" {
  let text = "BEGIN_SEQUENCE r1\nBEGIN_DNA\nA 10 1\nT 20 2\nEND_DNA\nEND_SEQUENCE\n"
  let file = @src.phd_parse(text)
  assert_eq(file.num_reads(), 1)
  let read = file.get(0).unwrap()
  assert_eq(read.name(), "r1")
  assert_eq(read.length(), 2)
  assert_eq(read.sequence(), "AT")
  assert_eq(read.quality(), [10, 20])
}

test "phd_parse_no_dna_block" {
  let text = "BEGIN_SEQUENCE r1\nBEGIN_COMMENT\nCHROMAT_FILE: r.ab1\nEND_COMMENT\nEND_SEQUENCE\n"
  let file = @src.phd_parse(text)
  assert_eq(file.num_reads(), 1)
  let read = file.get(0).unwrap()
  assert_eq(read.name(), "r1")
  assert_eq(read.length(), 0)
  assert_eq(read.comment().chromat_file(), "r.ab1")
}

test "phd_parse_trim_values" {
  let text = "BEGIN_SEQUENCE r1\nBEGIN_COMMENT\nTRIM: 5 50 0.02\nEND_COMMENT\nBEGIN_DNA\nA 10 1\nEND_DNA\nEND_SEQUENCE\n"
  let file = @src.phd_parse(text)
  let read = file.get(0).unwrap()
  assert_eq(read.comment().trim_start(), 5)
  assert_eq(read.comment().trim_end(), 50)
}

// ---------------------------------------------------------------------------
// Formatting
// ---------------------------------------------------------------------------

test "phd_read_to_string" {
  let text = @src.phd_sample_text()
  let file = @src.phd_parse(text)
  let read = file.get(0).unwrap()
  let s = read.to_string()
  assert_true(s.contains("PhdRead"))
  assert_true(s.contains("sample_read"))
  assert_true(s.contains("ACGTACGTAC"))
}

test "phd_file_to_string" {
  let text = @src.phd_sample_text()
  let file = @src.phd_parse(text)
  let s = file.to_string()
  assert_true(s.contains("PhdFile"))
  assert_true(s.contains("1 reads"))
}
