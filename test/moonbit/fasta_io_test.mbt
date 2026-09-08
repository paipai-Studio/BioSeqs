// Whitebox tests for FASTA/FASTQ low-level parsers and iterators:
//   - SimpleFastaParser (low-level FASTA parser)
//   - FastqGeneralIterator (low-level FASTQ parser)
//   - FastaTwoLineParser (compact two-line FASTA parser)
//   - FastqTwoLineParser (compact two-line FASTQ parser)
//   - FastaIterator (format-specific iterator)
//   - FastaTwoLineIterator (compact FASTA iterator)
//   - Empty line handling strategies
//   - Error handling with line numbers
// Run `moon test` to execute.

// --- SimpleFastaParser Tests ---

///|
test "simple_fasta_parser_basic" {
  let fasta = ">seq1\nACGT\n>seq2\nTGCA"
  let records = @src.simple_fasta_parser(fasta)
  assert_eq(records.length(), 2)
  assert_eq(records[0].0, "seq1")
  assert_eq(records[0].1, "ACGT")
  assert_eq(records[1].0, "seq2")
  assert_eq(records[1].1, "TGCA")
}

///|
test "simple_fasta_parser_with_spaces" {
  let fasta = ">seq1\nA C G T\n>seq2\nTG CA"
  let records = @src.simple_fasta_parser(fasta)
  assert_eq(records.length(), 2)
  assert_eq(records[0].1, "ACGT")
  assert_eq(records[1].1, "TGCA")
}

///|
test "simple_fasta_parser_empty_lines" {
  let fasta = ">seq1\n\nACGT\n\n>seq2\nTGCA\n"
  let records = @src.simple_fasta_parser(fasta)
  assert_eq(records.length(), 2)
  assert_eq(records[0].1, "ACGT")
  assert_eq(records[1].1, "TGCA")
}

///|
test "simple_fasta_parser_whitespace_only_lines" {
  let fasta = ">seq1\n   \n\t\nACGT\n>seq2\nTGCA"
  let records = @src.simple_fasta_parser(fasta)
  assert_eq(records.length(), 2)
  assert_eq(records[0].1, "ACGT")
  assert_eq(records[1].1, "TGCA")
}

///|
test "simple_fasta_parser_crlf" {
  let fasta = ">seq1\r\nACGT\r\n>seq2\r\nTGCA\r\n"
  let records = @src.simple_fasta_parser(fasta)
  assert_eq(records.length(), 2)
  assert_eq(records[0].0, "seq1")
  assert_eq(records[0].1, "ACGT")
  assert_eq(records[1].0, "seq2")
  assert_eq(records[1].1, "TGCA")
}

///|
test "simple_fasta_parser_empty_input" {
  let records = @src.simple_fasta_parser("")
  assert_eq(records.length(), 0)
}

// --- FastqGeneralIterator Tests ---

///|
test "fastq_general_iterator_basic" {
  let fastq = "@seq1\nACGT\n+\n!''*\n@seq2\nTGCA\n+\n''*!"
  let records = @src.fastq_general_parser(fastq)
  assert_eq(records.length(), 2)
  let (title1, seq1, quals1) = records[0]
  assert_eq(title1, "seq1")
  assert_eq(seq1, "ACGT")
  assert_eq(quals1.length(), 4)
  // '!' = 0, '\'' = 6, '*' = 9 (Phred+33)
  assert_eq(quals1[0], 0)
  assert_eq(quals1[1], 6)
  assert_eq(quals1[3], 9)
  let (title2, seq2, _) = records[1]
  assert_eq(title2, "seq2")
  assert_eq(seq2, "TGCA")
}

///|
test "fastq_general_iterator_raw_sequence_line" {
  let fastq = "@seq1\nA C G T\n+\n!''*\n"
  let records = @src.fastq_general_parser(fastq)
  assert_eq(records.length(), 1)
  // FastqGeneralIterator yields the raw sequence line (no space stripping)
  assert_eq(records[0].1, "A C G T")
}

///|
test "fastq_general_iterator_skips_blank_and_junk_lines" {
  let fastq = "\n@seq1\nACGT\n+\n!''*\n\n@seq2\nTGCA\n+\n''*!"
  let records = @src.fastq_general_parser(fastq)
  assert_eq(records.length(), 2)
  assert_eq(records[0].0, "seq1")
  assert_eq(records[1].0, "seq2")
}

///|
test "fastq_general_iterator_crlf" {
  let fastq = "@seq1\r\nACGT\r\n+\r\n!''*\r\n"
  let records = @src.fastq_general_parser(fastq)
  assert_eq(records.length(), 1)
  assert_eq(records[0].0, "seq1")
  assert_eq(records[0].1, "ACGT")
  assert_eq(records[0].2.length(), 4)
}

///|
test "fastq_general_iterator_incomplete_trailing_record_skipped" {
  let fastq = "@seq1\nACGT\n+\n!''*\n@trunc\nACGT\n+"
  let records = @src.fastq_general_parser(fastq)
  assert_eq(records.length(), 1)
  assert_eq(records[0].0, "seq1")
}

// --- FastaTwoLineParser Tests ---

///|
test "fasta_two_line_parser_basic" {
  let fasta = ">seq1 ACGT\n>seq2 GCTA\n"
  let records = @src.fasta_two_line_parser(fasta)
  assert_eq(records.length(), 2)
  assert_eq(records[0].0, "seq1")
  assert_eq(records[0].1, "ACGT")
  assert_eq(records[1].0, "seq2")
  assert_eq(records[1].1, "GCTA")
}

///|
test "fasta_two_line_parser_no_space" {
  let fasta = ">seq1ACGT\n"
  let records = @src.fasta_two_line_parser(fasta)
  assert_eq(records.length(), 1)
  assert_eq(records[0].0, "seq1ACGT")
  assert_eq(records[0].1, "")
}

///|
test "fasta_two_line_parser_skips_non_header_lines" {
  let fasta = "junk\n>seq1 ACGT\nmore junk\n>seq2 GCTA"
  let records = @src.fasta_two_line_parser(fasta)
  assert_eq(records.length(), 2)
  assert_eq(records[0].0, "seq1")
  assert_eq(records[1].0, "seq2")
}

///|
test "fasta_two_line_parser_crlf" {
  let fasta = ">seq1 ACGT\r\n>seq2 GCTA\r\n"
  let records = @src.fasta_two_line_parser(fasta)
  assert_eq(records.length(), 2)
  assert_eq(records[0].1, "ACGT")
  assert_eq(records[1].1, "GCTA")
}

// --- FastqTwoLineParser Tests ---

///|
test "fastq_two_line_parser_basic" {
  let fastq = "@seq1 ACGT\n!''*\n@seq2 TGCA\n''*!"
  let records = @src.fastq_two_line_parser(fastq)
  assert_eq(records.length(), 2)
  let (title1, seq1, quals1) = records[0]
  assert_eq(title1, "seq1")
  assert_eq(seq1, "ACGT")
  assert_eq(quals1.length(), 4)
  let (title2, seq2, _) = records[1]
  assert_eq(title2, "seq2")
  assert_eq(seq2, "TGCA")
}

///|
test "fastq_two_line_parser_no_space" {
  let fastq = "@seq1ACGT\n!''*\n"
  let records = @src.fastq_two_line_parser(fastq)
  assert_eq(records.length(), 1)
  // No space: whole text after '@' is the title, sequence is empty
  assert_eq(records[0].0, "seq1ACGT")
  assert_eq(records[0].1, "")
}

///|
test "fastq_two_line_parser_skips_non_at_lines" {
  let fastq = "junk\n@seq1 ACGT\n!''*"
  let records = @src.fastq_two_line_parser(fastq)
  assert_eq(records.length(), 1)
  assert_eq(records[0].0, "seq1")
  assert_eq(records[0].1, "ACGT")
}

///|
test "fastq_two_line_parser_crlf" {
  let fastq = "@seq1 ACGT\r\n!''*\r\n"
  let records = @src.fastq_two_line_parser(fastq)
  assert_eq(records.length(), 1)
  assert_eq(records[0].1, "ACGT")
  assert_eq(records[0].2.length(), 4)
}

// --- FastaIterator Tests ---

///|
test "fasta_iterator_basic" {
  let fasta = ">seq1\nACGT\n>seq2\nTGCA"
  let records = @src.fasta_iterator(fasta)
  assert_eq(records.length(), 2)
  assert_eq(records[0].id, "seq1")
  assert_eq(records[0].seq.to_string(), "ACGT")
  assert_eq(records[1].id, "seq2")
  assert_eq(records[1].seq.to_string(), "TGCA")
}

///|
test "fasta_iterator_empty_input" {
  let records = @src.fasta_iterator("")
  assert_eq(records.length(), 0)
}

///|
test "fasta_iterator_description_is_full_title" {
  let records = @src.fasta_iterator(">seq1 some description\nACGT")
  assert_eq(records[0].id, "seq1")
  assert_eq(records[0].description, "seq1 some description")
}

// --- FastaTwoLineIterator Tests ---

///|
test "fasta_two_line_iterator_basic" {
  let fasta = ">seq1 ACGT\n>seq2 GCTA"
  let records = @src.fasta_two_line_iterator(fasta)
  assert_eq(records.length(), 2)
  assert_eq(records[0].id, "seq1")
  assert_eq(records[0].seq.to_string(), "ACGT")
  assert_eq(records[1].id, "seq2")
  assert_eq(records[1].seq.to_string(), "GCTA")
}

// --- Empty Line Handling Strategy Tests ---

///|
test "empty_lines_at_start" {
  let fasta = "\n\n>seq1\nACGT"
  let records = @src.parse_fasta(fasta)
  assert_eq(records.length(), 1)
  assert_eq(records[0].id, "seq1")
}

///|
test "empty_lines_at_end" {
  let fasta = ">seq1\nACGT\n\n\n"
  let records = @src.parse_fasta(fasta)
  assert_eq(records.length(), 1)
}

///|
test "empty_lines_between_records" {
  let fasta = ">seq1\nACGT\n\n\n>seq2\nTGCA"
  let records = @src.parse_fasta(fasta)
  assert_eq(records.length(), 2)
}

///|
test "whitespace_only_lines" {
  let fasta = ">seq1\n   \n\t\n  \nACGT\n>seq2\nTGCA"
  let records = @src.parse_fasta(fasta)
  assert_eq(records.length(), 2)
}

///|
test "mixed_whitespace_lines" {
  let fasta = ">seq1\n\r\n \t \r\nACGT"
  let records = @src.parse_fasta(fasta)
  assert_eq(records.length(), 1)
  assert_eq(records[0].seq.to_string(), "ACGT")
}

// --- Error Handling with Line Numbers Tests ---

///|
test "seqio_error_with_line_basic" {
  let err = @src.seqio_error_with_line("Test error message", 42)
  guard err is @src.SeqIOError(msg)
  assert_true(msg.contains("Test error message"))
  assert_true(msg.contains("42"))
}

///|
test "seqio_error_with_line_zero" {
  let err = @src.seqio_error_with_line("Error at start", 0)
  guard err is @src.SeqIOError(msg)
  assert_true(msg.contains("0"))
}

// --- Roundtrip Tests ---

///|
test "fasta_roundtrip_basic" {
  let original = ">seq1\nACGT\n>seq2\nTGCA"
  let records = @src.parse_fasta(original)
  let written = @src.write_fasta(records)
  let parsed_again = @src.parse_fasta(written)
  assert_eq(parsed_again.length(), 2)
  assert_eq(parsed_again[0].id, "seq1")
  assert_eq(parsed_again[1].id, "seq2")
}

///|
test "fasta_roundtrip_with_description" {
  let original = ">seq1 First description\nACGT\n>seq2 Second\nTGCA"
  let records = @src.parse_fasta(original)
  let written = @src.write_fasta(records)
  let parsed_again = @src.parse_fasta(written)
  assert_eq(parsed_again[0].description, "seq1 First description")
  assert_eq(parsed_again[0].seq.to_string(), "ACGT")
}

///|
test "fasta_roundtrip_multiline_sequence" {
  let original = ">seq1\nACGTACGT\nACGTACGT\n>seq2\nTTTT"
  let records = @src.parse_fasta(original)
  let written = @src.write_fasta(records)
  let parsed_again = @src.parse_fasta(written)
  assert_eq(parsed_again[0].seq.to_string(), "ACGTACGTACGTACGT")
}

///|
test "write_fasta_record_custom_width" {
  let records = @src.parse_fasta(">seq1\nACGTACGTAC")
  let written = @src.write_fasta_record(records[0], width=5)
  assert_eq(written, ">seq1\nACGTA\nCGTAC\n")
}

///|
test "write_fasta_records_matches_write_fasta" {
  let records = @src.parse_fasta(">a\nACGT\n>b\nTGCA")
  assert_eq(@src.write_fasta(records), @src.write_fasta_records(records))
}

// --- Edge Cases ---

///|
test "empty_sequence" {
  let fasta = ">seq1\n\n>seq2\nACGT"
  let records = @src.parse_fasta(fasta)
  assert_eq(records.length(), 2)
  assert_eq(records[0].seq.to_string(), "")
}

///|
test "single_header_no_sequence" {
  let fasta = ">seq1"
  let records = @src.parse_fasta(fasta)
  assert_eq(records.length(), 1)
  assert_eq(records[0].seq.to_string(), "")
}

///|
test "special_characters_in_header" {
  let fasta = ">seq1|organism=Homo sapiens|gene=BRCA1\nACGT"
  let records = @src.parse_fasta(fasta)
  assert_eq(records.length(), 1)
  // First whitespace-delimited token becomes the id (Biopython behavior)
  assert_eq(records[0].id, "seq1|organism=Homo")
  assert_eq(records[0].description, "seq1|organism=Homo sapiens|gene=BRCA1")
}

///|
test "tab_in_sequence" {
  let fasta = ">seq1\nAC\tGT"
  let records = @src.parse_fasta(fasta)
  assert_eq(records[0].seq.to_string(), "ACGT")
}
