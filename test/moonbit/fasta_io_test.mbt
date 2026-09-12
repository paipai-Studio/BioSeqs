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

// ============================================================================
// QUAL format tests
// ============================================================================

///|
fn close_double(a : Double, b : Double, tol? : Double = 1.0e-6) -> Bool {
  (a - b).abs() < tol
}

///|
test "parse_qual_basic" {
  let qual = ">seq1\n10 20 30 40\n>seq2\n50 60 70 80\n"
  let records = @src.parse_qual(qual)
  assert_eq(records.length(), 2)
  assert_eq(records[0].id, "seq1")
  assert_eq(records[0].seq.length(), 4)
  assert_eq(records[0].seq.to_string(), "NNNN")
  let quals0 = records[0].letter_annotations["phred_quality"]
  assert_eq(quals0, [10, 20, 30, 40])
  let quals1 = records[1].letter_annotations["phred_quality"]
  assert_eq(quals1, [50, 60, 70, 80])
}

///|
test "parse_qual_multiline_scores" {
  let qual = ">seq1\n10 20 30 40\n50 60\n>seq2\n70 80\n"
  let records = @src.parse_qual(qual)
  assert_eq(records.length(), 2)
  assert_eq(records[0].letter_annotations["phred_quality"], [
    10, 20, 30, 40, 50, 60,
  ])
  assert_eq(records[1].letter_annotations["phred_quality"], [70, 80])
}

///|
test "parse_qual_empty" {
  let records = @src.parse_qual("")
  assert_eq(records.length(), 0)
}

///|
test "write_qual_roundtrip" {
  let qual = ">seq1\n10 20 30 40\n>seq2\n50 60 70 80\n"
  let records = @src.parse_qual(qual)
  let written = @src.write_qual(records)
  let parsed_again = @src.parse_qual(written)
  assert_eq(parsed_again.length(), 2)
  assert_eq(parsed_again[0].letter_annotations["phred_quality"], [
    10, 20, 30, 40,
  ])
  assert_eq(parsed_again[1].letter_annotations["phred_quality"], [
    50, 60, 70, 80,
  ])
}

///|
test "write_qual_empty_records" {
  let records = @src.parse_qual("")
  let written = @src.write_qual(records)
  assert_eq(written, "")
}

// ============================================================================
// FASTQ variant encoding tests
// ============================================================================

///|
test "parse_fastq_illumina_basic" {
  // Illumina 1.3-1.7: PHRED scores with ASCII offset 64
  // '@' = 64, so '@' -> 0
  let fastq = "@seq1\nACGT\n+\n@ABC\n"
  let records = @src.parse_fastq_illumina(fastq)
  assert_eq(records.length(), 1)
  assert_eq(records[0].id, "seq1")
  assert_eq(records[0].seq.to_string(), "ACGT")
  let quals = records[0].letter_annotations["phred_quality"]
  // '@'=64-64=0, 'A'=65-64=1, 'B'=66-64=2, 'C'=67-64=3
  assert_eq(quals, [0, 1, 2, 3])
}

///|
test "parse_fastq_solexa_basic" {
  // Solexa/early Illumina: Solexa scores with ASCII offset 64
  // ';' = 59, so ';' -> -5 (minimum Solexa score)
  let fastq = "@seq1\nACGT\n+\n;>?@\n"
  let records = @src.parse_fastq_solexa(fastq)
  assert_eq(records.length(), 1)
  assert_eq(records[0].id, "seq1")
  let quals = records[0].letter_annotations["phred_quality"]
  // ';'=59-64=-5, '>'=62-64=-2, '?'=63-64=-1, '@'=64-64=0
  assert_eq(quals, [-5, -2, -1, 0])
}

///|
test "write_fastq_illumina_roundtrip" {
  let fastq = "@seq1\nACGT\n+\n@ABC\n"
  let records = @src.parse_fastq_illumina(fastq)
  let written = @src.write_fastq_illumina(records)
  let parsed_again = @src.parse_fastq_illumina(written)
  assert_eq(parsed_again.length(), 1)
  assert_eq(parsed_again[0].seq.to_string(), "ACGT")
  assert_eq(parsed_again[0].letter_annotations["phred_quality"], [0, 1, 2, 3])
}

// ============================================================================
// Quality score conversion tests
// ============================================================================

///|
test "solexa_quality_from_phred_values" {
  // High scores: approximately equal (but not exactly, due to -1 in formula)
  assert_true(
    close_double(@src.solexa_quality_from_phred(40.0), 40.0, tol=1.0e-3),
  )
  assert_true(
    close_double(@src.solexa_quality_from_phred(30.0), 30.0, tol=5.0e-3),
  )
  assert_true(
    close_double(@src.solexa_quality_from_phred(20.0), 19.96, tol=1.0e-2),
  )
  // Low scores: significant differences
  assert_true(
    close_double(@src.solexa_quality_from_phred(10.0), 9.54, tol=1.0e-2),
  )
  assert_true(
    close_double(@src.solexa_quality_from_phred(5.0), 3.35, tol=1.0e-2),
  )
  assert_true(
    close_double(@src.solexa_quality_from_phred(3.0), -0.02, tol=1.0e-2),
  )
  assert_true(
    close_double(@src.solexa_quality_from_phred(2.0), -2.33, tol=1.0e-2),
  )
  assert_true(
    close_double(@src.solexa_quality_from_phred(1.0), -5.0, tol=1.0e-2),
  )
  // Zero and negative: map to -5.0
  assert_true(
    close_double(@src.solexa_quality_from_phred(0.0), -5.0, tol=1.0e-2),
  )
  assert_true(
    close_double(@src.solexa_quality_from_phred(-1.0), -5.0, tol=1.0e-2),
  )
}

///|
test "phred_quality_from_solexa_values" {
  // Phred scores are always >= corresponding Solexa scores
  assert_true(
    close_double(@src.phred_quality_from_solexa(40.0), 40.0, tol=1.0e-3),
  )
  assert_true(
    close_double(@src.phred_quality_from_solexa(30.0), 30.0, tol=5.0e-3),
  )
  assert_true(
    close_double(@src.phred_quality_from_solexa(20.0), 20.0, tol=5.0e-2),
  )
  assert_true(
    close_double(@src.phred_quality_from_solexa(10.0), 10.0, tol=5.0e-1),
  )
  // Negative Solexa: map to zero
  assert_true(
    close_double(@src.phred_quality_from_solexa(-5.0), 0.0, tol=1.0e-2),
  )
  assert_true(
    close_double(@src.phred_quality_from_solexa(-6.0), 0.0, tol=1.0e-2),
  )
}

///|
test "solexa_phred_roundtrip" {
  let phred = [40.0, 30.0, 20.0, 10.0, 5.0, 3.0, 2.0, 1.0, 0.0]
  for q in phred {
    let solexa = @src.solexa_quality_from_phred(q)
    let back = @src.phred_quality_from_solexa(solexa)
    // At phred=1.0, solexa hits floor (-5.0), roundtrip gives 0.0
    assert_true(close_double(back, q, tol=1.5))
  }
}

// ============================================================================
// SeqIO unified interface tests (QUAL, FASTQ variants)
// ============================================================================

///|
test "seqio_parse_qual" {
  let qual = ">seq1\n10 20 30\n"
  let records = @src.seqio_parse(qual, "qual")
  assert_eq(records.length(), 1)
  assert_eq(records[0].letter_annotations["phred_quality"], [10, 20, 30])
}

///|
test "seqio_write_qual" {
  let qual = ">seq1\n10 20 30\n"
  let records = @src.seqio_parse(qual, "qual")
  let written = @src.seqio_write(records, "qual")
  assert_true(written.contains(">seq1"))
  assert_true(written.contains("10 20 30"))
}

///|
test "seqio_parse_fastq_illumina" {
  let fastq = "@seq1\nACGT\n+\n@ABC\n"
  let records = @src.seqio_parse(fastq, "fastq-illumina")
  assert_eq(records.length(), 1)
  assert_eq(records[0].letter_annotations["phred_quality"], [0, 1, 2, 3])
}

///|
test "seqio_parse_fastq_solexa" {
  let fastq = "@seq1\nACGT\n+\n;>?@\n"
  let records = @src.seqio_parse(fastq, "fastq-solexa")
  assert_eq(records.length(), 1)
}

///|
test "seqio_write_fastq_illumina" {
  let fastq = "@seq1\nACGT\n+\n@ABC\n"
  let records = @src.seqio_parse(fastq, "fastq-illumina")
  let written = @src.seqio_write(records, "fastq-illumina")
  let parsed_again = @src.seqio_parse(written, "fastq-illumina")
  assert_eq(parsed_again[0].letter_annotations["phred_quality"], [0, 1, 2, 3])
}

///|
test "seqio_convert_fasta_to_qual_errors" {
  // FASTA has no quality scores; converting to QUAL should raise an error
  let fasta = ">seq1\nACGT\n"
  let raised = try {
    let _ = @src.seqio_convert(fasta, "fasta", "qual")
    false
  } catch {
    _ => true
  }
  assert_true(raised)
}

///|
test "seqio_convert_fastq_to_qual" {
  let fastq = "@seq1\nACGT\n+\n!\"#$\n"
  let qual = @src.seqio_convert(fastq, "fastq", "qual")
  let records = @src.seqio_parse(qual, "qual")
  assert_eq(records.length(), 1)
  assert_eq(records[0].letter_annotations["phred_quality"], [0, 1, 2, 3])
}

///|
test "seqio_convert_qual_to_fastq" {
  let qual = ">seq1\n10 20 30 40\n"
  let fastq = @src.seqio_convert(qual, "qual", "fastq")
  let records = @src.seqio_parse(fastq, "fastq")
  assert_eq(records.length(), 1)
  assert_eq(records[0].letter_annotations["phred_quality"], [10, 20, 30, 40])
}

// ============================================================================
// PairedFastaQualIterator tests
// ============================================================================

///|
test "paired_fasta_qual_iterator_basic" {
  let fasta = ">seq1\nACGT\n>seq2\nTGCA\n"
  let qual = ">seq1\n10 20 30 40\n>seq2\n50 60 70 80\n"
  let records = @src.paired_fasta_qual_iterator(fasta, qual)
  assert_eq(records.length(), 2)
  assert_eq(records[0].id, "seq1")
  assert_eq(records[0].seq.to_string(), "ACGT")
  assert_eq(records[0].letter_annotations["phred_quality"], [10, 20, 30, 40])
  assert_eq(records[1].letter_annotations["phred_quality"], [50, 60, 70, 80])
}

///|
test "paired_fasta_qual_mismatched_length" {
  // Sequence length and quality length don't match: record is skipped
  let fasta = ">seq1\nACGT\n>seq2\nACGT\n"
  let qual = ">seq1\n10 20 30 40 50\n" // 5 qualities for 4 bases
  // seq1 should be skipped due to length mismatch
  // seq2 should be in the result but without quality (no match)
  let records = @src.paired_fasta_qual_iterator(fasta, qual)
  assert_eq(records.length(), 1)
  assert_eq(records[0].id, "seq2")
}

///|
test "paired_fasta_qual_no_qual_match" {
  // seq1 only in FASTA, not in QUAL
  let fasta = ">seq1\nACGT\n>seq2\nTGCA\n"
  let qual = ">seq2\n10 20 30 40\n"
  let records = @src.paired_fasta_qual_iterator(fasta, qual)
  assert_eq(records.length(), 2)
  // seq1 has no quality
  assert_eq(records[0].letter_annotations.get("phred_quality"), None)
  // seq2 has quality
  assert_eq(records[1].letter_annotations["phred_quality"], [10, 20, 30, 40])
}
