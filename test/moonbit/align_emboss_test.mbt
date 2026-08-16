// Black-box tests for Biopython Bio.Align.emboss-compatible support.

///|
fn align_emboss_test_parse(text : String) -> @src.AlignEmbossDocument {
  @src.align_emboss_parse(text) catch {
    AlignEmbossError(message) => abort("valid EMBOSS input failed: " + message)
  }
}

///|
fn align_emboss_test_raises(text : String) -> Bool {
  try {
    ignore(@src.align_emboss_parse(text))
    false
  } catch {
    AlignEmbossError(_) => true
  }
}

///|
fn align_emboss_test_sample() -> @src.AlignEmbossDocument {
  align_emboss_test_parse(@src.align_emboss_example_text())
}

///|
fn align_emboss_test_reverse_text() -> String {
  "########################################\n" +
  "# Program: water\n" +
  "# Rundate: Sat 22 Oct 2022 22:56:15\n" +
  "# Commandline: water -sreverse2\n" +
  "# Align_format: srspair\n" +
  "# Report_file: water_reverse4.txt\n" +
  "########################################\n" +
  "#=======================================\n" +
  "#\n" +
  "# Aligned_sequences: 2\n" +
  "# 1: seqA\n" +
  "# 2: seqB\n" +
  "# Matrix: EDNAFULL\n" +
  "# Gap_penalty: 1.0\n" +
  "# Extend_penalty: 0.5\n" +
  "#\n" +
  "# Length: 19\n" +
  "# Identity: 16/19 (84.2%)\n" +
  "# Similarity: 16/19 (84.2%)\n" +
  "# Gaps: 3/19 (15.8%)\n" +
  "# Score: 77.5\n" +
  "#\n" +
  "#=======================================\n" +
  "seqA               3 GGGCCCGGTT-TAAAAAAA     20\n" +
  "                     ||||||||   ||||||||\n" +
  "seqB              19 GGGCCCGG--GTAAAAAAA      3\n" +
  "#---------------------------------------\n" +
  "#---------------------------------------\n"
}

///|
fn align_emboss_test_simple_text() -> String {
  "########################################\n" +
  "# Program: matcher\n" +
  "# Rundate: Tue  8 Dec 2009 11:48:35\n" +
  "# Commandline: matcher -aformat simple\n" +
  "# Align_format: simple\n" +
  "# Report_file: matcher_simple.txt\n" +
  "########################################\n" +
  "#=======================================\n" +
  "#\n" +
  "# Aligned_sequences: 2\n" +
  "# 1: AF069992_1\n" +
  "# 2: CAA85685.1\n" +
  "# Matrix: EBLOSUM62\n" +
  "# Gap_penalty: 14\n" +
  "# Extend_penalty: 4\n" +
  "#\n" +
  "# Length: 16\n" +
  "# Identity: 7/16 (43.8%)\n" +
  "# Similarity: 8/16 (50.0%)\n" +
  "# Gaps: 0/16 ( 0.0%)\n" +
  "# Score: 29\n" +
  "#\n" +
  "#=======================================\n" +
  "AF069992_1        73 GPPPQSPDENRAGESS     88\n" +
  "                     |.||:......|.|||\n" +
  "CAA85685.1        47 GVPPEEAGAAVAAESS     62\n" +
  "#---------------------------------------\n" +
  "#---------------------------------------\n"
}

///|
fn align_emboss_test_multiple_rows_text() -> String {
  "########################################\n" +
  "# Program: alignret\n" +
  "# Align_format: simple\n" +
  "# Report_file: stdout\n" +
  "########################################\n" +
  "#=======================================\n" +
  "#\n" +
  "# Aligned_sequences: 4\n" +
  "# 1: alpha\n" +
  "# 2: beta\n" +
  "# 3: gamma\n" +
  "# 4: delta\n" +
  "# Matrix: EBLOSUM62\n" +
  "# Length: 6\n" +
  "# Identity: 3/6\n" +
  "# Similarity: 3/6\n" +
  "# Gaps: 4/6\n" +
  "# Score: 10\n" +
  "#\n" +
  "#=======================================\n" +
  "alpha              1 AC-GTT      5\n" +
  "beta               1 ACCG-T      5\n" +
  "gamma              1 A--GTT      4\n" +
  "delta              1 AC-G-T      4\n" +
  "                     |  | |\n" +
  "#---------------------------------------\n" +
  "#---------------------------------------\n"
}

///|
fn align_emboss_test_gap_blocks_text() -> String {
  "########################################\n" +
  "# Program: needle\n" +
  "# Align_format: srspair\n" +
  "########################################\n" +
  "#=======================================\n" +
  "# Aligned_sequences: 2\n" +
  "# 1: full\n" +
  "# 2: overhang\n" +
  "# Length: 6\n" +
  "# Identity: 3/6\n" +
  "# Similarity: 3/6\n" +
  "# Gaps: 3/6\n" +
  "#=======================================\n" +
  "full                1 AAA      3\n" +
  "                       \n" +
  "overhang            0 ---      0\n" +
  "\n" +
  "full                4 CCC      6\n" +
  "                     |||\n" +
  "overhang            1 CCC      3\n" +
  "#---------------------------------------\n"
}

///|
test "Bio.Align.emboss parses file metadata" {
  let document = align_emboss_test_sample()
  assert_eq(document.metadata.program, "water")
  assert_eq(document.metadata.rundate, "Wed Jan 16 17:23:19 2002")
  assert_eq(document.metadata.report_file, "stdout")
  assert_eq(document.metadata.align_format, "srspair")
}

///|
test "Bio.Align.emboss joins command continuations" {
  let document = align_emboss_test_sample()
  assert_eq(
    document.metadata.command_line,
    "water -asequence reference.fa -bsequence query.fa",
  )
}

///|
test "Bio.Align.emboss defaults to srspair" {
  let text = @src.align_emboss_example_text().replace_all(
    old="# Align_format: srspair\n",
    new="",
  )
  assert_eq(align_emboss_test_parse(text).metadata.align_format, "srspair")
}

///|
test "Bio.Align.emboss report summary" {
  assert_eq(
    align_emboss_test_sample().summary(),
    "EMBOSS report(program=water, format=srspair, alignments=1)",
  )
}

///|
test "Bio.Align.emboss document accessors" {
  let document = align_emboss_test_sample()
  assert_eq(document.num_alignments(), 1)
  assert_true(document.get(0) is Some(_))
  assert_true(document.get(-1) is None)
  assert_true(document.get(1) is None)
}

///|
test "Bio.Align.emboss parses alignment annotations" {
  let annotations = align_emboss_test_sample().alignments[0].annotations
  assert_eq(annotations.matrix, "EDNAFULL")
  assert_eq(annotations.gap_penalty, Some(10.0))
  assert_eq(annotations.extend_penalty, Some(0.5))
  assert_eq(annotations.length, 18)
  assert_eq(annotations.identity, Some(16))
  assert_eq(annotations.similarity, Some(17))
  assert_eq(annotations.gaps, Some(1))
  assert_eq(annotations.score, Some(72.5))
}

///|
test "Bio.Align.emboss parses row identifiers" {
  let alignment = align_emboss_test_sample().alignments[0]
  assert_eq(alignment.num_sequences(), 2)
  assert_eq(alignment.sequences[0].id, "reference_sequence")
  assert_eq(alignment.sequences[1].id, "query_sequence")
}

///|
test "Bio.Align.emboss stores aligned and ungapped rows" {
  let alignment = align_emboss_test_sample().alignments[0]
  assert_eq(alignment.sequences[0].aligned_sequence, "ACGTTGAGT-CTGGGATG")
  assert_eq(alignment.sequences[0].sequence, "ACGTTGAGTCTGGGATG")
  assert_eq(alignment.sequences[1].aligned_sequence, "ACGTTGAGTACTGGAATG")
  assert_eq(alignment.sequences[1].sequence, "ACGTTGAGTACTGGAATG")
}

///|
test "Bio.Align.emboss normalizes forward boundaries" {
  let alignment = align_emboss_test_sample().alignments[0]
  assert_eq(alignment.sequences[0].start, 78)
  assert_eq(alignment.sequences[0].end, 95)
  assert_eq(alignment.sequences[1].start, 0)
  assert_eq(alignment.sequences[1].end, 18)
  assert_eq(alignment.sequences[0].strand, AlignEmbossForward)
}

///|
test "Bio.Align.emboss preserves consensus columns" {
  let alignment = align_emboss_test_sample().alignments[0]
  assert_eq(alignment.consensus, "||||||||| ||||:|||")
  assert_eq(alignment.consensus.length(), alignment.alignment_length())
}

///|
test "Bio.Align.emboss alignment accessors" {
  let alignment = align_emboss_test_sample().alignments[0]
  assert_eq(alignment.alignment_length(), 18)
  assert_eq(alignment.find_sequence("query_sequence"), Some(1))
  assert_true(alignment.find_sequence("missing") is None)
  assert_eq(alignment.column(9), Some("-A"))
  assert_true(alignment.column(-1) is None)
  assert_true(alignment.column(18) is None)
}

///|
test "Bio.Align.emboss alignment summary" {
  assert_eq(
    align_emboss_test_sample().alignments[0].summary(),
    "EMBOSS alignment(rows=2, columns=18, matrix=EDNAFULL)",
  )
}

///|
test "Bio.Align.emboss maps forward columns to positions" {
  let alignment = align_emboss_test_sample().alignments[0]
  assert_eq(alignment.column_to_sequence_position(0, 0), Some(78))
  assert_eq(alignment.column_to_sequence_position(0, 8), Some(86))
  assert_true(alignment.column_to_sequence_position(0, 9) is None)
  assert_eq(alignment.column_to_sequence_position(0, 10), Some(87))
  assert_eq(alignment.column_to_sequence_position(0, 17), Some(94))
}

///|
test "Bio.Align.emboss maps positions to columns" {
  let alignment = align_emboss_test_sample().alignments[0]
  assert_eq(alignment.sequence_position_to_column(0, 78), Some(0))
  assert_eq(alignment.sequence_position_to_column(0, 87), Some(10))
  assert_eq(alignment.sequence_position_to_column(0, 94), Some(17))
  assert_true(alignment.sequence_position_to_column(0, 77) is None)
  assert_true(alignment.sequence_position_to_column(0, 95) is None)
}

///|
test "Bio.Align.emboss maps positions between rows" {
  let alignment = align_emboss_test_sample().alignments[0]
  assert_eq(alignment.map_position(0, 1, 78), Some(0))
  assert_eq(alignment.map_position(0, 1, 87), Some(10))
  assert_true(alignment.map_position(1, 0, 9) is None)
  assert_true(alignment.map_position(0, 4, 78) is None)
}

///|
test "Bio.Align.emboss aligned pairs use absolute coordinates" {
  let pairs = align_emboss_test_sample().alignments[0].aligned_pairs(0, 1) catch {
    _ => abort("valid row indices should produce aligned pairs")
  }
  assert_eq(pairs.length(), 18)
  assert_eq(pairs[0], (Some(78), Some(0)))
  assert_eq(pairs[9], (None, Some(9)))
  assert_eq(pairs[10], (Some(87), Some(10)))
  assert_eq(pairs[17], (Some(94), Some(17)))
}

///|
test "Bio.Align.emboss rejects invalid aligned-pair row" {
  let alignment = align_emboss_test_sample().alignments[0]
  let raised = try {
    ignore(alignment.aligned_pairs(0, 2))
    false
  } catch {
    AlignEmbossError(_) => true
  }
  assert_true(raised)
}

///|
test "Bio.Align.emboss builds compact coordinate path" {
  let path = align_emboss_test_sample().alignments[0].coordinate_path()
  assert_eq(path.length(), 2)
  assert_eq(path[0], [78, 87, 87, 95])
  assert_eq(path[1], [0, 9, 10, 18])
}

///|
test "Bio.Align.emboss counts aligned residues and identities" {
  let counts = align_emboss_test_sample().alignments[0].pair_counts(0, 1) catch {
    _ => abort("valid rows should be counted")
  }
  assert_eq(counts.columns, 18)
  assert_eq(counts.aligned, 17)
  assert_eq(counts.identities, 16)
  assert_eq(counts.mismatches, 1)
}

///|
test "Bio.Align.emboss counts insertion gaps and opens" {
  let counts = align_emboss_test_sample().alignments[0].pair_counts(0, 1) catch {
    _ => abort("valid rows should be counted")
  }
  assert_eq(counts.insertions, 1)
  assert_eq(counts.deletions, 0)
  assert_eq(counts.insertion_opens, 1)
  assert_eq(counts.deletion_opens, 0)
  assert_eq(counts.double_gap_columns, 0)
}

///|
test "Bio.Align.emboss counts consensus positives" {
  let counts = align_emboss_test_sample().alignments[0].pair_counts(0, 1) catch {
    _ => abort("valid rows should be counted")
  }
  assert_eq(counts.positives, 17)
  assert_eq(counts.identity(), 16.0 / 17.0)
}

///|
test "Bio.Align.emboss reverses gap direction with row order" {
  let counts = align_emboss_test_sample().alignments[0].pair_counts(1, 0) catch {
    _ => abort("valid rows should be counted")
  }
  assert_eq(counts.insertions, 0)
  assert_eq(counts.deletions, 1)
  assert_eq(counts.deletion_opens, 1)
}

///|
test "Bio.Align.emboss parses simple format" {
  let document = align_emboss_test_parse(align_emboss_test_simple_text())
  let alignment = document.alignments[0]
  assert_eq(document.metadata.align_format, "simple")
  assert_eq(alignment.sequences[0].start, 72)
  assert_eq(alignment.sequences[0].end, 88)
  assert_eq(alignment.sequences[1].start, 46)
  assert_eq(alignment.sequences[1].end, 62)
  assert_eq(alignment.annotations.identity, Some(7))
  assert_eq(alignment.annotations.similarity, Some(8))
}

///|
test "Bio.Align.emboss parses multiple sequence simple alignment" {
  let alignment = align_emboss_test_parse(
      align_emboss_test_multiple_rows_text(),
    ).alignments[0]
  assert_eq(alignment.num_sequences(), 4)
  assert_eq(alignment.column(2), Some("-C--"))
  assert_eq(alignment.consensus, "|  | |")
  assert_eq(alignment.coordinate_path().length(), 4)
}

///|
test "Bio.Align.emboss parses gap-only leading block" {
  let alignment = align_emboss_test_parse(align_emboss_test_gap_blocks_text()).alignments[0]
  assert_eq(alignment.sequences[1].aligned_sequence, "---CCC")
  assert_eq(alignment.sequences[1].sequence, "CCC")
  assert_eq(alignment.sequences[1].start, 0)
  assert_eq(alignment.sequences[1].end, 3)
  assert_eq(alignment.coordinate_path()[1], [0, 0, 3])
}

///|
test "Bio.Align.emboss parses reverse strand row" {
  let alignment = align_emboss_test_parse(align_emboss_test_reverse_text()).alignments[0]
  let reverse = alignment.sequences[1]
  assert_eq(reverse.start, 19)
  assert_eq(reverse.end, 2)
  assert_eq(reverse.strand, AlignEmbossReverse)
  assert_eq(reverse.sequence, "GGGCCCGGGTAAAAAAA")
}

///|
test "Bio.Align.emboss maps reverse coordinates" {
  let alignment = align_emboss_test_parse(align_emboss_test_reverse_text()).alignments[0]
  assert_eq(alignment.column_to_sequence_position(1, 0), Some(18))
  assert_eq(alignment.column_to_sequence_position(1, 7), Some(11))
  assert_true(alignment.column_to_sequence_position(1, 8) is None)
  assert_true(alignment.column_to_sequence_position(1, 9) is None)
  assert_eq(alignment.column_to_sequence_position(1, 10), Some(10))
  assert_eq(alignment.column_to_sequence_position(1, 18), Some(2))
}

///|
test "Bio.Align.emboss reverse position lookup" {
  let alignment = align_emboss_test_parse(align_emboss_test_reverse_text()).alignments[0]
  assert_eq(alignment.sequence_position_to_column(1, 18), Some(0))
  assert_eq(alignment.sequence_position_to_column(1, 10), Some(10))
  assert_eq(alignment.sequence_position_to_column(1, 2), Some(18))
  assert_true(alignment.sequence_position_to_column(1, 19) is None)
}

///|
test "Bio.Align.emboss reverse coordinate path decreases" {
  let path = align_emboss_test_parse(align_emboss_test_reverse_text()).alignments[0].coordinate_path()
  assert_eq(path[0], [2, 10, 12, 12, 20])
  assert_eq(path[1], [19, 11, 11, 10, 2])
}

///|
test "Bio.Align.emboss reverse pair counts" {
  let alignment = align_emboss_test_parse(align_emboss_test_reverse_text()).alignments[0]
  let counts = alignment.pair_counts(0, 1) catch {
    _ => abort("reverse rows should be counted")
  }
  assert_eq(counts.aligned, 16)
  assert_eq(counts.identities, 16)
  assert_eq(counts.insertions, 1)
  assert_eq(counts.deletions, 2)
  assert_eq(counts.insertion_opens, 1)
  assert_eq(counts.deletion_opens, 1)
}

///|
test "Bio.Align.emboss accepts CRLF" {
  let crlf = @src.align_emboss_example_text().replace_all(old="\n", new="\r\n")
  assert_eq(align_emboss_test_parse(crlf), align_emboss_test_sample())
}

///|
test "Bio.Align.emboss accepts trailing spaces on separators" {
  let text = @src.align_emboss_example_text().replace_all(
    old="#---------------------------------------\n",
    new="#---------------------------------------       \n",
  )
  assert_eq(align_emboss_test_parse(text).alignments[0].annotations.length, 18)
}

///|
test "Bio.Align.emboss parses nobrief equal annotations" {
  let text = @src.align_emboss_example_text().replace_all(
    old="# Score: 72.5\n",
    new="# Score: 72.5\n# Longest_Identity = 9\n# Shortest_Similarity = 2\n",
  )
  let annotations = align_emboss_test_parse(text).alignments[0].annotations
  assert_eq(annotations.longest_identity, "9")
  assert_eq(annotations.shortest_similarity, "2")
}

///|
test "Bio.Align.emboss writes canonical report" {
  let output = @src.align_emboss_write(align_emboss_test_sample()) catch {
    _ => abort("valid EMBOSS document should serialize")
  }
  assert_true(output.has_prefix("########################################\n"))
  assert_true(output.contains("# Aligned_sequences: 2\n"))
  assert_true(output.contains("# Length: 18\n"))
  assert_true(output.contains("reference_seq"))
}

///|
test "Bio.Align.emboss writer round trip" {
  let document = align_emboss_test_sample()
  let output = @src.align_emboss_write(document) catch {
    _ => abort("valid EMBOSS document should serialize")
  }
  assert_eq(align_emboss_test_parse(output), document)
}

///|
test "Bio.Align.emboss writer wraps sequence blocks" {
  let document = align_emboss_test_sample()
  let output = @src.align_emboss_write(document, line_width=7) catch {
    _ => abort("valid EMBOSS document should wrap")
  }
  assert_eq(align_emboss_test_parse(output), document)
}

///|
test "Bio.Align.emboss reverse writer round trip" {
  let document = align_emboss_test_parse(align_emboss_test_reverse_text())
  let output = @src.align_emboss_write(document, line_width=1) catch {
    _ => abort("reverse EMBOSS document should serialize")
  }
  assert_eq(align_emboss_test_parse(output), document)
}

///|
test "Bio.Align.emboss multiple-row writer round trip" {
  let document = align_emboss_test_parse(align_emboss_test_multiple_rows_text())
  let output = @src.align_emboss_write(document, line_width=3) catch {
    _ => abort("multiple-row EMBOSS document should serialize")
  }
  assert_eq(align_emboss_test_parse(output), document)
}

///|
test "Bio.Align.emboss writer supports multiple alignments" {
  let document = align_emboss_test_sample()
  let doubled = @src.AlignEmbossDocument::create(document.metadata, [
    document.alignments[0],
    document.alignments[0],
  ]) catch {
    _ => abort("valid multi-alignment document should construct")
  }
  let output = @src.align_emboss_write(doubled, line_width=9) catch {
    _ => abort("multi-alignment EMBOSS document should serialize")
  }
  let parsed = align_emboss_test_parse(output)
  assert_eq(parsed.num_alignments(), 2)
  assert_eq(parsed, doubled)
}

///|
test "Bio.Align.emboss writer rejects zero width" {
  let raised = try {
    ignore(@src.align_emboss_write(align_emboss_test_sample(), line_width=0))
    false
  } catch {
    AlignEmbossError(_) => true
  }
  assert_true(raised)
}

///|
test "Bio.Align.emboss writer rejects width above format maximum" {
  let raised = try {
    ignore(@src.align_emboss_write(align_emboss_test_sample(), line_width=51))
    false
  } catch {
    AlignEmbossError(_) => true
  }
  assert_true(raised)
}

///|
test "Bio.Align.emboss constructs metadata" {
  let metadata = @src.AlignEmbossMetadata::create(
    program="needle",
    align_format="pair",
  ) catch {
    _ => abort("valid metadata should construct")
  }
  assert_eq(metadata.program, "needle")
  assert_eq(metadata.align_format, "pair")
}

///|
test "Bio.Align.emboss rejects invalid metadata format" {
  let raised = try {
    ignore(@src.AlignEmbossMetadata::create(align_format="json"))
    false
  } catch {
    AlignEmbossError(_) => true
  }
  assert_true(raised)
}

///|
test "Bio.Align.emboss constructs sequence row" {
  let sequence = @src.AlignEmbossSequence::create("row", "AC-GT", 5, 9) catch {
    _ => abort("valid row should construct")
  }
  assert_eq(sequence.sequence, "ACGT")
  assert_eq(sequence.strand, AlignEmbossForward)
}

///|
test "Bio.Align.emboss constructs reverse sequence row" {
  let sequence = @src.AlignEmbossSequence::create("row", "AC-GT", 9, 5) catch {
    _ => abort("valid reverse row should construct")
  }
  assert_eq(sequence.sequence, "ACGT")
  assert_eq(sequence.strand, AlignEmbossReverse)
}

///|
test "Bio.Align.emboss rejects row span mismatch" {
  let raised = try {
    ignore(@src.AlignEmbossSequence::create("row", "AC-GT", 5, 10))
    false
  } catch {
    AlignEmbossError(_) => true
  }
  assert_true(raised)
}

///|
test "Bio.Align.emboss rejects empty input" {
  assert_true(align_emboss_test_raises(""))
}

///|
test "Bio.Align.emboss rejects missing file divider" {
  assert_true(align_emboss_test_raises("# Program: water\n"))
}

///|
test "Bio.Align.emboss rejects truncated header" {
  assert_true(
    align_emboss_test_raises(
      "########################################\n# Program: water\n",
    ),
  )
}

///|
test "Bio.Align.emboss rejects malformed header line" {
  let text = @src.align_emboss_example_text().replace_all(
    old="# Program: water\n",
    new="Program: water\n",
  )
  assert_true(align_emboss_test_raises(text))
}

///|
test "Bio.Align.emboss rejects command continuation before commandline" {
  let text = @src.align_emboss_example_text().replace_all(
    old="# Program: water\n",
    new="# Program: water\n#    orphan\n",
  )
  assert_true(align_emboss_test_raises(text))
}

///|
test "Bio.Align.emboss rejects invalid align format" {
  let text = @src.align_emboss_example_text().replace_all(
    old="# Align_format: srspair\n",
    new="# Align_format: unknown\n",
  )
  assert_true(align_emboss_test_raises(text))
}

///|
test "Bio.Align.emboss rejects header-only report" {
  assert_true(
    align_emboss_test_raises(
      "########################################\n" +
      "# Program: water\n" +
      "########################################\n",
    ),
  )
}

///|
test "Bio.Align.emboss rejects content before alignment" {
  let text = @src.align_emboss_example_text().replace_all(
    old="#=======================================\n",
    new="unexpected\n#=======================================\n",
  )
  assert_true(align_emboss_test_raises(text))
}

///|
test "Bio.Align.emboss rejects duplicate sequence count" {
  let text = @src.align_emboss_example_text().replace_all(
    old="# Matrix: EDNAFULL\n",
    new="# Aligned_sequences: 2\n# 1: x\n# 2: y\n# Matrix: EDNAFULL\n",
  )
  assert_true(align_emboss_test_raises(text))
}

///|
test "Bio.Align.emboss rejects zero sequence count" {
  let text = @src.align_emboss_example_text().replace_all(
    old="# Aligned_sequences: 2\n",
    new="# Aligned_sequences: 0\n",
  )
  assert_true(align_emboss_test_raises(text))
}

///|
test "Bio.Align.emboss rejects non-contiguous sequence numbering" {
  let text = @src.align_emboss_example_text().replace_all(
    old="# 2: query_sequence\n",
    new="# 3: query_sequence\n",
  )
  assert_true(align_emboss_test_raises(text))
}

///|
test "Bio.Align.emboss rejects whitespace in declared identifier" {
  let text = @src.align_emboss_example_text().replace_all(
    old="# 1: reference_sequence\n",
    new="# 1: reference sequence\n",
  )
  assert_true(align_emboss_test_raises(text))
}

///|
test "Bio.Align.emboss rejects unknown alignment annotation" {
  let text = @src.align_emboss_example_text().replace_all(
    old="# Matrix: EDNAFULL\n",
    new="# Unknown: value\n# Matrix: EDNAFULL\n",
  )
  assert_true(align_emboss_test_raises(text))
}

///|
test "Bio.Align.emboss rejects missing alignment length" {
  let text = @src.align_emboss_example_text().replace_all(
    old="# Length: 18\n",
    new="",
  )
  assert_true(align_emboss_test_raises(text))
}

///|
test "Bio.Align.emboss rejects invalid alignment length" {
  let text = @src.align_emboss_example_text().replace_all(
    old="# Length: 18\n",
    new="# Length: eighteen\n",
  )
  assert_true(align_emboss_test_raises(text))
}

///|
test "Bio.Align.emboss rejects annotation count above length" {
  let text = @src.align_emboss_example_text().replace_all(
    old="# Identity: 16/18 (88.9%)\n",
    new="# Identity: 19/18 (105.6%)\n",
  )
  assert_true(align_emboss_test_raises(text))
}

///|
test "Bio.Align.emboss rejects malformed body line" {
  let text = @src.align_emboss_example_text().replace_all(
    old="reference_seq     79 ACGTTGAGT-CTGGGATG     95\n",
    new="short\n",
  )
  assert_true(align_emboss_test_raises(text))
}

///|
test "Bio.Align.emboss rejects unexpected row identifier" {
  let text = @src.align_emboss_example_text().replace_all(
    old="reference_seq     79",
    new="other_id          79",
  )
  assert_true(align_emboss_test_raises(text))
}

///|
test "Bio.Align.emboss rejects invalid sequence start" {
  let text = @src.align_emboss_example_text().replace_all(
    old="reference_seq     79",
    new="reference_seq     xx",
  )
  assert_true(align_emboss_test_raises(text))
}

///|
test "Bio.Align.emboss rejects invalid sequence end" {
  let text = @src.align_emboss_example_text().replace_all(
    old="ACGTTGAGT-CTGGGATG     95",
    new="ACGTTGAGT-CTGGGATG     xx",
  )
  assert_true(align_emboss_test_raises(text))
}

///|
test "Bio.Align.emboss rejects inconsistent block width" {
  let text = @src.align_emboss_example_text().replace_all(
    old="query_sequence     1 ACGTTGAGTACTGGAATG     18",
    new="query_sequence     1 ACGTTGAGTACTGGAAT      17",
  )
  assert_true(align_emboss_test_raises(text))
}

///|
test "Bio.Align.emboss rejects coordinate consumption mismatch" {
  let text = @src.align_emboss_example_text().replace_all(
    old="reference_seq     79 ACGTTGAGT-CTGGGATG     95",
    new="reference_seq     79 ACGTTGAGT-CTGGGATG     94",
  )
  assert_true(align_emboss_test_raises(text))
}

///|
test "Bio.Align.emboss rejects declared width mismatch" {
  let text = @src.align_emboss_example_text()
    .replace_all(old="# Length: 18\n", new="# Length: 19\n")
    .replace_all(old="# Identity: 16/18 (88.9%)\n", new="# Identity: 16/19\n")
    .replace_all(
      old="# Similarity: 17/18 (94.4%)\n",
      new="# Similarity: 17/19\n",
    )
    .replace_all(old="# Gaps: 1/18 ( 5.6%)\n", new="# Gaps: 1/19\n")
  assert_true(align_emboss_test_raises(text))
}

///|
test "Bio.Align.emboss rejects consensus before sequence block" {
  let text = @src.align_emboss_example_text().replace_all(
    old="reference_seq     79 ACGTTGAGT-CTGGGATG     95\n",
    new="                     ||||||||||||||||||\n" +
      "reference_seq     79 ACGTTGAGT-CTGGGATG     95\n",
  )
  assert_true(align_emboss_test_raises(text))
}

///|
test "Bio.Align.emboss rejects consensus wider than block" {
  let text = @src.align_emboss_example_text().replace_all(
    old="                     ||||||||| ||||:|||\n",
    new="                     ||||||||| ||||:|||X\n",
  )
  assert_true(align_emboss_test_raises(text))
}

///|
test "Bio.Align.emboss rejects identity disagreement" {
  let text = @src.align_emboss_example_text().replace_all(
    old="# Identity: 16/18 (88.9%)\n",
    new="# Identity: 15/18\n",
  )
  assert_true(align_emboss_test_raises(text))
}

///|
test "Bio.Align.emboss rejects similarity disagreement" {
  let text = @src.align_emboss_example_text().replace_all(
    old="# Similarity: 17/18 (94.4%)\n",
    new="# Similarity: 16/18\n",
  )
  assert_true(align_emboss_test_raises(text))
}

///|
test "Bio.Align.emboss rejects gap disagreement" {
  let text = @src.align_emboss_example_text().replace_all(
    old="# Gaps: 1/18 ( 5.6%)\n",
    new="# Gaps: 2/18\n",
  )
  assert_true(align_emboss_test_raises(text))
}

///|
test "Bio.Align.emboss rejects invalid reverse coordinates" {
  let text = align_emboss_test_reverse_text().replace_all(
    old="seqB              19 GGGCCCGG--GTAAAAAAA      3",
    new="seqB              19 GGGCCCGG--GTAAAAAAA      4",
  )
  assert_true(align_emboss_test_raises(text))
}

///|
test "Bio.Align.emboss rejects out-of-range count rows" {
  let alignment = align_emboss_test_sample().alignments[0]
  let raised = try {
    ignore(alignment.pair_counts(-1, 1))
    false
  } catch {
    AlignEmbossError(_) => true
  }
  assert_true(raised)
}

///|
test "Bio.Align.emboss constructor rejects empty document" {
  let metadata = @src.AlignEmbossMetadata::create() catch {
    _ => abort("default metadata should construct")
  }
  let raised = try {
    ignore(@src.AlignEmbossDocument::create(metadata, []))
    false
  } catch {
    AlignEmbossError(_) => true
  }
  assert_true(raised)
}

///|
test "Bio.Align.emboss constructor rejects invalid annotation length" {
  let raised = try {
    ignore(@src.AlignEmbossAnnotations::create(0))
    false
  } catch {
    AlignEmbossError(_) => true
  }
  assert_true(raised)
}

///|
test "Bio.Align.emboss constructor rejects negative gap penalty" {
  let raised = try {
    ignore(@src.AlignEmbossAnnotations::create(5, gap_penalty=Some(-1.0)))
    false
  } catch {
    AlignEmbossError(_) => true
  }
  assert_true(raised)
}

///|
test "Bio.Align.emboss constructor rejects empty row identifier" {
  let raised = try {
    ignore(@src.AlignEmbossSequence::create("", "AC", 0, 2))
    false
  } catch {
    AlignEmbossError(_) => true
  }
  assert_true(raised)
}

///|
test "Bio.Align.emboss constructor rejects empty aligned row" {
  let raised = try {
    ignore(@src.AlignEmbossSequence::create("row", "", 0, 0))
    false
  } catch {
    AlignEmbossError(_) => true
  }
  assert_true(raised)
}
