// Black-box tests for Biopython Bio.Align.phylip-compatible support.

///|
fn align_phylip_test_parse(text : String) -> @src.AlignPhylipAlignment {
  @src.align_phylip_parse(text) catch {
    AlignPhylipError(message) => abort("valid PHYLIP input failed: " + message)
  }
}

///|
fn align_phylip_test_rejects(text : String) -> Bool {
  try {
    ignore(@src.align_phylip_parse(text))
    false
  } catch {
    AlignPhylipError(_) => true
  }
}

///|
fn align_phylip_test_write(
  alignment : @src.AlignPhylipAlignment,
  layout : @src.AlignPhylipLayout,
  block_width : Int,
  group_width : Int,
) -> String {
  @src.align_phylip_write(alignment, layout~, block_width~, group_width~) catch {
    AlignPhylipError(message) => abort("PHYLIP writing failed: " + message)
  }
}

///|
fn align_phylip_test_write_rejects(
  alignment : @src.AlignPhylipAlignment,
  block_width : Int,
  group_width : Int,
) -> Bool {
  try {
    ignore(@src.align_phylip_write(alignment, block_width~, group_width~))
    false
  } catch {
    AlignPhylipError(_) => true
  }
}

///|
fn align_phylip_test_from_rejects(
  ids : Array[String],
  rows : Array[String],
) -> Bool {
  try {
    ignore(@src.align_phylip_from_aligned(ids, rows))
    false
  } catch {
    AlignPhylipError(_) => true
  }
}

///|
fn align_phylip_test_sequence_rejects(id : String, row : String) -> Bool {
  try {
    ignore(@src.AlignPhylipSequence::create(id, row))
    false
  } catch {
    AlignPhylipError(_) => true
  }
}

///|
fn align_phylip_test_sample() -> @src.AlignPhylipAlignment {
  align_phylip_test_parse(@src.align_phylip_example_text())
}

///|
fn align_phylip_test_sequential_text() -> String {
  "3 12\n" +
  "Reference MKT--AIL\n" +
  "          GHAA\n" +
  "Query_one M-TQQAIL\n" +
  "          G-AA\n" +
  "Query_two MKT--AVL\n" +
  "          GHAA\n"
}

///|
fn align_phylip_test_single_line_text() -> String {
  "3 12\n" +
  "Reference MKT--AILGHAA\n" +
  "Query_one M-TQQAILG-AA\n" +
  "Query_two MKT--AVLGHAA\n"
}

///|
fn align_phylip_test_two_dat_style() -> String {
  "5 20\n" +
  "Tax1      CCATC TCACG\n" +
  "Tax2      CCATC TCACG\n" +
  "Tax3      CCATC TCCCG\n" +
  "Tax4      TCATC TCATG\n" +
  "Tax5      CCATC TCACG\n" +
  "\n" +
  "GTCG GTACGA\n" +
  "GTCA GTAAGA\n" +
  "CTCA GTAAGA\n" +
  "GTCA ATAAGA\n" +
  "GTCG GTAAGA\n"
}

///|
test "Bio.Align.phylip parses declared shape" {
  let alignment = align_phylip_test_sample()
  assert_eq(alignment.num_sequences(), 3)
  assert_eq(alignment.alignment_length(), 12)
}

///|
test "Bio.Align.phylip detects interleaved layout" {
  assert_eq(align_phylip_test_sample().source_layout, @src.PhylipInterleaved)
}

///|
test "Bio.Align.phylip detects wrapped sequential layout" {
  let alignment = align_phylip_test_parse(align_phylip_test_sequential_text())
  assert_eq(alignment.source_layout, @src.PhylipSequential)
  assert_eq(alignment.num_sequences(), 3)
}

///|
test "Bio.Align.phylip treats complete one-line rows as sequential" {
  let alignment = align_phylip_test_parse(align_phylip_test_single_line_text())
  assert_eq(alignment.source_layout, @src.PhylipSequential)
  assert_eq(alignment.alignment_length(), 12)
}

///|
test "Bio.Align.phylip joins interleaved blocks" {
  let alignment = align_phylip_test_sample()
  assert_eq(alignment.sequences[0].aligned_sequence, "MKT--AILGHAA")
  assert_eq(alignment.sequences[1].aligned_sequence, "M-TQQAILG-AA")
  assert_eq(alignment.sequences[2].aligned_sequence, "MKT--AVLGHAA")
}

///|
test "Bio.Align.phylip joins wrapped sequential rows" {
  let alignment = align_phylip_test_parse(align_phylip_test_sequential_text())
  assert_eq(alignment.sequences[0].aligned_sequence, "MKT--AILGHAA")
  assert_eq(alignment.sequences[1].aligned_sequence, "M-TQQAILG-AA")
  assert_eq(alignment.sequences[2].aligned_sequence, "MKT--AVLGHAA")
}

///|
test "Bio.Align.phylip derives ungapped biological sequences" {
  let alignment = align_phylip_test_sample()
  assert_eq(alignment.sequences[0].sequence, "MKTAILGHAA")
  assert_eq(alignment.sequences[1].sequence, "MTQQAILGAA")
}

///|
test "Bio.Align.phylip parses grouped official-style blocks" {
  let alignment = align_phylip_test_parse(align_phylip_test_two_dat_style())
  assert_eq(alignment.source_layout, @src.PhylipInterleaved)
  assert_eq(alignment.num_sequences(), 5)
  assert_eq(alignment.alignment_length(), 20)
  assert_eq(alignment.sequences[0].aligned_sequence, "CCATCTCACGGTCGGTACGA")
}

///|
test "Bio.Align.phylip preserves fixed-field identifiers" {
  let alignment = align_phylip_test_parse(align_phylip_test_two_dat_style())
  assert_eq(alignment.sequences[0].id, "Tax1")
  assert_eq(alignment.sequences[4].id, "Tax5")
}

///|
test "Bio.Align.phylip accepts internal spaces in fixed identifiers" {
  let alignment = align_phylip_test_parse(
    "2 4\nTax one   ACGT\nTax two   AC-T\n",
  )
  assert_eq(alignment.sequences[0].id, "Tax one")
  assert_eq(alignment.sequences[1].id, "Tax two")
}

///|
test "Bio.Align.phylip accepts empty fixed identifiers" {
  let alignment = align_phylip_test_parse(
    "2 4\n          ACGT\nTaxon     AC-T\n",
  )
  assert_eq(alignment.sequences[0].id, "")
  assert_eq(alignment.sequences[0].sequence, "ACGT")
}

///|
test "Bio.Align.phylip retains duplicate identifiers" {
  let alignment = align_phylip_test_parse(
    "2 4\nduplicate ACGT\nduplicate AC-T\n",
  )
  assert_eq(alignment.find_sequences("duplicate"), [0, 1])
}

///|
test "Bio.Align.phylip accepts CRLF input" {
  let text = @src.align_phylip_example_text().replace(old="\n", new="\r\n")
  let alignment = align_phylip_test_parse(text)
  assert_eq(alignment.num_sequences(), 3)
  assert_eq(alignment.alignment_length(), 12)
}

///|
test "Bio.Align.phylip ignores trailing blank lines" {
  let alignment = align_phylip_test_parse(
    align_phylip_test_single_line_text() + "\n\n",
  )
  assert_eq(alignment.num_sequences(), 3)
}

///|
test "Bio.Align.phylip parses one sequence" {
  let alignment = align_phylip_test_parse("1 4\nTaxon     ACGT\n")
  assert_eq(alignment.num_sequences(), 1)
  assert_eq(alignment.sequences[0].sequence, "ACGT")
}

///|
test "Bio.Align.phylip reports compact summary" {
  assert_eq(
    align_phylip_test_sample().summary(),
    "AlignPhylipAlignment(layout=interleaved, sequences=3, columns=12, normalized_collisions=false)",
  )
}

///|
test "Bio.Align.phylip returns alignment columns" {
  let alignment = align_phylip_test_sample()
  assert_eq(alignment.column(0), Some("MMM"))
  assert_eq(alignment.column(3), Some("-Q-"))
  assert_eq(alignment.column(99), None)
}

///|
test "Bio.Align.phylip maps residue positions to columns" {
  let alignment = align_phylip_test_sample()
  assert_eq(alignment.sequence_position_to_column(0, 0), Some(0))
  assert_eq(alignment.sequence_position_to_column(0, 3), Some(5))
  assert_eq(alignment.sequence_position_to_column(1, 2), Some(3))
}

///|
test "Bio.Align.phylip maps columns to residue positions" {
  let alignment = align_phylip_test_sample()
  assert_eq(alignment.column_to_sequence_position(0, 3), None)
  assert_eq(alignment.column_to_sequence_position(1, 3), Some(2))
  assert_eq(alignment.column_to_sequence_position(0, 5), Some(3))
}

///|
test "Bio.Align.phylip rejects coordinate bounds with None" {
  let alignment = align_phylip_test_sample()
  assert_eq(alignment.sequence_position_to_column(-1, 0), None)
  assert_eq(alignment.sequence_position_to_column(0, 99), None)
  assert_eq(alignment.column_to_sequence_position(9, 0), None)
  assert_eq(alignment.column_to_sequence_position(0, -1), None)
}

///|
test "Bio.Align.phylip projects positions between rows" {
  let alignment = align_phylip_test_sample()
  assert_eq(alignment.map_position(0, 1, 0), Some(0))
  assert_eq(alignment.map_position(0, 1, 1), None)
  assert_eq(alignment.map_position(0, 1, 3), Some(4))
}

///|
test "Bio.Align.phylip returns aligned position pairs" {
  let pairs = align_phylip_test_sample().aligned_pairs(0, 1) catch {
    AlignPhylipError(message) => abort(message)
  }
  assert_eq(pairs.length(), 12)
  assert_eq(pairs[0], (Some(0), Some(0)))
  assert_eq(pairs[1], (Some(1), None))
  assert_eq(pairs[3], (None, Some(2)))
}

///|
test "Bio.Align.phylip emits compressed coordinate paths" {
  let paths = align_phylip_test_sample().coordinate_path()
  assert_eq(paths[0], [0, 1, 2, 3, 3, 7, 8, 10])
  assert_eq(paths[1], [0, 1, 1, 2, 4, 8, 8, 10])
  assert_eq(paths[2], [0, 1, 2, 3, 3, 7, 8, 10])
}

///|
test "Bio.Align.phylip computes pair counts" {
  let counts = align_phylip_test_sample().pair_counts(0, 1) catch {
    AlignPhylipError(message) => abort(message)
  }
  assert_eq(counts.pairs, 1)
  assert_eq(counts.aligned, 8)
  assert_eq(counts.identities, 8)
  assert_eq(counts.mismatches, 0)
  assert_eq(counts.gap_columns, 4)
  assert_eq(counts.double_gap_columns, 0)
  assert_eq(counts.gap_opens, 3)
}

///|
test "Bio.Align.phylip computes exact pair identity" {
  let counts = align_phylip_test_sample().pair_counts(0, 1) catch {
    AlignPhylipError(message) => abort(message)
  }
  assert_eq(counts.identity(), 1.0)
}

///|
test "Bio.Align.phylip counts mismatches case-sensitively" {
  let alignment = @src.align_phylip_from_aligned(["first", "second"], [
    "ACGT", "AcGA",
  ]) catch {
    AlignPhylipError(message) => abort(message)
  }
  let counts = alignment.pair_counts(0, 1) catch {
    AlignPhylipError(message) => abort(message)
  }
  assert_eq(counts.identities, 2)
  assert_eq(counts.mismatches, 2)
}

///|
test "Bio.Align.phylip aggregates all unordered pairs" {
  let counts = align_phylip_test_sample().counts()
  assert_eq(counts.pairs, 3)
  assert_eq(counts.aligned, 26)
  assert_eq(counts.identities, 24)
  assert_eq(counts.mismatches, 2)
  assert_eq(counts.gap_columns, 8)
  assert_eq(counts.gap_opens, 6)
}

///|
test "Bio.Align.phylip reports per-column occupancy" {
  let occupancy = align_phylip_test_sample().occupancy()
  assert_eq(occupancy.length(), 12)
  assert_eq(occupancy[0], 1.0)
  assert_eq(occupancy[1], 2.0 / 3.0)
  assert_eq(occupancy[3], 1.0 / 3.0)
}

///|
test "Bio.Align.phylip calculates majority consensus" {
  let consensus = align_phylip_test_sample().majority_consensus() catch {
    AlignPhylipError(message) => abort(message)
  }
  assert_eq(consensus, "MKTQQAILGHAA")
}

///|
test "Bio.Align.phylip applies consensus threshold" {
  let alignment = @src.align_phylip_from_aligned(["one", "two", "three"], [
    "ACGT", "ATGT", "AGGT",
  ]) catch {
    AlignPhylipError(message) => abort(message)
  }
  let consensus = alignment.majority_consensus(minimum_fraction=0.5) catch {
    AlignPhylipError(message) => abort(message)
  }
  assert_eq(consensus, "AXGT")
}

///|
test "Bio.Align.phylip normalizes punctuation in identifiers" {
  let name = @src.align_phylip_normalize_id("A:B;C(D),E") catch {
    AlignPhylipError(message) => abort(message)
  }
  assert_eq(name, "A|B|CDE")
}

///|
test "Bio.Align.phylip truncates normalized identifiers to ten columns" {
  let name = @src.align_phylip_normalize_id("Taxon[1]:extra") catch {
    AlignPhylipError(message) => abort(message)
  }
  assert_eq(name, "Taxon1|ext")
}

///|
test "Bio.Align.phylip detects normalization collisions" {
  let alignment = @src.align_phylip_from_aligned(
    ["B_subtilis_YXEM", "B_subtilis_GlnH"],
    ["ACGT", "AC-T"],
  ) catch {
    AlignPhylipError(message) => abort(message)
  }
  assert_true(alignment.has_normalized_id_collisions())
  assert_eq(alignment.normalized_ids(), ["B_subtilis", "B_subtilis"])
}

///|
test "Bio.Align.phylip normalized identifiers are defensive copies" {
  let alignment = align_phylip_test_sample()
  let names = alignment.normalized_ids()
  names[0] = "changed"
  assert_eq(alignment.normalized_ids()[0], "Reference")
}

///|
test "Bio.Align.phylip constructor copies input arrays" {
  let ids = ["one", "two"]
  let rows = ["ACGT", "AC-T"]
  let alignment = @src.align_phylip_from_aligned(ids, rows) catch {
    AlignPhylipError(message) => abort(message)
  }
  ids[0] = "changed"
  rows[0] = "TTTT"
  assert_eq(alignment.sequences[0].id, "one")
  assert_eq(alignment.sequences[0].aligned_sequence, "ACGT")
}

///|
test "Bio.Align.phylip coordinate paths are defensive copies" {
  let alignment = align_phylip_test_sample()
  let paths = alignment.coordinate_path()
  paths[0][0] = 99
  assert_eq(alignment.coordinate_path()[0][0], 0)
}

///|
test "Bio.Align.phylip occupancy arrays are defensive copies" {
  let alignment = align_phylip_test_sample()
  let occupancy = alignment.occupancy()
  occupancy[0] = 0.0
  assert_eq(alignment.occupancy()[0], 1.0)
}

///|
test "Bio.Align.phylip default writer matches modern sequential form" {
  let text = align_phylip_test_write(
    align_phylip_test_sample(),
    @src.PhylipSequential,
    0,
    0,
  )
  assert_true(text.has_prefix("3 12\nReference MKT--AILGHAA\n"))
  assert_true(text.contains("Query_one M-TQQAILG-AA\n"))
}

///|
test "Bio.Align.phylip writer normalizes long identifiers" {
  let alignment = @src.align_phylip_from_aligned(["Taxon[1]:extra", "short"], [
    "ACGT", "AC-T",
  ]) catch {
    AlignPhylipError(message) => abort(message)
  }
  let text = align_phylip_test_write(alignment, @src.PhylipSequential, 0, 0)
  assert_true(text.contains("Taxon1|extACGT\n"))
  assert_true(text.contains("short     AC-T\n"))
}

///|
test "Bio.Align.phylip writes wrapped sequential rows" {
  let text = align_phylip_test_write(
    align_phylip_test_sample(),
    @src.PhylipSequential,
    4,
    0,
  )
  assert_true(text.contains("Reference MKT-\n          -AIL\n          GHAA\n"))
  let reparsed = align_phylip_test_parse(text)
  assert_eq(reparsed.sequences[0].aligned_sequence, "MKT--AILGHAA")
}

///|
test "Bio.Align.phylip writes grouped sequential rows" {
  let text = align_phylip_test_write(
    align_phylip_test_sample(),
    @src.PhylipSequential,
    4,
    2,
  )
  assert_true(text.contains("Reference MK T-\n"))
  assert_true(text.contains("          -A IL\n"))
}

///|
test "Bio.Align.phylip writes interleaved blocks" {
  let text = align_phylip_test_write(
    align_phylip_test_sample(),
    @src.PhylipInterleaved,
    4,
    0,
  )
  assert_true(text.contains("Reference MKT-\nQuery_one M-TQ\nQuery_two MKT-\n"))
  assert_true(text.contains("\n          -AIL\n          QAIL\n"))
  let reparsed = align_phylip_test_parse(text)
  assert_eq(reparsed.source_layout, @src.PhylipInterleaved)
}

///|
test "Bio.Align.phylip writes grouped interleaved blocks" {
  let text = align_phylip_test_write(
    align_phylip_test_sample(),
    @src.PhylipInterleaved,
    6,
    3,
  )
  assert_true(text.contains("Reference MKT --A\n"))
  assert_true(text.contains("Query_one M-T QQA\n"))
}

///|
test "Bio.Align.phylip canonical sequential round trip is stable" {
  let first = align_phylip_test_write(
    align_phylip_test_sample(),
    @src.PhylipSequential,
    0,
    0,
  )
  let second = align_phylip_test_write(
    align_phylip_test_parse(first),
    @src.PhylipSequential,
    0,
    0,
  )
  assert_eq(second, first)
}

///|
test "Bio.Align.phylip canonical interleaved round trip is stable" {
  let first = align_phylip_test_write(
    align_phylip_test_sample(),
    @src.PhylipInterleaved,
    4,
    2,
  )
  let second = align_phylip_test_write(
    align_phylip_test_parse(first),
    @src.PhylipInterleaved,
    4,
    2,
  )
  assert_eq(second, first)
}

///|
test "Bio.Align.phylip writer retains legal duplicate names" {
  let alignment = align_phylip_test_parse(
    "2 4\nduplicate ACGT\nduplicate AC-T\n",
  )
  let text = align_phylip_test_write(alignment, @src.PhylipSequential, 0, 0)
  assert_eq(text, "2 4\nduplicate ACGT\nduplicate AC-T\n")
}

///|
test "Bio.Align.phylip rejects empty input" {
  assert_true(align_phylip_test_rejects(""))
}

///|
test "Bio.Align.phylip rejects whitespace-only input" {
  assert_true(align_phylip_test_rejects("\n"))
}

///|
test "Bio.Align.phylip rejects malformed headers" {
  assert_true(align_phylip_test_rejects("two 4\nTaxon     ACGT\n"))
  assert_true(align_phylip_test_rejects("2\nTaxon     ACGT\n"))
  assert_true(align_phylip_test_rejects("2 4 extra\nTaxon     ACGT\n"))
}

///|
test "Bio.Align.phylip rejects non-positive dimensions" {
  assert_true(align_phylip_test_rejects("0 4\n"))
  assert_true(align_phylip_test_rejects("1 0\nTaxon     ACGT\n"))
  assert_true(align_phylip_test_rejects("-1 4\nTaxon     ACGT\n"))
}

///|
test "Bio.Align.phylip rejects integer overflow" {
  assert_true(align_phylip_test_rejects("2147483648 4\nTaxon     ACGT\n"))
}

///|
test "Bio.Align.phylip rejects fewer rows than declared" {
  assert_true(align_phylip_test_rejects("2 4\nTaxon     ACGT\n"))
}

///|
test "Bio.Align.phylip rejects empty first-block rows" {
  assert_true(
    align_phylip_test_rejects("2 4\nTaxon     ACGT\n\nOther     AC-T\n"),
  )
}

///|
test "Bio.Align.phylip rejects rows shorter than identifier field" {
  assert_true(align_phylip_test_rejects("1 4\nTax ACGT\n"))
}

///|
test "Bio.Align.phylip rejects empty sequence segments" {
  assert_true(align_phylip_test_rejects("1 4\nTaxon     \n"))
  assert_true(align_phylip_test_sequence_rejects("Taxon", ""))
}

///|
test "Bio.Align.phylip rejects dot gap notation" {
  assert_true(align_phylip_test_rejects("1 4\nTaxon     AC.T\n"))
  assert_true(align_phylip_test_sequence_rejects("Taxon", "AC.T"))
}

///|
test "Bio.Align.phylip rejects short sequential rows" {
  assert_true(
    align_phylip_test_rejects("2 8\nTaxon     ACGT\nOther     AC-T\n"),
  )
}

///|
test "Bio.Align.phylip rejects overlong sequential rows" {
  assert_true(
    align_phylip_test_rejects("2 4\nTaxon     ACGTA\nOther     AC-TT\n"),
  )
}

///|
test "Bio.Align.phylip rejects blank lines inside sequential rows" {
  assert_true(
    align_phylip_test_rejects(
      "2 8\nTaxon     ACGT\n\n          ACGT\nOther     AC-TACGT\n",
    ),
  )
}

///|
test "Bio.Align.phylip rejects incomplete interleaved blocks" {
  assert_true(
    align_phylip_test_rejects(
      "3 8\nTaxon1    ACGT\nTaxon2    AC-T\nTaxon3    ACGT\n\nACGT\nACGT\n",
    ),
  )
}

///|
test "Bio.Align.phylip rejects unequal interleaved block widths" {
  assert_true(
    align_phylip_test_rejects(
      "2 8\nTaxon1    ACGT\nTaxon2    AC-T\n\nACGT\nACG\n",
    ),
  )
}

///|
test "Bio.Align.phylip rejects trailing non-whitespace data" {
  assert_true(
    align_phylip_test_rejects(
      align_phylip_test_single_line_text() + "unexpected\n",
    ),
  )
}

///|
test "Bio.Align.phylip rejects all-gap columns" {
  assert_true(align_phylip_test_rejects("2 2\nTaxon1    A-\nTaxon2    C-\n"))
  assert_true(align_phylip_test_from_rejects(["one", "two"], ["A-", "C-"]))
}

///|
test "Bio.Align.phylip rejects unequal constructor row widths" {
  assert_true(align_phylip_test_from_rejects(["one", "two"], ["ACGT", "ACG"]))
}

///|
test "Bio.Align.phylip rejects constructor dimension mismatch" {
  assert_true(align_phylip_test_from_rejects(["one"], ["ACGT", "AC-T"]))
}

///|
test "Bio.Align.phylip rejects empty constructor alignments" {
  assert_true(align_phylip_test_from_rejects([], []))
}

///|
test "Bio.Align.phylip rejects control characters in identifiers" {
  assert_true(align_phylip_test_sequence_rejects("bad\nid", "ACGT"))
}

///|
test "Bio.Align.phylip rejects invalid consensus thresholds" {
  let alignment = align_phylip_test_sample()
  let below = try {
    ignore(alignment.majority_consensus(minimum_fraction=-0.1))
    false
  } catch {
    AlignPhylipError(_) => true
  }
  let above = try {
    ignore(alignment.majority_consensus(minimum_fraction=1.1))
    false
  } catch {
    AlignPhylipError(_) => true
  }
  assert_true(below)
  assert_true(above)
}

///|
test "Bio.Align.phylip rejects invalid pair row indices" {
  let alignment = align_phylip_test_sample()
  let rejected = try {
    ignore(alignment.pair_counts(0, 9))
    false
  } catch {
    AlignPhylipError(_) => true
  }
  assert_true(rejected)
}

///|
test "Bio.Align.phylip rejects invalid aligned-pair row indices" {
  let alignment = align_phylip_test_sample()
  let rejected = try {
    ignore(alignment.aligned_pairs(-1, 0))
    false
  } catch {
    AlignPhylipError(_) => true
  }
  assert_true(rejected)
}

///|
test "Bio.Align.phylip writer rejects negative block widths" {
  assert_true(
    align_phylip_test_write_rejects(align_phylip_test_sample(), -1, 0),
  )
}

///|
test "Bio.Align.phylip writer rejects negative group widths" {
  assert_true(
    align_phylip_test_write_rejects(align_phylip_test_sample(), 4, -1),
  )
}
