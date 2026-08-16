///|
fn align_mauve_test_parse(text : String) -> @src.AlignMauveDocument {
  @src.align_mauve_parse(text) catch {
    AlignMauveError(message) => abort("valid XMFA input failed: " + message)
  }
}

///|
fn align_mauve_test_document() -> @src.AlignMauveDocument {
  align_mauve_test_parse(@src.align_mauve_example_text())
}

///|
fn align_mauve_test_rejects(text : String) -> Bool {
  try {
    ignore(@src.align_mauve_parse(text))
    false
  } catch {
    AlignMauveError(_) => true
  }
}

///|
fn align_mauve_test_one_source(body : String) -> String {
  "#FormatVersion Mauve1\n" +
  "#Sequence1File a.fa\n" +
  "#Sequence1Entry 1\n" +
  "#Sequence1Format FastA\n" +
  body
}

///|
fn align_mauve_test_two_sources(body : String) -> String {
  "#FormatVersion Mauve1\n" +
  "#Sequence1File a.fa\n" +
  "#Sequence1Format FastA\n" +
  "#Sequence2File b.fa\n" +
  "#Sequence2Format FastA\n" +
  body
}

///|
fn align_mauve_test_separate_text() -> String {
  "#FormatVersion Mauve1\n" +
  "#Sequence1File  equCab1.fa\n" +
  "#Sequence1Format        FastA\n" +
  "#Sequence2File  canFam2.fa\n" +
  "#Sequence2Format        FastA\n" +
  "#Sequence3File  mm9.fa\n" +
  "#Sequence3Format        FastA\n" +
  "#BackboneFile   separate.xmfa.bbcols\n" +
  "> 1:1-50 - equCab1.fa\n" +
  "TAAGCCCTCCTAGCACACACCCGGAGTGGCC-GGGCCGTAC-TTTCCTTTTC\n" +
  "> 2:1-49 + canFam2.fa\n" +
  "CAAGCCCTGC--GCGCTCAGCCGGAGTGTCCCGGGCCCTGC-TTTCCTTTTC\n" +
  "> 3:1-19 - mm9.fa\n" +
  "---------------------------------GGATCTACTTTTCCTCTTC\n" +
  "=\n" +
  "> 3:20-41 + mm9.fa\n" +
  "CTGGCGTCCGGAGCTGGGACGT\n" +
  "=\n"
}

///|
fn align_mauve_test_projection_text() -> String {
  align_mauve_test_two_sources(
    "> 1:1-5 + a.fa\n" + "AC-GTA\n" + "> 2:10-14 - b.fa\n" + "A-TGCA\n" + "=\n",
  )
}

///|
fn align_mauve_test_conflict_text() -> String {
  align_mauve_test_one_source(
    "> 1:1-3 + a.fa\n" + "AAA\n" + "=\n" + "> 1:2-4 + a.fa\n" + "CCC\n" + "=\n",
  )
}

///|
test "Bio.Align.mauve parses the official combined-file format version" {
  assert_eq(align_mauve_test_document().format_version(), "Mauve1")
}

///|
test "Bio.Align.mauve parses three combined-file sources" {
  assert_eq(align_mauve_test_document().source_count(), 3)
}

///|
test "Bio.Align.mauve derives zero-based combined-file identifiers" {
  assert_eq(align_mauve_test_document().identifiers(), ["0", "1", "2"])
}

///|
test "Bio.Align.mauve preserves combined-file source entries" {
  let source = align_mauve_test_document().sources()[2]
  assert_eq(source.number(), 3)
  assert_eq(source.file(), "combined.fa")
  assert_eq(source.entry(), Some(3))
  assert_eq(source.format(), "FastA")
  assert_eq(source.identifier(), "2")
}

///|
test "Bio.Align.mauve preserves non-sequence metadata" {
  let document = align_mauve_test_document()
  assert_eq(document.metadata_value("FormatVersion"), Some("Mauve1"))
  assert_eq(
    document.metadata_value("BackboneFile"),
    Some("combined.xmfa.bbcols"),
  )
}

///|
test "Bio.Align.mauve parses all official combined-file blocks" {
  assert_eq(align_mauve_test_document().block_count(), 6)
}

///|
test "Bio.Align.mauve parses the first official block shape" {
  let block = align_mauve_test_document().blocks()[0]
  assert_eq(block.row_count(), 3)
  assert_eq(block.width(), 49)
}

///|
test "Bio.Align.mauve normalizes reverse coordinates to half-open intervals" {
  let row = align_mauve_test_document().blocks()[0].rows()[0]
  assert_eq(row.start(), 1)
  assert_eq(row.end(), 49)
  assert_eq(row.strand(), "-")
  assert_eq(row.coordinate_at_boundary(0), 49)
  assert_eq(row.coordinate_at_boundary(49), 1)
}

///|
test "Bio.Align.mauve reconstructs a forward-oriented reverse row sequence" {
  let row = align_mauve_test_document().blocks()[0].rows()[0]
  assert_eq(
    row.source_sequence(),
    "AAAAGGAAAGTACGGCCCGGCCACTCCGGGTGTGTGCTAGGAGGGCTT",
  )
}

///|
test "Bio.Align.mauve preserves a zero-zero all-gap row" {
  let row = align_mauve_test_document().blocks()[0].rows()[1]
  assert_eq(row.start(), 0)
  assert_eq(row.end(), 0)
  assert_eq(row.source_sequence(), "")
  assert_eq(row.coordinates().length(), 50)
  assert_eq(row.coordinate_at_boundary(24), 0)
}

///|
test "Bio.Align.mauve reproduces Biopython compact coordinate paths" {
  let compact = align_mauve_test_document().blocks()[0].compact_coordinates()
  assert_eq(compact[0], [49, 40, 38, 21, 21, 1])
  assert_eq(compact[1], [0, 0, 0, 0, 0, 0])
  assert_eq(compact[2], [1, 10, 10, 27, 28, 48])
}

///|
test "Bio.Align.mauve exposes forward residue coordinates on reverse rows" {
  let row = align_mauve_test_document().blocks()[0].rows()[0]
  assert_eq(row.residue_coordinate(0), Some(48))
  assert_eq(row.residue_coordinate(28), None)
  assert_eq(row.column_for_coordinate(48), Some(0))
  assert_eq(row.column_for_coordinate(39), Some(9))
}

///|
test "Bio.Align.mauve calculates official first-block pairwise counts" {
  let counts = align_mauve_test_document().blocks()[0].counts()
  assert_eq(counts.pair_count(), 3)
  assert_eq(counts.columns(), 147)
  assert_eq(counts.aligned(), 46)
  assert_eq(counts.identities(), 39)
  assert_eq(counts.mismatches(), 7)
  assert_eq(counts.gaps(), 98)
  assert_eq(counts.gap_opens(), 4)
  assert_eq(counts.double_gap_columns(), 3)
}

///|
test "Bio.Align.mauve aggregates document statistics" {
  let summary = @src.align_mauve_summary(align_mauve_test_document())
  assert_eq(summary.sources(), 3)
  assert_eq(summary.blocks(), 6)
  assert_eq(summary.rows(), 8)
  assert_eq(summary.alignment_columns(), 94)
  assert_eq(summary.aligned_pairs(), 46)
  assert_eq(summary.identities(), 39)
  assert_eq(summary.mismatches(), 7)
  assert_eq(summary.gaps(), 98)
}

///|
test "Bio.Align.mauve parses official separate-file identifiers" {
  let document = align_mauve_test_parse(align_mauve_test_separate_text())
  assert_eq(document.identifiers(), ["equCab1.fa", "canFam2.fa", "mm9.fa"])
  assert_eq(document.sources()[0].entry(), None)
}

///|
test "Bio.Align.mauve parses official separate-file block structure" {
  let document = align_mauve_test_parse(align_mauve_test_separate_text())
  assert_eq(document.block_count(), 2)
  assert_eq(document.blocks()[0].row_count(), 3)
  assert_eq(document.blocks()[1].row_count(), 1)
}

///|
test "Bio.Align.mauve preserves separate-file reverse coordinates" {
  let row = align_mauve_test_parse(align_mauve_test_separate_text()).blocks()[0].rows()[0]
  assert_eq(row.identifier(), "equCab1.fa")
  assert_eq(row.coordinate_at_boundary(0), 50)
  assert_eq(row.coordinate_at_boundary(row.aligned_sequence().length()), 0)
}

///|
test "Bio.Align.mauve accepts wrapped alignment rows" {
  let text = align_mauve_test_one_source(
    "> 1:1-6 + a.fa\n" + "AC-\n" + "GTAC\n" + "=\n",
  )
  let row = align_mauve_test_parse(text).blocks()[0].rows()[0]
  assert_eq(row.aligned_sequence(), "AC-GTAC")
  assert_eq(row.source_sequence(), "ACGTAC")
}

///|
test "Bio.Align.mauve accepts CRLF input" {
  let text = "#FormatVersion Mauve1\r\n" +
    "#Sequence1File a.fa\r\n" +
    "#Sequence1Entry 1\r\n" +
    "#Sequence1Format FastA\r\n" +
    "> 1:1-2 + a.fa\r\n" +
    "AC\r\n" +
    "=\r\n"
  assert_eq(align_mauve_test_parse(text).block_count(), 1)
}

///|
test "Bio.Align.mauve accepts a header-only document" {
  let document = align_mauve_test_parse(align_mauve_test_one_source(""))
  assert_eq(document.block_count(), 0)
}

///|
test "Bio.Align.mauve canonical writer round-trips the official fixture" {
  let document = align_mauve_test_document()
  let written = @src.align_mauve_write(document) catch {
    AlignMauveError(message) => abort("valid XMFA write failed: " + message)
  }
  assert_eq(align_mauve_test_parse(written), document)
}

///|
test "Bio.Align.mauve writer wraps rows deterministically" {
  let written = @src.align_mauve_write(
    align_mauve_test_parse(align_mauve_test_projection_text()),
    line_width=2,
  ) catch {
    AlignMauveError(message) => abort("valid wrapped write failed: " + message)
  }
  assert_true(written.contains("AC\n-G\nTA\n"))
  assert_eq(
    align_mauve_test_parse(written),
    align_mauve_test_parse(align_mauve_test_projection_text()),
  )
}

///|
test "Bio.Align.mauve row coordinate arrays are defensive copies" {
  let row = align_mauve_test_document().blocks()[0].rows()[0]
  let coordinates = row.coordinates()
  coordinates[0] = 999
  assert_eq(row.coordinate_at_boundary(0), 49)
}

///|
test "Bio.Align.mauve nested block accessors are defensive copies" {
  let document = align_mauve_test_document()
  let blocks = document.blocks()
  blocks[0].rows[0].coordinates[0] = 999
  assert_eq(document.blocks()[0].rows()[0].coordinate_at_boundary(0), 49)
}

///|
test "Bio.Align.mauve custom block calculates mixed-gap counts" {
  let counts = align_mauve_test_parse(align_mauve_test_projection_text()).blocks()[0].counts()
  assert_eq(counts.pair_count(), 1)
  assert_eq(counts.columns(), 6)
  assert_eq(counts.aligned(), 4)
  assert_eq(counts.identities(), 3)
  assert_eq(counts.mismatches(), 1)
  assert_eq(counts.gaps(), 2)
  assert_eq(counts.gap_opens(), 2)
  assert_eq(counts.double_gap_columns(), 0)
}

///|
test "Bio.Align.mauve indexes all non-empty combined-file rows" {
  let index = @src.align_mauve_index(align_mauve_test_document())
  assert_eq(index.entries().length(), 7)
}

///|
test "Bio.Align.mauve interval index uses half-open overlap semantics" {
  let index = @src.align_mauve_index(align_mauve_test_document())
  let before = index.query("0", 0, 1) catch {
    AlignMauveError(message) => abort("valid index query failed: " + message)
  }
  let inside = index.query("0", 1, 2) catch {
    AlignMauveError(message) => abort("valid index query failed: " + message)
  }
  assert_eq(before.length(), 1)
  assert_eq(before[0].block_index(), 1)
  assert_eq(inside.length(), 1)
  assert_eq(inside[0].block_index(), 0)
  assert_eq(inside[0].strand(), "-")
}

///|
test "Bio.Align.mauve block position projection handles reverse targets" {
  let block = align_mauve_test_parse(align_mauve_test_projection_text()).blocks()[0]
  assert_eq(block.map_position("a.fa", "b.fa", 0), Some(13))
  assert_eq(block.map_position("a.fa", "b.fa", 1), None)
  assert_eq(block.map_position("a.fa", "b.fa", 2), Some(11))
}

///|
test "Bio.Align.mauve document position projection reports block and column" {
  let mappings = @src.align_mauve_map_position(
    align_mauve_test_parse(align_mauve_test_projection_text()),
    "a.fa",
    "b.fa",
    2,
  ) catch {
    AlignMauveError(message) =>
      abort("valid position mapping failed: " + message)
  }
  assert_eq(mappings.length(), 1)
  assert_eq(mappings[0].block_index(), 0)
  assert_eq(mappings[0].column(), 3)
  assert_eq(mappings[0].source_position(), 2)
  assert_eq(mappings[0].target_position(), 11)
}

///|
test "Bio.Align.mauve range projection splits runs at gaps" {
  let mappings = @src.align_mauve_map_range(
    align_mauve_test_parse(align_mauve_test_projection_text()),
    "a.fa",
    "b.fa",
    0,
    5,
  ) catch {
    AlignMauveError(message) => abort("valid range mapping failed: " + message)
  }
  assert_eq(mappings.length(), 2)
  assert_eq(mappings[0].column_start(), 0)
  assert_eq(mappings[0].column_end(), 1)
  assert_eq(mappings[0].source_start(), 0)
  assert_eq(mappings[0].source_end(), 1)
  assert_eq(mappings[0].target_start(), 13)
  assert_eq(mappings[0].target_end(), 14)
  assert_eq(mappings[1].column_start(), 3)
  assert_eq(mappings[1].column_end(), 6)
  assert_eq(mappings[1].source_start(), 2)
  assert_eq(mappings[1].source_end(), 5)
  assert_eq(mappings[1].target_start(), 9)
  assert_eq(mappings[1].target_end(), 12)
  assert_eq(mappings[1].source_strand(), "+")
  assert_eq(mappings[1].target_strand(), "-")
}

///|
test "Bio.Align.mauve reconstructs forward source sequences" {
  let document = align_mauve_test_parse(align_mauve_test_projection_text())
  assert_eq(@src.align_mauve_reconstruct_sequence(document, "a.fa"), "ACGTA")
  assert_eq(
    @src.align_mauve_reconstruct_sequence(document, "b.fa"),
    "NNNNNNNNNTGCAT",
  )
}

///|
test "Bio.Align.mauve reconstruction accepts an IUPAC fill character" {
  let sequence = @src.align_mauve_reconstruct_sequence(
    align_mauve_test_parse(align_mauve_test_projection_text()),
    "b.fa",
    fill="X",
  ) catch {
    AlignMauveError(message) =>
      abort("valid reconstruction fill failed: " + message)
  }
  assert_eq(sequence, "XXXXXXXXXTGCAT")
}

///|
test "Bio.Align.mauve rejects empty input" {
  assert_true(align_mauve_test_rejects(""))
}

///|
test "Bio.Align.mauve rejects input without metadata" {
  assert_true(align_mauve_test_rejects("> 1:1-1 + a.fa\nA\n=\n"))
}

///|
test "Bio.Align.mauve rejects a missing format version" {
  let text = "#Sequence1File a.fa\n" +
    "#Sequence1Entry 1\n" +
    "#Sequence1Format FastA\n"
  assert_true(align_mauve_test_rejects(text))
}

///|
test "Bio.Align.mauve rejects a header without sources" {
  assert_true(align_mauve_test_rejects("#FormatVersion Mauve1\n"))
}

///|
test "Bio.Align.mauve rejects duplicate header keys" {
  let text = "#FormatVersion Mauve1\n" +
    "#FormatVersion Mauve2\n" +
    "#Sequence1File a.fa\n" +
    "#Sequence1Entry 1\n" +
    "#Sequence1Format FastA\n"
  assert_true(align_mauve_test_rejects(text))
}

///|
test "Bio.Align.mauve rejects missing source file metadata" {
  let text = "#FormatVersion Mauve1\n" +
    "#Sequence1Entry 1\n" +
    "#Sequence1Format FastA\n"
  assert_true(align_mauve_test_rejects(text))
}

///|
test "Bio.Align.mauve rejects missing source format metadata" {
  let text = "#FormatVersion Mauve1\n" +
    "#Sequence1File a.fa\n" +
    "#Sequence1Entry 1\n"
  assert_true(align_mauve_test_rejects(text))
}

///|
test "Bio.Align.mauve rejects non-contiguous sequence declarations" {
  let text = "#FormatVersion Mauve1\n" +
    "#Sequence2File b.fa\n" +
    "#Sequence2Entry 1\n" +
    "#Sequence2Format FastA\n"
  assert_true(align_mauve_test_rejects(text))
}

///|
test "Bio.Align.mauve rejects a combined file without entries" {
  let text = "#FormatVersion Mauve1\n" +
    "#Sequence1File combined.fa\n" +
    "#Sequence1Format FastA\n" +
    "#Sequence2File combined.fa\n" +
    "#Sequence2Format FastA\n"
  assert_true(align_mauve_test_rejects(text))
}

///|
test "Bio.Align.mauve rejects duplicate normalized source identifiers" {
  let text = "#FormatVersion Mauve1\n" +
    "#Sequence1File combined.fa\n" +
    "#Sequence1Entry 1\n" +
    "#Sequence1Format FastA\n" +
    "#Sequence2File combined.fa\n" +
    "#Sequence2Entry 1\n" +
    "#Sequence2Format FastA\n"
  assert_true(align_mauve_test_rejects(text))
}

///|
test "Bio.Align.mauve rejects an unexpected sequence header keyword" {
  let text = "#FormatVersion Mauve1\n" +
    "#Sequence1File a.fa\n" +
    "#Sequence1Entry 1\n" +
    "#Sequence1Format FastA\n" +
    "#Sequence1Length 10\n"
  assert_true(align_mauve_test_rejects(text))
}

///|
test "Bio.Align.mauve rejects sequence data before a row description" {
  assert_true(align_mauve_test_rejects(align_mauve_test_one_source("AC\n=\n")))
}

///|
test "Bio.Align.mauve rejects an empty block terminator" {
  assert_true(align_mauve_test_rejects(align_mauve_test_one_source("=\n")))
}

///|
test "Bio.Align.mauve rejects a missing final block terminator" {
  let text = align_mauve_test_one_source("> 1:1-2 + a.fa\nAC\n")
  assert_true(align_mauve_test_rejects(text))
}

///|
test "Bio.Align.mauve rejects metadata after alignment rows" {
  let text = align_mauve_test_one_source(
    "> 1:1-2 + a.fa\n" + "AC\n" + "#BackboneFile late.bbcols\n" + "=\n",
  )
  assert_true(align_mauve_test_rejects(text))
}

///|
test "Bio.Align.mauve rejects rows for undeclared sources" {
  let text = align_mauve_test_one_source("> 2:1-1 + b.fa\nA\n=\n")
  assert_true(align_mauve_test_rejects(text))
}

///|
test "Bio.Align.mauve rejects duplicate source rows in one block" {
  let text = align_mauve_test_one_source(
    "> 1:1-1 + a.fa\n" + "A\n" + "> 1:2-2 + a.fa\n" + "C\n" + "=\n",
  )
  assert_true(align_mauve_test_rejects(text))
}

///|
test "Bio.Align.mauve rejects malformed row descriptions" {
  let text = align_mauve_test_one_source("> 1:1-1 +\nA\n=\n")
  assert_true(align_mauve_test_rejects(text))
}

///|
test "Bio.Align.mauve rejects invalid strands" {
  let text = align_mauve_test_one_source("> 1:1-1 ? a.fa\nA\n=\n")
  assert_true(align_mauve_test_rejects(text))
}

///|
test "Bio.Align.mauve rejects zero starts with nonzero ends" {
  let text = align_mauve_test_one_source("> 1:0-2 + a.fa\nAA\n=\n")
  assert_true(align_mauve_test_rejects(text))
}

///|
test "Bio.Align.mauve rejects descending file intervals" {
  let text = align_mauve_test_one_source("> 1:3-2 + a.fa\nAA\n=\n")
  assert_true(align_mauve_test_rejects(text))
}

///|
test "Bio.Align.mauve rejects coordinate span mismatches" {
  let text = align_mauve_test_one_source("> 1:1-3 + a.fa\nAA\n=\n")
  assert_true(align_mauve_test_rejects(text))
}

///|
test "Bio.Align.mauve rejects invalid sequence residues" {
  let text = align_mauve_test_one_source("> 1:1-2 + a.fa\nAZ\n=\n")
  assert_true(align_mauve_test_rejects(text))
}

///|
test "Bio.Align.mauve rejects internal sequence whitespace" {
  let text = align_mauve_test_one_source("> 1:1-2 + a.fa\nA C\n=\n")
  assert_true(align_mauve_test_rejects(text))
}

///|
test "Bio.Align.mauve rejects unequal row widths" {
  let text = align_mauve_test_two_sources(
    "> 1:1-2 + a.fa\n" + "AC\n" + "> 2:1-2 + b.fa\n" + "A-C\n" + "=\n",
  )
  assert_true(align_mauve_test_rejects(text))
}

///|
test "Bio.Align.mauve rejects all-gap-only blocks" {
  let text = align_mauve_test_one_source("> 1:0-0 + a.fa\n---\n=\n")
  assert_true(align_mauve_test_rejects(text))
}

///|
test "Bio.Align.mauve rejects invalid writer widths" {
  let rejected = try {
    ignore(@src.align_mauve_write(align_mauve_test_document(), line_width=-1))
    false
  } catch {
    AlignMauveError(_) => true
  }
  assert_true(rejected)
}

///|
test "Bio.Align.mauve rejects invalid interval-index queries" {
  let index = @src.align_mauve_index(align_mauve_test_document())
  let empty = try {
    ignore(index.query("0", 2, 2))
    false
  } catch {
    AlignMauveError(_) => true
  }
  let unknown = try {
    ignore(index.query("missing", 0, 1))
    false
  } catch {
    AlignMauveError(_) => true
  }
  assert_true(empty)
  assert_true(unknown)
}

///|
test "Bio.Align.mauve rejects unknown position-projection identifiers" {
  let rejected = try {
    ignore(
      @src.align_mauve_map_position(
        align_mauve_test_document(),
        "missing",
        "1",
        0,
      ),
    )
    false
  } catch {
    AlignMauveError(_) => true
  }
  assert_true(rejected)
}

///|
test "Bio.Align.mauve rejects invalid range projections" {
  let rejected = try {
    ignore(
      @src.align_mauve_map_range(align_mauve_test_document(), "0", "1", 3, 3),
    )
    false
  } catch {
    AlignMauveError(_) => true
  }
  assert_true(rejected)
}

///|
test "Bio.Align.mauve rejects unknown reconstruction identifiers" {
  let rejected = try {
    ignore(
      @src.align_mauve_reconstruct_sequence(
        align_mauve_test_document(),
        "missing",
      ),
    )
    false
  } catch {
    AlignMauveError(_) => true
  }
  assert_true(rejected)
}

///|
test "Bio.Align.mauve rejects invalid reconstruction fill characters" {
  let rejected = try {
    ignore(
      @src.align_mauve_reconstruct_sequence(
        align_mauve_test_document(),
        "0",
        fill="-",
      ),
    )
    false
  } catch {
    AlignMauveError(_) => true
  }
  assert_true(rejected)
}

///|
test "Bio.Align.mauve detects conflicting overlapping reconstruction rows" {
  let rejected = try {
    ignore(
      @src.align_mauve_reconstruct_sequence(
        align_mauve_test_parse(align_mauve_test_conflict_text()),
        "0",
      ),
    )
    false
  } catch {
    AlignMauveError(_) => true
  }
  assert_true(rejected)
}

///|
test "Bio.Align.mauve validates public source constructors" {
  let rejected = try {
    ignore(@src.AlignMauveSource::create(0, "a.fa"))
    false
  } catch {
    AlignMauveError(_) => true
  }
  assert_true(rejected)
}

///|
test "Bio.Align.mauve validates public row constructors" {
  let source = @src.AlignMauveSource::create(1, "a.fa") catch {
    AlignMauveError(message) => abort("valid source failed: " + message)
  }
  let rejected = try {
    ignore(@src.AlignMauveRow::create(source, 0, 2, "+", "a.fa", "A"))
    false
  } catch {
    AlignMauveError(_) => true
  }
  assert_true(rejected)
}
