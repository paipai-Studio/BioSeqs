// Black-box tests for Biopython Bio.Align.nexus-compatible support.

///|
fn align_nexus_test_parse(text : String) -> @src.AlignNexusAlignment {
  @src.align_nexus_parse(text) catch {
    AlignNexusError(message) => abort("valid NEXUS input failed: " + message)
  }
}

///|
fn align_nexus_test_raises(text : String) -> Bool {
  try {
    ignore(@src.align_nexus_parse(text))
    false
  } catch {
    AlignNexusError(_) => true
  }
}

///|
fn align_nexus_test_write(
  alignment : @src.AlignNexusAlignment,
  interleave : Bool?,
  block_width : Int,
) -> String {
  @src.align_nexus_write(alignment, interleave~, block_width~) catch {
    AlignNexusError(message) => abort("NEXUS writing failed: " + message)
  }
}

///|
fn align_nexus_test_write_raises(
  alignment : @src.AlignNexusAlignment,
  block_width : Int,
) -> Bool {
  try {
    ignore(@src.align_nexus_write(alignment, block_width~))
    false
  } catch {
    AlignNexusError(_) => true
  }
}

///|
fn align_nexus_test_from_raises(
  ids : Array[String],
  rows : Array[String],
) -> Bool {
  try {
    ignore(@src.align_nexus_from_aligned(ids, rows, @src.AlignNexusDna))
    false
  } catch {
    AlignNexusError(_) => true
  }
}

///|
fn align_nexus_test_sample() -> @src.AlignNexusAlignment {
  align_nexus_test_parse(@src.align_nexus_example_text())
}

///|
fn align_nexus_test_codon_text() -> String {
  "#NEXUS\n" +
  "[MacClade fixture]\n" +
  "BEGIN DATA;\n" +
  "DIMENSIONS NTAX=2 NCHAR=22;\n" +
  "FORMAT DATATYPE=DNA MISSING=? GAP=-;\n" +
  "MATRIX\n" +
  "[                  10        20]\n" +
  "Aegotheles         AAAAAGGCATTGTGGTGGGAAT [22]\n" +
  "Aerodramus         ?????????TTGTGGTGGGAAT [13]\n" +
  ";\n" +
  "END;\n" +
  "BEGIN CODONS;\n" +
  "CODONPOSSET * CodonPositions = N: 1-10, 1: 11-22\\3;\n" +
  "END;\n"
}

///|
fn align_nexus_test_official_text() -> String {
  "#NEXUS\n" +
  "[ NEXUS test file. Weird, but NEXUS standard compatible.]\n" +
  "begin data;\n" +
  "dimensions ntax=9 nchar=48;\n" +
  "format datatype=dna gap=- missing=? interleave=yes;\n" +
  "options gapmode=newstate;\n" +
  "taxlabels t1\n" +
  "'t2 the name'\n" +
  "'isn''that [a] strange name?' 'one should be punished, for (that)!'\n" +
  "t5 t6 t7 t8 t9;\n" +
  "charlabels 1 'a', 2 'b', 48 final;\n" +
  "matrix\n" +
  "[interleaved matrix]\n" +
  "t1                                      A-C-G-Tc-gtgtgtgctct-t-t\n" +
  "'t2 the name'                           A-C-GcTc-gtg-----tct-t-t\n" +
  "'isn''that [a] strange name?'           A-CcGcTc-gtgtgtgct------\n" +
  "'one should be punished, for (that)!'   A-C-G-T-----------------\n" +
  "t5                                      A-C?G-T?--acgt??-???-???\n" +
  "t6                                      AcCaGtTc-?--aaaaaaaa-a-a\n" +
  "t7                                      A?C-GgTg-ggggggggggg-g-g\n" +
  "t8                                      AtCtGtTt-ttttttttttt-?-?\n" +
  "[last taxon of the matrix]\n" +
  "t9                                      cccccccc-cccccccccccNc-c\n" +
  "t1                                      ------ac-gtgtgtgctct-t-t\n" +
  "t2_the_name                             ----acac-gtg-----tct-t-t\n" +
  "'isn''that [a] strange name?'           --acacac-gtgtgtgct------\n" +
  "'one should be punished, for (that)!'   ------------------------\n" +
  "t5                                      --??---?--acgt??-???-???\n" +
  "t6                                      acgactac-?--aaaaaaaa-a-a\n" +
  "t7                                      ??--?ggg-ggggggggggg-g-g\n" +
  "t8                                      tttttttt-ttttttttttt-?-?\n" +
  "[last taxon of the matrix]\n" +
  "t9                                      cccccccc-cccccccccccNc-c\n" +
  ";\n" +
  "end;\n" +
  "[ this is a [nested ] comment]\n" +
  "begin sets;\n" +
  "charset c1=1-8;\n" +
  "end;\n" +
  "begin spam;\n" +
  "this is an unknown block;\n" +
  "end;\n"
}

///|
fn align_nexus_test_official() -> @src.AlignNexusAlignment {
  align_nexus_test_parse(align_nexus_test_official_text())
}

///|
test "Bio.Align.nexus parses interleaved example metadata" {
  let alignment = align_nexus_test_sample()
  assert_eq(alignment.metadata.block_name, "data")
  assert_eq(alignment.metadata.declared_taxa, 3)
  assert_eq(alignment.metadata.declared_characters, 12)
  assert_eq(alignment.metadata.data_type, @src.AlignNexusDna)
  assert_true(alignment.metadata.interleaved)
  assert_eq(alignment.metadata.match_character, Some("."))
}

///|
test "Bio.Align.nexus resolves MATCHCHAR across interleaved blocks" {
  let alignment = align_nexus_test_sample()
  assert_eq(alignment.sequences[0].aligned_sequence, "ACGT--ACGTAC")
  assert_eq(alignment.sequences[1].aligned_sequence, "ACGT??ACG-AC")
  assert_eq(alignment.sequences[2].aligned_sequence, "AC-T--A-GTAC")
}

///|
test "Bio.Align.nexus preserves quoted and escaped identifiers" {
  let alignment = align_nexus_test_sample()
  assert_eq(alignment.sequences[1].id, "query one")
  assert_eq(alignment.sequences[2].id, "isn't_query")
}

///|
test "Bio.Align.nexus treats underscores and spaces as interleave equivalents" {
  let alignment = align_nexus_test_sample()
  assert_eq(alignment.sequences[1].id, "query one")
  assert_eq(alignment.sequences[1].sequence, "ACGT??ACGAC")
}

///|
test "Bio.Align.nexus reports sample shape" {
  let alignment = align_nexus_test_sample()
  assert_eq(alignment.num_sequences(), 3)
  assert_eq(alignment.alignment_length(), 12)
  assert_eq(alignment.source_alignment_length(), 12)
}

///|
test "Bio.Align.nexus reports sample summary" {
  assert_eq(
    align_nexus_test_sample().summary(),
    "AlignNexusAlignment(type=dna, sequences=3, columns=12, " +
    "source_columns=12, removed_all_gap=0)",
  )
}

///|
test "Bio.Align.nexus exposes molecule type annotations" {
  assert_eq(@src.AlignNexusDna.molecule_type(), "DNA")
  assert_eq(@src.AlignNexusRna.molecule_type(), "RNA")
  assert_eq(@src.AlignNexusProtein.molecule_type(), "protein")
  assert_eq(@src.AlignNexusStandard.molecule_type(), "")
}

///|
test "Bio.Align.nexus exposes canonical datatype codes" {
  assert_eq(@src.AlignNexusDna.code(), "dna")
  assert_eq(@src.AlignNexusRna.code(), "rna")
  assert_eq(@src.AlignNexusProtein.code(), "protein")
  assert_eq(@src.AlignNexusStandard.code(), "standard")
}

///|
test "Bio.Align.nexus finds duplicate-safe row identifiers" {
  let alignment = align_nexus_test_sample()
  assert_eq(alignment.find_sequence("reference"), Some(0))
  assert_eq(alignment.find_sequence("query one"), Some(1))
  assert_eq(alignment.find_sequence("missing"), None)
}

///|
test "Bio.Align.nexus returns printed columns" {
  let alignment = align_nexus_test_sample()
  assert_eq(alignment.column(0), Some("AAA"))
  assert_eq(alignment.column(2), Some("GG-"))
  assert_eq(alignment.column(4), Some("-?-"))
  assert_eq(alignment.column(-1), None)
  assert_eq(alignment.column(12), None)
}

///|
test "Bio.Align.nexus maps sequence positions to columns" {
  let alignment = align_nexus_test_sample()
  assert_eq(alignment.sequence_position_to_column(1, 0), Some(0))
  assert_eq(alignment.sequence_position_to_column(1, 7), Some(7))
  assert_eq(alignment.sequence_position_to_column(1, 8), Some(8))
  assert_eq(alignment.sequence_position_to_column(1, 11), None)
}

///|
test "Bio.Align.nexus maps columns to sequence positions" {
  let alignment = align_nexus_test_sample()
  assert_eq(alignment.column_to_sequence_position(1, 0), Some(0))
  assert_eq(alignment.column_to_sequence_position(1, 8), Some(8))
  assert_eq(alignment.column_to_sequence_position(1, 9), None)
  assert_eq(alignment.column_to_sequence_position(1, 12), None)
}

///|
test "Bio.Align.nexus maps residues between rows" {
  let alignment = align_nexus_test_sample()
  assert_eq(alignment.map_position(0, 1, 0), Some(0))
  assert_eq(alignment.map_position(0, 1, 5), Some(7))
  assert_eq(alignment.map_position(0, 1, 7), None)
  assert_eq(alignment.map_position(0, 9, 0), None)
}

///|
test "Bio.Align.nexus returns aligned coordinate pairs" {
  let pairs = align_nexus_test_sample().aligned_pairs(0, 1) catch {
    AlignNexusError(message) => abort(message)
  }
  assert_eq(pairs.length(), 12)
  assert_eq(pairs[0], (Some(0), Some(0)))
  assert_eq(pairs[4], (None, Some(4)))
  assert_eq(pairs[9], (Some(7), None))
  assert_eq(pairs[11], (Some(9), Some(10)))
}

///|
test "Bio.Align.nexus builds a compact coordinate path" {
  let path = align_nexus_test_sample().coordinate_path()
  assert_eq(path.length(), 3)
  assert_eq(path[0], [0, 2, 3, 4, 4, 5, 6, 7, 8, 10])
  assert_eq(path[1], [0, 2, 3, 4, 6, 7, 8, 9, 9, 11])
  assert_eq(path[2], [0, 2, 2, 3, 3, 4, 4, 5, 6, 8])
}

///|
test "Bio.Align.nexus computes pair counts" {
  let counts = align_nexus_test_sample().pair_counts(0, 1) catch {
    AlignNexusError(message) => abort(message)
  }
  assert_eq(counts.pairs, 1)
  assert_eq(counts.aligned, 9)
  assert_eq(counts.identities, 9)
  assert_eq(counts.mismatches, 0)
  assert_eq(counts.gap_columns, 3)
  assert_eq(counts.double_gap_columns, 0)
  assert_eq(counts.gap_opens, 2)
  assert_eq(counts.identity(), 1.0)
}

///|
test "Bio.Align.nexus aggregates all unordered row pairs" {
  let counts = align_nexus_test_sample().counts()
  assert_eq(counts.pairs, 3)
  assert_eq(counts.aligned, 24)
  assert_eq(counts.identities, 24)
  assert_eq(counts.mismatches, 0)
  assert_eq(counts.gap_columns, 10)
  assert_eq(counts.double_gap_columns, 2)
}

///|
test "Bio.Align.nexus calculates consensus" {
  let consensus = align_nexus_test_sample().consensus() catch {
    AlignNexusError(message) => abort(message)
  }
  assert_eq(consensus, "ACGT??ACGTAC")
}

///|
test "Bio.Align.nexus applies consensus threshold" {
  let consensus = align_nexus_test_sample().consensus(minimum_fraction=0.75) catch {
    AlignNexusError(message) => abort(message)
  }
  assert_eq(consensus, "ACGT??ACGTAC")
}

///|
test "Bio.Align.nexus computes occupancy" {
  let occupancy = align_nexus_test_sample().occupancy()
  assert_eq(occupancy.length(), 12)
  assert_eq(occupancy[0], 1.0)
  assert_eq(occupancy[4], 1.0 / 3.0)
  assert_eq(occupancy[6], 1.0)
  assert_eq(occupancy[9], 2.0 / 3.0)
}

///|
test "Bio.Align.nexus parses official codon fixture" {
  let alignment = align_nexus_test_parse(align_nexus_test_codon_text())
  assert_eq(alignment.num_sequences(), 2)
  assert_eq(alignment.alignment_length(), 22)
  assert_eq(alignment.sequences[0].id, "Aegotheles")
  assert_eq(alignment.sequences[1].id, "Aerodramus")
}

///|
test "Bio.Align.nexus preserves codon fixture sequences" {
  let alignment = align_nexus_test_parse(align_nexus_test_codon_text())
  assert_eq(alignment.sequences[0].sequence, "AAAAAGGCATTGTGGTGGGAAT")
  assert_eq(alignment.sequences[1].sequence, "?????????TTGTGGTGGGAAT")
  assert_eq(alignment.coordinate_path(), [[0, 22], [0, 22]])
}

///|
test "Bio.Align.nexus matches codon fixture counts" {
  let counts = align_nexus_test_parse(align_nexus_test_codon_text()).counts()
  assert_eq(counts.pairs, 1)
  assert_eq(counts.aligned, 22)
  assert_eq(counts.identities, 13)
  assert_eq(counts.mismatches, 9)
  assert_eq(counts.gap_columns, 0)
}

///|
test "Bio.Align.nexus parses official nine-taxon fixture" {
  let alignment = align_nexus_test_official()
  assert_eq(alignment.num_sequences(), 9)
  assert_eq(alignment.metadata.declared_characters, 48)
  assert_eq(alignment.alignment_length(), 46)
  assert_eq(alignment.metadata.removed_all_gap_columns, 2)
}

///|
test "Bio.Align.nexus preserves official quoted identifiers" {
  let alignment = align_nexus_test_official()
  assert_eq(alignment.sequences[0].id, "t1")
  assert_eq(alignment.sequences[1].id, "t2 the name")
  assert_eq(alignment.sequences[2].id, "isn'that [a] strange name?")
  assert_eq(alignment.sequences[3].id, "one should be punished, for (that)!")
  assert_eq(alignment.sequences[8].id, "t9")
}

///|
test "Bio.Align.nexus preserves official ungapped sequences" {
  let alignment = align_nexus_test_official()
  assert_eq(
    alignment.sequences[0].sequence,
    "ACGTcgtgtgtgctctttacgtgtgtgctcttt",
  )
  assert_eq(alignment.sequences[1].sequence, "ACGcTcgtgtctttacacgtgtcttt")
  assert_eq(alignment.sequences[2].sequence, "ACcGcTcgtgtgtgctacacacgtgtgtgct")
  assert_eq(alignment.sequences[3].sequence, "ACGT")
}

///|
test "Bio.Align.nexus preserves official missing data" {
  let alignment = align_nexus_test_official()
  assert_eq(
    alignment.sequences[4].sequence,
    "AC?GT?acgt???????????acgt????????",
  )
  assert_eq(
    alignment.sequences[8].sequence,
    "cccccccccccccccccccNcccccccccccccccccccccNcc",
  )
}

///|
test "Bio.Align.nexus removes official all-gap columns" {
  let alignment = align_nexus_test_official()
  assert_eq(
    alignment.sequences[0].aligned_sequence,
    "A-C-G-Tcgtgtgtgctct-t-t------acgtgtgtgctct-t-t",
  )
  assert_eq(
    alignment.sequences[3].aligned_sequence,
    "A-C-G-T---------------------------------------",
  )
}

///|
test "Bio.Align.nexus matches official first coordinate path" {
  let path = align_nexus_test_official().coordinate_path()
  assert_eq(path[0], [
    0, 1, 1, 2, 2, 3, 3, 4, 5, 6, 8, 12, 13, 14, 16, 16, 17, 17, 18, 18, 18, 18,
    19, 20, 21, 23, 27, 28, 29, 31, 31, 32, 32, 33,
  ])
}

///|
test "Bio.Align.nexus matches official sparse coordinate path" {
  let path = align_nexus_test_official().coordinate_path()
  assert_eq(path[3], [
    0, 1, 1, 2, 2, 3, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
    4, 4, 4, 4, 4, 4, 4, 4,
  ])
}

///|
test "Bio.Align.nexus matches official final coordinate path" {
  let path = align_nexus_test_official().coordinate_path()
  assert_eq(path[8], [
    0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 11, 15, 16, 17, 19, 20, 21, 21, 22, 24, 26, 28,
    29, 30, 31, 33, 37, 38, 39, 41, 42, 43, 43, 44,
  ])
}

///|
test "Bio.Align.nexus matches official aggregate counts" {
  let counts = align_nexus_test_official().counts()
  assert_eq(counts.pairs, 36)
  assert_eq(counts.aligned, 862)
  assert_eq(counts.identities, 256)
  assert_eq(counts.mismatches, 606)
  assert_eq(counts.gap_columns, 596)
}

///|
test "Bio.Align.nexus writes official canonical dimensions" {
  let output = align_nexus_test_write(
    align_nexus_test_official(),
    Some(false),
    70,
  )
  assert_true(output.has_prefix("#NEXUS\nbegin data;\n"))
  assert_true(output.contains("dimensions ntax=9 nchar=46;"))
  assert_true(output.contains("format datatype=dna missing=? gap=-;"))
}

///|
test "Bio.Align.nexus writer quotes official identifiers" {
  let output = align_nexus_test_write(
    align_nexus_test_official(),
    Some(false),
    70,
  )
  assert_true(output.contains("'t2 the name'"))
  assert_true(output.contains("'isn''that [a] strange name?'"))
  assert_true(output.contains("'one should be punished, for (that)!'"))
}

///|
test "Bio.Align.nexus sequential writer round-trips" {
  let alignment = align_nexus_test_official()
  let output = align_nexus_test_write(alignment, Some(false), 70)
  let reparsed = align_nexus_test_parse(output)
  assert_eq(reparsed.metadata.declared_characters, 46)
  assert_eq(reparsed.metadata.removed_all_gap_columns, 0)
  assert_eq(reparsed.sequences, alignment.sequences)
  assert_eq(reparsed.coordinate_path(), alignment.coordinate_path())
}

///|
test "Bio.Align.nexus interleaved writer round-trips" {
  let alignment = align_nexus_test_sample()
  let output = align_nexus_test_write(alignment, Some(true), 4)
  assert_true(output.contains(" interleave;"))
  let reparsed = align_nexus_test_parse(output)
  assert_eq(reparsed.sequences, alignment.sequences)
  assert_eq(reparsed.coordinate_path(), alignment.coordinate_path())
}

///|
test "Bio.Align.nexus sequential writer omits interleave" {
  let output = align_nexus_test_write(align_nexus_test_sample(), Some(false), 4)
  assert_false(output.contains(" interleave;"))
}

///|
test "Bio.Align.nexus writer uses canonical terminator" {
  let output = align_nexus_test_write(
    align_nexus_test_sample(),
    Some(false),
    70,
  )
  assert_true(output.has_suffix(";\nend;\n"))
}

///|
test "Bio.Align.nexus writer auto-interleaves above 1000 columns" {
  let short_builder = StringBuilder::new(size_hint=1000)
  let long_builder = StringBuilder::new(size_hint=1001)
  for _index = 0; _index < 1000; _index = _index + 1 {
    short_builder.write_char('A')
    long_builder.write_char('A')
  }
  long_builder.write_char('A')
  let short_alignment = @src.align_nexus_from_aligned(
    ["alpha"],
    [short_builder.to_string()],
    @src.AlignNexusDna,
  ) catch {
    AlignNexusError(message) => abort(message)
  }
  let long_alignment = @src.align_nexus_from_aligned(
    ["alpha"],
    [long_builder.to_string()],
    @src.AlignNexusDna,
  ) catch {
    AlignNexusError(message) => abort(message)
  }
  let short_output = @src.align_nexus_write(short_alignment) catch {
    AlignNexusError(message) => abort(message)
  }
  let long_output = @src.align_nexus_write(long_alignment) catch {
    AlignNexusError(message) => abort(message)
  }
  assert_false(short_output.contains(" interleave;"))
  assert_true(long_output.contains(" interleave;"))
}

///|
test "Bio.Align.nexus parses CHARACTERS with TAXA association" {
  let text = "#NEXUS\n" +
    "begin taxa;\n" +
    "dimensions ntax=2;\n" +
    "taxlabels 'taxon one' taxon_two;\n" +
    "end;\n" +
    "begin characters;\n" +
    "dimensions nchar=4;\n" +
    "format datatype=dna labels=no;\n" +
    "matrix\nACGT\nA-GT\n;\n" +
    "end;\n"
  let alignment = align_nexus_test_parse(text)
  assert_eq(alignment.metadata.block_name, "characters")
  assert_eq(alignment.sequences[0].id, "taxon one")
  assert_eq(alignment.sequences[1].id, "taxon_two")
}

///|
test "Bio.Align.nexus parses wrapped sequential rows" {
  let text = "#NEXUS\nbegin data;\n" +
    "dimensions ntax=2 nchar=8;\nformat datatype=dna;\nmatrix\n" +
    "alpha ACGT\nACGT\nbeta\nAC--\nACGT\n;\nend;\n"
  let alignment = align_nexus_test_parse(text)
  assert_eq(alignment.sequences[0].aligned_sequence, "ACGTACGT")
  assert_eq(alignment.sequences[1].aligned_sequence, "AC--ACGT")
}

///|
test "Bio.Align.nexus accepts CRLF input" {
  let text = "#NEXUS\r\nBEGIN DATA;\r\n" +
    "DIMENSIONS NTAX=1 NCHAR=4;\r\n" +
    "FORMAT DATATYPE=DNA;\r\nMATRIX\r\nalpha ACGT\r\n;\r\nEND;\r\n"
  let alignment = align_nexus_test_parse(text)
  assert_eq(alignment.sequences[0].sequence, "ACGT")
}

///|
test "Bio.Align.nexus preserves duplicate taxon names" {
  let text = "#NEXUS\nbegin data;\ndimensions ntax=2 nchar=4;\n" +
    "format datatype=dna;\nmatrix\nalpha ACGT\nalpha A-GT\n;\nend;\n"
  let alignment = align_nexus_test_parse(text)
  assert_eq(alignment.sequences[0].id, "alpha")
  assert_eq(alignment.sequences[1].id, "alpha")
  assert_eq(alignment.find_sequence("alpha"), Some(0))
}

///|
test "Bio.Align.nexus preserves lower-case residues" {
  let text = "#NEXUS\nbegin data;\ndimensions ntax=1 nchar=6;\n" +
    "format datatype=dna;\nmatrix\nalpha AcgTnN\n;\nend;\n"
  let alignment = align_nexus_test_parse(text)
  assert_eq(alignment.sequences[0].sequence, "AcgTnN")
}

///|
test "Bio.Align.nexus parses RNA data" {
  let text = "#NEXUS\nbegin data;\ndimensions ntax=2 nchar=5;\n" +
    "format datatype=rna;\nmatrix\nalpha ACGUN\nbeta A-GUN\n;\nend;\n"
  let alignment = align_nexus_test_parse(text)
  assert_eq(alignment.metadata.data_type, @src.AlignNexusRna)
  assert_eq(alignment.sequences[1].sequence, "AGUN")
}

///|
test "Bio.Align.nexus parses protein data" {
  let text = "#NEXUS\nbegin data;\ndimensions ntax=2 nchar=7;\n" +
    "format datatype=protein;\nmatrix\nalpha ACDE*XZ\nbeta A-DEBXZ\n;\nend;\n"
  let alignment = align_nexus_test_parse(text)
  assert_eq(alignment.metadata.data_type, @src.AlignNexusProtein)
  assert_eq(alignment.sequences[0].sequence, "ACDE*XZ")
}

///|
test "Bio.Align.nexus parses STANDARD symbols" {
  let text = "#NEXUS\nbegin characters;\ndimensions ntax=2 nchar=5;\n" +
    "format datatype=standard symbols=\"01\" missing=? gap=-;\n" +
    "matrix\nalpha 0101?\nbeta 0-11?\n;\nend;\n"
  let alignment = align_nexus_test_parse(text)
  assert_eq(alignment.metadata.data_type, @src.AlignNexusStandard)
  assert_eq(alignment.metadata.symbols, "01")
  assert_eq(alignment.sequences[1].sequence, "011?")
}

///|
test "Bio.Align.nexus STANDARD is case-insensitive by default" {
  let text = "#NEXUS\nbegin data;\ndimensions ntax=1 nchar=4;\n" +
    "format datatype=standard symbols=\"AB\";\nmatrix\nalpha AbBa\n;\nend;\n"
  let alignment = align_nexus_test_parse(text)
  assert_eq(alignment.sequences[0].sequence, "AbBa")
}

///|
test "Bio.Align.nexus supports custom gap and missing symbols" {
  let text = "#NEXUS\nbegin data;\ndimensions ntax=2 nchar=5;\n" +
    "format datatype=dna missing=X gap=~;\nmatrix\nalpha AC~TX\n" +
    "beta A~~TX\n;\nend;\n"
  let alignment = align_nexus_test_parse(text)
  assert_eq(alignment.metadata.missing_character, "X")
  assert_eq(alignment.metadata.gap_character, "~")
  assert_eq(alignment.metadata.removed_all_gap_columns, 1)
  assert_eq(alignment.sequences[0].aligned_sequence, "ACTX")
  let output = align_nexus_test_write(alignment, Some(false), 70)
  assert_true(output.contains("missing=X gap=~;"))
  assert_true(output.contains("ACTX"))
}

///|
test "Bio.Align.nexus removes constructor all-gap columns" {
  let alignment = @src.align_nexus_from_aligned(
    ["alpha", "beta"],
    ["A-C--G", "ATC--G"],
    @src.AlignNexusDna,
  ) catch {
    AlignNexusError(message) => abort(message)
  }
  assert_eq(alignment.metadata.declared_characters, 6)
  assert_eq(alignment.metadata.removed_all_gap_columns, 2)
  assert_eq(alignment.alignment_length(), 4)
  assert_eq(alignment.sequences[0].aligned_sequence, "A-CG")
}

///|
test "Bio.Align.nexus constructor supports all-gap row" {
  let alignment = @src.align_nexus_from_aligned(
    ["alpha", "empty"],
    ["ACGT", "----"],
    @src.AlignNexusDna,
  ) catch {
    AlignNexusError(message) => abort(message)
  }
  assert_eq(alignment.sequences[1].sequence, "")
  assert_eq(alignment.sequences[1].aligned_sequence, "----")
}

///|
test "Bio.Align.nexus row constructor normalizes custom gap" {
  let sequence = @src.AlignNexusSequence::create(
    "alpha",
    "AC~GT",
    @src.AlignNexusDna,
    gap_character="~",
  ) catch {
    AlignNexusError(message) => abort(message)
  }
  assert_eq(sequence.aligned_sequence, "AC-GT")
  assert_eq(sequence.sequence, "ACGT")
}

///|
test "Bio.Align.nexus metadata constructor preserves source facts" {
  let metadata = @src.AlignNexusMetadata::create(
    3,
    12,
    @src.AlignNexusDna,
    block_name="CHARACTERS",
    match_character=Some("."),
    interleaved=true,
    removed_all_gap_columns=2,
  ) catch {
    AlignNexusError(message) => abort(message)
  }
  assert_eq(metadata.block_name, "characters")
  assert_eq(metadata.match_character, Some("."))
  assert_eq(metadata.removed_all_gap_columns, 2)
}

///|
test "Bio.Align.nexus rejects empty input" {
  assert_true(align_nexus_test_raises(""))
}

///|
test "Bio.Align.nexus rejects whitespace-only input" {
  assert_true(align_nexus_test_raises(" \n\t"))
}

///|
test "Bio.Align.nexus rejects missing header" {
  assert_true(
    align_nexus_test_raises(
      "begin data; dimensions ntax=1 nchar=1; " +
      "format datatype=dna; matrix alpha A; end;",
    ),
  )
}

///|
test "Bio.Align.nexus rejects wrong-case header" {
  assert_true(
    align_nexus_test_raises(
      "#nexus\nbegin data; dimensions ntax=1 nchar=1; " +
      "format datatype=dna; matrix alpha A; end;",
    ),
  )
}

///|
test "Bio.Align.nexus rejects missing matrix" {
  assert_true(
    align_nexus_test_raises(
      "#NEXUS\nbegin data;\ndimensions ntax=1 nchar=1;\n" +
      "format datatype=dna;\nend;\n",
    ),
  )
}

///|
test "Bio.Align.nexus rejects missing semicolon" {
  assert_true(
    align_nexus_test_raises(
      "#NEXUS\nbegin data;\ndimensions ntax=1 nchar=1;\n" +
      "format datatype=dna;\nmatrix\nalpha A\n;\nend\n",
    ),
  )
}

///|
test "Bio.Align.nexus rejects nested blocks" {
  assert_true(
    align_nexus_test_raises("#NEXUS\nbegin data;\nbegin taxa;\nend;\nend;\n"),
  )
}

///|
test "Bio.Align.nexus rejects unmatched END" {
  assert_true(align_nexus_test_raises("#NEXUS\nend;\n"))
}

///|
test "Bio.Align.nexus rejects unterminated comment" {
  assert_true(align_nexus_test_raises("#NEXUS\n[broken\nbegin data;\nend;\n"))
}

///|
test "Bio.Align.nexus rejects unmatched comment close" {
  assert_true(align_nexus_test_raises("#NEXUS\n]\n"))
}

///|
test "Bio.Align.nexus rejects unterminated quote" {
  assert_true(
    align_nexus_test_raises(
      "#NEXUS\nbegin data;\ndimensions ntax=1 nchar=1;\n" +
      "format datatype=dna;\nmatrix\n'alpha A\n;\nend;\n",
    ),
  )
}

///|
test "Bio.Align.nexus rejects missing NTAX" {
  assert_true(
    align_nexus_test_raises(
      "#NEXUS\nbegin data;\ndimensions nchar=1;\n" +
      "format datatype=dna;\nmatrix\nalpha A\n;\nend;\n",
    ),
  )
}

///|
test "Bio.Align.nexus rejects missing NCHAR" {
  assert_true(
    align_nexus_test_raises(
      "#NEXUS\nbegin data;\ndimensions ntax=1;\n" +
      "format datatype=dna;\nmatrix\nalpha A\n;\nend;\n",
    ),
  )
}

///|
test "Bio.Align.nexus rejects zero dimensions" {
  assert_true(
    align_nexus_test_raises(
      "#NEXUS\nbegin data;\ndimensions ntax=0 nchar=1;\nend;\n",
    ),
  )
}

///|
test "Bio.Align.nexus rejects inconsistent TAXA dimensions" {
  let text = "#NEXUS\nbegin taxa;\ndimensions ntax=3;\n" +
    "taxlabels a b c;\nend;\nbegin characters;\n" +
    "dimensions ntax=2 nchar=1;\nformat datatype=dna;\n" +
    "matrix\na A\nb C\n;\nend;\n"
  assert_true(align_nexus_test_raises(text))
}

///|
test "Bio.Align.nexus rejects TAXLABELS count mismatch" {
  let text = "#NEXUS\nbegin taxa;\ndimensions ntax=2;\n" +
    "taxlabels alpha;\nend;\nbegin characters;\n" +
    "dimensions nchar=1;\nformat datatype=dna labels=no;\n" +
    "matrix\nA\nC\n;\nend;\n"
  assert_true(align_nexus_test_raises(text))
}

///|
test "Bio.Align.nexus rejects short sequential row" {
  let text = "#NEXUS\nbegin data;\ndimensions ntax=1 nchar=5;\n" +
    "format datatype=dna;\nmatrix\nalpha ACGT\n;\nend;\n"
  assert_true(align_nexus_test_raises(text))
}

///|
test "Bio.Align.nexus rejects long sequential row" {
  let text = "#NEXUS\nbegin data;\ndimensions ntax=1 nchar=3;\n" +
    "format datatype=dna;\nmatrix\nalpha ACGT\n;\nend;\n"
  assert_true(align_nexus_test_raises(text))
}

///|
test "Bio.Align.nexus rejects incomplete interleaved block" {
  let text = "#NEXUS\nbegin data;\ndimensions ntax=2 nchar=4;\n" +
    "format datatype=dna interleave;\nmatrix\n" +
    "alpha AC\nbeta AC\nalpha GT\n;\nend;\n"
  assert_true(align_nexus_test_raises(text))
}

///|
test "Bio.Align.nexus rejects unequal interleaved chunk widths" {
  let text = "#NEXUS\nbegin data;\ndimensions ntax=2 nchar=4;\n" +
    "format datatype=dna interleave;\nmatrix\n" +
    "alpha AC\nbeta A\nalpha GT\nbeta CGT\n;\nend;\n"
  assert_true(align_nexus_test_raises(text))
}

///|
test "Bio.Align.nexus rejects interleaved taxon reorder" {
  let text = "#NEXUS\nbegin data;\ndimensions ntax=2 nchar=4;\n" +
    "format datatype=dna interleave;\nmatrix\n" +
    "alpha AC\nbeta AC\nbeta GT\nalpha GT\n;\nend;\n"
  assert_true(align_nexus_test_raises(text))
}

///|
test "Bio.Align.nexus rejects DNA U" {
  let text = "#NEXUS\nbegin data;\ndimensions ntax=1 nchar=4;\n" +
    "format datatype=dna;\nmatrix\nalpha ACGU\n;\nend;\n"
  assert_true(align_nexus_test_raises(text))
}

///|
test "Bio.Align.nexus rejects RNA T" {
  let text = "#NEXUS\nbegin data;\ndimensions ntax=1 nchar=4;\n" +
    "format datatype=rna;\nmatrix\nalpha ACGT\n;\nend;\n"
  assert_true(align_nexus_test_raises(text))
}

///|
test "Bio.Align.nexus rejects illegal protein residue" {
  let text = "#NEXUS\nbegin data;\ndimensions ntax=1 nchar=4;\n" +
    "format datatype=protein;\nmatrix\nalpha ACOU\n;\nend;\n"
  assert_true(align_nexus_test_raises(text))
}

///|
test "Bio.Align.nexus rejects STANDARD without symbols" {
  let text = "#NEXUS\nbegin data;\ndimensions ntax=1 nchar=2;\n" +
    "format datatype=standard;\nmatrix\nalpha 01\n;\nend;\n"
  assert_true(align_nexus_test_raises(text))
}

///|
test "Bio.Align.nexus rejects illegal STANDARD state" {
  let text = "#NEXUS\nbegin data;\ndimensions ntax=1 nchar=2;\n" +
    "format datatype=standard symbols=\"01\";\nmatrix\nalpha 02\n;\nend;\n"
  assert_true(align_nexus_test_raises(text))
}

///|
test "Bio.Align.nexus respects case for STANDARD symbols" {
  let text = "#NEXUS\nbegin data;\ndimensions ntax=1 nchar=2;\n" +
    "format datatype=standard symbols=\"AB\" respectcase;\n" +
    "matrix\nalpha Ab\n;\nend;\n"
  assert_true(align_nexus_test_raises(text))
}

///|
test "Bio.Align.nexus rejects match character in reference" {
  let text = "#NEXUS\nbegin data;\ndimensions ntax=2 nchar=4;\n" +
    "format datatype=dna matchchar=.;\nmatrix\n" +
    "alpha A.GT\nbeta ....\n;\nend;\n"
  assert_true(align_nexus_test_raises(text))
}

///|
test "Bio.Align.nexus rejects all-gap matrix" {
  let text = "#NEXUS\nbegin data;\ndimensions ntax=2 nchar=4;\n" +
    "format datatype=dna;\nmatrix\nalpha ----\nbeta ----\n;\nend;\n"
  assert_true(align_nexus_test_raises(text))
}

///|
test "Bio.Align.nexus rejects duplicate FORMAT character roles" {
  let text = "#NEXUS\nbegin data;\ndimensions ntax=1 nchar=1;\n" +
    "format datatype=dna missing=- gap=-;\nmatrix\nalpha A\n;\nend;\n"
  assert_true(align_nexus_test_raises(text))
}

///|
test "Bio.Align.nexus rejects unsupported datatype" {
  let text = "#NEXUS\nbegin data;\ndimensions ntax=1 nchar=1;\n" +
    "format datatype=continuous;\nmatrix\nalpha A\n;\nend;\n"
  assert_true(align_nexus_test_raises(text))
}

///|
test "Bio.Align.nexus rejects transpose matrices" {
  let text = "#NEXUS\nbegin data;\ndimensions ntax=1 nchar=1;\n" +
    "format datatype=dna transpose;\nmatrix\nalpha A\n;\nend;\n"
  assert_true(align_nexus_test_raises(text))
}

///|
test "Bio.Align.nexus rejects tokenized matrices" {
  let text = "#NEXUS\nbegin data;\ndimensions ntax=1 nchar=1;\n" +
    "format datatype=dna tokens;\nmatrix\nalpha A\n;\nend;\n"
  assert_true(align_nexus_test_raises(text))
}

///|
test "Bio.Align.nexus rejects labels=no without TAXLABELS" {
  let text = "#NEXUS\nbegin data;\ndimensions ntax=1 nchar=1;\n" +
    "format datatype=dna labels=no;\nmatrix\nA\n;\nend;\n"
  assert_true(align_nexus_test_raises(text))
}

///|
test "Bio.Align.nexus rejects malformed DIMENSIONS" {
  let text = "#NEXUS\nbegin data;\ndimensions ntax 1 nchar=1;\nend;\n"
  assert_true(align_nexus_test_raises(text))
}

///|
test "Bio.Align.nexus rejects unknown FORMAT option" {
  let text = "#NEXUS\nbegin data;\ndimensions ntax=1 nchar=1;\n" +
    "format datatype=dna mystery=yes;\nmatrix\nalpha A\n;\nend;\n"
  assert_true(align_nexus_test_raises(text))
}

///|
test "Bio.Align.nexus rejects invalid INTERLEAVE value" {
  let text = "#NEXUS\nbegin data;\ndimensions ntax=1 nchar=1;\n" +
    "format datatype=dna interleave=maybe;\nmatrix\nalpha A\n;\nend;\n"
  assert_true(align_nexus_test_raises(text))
}

///|
test "Bio.Align.nexus rejects mismatched constructor arrays" {
  assert_true(align_nexus_test_from_raises(["alpha"], ["ACGT", "ACGT"]))
}

///|
test "Bio.Align.nexus rejects empty constructor arrays" {
  assert_true(align_nexus_test_from_raises([], []))
}

///|
test "Bio.Align.nexus rejects unequal constructor widths" {
  assert_true(align_nexus_test_from_raises(["alpha", "beta"], ["ACGT", "ACG"]))
}

///|
test "Bio.Align.nexus rejects all-gap constructor alignment" {
  assert_true(align_nexus_test_from_raises(["alpha", "beta"], ["----", "----"]))
}

///|
test "Bio.Align.nexus rejects empty row identifier" {
  let failed = try {
    ignore(@src.AlignNexusSequence::create("", "ACGT", @src.AlignNexusDna))
    false
  } catch {
    AlignNexusError(_) => true
  }
  assert_true(failed)
}

///|
test "Bio.Align.nexus rejects invalid metadata block" {
  let failed = try {
    ignore(
      @src.AlignNexusMetadata::create(
        1,
        4,
        @src.AlignNexusDna,
        block_name="trees",
      ),
    )
    false
  } catch {
    AlignNexusError(_) => true
  }
  assert_true(failed)
}

///|
test "Bio.Align.nexus rejects invalid metadata removed count" {
  let failed = try {
    ignore(
      @src.AlignNexusMetadata::create(
        1,
        4,
        @src.AlignNexusDna,
        removed_all_gap_columns=4,
      ),
    )
    false
  } catch {
    AlignNexusError(_) => true
  }
  assert_true(failed)
}

///|
test "Bio.Align.nexus rejects invalid writer block width" {
  assert_true(align_nexus_test_write_raises(align_nexus_test_sample(), 0))
}

///|
test "Bio.Align.nexus rejects invalid pair row" {
  let failed = try {
    ignore(align_nexus_test_sample().pair_counts(0, 9))
    false
  } catch {
    AlignNexusError(_) => true
  }
  assert_true(failed)
}

///|
test "Bio.Align.nexus rejects invalid consensus threshold" {
  let failed = try {
    ignore(align_nexus_test_sample().consensus(minimum_fraction=1.1))
    false
  } catch {
    AlignNexusError(_) => true
  }
  assert_true(failed)
}
