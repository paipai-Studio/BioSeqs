// Black-box tests for Biopython Bio.Align.stockholm-compatible support.

///|
fn align_stockholm_test_parse(text : String) -> @src.AlignStockholmAlignment {
  @src.align_stockholm_parse(text) catch {
    AlignStockholmError(message) =>
      abort("valid Stockholm input failed: " + message)
  }
}

///|
fn align_stockholm_test_parse_all(
  text : String,
) -> Array[@src.AlignStockholmAlignment] {
  @src.align_stockholm_parse_all(text) catch {
    AlignStockholmError(message) =>
      abort("valid Stockholm input failed: " + message)
  }
}

///|
fn align_stockholm_test_raises(text : String) -> Bool {
  try {
    ignore(@src.align_stockholm_parse(text))
    false
  } catch {
    AlignStockholmError(_) => true
  }
}

///|
fn align_stockholm_test_parse_all_raises(text : String) -> Bool {
  try {
    ignore(@src.align_stockholm_parse_all(text))
    false
  } catch {
    AlignStockholmError(_) => true
  }
}

///|
fn align_stockholm_test_write(
  alignment : @src.AlignStockholmAlignment,
) -> String {
  @src.align_stockholm_write(alignment) catch {
    AlignStockholmError(message) =>
      abort("Stockholm writing failed: " + message)
  }
}

///|
fn align_stockholm_test_write_raises(
  alignment : @src.AlignStockholmAlignment,
) -> Bool {
  try {
    ignore(@src.align_stockholm_write(alignment))
    false
  } catch {
    AlignStockholmError(_) => true
  }
}

///|
fn align_stockholm_test_from(
  ids : Array[String],
  rows : Array[String],
) -> @src.AlignStockholmAlignment {
  @src.align_stockholm_from_aligned(ids, rows) catch {
    AlignStockholmError(message) =>
      abort("Stockholm construction failed: " + message)
  }
}

///|
fn align_stockholm_test_from_raises(
  ids : Array[String],
  rows : Array[String],
) -> Bool {
  try {
    ignore(@src.align_stockholm_from_aligned(ids, rows))
    false
  } catch {
    AlignStockholmError(_) => true
  }
}

///|
fn align_stockholm_test_sequence(
  id : String,
  row : String,
) -> @src.AlignStockholmSequence {
  @src.AlignStockholmSequence::create(id, row) catch {
    AlignStockholmError(message) =>
      abort("Stockholm sequence construction failed: " + message)
  }
}

///|
fn align_stockholm_test_sample() -> @src.AlignStockholmAlignment {
  align_stockholm_test_parse(@src.align_stockholm_example_text())
}

///|
fn align_stockholm_test_custom_text() -> String {
  "# STOCKHOLM 1.0\n" +
  "# ordinary comment\n" +
  "#=GF AU Alice;\n" +
  "#=GF AU Bob;\n" +
  "#=GF WK Long/\n" +
  "#=GF WK Name\n" +
  "#=GF SM first command\n" +
  "#=GF SM second command\n" +
  "#=GF XX ignored alignment annotation\n" +
  "#=GF RC reference note\n" +
  "#=GF RN [2]\n" +
  "#=GF RM 12345\n" +
  "#=GF RT title part one\n" +
  "#=GF RT title part two\n" +
  "#=GF RA Example Author;\n" +
  "#=GF RL Example Journal.\n" +
  "#=GF DR DB; ONE;\n" +
  "#=GF DC database comment\n" +
  "#=GF NE PF00001\n" +
  "#=GF NL 1-4\n" +
  "#=GF SQ 2\n" +
  "#=GS alpha AC A1\n" +
  "#=GS alpha OS organism one\n" +
  "#=GS alpha ZZ custom sequence value\n" +
  "#=GS alpha DE Alpha description\n" +
  "#=GS alpha DR DBX; 1;\n" +
  "#=GS alpha DR DBX; 2;\n" +
  "alpha A.-C--G\n" +
  "#=GR alpha SS x..y..z\n" +
  "beta AT-C--G\n" +
  "#=GR beta ZZ wx.y..z\n" +
  "#=GC RF 1234567\n" +
  "#=GC ZZ abcdefg\n" +
  "//\n"
}

///|
fn align_stockholm_test_custom() -> @src.AlignStockholmAlignment {
  align_stockholm_test_parse(align_stockholm_test_custom_text())
}

///|
test "Bio.Align.stockholm parses official HAT shape" {
  let alignment = align_stockholm_test_sample()
  assert_eq(alignment.version, "1.0")
  assert_eq(alignment.num_sequences(), 3)
  assert_eq(alignment.alignment_length(), 33)
  assert_eq(alignment.source_alignment_length(), 33)
  assert_eq(alignment.removed_all_gap_columns, 0)
}

///|
test "Bio.Align.stockholm reports official HAT summary" {
  assert_eq(
    align_stockholm_test_sample().summary(),
    "AlignStockholmAlignment(version=1.0, sequences=3, columns=33, " +
    "source_columns=33, removed_all_gap=0)",
  )
}

///|
test "Bio.Align.stockholm preserves official HAT identifiers" {
  let alignment = align_stockholm_test_sample()
  assert_eq(alignment.sequences[0].id, "CRN_DROME/191-222")
  assert_eq(alignment.sequences[1].id, "CLF1_SCHPO/185-216")
  assert_eq(alignment.sequences[2].id, "O16376_CAEEL/201-233")
}

///|
test "Bio.Align.stockholm normalizes official insertion gaps" {
  let alignment = align_stockholm_test_sample()
  assert_eq(
    alignment.sequences[0].aligned_sequence,
    "KEIDRAREIYERFVYVH-PDVKNWIKFARFEES",
  )
  assert_eq(
    alignment.sequences[1].aligned_sequence,
    "HENERARGIYERFVVVH-PEVTNWLRWARFEEE",
  )
  assert_eq(
    alignment.sequences[2].aligned_sequence,
    "KEIDRARSVYQRFLHVHGINVQNWIKYAKFEER",
  )
}

///|
test "Bio.Align.stockholm stores official ungapped sequences" {
  let alignment = align_stockholm_test_sample()
  assert_eq(alignment.sequences[0].sequence, "KEIDRAREIYERFVYVHPDVKNWIKFARFEES")
  assert_eq(alignment.sequences[1].sequence, "HENERARGIYERFVVVHPEVTNWLRWARFEEE")
  assert_eq(
    alignment.sequences[2].sequence,
    "KEIDRARSVYQRFLHVHGINVQNWIKYAKFEER",
  )
}

///|
test "Bio.Align.stockholm records official insertion operations" {
  assert_eq(
    align_stockholm_test_sample().operations,
    "MMMMMMMMMMMMMMMMMIMMMMMMMMMMMMMMM",
  )
}

///|
test "Bio.Align.stockholm maps official identifier annotation" {
  let alignment = align_stockholm_test_sample()
  assert_eq(alignment.annotation("ID"), Some("HAT"))
  assert_eq(alignment.annotation("identifier"), Some("HAT"))
}

///|
test "Bio.Align.stockholm maps official accession annotation" {
  let alignment = align_stockholm_test_sample()
  assert_eq(alignment.annotation("AC"), Some("PF02184.18"))
  assert_eq(alignment.annotation("accession"), Some("PF02184.18"))
}

///|
test "Bio.Align.stockholm maps official definition annotation" {
  assert_eq(
    align_stockholm_test_sample().annotation("definition"),
    Some("HAT (Half-A-TPR) repeat"),
  )
}

///|
test "Bio.Align.stockholm preserves repeated author values" {
  assert_eq(align_stockholm_test_sample().annotation_values("author"), [
    "SMART;",
  ])
}

///|
test "Bio.Align.stockholm maps official seed source" {
  assert_eq(
    align_stockholm_test_sample().annotation("source of seed"),
    Some("Alignment kindly provided by SMART"),
  )
}

///|
test "Bio.Align.stockholm maps official cutoff annotations" {
  let alignment = align_stockholm_test_sample()
  assert_eq(alignment.annotation("gathering method"), Some("21.00 21.00;"))
  assert_eq(alignment.annotation("trusted cutoff"), Some("21.00 21.00;"))
  assert_eq(alignment.annotation("noise cutoff"), Some("20.90 20.90;"))
}

///|
test "Bio.Align.stockholm maps official build method" {
  assert_eq(
    align_stockholm_test_sample().annotation("build method"),
    Some("hmmbuild HMM.ann SEED.ann"),
  )
}

///|
test "Bio.Align.stockholm maps official search method" {
  assert_eq(
    align_stockholm_test_sample().annotation("search method"),
    Some("hmmsearch -Z 57096847 -E 1000 --cpu 4 HMM pfamseq"),
  )
}

///|
test "Bio.Align.stockholm maps official type and clan" {
  let alignment = align_stockholm_test_sample()
  assert_eq(alignment.annotation("type"), Some("Repeat"))
  assert_eq(alignment.annotation("clan"), Some("CL0020"))
}

///|
test "Bio.Align.stockholm joins official comments" {
  assert_eq(
    align_stockholm_test_sample().annotation("comment"),
    Some(
      "The HAT (Half A TPR) repeat is found in several RNA processing " +
      "proteins [1].",
    ),
  )
}

///|
test "Bio.Align.stockholm parses official reference number and medline" {
  let references = align_stockholm_test_sample().references
  assert_eq(references.length(), 1)
  assert_eq(references[0].number, 1)
  assert_eq(references[0].medline, "9478129")
  assert_eq(references[0].comment, "")
}

///|
test "Bio.Align.stockholm parses official reference citation" {
  let reference = align_stockholm_test_sample().references[0]
  assert_eq(
    reference.title,
    "The HAT helix, a repetitive motif implicated in RNA processing.",
  )
  assert_eq(reference.author, "Preker PJ, Keller W;")
  assert_eq(reference.location, "Trends Biochem Sci 1998;23:15-16.")
}

///|
test "Bio.Align.stockholm parses official database references" {
  let references = align_stockholm_test_sample().database_references
  assert_eq(references.length(), 3)
  assert_eq(references[0].reference, "INTERPRO; IPR003107;")
  assert_eq(references[1].reference, "SMART; HAT;")
  assert_eq(references[2].reference, "SO; 0001068; polypeptide_repeat;")
}

///|
test "Bio.Align.stockholm parses official sequence accessions" {
  let sequences = align_stockholm_test_sample().sequences
  assert_eq(sequences[0].annotation("accession"), Some("P17886.2"))
  assert_eq(sequences[1].annotation("accession"), Some("P87312.1"))
  assert_eq(sequences[2].annotation("accession"), Some("O16376.2"))
}

///|
test "Bio.Align.stockholm parses official sequence database reference" {
  assert_eq(align_stockholm_test_sample().sequences[1].dbxrefs, [
    "PDB; 3JB9 R; 185-216;",
  ])
}

///|
test "Bio.Align.stockholm maps official GR secondary structure" {
  let sequence = align_stockholm_test_sample().sequences[1]
  assert_eq(
    sequence.letter_annotation("secondary structure"),
    Some("--HHHHHHHHHHHHHHS--HHHHHHHHHHHHH"),
  )
}

///|
test "Bio.Align.stockholm preserves official aligned GR value" {
  let annotation = align_stockholm_test_sample().sequences[1].letter_annotations[0]
  assert_eq(annotation.code, "SS")
  assert_eq(annotation.feature, "secondary structure")
  assert_eq(annotation.aligned_value, "--HHHHHHHHHHHHHHS.--HHHHHHHHHHHHH")
}

///|
test "Bio.Align.stockholm maps official GC secondary structure" {
  assert_eq(
    align_stockholm_test_sample().column_annotation(
      "consensus secondary structure",
    ),
    Some("--HHHHHHHHHHHHHHS.--HHHHHHHHHHHHH"),
  )
}

///|
test "Bio.Align.stockholm maps official GC consensus sequence" {
  assert_eq(
    align_stockholm_test_sample().column_annotation("consensus sequence"),
    Some("KEIDRARuIYERFVaVH.P-VpNWIKaARFEEc"),
  )
}

///|
test "Bio.Align.stockholm finds official rows" {
  let alignment = align_stockholm_test_sample()
  assert_eq(alignment.find_sequence("CRN_DROME/191-222"), Some(0))
  assert_eq(alignment.find_sequence("CLF1_SCHPO/185-216"), Some(1))
  assert_eq(alignment.find_sequence("missing"), None)
}

///|
test "Bio.Align.stockholm returns official printed columns" {
  let alignment = align_stockholm_test_sample()
  assert_eq(alignment.column(0), Some("KHK"))
  assert_eq(alignment.column(17), Some("--G"))
  assert_eq(alignment.column(32), Some("SER"))
  assert_eq(alignment.column(-1), None)
  assert_eq(alignment.column(33), None)
}

///|
test "Bio.Align.stockholm maps official sequence positions to columns" {
  let alignment = align_stockholm_test_sample()
  assert_eq(alignment.sequence_position_to_column(0, 0), Some(0))
  assert_eq(alignment.sequence_position_to_column(0, 16), Some(16))
  assert_eq(alignment.sequence_position_to_column(0, 17), Some(18))
  assert_eq(alignment.sequence_position_to_column(0, 32), None)
}

///|
test "Bio.Align.stockholm maps official columns to sequence positions" {
  let alignment = align_stockholm_test_sample()
  assert_eq(alignment.column_to_sequence_position(0, 16), Some(16))
  assert_eq(alignment.column_to_sequence_position(0, 17), None)
  assert_eq(alignment.column_to_sequence_position(0, 18), Some(17))
  assert_eq(alignment.column_to_sequence_position(2, 17), Some(17))
}

///|
test "Bio.Align.stockholm maps official residues between rows" {
  let alignment = align_stockholm_test_sample()
  assert_eq(alignment.map_position(0, 1, 0), Some(0))
  assert_eq(alignment.map_position(0, 2, 16), Some(16))
  assert_eq(alignment.map_position(2, 0, 17), None)
  assert_eq(alignment.map_position(0, 9, 0), None)
}

///|
test "Bio.Align.stockholm returns official aligned pairs" {
  let pairs = align_stockholm_test_sample().aligned_pairs(0, 2) catch {
    AlignStockholmError(message) => abort(message)
  }
  assert_eq(pairs.length(), 33)
  assert_eq(pairs[0], (Some(0), Some(0)))
  assert_eq(pairs[17], (None, Some(17)))
  assert_eq(pairs[18], (Some(17), Some(18)))
  assert_eq(pairs[32], (Some(31), Some(32)))
}

///|
test "Bio.Align.stockholm builds official compact coordinate path" {
  assert_eq(align_stockholm_test_sample().coordinate_path(), [
    [0, 17, 17, 32],
    [0, 17, 17, 32],
    [0, 17, 18, 33],
  ])
}

///|
test "Bio.Align.stockholm computes official first pair counts" {
  let counts = align_stockholm_test_sample().pair_counts(0, 1) catch {
    AlignStockholmError(message) => abort(message)
  }
  assert_eq(counts.pairs, 1)
  assert_eq(counts.aligned, 32)
  assert_eq(counts.identities, 21)
  assert_eq(counts.mismatches, 11)
  assert_eq(counts.gap_columns, 0)
  assert_eq(counts.double_gap_columns, 1)
  assert_eq(counts.gap_opens, 1)
}

///|
test "Bio.Align.stockholm computes official insertion pair counts" {
  let counts = align_stockholm_test_sample().pair_counts(0, 2) catch {
    AlignStockholmError(message) => abort(message)
  }
  assert_eq(counts.aligned, 32)
  assert_eq(counts.identities, 21)
  assert_eq(counts.mismatches, 11)
  assert_eq(counts.gap_columns, 1)
  assert_eq(counts.double_gap_columns, 0)
  assert_eq(counts.gap_opens, 1)
}

///|
test "Bio.Align.stockholm aggregates official all-pair counts" {
  let counts = align_stockholm_test_sample().counts()
  assert_eq(counts.pairs, 3)
  assert_eq(counts.aligned, 96)
  assert_eq(counts.identities, 58)
  assert_eq(counts.mismatches, 38)
  assert_eq(counts.gap_columns, 2)
  assert_eq(counts.double_gap_columns, 1)
  assert_eq(counts.gap_opens, 3)
}

///|
test "Bio.Align.stockholm calculates official majority consensus" {
  let consensus = align_stockholm_test_sample().consensus() catch {
    AlignStockholmError(message) => abort(message)
  }
  assert_eq(consensus, "KEIDRARXIYERFVXVHGPXVXNWIKXARFEEX")
}

///|
test "Bio.Align.stockholm computes official occupancy" {
  let occupancy = align_stockholm_test_sample().occupancy()
  assert_eq(occupancy.length(), 33)
  assert_eq(occupancy[0], 1.0)
  assert_eq(occupancy[17], 1.0 / 3.0)
  assert_eq(occupancy[32], 1.0)
}

///|
test "Bio.Align.stockholm writes canonical header and terminator" {
  let output = align_stockholm_test_write(align_stockholm_test_sample())
  assert_true(output.has_prefix("# STOCKHOLM 1.0\n"))
  assert_true(output.has_suffix("//\n"))
  assert_true(output.contains("#=GF SQ   3\n"))
}

///|
test "Bio.Align.stockholm writer restores insertion dots" {
  let output = align_stockholm_test_write(align_stockholm_test_sample())
  assert_true(output.contains("KEIDRAREIYERFVYVH.PDVKNWIKFARFEES"))
  assert_true(output.contains("HENERARGIYERFVVVH.PEVTNWLRWARFEEE"))
}

///|
test "Bio.Align.stockholm writer expands residue GR annotations" {
  let output = align_stockholm_test_write(align_stockholm_test_sample())
  assert_true(output.contains("#=GR CLF1_SCHPO/185-216"))
  assert_true(output.contains("--HHHHHHHHHHHHHHS.--HHHHHHHHHHHHH"))
}

///|
test "Bio.Align.stockholm canonical writer round-trips rows" {
  let alignment = align_stockholm_test_sample()
  let reparsed = align_stockholm_test_parse(
    align_stockholm_test_write(alignment),
  )
  assert_eq(reparsed.sequences, alignment.sequences)
  assert_eq(reparsed.operations, alignment.operations)
  assert_eq(reparsed.coordinate_path(), alignment.coordinate_path())
}

///|
test "Bio.Align.stockholm canonical writer round-trips metadata" {
  let alignment = align_stockholm_test_sample()
  let reparsed = align_stockholm_test_parse(
    align_stockholm_test_write(alignment),
  )
  assert_eq(reparsed.annotations, alignment.annotations)
  assert_eq(reparsed.references, alignment.references)
  assert_eq(reparsed.database_references, alignment.database_references)
  assert_eq(reparsed.column_annotations, alignment.column_annotations)
}

///|
test "Bio.Align.stockholm writes and parses multiple records" {
  let alignments = [
    align_stockholm_test_sample(),
    align_stockholm_test_custom(),
  ]
  let output = @src.align_stockholm_write_all(alignments) catch {
    AlignStockholmError(message) => abort(message)
  }
  let reparsed = align_stockholm_test_parse_all(output)
  assert_eq(reparsed.length(), 2)
  assert_eq(reparsed[0].annotation("identifier"), Some("HAT"))
  assert_eq(reparsed[1].num_sequences(), 2)
}

///|
test "Bio.Align.stockholm parse rejects multiple records" {
  let text = @src.align_stockholm_example_text() +
    align_stockholm_test_custom_text()
  assert_true(align_stockholm_test_raises(text))
}

///|
test "Bio.Align.stockholm parse_all accepts empty input" {
  assert_eq(align_stockholm_test_parse_all("").length(), 0)
  assert_eq(align_stockholm_test_parse_all("\n \n").length(), 0)
}

///|
test "Bio.Align.stockholm write_all accepts an empty collection" {
  let output = @src.align_stockholm_write_all([]) catch {
    AlignStockholmError(message) => abort(message)
  }
  assert_eq(output, "")
}

///|
test "Bio.Align.stockholm accepts CRLF records" {
  let text = "# STOCKHOLM 1.0\r\n#=GF SQ 1\r\nalpha A-C\r\n//\r\n"
  let alignment = align_stockholm_test_parse(text)
  assert_eq(alignment.sequences[0].aligned_sequence, "AC")
  assert_eq(alignment.operations, "MM")
  assert_eq(alignment.source_alignment_length(), 3)
  assert_eq(alignment.removed_all_gap_columns, 1)
}

///|
test "Bio.Align.stockholm ignores ordinary comments" {
  let text = "# STOCKHOLM 1.0\n# first\nalpha ACG\n# second\n//\n"
  let alignment = align_stockholm_test_parse(text)
  assert_eq(alignment.sequences[0].sequence, "ACG")
}

///|
test "Bio.Align.stockholm ignores unknown GF annotations" {
  let alignment = align_stockholm_test_custom()
  assert_eq(alignment.annotation("XX"), None)
  assert_eq(alignment.annotation("ignored alignment annotation"), None)
}

///|
test "Bio.Align.stockholm preserves custom GS annotations" {
  let sequence = align_stockholm_test_custom().sequences[0]
  assert_eq(sequence.annotation("ZZ"), Some("custom sequence value"))
  assert_eq(sequence.annotations[2].feature, "ZZ")
}

///|
test "Bio.Align.stockholm preserves custom GR annotations" {
  let sequence = align_stockholm_test_custom().sequences[1]
  assert_eq(sequence.letter_annotation("ZZ"), Some("wxyz"))
  assert_eq(sequence.letter_annotations[0].aligned_value, "wxyz")
}

///|
test "Bio.Align.stockholm preserves custom GC annotations" {
  let alignment = align_stockholm_test_custom()
  assert_eq(alignment.column_annotation("ZZ"), Some("abdg"))
}

///|
test "Bio.Align.stockholm preserves repeated AU values in order" {
  assert_eq(align_stockholm_test_custom().annotation_values("author"), [
    "Alice;", "Bob;",
  ])
}

///|
test "Bio.Align.stockholm joins repeated search methods" {
  assert_eq(
    align_stockholm_test_custom().annotation("search method"),
    Some("first command second command"),
  )
}

///|
test "Bio.Align.stockholm joins continued wikipedia values" {
  assert_eq(align_stockholm_test_custom().annotation_values("wikipedia"), [
    "LongName",
  ])
}

///|
test "Bio.Align.stockholm parses GS description" {
  assert_eq(
    align_stockholm_test_custom().sequences[0].description,
    "Alpha description",
  )
}

///|
test "Bio.Align.stockholm preserves repeated GS database references" {
  assert_eq(align_stockholm_test_custom().sequences[0].dbxrefs, [
    "DBX; 1;", "DBX; 2;",
  ])
}

///|
test "Bio.Align.stockholm maps standard GS organism" {
  let sequence = align_stockholm_test_custom().sequences[0]
  assert_eq(sequence.annotation("organism"), Some("organism one"))
}

///|
test "Bio.Align.stockholm parses nested domains" {
  let domains = align_stockholm_test_custom().nested_domains
  assert_eq(domains.length(), 1)
  assert_eq(domains[0].accession, "PF00001")
  assert_eq(domains[0].location, "1-4")
}

///|
test "Bio.Align.stockholm parses database reference comments" {
  let references = align_stockholm_test_custom().database_references
  assert_eq(references.length(), 1)
  assert_eq(references[0].reference, "DB; ONE;")
  assert_eq(references[0].comment, "database comment")
}

///|
test "Bio.Align.stockholm attaches pending reference comments" {
  let reference = align_stockholm_test_custom().references[0]
  assert_eq(reference.number, 2)
  assert_eq(reference.comment, "reference note")
  assert_eq(reference.medline, "12345")
}

///|
test "Bio.Align.stockholm joins multiline reference fields" {
  let reference = align_stockholm_test_custom().references[0]
  assert_eq(reference.title, "title part one title part two")
  assert_eq(reference.author, "Example Author;")
  assert_eq(reference.location, "Example Journal.")
}

///|
test "Bio.Align.stockholm removes all-gap source columns" {
  let alignment = align_stockholm_test_custom()
  assert_eq(alignment.source_alignment_length(), 7)
  assert_eq(alignment.alignment_length(), 4)
  assert_eq(alignment.removed_all_gap_columns, 3)
  assert_eq(alignment.sequences[0].aligned_sequence, "A-CG")
  assert_eq(alignment.sequences[1].aligned_sequence, "ATCG")
}

///|
test "Bio.Align.stockholm preserves match and insertion operations" {
  assert_eq(align_stockholm_test_custom().operations, "MIMM")
}

///|
test "Bio.Align.stockholm compresses GC with removed columns" {
  let alignment = align_stockholm_test_custom()
  assert_eq(alignment.column_annotation("RF"), Some("1247"))
  assert_eq(alignment.column_annotation("ZZ"), Some("abdg"))
}

///|
test "Bio.Align.stockholm compresses aligned GR with removed columns" {
  let annotation = align_stockholm_test_custom().sequences[0].letter_annotations[0]
  assert_eq(annotation.value, "xyz")
  assert_eq(annotation.aligned_value, "x.yz")
}

///|
test "Bio.Align.stockholm maps custom compact coordinates" {
  let alignment = align_stockholm_test_custom()
  assert_eq(alignment.coordinate_path(), [[0, 1, 1, 3], [0, 1, 2, 4]])
  assert_eq(alignment.map_position(1, 0, 1), None)
  assert_eq(alignment.map_position(1, 0, 2), Some(1))
}

///|
test "Bio.Align.stockholm computes custom compact pair counts" {
  let counts = align_stockholm_test_custom().pair_counts(0, 1) catch {
    AlignStockholmError(message) => abort(message)
  }
  assert_eq(counts.aligned, 3)
  assert_eq(counts.identities, 3)
  assert_eq(counts.mismatches, 0)
  assert_eq(counts.gap_columns, 1)
  assert_eq(counts.gap_opens, 1)
  assert_eq(counts.identity(), 1.0)
}

///|
test "Bio.Align.stockholm computes custom compact occupancy" {
  assert_eq(align_stockholm_test_custom().occupancy(), [1.0, 0.5, 1.0, 1.0])
}

///|
test "Bio.Align.stockholm computes custom compact consensus" {
  let consensus = align_stockholm_test_custom().consensus() catch {
    AlignStockholmError(message) => abort(message)
  }
  assert_eq(consensus, "ATCG")
}

///|
test "Bio.Align.stockholm slices rows and operations" {
  let sliced = align_stockholm_test_custom().slice_columns(1, 4) catch {
    AlignStockholmError(message) => abort(message)
  }
  assert_eq(sliced.sequences[0].aligned_sequence, "-CG")
  assert_eq(sliced.sequences[1].aligned_sequence, "TCG")
  assert_eq(sliced.operations, "IMM")
  assert_eq(sliced.source_alignment_length(), 3)
}

///|
test "Bio.Align.stockholm slices GC annotations" {
  let sliced = align_stockholm_test_custom().slice_columns(1, 4) catch {
    AlignStockholmError(message) => abort(message)
  }
  assert_eq(sliced.column_annotation("RF"), Some("247"))
  assert_eq(sliced.column_annotation("ZZ"), Some("bdg"))
}

///|
test "Bio.Align.stockholm slices GR annotations" {
  let sliced = align_stockholm_test_custom().slice_columns(1, 4) catch {
    AlignStockholmError(message) => abort(message)
  }
  let first = sliced.sequences[0].letter_annotations[0]
  assert_eq(first.aligned_value, ".yz")
  assert_eq(first.value, "yz")
  assert_eq(sliced.sequences[1].letter_annotations[0].value, "xyz")
}

///|
test "Bio.Align.stockholm preserves metadata across column slicing" {
  let alignment = align_stockholm_test_custom()
  let sliced = alignment.slice_columns(1, 4) catch {
    AlignStockholmError(message) => abort(message)
  }
  assert_eq(sliced.annotations, alignment.annotations)
  assert_eq(sliced.references, alignment.references)
  assert_eq(sliced.database_references, alignment.database_references)
  assert_eq(sliced.nested_domains, alignment.nested_domains)
}

///|
test "Bio.Align.stockholm constructs compact rows directly" {
  let alignment = align_stockholm_test_from(["alpha", "beta"], [
    "A.-C--G", "AT-C--G",
  ])
  assert_eq(alignment.operations, "MIMM")
  assert_eq(alignment.source_alignment_length(), 7)
  assert_eq(alignment.removed_all_gap_columns, 3)
}

///|
test "Bio.Align.stockholm sequence constructor normalizes dots" {
  let sequence = align_stockholm_test_sequence("alpha", "A.C-")
  assert_eq(sequence.aligned_sequence, "A-C-")
  assert_eq(sequence.sequence, "AC")
}

///|
test "Bio.Align.stockholm annotation constructor defaults feature" {
  let annotation = @src.AlignStockholmAnnotation::create("ZZ", "value") catch {
    AlignStockholmError(message) => abort(message)
  }
  assert_eq(annotation.code, "ZZ")
  assert_eq(annotation.feature, "ZZ")
  assert_eq(annotation.value, "value")
}

///|
test "Bio.Align.stockholm annotation constructor accepts custom feature" {
  let annotation = @src.AlignStockholmAnnotation::create(
    "ZZ",
    "value",
    feature=Some("custom feature"),
  ) catch {
    AlignStockholmError(message) => abort(message)
  }
  assert_eq(annotation.feature, "custom feature")
}

///|
test "Bio.Align.stockholm column constructor maps standard feature" {
  let annotation = @src.AlignStockholmColumnAnnotation::create(
    "SS_cons", "<<>>",
  ) catch {
    AlignStockholmError(message) => abort(message)
  }
  assert_eq(annotation.feature, "consensus secondary structure")
}

///|
test "Bio.Align.stockholm letter constructor maps standard feature" {
  let annotation = @src.AlignStockholmLetterAnnotation::create(
    "PP", "999", "9.99",
  ) catch {
    AlignStockholmError(message) => abort(message)
  }
  assert_eq(annotation.feature, "posterior probability")
}

///|
test "Bio.Align.stockholm alignment constructor infers deletion operations" {
  let sequences = [
    align_stockholm_test_sequence("alpha", "A-C"),
    align_stockholm_test_sequence("beta", "ATC"),
  ]
  let alignment = @src.AlignStockholmAlignment::create(sequences) catch {
    AlignStockholmError(message) => abort(message)
  }
  assert_eq(alignment.operations, "MDM")
  assert_eq(alignment.source_alignment_length(), 3)
}

///|
test "Bio.Align.stockholm writer preserves custom GS GR and GC" {
  let output = align_stockholm_test_write(align_stockholm_test_custom())
  assert_true(output.contains("#=GS alpha  ZZ custom sequence value"))
  assert_true(output.contains("#=GR beta   ZZ "))
  assert_true(output.contains("#=GC ZZ "))
  let reparsed = align_stockholm_test_parse(output)
  assert_eq(
    reparsed.sequences[0].annotation("ZZ"),
    Some("custom sequence value"),
  )
  assert_eq(reparsed.sequences[1].letter_annotation("ZZ"), Some("wxyz"))
  assert_eq(reparsed.column_annotation("ZZ"), Some("abdg"))
}

///|
test "Bio.Align.stockholm writer rejects unknown GF annotations" {
  let annotation = @src.AlignStockholmAnnotation::create("ZZ", "value") catch {
    AlignStockholmError(message) => abort(message)
  }
  let sequence = align_stockholm_test_sequence("alpha", "ACG")
  let alignment = @src.AlignStockholmAlignment::create([sequence], annotations=[
    annotation,
  ]) catch {
    AlignStockholmError(message) => abort(message)
  }
  assert_true(align_stockholm_test_write_raises(alignment))
}

///|
test "Bio.Align.stockholm writer wraps long comments at 79 columns" {
  let comment = @src.AlignStockholmAnnotation::create(
    "CC",
    "This deliberately long Stockholm comment contains enough words to " +
    "require canonical wrapping across multiple GF comment lines.",
    feature=Some("comment"),
  ) catch {
    AlignStockholmError(message) => abort(message)
  }
  let alignment = @src.AlignStockholmAlignment::create(
    [align_stockholm_test_sequence("alpha", "ACG")],
    annotations=[comment],
  ) catch {
    AlignStockholmError(message) => abort(message)
  }
  let output = align_stockholm_test_write(alignment)
  let mut comment_lines = 0
  for line in output.split("\n") {
    if line.to_owned().has_prefix("#=GF CC") {
      comment_lines = comment_lines + 1
      assert_true(line.length() <= 79)
    }
  }
  assert_true(comment_lines >= 2)
}

///|
test "Bio.Align.stockholm compares residues case-insensitively" {
  let counts = align_stockholm_test_from(["alpha", "beta"], ["aCg", "AcG"]).counts()
  assert_eq(counts.identities, 3)
  assert_eq(counts.mismatches, 0)
}

///|
test "Bio.Align.stockholm returns zero identity without aligned residues" {
  let counts = align_stockholm_test_from(["alpha", "beta"], ["A-", "-B"]).counts()
  assert_eq(counts.aligned, 0)
  assert_eq(counts.gap_columns, 2)
  assert_eq(counts.identity(), 0.0)
}

///|
test "Bio.Align.stockholm consensus is case-insensitive" {
  let alignment = align_stockholm_test_from(["alpha", "beta", "gamma"], [
    "aCG", "ACG", "tCG",
  ])
  let consensus = alignment.consensus() catch {
    AlignStockholmError(message) => abort(message)
  }
  assert_eq(consensus, "ACG")
}

///|
test "Bio.Align.stockholm consensus honors strict thresholds" {
  let alignment = align_stockholm_test_from(["alpha", "beta", "gamma"], [
    "ACG", "TCG", "TCG",
  ])
  let consensus = alignment.consensus(minimum_fraction=0.75) catch {
    AlignStockholmError(message) => abort(message)
  }
  assert_eq(consensus, "XCG")
}

///|
test "Bio.Align.stockholm consensus accepts a custom ambiguous symbol" {
  let alignment = align_stockholm_test_from(["alpha", "beta"], ["AC", "TC"])
  let consensus = alignment.consensus(minimum_fraction=0.75, ambiguous="?") catch {
    AlignStockholmError(message) => abort(message)
  }
  assert_eq(consensus, "?C")
}

///|
test "Bio.Align.stockholm rejects invalid aligned-pair row indexes" {
  let alignment = align_stockholm_test_sample()
  let low = try {
    ignore(alignment.aligned_pairs(-1, 0))
    false
  } catch {
    AlignStockholmError(_) => true
  }
  let high = try {
    ignore(alignment.aligned_pairs(0, 3))
    false
  } catch {
    AlignStockholmError(_) => true
  }
  assert_true(low)
  assert_true(high)
}

///|
test "Bio.Align.stockholm rejects invalid pair-count row indexes" {
  let alignment = align_stockholm_test_sample()
  let low = try {
    ignore(alignment.pair_counts(-1, 0))
    false
  } catch {
    AlignStockholmError(_) => true
  }
  let high = try {
    ignore(alignment.pair_counts(0, 3))
    false
  } catch {
    AlignStockholmError(_) => true
  }
  assert_true(low)
  assert_true(high)
}

///|
test "Bio.Align.stockholm returns None for invalid mappings" {
  let alignment = align_stockholm_test_sample()
  assert_eq(alignment.sequence_position_to_column(-1, 0), None)
  assert_eq(alignment.sequence_position_to_column(3, 0), None)
  assert_eq(alignment.sequence_position_to_column(0, -1), None)
  assert_eq(alignment.column_to_sequence_position(-1, 0), None)
  assert_eq(alignment.column_to_sequence_position(3, 0), None)
  assert_eq(alignment.column_to_sequence_position(0, 33), None)
}

///|
test "Bio.Align.stockholm rejects invalid column slices" {
  let alignment = align_stockholm_test_sample()
  let first = try {
    ignore(alignment.slice_columns(-1, 2))
    false
  } catch {
    AlignStockholmError(_) => true
  }
  let second = try {
    ignore(alignment.slice_columns(2, 2))
    false
  } catch {
    AlignStockholmError(_) => true
  }
  let third = try {
    ignore(alignment.slice_columns(0, 34))
    false
  } catch {
    AlignStockholmError(_) => true
  }
  assert_true(first)
  assert_true(second)
  assert_true(third)
}

///|
test "Bio.Align.stockholm rejects a negative consensus threshold" {
  let raised = try {
    ignore(align_stockholm_test_sample().consensus(minimum_fraction=-0.1))
    false
  } catch {
    AlignStockholmError(_) => true
  }
  assert_true(raised)
}

///|
test "Bio.Align.stockholm rejects a consensus threshold above one" {
  let raised = try {
    ignore(align_stockholm_test_sample().consensus(minimum_fraction=1.1))
    false
  } catch {
    AlignStockholmError(_) => true
  }
  assert_true(raised)
}

///|
test "Bio.Align.stockholm rejects a multi-character ambiguous symbol" {
  let raised = try {
    ignore(align_stockholm_test_sample().consensus(ambiguous="NN"))
    false
  } catch {
    AlignStockholmError(_) => true
  }
  assert_true(raised)
}

///|
test "Bio.Align.stockholm rejects empty single-record input" {
  assert_true(align_stockholm_test_raises(""))
  assert_true(align_stockholm_test_raises("\n\n"))
}

///|
test "Bio.Align.stockholm rejects an invalid header" {
  assert_true(align_stockholm_test_raises("# STOCKHOLM 1.1\nalpha ACG\n//\n"))
  assert_true(align_stockholm_test_raises("STOCKHOLM 1.0\nalpha ACG\n//\n"))
}

///|
test "Bio.Align.stockholm rejects a missing terminator" {
  assert_true(
    align_stockholm_test_parse_all_raises("# STOCKHOLM 1.0\nalpha ACG\n"),
  )
}

///|
test "Bio.Align.stockholm rejects a nested header without terminator" {
  let text = "# STOCKHOLM 1.0\nalpha ACG\n" + "# STOCKHOLM 1.0\nbeta ACG\n//\n"
  assert_true(align_stockholm_test_parse_all_raises(text))
}

///|
test "Bio.Align.stockholm rejects an empty record" {
  assert_true(align_stockholm_test_raises("# STOCKHOLM 1.0\n//\n"))
  assert_true(
    align_stockholm_test_raises("# STOCKHOLM 1.0\n#=GF ID empty\n//\n"),
  )
}

///|
test "Bio.Align.stockholm rejects duplicate identifiers" {
  let text = "# STOCKHOLM 1.0\nalpha ACG\nalpha ACG\n//\n"
  assert_true(align_stockholm_test_raises(text))
}

///|
test "Bio.Align.stockholm rejects wrapped sequence rows" {
  let text = "# STOCKHOLM 1.0\nalpha ACG\nbeta ACG\nalpha TTT\n//\n"
  assert_true(align_stockholm_test_raises(text))
}

///|
test "Bio.Align.stockholm rejects unequal row widths" {
  let text = "# STOCKHOLM 1.0\nalpha ACG\nbeta AC\n//\n"
  assert_true(align_stockholm_test_raises(text))
}

///|
test "Bio.Align.stockholm rejects malformed sequence lines" {
  let text = "# STOCKHOLM 1.0\nalpha ACG extra\n//\n"
  assert_true(align_stockholm_test_raises(text))
}

///|
test "Bio.Align.stockholm rejects mixed insertion and deletion gaps" {
  let text = "# STOCKHOLM 1.0\nalpha A.C\nbeta A-C\n//\n"
  assert_true(align_stockholm_test_raises(text))
}

///|
test "Bio.Align.stockholm rejects all-gap-only records" {
  let text = "# STOCKHOLM 1.0\nalpha .--\nbeta .--\n//\n"
  assert_true(align_stockholm_test_raises(text))
}

///|
test "Bio.Align.stockholm rejects malformed GF annotations" {
  assert_true(
    align_stockholm_test_raises("# STOCKHOLM 1.0\n#=GF ID\nalpha ACG\n//\n"),
  )
  assert_true(
    align_stockholm_test_raises(
      "# STOCKHOLM 1.0\n#=GFID value\nalpha ACG\n//\n",
    ),
  )
}

///|
test "Bio.Align.stockholm rejects malformed GS annotations" {
  assert_true(
    align_stockholm_test_raises("# STOCKHOLM 1.0\n#=GS alpha\nalpha ACG\n//\n"),
  )
  assert_true(
    align_stockholm_test_raises(
      "# STOCKHOLM 1.0\n#=GSalpha AC A1\nalpha ACG\n//\n",
    ),
  )
}

///|
test "Bio.Align.stockholm rejects malformed GC annotations" {
  assert_true(
    align_stockholm_test_raises("# STOCKHOLM 1.0\nalpha ACG\n#=GC RF\n//\n"),
  )
  assert_true(
    align_stockholm_test_raises("# STOCKHOLM 1.0\nalpha ACG\n#=GCRF xxx\n//\n"),
  )
}

///|
test "Bio.Align.stockholm rejects malformed GR annotations" {
  assert_true(
    align_stockholm_test_raises(
      "# STOCKHOLM 1.0\nalpha ACG\n#=GR alpha SS\n//\n",
    ),
  )
  assert_true(
    align_stockholm_test_raises(
      "# STOCKHOLM 1.0\nalpha ACG\n#=GRalpha SS xxx\n//\n",
    ),
  )
}

///|
test "Bio.Align.stockholm requires GR after its row" {
  let before = "# STOCKHOLM 1.0\n#=GR alpha SS xxx\nalpha ACG\n//\n"
  let wrong = "# STOCKHOLM 1.0\nalpha ACG\nbeta ACG\n" +
    "#=GR alpha SS xxx\n//\n"
  assert_true(align_stockholm_test_raises(before))
  assert_true(align_stockholm_test_raises(wrong))
}

///|
test "Bio.Align.stockholm rejects duplicate GR annotations" {
  let text = "# STOCKHOLM 1.0\nalpha ACG\n" +
    "#=GR alpha SS xxx\n#=GR alpha SS yyy\n//\n"
  assert_true(align_stockholm_test_raises(text))
}

///|
test "Bio.Align.stockholm validates GF SQ row count" {
  let text = "# STOCKHOLM 1.0\n#=GF SQ 2\nalpha ACG\n//\n"
  assert_true(align_stockholm_test_raises(text))
}

///|
test "Bio.Align.stockholm rejects duplicate GF SQ" {
  let text = "# STOCKHOLM 1.0\n#=GF SQ 1\n#=GF SQ 1\nalpha ACG\n//\n"
  assert_true(align_stockholm_test_raises(text))
}

///|
test "Bio.Align.stockholm rejects malformed GF SQ integers" {
  let nonnumeric = "# STOCKHOLM 1.0\n#=GF SQ one\nalpha ACG\n//\n"
  let zero = "# STOCKHOLM 1.0\n#=GF SQ 0\nalpha ACG\n//\n"
  let overflow = "# STOCKHOLM 1.0\n#=GF SQ 2147483648\nalpha ACG\n//\n"
  assert_true(align_stockholm_test_raises(nonnumeric))
  assert_true(align_stockholm_test_raises(zero))
  assert_true(align_stockholm_test_raises(overflow))
}

///|
test "Bio.Align.stockholm rejects duplicate single GF annotations" {
  let text = "# STOCKHOLM 1.0\n#=GF ID one\n#=GF ID two\nalpha ACG\n//\n"
  assert_true(align_stockholm_test_raises(text))
}

///|
test "Bio.Align.stockholm rejects reference fields before RN" {
  let rm = "# STOCKHOLM 1.0\n#=GF RM 123\nalpha ACG\n//\n"
  let rt = "# STOCKHOLM 1.0\n#=GF RT title\nalpha ACG\n//\n"
  assert_true(align_stockholm_test_raises(rm))
  assert_true(align_stockholm_test_raises(rt))
}

///|
test "Bio.Align.stockholm rejects dangling reference comments" {
  let text = "# STOCKHOLM 1.0\n#=GF RC note\nalpha ACG\n//\n"
  assert_true(align_stockholm_test_raises(text))
}

///|
test "Bio.Align.stockholm rejects malformed reference numbers" {
  let brackets = "# STOCKHOLM 1.0\n#=GF RN 1\nalpha ACG\n//\n"
  let zero = "# STOCKHOLM 1.0\n#=GF RN [0]\nalpha ACG\n//\n"
  assert_true(align_stockholm_test_raises(brackets))
  assert_true(align_stockholm_test_raises(zero))
}

///|
test "Bio.Align.stockholm rejects duplicate reference medline" {
  let text = "# STOCKHOLM 1.0\n#=GF RN [1]\n" +
    "#=GF RM 123\n#=GF RM 456\nalpha ACG\n//\n"
  assert_true(align_stockholm_test_raises(text))
}

///|
test "Bio.Align.stockholm rejects database comments before DR" {
  let text = "# STOCKHOLM 1.0\n#=GF DC comment\nalpha ACG\n//\n"
  assert_true(align_stockholm_test_raises(text))
}

///|
test "Bio.Align.stockholm rejects duplicate database comments" {
  let text = "# STOCKHOLM 1.0\n#=GF DR DB; ONE;\n" +
    "#=GF DC first\n#=GF DC second\nalpha ACG\n//\n"
  assert_true(align_stockholm_test_raises(text))
}

///|
test "Bio.Align.stockholm rejects nested locations before NE" {
  let text = "# STOCKHOLM 1.0\n#=GF NL 1-3\nalpha ACG\n//\n"
  assert_true(align_stockholm_test_raises(text))
}

///|
test "Bio.Align.stockholm rejects duplicate nested locations" {
  let text = "# STOCKHOLM 1.0\n#=GF NE PF00001\n" +
    "#=GF NL 1-3\n#=GF NL 2-4\nalpha ACG\n//\n"
  assert_true(align_stockholm_test_raises(text))
}

///|
test "Bio.Align.stockholm rejects GS for unknown rows" {
  let text = "# STOCKHOLM 1.0\n#=GS missing AC A1\nalpha ACG\n//\n"
  assert_true(align_stockholm_test_raises(text))
}

///|
test "Bio.Align.stockholm rejects duplicate GS annotations" {
  let text = "# STOCKHOLM 1.0\n#=GS alpha AC A1\n" +
    "#=GS alpha AC A2\nalpha ACG\n//\n"
  assert_true(align_stockholm_test_raises(text))
}

///|
test "Bio.Align.stockholm rejects duplicate GS descriptions" {
  let text = "# STOCKHOLM 1.0\n#=GS alpha DE first\n" +
    "#=GS alpha DE second\nalpha ACG\n//\n"
  assert_true(align_stockholm_test_raises(text))
}

///|
test "Bio.Align.stockholm validates GC source width" {
  let short = "# STOCKHOLM 1.0\nalpha ACG\n#=GC RF xx\n//\n"
  let long = "# STOCKHOLM 1.0\nalpha ACG\n#=GC RF xxxx\n//\n"
  assert_true(align_stockholm_test_raises(short))
  assert_true(align_stockholm_test_raises(long))
}

///|
test "Bio.Align.stockholm validates GR source width" {
  let text = "# STOCKHOLM 1.0\nalpha ACG\n#=GR alpha SS xx\n//\n"
  assert_true(align_stockholm_test_raises(text))
}

///|
test "Bio.Align.stockholm validates GR residue count" {
  let text = "# STOCKHOLM 1.0\nalpha A-C\n#=GR alpha SS xyz\n//\n"
  assert_true(align_stockholm_test_raises(text))
}

///|
test "Bio.Align.stockholm parses CSA residue annotations" {
  let text = "# STOCKHOLM 1.0\nalpha A-C\n" +
    "#=GR alpha CSA x-y\nbeta ATC\n//\n"
  let alignment = align_stockholm_test_parse(text)
  let annotation = alignment.sequences[0].letter_annotations[0]
  assert_eq(annotation.feature, "Catalytic Site Atlas")
  assert_eq(annotation.value, "xy")
  assert_eq(annotation.aligned_value, "x-y")
}

///|
test "Bio.Align.stockholm constructors reject missing identifiers" {
  assert_true(align_stockholm_test_from_raises([], []))
  assert_true(align_stockholm_test_from_raises(["alpha"], []))
}

///|
test "Bio.Align.stockholm constructors reject unequal widths" {
  assert_true(
    align_stockholm_test_from_raises(["alpha", "beta"], ["ACG", "AC"]),
  )
}

///|
test "Bio.Align.stockholm constructors reject empty rows" {
  assert_true(align_stockholm_test_from_raises(["alpha"], [""]))
}

///|
test "Bio.Align.stockholm constructors reject all-gap columns only" {
  assert_true(
    align_stockholm_test_from_raises(["alpha", "beta"], ["---", "..."]),
  )
}

///|
test "Bio.Align.stockholm sequence constructor rejects empty rows" {
  let raised = try {
    ignore(@src.AlignStockholmSequence::create("alpha", ""))
    false
  } catch {
    AlignStockholmError(_) => true
  }
  assert_true(raised)
}

///|
test "Bio.Align.stockholm sequence constructor rejects whitespace IDs" {
  let raised = try {
    ignore(@src.AlignStockholmSequence::create("bad id", "ACG"))
    false
  } catch {
    AlignStockholmError(_) => true
  }
  assert_true(raised)
}

///|
test "Bio.Align.stockholm annotation constructors reject empty codes" {
  let generic = try {
    ignore(@src.AlignStockholmAnnotation::create("", "value"))
    false
  } catch {
    AlignStockholmError(_) => true
  }
  let column = try {
    ignore(@src.AlignStockholmColumnAnnotation::create("", "ACG"))
    false
  } catch {
    AlignStockholmError(_) => true
  }
  let letter = try {
    ignore(@src.AlignStockholmLetterAnnotation::create("", "xyz", "xyz"))
    false
  } catch {
    AlignStockholmError(_) => true
  }
  assert_true(generic)
  assert_true(column)
  assert_true(letter)
}

///|
test "Bio.Align.stockholm annotation constructors reject whitespace codes" {
  let raised = try {
    ignore(@src.AlignStockholmAnnotation::create("bad code", "value"))
    false
  } catch {
    AlignStockholmError(_) => true
  }
  assert_true(raised)
}

///|
test "Bio.Align.stockholm alignment constructor rejects duplicate IDs" {
  let sequences = [
    align_stockholm_test_sequence("alpha", "ACG"),
    align_stockholm_test_sequence("alpha", "ACG"),
  ]
  let raised = try {
    ignore(@src.AlignStockholmAlignment::create(sequences))
    false
  } catch {
    AlignStockholmError(_) => true
  }
  assert_true(raised)
}

///|
test "Bio.Align.stockholm alignment constructor validates operations length" {
  let sequence = align_stockholm_test_sequence("alpha", "ACG")
  let raised = try {
    ignore(@src.AlignStockholmAlignment::create([sequence], operations="MM"))
    false
  } catch {
    AlignStockholmError(_) => true
  }
  assert_true(raised)
}

///|
test "Bio.Align.stockholm alignment constructor validates operation symbols" {
  let sequence = align_stockholm_test_sequence("alpha", "ACG")
  let raised = try {
    ignore(@src.AlignStockholmAlignment::create([sequence], operations="MXM"))
    false
  } catch {
    AlignStockholmError(_) => true
  }
  assert_true(raised)
}

///|
test "Bio.Align.stockholm alignment constructor rejects match gaps" {
  let sequences = [
    align_stockholm_test_sequence("alpha", "A-C"),
    align_stockholm_test_sequence("beta", "ATC"),
  ]
  let raised = try {
    ignore(@src.AlignStockholmAlignment::create(sequences, operations="MMM"))
    false
  } catch {
    AlignStockholmError(_) => true
  }
  assert_true(raised)
}

///|
test "Bio.Align.stockholm alignment constructor validates source metadata" {
  let sequence = align_stockholm_test_sequence("alpha", "ACG")
  let first = try {
    ignore(
      @src.AlignStockholmAlignment::create(
        [sequence],
        source_columns=4,
        removed_all_gap_columns=0,
      ),
    )
    false
  } catch {
    AlignStockholmError(_) => true
  }
  let second = try {
    ignore(
      @src.AlignStockholmAlignment::create(
        [sequence],
        source_columns=3,
        removed_all_gap_columns=-1,
      ),
    )
    false
  } catch {
    AlignStockholmError(_) => true
  }
  assert_true(first)
  assert_true(second)
}

///|
test "Bio.Align.stockholm alignment constructor validates GC width" {
  let annotation = @src.AlignStockholmColumnAnnotation::create("RF", "xx") catch {
    AlignStockholmError(message) => abort(message)
  }
  let raised = try {
    ignore(
      @src.AlignStockholmAlignment::create(
        [align_stockholm_test_sequence("alpha", "ACG")],
        column_annotations=[annotation],
      ),
    )
    false
  } catch {
    AlignStockholmError(_) => true
  }
  assert_true(raised)
}

///|
test "Bio.Align.stockholm sequence constructor validates GR residue width" {
  let annotation = @src.AlignStockholmLetterAnnotation::create(
    "SS", "xx", "x.xx",
  ) catch {
    AlignStockholmError(message) => abort(message)
  }
  let raised = try {
    ignore(
      @src.AlignStockholmSequence::create("alpha", "A-CG", letter_annotations=[
        annotation,
      ]),
    )
    false
  } catch {
    AlignStockholmError(_) => true
  }
  assert_true(raised)
}

///|
test "Bio.Align.stockholm sequence constructor validates aligned GR width" {
  let annotation = @src.AlignStockholmLetterAnnotation::create(
    "SS", "xyz", "xyz",
  ) catch {
    AlignStockholmError(message) => abort(message)
  }
  let raised = try {
    ignore(
      @src.AlignStockholmSequence::create("alpha", "A-CG", letter_annotations=[
        annotation,
      ]),
    )
    false
  } catch {
    AlignStockholmError(_) => true
  }
  assert_true(raised)
}
