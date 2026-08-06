// Black-box tests for Biopython Bio.Align.clustal-compatible support.

///|
fn align_clustal_test_parse(text : String) -> @src.AlignClustalAlignment {
  @src.align_clustal_parse(text) catch {
    AlignClustalError(message) =>
      abort("valid CLUSTAL input failed: " + message)
  }
}

///|
fn align_clustal_test_rejects(text : String) -> Bool {
  try {
    ignore(@src.align_clustal_parse(text))
    false
  } catch {
    AlignClustalError(_) => true
  }
}

///|
fn align_clustal_test_write(
  alignment : @src.AlignClustalAlignment,
  block_width : Int,
  name_width : Int,
  include_counts : Bool,
  include_consensus : Bool,
) -> String {
  @src.align_clustal_write(
    alignment,
    block_width~,
    name_width~,
    include_counts~,
    include_consensus~,
  ) catch {
    AlignClustalError(message) => abort("CLUSTAL writing failed: " + message)
  }
}

///|
fn align_clustal_test_write_rejects(
  alignment : @src.AlignClustalAlignment,
  block_width : Int,
  name_width : Int,
) -> Bool {
  try {
    ignore(@src.align_clustal_write(alignment, block_width~, name_width~))
    false
  } catch {
    AlignClustalError(_) => true
  }
}

///|
fn align_clustal_test_from_rejects(
  ids : Array[String],
  rows : Array[String],
  consensus : String?,
) -> Bool {
  try {
    ignore(@src.align_clustal_from_aligned(ids, rows, consensus~))
    false
  } catch {
    AlignClustalError(_) => true
  }
}

///|
fn align_clustal_test_metadata_rejects(
  program : String,
  version : String,
) -> Bool {
  try {
    ignore(@src.AlignClustalMetadata::create(program, version~))
    false
  } catch {
    AlignClustalError(_) => true
  }
}

///|
fn align_clustal_test_sequence_rejects(id : String, row : String) -> Bool {
  try {
    ignore(@src.AlignClustalSequence::create(id, row))
    false
  } catch {
    AlignClustalError(_) => true
  }
}

///|
fn align_clustal_test_sample() -> @src.AlignClustalAlignment {
  align_clustal_test_parse(@src.align_clustal_example_text())
}

///|
fn align_clustal_test_no_consensus_text() -> String {
  "Kalign (2.0) alignment in ClustalW format\n" +
  "\n" +
  "\n" +
  "Test1seq      GCTGGGGATGGAGAGGGAACAGAGT-T\n" +
  "AT3G20900     GCTGGGGATGGAGAGGGAACAGAGTAG\n"
}

///|
fn align_clustal_test_sparse_consensus_text() -> String {
  "MUSCLE (3.8) multiple sequence alignment\n" +
  "\n" +
  "\n" +
  "alpha      AC-G 3\n" +
  "beta       ATCG 4\n" +
  "\n" +
  "alpha      TT 5\n" +
  "beta       TT 6\n" +
  "           **\n"
}

///|
fn align_clustal_test_probcons_text() -> String {
  "PROBCONS version 1.12 multiple sequence alignment\n" +
  "\n" +
  "plas_horvu    D-VLLGANGGVL\n" +
  "plas_chlre    -AVKLGADSGAL\n" +
  "plas_anava    --VKLGSDKGLL\n" +
  "              ::    .   \n"
}

///|
fn align_clustal_test_msaprobs_text() -> String {
  "MSAPROBS version 0.9.7 multiple sequence alignment\n" +
  "\n" +
  "V_Harveyi_PATH    MKNW--------IKV\n" +
  "B_subtilis_YXEM   MKMKKWAAAAAATVL\n" +
  "                  **         :\n"
}

///|
fn align_clustal_test_clustalw_text() -> String {
  "CLUSTAL W (1.81) multiple sequence alignment\n" +
  "\n" +
  "\n" +
  "gi|4959044|gb|AAD34209.1|AF069      MENSDSNDKG 10\n" +
  "gi|671626|emb|CAA85685.1|           ---------M 1\n" +
  "                                    *\n"
}

///|
test "Bio.Align.clustal parses CLUSTAL metadata" {
  let alignment = align_clustal_test_sample()
  assert_eq(alignment.metadata.program, "CLUSTAL")
  assert_eq(alignment.metadata.version, "2.1")
  assert_true(alignment.metadata.header.contains("CLUSTAL W (2.1)"))
}

///|
test "Bio.Align.clustal parses sample shape" {
  let alignment = align_clustal_test_sample()
  assert_eq(alignment.num_sequences(), 3)
  assert_eq(alignment.alignment_length(), 10)
}

///|
test "Bio.Align.clustal joins interleaved rows" {
  let alignment = align_clustal_test_sample()
  assert_eq(alignment.sequences[0].aligned_sequence, "MKT--AILGH")
  assert_eq(alignment.sequences[1].aligned_sequence, "M-TQQAILG-")
  assert_eq(alignment.sequences[2].aligned_sequence, "MKT--AVLGH")
}

///|
test "Bio.Align.clustal derives ungapped rows" {
  let alignment = align_clustal_test_sample()
  assert_eq(alignment.sequences[0].sequence, "MKTAILGH")
  assert_eq(alignment.sequences[1].sequence, "MTQQAILG")
  assert_eq(alignment.sequences[2].sequence, "MKTAVLGH")
}

///|
test "Bio.Align.clustal retains column consensus" {
  assert_eq(align_clustal_test_sample().consensus, Some("* *  *:** "))
}

///|
test "Bio.Align.clustal accepts CRLF input" {
  let text = @src.align_clustal_example_text().replace(old="\n", new="\r\n")
  let alignment = align_clustal_test_parse(text)
  assert_eq(alignment.num_sequences(), 3)
  assert_eq(alignment.metadata.version, "2.1")
}

///|
test "Bio.Align.clustal parses CLUSTAL W official-style header" {
  let alignment = align_clustal_test_parse(align_clustal_test_clustalw_text())
  assert_eq(alignment.metadata.program, "CLUSTAL")
  assert_eq(alignment.metadata.version, "1.81")
  assert_eq(alignment.alignment_length(), 10)
}

///|
test "Bio.Align.clustal parses MUSCLE header" {
  let alignment = align_clustal_test_parse(
    align_clustal_test_sparse_consensus_text(),
  )
  assert_eq(alignment.metadata.program, "MUSCLE")
  assert_eq(alignment.metadata.version, "3.8")
}

///|
test "Bio.Align.clustal parses PROBCONS header" {
  let alignment = align_clustal_test_parse(align_clustal_test_probcons_text())
  assert_eq(alignment.metadata.program, "PROBCONS")
  assert_eq(alignment.metadata.version, "1.12")
  assert_eq(alignment.num_sequences(), 3)
}

///|
test "Bio.Align.clustal parses MSAPROBS header" {
  let alignment = align_clustal_test_parse(align_clustal_test_msaprobs_text())
  assert_eq(alignment.metadata.program, "MSAPROBS")
  assert_eq(alignment.metadata.version, "0.9.7")
}

///|
test "Bio.Align.clustal parses Kalign header" {
  let alignment = align_clustal_test_parse(
    align_clustal_test_no_consensus_text(),
  )
  assert_eq(alignment.metadata.program, "Kalign")
  assert_eq(alignment.metadata.version, "2.0")
}

///|
test "Bio.Align.clustal parses Biopython header" {
  let alignment = align_clustal_test_parse(
    "Biopython 1.86 multiple sequence alignment\n\nalpha ACG\nbeta AC-\n",
  )
  assert_eq(alignment.metadata.program, "Biopython")
  assert_eq(alignment.metadata.version, "1.86")
}

///|
test "Bio.Align.clustal supports a header without version" {
  let alignment = align_clustal_test_parse(
    "CLUSTAL format multiple sequence alignment by PROMALS3D\n\nalpha ACG\nbeta ATG\n",
  )
  assert_eq(alignment.metadata.version, "")
}

///|
test "Bio.Align.clustal validates official cumulative counts" {
  let alignment = align_clustal_test_parse(align_clustal_test_clustalw_text())
  assert_eq(alignment.sequences[0].sequence.length(), 10)
  assert_eq(alignment.sequences[1].sequence.length(), 1)
}

///|
test "Bio.Align.clustal pads sparse consensus blocks" {
  let alignment = align_clustal_test_parse(
    align_clustal_test_sparse_consensus_text(),
  )
  assert_eq(alignment.consensus, Some("    **"))
}

///|
test "Bio.Align.clustal pads trailing consensus spaces" {
  let alignment = align_clustal_test_parse(
    "CLUSTAL 2.0 multiple sequence alignment\n\nalpha ACGT\nbeta ACGT\n      **\n",
  )
  assert_eq(alignment.consensus, Some("**  "))
}

///|
test "Bio.Align.clustal retains an all-space consensus" {
  let alignment = align_clustal_test_parse(
    "CLUSTAL 2.0 multiple sequence alignment\n\nalpha ACGT\nbeta TGCA\n          \n",
  )
  assert_eq(alignment.consensus, Some("    "))
}

///|
test "Bio.Align.clustal distinguishes absent consensus" {
  let alignment = align_clustal_test_parse(
    align_clustal_test_no_consensus_text(),
  )
  assert_eq(alignment.consensus, None)
}

///|
test "Bio.Align.clustal locates exact sequence identifiers" {
  let alignment = align_clustal_test_sample()
  assert_eq(alignment.find_sequence("query_one"), Some(1))
  assert_eq(alignment.find_sequence("missing"), None)
}

///|
test "Bio.Align.clustal returns printed columns" {
  let alignment = align_clustal_test_sample()
  assert_eq(alignment.column(0), Some("MMM"))
  assert_eq(alignment.column(3), Some("-Q-"))
  assert_eq(alignment.column(10), None)
}

///|
test "Bio.Align.clustal maps sequence positions to columns" {
  let alignment = align_clustal_test_sample()
  assert_eq(alignment.sequence_position_to_column(0, 3), Some(5))
  assert_eq(alignment.sequence_position_to_column(1, 2), Some(3))
  assert_eq(alignment.sequence_position_to_column(0, 8), None)
}

///|
test "Bio.Align.clustal maps columns to sequence positions" {
  let alignment = align_clustal_test_sample()
  assert_eq(alignment.column_to_sequence_position(0, 5), Some(3))
  assert_eq(alignment.column_to_sequence_position(0, 3), None)
  assert_eq(alignment.column_to_sequence_position(3, 0), None)
}

///|
test "Bio.Align.clustal projects residues between rows" {
  let alignment = align_clustal_test_sample()
  assert_eq(alignment.map_position(0, 1, 3), Some(4))
  assert_eq(alignment.map_position(0, 1, 1), None)
}

///|
test "Bio.Align.clustal returns aligned coordinate pairs" {
  let pairs = align_clustal_test_sample().aligned_pairs(0, 1) catch {
    AlignClustalError(message) => abort(message)
  }
  assert_eq(pairs.length(), 10)
  assert_eq(pairs[0], (Some(0), Some(0)))
  assert_eq(pairs[1], (Some(1), None))
  assert_eq(pairs[3], (None, Some(2)))
}

///|
test "Bio.Align.clustal compresses coordinate movement paths" {
  let paths = align_clustal_test_sample().coordinate_path()
  assert_eq(paths[0], [0, 1, 2, 3, 3, 7, 8])
  assert_eq(paths[1], [0, 1, 1, 2, 4, 8, 8])
  assert_eq(paths[2], [0, 1, 2, 3, 3, 7, 8])
}

///|
test "Bio.Align.clustal computes pair statistics" {
  let counts = align_clustal_test_sample().pair_counts(0, 1) catch {
    AlignClustalError(message) => abort(message)
  }
  assert_eq(counts.pairs, 1)
  assert_eq(counts.aligned, 6)
  assert_eq(counts.identities, 6)
  assert_eq(counts.mismatches, 0)
  assert_eq(counts.gap_columns, 4)
  assert_eq(counts.double_gap_columns, 0)
  assert_eq(counts.gap_opens, 3)
}

///|
test "Bio.Align.clustal computes exact pair identity" {
  let counts = align_clustal_test_sample().pair_counts(0, 1) catch {
    AlignClustalError(message) => abort(message)
  }
  assert_eq(counts.identity(), 1.0)
}

///|
test "Bio.Align.clustal aggregates all unordered pairs" {
  let counts = align_clustal_test_sample().counts()
  assert_eq(counts.pairs, 3)
  assert_true(counts.aligned > 0)
  assert_true(counts.identities > 0)
}

///|
test "Bio.Align.clustal reports occupancy" {
  let occupancy = align_clustal_test_sample().occupancy()
  assert_eq(occupancy.length(), 10)
  assert_eq(occupancy[0], 1.0)
  assert_eq(occupancy[1], 2.0 / 3.0)
  assert_eq(occupancy[3], 1.0 / 3.0)
}

///|
test "Bio.Align.clustal calculates majority residues" {
  let consensus = align_clustal_test_sample().majority_consensus() catch {
    AlignClustalError(message) => abort(message)
  }
  assert_eq(consensus, "MKTQQAILGH")
}

///|
test "Bio.Align.clustal calculates protein conservation symbols" {
  assert_eq(align_clustal_test_sample().calculated_consensus(), "* *  *:** ")
}

///|
test "Bio.Align.clustal limits nucleotide conservation to identity" {
  let alignment = @src.align_clustal_from_aligned(["alpha", "beta"], [
    "ACGT", "ACAT",
  ]) catch {
    AlignClustalError(message) => abort(message)
  }
  assert_eq(alignment.calculated_consensus(), "** *")
}

///|
test "Bio.Align.clustal reports summary" {
  assert_eq(
    align_clustal_test_sample().summary(),
    "AlignClustalAlignment(program=CLUSTAL, version=2.1, sequences=3, " +
    "columns=10, consensus=yes)",
  )
}

///|
test "Bio.Align.clustal writer emits canonical metadata" {
  let text = align_clustal_test_write(
    align_clustal_test_sample(),
    50,
    36,
    false,
    true,
  )
  assert_true(text.has_prefix("CLUSTAL 2.1 multiple sequence alignment"))
}

///|
test "Bio.Align.clustal writer splits blocks" {
  let text = align_clustal_test_write(
    align_clustal_test_sample(),
    4,
    20,
    false,
    true,
  )
  let prefix = "reference" + " ".repeat(11)
  assert_true(text.contains(prefix + "MKT-"))
  assert_true(text.contains(prefix + "-AIL"))
  assert_true(text.contains(prefix + "GH"))
}

///|
test "Bio.Align.clustal writer optionally emits residue counts" {
  let text = align_clustal_test_write(
    align_clustal_test_sample(),
    7,
    20,
    true,
    true,
  )
  assert_true(text.contains("MKT--AI 5"))
  assert_true(text.contains("LGH 8"))
}

///|
test "Bio.Align.clustal writer omits counts by default" {
  let text = align_clustal_test_write(
    align_clustal_test_sample(),
    10,
    20,
    false,
    true,
  )
  assert_false(text.contains("MKT--AILGH 8"))
}

///|
test "Bio.Align.clustal writer can omit consensus" {
  let text = align_clustal_test_write(
    align_clustal_test_sample(),
    10,
    20,
    false,
    false,
  )
  let parsed = align_clustal_test_parse(text)
  assert_eq(parsed.consensus, None)
}

///|
test "Bio.Align.clustal canonical round trip preserves data" {
  let original = align_clustal_test_sample()
  let text = align_clustal_test_write(original, 4, 20, true, true)
  let parsed = align_clustal_test_parse(text)
  assert_eq(parsed.metadata.program, original.metadata.program)
  assert_eq(parsed.metadata.version, original.metadata.version)
  assert_eq(parsed.sequences, original.sequences)
  assert_eq(parsed.consensus, original.consensus)
}

///|
test "Bio.Align.clustal canonical output is stable" {
  let first = align_clustal_test_write(
    align_clustal_test_sample(),
    4,
    20,
    true,
    true,
  )
  let second = align_clustal_test_write(
    align_clustal_test_parse(first),
    4,
    20,
    true,
    true,
  )
  assert_eq(second, first)
}

///|
test "Bio.Align.clustal constructor creates coordinate-aware rows" {
  let alignment = @src.align_clustal_from_aligned(
    ["alpha", "beta"],
    ["AC-G", "ACTG"],
    program="MUSCLE",
    version="5.1",
    consensus=Some("** *"),
  ) catch {
    AlignClustalError(message) => abort(message)
  }
  assert_eq(alignment.metadata.program, "MUSCLE")
  assert_eq(alignment.sequences[0].sequence, "ACG")
  assert_eq(alignment.consensus, Some("** *"))
}

///|
test "Bio.Align.clustal constructor defensively copies sequence arrays" {
  let rows = [
    @src.AlignClustalSequence::create("alpha", "AC-G") catch {
      AlignClustalError(message) => abort(message)
    },
    @src.AlignClustalSequence::create("beta", "ACTG") catch {
      AlignClustalError(message) => abort(message)
    },
  ]
  let metadata = @src.AlignClustalMetadata::create("Biopython") catch {
    AlignClustalError(message) => abort(message)
  }
  let alignment = @src.AlignClustalAlignment::create(metadata, rows) catch {
    AlignClustalError(message) => abort(message)
  }
  ignore(rows.pop())
  assert_eq(alignment.num_sequences(), 2)
}

///|
test "Bio.Align.clustal rejects empty input" {
  assert_true(align_clustal_test_rejects(""))
}

///|
test "Bio.Align.clustal rejects an unknown generator" {
  assert_true(
    align_clustal_test_rejects(
      "UNKNOWN 1.0 multiple sequence alignment\n\nalpha ACG\nbeta ACG\n",
    ),
  )
}

///|
test "Bio.Align.clustal rejects header-only input" {
  assert_true(
    align_clustal_test_rejects("CLUSTAL 2.0 multiple sequence alignment\n\n"),
  )
}

///|
test "Bio.Align.clustal rejects malformed sequence lines" {
  assert_true(
    align_clustal_test_rejects(
      "CLUSTAL 2.0 multiple sequence alignment\n\nalpha\n",
    ),
  )
}

///|
test "Bio.Align.clustal rejects inconsistent block counts" {
  assert_true(
    align_clustal_test_rejects(
      "CLUSTAL 2.0 multiple sequence alignment\n\nalpha AC\nbeta AC\n\nalpha GT\n",
    ),
  )
}

///|
test "Bio.Align.clustal rejects changed block identifiers" {
  assert_true(
    align_clustal_test_rejects(
      "CLUSTAL 2.0 multiple sequence alignment\n\nalpha AC\nbeta AC\n\nalpha GT\ngamma GT\n",
    ),
  )
}

///|
test "Bio.Align.clustal rejects reordered blocks" {
  assert_true(
    align_clustal_test_rejects(
      "CLUSTAL 2.0 multiple sequence alignment\n\nalpha AC\nbeta AC\n\nbeta GT\nalpha GT\n",
    ),
  )
}

///|
test "Bio.Align.clustal rejects incorrect cumulative counts" {
  assert_true(
    align_clustal_test_rejects(
      "CLUSTAL 2.0 multiple sequence alignment\n\nalpha AC- 3\nbeta ACT 3\n",
    ),
  )
}

///|
test "Bio.Align.clustal rejects nonnumeric residue counts" {
  assert_true(
    align_clustal_test_rejects(
      "CLUSTAL 2.0 multiple sequence alignment\n\nalpha ACG three\nbeta ACG 3\n",
    ),
  )
}

///|
test "Bio.Align.clustal rejects negative residue counts" {
  assert_true(
    align_clustal_test_rejects(
      "CLUSTAL 2.0 multiple sequence alignment\n\nalpha ACG -1\nbeta ACG 3\n",
    ),
  )
}

///|
test "Bio.Align.clustal rejects integer overflow" {
  assert_true(
    align_clustal_test_rejects(
      "CLUSTAL 2.0 multiple sequence alignment\n\nalpha ACG 2147483648\nbeta ACG 3\n",
    ),
  )
}

///|
test "Bio.Align.clustal rejects unequal block row widths" {
  assert_true(
    align_clustal_test_rejects(
      "CLUSTAL 2.0 multiple sequence alignment\n\nalpha ACG\nbeta AC\n",
    ),
  )
}

///|
test "Bio.Align.clustal rejects duplicate identifiers" {
  assert_true(
    align_clustal_test_rejects(
      "CLUSTAL 2.0 multiple sequence alignment\n\nalpha ACG\nalpha ACT\n",
    ),
  )
}

///|
test "Bio.Align.clustal rejects consensus before rows" {
  assert_true(
    align_clustal_test_rejects(
      "CLUSTAL 2.0 multiple sequence alignment\n\n      ***\nalpha ACG\nbeta ACG\n",
    ),
  )
}

///|
test "Bio.Align.clustal rejects multiple consensus lines" {
  assert_true(
    align_clustal_test_rejects(
      "CLUSTAL 2.0 multiple sequence alignment\n\nalpha ACG\nbeta ACG\n      ***\n      ***\n",
    ),
  )
}

///|
test "Bio.Align.clustal rejects rows after consensus" {
  assert_true(
    align_clustal_test_rejects(
      "CLUSTAL 2.0 multiple sequence alignment\n\nalpha ACG\n      ***\nbeta ACG\n",
    ),
  )
}

///|
test "Bio.Align.clustal rejects invalid consensus symbols" {
  assert_true(
    align_clustal_test_rejects(
      "CLUSTAL 2.0 multiple sequence alignment\n\nalpha ACG\nbeta ACG\n      *+*\n",
    ),
  )
}

///|
test "Bio.Align.clustal rejects invalid residue symbols" {
  assert_true(
    align_clustal_test_rejects(
      "CLUSTAL 2.0 multiple sequence alignment\n\nalpha AC1\nbeta ACG\n",
    ),
  )
}

///|
test "Bio.Align.clustal rejects all-gap rows" {
  assert_true(align_clustal_test_sequence_rejects("alpha", "---"))
}

///|
test "Bio.Align.clustal rejects all-gap columns" {
  assert_true(
    align_clustal_test_from_rejects(["alpha", "beta"], ["A-C", "T-G"], None),
  )
}

///|
test "Bio.Align.clustal rejects mismatched constructor dimensions" {
  assert_true(align_clustal_test_from_rejects(["alpha"], ["ACG", "ACT"], None))
}

///|
test "Bio.Align.clustal rejects an empty constructor" {
  assert_true(align_clustal_test_from_rejects([], [], None))
}

///|
test "Bio.Align.clustal rejects constructor consensus length mismatch" {
  assert_true(
    align_clustal_test_from_rejects(
      ["alpha", "beta"],
      ["ACG", "ACT"],
      Some("**"),
    ),
  )
}

///|
test "Bio.Align.clustal rejects constructor consensus symbols" {
  assert_true(
    align_clustal_test_from_rejects(
      ["alpha", "beta"],
      ["ACG", "ACT"],
      Some("*+*"),
    ),
  )
}

///|
test "Bio.Align.clustal rejects identifiers with whitespace" {
  assert_true(align_clustal_test_sequence_rejects("alpha one", "ACG"))
}

///|
test "Bio.Align.clustal rejects unknown metadata programs" {
  assert_true(align_clustal_test_metadata_rejects("PROMALS3D", "1.0"))
}

///|
test "Bio.Align.clustal rejects whitespace in versions" {
  assert_true(align_clustal_test_metadata_rejects("CLUSTAL", "2.0 beta"))
}

///|
test "Bio.Align.clustal rejects invalid writer block width" {
  assert_true(
    align_clustal_test_write_rejects(align_clustal_test_sample(), 0, 36),
  )
}

///|
test "Bio.Align.clustal rejects undersized writer name width" {
  assert_true(
    align_clustal_test_write_rejects(align_clustal_test_sample(), 50, 9),
  )
}

///|
test "Bio.Align.clustal does not silently truncate long identifiers" {
  let alignment = @src.align_clustal_from_aligned(
    ["abcdefghijklmnopqrstuvwxyzABCDE", "beta"],
    ["ACG", "ACT"],
  ) catch {
    AlignClustalError(message) => abort(message)
  }
  assert_true(align_clustal_test_write_rejects(alignment, 50, 36))
}

///|
test "Bio.Align.clustal accepts long identifiers with explicit width" {
  let alignment = @src.align_clustal_from_aligned(
    ["abcdefghijklmnopqrstuvwxyzABCDE", "beta"],
    ["ACG", "ACT"],
  ) catch {
    AlignClustalError(message) => abort(message)
  }
  let text = align_clustal_test_write(alignment, 50, 40, false, true)
  assert_eq(
    align_clustal_test_parse(text).sequences[0].id,
    "abcdefghijklmnopqrstuvwxyzABCDE",
  )
}

///|
test "Bio.Align.clustal rejects invalid majority threshold" {
  let low = try {
    ignore(
      align_clustal_test_sample().majority_consensus(minimum_fraction=-0.1),
    )
    false
  } catch {
    AlignClustalError(_) => true
  }
  let high = try {
    ignore(align_clustal_test_sample().majority_consensus(minimum_fraction=1.1))
    false
  } catch {
    AlignClustalError(_) => true
  }
  assert_true(low)
  assert_true(high)
}

///|
test "Bio.Align.clustal rejects out-of-range pair rows" {
  let rejected = try {
    ignore(align_clustal_test_sample().pair_counts(-1, 0))
    false
  } catch {
    AlignClustalError(_) => true
  }
  assert_true(rejected)
}
