// Black-box tests for Biopython Bio.Align.msf-compatible support.

///|
fn msf_test_parse(text : String) -> @src.MsfAlignment {
  @src.msf_parse(text) catch {
    MsfError(message) => abort("valid MSF input failed: " + message)
  }
}

///|
fn msf_test_parse_unchecked(text : String) -> @src.MsfAlignment {
  @src.msf_parse(text, verify_checksums=false) catch {
    MsfError(message) => abort("unchecked MSF input failed: " + message)
  }
}

///|
fn msf_test_raises(text : String) -> Bool {
  try {
    ignore(@src.msf_parse(text))
    false
  } catch {
    MsfError(_) => true
  }
}

///|
fn msf_test_write_raises(
  alignment : @src.MsfAlignment,
  block_width : Int,
  group_width : Int,
  gap_character : String,
) -> Bool {
  try {
    ignore(
      @src.msf_write(
        alignment,
        block_width=block_width,
        group_width=group_width,
        gap_character=gap_character,
      ),
    )
    false
  } catch {
    MsfError(_) => true
  }
}

///|
fn msf_test_from_aligned_raises(
  ids : Array[String],
  rows : Array[String],
  weights : Array[Double],
) -> Bool {
  try {
    ignore(
      @src.msf_from_aligned(
        ids,
        rows,
        @src.MsfNucleotide,
        weights=weights,
      ),
    )
    false
  } catch {
    MsfError(_) => true
  }
}

///|
fn msf_test_sequence_raises(
  id : String,
  row : String,
  sequence_type : @src.MsfSequenceType,
  weight : Double,
  checksum : Int,
) -> Bool {
  try {
    ignore(
      @src.MsfSequence::create(
        id,
        row,
        sequence_type,
        weight=weight,
        checksum=Some(checksum),
      ),
    )
    false
  } catch {
    MsfError(_) => true
  }
}

///|
fn msf_test_sample() -> @src.MsfAlignment {
  msf_test_parse(@src.msf_example_text())
}

///|
fn msf_test_single_text() -> String {
  "!!NA_MULTIPLE_ALIGNMENT 1.0\n" +
  "\n" +
  "Single MSF: 4 Type: N Check: 748 ..\n" +
  "\n" +
  " Name: alpha Len: 4 Check: 748 Weight: 1.0\n" +
  "//\n" +
  "\n" +
  " alpha ACGT\n"
}

///|
fn msf_test_protein_text() -> String {
  "!!AA_MULTIPLE_ALIGNMENT 1.0\n" +
  "\n" +
  "Protein fixture MSF: 12 Type: P Check: 8342 ..\n" +
  "\n" +
  " Name: alpha Len: 12 Check: 5761 Weight: 1.0\n" +
  " Name: beta Len: 8 Check: 2581 Weight: 0.5\n" +
  "//\n" +
  "\n" +
  " alpha ACDEFG HIKLMN\n" +
  " beta ACDFGH IK\n"
}

///|
fn msf_test_doa_length_mismatch_text() -> String {
  "!!AA_MULTIPLE_ALIGNMENT\n" +
  "\n" +
  "DOA-style MSF: 2 Type: P Check: 0 ..\n" +
  "\n" +
  " Name: full Len: 4 Check: 0 Weight: 1.0\n" +
  " Name: short Len: 2 Check: 0 Weight: 1.0\n" +
  "//\n" +
  "\n" +
  " full ACDE\n" +
  " short AC\n"
}

///|
fn msf_test_w_protein_text() -> String {
  "!!AA_MULTIPLE_ALIGNMENT\n" +
  "   MSF: 99  Type: P  Oct 18, 2017  11:35  Check: 0 ..\n" +
  " Name: W*01:01:01:01    Len:    99  Check: 7236  Weight:  1.00\n" +
  " Name: W*01:01:01:02    Len:    99  Check: 7236  Weight:  1.00\n" +
  " Name: W*01:01:01:03    Len:    99  Check: 7236  Weight:  1.00\n" +
  " Name: W*01:01:01:04    Len:    99  Check: 7236  Weight:  1.00\n" +
  " Name: W*01:01:01:05    Len:    99  Check: 7236  Weight:  1.00\n" +
  " Name: W*01:01:01:06    Len:    99  Check: 7236  Weight:  1.00\n" +
  " Name: W*02:01          Len:    93  Check: 9483  Weight:  1.00\n" +
  " Name: W*03:01:01:01    Len:    93  Check: 9974  Weight:  1.00\n" +
  " Name: W*03:01:01:02    Len:    93  Check: 9974  Weight:  1.00\n" +
  " Name: W*04:01          Len:    93  Check: 9169  Weight:  1.00\n" +
  " Name: W*05:01          Len:    99  Check: 7331  Weight:  1.00\n" +
  "//\n" +
  "\n" +
  " W*01:01:01:01 GLTPFNGYTA ATWTRTAVSS VGMNIPYHGA SYLVRNQELR SWTAADKAAQ\n" +
  " W*01:01:01:02 GLTPFNGYTA ATWTRTAVSS VGMNIPYHGA SYLVRNQELR SWTAADKAAQ\n" +
  " W*01:01:01:03 GLTPFNGYTA ATWTRTAVSS VGMNIPYHGA SYLVRNQELR SWTAADKAAQ\n" +
  " W*01:01:01:04 GLTPFNGYTA ATWTRTAVSS VGMNIPYHGA SYLVRNQELR SWTAADKAAQ\n" +
  " W*01:01:01:05 GLTPFNGYTA ATWTRTAVSS VGMNIPYHGA SYLVRNQELR SWTAADKAAQ\n" +
  " W*01:01:01:06 GLTPFNGYTA ATWTRTAVSS VGMNIPYHGA SYLVRNQELR SWTAADKAAQ\n" +
  " W*02:01 GLTPSNGYTA ATWTRTAASS VGMNIPYDGA SYLVRNQELR SWTAADKAAQ\n" +
  " W*03:01:01:01 GLTPSSGYTA ATWTRTAVSS VGMNIPYHGA SYLVRNQELR SWTAADKAAQ\n" +
  " W*03:01:01:02 GLTPSSGYTA ATWTRTAVSS VGMNIPYHGA SYLVRNQELR SWTAADKAAQ\n" +
  " W*04:01 GLTPSNGYTA ATWTRTAASS VGMNIPYDGA SYLVRNQELR SWTAADKAAQ\n" +
  " W*05:01 GLTPSSGYTA ATWTRTAVSS VGMNIPYHGA SYLVRNQELR SWTAADKAAQ\n" +
  " W*01:01:01:01 MPWRRNRQSC SKPTCREGGR SGSAKSLRMG RRGCSAQNPK DSHDPPPHL\n" +
  " W*01:01:01:02 MPWRRNRQSC SKPTCREGGR SGSAKSLRMG RRGCSAQNPK DSHDPPPHL\n" +
  " W*01:01:01:03 MPWRRNRQSC SKPTCREGGR SGSAKSLRMG RRGCSAQNPK DSHDPPPHL\n" +
  " W*01:01:01:04 MPWRRNRQSC SKPTCREGGR SGSAKSLRMG RRGCSAQNPK DSHDPPPHL\n" +
  " W*01:01:01:05 MPWRRNRQSC SKPTCREGGR SGSAKSLRMG RRGCSAQNPK DSHDPPPHL\n" +
  " W*01:01:01:06 MPWRRNRQSC SKPTCREGGR SGSAKSLRMG RRGCSAQNPK DSHDPPPHL\n" +
  " W*02:01 MPWRRNMQSC SKPTCREGGR SGSAKSLRMG RRRCTAQNPK RLT\n" +
  " W*03:01:01:01 MPWRRNRQSC SKPTCREGGR SGSAKSLRMG RRGCSAQNPK RLT\n" +
  " W*03:01:01:02 MPWRRNRQSC SKPTCREGGR SGSAKSLRMG RRGCSAQNPK RLT\n" +
  " W*04:01 MPWRRNMQSC SKPTCREGGR SGSAKSLRMG RRGCSAQNPK RLT\n" +
  " W*05:01 MPWRRNRQSC SKPTCREGGR SGSAKSLRMG RRGCSAQNPK DSHDPPPHL\n"
}

///|
test "Bio.Align.msf parses sample metadata" {
  let alignment = msf_test_sample()
  assert_eq(alignment.metadata.header, "!!NA_MULTIPLE_ALIGNMENT 1.0")
  assert_eq(alignment.metadata.title, "MoonBit example")
  assert_eq(alignment.metadata.declared_length, 12)
  assert_eq(alignment.metadata.sequence_type, @src.MsfNucleotide)
  assert_eq(alignment.metadata.checksum_label, "Check:")
}

///|
test "Bio.Align.msf parses interleaved sample rows" {
  let alignment = msf_test_sample()
  assert_eq(alignment.num_sequences(), 3)
  assert_eq(alignment.alignment_length(), 12)
  assert_eq(alignment.sequences[0].aligned_sequence, "ACGTACGTACGT")
  assert_eq(alignment.sequences[1].aligned_sequence, "ACG-ACGTAC-T")
  assert_eq(alignment.sequences[2].aligned_sequence, "A-GTTCGTACGT")
}

///|
test "Bio.Align.msf preserves ungapped sequences" {
  let alignment = msf_test_sample()
  assert_eq(alignment.sequences[0].sequence, "ACGTACGTACGT")
  assert_eq(alignment.sequences[1].sequence, "ACGACGTACT")
  assert_eq(alignment.sequences[2].sequence, "AGTTCGTACGT")
}

///|
test "Bio.Align.msf preserves lengths and weights" {
  let alignment = msf_test_sample()
  assert_eq(alignment.sequences[0].length, 12)
  assert_eq(alignment.sequences[1].length, 10)
  assert_eq(alignment.sequences[2].length, 11)
  assert_eq(alignment.sequences[1].weight, 0.75)
  assert_eq(alignment.sequences[2].weight, 1.25)
}

///|
test "Bio.Align.msf validates sample checksums" {
  let alignment = msf_test_sample()
  assert_true(alignment.sequences[0].checksum_valid())
  assert_true(alignment.sequences[1].checksum_valid())
  assert_true(alignment.sequences[2].checksum_valid())
  assert_true(alignment.checksums_valid())
  assert_eq(alignment.computed_checksum(), 4573)
}

///|
test "Bio.Align.msf reports summary" {
  assert_eq(
    msf_test_sample().summary(),
    "MsfAlignment(type=N, sequences=3, columns=12, declared=12, checksum=4573)",
  )
}

///|
test "Bio.Align.msf finds rows by identifier" {
  let alignment = msf_test_sample()
  assert_eq(alignment.find_sequence("reference"), Some(0))
  assert_eq(alignment.find_sequence("query_two"), Some(2))
  assert_eq(alignment.find_sequence("missing"), None)
}

///|
test "Bio.Align.msf returns alignment columns" {
  let alignment = msf_test_sample()
  assert_eq(alignment.column(0), Some("AAA"))
  assert_eq(alignment.column(1), Some("CC-"))
  assert_eq(alignment.column(3), Some("T-T"))
  assert_eq(alignment.column(-1), None)
  assert_eq(alignment.column(12), None)
}

///|
test "Bio.Align.msf maps sequence positions to columns" {
  let alignment = msf_test_sample()
  assert_eq(alignment.sequence_position_to_column(1, 0), Some(0))
  assert_eq(alignment.sequence_position_to_column(1, 2), Some(2))
  assert_eq(alignment.sequence_position_to_column(1, 3), Some(4))
  assert_eq(alignment.sequence_position_to_column(1, 9), Some(11))
  assert_eq(alignment.sequence_position_to_column(1, 10), None)
}

///|
test "Bio.Align.msf maps columns to sequence positions" {
  let alignment = msf_test_sample()
  assert_eq(alignment.column_to_sequence_position(1, 2), Some(2))
  assert_eq(alignment.column_to_sequence_position(1, 3), None)
  assert_eq(alignment.column_to_sequence_position(1, 4), Some(3))
  assert_eq(alignment.column_to_sequence_position(1, 11), Some(9))
}

///|
test "Bio.Align.msf maps residues between rows" {
  let alignment = msf_test_sample()
  assert_eq(alignment.map_position(0, 1, 2), Some(2))
  assert_eq(alignment.map_position(0, 1, 3), None)
  assert_eq(alignment.map_position(0, 1, 4), Some(3))
  assert_eq(alignment.map_position(1, 0, 3), Some(4))
  assert_eq(alignment.map_position(0, 3, 0), None)
}

///|
test "Bio.Align.msf emits aligned coordinate pairs" {
  let pairs = msf_test_sample().aligned_pairs(0, 1) catch {
    MsfError(message) => abort(message)
  }
  assert_eq(pairs.length(), 12)
  assert_eq(pairs[0], (Some(0), Some(0)))
  assert_eq(pairs[3], (Some(3), None))
  assert_eq(pairs[4], (Some(4), Some(3)))
  assert_eq(pairs[10], (Some(10), None))
  assert_eq(pairs[11], (Some(11), Some(9)))
}

///|
test "Bio.Align.msf builds compact coordinate path" {
  let path = msf_test_sample().coordinate_path()
  assert_eq(path.length(), 3)
  assert_eq(path[0], [0, 1, 2, 3, 4, 10, 11, 12])
  assert_eq(path[1], [0, 1, 2, 3, 3, 9, 9, 10])
  assert_eq(path[2], [0, 1, 1, 2, 3, 9, 10, 11])
}

///|
test "Bio.Align.msf counts identity and gaps" {
  let counts = msf_test_sample().pair_counts(0, 1) catch {
    MsfError(message) => abort(message)
  }
  assert_eq(counts.columns, 12)
  assert_eq(counts.aligned, 10)
  assert_eq(counts.identities, 10)
  assert_eq(counts.mismatches, 0)
  assert_eq(counts.gap_columns, 2)
  assert_eq(counts.double_gap_columns, 0)
  assert_eq(counts.gap_opens, 2)
  assert_eq(counts.identity(), 1.0)
}

///|
test "Bio.Align.msf counts mismatches" {
  let counts = msf_test_sample().pair_counts(0, 2) catch {
    MsfError(message) => abort(message)
  }
  assert_eq(counts.aligned, 11)
  assert_eq(counts.identities, 10)
  assert_eq(counts.mismatches, 1)
  assert_eq(counts.gap_columns, 1)
  assert_eq(counts.gap_opens, 1)
}

///|
test "Bio.Align.msf counts double-gap columns" {
  let alignment = @src.msf_from_aligned(
    ["reference", "first", "second"],
    ["ACGT", "A--T", "A--T"],
    @src.MsfNucleotide,
  ) catch {
    MsfError(message) => abort(message)
  }
  let counts = alignment.pair_counts(1, 2) catch {
    MsfError(message) => abort(message)
  }
  assert_eq(counts.aligned, 2)
  assert_eq(counts.double_gap_columns, 2)
  assert_eq(counts.gap_columns, 0)
}

///|
test "Bio.Align.msf calculates majority consensus" {
  let alignment = msf_test_sample()
  assert_eq(
    alignment.consensus() catch {
      MsfError(message) => abort(message)
    },
    "ACGTACGTACGT",
  )
}

///|
test "Bio.Align.msf applies consensus threshold" {
  let consensus = msf_test_sample().consensus(minimum_fraction=1.0) catch {
    MsfError(message) => abort(message)
  }
  assert_eq(consensus, "ACGTXCGTACGT")
}

///|
test "Bio.Align.msf calculates column occupancy" {
  let occupancy = msf_test_sample().occupancy()
  assert_eq(occupancy.length(), 12)
  assert_eq(occupancy[0], 1.0)
  assert_eq(occupancy[1], 2.0 / 3.0)
  assert_eq(occupancy[3], 2.0 / 3.0)
  assert_eq(occupancy[4], 1.0)
  assert_eq(occupancy[10], 2.0 / 3.0)
}

///|
test "Bio.Align.msf implements standard GCG checksum" {
  assert_eq(@src.msf_gcg_checksum("ACGTACGTACGT"), 5688)
  assert_eq(@src.msf_gcg_checksum("ACDEFGHIKLMN"), 5761)
  assert_eq(@src.msf_gcg_checksum("ACDFGHIK"), 2581)
}

///|
test "Bio.Align.msf checksum is case insensitive" {
  assert_eq(@src.msf_gcg_checksum("acgt"), 748)
  assert_eq(@src.msf_gcg_checksum("AcGt"), 748)
}

///|
test "Bio.Align.msf checksum cycles position after 57" {
  assert_eq(
    @src.msf_gcg_checksum(
      "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA",
    ),
    7510,
  )
}

///|
test "Bio.Align.msf parses official W protein shape" {
  let alignment = msf_test_parse(msf_test_w_protein_text())
  assert_eq(alignment.metadata.sequence_type, @src.MsfProtein)
  assert_eq(alignment.num_sequences(), 11)
  assert_eq(alignment.alignment_length(), 99)
  assert_false(alignment.length_mismatch())
}

///|
test "Bio.Align.msf parses official W identifiers" {
  let alignment = msf_test_parse(msf_test_w_protein_text())
  assert_eq(alignment.sequences[0].id, "W*01:01:01:01")
  assert_eq(alignment.sequences[6].id, "W*02:01")
  assert_eq(alignment.sequences[10].id, "W*05:01")
}

///|
test "Bio.Align.msf parses official W sequence" {
  let alignment = msf_test_parse(msf_test_w_protein_text())
  assert_eq(
    alignment.sequences[0].sequence,
    "GLTPFNGYTAATWTRTAVSSVGMNIPYHGASYLVRNQELRSWTAADKAAQMPWRRNRQSCSKPTCREGGRSGSAKSLRMGRRGCSAQNPKDSHDPPPHL",
  )
  assert_eq(alignment.sequences[0].checksum, 7236)
}

///|
test "Bio.Align.msf pads official W short rows" {
  let alignment = msf_test_parse(msf_test_w_protein_text())
  assert_eq(alignment.sequences[6].length, 93)
  assert_eq(alignment.sequences[6].aligned_sequence.length(), 99)
  assert_true(alignment.sequences[6].aligned_sequence.ends_with("RLT------"))
  assert_eq(alignment.sequences[6].checksum, 9483)
}

///|
test "Bio.Align.msf builds official W coordinate path" {
  let path = msf_test_parse(msf_test_w_protein_text()).coordinate_path()
  assert_eq(path[0], [0, 93, 99])
  assert_eq(path[6], [0, 93, 93])
  assert_eq(path[10], [0, 93, 99])
}

///|
test "Bio.Align.msf preserves DOA-style declared length mismatch" {
  let alignment = msf_test_parse(msf_test_doa_length_mismatch_text())
  assert_eq(alignment.metadata.declared_length, 2)
  assert_eq(alignment.alignment_length(), 4)
  assert_true(alignment.length_mismatch())
}

///|
test "Bio.Align.msf pads DOA-style completed rows" {
  let alignment = msf_test_parse(msf_test_doa_length_mismatch_text())
  assert_eq(alignment.sequences[1].sequence, "AC")
  assert_eq(alignment.sequences[1].aligned_sequence, "AC--")
  assert_eq(alignment.coordinate_path()[1], [0, 2, 2])
}

///|
test "Bio.Align.msf parses PileUp header" {
  let text = "PileUp\n\n" +
    "PileUp sample MSF: 4 Type: N Check: 0 ..\n\n" +
    " Name: alpha oo Len: 4 Check: 0 Weight: 1.0\n" +
    "//\n\nalpha ACGT\n"
  let alignment = msf_test_parse(text)
  assert_eq(alignment.metadata.header, "PileUp")
  assert_eq(alignment.metadata.sequence_type, @src.MsfNucleotide)
}

///|
test "Bio.Align.msf parses EMBOSS CompCheck header" {
  let text = "!!NA_MULTIPLE_ALIGNMENT 1.0\n\n" +
    "stdout MSF: 4 Type: N 01/08/19 CompCheck: 748 ..\n\n" +
    " Name: alpha Len: 4 Check: 748 Weight: 1.0\n" +
    "//\n\nalpha ACGT\n"
  let alignment = msf_test_parse(text)
  assert_eq(alignment.metadata.checksum_label, "CompCheck:")
  assert_eq(alignment.metadata.date_text, "01/08/19")
}

///|
test "Bio.Align.msf preserves title date and preamble" {
  let text = "PileUp\n" +
    "Generated by MoonBit\n" +
    "Project alpha MSF: 4 Type: N Jan 2 2026 Check: 0 ..\n\n" +
    " Name: alpha Len: 4 Check: 0 Weight: 1.0\n" +
    "//\n\nalpha ACGT\n"
  let alignment = msf_test_parse(text)
  assert_eq(alignment.metadata.preamble, ["Generated by MoonBit"])
  assert_eq(alignment.metadata.title, "Project alpha")
  assert_eq(alignment.metadata.date_text, "Jan 2 2026")
}

///|
test "Bio.Align.msf accepts CRLF and lowercase residues" {
  let text = "!!NA_MULTIPLE_ALIGNMENT\r\n\r\n" +
    "MSF: 4 Type: N Check: 0 ..\r\n\r\n" +
    "Name: alpha Len: 4 Check: 0 Weight: 1.0\r\n" +
    "//\r\n\r\nalpha acgt\r\n"
  let alignment = msf_test_parse(text)
  assert_eq(alignment.sequences[0].sequence, "ACGT")
}

///|
test "Bio.Align.msf normalizes dot tilde and hyphen gaps" {
  let text = "!!NA_MULTIPLE_ALIGNMENT\n" +
    "MSF: 4 Type: N Check: 0 ..\n" +
    "Name: dot Len: 3 Check: 0 Weight: 1.0\n" +
    "Name: tilde Len: 3 Check: 0 Weight: 1.0\n" +
    "Name: hyphen Len: 3 Check: 0 Weight: 1.0\n" +
    "//\n\n" +
    "dot A.CG\n" +
    "tilde AT~G\n" +
    "hyphen A-CG\n"
  let alignment = msf_test_parse(text)
  assert_eq(alignment.sequences[0].aligned_sequence, "A-CG")
  assert_eq(alignment.sequences[1].aligned_sequence, "AT-G")
  assert_eq(alignment.sequences[2].aligned_sequence, "A-CG")
}

///|
test "Bio.Align.msf ignores numeric coordinate lines" {
  let text = "!!NA_MULTIPLE_ALIGNMENT\n" +
    "MSF: 4 Type: N Check: 0 ..\n" +
    "Name: alpha Len: 4 Check: 0 Weight: 1.0\n" +
    "//\n\n" +
    "1 2\nalpha AC\n\n3 4\nalpha GT\n"
  assert_eq(msf_test_parse(text).sequences[0].sequence, "ACGT")
}

///|
test "Bio.Align.msf writes canonical default format" {
  let output = @src.msf_write(msf_test_sample()) catch {
    MsfError(message) => abort(message)
  }
  assert_true(output.contains("MSF: 12 Type: N"))
  assert_true(output.contains("Check: 4573 .."))
  assert_true(output.contains("Name: reference Len: 12 Check: 5688"))
  assert_true(output.contains("query_one  ACG.ACGTAC .T"))
}

///|
test "Bio.Align.msf writes tilde gaps" {
  let output = @src.msf_write(
    msf_test_sample(),
    block_width=12,
    group_width=4,
    gap_character="~",
  ) catch {
    MsfError(message) => abort(message)
  }
  assert_true(output.contains("ACG~ ACGT AC~T"))
}

///|
test "Bio.Align.msf writes hyphen gaps" {
  let output = @src.msf_write(
    msf_test_sample(),
    block_width=12,
    group_width=6,
    gap_character="-",
  ) catch {
    MsfError(message) => abort(message)
  }
  assert_true(output.contains("ACG-AC GTAC-T"))
}

///|
test "Bio.Align.msf writes custom interleaved blocks" {
  let output = @src.msf_write(
    msf_test_sample(),
    block_width=6,
    group_width=3,
  ) catch {
    MsfError(message) => abort(message)
  }
  assert_true(output.contains("reference  ACG TAC"))
  assert_true(output.contains("reference  GTA CGT"))
}

///|
test "Bio.Align.msf round trips canonical output" {
  let original = msf_test_sample()
  let output = @src.msf_write(original) catch {
    MsfError(message) => abort(message)
  }
  let parsed = msf_test_parse(output)
  assert_eq(parsed.num_sequences(), original.num_sequences())
  assert_eq(parsed.alignment_length(), original.alignment_length())
  for index = 0; index < parsed.num_sequences(); index = index + 1 {
    assert_eq(
      parsed.sequences[index].aligned_sequence,
      original.sequences[index].aligned_sequence,
    )
  }
}

///|
test "Bio.Align.msf round trip preserves weights" {
  let output = @src.msf_write(msf_test_sample()) catch {
    MsfError(message) => abort(message)
  }
  let parsed = msf_test_parse(output)
  assert_eq(parsed.sequences[0].weight, 1.0)
  assert_eq(parsed.sequences[1].weight, 0.75)
  assert_eq(parsed.sequences[2].weight, 1.25)
}

///|
test "Bio.Align.msf round trip recomputes checksums" {
  let output = @src.msf_write(msf_test_sample()) catch {
    MsfError(message) => abort(message)
  }
  let parsed = msf_test_parse(output)
  assert_eq(parsed.metadata.declared_checksum, 4573)
  assert_true(parsed.checksums_valid())
}

///|
test "Bio.Align.msf constructs alignment from rows" {
  let alignment = @src.msf_from_aligned(
    ["alpha", "beta"],
    ["ACGT", "A-GT"],
    @src.MsfNucleotide,
    weights=[1.0, 0.5],
  ) catch {
    MsfError(message) => abort(message)
  }
  assert_eq(alignment.num_sequences(), 2)
  assert_eq(alignment.sequences[1].sequence, "AGT")
  assert_eq(alignment.sequences[1].weight, 0.5)
  assert_eq(alignment.metadata.declared_length, 4)
}

///|
test "Bio.Align.msf sequence type emits canonical code" {
  assert_eq(@src.MsfProtein.code(), "P")
  assert_eq(@src.MsfNucleotide.code(), "N")
}

///|
test "Bio.Align.msf constructs normalized sequence" {
  let sequence = @src.MsfSequence::create(
    "alpha",
    "acg.t",
    @src.MsfNucleotide,
    weight=0.25,
  ) catch {
    MsfError(message) => abort(message)
  }
  assert_eq(sequence.aligned_sequence, "ACG-T")
  assert_eq(sequence.sequence, "ACGT")
  assert_eq(sequence.length, 4)
  assert_eq(sequence.checksum, 748)
  assert_eq(sequence.weight, 0.25)
}

///|
test "Bio.Align.msf constructs canonical metadata" {
  let metadata = @src.MsfMetadata::create(4, @src.MsfProtein) catch {
    MsfError(message) => abort(message)
  }
  assert_eq(metadata.header, "!!AA_MULTIPLE_ALIGNMENT 1.0")
  assert_eq(metadata.declared_length, 4)
  assert_eq(metadata.sequence_type, @src.MsfProtein)
}

///|
test "Bio.Align.msf accepts protein wildcard residues" {
  let sequence = @src.MsfSequence::create(
    "protein",
    "ACDX*?",
    @src.MsfProtein,
  ) catch {
    MsfError(message) => abort(message)
  }
  assert_eq(sequence.sequence, "ACDX*?")
}

///|
test "Bio.Align.msf accepts nucleotide ambiguity codes" {
  let sequence = @src.MsfSequence::create(
    "dna",
    "ACGTURYSWKMBDHVNX",
    @src.MsfNucleotide,
  ) catch {
    MsfError(message) => abort(message)
  }
  assert_eq(sequence.length, 17)
}

///|
test "Bio.Align.msf treats zero checksums as unspecified" {
  let alignment = msf_test_parse(msf_test_doa_length_mismatch_text())
  assert_true(alignment.sequences[0].checksum_valid())
  assert_true(alignment.checksums_valid())
}

///|
test "Bio.Align.msf rejects empty input" {
  assert_true(msf_test_raises(""))
}

///|
test "Bio.Align.msf rejects leading blank input" {
  assert_true(msf_test_raises("\n" + msf_test_single_text()))
}

///|
test "Bio.Align.msf rejects unknown header" {
  assert_true(
    msf_test_raises(
      "CLUSTAL\nMSF: 4 Type: N Check: 0 ..\n" +
      "Name: alpha Len: 4 Check: 0 Weight: 1.0\n//\n\nalpha ACGT\n",
    ),
  )
}

///|
test "Bio.Align.msf rejects missing alignment header" {
  assert_true(msf_test_raises("!!NA_MULTIPLE_ALIGNMENT\nName: alpha\n"))
}

///|
test "Bio.Align.msf rejects malformed alignment header" {
  assert_true(
    msf_test_raises(
      "!!NA_MULTIPLE_ALIGNMENT\nMSF: 4 Kind: N Check: 0 ..\n",
    ),
  )
}

///|
test "Bio.Align.msf rejects invalid sequence type" {
  assert_true(
    msf_test_raises(
      "PileUp\nMSF: 4 Type: X Check: 0 ..\n" +
      "Name: alpha Len: 4 Check: 0 Weight: 1.0\n//\n\nalpha ACGT\n",
    ),
  )
}

///|
test "Bio.Align.msf rejects header type conflict" {
  assert_true(
    msf_test_raises(
      "!!AA_MULTIPLE_ALIGNMENT\nMSF: 4 Type: N Check: 0 ..\n" +
      "Name: alpha Len: 4 Check: 0 Weight: 1.0\n//\n\nalpha ACGT\n",
    ),
  )
}

///|
test "Bio.Align.msf rejects zero declared width" {
  assert_true(
    msf_test_raises(
      "PileUp\nMSF: 0 Type: N Check: 0 ..\n" +
      "Name: alpha Len: 4 Check: 0 Weight: 1.0\n//\n\nalpha ACGT\n",
    ),
  )
}

///|
test "Bio.Align.msf rejects integer overflow" {
  assert_true(
    msf_test_raises(
      "PileUp\nMSF: 2147483648 Type: N Check: 0 ..\n",
    ),
  )
}

///|
test "Bio.Align.msf rejects file checksum outside range" {
  assert_true(
    msf_test_raises(
      "PileUp\nMSF: 4 Type: N Check: 10000 ..\n",
    ),
  )
}

///|
test "Bio.Align.msf rejects header without names" {
  assert_true(
    msf_test_raises("PileUp\nMSF: 4 Type: N Check: 0 ..\n//\n\nalpha ACGT\n"),
  )
}

///|
test "Bio.Align.msf rejects duplicate names" {
  assert_true(
    msf_test_raises(
      "PileUp\nMSF: 4 Type: N Check: 0 ..\n" +
      "Name: alpha Len: 4 Check: 0 Weight: 1.0\n" +
      "Name: alpha Len: 4 Check: 0 Weight: 1.0\n" +
      "//\n\nalpha ACGT\n",
    ),
  )
}

///|
test "Bio.Align.msf rejects missing header terminator" {
  assert_true(
    msf_test_raises(
      "PileUp\nMSF: 4 Type: N Check: 0 ..\n" +
      "Name: alpha Len: 4 Check: 0 Weight: 1.0\n",
    ),
  )
}

///|
test "Bio.Align.msf rejects missing blank after terminator" {
  assert_true(
    msf_test_raises(
      "PileUp\nMSF: 4 Type: N Check: 0 ..\n" +
      "Name: alpha Len: 4 Check: 0 Weight: 1.0\n//\nalpha ACGT\n",
    ),
  )
}

///|
test "Bio.Align.msf rejects empty sequence body" {
  assert_true(
    msf_test_raises(
      "PileUp\nMSF: 4 Type: N Check: 0 ..\n" +
      "Name: alpha Len: 4 Check: 0 Weight: 1.0\n//\n\n",
    ),
  )
}

///|
test "Bio.Align.msf rejects unknown body row" {
  assert_true(
    msf_test_raises(
      "PileUp\nMSF: 4 Type: N Check: 0 ..\n" +
      "Name: alpha Len: 4 Check: 0 Weight: 1.0\n//\n\nbeta ACGT\n",
    ),
  )
}

///|
test "Bio.Align.msf rejects row without residues" {
  assert_true(
    msf_test_raises(
      "PileUp\nMSF: 4 Type: N Check: 0 ..\n" +
      "Name: alpha Len: 4 Check: 0 Weight: 1.0\n//\n\nalpha\n",
    ),
  )
}

///|
test "Bio.Align.msf rejects excess residues" {
  assert_true(
    msf_test_raises(
      "PileUp\nMSF: 4 Type: N Check: 0 ..\n" +
      "Name: alpha Len: 3 Check: 0 Weight: 1.0\n//\n\nalpha ACGT\n",
    ),
  )
}

///|
test "Bio.Align.msf rejects missing residues" {
  assert_true(
    msf_test_raises(
      "PileUp\nMSF: 4 Type: N Check: 0 ..\n" +
      "Name: alpha Len: 4 Check: 0 Weight: 1.0\n//\n\nalpha ACG\n",
    ),
  )
}

///|
test "Bio.Align.msf rejects invalid nucleotide residue" {
  assert_true(
    msf_test_raises(
      "PileUp\nMSF: 4 Type: N Check: 0 ..\n" +
      "Name: alpha Len: 4 Check: 0 Weight: 1.0\n//\n\nalpha ACGZ\n",
    ),
  )
}

///|
test "Bio.Align.msf rejects invalid protein residue" {
  assert_true(
    msf_test_raises(
      "PileUp\nMSF: 4 Type: P Check: 0 ..\n" +
      "Name: alpha Len: 4 Check: 0 Weight: 1.0\n//\n\nalpha ACD%\n",
    ),
  )
}

///|
test "Bio.Align.msf rejects all-gap alignment column" {
  assert_true(
    msf_test_raises(
      "PileUp\nMSF: 3 Type: N Check: 0 ..\n" +
      "Name: alpha Len: 2 Check: 0 Weight: 1.0\n" +
      "Name: beta Len: 2 Check: 0 Weight: 1.0\n" +
      "//\n\nalpha A.C\nbeta A-C\n",
    ),
  )
}

///|
test "Bio.Align.msf verifies sequence checksum by default" {
  let text = "PileUp\nMSF: 4 Type: N Check: 0 ..\n" +
    "Name: alpha Len: 4 Check: 1 Weight: 1.0\n" +
    "//\n\nalpha ACGT\n"
  assert_true(msf_test_raises(text))
}

///|
test "Bio.Align.msf can disable checksum verification" {
  let text = "PileUp\nMSF: 4 Type: N Check: 1 ..\n" +
    "Name: alpha Len: 4 Check: 1 Weight: 1.0\n" +
    "//\n\nalpha ACGT\n"
  let alignment = msf_test_parse_unchecked(text)
  assert_false(alignment.checksums_valid())
  assert_eq(alignment.sequences[0].sequence, "ACGT")
}

///|
test "Bio.Align.msf verifies file checksum by default" {
  let text = "PileUp\nMSF: 4 Type: N Check: 1 ..\n" +
    "Name: alpha Len: 4 Check: 748 Weight: 1.0\n" +
    "//\n\nalpha ACGT\n"
  assert_true(msf_test_raises(text))
}

///|
test "Bio.Align.msf rejects malformed weight" {
  assert_true(
    msf_test_raises(
      "PileUp\nMSF: 4 Type: N Check: 0 ..\n" +
      "Name: alpha Len: 4 Check: 0 Weight: heavy\n//\n\nalpha ACGT\n",
    ),
  )
}

///|
test "Bio.Align.msf rejects negative weight" {
  assert_true(
    msf_test_raises(
      "PileUp\nMSF: 4 Type: N Check: 0 ..\n" +
      "Name: alpha Len: 4 Check: 0 Weight: -1.0\n//\n\nalpha ACGT\n",
    ),
  )
}

///|
test "Bio.Align.msf rejects malformed name descriptor" {
  assert_true(
    msf_test_raises(
      "PileUp\nMSF: 4 Type: N Check: 0 ..\n" +
      "Name: alpha Length: 4 Check: 0 Weight: 1.0\n//\n\nalpha ACGT\n",
    ),
  )
}

///|
test "Bio.Align.msf rejects unexpected descriptor token" {
  assert_true(
    msf_test_raises(
      "PileUp\nMSF: 4 Type: N Check: 0 ..\n" +
      "Name: alpha extra Len: 4 Check: 0 Weight: 1.0\n//\n\nalpha ACGT\n",
    ),
  )
}

///|
test "Bio.Align.msf rejects empty row constructor input" {
  assert_true(
    msf_test_from_aligned_raises([], [], []),
  )
}

///|
test "Bio.Align.msf rejects mismatched row constructor input" {
  assert_true(
    msf_test_from_aligned_raises(["alpha"], ["ACGT", "ACGT"], []),
  )
}

///|
test "Bio.Align.msf rejects mismatched weights" {
  assert_true(
    msf_test_from_aligned_raises(
      ["alpha", "beta"],
      ["ACGT", "ACGT"],
      [1.0],
    ),
  )
}

///|
test "Bio.Align.msf rejects unequal row widths" {
  assert_true(
    msf_test_from_aligned_raises(
      ["alpha", "beta"],
      ["ACGT", "ACG"],
      [],
    ),
  )
}

///|
test "Bio.Align.msf rejects duplicate constructor identifiers" {
  assert_true(
    msf_test_from_aligned_raises(
      ["alpha", "alpha"],
      ["ACGT", "ACGT"],
      [],
    ),
  )
}

///|
test "Bio.Align.msf rejects whitespace in sequence identifier" {
  assert_true(
    msf_test_sequence_raises(
      "bad id",
      "ACGT",
      @src.MsfNucleotide,
      1.0,
      748,
    ),
  )
}

///|
test "Bio.Align.msf rejects empty aligned sequence" {
  assert_true(
    msf_test_sequence_raises("", "", @src.MsfNucleotide, 1.0, 0),
  )
}

///|
test "Bio.Align.msf rejects all-gap sequence" {
  assert_true(
    msf_test_sequence_raises(
      "alpha",
      "...",
      @src.MsfNucleotide,
      1.0,
      0,
    ),
  )
}

///|
test "Bio.Align.msf rejects sequence checksum outside range" {
  assert_true(
    msf_test_sequence_raises(
      "alpha",
      "ACGT",
      @src.MsfNucleotide,
      1.0,
      10000,
    ),
  )
}

///|
test "Bio.Align.msf rejects negative sequence weight" {
  assert_true(
    msf_test_sequence_raises(
      "alpha",
      "ACGT",
      @src.MsfNucleotide,
      -1.0,
      748,
    ),
  )
}

///|
test "Bio.Align.msf rejects invalid metadata width" {
  let failed = try {
    ignore(@src.MsfMetadata::create(0, @src.MsfNucleotide))
    false
  } catch {
    MsfError(_) => true
  }
  assert_true(failed)
}

///|
test "Bio.Align.msf rejects invalid metadata checksum label" {
  let failed = try {
    ignore(
      @src.MsfMetadata::create(
        4,
        @src.MsfNucleotide,
        checksum_label="Checksum:",
      ),
    )
    false
  } catch {
    MsfError(_) => true
  }
  assert_true(failed)
}

///|
test "Bio.Align.msf rejects invalid writer block width" {
  assert_true(msf_test_write_raises(msf_test_sample(), 0, 1, "."))
}

///|
test "Bio.Align.msf rejects indivisible writer groups" {
  assert_true(msf_test_write_raises(msf_test_sample(), 10, 3, "."))
}

///|
test "Bio.Align.msf rejects invalid writer gap character" {
  assert_true(msf_test_write_raises(msf_test_sample(), 10, 5, "_"))
}

///|
test "Bio.Align.msf rejects invalid consensus threshold" {
  let failed = try {
    ignore(msf_test_sample().consensus(minimum_fraction=1.1))
    false
  } catch {
    MsfError(_) => true
  }
  assert_true(failed)
}

///|
test "Bio.Align.msf rejects invalid pairwise row" {
  let failed = try {
    ignore(msf_test_sample().pair_counts(0, 3))
    false
  } catch {
    MsfError(_) => true
  }
  assert_true(failed)
}
