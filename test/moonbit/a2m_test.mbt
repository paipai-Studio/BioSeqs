// Black-box tests for Biopython Bio.Align.a2m-compatible support.

///|
fn a2m_sample() -> @src.A2mAlignment {
  @src.a2m_parse(@src.a2m_example_text()) catch {
    _ => abort("valid A2M example should parse")
  }
}

///|
fn a2m_raises(text : String) -> Bool {
  try {
    ignore(@src.a2m_parse(text))
    false
  } catch {
    A2mError(_) => true
  }
}

///|
test "a2m parses sequence identifiers" {
  let alignment = a2m_sample()
  assert_eq(alignment.num_sequences(), 3)
  assert_eq(alignment.sequences[0].id, "reference")
  assert_eq(alignment.sequences[1].id, "query_one")
  assert_eq(alignment.sequences[2].id, "query_two")
}

///|
test "a2m parses descriptions" {
  let alignment = a2m_sample()
  assert_eq(alignment.sequences[0].description, "SAM model reference")
  assert_eq(alignment.sequences[1].description, "insertion and deletion")
  assert_eq(alignment.sequences[2].description, "divergent homolog")
}

///|
test "a2m derives states from first row" {
  let alignment = a2m_sample()
  assert_eq(alignment.states, "DDDIIDDDDIDD")
  assert_eq(alignment.alignment_length(), 12)
}

///|
test "a2m normalizes match and insertion residues" {
  let alignment = a2m_sample()
  assert_eq(alignment.sequences[0].aligned_sequence, "ACDEFG-HIKLM")
  assert_eq(alignment.sequences[1].aligned_sequence, "ACD--GTHI-LM")
  assert_eq(alignment.sequences[2].aligned_sequence, "ATDQRG-HI-L-")
}

///|
test "a2m stores ungapped sequences" {
  let alignment = a2m_sample()
  assert_eq(alignment.sequences[0].sequence, "ACDEFGHIKLM")
  assert_eq(alignment.sequences[1].sequence, "ACDGTHILM")
  assert_eq(alignment.sequences[2].sequence, "ATDQRGHIL")
}

///|
test "a2m accepts wrapped sequence lines" {
  let text = ">ref description\nACDe\nfG-H\nIkLM\n>query\nACD.\n.GTH\nI.LM\n"
  let alignment = @src.a2m_parse(text) catch {
    _ => abort("wrapped A2M should parse")
  }
  assert_eq(alignment.num_sequences(), 2)
  assert_eq(alignment.states, "DDDIIDDDDIDD")
  assert_eq(alignment.sequences[0].aligned_sequence, "ACDEFG-HIKLM")
  assert_eq(alignment.sequences[1].aligned_sequence, "ACD--GTHI-LM")
}

///|
test "a2m ignores blank lines" {
  let text = "\n>ref\nACdE\n\n>query\nAC.E\n\n"
  let alignment = @src.a2m_parse(text) catch {
    _ => abort("A2M with blank lines should parse")
  }
  assert_eq(alignment.states, "DDID")
  assert_eq(alignment.num_sequences(), 2)
}

///|
test "a2m accepts CRLF input" {
  let text = @src.a2m_example_text().replace_all(old="\n", new="\r\n")
  assert_eq(
    @src.a2m_parse(text) catch {
      _ => abort("CRLF A2M should parse")
    },
    a2m_sample(),
  )
}

///|
test "a2m writes canonical character classes" {
  let alignment = a2m_sample()
  assert_eq(
    @src.a2m_write(alignment) catch {
      _ => abort("valid A2M should serialize")
    },
    @src.a2m_example_text(),
  )
}

///|
test "a2m writer wraps rows" {
  let text = @src.a2m_write(a2m_sample(), line_width=5) catch {
    _ => abort("valid wrapped serialization should succeed")
  }
  assert_true(text.contains("ACDef\nG-HIk\nLM\n"))
  assert_true(text.contains("ACD..\nGTHI.\nLM\n"))
}

///|
test "a2m wrapped serialization round trips" {
  let original = a2m_sample()
  let text = @src.a2m_write(original, line_width=4) catch {
    _ => abort("valid wrapped serialization should succeed")
  }
  let reparsed = @src.a2m_parse(text) catch {
    _ => abort("wrapped serialization should reparse")
  }
  assert_eq(reparsed, original)
}

///|
test "a2m reports match and insertion column counts" {
  let alignment = a2m_sample()
  assert_eq(alignment.match_columns(), 9)
  assert_eq(alignment.insertion_columns(), 3)
}

///|
test "a2m returns column states" {
  let alignment = a2m_sample()
  assert_eq(alignment.column_state(0), Some(@src.A2mMatch))
  assert_eq(alignment.column_state(3), Some(@src.A2mInsertion))
  assert_eq(alignment.column_state(9), Some(@src.A2mInsertion))
  assert_true(alignment.column_state(-1) is None)
  assert_true(alignment.column_state(12) is None)
}

///|
test "a2m returns alignment columns" {
  let alignment = a2m_sample()
  assert_eq(alignment.column(0), Some("AAA"))
  assert_eq(alignment.column(3), Some("E-Q"))
  assert_eq(alignment.column(6), Some("-T-"))
  assert_true(alignment.column(-1) is None)
  assert_true(alignment.column(12) is None)
}

///|
test "a2m finds sequence identifiers" {
  let alignment = a2m_sample()
  assert_eq(alignment.find_sequence("query_one"), Some(1))
  assert_eq(alignment.find_sequence("reference"), Some(0))
  assert_true(alignment.find_sequence("missing") is None)
}

///|
test "a2m reconstructs encoded rows" {
  let alignment = a2m_sample()
  assert_eq(alignment.encoded_sequence(0), Some("ACDefG-HIkLM"))
  assert_eq(alignment.encoded_sequence(1), Some("ACD..GTHI.LM"))
  assert_eq(alignment.encoded_sequence(2), Some("ATDqrG-HI.L-"))
  assert_true(alignment.encoded_sequence(-1) is None)
  assert_true(alignment.encoded_sequence(3) is None)
}

///|
test "a2m maps sequence coordinates to columns" {
  let alignment = a2m_sample()
  assert_eq(alignment.sequence_position_to_column(0, 0), Some(0))
  assert_eq(alignment.sequence_position_to_column(0, 3), Some(3))
  assert_eq(alignment.sequence_position_to_column(0, 6), Some(7))
  assert_eq(alignment.sequence_position_to_column(0, 10), Some(11))
  assert_true(alignment.sequence_position_to_column(0, 11) is None)
}

///|
test "a2m maps columns to sequence coordinates" {
  let alignment = a2m_sample()
  assert_eq(alignment.column_to_sequence_position(1, 0), Some(0))
  assert_eq(alignment.column_to_sequence_position(1, 5), Some(3))
  assert_eq(alignment.column_to_sequence_position(1, 6), Some(4))
  assert_true(alignment.column_to_sequence_position(1, 3) is None)
  assert_eq(alignment.column_to_sequence_position(1, 11), Some(8))
}

///|
test "a2m maps positions between rows" {
  let alignment = a2m_sample()
  assert_eq(alignment.map_position(0, 1, 0), Some(0))
  assert_true(alignment.map_position(0, 1, 3) is None)
  assert_eq(alignment.map_position(0, 1, 5), Some(3))
  assert_eq(alignment.map_position(0, 1, 6), Some(5))
  assert_true(alignment.map_position(0, 1, 8) is None)
}

///|
test "a2m maps insertion residues to reference gaps" {
  let alignment = a2m_sample()
  assert_eq(alignment.sequence_position_to_column(1, 4), Some(6))
  assert_true(alignment.map_position(1, 0, 4) is None)
}

///|
test "a2m returns aligned coordinate pairs" {
  let pairs = a2m_sample().aligned_pairs(0, 1) catch {
    _ => abort("valid row pair should succeed")
  }
  assert_eq(pairs.length(), 12)
  assert_eq(pairs[0], (Some(0), Some(0)))
  assert_eq(pairs[3], (Some(3), None))
  assert_eq(pairs[6], (None, Some(4)))
  assert_eq(pairs[11], (Some(10), Some(8)))
}

///|
test "a2m counts identical pair" {
  let counts = a2m_sample().pair_counts(0, 1) catch {
    _ => abort("valid row pair should count")
  }
  assert_eq(counts.columns, 12)
  assert_eq(counts.aligned, 8)
  assert_eq(counts.identities, 8)
  assert_eq(counts.mismatches, 0)
  assert_eq(counts.gap_columns, 4)
  assert_eq(counts.double_gap_columns, 0)
  assert_eq(counts.match_aligned, 8)
  assert_eq(counts.insertion_aligned, 0)
  assert_true((counts.identity() - 1.0).abs() < 1.0e-12)
}

///|
test "a2m counts divergent pair and insertion states" {
  let counts = a2m_sample().pair_counts(0, 2) catch {
    _ => abort("valid row pair should count")
  }
  assert_eq(counts.aligned, 9)
  assert_eq(counts.identities, 6)
  assert_eq(counts.mismatches, 3)
  assert_eq(counts.gap_columns, 2)
  assert_eq(counts.double_gap_columns, 1)
  assert_eq(counts.match_aligned, 7)
  assert_eq(counts.insertion_aligned, 2)
  assert_true((counts.identity() - 2.0 / 3.0).abs() < 1.0e-12)
}

///|
test "a2m identity is zero without paired residues" {
  let alignment = @src.a2m_from_aligned(["left", "right"], ["A-", "-C"], "DD") catch {
    _ => abort("valid disjoint alignment should build")
  }
  let counts = alignment.pair_counts(0, 1) catch {
    _ => abort("valid disjoint pair should count")
  }
  assert_eq(counts.aligned, 0)
  assert_true(counts.identity().abs() < 1.0e-12)
}

///|
test "a2m locates insertion runs" {
  let runs = a2m_sample().insertion_runs()
  assert_eq(runs.length(), 2)
  assert_eq(runs[0].slot, 3)
  assert_eq(runs[0].start_column, 3)
  assert_eq(runs[0].end_column, 5)
  assert_eq(runs[0].width, 2)
  assert_eq(runs[1].slot, 7)
  assert_eq(runs[1].start_column, 9)
  assert_eq(runs[1].end_column, 10)
}

///|
test "a2m handles terminal insertion slots" {
  let alignment = @src.a2m_parse(">ref\naBCd\n>query\n.EF.\n") catch {
    _ => abort("terminal insertion alignment should parse")
  }
  let runs = alignment.insertion_runs()
  assert_eq(runs.length(), 2)
  assert_eq(runs[0].slot, 0)
  assert_eq(runs[0].width, 1)
  assert_eq(runs[1].slot, 2)
  assert_eq(runs[1].width, 1)
}

///|
test "a2m computes all-column consensus" {
  assert_eq(
    a2m_sample().consensus() catch {
      _ => abort("valid consensus should succeed")
    },
    "ACDEFGTHIKLM",
  )
}

///|
test "a2m computes match-state consensus" {
  assert_eq(
    a2m_sample().consensus(include_insertions=false) catch {
      _ => abort("valid match consensus should succeed")
    },
    "ACDGTHILM",
  )
}

///|
test "a2m applies consensus threshold" {
  assert_eq(
    a2m_sample().consensus(minimum_fraction=0.8) catch {
      _ => abort("valid threshold consensus should succeed")
    },
    "AXDXXGTHIKLM",
  )
}

///|
test "a2m consensus breaks ties lexicographically" {
  let alignment = @src.a2m_from_aligned(["z", "a"], ["T", "C"], "D") catch {
    _ => abort("valid tie alignment should build")
  }
  assert_eq(
    alignment.consensus() catch {
      _ => abort("valid tie consensus should succeed")
    },
    "C",
  )
}

///|
test "a2m consensus emits gap for empty column" {
  let alignment = @src.a2m_from_aligned(["one", "two"], ["-", "-"], "D") catch {
    _ => abort("all-gap alignment should build")
  }
  assert_eq(
    alignment.consensus() catch {
      _ => abort("all-gap consensus should succeed")
    },
    "-",
  )
}

///|
test "a2m reports occupancy" {
  let occupancy = a2m_sample().occupancy()
  assert_eq(occupancy.length(), 12)
  assert_true((occupancy[0] - 1.0).abs() < 1.0e-12)
  assert_true((occupancy[3] - 2.0 / 3.0).abs() < 1.0e-12)
  assert_true((occupancy[6] - 1.0 / 3.0).abs() < 1.0e-12)
  assert_true((occupancy[9] - 1.0 / 3.0).abs() < 1.0e-12)
  assert_true((occupancy[11] - 2.0 / 3.0).abs() < 1.0e-12)
}

///|
test "a2m projects match columns" {
  let projected = a2m_sample().match_projection() catch {
    _ => abort("match projection should succeed")
  }
  assert_eq(projected.states, "DDDDDDDDD")
  assert_eq(projected.sequences[0].aligned_sequence, "ACDG-HILM")
  assert_eq(projected.sequences[1].aligned_sequence, "ACDGTHILM")
  assert_eq(projected.sequences[2].aligned_sequence, "ATDG-HIL-")
  assert_eq(projected.sequences[0].description, "SAM model reference")
}

///|
test "a2m slices columns and preserves states" {
  let sliced = a2m_sample().slice_columns(2, 10) catch {
    _ => abort("valid A2M slice should succeed")
  }
  assert_eq(sliced.states, "DIIDDDDI")
  assert_eq(sliced.sequences[0].aligned_sequence, "DEFG-HIK")
  assert_eq(sliced.sequences[1].aligned_sequence, "D--GTHI-")
  assert_eq(sliced.sequences[2].aligned_sequence, "DQRG-HI-")
}

///|
test "a2m summarizes alignment" {
  assert_eq(
    a2m_sample().summary(),
    "A2mAlignment(sequences=3, columns=12, match=9, insertion=3)",
  )
}

///|
test "a2m sequence constructor normalizes case" {
  let sequence = @src.A2mSequence::create(
    "protein",
    "AcD-Ef",
    description="mixed case",
  ) catch {
    _ => abort("valid normalized row should build")
  }
  assert_eq(sequence.aligned_sequence, "ACD-EF")
  assert_eq(sequence.sequence, "ACDEF")
  assert_eq(sequence.description, "mixed case")
}

///|
test "a2m alignment constructor copies rows" {
  let rows = [
    @src.A2mSequence::create("one", "AC-") catch {
      _ => abort("valid row should build")
    },
    @src.A2mSequence::create("two", "A-C") catch {
      _ => abort("valid row should build")
    },
  ]
  let alignment = @src.A2mAlignment::create(rows, "DID") catch {
    _ => abort("valid alignment should build")
  }
  rows.push(
    @src.A2mSequence::create("three", "ACC") catch {
      _ => abort("valid row should build")
    },
  )
  assert_eq(alignment.num_sequences(), 2)
}

///|
test "a2m builds from aligned rows without descriptions" {
  let alignment = @src.a2m_from_aligned(["one", "two"], ["AC-", "A-C"], "DID") catch {
    _ => abort("valid direct alignment should build")
  }
  assert_eq(alignment.sequences[0].description, "")
  assert_eq(alignment.encoded_sequence(0), Some("Ac-"))
  assert_eq(alignment.encoded_sequence(1), Some("A.C"))
}

///|
test "a2m builds from aligned rows with descriptions" {
  let alignment = @src.a2m_from_aligned(["one", "two"], ["AC-", "A-C"], "DID", descriptions=[
    "first row", "second row",
  ]) catch {
    _ => abort("valid described alignment should build")
  }
  assert_eq(alignment.sequences[1].description, "second row")
}

///|
test "a2m accepts duplicate identifiers like Biopython" {
  let alignment = @src.a2m_parse(">same\nAC\n>same\nAT\n") catch {
    _ => abort("duplicate A2M identifiers should parse")
  }
  assert_eq(alignment.num_sequences(), 2)
  assert_eq(alignment.find_sequence("same"), Some(0))
}

///|
test "a2m supports a single sequence" {
  let alignment = @src.a2m_parse(">only\nACdE\n") catch {
    _ => abort("single-row A2M should parse")
  }
  assert_eq(alignment.num_sequences(), 1)
  assert_eq(alignment.states, "DDID")
  assert_eq(alignment.occupancy(), [1.0, 1.0, 1.0, 1.0])
}

///|
test "a2m rejects empty input" {
  assert_true(a2m_raises(""))
  assert_true(a2m_raises("\n\r\n"))
}

///|
test "a2m rejects sequence before header" {
  assert_true(a2m_raises("ACDE\n>query\nACDE\n"))
}

///|
test "a2m rejects empty header" {
  assert_true(a2m_raises(">\nACDE\n"))
  assert_true(a2m_raises(">   \nACDE\n"))
}

///|
test "a2m rejects empty first row" {
  assert_true(a2m_raises(">empty\n>filled\nACDE\n"))
}

///|
test "a2m rejects unequal row widths" {
  assert_true(a2m_raises(">one\nACDE\n>two\nACD\n"))
  assert_true(a2m_raises(">one\nACD\n>two\nACDE\n"))
}

///|
test "a2m rejects invalid sequence characters" {
  assert_true(a2m_raises(">one\nAC*D\n"))
  assert_true(a2m_raises(">one\nAC_D\n"))
  assert_true(a2m_raises(">one\nAC/D\n"))
}

///|
test "a2m rejects embedded sequence whitespace" {
  assert_true(a2m_raises(">one\nAC DE\n"))
  assert_true(a2m_raises(">one\nAC\tDE\n"))
}

///|
test "a2m rejects insertion character in match column" {
  assert_true(a2m_raises(">model\nACDE\n>query\nACdE\n"))
  assert_true(a2m_raises(">model\nAC-E\n>query\nAC.E\n"))
}

///|
test "a2m rejects match character in insertion column" {
  assert_true(a2m_raises(">model\nACdE\n>query\nACDE\n"))
  assert_true(a2m_raises(">model\nAC.E\n>query\nAC-E\n"))
}

///|
test "a2m rejects invalid constructor state" {
  let raised = try {
    ignore(@src.a2m_from_aligned(["one"], ["ACD"], "DXD"))
    false
  } catch {
    A2mError(_) => true
  }
  assert_true(raised)
}

///|
test "a2m rejects state width mismatch" {
  let raised = try {
    ignore(@src.a2m_from_aligned(["one"], ["ACD"], "DD"))
    false
  } catch {
    A2mError(_) => true
  }
  assert_true(raised)
}

///|
test "a2m rejects empty alignment constructor" {
  let rows : Array[@src.A2mSequence] = []
  let raised = try {
    ignore(@src.A2mAlignment::create(rows, "D"))
    false
  } catch {
    A2mError(_) => true
  }
  assert_true(raised)
}

///|
test "a2m rejects empty sequence constructor" {
  let raised = try {
    ignore(@src.A2mSequence::create("one", ""))
    false
  } catch {
    A2mError(_) => true
  }
  assert_true(raised)
}

///|
test "a2m rejects invalid normalized row characters" {
  let dot = try {
    ignore(@src.A2mSequence::create("one", "AC.D"))
    false
  } catch {
    A2mError(_) => true
  }
  let star = try {
    ignore(@src.A2mSequence::create("one", "AC*D"))
    false
  } catch {
    A2mError(_) => true
  }
  assert_true(dot)
  assert_true(star)
}

///|
test "a2m rejects empty sequence identifier" {
  let raised = try {
    ignore(@src.A2mSequence::create("", "ACD"))
    false
  } catch {
    A2mError(_) => true
  }
  assert_true(raised)
}

///|
test "a2m rejects whitespace in sequence identifier" {
  let space = try {
    ignore(@src.A2mSequence::create("bad id", "ACD"))
    false
  } catch {
    A2mError(_) => true
  }
  let tab = try {
    ignore(@src.A2mSequence::create("bad\tid", "ACD"))
    false
  } catch {
    A2mError(_) => true
  }
  assert_true(space)
  assert_true(tab)
}

///|
test "a2m rejects line breaks in description" {
  let raised = try {
    ignore(
      @src.A2mSequence::create("one", "ACD", description="bad\ndescription"),
    )
    false
  } catch {
    A2mError(_) => true
  }
  assert_true(raised)
}

///|
test "a2m rejects direct identifier row count mismatch" {
  let raised = try {
    ignore(@src.a2m_from_aligned(["one", "two"], ["ACD"], "DDD"))
    false
  } catch {
    A2mError(_) => true
  }
  assert_true(raised)
}

///|
test "a2m rejects direct description count mismatch" {
  let raised = try {
    ignore(
      @src.a2m_from_aligned(["one", "two"], ["ACD", "ACD"], "DDD", descriptions=[
        "only one",
      ]),
    )
    false
  } catch {
    A2mError(_) => true
  }
  assert_true(raised)
}

///|
test "a2m rejects negative writer width" {
  let raised = try {
    ignore(@src.a2m_write(a2m_sample(), line_width=-1))
    false
  } catch {
    A2mError(_) => true
  }
  assert_true(raised)
}

///|
test "a2m rejects match projection without match states" {
  let alignment = @src.a2m_from_aligned(["one", "two"], ["AC", "A-"], "II") catch {
    _ => abort("insertion-only alignment should build")
  }
  let raised = try {
    ignore(alignment.match_projection())
    false
  } catch {
    A2mError(_) => true
  }
  assert_true(raised)
}

///|
test "a2m rejects negative column slice" {
  let raised = try {
    ignore(a2m_sample().slice_columns(-1, 2))
    false
  } catch {
    A2mError(_) => true
  }
  assert_true(raised)
}

///|
test "a2m rejects empty column slice" {
  let raised = try {
    ignore(a2m_sample().slice_columns(2, 2))
    false
  } catch {
    A2mError(_) => true
  }
  assert_true(raised)
}

///|
test "a2m rejects column slice past end" {
  let raised = try {
    ignore(a2m_sample().slice_columns(0, 13))
    false
  } catch {
    A2mError(_) => true
  }
  assert_true(raised)
}

///|
test "a2m rejects invalid aligned pair row" {
  let low = try {
    ignore(a2m_sample().aligned_pairs(-1, 0))
    false
  } catch {
    A2mError(_) => true
  }
  let high = try {
    ignore(a2m_sample().aligned_pairs(0, 3))
    false
  } catch {
    A2mError(_) => true
  }
  assert_true(low)
  assert_true(high)
}

///|
test "a2m rejects invalid pair count row" {
  let raised = try {
    ignore(a2m_sample().pair_counts(3, 0))
    false
  } catch {
    A2mError(_) => true
  }
  assert_true(raised)
}

///|
test "a2m rejects invalid consensus threshold" {
  let low = try {
    ignore(a2m_sample().consensus(minimum_fraction=-0.1))
    false
  } catch {
    A2mError(_) => true
  }
  let high = try {
    ignore(a2m_sample().consensus(minimum_fraction=1.1))
    false
  } catch {
    A2mError(_) => true
  }
  assert_true(low)
  assert_true(high)
}

///|
test "a2m accessors reject invalid rows and positions" {
  let alignment = a2m_sample()
  assert_true(alignment.sequence_position_to_column(-1, 0) is None)
  assert_true(alignment.sequence_position_to_column(3, 0) is None)
  assert_true(alignment.sequence_position_to_column(0, -1) is None)
  assert_true(alignment.column_to_sequence_position(-1, 0) is None)
  assert_true(alignment.column_to_sequence_position(3, 0) is None)
  assert_true(alignment.column_to_sequence_position(0, -1) is None)
  assert_true(alignment.column_to_sequence_position(0, 12) is None)
}

///|
test "a2m map rejects invalid target row" {
  let alignment = a2m_sample()
  assert_true(alignment.map_position(0, -1, 0) is None)
  assert_true(alignment.map_position(0, 3, 0) is None)
  assert_true(alignment.map_position(3, 0, 0) is None)
  assert_true(alignment.map_position(0, 1, -1) is None)
}

///|
test "a2m insertion-only alignment retains coordinates" {
  let alignment = @src.a2m_from_aligned(["one", "two"], ["AC", "A-"], "II") catch {
    _ => abort("insertion-only alignment should build")
  }
  assert_eq(alignment.match_columns(), 0)
  assert_eq(alignment.insertion_columns(), 2)
  assert_eq(alignment.map_position(0, 1, 0), Some(0))
  assert_true(alignment.map_position(0, 1, 1) is None)
  assert_eq(alignment.insertion_runs()[0].slot, 0)
}
