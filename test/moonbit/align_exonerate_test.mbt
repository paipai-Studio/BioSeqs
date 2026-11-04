// Black-box tests for Biopython Bio.Align.exonerate-compatible support.

///|
fn align_exonerate_test_parse(text : String) -> @src.AlignExonerateDocument {
  @src.align_exonerate_parse(text) catch {
    AlignExonerateError(message) =>
      abort("valid alignment-aware Exonerate input failed: " + message)
  }
}

///|
fn align_exonerate_test_raises(text : String) -> Bool {
  try {
    ignore(@src.align_exonerate_parse(text))
    false
  } catch {
    AlignExonerateError(_) => true
  }
}

///|
fn align_exonerate_test_write_raises(
  document : @src.AlignExonerateDocument,
  format : String,
) -> Bool {
  try {
    ignore(@src.align_exonerate_write(document, format~))
    false
  } catch {
    AlignExonerateError(_) => true
  }
}

///|
fn align_exonerate_test_sample() -> @src.AlignExonerateDocument {
  align_exonerate_test_parse(@src.align_exonerate_example_text())
}

///|
fn align_exonerate_test_cigar_text() -> String {
  "Command line: [exonerate -m affine:local q.fa t.fa]\n" +
  "Hostname: [blackbriar]\n" +
  "cigar: query 0 8 + target 20 11 - 42 M 3 I 2 M 2 D 3 M 1\n" +
  "-- completed exonerate analysis\n"
}

///|
fn align_exonerate_test_all_operations_text() -> String {
  "Command line: [exonerate -m coding2genome q.fa t.fa]\n" +
  "Hostname: [worker]\n" +
  "vulgar: q 0 15 + t 0 31 + 7 M 3 3 5 0 2 I 0 10 3 0 2 C 3 3 G 2 0 G 0 4 N 5 6 S 1 1 F 1 0\n" +
  "-- completed exonerate analysis\n"
}

///|
fn align_exonerate_test_empty_text() -> String {
  "Command line: [exonerate -m est2genome q.fa t.fa]\n" +
  "Hostname: [worker]\n" +
  "-- completed exonerate analysis\n"
}

///|
test "Bio.Align.exonerate parses metadata" {
  let document = align_exonerate_test_sample()
  assert_eq(document.metadata.program, "exonerate")
  assert_eq(
    document.metadata.command_line,
    "exonerate -m est2genome query.fa target.fa --showvulgar yes",
  )
  assert_eq(document.metadata.hostname, "moonbit")
}

///|
test "Bio.Align.exonerate report summary" {
  assert_eq(
    align_exonerate_test_sample().summary(),
    "Exonerate report(host=moonbit, alignments=2)",
  )
}

///|
test "Bio.Align.exonerate document accessors" {
  let document = align_exonerate_test_sample()
  assert_eq(document.num_alignments(), 2)
  assert_true(document.get(0) is Some(_))
  assert_true(document.get(-1) is None)
  assert_true(document.get(2) is None)
}

///|
test "Bio.Align.exonerate query lookup" {
  let document = align_exonerate_test_sample()
  assert_eq(document.find_by_query("transcript").length(), 1)
  assert_eq(document.find_by_query("missing").length(), 0)
}

///|
test "Bio.Align.exonerate selects best score" {
  let best = align_exonerate_test_sample().best_alignment()
  match best {
    Some(alignment) => {
      assert_eq(alignment.query_id, "transcript")
      assert_eq(alignment.score, 250.0)
    }
    None => abort("sample must have a best alignment")
  }
}

///|
test "Bio.Align.exonerate parses identifiers and score" {
  let alignment = align_exonerate_test_sample().alignments[0]
  assert_eq(alignment.query_id, "transcript")
  assert_eq(alignment.target_id, "chromosome")
  assert_eq(alignment.score, 250.0)
}

///|
test "Bio.Align.exonerate normalizes spliced operations" {
  let alignment = align_exonerate_test_sample().alignments[0]
  assert_eq(alignment.operation_codes(), "M5N3M")
  assert_eq(alignment.num_operations(), 5)
  assert_eq(alignment.operations[2].query_step, 0)
  assert_eq(alignment.operations[2].target_step, 96)
}

///|
test "Bio.Align.exonerate builds spliced coordinate path" {
  let path = align_exonerate_test_sample().alignments[0].coordinate_path()
  assert_eq(path[0], [100, 106, 108, 204, 206, 218])
  assert_eq(path[1], [0, 6, 6, 6, 6, 18])
}

///|
test "Bio.Align.exonerate preserves forward bounds" {
  let alignment = align_exonerate_test_sample().alignments[0]
  assert_eq(alignment.query_start, 0)
  assert_eq(alignment.query_end, 18)
  assert_eq(alignment.target_start, 100)
  assert_eq(alignment.target_end, 218)
  assert_eq(alignment.query_strand, @src.AlignExonerateForward)
  assert_eq(alignment.target_strand, @src.AlignExonerateForward)
}

///|
test "Bio.Align.exonerate reports inferred sequence lengths" {
  let first = align_exonerate_test_sample().alignments[0]
  let second = align_exonerate_test_sample().alignments[1]
  assert_eq(first.query_length(), 18)
  assert_eq(first.target_length(), 218)
  assert_eq(second.query_length(), 4)
  assert_eq(second.target_length(), 500)
}

///|
test "Bio.Align.exonerate maps forward query positions" {
  let alignment = align_exonerate_test_sample().alignments[0]
  assert_eq(alignment.query_position_to_target(0), Some(100))
  assert_eq(alignment.query_position_to_target(5), Some(105))
  assert_eq(alignment.query_position_to_target(6), Some(206))
  assert_eq(alignment.query_position_to_target(17), Some(217))
}

///|
test "Bio.Align.exonerate maps forward target positions" {
  let alignment = align_exonerate_test_sample().alignments[0]
  assert_eq(alignment.target_position_to_query(100), Some(0))
  assert_eq(alignment.target_position_to_query(217), Some(17))
}

///|
test "Bio.Align.exonerate mapping returns none in splice and intron" {
  let alignment = align_exonerate_test_sample().alignments[0]
  assert_eq(alignment.target_position_to_query(106), None)
  assert_eq(alignment.target_position_to_query(108), None)
  assert_eq(alignment.target_position_to_query(204), None)
}

///|
test "Bio.Align.exonerate mapping rejects outside positions" {
  let alignment = align_exonerate_test_sample().alignments[0]
  assert_eq(alignment.query_position_to_target(-1), None)
  assert_eq(alignment.query_position_to_target(18), None)
  assert_eq(alignment.target_position_to_query(218), None)
}

///|
test "Bio.Align.exonerate parses protein strand" {
  let alignment = align_exonerate_test_sample().alignments[1]
  assert_eq(alignment.query_strand, @src.AlignExonerateProtein)
  assert_eq(alignment.target_strand, @src.AlignExonerateReverse)
  assert_eq(alignment.operation_codes(), "M")
}

///|
test "Bio.Align.exonerate builds translated reverse path" {
  let path = align_exonerate_test_sample().alignments[1].coordinate_path()
  assert_eq(path[0], [500, 488])
  assert_eq(path[1], [0, 4])
}

///|
test "Bio.Align.exonerate maps protein to reverse codons" {
  let alignment = align_exonerate_test_sample().alignments[1]
  assert_eq(alignment.query_position_to_target(0), Some(499))
  assert_eq(alignment.query_position_to_target(1), Some(496))
  assert_eq(alignment.query_position_to_target(3), Some(490))
}

///|
test "Bio.Align.exonerate maps reverse DNA to protein" {
  let alignment = align_exonerate_test_sample().alignments[1]
  assert_eq(alignment.target_position_to_query(499), Some(0))
  assert_eq(alignment.target_position_to_query(498), Some(0))
  assert_eq(alignment.target_position_to_query(497), Some(0))
  assert_eq(alignment.target_position_to_query(496), Some(1))
}

///|
test "Bio.Align.exonerate emits translated aligned pairs" {
  let pairs = align_exonerate_test_sample().alignments[1].aligned_pairs()
  assert_eq(pairs.length(), 4)
  assert_eq(pairs[0], (499, 0))
  assert_eq(pairs[1], (496, 1))
  assert_eq(pairs[3], (490, 3))
}

///|
test "Bio.Align.exonerate counts spliced alignment" {
  let counts = align_exonerate_test_sample().alignments[0].counts()
  assert_eq(counts.operations, 5)
  assert_eq(counts.matching_operations, 2)
  assert_eq(counts.aligned_query_units, 18)
  assert_eq(counts.aligned_target_units, 18)
  assert_eq(counts.introns, 1)
  assert_eq(counts.intron_units, 96)
}

///|
test "Bio.Align.exonerate counts contiguous target gaps" {
  let counts = align_exonerate_test_sample().alignments[0].counts()
  assert_eq(counts.insertions, 0)
  assert_eq(counts.deletions, 100)
  assert_eq(counts.insertion_opens, 0)
  assert_eq(counts.deletion_opens, 1)
}

///|
test "Bio.Align.exonerate alignment summary" {
  assert_eq(
    align_exonerate_test_sample().alignments[1].summary(),
    "Exonerate alignment(query=protein, target=genome, score=80, operations=1)",
  )
}

///|
test "Bio.Align.exonerate parses reverse cigar" {
  let alignment = align_exonerate_test_parse(align_exonerate_test_cigar_text()).alignments[0]
  assert_eq(alignment.query_strand, @src.AlignExonerateForward)
  assert_eq(alignment.target_strand, @src.AlignExonerateReverse)
  assert_eq(alignment.operation_codes(), "MIMDM")
}

///|
test "Bio.Align.exonerate cigar coordinate path" {
  let path = align_exonerate_test_parse(align_exonerate_test_cigar_text()).alignments[0].coordinate_path()
  assert_eq(path[0], [20, 17, 17, 15, 12, 11])
  assert_eq(path[1], [0, 3, 5, 7, 7, 8])
}

///|
test "Bio.Align.exonerate parses protein2dna cigar" {
  let text = "Command line: [protein2dna]\n" +
    "Hostname: [host]\n" +
    "cigar: protein 0 4 . genome 100 91 - 20 M 6 I 1 M 3\n" +
    "-- completed exonerate analysis\n"
  let alignment = align_exonerate_test_parse(text).alignments[0]
  assert_eq(alignment.operation_codes(), "MIM")
  assert_eq(alignment.operations[0].query_step, 2)
  assert_eq(alignment.operations[0].target_step, 6)
  assert_eq(alignment.operations[1].query_step, 1)
}

///|
test "Bio.Align.exonerate protein2dna cigar path" {
  let text = "Command line: [protein2dna]\n" +
    "Hostname: [host]\n" +
    "cigar: protein 0 4 . genome 100 91 - 20 M 6 I 1 M 3\n" +
    "-- completed exonerate analysis\n"
  let path = align_exonerate_test_parse(text).alignments[0].coordinate_path()
  assert_eq(path[0], [100, 94, 94, 91])
  assert_eq(path[1], [0, 2, 3, 4])
}

///|
test "Bio.Align.exonerate parses dna2protein cigar" {
  let text = "Command line: [dna2protein]\n" +
    "Hostname: [host]\n" +
    "cigar: dna 0 9 + protein 0 4 . 20 M 6 D 1 M 3\n" +
    "-- completed exonerate analysis\n"
  let alignment = align_exonerate_test_parse(text).alignments[0]
  assert_eq(alignment.operation_codes(), "MDM")
  assert_eq(alignment.operations[0].query_step, 6)
  assert_eq(alignment.operations[0].target_step, 2)
  assert_eq(alignment.operations[1].target_step, 1)
}

///|
test "Bio.Align.exonerate dna2protein mapping" {
  let text = "Command line: [dna2protein]\n" +
    "Hostname: [host]\n" +
    "cigar: dna 0 9 + protein 0 4 . 20 M 6 D 1 M 3\n" +
    "-- completed exonerate analysis\n"
  let alignment = align_exonerate_test_parse(text).alignments[0]
  assert_eq(alignment.query_position_to_target(0), Some(0))
  assert_eq(alignment.query_position_to_target(3), Some(1))
  assert_eq(alignment.target_position_to_query(2), None)
  assert_eq(alignment.query_position_to_target(6), Some(3))
}

///|
test "Bio.Align.exonerate parses every vulgar operation" {
  let alignment = align_exonerate_test_parse(
      align_exonerate_test_all_operations_text(),
    ).alignments[0]
  assert_eq(alignment.operation_codes(), "M5N3CIDUUSF")
  assert_eq(alignment.num_operations(), 11)
}

///|
test "Bio.Align.exonerate splits non-equivalenced operation" {
  let alignment = align_exonerate_test_parse(
      align_exonerate_test_all_operations_text(),
    ).alignments[0]
  assert_eq(alignment.operations[7].code, "U")
  assert_eq(alignment.operations[7].query_step, 0)
  assert_eq(alignment.operations[7].target_step, 6)
  assert_eq(alignment.operations[8].query_step, 5)
  assert_eq(alignment.operations[8].target_step, 0)
}

///|
test "Bio.Align.exonerate counts special operations" {
  let counts = align_exonerate_test_parse(
    align_exonerate_test_all_operations_text(),
  ).alignments[0].counts()
  assert_eq(counts.introns, 1)
  assert_eq(counts.intron_units, 10)
  assert_eq(counts.non_equivalenced_units, 11)
  assert_eq(counts.split_codons, 1)
  assert_eq(counts.frameshifts, 1)
}

///|
test "Bio.Align.exonerate counts insertions and deletions" {
  let counts = align_exonerate_test_parse(
    align_exonerate_test_all_operations_text(),
  ).alignments[0].counts()
  assert_eq(counts.insertions, 8)
  assert_eq(counts.deletions, 24)
  assert_eq(counts.insertion_opens, 3)
  assert_eq(counts.deletion_opens, 2)
}

///|
test "Bio.Align.exonerate vulgar canonical round trip" {
  let document = align_exonerate_test_parse(
    align_exonerate_test_all_operations_text(),
  )
  let output = @src.align_exonerate_write(document, format="vulgar") catch {
    AlignExonerateError(message) => abort("vulgar writing failed: " + message)
  }
  assert_eq(align_exonerate_test_parse(output), document)
}

///|
test "Bio.Align.exonerate sample vulgar round trip" {
  let document = align_exonerate_test_sample()
  let output = @src.align_exonerate_write(document) catch {
    AlignExonerateError(message) =>
      abort("default vulgar writing failed: " + message)
  }
  assert_eq(align_exonerate_test_parse(output), document)
}

///|
test "Bio.Align.exonerate cigar round trip preserves cigar paths" {
  let document = align_exonerate_test_parse(align_exonerate_test_cigar_text())
  let output = @src.align_exonerate_write(document, format="cigar") catch {
    AlignExonerateError(message) => abort("cigar writing failed: " + message)
  }
  assert_eq(align_exonerate_test_parse(output), document)
}

///|
test "Bio.Align.exonerate vulgar to cigar preserves paths" {
  let document = align_exonerate_test_sample()
  let output = @src.align_exonerate_write(document, format="cigar") catch {
    AlignExonerateError(message) =>
      abort("vulgar-to-cigar writing failed: " + message)
  }
  let converted = align_exonerate_test_parse(output)
  assert_eq(
    converted.alignments[0].coordinate_path(),
    document.alignments[0].coordinate_path(),
  )
  assert_eq(
    converted.alignments[1].coordinate_path(),
    document.alignments[1].coordinate_path(),
  )
}

///|
test "Bio.Align.exonerate translated cigar round trip" {
  let document = align_exonerate_test_sample()
  let output = @src.align_exonerate_write(document, format="cigar") catch {
    AlignExonerateError(message) =>
      abort("translated cigar writing failed: " + message)
  }
  let converted = align_exonerate_test_parse(output)
  assert_eq(converted.alignments[1].operations[0].query_step, 4)
  assert_eq(converted.alignments[1].operations[0].target_step, 12)
}

///|
test "Bio.Align.exonerate writes canonical metadata and footer" {
  let output = @src.align_exonerate_write(
    align_exonerate_test_sample(),
    format="vulgar",
  ) catch {
    AlignExonerateError(message) =>
      abort("canonical writing failed: " + message)
  }
  assert_true(output.has_prefix("Command line: ["))
  assert_true(output.contains("\nHostname: [moonbit]\n"))
  assert_true(output.has_suffix("-- completed exonerate analysis\n"))
}

///|
test "Bio.Align.exonerate supports valid zero-hit report" {
  let document = align_exonerate_test_parse(align_exonerate_test_empty_text())
  assert_eq(document.num_alignments(), 0)
  assert_eq(document.best_alignment(), None)
}

///|
test "Bio.Align.exonerate zero-hit round trip" {
  let document = align_exonerate_test_parse(align_exonerate_test_empty_text())
  let output = @src.align_exonerate_write(document) catch {
    AlignExonerateError(message) => abort("zero-hit writing failed: " + message)
  }
  assert_eq(align_exonerate_test_parse(output), document)
}

///|
test "Bio.Align.exonerate accepts CRLF" {
  let text = @src.align_exonerate_example_text().replace_all(
    old="\n",
    new="\r\n",
  )
  assert_eq(align_exonerate_test_parse(text).num_alignments(), 2)
}

///|
test "Bio.Align.exonerate accepts blank body lines" {
  let text = @src.align_exonerate_example_text().replace_all(
    old="Hostname: [moonbit]\n",
    new="Hostname: [moonbit]\n\n",
  )
  assert_eq(align_exonerate_test_parse(text).num_alignments(), 2)
}

///|
test "Bio.Align.exonerate parses mixed cigar and vulgar lines" {
  let text = "Command line: [mixed]\n" +
    "Hostname: [host]\n" +
    "vulgar: q 0 3 + t 0 3 + 1 M 3 3\n" +
    "cigar: q 0 3 + t 0 3 + 1 M 3\n" +
    "-- completed exonerate analysis\n"
  let document = align_exonerate_test_parse(text)
  assert_eq(document.num_alignments(), 2)
  assert_eq(
    document.alignments[0].coordinate_path(),
    document.alignments[1].coordinate_path(),
  )
}

///|
test "Bio.Align.exonerate metadata constructor" {
  let metadata = @src.AlignExonerateMetadata::create(
    command_line="exonerate query.fa target.fa",
    hostname="node",
  ) catch {
    AlignExonerateError(message) =>
      abort("metadata construction failed: " + message)
  }
  assert_eq(metadata.program, "exonerate")
  assert_eq(metadata.hostname, "node")
}

///|
test "Bio.Align.exonerate operation constructors" {
  let match_op = @src.AlignExonerateOperation::create("M", 3, 3) catch {
    AlignExonerateError(message) =>
      abort("match construction failed: " + message)
  }
  let intron = @src.AlignExonerateOperation::create("N", 0, 10) catch {
    AlignExonerateError(message) =>
      abort("intron construction failed: " + message)
  }
  assert_eq(match_op.code, "M")
  assert_eq(intron.target_step, 10)
}

///|
test "Bio.Align.exonerate alignment constructor" {
  let operation = @src.AlignExonerateOperation::create("M", 3, 3) catch {
    AlignExonerateError(message) =>
      abort("operation construction failed: " + message)
  }
  let alignment = @src.AlignExonerateAlignment::create(
    "q",
    0,
    3,
    @src.AlignExonerateForward,
    "t",
    10,
    7,
    @src.AlignExonerateReverse,
    5.0,
    [operation],
  ) catch {
    AlignExonerateError(message) =>
      abort("alignment construction failed: " + message)
  }
  assert_eq(alignment.coordinate_path()[0], [10, 7])
}

///|
test "Bio.Align.exonerate document constructor copies alignments" {
  let metadata = @src.AlignExonerateMetadata::create() catch {
    AlignExonerateError(message) =>
      abort("metadata construction failed: " + message)
  }
  let alignments : Array[@src.AlignExonerateAlignment] = []
  let document = @src.AlignExonerateDocument::create(metadata, alignments) catch {
    AlignExonerateError(message) =>
      abort("document construction failed: " + message)
  }
  alignments.push(align_exonerate_test_sample().alignments[0])
  assert_eq(document.num_alignments(), 0)
}

///|
test "Bio.Align.exonerate rejects empty input" {
  assert_true(align_exonerate_test_raises(""))
}

///|
test "Bio.Align.exonerate rejects missing command header" {
  assert_true(
    align_exonerate_test_raises(
      "Hostname: [host]\n-- completed exonerate analysis\n",
    ),
  )
}

///|
test "Bio.Align.exonerate rejects missing hostname header" {
  assert_true(
    align_exonerate_test_raises(
      "Command line: [cmd]\n-- completed exonerate analysis\n",
    ),
  )
}

///|
test "Bio.Align.exonerate rejects unbracketed command" {
  assert_true(
    align_exonerate_test_raises(
      "Command line: cmd\nHostname: [host]\n-- completed exonerate analysis\n",
    ),
  )
}

///|
test "Bio.Align.exonerate rejects unbracketed hostname" {
  assert_true(
    align_exonerate_test_raises(
      "Command line: [cmd]\nHostname: host\n-- completed exonerate analysis\n",
    ),
  )
}

///|
test "Bio.Align.exonerate rejects missing footer" {
  assert_true(
    align_exonerate_test_raises("Command line: [cmd]\nHostname: [host]\n"),
  )
}

///|
test "Bio.Align.exonerate rejects data after footer" {
  assert_true(
    align_exonerate_test_raises(
      align_exonerate_test_empty_text() + "vulgar: q 0 1 + t 0 1 + 1 M 1 1\n",
    ),
  )
}

///|
test "Bio.Align.exonerate rejects unknown body line" {
  assert_true(
    align_exonerate_test_raises(
      "Command line: [cmd]\nHostname: [host]\nunknown\n-- completed exonerate analysis\n",
    ),
  )
}

///|
test "Bio.Align.exonerate rejects malformed vulgar field count" {
  let text = "Command line: [cmd]\nHostname: [host]\n" +
    "vulgar: q 0 3 + t 0 3 + 1 M 3\n" +
    "-- completed exonerate analysis\n"
  assert_true(align_exonerate_test_raises(text))
}

///|
test "Bio.Align.exonerate rejects malformed cigar field count" {
  let text = "Command line: [cmd]\nHostname: [host]\n" +
    "cigar: q 0 3 + t 0 3 + 1 M\n" +
    "-- completed exonerate analysis\n"
  assert_true(align_exonerate_test_raises(text))
}

///|
test "Bio.Align.exonerate rejects negative coordinate" {
  let text = "Command line: [cmd]\nHostname: [host]\n" +
    "vulgar: q -1 2 + t 0 3 + 1 M 3 3\n" +
    "-- completed exonerate analysis\n"
  assert_true(align_exonerate_test_raises(text))
}

///|
test "Bio.Align.exonerate rejects coordinate overflow" {
  let text = "Command line: [cmd]\nHostname: [host]\n" +
    "vulgar: q 0 2147483648 + t 0 1 + 1 M 1 1\n" +
    "-- completed exonerate analysis\n"
  assert_true(align_exonerate_test_raises(text))
}

///|
test "Bio.Align.exonerate rejects invalid score" {
  let text = "Command line: [cmd]\nHostname: [host]\n" +
    "vulgar: q 0 3 + t 0 3 + score M 3 3\n" +
    "-- completed exonerate analysis\n"
  assert_true(align_exonerate_test_raises(text))
}

///|
test "Bio.Align.exonerate rejects non-finite score" {
  let text = "Command line: [cmd]\nHostname: [host]\n" +
    "vulgar: q 0 3 + t 0 3 + 1e301 M 3 3\n" +
    "-- completed exonerate analysis\n"
  assert_true(align_exonerate_test_raises(text))
}

///|
test "Bio.Align.exonerate rejects invalid query strand" {
  let text = "Command line: [cmd]\nHostname: [host]\n" +
    "vulgar: q 0 3 ? t 0 3 + 1 M 3 3\n" +
    "-- completed exonerate analysis\n"
  assert_true(align_exonerate_test_raises(text))
}

///|
test "Bio.Align.exonerate rejects invalid target strand" {
  let text = "Command line: [cmd]\nHostname: [host]\n" +
    "vulgar: q 0 3 + t 0 3 ? 1 M 3 3\n" +
    "-- completed exonerate analysis\n"
  assert_true(align_exonerate_test_raises(text))
}

///|
test "Bio.Align.exonerate rejects forward decreasing coordinates" {
  let text = "Command line: [cmd]\nHostname: [host]\n" +
    "vulgar: q 3 0 + t 0 3 + 1 M 3 3\n" +
    "-- completed exonerate analysis\n"
  assert_true(align_exonerate_test_raises(text))
}

///|
test "Bio.Align.exonerate rejects reverse increasing coordinates" {
  let text = "Command line: [cmd]\nHostname: [host]\n" +
    "vulgar: q 0 3 - t 0 3 + 1 M 3 3\n" +
    "-- completed exonerate analysis\n"
  assert_true(align_exonerate_test_raises(text))
}

///|
test "Bio.Align.exonerate rejects protein decreasing coordinates" {
  let text = "Command line: [cmd]\nHostname: [host]\n" +
    "vulgar: q 3 0 . t 0 3 + 1 M 3 3\n" +
    "-- completed exonerate analysis\n"
  assert_true(align_exonerate_test_raises(text))
}

///|
test "Bio.Align.exonerate rejects query span mismatch" {
  let text = "Command line: [cmd]\nHostname: [host]\n" +
    "vulgar: q 0 4 + t 0 3 + 1 M 3 3\n" +
    "-- completed exonerate analysis\n"
  assert_true(align_exonerate_test_raises(text))
}

///|
test "Bio.Align.exonerate rejects target span mismatch" {
  let text = "Command line: [cmd]\nHostname: [host]\n" +
    "vulgar: q 0 3 + t 0 4 + 1 M 3 3\n" +
    "-- completed exonerate analysis\n"
  assert_true(align_exonerate_test_raises(text))
}

///|
test "Bio.Align.exonerate rejects unknown vulgar operation" {
  let text = "Command line: [cmd]\nHostname: [host]\n" +
    "vulgar: q 0 3 + t 0 3 + 1 X 3 3\n" +
    "-- completed exonerate analysis\n"
  assert_true(align_exonerate_test_raises(text))
}

///|
test "Bio.Align.exonerate rejects empty vulgar operation" {
  let text = "Command line: [cmd]\nHostname: [host]\n" +
    "vulgar: q 0 1 + t 0 1 + 1 M 0 0\n" +
    "-- completed exonerate analysis\n"
  assert_true(align_exonerate_test_raises(text))
}

///|
test "Bio.Align.exonerate rejects two-axis vulgar gap" {
  let text = "Command line: [cmd]\nHostname: [host]\n" +
    "vulgar: q 0 1 + t 0 1 + 1 G 1 1\n" +
    "-- completed exonerate analysis\n"
  assert_true(align_exonerate_test_raises(text))
}

///|
test "Bio.Align.exonerate rejects empty vulgar gap" {
  let text = "Command line: [cmd]\nHostname: [host]\n" +
    "vulgar: q 0 1 + t 0 1 + 1 G 0 0\n" +
    "-- completed exonerate analysis\n"
  assert_true(align_exonerate_test_raises(text))
}

///|
test "Bio.Align.exonerate rejects invalid splice site" {
  let text = "Command line: [cmd]\nHostname: [host]\n" +
    "vulgar: q 0 1 + t 0 1 + 1 5 1 1\n" +
    "-- completed exonerate analysis\n"
  assert_true(align_exonerate_test_raises(text))
}

///|
test "Bio.Align.exonerate rejects invalid codon step" {
  let text = "Command line: [cmd]\nHostname: [host]\n" +
    "vulgar: q 0 2 + t 0 2 + 1 C 2 2\n" +
    "-- completed exonerate analysis\n"
  assert_true(align_exonerate_test_raises(text))
}

///|
test "Bio.Align.exonerate rejects two-axis intron" {
  let text = "Command line: [cmd]\nHostname: [host]\n" +
    "vulgar: q 0 1 + t 0 1 + 1 I 1 1\n" +
    "-- completed exonerate analysis\n"
  assert_true(align_exonerate_test_raises(text))
}

///|
test "Bio.Align.exonerate rejects empty non-equivalenced region" {
  let text = "Command line: [cmd]\nHostname: [host]\n" +
    "vulgar: q 0 1 + t 0 1 + 1 N 0 0\n" +
    "-- completed exonerate analysis\n"
  assert_true(align_exonerate_test_raises(text))
}

///|
test "Bio.Align.exonerate rejects unknown cigar operation" {
  let text = "Command line: [cmd]\nHostname: [host]\n" +
    "cigar: q 0 3 + t 0 3 + 1 X 3\n" +
    "-- completed exonerate analysis\n"
  assert_true(align_exonerate_test_raises(text))
}

///|
test "Bio.Align.exonerate rejects zero cigar step" {
  let text = "Command line: [cmd]\nHostname: [host]\n" +
    "cigar: q 0 1 + t 0 1 + 1 M 0\n" +
    "-- completed exonerate analysis\n"
  assert_true(align_exonerate_test_raises(text))
}

///|
test "Bio.Align.exonerate rejects header-only alignment line" {
  let text = "Command line: [cmd]\nHostname: [host]\n" +
    "cigar: q 0 1 + t 0 1 + 1\n" +
    "-- completed exonerate analysis\n"
  assert_true(align_exonerate_test_raises(text))
}

///|
test "Bio.Align.exonerate metadata rejects newline" {
  let raised = try {
    ignore(@src.AlignExonerateMetadata::create(hostname="bad\nhost"))
    false
  } catch {
    AlignExonerateError(_) => true
  }
  assert_true(raised)
}

///|
test "Bio.Align.exonerate operation rejects negative step" {
  let raised = try {
    ignore(@src.AlignExonerateOperation::create("M", -1, 1))
    false
  } catch {
    AlignExonerateError(_) => true
  }
  assert_true(raised)
}

///|
test "Bio.Align.exonerate operation rejects unknown code" {
  let raised = try {
    ignore(@src.AlignExonerateOperation::create("Q", 1, 1))
    false
  } catch {
    AlignExonerateError(_) => true
  }
  assert_true(raised)
}

///|
test "Bio.Align.exonerate alignment rejects whitespace identifier" {
  let operation = @src.AlignExonerateOperation::create("M", 1, 1) catch {
    AlignExonerateError(message) =>
      abort("operation construction failed: " + message)
  }
  let raised = try {
    ignore(
      @src.AlignExonerateAlignment::create(
        "bad query",
        0,
        1,
        @src.AlignExonerateForward,
        "target",
        0,
        1,
        @src.AlignExonerateForward,
        1.0,
        [operation],
      ),
    )
    false
  } catch {
    AlignExonerateError(_) => true
  }
  assert_true(raised)
}

///|
test "Bio.Align.exonerate alignment rejects no operations" {
  let raised = try {
    ignore(
      @src.AlignExonerateAlignment::create(
        "q",
        0,
        1,
        @src.AlignExonerateForward,
        "t",
        0,
        1,
        @src.AlignExonerateForward,
        1.0,
        [],
      ),
    )
    false
  } catch {
    AlignExonerateError(_) => true
  }
  assert_true(raised)
}

///|
test "Bio.Align.exonerate rejects unsupported writer format" {
  assert_true(
    align_exonerate_test_write_raises(align_exonerate_test_sample(), "text"),
  )
}

///|
test "Bio.Align.exonerate cigar writer rejects unequal plain match" {
  let operation = @src.AlignExonerateOperation::create("M", 1, 2) catch {
    AlignExonerateError(message) =>
      abort("operation construction failed: " + message)
  }
  let alignment = @src.AlignExonerateAlignment::create(
    "q",
    0,
    1,
    @src.AlignExonerateForward,
    "t",
    0,
    2,
    @src.AlignExonerateForward,
    1.0,
    [operation],
  ) catch {
    AlignExonerateError(message) =>
      abort("alignment construction failed: " + message)
  }
  let metadata = @src.AlignExonerateMetadata::create() catch {
    AlignExonerateError(message) =>
      abort("metadata construction failed: " + message)
  }
  let document = @src.AlignExonerateDocument::create(metadata, [alignment]) catch {
    AlignExonerateError(message) =>
      abort("document construction failed: " + message)
  }
  assert_true(align_exonerate_test_write_raises(document, "cigar"))
}
