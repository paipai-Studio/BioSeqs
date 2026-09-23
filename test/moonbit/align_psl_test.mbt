///|
fn psl_test_join(values : Array[String], separator : String) -> String {
  let output = StringBuilder()
  for index in 0..<values.length() {
    if index > 0 {
      output.write_string(separator)
    }
    output.write_string(values[index])
  }
  output.to_string()
}

///|
fn psl_test_fields() -> Array[String] {
  [
    "6", "1", "0", "1", "1", "2", "1", "3", "+", "query1", "20", "1", "11", "chr1",
    "30", "2", "13", "2", "4,4,", "1,7,", "2,9,",
  ]
}

///|
fn psl_test_line(fields : Array[String]) -> String {
  psl_test_join(fields, "\t") + "\n"
}

///|
fn psl_test_parse_raises(content : String) -> Bool {
  try {
    ignore(@src.psl_parse(content))
    false
  } catch {
    PslError(_) => true
  }
}

///|
fn psl_test_create_raises(
  target_coordinates : Array[Int],
  query_coordinates : Array[Int],
) -> Bool {
  try {
    ignore(
      @src.PslAlignment::create(
        "target", 30, "query", 20, target_coordinates, query_coordinates,
      ),
    )
    false
  } catch {
    PslError(_) => true
  }
}

///|
fn psl_test_forward() -> @src.PslAlignment {
  @src.PslAlignment::create(
    "chr1",
    30,
    "query1",
    20,
    [2, 6, 9, 9, 13],
    [1, 5, 5, 7, 11],
    matches=6,
    mismatches=1,
    repeat_matches=0,
    n_count=1,
  ) catch {
    PslError(message) => abort("valid forward PSL failed: " + message)
  }
}

///|
fn psl_test_reverse() -> @src.PslAlignment {
  @src.PslAlignment::create("chrR", 30, "queryR", 10, [4, 8, 10, 14], [
    8, 4, 4, 0,
  ]) catch {
    PslError(message) => abort("valid reverse PSL failed: " + message)
  }
}

///|
fn psl_test_protein(reverse : Bool) -> @src.PslAlignment {
  if reverse {
    @src.PslAlignment::create(
      "codingReverse",
      18,
      "peptide",
      4,
      [18, 12, 9, 3],
      [0, 2, 2, 4],
      sequence_type=@src.PslProtein,
      target_sequence="AAATCCAAAAAAAGCCAT",
      query_sequence="MAFG",
    ) catch {
      PslError(message) =>
        abort("valid reverse translated PSL failed: " + message)
    }
  } else {
    @src.PslAlignment::create(
      "codingForward",
      18,
      "peptide",
      4,
      [0, 6, 9, 15],
      [0, 2, 2, 4],
      sequence_type=@src.PslProtein,
      target_sequence="ATGGCTAAATTTGGATAG",
      query_sequence="MAFG",
    ) catch {
      PslError(message) =>
        abort("valid forward translated PSL failed: " + message)
    }
  }
}

///|
test "Bio.Align.psl constructor: forward path and metadata" {
  let alignment = psl_test_forward()
  assert_eq(alignment.target_name, "chr1")
  assert_eq(alignment.target_size, 30)
  assert_eq(alignment.query_name, "query1")
  assert_eq(alignment.query_size, 20)
  assert_eq(alignment.target_coordinates, [2, 6, 9, 9, 13])
  assert_eq(alignment.query_coordinates, [1, 5, 5, 7, 11])
  assert_true(alignment.sequence_type is @src.PslSequenceType::PslNucleotide)
}

///|
test "Bio.Align.psl constructor: omitted categories become matches" {
  let alignment = @src.PslAlignment::create("target", 10, "query", 10, [1, 5], [
    2, 6,
  ]) catch {
    _ => abort("valid default categories should construct")
  }
  assert_eq(alignment.matches, 4)
  assert_eq(alignment.mismatches, 0)
  assert_eq(alignment.repeat_matches, 0)
  assert_eq(alignment.n_count, 0)
}

///|
test "Bio.Align.psl constructor: defensive coordinate copy" {
  let target = [1, 5]
  let query = [2, 6]
  let alignment = @src.PslAlignment::create(
    "target", 10, "query", 10, target, query,
  ) catch {
    _ => abort("valid alignment should construct")
  }
  target[0] = 0
  query[0] = 0
  assert_eq(alignment.target_coordinates, [1, 5])
  assert_eq(alignment.query_coordinates, [2, 6])
}

///|
test "Bio.Align.psl constructor: decreasing target is normalized" {
  let alignment = @src.PslAlignment::create(
    "target",
    20,
    "query",
    20,
    [10, 6],
    [8, 4],
  ) catch {
    _ => abort("decreasing nucleotide target should normalize")
  }
  assert_eq(alignment.target_coordinates, [6, 10])
  assert_eq(alignment.query_coordinates, [4, 8])
  assert_false(alignment.is_reverse())
}

///|
test "Bio.Align.psl constructor: decreasing protein query is normalized" {
  let alignment = @src.PslAlignment::create(
    "target",
    12,
    "query",
    4,
    [0, 6],
    [2, 0],
    sequence_type=@src.PslProtein,
  ) catch {
    _ => abort("decreasing protein query should normalize")
  }
  assert_eq(alignment.target_coordinates, [6, 0])
  assert_eq(alignment.query_coordinates, [0, 2])
  assert_true(alignment.is_reverse())
}

///|
test "Bio.Align.psl constructor: empty target name is rejected" {
  let raised = try {
    ignore(@src.PslAlignment::create("", 10, "q", 10, [0, 2], [0, 2]))
    false
  } catch {
    PslError(_) => true
  }
  assert_true(raised)
}

///|
test "Bio.Align.psl constructor: control characters are rejected" {
  let raised = try {
    ignore(
      @src.PslAlignment::create("target", 10, "bad\tquery", 10, [0, 2], [0, 2]),
    )
    false
  } catch {
    PslError(_) => true
  }
  assert_true(raised)
}

///|
test "Bio.Align.psl constructor: non-positive sequence size is rejected" {
  let raised = try {
    ignore(@src.PslAlignment::create("t", 0, "q", 10, [0, 2], [0, 2]))
    false
  } catch {
    PslError(_) => true
  }
  assert_true(raised)
}

///|
test "Bio.Align.psl constructor: unequal coordinate arrays are rejected" {
  assert_true(psl_test_create_raises([0, 2], [0, 2, 3]))
}

///|
test "Bio.Align.psl constructor: one coordinate point is rejected" {
  assert_true(psl_test_create_raises([0], [0]))
}

///|
test "Bio.Align.psl constructor: target direction changes are rejected" {
  assert_true(psl_test_create_raises([0, 4, 2], [0, 4, 6]))
}

///|
test "Bio.Align.psl constructor: query direction changes are rejected" {
  assert_true(psl_test_create_raises([0, 4, 6], [0, 4, 2]))
}

///|
test "Bio.Align.psl constructor: out-of-bounds target is rejected" {
  assert_true(psl_test_create_raises([29, 31], [0, 2]))
}

///|
test "Bio.Align.psl constructor: out-of-bounds query is rejected" {
  assert_true(psl_test_create_raises([0, 2], [19, 21]))
}

///|
test "Bio.Align.psl constructor: zero-length path step is rejected" {
  assert_true(psl_test_create_raises([0, 0, 2], [0, 0, 2]))
}

///|
test "Bio.Align.psl constructor: nucleotide aligned ratio is enforced" {
  assert_true(psl_test_create_raises([0, 4], [0, 3]))
}

///|
test "Bio.Align.psl constructor: translated aligned ratio is enforced" {
  let raised = try {
    ignore(
      @src.PslAlignment::create(
        "t",
        20,
        "q",
        10,
        [0, 5],
        [0, 2],
        sequence_type=@src.PslProtein,
      ),
    )
    false
  } catch {
    PslError(_) => true
  }
  assert_true(raised)
}

///|
test "Bio.Align.psl constructor: gap-only path is rejected" {
  assert_true(psl_test_create_raises([0, 0, 4], [0, 2, 2]))
}

///|
test "Bio.Align.psl constructor: leading gap is rejected" {
  assert_true(psl_test_create_raises([0, 0, 4], [0, 2, 6]))
}

///|
test "Bio.Align.psl constructor: trailing gap is rejected" {
  assert_true(psl_test_create_raises([0, 4, 4], [0, 4, 6]))
}

///|
test "Bio.Align.psl constructor: sequence length is validated" {
  let raised = try {
    ignore(
      @src.PslAlignment::create(
        "t",
        5,
        "q",
        5,
        [0, 5],
        [0, 5],
        target_sequence="ACGT",
      ),
    )
    false
  } catch {
    PslError(_) => true
  }
  assert_true(raised)
}

///|
test "Bio.Align.psl constructor: gapped sequence is rejected" {
  let raised = try {
    ignore(
      @src.PslAlignment::create(
        "t",
        5,
        "q",
        5,
        [0, 5],
        [0, 5],
        query_sequence="AC-TG",
      ),
    )
    false
  } catch {
    PslError(_) => true
  }
  assert_true(raised)
}

///|
test "Bio.Align.psl constructor: both PSLX sequence arrays are required" {
  let raised = try {
    ignore(
      @src.PslAlignment::create("t", 5, "q", 5, [0, 5], [0, 5], query_block_sequences=[
        "ACGTA",
      ]),
    )
    false
  } catch {
    PslError(_) => true
  }
  assert_true(raised)
}

///|
test "Bio.Align.psl constructor: PSLX block count is validated" {
  let raised = try {
    ignore(
      @src.PslAlignment::create(
        "t",
        12,
        "q",
        12,
        [0, 4, 8, 12],
        [0, 4, 4, 8],
        query_block_sequences=["AAAA"],
        target_block_sequences=["AAAA"],
      ),
    )
    false
  } catch {
    PslError(_) => true
  }
  assert_true(raised)
}

///|
test "Bio.Align.psl constructor: PSLX fragment length is validated" {
  let raised = try {
    ignore(
      @src.PslAlignment::create(
        "t",
        5,
        "q",
        5,
        [0, 5],
        [0, 5],
        query_block_sequences=["ACGT"],
        target_block_sequences=["ACGTA"],
      ),
    )
    false
  } catch {
    PslError(_) => true
  }
  assert_true(raised)
}

///|
test "Bio.Align.psl constructor: PSLX fragment delimiter is rejected" {
  let raised = try {
    ignore(
      @src.PslAlignment::create(
        "t",
        5,
        "q",
        5,
        [0, 5],
        [0, 5],
        query_block_sequences=["AC,TA"],
        target_block_sequences=["ACGTA"],
      ),
    )
    false
  } catch {
    PslError(_) => true
  }
  assert_true(raised)
}

///|
test "Bio.Align.psl constructor: partial categories are rejected" {
  let raised = try {
    ignore(@src.PslAlignment::create("t", 5, "q", 5, [0, 5], [0, 5], matches=5))
    false
  } catch {
    PslError(_) => true
  }
  assert_true(raised)
}

///|
test "Bio.Align.psl constructor: category sum is validated" {
  let raised = try {
    ignore(
      @src.PslAlignment::create(
        "t",
        5,
        "q",
        5,
        [0, 5],
        [0, 5],
        matches=4,
        mismatches=0,
        repeat_matches=0,
        n_count=0,
      ),
    )
    false
  } catch {
    PslError(_) => true
  }
  assert_true(raised)
}

///|
test "Bio.Align.psl blocks: aligned segments preserve orientation" {
  let blocks = psl_test_reverse().blocks()
  assert_eq(blocks.length(), 2)
  assert_eq(blocks[0].target_start, 4)
  assert_eq(blocks[0].target_end, 8)
  assert_eq(blocks[0].query_start, 8)
  assert_eq(blocks[0].query_end, 4)
  assert_eq(blocks[0].target_size, 4)
  assert_eq(blocks[0].query_size, 4)
}

///|
test "Bio.Align.psl counts: query and target gaps are independent" {
  let counts = psl_test_forward().counts()
  assert_eq(counts.aligned_query_units, 8)
  assert_eq(counts.query_insert_count, 1)
  assert_eq(counts.query_insert_bases, 2)
  assert_eq(counts.target_insert_count, 1)
  assert_eq(counts.target_insert_bases, 3)
  assert_eq(counts.block_count, 2)
}

///|
test "Bio.Align.psl metrics: identity and score" {
  let alignment = psl_test_forward()
  assert_eq(alignment.identity(), 0.75)
  assert_eq(alignment.score(), 3)
}

///|
test "Bio.Align.psl mapping: forward target to query" {
  let alignment = psl_test_forward()
  assert_eq(alignment.target_to_query(2), Some(1))
  assert_eq(alignment.target_to_query(5), Some(4))
  assert_eq(alignment.target_to_query(6), None)
  assert_eq(alignment.target_to_query(9), Some(7))
  assert_eq(alignment.target_to_query(-1), None)
}

///|
test "Bio.Align.psl mapping: forward query to target" {
  let alignment = psl_test_forward()
  assert_eq(alignment.query_to_target(1), Some(2))
  assert_eq(alignment.query_to_target(5), None)
  assert_eq(alignment.query_to_target(7), Some(9))
  assert_eq(alignment.query_to_target_interval(10), Some((12, 13)))
  assert_eq(alignment.query_to_target(20), None)
}

///|
test "Bio.Align.psl mapping: reverse nucleotide coordinates" {
  let alignment = psl_test_reverse()
  assert_eq(alignment.target_to_query(4), Some(7))
  assert_eq(alignment.target_to_query(7), Some(4))
  assert_eq(alignment.target_to_query(8), None)
  assert_eq(alignment.query_to_target(7), Some(4))
  assert_eq(alignment.query_to_target(4), Some(7))
  assert_eq(alignment.query_to_target(3), Some(10))
}

///|
test "Bio.Align.psl mapping: translated codon intervals" {
  let alignment = psl_test_protein(false)
  assert_eq(alignment.query_to_target_interval(0), Some((0, 3)))
  assert_eq(alignment.query_to_target_interval(2), Some((9, 12)))
  assert_eq(alignment.target_to_query(10), Some(2))
  assert_eq(alignment.target_to_query(12), Some(3))
  assert_eq(alignment.target_to_query(7), None)
}

///|
test "Bio.Align.psl mapping: reverse translated codon intervals" {
  let alignment = psl_test_protein(true)
  assert_eq(alignment.query_to_target_interval(0), Some((15, 18)))
  assert_eq(alignment.query_to_target_interval(2), Some((6, 9)))
  assert_eq(alignment.target_to_query(17), Some(0))
  assert_eq(alignment.target_to_query(14), Some(1))
  assert_eq(alignment.target_to_query(10), None)
}

///|
test "Bio.Align.psl recount: mismatch and wildcard" {
  let alignment = @src.PslAlignment::create(
    "target",
    5,
    "query",
    5,
    [0, 5],
    [0, 5],
    target_sequence="ACGTN",
    query_sequence="ACGCA",
  ) catch {
    _ => abort("valid sequence alignment should construct")
  }
  let recounted = alignment.recount() catch {
    _ => abort("sequence-aware recount should succeed")
  }
  assert_eq(recounted.matches, 3)
  assert_eq(recounted.mismatches, 1)
  assert_eq(recounted.repeat_matches, 0)
  assert_eq(recounted.n_count, 1)
}

///|
test "Bio.Align.psl recount: lowercase target repeat masking" {
  let alignment = @src.PslAlignment::create(
    "target",
    5,
    "query",
    5,
    [0, 5],
    [0, 5],
    target_sequence="AcGTN",
    query_sequence="ACGTA",
  ) catch {
    _ => abort("valid sequence alignment should construct")
  }
  let recounted = alignment.recount(mask=@src.PslMaskLower) catch {
    _ => abort("lowercase repeat recount should succeed")
  }
  assert_eq(recounted.matches, 3)
  assert_eq(recounted.repeat_matches, 1)
  assert_eq(recounted.n_count, 1)
}

///|
test "Bio.Align.psl recount: uppercase target repeat masking" {
  let alignment = @src.PslAlignment::create(
    "target",
    5,
    "query",
    5,
    [0, 5],
    [0, 5],
    target_sequence="AcGTN",
    query_sequence="ACGTA",
  ) catch {
    _ => abort("valid sequence alignment should construct")
  }
  let recounted = alignment.recount(mask=@src.PslMaskUpper) catch {
    _ => abort("uppercase repeat recount should succeed")
  }
  assert_eq(recounted.matches, 1)
  assert_eq(recounted.repeat_matches, 3)
  assert_eq(recounted.n_count, 1)
}

///|
test "Bio.Align.psl recount: custom wildcard" {
  let alignment = @src.PslAlignment::create(
    "target",
    4,
    "query",
    4,
    [0, 4],
    [0, 4],
    target_sequence="ACXT",
    query_sequence="ACGT",
  ) catch {
    _ => abort("valid sequence alignment should construct")
  }
  let recounted = alignment.recount(wildcard='X') catch {
    _ => abort("custom wildcard recount should succeed")
  }
  assert_eq(recounted.matches, 3)
  assert_eq(recounted.n_count, 1)
}

///|
test "Bio.Align.psl recount: reverse complement query" {
  let alignment = @src.PslAlignment::create(
    "target",
    4,
    "query",
    4,
    [0, 4],
    [4, 0],
    target_sequence="AGTC",
    query_sequence="GACT",
  ) catch {
    _ => abort("valid reverse sequence alignment should construct")
  }
  let recounted = alignment.recount() catch {
    _ => abort("reverse-complement recount should succeed")
  }
  assert_eq(recounted.matches, 4)
  assert_eq(recounted.mismatches, 0)
}

///|
test "Bio.Align.psl recount: translated forward DNA" {
  let recounted = psl_test_protein(false).recount() catch {
    _ => abort("forward translated recount should succeed")
  }
  assert_eq(recounted.matches, 4)
  assert_eq(recounted.mismatches, 0)
}

///|
test "Bio.Align.psl recount: translated reverse-complement DNA" {
  let recounted = psl_test_protein(true).recount() catch {
    _ => abort("reverse translated recount should succeed")
  }
  assert_eq(recounted.matches, 4)
  assert_eq(recounted.mismatches, 0)
}

///|
test "Bio.Align.psl recount: concrete sequences are required" {
  let raised = try {
    ignore(psl_test_forward().recount())
    false
  } catch {
    PslError(_) => true
  }
  assert_true(raised)
}

///|
test "Bio.Align.psl recount: protein repeat mask is rejected" {
  let raised = try {
    ignore(psl_test_protein(false).recount(mask=@src.PslMaskLower))
    false
  } catch {
    PslError(_) => true
  }
  assert_true(raised)
}

///|
test "Bio.Align.psl writer config: defaults" {
  let config = @src.PslWriteConfig::default()
  assert_true(config.header)
  assert_false(config.pslx)
  assert_false(config.recount)
  assert_true(config.mask is @src.PslMaskMode::PslNoMask)
  assert_eq(config.wildcard, 'N')
}

///|
test "Bio.Align.psl writer config: custom options" {
  let config = @src.PslWriteConfig::create(
    header=false,
    pslx=true,
    recount=true,
    mask=@src.PslMaskLower,
    wildcard='X',
  ) catch {
    _ => abort("valid writer config should construct")
  }
  assert_false(config.header)
  assert_true(config.pslx)
  assert_true(config.recount)
  assert_true(config.mask is @src.PslMaskMode::PslMaskLower)
  assert_eq(config.wildcard, 'X')
}

///|
test "Bio.Align.psl writer config: delimiter wildcard is rejected" {
  let raised = try {
    ignore(@src.PslWriteConfig::create(wildcard=','))
    false
  } catch {
    PslError(_) => true
  }
  assert_true(raised)
}

///|
test "Bio.Align.psl header: version and column labels" {
  let header = @src.psl_header(version="4") catch {
    _ => abort("valid PSL header should format")
  }
  assert_true(header.has_prefix("psLayout version 4\n\n"))
  assert_true(header.contains("blockSizes"))
}

///|
test "Bio.Align.psl header: whitespace version is rejected" {
  let raised = try {
    ignore(@src.psl_header(version="3 draft"))
    false
  } catch {
    PslError(_) => true
  }
  assert_true(raised)
}

///|
test "Bio.Align.psl format: exact forward 21-column record" {
  let config = @src.PslWriteConfig::create(header=false) catch {
    _ => abort("valid writer config should construct")
  }
  let line = psl_test_forward().format(config~) catch {
    _ => abort("valid PSL should format")
  }
  assert_eq(line, psl_test_line(psl_test_fields()))
}

///|
test "Bio.Align.psl format: reverse nucleotide storage coordinates" {
  let config = @src.PslWriteConfig::create(header=false) catch {
    _ => abort("valid writer config should construct")
  }
  let line = psl_test_reverse().format(config~) catch {
    _ => abort("valid reverse PSL should format")
  }
  assert_true(line.contains("\t-\tqueryR\t10\t0\t8\tchrR\t30\t4\t14\t"))
  assert_true(line.has_suffix("\t4,4,\t2,6,\t4,10,\n"))
}

///|
test "Bio.Align.psl format: translated strand forms" {
  let config = @src.PslWriteConfig::create(header=false) catch {
    _ => abort("valid writer config should construct")
  }
  let forward = psl_test_protein(false).format(config~) catch {
    _ => abort("forward protein PSL should format")
  }
  let reverse = psl_test_protein(true).format(config~) catch {
    _ => abort("reverse protein PSL should format")
  }
  assert_true(forward.contains("\t++\tpeptide\t"))
  assert_true(reverse.contains("\t+-\tpeptide\t"))
  assert_true(reverse.has_suffix("\t2,2,\t0,2,\t0,9,\n"))
}

///|
test "Bio.Align.psl format: generated nucleotide PSLX blocks" {
  let alignment = @src.PslAlignment::create(
    "target",
    8,
    "query",
    8,
    [0, 4, 8],
    [8, 4, 0],
    target_sequence="ACGTTTAA",
    query_sequence="TTAAAACG",
  ) catch {
    _ => abort("valid reverse PSLX source should construct")
  }
  let config = @src.PslWriteConfig::create(header=false, pslx=true) catch {
    _ => abort("valid PSLX config should construct")
  }
  let line = alignment.format(config~) catch {
    _ => abort("PSLX output should derive block sequences")
  }
  assert_true(line.has_suffix("\tCGTT,TTAA,\tACGT,TTAA,\n"))
}

///|
test "Bio.Align.psl format: generated translated PSLX blocks" {
  let config = @src.PslWriteConfig::create(header=false, pslx=true) catch {
    _ => abort("valid PSLX config should construct")
  }
  let line = psl_test_protein(false).format(config~) catch {
    _ => abort("translated PSLX output should derive amino-acid blocks")
  }
  assert_true(line.has_suffix("\tMA,FG,\tMA,FG,\n"))
}

///|
test "Bio.Align.psl format: PSLX needs sequence data" {
  let config = @src.PslWriteConfig::create(header=false, pslx=true) catch {
    _ => abort("valid PSLX config should construct")
  }
  let raised = try {
    ignore(psl_test_forward().format(config~))
    false
  } catch {
    PslError(_) => true
  }
  assert_true(raised)
}

///|
test "Bio.Align.psl format: writer can recount sequences" {
  let alignment = @src.PslAlignment::create(
    "target",
    4,
    "query",
    4,
    [0, 4],
    [0, 4],
    target_sequence="ACGT",
    query_sequence="ACGA",
  ) catch {
    _ => abort("valid alignment should construct")
  }
  let config = @src.PslWriteConfig::create(header=false, recount=true) catch {
    _ => abort("valid recount config should construct")
  }
  let line = alignment.format(config~) catch {
    _ => abort("recounting writer should succeed")
  }
  assert_true(line.has_prefix("3\t1\t0\t0\t"))
}

///|
test "Bio.Align.psl parse: headerless 21-column record" {
  let document = @src.psl_parse(psl_test_line(psl_test_fields())) catch {
    _ => abort("valid headerless PSL should parse")
  }
  assert_false(document.has_header)
  assert_eq(document.version, "")
  assert_eq(document.alignments.length(), 1)
  assert_eq(document.alignments[0].target_coordinates, [2, 6, 9, 9, 13])
  assert_eq(document.alignments[0].query_coordinates, [1, 5, 5, 7, 11])
}

///|
test "Bio.Align.psl parse: generated header and CRLF" {
  let text = @src.psl_write([psl_test_forward()]) catch {
    _ => abort("valid PSL document should write")
  }
  let crlf = text.replace_all(old="\n", new="\r\n")
  let document = @src.psl_parse(crlf) catch {
    _ => abort("CRLF PSL should parse")
  }
  assert_true(document.has_header)
  assert_eq(document.version, "3")
  assert_eq(document.alignments.length(), 1)
}

///|
test "Bio.Align.psl parse: 23-column PSLX record" {
  let fields = psl_test_fields()
  fields.push("ACGT,NNAA,")
  fields.push("ACGT,NNTA,")
  let document = @src.psl_parse(psl_test_line(fields)) catch {
    _ => abort("valid PSLX should parse")
  }
  let alignment = document.alignments[0]
  assert_true(alignment.is_pslx())
  assert_eq(alignment.query_block_sequences, ["ACGT", "NNAA"])
  assert_eq(alignment.target_block_sequences, ["ACGT", "NNTA"])
  assert_eq(alignment.blocks()[1].query_sequence, "NNAA")
}

///|
test "Bio.Align.psl parse: reverse nucleotide path" {
  let config = @src.PslWriteConfig::create(header=false) catch {
    _ => abort("valid config should construct")
  }
  let line = psl_test_reverse().format(config~) catch {
    _ => abort("reverse PSL should format")
  }
  let parsed = @src.psl_parse(line) catch {
    _ => abort("reverse PSL should parse")
  }
  assert_eq(parsed.alignments[0].query_coordinates, [8, 4, 4, 0])
  assert_true(parsed.alignments[0].is_reverse())
  let biopython_fixture = "175\t0\t6\t0\t0\t0\t2\t1530\t-\tNR_046654.1\t181\t0\t181\tchr3\t198295559\t42530895\t42532606\t3\t63,75,43,\t0,63,138,\t42530895,42532020,42532563,\n"
  let fixture = @src.psl_parse(biopython_fixture) catch {
    _ => abort("Biopython negative-strand fixture should parse")
  }
  let fixture_alignment = fixture.alignments[0]
  assert_eq(fixture_alignment.target_coordinates, [
    42530895, 42530958, 42532020, 42532095, 42532563, 42532606,
  ])
  assert_eq(fixture_alignment.query_coordinates, [181, 118, 118, 43, 43, 0])
  assert_eq(fixture_alignment.counts().target_insert_bases, 1530)
  assert_eq(
    fixture_alignment.format(config~) catch {
      _ => abort("Biopython negative-strand fixture should format")
    },
    biopython_fixture,
  )
}

///|
test "Bio.Align.psl parse: translated forward and reverse paths" {
  let config = @src.PslWriteConfig::create(header=false) catch {
    _ => abort("valid config should construct")
  }
  let forward = psl_test_protein(false).format(config~) catch {
    _ => abort("forward protein PSL should format")
  }
  let reverse = psl_test_protein(true).format(config~) catch {
    _ => abort("reverse protein PSL should format")
  }
  let text = forward + reverse
  let document = @src.psl_parse(text) catch {
    _ => abort("translated PSL records should parse")
  }
  assert_eq(document.alignments.length(), 2)
  assert_eq(document.alignments[0].target_coordinates, [0, 6, 9, 15])
  assert_eq(document.alignments[1].target_coordinates, [18, 12, 9, 3])
}

///|
test "Bio.Align.psl parse: empty input is rejected" {
  assert_true(psl_test_parse_raises("\n\r\n"))
}

///|
test "Bio.Align.psl parse: malformed header is rejected" {
  assert_true(psl_test_parse_raises("psLayout bad 3\n-----\n"))
}

///|
test "Bio.Align.psl parse: missing header separator is rejected" {
  assert_true(psl_test_parse_raises("psLayout version 3\ncolumns only\n"))
}

///|
test "Bio.Align.psl parse: wrong column count is rejected" {
  assert_true(psl_test_parse_raises("1\t2\t3\n"))
}

///|
test "Bio.Align.psl parse: integer trailing garbage is rejected" {
  let fields = psl_test_fields()
  fields[0] = "6oops"
  assert_true(psl_test_parse_raises(psl_test_line(fields)))
}

///|
test "Bio.Align.psl parse: integer overflow is rejected" {
  let fields = psl_test_fields()
  fields[10] = "2147483648"
  assert_true(psl_test_parse_raises(psl_test_line(fields)))
}

///|
test "Bio.Align.psl parse: unsupported strand is rejected" {
  let fields = psl_test_fields()
  fields[8] = "--"
  assert_true(psl_test_parse_raises(psl_test_line(fields)))
}

///|
test "Bio.Align.psl parse: zero block count is rejected" {
  let fields = psl_test_fields()
  fields[17] = "0"
  assert_true(psl_test_parse_raises(psl_test_line(fields)))
}

///|
test "Bio.Align.psl parse: block array count is validated" {
  let fields = psl_test_fields()
  fields[18] = "4,"
  assert_true(psl_test_parse_raises(psl_test_line(fields)))
}

///|
test "Bio.Align.psl parse: empty CSV element is rejected" {
  let fields = psl_test_fields()
  fields[18] = "4,,4,"
  assert_true(psl_test_parse_raises(psl_test_line(fields)))
}

///|
test "Bio.Align.psl parse: zero block size is rejected" {
  let fields = psl_test_fields()
  fields[18] = "4,0,"
  assert_true(psl_test_parse_raises(psl_test_line(fields)))
}

///|
test "Bio.Align.psl parse: overlapping target blocks are rejected" {
  let fields = psl_test_fields()
  fields[20] = "2,5,"
  assert_true(psl_test_parse_raises(psl_test_line(fields)))
}

///|
test "Bio.Align.psl parse: overlapping query blocks are rejected" {
  let fields = psl_test_fields()
  fields[19] = "1,4,"
  assert_true(psl_test_parse_raises(psl_test_line(fields)))
}

///|
test "Bio.Align.psl parse: block beyond target size is rejected" {
  let fields = psl_test_fields()
  fields[14] = "12"
  assert_true(psl_test_parse_raises(psl_test_line(fields)))
}

///|
test "Bio.Align.psl parse: match category sum is validated" {
  let fields = psl_test_fields()
  fields[0] = "7"
  assert_true(psl_test_parse_raises(psl_test_line(fields)))
}

///|
test "Bio.Align.psl parse: query insertion statistics are validated" {
  let fields = psl_test_fields()
  fields[4] = "0"
  assert_true(psl_test_parse_raises(psl_test_line(fields)))
}

///|
test "Bio.Align.psl parse: target insertion statistics are validated" {
  let fields = psl_test_fields()
  fields[7] = "2"
  assert_true(psl_test_parse_raises(psl_test_line(fields)))
}

///|
test "Bio.Align.psl parse: query interval is validated" {
  let fields = psl_test_fields()
  fields[11] = "0"
  assert_true(psl_test_parse_raises(psl_test_line(fields)))
}

///|
test "Bio.Align.psl parse: target interval is validated" {
  let fields = psl_test_fields()
  fields[16] = "14"
  assert_true(psl_test_parse_raises(psl_test_line(fields)))
}

///|
test "Bio.Align.psl parse: PSLX sequence count is validated" {
  let fields = psl_test_fields()
  fields.push("ACGT,")
  fields.push("ACGT,NNTA,")
  assert_true(psl_test_parse_raises(psl_test_line(fields)))
}

///|
test "Bio.Align.psl parse: PSLX fragment length is validated" {
  let fields = psl_test_fields()
  fields.push("ACG,NNAA,")
  fields.push("ACGT,NNTA,")
  assert_true(psl_test_parse_raises(psl_test_line(fields)))
}

///|
test "Bio.Align.psl parse: PSLX invalid residue delimiter is rejected" {
  let fields = psl_test_fields()
  fields.push("AC-T,NNAA,")
  fields.push("ACGT,NNTA,")
  assert_true(psl_test_parse_raises(psl_test_line(fields)))
}

///|
test "Bio.Align.psl document: query and target filters" {
  let config = @src.PslWriteConfig::create(header=false) catch {
    _ => abort("valid config should construct")
  }
  let text = @src.psl_write(
    [psl_test_forward(), psl_test_reverse(), psl_test_protein(false)],
    config~,
  ) catch {
    _ => abort("valid multi-record PSL should write")
  }
  let document = @src.psl_parse(text) catch {
    _ => abort("valid multi-record PSL should parse")
  }
  assert_eq(document.query("peptide").length(), 1)
  assert_eq(document.query("missing").length(), 0)
  assert_eq(document.target("chr1").length(), 1)
}

///|
test "Bio.Align.psl document: summary aggregates path types" {
  let fields = psl_test_fields()
  fields.push("ACGT,NNAA,")
  fields.push("ACGT,NNTA,")
  let config = @src.PslWriteConfig::create(header=false) catch {
    _ => abort("valid config should construct")
  }
  let text = psl_test_line(fields) +
    (psl_test_protein(false).format(config~) catch {
      _ => abort("valid translated PSL should format")
    })
  let summary = (@src.psl_parse(text) catch {
    _ => abort("valid mixed PSL should parse")
  }).summary()
  assert_eq(summary.alignment_count, 2)
  assert_eq(summary.pslx_count, 1)
  assert_eq(summary.nucleotide_count, 1)
  assert_eq(summary.protein_count, 1)
  assert_eq(summary.aligned_query_units, 12)
  assert_eq(summary.query_insert_bases, 2)
  assert_eq(summary.target_insert_bases, 6)
}

///|
test "Bio.Align.psl document: parse-write round trip" {
  let config = @src.PslWriteConfig::create(header=false) catch {
    _ => abort("valid config should construct")
  }
  let initial = @src.psl_write(
    [psl_test_forward(), psl_test_reverse(), psl_test_protein(true)],
    config~,
  ) catch {
    _ => abort("valid PSL records should write")
  }
  let document = @src.psl_parse(initial) catch {
    _ => abort("valid PSL records should parse")
  }
  let rewritten = document.write(config~) catch {
    _ => abort("parsed PSL document should write")
  }
  assert_eq(rewritten, initial)
}

///|
test "Bio.Align.psl document: header-only input is valid" {
  let header = @src.psl_header() catch {
    _ => abort("default PSL header should format")
  }
  let document = @src.psl_parse(header) catch {
    _ => abort("header-only PSL should parse")
  }
  assert_true(document.has_header)
  assert_eq(document.alignments.length(), 0)
}

///|
test "Bio.Align.psl coordinate adapter: forward conversion" {
  let coordinate = psl_test_forward().to_coordinate_alignment() catch {
    _ => abort("forward PSL should adapt to coordinate alignment")
  }
  assert_eq(coordinate.target_name, "chr1")
  assert_eq(coordinate.query_name, "query1")
  assert_eq(coordinate.target_coordinates, [2, 6, 9, 9, 13])
  assert_eq(coordinate.query_coordinates, [1, 5, 5, 7, 11])
}

///|
test "Bio.Align.psl coordinate adapter: reverse nucleotide conversion" {
  let coordinate = psl_test_reverse().to_coordinate_alignment() catch {
    _ => abort("reverse-query PSL should adapt")
  }
  assert_eq(coordinate.query_coordinates, [8, 4, 4, 0])
}

///|
test "Bio.Align.psl coordinate adapter: reverse translated target is rejected" {
  let raised = try {
    ignore(psl_test_protein(true).to_coordinate_alignment())
    false
  } catch {
    PslError(_) => true
  }
  assert_true(raised)
}

///|
test "Bio.Align.psl coordinate adapter: build from coordinate alignment" {
  let coordinate = @src.coordinate_pairwise_alignment_with_lengths(
    "target",
    20,
    "query",
    20,
    [2, 6, 9, 13],
    [1, 5, 5, 9],
  ) catch {
    _ => abort("valid coordinate alignment should construct")
  }
  let alignment = @src.psl_from_coordinate_alignment(coordinate) catch {
    _ => abort("coordinate alignment should convert to PSL")
  }
  assert_eq(alignment.target_coordinates, [2, 6, 9, 13])
  assert_eq(alignment.query_coordinates, [1, 5, 5, 9])
  assert_eq(alignment.counts().target_insert_bases, 3)
}

///|
test "Bio.Align.psl example data: all records construct and summarize" {
  let alignments = @src.psl_example_data() catch {
    PslError(message) => abort("PSL example data failed: " + message)
  }
  assert_eq(alignments.length(), 3)
  assert_false(alignments[0].is_reverse())
  assert_true(alignments[1].is_reverse())
  assert_true(alignments[2].sequence_type is @src.PslSequenceType::PslProtein)
}
