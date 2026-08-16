///|
fn align_tabular_test_join(
  values : Array[String],
  separator : String,
) -> String {
  let mut result = ""
  for index in 0..<values.length() {
    if index > 0 {
      result = result + separator
    }
    result = result + values[index]
  }
  result
}

///|
fn align_tabular_test_block(
  program : String,
  query : String,
  database : String,
  fields : Array[String],
  rows : Array[Array[String]],
  hits : Int,
) -> String {
  let mut result = "# " +
    program +
    "\n# Query: " +
    query +
    "\n# Database: " +
    database +
    "\n"
  if fields.length() > 0 {
    result = result +
      "# Fields: " +
      align_tabular_test_join(fields, ", ") +
      "\n"
  }
  result = result + "# " + hits.to_string() + " hits found\n"
  for row in rows {
    result = result + align_tabular_test_join(row, "\t") + "\n"
  }
  result
}

///|
fn align_tabular_test_fields() -> Array[String] {
  [
    "query id", "subject id", "% identity", "alignment length", "mismatches", "gap opens",
    "q. start", "q. end", "s. start", "s. end", "evalue", "bit score", "BTOP", "query length",
    "subject length", "query seq", "subject seq",
  ]
}

///|
fn align_tabular_test_values() -> Array[String] {
  [
    "q1", "s1", "75.0", "8", "2", "2", "3", "9", "10", "16", "1e-5", "42.0", "2AC1-GT-2",
    "20", "30", "AAAT-TGG", "AACTG-GG",
  ]
}

///|
fn align_tabular_test_document() -> @src.AlignTabularDocument {
  let text = align_tabular_test_block(
    "BLASTP 2.14.1+",
    "q1 protein query",
    "proteins",
    align_tabular_test_fields(),
    [align_tabular_test_values()],
    1,
  )
  @src.align_tabular_parse(text) catch {
    _ => abort("valid alignment-aware BLAST document should parse")
  }
}

///|
fn align_tabular_test_raises(content : String) -> Bool {
  try {
    ignore(@src.align_tabular_parse(content))
    false
  } catch {
    AlignTabularError(_) => true
  }
}

///|
fn align_tabular_test_btop_raises(value : String) -> Bool {
  try {
    ignore(@src.align_tabular_parse_btop(value))
    false
  } catch {
    AlignTabularError(_) => true
  }
}

///|
fn align_tabular_test_cigar_raises(value : String) -> Bool {
  try {
    ignore(@src.align_tabular_parse_cigar(value))
    false
  } catch {
    AlignTabularError(_) => true
  }
}

///|
test "Bio.Align.tabular BTOP: one match run" {
  let (target, query, operations) = @src.align_tabular_parse_btop("12") catch {
    _ => abort("valid BTOP should parse")
  }
  assert_eq(target, [0, 12])
  assert_eq(query, [0, 12])
  assert_eq(operations.length(), 1)
  assert_true(operations[0] is @src.AlignTabularOperation::Aligned)
}

///|
test "Bio.Align.tabular BTOP: mismatches merge with matches" {
  let (target, query, operations) = @src.align_tabular_parse_btop("2AC3GT1") catch {
    _ => abort("valid BTOP should parse")
  }
  assert_eq(target, [0, 8])
  assert_eq(query, [0, 8])
  assert_eq(operations.length(), 1)
}

///|
test "Bio.Align.tabular BTOP: query gap consumes target" {
  let (target, query, operations) = @src.align_tabular_parse_btop("2-G3") catch {
    _ => abort("valid BTOP should parse")
  }
  assert_eq(target, [0, 2, 3, 6])
  assert_eq(query, [0, 2, 2, 5])
  assert_true(operations[1] is @src.AlignTabularOperation::QueryGap)
}

///|
test "Bio.Align.tabular BTOP: target gap consumes query" {
  let (target, query, operations) = @src.align_tabular_parse_btop("2G-3") catch {
    _ => abort("valid BTOP should parse")
  }
  assert_eq(target, [0, 2, 2, 5])
  assert_eq(query, [0, 2, 3, 6])
  assert_true(operations[1] is @src.AlignTabularOperation::TargetGap)
}

///|
test "Bio.Align.tabular BTOP: adjacent gaps form one run" {
  let (target, query, operations) = @src.align_tabular_parse_btop("1-A-C2") catch {
    _ => abort("valid BTOP should parse")
  }
  assert_eq(target, [0, 1, 3, 5])
  assert_eq(query, [0, 1, 1, 3])
  assert_eq(operations.length(), 3)
}

///|
test "Bio.Align.tabular BTOP: stop residues are accepted" {
  let (target, query, operations) = @src.align_tabular_parse_btop("1*A2") catch {
    _ => abort("protein stop mismatch should parse")
  }
  assert_eq(target, [0, 4])
  assert_eq(query, [0, 4])
  assert_eq(operations.length(), 1)
}

///|
test "Bio.Align.tabular CIGAR: one match run" {
  let (target, query, operations) = @src.align_tabular_parse_cigar("218M") catch {
    _ => abort("valid aln_code should parse")
  }
  assert_eq(target, [0, 218])
  assert_eq(query, [0, 218])
  assert_eq(operations.length(), 1)
  assert_true(operations[0] is @src.AlignTabularOperation::Aligned)
}

///|
test "Bio.Align.tabular CIGAR: FASTA insertion consumes target" {
  let (target, query, operations) = @src.align_tabular_parse_cigar("3M2I4M") catch {
    _ => abort("valid aln_code should parse")
  }
  assert_eq(target, [0, 3, 5, 9])
  assert_eq(query, [0, 3, 3, 7])
  assert_true(operations[1] is @src.AlignTabularOperation::QueryGap)
}

///|
test "Bio.Align.tabular CIGAR: FASTA deletion consumes query" {
  let (target, query, operations) = @src.align_tabular_parse_cigar("3M2D4M") catch {
    _ => abort("valid aln_code should parse")
  }
  assert_eq(target, [0, 3, 3, 7])
  assert_eq(query, [0, 3, 5, 9])
  assert_true(operations[1] is @src.AlignTabularOperation::TargetGap)
}

///|
test "Bio.Align.tabular CIGAR: adjacent operations are normalized" {
  let (target, query, operations) = @src.align_tabular_parse_cigar("2M3M1I2I") catch {
    _ => abort("valid aln_code should parse")
  }
  assert_eq(target, [0, 5, 8])
  assert_eq(query, [0, 5, 5])
  assert_eq(operations.length(), 2)
}

///|
test "Bio.Align.tabular BLAST: metadata and descriptions" {
  let document = align_tabular_test_document()
  let query = document.queries[0]
  assert_eq(query.metadata.program, "BLASTP")
  assert_eq(query.metadata.version, "2.14.1+")
  assert_eq(query.metadata.database, "proteins")
  assert_eq(query.query_id, "q1")
  assert_eq(query.query_description, "protein query")
  assert_eq(query.declared_hits, 1)
}

///|
test "Bio.Align.tabular BLAST: intervals are zero-based half-open" {
  let alignment = align_tabular_test_document().queries[0].alignments[0]
  assert_eq(alignment.query_start, Some(2))
  assert_eq(alignment.query_end, Some(9))
  assert_eq(alignment.target_start, Some(9))
  assert_eq(alignment.target_end, Some(16))
  assert_true(alignment.query_strand is @src.Strand::Plus)
  assert_true(alignment.target_strand is @src.Strand::Plus)
}

///|
test "Bio.Align.tabular BLAST: BTOP creates absolute coordinate path" {
  let alignment = align_tabular_test_document().queries[0].alignments[0]
  assert_eq(alignment.target_coordinates, [9, 13, 14, 14, 16])
  assert_eq(alignment.query_coordinates, [2, 6, 6, 7, 9])
  assert_true(alignment.path_is_absolute)
  assert_true(alignment.trace is @src.AlignTabularTrace::Btop)
}

///|
test "Bio.Align.tabular BLAST: path operations preserve gap direction" {
  let alignment = align_tabular_test_document().queries[0].alignments[0]
  assert_eq(alignment.operations.length(), 4)
  assert_true(alignment.operations[0] is @src.AlignTabularOperation::Aligned)
  assert_true(alignment.operations[1] is @src.AlignTabularOperation::QueryGap)
  assert_true(alignment.operations[2] is @src.AlignTabularOperation::TargetGap)
  assert_true(alignment.operations[3] is @src.AlignTabularOperation::Aligned)
}

///|
test "Bio.Align.tabular BLAST: path statistics" {
  let alignment = align_tabular_test_document().queries[0].alignments[0]
  assert_eq(alignment.alignment_columns(), 8)
  assert_eq(alignment.aligned_residues(), 6)
  assert_eq(alignment.query_gap_residues(), 1)
  assert_eq(alignment.target_gap_residues(), 1)
  assert_eq(alignment.gap_events(), 2)
}

///|
test "Bio.Align.tabular BLAST: typed and raw field lookup" {
  let alignment = align_tabular_test_document().queries[0].alignments[0]
  assert_eq(alignment.field("subject id"), Some("s1"))
  assert_eq(
    alignment.integer_field("alignment length") catch {
      _ => abort("integer field should parse")
    },
    Some(8),
  )
  let evalue = alignment.number_field("evalue") catch {
    _ => abort("numeric field should parse")
  }
  match evalue {
    Some(value) => assert_true((value - 1.0e-5).abs() < 1.0e-15)
    None => abort("numeric field should exist")
  }
  assert_eq(alignment.field("missing"), None)
}

///|
test "Bio.Align.tabular BLAST: summary includes score and path statistics" {
  let summary = align_tabular_test_document().queries[0].alignments[0].summary()
  assert_true(summary.contains("q1 -> s1"))
  assert_true(summary.contains("E=1e-5"))
  assert_true(summary.contains("bits=42.0"))
  assert_true(summary.contains("columns=8"))
  assert_true(summary.contains("gaps=2"))
}

///|
test "Bio.Align.tabular BLAST: coordinate alignment conversion" {
  let converted = align_tabular_test_document().queries[0].alignments[0].coordinate_alignment() catch {
    _ => abort("absolute path should convert")
  }
  assert_eq(converted.target_name, "s1")
  assert_eq(converted.query_name, "q1")
  assert_eq(converted.target_length, 30)
  assert_eq(converted.query_length, 20)
  assert_eq(converted.target_coordinates, [9, 13, 14, 14, 16])
  assert_eq(converted.query_coordinates, [2, 6, 6, 7, 9])
}

///|
test "Bio.Align.tabular BLAST: reverse query coordinates" {
  let values = align_tabular_test_values()
  values[6] = "9"
  values[7] = "3"
  let document = @src.align_tabular_parse(
    align_tabular_test_block(
      "BLASTN 2.14.1+",
      "q1 reverse query",
      "nucleotide",
      align_tabular_test_fields(),
      [values],
      1,
    ),
  ) catch {
    _ => abort("reverse query should parse")
  }
  let alignment = document.queries[0].alignments[0]
  assert_true(alignment.query_is_reverse())
  assert_eq(alignment.query_coordinates, [9, 5, 5, 4, 2])
}

///|
test "Bio.Align.tabular BLAST: reverse target conversion normalizes target" {
  let values = align_tabular_test_values()
  values[8] = "16"
  values[9] = "10"
  let document = @src.align_tabular_parse(
    align_tabular_test_block(
      "BLASTN 2.14.1+",
      "q1 reverse target",
      "nucleotide",
      align_tabular_test_fields(),
      [values],
      1,
    ),
  ) catch {
    _ => abort("reverse target should parse")
  }
  let alignment = document.queries[0].alignments[0]
  assert_true(alignment.target_is_reverse())
  assert_eq(alignment.target_coordinates, [16, 12, 11, 11, 9])
  let converted = alignment.coordinate_alignment() catch {
    _ => abort("reverse target should normalize")
  }
  assert_eq(converted.target_coordinates, [9, 11, 11, 12, 16])
  assert_eq(converted.query_coordinates, [9, 7, 6, 6, 2])
}

///|
test "Bio.Align.tabular FASTA 8CC: command line and query size" {
  let text = "# fasta36 -q -m 8CC query.aa database.aa\n" +
    "# FASTA 36.3.8h May, 2020\n" +
    "# Query: qf protein query - 218 aa\n" +
    "# Database: database.aa\n" +
    "# Fields: query id, subject id, % identity, alignment length, mismatches, gap opens, q. start, q. end, s. start, s. end, evalue, bit score, aln_code\n" +
    "# 1 hits found\n" +
    "qf\tsf\t50.0\t9\t4\t2\t4\t10\t6\t13\t1e-4\t30.0\t3M2I1D3M\n" +
    "# FASTA processed 1 queries\n"
  let document = @src.align_tabular_parse(text) catch {
    _ => abort("valid FASTA 8CC should parse")
  }
  let query = document.queries[0]
  assert_eq(
    query.metadata.command_line,
    "fasta36 -q -m 8CC query.aa database.aa",
  )
  assert_eq(query.metadata.program, "FASTA")
  assert_eq(query.query_length, Some(218))
  assert_eq(query.query_unit, "aa")
  assert_eq(document.processed_queries, Some(1))
}

///|
test "Bio.Align.tabular FASTA 8CC: CIGAR path" {
  let text = "# fasta36 -m 8CC q db\n" +
    "# FASTA 36.3.8h\n" +
    "# Query: qf protein - 218 aa\n" +
    "# Database: db\n" +
    "# Fields: query id, subject id, alignment length, q. start, q. end, s. start, s. end, aln_code\n" +
    "# 1 hits found\n" +
    "qf\tsf\t9\t4\t10\t6\t13\t3M2I1D3M\n"
  let alignment = (@src.align_tabular_parse(text) catch {
      _ => abort("valid FASTA CIGAR should parse")
    }).queries[0].alignments[0]
  assert_true(alignment.trace is @src.AlignTabularTrace::Cigar)
  assert_eq(alignment.target_coordinates, [5, 8, 10, 10, 13])
  assert_eq(alignment.query_coordinates, [3, 6, 6, 7, 10])
  assert_eq(alignment.alignment_columns(), 9)
}

///|
test "Bio.Align.tabular BLAST: multiple queries include zero-hit blocks" {
  let first = align_tabular_test_block(
    "BLASTP 2.14.1+",
    "empty no matches",
    "proteins",
    [],
    [],
    0,
  )
  let second = align_tabular_test_block(
    "BLASTP 2.14.1+",
    "q2 one match",
    "proteins",
    ["query id", "subject id", "evalue"],
    [["q2", "s2", "0.01"]],
    1,
  )
  let document = @src.align_tabular_parse(
    first + second + "# BLAST processed 2 queries\n",
  ) catch {
    _ => abort("multi-query BLAST should parse")
  }
  assert_eq(document.queries.length(), 2)
  assert_eq(document.queries[0].alignments.length(), 0)
  assert_eq(document.queries[1].alignments.length(), 1)
  assert_eq(document.alignments().length(), 1)
}

///|
test "Bio.Align.tabular BLAST: processed summary is optional" {
  let document = align_tabular_test_document()
  assert_eq(document.processed_queries, None)
}

///|
test "Bio.Align.tabular BLAST: query and subject accession fallback" {
  let document = @src.align_tabular_parse(
    align_tabular_test_block(
      "BLASTP 2.14.1+",
      "QACC.1 accession query",
      "proteins",
      ["query acc.ver", "subject acc.ver", "evalue"],
      [["QACC.1", "SACC.2", "0.2"]],
      1,
    ),
  ) catch {
    _ => abort("accession-only identifiers should parse")
  }
  let alignment = document.queries[0].alignments[0]
  assert_eq(alignment.query_id, "QACC.1")
  assert_eq(alignment.target_id, "SACC.2")
}

///|
test "Bio.Align.tabular BLAST: rows without traceback retain annotations" {
  let document = @src.align_tabular_parse(
    align_tabular_test_block(
      "BLASTP 2.14.1+",
      "q1 query",
      "proteins",
      ["query id", "subject id", "alignment length", "evalue"],
      [["q1", "s1", "12", "0.4"]],
      1,
    ),
  ) catch {
    _ => abort("hit-table row should parse without traceback")
  }
  let alignment = document.queries[0].alignments[0]
  assert_true(alignment.trace is @src.AlignTabularTrace::NoTrace)
  assert_eq(alignment.operations.length(), 0)
  assert_true(!alignment.path_is_absolute)
}

///|
test "Bio.Align.tabular BLAST: complete upstream field vocabulary" {
  let fields = [
    "query id", "subject id", "% identity", "alignment length", "mismatches", "gap opens",
    "q. start", "q. end", "s. start", "s. end", "evalue", "bit score", "BTOP", "query gi",
    "query acc.", "query acc.ver", "query length", "subject ids", "subject gi", "subject gis",
    "subject acc.", "subject acc.ver", "subject accs.", "subject length", "query seq",
    "subject seq", "score", "identical", "positives", "gaps", "% positives", "% hsp coverage",
    "query/sbjct frames", "query frame", "sbjct frame", "subject tax ids", "subject sci names",
    "subject com names", "subject blast names", "subject super kingdoms", "subject title",
    "subject titles", "subject strand", "% subject coverage",
  ]
  let values = [
    "q1", "s1", "100", "3", "0", "0", "1", "3", "2", "4", "1e-20", "60.5", "3", "11",
    "QACC", "q1", "10", "s1;s2", "22", "22;23", "SACC", "SACC.1", "SACC;ALT", "20",
    "AAA", "AAA", "120", "3", "3", "0", "100", "30", "1/1", "1", "1", "9606", "Homo sapiens",
    "human", "primates", "Eukaryota", "subject title with spaces", "title one<>title two",
    "plus", "75",
  ]
  let alignment = (@src.align_tabular_parse(
      align_tabular_test_block(
        "BLASTP 2.14.1+",
        "q1 full fields",
        "database",
        fields,
        [values],
        1,
      ),
    ) catch {
      _ => abort("complete BLAST field vocabulary should parse")
    }).queries[0].alignments[0]
  assert_eq(alignment.fields.length(), 44)
  assert_eq(alignment.field("subject title"), Some("subject title with spaces"))
  assert_eq(alignment.field("subject tax ids"), Some("9606"))
}

///|
test "Bio.Align.tabular BLASTX: query nucleotide axis is tripled" {
  let alignment = (@src.align_tabular_parse(
      align_tabular_test_block(
        "BLASTX 2.14.1+",
        "qx translated query",
        "proteins",
        [
          "query id", "subject id", "query length", "subject length", "q. start",
          "q. end", "s. start", "s. end", "alignment length", "BTOP",
        ],
        [["qx", "sp", "30", "20", "1", "12", "5", "8", "4", "4"]],
        1,
      ),
    ) catch {
      _ => abort("BLASTX coordinates should parse")
    }).queries[0].alignments[0]
  assert_eq(alignment.query_units_per_residue, 3)
  assert_eq(alignment.target_units_per_residue, 1)
  assert_eq(alignment.query_coordinates, [0, 12])
  assert_eq(alignment.target_coordinates, [4, 8])
  assert_eq(alignment.alignment_columns(), 4)
}

///|
test "Bio.Align.tabular TBLASTN: target nucleotide axis is tripled" {
  let alignment = (@src.align_tabular_parse(
      align_tabular_test_block(
        "TBLASTN 2.14.1+",
        "qp protein query",
        "nucleotides",
        [
          "query id", "subject id", "query length", "subject length", "q. start",
          "q. end", "s. start", "s. end", "alignment length", "BTOP",
        ],
        [["qp", "sn", "10", "60", "2", "5", "7", "18", "4", "4"]],
        1,
      ),
    ) catch {
      _ => abort("TBLASTN coordinates should parse")
    }).queries[0].alignments[0]
  assert_eq(alignment.query_units_per_residue, 1)
  assert_eq(alignment.target_units_per_residue, 3)
  assert_eq(alignment.query_coordinates, [1, 5])
  assert_eq(alignment.target_coordinates, [6, 18])
}

///|
test "Bio.Align.tabular TBLASTX: both reverse nucleotide axes are tripled" {
  let alignment = (@src.align_tabular_parse(
      align_tabular_test_block(
        "TBLASTX 2.14.1+",
        "qn translated query",
        "nucleotides",
        [
          "query id", "subject id", "query length", "subject length", "q. start",
          "q. end", "s. start", "s. end", "alignment length", "BTOP",
        ],
        [["qn", "sn", "60", "90", "30", "19", "45", "34", "4", "4"]],
        1,
      ),
    ) catch {
      _ => abort("TBLASTX reverse coordinates should parse")
    }).queries[0].alignments[0]
  assert_eq(alignment.query_coordinates, [30, 18])
  assert_eq(alignment.target_coordinates, [45, 33])
  assert_true(alignment.query_is_reverse())
  assert_true(alignment.target_is_reverse())
}

///|
test "Bio.Align.tabular FASTX: query translated axis is tripled" {
  let text = "# fastx36 -m 8CB q.nt db.aa\n" +
    "# FASTX 36.3.8h\n" +
    "# Query: qx translated - 30 nt\n" +
    "# Database: db.aa\n" +
    "# Fields: query id, subject id, alignment length, q. start, q. end, s. start, s. end, BTOP\n" +
    "# 1 hits found\n" +
    "qx\tsp\t4\t1\t12\t5\t8\t4\n"
  let alignment = (@src.align_tabular_parse(text) catch {
      _ => abort("FASTX coordinates should parse")
    }).queries[0].alignments[0]
  assert_eq(alignment.query_units_per_residue, 3)
  assert_eq(alignment.query_coordinates, [0, 12])
}

///|
test "Bio.Align.tabular TFASTX: target translated axis is tripled" {
  let text = "# tfastx36 -m 8CB q.aa db.nt\n" +
    "# TFASTX 36.3.8h\n" +
    "# Query: qp protein - 10 aa\n" +
    "# Database: db.nt\n" +
    "# Fields: query id, subject id, alignment length, q. start, q. end, s. start, s. end, BTOP\n" +
    "# 1 hits found\n" +
    "qp\tsn\t4\t2\t5\t7\t18\t4\n"
  let alignment = (@src.align_tabular_parse(text) catch {
      _ => abort("TFASTX coordinates should parse")
    }).queries[0].alignments[0]
  assert_eq(alignment.target_units_per_residue, 3)
  assert_eq(alignment.target_coordinates, [6, 18])
}

///|
test "Bio.Align.tabular translated relative paths preserve axis units" {
  let alignment = (@src.align_tabular_parse(
      align_tabular_test_block(
        "BLASTX 2.14.1+",
        "qx translated query",
        "proteins",
        ["query id", "subject id", "alignment length", "BTOP"],
        [["qx", "sp", "4", "4"]],
        1,
      ),
    ) catch {
      _ => abort("relative translated path should parse")
    }).queries[0].alignments[0]
  assert_eq(alignment.query_coordinates, [0, 12])
  assert_eq(alignment.target_coordinates, [0, 4])
  assert_eq(alignment.alignment_columns(), 4)
  assert_true(!alignment.path_is_absolute)
}

///|
test "Bio.Align.tabular query lookup and flattening" {
  let document = align_tabular_test_document()
  match document.find_query("q1") {
    Some(query) => assert_eq(query.alignments.length(), 1)
    None => abort("query should be found")
  }
  assert_true(document.find_query("missing") is None)
  assert_eq(document.alignments().length(), 1)
}

///|
test "Bio.Align.tabular best hit uses minimum E-value" {
  let document = @src.align_tabular_parse(
    align_tabular_test_block(
      "BLASTP 2.14.1+",
      "q1 query",
      "proteins",
      ["query id", "subject id", "evalue"],
      [["q1", "late", "0.1"], ["q1", "best", "1e-20"], ["q1", "middle", "1e-5"]],
      3,
    ),
  ) catch {
    _ => abort("three hits should parse")
  }
  let best = document.queries[0].best_by_evalue() catch {
    _ => abort("best E-value should be computable")
  }
  match best {
    Some(alignment) => assert_eq(alignment.target_id, "best")
    None => abort("best hit should exist")
  }
}

///|
test "Bio.Align.tabular E-value filtering is inclusive" {
  let document = @src.align_tabular_parse(
    align_tabular_test_block(
      "BLASTP 2.14.1+",
      "q1 query",
      "proteins",
      ["query id", "subject id", "evalue"],
      [["q1", "a", "0.01"], ["q1", "b", "0.1"]],
      2,
    ),
  ) catch {
    _ => abort("two hits should parse")
  }
  let filtered = document.queries[0].filter_evalue(0.01) catch {
    _ => abort("finite threshold should filter")
  }
  assert_eq(filtered.length(), 1)
  assert_eq(filtered[0].target_id, "a")
}

///|
test "Bio.Align.tabular empty document is rejected" {
  assert_true(align_tabular_test_raises(" \n\t"))
}

///|
test "Bio.Align.tabular missing hash header is rejected" {
  assert_true(align_tabular_test_raises("BLASTP 2.14.1+\n"))
}

///|
test "Bio.Align.tabular missing FASTA program header is rejected" {
  assert_true(align_tabular_test_raises("# fasta36 -m 8CB q db\n"))
}

///|
test "Bio.Align.tabular FASTA command followed by BLAST is rejected" {
  assert_true(
    align_tabular_test_raises("# fasta36 -m 8CB q db\n# BLASTP 2.14.1+\n"),
  )
}

///|
test "Bio.Align.tabular missing Query header is rejected" {
  let text = align_tabular_test_block(
    "BLASTP 2.14.1+",
    "q1 query",
    "db",
    [],
    [],
    0,
  ).replace_all(old="# Query: q1 query\n", new="")
  assert_true(align_tabular_test_raises(text))
}

///|
test "Bio.Align.tabular missing Database header is rejected" {
  let text = align_tabular_test_block(
    "BLASTP 2.14.1+",
    "q1 query",
    "db",
    [],
    [],
    0,
  ).replace_all(old="# Database: db\n", new="")
  assert_true(align_tabular_test_raises(text))
}

///|
test "Bio.Align.tabular missing hit count is rejected" {
  let text = align_tabular_test_block(
    "BLASTP 2.14.1+",
    "q1 query",
    "db",
    [],
    [],
    0,
  ).replace_all(old="# 0 hits found\n", new="")
  assert_true(align_tabular_test_raises(text))
}

///|
test "Bio.Align.tabular non-empty block requires Fields" {
  let text = align_tabular_test_block(
    "BLASTP 2.14.1+",
    "q1 query",
    "db",
    [],
    [["q1", "s1"]],
    1,
  )
  assert_true(align_tabular_test_raises(text))
}

///|
test "Bio.Align.tabular unknown field is rejected" {
  let text = align_tabular_test_block(
    "BLASTP 2.14.1+",
    "q1 query",
    "db",
    ["query id", "subject id", "unknown score"],
    [["q1", "s1", "1"]],
    1,
  )
  assert_true(align_tabular_test_raises(text))
}

///|
test "Bio.Align.tabular duplicate field is rejected" {
  let text = align_tabular_test_block(
    "BLASTP 2.14.1+",
    "q1 query",
    "db",
    ["query id", "subject id", "subject id"],
    [["q1", "s1", "s1"]],
    1,
  )
  assert_true(align_tabular_test_raises(text))
}

///|
test "Bio.Align.tabular empty field name is rejected" {
  let text = align_tabular_test_block(
    "BLASTP 2.14.1+",
    "q1 query",
    "db",
    ["query id", "", "subject id"],
    [["q1", "", "s1"]],
    1,
  )
  assert_true(align_tabular_test_raises(text))
}

///|
test "Bio.Align.tabular duplicate Query header is rejected" {
  let text = align_tabular_test_block(
    "BLASTP 2.14.1+",
    "q1 query",
    "db",
    [],
    [],
    0,
  ).replace_all(
    old="# Query: q1 query\n",
    new="# Query: q1 query\n# Query: q1 query\n",
  )
  assert_true(align_tabular_test_raises(text))
}

///|
test "Bio.Align.tabular duplicate Database header is rejected" {
  let text = align_tabular_test_block(
    "BLASTP 2.14.1+",
    "q1 query",
    "db",
    [],
    [],
    0,
  ).replace_all(old="# Database: db\n", new="# Database: db\n# Database: db\n")
  assert_true(align_tabular_test_raises(text))
}

///|
test "Bio.Align.tabular duplicate RID header is rejected" {
  let text = "# BLASTP 2.14.1+\n" +
    "# Query: q1 query\n" +
    "# Database: db\n" +
    "# RID: one\n" +
    "# RID: two\n" +
    "# 0 hits found\n"
  assert_true(align_tabular_test_raises(text))
}

///|
test "Bio.Align.tabular duplicate Fields header is rejected" {
  let text = "# BLASTP 2.14.1+\n" +
    "# Query: q1 query\n" +
    "# Database: db\n" +
    "# Fields: query id, subject id\n" +
    "# Fields: query id, subject id\n" +
    "# 0 hits found\n"
  assert_true(align_tabular_test_raises(text))
}

///|
test "Bio.Align.tabular negative hit count is rejected" {
  let text = align_tabular_test_block(
    "BLASTP 2.14.1+",
    "q1 query",
    "db",
    [],
    [],
    -1,
  )
  assert_true(align_tabular_test_raises(text))
}

///|
test "Bio.Align.tabular declared hit count must match rows" {
  let text = align_tabular_test_block(
    "BLASTP 2.14.1+",
    "q1 query",
    "db",
    ["query id", "subject id"],
    [["q1", "s1"]],
    2,
  )
  assert_true(align_tabular_test_raises(text))
}

///|
test "Bio.Align.tabular processed count must match query blocks" {
  let text = align_tabular_test_block(
      "BLASTP 2.14.1+",
      "q1 query",
      "db",
      [],
      [],
      0,
    ) +
    "# BLAST processed 2 queries\n"
  assert_true(align_tabular_test_raises(text))
}

///|
test "Bio.Align.tabular duplicate processed summary is rejected" {
  let text = align_tabular_test_block(
      "BLASTP 2.14.1+",
      "q1 query",
      "db",
      [],
      [],
      0,
    ) +
    "# BLAST processed 1 queries\n" +
    "# BLAST processed 1 queries\n"
  assert_true(align_tabular_test_raises(text))
}

///|
test "Bio.Align.tabular negative processed count is rejected" {
  let text = align_tabular_test_block(
      "BLASTP 2.14.1+",
      "q1 query",
      "db",
      [],
      [],
      0,
    ) +
    "# BLAST processed -1 queries\n"
  assert_true(align_tabular_test_raises(text))
}

///|
test "Bio.Align.tabular row column count is validated" {
  let text = align_tabular_test_block(
    "BLASTP 2.14.1+",
    "q1 query",
    "db",
    ["query id", "subject id", "evalue"],
    [["q1", "s1"]],
    1,
  )
  assert_true(align_tabular_test_raises(text))
}

///|
test "Bio.Align.tabular invalid integer field is rejected" {
  let values = align_tabular_test_values()
  values[3] = "8oops"
  assert_true(
    align_tabular_test_raises(
      align_tabular_test_block(
        "BLASTP 2.14.1+",
        "q1 query",
        "db",
        align_tabular_test_fields(),
        [values],
        1,
      ),
    ),
  )
}

///|
test "Bio.Align.tabular invalid floating field is rejected" {
  let values = align_tabular_test_values()
  values[10] = "1e-5oops"
  assert_true(
    align_tabular_test_raises(
      align_tabular_test_block(
        "BLASTP 2.14.1+",
        "q1 query",
        "db",
        align_tabular_test_fields(),
        [values],
        1,
      ),
    ),
  )
}

///|
test "Bio.Align.tabular NaN is rejected" {
  let values = align_tabular_test_values()
  values[10] = "NaN"
  assert_true(
    align_tabular_test_raises(
      align_tabular_test_block(
        "BLASTP 2.14.1+",
        "q1 query",
        "db",
        align_tabular_test_fields(),
        [values],
        1,
      ),
    ),
  )
}

///|
test "Bio.Align.tabular percentage above one hundred is rejected" {
  let values = align_tabular_test_values()
  values[2] = "100.1"
  assert_true(
    align_tabular_test_raises(
      align_tabular_test_block(
        "BLASTP 2.14.1+",
        "q1 query",
        "db",
        align_tabular_test_fields(),
        [values],
        1,
      ),
    ),
  )
}

///|
test "Bio.Align.tabular negative E-value is rejected" {
  let values = align_tabular_test_values()
  values[10] = "-0.1"
  assert_true(
    align_tabular_test_raises(
      align_tabular_test_block(
        "BLASTP 2.14.1+",
        "q1 query",
        "db",
        align_tabular_test_fields(),
        [values],
        1,
      ),
    ),
  )
}

///|
test "Bio.Align.tabular zero coordinate is rejected" {
  let values = align_tabular_test_values()
  values[6] = "0"
  assert_true(
    align_tabular_test_raises(
      align_tabular_test_block(
        "BLASTP 2.14.1+",
        "q1 query",
        "db",
        align_tabular_test_fields(),
        [values],
        1,
      ),
    ),
  )
}

///|
test "Bio.Align.tabular coordinate beyond sequence length is rejected" {
  let values = align_tabular_test_values()
  values[7] = "21"
  assert_true(
    align_tabular_test_raises(
      align_tabular_test_block(
        "BLASTP 2.14.1+",
        "q1 query",
        "db",
        align_tabular_test_fields(),
        [values],
        1,
      ),
    ),
  )
}

///|
test "Bio.Align.tabular paired coordinates are required" {
  let text = align_tabular_test_block(
    "BLASTP 2.14.1+",
    "q1 query",
    "db",
    ["query id", "subject id", "q. start"],
    [["q1", "s1", "1"]],
    1,
  )
  assert_true(align_tabular_test_raises(text))
}

///|
test "Bio.Align.tabular traceback requires all four coordinates" {
  let text = align_tabular_test_block(
    "BLASTP 2.14.1+",
    "q1 query",
    "db",
    ["query id", "subject id", "q. start", "q. end", "alignment length", "BTOP"],
    [["q1", "s1", "1", "4", "4", "4"]],
    1,
  )
  assert_true(align_tabular_test_raises(text))
}

///|
test "Bio.Align.tabular traceback span must match reported coordinates" {
  let values = align_tabular_test_values()
  values[7] = "10"
  assert_true(
    align_tabular_test_raises(
      align_tabular_test_block(
        "BLASTP 2.14.1+",
        "q1 query",
        "db",
        align_tabular_test_fields(),
        [values],
        1,
      ),
    ),
  )
}

///|
test "Bio.Align.tabular traceback columns must match alignment length" {
  let values = align_tabular_test_values()
  values[3] = "9"
  assert_true(
    align_tabular_test_raises(
      align_tabular_test_block(
        "BLASTP 2.14.1+",
        "q1 query",
        "db",
        align_tabular_test_fields(),
        [values],
        1,
      ),
    ),
  )
}

///|
test "Bio.Align.tabular aligned sequence lengths must agree" {
  let values = align_tabular_test_values()
  values[16] = "AACTG-G"
  assert_true(
    align_tabular_test_raises(
      align_tabular_test_block(
        "BLASTP 2.14.1+",
        "q1 query",
        "db",
        align_tabular_test_fields(),
        [values],
        1,
      ),
    ),
  )
}

///|
test "Bio.Align.tabular aligned sequence columns must match traceback" {
  let values = align_tabular_test_values()
  values[15] = "AAAT-TG"
  values[16] = "AACTG-G"
  assert_true(
    align_tabular_test_raises(
      align_tabular_test_block(
        "BLASTP 2.14.1+",
        "q1 query",
        "db",
        align_tabular_test_fields(),
        [values],
        1,
      ),
    ),
  )
}

///|
test "Bio.Align.tabular sequence consumption must match traceback" {
  let values = align_tabular_test_values()
  values[15] = "AAAATTGG"
  assert_true(
    align_tabular_test_raises(
      align_tabular_test_block(
        "BLASTP 2.14.1+",
        "q1 query",
        "db",
        align_tabular_test_fields(),
        [values],
        1,
      ),
    ),
  )
}

///|
test "Bio.Align.tabular sequence whitespace is rejected" {
  let fields = ["query id", "subject id", "query seq", "subject seq"]
  let text = align_tabular_test_block(
    "BLASTP 2.14.1+",
    "q1 query",
    "db",
    fields,
    [["q1", "s1", "AA AA", "AAAAA"]],
    1,
  )
  assert_true(align_tabular_test_raises(text))
}

///|
test "Bio.Align.tabular BTOP and CIGAR cannot coexist" {
  let fields = ["query id", "subject id", "BTOP", "aln_code"]
  let text = align_tabular_test_block(
    "BLASTP 2.14.1+",
    "q1 query",
    "db",
    fields,
    [["q1", "s1", "4", "4M"]],
    1,
  )
  assert_true(align_tabular_test_raises(text))
}

///|
test "Bio.Align.tabular query identifier must match header" {
  let text = align_tabular_test_block(
    "BLASTP 2.14.1+",
    "q1 query",
    "db",
    ["query id", "subject id"],
    [["other", "s1"]],
    1,
  )
  assert_true(align_tabular_test_raises(text))
}

///|
test "Bio.Align.tabular subject identifier is required" {
  let text = align_tabular_test_block(
    "BLASTP 2.14.1+",
    "q1 query",
    "db",
    ["query id", "evalue"],
    [["q1", "0.1"]],
    1,
  )
  assert_true(align_tabular_test_raises(text))
}

///|
test "Bio.Align.tabular query length must be stable in block" {
  let text = align_tabular_test_block(
    "BLASTP 2.14.1+",
    "q1 query",
    "db",
    ["query id", "subject id", "query length"],
    [["q1", "s1", "10"], ["q1", "s2", "11"]],
    2,
  )
  assert_true(align_tabular_test_raises(text))
}

///|
test "Bio.Align.tabular FASTA query size must agree with row" {
  let text = "# fasta36 -m 8CB q db\n" +
    "# FASTA 36.3.8h\n" +
    "# Query: q1 query - 10 aa\n" +
    "# Database: db\n" +
    "# Fields: query id, subject id, query length\n" +
    "# 1 hits found\n" +
    "q1\ts1\t11\n"
  assert_true(align_tabular_test_raises(text))
}

///|
test "Bio.Align.tabular FASTA query size must be positive" {
  let text = "# fasta36 -m 8CB q db\n" +
    "# FASTA 36.3.8h\n" +
    "# Query: q1 query - 0 aa\n" +
    "# Database: db\n" +
    "# 0 hits found\n"
  assert_true(align_tabular_test_raises(text))
}

///|
test "Bio.Align.tabular BTOP rejects empty value" {
  assert_true(align_tabular_test_btop_raises(""))
}

///|
test "Bio.Align.tabular BTOP rejects zero match run" {
  assert_true(align_tabular_test_btop_raises("0"))
}

///|
test "Bio.Align.tabular BTOP rejects truncated pair" {
  assert_true(align_tabular_test_btop_raises("2A"))
}

///|
test "Bio.Align.tabular BTOP rejects invalid residue pair" {
  assert_true(align_tabular_test_btop_raises("A?"))
}

///|
test "Bio.Align.tabular BTOP rejects two gaps" {
  assert_true(align_tabular_test_btop_raises("--"))
}

///|
test "Bio.Align.tabular CIGAR rejects empty value" {
  assert_true(align_tabular_test_cigar_raises(""))
}

///|
test "Bio.Align.tabular CIGAR requires leading length" {
  assert_true(align_tabular_test_cigar_raises("M3"))
}

///|
test "Bio.Align.tabular CIGAR requires trailing operation" {
  assert_true(align_tabular_test_cigar_raises("3"))
}

///|
test "Bio.Align.tabular CIGAR rejects zero length" {
  assert_true(align_tabular_test_cigar_raises("0M"))
}

///|
test "Bio.Align.tabular CIGAR rejects unsupported operation" {
  assert_true(align_tabular_test_cigar_raises("3S"))
}

///|
test "Bio.Align.tabular relative path cannot convert to coordinate alignment" {
  let alignment = (@src.align_tabular_parse(
      align_tabular_test_block(
        "BLASTP 2.14.1+",
        "q1 query",
        "db",
        ["query id", "subject id", "alignment length", "BTOP"],
        [["q1", "s1", "4", "4"]],
        1,
      ),
    ) catch {
      _ => abort("relative path should parse")
    }).queries[0].alignments[0]
  let rejected = try {
    ignore(alignment.coordinate_alignment())
    false
  } catch {
    AlignTabularError(_) => true
  }
  assert_true(rejected)
}

///|
test "Bio.Align.tabular E-value filter rejects negative threshold" {
  let rejected = try {
    ignore(align_tabular_test_document().queries[0].filter_evalue(-0.1))
    false
  } catch {
    AlignTabularError(_) => true
  }
  assert_true(rejected)
}
