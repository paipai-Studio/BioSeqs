///|
/// Tests for Bio.SearchIO.FastaIO - FASTA -m8/-m9 search output parsing.

///|
test "fasta_sio_parse_m8_basic" {
  let sample = @src.fasta_m8_sample()
  let results = @src.parse_fasta_m8(sample)
  // 3 rows: query1 has 2 hits (subject1, subject2), query2 has 1 hit.
  assert_eq(results.length(), 2)
  assert_eq(results[0].id, "query1")
  assert_eq(results[0].hits.length(), 2)
  assert_eq(results[1].id, "query2")
  assert_eq(results[1].hits.length(), 1)
}

///|
test "fasta_sio_parse_m8_empty" {
  let results = @src.parse_fasta_m8("")
  assert_eq(results.length(), 0)
}

///|
test "fasta_sio_parse_m8_comment_lines_skipped" {
  let content = "# comment line\n# another comment\nquery1\tsubject1\t95.0\t100\t5\t0\t1\t100\t1\t100\t1e-50\t200.3\n"
  let results = @src.parse_fasta_m8(content)
  assert_eq(results.length(), 1)
  assert_eq(results[0].hits.length(), 1)
}

///|
test "fasta_sio_parse_m8_multiple_hsp_per_target" {
  let content = "query1\tsubject1\t95.0\t100\t5\t0\t1\t100\t1\t100\t1e-50\t200.3\n" +
    "query1\tsubject1\t80.0\t80\t16\t2\t1\t80\t10\t90\t1e-30\t150.1\n"
  let results = @src.parse_fasta_m8(content)
  assert_eq(results.length(), 1)
  assert_eq(results[0].hits.length(), 1)
  // Two rows for the same target -> two HSPs.
  assert_eq(results[0].hits[0].hsps.length(), 2)
}

///|
test "fasta_sio_parse_m8_fields" {
  let content = "q1\ts1\t95.0\t100\t5\t0\t1\t100\t1\t100\t1e-50\t200.3\n"
  let results = @src.parse_fasta_m8(content)
  let hsp = results[0].hits[0].hsps[0]
  let frag = hsp.fragments[0]
  assert_eq(frag.query_start, 1)
  assert_eq(frag.query_end, 100)
  assert_eq(frag.hit_start, 1)
  assert_eq(frag.hit_end, 100)
  assert_eq(frag.aln_span, 100)
  assert_eq(frag.gap_num, 0)
  assert_eq(frag.bitscore > 200.0, true)
  assert_eq(frag.evalue < 1.0e-40, true)
}

///|
test "fasta_sio_parse_m8_whitespace_separated" {
  // Lines without tabs should fall back to whitespace splitting.
  let content = "q1 s1 95.0 100 5 0 1 100 1 100 1e-50 200.3\n"
  let results = @src.parse_fasta_m8(content)
  assert_eq(results.length(), 1)
  assert_eq(results[0].hits.length(), 1)
}

///|
test "fasta_sio_parse_m8_short_row_padded" {
  // Only 6 columns - should be padded with defaults.
  let content = "q1\ts1\t95.0\t100\t5\t0\n"
  let results = @src.parse_fasta_m8(content)
  assert_eq(results.length(), 1)
  let frag = results[0].hits[0].hsps[0].fragments[0]
  assert_eq(frag.query_start, 0)
  assert_eq(frag.query_end, 0)
}

///|
test "fasta_sio_parse_m9_basic" {
  let sample = @src.fasta_m9_sample()
  let result = @src.parse_fasta_m9(sample)
  assert_eq(result.queries.length(), 2)
  assert_eq(result.queries[0].id, "query1")
  assert_eq(result.queries[0].hits.length(), 2)
  assert_eq(result.queries[1].id, "query2")
  assert_eq(result.queries[1].hits.length(), 1)
}

///|
test "fasta_sio_parse_m9_meta_program" {
  let sample = @src.fasta_m9_sample()
  let result = @src.parse_fasta_m9(sample)
  assert_eq(result.meta.program, "FASTA")
  assert_eq(result.meta.version, "36.3.8h")
}

///|
test "fasta_sio_parse_m9_meta_database" {
  let sample = @src.fasta_m9_sample()
  let result = @src.parse_fasta_m9(sample)
  assert_eq(result.meta.database, "uniprot_sprot.fasta")
}

///|
test "fasta_sio_parse_m9_comments_collected" {
  let sample = @src.fasta_m9_sample()
  let result = @src.parse_fasta_m9(sample)
  // The sample has 4 comment lines.
  assert_eq(result.meta.comments.length() >= 4, true)
}

///|
test "fasta_sio_parse_m9_empty" {
  let result = @src.parse_fasta_m9("")
  assert_eq(result.queries.length(), 0)
  assert_eq(result.meta.comments.length(), 0)
}

///|
test "fasta_sio_parse_m9_comments_only" {
  let content = "# FASTA 36.3.8h\n# Database: test.fasta\n"
  let result = @src.parse_fasta_m9(content)
  assert_eq(result.queries.length(), 0)
  assert_eq(result.meta.program, "FASTA")
  assert_eq(result.meta.database, "test.fasta")
}

///|
test "fasta_sio_parse_m8_multiple_queries" {
  let content = "q1\ts1\t95.0\t100\t5\t0\t1\t100\t1\t100\t1e-50\t200.3\n" +
    "q2\ts2\t80.0\t80\t16\t2\t1\t80\t10\t90\t1e-30\t150.1\n" +
    "q3\ts1\t70.0\t60\t18\t3\t1\t60\t5\t65\t1e-20\t100.5\n"
  let results = @src.parse_fasta_m8(content)
  assert_eq(results.length(), 3)
  assert_eq(results[0].id, "q1")
  assert_eq(results[1].id, "q2")
  assert_eq(results[2].id, "q3")
}

///|
test "fasta_sio_sample_m8_not_empty" {
  let sample = @src.fasta_m8_sample()
  assert_eq(sample.length() > 0, true)
}

///|
test "fasta_sio_sample_m9_not_empty" {
  let sample = @src.fasta_m9_sample()
  assert_eq(sample.length() > 0, true)
}

///|
test "fasta_sio_parse_m8_evalue_bitscore" {
  let content = "q1\ts1\t95.0\t100\t5\t0\t1\t100\t1\t100\t1e-50\t200.3\n"
  let results = @src.parse_fasta_m8(content)
  let hit = results[0].hits[0]
  assert_eq(hit.bitscore > 200.0, true)
  assert_eq(hit.evalue < 1.0e-40, true)
}

///|
test "fasta_sio_parse_m8_reverse_strand" {
  // query_end < query_start indicates reverse strand alignment.
  let content = "q1\ts1\t95.0\t100\t5\t0\t100\t1\t1\t100\t1e-50\t200.3\n"
  let results = @src.parse_fasta_m8(content)
  let frag = results[0].hits[0].hsps[0].fragments[0]
  assert_eq(frag.query_strand == @src.strand_minus(), true)
}

///|
test "fasta_sio_parse_m9_data_after_comments" {
  let content = "# FASTA 36.3\n# Database: db\n" +
    "q1\ts1\t95.0\t100\t5\t0\t1\t100\t1\t100\t1e-50\t200.3\n"
  let result = @src.parse_fasta_m9(content)
  assert_eq(result.queries.length(), 1)
  assert_eq(result.queries[0].hits.length(), 1)
}
