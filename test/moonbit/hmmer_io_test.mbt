///|
/// Tests for Bio.SearchIO.HmmerIO - HMMER3 domtblout and text output parsing.

///|
test "hmmer_io_parse_domtab_basic" {
  let sample = @src.hmmer3_domtab_sample()
  let results = @src.parse_hmmer3_domtab(sample)
  // 3 data rows: 2 targets for PF00001.28, 1 target for PF00002.10
  assert_eq(results.length(), 1)
  let qr = results[0]
  assert_eq(qr.id, "PF00001.28")
  assert_eq(qr.hits.length(), 2)
}

///|
test "hmmer_io_parse_domtab_multiple_queries" {
  let content = "# header\n" +
    "target1 acc1 queryA qacc1 1e-50 200.0 1.0 1 1 1e-50 1e-53 200.0 1.0 1 150 1 150 1 150 1.0 1 0 0 1 1 1 1 desc1\n" +
    "target2 acc2 queryB qacc2 1e-30 150.0 2.0 1 1 1e-30 1e-33 150.0 2.0 1 120 1 120 1 120 1.0 1 0 0 1 1 1 1 desc2\n"
  let results = @src.parse_hmmer3_domtab(content)
  assert_eq(results.length(), 2)
  assert_eq(results[0].id, "queryA")
  assert_eq(results[1].id, "queryB")
}

///|
test "hmmer_io_parse_domtab_multi_domain" {
  let content = "# header\n" +
    "target1 acc1 queryA qacc1 1e-50 200.0 1.0 1 2 1e-50 1e-53 100.0 1.0 1 80 1 80 1 80 1.0 1 0 0 1 1 1 1 desc1\n" +
    "target1 acc1 queryA qacc1 1e-50 200.0 1.0 2 2 1e-40 1e-43 100.0 1.0 81 150 81 150 81 150 1.0 1 0 0 1 1 1 1 desc1\n"
  let results = @src.parse_hmmer3_domtab(content)
  assert_eq(results.length(), 1)
  assert_eq(results[0].hits.length(), 1)
  // Two domains -> two HSPs.
  assert_eq(results[0].hits[0].hsps.length(), 2)
}

///|
test "hmmer_io_parse_domtab_empty" {
  let results = @src.parse_hmmer3_domtab("")
  assert_eq(results.length(), 0)
}

///|
test "hmmer_io_parse_domtab_comments_only" {
  let content = "# hmmsearch\n# HMMER 3.3.2\n# header line\n"
  let results = @src.parse_hmmer3_domtab(content)
  assert_eq(results.length(), 0)
}

///|
test "hmmer_io_parse_domtab_short_line_skipped" {
  let content = "# header\nshort line\n"
  let results = @src.parse_hmmer3_domtab(content)
  assert_eq(results.length(), 0)
}

///|
test "hmmer_io_parse_domtab_evalue_bitscore" {
  let sample = @src.hmmer3_domtab_sample()
  let results = @src.parse_hmmer3_domtab(sample)
  let qr = results[0]
  let hit0 = qr.hits[0]
  // First hit is sp|Q9Y2W8|A1BG_HUMAN with i-Evalue 1.5e-53.
  assert_eq(hit0.id, "sp|Q9Y2W8|A1BG_HUMAN")
  assert_eq(hit0.evalue < 1.0e-50, true)
  assert_eq(hit0.bitscore > 100.0, true)
}

///|
test "hmmer_io_parse_domtab_count_hsps" {
  let sample = @src.hmmer3_domtab_sample()
  let results = @src.parse_hmmer3_domtab(sample)
  let count = @src.hmmer_io_count_hsps(results[0])
  // 3 data rows -> 3 HSPs total.
  assert_eq(count, 3)
}

///|
test "hmmer_io_parse_text_basic" {
  let sample = @src.hmmer3_text_sample()
  let results = @src.parse_hmmer3_text(sample)
  assert_eq(results.length(), 1)
  let qr = results[0]
  assert_eq(qr.id, "PF00001.28")
  assert_eq(qr.seq_len, 150)
  assert_eq(qr.hits.length(), 2)
}

///|
test "hmmer_io_parse_text_empty" {
  let results = @src.parse_hmmer3_text("")
  assert_eq(results.length(), 0)
}

///|
test "hmmer_io_parse_text_query_length" {
  let content = "# header\n" +
    "Query:       PF00001  7tm_1\n" +
    "Query length: 245\n" +
    ">> target1  description1\n" +
    "  == domain 1   score: 200.3 bits;  conditional E-value: 1.5e-53\n"
  let results = @src.parse_hmmer3_text(content)
  assert_eq(results.length(), 1)
  assert_eq(results[0].seq_len, 245)
  assert_eq(results[0].hits.length(), 1)
  assert_eq(results[0].hits[0].hsps.length(), 1)
}

///|
test "hmmer_io_parse_text_multi_query" {
  let content = "# header\n" +
    "Query:       q1  desc1\n" +
    "Query length: 100\n" +
    ">> t1  d1\n" +
    "  == domain 1   score: 100.0 bits;  conditional E-value: 1e-20\n" +
    "Query:       q2  desc2\n" +
    "Query length: 200\n" +
    ">> t2  d2\n" +
    "  == domain 1   score: 50.0 bits;  conditional E-value: 1e-10\n"
  let results = @src.parse_hmmer3_text(content)
  assert_eq(results.length(), 2)
  assert_eq(results[0].id, "q1")
  assert_eq(results[0].seq_len, 100)
  assert_eq(results[1].id, "q2")
  assert_eq(results[1].seq_len, 200)
}

///|
test "hmmer_io_parse_text_multi_domain" {
  let content = "# header\n" +
    "Query:       q1  desc1\n" +
    "Query length: 100\n" +
    ">> t1  d1\n" +
    "  == domain 1   score: 100.0 bits;  conditional E-value: 1e-20\n" +
    "  == domain 2   score: 80.0 bits;  conditional E-value: 1e-15\n"
  let results = @src.parse_hmmer3_text(content)
  assert_eq(results.length(), 1)
  assert_eq(results[0].hits.length(), 1)
  assert_eq(results[0].hits[0].hsps.length(), 2)
}

///|
test "hmmer_io_parse_text_score_extraction" {
  let content = "Query: q1\n" +
    ">> t1  d1\n" +
    "  == domain 1   score: 200.3 bits;  conditional E-value: 1.5e-53\n"
  let results = @src.parse_hmmer3_text(content)
  let hsp = results[0].hits[0].hsps[0]
  assert_eq(hsp.bitscore > 200.0, true)
  assert_eq(hsp.evalue < 1.0e-50, true)
}

///|
test "hmmer_io_parse_text_multi_target" {
  let content = "Query: q1\n" +
    ">> t1  d1\n" +
    "  == domain 1   score: 100.0 bits;  conditional E-value: 1e-20\n" +
    ">> t2  d2\n" +
    "  == domain 1   score: 80.0 bits;  conditional E-value: 1e-15\n"
  let results = @src.parse_hmmer3_text(content)
  assert_eq(results[0].hits.length(), 2)
  assert_eq(results[0].hits[0].id, "t1")
  assert_eq(results[0].hits[1].id, "t2")
}

///|
test "hmmer_io_parse_text_description" {
  let content = "Query: q1  some query description\n" +
    ">> t1  Alpha-1-B glycoprotein\n" +
    "  == domain 1   score: 100.0 bits;  conditional E-value: 1e-20\n"
  let results = @src.parse_hmmer3_text(content)
  assert_eq(results[0].description, "some query description")
  assert_eq(results[0].hits[0].description, "Alpha-1-B glycoprotein")
}

///|
test "hmmer_io_sample_domtab_not_empty" {
  let sample = @src.hmmer3_domtab_sample()
  assert_eq(sample.length() > 0, true)
}

///|
test "hmmer_io_sample_text_not_empty" {
  let sample = @src.hmmer3_text_sample()
  assert_eq(sample.length() > 0, true)
}

///|
test "hmmer_io_parse_text_skip_header_lines" {
  let content = "# hmmsearch :: search profile(s) against a sequence database\n" +
    "# HMMER 3.3.2 (Nov 2020)\n" +
    "#\n" +
    "Query: q1\n" +
    "Scores for complete sequences (score includes all domains):\n" +
    "   E-value  score  bias    E-value  score  bias    exp  N  Domain  Description\n" +
    "   ------- ------ -----    ------- ------ -----    ---- ---  ------  -----------\n" +
    ">> t1  d1\n" +
    "  == domain 1   score: 100.0 bits;  conditional E-value: 1e-20\n"
  let results = @src.parse_hmmer3_text(content)
  assert_eq(results.length(), 1)
  assert_eq(results[0].hits.length(), 1)
}
