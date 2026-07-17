///|
test "search_io_parse_hmmer3_tab_basic" {
  let hmmer3_tab = "# hmmscan 3.3.2 (March 2023)\n" +
    "# query name           target name        accession   E-value  score  bias\n" +
    "PF00001.28            sp|Q9Y2W8|A1BG_HUMAN  Q9Y2W8     1e-50    200    5\n" +
    "PF00001.28            sp|P12345|ABC_HUMAN   P12345     1e-30    150    2\n"
  let results = @src.parse_hmmer3_tab(hmmer3_tab)
  assert_eq(results.length(), 1)
  let qr = results[0]
  assert_eq(qr.id, "PF00001.28")
  assert_eq(qr.hits.length(), 2)
}

///|
test "search_io_parse_hmmer3_tab_empty" {
  let results = @src.parse_hmmer3_tab("")
  assert_eq(results.length(), 0)
}

///|
test "search_io_parse_blat_psl_basic" {
  let psl = "100\t0\t0\t0\t0\t0\t0\t0\t+\tquery1\t1000\t0\t100\ttarget1\t2000\t0\t100\t1\t100,\t0,\t0,\n" +
    "80\t0\t0\t0\t0\t0\t0\t0\t-\tquery1\t1000\t100\t180\ttarget2\t1500\t50\t130\t1\t80,\t100,\t50,\n"
  let results = @src.parse_blat_psl(psl)
  assert_eq(results.length(), 1)
  let qr = results[0]
  assert_eq(qr.id, "query1")
  assert_eq(qr.seq_len, 1000)
  assert_eq(qr.hits.length(), 2)
}

///|
test "search_io_parse_blat_psl_empty" {
  let results = @src.parse_blat_psl("")
  assert_eq(results.length(), 0)
}

///|
test "search_io_count_hsps" {
  let hmmer3_tab = "# hmmscan 3.3.2\n" +
    "PF00001.28            sp|Q9Y2W8|A1BG_HUMAN  Q9Y2W8     1e-50    200    5\n" +
    "PF00001.28            sp|P12345|ABC_HUMAN   P12345     1e-30    150    2\n"
  let results = @src.parse_hmmer3_tab(hmmer3_tab)
  assert_eq(results.length(), 1)
  let count = @src.count_hsps(results[0])
  assert_eq(count, 2)
}

///|
test "search_io_top_hits" {
  let hmmer3_tab = "# hmmscan 3.3.2\n" +
    "PF00001.28            sp|Q9Y2W8|A1BG_HUMAN  Q9Y2W8     1e-50    200    5\n" +
    "PF00001.28            sp|P12345|ABC_HUMAN   P12345     1e-30    150    2\n" +
    "PF00001.28            sp|X99999|XYZ_HUMAN   X99999     0.5      50     1\n"
  let results = @src.parse_hmmer3_tab(hmmer3_tab)
  assert_eq(results.length(), 1)
  let tops = @src.top_hits(results[0], 2, 0.01)
  assert_eq(tops.length(), 2)
}

///|
test "search_io_top_hits_threshold" {
  let hmmer3_tab = "# hmmscan 3.3.2\n" +
    "PF00001.28            sp|Q9Y2W8|A1BG_HUMAN  Q9Y2W8     1e-50    200    5\n" +
    "PF00001.28            sp|P12345|ABC_HUMAN   P12345     0.02     150    2\n" +
    "PF00001.28            sp|X99999|XYZ_HUMAN   X99999     0.5      50     1\n"
  let results = @src.parse_hmmer3_tab(hmmer3_tab)
  assert_eq(results.length(), 1)
  let tops = @src.top_hits(results[0], 10, 0.01)
  assert_eq(tops.length(), 1)
}
