///|
/// Tests for biotite_sequence_search.mbt.
test "biotite_blast_local_command" {
  let query = @src.biotite_blast_query("blastn", "nt")
  let command = @src.biotite_blast_local_command(query, "query.fa", "out.txt")
  assert_true(command.contains("blastn"))
  assert_true(command.contains("-query query.fa"))
  assert_true(command.contains("-db nt"))
  assert_true(command.contains("-outfmt 6"))
}

///|
test "biotite_blast_remote_url" {
  let query = @src.biotite_blast_query("blastp", "nr")
  let url = @src.biotite_blast_remote_url(query, "ACDE")
  assert_true(url.contains("PROGRAM=blastp"))
  assert_true(url.contains("DATABASE=nr"))
  assert_true(url.contains("QUERY=ACDE"))
}

///|
test "biotite_blast_parse_tabular" {
  let text = "# BLASTP\nq1\ts1\t99.5\t100\t0\t0\t1\t100\t1\t100\t0.001\t150\n"
  let rows = @src.biotite_blast_parse_tabular(text)
  assert_eq(rows.length(), 1)
  assert_eq(rows[0].0, "q1")
  assert_eq(rows[0].1, "s1")
  assert_true(@moonbit.assert_num_close(rows[0].2, 99.5, 1.0e-12))
  assert_true(@moonbit.assert_num_close(rows[0].3, 0.001, 1.0e-12))
}

///|
test "biotite_blast_parse_tabular_extended" {
  let text = "q1\ts1\t95.0\t120\t0\t0\t1\t120\t1\t120\t1.5\t200\n"
  let rows = @src.biotite_blast_parse_tabular_extended(text)
  assert_eq(rows.length(), 1)
  assert_eq(rows[0].3, 120)
  assert_true(@moonbit.assert_num_close(rows[0].4, 1.5, 1.0e-12))
}

///|
test "biotite_blast_parse_skips_bad_lines" {
  assert_eq(@src.biotite_blast_parse_tabular("").length(), 0)
  assert_eq(@src.biotite_blast_parse_tabular("too\tshort\n").length(), 0)
}
