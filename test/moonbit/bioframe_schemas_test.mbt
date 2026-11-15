///|
/// Tests for bioframe_schemas.mbt — schema registry + tabular I/O.

// ---- bf_schemas_available -------------------------------------------------

test "bf_schemas_available: 18 schemas" {
  let all = @src.bf_schemas_available()
  assert_eq(all.length(), 18)
  assert_true(all.contains("bed3"))
  assert_true(all.contains("bed12"))
  assert_true(all.contains("narrowPeak"))
  assert_true(all.contains("pairs"))
  assert_true(all.contains("refSeq"))
}

// ---- bf_schema_cols -------------------------------------------------------

test "bf_schema_cols: bed3" {
  match @src.bf_schema_cols("bed3") {
    Some(cols) => {
      assert_eq(cols.length(), 3)
      assert_eq(cols[0], "chrom")
      assert_eq(cols[1], "start")
      assert_eq(cols[2], "end")
    }
    None => assert_true(false)
  }
}

test "bf_schema_cols: narrowPeak" {
  match @src.bf_schema_cols("narrowPeak") {
    Some(cols) => assert_eq(cols.length(), 10)
    None => assert_true(false)
  }
}

test "bf_schema_cols: unknown → None" {
  match @src.bf_schema_cols("notASchema") {
    Some(_) => assert_true(false)
    None => assert_true(true)
  }
}

// ---- bf_read_table: with schema ------------------------------------------

test "bf_read_table: no header + schema=bed3" {
  let text = "chr1\t0\t100\nchr1\t200\t300\n"
  let (header, rows) = @src.bf_read_table(text, schema=Some("bed3"))
  assert_eq(header.length(), 3)
  assert_eq(header[0], "chrom")
  assert_eq(rows.length(), 2)
}

test "bf_read_table: # comment lines skipped" {
  let text = "# this is a comment\nchr1\t0\t100\n# another comment\nchr1\t200\t300\n"
  let (_, rows) = @src.bf_read_table(text, schema=Some("bed3"))
  assert_eq(rows.length(), 2)
}

test "bf_read_table: track/browser UCSC lines skipped" {
  let text = "track name=myTrack\nchr1\t0\t100\nbrowser position chr1\nchr1\t200\t300\n"
  let (_, rows) = @src.bf_read_table(text, schema=Some("bed3"))
  assert_eq(rows.length(), 2)
}

test "bf_read_table: with explicit header (auto-detected)" {
  let text = "chrom\tstart\tend\nchr1\t0\t100\nchr1\t200\t300\n"
  let (header, rows) = @src.bf_read_table(text)
  assert_eq(header[0], "chrom")
  assert_eq(header[1], "start")
  assert_eq(header[2], "end")
  assert_eq(rows.length(), 2)
}

test "bf_read_table: all-numeric first row → default header" {
  let text = "1\t2\t3\n10\t20\t30\n"
  let (header, rows) = @src.bf_read_table(text)
  // first row looks numeric → default col0/col1/col2
  assert_eq(header.length(), 3)
  assert_true(header[0].has_prefix("col"))
  assert_eq(rows.length(), 2)
}

// ---- bf_row_get ----------------------------------------------------------

test "bf_row_get: retrieve by column name" {
  let text = "chrom\tstart\tend\nchr1\t0\t100\n"
  let (header, rows) = @src.bf_read_table(text)
  match rows[0].field(0) {
    Some(v) => assert_eq(v, "chr1")
    None => assert_true(false)
  }
}

// ---- bf_rows_to_intervals -------------------------------------------------

test "bf_rows_to_intervals: bed3 rows → BedInterval array" {
  let text = "chr1\t0\t100\nchr1\t200\t300\nchr2\t50\t150\n"
  let (header, rows) = @src.bf_read_table(text, schema=Some("bed3"))
  let ivs = @src.bf_rows_to_intervals(header, rows)
  assert_eq(ivs.length(), 3)
  assert_eq(ivs[0].chrom, "chr1")
  assert_eq(ivs[0].start, 0)
  assert_eq(ivs[0].end, 100)
  assert_eq(ivs[2].chrom, "chr2")
  assert_eq(ivs[2].start, 50)
  assert_eq(ivs[2].end, 150)
}

test "bf_rows_to_intervals: empty → empty" {
  let ivs = @src.bf_rows_to_intervals([], [])
  assert_eq(ivs.length(), 0)
}
