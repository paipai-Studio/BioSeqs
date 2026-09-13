///|
/// Tests for the Biopython Bio.Affy.CelFile-compatible parser.

///|
test "cel parses header and grid dimensions" {
  let record = @src.cel_read(@src.cel_sample_text()) catch {
    _ => abort("failed to parse sample CEL")
  }
  assert_eq(record.version, "3.4")
  assert_eq(record.n_rows, 2)
  assert_eq(record.n_cols, 3)
  assert_eq(record.n_cells, 6)
  assert_eq(record.cells_read, 6)
}

///|
test "cel parses intensity matrix" {
  let record = @src.cel_read(@src.cel_sample_text()) catch {
    _ => abort("failed to parse sample CEL")
  }
  // Row 0: x=0 -> 91.5, x=1 -> 88.0, x=2 -> 102.3
  assert_true((record.intensity_at(0, 0) - 91.5).abs() < 1.0e-9)
  assert_true((record.intensity_at(1, 0) - 88.0).abs() < 1.0e-9)
  assert_true((record.intensity_at(2, 0) - 102.3).abs() < 1.0e-9)
  // Row 1: x=0 -> 141.2, x=1 -> 130.0, x=2 -> 120.9
  assert_true((record.intensity_at(0, 1) - 141.2).abs() < 1.0e-9)
  assert_true((record.intensity_at(1, 1) - 130.0).abs() < 1.0e-9)
  assert_true((record.intensity_at(2, 1) - 120.9).abs() < 1.0e-9)
}

///|
test "cel parses stdev and npix" {
  let record = @src.cel_read(@src.cel_sample_text()) catch {
    _ => abort("failed to parse sample CEL")
  }
  assert_true((record.stdev_at(0, 0) - 4.5).abs() < 1.0e-9)
  assert_eq(record.npix[0][0], 25)
}

///|
test "cel parses mask and outliers" {
  let record = @src.cel_read(@src.cel_sample_text()) catch {
    _ => abort("failed to parse sample CEL")
  }
  assert_eq(record.mask, [3])
  assert_eq(record.outliers, [5])
  assert_eq(record.modified, [])
}

///|
test "cel header lookup" {
  let record = @src.cel_read(@src.cel_sample_text()) catch {
    _ => abort("failed to parse sample CEL")
  }
  assert_eq(record.get_header("Algorithm"), Some("Percentile"))
  assert_eq(record.get_header("MissingKey"), None)
}

///|
test "cel statistics" {
  let record = @src.cel_read(@src.cel_sample_text()) catch {
    _ => abort("failed to parse sample CEL")
  }
  // Mean of [91.5, 88.0, 102.3, 141.2, 130.0, 120.9] = 673.9 / 6
  let mean = record.mean_intensity()
  assert_true((mean - 112.31666666666667).abs() < 1.0e-6)
  let (lo, hi) = record.intensity_range()
  assert_true((lo - 88.0).abs() < 1.0e-9)
  assert_true((hi - 141.2).abs() < 1.0e-9)
}

///|
test "cel percent masked and outliers" {
  let record = @src.cel_read(@src.cel_sample_text()) catch {
    _ => abort("failed to parse sample CEL")
  }
  assert_true((record.percent_masked() - 1.0 / 6.0).abs() < 1.0e-9)
  assert_true((record.percent_outliers() - 1.0 / 6.0).abs() < 1.0e-9)
}

///|
test "cel percentile and background" {
  let record = @src.cel_read(@src.cel_sample_text()) catch {
    _ => abort("failed to parse sample CEL")
  }
  let p50 = record.intensity_percentile(50.0)
  // Sorted: 88.0, 91.5, 102.3, 120.9, 130.0, 141.2 -> median ~111.6
  assert_true(p50 > 100.0 && p50 < 125.0)
  let bg = record.background()
  assert_true(bg >= 88.0 && bg <= 95.0)
}

///|
test "cel dat header grid parsing" {
  let (cols, rows) = @src.cel_parse_dat_header(
    "[0..3]  712:712  712  712  Expressions",
  )
  assert_eq(cols, 712)
  assert_eq(rows, 712)
}

///|
test "cel empty input raises" {
  let raised = try {
    let _ = @src.cel_read("")
    false
  } catch {
    _ => true
  }
  assert_true(raised)
}