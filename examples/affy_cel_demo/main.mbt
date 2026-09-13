///|
fn main {
  println("=== Bio.Affy.CelFile Demo ===")

  let record = @src.cel_read(@src.cel_sample_text()) catch {
    _ => abort("failed to parse sample CEL")
  }

  println("\n1. Header")
  println("  Version: " + record.version)
  println("  Grid: " + record.n_cols.to_string() + " cols x " + record.n_rows.to_string() + " rows")
  println("  Cells declared: " + record.n_cells.to_string() + ", read: " + record.cells_read.to_string())

  println("\n2. Intensity matrix")
  for y in 0..<record.n_rows {
    let mut row = ""
    for x in 0..<record.n_cols {
      row = row + record.intensity_at(x, y).to_string() + "  "
    }
    println("  row " + y.to_string() + ": " + row)
  }

  println("\n3. Probe-level statistics")
  println("  Mean intensity:   " + record.mean_intensity().to_string())
  println("  Median intensity: " + record.median_intensity().to_string())
  let (lo, hi) = record.intensity_range()
  println("  Range:            [" + lo.to_string() + ", " + hi.to_string() + "]")
  println("  Background (P5):  " + record.background().to_string())
  println("  IQR:              " + record.intensity_iqr().to_string())

  println("\n4. Quality flags")
  println("  Masked cells:    " + record.mask.length().to_string() + " (" + record.percent_masked().to_string() + ")")
  println("  Outlier cells:   " + record.outliers.length().to_string() + " (" + record.percent_outliers().to_string() + ")")
}