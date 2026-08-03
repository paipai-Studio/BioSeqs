///|
/// MatrixGenerics module demo - Bioconductor-style row/col summary statistics.
///
/// Run: moon run examples/matrix_generics_demo
fn main {
  println("=== MatrixGenerics Module Demo ===")
  println(" efficient row and column summary statistics for matrix-like objects")
  println(" Reference: Bioconductor MatrixGenerics package")
  println("")

  // ------------------------------------------------------------------
  // 1. Create a sample matrix
  // ------------------------------------------------------------------
  println("=== 1. Create sample matrix ===")
  let mat = @bio.MgMatrix::new([[1.0, 2.0, 3.0], [4.0, 5.0, 6.0], [7.0, 8.0, 9.0]])
  println("Dimensions: \{@bio.MgMatrix::dim_rows(mat)} rows x \{@bio.MgMatrix::dim_cols(mat)} cols")
  println("Matrix data:")
  for i in 0..<@bio.MgMatrix::dim_rows(mat) {
    let row = @bio.MgMatrix::get_row(mat, i)
    println("  row \{i}: \{row}")
  }
  println("")

  // ------------------------------------------------------------------
  // 2. Row summary statistics
  // ------------------------------------------------------------------
  println("=== 2. Row summary statistics ===")
  let row_means = @bio.mg_row_means(mat)
  let row_sums = @bio.mg_row_sums(mat)
  let row_medians = @bio.mg_row_medians(mat)
  let row_mins = @bio.mg_row_mins(mat)
  let row_maxs = @bio.mg_row_maxs(mat)
  let row_vars = @bio.mg_row_vars(mat)
  let row_sds = @bio.mg_row_sds(mat)
  println("Row means:    \{row_means}")
  println("Row sums:     \{row_sums}")
  println("Row medians:  \{row_medians}")
  println("Row mins:     \{row_mins}")
  println("Row maxs:     \{row_maxs}")
  println("Row vars:     \{row_vars}")
  println("Row sds:      \{row_sds}")
  println("")

  // ------------------------------------------------------------------
  // 3. Column summary statistics
  // ------------------------------------------------------------------
  println("=== 3. Column summary statistics ===")
  let col_means = @bio.mg_col_means(mat)
  let col_sums = @bio.mg_col_sums(mat)
  let col_medians = @bio.mg_col_medians(mat)
  let col_mins = @bio.mg_col_mins(mat)
  let col_maxs = @bio.mg_col_maxs(mat)
  println("Col means:    \{col_means}")
  println("Col sums:     \{col_sums}")
  println("Col medians:  \{col_medians}")
  println("Col mins:     \{col_mins}")
  println("Col maxs:     \{col_maxs}")
  println("")

  // ------------------------------------------------------------------
  // 4. Row anys and alls with a threshold
  // ------------------------------------------------------------------
  println("=== 4. Row anys and row_alls (threshold = 5.0) ===")
  let row_anys = @bio.mg_row_anys(mat, 5.0)
  let row_alls = @bio.mg_row_alls(mat, 5.0)
  println("Row anys (>5.0): \{row_anys}")
  println("Row alls (>5.0): \{row_alls}")
  println("")

  // ------------------------------------------------------------------
  // 5. Block processing with mg_block_apply_rows
  // ------------------------------------------------------------------
  println("=== 5. Block processing (block_size=2, row_means) ===")
  let block_result = @bio.mg_block_apply_rows(mat, 2, (block) => {
    @bio.mg_row_means(block)
  })
  println("Block-applied row means: \{block_result}")
  println("")

  // ------------------------------------------------------------------
  // 6. Standalone statistical functions
  // ------------------------------------------------------------------
  println("=== 6. Standalone statistical functions ===")
  let sample = [3.0, 1.0, 4.0, 1.0, 5.0, 9.0, 2.0, 6.0]
  println("Sample data: \{sample}")
  println("Mean:   \{@bio.mg_mean(sample)}")
  println("Median: \{@bio.mg_median(sample)}")
  println("Min:    \{@bio.mg_min(sample)}")
  println("Max:    \{@bio.mg_max(sample)}")
  println("Var:    \{@bio.mg_variance(sample)}")

  println("=== Demo complete ===")
}