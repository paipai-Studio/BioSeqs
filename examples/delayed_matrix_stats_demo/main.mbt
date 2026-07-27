///|
/// DelayedMatrixStats Demo - Comprehensive demonstration of the
/// DelayedMatrixStats module for memory-efficient row/column statistics
/// on DelayedArray matrices, inspired by Bioconductor's DelayedMatrixStats.

///|
fn main {
  println("=== DelayedMatrixStats Demo ===\n")

  let row_names = ["gene_a", "gene_b", "gene_c", "gene_d"]
  let col_names = ["sample_1", "sample_2", "sample_3", "sample_4"]

  println("=== 1. Creating a DelayedMatrix ===")
  let data : Array[Array[Double]] = [
    [1.0, 2.0, 3.0, 4.0],
    [2.5, 3.5, 4.5, 5.5],
    [10.0, 20.0, 30.0, 40.0],
    [5.0, 7.0, 9.0, 11.0],
  ]
  let matrix = @bio.create_delayed_matrix(data, row_names, col_names)
  println("Dimensions: " + matrix.nrow().to_string() +
    " x " + matrix.ncol().to_string())
  let dims = matrix.dim()
  println("Dim vector: [" + dims[0].to_string() +
    ", " + dims[1].to_string() + "]")
  println("Element [0,0]: " + matrix.get(0, 0).to_string())
  println("Element [2,3]: " + matrix.get(2, 3).to_string())

  println("\n=== 2. Row and Column Means ===")
  let row_means = @bio.dms_row_means(matrix)
  let mut ri = 0
  while ri < row_means.length() {
    println("  row mean [" + row_names[ri] + "]: " +
      row_means[ri].to_string())
    ri = ri + 1
  }
  let col_means = @bio.dms_col_means(matrix)
  let mut ci = 0
  while ci < col_means.length() {
    println("  col mean [" + col_names[ci] + "]: " +
      col_means[ci].to_string())
    ci = ci + 1
  }

  println("\n=== 3. Row and Column Medians ===")
  let row_medians = @bio.row_medians(matrix)
  ri = 0
  while ri < row_medians.length() {
    println("  row median [" + row_names[ri] + "]: " +
      row_medians[ri].to_string())
    ri = ri + 1
  }
  let col_medians = @bio.col_medians(matrix)
  ci = 0
  while ci < col_medians.length() {
    println("  col median [" + col_names[ci] + "]: " +
      col_medians[ci].to_string())
    ci = ci + 1
  }

  println("\n=== 4. Row and Column Variances ===")
  let row_vars = @bio.dms_row_vars(matrix)
  ri = 0
  while ri < row_vars.length() {
    println("  row var [" + row_names[ri] + "]: " +
      row_vars[ri].to_string())
    ri = ri + 1
  }
  let col_vars = @bio.col_vars(matrix)
  ci = 0
  while ci < col_vars.length() {
    println("  col var [" + col_names[ci] + "]: " +
      col_vars[ci].to_string())
    ci = ci + 1
  }

  println("\n=== 5. Handling NA Values ===")
  let na_data : Array[Array[Double]] = [
    [1.0, Double::nan(), 3.0, 4.0],
    [Double::nan(), 5.5, Double::nan(), 7.5],
    [10.0, 20.0, Double::nan(), 40.0],
    [5.0, 7.0, 9.0, 11.0],
  ]
  let na_matrix = @bio.create_delayed_matrix(
    na_data, row_names, col_names,
  )
  println("Dimensions with NAs: " +
    na_matrix.nrow().to_string() + " x " +
    na_matrix.ncol().to_string())

  let na_row_means = @bio.dms_row_means(na_matrix)
  ri = 0
  while ri < na_row_means.length() {
    println("  row mean [" + row_names[ri] + "]: " +
      na_row_means[ri].to_string())
    ri = ri + 1
  }

  let na_col_means = @bio.dms_col_means(na_matrix)
  ci = 0
  while ci < na_col_means.length() {
    println("  col mean [" + col_names[ci] + "]: " +
      na_col_means[ci].to_string())
    ci = ci + 1
  }

  println("\n=== 6. Subsetting the Matrix ===")
  let subset_rows = [0, 2]
  let subset_cols = [1, 3]
  let sub_matrix = @bio.delayed_matrix_subset(
    na_matrix, subset_rows, subset_cols,
  )
  println("Subset dimensions: " +
    sub_matrix.nrow().to_string() + " x " +
    sub_matrix.ncol().to_string())
  let sub_dims = sub_matrix.dim()
  println("Subset dim vector: [" +
    sub_dims[0].to_string() + ", " +
    sub_dims[1].to_string() + "]")

  let sub_row_means = @bio.dms_row_means(sub_matrix)
  ri = 0
  while ri < sub_row_means.length() {
    println("  subset row mean: " + sub_row_means[ri].to_string())
    ri = ri + 1
  }

  println("\n=== 7. Row n (non-NA count) ===")
  let row_counts = @bio.row_n(na_matrix)
  ri = 0
  while ri < row_counts.length() {
    println("  row n [" + row_names[ri] + "]: " +
      row_counts[ri].to_string())
    ri = ri + 1
  }
  let col_counts = @bio.col_n(na_matrix)
  ci = 0
  while ci < col_counts.length() {
    println("  col n [" + col_names[ci] + "]: " +
      col_counts[ci].to_string())
    ci = ci + 1
  }

  println("\n=== 8. Additional Statistics ===")
  let row_sums = @bio.row_allsums(matrix)
  ri = 0
  while ri < row_sums.length() {
    println("  row sum [" + row_names[ri] + "]: " +
      row_sums[ri].to_string())
    ri = ri + 1
  }

  println("\n=== DelayedMatrixStats Demo Done ===")
}
