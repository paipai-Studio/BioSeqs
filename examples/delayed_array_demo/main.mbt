///|
fn main {
  let darr = @bio.da_create_example_data()

  println(
    "DelayedArray created: " +
    darr.nrow().to_string() +
    " x " +
    darr.ncol().to_string(),
  )

  let row_sums = @bio.da_row_sums(darr)
  println("Row sums computed: " + row_sums.length().to_string() + " values")

  let col_means = @bio.da_col_means(darr)
  println(
    "Column means computed: " + col_means.length().to_string() + " values",
  )

  let transposed = @bio.da_transpose(darr)
  println(
    "Transposed array: " +
    transposed.nrow().to_string() +
    " x " +
    transposed.ncol().to_string(),
  )

  let subset = @bio.da_subset_rows(darr, [0, 1, 2, 3, 4])
  println("Row subset: " + subset.nrow().to_string() + " rows")

  let summary = @bio.da_summary(darr)
  println(summary)
}
