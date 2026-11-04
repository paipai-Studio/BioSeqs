///|
/// SummarizedExperiment Demo.
fn main {
  println("=== SummarizedExperiment Demo ===\n")

  let assays = Map([
    ("counts", [[1.0, 2.0, 3.0], [4.0, 5.0, 6.0], [7.0, 8.0, 9.0]]),
  ])
  let row_ranges = [("chr1", 1, 100), ("chr1", 200, 300), ("chr2", 100, 200)]
  let col_data = [
    Map([("sample", "s1"), ("group", "A")]),
    Map([("sample", "s2"), ("group", "A")]),
    Map([("sample", "s3"), ("group", "B")]),
  ]
  let metadata = Map([("experiment", "RNA-seq"), ("platform", "Illumina")])

  println("1. Creating SummarizedExperiment...")
  let se = @src.summarized_experiment(assays, row_ranges, col_data, metadata)
  println("   Rows: \{@src.se_nrow(se)}")
  println("   Cols: \{@src.se_ncol(se)}")
  println("   Assays: \{@src.se_assay_names(se).join(", ")}")

  println("\n2. Getting row and column means...")
  let row_means = @src.se_row_means(se, "counts")
  let col_means = @src.se_col_means(se, "counts")
  let row_str = row_means.map(fn(x) { x.to_string() }).join(", ")
  let col_str = col_means.map(fn(x) { x.to_string() }).join(", ")
  println("   Row means: \{row_str}")
  println("   Col means: \{col_str}")

  println("\n3. Subsetting rows [0, 2]...")
  let se_sub = @src.se_subset_rows(se, [0, 2])
  println("   New rows: \{@src.se_nrow(se_sub)}")

  println("\n4. Subsetting columns [0, 1]...")
  let se_sub2 = @src.se_subset_cols(se, [0, 1])
  println("   New cols: \{@src.se_ncol(se_sub2)}")

  println("\n=== SummarizedExperiment Demo Complete ===")
}
