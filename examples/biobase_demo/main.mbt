///|
fn main {
  println("=== Biobase ExpressionSet Example ===")
  println("Inspired by Bioconductor Biobase package")
  println("")

  let assay_data = [
    [1.2, 2.3, 3.4, 4.5],
    [2.1, 3.2, 4.3, 5.4],
    [3.0, 4.1, 5.2, 6.3],
    [4.5, 5.6, 6.7, 7.8],
    [5.4, 6.5, 7.6, 8.7],
  ]

  let pheno_data : Map[String, Array[String]] = Map([], capacity=3)
  pheno_data["sample_name"] = ["sample1", "sample2", "sample3", "sample4"]
  pheno_data["group"] = ["control", "control", "treatment", "treatment"]
  pheno_data["sex"] = ["male", "female", "male", "female"]

  let feature_data : Map[String, Array[String]] = Map([], capacity=2)
  feature_data["feature_name"] = ["TP53", "BRCA1", "BRCA2", "ATM", "CDKN1A"]
  feature_data["chromosome"] = ["chr17", "chr17", "chr13", "chr11", "chr6"]

  let eset = @bio.ExpressionSet::new(
    assay_data,
    pheno_data,
    feature_data,
    "Homo sapiens",
    Map([]),
  )

  println("1. ExpressionSet Summary")
  println(@bio.expression_set_summary(eset))
  println("")

  println("2. Dimensions")
  println("   Rows (features): " + @bio.expression_set_nrow(eset).to_string())
  println("   Columns (samples): " + @bio.expression_set_ncol(eset).to_string())
  println("")

  println("3. Sample Names")
  let sample_names = @bio.expression_set_sample_names(eset)
  let mut i = 0
  while i < sample_names.length() {
    println("   - " + sample_names[i])
    i = i + 1
  }
  println("")

  println("4. Feature Names")
  let feature_names = @bio.expression_set_feature_names(eset)
  i = 0
  while i < feature_names.length() {
    println("   - " + feature_names[i])
    i = i + 1
  }
  println("")

  println("5. Subset (rows 0-1, columns 0-1)")
  let subset = @bio.expression_set_subset(eset, [0, 1], [0, 1])
  println("   Subset rows: " + @bio.expression_set_nrow(subset).to_string())
  println("   Subset cols: " + @bio.expression_set_ncol(subset).to_string())
  println("")

  println("6. Filter by row (max expression >= 5.0)")
  let filtered = @bio.expression_set_filter_by_row(eset, 5.0)
  println("   Filtered rows: " + @bio.expression_set_nrow(filtered).to_string())
  println("")

  println("7. Normalize (z-score)")
  let normalized = @bio.expression_set_normalize(eset)
  println("   Normalized assay data computed")
  println("")

  println("8. Log2 Transform")
  let log2_transformed = @bio.expression_set_log2(eset)
  println("   Log2 transformed assay data computed")
  println("")

  println("=== AnnotatedDataFrame Example ===")
  println("")

  let adf_data : Map[String, Array[String]] = Map([], capacity=3)
  adf_data["ID"] = ["1", "2", "3", "4"]
  adf_data["Group"] = ["A", "A", "B", "B"]
  adf_data["Treatment"] = ["Control", "Control", "Drug", "Drug"]

  let var_metadata : Map[String, String] = Map([], capacity=3)
  var_metadata["ID"] = "Sample identifier"
  var_metadata["Group"] = "Experimental group"
  var_metadata["Treatment"] = "Treatment type"

  let adf = @bio.AnnotatedDataFrame::new(adf_data, var_metadata)

  println("1. AnnotatedDataFrame Summary")
  println(@bio.annotated_data_frame_summary(adf))
  println("")

  println("2. Dimensions")
  println("   Rows: " + @bio.annotated_data_frame_nrow(adf).to_string())
  println("   Columns: " + @bio.annotated_data_frame_ncol(adf).to_string())
}
