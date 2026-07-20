///|
fn main {
  println("=== GEOquery GEO Series Matrix Parsing Example ===")
  println("Inspired by Bioconductor GEOquery package")
  println("")

  let geo_content = "!Series_geo_accession = GSE12345\n" +
    "!Platform_geo_accession = GPL96\n" +
    "!Sample_geo_accession = GSM123456\n" +
    "!Sample_geo_accession = GSM123457\n" +
    "!Sample_title = Control sample 1\n" +
    "!Sample_title = Treatment sample 1\n" +
    "!series_matrix_table_begin\n" +
    "ID_REF\tGSM123456\tGSM123457\n" +
    "TP53\t10.5\t15.2\n" +
    "BRCA1\t8.3\t12.1\n" +
    "BRCA2\t6.7\t9.8\n" +
    "!series_matrix_table_end\n"

  let geo_data = @bio.parse_geo_series_matrix(geo_content)

  println("1. GEOData Summary")
  println(@bio.geo_data_summary(geo_data))
  println("")

  println("2. Parse Results")
  println("   Series ID: " + geo_data.series_id)
  println("   Platform: " + geo_data.platform_id)
  println("")

  println("3. Sample IDs")
  let mut i = 0
  while i < geo_data.sample_ids.length() {
    println("   - " + geo_data.sample_ids[i])
    i = i + 1
  }
  println("")

  println("4. Gene Names")
  i = 0
  while i < geo_data.gene_names.length() {
    println("   - " + geo_data.gene_names[i])
    i = i + 1
  }
  println("")

  println("5. Expression Data")
  i = 0
  while i < geo_data.data.length() && i < 3 {
    let row = geo_data.data[i]
    let mut row_str = "   " + geo_data.gene_names[i] + ": "
    let mut j = 0
    while j < row.length() {
      if j > 0 {
        row_str = row_str + ", "
      }
      row_str = row_str + row[j].to_string()
      j = j + 1
    }
    println(row_str)
    i = i + 1
  }
  println("")

  println("6. Convert to ExpressionSet")
  let eset = @bio.geo_data_to_expression_set(geo_data)
  println("   Features: " + @bio.expression_set_nrow(eset).to_string())
  println("   Samples: " + @bio.expression_set_ncol(eset).to_string())
  println("")

  println("7. Filter Genes (expression 5.0 - 20.0)")
  let filtered = @bio.geo_data_filter_genes(geo_data, 5.0, 20.0)
  println("   Filtered genes: " + filtered.gene_names.length().to_string())
  println("")

  println("8. Log2 Transform")
  let log2_data = @bio.geo_data_log2_transform(geo_data)
  println("   Log2 transformed data computed")
  println("")

  println("=== GEOquery SOFT Format Parsing Example ===")
  println("")

  let soft_content = "^SERIES = GSE12345\n" +
    "!Series_geo_accession = GSE12345\n" +
    "!Series_title = Example dataset\n" +
    "!Series_summary = A test dataset for GEOquery\n" +
    "^SAMPLE = GSM123456\n" +
    "!Sample_platform_id = GPL96\n" +
    "^SAMPLE = GSM123457\n" +
    "!Sample_platform_id = GPL96\n"

  let series = @bio.parse_soft_format(soft_content)

  println("1. GEOSeries Summary")
  println(@bio.geo_series_summary(series))
}
