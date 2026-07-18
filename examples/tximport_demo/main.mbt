///|
fn main {
  println("=== tximport Transcript Quantification Import Example ===")
  println("Inspired by Bioconductor tximport package")
  println("")

  let quant_sf_content1 = "# Name\tLength\tEffectiveLength\tTPM\tNumReads\n" +
    "ENST00000357654.5\t1500\t1480\t10.5\t150\n" +
    "ENST00000357655.3\t1200\t1180\t8.3\t120\n" +
    "ENST00000357656.2\t900\t880\t6.7\t95\n"

  let quant_sf_content2 = "# Name\tLength\tEffectiveLength\tTPM\tNumReads\n" +
    "ENST00000357654.5\t1500\t1480\t15.2\t220\n" +
    "ENST00000357655.3\t1200\t1180\t12.1\t175\n" +
    "ENST00000357656.2\t900\t880\t9.8\t140\n"

  let tx2gene = [
    @bio.tx2gene_new("ENST00000357654", "ENSG00000130203", "TP53"),
    @bio.tx2gene_new("ENST00000357655", "ENSG00000012048", "BRCA1"),
    @bio.tx2gene_new("ENST00000357656", "ENSG00000139618", "BRCA2"),
  ]

  println("1. Parse quant.sf File")
  let (tx_names, abundance, counts, length) = @bio.parse_quant_sf(
    quant_sf_content1,
  )
  println("   Transcripts: " + tx_names.length().to_string())
  println("")

  println("2. tx_import (multi-sample)")
  let txi = @bio.tx_import(
    [quant_sf_content1, quant_sf_content2],
    tx2gene,
    true,
  )
  println(@bio.tx_import_data_summary(txi))
  println("")

  println("3. Transcript-level Abundance")
  let mut i = 0
  while i < txi.tx_names.length() && i < 3 {
    println("   - " + txi.tx_names[i] + " -> " + txi.gene_names[i])
    i = i + 1
  }
  println("")

  println("4. Summarize to Gene Level")
  let gene_level = @bio.summarize_to_gene(txi, tx2gene, true)
  println(
    "   Gene-level data: " +
    gene_level.gene_names.length().to_string() +
    " genes",
  )
  println("")

  println("5. Gene Names")
  i = 0
  while i < gene_level.gene_names.length() {
    println("   - " + gene_level.gene_names[i])
    i = i + 1
  }
  println("")

  println("6. Convert to ExpressionSet")
  let eset_counts = @bio.tx_import_data_to_expression_set(gene_level, "counts")
  println("   Counts ExpressionSet:")
  println("     Features: " + @bio.expression_set_nrow(eset_counts).to_string())
  println("     Samples: " + @bio.expression_set_ncol(eset_counts).to_string())
  println("")

  let eset_tpm = @bio.tx_import_data_to_expression_set(gene_level, "abundance")
  println("   TPM ExpressionSet:")
  println("     Features: " + @bio.expression_set_nrow(eset_tpm).to_string())
  println("     Samples: " + @bio.expression_set_ncol(eset_tpm).to_string())
  println("")

  println("7. Filter Low Expression (min 100 counts)")
  let filtered = @bio.tx_import_filter_low_expression(txi, 100.0)
  println(
    "   After filtering: " +
    filtered.tx_names.length().to_string() +
    " transcripts",
  )
  println("")

  println("8. Sample Names")
  i = 0
  while i < txi.sample_names.length() {
    println("   - " + txi.sample_names[i])
    i = i + 1
  }
}
