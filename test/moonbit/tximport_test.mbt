///|
test "tximport_parse_quant_sf" {
  let content = "Name\tLength\tEffectiveLength\tTPM\tNumReads\n" +
    "tx1\t1000\t900\t10.5\t100\n" +
    "tx2\t2000\t1800\t20.3\t200\n"

  let (tx_names, abundance, counts, length) = @src.parse_quant_sf(content)

  assert_eq(tx_names.length(), 2)
  assert_eq(tx_names[0], "tx1")
  assert_eq(abundance[0], 10.5)
  assert_eq(counts[0], 100.0)
}

///|
test "tximport_tx_import" {
  let file1 = "Name\tLength\tEffectiveLength\tTPM\tNumReads\n" +
    "tx1\t1000\t900\t10.5\t100\n" +
    "tx2\t2000\t1800\t20.3\t200\n"

  let tx2gene = [
    @src.tx2gene_new("tx1", "gene1", "Gene1"),
    @src.tx2gene_new("tx2", "gene2", "Gene2"),
  ]

  let txi = @src.tx_import([file1], tx2gene, false)

  assert_eq(txi.tx_names.length(), 2)
  assert_eq(txi.sample_names.length(), 1)
}

///|
test "tximport_summarize_to_gene" {
  let abundance = [[10.5, 20.3], [15.2, 25.8]]
  let counts = [[100.0, 200.0], [150.0, 250.0]]
  let length = [[900.0, 1800.0], [950.0, 1850.0]]

  let txi = @src.tx_import_data_new(
    abundance,
    counts,
    length,
    ["tx1", "tx2", "tx1a", "tx2a"],
    ["s1", "s2"],
    ["Gene1", "Gene2", "Gene1", "Gene2"],
  )

  let tx2gene = [
    @src.tx2gene_new("tx1", "gene1", "Gene1"),
    @src.tx2gene_new("tx2", "gene2", "Gene2"),
  ]

  let gene_txi = @src.summarize_to_gene(txi, tx2gene, false)

  assert_eq(gene_txi.gene_names.length(), 2)
}

///|
test "tximport_filter_low_expression" {
  let counts = [[0.5, 0.0, 0.0], [5.0, 200.0, 300.0]]
  let abundance = [[0.5, 0.0, 0.0], [50.0, 100.0, 150.0]]
  let length = [[1000.0, 2000.0, 3000.0], [1000.0, 2000.0, 3000.0]]

  let txi = @src.tx_import_data_new(
    abundance,
    counts,
    length,
    ["tx1", "tx2", "tx3"],
    ["s1", "s2"],
    ["Gene1", "Gene2", "Gene3"],
  )

  let filtered = @src.tx_import_filter_low_expression(txi, 10.0)

  assert_eq(filtered.tx_names.length(), 2)
}

///|
test "tximport_tx_import_data_summary" {
  let abundance = [[10.5, 20.3]]
  let counts = [[100.0, 200.0]]
  let length = [[900.0, 1800.0]]

  let txi = @src.tx_import_data_new(
    abundance,
    counts,
    length,
    ["tx1", "tx2"],
    ["s1"],
    ["Gene1", "Gene2"],
  )

  let summary = @src.tx_import_data_summary(txi)

  assert_true(summary.contains("TxImportData Summary"))
  assert_true(summary.contains("Transcripts: 2"))
}
