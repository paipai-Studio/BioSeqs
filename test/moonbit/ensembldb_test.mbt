///|
/// Tests for ensembldb module.

test "EnsDb creation" {
  let db = @src.create_example_ensdb()
  assert_eq(@src.edb_num_genes(db), 3)
  assert_eq(@src.edb_num_transcripts(db), 4)
}

test "edb_get_gene_by_id" {
  let db = @src.create_example_ensdb()
  let gene = @src.edb_get_gene_by_id(db, "ENSG000001")
  assert_true(gene is Some(_))
}

test "edb_get_gene_by_name" {
  let db = @src.create_example_ensdb()
  let genes = @src.edb_get_gene_by_name(db, "ACTB")
  assert_eq(genes.length(), 1)
}

test "edb_get_transcripts_by_gene" {
  let db = @src.create_example_ensdb()
  let txs = @src.edb_get_transcripts_by_gene(db, "ENSG000001")
  assert_eq(txs.length(), 2)
}

test "edb_filter_by_chromosome" {
  let db = @src.create_example_ensdb()
  let filtered = @src.edb_filter_by_chromosome(db, "17")
  assert_eq(@src.edb_num_genes(filtered), 2)
}

test "edb_get_gene_length" {
  let db = @src.create_example_ensdb()
  let length = @src.edb_get_gene_length(db, "ENSG000001")
  assert_eq(length, 6257)
}