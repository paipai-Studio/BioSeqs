///|
test "biomart_create_human" {
  let mart = @src.create_human_mart()
  assert_eq(mart.name, "hsapiens_gene_ensembl")
}

///|
test "biomart_get_mapping" {
  let mart = @src.create_human_mart()
  let mapping = mart.get_mapping("ENSG00000139618")
  assert_true(mapping is Some(_))
}

///|
test "biomart_get_annotation" {
  let mart = @src.create_human_mart()
  let desc = mart.get_annotation("ENSG00000139618", "description")
  assert_true(desc is Some(_))
}

///|
test "biomart_query" {
  let mart = @src.create_human_mart()
  let gene_ids = Array::new()
  gene_ids.push("ENSG00000139618")
  let attributes = Array::new()
  attributes.push("external_gene_name")
  attributes.push("description")
  let results = mart.query(gene_ids, attributes)
  assert_eq(results.length(), 2)
}
