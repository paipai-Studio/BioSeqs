///|
test "genomic_features_create_example" {
  let annotation = @src.create_example_annotation()
  assert_eq(@src.genome_annotation_count_genes(annotation), 5)
}

///|
test "genomic_features_get_gene_by_id" {
  let annotation = @src.create_example_annotation()
  let gene = @src.genome_annotation_get_gene_by_id(annotation, "ENSG00000130203")
  assert_true(gene.is_some())
  if gene.is_some() {
    assert_eq(gene.unwrap().gene_name, "TP53")
  }
}

///|
test "genomic_features_get_gene_by_name" {
  let annotation = @src.create_example_annotation()
  let genes = @src.genome_annotation_get_gene_by_name(annotation, "TP53")
  assert_eq(genes.length(), 1)
}

///|
test "genomic_features_genes_on_chromosome" {
  let annotation = @src.create_example_annotation()
  let genes_chr = @src.genome_annotation_genes_on_chromosome(annotation, "chr17")
  assert_eq(genes_chr.length(), 2)
}

///|
test "genomic_features_get_transcripts_by_gene" {
  let annotation = @src.create_example_annotation()
  let transcripts = @src.genome_annotation_get_transcripts_by_gene(annotation, "ENSG00000130203")
  assert_eq(transcripts.length(), 1)
}

///|
test "genomic_features_genes_in_region" {
  let annotation = @src.create_example_annotation()
  let genes_region = @src.genome_annotation_genes_in_region(annotation, "chr17", 7000000, 8000000)
  assert_eq(genes_region.length(), 1)
}

///|
test "genomic_features_count_transcripts" {
  let annotation = @src.create_example_annotation()
  assert_eq(@src.genome_annotation_count_transcripts(annotation), 4)
}

///|
test "genomic_features_summarize" {
  let annotation = @src.create_example_annotation()
  let summary = @src.genome_annotation_summarize(annotation)
  assert_true(@src.contains(summary, "Genes"))
}

///|
test "genomic_features_parse_gtf" {
  let gtf_content = 
    "chr17\tENSEMBL\tgene\t7661779\t7687550\t.\t-\t.\tgene_id \"ENSG00000130203\"; gene_name \"TP53\"; gene_biotype \"protein_coding\"; gene_description \"Tumor protein p53\";\n"
  
  let annotation = @src.gf_parse_gtf(gtf_content)
  assert_eq(@src.genome_annotation_count_genes(annotation), 1)
}

///|
test "genomic_features_gene_to_string" {
  let annotation = @src.create_example_annotation()
  let gene = @src.genome_annotation_get_gene_by_id(annotation, "ENSG00000130203")
  if gene.is_some() {
    let str = gene.unwrap().to_string()
    assert_true(@src.contains(str, "TP53"))
  }
}
