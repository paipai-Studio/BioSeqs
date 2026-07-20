///|
test "annotation_dbi_create_db" {
  let db = @src.create_annotation_db("org.Hs.eg.db")

  assert_eq(db.org_db, "org.Hs.eg.db")
  assert_eq(@src.nrow(db), 0)
}

///|
test "annotation_dbi_add_gene_info" {
  let db = @src.create_annotation_db("org.Hs.eg.db")
  let gene_info = @src.GeneInfo::new(
    "ENSG00000130203", "TP53", "Tumor protein p53", "chr17", 7661779, 7687550, "+",
    "protein_coding",
  )

  @src.add_gene_info(db, gene_info)

  assert_eq(@src.nrow(db), 1)
}

///|
test "annotation_dbi_get_gene_info" {
  let db = @src.create_annotation_db("org.Hs.eg.db")
  let gene_info = @src.GeneInfo::new(
    "ENSG00000130203", "TP53", "Tumor protein p53", "chr17", 7661779, 7687550, "+",
    "protein_coding",
  )

  @src.add_gene_info(db, gene_info)

  let result = @src.annotation_dbi_get_gene_info(db, "ENSG00000130203")
  assert_true(result.is_some())
}

///|
test "annotation_dbi_map_ids" {
  let db = @src.create_annotation_db("org.Hs.eg.db")
  let mappings = Map(
    [("TP53", "ENSG00000130203"), ("BRCA1", "ENSG00000012048")],
    capacity=2,
  )

  @src.add_id_map(db, "symbol", mappings)

  let result = @src.map_ids(db, ["TP53", "BRCA1"], "symbol", "ensembl")
  assert_eq(result.length(), 2)
}

///|
test "annotation_dbi_select" {
  let db = @src.create_annotation_db("org.Hs.eg.db")
  let gene_info = @src.GeneInfo::new(
    "ENSG00000130203", "TP53", "Tumor protein p53", "chr17", 7661779, 7687550, "+",
    "protein_coding",
  )

  @src.add_gene_info(db, gene_info)

  let result = @src.select(db, ["ENSG00000130203"], ["gene_id", "symbol"])
  assert_eq(result.length(), 1)
}

///|
test "annotation_dbi_add_go_terms" {
  let db = @src.create_annotation_db("org.Hs.eg.db")
  let go_terms = [
    @src.AnnotationGOTerm::new("GO:0005515", "protein binding", "IEA", "MF"),
  ]

  @src.add_go_terms(db, "ENSG00000130203", go_terms)

  let result = @src.get_go_annotations(db, ["ENSG00000130203"])
  assert_eq(result.length(), 1)
}

///|
test "annotation_dbi_add_kegg_terms" {
  let db = @src.create_annotation_db("org.Hs.eg.db")
  let kegg_terms = [
    @src.AnnotationKEGGTerm::new("hsa04115", "p53 signaling pathway", "Pathway"),
  ]

  @src.add_kegg_terms(db, "ENSG00000130203", kegg_terms)

  let result = @src.get_kegg_annotations(db, ["ENSG00000130203"])
  assert_eq(result.length(), 1)
}

///|
test "annotation_dbi_head" {
  let db = @src.create_annotation_db("org.Hs.eg.db")
  let gene_info1 = @src.GeneInfo::new(
    "ENSG00000130203", "TP53", "Tumor protein p53", "chr17", 7661779, 7687550, "+",
    "protein_coding",
  )
  let gene_info2 = @src.GeneInfo::new(
    "ENSG00000012048", "BRCA1", "Breast cancer type 1 susceptibility protein", "chr17",
    43044295, 43125370, "-", "protein_coding",
  )

  @src.add_gene_info(db, gene_info1)
  @src.add_gene_info(db, gene_info2)

  let head_result = @src.head(db, 1)
  assert_eq(head_result.length(), 1)
}

///|
test "annotation_dbi_summarize" {
  let db = @src.create_annotation_db("org.Hs.eg.db")
  let gene_info = @src.GeneInfo::new(
    "ENSG00000130203", "TP53", "Tumor protein p53", "chr17", 7661779, 7687550, "+",
    "protein_coding",
  )

  @src.add_gene_info(db, gene_info)

  let summary = @src.summarize_db(db)
  assert_true(summary.length() > 0)
}
