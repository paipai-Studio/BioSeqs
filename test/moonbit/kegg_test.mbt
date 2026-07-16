///|
test "kegg_gene_new" {
  let pathways = ["path:hsa00010", "path:hsa00020"]
  let dblinks = [("NCBI-GeneID", "12345"), ("UniProt", "P12345")]

  let gene = @src.KEGGGene::new(
    "hsa:12345", "GeneA", "Description of gene A", "K00001", "Homo sapiens", pathways,
    "1p36.11", dblinks,
  )

  assert_eq(gene.entry, "hsa:12345")
  assert_eq(gene.name, "GeneA")
  assert_eq(gene.pathways.length(), 2)
  assert_eq(gene.dblinks.length(), 2)
}

///|
test "kegg_parse_gene" {
  let text = "ENTRY       hsa:12345\nNAME        GeneA\nDEFINITION  Description of gene A\nORTHOLOGY   K00001\nORGANISM    Homo sapiens\nPATHWAY     path:hsa00010 Glycolysis / Gluconeogenesis\nPATHWAY     path:hsa00020 Citrate cycle\nPOSITION    1p36.11\nDBLINKS     NCBI-GeneID: 12345"

  let gene = @src.parse_kegg_gene(text)

  assert_true(gene.is_some())
  assert_eq(gene.unwrap().entry, "hsa:12345")
  assert_eq(gene.unwrap().name, "GeneA")
  assert_eq(gene.unwrap().pathways.length(), 2)
}

///|
test "kegg_parse_pathway" {
  let text = "ENTRY       path:hsa00010\nNAME        Glycolysis / Gluconeogenesis\nDEFINITION  Glycolysis is the metabolic pathway...\nCLASS       Carbohydrate metabolism\nGENE        hsa:12345  GeneA\nGENE        hsa:67890  GeneB\nCOMPOUND    C00033 Glucose\nENZYME      EC:2.7.1.1\nORGANISM    Homo sapiens"

  let pathway = @src.parse_kegg_pathway(text)

  assert_true(pathway.is_some())
  assert_eq(pathway.unwrap().entry, "path:hsa00010")
  assert_eq(pathway.unwrap().name, "Glycolysis / Gluconeogenesis")
  assert_eq(pathway.unwrap().genes.length(), 2)
}

///|
test "kegg_parse_compound" {
  let text = "ENTRY       C00033\nNAME        Glucose\nFORMULA     C6H12O6\nMASS        180.156\nPATHWAY     path:hsa00010\nENZYME      EC:2.7.1.1\nREACTION    R01786"

  let compound = @src.parse_kegg_compound(text)

  assert_true(compound.is_some())
  assert_eq(compound.unwrap().entry, "C00033")
  assert_eq(compound.unwrap().name, "Glucose")
  assert_true(compound.unwrap().mass == 180.156)
}

///|
test "kegg_parse_enzyme" {
  let text = "ENTRY       EC:2.7.1.1\nNAME        Hexokinase\nCLASS       Transferases; Transferring phosphorus-containing groups\nGENE        hsa:1234  HK1\nGENE        hsa:1235  HK2\nPATHWAY     path:hsa00010\nREACTION    R01786\nCOFACTOR    ATP"

  let enzyme = @src.parse_kegg_enzyme(text)

  assert_true(enzyme.is_some())
  assert_eq(enzyme.unwrap().entry, "EC:2.7.1.1")
  assert_eq(enzyme.unwrap().name, "Hexokinase")
  assert_eq(enzyme.unwrap().genes.length(), 2)
}

///|
test "kegg_pathway_gene_ids" {
  let genes = [("hsa:12345", "GeneA"), ("hsa:67890", "GeneB")]
  let pathway = @src.KEGGPathway::new(
    "path:hsa00010",
    "Glycolysis",
    "desc",
    "Carbohydrate",
    genes,
    [],
    [],
    [],
  )

  let ids = @src.kegg_pathway_gene_ids(pathway)

  assert_eq(ids.length(), 2)
  assert_eq(ids[0], "hsa:12345")
}

///|
test "kegg_gene_pathway_count" {
  let gene = @src.KEGGGene::new(
    "hsa:12345",
    "GeneA",
    "desc",
    "K00001",
    "Homo sapiens",
    ["path:hsa00010", "path:hsa00020"],
    "",
    [],
  )

  let count = @src.kegg_gene_pathway_count(gene)

  assert_eq(count, 2)
}

///|
test "kegg_compound_in_pathway" {
  let compound = @src.KEGGCompound::new(
    "C00033",
    "Glucose",
    "C6H12O6",
    180.156,
    ["path:hsa00010 Glycolysis", "path:hsa00020 Citrate"],
    [],
    [],
  )

  assert_true(@src.kegg_compound_in_pathway(compound, "path:hsa00010"))
  assert_true(@src.kegg_compound_in_pathway(compound, "path:hsa00020"))
  assert_true(!@src.kegg_compound_in_pathway(compound, "path:hsa00030"))
}
