///|
test "gmt_gene_set_new" {
  let genes = ["TP53", "BRCA1", "EGFR", "MYC"]
  let gs = @src.GmtGeneSet::new("oncogenes", "Cancer genes", genes)
  assert_eq(gs.name, "oncogenes")
  assert_eq(gs.description, "Cancer genes")
  assert_eq(gs.gene_ids.length(), 4)
}

///|
test "gmt_gene_set_size" {
  let genes = ["A", "B", "C"]
  let gs = @src.GmtGeneSet::new("t", "d", genes)
  assert_eq(@src.GmtGeneSet::size(gs), 3)
}

///|
test "gmt_gene_set_has_gene" {
  let genes = ["TP53", "BRCA1"]
  let gs = @src.GmtGeneSet::new("g", "d", genes)
  assert_true(@src.GmtGeneSet::has_gene(gs, "TP53"))
  assert_true(@src.GmtGeneSet::has_gene(gs, "brca1")) // case insensitive
  assert_true(!@src.GmtGeneSet::has_gene(gs, "EGFR"))
}

///|
test "gmt_gene_set_with_annotation" {
  let genes = ["G1"]
  let ct = @src.GeneSetCollectionType::from_string("canonical")
  let gs = @src.GmtGeneSet::with_annotation(
    "g",
    "desc",
    genes,
    ct,
    "human",
    "GS123",
  )
  assert_eq(gs.name, "g")
  assert_eq(gs.organism, "human")
  assert_eq(gs.id, "GS123")
}

///|
test "gene_set_overlap_union_intersect" {
  let a = @src.GmtGeneSet::new("a", "d", ["1", "2", "3"])
  let b = @src.GmtGeneSet::new("b", "d", ["2", "3", "4"])
  let ov = @src.gene_set_overlap(a, b)
  assert_eq(ov, 2)

  let u = @src.gene_set_union_arr(a, b)
  assert_eq(u.length(), 4)

  let i_arr = @src.gene_set_intersect_arr(a, b)
  assert_eq(i_arr.length(), 2)

  let j = @src.gene_set_jaccard(a, b)
  assert_true((j - 0.5).abs() < 0.001)
}

///|
test "gene_set_overlap_coef" {
  let a = @src.GmtGeneSet::new("a", "d", ["1", "2", "3"])
  let b = @src.GmtGeneSet::new("b", "d", ["1", "2", "3", "4", "5"])
  let oc = @src.gene_set_overlap_coef(a, b)
  // overlap = 3, min(|a|, |b|) = 3, overlap_coef = 1
  assert_true((oc - 1.0).abs() < 0.001)
}

///|
test "gene_set_setdiff" {
  let a = @src.GmtGeneSet::new("a", "d", ["1", "2", "3", "4"])
  let b = @src.GmtGeneSet::new("b", "d", ["3", "4", "5"])
  let diff = @src.gene_set_setdiff(a, b)
  assert_eq(diff.length(), 2)
}

///|
test "gmt_gene_set_collection_new" {
  let g1 = @src.GmtGeneSet::new("s1", "d1", ["A", "B"])
  let g2 = @src.GmtGeneSet::new("s2", "d2", ["C", "D", "E"])
  let col = @src.GmtGeneSetCollection::new([g1, g2])
  assert_eq(col.size, 2)
  assert_eq(col.sets.length(), 2)
  assert_eq(col.names.length(), 2)
  assert_eq(col.names[0], "s1")
}

///|
test "gmt_gene_set_collection_push_s" {
  let col = @src.GmtGeneSetCollection::new([])
  assert_eq(col.size, 0)
  let gs = @src.GmtGeneSet::new("s1", "d", ["X"])
  let col2 = @src.GmtGeneSetCollection::push_s(col, gs)
  assert_eq(col2.size, 1)
  assert_eq(col2.names[0], "s1")
}

///|
test "gmt_gene_set_collection_by_name_and_at" {
  let g1 = @src.GmtGeneSet::new("s1", "d1", ["A", "B"])
  let g2 = @src.GmtGeneSet::new("s2", "d2", ["C"])
  let col = @src.GmtGeneSetCollection::new([g1, g2])
  let found = @src.GmtGeneSetCollection::by_name(col, "s1")
  match found {
    Some(gs) => assert_eq(gs.name, "s1")
    None => assert_true(false)
  }
  let at_gs = @src.GmtGeneSetCollection::at(col, 1)
  match at_gs {
    Some(gs) => assert_eq(gs.name, "s2")
    None => assert_true(false)
  }
}

///|
test "gmt_gene_set_collection_summary_and_sizes" {
  let g1 = @src.GmtGeneSet::new("s1", "d1", ["A", "B", "C"])
  let g2 = @src.GmtGeneSet::new("s2", "d2", ["D", "E"])
  let col = @src.GmtGeneSetCollection::new([g1, g2])
  let sz = @src.GmtGeneSetCollection::all_sizes(col)
  assert_eq(sz.length(), 2)
  assert_eq(sz[0], 3)
  assert_eq(sz[1], 2)
  let s = @src.GmtGeneSetCollection::summary(col)
  assert_true(s.contains("2 sets"))
}

///|
test "parse_gmt_basic" {
  let gmt = "set1\tdesc1\tGeneA\tGeneB\nset2\tdesc2\tGeneC\n"
  let col = @src.parse_gmt(gmt)
  assert_eq(col.size, 2)
  assert_eq(col.sets[0].name, "set1")
  assert_eq(col.sets[0].gene_ids.length(), 2)
  assert_eq(col.sets[1].name, "set2")
  assert_eq(col.sets[1].gene_ids.length(), 1)
}

///|
test "parse_gmx_basic" {
  let gmx = "set1\tset2\ndesc1\tdesc2\nGeneA\tGeneC\nGeneB\t\n"
  let col = @src.parse_gmx(gmx)
  assert_eq(col.size, 2)
  assert_eq(col.sets[0].name, "set1")
  assert_eq(col.sets[0].gene_ids.length(), 2)
}

///|
test "write_gmt_roundtrip" {
  let g1 = @src.GmtGeneSet::new("s1", "d1", ["A", "B"])
  let g2 = @src.GmtGeneSet::new("s2", "d2", ["C"])
  let col = @src.GmtGeneSetCollection::new([g1, g2])
  let gmt_str = @src.write_gmt(col)
  let col2 = @src.parse_gmt(gmt_str)
  assert_eq(col2.size, 2)
  assert_eq(col2.sets[0].name, "s1")
  assert_eq(col2.sets[1].name, "s2")
}

///|
test "write_gmx_roundtrip" {
  let g1 = @src.GmtGeneSet::new("s1", "d1", ["A", "B", "C"])
  let g2 = @src.GmtGeneSet::new("s2", "d2", ["D", "E"])
  let col = @src.GmtGeneSetCollection::new([g1, g2])
  let gmx_str = @src.write_gmx(col)
  assert_true(gmx_str.starts_with("s1\ts2"))
}

///|
test "filter_collection_by_gene_and_size" {
  let g1 = @src.GmtGeneSet::new("s1", "d", ["A", "B", "C"])
  let g2 = @src.GmtGeneSet::new("s2", "d", ["D", "E"])
  let g3 = @src.GmtGeneSet::new("s3", "d", ["B", "D", "F"])
  let col = @src.GmtGeneSetCollection::new([g1, g2, g3])
  let filtered = @src.collection_filter_by_gene(col, ["A", "B", "D", "F"])
  assert_eq(filtered.size, 3)
  assert_eq(filtered.sets[0].gene_ids.length(), 2)

  let by_size = @src.collection_filter_by_size(col, 2, 2)
  assert_eq(by_size.size, 1)
  assert_eq(by_size.sets[0].name, "s2")
}

///|
test "collection_gene_union" {
  let g1 = @src.GmtGeneSet::new("s1", "d", ["A", "B", "C"])
  let g2 = @src.GmtGeneSet::new("s2", "d", ["B", "C", "D"])
  let col = @src.GmtGeneSetCollection::new([g1, g2])
  let union = @src.collection_gene_union(col)
  assert_eq(union.length(), 4)
}

///|
test "collection_pairwise_jaccard" {
  let g1 = @src.GmtGeneSet::new("s1", "d", ["A", "B", "C"])
  let g2 = @src.GmtGeneSet::new("s2", "d", ["B", "C", "D"])
  let g3 = @src.GmtGeneSet::new("s3", "d", ["X", "Y"])
  let col = @src.GmtGeneSetCollection::new([g1, g2, g3])
  let jm = @src.collection_pairwise_jaccard(col)
  assert_eq(jm.length(), 3)
  assert_eq(jm[0].length(), 3)
  // jm[0][0] = self jaccard = 1
  assert_true((jm[0][0] - 1.0).abs() < 0.001)
  // jm[0][1] = overlap 2/4 = 0.5
  assert_true((jm[0][1] - 0.5).abs() < 0.01)
  // jm[0][2] = no overlap = 0
  assert_true(jm[0][2] < 0.001)
}

///|
test "gsea_as_plain_name_genes" {
  let gs = @src.GmtGeneSet::new("myset", "mydesc", ["G1", "G2"])
  let (name, genes) = @src.gsea_as_plain_name_genes(gs)
  assert_eq(name, "myset")
  assert_eq(genes.length(), 2)
}

///|
test "gene_set_collection_type_show_and_from" {
  let ct = @src.GeneSetCollectionType::from_string("canonical")
  assert_eq(@src.GeneSetCollectionType::show(ct), "Canonical")
  let parsed = @src.GeneSetCollectionType::from_string("kegg")
  assert_eq(@src.GeneSetCollectionType::show(parsed), "KEGG")
}
