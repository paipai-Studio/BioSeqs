fn main {
  println("=== Bioconductor-style GSEABase Demo ===\n")

  println("--- 1. Create Gene Sets ---")
  let ct = @bio.GeneSetCollectionType::from_string("hallmark")
  let gs1 = @bio.GmtGeneSet::with_annotation(
    "HALLMARK_APOPTOSIS",
    "Genes involved in apoptosis",
    ["TP53".to_string(), "BAX".to_string(), "BCL2".to_string(), "CASP3".to_string(), "CASP8".to_string()],
    ct,
    "human".to_string(),
    "H00001".to_string(),
  )
  let size1 = @bio.GmtGeneSet::size(gs1)
  println("  gs1: \{gs1.name} [genes: \{size1}]")
  println("       desc: \{gs1.description}")
  println("       organism: \{gs1.organism}")
  let has_tp53 = @bio.GmtGeneSet::has_gene(gs1, "TP53")
  let has_kras = @bio.GmtGeneSet::has_gene(gs1, "KRAS")
  println("       has TP53: \{has_tp53}")
  println("       has KRAS: \{has_kras}")

  let gs2 = @bio.GmtGeneSet::new(
    "HALLMARK_CELL_CYCLE",
    "Cell cycle related genes",
    ["CDK1".to_string(), "CCNB1".to_string(), "TP53".to_string(), "BCL2".to_string(), "MYC".to_string()],
  )
  let size2 = @bio.GmtGeneSet::size(gs2)
  println("\n  gs2: \{gs2.name} [genes: \{size2}]")

  println("\n--- 2. Gene Set Operations ---")
  let ov = @bio.gene_set_overlap(gs1, gs2)
  println("  overlap count (shared genes): \{ov}")
  let oc = @bio.gene_set_overlap_coef(gs1, gs2)
  println("  overlap coefficient: \{oc}")
  let jacc = @bio.gene_set_jaccard(gs1, gs2)
  println("  Jaccard index: \{jacc}")
  let diff_arr = @bio.gene_set_setdiff(gs1, gs2)
  println("  setdiff (gs1 - gs2) size: \{diff_arr.length()}")
  let union_arr = @bio.gene_set_union_arr(gs1, gs2)
  println("  union size: \{union_arr.length()}")
  let intersect_arr = @bio.gene_set_intersect_arr(gs1, gs2)
  println("  intersection size: \{intersect_arr.length()}")

  println("\n--- 3. Gene Set Collection ---")
  let col = @bio.GmtGeneSetCollection::new([gs1, gs2])
  println("  Collection size: \{col.size}")
  let sizes = @bio.GmtGeneSetCollection::all_sizes(col)
  println("  All sizes: [\{sizes[0]}, \{sizes[1]}]")
  println("  Collection summary:")
  println(@bio.GmtGeneSetCollection::summary(col))

  println("\n--- 4. Query Collection ---")
  let found = @bio.GmtGeneSetCollection::by_name(col, "HALLMARK_APOPTOSIS")
  match found {
    Some(g) => {
      let sg = @bio.GmtGeneSet::size(g)
      println("  Found by name: \{g.name} (genes: \{sg})")
    }
    None => println("  Not found")
  }
  let at_gs = @bio.GmtGeneSetCollection::at(col, 1)
  match at_gs {
    Some(g) => println("  At index 1: \{g.name}")
    None => println("  Out of range")
  }

  println("\n--- 5. Push additional sets ---")
  let gs3 = @bio.GmtGeneSet::new("HALLMARK_DNA_REPAIR", "DNA repair genes", ["BRCA1".to_string(), "BRCA2".to_string(), "ATM".to_string()])
  let col2 = @bio.GmtGeneSetCollection::push_s(col, gs3)
  println("  After push: \{col2.size} sets total")
  let mut names_s = "  All names: ".to_string()
  let mut ni = 0
  while ni < col2.names.length() {
    if ni > 0 { names_s = names_s + ", " }
    names_s = names_s + col2.names[ni]
    ni = ni + 1
  }
  println(names_s)

  println("\n--- 6. GMT Format Roundtrip ---")
  let gmt_str = @bio.write_gmt(col2)
  println("  GMT output:")
  println("    \{gmt_str}")
  let parsed = @bio.parse_gmt(gmt_str)
  println("  Re-parsed collection size: \{parsed.size}")

  println("\n--- 7. GeneSetCollectionType ---")
  let ct_kegg = @bio.GeneSetCollectionType::from_string("kegg")
  let show_kegg = @bio.GeneSetCollectionType::show(ct_kegg)
  println("  from_string('kegg') -> show: \{show_kegg}")
  let ct_other = @bio.GeneSetCollectionType::from_string("gmt_custom")
  let show_other = @bio.GeneSetCollectionType::show(ct_other)
  println("  from_string('gmt_custom') -> show: \{show_other}")

  println("\nGSEABase demo completed successfully!")
}
