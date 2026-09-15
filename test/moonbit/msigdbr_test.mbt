///|
/// Tests for msigdbr module (MSigDB gene set collections).

test "msig_category_code" {
  assert_eq(@src.MsigCategory::code(@src.MsigHallmark), "H")
  assert_eq(@src.MsigCategory::code(@src.MsigCurated), "C2")
  assert_eq(@src.MsigCategory::code(@src.MsigOntology), "C5")
}

test "msig_category_label" {
  let l = @src.MsigCategory::label(@src.MsigHallmark)
  assert_true(l.contains("Hallmark"))
}

test "msig_category_from_string" {
  let h = @src.msig_category_from_string("hallmark")
  assert_eq(@src.MsigCategory::code(h), "H")
  let c2 = @src.msig_category_from_string("C2")
  assert_eq(@src.MsigCategory::code(c2), "C2")
  let go = @src.msig_category_from_string("ontology")
  assert_eq(@src.MsigCategory::code(go), "C5")
  let other = @src.msig_category_from_string("unknown")
  assert_eq(@src.MsigCategory::code(other), "other")
}

test "msig_hallmark" {
  let col = @src.msig_hallmark()
  assert_true(col.size() >= 5)
  assert_true(col.category is @src.MsigHallmark)
  for s in col.sets {
    assert_true(s.size() >= 5)
    assert_true(s.name.has_prefix("HALLMARK_"))
  }
}

test "msig_curated_canonical" {
  let col = @src.msig_curated_canonical()
  assert_true(col.size() >= 3)
  assert_true(col.category is @src.MsigCurated)
}

test "msig_ontology_go_bp" {
  let col = @src.msig_ontology_go_bp()
  assert_true(col.size() >= 3)
  assert_true(col.category is @src.MsigOntology)
}

test "msig_catalog" {
  let catalog = @src.msig_catalog()
  assert_eq(catalog.length(), 3)
  let total = @src.msig_catalog_size(catalog)
  assert_true(total >= 11)
}

test "msig_lookup" {
  let catalog = @src.msig_catalog()
  let gs = @src.msig_lookup(catalog, "HALLMARK_TNFA_SIGNALING_VIA_NFKB")
  assert_true(gs is Some(_))
  match gs {
    Some(s) => assert_true(s.size() >= 10)
    None => fail("expected gene set")
  }
  assert_true(@src.msig_lookup(catalog, "NONEXISTENT_SET") is None)
}

test "msig_filter_by_size" {
  let col = @src.msig_hallmark()
  let big = @src.msig_filter_by_size(col, 12, 1000)
  for s in big.sets {
    assert_true(s.size() >= 12)
  }
  let all = @src.msig_filter_by_size(col, 1, 1000)
  assert_eq(all.size(), col.size())
}

test "msig_to_gmt_roundtrip" {
  let col = @src.msig_hallmark()
  let gmt = @src.msig_to_gmt(col)
  assert_true(gmt.contains("HALLMARK_TNFA_SIGNALING_VIA_NFKB"))
  let parsed = @src.msig_from_gmt(gmt, @src.MsigHallmark)
  assert_eq(parsed.size(), col.size())
  match parsed.sets[0] {
    s => assert_true(s.name.has_prefix("HALLMARK_"))
  }
}

test "msig_id_table" {
  let pairs = [("TP53", "7157"), ("BRCA1", "672")]
  let table = @src.msig_id_table(pairs)
  match table.get("TP53") {
    Some(v) => assert_eq(v, "7157")
    None => fail("expected TP53 mapping")
  }
  match table.get("BRCA1") {
    Some(v) => assert_eq(v, "672")
    None => fail("expected BRCA1 mapping")
  }
  assert_true(table.get("EGFR") is None)
}

test "msig_map_ids" {
  let gs = @src.GmtGeneSet::new("test", "test set", ["TP53", "BRCA1", "EGFR"])
  let table = @src.msig_id_table([("TP53", "7157"), ("BRCA1", "672")])
  let mapped = @src.msig_map_ids(gs, table)
  assert_eq(mapped.gene_ids[0], "7157")
  assert_eq(mapped.gene_ids[1], "672")
  assert_eq(mapped.gene_ids[2], "EGFR")
}

test "msig_map_collection" {
  let col = @src.msig_hallmark()
  let table = @src.msig_id_table([("TP53", "7157")])
  let mapped = @src.msig_map_collection(col, table)
  assert_eq(mapped.size(), col.size())
  let mut found = false
  for s in mapped.sets {
    for g in s.gene_ids {
      if g == "7157" {
        found = true
      }
    }
  }
  assert_true(found)
}

test "msig_gene_frequency" {
  let col = @src.msig_hallmark()
  let freq = @src.msig_gene_frequency(col)
  assert_true(freq.size() > 0)
  let tp53 = freq.get("TP53")
  assert_true(tp53 is Some(_))
}

test "msig_shared_genes" {
  let a = @src.GmtGeneSet::new("a", "a", ["TP53", "BRCA1", "EGFR"])
  let b = @src.GmtGeneSet::new("b", "b", ["TP53", "MYC", "KRAS"])
  let shared = @src.msig_shared_genes(a, b)
  assert_eq(shared.length(), 1)
  assert_eq(shared[0], "TP53")
}

test "msig_jaccard" {
  let a = @src.GmtGeneSet::new("a", "a", ["TP53", "BRCA1", "EGFR"])
  let b = @src.GmtGeneSet::new("b", "b", ["TP53", "MYC", "KRAS"])
  let j = @src.msig_jaccard(a, b)
  assert_true(j > 0.0 && j < 1.0)
  let same = @src.msig_jaccard(a, a)
  assert_true((same - 1.0).abs() < 1.0e-9)
}

test "msig_summary" {
  let col = @src.msig_hallmark()
  let s = @src.msig_summary(col)
  assert_true(s.contains("MsigCollection[H]"))
  assert_true(s.contains("HALLMARK_TNFA_SIGNALING_VIA_NFKB"))
}
