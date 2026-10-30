// ===========================================================================
// kraken_lite_test.mbt — 宏基因组 k-mer LCA 分类与丰度定量测试.
// ===========================================================================

///|
fn kraken_test_taxa() -> Array[@src.KrakenTaxon] {
  [
    @src.kraken_taxon_new("1", "root", "root", None),
    @src.kraken_taxon_new("10", "Alpha", "genus", Some("1")),
    @src.kraken_taxon_new("11", "Alpha one", "species", Some("10")),
    @src.kraken_taxon_new("12", "Alpha two", "species", Some("10")),
    @src.kraken_taxon_new("20", "Beta", "genus", Some("1")),
    @src.kraken_taxon_new("21", "Beta one", "species", Some("20")),
  ]
}

///|
fn kraken_test_references() -> Array[@src.KrakenReference] {
  [
    @src.kraken_reference_new("11", "ACGTAACCGGTT"),
    @src.kraken_reference_new("12", "ACGTAACCGGAA"),
    @src.kraken_reference_new("21", "TTTTGGGGCCCCAAAA"),
  ]
}

///|
fn kraken_test_database() -> @src.KrakenDatabase {
  @src.kraken_build_database(5, kraken_test_taxa(), kraken_test_references()) catch {
    KrakenLiteError(message) => abort("database build failed: " + message)
  }
}

///|
fn kraken_test_reads() -> Map[String, String] {
  Map([
    ("r1", "ACGTAACCGGTT"),
    ("r2", "ACGTAACCGGTT"),
    ("r3", "ACGTAACCGGTT"),
    ("r4", "ACGTAACCGG"),
    ("r5", "TTTTGGGGCCCCAAAA"),
    ("r6", "TTTTGGGGCCCCAAAA"),
  ])
}

///|
test "kraken_build_indexes_kmer_lca" {
  let db = kraken_test_database()
  assert_eq(db.k, 5)
  assert_eq(db.root, "1")
  // 两个物种共有的 k-mer -> 属级 LCA.
  assert_eq(@src.kraken_lookup_kmer(db, "ACGTA"), Some("10"))
  assert_eq(@src.kraken_lookup_kmer(db, "CGTAA"), Some("10"))
  // 物种特异的 k-mer -> 种级.
  assert_eq(@src.kraken_lookup_kmer(db, "CCGGT"), Some("11"))
  assert_eq(@src.kraken_lookup_kmer(db, "CGGAA"), Some("12"))
  assert_eq(@src.kraken_lookup_kmer(db, "TTTTG"), Some("21"))
  assert_eq(@src.kraken_lookup_kmer(db, "AAAAA"), None)
}

///|
test "kraken_build_computes_subtree_support" {
  let db = kraken_test_database()
  assert_eq(@src.kraken_kmer_support(db, "11"), 2)
  assert_eq(@src.kraken_kmer_support(db, "12"), 2)
  assert_eq(@src.kraken_kmer_support(db, "21"), 12)
  assert_eq(@src.kraken_kmer_support(db, "10"), 10)
  assert_eq(@src.kraken_kmer_support(db, "20"), 12)
  assert_eq(@src.kraken_kmer_support(db, "1"), 22)
}

///|
test "kraken_classify_unique_read" {
  let db = kraken_test_database()
  let result = @src.kraken_classify(db, "read11", "ACGTAACCGGTT")
  assert_eq(result.taxid, Some("11"))
  assert_eq(result.rank, "species")
  assert_eq(result.total_kmers, 8)
  assert_eq(result.hits, 8)
  assert_true(result.confidence > 0.999)
}

///|
test "kraken_classify_shared_read_stops_at_lca" {
  let db = kraken_test_database()
  let result = @src.kraken_classify(db, "shared", "ACGTAACCGG")
  assert_eq(result.taxid, Some("10"))
  assert_eq(result.rank, "genus")
  assert_true(result.confidence > 0.999)
}

///|
test "kraken_classify_ignores_invalid_or_short_reads" {
  let db = kraken_test_database()
  let ambiguous = @src.kraken_classify(db, "n", "NNNNNN")
  assert_eq(ambiguous.taxid, None)
  assert_eq(ambiguous.rank, "unclassified")
  let short = @src.kraken_classify(db, "short", "ACGT")
  assert_eq(short.taxid, None)
}

///|
test "kraken_classify_many_is_sorted_and_consistent" {
  let db = kraken_test_database()
  let results = @src.kraken_classify_many(db, kraken_test_reads())
  assert_eq(results.length(), 6)
  assert_eq(results[0].name, "r1")
  assert_eq(results[2].taxid, Some("11"))
  assert_eq(results[3].taxid, Some("10"))
  assert_eq(results[5].taxid, Some("21"))
}

///|
test "kraken_report_rolls_up_counts" {
  let db = kraken_test_database()
  let results = @src.kraken_classify_many(db, kraken_test_reads())
  let species = @src.kraken_report(db, results, "species")
  assert_eq(species.length(), 2)
  assert_eq(species[0].taxid, "11")
  assert_eq(species[0].count, 3.0)
  assert_eq(species[1].taxid, "21")
  assert_eq(species[1].count, 2.0)
  assert_true((species[0].relative_abundance - 0.5).abs() < 1.0e-9)

  let genera = @src.kraken_report(db, results, "genus")
  assert_eq(genera.length(), 2)
  assert_eq(genera[0].taxid, "10")
  assert_eq(genera[0].count, 4.0)
  assert_eq(genera[1].taxid, "20")
  assert_eq(genera[1].count, 2.0)
}

///|
test "kraken_bracken_redistributes_genus_reads" {
  let db = kraken_test_database()
  let results = @src.kraken_classify_many(db, kraken_test_reads())
  let species = @src.kraken_bracken(db, results, "species")
  assert_eq(species.length(), 3)
  assert_eq(species[0].taxid, "11")
  assert_true((species[0].count - 3.5).abs() < 1.0e-9)
  assert_eq(species[1].taxid, "21")
  assert_true((species[1].count - 2.0).abs() < 1.0e-9)
  assert_eq(species[2].taxid, "12")
  assert_true((species[2].count - 0.5).abs() < 1.0e-9)
  // 重分摊后相对丰度总和为 1.
  let mut total = 0.0
  for entry in species {
    total = total + entry.relative_abundance
  }
  assert_true((total - 1.0).abs() < 1.0e-9)
}

///|
test "kraken_json_and_tsv_outputs" {
  let db = kraken_test_database()
  let results = @src.kraken_classify_many(db, kraken_test_reads())
  let species = @src.kraken_report(db, results, "species")
  let json = @src.kraken_abundance_json(db, species)
  let parsed = @src.parse_json(json)
  assert_true(parsed is Some(_))
  assert_true(json.contains("\"taxid\""))
  assert_true(json.contains("\"relative_abundance\""))

  let tsv = @src.kraken_abundance_tsv(species)
  let lines = tsv.split("\n").to_array()
  assert_eq(lines[0].to_owned(), "taxid\tname\trank\tcount\trelative_abundance")
  assert_eq(lines.length(), 4)
}

///|
test "kraken_build_rejects_invalid_taxonomy" {
  let bad_parent = [
    @src.kraken_taxon_new("1", "root", "root", None),
    @src.kraken_taxon_new("2", "orphan", "genus", Some("999")),
  ]
  let raised = try {
    ignore(@src.kraken_build_database(5, bad_parent, []))
    false
  } catch {
    @src.KrakenLiteError(_) => true
  }
  assert_true(raised)
}

// ===========================================================================
// H3: minimizer 建库 / 数据库 JSON / HTML 报告 / MultiQC 采集器
// ===========================================================================

///|
fn kraken_test_minimizer_database(w : Int32) -> @src.KrakenDatabase {
  @src.kraken_build_database_minimizer(
    5,
    w,
    kraken_test_taxa(),
    kraken_test_references(),
  ) catch {
    KrakenLiteError(message) => abort("minimizer build failed: " + message)
  }
}

///|
test "kraken_minimizer_database_is_compacter" {
  let exact = kraken_test_database()
  let mini = kraken_test_minimizer_database(4)
  assert_eq(mini.minimizer_w, 4)
  // minimizer 建库保留的索引项应严格少于精确建库.
  assert_true(mini.kmer_lca.length() < exact.kmer_lca.length())
}

///|
test "kraken_minimizer_classifies_unique_read" {
  let db = kraken_test_minimizer_database(4)
  let result = @src.kraken_classify(db, "read21", "TTTTGGGGCCCCAAAA")
  assert_eq(result.taxid, Some("21"))
  assert_eq(result.rank, "species")
}

///|
test "kraken_minimizer_window_must_be_positive" {
  let raised = try {
    ignore(
      @src.kraken_build_database_minimizer(
        5,
        0,
        kraken_test_taxa(),
        kraken_test_references(),
      ),
    )
    false
  } catch {
    @src.KrakenLiteError(_) => true
  }
  assert_true(raised)
}

///|
test "kraken_database_json_roundtrip" {
  let db = kraken_test_minimizer_database(4)
  let json = @src.kraken_database_to_json(db)
  assert_true(json.contains("\"minimizer_w\":4"))
  match @src.kraken_database_from_json(json) {
    Some(restored) => {
      assert_eq(restored.k, db.k)
      assert_eq(restored.root, db.root)
      assert_eq(restored.minimizer_w, db.minimizer_w)
      assert_eq(restored.taxa.length(), db.taxa.length())
      assert_eq(restored.kmer_lca.length(), db.kmer_lca.length())
      assert_eq(restored.subtree_kmers.length(), db.subtree_kmers.length())
      // 还原后分类结果一致.
      let before = @src.kraken_classify(db, "r", "TTTTGGGGCCCCAAAA")
      let after = @src.kraken_classify(restored, "r", "TTTTGGGGCCCCAAAA")
      assert_eq(before.taxid, after.taxid)
    }
    None => assert_true(false)
  }
}

///|
test "kraken_database_json_roundtrip_exact_mode" {
  let db = kraken_test_database()
  let json = @src.kraken_database_to_json(db)
  match @src.kraken_database_from_json(json) {
    Some(restored) => {
      assert_eq(restored.minimizer_w, 0)
      assert_eq(restored.kmer_lca.length(), db.kmer_lca.length())
      assert_eq(@src.kraken_lookup_kmer(restored, "ACGTA"), Some("10"))
    }
    None => assert_true(false)
  }
}

///|
test "kraken_database_from_json_rejects_garbage" {
  assert_true(@src.kraken_database_from_json("not json") is None)
  assert_true(@src.kraken_database_from_json("{}") is None)
}

///|
test "kraken_report_html_renders_svg_and_taxa" {
  let db = kraken_test_database()
  let results = @src.kraken_classify_many(db, kraken_test_reads())
  let species = @src.kraken_report(db, results, "species")
  let html = @src.kraken_report_html(species, title="Demo Report")
  assert_true(html.contains("<svg"))
  assert_true(html.contains("Demo Report"))
  assert_true(html.contains("Alpha one"))
  assert_true(html.contains("</html>"))
}

///|
test "multiqc_collect_kraken_reports_metrics" {
  let db = kraken_test_database()
  let results = @src.kraken_classify_many(db, kraken_test_reads())
  let species = @src.kraken_report(db, results, "species")
  let module_ = @src.multiqc_collect_kraken(results, species, "sample1")
  assert_eq(module_.name, "kraken-lite")
  let metrics = module_.samples.get("sample1").unwrap()
  assert_eq(metrics.get("reads").unwrap(), 6.0)
  assert_eq(metrics.get("classified").unwrap(), 6.0)
  assert_eq(metrics.get("unclassified").unwrap(), 0.0)
  assert_eq(metrics.get("taxa_detected").unwrap(), 2.0)
  assert_true(metrics.get("mean_confidence").unwrap() > 0.9)
}
