///|
/// 宏基因组物种丰度定量演示(Kraken/Bracken 风格):
/// 分类树 + 参考基因组 -> k-mer LCA 库 -> 读段分类 -> 丰度报告.
fn main {
  println("=== 宏基因组物种丰度定量 (Kraken/Bracken-lite) ===")

  // 分类树.
  let taxa : Array[@src.KrakenTaxon] = [
    @src.kraken_taxon_new("1", "root", "root", None),
    @src.kraken_taxon_new("10", "Alpha", "genus", Some("1")),
    @src.kraken_taxon_new("11", "Alpha one", "species", Some("10")),
    @src.kraken_taxon_new("12", "Alpha two", "species", Some("10")),
    @src.kraken_taxon_new("20", "Beta", "genus", Some("1")),
    @src.kraken_taxon_new("21", "Beta one", "species", Some("20")),
  ]

  // 参考基因组(k=5).
  let references : Array[@src.KrakenReference] = [
    @src.kraken_reference_new("11", "ACGTAACCGGTT"),
    @src.kraken_reference_new("12", "ACGTAACCGGAA"),
    @src.kraken_reference_new("21", "TTTTGGGGCCCCAAAA"),
  ]
  let db = @src.kraken_build_database(5, taxa, references) catch {
    KrakenLiteError(message) => {
      println("建库失败: " + message)
      return
    }
  }
  println("\nk = \{db.k}, 分类单元数 = \{db.taxa.length()}")

  // 模拟读段: Alpha one x3, 属级共享 x1, Beta one x2.
  let reads : Map[String, String] = Map([
    ("r1", "ACGTAACCGGTT"),
    ("r2", "ACGTAACCGGTT"),
    ("r3", "ACGTAACCGGTT"),
    ("r4", "ACGTAACCGG"),
    ("r5", "TTTTGGGGCCCCAAAA"),
    ("r6", "TTTTGGGGCCCCAAAA"),
  ])
  let results = @src.kraken_classify_many(db, reads)
  println("\n1. 逐读段分类")
  for result in results {
    let taxid = match result.taxid {
      Some(value) => value
      None => "-"
    }
    println(
      "  \{result.name}: taxid=\{taxid}, rank=\{result.rank}, 支持=\{result.hits}/\{result.total_kmers}, 置信度=\{result.confidence}",
    )
  }

  println("\n2. Kraken 风格种级报告(累计计数)")
  let species = @src.kraken_report(db, results, "species")
  println(@src.kraken_abundance_tsv(species))

  println(
    "3. Bracken 风格种级重分摊(属级读段按 k-mer 支持度下推)",
  )
  let bracken = @src.kraken_bracken(db, results, "species")
  println(@src.kraken_abundance_tsv(bracken))

  println("4. JSON 输出")
  println(@src.kraken_abundance_json(db, bracken))

  println("\n5. minimizer 建库(内存压缩)")
  let mini = @src.kraken_build_database_minimizer(5, 4, taxa, references) catch {
    KrakenLiteError(message) => {
      println("minimizer 建库失败: " + message)
      return
    }
  }
  println(
    "  精确库索引项 = \{db.kmer_lca.length()}, minimizer(w=4) 索引项 = \{mini.kmer_lca.length()}",
  )
  let mini_read = @src.kraken_classify(mini, "r5", "TTTTGGGGCCCCAAAA")
  let mini_taxid = match mini_read.taxid {
    Some(value) => value
    None => "-"
  }
  println(
    "  minimizer 分类 r5 -> taxid=\{mini_taxid}, rank=\{mini_read.rank}",
  )

  println("\n6. 数据库 JSON 序列化 / roundtrip")
  let db_json = @src.kraken_database_to_json(mini)
  println("  序列化长度 = \{db_json.length()} 字节")
  match @src.kraken_database_from_json(db_json) {
    Some(restored) =>
      println(
        "  还原: k=\{restored.k}, minimizer_w=\{restored.minimizer_w}, 索引项=\{restored.kmer_lca.length()}",
      )
    None => println("  roundtrip 失败")
  }

  println("\n7. HTML/SVG 报告")
  let html = @src.kraken_report_html(bracken, title="Kraken-lite Demo")
  println(
    "  HTML 长度 = \{html.length()} 字节, 含 <svg> = \{html.contains("<svg")}",
  )

  println("\n8. MultiQC 采集器")
  let module_ = @src.multiqc_collect_kraken(results, species, "demo_sample")
  let report = @src.multiqc_report(title="Metagenomics")
  ignore(report.add_module(module_))
  println("  " + report.summary())
}
