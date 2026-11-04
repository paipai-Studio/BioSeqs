///|
/// 免疫组库分析演示: V(D)J 指认 -> 克隆型聚类 -> 多样性 -> AIRR TSV -> TCRdist.
fn main raise {
  println("=== 免疫组库分析 (TRB, 教学级 germline 库) ===")

  // 构建自定义 germline 库(生产可用 vdj_build_gene_db 导入 IMGT 完整库)
  let v : Array[@src.VdjGene] = [
    {
      name: "TRBV01",
      locus: "TRB",
      kind: @src.VGene,
      sequence: "AAAACTGGGGTTTACCTGT",
    },
    {
      name: "TRBV02",
      locus: "TRB",
      kind: @src.VGene,
      sequence: "CCCCATCGGGGAACGATGT",
    },
    {
      name: "TRBV03",
      locus: "TRB",
      kind: @src.VGene,
      sequence: "TTTACGGCCTTAAGGTGT",
    },
  ]
  let d : Array[@src.VdjGene] = [
    {
      name: "TRBD1",
      locus: "TRB",
      kind: @src.DGene,
      sequence: "GCGCGATTATCAGGCA",
    },
    {
      name: "TRBD2",
      locus: "TRB",
      kind: @src.DGene,
      sequence: "TTATCGGATTAAGCCT",
    },
  ]
  let j : Array[@src.VdjGene] = [
    {
      name: "TRBJ01",
      locus: "TRB",
      kind: @src.JGene,
      sequence: "TGGGGCCATTTGAAAGTGG",
    },
    {
      name: "TRBJ02",
      locus: "TRB",
      kind: @src.JGene,
      sequence: "TGGGGTTAAACCGTTAGG",
    },
  ]
  let db = @src.vdj_build_gene_db(v, d, j, k=6)

  // 模拟一条样本的 read 集(含克隆型扩张)
  let v1 = "AAAACTGGGGTTTACCTGT"
  let v2 = "CCCCATCGGGGAACGATGT"
  let v3 = "TTTACGGCCTTAAGGTGT"
  let j1 = "TGGGGCCATTTGAAAGTGG"
  let j2 = "TGGGGTTAAACCGTTAGG"
  let reads : Array[(String, String)] = [
    (v1 + "ACG" + j1, "read01"),
    (v1 + "ACG" + j1, "read02"),
    (v1 + "ACG" + j1, "read03"),
    (v1 + "GCA" + j1, "read04"),
    (v2 + "TTA" + j2, "read05"),
    (v2 + "TTA" + j2, "read06"),
    (v3 + "GCC" + j1, "read07"),
    (v1 + "ACGT" + j1, "read08"),
  ]

  println("\n1. V(D)J 指认 (8 条 read)")
  let chains = @src.vdj_assign_all(reads, db)
  for ch in chains {
    println(
      "  " +
      ch.sequence_id +
      ": " +
      ch.v_call +
      "/" +
      ch.d_call +
      "/" +
      ch.j_call +
      "  CDR3=" +
      ch.cdr3_aa +
      " (" +
      ch.cdr3_nt.length().to_string() +
      "nt)  productive=" +
      ch.productive.to_string(),
    )
  }

  println("\n2. 克隆型聚类 (top)")
  let clonotypes = @src.vdj_clonotypes(chains)
  for c in clonotypes {
    println(
      "  #" +
      c.clonotype_id.to_string() +
      " " +
      c.cdr3_aa +
      " " +
      c.v_family +
      "/" +
      c.j_family +
      "  count=" +
      c.count.to_string() +
      " freq=" +
      (c.frequency * 100.0).to_string() +
      "%",
    )
  }

  println("\n3. 多样性指标")
  let div = @src.vdj_diversity(chains)
  println("  total chains:      " + div.total_chains.to_string())
  println("  total clonotypes:  " + div.total_clonotypes.to_string())
  println("  shannon:           " + div.shannon.to_string())
  println("  simpson:           " + div.simpson.to_string())
  println("  gini:              " + div.gini.to_string())
  println("  clonality:         " + div.clonality.to_string())
  println("  d50:               " + div.d50.to_string())

  println("\n4. AIRR TSV 输出")
  let records = @src.airr_from_vdj_chains(chains, locus="TRB")
  let tsv = @src.airr_write(records)
  let lines = tsv.split("\n").to_array()
  println("  [header] " + lines[0].to_owned())
  for li = 1; li < lines.length() && li <= 4; li = li + 1 {
    println("  [row]    " + lines[li].to_owned())
  }
  // roundtrip 验证
  let reparsed = @src.airr_parse(tsv)
  println("  roundtrip: " + reparsed.length().to_string() + " records reparsed")

  println("\n5. TCRdist-lite 克隆型间距离矩阵 (top 3)")
  let top_aas : Array[String] = []
  for i = 0; i < 3 && i < clonotypes.length(); i = i + 1 {
    top_aas.push(clonotypes[i].cdr3_aa)
  }
  let matrix = @src.vdj_tcrdist_matrix(top_aas)
  println("        " + top_aas.join("   "))
  for i = 0; i < matrix.length(); i = i + 1 {
    let cells : Array[String] = []
    for x in matrix[i] {
      cells.push(x.to_string())
    }
    println("  " + top_aas[i] + "  " + cells.join("  "))
  }
  println(
    "\n完成: 指认 -> 克隆型 -> 多样性 -> AIRR -> TCRdist 全链路本地完成.",
  )
}
