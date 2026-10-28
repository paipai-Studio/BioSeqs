///|
/// MHC 结合预测演示: 新抗原候选筛查(锚定规则框架).
fn main raise {
  println("=== MHC-I: HLA-A*02:01 九肽筛查 ===")
  let a0201 = @src.mhc_allele("HLA-A*02:01").unwrap()
  // 含经典流感表位 GILGFVFTL 的蛋白片段
  let protein1 = "MDSKGSCWGPSGDWEMMKVLDKLSEGKQGILGFVFTLTVPSERGVSILDSIRSLYNTVATL"
  let hits1 = @src.mhc_predict(protein1, a0201)
  let binders = @src.mhc_binders(hits1)
  println(
    "  蛋白 " +
    protein1.length().to_string() +
    "aa, 窗口 " +
    hits1.length().to_string() +
    " 个, 弱结合以上 " +
    binders.length().to_string() +
    " 个",
  )
  println("  Top 5 (锚匹配/得分/分级):")
  for h in hits1 {
    if h.binding_category() != "none" {
      println(
        "    " +
        h.peptide +
        " @" +
        h.start.to_string() +
        "  " +
        h.anchor_matches.to_string() +
        "/" +
        h.anchor_total.to_string() +
        "  " +
        h.score.to_string() +
        "  " +
        h.binding_category(),
      )
      if h.is_weak() {
        break
      }
    }
  }
  println("=== MHC-I: 小鼠 H-2Kb 八肽(OVA 表位) ===")
  let kb = @src.mhc_allele("H-2Kb").unwrap()
  let ova = "MQIERHVDGRTQIINFEKLTEWTSSNVMEERKIHQDDRHGLDLVQSLDLMKKLQQIV"
  let hits_kb = @src.mhc_predict(ova, kb)
  let top_kb = hits_kb[0]
  println(
    "  最优: " +
    top_kb.peptide +
    " @" +
    top_kb.start.to_string() +
    "  锚 " +
    top_kb.anchor_matches.to_string() +
    "/" +
    top_kb.anchor_total.to_string() +
    "  " +
    top_kb.binding_category(),
  )
  println("=== MHC-II: HLA-DRB1*01:01 核心 9-mer 扫描 ===")
  let dr = @src.mhc_allele("HLA-DRB1*01:01").unwrap()
  let protein2 = "MKTFVAAWTLKSAAVDQAKFIILCASTVAKFISAAWTLKSAAVDQAKFVAAWTLKSAAVD"
  let hits2 = @src.mhc_predict(protein2, dr)
  let uniq = @src.mhc_unique_epitopes(@src.mhc_binders(hits2))
  println(
    "  15-mer 窗口 " +
    hits2.length().to_string() +
    " 个, 去重后表位 " +
    uniq.length().to_string() +
    " 个(核心 9-mer 身份)",
  )
  for h in uniq {
    let core = h.peptide[h.core_start():h.core_start() + 9].to_owned()
    println(
      "    " +
      core +
      " (窗口 @" +
      h.start.to_string() +
      ", 框架 +" +
      h.core_start.to_string() +
      ")  " +
      h.score.to_string(),
    )
  }
  println("=== 内置等位基因 ===")
  println("  " + @src.mhc_allele_names().join(", "))
  println("")
  println(
    "注: 锚定规则筛查级工具, 精度低于 netMHCpan 神经网络方法(权重未公开无法移植);",
  )
  println("    候选集应交由实验或成熟工具验证.")
}
