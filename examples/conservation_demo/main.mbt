///|
/// conservation_demo — GERP/phastCons-lite 进化约束打分演示:
/// 合成 5 物种 MSA(保守-可变-保守三段) -> GERP-lite 逐列拒绝打分 ->
/// 保守后验 -> 两状态 Viterbi 保守区段 -> GenomicScores 容器闭环.

///|
fn csd_msa() -> Array[String] {
  [
    "ACGTACGT" + "TTTTGGGG" + "TGCA",
    "ACGTACGT" + "CCCCAAAA" + "TGCA",
    "ACGTACGT" + "TTTTGGGG" + "TGCA",
    "ACGTACGT" + "CCCCAAAA" + "TGCA",
    "ACGTACGT" + "ACGTACGT" + "TGCA",
  ]
}

///|
fn main raise {
  println(
    "1. GERP/phastCons-lite 进化约束打分演示(熵代理 + 两状态 HMM)",
  )
  let msa = csd_msa()
  println(
    "2. 合成 MSA: \{msa.length()} 物种 x \{msa[0].length()} 列(保守 0-7 / 可变 8-15 / 保守 16-19)",
  )
  // 3. GERP-lite
  let gerp = @src.conservation_gerp(msa, tree_length=2.0)
  println("3. GERP-lite 逐列拒绝打分(树长 2.0):")
  let mut row = "  "
  for i = 0; i < 20; i = i + 1 {
    let v = gerp[i]
    let tag = if v > 1.5 { "H" } else if v > 0.5 { "m" } else { "." }
    row = row + tag
  }
  println("\{row}  (H=高约束 m=中 .=低)")
  println("  列 0 = \{gerp[0]}, 列 8 = \{gerp[8]}, 列 16 = \{gerp[16]}")
  // 4. 后验 + Viterbi
  let post = @src.conservation_posterior(msa)
  let segs = @src.conservation_viterbi_segments(msa, posterior=Some(post))
  println("4. 保守后验与 Viterbi 区段:")
  for s in segs {
    println(
      "  保守区段 [\{s.start}, \{s.end_}) 平均后验=\{s.mean_posterior}",
    )
  }
  // 5. GenomicScores 闭环
  let gs = @src.GScores::new("GERP-lite", "demo")
  @src.conservation_to_gscores(gs, "chr1", gerp)
  gs.set_default_population("gerp_lite")
  let v = gs.score_at("chr1", 1)
  println("5. GenomicScores 容器闭环: chr1:1 分数 = \{v}")
  // 6. 导出
  let tsv = @src.conservation_tsv("chr1", gerp, post)
  let json = @src.conservation_json(gerp, post, segs)
  println("6. TSV \{tsv.length()} 字节, JSON \{json.length()} 字节")
  println(@src.conservation_summary_line(gerp, segs))
  println(
    "演示完成: MSA -> GERP-lite -> 后验 -> 保守区段 -> 打分容器闭环",
  )
}
