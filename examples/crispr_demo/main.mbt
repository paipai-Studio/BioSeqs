///|
/// CRISPR gRNA 设计演示: PAM 扫描 -> 打分 -> off-target -> JSON.
fn main {
  println("=== CRISPR gRNA 设计器 (SpCas9 NGG) ===")
  // 目标序列(含一处重复片段制造 off-target)
  let target = "ATGCGTACCTGAGGCTAAGCTTACGGATCCGGGATCTAGAATTCCGGAAGCTTGGATCCCTAGCGGATC"
  let reference = "GCTTAAGCTTACGGATCCGGGATCTAGAATTCCGGGATC" + target
  println(
    "目标 \{target.length()}bp / 参考基因组 \{reference.length()}bp",
  )

  // 1. 扫描全部候选
  println("\n1. 候选扫描")
  let all = @src.crispr_scan_guides(target)
  println("  共 \{all.length()} 条候选(含正/反链)")

  // 2. 完整设计(打分+off-target+排序)
  println("\n2. Top 5 设计结果")
  let guides = @src.crispr_design_guides(target, reference~, top_n=5)
  println("  spacer              PAM  strand  on-target  spec     total")
  for g in guides {
    println(
      "  \{g.spacer}  \{g.pam}  \{g.strand}       \{g.on_target}  \{g.specificity}  \{g.total_score}",
    )
    for w in g.warnings {
      println("    [警告] \{w}")
    }
  }

  // 3. off-target 明细
  println("\n3. Top 1 候选 off-target 明细")
  if guides.length() > 0 {
    let top = guides[0]
    println(
      "  spacer=\{top.spacer} 共 \{top.off_targets.length()} 个脱靶位点:",
    )
    for o in top.off_targets {
      println(
        "    pos=\{o.pos} mismatches=\{o.mismatches} score=\{o.site_score}",
      )
    }
  }

  // 4. Cas12a 演示
  println("\n4. Cas12a(TTTV PAM) 候选")
  let cas12a = @src.crispr_scan_guides(
    "AAAATTTGAAACCCCGGGGTTTTAAAACCCCGGGGTTTTAAA",
    pam=@src.crispr_pam_cas12a(),
    spacer_len=23,
  )
  println("  共 \{cas12a.length()} 条候选")

  // 5. JSON 输出
  println("\n5. JSON 输出(parse_json 闭环)")
  let json = @src.crispr_guides_json(guides)
  let v = @src.parse_json(json)
  println("  \{json.length()} 字符, 解析: \{v is Some(_)}")

  println(
    "\n完成: PAM 扫描 -> Doench 简化打分 -> FM-index 脱靶评估 -> JSON 全链路本地完成.",
  )
}
