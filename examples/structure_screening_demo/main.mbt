///|
/// 蛋白结构筛查演示: 二级结构 -> 折叠类 -> 无序/卷曲螺旋 -> 全景报告.
fn main raise {
  println("=== 蛋白结构筛查 (screening-grade, 非深度学习预测) ===")

  // 1. 全 alpha 型(螺旋富集, 血红蛋白样)
  println("\n1. 全 alpha 型序列")
  let alpha_seq = "AEAALKREAAALKREAAALKREAAALKREAAALKREAAALKREAAALK"
  let alpha_report = @src.structure_screening(alpha_seq)
  println(@src.structure_screening_report_text(alpha_report))
  println("  SS: " + alpha_report.ss_string())

  // 2. 全 beta 型(免疫球蛋白样)
  println("\n2. 全 beta 型序列")
  let beta_seq = "VVIVTYYIVFVVIVTYYIVFVVIVTYYIVFVVIVTYYIVFVVIVTYYIVF"
  let beta_report = @src.structure_screening(beta_seq)
  println(@src.structure_screening_report_text(beta_report))
  println("  SS: " + beta_report.ss_string())

  // 3. 无序倾向序列
  println("\n3. 无序倾向序列")
  let coil_seq = "GSGNQGSGNQGSGNQGSGNQGSGNQGSGNQGSGNQGSGNQGSGNQ"
  let coil_report = @src.structure_screening(coil_seq)
  println(@src.structure_screening_report_text(coil_report))

  // 4. JSON 报告
  println("\n4. JSON 报告(parse_json 闭环)")
  let json = @src.structure_screening_report_json(alpha_report)
  let lines = json.split("\n").to_array()
  println("  长度: \{json.length()} 字符")
  let v = @src.parse_json(json)
  println("  parse_json 解析: \{v is Some(_)}")

  // 5. 与 multiqc-lite 集成
  println("\n5. MultiQC-lite 集成")
  let m = @src.multiqc_module("structure-screening")
  let _ = m.add_sample("protein_alpha", [
    ("helix_percent", alpha_report.helix_percent),
    ("sheet_percent", alpha_report.sheet_percent),
    ("coil_percent", alpha_report.coil_percent),
  ])
  let _ = m.add_sample("protein_beta", [
    ("helix_percent", beta_report.helix_percent),
    ("sheet_percent", beta_report.sheet_percent),
    ("coil_percent", beta_report.coil_percent),
  ])
  let report = @src.multiqc_report(title="Structure Screening")
  let _ = report.add_module(m)
  println(@src.multiqc_write_markdown(report))
  println(
    "完成: 二级结构 -> 折叠类 -> 全景报告 -> 聚合 全链路本地完成.",
  )
}
