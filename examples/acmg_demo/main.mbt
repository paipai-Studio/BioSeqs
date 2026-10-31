///|
/// acmg_demo — ACMG 2015 自动化分级 + ClinVar 摄取演示:
/// 合成 ClinVar 库 + 6 类典型变异(LoF/低频错义/高频多态/框内缺失/
/// 同义/ClinVar 致病位点) -> 证据推导 -> 五级分类 -> TSV/JSON 报告.

///|
fn acd_clinvar_tsv() -> String {
  "#ClinVar simplified TSV\n" +
  "chrom\tpos\tref\talt\tclnsig\tstar\n" +
  "chr1\t101\tA\tG\tPathogenic\t3\n" +
  "chr1\t805\tC\tT\tBenign\t2\n"
}

///|
fn main {
  println("1. ACMG 2015 自动化分级演示(InterVar 风格)")
  let clinvar = @src.clinvar_parse(acd_clinvar_tsv())
  println("2. ClinVar 库摄取: \{clinvar.length()} 条记录")
  // 3. 构造 6 类典型变异
  let inputs : Array[@src.AcmgInput] = []
  // (a) LoF + 低频 + 预测有害 + 功能阳性 -> PVS1+PM2+PP3+PS3
  inputs.push(
    @src.acmg_make_input(
      "chr1",
      701,
      "C",
      "T",
      "stop_gained",
      popmax_af=0.0001,
      in_silico_score=Some(28.0),
      functional_score=Some(1.5),
    ),
  )
  // (b) LoF + 低频 + de novo -> PVS1+PM2+PM6
  inputs.push(
    @src.acmg_make_input(
      "chr1",
      702,
      "CA",
      "C",
      "frameshift_variant",
      popmax_af=0.0002,
      de_novo=true,
    ),
  )
  // (c) 高频错义 -> BA1
  inputs.push(
    @src.acmg_make_input(
      "chr1",
      703,
      "A",
      "G",
      "missense_variant",
      popmax_af=0.12,
    ),
  )
  // (d) 低频错义 missense-only 基因 + 预测有害 -> PM2+PP2+PP3
  inputs.push(
    @src.acmg_make_input(
      "chr1",
      704,
      "G",
      "A",
      "missense_variant",
      popmax_af=0.0005,
      gene_missense_only=true,
      in_silico_score=Some(25.0),
    ),
  )
  // (e) 同义变异 -> BP7
  inputs.push(@src.acmg_make_input("chr1", 705, "T", "C", "synonymous_variant"))
  // (f) ClinVar 已知致病位点(>=2星) -> PS1+PP5+PM2
  inputs.push(
    @src.acmg_make_input(
      "chr1",
      101,
      "A",
      "G",
      "missense_variant",
      popmax_af=0.0005,
    ),
  )
  println(
    "3. 输入 6 类典型变异(LoF x2 / 高频 / miss-only / 同义 / ClinVar 位点)",
  )
  // 4. 批量解读
  let results = @src.acmg_annotate(inputs, clinvar~)
  println("4. 逐条解读:")
  for r in results {
    let mut ev = ""
    for i = 0; i < r.evidence.length(); i = i + 1 {
      if i > 0 {
        ev = ev + "&"
      }
      ev = ev + r.evidence[i]
    }
    println(
      "  chr1:\{r.pos} \{r.ref_allele}>\{r.alt} [\{r.consequence}] 证据=\{ev}",
    )
    let cv_note = if r.clinvar_sig == "" {
      ""
    } else {
      " (ClinVar: " + r.clinvar_sig + ")"
    }
    println("    分级=\{r.classification}\{cv_note}")
    for reason in r.reasons {
      println("    理由: \{reason}")
    }
  }
  // 5. 汇总与导出
  println("5. \{results.length()} 条变异")
  println(@src.acmg_summary_line(results))
  let tsv = @src.acmg_report_tsv(results)
  let json = @src.acmg_json(results)
  println("6. TSV \{tsv.length()} 字节, JSON \{json.length()} 字节")
  println(
    "演示完成: ClinVar 摄取 -> 证据推导 -> ACMG 五级分类 -> 报告闭环",
  )
}
