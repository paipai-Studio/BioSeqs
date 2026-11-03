///|
/// methylation_caller_demo — Bismark 风格甲基化调用器演示:
/// 合成参考 + 正/反链转化读段 -> 逐位点 meth/unmeth 计数 ->
/// Bismark TSV -> bsseq_parse_bismark 闭环 -> 汇总/JSON.

///|
/// 20bp 参考(含 CpG/CHG/CHH 位点).
fn mcd_reference() -> String {
  "ACGTCAGCAATCGTTCAACA"
}

///|
/// 正链读段(部分甲基化).
fn mcd_forward_read() -> String {
  "ACGTTAGCAATTGTTTAACA"
}

///|
/// 反链读段(正链方向存储).
fn mcd_reverse_read() -> String {
  "ACGTCAACAATCGTTCAACA"
}

///|
fn main {
  let reference = mcd_reference()
  println("1. 合成参考: 20bp (CpG/CHG/CHH 位点混布)")
  // 2. 合成 SAM: 4 条正链 + 3 条反链
  let sb = StringBuilder()
  sb.write_string("@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:1000\n")
  for i = 0; i < 4; i = i + 1 {
    sb.write_string(
      "f\{i}\t0\tchr1\t1\t60\t20M\t*\t0\t0\t\{mcd_forward_read()}\t*\n",
    )
  }
  for i = 0; i < 3; i = i + 1 {
    sb.write_string(
      "r\{i}\t16\tchr1\t1\t60\t20M\t*\t0\t0\t\{mcd_reverse_read()}\t*\n",
    )
  }
  let sam = sb.to_string()
  println("2. 合成 SAM: 4 条正链(C->T 转化) + 3 条反链(G->A 转化)")
  // 3. 提取
  let calls = @src.mc_extract("chr1", reference, sam)
  println("3. 逐位点调用 (\{calls.length()} 个位点):")
  println("  pos strand ctx meth unmeth rate")
  for c in calls {
    let total = c.meth + c.unmeth
    let rate = if total > 0 {
      c.meth.to_double() / total.to_double()
    } else {
      0.0
    }
    println("  \{c.pos} \{c.strand} \{c.context} \{c.meth} \{c.unmeth} \{rate}")
  }
  // 4. Bismark TSV + bsseq 闭环
  let tsv = @src.mc_bismark_tsv(calls)
  let bsdata = @src.bsseq_parse_bismark(tsv)
  let n_bsrecords = bsdata.records.length()
  println("4. Bismark TSV -> bsseq 闭环: \{n_bsrecords} 条 BSRecord")
  // 5. 汇总
  let summary = @src.mc_summary(sam, calls)
  println(
    "5. CpG=\{summary.n_cpg_sites}, CHG=\{summary.n_chg_sites}, CHH=\{summary.n_chh_sites}, CpG 率=\{summary.cpg_methylation_rate}",
  )
  println(@src.mc_summary_line(summary))
  // 6. JSON
  let json = @src.mc_json(calls)
  println("6. JSON 导出: \{json.length()} 字节")
  println(
    "演示完成: 转化读段 -> 逐位点计数 -> Bismark TSV -> bsseq 统计链闭环",
  )
}
