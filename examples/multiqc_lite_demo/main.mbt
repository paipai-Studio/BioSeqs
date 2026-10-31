///|
/// MultiQC-lite 聚合报告演示: 多样本 x 多模块 -> HTML/JSON/Markdown.
fn main raise {
  println("=== MultiQC-lite 多样本聚合报告 ===")

  // 模拟两个样本的 FASTQ 与定量/变异结果
  let fastq_a = "@a1\nACGTACGTACGT\n+\nIIIIIIIIIIII\n@a2\nGGCCTTAAGGCC\n+\nJJJJJJJJJJJJ\n"
  let fastq_b = "@b1\nTTAACCGGTTAA\n+\nIIIIIIIIIIII\n@b2\nCCCCGGGGTTTT\n+\nHHHHHHHHHHHH\n@b3\nACGTACGTNNNN\n+\nFFFFFFFFFFFF\n"
  let reference = "ATCGGATTCAGCCTAGGATCAAGTTCGGATACCTAGGCATTGCAAGGTCCTAGAAGT"
  let reads : Array[@src.GatkRead] = []
  for i = 0; i < 10; i = i + 1 {
    let sliced = reference[10:30].to_owned()
    let mutated = if i < 5 {
      sliced[0:5].to_owned() + "C" + sliced[6:].to_owned()
    } else {
      sliced
    }
    reads.push({ name: "r" + i.to_string(), seq: mutated, pos: 10, })
  }
  let variants = @src.gatk_call_variants("chr1", reference, reads)
  let quant = @src.kallisto_run(
    ["t1", "t2"],
    [
      "AAAAACCCCCGGGGGTTTTTATATGCGCGCTATAT", "TTTTTGGGGGCCCCCAAAAAATATGCGCGCTATAT",
    ],
    [
      "AAAAACCCCCGGGGGTTTTTATAT", "TTTTTGGGGGCCCCCAAAAA", "AAAAACCCCCGGGGGTTTTTATAT",
    ],
    k=9,
  )

  // 聚合
  let fqqc_a = @src.multiqc_collect_fqqc(fastq_a, "sampleA")
  let fqqc_b = @src.multiqc_collect_fqqc(fastq_b, "sampleB")
  let report = @src.multiqc_report(title="Demo Run")
  let _ = report.add_module(fqqc_a)
  let _ = report.add_module(fqqc_b)
  let _ = report.add_module(@src.multiqc_collect_kallisto(quant, "sampleA"))
  let _ = report.add_module(@src.multiqc_collect_gatk(variants, "sampleB"))
  // 合并同名模块(多样本聚合, 对齐 MultiQC 语义)
  let merged = report.merge_modules()
  println(merged.summary())

  // Markdown 输出
  println("\n1. Markdown 表格")
  println(@src.multiqc_write_markdown(merged))

  // JSON 输出与 roundtrip
  println("2. JSON 输出(roundtrip)")
  let json = @src.multiqc_write_json(merged)
  println("  " + json.substring(start=0, end=120) + "...")
  let parsed = @src.multiqc_parse_json(json).unwrap()
  println("  reparsed: " + parsed.summary())

  // HTML 仪表盘
  println("\n3. HTML 仪表盘(CSS 内嵌, 无外部依赖)")
  let html = @src.multiqc_write_html(merged)
  println("  长度: \{html.length()} 字符")
  println(
    "  包含着色单元格: \{html.contains("class=\"good\"") || html.contains("class=\"bad\"")}",
  )
  println("\n完成: 多样本多模块聚合报告全链路本地完成.")
}
