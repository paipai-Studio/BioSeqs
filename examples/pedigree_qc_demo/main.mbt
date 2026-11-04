///|
/// 家系与基因型一致性 QC 演示:
/// .fam 家谱 -> trio 提取 -> 孟德尔错误 -> 样本一致性 -> ROH -> JSON.

///|
fn pq_demo_fam() -> String {
  "F1 F 0 0 1 1\n" + "F1 M 0 0 2 1\n" + "F1 C F M 1 2\n"
}

///|
fn pq_demo_vcf() -> String {
  "##fileformat=VCFv4.2\n" +
  "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tF\tM\tC\n" +
  "chr1\t10\t.\tA\tG\t.\t.\t.\tGT\t0/1\t0/1\t0/1\n" +
  "chr1\t20\t.\tC\tT\t.\t.\t.\tGT\t0/0\t1/1\t0/1\n" +
  "chr1\t30\t.\tG\tA\t.\t.\t.\tGT\t0/1\t0/0\t0/0\n" +
  "chr1\t40\t.\tT\tC\t.\t.\t.\tGT\t0/0\t0/0\t0/1\n" +
  "chr1\t50\t.\tA\tG\t.\t.\t.\tGT\t0/0\t0/0\t0/0\n"
}

///|
fn pq_demo_roh_vcf() -> String {
  "##fileformat=VCFv4.2\n" +
  "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS\n" +
  "chr1\t100\t.\tA\tG\t.\t.\t.\tGT\t0/0\n" +
  "chr1\t200\t.\tC\tT\t.\t.\t.\tGT\t0/0\n" +
  "chr1\t300\t.\tG\tA\t.\t.\t.\tGT\t1/1\n" +
  "chr1\t400\t.\tT\tC\t.\t.\t.\tGT\t0/0\n" +
  "chr1\t500\t.\tA\tG\t.\t.\t.\tGT\t0/1\n" +
  "chr1\t600\t.\tC\tT\t.\t.\t.\tGT\t0/0\n"
}

///|
fn main {
  println("=== 家系与基因型一致性 QC ===")
  let pedigree = @src.pedigree_parse(pq_demo_fam())
  let vcf = @src.parse_vcf(pq_demo_vcf())
  println(
    "个体数: \{pedigree.individuals.length()}, VCF 位点: \{vcf.records.length()}",
  )

  println("\n1. 核心家系(trio)")
  let trios = @src.pedigree_trios(pedigree)
  for trio in trios {
    println("  child=\{trio.child} father=\{trio.father} mother=\{trio.mother}")
  }

  println("\n2. 孟德尔错误检测")
  for trio in trios {
    let qc = @src.pedigree_mendelian_check(vcf, trio)
    println(
      "  trio \{trio.child}: 检查 \{qc.sites_checked} 位点, 错误 \{qc.errors.length()}",
    )
    for err in qc.errors {
      println(
        "    ERROR \{err.chrom}:\{err.pos} child=\{err.child}(\{err.child_gt}) 父=\{err.father_gt} 母=\{err.mother_gt}",
      )
    }
  }
  println("  摘要: \{@src.pedigree_qc_summary(vcf, pedigree)}")

  println("\n3. 样本间基因型一致性")
  let fm = @src.pedigree_concordance(vcf, "F", "M")
  println(
    "  F vs M: 比较 \{fm.sites_compared}, 一致 \{fm.concordant}, 不一致 \{fm.discordant}, concordance=\{fm.concordance}",
  )
  let cf = @src.pedigree_concordance(vcf, "C", "F")
  println("  C vs F: concordance=\{cf.concordance}")

  println("\n4. 纯合区段(ROH)检测")
  let roh_vcf = @src.parse_vcf(pq_demo_roh_vcf())
  let roh = @src.pedigree_roh(roh_vcf, sample="S", min_markers=3, max_gap=1000)
  println(
    "  样本 S: 标记 \{roh.total_markers}, ROH 区段 \{roh.segments.length()}",
  )
  for seg in roh.segments {
    println(
      "    \{seg.chrom}:\{seg.start}-\{seg.end} (\{seg.n_markers} 标记)",
    )
  }

  println("\n5. JSON 汇总")
  println(@src.pedigree_qc_report_json(pedigree, vcf))

  println(
    "\n演示完成: 家系结构 + 孟德尔错误 + 一致性 + ROH 闭环",
  )
}
