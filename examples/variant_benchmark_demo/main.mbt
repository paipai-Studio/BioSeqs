///|
/// 变异检出基准评估演示: truth vs query.
fn main raise {
  println("=== 变异检出基准评估 (vcfeval 风格, 本地完成) ===")

  let truth_text = "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\nchr1\t100\t.\tA\tT\t50\tPASS\t.\tGT\t0/1\nchr1\t200\t.\tAT\tA\t40\tPASS\t.\tGT\t1/1\nchr2\t300\t.\tG\tC\t30\tPASS\t.\tGT\t0/0\n"
  let query_text = "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\nchr1\t100\t.\tA\tT\t60\tPASS\t.\tGT\t0/1\nchr1\t200\t.\tAT\tA\t20\tPASS\t.\tGT\t0/1\nchr1\t400\t.\tC\tG\t10\tPASS\t.\tGT\t0/1\n"

  let truth = @src.load_variants(@src.parse_vcf(truth_text), sample="S1")
  let query = @src.load_variants(@src.parse_vcf(query_text), sample="S1")
  let stats = @src.compare_variants(truth, query)

  println("\ntruth = \{stats.truth_total}, query = \{stats.query_total}")
  println("TP = \{stats.tp}, FP = \{stats.fp}, FN = \{stats.fn_}")
  println(
    "Precision = \{stats.precision}, Recall = \{stats.recall}, F1 = \{stats.f1}",
  )
  println("基因型不一致 = \{stats.genotype_mismatches}")

  println("\n按变异类型:")
  println(
    "  SNP: truth=\{stats.snp_truth}, query=\{stats.snp_query}, TP=\{stats.snp_tp}",
  )
  println(
    "  Indel: truth=\{stats.indel_truth}, query=\{stats.indel_query}, TP=\{stats.indel_tp}",
  )

  println("\nPR 曲线(按 query QUAL 阈值):")
  let points = @src.benchmark_pr_curve(truth, query)
  for p in points {
    println(
      "  QUAL >= \{p.threshold}: TP=\{p.tp} FP=\{p.fp} FN=\{p.fn_} P=\{p.precision} R=\{p.recall}",
    )
  }

  println("\nJSON 输出:")
  println(@src.benchmark_to_json(stats))

  println("\n完成: 基准评估全程本地完成, 数据不出本机.")
}
