// ===========================================================================
// variant_benchmark_test.mbt — truth vs query 变异基准评估测试
// ===========================================================================

///|
fn variant_bench_test_truth() -> String {
  "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\nchr1\t100\t.\tA\tT\t50\tPASS\t.\tGT\t0/1\nchr1\t200\t.\tAT\tA\t40\tPASS\t.\tGT\t1/1\nchr2\t300\t.\tG\tC\t30\tPASS\t.\tGT\t0/0\n"
}

///|
fn variant_bench_test_query() -> String {
  "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\nchr1\t100\t.\tA\tT\t60\tPASS\t.\tGT\t0/1\nchr1\t200\t.\tAT\tA\t20\tPASS\t.\tGT\t0/1\nchr1\t400\t.\tC\tG\t10\tPASS\t.\tGT\t0/1\n"
}

///|
test "variant benchmark normalizes ref alt by trimming prefix suffix" {
  let vcf = @src.parse_vcf(
    "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\nchr1\t100\t.\tAGC\tATG\t50\tPASS\t.\n",
  )
  let variant = @src.load_variants(vcf)[0]
  assert_eq(variant.pos, 101)
  assert_eq(variant.ref_allele, "GC")
  assert_eq(variant.alt_allele, "TG")
  assert_true(variant.kind == @src.Mnv)
}

///|
test "variant benchmark parses genotype dosage" {
  let truth = @src.parse_vcf(variant_bench_test_truth())
  let variants = @src.load_variants(truth, sample="S1")
  assert_eq(variants.length(), 3)
  assert_eq(variants[0].gt.unwrap(), "0/1")
  assert_eq(variants[0].alt_dosage.unwrap(), 1)
  assert_eq(variants[1].alt_dosage.unwrap(), 2)
  assert_eq(variants[2].alt_dosage.unwrap(), 0)
}

///|
test "variant benchmark compares tp fp fn and genotype mismatches" {
  let truth = @src.load_variants(
    @src.parse_vcf(variant_bench_test_truth()),
    sample="S1",
  )
  let query = @src.load_variants(
    @src.parse_vcf(variant_bench_test_query()),
    sample="S1",
  )
  let stats = @src.compare_variants(truth, query)
  assert_eq(stats.truth_total, 3)
  assert_eq(stats.query_total, 3)
  assert_eq(stats.tp, 2)
  assert_eq(stats.fp, 1)
  assert_eq(stats.fn_, 1)
  assert_eq(stats.genotype_mismatches, 1)
  assert_true(stats.precision > 0.66 && stats.precision < 0.67)
  assert_true(stats.recall > 0.66 && stats.recall < 0.67)
  assert_true(stats.f1 > 0.66 && stats.f1 < 0.67)
  assert_eq(stats.snp_truth, 2)
  assert_eq(stats.snp_query, 2)
  assert_eq(stats.snp_tp, 1)
  assert_eq(stats.indel_truth, 1)
  assert_eq(stats.indel_query, 1)
  assert_eq(stats.indel_tp, 1)
}

///|
test "variant benchmark pr curve respects quality thresholds" {
  let truth = @src.load_variants(
    @src.parse_vcf(variant_bench_test_truth()),
    sample="S1",
  )
  let query = @src.load_variants(
    @src.parse_vcf(variant_bench_test_query()),
    sample="S1",
  )
  let points = @src.benchmark_pr_curve(truth, query, thresholds=[
    0.0, 30.0, 55.0,
  ])
  assert_eq(points.length(), 3)
  assert_eq(points[0].tp, 2)
  assert_eq(points[1].tp, 1)
  assert_eq(points[1].fp, 0)
  assert_eq(points[2].tp, 1)
  assert_true(points[1].precision == 1.0)
}

///|
test "variant benchmark pass only keeps pass records" {
  let truth = @src.load_variants(
    @src.parse_vcf(variant_bench_test_truth()),
    sample="S1",
  )
  let query = @src.load_variants(
    @src.parse_vcf(variant_bench_test_query()),
    sample="S1",
  )
  let stats = @src.compare_variants(truth, query, pass_only=true)
  assert_eq(stats.fp, 1)

  let json = @src.benchmark_to_json(stats)
  assert_true(json.contains("precision"))
  assert_true(json.contains("\"tp\":2") || json.contains("\"tp\": 2"))
}
