///|
/// P3: pysam_vcf_filter 单元测试
let filter_vcf_content : String = "##fileformat=VCFv4.2\n" +
  "##contig=<ID=chr1,length=1000>\n" +
  "##INFO=<ID=DP,Number=1,Type=Integer,Description=\"depth\">\n" +
  "##INFO=<ID=CSQ,Number=A,Type=String,Description=\"consequence\">\n" +
  "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n" +
  "chr1\t100\trs1\tA\tT\t60.0\tPASS\tDP=10;CSQ=missense\n" +
  "chr1\t200\trs2\tC\tG\t5.0\tLowQual\tDP=2\n" +
  "chr1\t300\t.\tG\tA\t30.0\t.\tDP=8;CSQ=synonymous\n" +
  "chr1\t400\trs4\tT\tC\t.\tPASS\t.\n" +
  "chr1\t500\t.\tA\tG\t40.0\tFAIL\tDP=15\n"

///|
test "filter_qual_min" {
  let vcf = @src.parse_vcf(filter_vcf_content)
  let result = @src.vcf_filter_qual_min(vcf, 10.0)
  // rs1(60), rs3(30), rs5(40) 通过；rs2(5) 和 rs4(无 qual) 被过滤
  assert_eq(result.len(), 3)
}

///|
test "filter_qual_range" {
  let vcf = @src.parse_vcf(filter_vcf_content)
  let cfg = @src.vcf_filter_qual(20.0, max=50.0)
  let result = @src.vcf_filter(vcf, cfg)
  // 只有 rs3(30.0) 和 rs5(40.0) 在 [20, 50] 内
  assert_eq(result.len(), 2)
}

///|
test "filter_pass_only" {
  let vcf = @src.parse_vcf(filter_vcf_content)
  let result = @src.vcf_filter_pass_only(vcf)
  // rs1 PASS, rs2 LowQual, rs3 空 filter [], rs4 PASS, rs5 FAIL
  // pass_only 逻辑：length == 1 且值为 PASS 或 "."
  // rs3 filter 是空数组 length=0 → 不满足 length==1 → 被过滤
  assert_eq(result.len(), 2)
}

///|
test "filter_info_required" {
  let vcf = @src.parse_vcf(filter_vcf_content)
  let cfg = @src.VcfFilterConfig::{
    qual_min: None,
    qual_max: None,
    pass_only: false,
    info_required_keys: ["CSQ"],
    min_alt_count: None,
    max_alt_count: None,
    id_known: None,
  }
  let result = @src.vcf_filter(vcf, cfg)
  // 只有 rs1 和 rs3 有 CSQ
  assert_eq(result.len(), 2)
}

///|
test "filter_record_passes_qual" {
  let vcf = @src.parse_vcf(filter_vcf_content)
  let cfg = @src.vcf_filter_qual(10.0)
  let r1 = vcf.records[0] // qual=60
  let r2 = vcf.records[1] // qual=5
  let r4 = vcf.records[3] // qual=None
  assert_true(@src.vcf_record_passes(r1, cfg))
  assert_true(!@src.vcf_record_passes(r2, cfg))
  assert_true(!@src.vcf_record_passes(r4, cfg))
}

///|
test "filter_pass_rate_summary" {
  let vcf = @src.parse_vcf(filter_vcf_content)
  let cfg = @src.vcf_filter_qual(10.0)
  let (passes, fails) = @src.vcf_filter_summary(vcf, cfg)
  assert_eq(passes, 3)
  assert_eq(fails, 2)
  assert_true(
    @moonbit.assert_close_e6(@src.vcf_filter_pass_rate(vcf, cfg), 0.6),
  )
}

///|
test "filter_set_filter_mode" {
  let vcf = @src.parse_vcf(filter_vcf_content)
  let cfg = @src.vcf_filter_qual(10.0)
  let result = @src.vcf_filter_set_filter(vcf, cfg, filter_tag="LowQual")
  // 总条数不变
  assert_eq(result.len(), 5)
  // rs1 qual=60 通过，filter 保持 PASS
  assert_eq(result.records[0].filter.length(), 1)
  assert_eq(result.records[0].filter[0], "PASS")
  // rs2 qual=5 失败，filter 应该追加 LowQual
  assert_eq(result.records[1].filter.length(), 2)
  assert_eq(result.records[1].filter[1], "LowQual")
}

///|
test "filter_default_config_all_pass" {
  let vcf = @src.parse_vcf(filter_vcf_content)
  let cfg = @src.vcf_filter_config_default()
  let result = @src.vcf_filter(vcf, cfg)
  assert_eq(result.len(), 5)
}
