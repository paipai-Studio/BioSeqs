///|
/// P3: pysam_vcf_isec 单元测试

let isec_vcf_a : String =
  "##fileformat=VCFv4.2\n" +
  "##contig=<ID=chr1,length=1000>\n" +
  "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n" +
  "chr1\t100\ta1\tA\tT\t60\tPASS\t.\n" +
  "chr1\t200\ta2\tC\tG\t50\tPASS\t.\n" +
  "chr1\t300\ta3\tG\tA\t40\tPASS\t.\n" +
  "chr1\t400\ta4\tT\tC\t30\tPASS\t.\n"

let isec_vcf_b : String =
  "##fileformat=VCFv4.2\n" +
  "##contig=<ID=chr1,length=1000>\n" +
  "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n" +
  "chr1\t200\tb1\tC\tG\t55\tPASS\t.\n" +   // matches a2
  "chr1\t300\tb2\tG\tA\t35\tPASS\t.\n" +   // matches a3 (ref same)
  "chr1\t500\tb3\tA\tG\t45\tPASS\t.\n"    // only in b

///|
test "isec_intersection_exact" {
  let a = @src.parse_vcf(isec_vcf_a)
  let b = @src.parse_vcf(isec_vcf_b)
  let result = @src.vcf_isec_intersection(a, b)
  // pos=200, 300 匹配 → 2 条
  assert_eq(result.len(), 2)
  assert_eq(result.records[0].pos, 200)
  assert_eq(result.records[1].pos, 300)
}

///|
test "isec_intersection_fast_same" {
  let a = @src.parse_vcf(isec_vcf_a)
  let b = @src.parse_vcf(isec_vcf_b)
  let slow = @src.vcf_isec_intersection(a, b)
  let fast = @src.vcf_isec_intersection_fast(a, b)
  assert_eq(slow.len(), fast.len())
}

///|
test "isec_unique_first" {
  let a = @src.parse_vcf(isec_vcf_a)
  let b = @src.parse_vcf(isec_vcf_b)
  let only_a = @src.vcf_isec_unique_first(a, b)
  // a 中 pos 100, 400 不在 b 里
  assert_eq(only_a.len(), 2)
  assert_eq(only_a.records[0].pos, 100)
  assert_eq(only_a.records[1].pos, 400)
}

///|
test "isec_unique_second" {
  let a = @src.parse_vcf(isec_vcf_a)
  let b = @src.parse_vcf(isec_vcf_b)
  let only_b = @src.vcf_isec_unique_second(a, b)
  assert_eq(only_b.len(), 1)
  assert_eq(only_b.records[0].pos, 500)
}

///|
test "isec_union" {
  let a = @src.parse_vcf(isec_vcf_a)
  let b = @src.parse_vcf(isec_vcf_b)
  let u = @src.vcf_isec_union(a, b)
  // a 有 4, b 有 3, 共享 2 → 4+3-2 = 5
  assert_eq(u.len(), 5)
}

///|
test "isec_summary_count" {
  let a = @src.parse_vcf(isec_vcf_a)
  let b = @src.parse_vcf(isec_vcf_b)
  let (only_a, only_b, shared) = @src.vcf_isec_summary(a, b)
  assert_eq(only_a, 2)
  assert_eq(only_b, 1)
  assert_eq(shared, 2)
}

///|
test "isec_jaccard" {
  let a = @src.parse_vcf(isec_vcf_a)
  let b = @src.parse_vcf(isec_vcf_b)
  let j = @src.vcf_isec_jaccard(a, b)
  // J = |∩| / |∪| = 2 / 5 = 0.4
  assert_true(@moonbit.assert_close_e6(j, 0.4))
}

///|
test "isec_intersect_all_n_way" {
  let a = @src.parse_vcf(isec_vcf_a)
  let b_text : String =
    "##fileformat=VCFv4.2\n##contig=<ID=chr1,length=1000>\n" +
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n" +
    "chr1\t200\t.\tC\tG\t55\tPASS\t.\n" +
    "chr1\t300\t.\tG\tA\t35\tPASS\t.\n" +
    "chr1\t999\t.\tA\tG\t45\tPASS\t.\n"
  let c_text : String =
    "##fileformat=VCFv4.2\n##contig=<ID=chr1,length=1000>\n" +
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n" +
    "chr1\t200\t.\tC\tG\t55\tPASS\t.\n" +
    "chr1\t300\t.\tG\tA\t35\tPASS\t.\n" +
    "chr1\t1000\t.\tA\tG\t45\tPASS\t.\n"
  let results = @src.vcf_isec_intersect_all([a, @src.parse_vcf(b_text), @src.parse_vcf(c_text)])
  // 三个 VCF 都有 pos=200, 300
  assert_eq(results.len(), 2)
}

///|
test "isec_empty_input" {
  let a = @src.parse_vcf(isec_vcf_a)
  let empty_vcf = @src.parse_vcf("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
  let r = @src.vcf_isec_intersection(a, empty_vcf)
  assert_eq(r.len(), 0)
  let r2 = @src.vcf_isec_intersect_all([])
  assert_eq(r2.len(), 0)
}

///|
test "isec_match_mode_exact_with_ref" {
  let a = @src.parse_vcf(isec_vcf_a)
  let b = @src.parse_vcf(isec_vcf_b)
  let result = @src.vcf_isec_intersection_fast(a, b, mode=@src.isec_mode_exact_with_ref())
  assert_eq(result.len(), 2)
}
