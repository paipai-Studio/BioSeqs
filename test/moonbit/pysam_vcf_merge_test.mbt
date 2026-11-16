///|
/// P3 Batch 5: pysam_vcf_merge 单元测试
let merge_vcf_a : String = "##fileformat=VCFv4.2\n" +
  "##contig=<ID=chr1,length=1000>\n" +
  "##contig=<ID=chr2,length=1000>\n" +
  "##INFO=<ID=DP,Number=1,Type=Integer,Description=\"depth\">\n" +
  "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n" +
  "chr1\t100\ta1\tA\tT\t60\tPASS\tDP=10\n" +
  "chr1\t300\ta2\tC\tG\t50\tPASS\tDP=20\n" +
  "chr2\t50\ta3\tG\tA\t40\tPASS\tDP=15\n"

///|
let merge_vcf_b : String = "##fileformat=VCFv4.2\n" +
  "##contig=<ID=chr1,length=1000>\n" +
  "##contig=<ID=chr2,length=800>\n" +
  "##contig=<ID=chr3,length=500>\n" +
  "##INFO=<ID=AF,Number=A,Type=Float,Description=\"allele freq\">\n" +
  "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n" +
  "chr1\t200\tb1\tC\tT\t55\tPASS\tAF=0.3\n" +
  "chr1\t400\tb2\tT\tC\t35\tPASS\tAF=0.7\n" +
  "chr3\t100\tb3\tA\tG\t45\tPASS\tAF=0.1\n"

///|
test "vcf_merge_sorted_two_files" {
  let a = @src.parse_vcf(merge_vcf_a)
  let b = @src.parse_vcf(merge_vcf_b)
  let merged = @src.vcf_merge_sorted([a, b])
  assert_eq(merged.records.length(), 6)
  // (chrom, pos) 排序：chr1:100, chr1:200, chr1:300, chr1:400, chr2:50, chr3:100
  assert_eq(merged.records[0].pos, 100)
  assert_eq(merged.records[1].pos, 200)
  assert_eq(merged.records[2].pos, 300)
  assert_eq(merged.records[3].pos, 400)
  assert_eq(merged.records[4].chrom, "chr2")
  assert_eq(merged.records[4].pos, 50)
  assert_eq(merged.records[5].chrom, "chr3")
  assert_eq(merged.records[5].pos, 100)
}

///|
test "vcf_merge_unsorted_preserves_input_order" {
  let a = @src.parse_vcf(merge_vcf_a)
  let b = @src.parse_vcf(merge_vcf_b)
  let merged = @src.vcf_merge_unsorted([a, b])
  assert_eq(merged.records.length(), 6)
  // a 的 3 条先出现，然后 b 的 3 条
  assert_eq(merged.records[0].id[0], "a1")
  assert_eq(merged.records[1].id[0], "a2")
  assert_eq(merged.records[2].id[0], "a3")
  assert_eq(merged.records[3].id[0], "b1")
}

///|
test "vcf_merge_meta_headers_dedup" {
  let a = @src.parse_vcf(merge_vcf_a)
  let b = @src.parse_vcf(merge_vcf_b)
  let merged = @src.merge_vcf_meta([a, b])
  // contig chr1 在两个 VCF 都有，但 value 不同 (1000 vs 1000) → 相同 dedup
  // contig chr2 也重复出现（a 有 chr2 length=1000，b 有 chr2 length=800）→ value 不同，都保留
  // 注意：merge_vcf_meta 按 key+value 去重，chr2 length 不同所以两个都保留
  // contig chr3 只在 b 里
  let mut chr1_count = 0
  let mut chr2_count = 0
  let mut chr3_count = 0
  let mut dp_count = 0
  let mut af_count = 0
  for mh in merged {
    if mh.key == "contig" && mh.value.contains("chr1") {
      chr1_count = chr1_count + 1
    }
    if mh.key == "contig" && mh.value.contains("chr2") {
      chr2_count = chr2_count + 1
    }
    if mh.key == "contig" && mh.value.contains("chr3") {
      chr3_count = chr3_count + 1
    }
    if mh.key == "INFO" && mh.value.contains("ID=DP") {
      dp_count = dp_count + 1
    }
    if mh.key == "INFO" && mh.value.contains("ID=AF") {
      af_count = af_count + 1
    }
  }
  assert_eq(chr1_count, 1)
  assert_eq(chr2_count, 2) // length 不同，都保留
  assert_eq(chr3_count, 1)
  assert_eq(dp_count, 1)
  assert_eq(af_count, 1)
}

///|
test "vcf_merge_single_file_passthrough" {
  let a = @src.parse_vcf(merge_vcf_a)
  let merged = @src.vcf_merge_sorted([a])
  assert_eq(merged.records.length(), 3)
  assert_eq(merged.records[0].id[0], "a1")
}

///|
test "vcf_merge_empty_input" {
  let merged = @src.vcf_merge_sorted([])
  assert_eq(merged.records.length(), 0)
  match merged.fixed_header {
    Some(fh) => assert_true(fh.columns.length() >= 8)
    None => abort("expected fixed header")
  }
}

///|
test "vcf_merge_stats_totals" {
  let a = @src.parse_vcf(merge_vcf_a)
  let b = @src.parse_vcf(merge_vcf_b)
  let merged = @src.vcf_merge_sorted([a, b])
  let (total_in, total_out, unique_contigs) = @src.vcf_merge_stats(merged, [
    a, b,
  ])
  assert_eq(total_in, 6)
  assert_eq(total_out, 6)
  assert_eq(unique_contigs, 3)
}

///|
test "vcf_merge_three_files" {
  let a = @src.parse_vcf(merge_vcf_a)
  let b = @src.parse_vcf(merge_vcf_b)
  let c_text : String = "##fileformat=VCFv4.2\n##contig=<ID=chr1,length=1000>\n" +
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n" +
    "chr1\t50\tc1\tA\tC\t60\tPASS\t.\n" +
    "chr3\t200\tc2\tG\tT\t45\tPASS\t.\n"
  let c = @src.parse_vcf(c_text)
  let merged = @src.vcf_merge_sorted([a, b, c])
  assert_eq(merged.records.length(), 8)
  assert_eq(merged.records[0].pos, 50) // c1 最小
  assert_eq(merged.records[1].pos, 100)
}

///|
test "vcf_merge_sort_records_only_preserves_header" {
  let a = @src.parse_vcf(merge_vcf_a)
  // 乱序 records
  let b_text : String = "##fileformat=VCFv4.2\n##contig=<ID=chr1,length=1000>\n" +
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n" +
    "chr1\t300\tb1\tC\tG\t50\tPASS\t.\n" +
    "chr1\t100\tb2\tA\tT\t60\tPASS\t.\n" +
    "chr1\t200\tb3\tG\tA\t40\tPASS\t.\n"
  let b = @src.parse_vcf(b_text)
  let merged = @src.vcf_merge_sort_records_only(b)
  assert_eq(merged.records[0].pos, 100)
  assert_eq(merged.records[1].pos, 200)
  assert_eq(merged.records[2].pos, 300)
  // 原 header 不变
  match merged.fixed_header {
    Some(fh) => assert_true(fh.columns.length() >= 8)
    None => abort("expected fixed header")
  }
}

///|
test "vcf_merge_duplicate_positions_keeps_both" {
  // bcftools merge 默认不去重同 pos 的 variants，两个都保留
  let a_text : String = "##fileformat=VCFv4.2\n##contig=<ID=chr1,length=1000>\n" +
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n" +
    "chr1\t100\t.\tA\tT\t60\tPASS\t.\n"
  let b_text : String = "##fileformat=VCFv4.2\n##contig=<ID=chr1,length=1000>\n" +
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n" +
    "chr1\t100\t.\tA\tC\t50\tPASS\t.\n"
  let merged = @src.vcf_merge_sorted([
    @src.parse_vcf(a_text),
    @src.parse_vcf(b_text),
  ])
  assert_eq(merged.records.length(), 2)
  assert_eq(merged.records[0].pos, 100)
  assert_eq(merged.records[1].pos, 100)
  // 按输入顺序 a 先 b 后
  assert_eq(merged.records[0].alt[0], "T")
  assert_eq(merged.records[1].alt[0], "C")
}
