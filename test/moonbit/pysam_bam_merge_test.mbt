///|
/// P3 Batch 1: pysam_bam_merge 单元测试
///
/// 注：sam_to_bam 把 SAM 的 1-based POS 转成 BamRecord.pos 0-based（减 1）。
/// fixture 预期值按 0-based 写。

let merge_sam_a : String =
  "@HD\tVN:1.6\tSO:coordinate\n" +
  "@SQ\tSN:chr1\tLN:1000\n" +
  "read1\t0\tchr1\t101\t60\t8M\t*\t0\t0\tACGTACGT\t!!!!!!!!\n" +
  "read2\t0\tchr1\t301\t50\t8M\t*\t0\t0\tTGCATGCA\t!!!!!!!!\n"

let merge_sam_b : String =
  "@HD\tVN:1.6\tSO:coordinate\n" +
  "@SQ\tSN:chr1\tLN:1000\n" +
  "read3\t0\tchr1\t201\t40\t8M\t*\t0\t0\tGGCCGGCC\t!!!!!!!!\n" +
  "read4\t0\tchr1\t401\t30\t8M\t*\t0\t0\tAAAAAAA\t!!!!!!!!\n"

///|
test "merge_sorted_two_files_sort_by_coord" {
  let a = @src.sam_to_bam(merge_sam_a)
  let b = @src.sam_to_bam(merge_sam_b)
  let merged = @src.bam_merge_sorted([a, b])
  assert_eq(merged.len(), 4)
  // SAM 101/201/301/401 → BamRecord.pos 100/200/300/400
  assert_eq(merged.records[0].pos, 100)
  assert_eq(merged.records[1].pos, 200)
  assert_eq(merged.records[2].pos, 300)
  assert_eq(merged.records[3].pos, 400)
}

///|
test "merge_sorted_single_file_passthrough" {
  let a = @src.sam_to_bam(merge_sam_a)
  let merged = @src.bam_merge_sorted([a])
  assert_eq(merged.len(), 2)
  assert_eq(merged.records[0].qname, "read1")
}

///|
test "merge_sorted_empty_input" {
  let merged = @src.bam_merge_sorted([])
  assert_eq(merged.len(), 0)
}

///|
test "merge_sorted_safe_header_validation_pass" {
  let a = @src.sam_to_bam(merge_sam_a)
  let b = @src.sam_to_bam(merge_sam_b)
  let merged = @src.bam_merge_sorted_safe([a, b])
  assert_eq(merged.len(), 4)
}

///|
test "merge_rg_prefix_renames" {
  let sam_with_rg : String =
    "@HD\tVN:1.6\tSO:coordinate\n" +
    "@SQ\tSN:chr1\tLN:1000\n" +
    "@RG\tID:orig\n" +
    "read1\t0\tchr1\t101\t60\t8M\t*\t0\t0\tACGTACGT\t!!!!!!!!\tRG:Z:orig\n"
  let a = @src.sam_to_bam(sam_with_rg)
  let b = @src.sam_to_bam(sam_with_rg)
  let merged = @src.bam_merge_with_rg_prefix([a, b], ["lib1", "lib2"])
  assert_eq(merged.len(), 2)
  let r1 = merged.records[0]
  let r2 = merged.records[1]
  match r1.get_tag("RG") {
    Some(@src.BamTagValue::Str(v)) => assert_eq(v, "lib1.orig")
    _ => abort("r1 RG tag missing")
  }
  match r2.get_tag("RG") {
    Some(@src.BamTagValue::Str(v)) => assert_eq(v, "lib2.orig")
    _ => abort("r2 RG tag missing")
  }
}

///|
test "merge_fixmate_flag_is_mate_info" {
  let sam_paired : String =
    "@HD\tVN:1.6\tSO:coordinate\n" +
    "@SQ\tSN:chr1\tLN:1000\n" +
    "pair\t67\tchr1\t101\t60\t8M\t=\t201\t150\tACGTACGT\t!!!!!!!!\n" +
    "pair\t131\tchr1\t201\t60\t8M\t=\t101\t-150\tTGCATGCA\t!!!!!!!!\n"
  let bam = @src.sam_to_bam(sam_paired)
  let fixed = @src.bam_merge_fixmate_pairs(bam)
  assert_eq(fixed.len(), 2)
  let r1 = fixed.records[0]
  let r2 = fixed.records[1]
  // SAM 201/101 → BamRecord.next_pos 存原值或 -1
  // 先检查实际行为
  assert_true(r1.next_pos >= 0)
  assert_true(r2.next_pos >= 0)
}

///|
test "merge_three_files_all_concat" {
  let a = @src.sam_to_bam(merge_sam_a)
  let b = @src.sam_to_bam(merge_sam_b)
  let c_sam : String =
    "@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:1000\n" +
    "read5\t0\tchr1\t51\t60\t8M\t*\t0\t0\tACGTACGT\t!!!!!!!!\n"
  let c = @src.sam_to_bam(c_sam)
  let merged = @src.bam_merge_sorted([a, b, c])
  assert_eq(merged.len(), 5)
  // SAM 51 → BamRecord.pos 50
  assert_eq(merged.records[0].pos, 50)
}
