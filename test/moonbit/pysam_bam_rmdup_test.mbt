///|
/// P3: pysam_bam_rmdup 单元测试

// SAM fixture: 4 条 reads
//   r1: qname=readA, chr1, pos=100, mapq=60
//   r2: qname=readA, chr1, pos=100, mapq=30  ← dup of r1
//   r3: qname=readB, chr1, pos=200, mapq=50
//   r4: qname=readC, *,   pos=0,   mapq=10  ← unmapped (always kept)
let rmdup_sam : String = "@HD\tVN:1.6\tSO:coordinate\n" +
  "@SQ\tSN:chr1\tLN:1000\n" +
  "readA\t0\tchr1\t100\t60\t8M\t*\t0\t0\tACGTACGT\t!!!!!!!!\n" +
  "readA\t0\tchr1\t100\t30\t8M\t*\t0\t0\tACGTACGT\t!!!!!!!!\n" +
  "readB\t0\tchr1\t200\t50\t8M\t*\t0\t0\tTGCATGCA\t!!!!!!!!\n" +
  "readC\t4\t*\t0\t10\t*\t*\t0\t0\tNNNNNNNN\t!!!!!!!!\n"

///|
test "rmdup_basic_drop_duplicates" {
  let bam = @src.sam_to_bam(rmdup_sam)
  assert_eq(bam.len(), 4)
  let result = @src.bam_rmdup(bam)
  // r1 (mapq=60) kept, r2 dropped (dup), r3 kept, r4 unmapped kept → 3 条
  assert_eq(result.len(), 3)
  // 保留的应该是 mapq 更高的那条
  assert_eq(result.records[0].mapq, 60)
  assert_eq(result.records[0].qname, "readA")
  assert_eq(result.records[1].qname, "readB")
  assert_eq(result.records[2].qname, "readC")
}

///|
test "rmdup_count_duplicates" {
  let bam = @src.sam_to_bam(rmdup_sam)
  assert_eq(@src.bam_count_duplicates(bam), 1)
  assert_true(@moonbit.assert_close_e6(@src.bam_duplicate_rate(bam), 0.25))
}

///|
test "rmdup_no_duplicates" {
  // 每条 read 都在不同 pos，没 dup
  let sam_no_dup : String = "@HD\tVN:1.6\tSO:coordinate\n" +
    "@SQ\tSN:chr1\tLN:1000\n" +
    "r1\t0\tchr1\t100\t60\t8M\t*\t0\t0\tACGTACGT\t!!!!!!!!\n" +
    "r2\t0\tchr1\t200\t30\t8M\t*\t0\t0\tACGTACGT\t!!!!!!!!\n" +
    "r3\t0\tchr1\t300\t50\t8M\t*\t0\t0\tACGTACGT\t!!!!!!!!\n"
  let bam = @src.sam_to_bam(sam_no_dup)
  let result = @src.bam_rmdup(bam)
  assert_eq(result.len(), 3)
  assert_eq(@src.bam_count_duplicates(bam), 0)
}

///|
test "rmdup_all_duplicates_same_group" {
  // 4 条完全相同 qname+chr+pos 的 reads → rmdup 只保留 mapq 最高的那条
  let sam_multi : String = "@HD\tVN:1.6\tSO:coordinate\n" +
    "@SQ\tSN:chr1\tLN:1000\n" +
    "dup_read\t0\tchr1\t100\t10\t8M\t*\t0\t0\tACGTACGT\t!!!!!!!!\n" +
    "dup_read\t0\tchr1\t100\t60\t8M\t*\t0\t0\tACGTACGT\t!!!!!!!!\n" +
    "dup_read\t0\tchr1\t100\t30\t8M\t*\t0\t0\tACGTACGT\t!!!!!!!!\n" +
    "dup_read\t0\tchr1\t100\t20\t8M\t*\t0\t0\tACGTACGT\t!!!!!!!!\n"
  let bam = @src.sam_to_bam(sam_multi)
  let result = @src.bam_rmdup(bam)
  assert_eq(result.len(), 1)
  assert_eq(result.records[0].mapq, 60)
  assert_eq(@src.bam_count_duplicates(bam), 3)
}

///|
test "rmdup_markdup_mode_sets_flag_no_drop" {
  let bam = @src.sam_to_bam(rmdup_sam)
  let result = @src.bam_markdup_rmdup_mode(bam)
  // 4 条全保留，但 r2 应该被打 duplicate flag (0x400 = 1024)
  assert_eq(result.len(), 4)
  let r1 = result.records[0]
  let r2 = result.records[1]
  // r1 是 mapq 高的那条，不被标 dup
  assert_true(!r1.is_duplicate())
  // r2 是 dup
  assert_true(r2.is_duplicate())
}

///|
test "rmdup_duplicate_rate_empty" {
  let sam_empty : String = "@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:1000\n"
  let bam = @src.sam_to_bam(sam_empty)
  assert_eq(@src.bam_duplicate_rate(bam), 0.0)
}

///|
test "rmdup_paired_mode" {
  // qname=pair1 出现两次（一个 pair），pair2 出现一次
  let sam_pe : String = "@HD\tVN:1.6\tSO:coordinate\n" +
    "@SQ\tSN:chr1\tLN:1000\n" +
    "pair1\t67\tchr1\t100\t60\t8M\t=\t200\t150\tACGTACGT\t!!!!!!!!\n" +
    "pair1\t131\tchr1\t200\t60\t8M\t=\t100\t-150\tTGCATGCA\t!!!!!!!!\n" +
    "pair1\t0\tchr1\t100\t30\t8M\t*\t0\t0\tACGTACGT\t!!!!!!!!\n" +
    "pair2\t0\tchr1\t300\t50\t8M\t*\t0\t0\tACGTACGT\t!!!!!!!!\n"
  let bam = @src.sam_to_bam(sam_pe)
  let result = @src.bam_rmdup_paired(bam)
  // pair1 有 3 条（2+1），pair2 有 1 条
  // paired mode 下 pair1 的前两条（flag 67/131）是 pair，都保留
  // 第三条 qname=pair1 是单端 dup，不属 pair，用 single-end 规则 → drop
  // pair2 保留
  assert_eq(result.len(), 3)
}
