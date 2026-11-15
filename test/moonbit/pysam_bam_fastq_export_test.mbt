///|
/// P3 Batch 1: pysam_bam_fastq_export 单元测试

let fq_sam : String =
  "@HD\tVN:1.6\tSO:coordinate\n" +
  "@SQ\tSN:chr1\tLN:1000\n" +
  "read1\t0\tchr1\t100\t60\t8M\t*\t0\t0\tACGTACGT\t!!!!!!!!\n" +
  "read2\t16\tchr1\t200\t50\t8M\t*\t0\t0\tTGCATGCA\t!!!!!!!!\n" +
  "pair_r1\t67\tchr1\t300\t40\t8M\t=\t400\t150\tGGCCGGCC\t!!!!!!!!\n" +
  "pair_r2\t131\tchr1\t400\t40\t8M\t=\t300\t-150\tCCGGCCGG\t!!!!!!!!\n"

fn count_nonempty_lines(s : String) -> Int32 {
  let lines = s.split("\n").to_array()
  let mut count = 0
  for l in lines { if l.length() > 0 { count = count + 1 } }
  count
}

///|
test "fastq_default_emit_all" {
  // bam_to_fastq_default 用 mate_of:0 —— 只 emit 非 paired 或 paired-but-unmapped。
  // fixture: read1/unpaired, read2/unpaired, pair_r1/pair_r2 (paired+mapped → 被过滤)
  let bam = @src.sam_to_bam(fq_sam)
  let fq = @src.bam_to_fastq_default(bam)
  // 2 条 unpaired × 4 行 = 8 行
  assert_eq(count_nonempty_lines(fq), 8)
  assert_true(fq.contains("@read1"))
  assert_true(fq.contains("@read2"))
  // pair_r1/pair_r2 因为 paired+mapped 默认被过滤
  assert_true(!fq.contains("pair_r1"))
  assert_true(!fq.contains("pair_r2"))
}

///|
test "fastq_record_header" {
  let bam = @src.sam_to_bam(fq_sam)
  let r1 = bam.records[0]
  let header = @src.bam_record_fastq_header(r1, "@", false)
  assert_eq(header[0].unsafe_to_char(), '@')
  assert_true(header.contains("read1"))
}

///|
test "fastq_paired_split_outputs" {
  // bam_to_fastq_paired 不分 mate_of —— 所有 mapped paired 进 r1/r2，其余进 singleton
  let bam = @src.sam_to_bam(fq_sam)
  let (r1_text, r2_text, singleton) = @src.bam_to_fastq_paired(bam)
  assert_true(r1_text.contains("@pair_r1"))
  assert_true(r2_text.contains("@pair_r2"))
  assert_true(singleton.contains("@read1"))
  assert_true(singleton.contains("@read2"))
  // pair_r1/pair_r2 各 4 行
  assert_eq(count_nonempty_lines(r1_text), 4)
  assert_eq(count_nonempty_lines(r2_text), 4)
  // singleton 里两条 unpaired
  assert_eq(count_nonempty_lines(singleton), 8)
}

///|
test "fastq_fasta_format" {
  let bam = @src.sam_to_bam(fq_sam)
  let fa = @src.bam_to_fasta(bam)
  // fasta 默认 emit 全部（不过滤 mate）：4 条 × 2 行 = 8 行
  assert_eq(count_nonempty_lines(fa), 8)
  assert_true(fa.contains(">read1"))
  assert_true(fa.contains("ACGTACGT"))
  assert_true(fa.contains(">pair_r1"))
}

///|
test "fastq_config_mate_filter_read1_only" {
  let bam = @src.sam_to_bam(fq_sam)
  let cfg = @src.FastqExportConfig::{
    mate_of: 1, name_prefix: "@",
    original_sequence: false, rg_filter: "",
    emit_qualities: true,
  }
  let fq = @src.bam_to_fastq(bam, cfg)
  // mate_of: 1 → 只 emit paired + mapped + read1
  assert_true(fq.contains("@pair_r1"))
  assert_true(!fq.contains("pair_r2"))
  // read1/read2 flag 0/16 不是 paired → 被 mate_of:1 过滤
  assert_true(!fq.contains("@read1"))
  assert_eq(count_nonempty_lines(fq), 4)
}

///|
test "fastq_config_mate_filter_all_mates" {
  // mate_of:0 应该 emit unpaired + paired singleton（mate-unmapped）
  let sam_singleton : String =
    "@HD\tVN:1.6\tSO:coordinate\n" +
    "@SQ\tSN:chr1\tLN:1000\n" +
    "sgl_pair\t75\tchr1\t100\t60\t8M\t*\t0\t0\tACGTACGT\t!!!!!!!!\n" +  // 75=PAIRED+PROPER+READ1+MATE_UNMAPPED
    "unp_read\t0\tchr1\t200\t50\t8M\t*\t0\t0\tTGCATGCA\t!!!!!!!!\n" +
    "proper_pair\t67\tchr1\t300\t50\t8M\t=\t400\t100\tGGCCGGCC\t!!!!!!!!\n"  // PROPER_PAIR, mate mapped → 不 emit
  let bam = @src.sam_to_bam(sam_singleton)
  let fq = @src.bam_to_fastq_default(bam)
  assert_true(fq.contains("sgl_pair"))    // singleton
  assert_true(fq.contains("unp_read"))    // unpaired
  assert_true(!fq.contains("proper_pair")) // properly paired → 被过滤
}

///|
test "fastq_config_default_roundtrip" {
  let bam = @src.sam_to_bam(fq_sam)
  let cfg = @src.fastq_config_default()
  let fq = @src.bam_to_fastq(bam, cfg)
  assert_eq(count_nonempty_lines(fq), 8)
  assert_true(fq.contains("@"))
  assert_true(fq.contains("+"))
}

///|
test "fastq_empty_bam" {
  let empty = @src.sam_to_bam("@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:1000\n")
  let fq = @src.bam_to_fastq_default(empty)
  assert_eq(fq.length(), 0)
  let fa = @src.bam_to_fasta(empty)
  assert_eq(fa.length(), 0)
}

///|
test "fastq_secondary_supplementary_skipped" {
  // flag 256 (0x100) = secondary；flag 2048 (0x800) = supplementary
  let sam : String =
    "@HD\tVN:1.6\tSO:coordinate\n" +
    "@SQ\tSN:chr1\tLN:1000\n" +
    "primary\t0\tchr1\t100\t60\t8M\t*\t0\t0\tACGTACGT\t!!!!!!!!\n" +
    "sec\t256\tchr1\t100\t30\t8M\t*\t0\t0\tACGTACGT\t!!!!!!!!\n" +
    "supp\t2048\tchr1\t200\t30\t8M\t*\t0\t0\tACGTACGT\t!!!!!!!!\n"
  let bam = @src.sam_to_bam(sam)
  let fq = @src.bam_to_fastq_default(bam)
  assert_eq(count_nonempty_lines(fq), 4)  // 只有 primary
  assert_true(fq.contains("primary"))
  assert_true(!fq.contains("sec"))
  assert_true(!fq.contains("supp"))
}
