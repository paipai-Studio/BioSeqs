///|
/// bam_postprocess_demo — BAM 比对后处理演示:
/// SAM 文本 -> sam_to_bam -> fixmate(mate 坐标/TLEN) ->
/// markdup(配对 mate-key + 光学重复) -> BQSR(碱基质量重校正).

///|
fn bp_demo_sam() -> String {
  "@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:1000\n" +
  // 两对 PCR 重复(同 5' 坐标, mapq 60 vs 30), 光学坐标距离 3/2
  "dup:1:1:100:200\t99\tchr1\t101\t60\t8M\t=\t201\t0\tACGTACGT\tIIIIIIII\n" +
  "dup:1:1:100:200\t147\tchr1\t201\t60\t8M\t=\t101\t0\tTTTTGGGG\tJJJJJJJJ\n" +
  "dup:1:1:103:202\t99\tchr1\t101\t30\t8M\t=\t201\t0\tACGTACGT\tIIIIIIII\n" +
  "dup:1:1:103:202\t147\tchr1\t201\t30\t8M\t=\t101\t0\tTTTTGGGG\tJJJJJJJJ\n" +
  // 一对正常读段
  "norm\t99\tchr1\t401\t60\t8M\t=\t501\t0\tCCGGTTAA\tIIIIIIII\n" +
  "norm\t147\tchr1\t501\t60\t8M\t=\t401\t0\tAACCGGTT\tJJJJJJJJ\n" +
  // 两条单端重复
  "s1\t0\tchr1\t601\t25\t6M\t*\t0\t0\tACGTAC\tIIIIII\n" +
  "s2\t0\tchr1\t601\t20\t6M\t*\t0\t0\tACGTAC\tIIIIII\n"
}

///|
fn main {
  // 1. SAM -> BAM
  let bam = @src.sam_to_bam(bp_demo_sam())
  println(
    "1. sam_to_bam: \{bam.records.length()} 条读段, 参考数 \{bam.num_references()}",
  )
  // 2. fixmate
  let fixed = @src.bam_fixmate(bam)
  println(
    "2. fixmate: norm 对 TLEN = \{fixed.records[4].tlen} / \{fixed.records[5].tlen}, mate 坐标已互填",
  )
  // 3. markdup
  let (marked, stats) = @src.bam_markdup_paired(fixed)
  println("3. markdup:")
  for r in marked.records {
    let dup = if (r.flag & 0x400) != 0 { "DUP" } else { "OK " }
    println("  \{dup} \{r.qname} pos=\{r.pos} mapq=\{r.mapq} tlen=\{r.tlen}")
  }
  println(
    "  组数=\{stats.n_pair_groups}, 重复对=\{stats.n_dup_pairs}, 光学=\{stats.n_optical}, 重复率=\{stats.dup_rate}",
  )
  println(@src.bam_markdup_stats_line(stats))
  // 4. BQSR
  let sb = StringBuilder()
  for _i = 0; _i < 10; _i = _i + 1 {
    sb.write_string("ACGT")
  }
  let reference = sb.to_string()
  let bqsr_sam = "@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:1000\n" +
    "b1\t0\tchr1\t1\t60\t8M\t*\t0\t0\tACGTACGT\t????????\n" +
    "b2\t0\tchr1\t1\t60\t8M\t*\t0\t0\tACGTACGT\t????????\n" +
    "b3\t0\tchr1\t1\t60\t8M\t*\t0\t0\tACGAACGT\t????????\n" +
    "b4\t0\tchr1\t1\t60\t8M\t*\t0\t0\tACGAACGT\t????????\n"
  let bqsr_bam = @src.sam_to_bam(bqsr_sam)
  // 对照: 无 known-sites 时错配位显著降质
  let (recal0, report0) = @src.bam_bqsr(bqsr_bam, reference)
  println("4. BQSR 无 known-sites 对照:")
  println("  b3(错配) qual: \{recal0.records[2].qual} (i=3 显著降质)")
  println("  观察=\{report0.n_observations}, 错配=\{report0.n_errors}")
  // known-sites(0-based [3,4)) 豁免该坐标上的全部读段位 -> 错配不计入训练
  let (recal, report) = @src.bam_bqsr(bqsr_bam, reference, known_sites=[
    (0, 3, 4),
  ])
  println("  BQSR 带 known-sites:")
  println("  b1(无错) qual: \{recal.records[0].qual}")
  println("  b3(错配) qual: \{recal.records[2].qual} (i=3 保持 30)")
  println(
    "  观察=\{report.n_observations}, 错配=\{report.n_errors}, 表项=\{report.n_table_entries}",
  )
  println(@src.bam_bqsr_report_line(report))
  println("演示完成: SAM -> BAM -> fixmate -> markdup -> BQSR 后处理链")
}
