// ===========================================================================
// qtlseq_test.mbt — QTL-seq / BSA 混池分离分析测试.
// ===========================================================================

///|
/// 双池 VCF 夹具: chr1 中性区域(两池 0.5/0.5) + 峰值区域 [4100,4300)
/// (mut 0.9 / wt 0.05) + 低深度位点 + 缺 AD 位点.
fn qs_test_vcf() -> String {
  "##fileformat=VCFv4.2\n" +
  "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tMUT\tWT\n" +
  "chr1\t1000\t.\tA\tG\t.\tPASS\t.\tGT:AD\t0/1:10,10\t0/1:10,10\n" +
  "chr1\t1100\t.\tC\tT\t.\tPASS\t.\tGT:AD\t0/1:12,12\t0/1:10,10\n" +
  "chr1\t4100\t.\tA\tG\t.\tPASS\t.\tGT:AD\t0/1:1,9\t0/1:19,1\n" +
  "chr1\t4200\t.\tC\tT\t.\tPASS\t.\tGT:AD\t0/1:1,9\t0/1:19,1\n" +
  "chr1\t4300\t.\tA\tG\t.\tPASS\t.\tGT:AD\t0/1:1,1\t0/1:19,1\n" +
  "chr1\t4400\t.\tC\tT\t.\tPASS\t.\tGT:AD\t0/1:9,9\t0/1:.\n" +
  "chr1\t9100\t.\tA\tG\t.\tPASS\t.\tGT:AD\t0/1:10,10\t0/1:10,10\n" +
  "chr1\t9200\t.\tC\tT\t.\tPASS\t.\tGT:AD\t0/1:12,12\t0/1:10,10\n"
}

///|
fn qs_small_params() -> @src.QtlseqParams {
  { min_depth: 3, window: 1000, step: 500, delta_threshold: 0.5, }
}

///|
test "qtlseq: SNP-index 计算(双池)" {
  let vcf = @src.parse_vcf(qs_test_vcf())
  let sites = @src.qtlseq_snp_index(
    vcf.records,
    "MUT",
    wt_bulk="WT",
    params=qs_small_params(),
  )
  // 8 条记录: 4300(深度 2)与 4400(WT 无 AD)被跳过 -> 6 位点
  inspect(sites.length(), content="6")
  let s0 = sites[0]
  inspect(s0.pos, content="1000")
  assert_eq(s0.mut_index, 0.5)
  assert_eq(s0.wt_index, 0.5)
  assert_eq(s0.delta, 0.0)
  // 峰值位点 4100: mut 9/10=0.9, wt 1/20=0.05
  let peak = sites[2]
  inspect(peak.pos, content="4100")
  assert_eq(peak.mut_index, 0.9)
  assert_eq(peak.wt_index, 0.05)
  assert_eq(peak.delta, 0.85)
  inspect(peak.mut_depth, content="10")
  inspect(peak.wt_depth, content="20")
}

///|
test "qtlseq: MutMap 单池模式" {
  let vcf = @src.parse_vcf(qs_test_vcf())
  let sites = @src.qtlseq_snp_index(
    vcf.records,
    "MUT",
    params=qs_small_params(),
  )
  // 单池: WT 不检查 -> 4300 深度 2 仍跳过, 4400 MUT 9,9 保留
  inspect(sites.length(), content="7")
  let s0 = sites[0]
  assert_eq(s0.wt_index, 0.0)
  assert_eq(s0.delta, 0.5)
  inspect(s0.wt_depth, content="0")
}

///|
test "qtlseq: 滑窗聚合" {
  let vcf = @src.parse_vcf(qs_test_vcf())
  let sites = @src.qtlseq_snp_index(
    vcf.records,
    "MUT",
    wt_bulk="WT",
    params=qs_small_params(),
  )
  let windows = @src.qtlseq_windows(sites, params=qs_small_params())
  // 网格起点 = 1000 - ((1000-1)/500)*500 = 500; 覆盖到 9200
  // 窗口: 500-1500(1 位点), 4000-5000(2 位点), 9000-10000(2 位点) 等
  let mut peak_win : @src.QtlseqWindow? = None
  for w in windows {
    if w.start <= 4100 && w.end_ >= 4200 && w.n_snps == 2 {
      peak_win = Some(w)
    }
  }
  let pw = peak_win.unwrap()
  assert_eq(pw.mean_mut_index, 0.9)
  assert_eq(pw.mean_wt_index, 0.05)
  assert_eq(pw.mean_delta, 0.85)
  // 中性窗口(500-1499 覆盖 1000/1100 两位点, delta 0)
  let mut neutral : @src.QtlseqWindow? = None
  for w in windows {
    if w.start == 500 && w.end_ == 1499 {
      neutral = Some(w)
    }
  }
  let nw = neutral.unwrap()
  inspect(nw.n_snps, content="2")
  assert_eq(nw.mean_delta, 0.0)
}

///|
test "qtlseq: 候选区间(峰值合并)" {
  let vcf = @src.parse_vcf(qs_test_vcf())
  let sites = @src.qtlseq_snp_index(
    vcf.records,
    "MUT",
    wt_bulk="WT",
    params=qs_small_params(),
  )
  let windows = @src.qtlseq_windows(sites, params=qs_small_params())
  let intervals = @src.qtlseq_intervals(windows, params=qs_small_params())
  // 阈值 0.5: 仅峰值窗口(|0.8|)入区间; 相邻显著窗合并为 1 个
  inspect(intervals.length(), content="1")
  let iv = intervals[0]
  inspect(iv.chrom, content="chr1")
  assert_eq(iv.peak_delta > 0.7, true)
  assert_eq(iv.direction, "+")
  assert_eq(iv.n_snps >= 2, true)
  // 中性区(|0|)不入区间
}

///|
test "qtlseq: 野生型富集方向" {
  // 构造 wt 富集峰(mut 0.1 / wt 0.9)
  let vcf_text = "##fileformat=VCFv4.2\n" +
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tMUT\tWT\n" +
    "chr1\t2100\t.\tA\tG\t.\tPASS\t.\tGT:AD\t0/1:9,1\t0/1:1,9\n" +
    "chr1\t2200\t.\tC\tT\t.\tPASS\t.\tGT:AD\t0/1:9,1\t0/1:1,9\n"
  let vcf = @src.parse_vcf(vcf_text)
  let sites = @src.qtlseq_snp_index(
    vcf.records,
    "MUT",
    wt_bulk="WT",
    params=qs_small_params(),
  )
  let windows = @src.qtlseq_windows(sites, params=qs_small_params())
  let intervals = @src.qtlseq_intervals(windows, params=qs_small_params())
  inspect(intervals.length(), content="1")
  inspect(intervals[0].direction, content="-")
}

///|
test "qtlseq: TSV/BED/JSON/摘要导出" {
  let vcf = @src.parse_vcf(qs_test_vcf())
  let sites = @src.qtlseq_snp_index(
    vcf.records,
    "MUT",
    wt_bulk="WT",
    params=qs_small_params(),
  )
  let windows = @src.qtlseq_windows(sites, params=qs_small_params())
  let intervals = @src.qtlseq_intervals(windows, params=qs_small_params())
  let tsv = @src.qtlseq_windows_tsv(windows)
  assert_eq(tsv.contains("chrom\tstart\tend\tn_snps"), true)
  assert_eq(tsv.contains("0.9\t0.05\t0.85"), true)
  let bed = @src.qtlseq_intervals_bed(intervals)
  // BED 0-based: start-1
  assert_eq(bed.contains("chr1\t"), true)
  assert_eq(bed.contains("dir=+"), true)
  assert_eq(bed.contains("peak="), true)
  let json = @src.qtlseq_json(sites, windows, intervals)
  assert_eq(json.contains("\"n_intervals\""), true)
  assert_eq(json.contains("\"mean_delta\""), true)
  assert_eq(json.contains("\"direction\""), true)
  let line = @src.qtlseq_summary_line(sites, windows, intervals)
  assert_eq(line.contains("sites=6"), true)
  assert_eq(line.contains("intervals=1"), true)
  assert_eq(line.contains("mutant-enriched=1"), true)
}

///|
test "qtlseq: 空输入与非 SNP 过滤" {
  let sites = @src.qtlseq_snp_index([], "MUT", params=qs_small_params())
  inspect(sites.length(), content="0")
  let windows = @src.qtlseq_windows(sites, params=qs_small_params())
  inspect(windows.length(), content="0")
  let intervals = @src.qtlseq_intervals(windows, params=qs_small_params())
  inspect(intervals.length(), content="0")
  // indel 记录被过滤
  let vcf_text = "##fileformat=VCFv4.2\n" +
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tMUT\tWT\n" +
    "chr1\t1000\t.\tAT\tA\t.\tPASS\t.\tGT:AD\t0/1:1,9\t0/1:9,1\n"
  let vcf = @src.parse_vcf(vcf_text)
  let sites2 = @src.qtlseq_snp_index(
    vcf.records,
    "MUT",
    wt_bulk="WT",
    params=qs_small_params(),
  )
  inspect(sites2.length(), content="0")
}
