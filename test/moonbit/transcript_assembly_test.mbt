// ===========================================================================
// transcript_assembly_test.mbt — StringTie 风格转录本组装测试.
// ===========================================================================

///|
/// SAM 头(单染色体 chr1).
fn ta_sam_header() -> String {
  "@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:2000\n"
}

///|
/// 三 isoform 场景读段(0-based 坐标语义):
/// E1 [10,100), E2 [200,300), E3 [400,440), E4 [600,640);
/// junction (100,200) x12, (300,400) x4, (300,600) x2;
/// 另加 20 条 E1 全覆盖未剪接读段.
fn ta_sam_three_isoforms() -> String {
  let sb = StringBuilder()
  sb.write_string(ta_sam_header())
  // 12 条跨 junction (100,200): pos 81, 20M100N40M
  for i = 0; i < 12; i = i + 1 {
    sb.write_string(
      "j1_\{i}\t0\tchr1\t81\t60\t20M100N40M\t*\t0\t0\tACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACG\t*\n",
    )
  }
  // 4 条跨 junction (300,400): pos 281, 20M100N40M
  for i = 0; i < 4; i = i + 1 {
    sb.write_string(
      "j2_\{i}\t0\tchr1\t281\t60\t20M100N40M\t*\t0\t0\tACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACG\t*\n",
    )
  }
  // 2 条跨 junction (300,600): pos 281, 20M300N40M
  for i = 0; i < 2; i = i + 1 {
    sb.write_string(
      "j3_\{i}\t0\tchr1\t281\t60\t20M300N40M\t*\t0\t0\tACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACG\t*\n",
    )
  }
  // 20 条 E1 全覆盖未剪接: pos 11, 90M
  for i = 0; i < 20; i = i + 1 {
    sb.write_string(
      "u_\{i}\t0\tchr1\t11\t60\t90M\t*\t0\t0\tACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACG\t*\n",
    )
  }
  sb.to_string()
}

///|
test "transcript_assembly: junction 提取" {
  let junctions = @src.ta_junctions_from_sam("chr1", ta_sam_three_isoforms())
  inspect(junctions.length(), content="3")
  // 排序后: (100,200) (300,400) (300,600)
  inspect(junctions[0].donor, content="100")
  inspect(junctions[0].acceptor, content="200")
  inspect(junctions[0].count, content="12")
  inspect(junctions[1].donor, content="300")
  inspect(junctions[1].acceptor, content="400")
  inspect(junctions[1].count, content="4")
  inspect(junctions[2].donor, content="300")
  inspect(junctions[2].acceptor, content="600")
  inspect(junctions[2].count, content="2")
}

///|
test "transcript_assembly: 三 isoform 组装" {
  let ts = @src.ta_assemble_transcripts(
    "chr1",
    ta_sam_three_isoforms(),
    strand="+",
  )
  // 贪心流: T1 = E1-E2-E3(瓶颈 4), T2 = E1-E2-E4(瓶颈 2),
  // 剩余 E1 深度 >= 2 -> T3 = E1 单外显子
  inspect(ts.length(), content="3")
  // T1: 三外显子
  let t1 = ts[0]
  inspect(t1.exons.length(), content="3")
  inspect(t1.exons[0].0, content="10")
  inspect(t1.exons[0].1, content="100")
  inspect(t1.exons[1].0, content="200")
  inspect(t1.exons[1].1, content="300")
  inspect(t1.exons[2].0, content="400")
  inspect(t1.exons[2].1, content="440")
  assert_true((t1.coverage - 4.0).abs() < 0.01)
  // T2: E1-E2-E4
  let t2 = ts[1]
  inspect(t2.exons.length(), content="3")
  inspect(t2.exons[2].0, content="600")
  inspect(t2.exons[2].1, content="640")
  assert_true((t2.coverage - 2.0).abs() < 0.01)
  // T3: E1 单外显子(剩余深度 20+12*20/90-4-2 = 20+2.667-6 = 16.67)
  let t3 = ts[2]
  inspect(t3.exons.length(), content="1")
  inspect(t3.exons[0].0, content="10")
  inspect(t3.exons[0].1, content="100")
  assert_true(t3.coverage > 15.0)
  assert_true(t3.coverage < 18.0)
}

///|
test "transcript_assembly: 表达量与 FPKM" {
  let ts = @src.ta_assemble_transcripts(
    "chr1",
    ta_sam_three_isoforms(),
    strand="+",
  )
  let t1 = ts[0]
  // avg_read_len = 2880/38 = 75.789(含 90bp 未剪接读段);
  // read_count = cov * eff_len / avg = 4*230/75.789 = 12.138
  inspect(t1.effective_length, content="230")
  assert_true((t1.read_count - 12.1379).abs() < 0.01)
  // n_reads = 38; fpkm 公式一致性(用库的 read_count 代入期望公式)
  let expect_fpkm = t1.read_count * 1.0e9 / (38.0 * 230.0)
  assert_true((t1.fpkm - expect_fpkm).abs() < 1.0)
  assert_true(t1.fpkm > 0.0)
}

///|
test "transcript_assembly: 汇总与 GTF" {
  let sam = ta_sam_three_isoforms()
  let ts = @src.ta_assemble_transcripts("chr1", sam, strand="+")
  let stats = @src.ta_assembly_stats("chr1", sam, ts)
  inspect(stats.n_reads, content="38")
  inspect(stats.n_junctions, content="3")
  inspect(stats.n_transcripts, content="3")
  let line = @src.ta_assembly_stats_line(stats)
  assert_eq(line.contains("junctions=3"), true)
  assert_eq(line.contains("transcripts=3"), true)
  let gtf = @src.ta_transcripts_to_gtf(ts)
  assert_eq(gtf.contains("##gff-version 2.2"), true)
  assert_eq(gtf.contains("transcript_assembly\ttranscript\t11\t440"), true)
  assert_eq(gtf.contains("exon_number \"1\""), true)
  assert_eq(gtf.contains("transcript_id \"TCONS_00000001\""), true)
  // GTF 1-based: E1 [10,100) -> 11..100
  assert_eq(gtf.contains("\texon\t11\t100\t"), true)
}

///|
test "transcript_assembly: JSON 导出" {
  let ts = @src.ta_assemble_transcripts(
    "chr1",
    ta_sam_three_isoforms(),
    strand="+",
  )
  let json = @src.ta_transcripts_json(ts)
  assert_eq(json.contains("\"TCONS_00000001\""), true)
  assert_eq(json.contains("\"fpkm\""), true)
  assert_eq(json.contains("\"exons\""), true)
}

///|
test "transcript_assembly: 无 junction 单外显子" {
  // 3 条读段: 平均深度 (80*2+80)/100 = 2.4 >= min_flow 2.0
  let sam = ta_sam_header() +
    "r1\t0\tchr1\t101\t60\t80M\t*\t0\t0\tACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACG\t*\n" +
    "r2\t0\tchr1\t121\t60\t80M\t*\t0\t0\tACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACG\t*\n" +
    "r3\t0\tchr1\t101\t60\t80M\t*\t0\t0\tACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACG\t*\n"
  let ts = @src.ta_assemble_transcripts("chr1", sam, strand="+")
  inspect(ts.length(), content="1")
  inspect(ts[0].exons.length(), content="1")
  inspect(ts[0].exons[0].0, content="100")
  inspect(ts[0].exons[0].1, content="200")
  assert_true((ts[0].coverage - 2.4).abs() < 0.01)
}

///|
test "transcript_assembly: 空输入" {
  let ts = @src.ta_assemble_transcripts("chr1", ta_sam_header(), strand="+")
  inspect(ts.length(), content="0")
  let stats = @src.ta_assembly_stats("chr1", ta_sam_header(), ts)
  inspect(stats.n_reads, content="0")
  inspect(stats.n_junctions, content="0")
}
