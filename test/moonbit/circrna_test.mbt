// ===========================================================================
// circrna_test.mbt — find_circ 风格环状 RNA 检测测试.
// ===========================================================================

///|
/// SAM 头.
fn cr_sam_header() -> String {
  "@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:10000\n"
}

///|
/// 构造一条 BSJ split-read:
/// 主定位 pos(1-based) 60M flag, SA = "chr1,sa_pos,+,60M,60,0".
fn cr_bsj_read(
  qname : String,
  main_pos : Int32,
  sa_pos : Int32,
  flag~ : Int32 = 0,
  sa_strand~ : String = "+",
) -> String {
  "\{qname}\t\{flag}\tchr1\t\{main_pos}\t60\t60M\t*\t0\t0\tACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTAC\t*\tSA:Z:chr1,\{sa_pos},\{sa_strand},60M,60,0\n"
}

///|
test "circrna: BSJ split-read 检测与聚合" {
  // circRNA [1000, 2000): 主 pos 1941(60M -> [1940,2000)),
  // SA pos 1001(60M -> [1000,1060)), SA 整体在主上游
  let sam = cr_sam_header() +
    cr_bsj_read("r1", 1941, 1001) +
    cr_bsj_read("r2", 1941, 1001) +
    cr_bsj_read("r3", 1941, 1001)
  let records = @src.circrna_detect(sam)
  inspect(records.length(), content="1")
  let r = records[0]
  inspect(r.chrom, content="chr1")
  inspect(r.start, content="1000")
  inspect(r.end_, content="2000")
  inspect(r.length, content="1000")
  inspect(r.support, content="3")
  inspect(r.strand, content="+")
  inspect(r.name, content="circ_chr1_1000_2000")
}

///|
test "circrna: 线性融合/剪接读段不计" {
  // SA 在主定位下游(正向拼接) = 线性, 不构成 BSJ
  let sam = cr_sam_header() +
    cr_bsj_read("f1", 501, 601) +
    cr_bsj_read("f2", 501, 601) +
    cr_bsj_read("f3", 501, 601)
  let records = @src.circrna_detect(sam)
  inspect(records.length(), content="0")
}

///|
test "circrna: 最小环长与最小支持过滤" {
  // 微小环(长 40 < 50): 主 pos 1031 10M -> [1030,1040),
  // SA pos 1001 10M -> [1000,1010), circ = [1000,1040) 长 40
  let tiny_sam = cr_sam_header() +
    "t1\t0\tchr1\t1031\t60\t10M\t*\t0\t0\tACGTACGTAC\t*\tSA:Z:chr1,1001,+,10M,60,0\n" +
    "t2\t0\tchr1\t1031\t60\t10M\t*\t0\t0\tACGTACGTAC\t*\tSA:Z:chr1,1001,+,10M,60,0\n"
  let tiny_records = @src.circrna_detect(tiny_sam)
  inspect(tiny_records.length(), content="0")
  // 单读支持(support 1 < 2)
  let single = cr_sam_header() + cr_bsj_read("s1", 3141, 3001)
  let single_records = @src.circrna_detect(single)
  inspect(single_records.length(), content="0")
  // 放宽参数后均通过
  let params : @src.CircRnaParams = { min_support: 1, min_length: 30, }
  let relaxed = @src.circrna_detect(tiny_sam, params~)
  inspect(relaxed.length(), content="1")
  inspect(relaxed[0].length, content="40")
}

///|
test "circrna: 链方向与跨染色体过滤" {
  // 反链 circRNA(主 flag 16 + SA "-")
  let sam = cr_sam_header() +
    cr_bsj_read("m1", 1941, 1001, flag=16, sa_strand="-") +
    cr_bsj_read("m2", 1941, 1001, flag=16, sa_strand="-")
  let records = @src.circrna_detect(sam)
  inspect(records.length(), content="1")
  inspect(records[0].strand, content="-")
  // 链不一致 -> "."
  let sam2 = cr_sam_header() +
    cr_bsj_read("x1", 1941, 1001, flag=0, sa_strand="-") +
    cr_bsj_read("x2", 1941, 1001, flag=0, sa_strand="-")
  let records2 = @src.circrna_detect(sam2)
  inspect(records2[0].strand, content=".")
  // 跨染色体 SA 不计
  let sam3 = cr_sam_header() +
    "y1\t0\tchr1\t1941\t60\t60M\t*\t0\t0\tACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTAC\t*\tSA:Z:chr2,1001,+,60M,60,0\n" +
    "y2\t0\tchr1\t1941\t60\t60M\t*\t0\t0\tACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTAC\t*\tSA:Z:chr2,1001,+,60M,60,0\n"
  let records3 = @src.circrna_detect(sam3)
  inspect(records3.length(), content="0")
}

///|
test "circrna: 多环与排序" {
  // 两个独立 circRNA: [1000,2000) support 3 与 [3000,3200) support 2
  let sam = cr_sam_header() +
    cr_bsj_read("a1", 1941, 1001) +
    cr_bsj_read("a2", 1941, 1001) +
    cr_bsj_read("a3", 1941, 1001) +
    cr_bsj_read("b1", 3141, 3001) +
    cr_bsj_read("b2", 3141, 3001)
  let records = @src.circrna_detect(sam)
  inspect(records.length(), content="2")
  // 按 start 升序
  inspect(records[0].start, content="1000")
  inspect(records[0].support, content="3")
  inspect(records[1].start, content="3000")
  inspect(records[1].support, content="2")
  inspect(records[1].length, content="200")
}

///|
test "circrna: secondary/supplementary 双计过滤" {
  // secondary(0x100) 与 supplementary(0x800) 记录不参与(避免双计)
  let sam = cr_sam_header() +
    cr_bsj_read("p1", 1941, 1001) +
    cr_bsj_read("p2", 1941, 1001) +
    cr_bsj_read("p3", 1941, 1001, flag=256) +
    cr_bsj_read("p4", 1941, 1001, flag=2048)
  let records = @src.circrna_detect(sam)
  inspect(records.length(), content="1")
  inspect(records[0].support, content="2")
}

///|
test "circrna: TSV/JSON/摘要导出" {
  let sam = cr_sam_header() +
    cr_bsj_read("a1", 1941, 1001) +
    cr_bsj_read("a2", 1941, 1001) +
    cr_bsj_read("b1", 3141, 3001) +
    cr_bsj_read("b2", 3141, 3001)
  let records = @src.circrna_detect(sam)
  let tsv = @src.circrna_tsv(records)
  assert_eq(
    tsv.contains("chrom\tstart\tend\tname\tsupport\tstrand\tlength"),
    true,
  )
  assert_eq(
    tsv.contains("chr1\t1000\t2000\tcirc_chr1_1000_2000\t2\t+\t1000"),
    true,
  )
  assert_eq(tsv.contains("chr1\t3000\t3200\tcirc_chr1_3000_3200\t2"), true)
  let json = @src.circrna_json(records)
  assert_eq(json.contains("\"name\":\"circ_chr1_1000_2000\""), true)
  assert_eq(json.contains("\"support\""), true)
  assert_eq(json.contains("\"strand\""), true)
  let line = @src.circrna_summary_line(records)
  assert_eq(line.contains("records=2"), true)
  assert_eq(line.contains("total_bsj_support=4"), true)
  let s = records[0].summary()
  assert_eq(s.contains("support=2"), true)
  assert_eq(s.contains("len=1000"), true)
}

///|
test "circrna: 空输入与无 SA 读段" {
  let records = @src.circrna_detect(cr_sam_header())
  inspect(records.length(), content="0")
  // 无 SA 标签的普通读段
  let sam = cr_sam_header() +
    "n1\t0\tchr1\t100\t60\t60M\t*\t0\t0\tACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTAC\t*\n"
  let records2 = @src.circrna_detect(sam)
  inspect(records2.length(), content="0")
}
