// ===========================================================================
// methylation_caller_test.mbt — Bismark 风格甲基化调用器测试.
// ===========================================================================

///|
/// 20bp 参考, 手工标注位点:
/// 正链: 1=CpG 4=CHG 7=CHH 11=CpG 15=CHH 18=CHH (0-based)
/// 反链 G: 2=CpG 6=CHG 12=CpG (0-based)
fn mc_reference() -> String {
  "ACGTCAGCAATCGTTCAACA"
}

///|
fn mc_sam_header() -> String {
  "@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:1000\n"
}

///|
/// 正链读段: CpG(1)meth / CHG(4)unmeth / CHH(7)meth / CpG(11)unmeth /
/// CHH(15)unmeth / CHH(18)meth 的 C->T 转化序列.
fn mc_forward_read() -> String {
  "ACGTTAGCAATTGTTTAACA"
}

///|
/// 反链读段(正链方向存储): G(2)meth / G(6)unmeth / G(12)meth.
fn mc_reverse_read() -> String {
  "ACGTCAACAATCGTTCAACA"
}

///|
test "methylation: 正链 C->T 计数与上下文" {
  let sam = mc_sam_header() +
    "r1\t0\tchr1\t1\t60\t20M\t*\t0\t0\t\{mc_forward_read()}\t*\n"
  let calls = @src.mc_extract("chr1", mc_reference(), sam)
  inspect(calls.length(), content="6")
  // 按 (pos, strand) 排序: 全部 '+'
  let c0 = calls[0]
  inspect(c0.pos, content="2")
  inspect(c0.strand, content="+")
  inspect(c0.context, content="CpG")
  inspect(c0.meth, content="1")
  inspect(c0.unmeth, content="0")
  let c1 = calls[1]
  inspect(c1.pos, content="5")
  inspect(c1.context, content="CHG")
  inspect(c1.meth, content="0")
  inspect(c1.unmeth, content="1")
  let c2 = calls[2]
  inspect(c2.pos, content="8")
  inspect(c2.context, content="CHH")
  inspect(c2.meth, content="1")
  let c3 = calls[3]
  inspect(c3.pos, content="12")
  inspect(c3.context, content="CpG")
  inspect(c3.unmeth, content="1")
  let c4 = calls[4]
  inspect(c4.pos, content="16")
  inspect(c4.context, content="CHH")
  let c5 = calls[5]
  inspect(c5.pos, content="19")
  inspect(c5.context, content="CHH")
  inspect(c5.meth, content="1")
}

///|
test "methylation: 反链 G->A 计数" {
  let sam = mc_sam_header() +
    "r1\t16\tchr1\t1\t60\t20M\t*\t0\t0\t\{mc_reverse_read()}\t*\n"
  let calls = @src.mc_extract("chr1", mc_reference(), sam)
  inspect(calls.length(), content="3")
  let c0 = calls[0]
  inspect(c0.pos, content="3")
  inspect(c0.strand, content="-")
  inspect(c0.context, content="CpG")
  inspect(c0.meth, content="1")
  inspect(c0.unmeth, content="0")
  let c1 = calls[1]
  inspect(c1.pos, content="7")
  inspect(c1.context, content="CHG")
  inspect(c1.meth, content="0")
  inspect(c1.unmeth, content="1")
  let c2 = calls[2]
  inspect(c2.pos, content="13")
  inspect(c2.context, content="CpG")
  inspect(c2.meth, content="1")
}

///|
test "methylation: 多读段聚合" {
  let sam = mc_sam_header() +
    "a\t0\tchr1\t1\t60\t20M\t*\t0\t0\t\{mc_forward_read()}\t*\n" +
    "b\t0\tchr1\t1\t60\t20M\t*\t0\t0\t\{mc_forward_read()}\t*\n" +
    "c\t0\tchr1\t1\t60\t20M\t*\t0\t0\t\{mc_forward_read()}\t*\n"
  let calls = @src.mc_extract("chr1", mc_reference(), sam)
  inspect(calls.length(), content="6")
  // pos 2 CpG: 3 条读段均 methylated
  inspect(calls[0].meth, content="3")
  inspect(calls[0].unmeth, content="0")
  // pos 12 CpG: 3 条均 unmethylated
  inspect(calls[3].unmeth, content="3")
}

///|
test "methylation: 正反链同位独立聚合" {
  // 正链 C 位点(0-based 1, 1-based 2)与反链 G 位点(0-based 2,
  // 1-based 3)属不同 record
  let sam = mc_sam_header() +
    "f\t0\tchr1\t1\t60\t20M\t*\t0\t0\t\{mc_forward_read()}\t*\n" +
    "r\t16\tchr1\t1\t60\t20M\t*\t0\t0\t\{mc_reverse_read()}\t*\n"
  let calls = @src.mc_extract("chr1", mc_reference(), sam)
  // 6 正链 + 3 反链
  inspect(calls.length(), content="9")
  let mut plus = 0
  let mut minus = 0
  for c in calls {
    if c.strand == "+" {
      plus = plus + 1
    } else {
      minus = minus + 1
    }
  }
  inspect(plus, content="6")
  inspect(minus, content="3")
}

///|
test "methylation: CIGAR indel 读段" {
  // 读段 pos 1: 4M1I15M — 插入 1 碱基不推进参考
  // 参考 C 位点 0-based 1 在 4M 内; 0-based 4 在 4M 末;
  // 插入后 15M 覆盖 0-based 4..19
  // 构造: 前 4 碱基 ACGT + 插入 G + 后 15 碱基为参考 4..19 的转化
  // 参考[4..19] = CAGCAATCGTTCAACA
  // 转化: C(4) unmeth->T, C(7) meth->C, C(11) unmeth->T,
  //       C(15) unmeth->T, C(18) meth->C
  let tail = "TAGCAATTGTTTAACA"
  let read = "ACGT" + "G" + tail
  let sam = mc_sam_header() +
    "r1\t0\tchr1\t1\t60\t4M1I15M\t*\t0\t0\t\{read}\t*\n"
  let calls = @src.mc_extract("chr1", mc_reference(), sam)
  // C 位点: 1(4M 内), 4, 7, 11, 15, 18(15M 内)
  inspect(calls.length(), content="6")
  // 0-based 4 CHG: 读 T(unmeth)
  let c4 = calls[1]
  inspect(c4.pos, content="5")
  inspect(c4.unmeth, content="1")
  // 0-based 18 CHH meth(读 C)
  let c18 = calls[5]
  inspect(c18.pos, content="19")
  inspect(c18.meth, content="1")
}

///|
test "methylation: Bismark TSV 与 bsseq_parse_bismark 闭环" {
  let sam = mc_sam_header() +
    "f\t0\tchr1\t1\t60\t20M\t*\t0\t0\t\{mc_forward_read()}\t*\n" +
    "r\t16\tchr1\t1\t60\t20M\t*\t0\t0\t\{mc_reverse_read()}\t*\n"
  let calls = @src.mc_extract("chr1", mc_reference(), sam)
  let tsv = @src.mc_bismark_tsv(calls)
  assert_eq(tsv.contains("chr1\t2\t+\t1\t0\tCpG"), true)
  assert_eq(tsv.contains("chr1\t3\t-\t1\t0\tCpG"), true)
  // 闭环: TSV -> bsseq_parse_bismark -> BSData
  let bsdata = @src.bsseq_parse_bismark(tsv)
  inspect(bsdata.records.length(), content="9")
  let mut cpg_records = 0
  for rec in bsdata.records {
    if rec.context == "CpG" {
      cpg_records = cpg_records + 1
    }
  }
  inspect(cpg_records, content="4")
}

///|
test "methylation: 衔接 bsseq BSRecord 统计链" {
  let sam = mc_sam_header() +
    "f\t0\tchr1\t1\t60\t20M\t*\t0\t0\t\{mc_forward_read()}\t*\n"
  let calls = @src.mc_extract("chr1", mc_reference(), sam)
  let records = @src.mc_to_bsrecords(calls)
  inspect(records.length(), content="6")
  // pos 2 CpG: meth=1 unmeth=0 -> rate 1.0
  let r0 = records[0]
  inspect(r0.pos, content="2")
  assert_eq(r0.rate(), 1.0)
  inspect(r0.coverage(), content="1")
  // pos 12 CpG: unmeth -> rate 0
  assert_eq(records[3].rate(), 0.0)
}

///|
test "methylation: 汇总与 JSON" {
  let sam = mc_sam_header() +
    "f\t0\tchr1\t1\t60\t20M\t*\t0\t0\t\{mc_forward_read()}\t*\n" +
    "r\t16\tchr1\t1\t60\t20M\t*\t0\t0\t\{mc_reverse_read()}\t*\n"
  let calls = @src.mc_extract("chr1", mc_reference(), sam)
  let summary = @src.mc_summary(sam, calls)
  inspect(summary.n_reads, content="2")
  inspect(summary.n_cpg_sites, content="4")
  inspect(summary.n_chg_sites, content="2")
  inspect(summary.n_chh_sites, content="3")
  // CpG: 正链 2(pos2 meth, pos12 unmeth) + 反链 2(pos3 meth, pos13 meth)
  // rate = 3/4
  assert_eq(summary.cpg_methylation_rate, 0.75)
  let line = @src.mc_summary_line(summary)
  assert_eq(line.contains("cpg=4"), true)
  assert_eq(line.contains("cpg_rate=0.75"), true)
  let json = @src.mc_json(calls)
  assert_eq(json.contains("\"context\""), true)
  assert_eq(json.contains("\"unmeth\""), true)
}

///|
test "methylation: 空输入" {
  let calls = @src.mc_extract("chr1", mc_reference(), mc_sam_header())
  inspect(calls.length(), content="0")
  let summary = @src.mc_summary(mc_sam_header(), calls)
  inspect(summary.n_reads, content="0")
  inspect(summary.cpg_methylation_rate, content="0")
}
