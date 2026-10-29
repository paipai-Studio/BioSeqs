// ===========================================================================
// fusion_test.mbt — 融合基因 caller 测试(split-read 断点簇)
// ===========================================================================

///|
/// 构造 SAM 记录行(带 SA 标签).
fn fusion_sam_line(
  qname : String,
  chrom : String,
  pos : Int32,
  cigar : String,
  sa : String,
) -> String {
  let seq = ""
  let seq_filled = "ACGTACGTACGTACGTACGTACGTACGTACGT"
  ignore(seq)
  let qual = "IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII"
  qname +
  "\t0\t" +
  chrom +
  "\t" +
  pos.to_string() +
  "\t60\t" +
  cigar +
  "\t*\t0\t0\t" +
  seq_filled +
  "\t" +
  qual +
  "\tSA:Z:" +
  sa
}

///|
fn fusion_test_sam() -> String {
  let header = "@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:1000000\n@SQ\tSN:chr2\tLN:1000000\n"
  // read1-3: chr1:1000 -- chr2:5000 (同一融合断点, 支持数 3)
  // read4-5: chr1:2000 -- chr2:8000 (另一断点, 支持数 2)
  // read6: chr1:5000 -- chr2:9000 (单读, 低于 min_support)
  let lines : Array[String] = [
    header,
    fusion_sam_line(
      "split1", "chr1", 1000, "16M16S", "chr2,5000,+,16S16M,60,0;",
    ),
    fusion_sam_line(
      "split2", "chr1", 1000, "16M16S", "chr2,5000,+,16S16M,60,0;",
    ),
    fusion_sam_line(
      "split3", "chr1", 1005, "16M16S", "chr2,5005,+,16S16M,60,0;",
    ),
    fusion_sam_line(
      "split4", "chr1", 2000, "16M16S", "chr2,8000,+,16S16M,60,0;",
    ),
    fusion_sam_line(
      "split5", "chr1", 2000, "16M16S", "chr2,8000,+,16S16M,60,0;",
    ),
    fusion_sam_line(
      "split6", "chr1", 5000, "16M16S", "chr2,9000,+,16S16M,60,0;",
    ),
  ]
  let sb = StringBuilder()
  for l in lines {
    sb.write_string(l)
    sb.write_char('\n')
  }
  sb.to_string()
}

///|
test "fusion junctions from sam" {
  let sam = fusion_test_sam()
  let js = @src.fusion_junctions_from_sam(sam)
  // 6 条 split read 各产 1 个 junction
  assert_eq(js.length(), 6)
  // 有序对: chr1 < chr2, chrom_a = chr1
  for j in js {
    assert_eq(j.chrom_a, "chr1")
    assert_eq(j.chrom_b, "chr2")
  }
  // split1: 主 16M -> 断点 1000+16=1016; SA 16M -> 断点 5000+16=5016
  let j1 = js[0]
  assert_eq(j1.break_a, 1016)
  assert_eq(j1.break_b, 5016)
  assert_eq(j1.qname, "split1")
}

///|
test "fusion junctions no sa skipped" {
  let sam = "@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:1000\nread1\t0\tchr1\t100\t60\t32M\t*\t0\t0\tACGTACGTACGTACGTACGTACGTACGTACGT\tIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII\n"
  let js = @src.fusion_junctions_from_sam(sam)
  assert_eq(js.length(), 0)
}

///|
test "fusion call clusters and support" {
  let sam = fusion_test_sam()
  let js = @src.fusion_junctions_from_sam(sam)
  let candidates = @src.fusion_call(js, min_support=2)
  // 两个候选(支持 3 与 2); 单读 split6 被过滤
  assert_eq(candidates.length(), 2)
  // 支持数降序: 第一候选 support=3
  let top = candidates[0]
  assert_eq(top.support, 3)
  assert_eq(top.chrom_a, "chr1")
  assert_eq(top.chrom_b, "chr2")
  assert_true(top.is_interchromosomal)
  // 断点簇中位: (1016,1016,1021) -> 1016; (5016,5016,5021) -> 5016
  assert_eq(top.break_a, 1016)
  assert_eq(top.break_b, 5016)
  assert_eq(top.read_names.length(), 3)
  // 第二候选
  let second = candidates[1]
  assert_eq(second.support, 2)
  assert_eq(second.break_a, 2016)
  assert_eq(second.break_b, 8016)
}

///|
test "fusion call single read filtered" {
  let js : Array[@src.FusionJunction] = [
    {
      qname: "only1",
      chrom_a: "chr1",
      break_a: 1000,
      strand_a: 0,
      chrom_b: "chr2",
      break_b: 2000,
      strand_b: 0,
    },
  ]
  let candidates = @src.fusion_call(js, min_support=2)
  assert_eq(candidates.length(), 0)
}

///|
test "fusion call readthrough filtered" {
  // 同链同染色体近距离 -> read-through
  let js : Array[@src.FusionJunction] = []
  for i = 0; i < 5; i = i + 1 {
    js.push({
      qname: "rt" + i.to_string(),
      chrom_a: "chr1",
      break_a: 100000,
      strand_a: 0,
      chrom_b: "chr1",
      break_b: 150000,
      strand_b: 0,
    })
  }
  let candidates = @src.fusion_call(js, min_support=2, readthrough_max=200000)
  assert_eq(candidates.length(), 1)
  assert_eq(candidates[0].filter, "readthrough")
  assert_eq(candidates[0].is_interchromosomal, false)
  assert_eq(candidates[0].span, 50000)
}

///|
test "fusion call large span not readthrough" {
  // 同链但大跨距 -> PASS
  let js : Array[@src.FusionJunction] = []
  for i = 0; i < 5; i = i + 1 {
    js.push({
      qname: "lg" + i.to_string(),
      chrom_a: "chr1",
      break_a: 100000,
      strand_a: 0,
      chrom_b: "chr1",
      break_b: 900000,
      strand_b: 0,
    })
  }
  let candidates = @src.fusion_call(js, min_support=2, readthrough_max=200000)
  assert_eq(candidates.length(), 1)
  assert_eq(candidates[0].filter, "PASS")
}

///|
test "fusion call opposite strands not readthrough" {
  // 反式剪接方向(链相反) -> 不算 read-through
  let js : Array[@src.FusionJunction] = []
  for i = 0; i < 5; i = i + 1 {
    js.push({
      qname: "os" + i.to_string(),
      chrom_a: "chr1",
      break_a: 100000,
      strand_a: 0,
      chrom_b: "chr1",
      break_b: 150000,
      strand_b: 1,
    })
  }
  let candidates = @src.fusion_call(js, min_support=2, readthrough_max=200000)
  assert_eq(candidates.length(), 1)
  assert_eq(candidates[0].filter, "PASS")
}

///|
test "fusion call junction order normalized" {
  // 输入乱序 (chr2 先): 应归一为 chr1--chr2
  let js : Array[@src.FusionJunction] = [
    {
      qname: "n1",
      chrom_a: "chr2",
      break_a: 5000,
      strand_a: 0,
      chrom_b: "chr1",
      break_b: 1000,
      strand_b: 0,
    },
    {
      qname: "n2",
      chrom_a: "chr1",
      break_a: 1000,
      strand_a: 0,
      chrom_b: "chr2",
      break_b: 5000,
      strand_b: 0,
    },
  ]
  let candidates = @src.fusion_call(js, min_support=2)
  assert_eq(candidates.length(), 1)
  assert_eq(candidates[0].chrom_a, "chr1")
  assert_eq(candidates[0].fusion_name, "chr1--chr2")
}

///|
test "fusion write tsv format" {
  let sam = fusion_test_sam()
  let js = @src.fusion_junctions_from_sam(sam)
  let candidates = @src.fusion_call(js, min_support=2)
  let tsv = @src.fusion_write_tsv(candidates)
  let lines = tsv.split("\n").to_array()
  assert_eq(
    lines[0].to_owned(),
    "#FusionName\tJunctionReadCount\tChromA\tBreakA\tStrandA\tChromB\tBreakB\tStrandB\tInterchromosomal\tSpan\tFilter",
  )
  assert_true(lines[1].to_owned().contains("chr1--chr2"))
  assert_true(lines[1].to_owned().contains("\t3\t"))
  assert_true(lines[1].to_owned().contains("\tY\t"))
}

///|
test "fusion write bedpe and parse roundtrip" {
  let sam = fusion_test_sam()
  let js = @src.fusion_junctions_from_sam(sam)
  let candidates = @src.fusion_call(js, min_support=2)
  let bedpe = @src.fusion_write_bedpe(candidates)
  // 与 parse_bedpe 闭环(去掉注释行)
  let content = bedpe.split("\n").to_array()
  let filtered : Array[String] = []
  for l in content {
    let s = l.to_owned()
    if s.length() > 0 && !s.has_prefix("#") {
      filtered.push(s)
    }
  }
  let joined : String = filtered.join("\n") + "\n"
  let pairs = @src.parse_bedpe(joined)
  assert_eq(pairs.length(), 2)
  assert_eq(pairs[0].chrom1, "chr1")
  assert_eq(pairs[0].start1, 1015)
  assert_eq(pairs[0].end1, 1016)
  assert_eq(pairs[0].chrom2, "chr2")
}

///|
test "fusion candidates summary" {
  let sam = fusion_test_sam()
  let js = @src.fusion_junctions_from_sam(sam)
  let candidates = @src.fusion_call(js, min_support=2)
  let s = @src.fusion_candidates_summary(candidates)
  assert_true(s.contains("fusion"))
  assert_true(s.contains("2 candidates"))
  assert_true(s.contains("2 PASS"))
  assert_true(s.contains("2 interchromosomal"))
}

///|
test "fusion candidate summary line" {
  let c : @src.FusionCandidate = {
    chrom_a: "chr1",
    break_a: 1000,
    chrom_b: "chr2",
    break_b: 2000,
    strand_a: 0,
    strand_b: 1,
    support: 5,
    read_names: ["r1", "r2", "r3", "r4", "r5"],
    is_interchromosomal: true,
    fusion_name: "chr1--chr2",
    span: -1,
    filter: "PASS",
  }
  let s = c.summary()
  assert_true(s.contains("chr1--chr2"))
  assert_true(s.contains("support=5"))
  assert_true(s.contains("interchromosomal"))
}

///|
test "fusion call cluster tolerance" {
  // 断点相距 40bp(<=50 容差)应聚同簇
  let js : Array[@src.FusionJunction] = []
  for i = 0; i < 4; i = i + 1 {
    js.push({
      qname: "c" + i.to_string(),
      chrom_a: "chr1",
      break_a: 1000 + i * 20,
      strand_a: 0,
      chrom_b: "chr2",
      break_b: 5000 + i * 15,
      strand_b: 0,
    })
  }
  let candidates = @src.fusion_call(js, cluster_bp=50, min_support=2)
  assert_eq(candidates.length(), 1)
  assert_eq(candidates[0].support, 4)
  // 断点相距 200bp 应分裂
  let js2 : Array[@src.FusionJunction] = []
  for i = 0; i < 4; i = i + 1 {
    js2.push({
      qname: "d" + i.to_string(),
      chrom_a: "chr1",
      break_a: 1000 + i * 200,
      strand_a: 0,
      chrom_b: "chr2",
      break_b: 5000 + i * 200,
      strand_b: 0,
    })
  }
  let candidates2 = @src.fusion_call(js2, cluster_bp=50, min_support=1)
  assert_eq(candidates2.length(), 4)
}
