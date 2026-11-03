// ===========================================================================
// bam_postprocess_test.mbt — SAM->BAM 桥 / fixmate / markdup / BQSR 测试.
// ===========================================================================

///|
fn bp_sam_header() -> String {
  "@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:1000\n"
}

///|
test "bam_postprocess: sam_to_bam 转换" {
  let sam = bp_sam_header() +
    "r1\t0\tchr1\t101\t60\t8M\t=\t109\t8\tACGTACGT\tIIIIIIII\tNM:i:0\n" +
    "r2\t16\tchr1\t51\t40\t4M\t*\t0\t0\tTTTT\tJJJJ\n"
  let bam = @src.sam_to_bam(sam)
  assert_eq(bam.num_references(), 1)
  assert_eq(bam.get_reference_names()[0], "chr1")
  assert_eq(bam.records.length(), 2)
  let r1 = bam.records[0]
  assert_eq(r1.qname, "r1")
  assert_eq(r1.ref_id, 0)
  assert_eq(r1.pos, 100)
  assert_eq(r1.next_ref_id, 0)
  assert_eq(r1.next_pos, 108)
  assert_eq(r1.tlen, 8)
  assert_eq(r1.cigar.length(), 1)
  assert_eq(r1.cigar[0].length, 8)
  assert_eq(r1.seq, "ACGTACGT")
  assert_eq(r1.qual, "IIIIIIII")
  let has_nm = match r1.get_tag("NM") {
    Some(_) => true
    None => false
  }
  assert_eq(has_nm, true)
  let r2 = bam.records[1]
  assert_eq(r2.ref_id, 0)
  assert_eq(r2.pos, 50)
  assert_eq(r2.next_ref_id, -1)
  assert_eq(r2.flag, 16)
  assert_eq(r2.qual, "JJJJ")
}

///|
test "bam_postprocess: fixmate 填充 mate 坐标与 TLEN" {
  // qname 相邻的配对(输入 name-sorted), TLEN 故意为 0
  let sam = bp_sam_header() +
    "p1\t99\tchr1\t101\t60\t8M\t=\t201\t0\tACGTACGT\tIIIIIIII\n" +
    "p1\t147\tchr1\t201\t60\t8M\t=\t101\t0\tTTTTGGGG\tJJJJJJJJ\n"
  let bam = @src.sam_to_bam(sam)
  let fixed = @src.bam_fixmate(bam)
  assert_eq(fixed.records.length(), 2)
  let first = fixed.records[0]
  let second = fixed.records[1]
  // 外坐标: 0-based [100, 208) -> span 108, 左端为正
  assert_eq(first.tlen, 108)
  assert_eq(second.tlen, -108)
  // mate 坐标互填(0-based)
  assert_eq(first.next_pos, 200)
  assert_eq(second.next_pos, 100)
  // flag 镜像: 99(0x63) 含 0x20(mate reverse), 147(0x93) 不含
  assert_eq((first.flag & 0x20) != 0, true)
  assert_eq((second.flag & 0x20) == 0, true)
  assert_eq((first.flag & 0x8) == 0, true)
}

///|
test "bam_postprocess: fixmate 非配对读段 tlen 置 0" {
  let sam = bp_sam_header() +
    "s1\t0\tchr1\t101\t60\t8M\t*\t0\t99\tACGTACGT\tIIIIIIII\n"
  let bam = @src.sam_to_bam(sam)
  let fixed = @src.bam_fixmate(bam)
  assert_eq(fixed.records[0].tlen, 0)
}

///|
test "bam_postprocess: 配对感知 markdup 与光学重复" {
  // 两对读段同 5'-coordinate key(mapq 120 vs 60), 光学坐标距离 3/2;
  // 单端 s1/s2 同 key(mapq 25 vs 20)
  let sam = bp_sam_header() +
    "dup:1:1:100:200\t99\tchr1\t101\t60\t8M\t=\t201\t108\tACGTACGT\tIIIIIIII\n" +
    "dup:1:1:100:200\t147\tchr1\t201\t60\t8M\t=\t101\t-108\tTTTTGGGG\tJJJJJJJJ\n" +
    "dup:1:1:103:202\t99\tchr1\t101\t30\t8M\t=\t201\t108\tACGTACGT\tIIIIIIII\n" +
    "dup:1:1:103:202\t147\tchr1\t201\t30\t8M\t=\t101\t-108\tTTTTGGGG\tJJJJJJJJ\n" +
    "s1\t0\tchr1\t301\t25\t6M\t*\t0\t0\tACGTAC\tIIIIII\n" +
    "s2\t0\tchr1\t301\t20\t6M\t*\t0\t0\tACGTAC\tIIIIII\n"
  let bam = @src.sam_to_bam(sam)
  let (out, stats) = @src.bam_markdup_paired(bam)
  // 高 mapq 对保留, 低 mapq 对两条打标
  assert_eq((out.records[0].flag & 0x400) == 0, true)
  assert_eq((out.records[1].flag & 0x400) == 0, true)
  assert_eq((out.records[2].flag & 0x400) != 0, true)
  assert_eq((out.records[3].flag & 0x400) != 0, true)
  // 单端 s2 打标
  assert_eq((out.records[4].flag & 0x400) == 0, true)
  assert_eq((out.records[5].flag & 0x400) != 0, true)
  // 统计
  assert_eq(stats.n_pair_groups, 1)
  assert_eq(stats.n_dup_pairs, 1)
  assert_eq(stats.n_optical, 1)
  assert_eq(stats.n_single_groups, 1)
  assert_eq(stats.n_dup_single, 1)
  assert_eq(stats.dup_rate, 0.5)
  let line = @src.bam_markdup_stats_line(stats)
  assert_eq(line.contains("dup_pairs=1"), true)
  assert_eq(line.contains("optical=1"), true)
}

///|
test "bam_postprocess: markdup mapq 平局取 qname 最小" {
  let sam = bp_sam_header() +
    "bb\t99\tchr1\t101\t60\t8M\t=\t201\t108\tACGTACGT\tIIIIIIII\n" +
    "bb\t147\tchr1\t201\t60\t8M\t=\t101\t-108\tTTTTGGGG\tJJJJJJJJ\n" +
    "aa\t99\tchr1\t101\t60\t8M\t=\t201\t108\tACGTACGT\tIIIIIIII\n" +
    "aa\t147\tchr1\t201\t60\t8M\t=\t101\t-108\tTTTTGGGG\tJJJJJJJJ\n"
  let bam = @src.sam_to_bam(sam)
  let (out, stats) = @src.bam_markdup_paired(bam)
  // aa(字典序小)保留
  let mut aa_marked = 0
  let mut bb_marked = 0
  for r in out.records {
    if r.qname == "aa" && (r.flag & 0x400) != 0 {
      aa_marked = aa_marked + 1
    }
    if r.qname == "bb" && (r.flag & 0x400) != 0 {
      bb_marked = bb_marked + 1
    }
  }
  assert_eq(aa_marked, 0)
  assert_eq(bb_marked, 2)
  assert_eq(stats.n_dup_pairs, 1)
  assert_eq(stats.n_optical, 0)
}

///|
test "bam_postprocess: BQSR 错配重估" {
  // 40bp 参考 "ACGT" 循环; A: 10 条全 match(q30); B: 1 条 i=3 错配
  let sb = StringBuilder()
  for _i = 0; _i < 10; _i = _i + 1 {
    sb.write_string("ACGT")
  }
  let reference = sb.to_string()
  let sam_sb = StringBuilder()
  sam_sb.write_string("@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:1000\n")
  for i = 0; i < 10; i = i + 1 {
    sam_sb.write_string(
      "a\{i}\t0\tchr1\t1\t60\t8M\t*\t0\t0\tACGTACGT\t????????\n",
    )
  }
  sam_sb.write_string("bb\t0\tchr1\t1\t60\t8M\t*\t0\t0\tACGAACGT\t????????\n")
  let bam = @src.sam_to_bam(sam_sb.to_string())
  let (out, report) = @src.bam_bqsr(bam, reference)
  // 11 条读段 x 8 位 = 88 观察数, 其中 1 个错配(B 的 i=3)
  assert_eq(report.n_observations, 88)
  assert_eq(report.n_errors, 1)
  assert_eq(report.mean_qual_before, 30.0)
  // A 读段: 无错组按经验重估(小样本下经验上界低于标称, 属预期行为)
  let a3q = out.records[0].qual.unsafe_get(3).to_int() - 33
  // B 读段: i=3 的 ctx "GAA" 有错 -> 下调到 min_qual 2(ASCII 35 = '#'),
  // 且显著低于无错同位 A i=3
  let mut b_rec : @src.BamRecord? = None
  for r in out.records {
    if r.qname == "bb" {
      b_rec = Some(r)
    }
  }
  match b_rec {
    Some(r) => {
      assert_eq(r.qual.unsafe_get(3).unsafe_to_char(), '#')
      let b3q = r.qual.unsafe_get(3).to_int() - 33
      assert_true(b3q < a3q)
    }
    None => assert_true(false)
  }
  let line = @src.bam_bqsr_report_line(report)
  assert_eq(line.contains("observations=88"), true)
  assert_eq(line.contains("errors=1"), true)
  let json = @src.json_to_string(@src.bam_bqsr_value(report))
  assert_eq(json.contains("n_errors"), true)
}

///|
test "bam_postprocess: BQSR known-sites 豁免" {
  // 同上但注入 known-sites(0-based [3,4)): 覆盖该参考坐标的全部读段位
  // 均不参与训练 -> 错配不计 err, 该 covariate 无表项时应用保持原值
  let sb = StringBuilder()
  for _i = 0; _i < 10; _i = _i + 1 {
    sb.write_string("ACGT")
  }
  let reference = sb.to_string()
  let sam_sb = StringBuilder()
  sam_sb.write_string("@HD\tVN:1.6\n@SQ\tSN:chr1\tLN:1000\n")
  for i = 0; i < 10; i = i + 1 {
    sam_sb.write_string(
      "a\{i}\t0\tchr1\t1\t60\t8M\t*\t0\t0\tACGTACGT\t????????\n",
    )
  }
  sam_sb.write_string("bb\t0\tchr1\t1\t60\t8M\t*\t0\t0\tACGAACGT\t????????\n")
  let bam = @src.sam_to_bam(sam_sb.to_string())
  let (out, report) = @src.bam_bqsr(bam, reference, known_sites=[(0, 3, 4)])
  // 88 位中每条读段的 ref_pos=3 位被豁免(11 条) -> 77 观察数, 0 错配
  assert_eq(report.n_observations, 77)
  assert_eq(report.n_errors, 0)
  // B 读段 i=3 的 covariate 无表项 -> 质量保持 30('?')
  let mut b_rec : @src.BamRecord? = None
  for r in out.records {
    if r.qname == "bb" {
      b_rec = Some(r)
    }
  }
  match b_rec {
    Some(r) => assert_eq(r.qual.unsafe_get(3).unsafe_to_char(), '?')
    None => assert_true(false)
  }
}

///|
test "bam_postprocess: BQSR 空读段与空 CIGAR 保留原样" {
  let sam = bp_sam_header() + "q1\t4\t*\t0\t0\t*\t*\t0\t0\tACGT\t****\n"
  let bam = @src.sam_to_bam(sam)
  let reference = "ACGTACGT"
  let (out, report) = @src.bam_bqsr(bam, reference)
  assert_eq(out.records.length(), 1)
  assert_eq(report.n_observations, 0)
  assert_eq(report.mean_qual_before, 0.0)
}
