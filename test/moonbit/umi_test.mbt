// ===========================================================================
// umi_test.mbt — UMI 提取 / 纠错 / 去重 / 共识测试
// ===========================================================================

///|
test "umi parse header underscore format" {
  let h = @src.umi_parse_header(
    "@read1_ACGTACGT some description",
    @src.Underscore,
  )
  assert_eq(h.name, "read1")
  assert_eq(h.umi, "ACGTACGT")
  assert_eq(@src.umi_extract_umi("read1_acgt", @src.Underscore), "ACGT")
  assert_eq(@src.umi_strip_header("read1_acgt", @src.Underscore), "read1")
}

///|
test "umi parse header keeps internal separators" {
  let h = @src.umi_parse_header("@sample_1_flowcell_ACGT", @src.Underscore)
  assert_eq(h.name, "sample_1_flowcell")
  assert_eq(h.umi, "ACGT")
}

///|
test "umi parse header hash illumina pipe formats" {
  assert_eq(@src.umi_extract_umi("read1#TTTT", @src.Hash), "TTTT")
  assert_eq(
    @src.umi_extract_umi("M1:R1:FC1:1:100:200:300:GGGG", @src.Illumina),
    "GGGG",
  )
  assert_eq(@src.umi_strip_header("read1|CCCC", @src.Pipe), "read1")
}

///|
test "umi parse header rejects invalid input" {
  let mut raised = false
  try {
    let h = @src.umi_parse_header("read1", @src.Underscore)
    ignore(h)
  } catch {
    _ => raised = true
  }
  assert_true(raised)
  raised = false
  try {
    let h = @src.umi_parse_header("read1_ACXG", @src.Underscore)
    ignore(h)
  } catch {
    _ => raised = true
  }
  assert_true(raised)
}

///|
test "umi quality ascii roundtrip" {
  let q = @src.umi_qualities_from_ascii("IIII#!")
  assert_eq(q.length(), 6)
  assert_eq(q[0], 40)
  assert_eq(q[4], 2)
  assert_eq(@src.umi_qualities_to_ascii(q), "IIII#!")
}

///|
test "umi clip read removes 5prime bases" {
  let (seq, qual) = @src.umi_clip_read(
    "ACGTACGTAAAA",
    [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12],
    8,
  )
  assert_eq(seq, "AAAA")
  assert_eq(qual, [9, 10, 11, 12])
}

///|
test "umi clip read validates range" {
  let mut raised = false
  try {
    let (seq, qual) = @src.umi_clip_read("ACGT", [1, 2, 3, 4], 5)
    ignore(seq)
    ignore(qual)
  } catch {
    _ => raised = true
  }
  assert_true(raised)
}

///|
test "umi consensus prefers high quality base" {
  // 第 2 位: 两条低质量 T vs 一条高质量 C
  let seqs = ["AAT", "AAT", "AAC"]
  let quals = [[40, 1, 1], [40, 1, 1], [40, 40, 40]]
  assert_eq(@src.umi_consensus(seqs, quals), "AAC")
}

///|
test "umi consensus ties choose fixed order" {
  let seqs = ["AC", "CA"]
  let quals = [[30, 30], [30, 30]]
  // 平票时按 A/C/G/T/N 顺序
  assert_eq(@src.umi_consensus(seqs, quals), "AA")
}

///|
test "umi edit distance hamming" {
  assert_eq(@src.umi_edit_distance("ACGT", "ACGT"), 0)
  assert_eq(@src.umi_edit_distance("ACGT", "ACGA"), 1)
  assert_eq(@src.umi_edit_distance("ACGT", "AC"), 99)
}

///|
test "umi dedup groups reads by umi and locus" {
  let reads : Array[@src.UmiRead] = [
    @src.umi_read(
      "r1",
      "AAAA",
      @src.umi_qualities_from_ascii("IIII"),
      "ACGT",
      "chr1",
      100,
    ),
    @src.umi_read(
      "r2",
      "AAAT",
      @src.umi_qualities_from_ascii("IIII"),
      "ACGT",
      "chr1",
      100,
    ),
    @src.umi_read(
      "r3",
      "CCCC",
      @src.umi_qualities_from_ascii("IIII"),
      "TTTT",
      "chr1",
      100,
    ),
    @src.umi_read(
      "r4",
      "GGGG",
      @src.umi_qualities_from_ascii("IIII"),
      "ACGT",
      "chr1",
      200,
    ),
  ]
  let families = @src.umi_dedup(reads)
  assert_eq(families.length(), 3)
  // 按 chrom/pos/strand/umi 排序: pos100 ACGT, pos100 TTTT, pos200 ACGT
  assert_eq(families[0].umi, "ACGT")
  assert_eq(families[0].pos, 100)
  assert_eq(families[0].size, 2)
  assert_eq(families[0].consensus, "AAAA")
  assert_eq(families[1].umi, "TTTT")
  assert_eq(families[1].size, 1)
  assert_eq(families[2].pos, 200)
}

///|
test "umi dedup applies minimum family size" {
  let reads : Array[@src.UmiRead] = [
    @src.umi_read(
      "r1",
      "AAAA",
      @src.umi_qualities_from_ascii("IIII"),
      "ACGT",
      "chr1",
      100,
    ),
    @src.umi_read(
      "r2",
      "CCCC",
      @src.umi_qualities_from_ascii("IIII"),
      "TTTT",
      "chr1",
      100,
    ),
  ]
  let families = @src.umi_dedup(reads, min_family_size=2)
  assert_eq(families.length(), 0)
}

///|
test "umi directional correction merges adjacent umi" {
  let families : Array[@src.UmiFamily] = [
    {
      umi: "ACGT",
      chrom: "chr1",
      pos: 100,
      strand: 1,
      size: 10,
      consensus: "AAAA",
      mean_quality: 30.0,
    },
    {
      umi: "ACGA",
      chrom: "chr1",
      pos: 100,
      strand: 1,
      size: 2,
      consensus: "AAAT",
      mean_quality: 20.0,
    },
    {
      umi: "ACGA",
      chrom: "chr1",
      pos: 200,
      strand: 1,
      size: 5,
      consensus: "CCCC",
      mean_quality: 25.0,
    },
  ]
  let corrected = @src.umi_directional_correct(families)
  assert_eq(corrected.length(), 2)
  assert_eq(corrected[0].umi, "ACGT")
  assert_eq(corrected[0].size, 12)
  assert_true(corrected[0].mean_quality > 28.0)
  assert_eq(corrected[1].pos, 200)
}

///|
test "stats and json serialization" {
  let families : Array[@src.UmiFamily] = [
    {
      umi: "ACGT",
      chrom: "chr1",
      pos: 100,
      strand: 1,
      size: 3,
      consensus: "AAAA",
      mean_quality: 30.0,
    },
    {
      umi: "TTTT",
      chrom: "chr1",
      pos: 100,
      strand: 1,
      size: 1,
      consensus: "CCCC",
      mean_quality: 20.0,
    },
  ]
  let stats = @src.umi_stats(families)
  assert_eq(stats.n_families, 2)
  assert_eq(stats.n_reads, 4)
  assert_eq(stats.singleton_families, 1)
  assert_eq(stats.max_family_size, 3)
  assert_eq(stats.mean_family_size, 2.0)

  let json = @src.umi_families_to_json(families)
  assert_true(json.contains("ACGT"))
}
