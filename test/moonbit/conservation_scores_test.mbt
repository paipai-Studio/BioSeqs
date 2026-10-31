// ===========================================================================
// conservation_scores_test.mbt — GERP/phastCons-lite 进化约束打分测试.
// ===========================================================================

///|
/// 5 序列 MSA: 列 0-7 与 16-19 全同(保守), 列 8-15 混合(可变).
fn cs_test_msa() -> Array[String] {
  [
    "ACGTACGT" + "TTTTGGGG" + "TGCA",
    "ACGTACGT" + "CCCCAAAA" + "TGCA",
    "ACGTACGT" + "TTTTGGGG" + "TGCA",
    "ACGTACGT" + "CCCCAAAA" + "TGCA",
    "ACGTACGT" + "ACGTACGT" + "TGCA",
  ]
}

///|
test "conservation: GERP-lite 保守/可变列区分" {
  let gerp = @src.conservation_gerp(cs_test_msa(), tree_length=2.0)
  inspect(gerp.length(), content="20")
  // 全同列(熵 0) -> tree_length
  assert_eq(gerp[0], 2.0)
  assert_eq(gerp[4], 2.0)
  assert_eq(gerp[16], 2.0)
  assert_eq(gerp[19], 2.0)
  // 混合列(列 8: 2T+2C+1A -> 熵>0) -> 低于 tree_length
  assert_true(gerp[8] < 2.0)
  assert_true(gerp[8] > 0.0)
  // 更均匀的列(majority 更小)分数更低: 列 8 majority=2,
  // 列 12(G2A3 majority=3)分数高于列 8
  assert_true(gerp[12] > gerp[8])
}

///|
test "conservation: GERP-lite gap 处理与 min_effective" {
  // 近全 gap 列(仅 1 有效碱基 < 默认 min_effective=2) -> 0
  let msa = ["A-", "--", "AC"]
  let gerp = @src.conservation_gerp(msa, tree_length=1.0)
  inspect(gerp.length(), content="2")
  assert_eq(gerp[0], 1.0)
  assert_eq(gerp[1], 0.0)
  // 有效碱基数 < min_effective 的列 -> 0
  let msa2 = ["AC", "AC", "--"]
  let gerp2 = @src.conservation_gerp(msa2, tree_length=1.0, min_effective=3)
  inspect(gerp2[0], content="0")
  inspect(gerp2[1], content="0")
  // N 同 gap(列 1 仅剩 1 个有效碱基 < 2)
  let msa3 = ["AN", "AN", "AC"]
  let gerp3 = @src.conservation_gerp(msa3, tree_length=1.0)
  assert_eq(gerp3[0], 1.0)
  assert_eq(gerp3[1], 0.0)
}

///|
test "conservation: 保守后验区分" {
  let post = @src.conservation_posterior(cs_test_msa())
  inspect(post.length(), content="20")
  // 全同列(5/5 majority) -> 高后验
  assert_true(post[0] > 0.9)
  assert_true(post[16] > 0.9)
  // 混合列(majority 2/5) -> 低后验
  assert_true(post[8] < 0.5)
  assert_true(post[12] < 0.5)
}

///|
test "conservation: Viterbi 保守区段" {
  let msa = cs_test_msa()
  let post = @src.conservation_posterior(msa)
  let segs = @src.conservation_viterbi_segments(msa, posterior=Some(post))
  // 两个保守区段: [0,8) 与 [16,20) 附近
  inspect(segs.length(), content="2")
  inspect(segs[0].start, content="0")
  assert_true(segs[0].end_ >= 7)
  assert_true(segs[0].end_ <= 9)
  inspect(segs[1].start, content="16")
  inspect(segs[1].end_, content="20")
  // 区段平均后验为高
  assert_true(segs[0].mean_posterior > 0.9)
  assert_true(segs[1].mean_posterior > 0.9)
}

///|
test "conservation: 全可变 MSA 无保守区段" {
  let msa = ["ACGT", "CGTA", "GTAC", "TACG", "ACGT"]
  // 每列 majority 2/5, 后验低 -> Viterbi 全 neutral
  let segs = @src.conservation_viterbi_segments(msa)
  inspect(segs.length(), content="0")
}

///|
test "conservation: GenomicScores 容器闭环" {
  let gs = @src.GScores::new("GERP-lite", "test")
  let msa = cs_test_msa()
  let gerp = @src.conservation_gerp(msa, tree_length=2.0)
  @src.conservation_to_gscores(gs, "chr1", gerp)
  gs.set_default_population("gerp_lite")
  // 逐位检索闭环(score_at 直接返回数值, 缺失时按插值策略)
  let v0 = gs.score_at("chr1", 1)
  assert_eq(v0, 2.0)
  // 区间均值(列 1-8, 全保守 -> 2.0)
  let r = gs.range_score("chr1", 1, 8)
  assert_true((r - 2.0).abs() < 1.0e-9)
}

///|
test "conservation: TSV/JSON/摘要导出" {
  let msa = cs_test_msa()
  let gerp = @src.conservation_gerp(msa, tree_length=2.0)
  let post = @src.conservation_posterior(msa)
  let segs = @src.conservation_viterbi_segments(msa, posterior=Some(post))
  let tsv = @src.conservation_tsv("chr1", gerp, post)
  assert_eq(tsv.contains("chrom\tpos\tgerp_lite\tposterior"), true)
  assert_eq(tsv.contains("chr1\t1\t2"), true)
  let json = @src.conservation_json(gerp, post, segs)
  assert_eq(json.contains("\"gerp_lite\""), true)
  assert_eq(json.contains("\"posterior\""), true)
  assert_eq(json.contains("\"conserved_segments\""), true)
  assert_eq(json.contains("\"mean_posterior\""), true)
  let line = @src.conservation_summary_line(gerp, segs)
  assert_eq(line.contains("cols=20"), true)
  assert_eq(line.contains("max_gerp=2"), true)
  assert_eq(line.contains("conserved_segments=2"), true)
}

///|
test "conservation: 空输入" {
  let empty_gerp = @src.conservation_gerp([])
  assert_eq(empty_gerp.length(), 0)
  let empty_post = @src.conservation_posterior([])
  assert_eq(empty_post.length(), 0)
  let empty_segs = @src.conservation_viterbi_segments([])
  assert_eq(empty_segs.length(), 0)
  let line = @src.conservation_summary_line([], [])
  assert_eq(line.contains("cols=0"), true)
}
