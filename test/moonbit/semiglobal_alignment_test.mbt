///|
/// Tests for exact semiglobal (glocal) alignment.

///|
test "sg_identical_sequences" {
  // 完全相同序列：无空位，得分 = len * match
  let r = @src.semiglobal_align(
    "ACGTACGT",
    "ACGTACGT",
    @src.SGScoring::new(2.0, -1.0, -2.0),
  )
  assert_eq(r.aligned1(), "ACGTACGT")
  assert_eq(r.aligned2(), "ACGTACGT")
  assert_eq(r.score(), 16.0)
  assert_eq(r.num_matches(), 8)
  assert_eq(r.num_mismatches(), 0)
  assert_eq(r.percent_identity(), 100.0)
}

///|
test "sg_free_terminal_gaps_core_property" {
  // 半全局核心性质：首尾自由空位不计罚。
  // read "GATTACA" 嵌入 reference "TTGATTACATTT" 中间,
  // 全局比对会对首尾空位罚 5+4=9 分, 半全局不应罚。
  let r = @src.semiglobal_align(
    "TTGATTACATTT",
    "GATTACA",
    @src.SGScoring::new(1.0, -1.0, -2.0),
  )
  // 7 个匹配, 0 个空位/错配 → score = 7
  assert_eq(r.score(), 7.0)
  assert_eq(r.aligned1(), "GATTACA")
  assert_eq(r.aligned2(), "GATTACA")
  assert_eq(r.start1(), 2)
  assert_eq(r.end1(), 9)
  assert_eq(r.start2(), 0)
  assert_eq(r.end2(), 7)
}

///|
test "sg_differs_from_global_on_terminal_gaps" {
  // 全局比对 (NW) 对同样的嵌入会给 7 - 9 = -2 分, 半全局给 7 分
  let seq1 = "TTGATTACATTT"
  let seq2 = "GATTACA"
  let sg = @src.semiglobal_align(
    seq1,
    seq2,
    @src.SGScoring::new(1.0, -1.0, -2.0),
  )
  let nw = @src.needleman_wunsch(
    seq1,
    seq2,
    @src.NWScoring::new(1.0, -1.0, -2.0),
  )
  assert_eq(sg.score(), 7.0)
  // 全局: 7 match - 5 gap * 2 = -3 (12 与 7 长度差至少 5 个空位)
  assert_eq(nw.score, -3.0)
}

///|
test "sg_trailing_free_gaps" {
  // 末端自由空位: seq1 比 seq2 长, 多出的尾巴不计罚
  let r = @src.semiglobal_align(
    "ACGTAAAA",
    "ACGT",
    @src.SGScoring::new(1.0, -1.0, -2.0),
  )
  assert_eq(r.score(), 4.0)
  assert_eq(r.aligned1(), "ACGT")
  assert_eq(r.aligned2(), "ACGT")
  assert_eq(r.end1(), 4)
}

///|
test "sg_internal_gaps_still_penalized" {
  // 中间空位仍要罚: gap(-1) 优于 mismatch(-3), 且空位两侧
  // 均有匹配 (ACG..C), 无法靠末端自由空位逃避
  let r = @src.semiglobal_align(
    "ACGTC",
    "ACGC",
    @src.SGScoring::new(2.0, -3.0, -1.0),
  )
  // 4 match * 2 - 1 gap * 1 = 7
  assert_eq(r.score(), 7.0)
  assert_eq(r.aligned1(), "ACGTC")
  assert_eq(r.aligned2(), "ACG-C")
  assert_eq(r.num_gaps_seq2(), 1)
  assert_eq(r.num_gaps_seq1(), 0)
}

///|
test "sg_mismatch_scoring" {
  // 错配按 mismatch_penalty 计
  let r = @src.semiglobal_align(
    "ACGT",
    "AGGT",
    @src.SGScoring::new(2.0, -1.0, -3.0),
  )
  assert_eq(r.aligned1(), "ACGT")
  assert_eq(r.aligned2(), "AGGT")
  // 3 match * 2 + 1 mismatch * (-1) = 5
  assert_eq(r.score(), 5.0)
  assert_eq(r.num_matches(), 3)
  assert_eq(r.num_mismatches(), 1)
  // 逐列: A|A, C.G, G|G, T|T
  assert_eq(r.match_line(), "|.||")
}

///|
test "sg_empty_sequences" {
  // 空序列对空序列: 得分 0
  let r = @src.semiglobal_align("", "", @src.SGScoring::new(1.0, -1.0, -1.0))
  assert_eq(r.score(), 0.0)
  assert_eq(r.aligned1(), "")
  assert_eq(r.aligned2(), "")
  // 空序列对非空: 末端自由空位使 score 仍为 0
  let r2 = @src.semiglobal_align(
    "ACGT",
    "",
    @src.SGScoring::new(1.0, -1.0, -1.0),
  )
  assert_eq(r2.score(), 0.0)
}

///|
test "sg_result_coordinates" {
  // start/end 应排除自由空位的侧翼
  let r = @src.semiglobal_align(
    "AAAACGTCCCGGG",
    "ACGT",
    @src.SGScoring::new(1.0, -1.0, -2.0),
  )
  assert_eq(r.start1(), 3)
  assert_eq(r.end1(), 7)
  assert_eq(r.start2(), 0)
  assert_eq(r.end2(), 4)
  assert_eq(r.score(), 4.0)
}

///|
test "sg_align_dna_convenience" {
  let r = @src.semiglobal_align_dna("GATTACA", "TTGATTACATTT")
  // default_dna: match=2, mismatch=-1, gap=-1
  assert_eq(r.score(), 14.0)
}

///|
test "sg_matrix_semiglobal_protein" {
  // 蛋白质序列 + BLOSUM62
  let r = @src.semiglobal_align_blosum62("HEAGAWGHEE", "PAWHEAE", -8.0)
  // PAWHEAE 应能通过自由空位对上 HEAGAWGHEE 的 HE..AWGHE 段
  assert_true(r.score() > 0.0)
  assert_true(r.aligned1().length() == r.aligned2().length())
}

///|
test "sg_pam250_convenience" {
  let r = @src.semiglobal_align_pam250("HEAGAWGHEE", "PAWHEAE", -8.0)
  assert_true(r.aligned1().length() == r.aligned2().length())
}

///|
test "sg_protein_scoring_convenience" {
  let r = @src.semiglobal_align_protein("MVLSPADKTNVKAAW", "MVLSPADKSNVKAAW")
  // default_protein: match=5, mismatch=-4, gap=-3 → 14*5 + 1*(-4) = 66
  assert_eq(r.score(), 66.0)
}

///|
test "sg_align_banded_semiglobal_now_exact" {
  // align_banded 的 Semiglobal 模式现已接入真实现:
  // 末端空位不再照罚 (此前是"用全局比对近似"的假实现)
  let r = @src.align_banded_semiglobal("TTGATTACATTT", "GATTACA", 3)
  assert_eq(r.get_score(), 7.0)
  assert_eq(r.get_aligned_seq1(), "GATTACA")
  assert_eq(r.get_aligned_seq2(), "GATTACA")
}

///|
test "sg_align_banded_mode_semiglobal" {
  let r = @src.align_banded(
    "ACGTAAAA",
    "ACGT",
    k=2,
    mode=@src.BandedAlignMode::Semiglobal,
    match_reward=1.0,
    mismatch_penalty=-1.0,
    gap_penalty=-2.0,
  )
  assert_eq(r.get_score(), 4.0)
  assert_eq(r.get_aligned_seq1(), "ACGT")
}

///|
test "sg_dp_matrix_dimensions" {
  // DP 矩阵应为 (n+1) x (m+1)
  let r = @src.semiglobal_align(
    "ACG",
    "AT",
    @src.SGScoring::new(1.0, -1.0, -1.0),
  )
  let dp = r.score_matrix()
  assert_eq(dp.length(), 4)
  assert_eq(dp[0].length(), 3)
  // 首行首列应全为 0 (自由空位初始化)
  for j = 0; j < 3; j = j + 1 {
    assert_eq(dp[0][j], 0.0)
  }
  for i = 0; i < 4; i = i + 1 {
    assert_eq(dp[i][0], 0.0)
  }
}

///|
test "sg_format_output" {
  let r = @src.semiglobal_align(
    "ACGT",
    "AGGT",
    @src.SGScoring::new(2.0, -1.0, -3.0),
  )
  let f = r.format()
  let lines = f.split("\n").to_array()
  assert_eq(lines[0].to_owned(), "ACGT")
  assert_eq(lines[1].to_owned(), "|.||")
  assert_eq(lines[2].to_owned(), "AGGT")
}

///|
test "sg_alignment_length_consistency" {
  let r = @src.semiglobal_align(
    "TTGATTACATTT",
    "GATTACA",
    @src.SGScoring::new(1.0, -1.0, -2.0),
  )
  assert_eq(r.alignment_length(), r.aligned1().length())
  assert_eq(r.alignment_length(), r.match_line().length())
}

///|
test "sg_format_wrapped" {
  let r = @src.semiglobal_align(
    "ACGTACGTAC",
    "ACGTACGTAC",
    @src.SGScoring::new(1.0, -1.0, -1.0),
  )
  let wrapped = r.format_wrapped(4)
  assert_true(wrapped.contains("Block 1"))
  assert_true(wrapped.contains("Block 3"))
}
