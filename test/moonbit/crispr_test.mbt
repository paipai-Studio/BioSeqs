// ===========================================================================
// crispr_test.mbt — CRISPR gRNA 设计器测试
// ===========================================================================

///|
fn crispr_test_seq() -> String {
  // 含多个 NGG 位点的序列
  "ATGCGTACCTGAGGCTAAGCTTACGGATCCGGGATCTAGAATTCCGGAAGCTTGGATCCCTAGCGGATC"
}

///|
test "crispr pam char matches" {
  assert_true(@src.crispr_pam_char_matches("N", 'A'))
  assert_true(@src.crispr_pam_char_matches("G", 'G'))
  assert_false(@src.crispr_pam_char_matches("G", 'A'))
  assert_true(@src.crispr_pam_char_matches("S", 'C'))
  assert_true(@src.crispr_pam_char_matches("S", 'G'))
  assert_false(@src.crispr_pam_char_matches("S", 'A'))
}

///|
test "crispr scan guides forward strand" {
  let seq = "AAAACCCCGGGGTTTTAGGAAAACCCCGGGGTTTTAGG"
  let guides = @src.crispr_scan_guides(seq)
  // 至少发现正向 NGG 位点(位置 16 的 AGG)
  assert_true(guides.length() > 0)
  let mut found_forward = false
  for g in guides {
    if g.strand == 0 {
      found_forward = true
      // 正向: PAM 前的 20nt 为 spacer
      assert_eq(seq[g.start:g.end].to_owned(), g.spacer)
    }
  }
  assert_true(found_forward)
}

///|
test "crispr scan guides reverse strand" {
  let seq = "AAAACCCCGGGGTTTTCCAAAAACCCCGGGGTTTTCC"
  let guides = @src.crispr_scan_guides(seq)
  // 反向链: 原链 CCN -> 反向互补链 NGG
  let mut found_reverse = false
  for g in guides {
    if g.strand == 1 {
      found_reverse = true
      // spacer 为原链区段的反向互补
      let rc = @src.reverse_complement(seq[g.start:g.end].to_owned())
      assert_eq(rc, g.spacer)
    }
  }
  assert_true(found_reverse)
}

///|
test "crispr scan guides cas12a" {
  // Cas12a: TTTV PAM 在 5' 端
  let seq = "AAAATTTGAAACCCCGGGGTTTTAAAACCCCGGGGTTTTAAA"
  let guides = @src.crispr_scan_guides(
    seq,
    pam=@src.crispr_pam_cas12a(),
    spacer_len=23,
  )
  assert_true(guides.length() > 0)
  for g in guides {
    // Cas12a spacer 紧跟 PAM 之后(正向)
    if g.strand == 0 {
      assert_eq(g.end - g.start, 23)
    }
  }
}

///|
test "crispr on target score gc optimal" {
  // GC 40-60% 无不良特征
  let good = "GATCGATCGATCGATCGATC"
  let s_good = @src.crispr_on_target_score(good)
  // 极端 GC
  let bad = "ATATATATATATATATATAT"
  let s_bad = @src.crispr_on_target_score(bad)
  assert_true(s_good > s_bad)
}

///|
test "crispr on target score tttt penalty" {
  let normal = "GATCGATCGATCGATCGATC"
  let term = "GATCGATTTTTGATCGATCG"
  assert_true(
    @src.crispr_on_target_score(normal) > @src.crispr_on_target_score(term),
  )
}

///|
test "crispr on target score homopolymer penalty" {
  let normal = "GATCGATCGATCGATCGATC"
  let homo = "GGGGGGGGGGGGGGGGGGGG"
  assert_true(
    @src.crispr_on_target_score(normal) > @src.crispr_on_target_score(homo),
  )
}

///|
test "crispr on target score range" {
  let s1 = @src.crispr_on_target_score("GATCGATCGATCGATCGATC")
  let s2 = @src.crispr_on_target_score("TTTTTTTTTTTTTTTTTTTT")
  assert_true(s1 >= 0.0 && s1 <= 1.0)
  assert_true(s2 >= 0.0 && s2 <= 1.0)
  assert_true(s2 <= 0.1)
}

///|
test "crispr off target site score perfect match" {
  let spacer = "GATCGATCGATCGATCGATC"
  let site = "GATCGATCGATCGATCGATC"
  assert_true(@src.crispr_off_target_site_score(spacer, site) > 0.99)
}

///|
test "crispr off target site score mismatches" {
  let spacer = "GATCGATCGATCGATCGATC"
  // 3' 端(PAM 近端)错配比 5' 端惩罚更重
  let mism_5p = "TATCGATCGATCGATCGATC" // 位置 0 错配
  let mism_3p = "GATCGATCGATCGATCGATT" // 位置 19 错配
  let s_5p = @src.crispr_off_target_site_score(spacer, mism_5p)
  let s_3p = @src.crispr_off_target_site_score(spacer, mism_3p)
  assert_true(s_5p > s_3p)
}

///|
test "crispr eval off targets finds exact duplicate" {
  // 参考含两个相同 spacer 拷贝(非周期序列避免移位错配命中)
  let spacer = "ACGTTGCAATCGGTACCATG"
  let reference = "AAAA" + spacer + "GGAA" + spacer + "GGAAAA"
  let guide : @src.CrisprGuide = {
    spacer,
    pam: "AGG",
    strand: 0,
    start: 4,
    end: 24,
    gc_percent: 50.0,
    on_target: 0.0,
    specificity: 1.0,
    total_score: 0.0,
    off_targets: [],
    warnings: [],
  }
  let (offs, spec) = @src.crispr_eval_off_targets(guide, reference)
  // 另一份拷贝(位置 4+20+4=28)应为 off-target, 原位点(位置 4)排除
  assert_eq(offs.length(), 1)
  assert_eq(offs[0].pos, 28)
  assert_eq(offs[0].mismatches, 0)
  assert_true(offs[0].site_score > 0.99)
  // 一个完全匹配 off-target -> 特异性接近 0
  assert_true(spec < 0.05)
}

///|
test "crispr eval off targets specific guide" {
  // 唯一非周期 spacer -> 无 off-target, 特异性 1
  let spacer = "TGCAATCGGTACCATGGCTA"
  let reference = "TTTT" + spacer + "GG" + "CCCCTAGGCATTGCATAGCTAGGCATCG"
  let guide : @src.CrisprGuide = {
    spacer,
    pam: "AGG",
    strand: 0,
    start: 4,
    end: 24,
    gc_percent: 45.0,
    on_target: 0.0,
    specificity: 1.0,
    total_score: 0.0,
    off_targets: [],
    warnings: [],
  }
  let (offs, spec) = @src.crispr_eval_off_targets(guide, reference)
  assert_eq(offs.length(), 0)
  assert_true(spec > 0.99)
}

///|
test "crispr design guides full flow" {
  let seq = "ATGCGTACCTGAGGCTAAGCTTACGGATCCGGGATCTAGAATTCCGGAAGCTTGGATCCCTAGCGGATC"
  let guides = @src.crispr_design_guides(seq, reference=seq, top_n=10)
  assert_true(guides.length() > 0)
  assert_true(guides.length() <= 10)
  // 排序: 总分降序
  for i = 1; i < guides.length(); i = i + 1 {
    assert_true(guides[i - 1].total_score >= guides[i].total_score)
  }
  // 全部字段已打分(specificity 可为 0: 序列含重复片段)
  for g in guides {
    assert_true(g.on_target > 0.0)
    assert_true(g.specificity >= 0.0)
    assert_true(g.total_score > 0.0)
  }
}

///|
test "crispr design guides no reference" {
  // 无参考序列: 跳过 off-target, 特异性默认 1
  let seq = "ATGCGTACCTGAGGCTAAGCTTACGGATCCGGGATCTAGAATTCCGGAAGCTTGGATCCCTAGCGGATC"
  let guides = @src.crispr_design_guides(seq, top_n=5)
  for g in guides {
    assert_eq(g.off_targets.length(), 0)
    assert_true(g.specificity > 0.99)
  }
}

///|
test "crispr design guides tttt warning" {
  let seq = "AAAAGATCGATTTTTGATCGATAGGAAAA"
  let guides = @src.crispr_design_guides(seq, top_n=50)
  let mut has_tttt_warn = false
  for g in guides {
    if g.spacer.contains("TTTT") {
      for w in g.warnings {
        if w.contains("TTTT") {
          has_tttt_warn = true
        }
      }
    }
  }
  assert_true(has_tttt_warn)
}

///|
test "crispr guides summary" {
  let seq = crispr_test_seq()
  let guides = @src.crispr_design_guides(seq, top_n=5)
  let s = @src.crispr_guides_summary(guides)
  assert_true(s.contains("crispr"))
  assert_true(s.contains("guides"))
}

///|
test "crispr guides json output" {
  let seq = crispr_test_seq()
  let guides = @src.crispr_design_guides(seq, top_n=3)
  let json = @src.crispr_guides_json(guides)
  assert_true(json.contains("\"spacer\""))
  assert_true(json.contains("\"pam\""))
  assert_true(json.contains("\"on_target\""))
  // parse_json 闭环
  let v = @src.parse_json(json)
  assert_true(v is Some(_))
}

///|
test "crispr spg pam" {
  let spg = @src.crispr_pam_spg()
  assert_eq(spg.pattern, ["S", "G", "G"])
  // SGG 匹配 GGG 但不匹配 AGG
  assert_true(@src.crispr_pam_char_matches("S", 'G'))
  assert_false(@src.crispr_pam_char_matches("S", 'A'))
}

///|
test "crispr pam reverse matching" {
  // NGG 反向互补 = CCN: 原链 CCT 应匹配
  let pam = @src.crispr_pam_spcas9()
  assert_true(@src.crispr_pam_matches_reverse(pam, "AACCTGGAA", 2))
  assert_false(@src.crispr_pam_matches_reverse(pam, "AAGGGGGAA", 2))
}
