// ===========================================================================
// acmg_test.mbt — ACMG 2015 自动化分级与 ClinVar 摄取测试.
// ===========================================================================

///|
fn ac_test_clinvar() -> String {
  "#ClinVar simplified TSV\n" +
  "chrom\tpos\tref\talt\tclnsig\tstar\n" +
  "chr1\t101\tA\tG\tPathogenic\t3\n" +
  "chr1\t202\tC\tT\tBenign\t2\n" +
  "chr1\t303\tG\tA\tConflicting_classifications\t1\n" +
  "chr1\t404\tT\tC\tLikely_pathogenic\t2\n"
}

///|
test "acmg: ClinVar TSV 摄取与查找" {
  let records = @src.clinvar_parse(ac_test_clinvar())
  inspect(records.length(), content="4")
  let r0 = records[0]
  inspect(r0.chrom, content="chr1")
  inspect(r0.pos, content="101")
  inspect(r0.ref_allele, content="A")
  inspect(r0.alt, content="G")
  inspect(r0.clnsig, content="Pathogenic")
  inspect(r0.review_star, content="3")
  // 查找命中
  let hit = @src.clinvar_lookup(records, "chr1", 101, "A", "G")
  assert_true(hit is Some(_))
  // 查找未命中(alt 不同)
  let miss = @src.clinvar_lookup(records, "chr1", 101, "A", "T")
  assert_true(miss is None)
}

///|
test "acmg: 证据推导(LoF/频率/后果)" {
  // LoF + 频率缺失
  let input = @src.acmg_default_input("chr1", 501, "C", "T", "stop_gained")
  let ev = @src.acmg_derive_evidence(input)
  inspect(ev.length(), content="1")
  inspect(ev[0], content="PVS1")
  // 移码 + 低频
  let input2 = @src.acmg_make_input(
    "chr1",
    502,
    "CA",
    "C",
    "frameshift_variant",
    popmax_af=0.0001,
  )
  let ev2 = @src.acmg_derive_evidence(input2)
  inspect(ev2.length(), content="2")
  assert_eq(ev2.contains("PVS1"), true)
  assert_eq(ev2.contains("PM2"), true)
  // 框内缺失
  let input3 = @src.acmg_default_input(
    "chr1", 503, "CAT", "C", "inframe_deletion",
  )
  let ev3 = @src.acmg_derive_evidence(input3)
  assert_eq(ev3.contains("PM4"), true)
  // 同义
  let input4 = @src.acmg_default_input(
    "chr1", 504, "A", "G", "synonymous_variant",
  )
  let ev4 = @src.acmg_derive_evidence(input4)
  assert_eq(ev4.contains("BP7"), true)
  // 高频良性
  let input5 = @src.acmg_make_input(
    "chr1",
    505,
    "A",
    "G",
    "missense_variant",
    popmax_af=0.1,
  )
  let ev5 = @src.acmg_derive_evidence(input5)
  assert_eq(ev5.contains("BA1"), true)
  // 中频 BS1
  let input5b = @src.acmg_make_input(
    "chr1",
    505,
    "A",
    "G",
    "missense_variant",
    popmax_af=0.02,
  )
  let ev5b = @src.acmg_derive_evidence(input5b)
  assert_eq(ev5b.contains("BS1"), true)
  assert_eq(ev5b.contains("BA1"), false)
}

///|
test "acmg: 证据推导(基因属性/de novo/预测/功能)" {
  // missense-only 基因的错义
  let input = @src.acmg_make_input(
    "chr1",
    601,
    "A",
    "G",
    "missense_variant",
    gene_missense_only=true,
  )
  let ev = @src.acmg_derive_evidence(input)
  assert_eq(ev.contains("PP2"), true)
  // LoF 基因的错义 -> BP1
  let input2 = @src.acmg_make_input(
    "chr1",
    601,
    "A",
    "G",
    "missense_variant",
    gene_lof_major=true,
  )
  let ev2 = @src.acmg_derive_evidence(input2)
  assert_eq(ev2.contains("BP1"), true)
  // de novo
  let input3 = @src.acmg_make_input(
    "chr1",
    601,
    "A",
    "G",
    "missense_variant",
    de_novo=true,
  )
  let ev3 = @src.acmg_derive_evidence(input3)
  assert_eq(ev3.contains("PM6"), true)
  // 计算预测有害/良性
  let input4 = @src.acmg_make_input(
    "chr1",
    601,
    "A",
    "G",
    "missense_variant",
    in_silico_score=Some(28.0),
  )
  let ev4 = @src.acmg_derive_evidence(input4)
  assert_eq(ev4.contains("PP3"), true)
  let input5 = @src.acmg_make_input(
    "chr1",
    601,
    "A",
    "G",
    "missense_variant",
    in_silico_score=Some(2.0),
  )
  let ev5 = @src.acmg_derive_evidence(input5)
  assert_eq(ev5.contains("BP4"), true)
  // 功能实验
  let input6 = @src.acmg_make_input(
    "chr1",
    601,
    "A",
    "G",
    "missense_variant",
    functional_score=Some(1.5),
  )
  let ev6 = @src.acmg_derive_evidence(input6)
  assert_eq(ev6.contains("PS3"), true)
  let input7 = @src.acmg_make_input(
    "chr1",
    601,
    "A",
    "G",
    "missense_variant",
    functional_score=Some(-1.0),
  )
  let ev7 = @src.acmg_derive_evidence(input7)
  assert_eq(ev7.contains("BS3"), true)
}

///|
test "acmg: ClinVar 证据补充" {
  let records = @src.clinvar_parse(ac_test_clinvar())
  // 同位点致病变异(>=2星) -> PS1 + PP5
  let hit = @src.clinvar_lookup(records, "chr1", 101, "A", "G")
  let ev = @src.acmg_clinvar_evidence(["PM2"], hit)
  inspect(ev.length(), content="3")
  assert_eq(ev.contains("PS1"), true)
  assert_eq(ev.contains("PP5"), true)
  // 良性收录 -> BP6
  let hit2 = @src.clinvar_lookup(records, "chr1", 202, "C", "T")
  let ev2 = @src.acmg_clinvar_evidence([], hit2)
  inspect(ev2.length(), content="1")
  assert_eq(ev2.contains("BP6"), true)
  // 冲突收录不计证据
  let hit3 = @src.clinvar_lookup(records, "chr1", 303, "G", "A")
  let ev3 = @src.acmg_clinvar_evidence(["PM2"], hit3)
  inspect(ev3.length(), content="1")
  // 1 星致病记录不计(可信度不足): 库内新增 1 星 Pathogenic 行
  let one_star = @src.clinvar_parse("chr1\t909\tA\tT\tPathogenic\t1\n")
  let hit5 = @src.clinvar_lookup(one_star, "chr1", 909, "A", "T")
  let ev5 = @src.acmg_clinvar_evidence([], hit5)
  inspect(ev5.length(), content="0")
}

///|
test "acmg: 组合规则 Pathogenic" {
  // PVS1 + 1PS(功能实验) -> Pathogenic
  let cls = @src.acmg_classify(["PVS1", "PS3", "PM2"])
  inspect(cls.classification, content="Pathogenic")
  assert_eq(cls.reasons.length() > 0, true)
  // 2PS -> Pathogenic
  let cls2 = @src.acmg_classify(["PS1", "PS3"])
  inspect(cls2.classification, content="Pathogenic")
  // PVS1 + 2PM -> Pathogenic
  let cls3 = @src.acmg_classify(["PVS1", "PM2", "PM4"])
  inspect(cls3.classification, content="Pathogenic")
  // 1PS + 3PM -> Pathogenic
  let cls4 = @src.acmg_classify(["PS3", "PM2", "PM4", "PM6"])
  inspect(cls4.classification, content="Pathogenic")
}

///|
test "acmg: 组合规则 Likely_pathogenic 与 VUS" {
  // PVS1 + 1PM -> Likely_pathogenic
  let cls = @src.acmg_classify(["PVS1", "PM2"])
  inspect(cls.classification, content="Likely_pathogenic")
  // 1PS + 1PM -> Likely_pathogenic
  let cls2 = @src.acmg_classify(["PS3", "PM2"])
  inspect(cls2.classification, content="Likely_pathogenic")
  // 3PM -> Likely_pathogenic
  let cls3 = @src.acmg_classify(["PM2", "PM4", "PM6"])
  inspect(cls3.classification, content="Likely_pathogenic")
  // 证据不足 -> VUS
  let cls4 = @src.acmg_classify(["PM2", "PP3"])
  inspect(cls4.classification, content="VUS")
  let cls5 = @src.acmg_classify([])
  inspect(cls5.classification, content="VUS")
}

///|
test "acmg: 组合规则 Benign/Likely_benign 与冲突降级" {
  // BA1 一票良性
  let cls = @src.acmg_classify(["BA1", "PM2", "PP3"])
  inspect(cls.classification, content="Benign")
  // 2BS -> Benign
  let cls2 = @src.acmg_classify(["BS1", "BS3"])
  inspect(cls2.classification, content="Benign")
  // 1BS + 1BP -> Likely_benign
  let cls3 = @src.acmg_classify(["BS1", "BP4"])
  inspect(cls3.classification, content="Likely_benign")
  // 2BP -> Likely_benign
  let cls4 = @src.acmg_classify(["BP4", "BP7"])
  inspect(cls4.classification, content="Likely_benign")
  // 致病与良性同时满足 -> 冲突降级 VUS
  let cls5 = @src.acmg_classify(["BA1", "PVS1", "PS3"])
  inspect(cls5.classification, content="VUS")
  assert_eq(cls5.reasons.length() > 0, true)
}

///|
test "acmg: 批量注释与报告" {
  let clinvar = @src.clinvar_parse(ac_test_clinvar())
  let inputs : Array[@src.AcmgInput] = []
  // LoF 低频 -> PVS1+PM2 -> Likely_pathogenic
  inputs.push(
    @src.acmg_make_input("chr1", 701, "C", "T", "stop_gained", popmax_af=0.0001),
  )
  // 高频错义 -> BA1 -> Benign
  inputs.push(
    @src.acmg_make_input(
      "chr1",
      702,
      "A",
      "G",
      "missense_variant",
      popmax_af=0.2,
    ),
  )
  // ClinVar 已知致病(chr1:101 A>G)同位点 -> PS1+PP5+PM2
  inputs.push(
    @src.acmg_make_input(
      "chr1",
      101,
      "A",
      "G",
      "missense_variant",
      popmax_af=0.0005,
    ),
  )
  let results = @src.acmg_annotate(inputs, clinvar~)
  inspect(results.length(), content="3")
  inspect(results[0].classification, content="Likely_pathogenic")
  inspect(results[1].classification, content="Benign")
  inspect(results[2].classification, content="Likely_pathogenic")
  inspect(results[2].clinvar_sig, content="Pathogenic(3star)")
  // TSV / JSON / 摘要
  let tsv = @src.acmg_report_tsv(results)
  assert_eq(tsv.contains("chrom\tpos\tref\talt"), true)
  assert_eq(tsv.contains("chr1\t701\tC\tT"), true)
  // 证据顺序: 频率类先于后果类
  assert_eq(tsv.contains("PM2&PVS1"), true)
  assert_eq(tsv.contains("Pathogenic(3star)"), true)
  let json = @src.acmg_json(results)
  assert_eq(json.contains("\"classification\""), true)
  assert_eq(json.contains("\"evidence\""), true)
  let line = @src.acmg_summary_line(results)
  assert_eq(line.contains("total=3"), true)
  assert_eq(line.contains("likely_pathogenic=2"), true)
  assert_eq(line.contains("benign=1"), true)
}
