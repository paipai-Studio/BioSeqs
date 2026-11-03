// ===========================================================================
// structure_screening_test.mbt — 蛋白结构筛查管线测试
// ===========================================================================

///|
test "structure screening empty raises" {
  let mut raised = false
  try {
    let r = @src.structure_screening("")
    ignore(r)
  } catch {
    _ => raised = true
  }
  assert_true(raised)
}

///|
test "structure screening alpha rich sequence" {
  // 富含 Ala/Leu/Glu/Lys 的序列 -> 高螺旋
  let seq = "AEAALKREAAALKREAAALKREAAALKREAAALK"
  let report = @src.structure_screening(seq)
  assert_eq(report.length, seq.length())
  assert_true(report.helix_percent > 30.0)
  assert_true(report.helix_regions.length() > 0 || report.helix_percent > 25.0)
  assert_eq(report.fold_class, "all-alpha")
}

///|
test "structure screening beta rich sequence" {
  // 富含 Val/Ile/Tyr/Phe 的序列 -> 高折叠
  let seq = "VVIVTYYIVFVVIVTYYIVFVVIVTYYIVFVVIVTYYIVF"
  let report = @src.structure_screening(seq)
  assert_true(report.sheet_percent > 30.0)
  assert_eq(report.fold_class, "all-beta")
}

///|
test "structure screening gly pro coil sequence" {
  // Gly/Pro/Ser 重复 -> 低规则结构
  let seq = "GPGPSPGPGPSPGPGPSPGPGPSPGPGPSPGPGPSPGPGPS"
  let report = @src.structure_screening(seq)
  assert_true(report.helix_percent + report.sheet_percent < 20.0)
  let fc = report.fold_class
  assert_true(fc == "coil" || fc == "mixed")
}

///|
test "structure screening percentages sum" {
  let seq = "AEAAALKREVVIVTYYIVFAEAAALKREVVIVTYYIVFAEAAALKREVVIVTYYIVF"
  let report = @src.structure_screening(seq)
  let total = report.helix_percent +
    report.sheet_percent +
    report.turn_percent +
    report.coil_percent
  assert_true(total > 99.9 && total < 100.1)
}

///|
test "structure screening report fields" {
  let seq = "MKVLAAALVKVLAAALGVLAAAL"
  let report = @src.structure_screening(seq)
  assert_eq(report.ss_prediction.length(), seq.length())
  for s in report.ss_prediction {
    assert_true(s == "H" || s == "E" || s == "T" || s == "C")
  }
}

///|
test "structure screening ss string" {
  let seq = "MKVLAAALVKVLAAALGVLAAAL"
  let report = @src.structure_screening(seq)
  let ss = report.ss_string()
  assert_eq(ss.length(), seq.length())
}

///|
test "structure screening json output" {
  let seq = "MKVLAAALVKVLAAALGVLAAAL"
  let report = @src.structure_screening(seq)
  let json = @src.structure_screening_report_json(report)
  assert_true(json.contains("\"module\":\"structure-screening\""))
  assert_true(json.contains("\"fold_class\""))
  assert_true(json.contains("\"helix_percent\""))
  // 与 parse_json 闭环
  let v = @src.parse_json(json)
  assert_true(v is Some(_))
}

///|
test "structure screening text output" {
  let seq = "MKVLAAALVKVLAAALGVLAAAL"
  let report = @src.structure_screening(seq)
  let text = @src.structure_screening_report_text(report)
  assert_true(text.contains("Structure Screening Report"))
  assert_true(text.contains("fold class"))
  assert_true(text.contains("screening-grade"))
}

///|
test "structure fold class rules" {
  // 直接测试分类规则
  let pred_h : Array[String] = []
  for i = 0; i < 30; i = i + 1 {
    pred_h.push("H")
  }
  assert_eq(@src.structure_fold_class(60.0, 5.0, pred_h), "all-alpha")
  let pred_e : Array[String] = []
  for i = 0; i < 30; i = i + 1 {
    pred_e.push("E")
  }
  assert_eq(@src.structure_fold_class(5.0, 60.0, pred_e), "all-beta")
  // 交替 -> alpha-beta
  let pred_alt : Array[String] = []
  for i = 0; i < 24; i = i + 1 {
    pred_alt.push(if i % 4 < 2 { "H" } else { "E" })
  }
  assert_eq(@src.structure_fold_class(30.0, 30.0, pred_alt), "alpha-beta")
  // 分离 -> alpha+beta
  let pred_sep : Array[String] = []
  for i = 0; i < 30; i = i + 1 {
    pred_sep.push(if i < 12 { "H" } else { "E" })
  }
  assert_eq(@src.structure_fold_class(30.0, 30.0, pred_sep), "alpha+beta")
  assert_eq(@src.structure_fold_class(10.0, 8.0, []), "coil")
}

///|
test "structure fold class description" {
  assert_true(
    @src.structure_fold_class_description("all-alpha").contains("alpha"),
  )
  assert_true(
    @src.structure_fold_class_description("all-beta").contains("beta"),
  )
  assert_true(
    @src.structure_fold_class_description("alpha-beta").contains("交替"),
  )
  assert_true(@src.structure_fold_class_description("coil").contains("无序"))
}

///|
test "structure screening high disorder warning" {
  // Gly/Ser/Asn/Gln 富集 -> IUPred 高无序 -> 警告
  let seq = "GSGNQGSGNQGSGNQGSGNQGSGNQGSGNQGSGNQGSGNQGSGNQGSGNQ"
  let report = @src.structure_screening(seq)
  // 至少有一条警告(无序或低规则结构)
  assert_true(report.warnings.length() > 0)
}

///|
test "structure screening multiqc integration" {
  // 与 E4 聚合报告集成: 结构筛查 -> multiqc 模块
  let report = @src.structure_screening("MKVLAAALVKVLAAALGVLAAAL")
  let m = @src.multiqc_module("structure-screening")
  let _ = m.add_sample("protein1", [
    ("helix_percent", report.helix_percent),
    ("sheet_percent", report.sheet_percent),
    ("coil_percent", report.coil_percent),
  ])
  let mr = @src.multiqc_report()
  let _ = mr.add_module(m)
  let md = @src.multiqc_write_markdown(mr)
  assert_true(md.contains("structure-screening"))
  assert_true(md.contains("protein1"))
}
