///|
/// structure_screening.mbt — 蛋白结构筛查管线(二级结构 -> 折叠类 -> 全景报告).
///
/// 定位: 结构预测**筛查级**(screening), 非深度学习预测器:
/// - 二级结构: 复用 Chou-Fasman(1974) 倾向法(chou_fasman_predict)
/// - 折叠类: SCOP 一级分类规则(all-alpha/all-beta/alpha-beta/alpha+beta/coil),
///   按二级结构含量与 H/E 段交替模式判定
/// - 全景整合: 无序区(IUPred) + 卷曲螺旋(COILS) + 溶解性警告
///
/// 输出 JSON/文本报告, 衔接结构分析栈(dssp/tm_score/cealign).
/// 明确声明: 筛查级方法勿用于替代 AF2/ESM 级精度预测. 数据不出本机.
pub suberror StructureScreeningError {
  StructureScreeningError(String)
}

///|
/// 结构筛查全景报告.
pub struct StructureScreeningReport {
  sequence : String
  length : Int32
  ss_prediction : Array[String]
  helix_percent : Double
  sheet_percent : Double
  turn_percent : Double
  coil_percent : Double
  fold_class : String
  helix_regions : Array[(Int32, Int32)]
  sheet_regions : Array[(Int32, Int32)]
  disordered_regions : Array[(Int32, Int32)]
  coiled_coil_regions : Array[(Int32, Int32)]
  warnings : Array[String]
}

///|
/// 一次性完成结构筛查(二级结构 + 折叠类 + 无序 + 卷曲螺旋).
pub fn structure_screening(
  sequence : String,
) -> StructureScreeningReport raise StructureScreeningError {
  let n = sequence.length()
  if n == 0 {
    raise StructureScreeningError("empty sequence")
  }
  let warnings : Array[String] = Array(capacity=0)
  // 1. 二级结构(Chou-Fasman)
  let cf : ChouFasmanResult = chou_fasman_predict(sequence) catch {
    e => raise StructureScreeningError("chou-fasman failed")
  }
  let pred = cf.prediction()
  if pred.length() != n {
    raise StructureScreeningError(
      "prediction length mismatch: \{pred.length()} != \{n}",
    )
  }
  // 2. 含量统计
  let mut n_h = 0
  let mut n_e = 0
  let mut n_t = 0
  for s in pred {
    if s == "H" {
      n_h = n_h + 1
    } else if s == "E" {
      n_e = n_e + 1
    } else if s == "T" {
      n_t = n_t + 1
    }
  }
  let helix_percent = n_h.to_double() / n.to_double() * 100.0
  let sheet_percent = n_e.to_double() / n.to_double() * 100.0
  let turn_percent = n_t.to_double() / n.to_double() * 100.0
  let coil_percent = 100.0 - helix_percent - sheet_percent - turn_percent
  // 3. 无序(IUPred)与卷曲螺旋(COILS)
  let disorder : Array[(Int32, Int32)] = iupred(sequence).disordered_regions() catch {
    _ => {
      warnings.push("iupred failed; disorder regions omitted")
      []
    }
  }
  let coils : Array[(Int32, Int32)] = predict_coiled_coils(sequence).coiled_coil_regions() catch {
    _ => {
      warnings.push("coils failed; coiled-coil regions omitted")
      []
    }
  }
  // 4. 折叠类(SCOP 一级分类简化规则)
  let fold_class = structure_fold_class(helix_percent, sheet_percent, pred)
  // 5. 警告: 无序占比高 -> 结构化风险
  let disorder_len : Int32 = {
    let mut acc = 0
    for region in disorder {
      acc = acc + (region.1 - region.0 + 1)
    }
    acc
  }
  if disorder_len.to_double() / n.to_double() > 0.4 {
    warnings.push(
      "high disorder content (\{disorder_len}/\{n} residues); structure prediction unreliable",
    )
  }
  if coils.length() > 0 {
    warnings.push(
      "\{coils.length()} coiled-coil segment(s) detected; verify oligomerization",
    )
  }
  if helix_percent + sheet_percent < 20.0 {
    warnings.push(
      "low regular secondary structure (<20%); likely intrinsically disordered",
    )
  }
  {
    sequence,
    length: n,
    ss_prediction: pred,
    helix_percent,
    sheet_percent,
    turn_percent,
    coil_percent,
    fold_class,
    helix_regions: cf.helices(),
    sheet_regions: cf.sheets(),
    disordered_regions: disorder,
    coiled_coil_regions: coils,
    warnings,
  }
}

///|
/// SCOP 一级折叠类判定(简化规则):
/// - all-alpha: 螺旋 >= 40% 且折叠 < 15%
/// - all-beta: 折叠 >= 40% 且螺旋 < 15%
/// - alpha-beta: 双 >= 25% 且 H/E 段交替(切换 >= 4 次, beta-alpha-beta 交叠模式)
/// - alpha+beta: 双 >= 25% 且非交替(空间分离的 alpha 与 beta 区)
/// - coil: 其余
pub fn structure_fold_class(
  helix_percent : Double,
  sheet_percent : Double,
  ss_prediction : Array[String],
) -> String {
  if helix_percent >= 40.0 && sheet_percent < 15.0 {
    return "all-alpha"
  }
  if sheet_percent >= 40.0 && helix_percent < 15.0 {
    return "all-beta"
  }
  if helix_percent >= 25.0 && sheet_percent >= 25.0 {
    // H/E 段交替统计
    let mut switches = 0
    let mut prev = ""
    for s in ss_prediction {
      if s == "H" || s == "E" {
        if prev.length() > 0 && s != prev {
          switches = switches + 1
        }
        prev = s
      }
    }
    if switches >= 4 {
      return "alpha-beta"
    }
    return "alpha+beta"
  }
  if helix_percent + sheet_percent < 20.0 {
    return "coil"
  }
  "mixed"
}

///|
/// 折叠类中文说明.
pub fn structure_fold_class_description(fold_class : String) -> String {
  match fold_class {
    "all-alpha" => "全 alpha(螺旋主导, 如球状血红蛋白超家族)"
    "all-beta" => "全 beta(折叠主导, 如免疫球蛋白超家族)"
    "alpha-beta" => "alpha/beta 交替(交叠 beta-alpha-beta 桶/折叠)"
    "alpha+beta" => "alpha+beta(空间分离的螺旋与折叠区)"
    "coil" => "卷曲主导(低规则结构, 可能内在无序)"
    _ => "混合"
  }
}

///|
/// JSON 报告值.
pub fn structure_screening_report_value(
  report : StructureScreeningReport,
) -> JsonValue {
  let ss : Array[JsonValue] = Array(capacity=0)
  for i = 0; i < report.ss_prediction.length(); i = i + 1 {
    ss.push(JsonStr(report.ss_prediction[i]))
  }
  let helix_regions : Array[JsonValue] = Array(capacity=0)
  for r in report.helix_regions {
    helix_regions.push(
      JsonArr([JsonNum(r.0.to_double()), JsonNum(r.1.to_double())]),
    )
  }
  let sheet_regions : Array[JsonValue] = Array(capacity=0)
  for r in report.sheet_regions {
    sheet_regions.push(
      JsonArr([JsonNum(r.0.to_double()), JsonNum(r.1.to_double())]),
    )
  }
  let disordered : Array[JsonValue] = Array(capacity=0)
  for r in report.disordered_regions {
    disordered.push(
      JsonArr([JsonNum(r.0.to_double()), JsonNum(r.1.to_double())]),
    )
  }
  let coiled : Array[JsonValue] = Array(capacity=0)
  for r in report.coiled_coil_regions {
    coiled.push(JsonArr([JsonNum(r.0.to_double()), JsonNum(r.1.to_double())]))
  }
  let warnings : Array[JsonValue] = Array(capacity=0)
  for w in report.warnings {
    warnings.push(JsonStr(w))
  }
  let fields : Array[(String, JsonValue)] = [
    ("module", JsonStr("structure-screening")),
    ("length", JsonNum(report.length.to_double())),
    ("fold_class", JsonStr(report.fold_class)),
    ("helix_percent", JsonNum(report.helix_percent)),
    ("sheet_percent", JsonNum(report.sheet_percent)),
    ("turn_percent", JsonNum(report.turn_percent)),
    ("coil_percent", JsonNum(report.coil_percent)),
    ("ss_prediction", JsonArr(ss)),
    ("helix_regions", JsonArr(helix_regions)),
    ("sheet_regions", JsonArr(sheet_regions)),
    ("disordered_regions", JsonArr(disordered)),
    ("coiled_coil_regions", JsonArr(coiled)),
    ("warnings", JsonArr(warnings)),
  ]
  let m : Map[String, JsonValue] = Map([], capacity=fields.length())
  for f in fields {
    m.set(f.0, f.1)
  }
  JsonObj(m)
}

///|
/// JSON 文本报告.
pub fn structure_screening_report_json(
  report : StructureScreeningReport,
) -> String {
  json_to_string(structure_screening_report_value(report))
}

///|
/// 人类可读文本报告.
pub fn structure_screening_report_text(
  report : StructureScreeningReport,
) -> String {
  let sb = StringBuilder()
  sb.write_string(
    "Structure Screening Report (screening-grade, not a deep-learning predictor)\n",
  )
  sb.write_string("  length: \{report.length}\n")
  sb.write_string(
    "  fold class: \{report.fold_class} — \{structure_fold_class_description(report.fold_class)}\n",
  )
  sb.write_string(
    "  secondary structure: helix \{report.helix_percent}% / sheet \{report.sheet_percent}% / turn \{report.turn_percent}% / coil \{report.coil_percent}%\n",
  )
  sb.write_string("  helix regions: \{report.helix_regions.length()}\n")
  sb.write_string("  sheet regions: \{report.sheet_regions.length()}\n")
  sb.write_string(
    "  disordered regions: \{report.disordered_regions.length()}\n",
  )
  sb.write_string(
    "  coiled-coil regions: \{report.coiled_coil_regions.length()}\n",
  )
  for w in report.warnings {
    sb.write_string("  warning: \{w}\n")
  }
  sb.to_string()
}

///|
/// 预测序列压缩为字符串("HEEETT.." 形式).
pub fn StructureScreeningReport::ss_string(
  self : StructureScreeningReport,
) -> String {
  let sb = StringBuilder()
  for s in self.ss_prediction {
    sb.write_string(s)
  }
  sb.to_string()
}
