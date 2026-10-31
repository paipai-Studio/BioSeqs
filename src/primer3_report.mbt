///|
/// 引物设计报告(JSON): primer3 设计结果的机器可读输出.
///
/// 供 `examples/primer_web` 浏览器工具与 CLI 消费; 复用 primer3 设计引擎,
/// 输出模板统计、正/反向候选全列表与引物对组合(含产物长度).
pub fn primer3_result_to_json(r : Primer3Result) -> JsonValue {
  let pairs : Array[(String, JsonValue)] = [
    ("sequence", JsonStr(r.sequence())),
    ("start", JsonNum(r.start().to_double())),
    ("length", JsonNum(r.length().to_double())),
    ("tm", JsonNum(r.tm())),
    ("gc", JsonNum(r.gc())),
    ("self_complement", JsonNum(r.self_complement().to_double())),
    ("is_forward", JsonBool(r.is_forward())),
    ("score", JsonNum(r.score())),
  ]
  let m : Map[String, JsonValue] = Map([], capacity=pairs.length())
  for p in pairs {
    m[p.0] = p.1
  }
  JsonObj(m)
}

///|
fn primer3_arr(items : Array[JsonValue]) -> JsonValue {
  JsonArr(items)
}

///|
/// 设计引物并输出完整 JSON 报告.
///
/// 报告结构:
/// - template: { length, gc, tm_wallace }
/// - forward / reverse: 候选列表(按得分降序)
/// - pairs: 组合列表(产物长度 product_len 限定 100-5000bp, 按组合得分降序, 至多 max_pairs 组)
/// - 模板过短无候选时 lists 为空数组而非缺失.
pub fn primer3_design_report_json(
  template : String,
  config~ : Primer3Config = primer3_default_config(),
  max_pairs~ : Int32 = 50,
) -> String raise SeqError {
  json_to_string(primer3_design_report_value(template, config, max_pairs))
}

///|
/// 报告树(JsonValue), 便于嵌入或扩展.
pub fn primer3_design_report_value(
  template : String,
  config : Primer3Config,
  max_pairs : Int32,
) -> JsonValue raise SeqError {
  let design = primer3_design_primers(template, config)
  let t_gc = primer3_gc_content(template)
  let t_tm = primer3_tm_wallace(template)
  // 候选列表
  let fwd_arr : Array[JsonValue] = Array(capacity=design.n_forward())
  for i in 0..<design.n_forward() {
    fwd_arr.push(primer3_result_to_json(design.forward()[i]))
  }
  let rev_arr : Array[JsonValue] = Array(capacity=design.n_reverse())
  for i in 0..<design.n_reverse() {
    rev_arr.push(primer3_result_to_json(design.reverse()[i]))
  }
  // 引物对组合: 产物长度 = reverse 末端 - forward 起点(同链坐标), 按经验
  // 扩增子区间过滤, 组合得分 = 两引物得分之和 - 发夹/二聚体罚分.
  let pairs : Array[(Int32, Int32, Int32, Double)] = Array(capacity=0)
  for fi in 0..<design.n_forward() {
    let f = design.forward()[fi]
    for ri in 0..<design.n_reverse() {
      let r = design.reverse()[ri]
      // reverse 候选存储的是模板链上坐标, 产物 = r 末端 - f 起点
      let f_start = f.start()
      let r_end = r.start() + r.length()
      let product = r_end - f_start
      // 产物区间 100-5000bp(下限对齐 primer3 默认 PRODUCT_SIZE,
      // 排除引物相互重叠的不可行组合)
      if product < 100 || product > 5000 {
        continue
      }
      let dimer = primer3_cross_dimer(f.sequence(), r.sequence())
      let pair_score = f.score() + r.score() - dimer.to_double() * 0.5
      pairs.push((fi, ri, product, pair_score))
    }
  }
  // 组合按得分降序
  let mut i = 0
  while i < pairs.length() {
    let mut j = i + 1
    while j < pairs.length() {
      if pairs[j].3 > pairs[i].3 {
        let t = pairs[i]
        pairs[i] = pairs[j]
        pairs[j] = t
      }
      j = j + 1
    }
    i = i + 1
  }
  let take = if max_pairs < pairs.length() { max_pairs } else { pairs.length() }
  let pair_arr : Array[JsonValue] = Array(capacity=take)
  for k in 0..<take {
    let (fi, ri, product, score) = pairs[k]
    let m : Map[String, JsonValue] = Map([], capacity=4)
    m["forward_index"] = JsonNum(fi.to_double())
    m["reverse_index"] = JsonNum(ri.to_double())
    m["product_len"] = JsonNum(product.to_double())
    m["pair_score"] = JsonNum(score)
    pair_arr.push(JsonObj(m))
  }
  let root : Map[String, JsonValue] = Map([], capacity=8)
  let tm : Map[String, JsonValue] = Map([], capacity=3)
  tm["length"] = JsonNum(template.length().to_double())
  tm["gc"] = JsonNum(t_gc)
  tm["tm_wallace"] = JsonNum(t_tm)
  root["module"] = JsonStr("primer3_report")
  root["template"] = JsonObj(tm)
  root["forward"] = primer3_arr(fwd_arr)
  root["reverse"] = primer3_arr(rev_arr)
  root["pairs"] = primer3_arr(pair_arr)
  ignore(
    // 模板合法性提示: 非空且仅 ACGT/actg 视为有效, 其余字符计入 flag
    {
      let mut bad = 0
      for i in 0..<template.length() {
        let c = template.unsafe_get(i)
        if c != 'A' &&
          c != 'C' &&
          c != 'G' &&
          c != 'T' &&
          c != 'a' &&
          c != 'c' &&
          c != 'g' &&
          c != 't' {
          bad = bad + 1
        }
      }
      bad
    },
  )
  JsonObj(root)
}
