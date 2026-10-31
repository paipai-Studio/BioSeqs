///|
/// multiqc_lite.mbt — MultiQC 风格多样本聚合报告(模块化指标 -> JSON/HTML/Markdown).
///
/// 对齐 MultiQC 的核心语义(本地无依赖实现):
/// - 模块化采集: 每个分析工具产出 MultiqcModule(样本 -> 指标键值)
/// - 聚合仪表盘: 样本 x 指标矩阵, 最优/最差单元格着色(方向感知)
/// - 三种输出: 单页 HTML(CSS 内嵌, 无外部依赖) / JSON(交换与 roundtrip) /
///   Markdown 表格(README 友好)
///
/// 内置采集器: fqqc(FASTQ QC), kallisto-lite(定量), gatk-lite(变异检出),
/// 并提供通用 multiqc_module_from_metrics 供任意模块接入. 数据不出本机.
pub(all) struct MultiqcModule {
  name : String
  /// 样本名 -> (指标名 -> 值)
  samples : Map[String, Map[String, Double]]
  /// 越低越好的指标(着色方向), 其余默认越高越好
  lower_better : Array[String]
}

///|
/// 聚合报告.
pub struct MultiqcReport {
  title : String
  modules : Array[MultiqcModule]
}

///|
/// 新建空模块.
pub fn multiqc_module(name : String) -> MultiqcModule {
  { name, samples: Map([], capacity=0), lower_better: [], }
}

///|
/// 添加样本指标(通用采集器入口).
pub fn MultiqcModule::add_sample(
  self : MultiqcModule,
  sample : String,
  metrics : Array[(String, Double)],
) -> MultiqcModule {
  let m : Map[String, Double] = Map([], capacity=0)
  for pair in metrics {
    m.set(pair.0, pair.1)
  }
  self.samples.set(sample, m)
  self
}

///|
/// 标记越低越好的指标.
pub fn MultiqcModule::with_lower_better(
  self : MultiqcModule,
  metrics : Array[String],
) -> MultiqcModule {
  { name: self.name, samples: self.samples, lower_better: metrics, }
}

///|
/// 新建空报告.
pub fn multiqc_report(title~ : String = "MultiQC-lite Report") -> MultiqcReport {
  { title, modules: Array(capacity=0), }
}

///|
/// 追加模块.
pub fn MultiqcReport::add_module(
  self : MultiqcReport,
  module_ : MultiqcModule,
) -> MultiqcReport {
  self.modules.push(module_)
  self
}

// ---------------------------------------------------------------------------
// 采集器
// ---------------------------------------------------------------------------

///|
/// fqqc 采集器: FASTQ QC 指标(reads/平均长度/平均质量/GC%).
pub fn multiqc_collect_fqqc(
  fastq : String,
  sample : String,
) -> MultiqcModule raise SeqIOError {
  let st = fqqc_compute(fastq)
  let mean_len = if st.n_reads > 0 {
    st.total_bases.to_double() / st.n_reads.to_double()
  } else {
    0.0
  }
  let mut q_sum = 0
  for q in st.mean_quals {
    q_sum = q_sum + q
  }
  let mean_q = if st.n_reads > 0 {
    q_sum.to_double() / st.n_reads.to_double()
  } else {
    0.0
  }
  let mut gc_sum = 0
  for g in st.gc_percents {
    gc_sum = gc_sum + g
  }
  let gc_pct = if st.n_reads > 0 {
    gc_sum.to_double() / st.n_reads.to_double()
  } else {
    0.0
  }
  let module_ = multiqc_module("fqqc")
  module_
  |> MultiqcModule::add_sample(sample, [
    ("reads", st.n_reads.to_double()),
    ("mean_length", mean_len),
    ("mean_quality", mean_q),
    ("gc_percent", gc_pct),
  ])
  |> MultiqcModule::with_lower_better([])
}

///|
/// kallisto-lite 采集器: 定量指标(reads/伪比对率/等价类数/检出转录本数).
pub fn multiqc_collect_kallisto(
  quant : KallistoQuant,
  sample : String,
) -> MultiqcModule {
  let pseudo_rate = if quant.n_reads > 0 {
    quant.n_pseudoaligned.to_double() / quant.n_reads.to_double() * 100.0
  } else {
    0.0
  }
  let mut n_expr = 0
  for c in quant.est_counts {
    if c > 0.5 {
      n_expr = n_expr + 1
    }
  }
  multiqc_module("kallisto-lite")
  |> MultiqcModule::add_sample(sample, [
    ("reads", quant.n_reads.to_double()),
    ("pseudoaligned_percent", pseudo_rate),
    ("equivalence_classes", quant.n_ec.to_double()),
    ("expressed_transcripts", n_expr.to_double()),
    ("em_iterations", quant.iterations.to_double()),
  ])
  |> MultiqcModule::with_lower_better(["em_iterations"])
}

///|
/// gatk-lite 采集器: 变异检出指标(SNP/插入/缺失/平均质量).
pub fn multiqc_collect_gatk(
  variants : Array[GatkVariant],
  sample : String,
) -> MultiqcModule {
  let mut n_snp = 0
  let mut n_ins = 0
  let mut n_del = 0
  let mut qual_sum = 0.0
  for v in variants {
    if v.ref.length() == 1 && v.alt.length() == 1 {
      n_snp = n_snp + 1
    } else if v.alt.length() > v.ref.length() {
      n_ins = n_ins + 1
    } else {
      n_del = n_del + 1
    }
    qual_sum = qual_sum + v.qual
  }
  let mean_qual = if variants.length() > 0 {
    qual_sum / variants.length().to_double()
  } else {
    0.0
  }
  multiqc_module("gatk-lite")
  |> MultiqcModule::add_sample(sample, [
    ("variants", variants.length().to_double()),
    ("snps", n_snp.to_double()),
    ("insertions", n_ins.to_double()),
    ("deletions", n_del.to_double()),
    ("mean_qual", mean_qual),
  ])
  |> MultiqcModule::with_lower_better([])
}

///|
/// kraken-lite 采集器: 宏基因组分类指标(读取/分类率/检出分类单元/平均置信度).
pub fn multiqc_collect_kraken(
  classifications : Array[KrakenClassification],
  abundances : Array[KrakenAbundance],
  sample : String,
) -> MultiqcModule {
  let total = classifications.length()
  let mut classified = 0
  let mut confidence_sum = 0.0
  for result in classifications {
    match result.taxid {
      Some(_) => classified = classified + 1
      None => ()
    }
    confidence_sum = confidence_sum + result.confidence
  }
  let unclassified = total - classified
  let classified_pct = if total > 0 {
    classified.to_double() / total.to_double() * 100.0
  } else {
    0.0
  }
  let mean_confidence = if total > 0 {
    confidence_sum / total.to_double()
  } else {
    0.0
  }
  let mut top_rel = 0.0
  for entry in abundances {
    if entry.relative_abundance > top_rel {
      top_rel = entry.relative_abundance
    }
  }
  multiqc_module("kraken-lite")
  |> MultiqcModule::add_sample(sample, [
    ("reads", total.to_double()),
    ("classified", classified.to_double()),
    ("unclassified", unclassified.to_double()),
    ("classified_percent", classified_pct),
    ("taxa_detected", abundances.length().to_double()),
    ("top_relative_percent", top_rel * 100.0),
    ("mean_confidence", mean_confidence),
  ])
  |> MultiqcModule::with_lower_better(["unclassified"])
}

// ---------------------------------------------------------------------------
// JSON 序列化 / 反序列化
// ---------------------------------------------------------------------------

///|
/// 报告 -> JSON 值.
fn mq_obj(pairs : Array[(String, JsonValue)]) -> JsonValue {
  let m : Map[String, JsonValue] = Map([], capacity=pairs.length())
  for p in pairs {
    m.set(p.0, p.1)
  }
  JsonObj(m)
}

///|
/// 报告 -> JSON 值.
pub fn multiqc_report_value(report : MultiqcReport) -> JsonValue {
  let modules : Array[JsonValue] = Array(capacity=0)
  for m in report.modules {
    let samples : Array[JsonValue] = Array(capacity=0)
    for sname in m.samples.keys() {
      let metrics = m.samples.get(sname).unwrap()
      let fields : Array[(String, JsonValue)] = Array(capacity=0)
      fields.push(("sample", JsonStr(sname)))
      for k in metrics.keys() {
        fields.push((k, JsonNum(metrics.get(k).unwrap())))
      }
      samples.push(mq_obj(fields))
    }
    let lb : Array[JsonValue] = Array(capacity=0)
    for x in m.lower_better {
      lb.push(JsonStr(x))
    }
    modules.push(
      mq_obj([
        ("module", JsonStr(m.name)),
        ("samples", JsonArr(samples)),
        ("lower_better", JsonArr(lb)),
      ]),
    )
  }
  mq_obj([("title", JsonStr(report.title)), ("modules", JsonArr(modules))])
}

///|
/// 报告 -> JSON 文本.
pub fn multiqc_write_json(report : MultiqcReport) -> String {
  json_to_string(multiqc_report_value(report))
}

///|
/// JSON 文本 -> 报告(roundtrip).
pub fn multiqc_parse_json(text : String) -> MultiqcReport? {
  let v = parse_json(text)
  match v {
    None => None
    Some(JsonObj(fields)) => {
      let mut title = "MultiQC-lite Report"
      let modules : Array[MultiqcModule] = Array(capacity=0)
      for f in fields {
        if f.0 == "title" {
          match f.1 {
            JsonStr(s) => title = s
            _ => ()
          }
        } else if f.0 == "modules" {
          match f.1 {
            JsonArr(arr) =>
              for mv in arr {
                match mv {
                  JsonObj(mfields) => {
                    let mut name = ""
                    let samples : Map[String, Map[String, Double]] = Map(
                      [],
                      capacity=0,
                    )
                    let lower_better : Array[String] = Array(capacity=0)
                    for mf in mfields {
                      if mf.0 == "module" {
                        match mf.1 {
                          JsonStr(s) => name = s
                          _ => ()
                        }
                      } else if mf.0 == "lower_better" {
                        match mf.1 {
                          JsonArr(lb) =>
                            for x in lb {
                              match x {
                                JsonStr(s) => lower_better.push(s)
                                _ => ()
                              }
                            }
                          _ => ()
                        }
                      } else if mf.0 == "samples" {
                        match mf.1 {
                          JsonArr(sarr) =>
                            for sv in sarr {
                              match sv {
                                JsonObj(sfields) => {
                                  let mut sname = ""
                                  let metrics : Map[String, Double] = Map(
                                    [],
                                    capacity=0,
                                  )
                                  for sf in sfields {
                                    if sf.0 == "sample" {
                                      match sf.1 {
                                        JsonStr(s) => sname = s
                                        _ => ()
                                      }
                                    } else {
                                      match sf.1 {
                                        JsonNum(x) => metrics.set(sf.0, x)
                                        _ => ()
                                      }
                                    }
                                  }
                                  if sname.length() > 0 {
                                    samples.set(sname, metrics)
                                  }
                                }
                                _ => ()
                              }
                            }
                          _ => ()
                        }
                      }
                    }
                    if name.length() > 0 {
                      modules.push({ name, samples, lower_better, })
                    }
                  }
                  _ => ()
                }
              }
            _ => ()
          }
        }
      }
      Some({ title, modules, })
    }
    _ => None
  }
}

// ---------------------------------------------------------------------------
// HTML 仪表盘
// ---------------------------------------------------------------------------

///|
fn mq_html_escape(s : String) -> String {
  s
  .replace(old="&", new="&amp;")
  .replace(old="<", new="&lt")
  .replace(old=">", new="&gt")
}

///|
fn mq_fmt(x : Double) -> String {
  if x == x.trunc() && x > -1.0e15 && x < 1.0e15 {
    x.to_int().to_string()
  } else {
    let r = x * 100.0
    if r == r.trunc() {
      x.to_string()
    } else {
      mq_fmt_small((x * 100.0).round() / 100.0)
    }
  }
}

///|
fn mq_fmt_small(x : Double) -> String {
  let s = x.to_string()
  if s.contains(".") || s.contains("e") {
    s
  } else {
    s + ".0"
  }
}

///|
/// 单页 HTML 仪表盘: 总览 + 每模块样本 x 指标矩阵(最优绿/最差红).
pub fn multiqc_write_html(report : MultiqcReport) -> String {
  let sb = StringBuilder()
  sb.write_string(
    "<!DOCTYPE html>\n<html lang=\"en\">\n<head>\n<meta charset=\"utf-8\">\n<title>",
  )
  sb.write_string(mq_html_escape(report.title))
  sb.write_string("</title>\n<style>\n")
  sb.write_string(
    "body{font-family:-apple-system,'Segoe UI',Roboto,sans-serif;margin:0;background:#f5f6fa;color:#2d3436;}\n",
  )
  sb.write_string("h1{margin:0;padding:24px 0 8px;font-size:22px;}\n")
  sb.write_string(
    ".wrap{max-width:1080px;margin:0 auto;padding:0 20px 40px;}\n",
  )
  sb.write_string(
    ".cards{display:flex;flex-wrap:wrap;gap:14px;margin:18px 0;}\n",
  )
  sb.write_string(
    ".card{background:#fff;border-radius:8px;padding:14px 18px;box-shadow:0 1px 3px rgba(0,0,0,.08);min-width:150px;}\n",
  )
  sb.write_string(".card .n{font-size:24px;font-weight:700;}\n")
  sb.write_string(".card .l{font-size:12px;color:#7f8c8d;margin-top:2px;}\n")
  sb.write_string(
    "h2{font-size:16px;margin:28px 0 10px;border-left:4px solid #0984e3;padding-left:10px;}\n",
  )
  sb.write_string(
    "table{border-collapse:collapse;width:100%;background:#fff;border-radius:8px;overflow:hidden;box-shadow:0 1px 3px rgba(0,0,0,.08);}\n",
  )
  sb.write_string("th,td{padding:8px 12px;text-align:right;font-size:13px;}\n")
  sb.write_string("th:first-child,td:first-child{text-align:left;}\n")
  sb.write_string("th{background:#dfe6e9;font-weight:600;}\n")
  sb.write_string("td.good{background:#d7f5dd;}\n")
  sb.write_string("td.bad{background:#fadbd8;}\n")
  sb.write_string(".foot{margin-top:24px;font-size:12px;color:#95a5a6;}\n")
  sb.write_string("</style>\n</head>\n<body>\n<div class=\"wrap\">\n")
  sb.write_string("<h1>")
  sb.write_string(mq_html_escape(report.title))
  sb.write_string("</h1>\n")
  // 总览卡片: 样本数 / 模块数
  let all_samples : Array[String] = []
  let sample_set : Map[String, Bool] = Map([], capacity=0)
  for m in report.modules {
    for s in m.samples.keys() {
      if sample_set.get(s) is None {
        sample_set.set(s, true)
        all_samples.push(s)
      }
    }
  }
  sb.write_string("<div class=\"cards\">\n")
  sb.write_string(
    "<div class=\"card\"><div class=\"n\">\{all_samples.length()}</div><div class=\"l\">samples</div></div>\n",
  )
  sb.write_string(
    "<div class=\"card\"><div class=\"n\">\{report.modules.length()}</div><div class=\"l\">modules</div></div>\n",
  )
  sb.write_string("</div>\n")
  // 每模块矩阵
  for m in report.modules {
    sb.write_string("<h2>")
    sb.write_string(mq_html_escape(m.name))
    sb.write_string("</h2>\n")
    // 收集列(所有样本指标名并集, 稳定排序)
    let col_set : Map[String, Bool] = Map([], capacity=0)
    for s in m.samples.keys() {
      let metrics = m.samples.get(s).unwrap()
      for k in metrics.keys() {
        col_set.set(k, true)
      }
    }
    let cols : Array[String] = []
    for k in col_set.keys() {
      cols.push(k)
    }
    cols.sort()
    let samples : Array[String] = []
    for s in m.samples.keys() {
      samples.push(s)
    }
    samples.sort()
    if cols.length() == 0 || samples.length() == 0 {
      sb.write_string("<p>no data</p>\n")
      continue
    }
    // 每列最优/最差
    let is_lower : Map[String, Bool] = Map([], capacity=0)
    for lb in m.lower_better {
      is_lower.set(lb, true)
    }
    sb.write_string("<table>\n<tr><th>sample</th>")
    for c in cols {
      sb.write_string("<th>")
      sb.write_string(mq_html_escape(c))
      sb.write_string("</th>")
    }
    sb.write_string("</tr>\n")
    for s in samples {
      let metrics = m.samples.get(s).unwrap()
      sb.write_string("<tr><td>")
      sb.write_string(mq_html_escape(s))
      sb.write_string("</td>")
      for c in cols {
        let v = metrics.get(c)
        sb.write_string("<td")
        match v {
          Some(x) => {
            // 列方向统计(样本数>1 才着色)
            if samples.length() > 1 {
              let mut mn = x
              let mut mx = x
              for s2 in samples {
                let m2 = m.samples.get(s2).unwrap()
                match m2.get(c) {
                  Some(y) => {
                    if y < mn {
                      mn = y
                    }
                    if y > mx {
                      mx = y
                    }
                  }
                  None => ()
                }
              }
              let lower = is_lower.get(c) is Some(_)
              let is_best = if lower { x == mn } else { x == mx }
              let is_worst = if lower { x == mx } else { x == mn }
              if is_best && mn != mx {
                sb.write_string(" class=\"good\"")
              } else if is_worst && mn != mx {
                sb.write_string(" class=\"bad\"")
              }
            }
            sb.write_string(">")
            sb.write_string(mq_fmt(x))
          }
          None => sb.write_string(">-")
        }
        sb.write_string("</td>")
      }
      sb.write_string("</tr>\n")
    }
    sb.write_string("</table>\n")
  }
  sb.write_string(
    "<div class=\"foot\">generated by BioSeqs multiqc-lite (local, no data leaves this machine)</div>\n",
  )
  sb.write_string("</div>\n</body>\n</html>\n")
  sb.to_string()
}

// ---------------------------------------------------------------------------
// Markdown 表格
// ---------------------------------------------------------------------------

///|
/// Markdown 表格输出(每模块一节).
pub fn multiqc_write_markdown(report : MultiqcReport) -> String {
  let sb = StringBuilder()
  sb.write_string("# ")
  sb.write_string(report.title)
  sb.write_string("\n\n")
  for m in report.modules {
    sb.write_string("## ")
    sb.write_string(m.name)
    sb.write_string("\n\n")
    let col_set : Map[String, Bool] = Map([], capacity=0)
    for s in m.samples.keys() {
      let metrics = m.samples.get(s).unwrap()
      for k in metrics.keys() {
        col_set.set(k, true)
      }
    }
    let cols : Array[String] = []
    for k in col_set.keys() {
      cols.push(k)
    }
    cols.sort()
    let samples : Array[String] = []
    for s in m.samples.keys() {
      samples.push(s)
    }
    samples.sort()
    if cols.length() == 0 || samples.length() == 0 {
      sb.write_string("_no data_\n\n")
      continue
    }
    sb.write_string("| sample |")
    for c in cols {
      sb.write_string(" ")
      sb.write_string(c)
      sb.write_string(" |")
    }
    sb.write_string("\n|---|")
    for c in cols {
      ignore(c)
      sb.write_string("---|")
    }
    sb.write_string("\n")
    for s in samples {
      let metrics = m.samples.get(s).unwrap()
      sb.write_string("| ")
      sb.write_string(s)
      sb.write_string(" |")
      for c in cols {
        match metrics.get(c) {
          Some(x) => {
            sb.write_string(" ")
            sb.write_string(mq_fmt(x))
            sb.write_string(" |")
          }
          None => sb.write_string(" - |")
        }
      }
      sb.write_string("\n")
    }
    sb.write_string("\n")
  }
  sb.to_string()
}

///|
/// 合并同名模块(多样本聚合语义, 对齐 MultiQC 模块合并).
pub fn MultiqcReport::merge_modules(self : MultiqcReport) -> MultiqcReport {
  let merged : Map[String, MultiqcModule] = Map([], capacity=0)
  let order : Array[String] = Array(capacity=0)
  for m in self.modules {
    match merged.get(m.name) {
      Some(acc) =>
        for s in m.samples.keys() {
          let metrics = m.samples.get(s).unwrap()
          let pairs : Array[(String, Double)] = Array(capacity=0)
          for k in metrics.keys() {
            pairs.push((k, metrics.get(k).unwrap()))
          }
          let _ = acc.add_sample(s, pairs)
        }
      None => {
        let fresh = multiqc_module(m.name)
        for s in m.samples.keys() {
          let metrics = m.samples.get(s).unwrap()
          let pairs : Array[(String, Double)] = Array(capacity=0)
          for k in metrics.keys() {
            pairs.push((k, metrics.get(k).unwrap()))
          }
          let _ = fresh.add_sample(s, pairs)
        }
        let _ = fresh.with_lower_better(m.lower_better)
        merged.set(m.name, fresh)
        order.push(m.name)
      }
    }
  }
  let modules : Array[MultiqcModule] = Array(capacity=0)
  for name in order {
    modules.push(merged.get(name).unwrap())
  }
  { title: self.title, modules, }
}

///|
/// 报告摘要行.
pub fn MultiqcReport::summary(self : MultiqcReport) -> String {
  let mut n_samples : Int32 = 0
  let sample_set : Map[String, Bool] = Map([], capacity=0)
  for m in self.modules {
    for s in m.samples.keys() {
      if sample_set.get(s) is None {
        sample_set.set(s, true)
        n_samples = n_samples + 1
      }
    }
  }
  "multiqc-lite: \{self.modules.length()} modules, \{n_samples} samples, title=\{self.title}"
}
