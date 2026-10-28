///|
/// primer_web - 浏览器端引物设计工具(MoonBit 计算侧).
///
/// 对模板序列运行 primer3 设计引擎并输出 JSON(report.json),
/// 由同目录 index.html 渲染候选表/引物对产物视图.
/// 生成: moon run examples/primer_web -- "$(cat template.txt)" > report.json
/// 演示: moon run examples/primer_web -- --demo > report.json
fn main raise {
  let argv : Array[String] = []
  for i, a in @env.args() {
    let is_prog = i == 0 && a.find("/") is Some(_)
    if a.length() > 0 && !is_prog {
      argv.push(a)
    }
  }
  let template = if argv.length() == 0 || argv[0] == "--demo" {
    @src.primer3_sample_template()
  } else {
    argv.join("")
  }
  // 盐校正 Tm 公式(81.5 + 16.6·log10[Na⁺] + 0.41·GC% - 600/n)在 50mM
  // 下基线约 59.8, 默认 55-65 区间过窄; 演示采用与库测试一致的宽容配置.
  let cfg = @src.Primer3Config::new(
    18, 22, 40.0, 75.0, 0.3, 0.8, 8, 8, 0.05, 50,
  )
  println(@src.primer3_design_report_json(template, config=cfg))
}
