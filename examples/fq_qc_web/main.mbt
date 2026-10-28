///|
/// fq_qc_web - 浏览器端 FASTQ 质控仪表盘(MoonBit 计算侧).
///
/// 计算全部 FastQC 风格指标并输出 JSON(report.json),
/// 由同目录 index.html 渲染为纯 SVG 仪表盘(无外部依赖).
/// 生成方式: moon run examples/fq_qc_web -- "$(cat sample.fq)" > report.json
/// 演示:    moon run examples/fq_qc_web -- --demo > report.json
fn main raise {
  let argv : Array[String] = []
  for i, a in @env.args() {
    let is_prog = i == 0 && a.find("/") is Some(_)
    if a.length() > 0 && !is_prog {
      argv.push(a)
    }
  }
  let fastq = if argv.length() == 0 || argv[0] == "--demo" {
    fq_qc_web_demo()
  } else {
    argv.join("\n")
  }
  let report = @src.fqqc_report(fastq)
  println(report)
}

///|
/// 生成演示数据: 300 reads, 质量随位置衰减, 含 adapter 污染与重复序列.
fn fq_qc_web_demo() -> String {
  let rng = @src.Rng::from_int(42)
  let bases = ['A', 'C', 'G', 'T']
  let adapter = "AGATCGGAAGAGC"
  let sb = StringBuilder()
  for i in 0..<300 {
    let len = 75 + rng.next_int(26)
    let with_adapter = i % 37 == 0
    let body_len = if with_adapter && len > adapter.length() {
      len - adapter.length()
    } else {
      len
    }
    let seq = StringBuilder()
    // GC 偏置: 前段 GC 富集
    for j in 0..<body_len {
      let bias = if j < body_len / 2 {
        rng.next_int(10) < 6
      } else {
        rng.next_int(10) < 4
      }
      let c = if bias {
        bases[2 + rng.next_int(2)]
      } else {
        bases[rng.next_int(2)]
      }
      seq.write_char(c)
    }
    if with_adapter {
      seq.write_string(adapter)
    }
    // 质量: 均值 36 -> 22 随位置指数衰减 + 噪声
    let qual = StringBuilder()
    for j in 0..<len {
      let decay = 36.0 - 14.0 * (j.to_double() / len.to_double()).sqrt()
      let q = rng.normal(decay, 3.0)
      let q = if q > 40.0 { 40.0 } else if q < 2.0 { 2.0 } else { q }
      qual.write_char(Char::from_int(33 + q.to_int()))
    }
    sb.write_string("@demo_read_")
    sb.write_string(i.to_string())
    if with_adapter {
      sb.write_string(" adapter-contaminated")
    }
    sb.write_char('\n')
    sb.write_string(seq.to_string())
    sb.write_char('\n')
    sb.write_string("+\n")
    sb.write_string(qual.to_string())
    sb.write_char('\n')
  }
  sb.to_string()
}
