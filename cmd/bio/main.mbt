///|
/// bio - BioSeqs 命令行工具(seqkit 风格).
///
/// 用法: bio <子命令> [FASTA/FASTQ 文本 | --demo]
/// 注: MoonBit core 无文件系统接口, 序列文本经参数传入;
/// 嵌入宿主程序读取文件后调用库 API 可获得完整文件处理能力.
/// 示例: moon run cmd/bio -- stats "$(cat seqs.fasta)"
///       moon run cmd/bio -- fq-qc --demo
fn bio_demo_fasta() -> String {
  ">seq1 kinase-like\nMKTLLLTLVVMKTAYIAKQRQISFVKSHFSRQLEER\n>seq2 promoter-region\nGCGCGTTTAAACGCGCGTTTAAACGGCCGGTT\n>seq3 at-rich\nAATTAATTAATTGCAAATTAA\n"
}

///|
fn bio_demo_fastq() -> String {
  "@r1 first read\nACGTACGTACGT\n+\nIIIIIIIIIIII\n@r2\nGGGGCCCCGGGG\n+\nFFFFFFFFFFFF\n@r3\nTTTTAAAATTTTAAAA\n+\n!!!!!!!!!!!!!!!!\n"
}

///|
fn bio_print_help() -> Unit {
  println("bio - BioSeqs 命令行工具 (seqkit 风格)")
  println("")
  println("用法: bio <子命令> [FASTA/FASTQ 文本 | --demo]")
  println("")
  println("子命令:")
  println("  gc      序列 GC 含量 (FASTA)")
  println("  stats   FASTA 统计: 条数/总长/均值/N50/最长")
  println(
    "  fq-qc   FASTQ 基础质控: reads 数/读长分布/质量均值/GC",
  )
  println("  rc      反向互补 (FASTA 输入输出)")
  println("  trans   翻译为蛋白 (FASTA -> FASTA, 框 1, NCBI 表 1)")
  println("  info    每条记录摘要: id/长度/GC%/分子量")
  println("")
  println("选项:")
  println("  --demo  使用内置示例数据")
  println("  -h      显示本帮助")
}

///|
fn bio_input_text(args : Array[String]) -> String? {
  let parts : Array[String] = []
  for a in args {
    if a == "--demo" {
      return Some(bio_demo_fasta())
    }
    if a.length() > 0 {
      parts.push(a)
    }
  }
  if parts.length() == 0 {
    return None
  }
  Some(parts.join("\n"))
}

///|
fn bio_require_input(cmd : String, args : Array[String]) -> String? {
  match bio_input_text(args) {
    Some(text) => Some(text)
    None => {
      println("错误: 子命令 '" + cmd + "' 需要序列文本或 --demo")
      None
    }
  }
}

///|
fn bio_is_fastq(cmd : String, args : Array[String]) -> Bool {
  for a in args {
    if a == "--demo" {
      return cmd == "fq-qc"
    }
  }
  false
}

///|
fn bio_pct(v : Double) -> String {
  (v * 100.0).to_string() + "%"
}

///|
fn bio_n50(lengths : Array[Int32]) -> Int32 {
  if lengths.length() == 0 {
    return 0
  }
  let sorted = lengths.copy()
  let mut i = 0
  while i < sorted.length() {
    let mut j = i + 1
    while j < sorted.length() {
      if sorted[j] > sorted[i] {
        let t = sorted[i]
        sorted[i] = sorted[j]
        sorted[j] = t
      }
      j = j + 1
    }
    i = i + 1
  }
  let total = {
    let mut s = 0
    for l in sorted {
      s = s + l
    }
    s
  }
  let mut acc = 0
  for l in sorted {
    acc = acc + l
    if acc * 2 >= total {
      return l
    }
  }
  0
}

///|
fn bio_cmd_gc(args : Array[String]) -> Unit {
  let input = match bio_require_input("gc", args) {
    Some(t) => t
    None => return
  }
  let records = @lib.parse_fasta(input)
  if records.length() == 0 {
    println("(无记录)")
    return
  }
  let mut total_len = 0
  let mut total_gc = 0.0
  for rec in records {
    let s = rec.seq.to_string()
    let gc = @lib.gc_content(s)
    total_len = total_len + s.length()
    total_gc = total_gc + gc * s.length().to_double()
    println("  " + rec.id + "\t" + bio_pct(gc))
  }
  let overall = if total_len > 0 {
    total_gc / total_len.to_double()
  } else {
    0.0
  }
  println("  合计\t" + bio_pct(overall))
}

///|
fn bio_cmd_stats(args : Array[String]) -> Unit raise @lib.SeqIOError {
  let input = match bio_require_input("stats", args) {
    Some(t) => t
    None => return
  }
  let lengths : Array[Int32] = []
  @lib.fasta_each_record(input, fn(rec) { lengths.push(rec.seq.length()) })
  if lengths.length() == 0 {
    println("(无记录)")
    return
  }
  let mut total = 0
  let mut max_len = 0
  let mut min_len = lengths[0]
  for l in lengths {
    total = total + l
    if l > max_len {
      max_len = l
    }
    if l < min_len {
      min_len = l
    }
  }
  let mean = total.to_double() / lengths.length().to_double()
  println("  条数:    " + lengths.length().to_string())
  println("  总长:    " + total.to_string())
  println("  均长:    " + mean.to_string())
  println("  最短:    " + min_len.to_string())
  println("  最长:    " + max_len.to_string())
  println("  N50:     " + bio_n50(lengths).to_string())
}

///|
fn bio_cmd_fq_qc(args : Array[String]) -> Unit raise @lib.SeqIOError {
  let input = if bio_is_fastq("fq-qc", args) {
    bio_demo_fastq()
  } else {
    match bio_require_input("fq-qc", args) {
      Some(t) => t
      None => return
    }
  }
  let mut n_reads = 0
  let mut total_bases = 0
  let mut q_sum = 0
  let mut q_count = 0
  let mut gc_bases = 0
  let mut min_len = -1
  let mut max_len = -1
  @lib.fastq_each_record(input, fn(rec) {
    let s = rec.seq.to_string()
    n_reads = n_reads + 1
    total_bases = total_bases + s.length()
    for i in 0..<s.length() {
      let c = s.unsafe_get(i)
      if c == 'G' || c == 'C' || c == 'g' || c == 'c' {
        gc_bases = gc_bases + 1
      }
    }
    match rec.letter_annotations.get("phred_quality") {
      Some(quals) =>
        for q in quals {
          q_sum = q_sum + q
          q_count = q_count + 1
        }
      None => ()
    }
    if min_len < 0 || s.length() < min_len {
      min_len = s.length()
    }
    if max_len < 0 || s.length() > max_len {
      max_len = s.length()
    }
  })
  if n_reads == 0 {
    println("(无 reads)")
    return
  }
  let mean_len = total_bases.to_double() / n_reads.to_double()
  let mean_q = if q_count > 0 {
    q_sum.to_double() / q_count.to_double()
  } else {
    0.0
  }
  let gc = if total_bases > 0 {
    gc_bases.to_double() / total_bases.to_double()
  } else {
    0.0
  }
  println("  reads:   " + n_reads.to_string())
  println("  碱基:    " + total_bases.to_string())
  println(
    "  读长:    min=" +
    min_len.to_string() +
    " mean=" +
    mean_len.to_string() +
    " max=" +
    max_len.to_string(),
  )
  println("  平均Q:   " + mean_q.to_string())
  println("  GC:      " + bio_pct(gc))
}

///|
fn bio_cmd_rc(args : Array[String]) -> Unit {
  let input = match bio_require_input("rc", args) {
    Some(t) => t
    None => return
  }
  let records = @lib.parse_fasta(input)
  for rec in records {
    let rc = @lib.Seq::new(rec.seq.to_string()).reverse_complement().to_string()
    println(">" + rec.description + "_rc")
    println(rc)
  }
}

///|
fn bio_cmd_trans(args : Array[String]) -> Unit {
  let input = match bio_require_input("trans", args) {
    Some(t) => t
    None => return
  }
  let records = @lib.parse_fasta(input)
  for rec in records {
    let protein : String? = Some(
      @lib.Seq::new(rec.seq.to_string()).translate().to_string(),
    ) catch {
      SeqError(msg) => {
        println(">" + rec.description + "_protein [翻译失败: " + msg + "]")
        None
      }
    }
    match protein {
      Some(p) => {
        println(">" + rec.description + "_protein")
        println(p)
      }
      None => ()
    }
  }
}

///|
fn bio_cmd_info(args : Array[String]) -> Unit {
  let input = match bio_require_input("info", args) {
    Some(t) => t
    None => return
  }
  let records = @lib.parse_fasta(input)
  println("  id\t长度\tGC%\t分子量(Da)")
  for rec in records {
    let s = rec.seq.to_string()
    let gc = @lib.gc_content(s)
    let mw = @lib.sequtils_molecular_weight(s, "DNA")
    println(
      "  " +
      rec.id +
      "\t" +
      s.length().to_string() +
      "\t" +
      bio_pct(gc) +
      "\t" +
      mw.to_string(),
    )
  }
}

///|
fn main raise {
  let argv = @env.args()
  let rest : Array[String] = []
  for i, a in argv {
    // argv[0] 为程序路径/程序名(moon run 下为 .wasm 路径), 跳过
    let is_prog = i == 0 && (a == "bio" || a.find("/") is Some(_))
    if a.length() > 0 && !is_prog {
      rest.push(a)
    }
  }
  if rest.length() == 0 || rest[0] == "-h" || rest[0] == "--help" {
    bio_print_help()
    return
  }
  let cmd = rest[0]
  let args : Array[String] = []
  for i in 1..<rest.length() {
    args.push(rest[i])
  }
  try {
    match cmd {
      "gc" => bio_cmd_gc(args)
      "stats" => bio_cmd_stats(args)
      "fq-qc" => bio_cmd_fq_qc(args)
      "rc" => bio_cmd_rc(args)
      "trans" => bio_cmd_trans(args)
      "info" => bio_cmd_info(args)
      _ => {
        println("未知子命令: " + cmd)
        bio_print_help()
      }
    }
  } catch {
    SeqIOError(msg) => println("解析错误: " + msg)
  }
}
