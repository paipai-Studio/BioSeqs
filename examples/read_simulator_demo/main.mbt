///|
/// 读段模拟器演示: 参考序列 -> 双端 FASTQ.
fn main raise {
  println(
    "=== 本地读段模拟器 (参考 -> 双端 FASTQ, wgsim 风格) ===",
  )

  // 构造 300bp 模拟参考.
  let sb = StringBuilder()
  for i = 0; i < 300; i = i + 1 {
    sb.write_char(['A', 'C', 'G', 'T'][(i * 5 + 1) % 4])
  }
  let reference = sb.to_string()

  let config = @src.read_simulator_config(
    read_length=100,
    coverage=20.0,
    paired_end=true,
    insert_mean=250,
    insert_sd=30,
    error_rate=0.005,
    quality_max=40,
    quality_min=2,
    quality_decay=12.0,
    snp_rate=0.005,
    indel_rate=0.001,
    adapter="AGATCGGAAGAG",
    adapter_rate=0.0,
    seed=2026UL,
  )
  let result = @src.simulate_reads(reference, config~)

  println(
    "\n参考长度: \{reference.length()}bp, 目标覆盖度: \{config.coverage}x",
  )
  println("生成 pairs = \{result.n_pairs}, reads = \{result.n_reads}")
  println("实际覆盖度 = \{result.achieved_coverage}x")
  println(
    "注入变异事件 = \{result.n_variants}, 碱基错误 = \{result.n_errors}",
  )

  println("\nFASTQ 前 2 条记录:")
  let lines = result.fastq.split("\n").to_array()
  for i = 0; i < 8 && i < lines.length(); i = i + 1 {
    println("  \{lines[i].to_owned()}")
  }

  println("\n完成: 读段模拟全程本地完成, 数据不出本机.")
}
