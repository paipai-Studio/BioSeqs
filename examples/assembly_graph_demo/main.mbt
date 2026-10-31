///|
/// 组装图简化与共识抛光演示:
/// GFA 图 -> tip 裁剪 / containment 坍缩 / bubble 解析 -> 统计,
/// 参考序列 + 模拟读段 -> racon-lite 抛光 -> 修正报告.
fn demo_repeat(base : Char, len : Int32) -> String {
  let sb = StringBuilder()
  for _i in 0..<len {
    sb.write_char(base)
  }
  sb.to_string()
}

///|
fn demo_reference(len : Int32) -> String {
  let bases = "ACGT"
  let sb = StringBuilder()
  let rng = @src.Rng::new(20261030UL)
  for _i in 0..<len {
    let idx = rng.randint(4)
    sb.write_char(bases.unsafe_get(idx).unsafe_to_char())
  }
  sb.to_string()
}

///|
fn demo_diff_count(a : String, b : String) -> Int32 {
  let n = if a.length() < b.length() { a.length() } else { b.length() }
  let mut diff = 0
  let mut i = 0
  while i < n {
    if a.unsafe_get(i) != b.unsafe_get(i) {
      diff = diff + 1
    }
    i = i + 1
  }
  diff +
  (if a.length() > b.length() {
    a.length() - b.length()
  } else {
    b.length() - a.length()
  })
}

///|
fn demo_has_segment(graph : @src.GfaGraph, name : String) -> Bool {
  match @src.gfa_get_segment(graph, name) {
    Some(_) => true
    None => false
  }
}

///|
fn main {
  println("=== 组装图简化与共识抛光 (racon-lite) ===")

  // ---- 1. 构造含短死端与 bubble 的 GFA 图 ----
  let unit = demo_repeat('A', 120)
  let gfa_text = "H\tVN:Z:1.0\n" +
    "S\tentry\t" +
    unit +
    "\n" +
    "S\ta1\t" +
    unit +
    "\tDP:i:10\n" +
    "S\ta2\t" +
    unit +
    "\tDP:i:10\n" +
    "S\tb1\t" +
    unit +
    "\tDP:i:2\n" +
    "S\texit\t" +
    unit +
    "\n" +
    "S\ttip\t" +
    demo_repeat('T', 20) +
    "\n" +
    "L\tentry\t+\ta1\t+\t0M\n" +
    "L\ta1\t+\ta2\t+\t0M\n" +
    "L\ta2\t+\texit\t+\t0M\n" +
    "L\tentry\t+\tb1\t+\t0M\n" +
    "L\tb1\t+\texit\t+\t0M\n" +
    "L\tentry\t+\ttip\t+\t0M\n" +
    "P\tp_tip\tentry+,a1+,a2+,exit+,tip+\t0M,0M,0M,0M\n" +
    "P\tp_main\tentry+,a1+,a2+,exit+\t0M,0M,0M\n"
  let graph = @src.gfa_parse(gfa_text)
  println("\n原始图: \{@src.gfa_graph_summary(graph)}")
  let bubbles = @src.assembly_detect_bubbles(graph)
  println(
    "原始图直接检测 bubble: \{bubbles.length()} (tip 使 entry 出度为 3, 需先裁剪 tip)",
  )

  // ---- 2. 简化 ----
  let simplified = @src.assembly_simplify(graph)
  let stats = simplified.stats()
  println("\n简化统计 JSON: \{@src.assembly_simplify_stats_to_json(stats)}")
  println("简化后图: \{@src.gfa_graph_summary(simplified.graph())}")
  let keep_a1 = demo_has_segment(simplified.graph(), "a1")
  let keep_b1 = demo_has_segment(simplified.graph(), "b1")
  let keep_tip = demo_has_segment(simplified.graph(), "tip")
  println("保留 a1: \{keep_a1}")
  println("保留 b1: \{keep_b1}")
  println("保留 tip: \{keep_tip}")
  println("保留 path 数: \{@src.gfa_graph_n_paths(simplified.graph())}")

  // ---- 3. 模拟读段并抛光 ----
  let reference = demo_reference(400)
  println("\n参考序列长度: \{reference.length()}")
  let mut sim_fastq = ""
  let mut n_reads = 0
  try {
    let sim_config = @src.read_simulator_config(
      read_length=100,
      coverage=30.0,
      paired_end=false,
      error_rate=0.02,
      snp_rate=0.01,
      indel_rate=0.003,
      seed=7UL,
    )
    let sim = @src.simulate_reads(reference, config=sim_config)
    sim_fastq = sim.fastq
    n_reads = sim.n_reads
    println(
      "模拟读段: \{sim.n_reads} 条, 实际覆盖度=\{sim.achieved_coverage}, 错误碱基数=\{sim.n_errors}",
    )
  } catch {
    @src.ReadSimulatorError(message) =>
      println("读段模拟失败: \{message}")
  }
  if sim_fastq == "" {
    return
  }
  let reads : Array[String] = Array(capacity=n_reads)
  try {
    let records = @src.parse_fastq(sim_fastq)
    for record in records {
      reads.push(record.seq.to_string())
    }
  } catch {
    _ => println("FASTQ 解析失败")
  }
  let mut polished_done = false
  try {
    let config = @src.polish_config(k=15, min_coverage=3, min_fraction=0.6)
    let result = @src.racon_polish(reference, reads, config~)
    println(
      "\n抛光结果: reads_used=\{result.reads_used()}/\{result.reads_total()}",
    )
    println("平均覆盖度: \{result.mean_coverage()}")
    println("采纳修正数: \{result.edits().length()}")
    let edits = result.edits()
    let shown = if edits.length() < 8 { edits.length() } else { 8 }
    for i in 0..<shown {
      let edit = edits[i]
      println(
        "  \{edit.kind()} @\{edit.position()} \{edit.ref_base()} -> \{edit.new_base()} (support=\{edit.support()}/\{edit.coverage()})",
      )
    }
    let mut n_sub = 0
    let mut n_del = 0
    let mut n_ins = 0
    for edit in result.edits() {
      if edit.kind() == "substitution" {
        n_sub = n_sub + 1
      } else if edit.kind() == "deletion" {
        n_del = n_del + 1
      } else {
        n_ins = n_ins + 1
      }
    }
    println(
      "修正分类: 替换=\{n_sub} 删除=\{n_del} 插入=\{n_ins}; consensus 长度=\{result.consensus().length()} (参考 \{reference.length()})",
    )
    ignore(demo_diff_count(result.consensus(), reference))
    println("JSON: \{@src.polish_result_to_json(result)}")
    polished_done = true
  } catch {
    @src.AssemblyGraphError(message) => println("抛光失败: \{message}")
  }
  if polished_done {
    println("\n演示完成: 图简化与读段共识抛光闭环")
  }
}
