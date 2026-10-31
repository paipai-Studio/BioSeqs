///|
/// msprime_demo — msprime 风格群体谱系模拟演示:
/// 时变 Ne 扩张模型 + 重组多块 -> SFS -> Tajima's D / Watterson theta ->
/// JSON / tree sequence.

///|
fn msd_demes() -> String {
  "# demes model: bottleneck\n" +
  "- start_time: 0\n  end_time: 500\n  ne: 1000\n" +
  "- start_time: 500\n  end_time: 10000\n  ne: 10000\n"
}

///|
fn main {
  println("1. msprime-lite 群体谱系模拟演示(时变 Ne + 重组 + SFS)")
  // 2. demes 解析
  let epochs = @src.msprime_parse_demes(msd_demes())
  println("2. demes 解析: \{epochs.length()} 个 epoch")
  for e in epochs {
    println("  [\{e.start_time}, \{e.end_time}) Ne=\{e.ne}")
  }
  // 3. 模拟
  let params : @src.MsprimeParams = {
    n_samples: 20,
    n_pops: 1,
    epochs,
    migration: None,
    recombination_rate: 1.0e-8,
    sequence_length: 10000,
    mutation_rate: 1.0e-8,
    n_blocks: 4,
  }
  let rng = @src.Rng::new(20261031UL)
  let result = @src.msprime_simulate(params, rng)
  println(
    "3. 模拟: \{result.n_mutations}) \{result.sequence_length}bp, \{result.blocks.length()} 块",
  )
  // 4. SFS 统计
  println(
    "4. 统计: Tajima's D = \{result.tajima_d}, Watterson theta = \{result.watterson_theta}",
  )
  // 5. 折叠 SFS 前 5 项
  println("5. 折叠 SFS (前 5 项):")
  for i = 0; i < 5 && i < result.folded_sfs.length(); i = i + 1 {
    println("  k=\{i}: \{result.folded_sfs[i]}")
  }
  // 6. 导出
  let json = @src.msprime_json(result)
  let ts = @src.msprime_tree_sequence(result)
  println(
    "6. JSON \{json.length()} 字节, tree sequence \{ts.length()} 字节",
  )
  println(@src.msprime_summary_line(result))
  println(
    "演示完成: demes 解析 -> 时变 Ne 模拟 -> 重组多块 -> SFS -> 统计量闭环",
  )
}
