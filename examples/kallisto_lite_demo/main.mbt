///|
/// kallisto-lite 伪比对定量演示: 索引 -> 等价类 -> EM -> TPM -> quant.sf -> tximport 闭环.
fn main raise {
  println("=== kallisto-lite 伪比对转录本定量 ===")

  // 模拟转录组: t1/t2 共享中段, t3 独特
  let tx_names : Array[String] = ["t1", "t2", "t3"]
  let tx_seqs : Array[String] = [
    "AAAAACCCCCGGGGGTTTTT" + "ATATGCGCGCTATAT" + "AACCTTGGAACCTTGGAA",
    "TTTTTGGGGGCCCCCAAAAA" + "ATATGCGCGCTATAT" + "GGTTCCAAGGTTCCAAGG",
    "ACACACGTGTGTATATATCGCGCGTATATCGCGCGTATAT",
  ]

  // 模拟 read 集: t1 独特区 x8, 共享区 x4, t3 独特区 x6
  let reads : Array[String] = []
  for i = 0; i < 8; i = i + 1 {
    reads.push("AAAAACCCCCGGGGGTTTTTATAT")
  }
  for i = 0; i < 4; i = i + 1 {
    reads.push("ATATGCGCGCTATAT")
  }
  for i = 0; i < 6; i = i + 1 {
    reads.push("ACACACGTGTGTATATATCGCGCG")
  }
  println(
    "转录本 \{tx_names.length()} 条, reads \{reads.length()} 条 (k=9)",
  )

  // 1. 建索引
  let index = @src.kallisto_build_index(tx_names, tx_seqs, k=9)
  println("\n1. 索引: 有效长度 " + repr(index.eff_lens))

  // 2. 等价类
  let ecs = @src.kallisto_pseudoalign(index, reads)
  println("\n2. 等价类 (\{ecs.length()} 个):")
  for ec in ecs {
    let names : Array[String] = []
    for id in ec.tx_ids {
      names.push(tx_names[id])
    }
    println("  {" + names.join(",") + "} count=" + ec.count.to_string())
  }

  // 3. EM 定量
  let quant = @src.kallisto_quantify(index, reads)
  println("\n3. EM 定量:")
  println("  tx   est_counts   TPM")
  for i = 0; i < tx_names.length(); i = i + 1 {
    println(
      "  " +
      tx_names[i] +
      "  " +
      quant.est_counts[i].to_string() +
      "        " +
      quant.tpm[i].to_string(),
    )
  }
  println("  " + quant.summary())

  // 4. quant.sf 输出与 tximport 闭环
  let sf = quant.to_quant_sf()
  println("\n4. quant.sf:")
  let lines = sf.split("\n").to_array()
  for l in lines {
    if l.to_owned().length() > 0 {
      println("  " + l.to_owned())
    }
  }
  let (parsed_names, abundance, counts, lengths) = @src.parse_quant_sf(sf)
  println(
    "\n5. tximport.parse_quant_sf 闭环: \{parsed_names.length()} 条转录本, " +
    "counts 总和 \{counts.fold(init=0.0, fn(acc, x) { acc + x }).to_string()}, " +
    "abundance(TPM) 总和 \{abundance.fold(init=0.0, fn(acc, x) { acc + x }).to_string()}, " +
    "lengths " +
    repr(lengths),
  )
  println(
    "\n完成: 伪比对 -> EM -> TPM -> quant.sf -> tximport 全链路本地完成.",
  )
}
