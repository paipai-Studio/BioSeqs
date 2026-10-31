///|
/// vg_lite_demo — vg 风格泛基因组变异图基因分型演示:
/// 参考 + VCF 变异集 -> 变异图(节点/边) -> 单倍型路径枚举 ->
/// 读段滑窗 Hamming 分派 -> 二倍体 GT -> GFA/JSON.

///|
fn vgd_reference() -> String {
  let sb = StringBuilder()
  for _i = 0; _i < 10; _i = _i + 1 {
    sb.write_string("ACGT")
  }
  sb.to_string()
}

///|
fn vgd_vcf() -> String {
  "##fileformat=VCFv4.2\n" +
  "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n" +
  "chr1\t5\t.\tA\tG\t.\tPASS\t.\n" +
  "chr1\t21\t.\tA\tC\t.\tPASS\t.\n" +
  "chr1\t33\t.\tT\tTA\t.\tPASS\t.\n"
}

///|
fn main {
  let reference = vgd_reference()
  let vcf = @src.parse_vcf(vgd_vcf())
  println("1. 输入: 40bp 参考 + 3 双等位变异(SNV x2 + 插入 x1)")
  // 2. 构图
  let graph = @src.vg_build_graph("chr1", reference, vcf.records)
  println(
    "2. 变异图: \{graph.nodes.length()} 节点 / \{graph.edges.length()} 边",
  )
  for n in graph.nodes {
    println("  节点 \{n.id} [\{n.kind}] \{n.seq.length()}bp")
  }
  // 3. 路径枚举
  let paths = @src.vg_paths(graph)
  println("3. 单倍型路径 \{paths.length()} 条:")
  for p in paths {
    let mut tag = ""
    for id in p.node_ids {
      tag = tag + "\{id} "
    }
    println("  [\{tag}] \{p.seq.length()}bp")
  }
  // 4. 分型(杂合样本: 3 ref + 3 变异1 读段)
  let reads = [
    "ACGTACGTAC", "ACGTACGTAC", "ACGTACGTAC", "ACGTGCGTAC", "ACGTGCGTAC", "ACGTGCGTAC",
  ]
  let gt = @src.vg_genotype(paths, reads)
  let g = gt.unwrap()
  println("4. 基因分型(6 条读段):")
  println("  GT = \{g.gt}")
  println("  等位支持 = \{g.counts[0]} / \{g.counts[1]}")
  let mut pc = ""
  for c in g.path_counts {
    pc = pc + "\{c} "
  }
  println("  全部路径计数: \{pc}")
  // 5. GFA / JSON
  let gfa = @src.vg_to_gfa(graph)
  let json = @src.json_to_string(@src.vg_value(graph, paths, gt))
  println("5. GFA \{gfa.length()} 字节, JSON \{json.length()} 字节")
  println(@src.vg_summary_line(graph, paths))
  println(
    "演示完成: VCF -> 变异图 -> 单倍型路径 -> 图上分型闭环",
  )
}
