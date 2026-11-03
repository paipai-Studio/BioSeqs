// ===========================================================================
// vg_lite_test.mbt — vg 风格变异图基因分型测试.
// ===========================================================================

///|
/// 40bp 参考("ACGT" x10): pos 5 与 pos 21 均为 'A'.
fn vg_test_reference() -> String {
  let sb = StringBuilder()
  for _i = 0; _i < 10; _i = _i + 1 {
    sb.write_string("ACGT")
  }
  sb.to_string()
}

///|
fn vg_test_vcf() -> String {
  "##fileformat=VCFv4.2\n" +
  "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n" +
  "chr1\t5\t.\tA\tG\t.\tPASS\t.\n" +
  "chr1\t21\t.\tA\tC\t.\tPASS\t.\n"
}

///|
test "vg: 变异图构造" {
  let vcf = @src.parse_vcf(vg_test_vcf())
  let graph = @src.vg_build_graph("chr1", vg_test_reference(), vcf.records)
  // 节点: [0,4)ref + A(ref) + G(alt) + [5,20)ref + A(ref) + C(alt) + [21,40)ref
  inspect(graph.nodes.length(), content="7")
  inspect(graph.nodes[0].seq, content="ACGT")
  inspect(graph.nodes[0].kind, content="ref")
  inspect(graph.nodes[1].seq, content="A")
  inspect(graph.nodes[1].kind, content="ref")
  inspect(graph.nodes[2].seq, content="G")
  inspect(graph.nodes[2].kind, content="alt")
  inspect(graph.nodes[3].seq.length(), content="15")
  inspect(graph.nodes[4].seq, content="A")
  inspect(graph.nodes[5].seq, content="C")
  inspect(graph.nodes[5].kind, content="alt")
  inspect(graph.nodes[6].seq.length(), content="19")
  // 边: 0->1, 0->2, 1->3, 2->3, 3->4, 3->5, 4->6, 5->6
  inspect(graph.edges.length(), content="8")
  let line = @src.vg_summary_line(graph, [])
  assert_eq(line.contains("nodes=7"), true)
  assert_eq(line.contains("edges=8"), true)
}

///|
test "vg: 路径枚举(4 单倍型)" {
  let vcf = @src.parse_vcf(vg_test_vcf())
  let graph = @src.vg_build_graph("chr1", vg_test_reference(), vcf.records)
  let paths = @src.vg_paths(graph)
  inspect(paths.length(), content="4")
  // 每条路径 = 40bp 单倍型
  for p in paths {
    inspect(p.seq.length(), content="40")
  }
  // 路径拼接: ref/ref 路径 0 与参考一致
  assert_eq(paths[0].seq, vg_test_reference())
  // alt(变异1) 路径: pos5 G
  assert_eq(paths[2].seq.unsafe_get(4).unsafe_to_char(), 'G')
  // alt(变异2) 路径: pos21 C
  assert_eq(paths[1].seq.unsafe_get(20).unsafe_to_char(), 'C')
  // 双 alt 路径
  assert_eq(paths[3].seq.unsafe_get(4).unsafe_to_char(), 'G')
  assert_eq(paths[3].seq.unsafe_get(20).unsafe_to_char(), 'C')
}

///|
test "vg: 杂合基因分型(0/2)" {
  let vcf = @src.parse_vcf(vg_test_vcf())
  let graph = @src.vg_build_graph("chr1", vg_test_reference(), vcf.records)
  let paths = @src.vg_paths(graph)
  // 2 条 ref 读段(参考 [0,10)) + 2 条变异1 读段(pos5 G)
  let reads = ["ACGTACGTAC", "ACGTACGTAC", "ACGTGCGTAC", "ACGTGCGTAC"]
  let gt = @src.vg_genotype(paths, reads)
  assert_true(gt is Some(_))
  let g = gt.unwrap()
  // 计数: 路径 0 = 2, 路径 2 = 2 -> 前二(降序, tie 索引小)
  inspect(g.gt, content="0/2")
  assert_eq(g.counts[0], 2)
  assert_eq(g.counts[1], 2)
  assert_eq(g.n_reads, 4)
}

///|
test "vg: 纯合基因分型(1/1)" {
  let vcf = @src.parse_vcf(vg_test_vcf())
  let graph = @src.vg_build_graph("chr1", vg_test_reference(), vcf.records)
  let paths = @src.vg_paths(graph)
  // 4 条变异2 读段(pos21 C): 参考 [16,26) 替换 [20]
  let reads = ["CGTACCGTAC", "CGTACCGTAC", "CGTACCGTAC", "CGTACCGTAC"]
  let gt = @src.vg_genotype(paths, reads)
  let g = gt.unwrap()
  inspect(g.gt, content="1/1")
  assert_eq(g.counts[0], 4)
  assert_eq(g.counts[1], 0)
}

///|
test "vg: 双杂合(3/3 纯合双 alt)" {
  let vcf = @src.parse_vcf(vg_test_vcf())
  let graph = @src.vg_build_graph("chr1", vg_test_reference(), vcf.records)
  let paths = @src.vg_paths(graph)
  // 双变异读段: [0,10) 变异1 + [16,26) 变异1
  let reads = ["ACGTGCGTAC", "CGTACCGTAC", "ACGTGCGTAC", "CGTACCGTAC"]
  let gt = @src.vg_genotype(paths, reads)
  let g = gt.unwrap()
  // 路径 2(变异1) 与路径 1(变异2) 各 2 -> 前二 tie 取索引小: 1 与 2
  inspect(g.gt, content="1/2")
}

///|
test "vg: indel 变异图与长读段分型" {
  let vcf = @src.parse_vcf(
    "##fileformat=VCFv4.2\n" +
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n" +
    "chr1\t5\t.\tA\tAGG\t.\tPASS\t.\n",
  )
  let graph = @src.vg_build_graph("chr1", vg_test_reference(), vcf.records)
  let paths = @src.vg_paths(graph)
  inspect(paths.length(), content="2")
  // ref 路径 40bp, 插入路径 40 - 1 + 3 = 42bp
  inspect(paths[0].seq.length(), content="40")
  inspect(paths[1].seq.length(), content="42")
  // 插入读段: "ACGT" + "AGG" + "CGTAC"
  let reads = ["ACGTAGGCGTAC", "ACGTAGGCGTAC", "ACGTAGGCGTAC"]
  let gt = @src.vg_genotype(paths, reads)
  let g = gt.unwrap()
  inspect(g.gt, content="1/1")
}

///|
test "vg: 重叠变异过滤" {
  let vcf = @src.parse_vcf(
    "##fileformat=VCFv4.2\n" +
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n" +
    "chr1\t5\t.\tAC\tGGG\t.\tPASS\t.\n" +
    "chr1\t6\t.\tC\tT\t.\tPASS\t.\n",
  )
  let graph = @src.vg_build_graph("chr1", vg_test_reference(), vcf.records)
  // pos5 ref "AC" 覆盖 [5,6], pos6 落在其内 -> 被过滤, 仅 1 个变异入图
  let paths = @src.vg_paths(graph)
  inspect(paths.length(), content="2")
}

///|
test "vg: GFA 与 JSON 导出" {
  let vcf = @src.parse_vcf(vg_test_vcf())
  let graph = @src.vg_build_graph("chr1", vg_test_reference(), vcf.records)
  let paths = @src.vg_paths(graph)
  let reads = ["ACGTACGTAC", "ACGTGCGTAC"]
  let gt = @src.vg_genotype(paths, reads)
  let gfa = @src.vg_to_gfa(graph)
  assert_eq(gfa.contains("H\tVN:Z:1.0"), true)
  assert_eq(gfa.contains("S\t0\tACGT"), true)
  assert_eq(gfa.contains("L\t0\t+\t1\t+\t*"), true)
  let json = @src.json_to_string(@src.vg_value(graph, paths, gt))
  assert_eq(json.contains("\"chrom\":\"chr1\""), true)
  assert_eq(json.contains("\"kind\":\"alt\""), true)
  assert_eq(json.contains("\"gt\":\"0/2\""), true)
  assert_eq(json.contains("\"paths\""), true)
}

///|
test "vg: 无变异输入" {
  let graph = @src.vg_build_graph("chr1", vg_test_reference(), [])
  inspect(graph.nodes.length(), content="1")
  inspect(graph.nodes[0].seq.length(), content="40")
  let paths = @src.vg_paths(graph)
  inspect(paths.length(), content="1")
  assert_eq(paths[0].seq, vg_test_reference())
  let gt = @src.vg_genotype(paths, ["ACGTACGTAC"])
  inspect(gt.unwrap().gt, content="0/0")
  // 空图
  let empty = @src.vg_build_graph("chr1", "", [])
  inspect(empty.nodes.length(), content="0")
  assert_true(@src.vg_genotype(@src.vg_paths(empty), []) is None)
}
