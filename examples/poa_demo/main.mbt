///|
/// POA 图多序列比对演示: 多样本共识与 MSA.
fn main {
  println("=== 场景 1: 测序误差下的共识 ===")
  let g1 = @src.PoaGraph::new("ACGTACGTAC")
  ignore(@src.poa_add(g1, "ACGTACGTAC"))
  ignore(@src.poa_add(g1, "ACCTACGTAC"))
  ignore(@src.poa_add(g1, "ACGTACGTTC"))
  println("  节点数: " + @src.poa_node_count(g1).to_string())
  println("  共识:   " + @src.poa_consensus(g1))
  println("  MSA:")
  for row in @src.poa_msa(g1) {
    println("    " + row)
  }
  println("=== 场景 2: 含 indel 的簇 ===")
  let g2 = @src.PoaGraph::new("ATGGCCTAGC")
  ignore(@src.poa_add(g2, "ATGGCCTTAGC"))
  ignore(@src.poa_add(g2, "ATGGCCTAGC"))
  println("  节点数: " + @src.poa_node_count(g2).to_string())
  println("  共识:   " + @src.poa_consensus(g2))
  println("  MSA:")
  for row in @src.poa_msa(g2) {
    println("    " + row)
  }
  println("=== 场景 3: 自定义打分 ===")
  let g3 = @src.PoaGraph::new("AAAA", match_score=5, mismatch_score=-4, gap=-3)
  let s = @src.poa_add(g3, "AAAT", match_score=5, mismatch_score=-4, gap=-3)
  println("  得分: " + s.to_string() + " (3 match + 1 mismatch)")
}
