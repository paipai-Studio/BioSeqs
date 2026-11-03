///|
test "poa single sequence consensus identity" {
  let g = @src.PoaGraph::new("ACGTACGT")
  assert_eq(@src.poa_consensus(g), "ACGTACGT")
  assert_eq(@src.poa_node_count(g), 8)
  assert_eq(@src.poa_msa(g), ["ACGTACGT"])
}

///|
test "poa identical sequences merge weights" {
  let g = @src.PoaGraph::new("ACGT")
  let score = @src.poa_add(g, "ACGT")
  assert_eq(score, 8)
  assert_eq(@src.poa_node_count(g), 4)
  assert_eq(@src.poa_consensus(g), "ACGT")
  let msa = @src.poa_msa(g)
  assert_eq(msa.length(), 2)
  assert_eq(msa[0], "ACGT")
  assert_eq(msa[1], "ACGT")
}

///|
test "poa majority vote over mismatch" {
  let g = @src.PoaGraph::new("ACGT")
  ignore(@src.poa_add(g, "ACGT"))
  ignore(@src.poa_add(g, "ACCT"))
  // 3 票中 A=3 C=3 G=1(ACCT 无 G) T=3 -> 共识 ACGT
  assert_eq(@src.poa_consensus(g), "ACGT")
}

///|
test "poa insertion handled" {
  let g = @src.PoaGraph::new("ACGT")
  ignore(@src.poa_add(g, "ACGGT"))
  // 多数列: A C G(ins G 为少数列) T
  let cons = @src.poa_consensus(g)
  assert_true(cons == "ACGT" || cons == "ACGGT")
  assert_true(cons.length() >= 4)
}

///|
test "poa deletion handled" {
  let g = @src.PoaGraph::new("ACGTACGT")
  ignore(@src.poa_add(g, "ACGACGT"))
  let cons = @src.poa_consensus(g)
  // 8 列多数: 第二个 T 为 2:1 -> ACGTACGT
  assert_eq(cons, "ACGTACGT")
}

///|
test "poa empty sequence" {
  let g = @src.PoaGraph::new("A")
  let s = @src.poa_add(g, "")
  assert_eq(s, 0)
  assert_eq(@src.poa_consensus(g), "A")
}

///|
test "poa score positive for similar" {
  let g = @src.PoaGraph::new("ACGT")
  let s = @src.poa_add(g, "ACGT")
  assert_true(s > 0)
}

///|
test "poa distant sequences lower score" {
  let g = @src.PoaGraph::new("AAAA")
  let near = @src.poa_add(g, "AAAT")
  let g2 = @src.PoaGraph::new("AAAA")
  let far = @src.poa_add(g2, "TTTT")
  assert_true(near > far)
}

///|
test "poa msa rows same length" {
  let g = @src.PoaGraph::new("ACGTACGT")
  ignore(@src.poa_add(g, "ACGT"))
  ignore(@src.poa_add(g, "ACGTTACGT"))
  let msa = @src.poa_msa(g)
  assert_eq(msa.length(), 3)
  let n0 = msa[0].length()
  for row in msa {
    assert_eq(row.length(), n0)
  }
}

///|
test "poa consensus over cluster with repeats" {
  let g = @src.PoaGraph::new("ATATATA")
  ignore(@src.poa_add(g, "ATATATA"))
  ignore(@src.poa_add(g, "ATATATA"))
  assert_eq(@src.poa_consensus(g), "ATATATA")
  assert_eq(@src.poa_node_count(g), 7)
}

///|
test "poa custom scoring" {
  let g = @src.PoaGraph::new("AAAA", match_score=5, mismatch_score=-4, gap=-3)
  let s = @src.poa_add(g, "AAAA", match_score=5, mismatch_score=-4, gap=-3)
  assert_eq(s, 20)
}

///|
test "poa last score tracking" {
  let g = @src.PoaGraph::new("ACGT")
  ignore(@src.poa_add(g, "ACGT"))
  assert_eq(@src.poa_last_score(g), 8)
}
