///|
/// Tests for scanpy tools: rank_genes_groups (t-test differential expression).

///|
test "scanpy_rank_genes_groups" {
  // 6 cells × 2 genes; 3 cells in group A, 3 in group B
  let x = [
    [1.0, 2.0],
    [2.0, 2.0],
    [3.0, 2.0],
    [4.0, 2.0],
    [5.0, 2.0],
    [6.0, 2.0],
  ]
  let groups = ["A", "A", "A", "B", "B", "B"]
  let res = @src.scanpy_rank_genes_groups(x, groups, "A", "B")
  assert_eq(res.length(), 2)
  // gene0: A mean=2 var=1, B mean=5 var=1 -> Welch t = 3/sqrt(2/3) > 3, p < 0.05
  assert_eq(res[0].mean_ref, 2.0)
  assert_eq(res[0].mean_target, 5.0)
  assert_eq(res[0].t_stat > 3.0, true)
  // logfoldchange = log2((5+1)/(2+1)) = 1.0
  assert_eq(sp_close(res[0].logfoldchange, 1.0), true)
  // gene1: constant across cells -> zero variance -> t=0, p=1
  assert_eq(res[1].t_stat, 0.0)
  assert_eq(res[1].p_value, 1.0)
}

// ============================================================
// tsne / dendrogram
// ============================================================

test "scanpy_tsne" {
  let x = [[1.0, 0.0], [0.0, 1.0], [1.0, 1.0], [0.0, 0.0]]
  let emb = @src.scanpy_tsne(x, perplexity=2.0, dims=2)
  assert_eq(emb.length(), 4)
  assert_eq(emb[0].length(), 2)
}

test "scanpy_dendrogram" {
  // 2 cells × 3 genes → 3 个基因的聚类树，merge 长度 = 2
  let x = [[1.0, 2.0, 3.0], [3.0, 2.0, 1.0]]
  let tree = @src.scanpy_dendrogram(x)
  assert_eq(tree.n, 3)
  assert_eq(tree.merge.length(), 2)
}
