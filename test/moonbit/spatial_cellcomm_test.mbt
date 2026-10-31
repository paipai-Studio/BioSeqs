///|
/// 构造 sender/receiver 相邻的四细胞布局(2x2 网格).
fn comm_test_coords() -> Array[(Double, Double)] {
  [(0.0, 0.0), (1.0, 0.0), (0.0, 1.0), (1.0, 1.0)]
}

///|
/// 空间远离布局.
fn comm_test_far_coords() -> Array[(Double, Double)] {
  [(0.0, 0.0), (1.0, 0.0), (100.0, 0.0), (101.0, 0.0)]
}

///|
test "spatial cellchat detects neighbor signal" {
  let coords = comm_test_coords()
  let labels : Array[String] = ["S", "S", "R", "R"]
  let expr : Array[Array[Double]] = [
    [5.0, 0.0],
    [5.0, 0.0],
    [0.0, 3.0],
    [0.0, 3.0],
  ]
  let genes : Array[String] = ["LIG", "REC"]
  let pairs : Array[@src.LRPair] = [@src.lr_pair("LIG", "REC")]
  let res = @src.spatial_cellchat(coords, labels, expr, genes, pairs, k=2)
  assert_eq(res.length(), 1)
  assert_eq(res[0].source(), "S")
  assert_eq(res[0].target(), "R")
  assert_eq(res[0].ligand(), "LIG")
  assert_eq(res[0].receptor(), "REC")
  assert_eq(res[0].n_pairs(), 2)
  assert_true(res[0].score() > 0.0)
}

///|
test "spatial cellchat requires spatial proximity" {
  let coords = comm_test_far_coords()
  let labels : Array[String] = ["S", "S", "R", "R"]
  let expr : Array[Array[Double]] = [
    [5.0, 0.0],
    [5.0, 0.0],
    [0.0, 3.0],
    [0.0, 3.0],
  ]
  let genes : Array[String] = ["LIG", "REC"]
  let pairs : Array[@src.LRPair] = [@src.lr_pair("LIG", "REC")]
  let res = @src.spatial_cellchat(coords, labels, expr, genes, pairs, k=1)
  assert_eq(res.length(), 0)
}

///|
test "spatial cellchat skips missing genes" {
  let coords = comm_test_coords()
  let labels : Array[String] = ["S", "S", "R", "R"]
  let expr : Array[Array[Double]] = [
    [5.0, 0.0],
    [5.0, 0.0],
    [0.0, 3.0],
    [0.0, 3.0],
  ]
  let genes : Array[String] = ["LIG", "REC"]
  let pairs : Array[@src.LRPair] = [@src.lr_pair("LIG", "NOT_A_GENE")]
  let res = @src.spatial_cellchat(coords, labels, expr, genes, pairs, k=2)
  assert_eq(res.length(), 0)
}

///|
test "spatial cellchat shape mismatch raises" {
  let coords = comm_test_coords()
  let labels : Array[String] = ["S", "S"]
  let expr : Array[Array[Double]] = [[5.0, 0.0]]
  let genes : Array[String] = ["LIG", "REC"]
  let pairs : Array[@src.LRPair] = [@src.lr_pair("LIG", "REC")]
  let msg = try {
    ignore(@src.spatial_cellchat(coords, labels, expr, genes, pairs))
    "NO_ERROR"
  } catch {
    @src.SpatialCommError(m) => m
  }
  assert_true(msg.contains("不一致"))
}

///|
test "spatial cellchat autocrine same type" {
  let coords : Array[(Double, Double)] = [(0.0, 0.0), (1.0, 0.0)]
  let labels : Array[String] = ["S", "S"]
  let expr : Array[Array[Double]] = [[2.0, 2.0], [2.0, 2.0]]
  let genes : Array[String] = ["LIG", "REC"]
  let pairs : Array[@src.LRPair] = [@src.lr_pair("LIG", "REC")]
  let res = @src.spatial_cellchat(coords, labels, expr, genes, pairs, k=1)
  assert_eq(res.length(), 1)
  assert_eq(res[0].source(), "S")
  assert_eq(res[0].target(), "S")
  assert_eq(res[0].n_pairs(), 2)
}

///|
test "spatial cellchat min pairs filter" {
  let coords = comm_test_coords()
  let labels : Array[String] = ["S", "S", "R", "R"]
  let expr : Array[Array[Double]] = [
    [5.0, 0.0],
    [5.0, 0.0],
    [0.0, 3.0],
    [0.0, 3.0],
  ]
  let genes : Array[String] = ["LIG", "REC"]
  let pairs : Array[@src.LRPair] = [@src.lr_pair("LIG", "REC")]
  let res = @src.spatial_cellchat(
    coords,
    labels,
    expr,
    genes,
    pairs,
    k=2,
    min_pairs=10,
  )
  assert_eq(res.length(), 0)
}

///|
test "spatial cellchat sorted by score" {
  let coords = comm_test_coords()
  let labels : Array[String] = ["S", "S", "R", "R"]
  let expr : Array[Array[Double]] = [
    [5.0, 0.0, 1.0],
    [5.0, 0.0, 1.0],
    [0.0, 3.0, 0.1],
    [0.0, 3.0, 0.1],
  ]
  let genes : Array[String] = ["LIG", "REC", "WEAK"]
  let pairs : Array[@src.LRPair] = [
    @src.lr_pair("LIG", "REC"),
    @src.lr_pair("WEAK", "REC"),
  ]
  let res = @src.spatial_cellchat(coords, labels, expr, genes, pairs, k=2)
  assert_true(res.length() >= 2)
  assert_true(res[0].score() >= res[1].score())
}

///|
test "spatial cellchat aggregate" {
  let coords = comm_test_coords()
  let labels : Array[String] = ["S", "S", "R", "R"]
  let expr : Array[Array[Double]] = [
    [5.0, 0.0],
    [5.0, 0.0],
    [0.0, 3.0],
    [0.0, 3.0],
  ]
  let genes : Array[String] = ["LIG", "REC"]
  let pairs : Array[@src.LRPair] = [@src.lr_pair("LIG", "REC")]
  let res = @src.spatial_cellchat(coords, labels, expr, genes, pairs, k=2)
  let agg = @src.spatial_cellchat_aggregate(res)
  assert_eq(agg.length(), 1)
  assert_true(agg.get("S->R").unwrap() > 0.0)
}
