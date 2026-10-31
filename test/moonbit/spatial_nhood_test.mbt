///|
/// 构造 n x n 网格坐标.
fn nhood_test_grid(n : Int32) -> Array[(Double, Double)] {
  let out : Array[(Double, Double)] = Array(capacity=0)
  let mut i = 0
  while i < n {
    let mut j = 0
    while j < n {
      out.push((i.to_double(), j.to_double()))
      j = j + 1
    }
    i = i + 1
  }
  out
}

///|
/// 左半 A / 右半 B 的类别标签.
fn nhood_test_split_labels(
  coords : Array[(Double, Double)],
  cut : Double,
) -> Array[String] {
  let out : Array[String] = Array(capacity=0)
  for c in coords {
    let (x, _) = c
    if x < cut {
      out.push("A")
    } else {
      out.push("B")
    }
  }
  out
}

///|
test "spatial knn edges nearest first" {
  let coords : Array[(Double, Double)] = [(0.0, 0.0), (1.0, 0.0), (5.0, 5.0)]
  let edges = @src.spatial_knn_edges(coords, 1)
  assert_eq(edges.length(), 3)
  assert_eq(edges[0][0], 1)
  assert_eq(edges[1][0], 0)
  assert_eq(edges[2][0], 1)
}

///|
test "spatial knn edges caps at n minus one" {
  let coords : Array[(Double, Double)] = [(0.0, 0.0), (1.0, 0.0)]
  let edges = @src.spatial_knn_edges(coords, 5)
  assert_eq(edges.length(), 2)
  assert_eq(edges[0].length(), 1)
}

///|
test "nhood label length mismatch raises" {
  let coords : Array[(Double, Double)] = [(0.0, 0.0)]
  let labels : Array[String] = []
  let msg = try {
    ignore(@src.nhood_enrichment(coords, labels))
    "NO_ERROR"
  } catch {
    @src.NhoodError(m) => m
  }
  assert_true(msg.contains("不一致"))
}

///|
test "nhood observed matrix counts neighbors" {
  let coords : Array[(Double, Double)] = [
    (0.0, 0.0),
    (1.0, 0.0),
    (0.0, 1.0),
    (1.0, 1.0),
  ]
  let labels : Array[String] = ["A", "B", "B", "A"]
  let r = @src.nhood_enrichment(coords, labels, k=1, n_perms=0)
  assert_eq(r.labels().length(), 2)
  // 每个点的 1 个邻居必为异类, 观察到 2 次 A->B 与 2 次 B->A
  assert_eq(r.observed()[0][1], 2.0)
  assert_eq(r.observed()[1][0], 2.0)
  assert_eq(r.observed()[0][0], 0.0)
}

///|
test "nhood zscore zero without permutations" {
  let coords = nhood_test_grid(4)
  let labels = nhood_test_split_labels(coords, 2.0)
  let r = @src.nhood_enrichment(coords, labels, k=3, n_perms=0)
  assert_eq(r.zscore_of("A", "B"), 0.0)
  assert_eq(r.n_perms(), 0)
}

///|
test "nhood clustered labels enriched within type" {
  let coords = nhood_test_grid(6)
  let labels = nhood_test_split_labels(coords, 3.0)
  let r = @src.nhood_enrichment(coords, labels, k=6, n_perms=50, seed=1)
  // 空间聚集: 同类相邻远多于随机
  assert_true(r.zscore_of("A", "A") > 1.0)
  assert_true(r.zscore_of("B", "B") > 1.0)
  // 异类相邻少于随机
  assert_true(r.zscore_of("A", "B") < 0.0)
}

///|
test "nhood label index and enrichment helper" {
  let coords = nhood_test_grid(6)
  let labels = nhood_test_split_labels(coords, 3.0)
  let r = @src.nhood_enrichment(coords, labels, k=6, n_perms=30, seed=2)
  assert_eq(r.label_index("A"), 0)
  assert_eq(r.label_index("B"), 1)
  assert_eq(r.label_index("C"), -1)
  assert_true(@src.nhood_is_enriched(r, "A", "A", threshold=0.5))
  assert_true(r.zscore_of("A", "B") < 0.0)
}

///|
test "nhood k and accessors" {
  let coords = nhood_test_grid(3)
  let labels = nhood_test_split_labels(coords, 1.5)
  let k : Int32 = 2
  let r = @src.nhood_enrichment(coords, labels, k~, n_perms=10, seed=3)
  assert_eq(r.k(), 2)
  assert_eq(r.n_perms(), 10)
  assert_eq(r.perm_mean().length(), 2)
  assert_eq(r.perm_std().length(), 2)
}
