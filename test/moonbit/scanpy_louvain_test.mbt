///|
/// Tests for scanpy `tl.louvain`.

///|
test "scanpy_louvain_two_clusters" {
  // two disconnected cliques of three nodes
  let a = [
    [0.0, 1.0, 1.0, 0.0, 0.0, 0.0],
    [1.0, 0.0, 1.0, 0.0, 0.0, 0.0],
    [1.0, 1.0, 0.0, 0.0, 0.0, 0.0],
    [0.0, 0.0, 0.0, 0.0, 1.0, 1.0],
    [0.0, 0.0, 0.0, 1.0, 0.0, 1.0],
    [0.0, 0.0, 0.0, 1.0, 1.0, 0.0],
  ]
  let labels = @src.scanpy_louvain(a)
  assert_eq(labels.length(), 6)
  // nodes 0,1,2 form one community; 3,4,5 another
  assert_eq((labels[0], labels[1], labels[2]), (labels[0], labels[0], labels[0]))
  assert_eq((labels[3], labels[4], labels[5]), (labels[3], labels[3], labels[3]))
  assert_true(labels[0] != labels[3])
}

///|
test "scanpy_louvain_two_cliques_with_bridge" {
  // two cliques joined by a weak bridge edge (3-4)
  let a = [
    [0.0, 1.0, 1.0, 0.0, 0.0, 0.0],
    [1.0, 0.0, 1.0, 0.0, 0.0, 0.0],
    [1.0, 1.0, 0.0, 1.0, 0.0, 0.0],
    [0.0, 0.0, 1.0, 0.0, 1.0, 1.0],
    [0.0, 0.0, 0.0, 1.0, 0.0, 1.0],
    [0.0, 0.0, 0.0, 1.0, 1.0, 0.0],
  ]
  let labels = @src.scanpy_louvain(a)
  // still splits into two communities
  assert_eq((labels[0], labels[1]), (labels[0], labels[0]))
  assert_eq((labels[4], labels[5]), (labels[4], labels[4]))
  assert_true(labels[0] != labels[4])
}