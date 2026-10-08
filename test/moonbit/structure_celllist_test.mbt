test "celllist_basic" {
  // 3 atoms in a line: 0, 1, 2 at x=0, 5, 10
  let coords = [
    [0.0, 0.0, 0.0],
    [5.0, 0.0, 0.0],
    [10.0, 0.0, 0.0],
  ]
  let cl = @src.celllist_new(coords, 3.0)
  // Neighbors of atom 0 within radius 6: should find atom 1
  let n0 = @src.celllist_neighbors(cl, 0, 6.0)
  assert_true(n0.length() == 1)
  assert_true(n0[0] == 1)
  // Neighbors of atom 1 within radius 6: should find atoms 0 and 2
  let n1 = @src.celllist_neighbors(cl, 1, 6.0)
  assert_true(n1.length() == 2)
}

test "celllist_all_pairs" {
  let coords = [
    [0.0, 0.0, 0.0],
    [1.0, 0.0, 0.0],
    [5.0, 0.0, 0.0],
    [6.0, 0.0, 0.0],
  ]
  let cl = @src.celllist_new(coords, 3.0)
  let pairs = @src.celllist_all_pairs(cl, 2.0)
  // Should find (0,1) and (2,3)
  assert_true(pairs.length() == 2)
}

test "celllist_point_query" {
  let coords = [
    [0.0, 0.0, 0.0],
    [1.0, 1.0, 1.0],
    [10.0, 10.0, 10.0],
  ]
  let cl = @src.celllist_new(coords, 2.0)
  // Query point near atom 0
  let neighbors = @src.celllist_neighbors_point(cl, 0.5, 0.5, 0.5, 2.0, -1)
  // Should find atoms 0 and 1
  assert_true(neighbors.length() == 2)
}

test "celllist_empty" {
  let coords : Array[Array[Double]] = []
  let cl = @src.celllist_new(coords, 1.0)
  let pairs = @src.celllist_all_pairs(cl, 1.0)
  assert_true(pairs.length() == 0)
}
