///|
test "tree_construction_distance_calculator" {
  let calc = @src.DistanceCalculator::new()

  let seq1 = "ACGTACGTACGT"
  let seq2 = "ACGTACGTACGT"
  let dist_identical = calc.pairwise(seq1, seq2)
  assert_true(dist_identical >= 0.0 && dist_identical < 0.001)

  let seq3 = "ACGTACGTACGT"
  let seq4 = "ACGTACGTATGT"
  let dist_diff = calc.pairwise(seq3, seq4)
  assert_true(dist_diff > 0.0)
}

///|
test "tree_construction_upgma" {
  let tree_cons = @src.TreeConstructor::new(tree_method="upgma")

  let dist_matrix = [
    [0.0, 5.0, 9.0, 9.0, 8.0],
    [5.0, 0.0, 10.0, 10.0, 9.0],
    [9.0, 10.0, 0.0, 7.0, 3.0],
    [9.0, 10.0, 7.0, 0.0, 4.0],
    [8.0, 9.0, 3.0, 4.0, 0.0],
  ]
  let names = ["A", "B", "C", "D", "E"]

  let tree = tree_cons.build_tree(dist_matrix, names)
  assert_true(tree.clades.length() == 2)
}

///|
test "tree_construction_nj" {
  let tree_cons = @src.TreeConstructor::new(tree_method="nj")

  let dist_matrix = [
    [0.0, 5.0, 9.0, 9.0, 8.0],
    [5.0, 0.0, 10.0, 10.0, 9.0],
    [9.0, 10.0, 0.0, 7.0, 3.0],
    [9.0, 10.0, 7.0, 0.0, 4.0],
    [8.0, 9.0, 3.0, 4.0, 0.0],
  ]
  let names = ["A", "B", "C", "D", "E"]

  let tree = tree_cons.build_tree(dist_matrix, names)
  assert_true(tree.clades.length() >= 2)
}

///|
test "tree_construction_wpgma" {
  let tree_cons = @src.TreeConstructor::new(tree_method="wpgma")

  let dist_matrix = [
    [0.0, 2.0, 6.0, 6.0],
    [2.0, 0.0, 6.0, 6.0],
    [6.0, 6.0, 0.0, 2.0],
    [6.0, 6.0, 2.0, 0.0],
  ]
  let names = ["W", "X", "Y", "Z"]

  let tree = tree_cons.build_tree(dist_matrix, names)
  assert_true(tree.clades.length() == 2)
}

///|
test "tree_construction_jukes_cantor" {
  let calc = @src.DistanceCalculator::new(model="jukes_cantor")

  let seq1 = "ACGTACGTACGTACGTACGT"
  let seq2 = "ACGTACGTACGTACGTACGT"
  let dist = calc.pairwise(seq1, seq2)
  assert_true(dist >= 0.0 && dist < 0.01)
}

///|
test "tree_construction_kimura" {
  let calc = @src.DistanceCalculator::new(model="kimura")

  let seq1 = "AAAAAAAAAA"
  let seq2 = "AGAAAAGAAA"
  let dist = calc.pairwise(seq1, seq2)
  assert_true(dist >= 0.0)
}

///|
test "tree_construction_distance_matrix" {
  let calc = @src.DistanceCalculator::new()
  let seqs = [("seq1", "AAAA"), ("seq2", "AAAT"), ("seq3", "AATT")]
  let matrix = calc.get_distance(seqs)
  assert_eq(matrix.length(), 3)
  assert_eq(matrix[0].length(), 3)
  assert_true(matrix[0][0] == 0.0)
  assert_true(matrix[0][1] >= 0.0)
}
