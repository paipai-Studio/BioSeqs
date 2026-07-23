///|
/// Tests for Matrix module.

test "BiocMatrix dense_from_array" {
  let mat = @src.BiocMatrix::dense_from_array([[1.0, 2.0], [3.0, 4.0]])
  assert_eq(mat.nrow(), 2)
  assert_eq(mat.ncol(), 2)
}

test "BiocMatrix get element" {
  let mat = @src.BiocMatrix::dense_from_array([[1.0, 2.0], [3.0, 4.0]])
  assert_eq(mat.get(0, 0), 1.0)
  assert_eq(mat.get(1, 1), 4.0)
}

test "BiocMatrix set element" {
  let mat = @src.BiocMatrix::dense_from_array([[1.0, 2.0], [3.0, 4.0]])
  let mat2 = mat.set(0, 0, 10.0)
  assert_eq(mat2.get(0, 0), 10.0)
}

test "BiocMatrix csc_from_triplets" {
  let mat = @src.BiocMatrix::csc_from_triplets([0, 1, 0, 1], [0, 0, 1, 1], [1.0, 3.0, 2.0, 4.0], 2, 2)
  assert_eq(mat.nrow(), 2)
  assert_eq(mat.ncol(), 2)
  assert_eq(mat.get(0, 0), 1.0)
}

test "BiocMatrix csr_from_triplets" {
  let mat = @src.BiocMatrix::csr_from_triplets([0, 0, 1, 1], [0, 1, 0, 1], [1.0, 2.0, 3.0, 4.0], 2, 2)
  assert_eq(mat.nrow(), 2)
  assert_eq(mat.ncol(), 2)
  assert_eq(mat.get(0, 1), 2.0)
}

test "BiocMatrix transpose" {
  let mat = @src.BiocMatrix::dense_from_array([[1.0, 2.0], [3.0, 4.0]])
  let t = mat.transpose()
  assert_eq(t.get(0, 1), 3.0)
  assert_eq(t.get(1, 0), 2.0)
}

test "BiocMatrix add" {
  let mat1 = @src.BiocMatrix::dense_from_array([[1.0, 2.0], [3.0, 4.0]])
  let mat2 = @src.BiocMatrix::dense_from_array([[1.0, 0.0], [0.0, 1.0]])
  let result = mat1.add(mat2)
  assert_true(result is Some(_))
  assert_eq(result.unwrap().get(0, 0), 2.0)
}

test "BiocMatrix multiply" {
  let mat1 = @src.BiocMatrix::dense_from_array([[1.0, 2.0], [3.0, 4.0]])
  let mat2 = @src.BiocMatrix::dense_from_array([[2.0, 0.0], [0.0, 2.0]])
  let result = mat1.multiply(mat2)
  assert_true(result is Some(_))
  assert_eq(result.unwrap().get(0, 0), 2.0)
}

test "BiocMatrix row_sums" {
  let mat = @src.BiocMatrix::dense_from_array([[1.0, 2.0], [3.0, 4.0]])
  let sums = mat.row_sums()
  assert_eq(sums.length(), 2)
  assert_eq(sums[0], 3.0)
}

test "BiocMatrix col_sums" {
  let mat = @src.BiocMatrix::dense_from_array([[1.0, 2.0], [3.0, 4.0]])
  let sums = mat.col_sums()
  assert_eq(sums.length(), 2)
  assert_eq(sums[0], 4.0)
}

test "BiocMatrix row_means" {
  let mat = @src.BiocMatrix::dense_from_array([[1.0, 3.0], [2.0, 4.0]])
  let means = mat.row_means()
  assert_eq(means[0], 2.0)
}

test "BiocMatrix col_means" {
  let mat = @src.BiocMatrix::dense_from_array([[1.0, 2.0], [3.0, 4.0]])
  let means = mat.col_means()
  assert_eq(means[0], 2.0)
}

test "BiocMatrix norm" {
  let mat = @src.BiocMatrix::dense_from_array([[3.0, 4.0]])
  let norm = mat.norm(2.0)
  assert_eq(norm, 5.0)
}

test "BiocMatrix dim" {
  let mat = @src.BiocMatrix::dense_from_array([[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]])
  let (nrow, ncol) = mat.dim()
  assert_eq(nrow, 2)
  assert_eq(ncol, 3)
}

test "BiocMatrix nnz" {
  let mat = @src.BiocMatrix::csc_from_triplets([0, 1], [0, 1], [1.0, 2.0], 2, 2)
  assert_eq(mat.nnz(), 2)
}

test "create_example_dense_matrix" {
  let mat = @src.create_example_dense_matrix()
  assert_eq(mat.nrow(), 3)
  assert_eq(mat.ncol(), 3)
}

test "create_example_sparse_matrix" {
  let mat = @src.create_example_sparse_matrix()
  assert_eq(mat.nrow(), 3)
  assert_eq(mat.ncol(), 2)
}
