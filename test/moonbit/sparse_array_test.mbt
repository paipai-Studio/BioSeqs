///|
fn sparse_array_matrix() -> @src.SparseArray {
  @src.SparseArray::from_matrix([[1.0, 0.0, 2.0], [0.0, 3.0, 0.0]]) catch {
    _ => abort("failed to create sparse matrix fixture")
  }
}

///|
test "sparse_array: canonical COO merges duplicates and removes zeros" {
  let array = @src.SparseArray::from_coo(
    [2, 3],
    [[1, 2], [0, 1], [1, 2], [0, 0]],
    [2.0, 3.0, -2.0, 0.0],
  ) catch {
    _ => abort("valid COO input should be accepted")
  }
  assert_eq(array.nnzero(), 1)
  assert_eq(array.nzcoordinates(), [[0, 1]])
  assert_eq(array.nzvalues(), [3.0])
}

///|
test "sparse_array: dimensions length density and sparse marker" {
  let array = @src.sparse_array_sample() catch {
    _ => abort("sample should be valid")
  }
  assert_eq(array.ndim(), 3)
  assert_eq(array.dim(), [4, 5, 2])
  assert_eq(array.length(), 40)
  assert_eq(array.nnzero(), 5)
  assert_eq(array.density(), 0.125)
  assert_true(array.is_sparse())
}

///|
test "sparse_array: zero extent array is valid and empty" {
  let array = @src.SparseArray::zeros([2, 0, 4]) catch {
    _ => abort("zero extent should be valid")
  }
  assert_eq(array.length(), 0)
  assert_eq(array.nnzero(), 0)
  assert_eq(array.density(), 0.0)
  assert_true(array.mean() is None)
  assert_true(array.minimum() is None)
  assert_true(array.maximum() is None)
}

///|
test "sparse_array: rejects invalid dimensions" {
  let empty_rank = try {
    ignore(@src.SparseArray::zeros([]))
    false
  } catch {
    SparseArrayError(_) => true
  }
  let negative_extent = try {
    ignore(@src.SparseArray::zeros([2, -1]))
    false
  } catch {
    SparseArrayError(_) => true
  }
  assert_true(empty_rank)
  assert_true(negative_extent)
}

///|
test "sparse_array: rejects malformed COO input" {
  let length_mismatch = try {
    ignore(@src.SparseArray::from_coo([2, 2], [[0, 0]], []))
    false
  } catch {
    SparseArrayError(_) => true
  }
  let rank_mismatch = try {
    ignore(@src.SparseArray::from_coo([2, 2], [[0]], [1.0]))
    false
  } catch {
    SparseArrayError(_) => true
  }
  let out_of_bounds = try {
    ignore(@src.SparseArray::from_coo([2, 2], [[0, 2]], [1.0]))
    false
  } catch {
    SparseArrayError(_) => true
  }
  assert_true(length_mismatch)
  assert_true(rank_mismatch)
  assert_true(out_of_bounds)
}

///|
test "sparse_array: flat data uses R column-major ordering" {
  let array = @src.SparseArray::from_flat([2, 3], [1.0, 2.0, 0.0, 3.0, 4.0, 0.0]) catch {
    _ => abort("flat input should be valid")
  }
  assert_eq(array.to_dense_matrix(), [[1.0, 0.0, 4.0], [2.0, 3.0, 0.0]])
  assert_eq(array.to_flat(), [1.0, 2.0, 0.0, 3.0, 4.0, 0.0])
}

///|
test "sparse_array: rejects flat length mismatch" {
  let raised = try {
    ignore(@src.SparseArray::from_flat([2, 3], [1.0, 2.0]))
    false
  } catch {
    SparseArrayError(_) => true
  }
  assert_true(raised)
}

///|
test "sparse_array: matrix construction and dense conversion" {
  let array = sparse_array_matrix()
  assert_eq(array.dim(), [2, 3])
  assert_eq(array.nnzero(), 3)
  assert_eq(array.to_dense_matrix(), [[1.0, 0.0, 2.0], [0.0, 3.0, 0.0]])
}

///|
test "sparse_array: rejects ragged matrix and multidimensional dense conversion" {
  let ragged = try {
    ignore(@src.SparseArray::from_matrix([[1.0], [2.0, 3.0]]))
    false
  } catch {
    SparseArrayError(_) => true
  }
  let sample = @src.sparse_array_sample() catch {
    _ => abort("sample should be valid")
  }
  let not_matrix = try {
    ignore(sample.to_dense_matrix())
    false
  } catch {
    SparseArrayError(_) => true
  }
  assert_true(ragged)
  assert_true(not_matrix)
}

///|
test "sparse_array: coordinate and linear access" {
  let array = sparse_array_matrix()
  assert_eq(array.get([0, 0]), 1.0)
  assert_eq(array.get([0, 1]), 0.0)
  assert_eq(array.get([1, 1]), 3.0)
  assert_eq(array.get_linear(0), 1.0)
  assert_eq(array.get_linear(2), 0.0)
  assert_eq(array.get_linear(3), 3.0)
  assert_eq(array.get_linear(4), 2.0)
}

///|
test "sparse_array: rejects invalid access" {
  let array = sparse_array_matrix()
  let bad_coordinate = try {
    ignore(array.get([2, 0]))
    false
  } catch {
    SparseArrayError(_) => true
  }
  let bad_rank = try {
    ignore(array.get([0]))
    false
  } catch {
    SparseArrayError(_) => true
  }
  let bad_linear = try {
    ignore(array.get_linear(6))
    false
  } catch {
    SparseArrayError(_) => true
  }
  assert_true(bad_coordinate)
  assert_true(bad_rank)
  assert_true(bad_linear)
}

///|
test "sparse_array: accessors return defensive coordinate copies" {
  let array = sparse_array_matrix()
  let dimensions = array.dim()
  dimensions[0] = 99
  let coordinates = array.nzcoordinates()
  coordinates[0][0] = 1
  let entries = array.nonzero_entries()
  let first_coordinates = entries[0].coordinates()
  first_coordinates[0] = 1
  assert_eq(array.dim(), [2, 3])
  assert_eq(array.get([0, 0]), 1.0)
}

///|
test "sparse_array: immutable single-value assignment inserts replaces and removes" {
  let original = sparse_array_matrix()
  let inserted = original.with_value([1, 2], 4.0)
  let replaced = inserted.with_value([0, 2], 5.0)
  let removed = replaced.with_value([1, 1], 0.0)
  assert_eq(original.get([1, 2]), 0.0)
  assert_eq(inserted.get([1, 2]), 4.0)
  assert_eq(replaced.get([0, 2]), 5.0)
  assert_eq(removed.get([1, 1]), 0.0)
  assert_eq(removed.nnzero(), 3)
}

///|
test "sparse_array: ordered multiple assignment uses last replacement" {
  let array = sparse_array_matrix().with_values([[0, 0], [1, 2], [0, 0]], [
    8.0, 4.0, 9.0,
  ])
  assert_eq(array.get([0, 0]), 9.0)
  assert_eq(array.get([1, 2]), 4.0)
}

///|
test "sparse_array: multiple assignment validates parallel lengths" {
  let raised = try {
    ignore(sparse_array_matrix().with_values([[0, 0]], []))
    false
  } catch {
    SparseArrayError(_) => true
  }
  assert_true(raised)
}

///|
test "sparse_array: subset preserves order and repeated indices" {
  let subset = sparse_array_matrix().subset([[1, 0, 1], [2, 1]])
  assert_eq(subset.dim(), [3, 2])
  assert_eq(subset.to_dense_matrix(), [[0.0, 3.0], [2.0, 0.0], [0.0, 3.0]])
  assert_eq(subset.nnzero(), 3)
}

///|
test "sparse_array: empty subset remains sparse" {
  let subset = sparse_array_matrix().subset([[], [0, 1]])
  assert_eq(subset.dim(), [0, 2])
  assert_eq(subset.length(), 0)
  assert_eq(subset.nnzero(), 0)
}

///|
test "sparse_array: subset validates rank and bounds" {
  let array = sparse_array_matrix()
  let bad_rank = try {
    ignore(array.subset([[0, 1]]))
    false
  } catch {
    SparseArrayError(_) => true
  }
  let bad_index = try {
    ignore(array.subset([[0, 2], [0]]))
    false
  } catch {
    SparseArrayError(_) => true
  }
  assert_true(bad_rank)
  assert_true(bad_index)
}

///|
test "sparse_array: end-exclusive multidimensional slice" {
  let sliced = sparse_array_matrix().slice([0, 1], [2, 3])
  assert_eq(sliced.dim(), [2, 2])
  assert_eq(sliced.to_dense_matrix(), [[0.0, 2.0], [3.0, 0.0]])
}

///|
test "sparse_array: slice validates bounds" {
  let array = sparse_array_matrix()
  let bad_rank = try {
    ignore(array.slice([0], [1]))
    false
  } catch {
    SparseArrayError(_) => true
  }
  let reversed = try {
    ignore(array.slice([1, 0], [0, 2]))
    false
  } catch {
    SparseArrayError(_) => true
  }
  let too_large = try {
    ignore(array.slice([0, 0], [3, 2]))
    false
  } catch {
    SparseArrayError(_) => true
  }
  assert_true(bad_rank)
  assert_true(reversed)
  assert_true(too_large)
}

///|
test "sparse_array: aperm reorders multidimensional coordinates" {
  let sample = @src.sparse_array_sample() catch {
    _ => abort("sample should be valid")
  }
  let permuted = sample.aperm([2, 0, 1])
  assert_eq(permuted.dim(), [2, 4, 5])
  assert_eq(permuted.get([1, 0, 4]), 40.0)
  assert_eq(permuted.get([1, 2, 1]), 50.0)
}

///|
test "sparse_array: aperm validates complete permutation" {
  let sample = @src.sparse_array_sample() catch {
    _ => abort("sample should be valid")
  }
  let wrong_length = try {
    ignore(sample.aperm([1, 0]))
    false
  } catch {
    SparseArrayError(_) => true
  }
  let duplicate = try {
    ignore(sample.aperm([0, 0, 2]))
    false
  } catch {
    SparseArrayError(_) => true
  }
  assert_true(wrong_length)
  assert_true(duplicate)
}

///|
test "sparse_array: matrix transpose remains sparse" {
  let transposed = sparse_array_matrix().transpose()
  assert_eq(transposed.dim(), [3, 2])
  assert_eq(transposed.to_dense_matrix(), [[1.0, 0.0], [0.0, 3.0], [2.0, 0.0]])
  assert_eq(transposed.nnzero(), 3)
}

///|
test "sparse_array: transpose rejects non-matrix" {
  let sample = @src.sparse_array_sample() catch {
    _ => abort("sample should be valid")
  }
  let raised = try {
    ignore(sample.transpose())
    false
  } catch {
    SparseArrayError(_) => true
  }
  assert_true(raised)
}

///|
test "sparse_array: bind offsets coordinates along selected dimension" {
  let left = @src.SparseArray::from_matrix([[1.0, 0.0], [0.0, 2.0]])
  let right = @src.SparseArray::from_matrix([[3.0, 4.0]])
  let bound = left.bind(right, 0)
  assert_eq(bound.dim(), [3, 2])
  assert_eq(bound.to_dense_matrix(), [[1.0, 0.0], [0.0, 2.0], [3.0, 4.0]])
}

///|
test "sparse_array: bind validates rank margin and dimensions" {
  let matrix = sparse_array_matrix()
  let sample = @src.sparse_array_sample() catch {
    _ => abort("sample should be valid")
  }
  let rank_mismatch = try {
    ignore(matrix.bind(sample, 0))
    false
  } catch {
    SparseArrayError(_) => true
  }
  let bad_margin = try {
    ignore(matrix.bind(matrix, 2))
    false
  } catch {
    SparseArrayError(_) => true
  }
  let other = @src.SparseArray::from_matrix([[1.0, 2.0]])
  let dimension_mismatch = try {
    ignore(matrix.bind(other, 1))
    false
  } catch {
    SparseArrayError(_) => true
  }
  assert_true(rank_mismatch)
  assert_true(bad_margin)
  assert_true(dimension_mismatch)
}

///|
test "sparse_array: addition merges coordinates and drops cancellation" {
  let left = @src.SparseArray::from_matrix([[1.0, 0.0], [2.0, -3.0]])
  let right = @src.SparseArray::from_matrix([[-1.0, 4.0], [0.0, 3.0]])
  let result = left.add(right)
  assert_eq(result.to_dense_matrix(), [[0.0, 4.0], [2.0, 0.0]])
  assert_eq(result.nnzero(), 2)
}

///|
test "sparse_array: subtraction preserves sparse union" {
  let left = @src.SparseArray::from_matrix([[1.0, 0.0], [2.0, -3.0]])
  let right = @src.SparseArray::from_matrix([[-1.0, 4.0], [0.0, 3.0]])
  assert_eq(left.subtract(right).to_dense_matrix(), [[2.0, -4.0], [2.0, -6.0]])
}

///|
test "sparse_array: Hadamard product only visits coordinate intersection" {
  let left = @src.SparseArray::from_matrix([[1.0, 0.0], [2.0, -3.0]])
  let right = @src.SparseArray::from_matrix([[-1.0, 4.0], [0.0, 3.0]])
  let result = left.hadamard(right)
  assert_eq(result.to_dense_matrix(), [[-1.0, 0.0], [0.0, -9.0]])
  assert_eq(result.nnzero(), 2)
}

///|
test "sparse_array: arithmetic validates dimensions" {
  let left = sparse_array_matrix()
  let right = @src.SparseArray::from_matrix([[1.0, 2.0]])
  let add_raised = try {
    ignore(left.add(right))
    false
  } catch {
    SparseArrayError(_) => true
  }
  let product_raised = try {
    ignore(left.hadamard(right))
    false
  } catch {
    SparseArrayError(_) => true
  }
  assert_true(add_raised)
  assert_true(product_raised)
}

///|
test "sparse_array: scale and zero-preserving map remove generated zeros" {
  let array = sparse_array_matrix()
  assert_eq(array.scale(0.0).nnzero(), 0)
  let mapped = array.map_nonzero(fn(value : Double) -> Double {
    if value == 2.0 {
      0.0
    } else {
      value * value
    }
  })
  assert_eq(mapped.to_dense_matrix(), [[1.0, 0.0, 0.0], [0.0, 9.0, 0.0]])
  assert_eq(mapped.nnzero(), 2)
}

///|
test "sparse_array: whole-array summaries include implicit zeros" {
  let array = @src.SparseArray::from_matrix([[-2.0, 0.0, 5.0], [0.0, 0.0, 1.0]])
  assert_eq(array.sum(), 4.0)
  match array.mean() {
    Some(value) => assert_true((value - 4.0 / 6.0).abs() < 1.0e-12)
    None => assert_true(false)
  }
  match array.minimum() {
    Some(value) => assert_eq(value, -2.0)
    None => assert_true(false)
  }
  match array.maximum() {
    Some(value) => assert_eq(value, 5.0)
    None => assert_true(false)
  }
}

///|
test "sparse_array: extrema do not inject zero into fully dense arrays" {
  let array = @src.SparseArray::from_matrix([[-2.0, -1.0]])
  match array.minimum() {
    Some(value) => assert_eq(value, -2.0)
    None => assert_true(false)
  }
  match array.maximum() {
    Some(value) => assert_eq(value, -1.0)
    None => assert_true(false)
  }
}

///|
test "sparse_array: row and column sums and means" {
  let array = @src.SparseArray::from_matrix([[-2.0, 0.0, 5.0], [0.0, 0.0, 1.0]])
  assert_eq(array.row_sums(), [3.0, 1.0])
  assert_eq(array.column_sums(), [-2.0, 0.0, 6.0])
  let row_means = array.row_means()
  assert_eq(row_means[0], 1.0)
  assert_true((row_means[1] - 1.0 / 3.0).abs() < 1.0e-12)
  assert_eq(array.column_means(), [-1.0, 0.0, 3.0])
}

///|
test "sparse_array: row and column nonzero counts" {
  let array = sparse_array_matrix()
  assert_eq(array.row_nonzero_counts(), [2, 1])
  assert_eq(array.column_nonzero_counts(), [1, 1, 1])
}

///|
test "sparse_array: matrix summaries reject rank and empty means" {
  let sample = @src.sparse_array_sample() catch {
    _ => abort("sample should be valid")
  }
  let rank_raised = try {
    ignore(sample.row_sums())
    false
  } catch {
    SparseArrayError(_) => true
  }
  let no_columns = @src.SparseArray::zeros([2, 0]) catch {
    _ => abort("zero-column matrix should be valid")
  }
  let row_mean_raised = try {
    ignore(no_columns.row_means())
    false
  } catch {
    SparseArrayError(_) => true
  }
  let no_rows = @src.SparseArray::zeros([0, 2]) catch {
    _ => abort("zero-row matrix should be valid")
  }
  let column_mean_raised = try {
    ignore(no_rows.column_means())
    false
  } catch {
    SparseArrayError(_) => true
  }
  assert_true(rank_raised)
  assert_true(row_mean_raised)
  assert_true(column_mean_raised)
}

///|
test "sparse_array: sparse matrix multiplication returns dense result" {
  let left = sparse_array_matrix()
  let right = @src.SparseArray::from_matrix([[0.0, 4.0], [5.0, 0.0], [6.0, 7.0]])
  assert_eq(left.matmul(right), [[12.0, 18.0], [15.0, 0.0]])
}

///|
test "sparse_array: crossprod and tcrossprod" {
  let array = sparse_array_matrix()
  assert_eq(array.crossprod(), [
    [1.0, 0.0, 2.0],
    [0.0, 9.0, 0.0],
    [2.0, 0.0, 4.0],
  ])
  assert_eq(array.tcrossprod(), [[5.0, 0.0], [0.0, 9.0]])
}

///|
test "sparse_array: matrix multiplication validates rank and dimensions" {
  let matrix = sparse_array_matrix()
  let incompatible = @src.SparseArray::from_matrix([[1.0, 2.0]])
  let dimension_raised = try {
    ignore(matrix.matmul(incompatible))
    false
  } catch {
    SparseArrayError(_) => true
  }
  let sample = @src.sparse_array_sample() catch {
    _ => abort("sample should be valid")
  }
  let rank_raised = try {
    ignore(sample.matmul(matrix))
    false
  } catch {
    SparseArrayError(_) => true
  }
  assert_true(dimension_raised)
  assert_true(rank_raised)
}

///|
test "sparse_array: three-dimensional flat round trip" {
  let original = @src.sparse_array_sample() catch {
    _ => abort("sample should be valid")
  }
  let restored = @src.SparseArray::from_flat(original.dim(), original.to_flat()) catch {
    _ => abort("round trip should be valid")
  }
  assert_eq(restored.dim(), original.dim())
  assert_eq(restored.nzcoordinates(), original.nzcoordinates())
  assert_eq(restored.nzvalues(), original.nzvalues())
}

///|
test "sparse_array: sample and summary" {
  let sample = @src.sparse_array_sample() catch {
    _ => abort("sample should be valid")
  }
  let summary = sample.summary()
  assert_true(summary.contains("SparseArray(4 x 5 x 2"))
  assert_true(summary.contains("nnzero=5"))
  assert_true(summary.contains("density=0.125"))
}
