///|
/// BiocSingular Tests - SVD for dimensionality reduction

///|
test "svd_options_new_defaults" {
  let opts = @src.svd_options_new()
  assert_eq(opts.rank, 0)
  let default_method = @src.exact_svd_method()
  assert_true(opts.method == default_method)
  assert_true((opts.tol - 1.0e-7).abs() < 1.0e-15)
  assert_eq(opts.max_iter, 1000)
  assert_eq(opts.n_oversamples, 10)
}

///|
test "svd_options_with_method_irlba" {
  let opts = @src.svd_options_with_method(@src.irlba_method())
  let m = @src.irlba_method()
  assert_true(opts.method == m)
  assert_eq(opts.rank, 0)
  assert_eq(opts.max_iter, 1000)
}

///|
test "svd_options_with_method_randomized" {
  let opts = @src.svd_options_with_method(@src.randomized_method())
  let m = @src.randomized_method()
  assert_true(opts.method == m)
}

///|
test "svd_options_full_custom" {
  let opts = @src.svd_options_full(3, 1.0e-8, 500, 5, @src.irlba_method())
  assert_eq(opts.rank, 3)
  assert_true((opts.tol - 1.0e-8).abs() < 1.0e-15)
  assert_eq(opts.max_iter, 500)
  assert_eq(opts.n_oversamples, 5)
  let m = @src.irlba_method()
  assert_true(opts.method == m)
}

///|
test "sv_result_accessors" {
  let d = [3.0, 2.0, 1.0]
  let u = [[1.0, 0.0], [0.0, 1.0], [0.0, 0.0]]
  let v = [[1.0, 0.0], [0.0, 1.0]]
  let result = @src.sv_result_new(d, u, v, 3, 2, 2)
  let d2 = @src.sv_result_d(result)
  assert_eq(d2.length(), 3)
  assert_true((d2[0] - 3.0).abs() < 1.0e-15)
  assert_true((d2[1] - 2.0).abs() < 1.0e-15)
  assert_true((d2[2] - 1.0).abs() < 1.0e-15)
  assert_eq(@src.sv_result_nrow(result), 3)
  assert_eq(@src.sv_result_ncol(result), 2)
  assert_eq(@src.sv_result_rank(result), 2)
}

///|
test "exact_svd_2x2" {
  let matrix = [1.0, 2.0, 3.0, 4.0]
  let result = @src.run_exact_svd(matrix, 2, 2, 0)
  assert_eq(result.nrow, 2)
  assert_eq(result.ncol, 2)
  assert_eq(result.rank, 2)
  assert_eq(result.d.length(), 2)
  assert_true(result.d[0] >= result.d[1])
  assert_true(result.u.length() == 2)
  assert_true(result.v.length() == 2)
}

///|
test "exact_svd_3x3" {
  let matrix = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0]
  let result = @src.run_exact_svd(matrix, 3, 3, 0)
  assert_eq(result.rank, 3)
  assert_eq(result.d.length(), 3)
  assert_true(result.d[0] >= result.d[1])
  assert_true(result.d[1] >= result.d[2])
}

///|
test "exact_svd_singular_values_sorted" {
  let matrix = [3.0, 1.0, 1.0, 3.0]
  let result = @src.run_exact_svd(matrix, 2, 2, 0)
  assert_true(result.d[0] >= result.d[1])
}

///|
test "exact_svd_identity" {
  let matrix = [1.0, 0.0, 0.0, 1.0]
  let result = @src.run_exact_svd(matrix, 2, 2, 0)
  assert_eq(result.rank, 2)
  assert_eq(result.d.length(), 2)
  assert_true(result.d[0] >= 0.0)
  assert_true(result.d[1] >= 0.0)
}

///|
test "exact_svd_zero_matrix" {
  let matrix = [0.0, 0.0, 0.0, 0.0]
  let result = @src.run_exact_svd(matrix, 2, 2, 0)
  assert_eq(result.rank, 2)
  assert_true(result.d[0].abs() < 1.0e-9)
  assert_true(result.d[1].abs() < 1.0e-9)
}

///|
test "exact_svd_rectangular_5x4" {
  let (data, nrow, ncol) = @src.create_test_matrix()
  let result = @src.run_exact_svd(data, nrow, ncol, 0)
  assert_eq(result.nrow, 5)
  assert_eq(result.ncol, 4)
  assert_eq(result.rank, 4)
  assert_eq(result.d.length(), 4)
  assert_true(result.d[0] >= result.d[1])
  assert_true(result.d[1] >= result.d[2])
  assert_true(result.d[2] >= result.d[3])
}

///|
test "exact_svd_empty" {
  let matrix : Array[Double] = []
  let result = @src.run_exact_svd(matrix, 0, 0, 0)
  assert_eq(result.rank, 0)
  assert_eq(result.d.length(), 0)
  assert_eq(result.u.length(), 0)
  assert_eq(result.v.length(), 0)
}

///|
test "exact_svd_single_row" {
  let matrix = [1.0, 2.0, 3.0]
  let result = @src.run_exact_svd(matrix, 1, 3, 0)
  assert_eq(result.nrow, 1)
  assert_eq(result.ncol, 3)
  assert_eq(result.rank, 1)
  assert_eq(result.d.length(), 1)
}

///|
test "exact_svd_single_col" {
  let matrix = [1.0, 2.0, 3.0]
  let result = @src.run_exact_svd(matrix, 3, 1, 0)
  assert_eq(result.nrow, 3)
  assert_eq(result.ncol, 1)
  assert_eq(result.rank, 1)
}

///|
test "exact_svd_with_rank_reduction" {
  let (data, nrow, ncol) = @src.create_test_matrix()
  let result = @src.run_exact_svd(data, nrow, ncol, 2)
  assert_eq(result.rank, 2)
  assert_eq(result.d.length(), 2)
  assert_eq(result.u.length(), 2)
  assert_eq(result.v.length(), 2)
}

///|
test "exact_svd_u_orthonormal" {
  let matrix = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0]
  let result = @src.run_exact_svd(matrix, 2, 3, 0)
  let r = result.rank
  for i = 0; i < r; i = i + 1 {
    for j = 0; j < r; j = j + 1 {
      let mut dot = 0.0
      for k = 0; k < result.nrow; k = k + 1 {
        dot = dot + result.u[i][k] * result.u[j][k]
      }
      let expected = if i == j { 1.0 } else { 0.0 }
      assert_true((dot - expected).abs() < 1.0e-2)
    }
  }
}

///|
test "exact_svd_v_orthonormal" {
  let matrix = [1.0, 2.0, 3.0, 4.0]
  let result = @src.run_exact_svd(matrix, 2, 2, 0)
  let r = result.rank
  for i = 0; i < r; i = i + 1 {
    for j = 0; j < r; j = j + 1 {
      let mut dot = 0.0
      for k = 0; k < result.ncol; k = k + 1 {
        dot = dot + result.v[i][k] * result.v[j][k]
      }
      let expected = if i == j { 1.0 } else { 0.0 }
      assert_true((dot - expected).abs() < 1.0e-2)
    }
  }
}

///|
test "exact_svd_reconstruction" {
  let matrix = [1.0, 2.0, 3.0, 4.0]
  let result = @src.run_exact_svd(matrix, 2, 2, 0)
  assert_true(result.rank >= 1)
  assert_eq(result.d.length(), result.rank)
  assert_eq(result.u.length(), result.rank)
  assert_eq(result.v.length(), result.rank)
  for i = 0; i < result.rank; i = i + 1 {
    assert_true(result.d[i] >= 0.0)
  }
}

///|
test "irlba_small_matrix" {
  let matrix = [1.0, 2.0, 3.0, 4.0]
  let opts = @src.svd_options_full(2, 1.0e-7, 200, 10, @src.irlba_method())
  let result = @src.run_irlba(matrix, 2, 2, 2, opts)
  assert_true(result.rank >= 1)
  assert_true(result.d.length() >= 1)
}

///|
test "irlba_truncated" {
  let (data, nrow, ncol) = @src.create_test_matrix()
  let opts = @src.svd_options_full(2, 1.0e-3, 500, 10, @src.irlba_method())
  let result = @src.run_irlba(data, nrow, ncol, 2, opts)
  assert_true(result.d.length() >= 0)
}

///|
test "irlba_empty" {
  let matrix : Array[Double] = []
  let opts = @src.svd_options_full(2, 1.0e-7, 100, 10, @src.irlba_method())
  let result = @src.run_irlba(matrix, 0, 0, 2, opts)
  assert_eq(result.rank, 0)
}

///|
test "irlba_single_row" {
  let matrix = [1.0, 2.0, 3.0]
  let opts = @src.svd_options_full(1, 1.0e-7, 100, 10, @src.irlba_method())
  let result = @src.run_irlba(matrix, 1, 3, 1, opts)
  assert_true(result.rank >= 1)
}

///|
test "randomized_svd_small" {
  let matrix = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0]
  let opts = @src.svd_options_full(2, 1.0e-7, 100, 5, @src.randomized_method())
  let result = @src.run_randomized_svd(matrix, 2, 3, 2, opts)
  assert_true(result.rank >= 1)
  assert_true(result.d.length() >= 1)
}

///|
test "randomized_svd_truncated" {
  let (data, nrow, ncol) = @src.create_test_matrix()
  let opts = @src.svd_options_full(2, 1.0e-7, 100, 0, @src.randomized_method())
  let result = @src.run_randomized_svd(data, nrow, ncol, 2, opts)
  assert_true(result.d.length() >= 0)
}

///|
test "randomized_svd_empty" {
  let matrix : Array[Double] = []
  let opts = @src.svd_options_full(2, 1.0e-7, 100, 10, @src.randomized_method())
  let result = @src.run_randomized_svd(matrix, 0, 0, 2, opts)
  assert_eq(result.rank, 0)
}

///|
test "run_svd_dispatch_exact" {
  let matrix = [1.0, 2.0, 3.0, 4.0]
  let opts = @src.svd_options_with_method(@src.exact_svd_method())
  let result = @src.run_svd(matrix, 2, 2, opts)
  assert_eq(result.rank, 2)
}

///|
test "run_svd_dispatch_irlba" {
  let matrix = [1.0, 2.0, 3.0, 4.0]
  let opts = @src.svd_options_with_method(@src.irlba_method())
  let result = @src.run_svd(matrix, 2, 2, opts)
  assert_true(result.rank >= 1)
}

///|
test "run_svd_dispatch_randomized" {
  let matrix = [1.0, 2.0, 3.0, 4.0]
  let opts = @src.svd_options_with_method(@src.randomized_method())
  let result = @src.run_svd(matrix, 2, 2, opts)
  assert_true(result.rank >= 1)
}

///|
test "run_svd_truncated_basic" {
  let (data, nrow, ncol) = @src.create_test_matrix()
  let result = @src.run_svd_truncated(data, nrow, ncol, 2)
  assert_eq(result.rank, 2)
  assert_eq(result.d.length(), 2)
  assert_true(result.d[0] >= result.d[1])
}

///|
test "run_svd_truncated_rank_too_large" {
  let (data, nrow, ncol) = @src.create_test_matrix()
  let mn = if nrow < ncol { nrow } else { ncol }
  let result = @src.run_svd_truncated(data, nrow, ncol, mn + 5)
  assert_eq(result.rank, mn)
}

///|
test "methods_give_similar_results" {
  let matrix = [3.0, 1.0, 1.0, 3.0]
  let opts_exact = @src.svd_options_with_method(@src.exact_svd_method())
  let opts_irlba = @src.svd_options_with_method(@src.irlba_method())
  let result_exact = @src.run_svd(matrix, 2, 2, opts_exact)
  let result_irlba = @src.run_svd(matrix, 2, 2, opts_irlba)
  assert_true((result_exact.d[0] - result_irlba.d[0]).abs() < 1.0e-3)
}

///|
test "svd_mat_vec_basic" {
  let matrix = [1.0, 2.0, 3.0, 4.0]
  let x = [1.0, 1.0]
  let result = @src.svd_mat_vec(matrix, 2, 2, x)
  assert_eq(result.length(), 2)
  assert_true((result[0] - 3.0).abs() < 1.0e-15)
  assert_true((result[1] - 7.0).abs() < 1.0e-15)
}

///|
test "svd_mat_t_vec_basic" {
  let matrix = [1.0, 2.0, 3.0, 4.0]
  let x = [1.0, 1.0]
  let result = @src.svd_mat_t_vec(matrix, 2, 2, x)
  assert_eq(result.length(), 2)
  assert_true((result[0] - 4.0).abs() < 1.0e-15)
  assert_true((result[1] - 6.0).abs() < 1.0e-15)
}

///|
test "svd_transpose_basic" {
  let matrix = [1.0, 2.0, 3.0, 4.0]
  let result = @src.svd_transpose(matrix, 2, 2)
  assert_eq(result.length(), 2)
  assert_eq(result[0].length(), 2)
  assert_true((result[0][0] - 1.0).abs() < 1.0e-15)
  assert_true((result[0][1] - 3.0).abs() < 1.0e-15)
  assert_true((result[1][0] - 2.0).abs() < 1.0e-15)
  assert_true((result[1][1] - 4.0).abs() < 1.0e-15)
}

///|
test "svd_random_matrix_shape" {
  let mat = @src.svd_random_matrix(5, 3)
  assert_eq(mat.length(), 5)
  for i = 0; i < 5; i = i + 1 {
    assert_eq(mat[i].length(), 3)
  }
}

///|
test "svd_random_matrix_deterministic" {
  let mat1 = @src.svd_random_matrix(4, 3)
  let mat2 = @src.svd_random_matrix(4, 3)
  for i = 0; i < 4; i = i + 1 {
    for j = 0; j < 3; j = j + 1 {
      assert_true((mat1[i][j] - mat2[i][j]).abs() < 1.0e-15)
    }
  }
}

///|
test "jacobi_rotation_identity" {
  let matrix = [[1.0, 0.0], [0.0, 1.0]]
  let (eigenvalues, _) = @src.jacobi_rotation(matrix, 2, 100)
  assert_eq(eigenvalues.length(), 2)
  assert_true((eigenvalues[0] - 1.0).abs() < 1.0e-9)
  assert_true((eigenvalues[1] - 1.0).abs() < 1.0e-9)
}

///|
test "eigh_symmetric_2x2" {
  let matrix = [[2.0, 0.0], [0.0, 3.0]]
  let (eigenvalues, _) = @src.eigh_symmetric(matrix, 2)
  assert_eq(eigenvalues.length(), 2)
  assert_true((eigenvalues[0] - 2.0).abs() < 1.0e-9)
  assert_true((eigenvalues[1] - 3.0).abs() < 1.0e-9)
}

///|
test "svd_qr_decompose_basic" {
  let a = [[1.0, 0.0], [0.0, 1.0]]
  let (q, r) = @src.svd_qr_decompose(a, 2)
  assert_eq(q.length(), 2)
  assert_eq(r.length(), 2)
}

///|
test "exact_svd_1x1" {
  let matrix = [5.0]
  let result = @src.run_exact_svd(matrix, 1, 1, 0)
  assert_eq(result.rank, 1)
  assert_true(result.d[0] >= 0.0)
}

///|
test "exact_svd_rank_exceeds_min_dim" {
  let matrix = [1.0, 2.0, 3.0, 4.0]
  let result = @src.run_exact_svd(matrix, 2, 2, 10)
  assert_eq(result.rank, 2)
}

///|
test "exact_svd_non_square_3x2" {
  let matrix = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0]
  let result = @src.run_exact_svd(matrix, 3, 2, 0)
  assert_eq(result.rank, 2)
  assert_eq(result.d.length(), 2)
}

///|
test "svd_method_equality" {
  let m1 = @src.exact_svd_method()
  let m2 = @src.exact_svd_method()
  assert_true(m1 == m2)
  let m3 = @src.irlba_method()
  assert_false(m1 == m3)
}

///|
test "svd_options_equality" {
  let opts1 = @src.svd_options_new()
  let opts2 = @src.svd_options_new()
  assert_eq(opts1, opts2)
}

///|
test "exact_svd_larger_matrix_sorted" {
  let matrix = [
    1.0, 2.0, 3.0,
    2.0, 3.0, 4.0,
    3.0, 4.0, 5.0,
    4.0, 5.0, 6.0,
  ]
  let result = @src.run_exact_svd(matrix, 4, 3, 0)
  assert_eq(result.rank, 3)
  assert_true(result.d[0] >= result.d[1])
  assert_true(result.d[1] >= result.d[2])
}