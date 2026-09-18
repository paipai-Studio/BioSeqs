///|
/// Tests for the shared linear-algebra primitives (linalg.mbt).

///|
fn linalg_test_close(a : Double, b : Double, tol : Double) -> Bool {
  (a - b).abs() < tol
}

///|
test "linalg_jacobi_eigenvalues" {
  // tridiagonal [[2,1,0],[1,2,1],[0,1,2]] → eigenvalues 2±√2 and 2
  let m = [[2.0, 1.0, 0.0], [1.0, 2.0, 1.0], [0.0, 1.0, 2.0]]
  let (evals, evecs) = @src.linalg_jacobi(m, 3, 100)
  let (sorted, _) = @src.linalg_sort_eigen(evals, evecs, 3)
  assert_true(linalg_test_close(sorted[0], 3.41421356, 1.0e-6))
  assert_true(linalg_test_close(sorted[1], 2.0, 1.0e-9))
  assert_true(linalg_test_close(sorted[2], 0.58578644, 1.0e-6))
}

///|
test "linalg_inverse_2x2" {
  let a = [[4.0, 2.0], [3.0, 1.0]]
  let inv = @src.linalg_inverse(a)
  // det = -2, inverse = [[-0.5, 1], [1.5, -2]]
  assert_true(linalg_test_close(inv[0][0], -0.5, 1.0e-9))
  assert_true(linalg_test_close(inv[0][1], 1.0, 1.0e-9))
  assert_true(linalg_test_close(inv[1][0], 1.5, 1.0e-9))
  assert_true(linalg_test_close(inv[1][1], -2.0, 1.0e-9))
}

///|
test "linalg_inverse_3x3_recovers_identity" {
  let a = [[2.0, 1.0, 0.0], [1.0, 2.0, 1.0], [0.0, 1.0, 2.0]]
  let inv = @src.linalg_inverse(a)
  let prod = @src.linalg_mat_mul(a, inv)
  for i = 0; i < 3; i = i + 1 {
    for j = 0; j < 3; j = j + 1 {
      let expect = if i == j { 1.0 } else { 0.0 }
      assert_true(linalg_test_close(prod[i][j], expect, 1.0e-9))
    }
  }
}

///|
test "linalg_mat_mul_transpose" {
  let a = [[1.0, 2.0], [3.0, 4.0]]
  let b = [[5.0, 6.0], [7.0, 8.0]]
  let c = @src.linalg_mat_mul(a, b)
  assert_true(linalg_test_close(c[0][0], 19.0, 1.0e-9))
  assert_true(linalg_test_close(c[0][1], 22.0, 1.0e-9))
  assert_true(linalg_test_close(c[1][0], 43.0, 1.0e-9))
  assert_true(linalg_test_close(c[1][1], 50.0, 1.0e-9))
  let t = @src.linalg_transpose(a)
  assert_true(linalg_test_close(t[0][1], 3.0, 1.0e-9))
  assert_true(linalg_test_close(t[1][0], 2.0, 1.0e-9))
}