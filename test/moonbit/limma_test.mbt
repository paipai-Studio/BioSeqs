///|
/// Tests for limma module.

///|
test "limma_matrix_mult" {
  let a = @src.Matrix::from_2d(2, 2, [[1.0, 2.0], [3.0, 4.0]])
  let b = @src.Matrix::from_2d(2, 2, [[5.0, 6.0], [7.0, 8.0]])
  let c = @src.mat_mult(a, b)
  
  assert_true(c.get(0, 0) == 19.0)
  assert_true(c.get(0, 1) == 22.0)
  assert_true(c.get(1, 0) == 43.0)
  assert_true(c.get(1, 1) == 50.0)
}

///|
test "limma_matrix_transpose" {
  let a = @src.Matrix::from_2d(2, 3, [[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]])
  let at = @src.mat_transpose(a)
  
  assert_eq(at.rows, 3)
  assert_eq(at.cols, 2)
  assert_true(at.get(0, 0) == 1.0)
  assert_true(at.get(0, 1) == 4.0)
  assert_true(at.get(1, 0) == 2.0)
  assert_true(at.get(1, 1) == 5.0)
}

///|
test "limma_matrix_inverse" {
  let a = @src.Matrix::from_2d(2, 2, [[4.0, 7.0], [2.0, 6.0]])
  let a_inv = @src.mat_inverse(a)
  
  let identity = @src.mat_mult(a, a_inv)
  assert_true(identity.get(0, 0).abs() - 1.0 < 0.0001)
  assert_true(identity.get(0, 1).abs() < 0.0001)
  assert_true(identity.get(1, 0).abs() < 0.0001)
  assert_true(identity.get(1, 1).abs() - 1.0 < 0.0001)
}

///|
test "limma_lm_fit" {
  let y = @src.Matrix::from_2d(3, 4, [
    [1.0, 2.0, 3.0, 4.0],
    [2.0, 3.0, 4.0, 5.0],
    [3.0, 4.0, 5.0, 6.0]
  ])
  let design = @src.Matrix::from_2d(4, 2, [
    [1.0, 0.0],
    [1.0, 0.0],
    [1.0, 1.0],
    [1.0, 1.0]
  ])
  
  let fit = @src.lm_fit(y, design)
  
  assert_eq(fit.coefficients.rows, 3)
  assert_eq(fit.coefficients.cols, 2)
  assert_true(fit.sigma[0] > 0.0)
}

///|
test "limma_ebayes" {
  let y = @src.Matrix::from_2d(10, 4, [
    [1.0, 2.0, 3.0, 4.0],
    [2.0, 3.0, 4.0, 5.0],
    [3.0, 4.0, 5.0, 6.0],
    [4.0, 5.0, 6.0, 7.0],
    [5.0, 6.0, 7.0, 8.0],
    [6.0, 7.0, 8.0, 9.0],
    [7.0, 8.0, 9.0, 10.0],
    [8.0, 9.0, 10.0, 11.0],
    [9.0, 10.0, 11.0, 12.0],
    [10.0, 11.0, 12.0, 13.0]
  ])
  let design = @src.Matrix::from_2d(4, 2, [
    [1.0, 0.0],
    [1.0, 0.0],
    [1.0, 1.0],
    [1.0, 1.0]
  ])
  
  let fit = @src.lm_fit(y, design)
  let eb_result = @src.ebayes(fit, 0.01)
  
  assert_eq(eb_result.t.rows, 10)
  assert_eq(eb_result.t.cols, 2)
  assert_eq(eb_result.p_value.rows, 10)
}

///|
test "limma_top_table" {
  let y = @src.Matrix::from_2d(10, 4, [
    [1.0, 2.0, 10.0, 11.0],
    [2.0, 3.0, 4.0, 5.0],
    [3.0, 4.0, 5.0, 6.0],
    [4.0, 5.0, 6.0, 7.0],
    [5.0, 6.0, 7.0, 8.0],
    [6.0, 7.0, 8.0, 9.0],
    [7.0, 8.0, 9.0, 10.0],
    [8.0, 9.0, 10.0, 11.0],
    [9.0, 10.0, 11.0, 12.0],
    [10.0, 11.0, 12.0, 13.0]
  ])
  let design = @src.Matrix::from_2d(4, 2, [
    [1.0, 0.0],
    [1.0, 0.0],
    [1.0, 1.0],
    [1.0, 1.0]
  ])
  
  let fit = @src.lm_fit(y, design)
  let eb_result = @src.ebayes(fit, 0.01)
  let top = @src.top_table(eb_result, 1, 5, "BH")
  
  assert_eq(top.length(), 5)
  assert_true(top[0].adj_p_value <= top[1].adj_p_value)
}

///|
test "limma_adjust_p_values_bh" {
  let p_values = [0.001, 0.01, 0.05, 0.1, 0.2, 0.5, 0.8, 1.0]
  let adj = @src.adjust_p_values(p_values, "BH")
  
  assert_eq(adj.length(), 8)
  assert_true(adj[0] <= adj[1])
  assert_true(adj[0] <= 1.0)
}

///|
test "limma_voom" {
  let counts = @src.Matrix::from_2d(3, 4, [
    [100.0, 120.0, 200.0, 220.0],
    [50.0, 60.0, 55.0, 45.0],
    [10.0, 15.0, 80.0, 90.0]
  ])
  let design = @src.Matrix::from_2d(4, 2, [
    [1.0, 0.0],
    [1.0, 0.0],
    [1.0, 1.0],
    [1.0, 1.0]
  ])
  
  let log_cpm = @src.voom(counts, design)
  
  assert_eq(log_cpm.rows, 3)
  assert_eq(log_cpm.cols, 4)
}

///|
test "limma_contrasts_fit" {
  let y = @src.Matrix::from_2d(5, 4, [
    [1.0, 2.0, 3.0, 4.0],
    [2.0, 3.0, 4.0, 5.0],
    [3.0, 4.0, 5.0, 6.0],
    [4.0, 5.0, 6.0, 7.0],
    [5.0, 6.0, 7.0, 8.0]
  ])
  let design = @src.Matrix::from_2d(4, 3, [
    [1.0, 0.0, 0.0],
    [1.0, 1.0, 0.0],
    [1.0, 0.0, 1.0],
    [1.0, 1.0, 1.0]
  ])
  
  let fit = @src.lm_fit(y, design)
  let eb_result = @src.ebayes(fit, 0.01)
  
  let contrasts = @src.Matrix::from_2d(3, 1, [[0.0], [-1.0], [1.0]])
  let contrast_result = @src.contrasts_fit(eb_result, contrasts)
  
  assert_eq(contrast_result.coefficients.cols, 1)
}

///|
test "limma_make_contrasts" {
  let design = @src.Matrix::from_2d(4, 3, [
    [1.0, 0.0, 0.0],
    [1.0, 1.0, 0.0],
    [1.0, 0.0, 1.0],
    [1.0, 1.0, 1.0]
  ])
  
  let contrast_list = [[0.0, -1.0, 1.0]]
  let contrasts = @src.make_contrasts(contrast_list, design)
  
  assert_eq(contrasts.rows, 3)
  assert_eq(contrasts.cols, 1)
  assert_true(contrasts.get(1, 0) == -1.0)
  assert_true(contrasts.get(2, 0) == 1.0)
}
