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
  let at = @src.matrix_transpose(a)

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
    [3.0, 4.0, 5.0, 6.0],
  ])
  let design = @src.Matrix::from_2d(4, 2, [
    [1.0, 0.0],
    [1.0, 0.0],
    [1.0, 1.0],
    [1.0, 1.0],
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
    [10.0, 11.0, 12.0, 13.0],
  ])
  let design = @src.Matrix::from_2d(4, 2, [
    [1.0, 0.0],
    [1.0, 0.0],
    [1.0, 1.0],
    [1.0, 1.0],
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
    [10.0, 11.0, 12.0, 13.0],
  ])
  let design = @src.Matrix::from_2d(4, 2, [
    [1.0, 0.0],
    [1.0, 0.0],
    [1.0, 1.0],
    [1.0, 1.0],
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
    [10.0, 15.0, 80.0, 90.0],
  ])
  let design = @src.Matrix::from_2d(4, 2, [
    [1.0, 0.0],
    [1.0, 0.0],
    [1.0, 1.0],
    [1.0, 1.0],
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
    [5.0, 6.0, 7.0, 8.0],
  ])
  let design = @src.Matrix::from_2d(4, 3, [
    [1.0, 0.0, 0.0],
    [1.0, 1.0, 0.0],
    [1.0, 0.0, 1.0],
    [1.0, 1.0, 1.0],
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
    [1.0, 1.0, 1.0],
  ])

  let contrast_list = [[0.0, -1.0, 1.0]]
  let contrasts = @src.make_contrasts(contrast_list, design)

  assert_eq(contrasts.rows, 3)
  assert_eq(contrasts.cols, 1)
  assert_true(contrasts.get(1, 0) == -1.0)
  assert_true(contrasts.get(2, 0) == 1.0)
}

///|
test "limma_remove_batch_effect" {
  let y = @src.Matrix::from_2d(5, 6, [
    [1.0, 2.0, 3.0, 10.0, 11.0, 12.0],
    [2.0, 3.0, 4.0, 11.0, 12.0, 13.0],
    [3.0, 4.0, 5.0, 12.0, 13.0, 14.0],
    [4.0, 5.0, 6.0, 13.0, 14.0, 15.0],
    [5.0, 6.0, 7.0, 14.0, 15.0, 16.0],
  ])
  let batch = ["B1", "B1", "B1", "B2", "B2", "B2"]
  let group = ["Control", "Control", "Control", "Treatment", "Treatment", "Treatment"]
  let result = @src.limma_remove_batch_effect(y, batch, group = group)
  assert_eq(result.rows, 5)
  assert_eq(result.cols, 6)
}

///|
test "limma_treat" {
  let y = @src.Matrix::from_2d(5, 4, [
    [1.0, 2.0, 3.0, 4.0],
    [2.0, 3.0, 4.0, 5.0],
    [3.0, 4.0, 5.0, 6.0],
    [4.0, 5.0, 6.0, 7.0],
    [5.0, 6.0, 7.0, 8.0],
  ])
  let design = @src.Matrix::from_2d(4, 2, [
    [1.0, 0.0],
    [1.0, 0.0],
    [1.0, 1.0],
    [1.0, 1.0],
  ])
  let fit = @src.lm_fit(y, design)
  let eb_result = @src.ebayes(fit, 0.01)
  let result = @src.treat(eb_result, 0.05, "both")
  assert_eq(result.t.rows, 5)
  assert_eq(result.p_value.rows, 5)
}

///|
test "limma_plot_mds_distances" {
  let y = @src.Matrix::from_2d(5, 6, [
    [1.0, 2.0, 3.0, 10.0, 11.0, 12.0],
    [2.0, 3.0, 4.0, 11.0, 12.0, 13.0],
    [3.0, 4.0, 5.0, 12.0, 13.0, 14.0],
    [4.0, 5.0, 6.0, 13.0, 14.0, 15.0],
    [5.0, 6.0, 7.0, 14.0, 15.0, 16.0],
  ])
  let result = @src.plot_mds_distances(y, 50)
  assert_eq(result.rows, 6)
  assert_eq(result.cols, 6)
}

///|
test "limma_rpkm_basic" {
  let counts = @src.Matrix::from_2d(3, 2, [
    [10.0, 20.0],
    [30.0, 40.0],
    [50.0, 60.0],
  ])
  let gene_lengths = [1000.0, 2000.0, 3000.0]
  let result = @src.limma_rpkm(counts, gene_lengths)
  assert_eq(result.rows, 3)
  assert_eq(result.cols, 2)
  // First gene, first sample: 10 * 1e9 / (150 * 1) = 66666666.67
  assert_true(result.get(0, 0) > 1.0e7)
}

///|
test "limma_rpkm_empty_lengths" {
  let counts = @src.Matrix::from_2d(2, 2, [
    [10.0, 20.0],
    [30.0, 40.0],
  ])
  let gene_lengths : Array[Double] = Array::new()
  let result = @src.limma_rpkm(counts, gene_lengths)
  assert_eq(result.rows, 2)
  assert_eq(result.cols, 2)
}

///|
test "limma_cpm_basic" {
  let counts = @src.Matrix::from_2d(3, 2, [
    [10.0, 20.0],
    [30.0, 40.0],
    [50.0, 60.0],
  ])
  let result = @src.limma_cpm(counts)
  assert_eq(result.rows, 3)
  assert_eq(result.cols, 2)
}

///|
test "limma_cpm_no_log" {
  let counts = @src.Matrix::from_2d(3, 2, [
    [10.0, 20.0],
    [30.0, 40.0],
    [50.0, 60.0],
  ])
  let result = @src.limma_cpm(counts, log = false)
  assert_eq(result.rows, 3)
  // No log transform - values should be counts per million
  // Sample 0 total = 90, so 10/90 * 1e6 ~ 111111.11
  assert_true(result.get(0, 0) > 100000.0)
}

///|
test "limma_quantile_normalize_basic" {
  let y = @src.Matrix::from_2d(5, 3, [
    [1.0, 2.0, 3.0],
    [2.0, 4.0, 6.0],
    [3.0, 6.0, 9.0],
    [4.0, 8.0, 12.0],
    [5.0, 10.0, 15.0],
  ])
  let result = @src.limma_quantile_normalize(y)
  assert_eq(result.rows, 5)
  assert_eq(result.cols, 3)
}

///|
test "limma_quantile_normalize_single" {
  let y = @src.Matrix::from_2d(3, 1, [
    [1.0],
    [2.0],
    [3.0],
  ])
  let result = @src.limma_quantile_normalize(y)
  assert_eq(result.rows, 3)
  assert_eq(result.cols, 1)
}

///|
test "limma_combat_basic" {
  let y = @src.Matrix::from_2d(5, 6, [
    [1.0, 2.0, 3.0, 10.0, 11.0, 12.0],
    [2.0, 3.0, 4.0, 11.0, 12.0, 13.0],
    [3.0, 4.0, 5.0, 12.0, 13.0, 14.0],
    [4.0, 5.0, 6.0, 13.0, 14.0, 15.0],
    [5.0, 6.0, 7.0, 14.0, 15.0, 16.0],
  ])
  let batch = ["B1", "B1", "B1", "B2", "B2", "B2"]
  let result = @src.limma_combat(y, batch)
  assert_eq(result.rows, 5)
  assert_eq(result.cols, 6)
}

///|
test "limma_combat_single_batch" {
  let y = @src.Matrix::from_2d(3, 4, [
    [1.0, 2.0, 3.0, 4.0],
    [5.0, 6.0, 7.0, 8.0],
    [9.0, 10.0, 11.0, 12.0],
  ])
  let batch = ["B1", "B1", "B1", "B1"]
  let result = @src.limma_combat(y, batch)
  // With single batch, should return unchanged
  assert_eq(result.rows, 3)
  assert_eq(result.cols, 4)
  assert_true(result.get(0, 0) == 1.0)
}

///|
test "limma_combat_with_group" {
  let y = @src.Matrix::from_2d(3, 6, [
    [1.0, 2.0, 3.0, 10.0, 11.0, 12.0],
    [2.0, 3.0, 4.0, 11.0, 12.0, 13.0],
    [3.0, 4.0, 5.0, 12.0, 13.0, 14.0],
  ])
  let batch = ["B1", "B1", "B1", "B2", "B2", "B2"]
  let group = ["C", "C", "C", "T", "T", "T"]
  let result = @src.limma_combat(y, batch, group = group)
  assert_eq(result.rows, 3)
  assert_eq(result.cols, 6)
}
