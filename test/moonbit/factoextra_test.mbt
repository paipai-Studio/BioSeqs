///|
/// Tests for factoextra - PCA/factor analysis

///|
test "factoextra_pca_basic" {
  let data = [
    [1.0, 2.0, 3.0],
    [4.0, 5.0, 6.0],
    [7.0, 8.0, 9.0],
    [10.0, 11.0, 12.0],
  ]
  let result = @src.facto_pca(data)
  assert_eq(result.n_ind, 4)
  assert_eq(result.n_var, 3)
  assert_eq(result.n_dims, 3)
}

///|
test "factoextra_pca_eigenvalues" {
  let data = [
    [2.0, 3.0],
    [5.0, 6.0],
    [8.0, 9.0],
  ]
  let result = @src.facto_pca(data)
  let eigenvalues = @src.facto_get_eigenvalue(result)
  assert_eq(eigenvalues.length(), 2)
  assert_true(eigenvalues[0].eigenvalue >= eigenvalues[1].eigenvalue)
}

///|
test "factoextra_pca_variance" {
  let data = [
    [1.0, 2.0],
    [4.0, 5.0],
    [7.0, 8.0],
    [10.0, 11.0],
  ]
  let result = @src.facto_pca(data)
  let eigenvalues = @src.facto_get_eigenvalue(result)
  let mut i = 1
  while i < eigenvalues.length() {
    assert_true(eigenvalues[i].cumulative_variance >= eigenvalues[i - 1].cumulative_variance)
    i = i + 1
  }
}

///|
test "factoextra_pca_individual_coords" {
  let data = [
    [1.0, 2.0, 3.0],
    [4.0, 5.0, 6.0],
    [7.0, 8.0, 9.0],
  ]
  let result = @src.facto_pca(data)
  let ind = @src.facto_get_pca_ind(result)
  assert_eq(ind.coord.length(), 3)
  assert_eq(ind.coord[0].length(), 3)
}

///|
test "factoextra_pca_variable_coords" {
  let data = [
    [1.0, 2.0, 3.0],
    [4.0, 5.0, 6.0],
    [7.0, 8.0, 9.0],
  ]
  let result = @src.facto_pca(data)
  let pca_var = @src.facto_get_pca_var(result)
  assert_eq(pca_var.coord.length(), 3)
  assert_eq(pca_var.coord[0].length(), 3)
}

///|
test "factoextra_pca_ncp" {
  let data = [
    [1.0, 2.0, 3.0, 4.0],
    [5.0, 6.0, 7.0, 8.0],
    [9.0, 10.0, 11.0, 12.0],
  ]
  let result = @src.facto_pca(data, ncp = 2)
  assert_eq(result.n_dims, 2)
  assert_eq(result.eigenvalues.length(), 2)
}

///|
test "factoextra_pca_cos2_individuals" {
  let data = [
    [1.0, 2.0],
    [4.0, 5.0],
    [7.0, 8.0],
  ]
  let result = @src.facto_pca(data)
  let ind = @src.facto_get_pca_ind(result)
  assert_eq(ind.cos2.length(), 3)
  let mut i = 0
  while i < 3 {
    let mut j = 0
    while j < 2 {
      assert_true(ind.cos2[i][j] >= 0.0)
      assert_true(ind.cos2[i][j] <= 1.0)
      j = j + 1
    }
    i = i + 1
  }
}

///|
test "factoextra_pca_contrib_variables" {
  let data = [
    [1.0, 2.0],
    [4.0, 5.0],
    [7.0, 8.0],
  ]
  let result = @src.facto_pca(data)
  let pca_var = @src.facto_get_pca_var(result)
  assert_eq(pca_var.contrib.length(), 2)
  assert_eq(pca_var.contrib[0].length(), 2)
}

///|
test "factoextra_pca_total_inertia" {
  let data = [
    [1.0, 2.0],
    [4.0, 5.0],
    [7.0, 8.0],
  ]
  let result = @src.facto_pca(data)
  let total = @src.facto_total_inertia(result)
  assert_true(total > 0.0)
}

///|
test "factoextra_pca_nb_dim" {
  let data = [
    [1.0, 2.0],
    [4.0, 5.0],
    [7.0, 8.0],
    [10.0, 11.0],
  ]
  let result = @src.facto_pca(data)
  let nb = @src.facto_nb_dim(result, 80.0)
  assert_true(nb >= 1)
  assert_true(nb <= 2)
}

///|
test "factoextra_pca_summary" {
  let data = [
    [1.0, 2.0, 3.0],
    [4.0, 5.0, 6.0],
  ]
  let result = @src.facto_pca(data)
  let summary = @src.facto_summary(result)
  assert_eq(summary.n_individuals, 2)
  assert_eq(summary.n_variables, 3)
  assert_true(summary.total_inertia > 0.0)
}

///|
test "factoextra_pca_empty" {
  let data : Array[Array[Double]] = []
  let result = @src.facto_pca(data)
  assert_eq(result.n_ind, 0)
  assert_eq(result.n_var, 0)
  assert_eq(result.n_dims, 0)
}

///|
test "factoextra_pca_contrib_dim" {
  let data = [
    [1.0, 2.0],
    [4.0, 5.0],
    [7.0, 8.0],
  ]
  let result = @src.facto_pca(data)
  let contribs = @src.facto_var_contrib_dim(result, 1)
  assert_eq(contribs.length(), 2)
  assert_true(contribs[0].1 >= contribs[1].1)
}

///|
test "factoextra_pca_cos2_dim" {
  let data = [
    [1.0, 2.0],
    [4.0, 5.0],
    [7.0, 8.0],
  ]
  let result = @src.facto_pca(data)
  let cos2s = @src.facto_var_cos2_dim(result, 1)
  assert_eq(cos2s.length(), 2)
}

///|
test "factoextra_pca_dimdesc" {
  let data = [
    [1.0, 2.0],
    [4.0, 5.0],
    [7.0, 8.0],
  ]
  let result = @src.facto_pca(data)
  let desc = @src.facto_dimdesc(data, result, 1)
  assert_eq(desc.length(), 2)
}

///|
test "factoextra_pca_loadings" {
  let data = [
    [1.0, 2.0],
    [4.0, 5.0],
    [7.0, 8.0],
  ]
  let result = @src.facto_pca(data)
  let loadings = @src.facto_loadings(result)
  assert_eq(loadings.length(), 2)
}

///|
test "factoextra_pca_print_summary" {
  let data = [
    [1.0, 2.0],
    [4.0, 5.0],
  ]
  let result = @src.facto_pca(data)
  let summary = @src.facto_print_summary(result)
  assert_true(summary.length() > 0)
}
