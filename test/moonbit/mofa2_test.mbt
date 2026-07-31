///|
/// Tests for MOFA2 multi-omics factor analysis module.

test "mofa_create_example" {
  let views = @src.mofa_create_example(30, 3)
  assert_true(views.length() == 2)
  assert_true(views[0].view_name == "rna")
  assert_true(views[1].view_name == "atac")
  assert_true(views[0].data.length() == 30)
  assert_true(views[0].data[0].length() == 50)
  assert_true(views[1].data[0].length() == 30)
}

test "mofa_params_default" {
  let params = @src.MofaParams::new()
  assert_true(params.n_factors == 10)
  assert_true(params.max_iterations == 100)
  assert_true(params.tolerance > 0.0)
  assert_true(params.seed == 42)
}

test "mofa_view_creation" {
  let data : Array[Array[Double]] = Array::new()
  data.push([1.0, 2.0, 3.0])
  data.push([4.0, 5.0, 6.0])
  let view = @src.MofaView::new("test", data, ["f1", "f2", "f3"], ["s1", "s2"], view_type = "gaussian")
  assert_true(view.view_name == "test")
  assert_true(view.data.length() == 2)
  assert_true(view.feature_names.length() == 3)
  assert_true(view.view_type == "gaussian")
}

test "mofa_run_basic" {
  let views = @src.mofa_create_example(20, 3)
  let params = @src.MofaParams::create(
    n_factors = 3,
    max_iterations = 20,
  )
  let result = @src.mofa_run(views, params)
  assert_true(result.factors.length() == 20)
  assert_true(result.loadings.length() == 2)
  assert_true(result.factor_vars.length() == 3)
  assert_true(result.n_iterations > 0)
}

test "mofa_convergence" {
  let views = @src.mofa_create_example(25, 2)
  let params = @src.MofaParams::create(
    n_factors = 2,
    max_iterations = 50,
  )
  let result = @src.mofa_run(views, params)
  assert_true(result.n_iterations <= 50)
}

test "mofa_variance_explained" {
  let views = @src.mofa_create_example(30, 4)
  let params = @src.MofaParams::create(
    n_factors = 4,
    max_iterations = 30,
  )
  let result = @src.mofa_run(views, params)
  // Variance should be non-negative
  for v in result.factor_vars {
    assert_true(v >= 0.0)
  }
  // Each view's variance explained should be non-negative
  let mut vi = 0
  while vi < result.view_vars.length() {
    let mut sum_v = 0.0
    for v in result.view_vars[vi] {
      assert_true(v >= 0.0)
      sum_v = sum_v + v
    }
    assert_true(sum_v >= 0.0)
    vi = vi + 1
  }
}

test "mofa_active_factors" {
  let views = @src.mofa_create_example(20, 2)
  let params = @src.MofaParams::create(
    n_factors = 5,
    max_iterations = 30,
  )
  let result = @src.mofa_run(views, params)
  assert_true(result.active_factors >= 0)
  assert_true(result.active_factors <= 5)
}

test "mofa_single_view" {
  // Should work with a single view
  let data : Array[Array[Double]] = Array::new()
  let mut s = 0
  while s < 15 {
    let row : Array[Double] = Array::make(10, 0.0)
    let mut f = 0
    while f < 10 {
      row[f] = (s.to_double() + f.to_double()) * 0.1
      f = f + 1
    }
    data.push(row)
    s = s + 1
  }
  let view = @src.MofaView::new("single", data, [], [])
  let views : Array[@src.MofaView] = Array::new()
  views.push(view)
  let params = @src.MofaParams::create(
    n_factors = 2,
    max_iterations = 15,
  )
  let result = @src.mofa_run(views, params)
  assert_true(result.factors.length() == 15)
  assert_true(result.loadings.length() == 1)
}

test "mofa_empty_views" {
  let views : Array[@src.MofaView] = Array::new()
  let params = @src.MofaParams::create(
    n_factors = 2,
  )
  let result = @src.mofa_run(views, params)
  assert_true(result.factors.length() == 0)
  assert_true(result.converged == false)
}
