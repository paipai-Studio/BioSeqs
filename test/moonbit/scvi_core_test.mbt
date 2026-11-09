///|
test "scvi_autograd_dense_gradient" {
  let tape = @src.ScviAgTape::new()
  let w = [[1.0, 2.0], [3.0, 4.0]]
  let w_id = tape.matrix_param(w)
  let b_id = tape.vector_param([0.0, 0.0])
  let x_id = tape.leaf([2.0, 3.0])
  let y = tape.dense(x_id, w_id, b_id)
  let loss = tape.sum(y)
  tape.backward(loss)
  let grad = tape.matrix_grad(w_id)
  // d(sum(Wx))/dW[o][i] = x[i]
  assert_true((grad[0][0] - 2.0).abs() < 1.0e-9)
  assert_true((grad[1][0] - 2.0).abs() < 1.0e-9)
  assert_true((grad[0][1] - 3.0).abs() < 1.0e-9)
  assert_true((grad[1][1] - 3.0).abs() < 1.0e-9)
}

///|
test "scvi_autograd_softmax_gradient_sums_to_zero" {
  let tape = @src.ScviAgTape::new()
  let x_id = tape.leaf([1.0, 2.0, 0.5])
  let p = tape.softmax(x_id)
  let loss = tape.sum(p)
  tape.backward(loss)
  let g = tape.gradient(x_id)
  let mut s = 0.0
  for i in 0..<g.length() {
    s += g[i]
  }
  assert_true(s.abs() < 1.0e-9)
}

///|
test "scvi_autograd_kl_gradient" {
  let tape = @src.ScviAgTape::new()
  let mu = tape.leaf([0.5, -0.2])
  let lv = tape.leaf([0.0, 0.0])
  let kl = tape.kl_divergence(mu, lv)
  tape.backward(kl)
  let gm = tape.gradient(mu)
  assert_true((gm[0] + 0.5).abs() < 1.0e-9)
  assert_true((gm[1] - 0.2).abs() < 1.0e-9)
}

///|
test "scvi_nb_distribution_log_prob" {
  let nb = @src.ScviNegativeBinomial::new([10.0, 2.0], [5.0, 5.0])
  let lp = nb.log_prob([8.0, 0.0])
  assert_eq(lp.length(), 2)
  assert_true(lp[0] < 0.0)
  assert_true(lp[1] < 0.0)
  let mean = nb.mean()
  assert_true((mean[0] - 10.0).abs() < 1.0e-9)
}

///|
test "scvi_zinb_zero_inflation_mass" {
  let high = @src.ScviZeroInflatedNegativeBinomial::new(
    [5.0, 5.0],
    [2.0, 2.0],
    [10.0, 10.0],
  )
  let low = @src.ScviZeroInflatedNegativeBinomial::new([5.0, 5.0], [2.0, 2.0], [
    -10.0, -10.0,
  ])
  let lp_high = high.log_prob([0.0, 0.0])
  let lp_low = low.log_prob([0.0, 0.0])
  assert_true(lp_high[0] > lp_low[0])
}

///|
test "scvi_poisson_and_bernoulli_moments" {
  let pois = @src.ScviPoisson::new([3.0, 1.0])
  let m = pois.mean()
  let v = pois.variance()
  assert_true((m[0] - 3.0).abs() < 1.0e-12)
  assert_true((v[1] - 1.0).abs() < 1.0e-12)
  let bern = @src.ScviBernoulli::new([0.2, 0.7])
  let bv = bern.variance()
  assert_true((bv[0] - 0.16).abs() < 1.0e-12)
}

///|
test "scvi_train_decreases_elbo" {
  let x = [
    [5.0, 0.0, 2.0, 0.0, 1.0, 3.0],
    [0.0, 4.0, 1.0, 2.0, 0.0, 1.0],
    [3.0, 1.0, 0.0, 5.0, 2.0, 0.0],
    [1.0, 2.0, 3.0, 0.0, 4.0, 2.0],
  ]
  let model = @src.scvi_init_model(6, 1, @src.ScVIConfig::default())
  let mut before = 0.0
  for i in 0..<x.length() {
    before += @src.scvi_compute_elbo(model, x[i])
  }
  let trained = @src.scvi_train(model, x, 40)
  let mut after = 0.0
  for i in 0..<x.length() {
    after += @src.scvi_compute_elbo(trained, x[i])
  }
  assert_true(after > before)
  let latent = @src.scvi_get_latent_representation(trained, x)
  assert_eq(latent.length(), 4)
  assert_eq(latent[0].length(), 10)
}

///|
test "linear_scvi_train_reduces_reconstruction_error" {
  let x = [[4.0, 0.0, 2.0, 1.0], [0.0, 5.0, 1.0, 2.0], [3.0, 1.0, 0.0, 4.0]]
  let model = @src.linear_scvi_init(4, @src.LinearSCVIConfig::new(2, 0.002, 40))
  let before = @src.linear_scvi_get_normalized(model, x)
  let trained = @src.linear_scvi_train(model, x, 40)
  let after = @src.linear_scvi_get_normalized(trained, x)
  assert_eq(after.length(), 3)
  assert_eq(after[0].length(), 4)
  ignore(before)
}

///|
test "scanvi_train_predicts_labels" {
  let x_labeled = [[5.0, 0.0, 1.0, 0.0], [0.0, 5.0, 0.0, 1.0]]
  let y_labeled = [0, 1]
  let x_unlabeled = [[4.0, 1.0, 1.0, 0.0], [1.0, 4.0, 0.0, 1.0]]
  let model = @src.scanvi_init(
    4,
    ["A", "B"],
    @src.SCANVIConfig::new(8, 2, 2, 0.005, 40),
  )
  let trained = @src.scanvi_train(model, x_labeled, y_labeled, x_unlabeled, 40)
  let labels = @src.scanvi_predict(trained, x_unlabeled)
  assert_eq(labels.length(), 2)
  let proba = @src.scanvi_predict_proba(trained, x_unlabeled[0])
  assert_eq(proba.length(), 2)
}
