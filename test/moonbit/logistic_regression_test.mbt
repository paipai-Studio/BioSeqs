///|
/// Tests for the LogisticRegression module (Biopython Bio.LogisticRegression).
///
/// Covers model construction, Newton-Raphson training, classification,
/// probability calculation, edge cases, and sample data validation.

// ===========================================================================
// 1. Model construction
// ===========================================================================

///|
test "lr_model_creation" {
  let model = @src.LogisticRegressionModel::new([0.0, 1.0, -1.0])
  assert_eq(model.n_weights(), 3)
  assert_true(model.beta()[0] == 0.0)
  assert_true(model.beta()[1] == 1.0)
  assert_true(model.beta()[2] == -1.0)
}

///|
test "lr_model_empty" {
  let model = @src.LogisticRegressionModel::new([])
  assert_eq(model.n_weights(), 0)
}

// ===========================================================================
// 2. Sample data
// ===========================================================================

///|
test "lr_sample_data_length" {
  let (xs, ys) = @src.logistic_regression_sample_data()
  assert_eq(xs.length(), 17)
  assert_eq(ys.length(), 17)
}

///|
test "lr_sample_data_dimensions" {
  let (xs, _) = @src.logistic_regression_sample_data()
  for x in xs {
    assert_eq(x.length(), 2)
  }
}

///|
test "lr_sample_data_labels" {
  let (_, ys) = @src.logistic_regression_sample_data()
  let mut n_op = 0
  let mut n_nop = 0
  for y in ys {
    if y == 1 {
      n_op = n_op + 1
    } else if y == 0 {
      n_nop = n_nop + 1
    }
  }
  assert_eq(n_op, 10)
  assert_eq(n_nop, 7)
}

// ===========================================================================
// 3. Training
// ===========================================================================

///|
test "lr_train_returns_model" {
  let (xs, ys) = @src.logistic_regression_sample_data()
  let model = @src.logistic_regression_train(xs, ys)
  // beta should have 3 weights: bias + 2 features.
  assert_eq(model.n_weights(), 3)
}

///|
test "lr_train_beta_nonzero" {
  let (xs, ys) = @src.logistic_regression_sample_data()
  let model = @src.logistic_regression_train(xs, ys)
  // At least one weight should be non-zero after training.
  let mut has_nonzero = false
  for b in model.beta() {
    if b != 0.0 {
      has_nonzero = true
    }
  }
  assert_true(has_nonzero)
}

///|
test "lr_train_simple_separable" {
  // Simple linearly separable data.
  let xs = [
    [0.0, 0.0],
    [0.0, 1.0],
    [1.0, 0.0],
    [1.0, 1.0],
    [5.0, 5.0],
    [6.0, 5.0],
    [5.0, 6.0],
    [6.0, 6.0],
  ]
  let ys = [0, 0, 0, 0, 1, 1, 1, 1]
  let model = @src.logistic_regression_train(xs, ys)
  // Classify training data - should mostly match.
  let mut correct = 0
  for i in 0..<xs.length() {
    let pred = @src.logistic_regression_classify(model, xs[i])
    if pred == ys[i] {
      correct = correct + 1
    }
  }
  // At least 6 out of 8 should be correct.
  assert_true(correct >= 6)
}

// ===========================================================================
// 4. Classification
// ===========================================================================

///|
test "lr_classify_returns_0_or_1" {
  let (xs, ys) = @src.logistic_regression_sample_data()
  let model = @src.logistic_regression_train(xs, ys)
  for x in xs {
    let pred = @src.logistic_regression_classify(model, x)
    assert_true(pred == 0 || pred == 1)
  }
}

///|
test "lr_classify_operon_data" {
  let (xs, ys) = @src.logistic_regression_sample_data()
  let model = @src.logistic_regression_train(xs, ys)
  // Check that training data is classified reasonably well.
  let mut correct = 0
  for i in 0..<xs.length() {
    let pred = @src.logistic_regression_classify(model, xs[i])
    if pred == ys[i] {
      correct = correct + 1
    }
  }
  // At least 12 out of 17 correct.
  assert_true(correct >= 12)
}

// ===========================================================================
// 5. Probability calculation
// ===========================================================================

///|
test "lr_calculate_returns_two_probs" {
  let (xs, _) = @src.logistic_regression_sample_data()
  let ys = [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0]
  let model = @src.logistic_regression_train(xs, ys)
  let probs = @src.logistic_regression_calculate(model, xs[0])
  assert_eq(probs.length(), 2)
}

///|
test "lr_calculate_probs_sum_to_one" {
  let (xs, _) = @src.logistic_regression_sample_data()
  let ys = [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0]
  let model = @src.logistic_regression_train(xs, ys)
  for x in xs {
    let probs = @src.logistic_regression_calculate(model, x)
    let sum = probs[0] + probs[1]
    assert_true(sum > 0.99)
    assert_true(sum < 1.01)
  }
}

///|
test "lr_calculate_probs_in_range" {
  let (xs, _) = @src.logistic_regression_sample_data()
  let ys = [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0]
  let model = @src.logistic_regression_train(xs, ys)
  for x in xs {
    let probs = @src.logistic_regression_calculate(model, x)
    for p in probs {
      assert_true(p >= 0.0)
      assert_true(p <= 1.0)
    }
  }
}

// ===========================================================================
// 6. Edge cases
// ===========================================================================

///|
test "lr_train_single_feature" {
  let xs = [[1.0], [2.0], [3.0], [10.0], [11.0], [12.0]]
  let ys = [0, 0, 0, 1, 1, 1]
  let model = @src.logistic_regression_train(xs, ys)
  // Should classify correctly.
  assert_eq(@src.logistic_regression_classify(model, [2.0]), 0)
  assert_eq(@src.logistic_regression_classify(model, [11.0]), 1)
}

///|
test "lr_train_all_same_class" {
  let xs = [[1.0, 2.0], [3.0, 4.0], [5.0, 6.0]]
  let ys = [1, 1, 1]
  let model = @src.logistic_regression_train(xs, ys)
  // All predictions should be 1.
  for x in xs {
    assert_eq(@src.logistic_regression_classify(model, x), 1)
  }
}

// ===========================================================================
// 7. Model accessors
// ===========================================================================

///|
test "lr_beta_accessor" {
  let beta = [1.0, -2.0, 3.0]
  let model = @src.LogisticRegressionModel::new(beta)
  let b = model.beta()
  assert_eq(b.length(), 3)
  assert_true(b[0] == 1.0)
  assert_true(b[1] == -2.0)
  assert_true(b[2] == 3.0)
}

///|
test "lr_n_weights_accessor" {
  let model = @src.LogisticRegressionModel::new([0.0, 0.0, 0.0, 0.0])
  assert_eq(model.n_weights(), 4)
}
