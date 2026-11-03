///|
/// Tests for the MaxEntropy module (Biopython Bio.MaxEntropy).
///
/// Covers feature construction, evaluation, model training, classification,
/// probability calculation, edge cases, and sample data validation.

// ===========================================================================
// 1. Feature construction
// ===========================================================================

///|
test "me_feature_creation" {
  let f = @src.MaxEntFeature::new(0, "Red", "Yes")
  assert_eq(f.attribute_index(), 0)
  assert_eq(f.attribute_value(), "Red")
  assert_eq(f.class_label(), "Yes")
}

///|
test "me_feature_eval_match" {
  let f = @src.MaxEntFeature::new(0, "Red", "Yes")
  let obs = ["Red", "Sports", "Domestic"]
  let val = @src.max_ent_eval_feature(f, obs, "Yes")
  assert_eq(val, 1)
}

///|
test "me_feature_eval_no_match_value" {
  let f = @src.MaxEntFeature::new(0, "Red", "Yes")
  let obs = ["Yellow", "Sports", "Domestic"]
  let val = @src.max_ent_eval_feature(f, obs, "Yes")
  assert_eq(val, 0)
}

///|
test "me_feature_eval_no_match_class" {
  let f = @src.MaxEntFeature::new(0, "Red", "Yes")
  let obs = ["Red", "Sports", "Domestic"]
  let val = @src.max_ent_eval_feature(f, obs, "No")
  assert_eq(val, 0)
}

///|
test "me_feature_eval_index_out_of_range" {
  let f = @src.MaxEntFeature::new(5, "Red", "Yes")
  let obs = ["Red", "Sports"]
  let val = @src.max_ent_eval_feature(f, obs, "Yes")
  assert_eq(val, 0)
}

// ===========================================================================
// 2. Sample data
// ===========================================================================

///|
test "me_sample_data_length" {
  let (xs, ys, features) = @src.max_entropy_sample_data()
  assert_eq(xs.length(), 10)
  assert_eq(ys.length(), 10)
  assert_true(features.length() > 0)
}

///|
test "me_sample_data_observation_dim" {
  let (xs, _, _) = @src.max_entropy_sample_data()
  for x in xs {
    assert_eq(x.length(), 3)
  }
}

///|
test "me_sample_data_classes" {
  let (_, ys, _) = @src.max_entropy_sample_data()
  let mut has_yes = false
  let mut has_no = false
  for y in ys {
    if y == "Yes" {
      has_yes = true
    } else if y == "No" {
      has_no = true
    }
  }
  assert_true(has_yes)
  assert_true(has_no)
}

///|
test "me_sample_data_features" {
  let (_, _, features) = @src.max_entropy_sample_data()
  // 3 attributes × 2 values × 2 classes = 12 features.
  assert_eq(features.length(), 12)
}

// ===========================================================================
// 3. Model training
// ===========================================================================

///|
test "me_train_returns_model" {
  let (xs, ys, features) = @src.max_entropy_sample_data()
  let model = @src.max_entropy_train(xs, ys, features)
  assert_true(model.classes().length() >= 2)
  assert_eq(model.alphas().length(), features.length())
  assert_eq(model.features().length(), features.length())
}

///|
test "me_train_alphas_modified" {
  let (xs, ys, features) = @src.max_entropy_sample_data()
  let model = @src.max_entropy_train(xs, ys, features)
  // At least one alpha should be non-zero after training.
  let mut has_nonzero = false
  for a in model.alphas() {
    if a != 0.0 {
      has_nonzero = true
    }
  }
  assert_true(has_nonzero)
}

///|
test "me_train_classes_contain_labels" {
  let (xs, ys, features) = @src.max_entropy_sample_data()
  let model = @src.max_entropy_train(xs, ys, features)
  let classes = model.classes()
  assert_true(classes.contains("Yes"))
  assert_true(classes.contains("No"))
}

// ===========================================================================
// 4. Classification
// ===========================================================================

///|
test "me_classify_returns_valid_class" {
  let (xs, ys, features) = @src.max_entropy_sample_data()
  let model = @src.max_entropy_train(xs, ys, features)
  for x in xs {
    let pred = @src.max_entropy_classify(model, x)
    assert_true(pred == "Yes" || pred == "No")
  }
}

///|
test "me_classify_new_observation" {
  let (xs, ys, features) = @src.max_entropy_sample_data()
  let model = @src.max_entropy_train(xs, ys, features)
  let obs = ["Red", "Sports", "Domestic"]
  let pred = @src.max_entropy_classify(model, obs)
  assert_true(pred == "Yes" || pred == "No")
}

// ===========================================================================
// 5. Probability calculation
// ===========================================================================

///|
test "me_calculate_returns_scores" {
  let (xs, ys, features) = @src.max_entropy_sample_data()
  let model = @src.max_entropy_train(xs, ys, features)
  let scores = @src.max_entropy_calculate(model, xs[0])
  assert_eq(scores.length(), model.classes().length())
}

///|
test "me_n_features" {
  let (xs, ys, features) = @src.max_entropy_sample_data()
  let model = @src.max_entropy_train(xs, ys, features)
  assert_eq(model.n_features(), features.length())
}

// ===========================================================================
// 6. Feature accessors
// ===========================================================================

///|
test "me_feature_accessors" {
  let f = @src.MaxEntFeature::new(2, "SUV", "No")
  assert_eq(f.attribute_index(), 2)
  assert_eq(f.attribute_value(), "SUV")
  assert_eq(f.class_label(), "No")
}

// ===========================================================================
// 7. Edge cases
// ===========================================================================

///|
test "me_train_single_class" {
  let xs = [["A"], ["A"], ["A"]]
  let ys = ["Yes", "Yes", "Yes"]
  let features = [@src.MaxEntFeature::new(0, "A", "Yes")]
  let model = @src.max_entropy_train(xs, ys, features)
  assert_eq(model.classes().length(), 1)
  assert_eq(@src.max_entropy_classify(model, ["A"]), "Yes")
}

///|
test "me_train_empty_features" {
  let xs = [["A"], ["B"]]
  let ys = ["Yes", "No"]
  let features : Array[@src.MaxEntFeature] = []
  let model = @src.max_entropy_train(xs, ys, features)
  assert_eq(model.n_features(), 0)
}

// ===========================================================================
// 8. Model accessors
// ===========================================================================

///|
test "me_model_accessors" {
  let (xs, ys, features) = @src.max_entropy_sample_data()
  let model = @src.max_entropy_train(xs, ys, features)
  let classes = model.classes()
  let alphas = model.alphas()
  let feats = model.features()
  assert_true(classes.length() >= 2)
  assert_true(alphas.length() > 0)
  assert_true(feats.length() > 0)
}
