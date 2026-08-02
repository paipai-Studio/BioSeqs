///|
/// Bio.MaxEntropy - Maximum entropy classifier using Improved Iterative Scaling
/// (IIS), inspired by Biopython Bio.MaxEntropy.
/// Reference: Chang J (2001) Biopython MaxEntropy module; Della Pietra et al.
/// (1997) IEEE Trans PAMI 19:380-393.

// ===========================================================================
// Feature definition
// ===========================================================================

///|
/// An indicator feature for the maximum entropy classifier.
/// The feature evaluates to 1 when the observation at `attribute_index`
/// equals `attribute_value` AND the class label equals `class_label`;
/// otherwise 0.
pub struct MaxEntFeature {
  attribute_index : Int
  attribute_value : String
  class_label : String
}

///|
/// Construct a new indicator feature.
pub fn MaxEntFeature::new(
  attribute_index : Int,
  attribute_value : String,
  class_label : String
) -> MaxEntFeature {
  MaxEntFeature::{ attribute_index, attribute_value, class_label }
}

///|
/// Accessor for the attribute index.
pub fn MaxEntFeature::attribute_index(self : MaxEntFeature) -> Int {
  self.attribute_index
}

///|
/// Accessor for the attribute value.
pub fn MaxEntFeature::attribute_value(self : MaxEntFeature) -> String {
  self.attribute_value
}

///|
/// Accessor for the class label.
pub fn MaxEntFeature::class_label(self : MaxEntFeature) -> String {
  self.class_label
}

///|
/// Evaluate the feature on an observation (array of attribute strings) and a
/// class label. Returns 1 if the feature fires, 0 otherwise.
pub fn max_ent_eval_feature(
  feature : MaxEntFeature,
  observation : Array[String],
  class_label : String
) -> Int {
  if feature.attribute_index >= 0 &&
    feature.attribute_index < observation.length() &&
    observation[feature.attribute_index] == feature.attribute_value &&
    class_label == feature.class_label {
    1
  } else {
    0
  }
}

// ===========================================================================
// Model struct
// ===========================================================================

///|
/// Maximum entropy classifier model.
pub struct MaxEntropyModel {
  classes : Array[String]
  alphas : Array[Double]
  features : Array[MaxEntFeature]
}

///|
/// Accessor for the class labels.
pub fn MaxEntropyModel::classes(self : MaxEntropyModel) -> Array[String] {
  self.classes
}

///|
/// Accessor for the feature weights (alphas).
pub fn MaxEntropyModel::alphas(self : MaxEntropyModel) -> Array[Double] {
  self.alphas
}

///|
/// Accessor for the features.
pub fn MaxEntropyModel::features(
  self : MaxEntropyModel
) -> Array[MaxEntFeature] {
  self.features
}

///|
/// Number of features.
pub fn MaxEntropyModel::n_features(self : MaxEntropyModel) -> Int {
  self.features.length()
}

// ===========================================================================
// Internal helpers
// ===========================================================================

///|
/// Safe exponential with overflow protection.
fn me_safe_exp(x : Double) -> Double {
  if x >= 40.0 {
    return 10000000000000.0 // exp(40) ~ 2.35e17, but clamp
  }
  if x <= -40.0 {
    return 0.0000000001
  }
  @math.exp(x)
}

///|
/// Safe natural log.
fn me_safe_log(x : Double) -> Double {
  if x <= 0.0 {
    return -40.0
  }
  @math.ln(x)
}

///|
/// Compute f^#(x) = sum over all classes of f(x, y), for a given feature
/// and observation. Since our features are indicator functions tied to a
/// specific class, this is either 0 or 1 (it's 1 only if the observation
/// matches and we're summing over the class that the feature refers to).
fn me_f_correction(
  feature : MaxEntFeature,
  observation : Array[String]
) -> Double {
  // f^#(x) = sum_y f(x, y). Since f is 1 only for one class, f^# is at most 1.
  if feature.attribute_index >= 0 &&
    feature.attribute_index < observation.length() &&
    observation[feature.attribute_index] == feature.attribute_value {
    1.0
  } else {
    0.0
  }
}

// ===========================================================================
// Public API
// ===========================================================================

///|
/// Calculate the unnormalized log-probability for each class given an
/// observation. Returns an array of log-probabilities parallel to
/// `model.classes`.
pub fn max_entropy_calculate(
  model : MaxEntropyModel,
  observation : Array[String]
) -> Array[Double] {
  let scores = Array::make(model.classes.length(), 0.0)
  for ci in 0..<model.classes.length() {
    let klass = model.classes[ci]
    let mut lprob = 0.0
    for fi in 0..<model.features.length() {
      let f_val = max_ent_eval_feature(
        model.features[fi], observation, klass
      )
      lprob = lprob + f_val.to_double() * model.alphas[fi]
    }
    scores[ci] = lprob
  }
  scores
}

///|
/// Classify an observation into the class with the highest probability.
pub fn max_entropy_classify(
  model : MaxEntropyModel,
  observation : Array[String]
) -> String {
  let scores = max_entropy_calculate(model, observation)
  let mut best_idx = 0
  let mut best_score = scores[0]
  for i in 1..<scores.length() {
    if scores[i] > best_score {
      best_score = scores[i]
      best_idx = i
    }
  }
  model.classes[best_idx]
}

///|
/// Train a maximum entropy classifier using Improved Iterative Scaling (IIS).
///
/// `training_set` is an array of observations (each observation is an array of
/// categorical attribute strings). `results` is an array of class labels.
/// `features` is an array of `MaxEntFeature` indicator features.
///
/// Returns a trained `MaxEntropyModel`.
pub fn max_entropy_train(
  training_set : Array[Array[String]],
  results : Array[String],
  features : Array[MaxEntFeature]
) -> MaxEntropyModel {
  let n = training_set.length()
  // Collect unique classes.
  let classes : Array[String] = Array::new()
  for r in results {
    if !classes.contains(r) {
      classes.push(r)
    }
  }
  let n_features = features.length()
  let n_classes = classes.length()
  // Initialize alphas to 0.
  let alphas = Array::make(n_features, 0.0)
  // Precompute empirical expectations: E_emp[f_i] = (1/N) * sum f_i(x, y).
  let emp_expects = Array::make(n_features, 0.0)
  for fi in 0..<n_features {
    let mut s = 0.0
    for i in 0..<n {
      s = s + max_ent_eval_feature(
        features[fi], training_set[i], results[i]
      ).to_double()
    }
    emp_expects[fi] = s / n.to_double()
  }
  // IIS iteration.
  let max_iterations = 500
  let converge_threshold = 0.00001
  for _ in 0..<max_iterations {
    let mut max_delta = 0.0
    for fi in 0..<n_features {
      // Skip features with zero empirical expectation.
      if emp_expects[fi] < 0.0000001 {
        continue
      }
      // Compute model expectation of f_i:
      // E_model[f_i] = (1/N) * sum_x sum_y P(y|x) f_i(x, y)
      let mut model_expect = 0.0
      for i in 0..<n {
        // Compute P(y|x) for each class.
        let scores = Array::make(n_classes, 0.0)
        for ci in 0..<n_classes {
          let mut s = 0.0
          for fj in 0..<n_features {
            let f_val = max_ent_eval_feature(
              features[fj], training_set[i], classes[ci]
            )
            s = s + f_val.to_double() * alphas[fj]
          }
          scores[ci] = s
        }
        // Normalize via softmax.
        let mut max_score = scores[0]
        for ci in 1..<n_classes {
          if scores[ci] > max_score {
            max_score = scores[ci]
          }
        }
        let mut sum_exp = 0.0
        for ci in 0..<n_classes {
          scores[ci] = me_safe_exp(scores[ci] - max_score)
          sum_exp = sum_exp + scores[ci]
        }
        // Add contribution of f_i to model expectation.
        for ci in 0..<n_classes {
          let prob = scores[ci] / sum_exp
          let f_val = max_ent_eval_feature(
            features[fi], training_set[i], classes[ci]
          )
          model_expect = model_expect + prob * f_val.to_double()
        }
      }
      model_expect = model_expect / n.to_double()
      if model_expect < 0.0000001 {
        continue
      }
      // Compute f^# for each training instance (correction term).
      // Since f^# is at most 1, use average f^# across training set.
      let mut avg_f_cor = 0.0
      for i in 0..<n {
        avg_f_cor = avg_f_cor + me_f_correction(features[fi], training_set[i])
      }
      avg_f_cor = avg_f_cor / n.to_double()
      if avg_f_cor < 0.0000001 {
        continue
      }
      // Solve for delta_i using Newton's method:
      // emp_expect = model_expect * exp(delta * f^#)
      // => delta = ln(emp_expect / model_expect) / f^#
      let ratio = emp_expects[fi] / model_expect
      let delta = me_safe_log(ratio) / avg_f_cor
      alphas[fi] = alphas[fi] + delta
      let abs_delta = if delta < 0.0 { -delta } else { delta }
      if abs_delta > max_delta {
        max_delta = abs_delta
      }
    }
    if max_delta < converge_threshold {
      break
    }
  }
  MaxEntropyModel::{ classes, alphas, features }
}

// ===========================================================================
// Sample data: car classification (Biopython tutorial example)
// ===========================================================================

///|
/// Sample training data for the car classification example.
/// Observations: [color, type, origin]. Labels: "Yes" / "No".
pub fn max_entropy_sample_data() -> (Array[
  Array[String]
], Array[String], Array[MaxEntFeature]) {
  let xs = [
    ["Red", "Sports", "Domestic"],
    ["Red", "Sports", "Domestic"],
    ["Red", "Sports", "Domestic"],
    ["Yellow", "Sports", "Domestic"],
    ["Yellow", "Sports", "Imported"],
    ["Yellow", "SUV", "Imported"],
    ["Yellow", "SUV", "Imported"],
    ["Yellow", "SUV", "Domestic"],
    ["Red", "SUV", "Imported"],
    ["Red", "Sports", "Imported"],
  ]
  let ys = [
    "Yes", "No", "Yes", "No", "Yes", "No", "Yes", "No", "No", "Yes",
  ]
  let features = [
    MaxEntFeature::new(0, "Red", "Yes"),
    MaxEntFeature::new(0, "Red", "No"),
    MaxEntFeature::new(0, "Yellow", "Yes"),
    MaxEntFeature::new(0, "Yellow", "No"),
    MaxEntFeature::new(1, "Sports", "Yes"),
    MaxEntFeature::new(1, "Sports", "No"),
    MaxEntFeature::new(1, "SUV", "Yes"),
    MaxEntFeature::new(1, "SUV", "No"),
    MaxEntFeature::new(2, "Domestic", "Yes"),
    MaxEntFeature::new(2, "Domestic", "No"),
    MaxEntFeature::new(2, "Imported", "Yes"),
    MaxEntFeature::new(2, "Imported", "No"),
  ]
  (xs, ys, features)
}
