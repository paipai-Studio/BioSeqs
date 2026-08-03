///|
/// LogisticRegression module demo - Binary logistic regression classifier.
///
/// Demonstrates Newton-Raphson training, probability calculation, and
/// classification on the operon prediction dataset.
///
/// Reference: Chang J (2002) Biopython Bio.LogisticRegression module.
/// Run: moon run examples/logistic_regression_demo
fn main {
  println("=== LogisticRegression Module Demo ===")
  println("Binary logistic regression with Newton-Raphson optimization")
  println("Reference: Biopython Bio.LogisticRegression")
  println("")

  // ------------------------------------------------------------------
  // 1. Load sample data (operon prediction in Bacillus subtilis)
  // ------------------------------------------------------------------
  println("=== 1. Load Sample Data ===")
  let (xs, ys) = @bio.logistic_regression_sample_data()
  println("Training samples: " + xs.length().to_string())
  let mut n_op = 0
  let mut n_nop = 0
  for y in ys {
    if y == 1 {
      n_op = n_op + 1
    } else {
      n_nop = n_nop + 1
    }
  }
  println("  Operon pairs (OP, class=1): " + n_op.to_string())
  println("  Non-operon pairs (NOP, class=0): " + n_nop.to_string())
  println("  Features: [intergenic_distance, expression_score]")
  println("")

  // ------------------------------------------------------------------
  // 2. Train the model
  // ------------------------------------------------------------------
  println("=== 2. Train Model ===")
  let model = @bio.logistic_regression_train(xs, ys)
  let beta = model.beta()
  println("Beta weights (bias + 2 features):")
  for i in 0..<beta.length() {
    println("  beta[" + i.to_string() + "] = " + beta[i].to_string())
  }
  println("")

  // ------------------------------------------------------------------
  // 3. Classify training data
  // ------------------------------------------------------------------
  println("=== 3. Classify Training Data ===")
  let mut correct = 0
  for i in 0..<xs.length() {
    let pred = @bio.logistic_regression_classify(model, xs[i])
    let probs = @bio.logistic_regression_calculate(model, xs[i])
    let label = if ys[i] == 1 { "OP" } else { "NOP" }
    let pred_label = if pred == 1 { "OP" } else { "NOP" }
    let match_str = if pred == ys[i] {
      correct = correct + 1
      "✓"
    } else {
      "✗"
    }
    println(
      "  " +
        match_str +
        " dist=" +
        xs[i][0].to_string() +
        " expr=" +
        xs[i][1].to_string() +
        " => P(OP)=" +
        probs[1].to_string() +
        " pred=" +
        pred_label +
        " actual=" +
        label,
    )
  }
  println("Accuracy: " + correct.to_string() + "/" + xs.length().to_string())
  println("")

  // ------------------------------------------------------------------
  // 4. Predict on new observations
  // ------------------------------------------------------------------
  println("=== 4. Predict on New Observations ===")
  let test_cases = [
    [20.0, -200.0], // Short distance, high similarity -> likely OP
    [200.0, -300.0], // Long distance, low similarity -> likely NOP
    [-30.0, -180.0], // Overlapping genes -> very likely OP
    [100.0, -250.0], // Borderline case
  ]
  for tc in test_cases {
    let pred = @bio.logistic_regression_classify(model, tc)
    let probs = @bio.logistic_regression_calculate(model, tc)
    let label = if pred == 1 { "OP (same operon)" } else { "NOP (different operon)" }
    println(
      "  dist=" +
        tc[0].to_string() +
        " expr=" +
        tc[1].to_string() +
        " => " +
        label +
        " (P=" +
        probs[1].to_string() +
        ")",
    )
  }
  println("")

  println("=== Demo complete ===")
}
