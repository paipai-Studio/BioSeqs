///|
/// MaxEntropy module demo - Maximum entropy classifier with IIS training.
///
/// Demonstrates feature construction, IIS training, and classification
/// on the car preference dataset.
///
/// Reference: Chang J (2001) Biopython Bio.MaxEntropy module.
/// Run: moon run examples/max_entropy_demo
fn main {
  println("=== MaxEntropy Module Demo ===")
  println("Maximum entropy classifier with Improved Iterative Scaling (IIS)")
  println("Reference: Biopython Bio.MaxEntropy")
  println("")

  // ------------------------------------------------------------------
  // 1. Load sample data (car preference classification)
  // ------------------------------------------------------------------
  println("=== 1. Load Sample Data ===")
  let (xs, ys, features) = @bio.max_entropy_sample_data()
  println("Training samples: " + xs.length().to_string())
  println("Features: " + features.length().to_string())
  println("Attributes: [color, type, origin]")
  println("")
  println("Training data:")
  for i in 0..<xs.length() {
    let color = xs[i][0]
    let typ = xs[i][1]
    let origin = xs[i][2]
    println("  " + color + " " + typ + " " + origin + " => " + ys[i])
  }
  println("")

  // ------------------------------------------------------------------
  // 2. Train the model
  // ------------------------------------------------------------------
  println("=== 2. Train Model (IIS) ===")
  let model = @bio.max_entropy_train(xs, ys, features)
  println("Classes: " + model.classes().to_string())
  println("Feature weights (alphas):")
  let alphas = model.alphas()
  let feats = model.features()
  for i in 0..<alphas.length() {
    let f = feats[i]
    println(
      "  alpha=" +
        alphas[i].to_string() +
        "  attr[" +
        f.attribute_index().to_string() +
        "]=" +
        f.attribute_value() +
        " class=" +
        f.class_label(),
    )
  }
  println("")

  // ------------------------------------------------------------------
  // 3. Classify training data
  // ------------------------------------------------------------------
  println("=== 3. Classify Training Data ===")
  let mut correct = 0
  for i in 0..<xs.length() {
    let pred = @bio.max_entropy_classify(model, xs[i])
    let match_str = if pred == ys[i] {
      correct = correct + 1
      "✓"
    } else {
      "✗"
    }
    println(
      "  " +
        match_str +
        " " +
        xs[i][0] +
        " " +
        xs[i][1] +
        " " +
        xs[i][2] +
        " => pred=" +
        pred +
        " actual=" +
        ys[i],
    )
  }
  println("Accuracy: " + correct.to_string() + "/" + xs.length().to_string())
  println("")

  // ------------------------------------------------------------------
  // 4. Predict on new observations
  // ------------------------------------------------------------------
  println("=== 4. Predict on New Observations ===")
  let test_cases = [
    ["Red", "Sports", "Domestic"],
    ["Yellow", "SUV", "Imported"],
    ["Red", "SUV", "Domestic"],
    ["Yellow", "Sports", "Imported"],
  ]
  for tc in test_cases {
    let pred = @bio.max_entropy_classify(model, tc)
    let scores = @bio.max_entropy_calculate(model, tc)
    println(
      "  " +
        tc[0] +
        " " +
        tc[1] +
        " " +
        tc[2] +
        " => " +
        pred +
        " (scores: " +
        scores.to_string() +
        ")",
    )
  }
  println("")

  println("=== Demo complete ===")
}
