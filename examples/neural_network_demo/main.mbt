///|
/// NeuralNetwork module demo - Feedforward neural network with backpropagation.
///
/// Demonstrates network construction, forward pass, backpropagation training,
/// and prediction on XOR and iris-like classification problems.
///
/// Reference: Biopython Bio.NeuralNetwork.BackPropagation module.
/// Run: moon run examples/neural_network_demo
fn main {
  println("=== NeuralNetwork Module Demo ===")
  println("Feedforward neural network with backpropagation training")
  println("Reference: Biopython Bio.NeuralNetwork")
  println("")

  // ------------------------------------------------------------------
  // 1. XOR problem
  // ------------------------------------------------------------------
  println("=== 1. XOR Problem ===")
  let (xor_inputs, xor_targets) = @bio.nn_sample_xor_data()
  println("Training data:")
  for i in 0..<xor_inputs.length() {
    println(
      "  [" +
        xor_inputs[i][0].to_string() +
        ", " +
        xor_inputs[i][1].to_string() +
        "] => " +
        xor_targets[i][0].to_string(),
    )
  }
  println("")

  let xor_net = @bio.NeuralNetwork::new(2, 4, 1)
  let xor_net = xor_net.set_learning_rate(0.5)
  let error = @bio.nn_train(xor_net, xor_inputs, xor_targets, 1000)
  println("Final MSE after 1000 epochs: " + error.to_string())
  println("")

  println("Predictions:")
  for i in 0..<xor_inputs.length() {
    let output = @bio.nn_forward(xor_net, xor_inputs[i])
    let pred = if output[0] > 0.5 { 1 } else { 0 }
    let target = xor_targets[i][0].to_int()
    let match_str = if pred == target { "✓" } else { "✗" }
    println(
      "  " +
        match_str +
        " [" +
        xor_inputs[i][0].to_string() +
        ", " +
        xor_inputs[i][1].to_string() +
        "] => output=" +
        output[0].to_string() +
        " pred=" +
        pred.to_string() +
        " target=" +
        target.to_string(),
    )
  }
  println("")

  // ------------------------------------------------------------------
  // 2. Iris-like classification
  // ------------------------------------------------------------------
  println("=== 2. Iris-like Classification ===")
  let (iris_inputs, iris_targets) = @bio.nn_sample_iris_data()
  println("Training data: " + iris_inputs.length().to_string() + " samples, 3 classes")
  let class_names = ["setosa", "versicolor", "virginica"]
  for i in 0..<class_names.length() {
    println("  Class " + i.to_string() + ": " + class_names[i])
  }
  println("")

  let iris_net = @bio.NeuralNetwork::new(2, 6, 3)
  let iris_net = iris_net.set_learning_rate(0.3)
  let iris_error = @bio.nn_train(iris_net, iris_inputs, iris_targets, 1000)
  println("Final MSE after 1000 epochs: " + iris_error.to_string())
  println("")

  println("Predictions:")
  let mut correct = 0
  for i in 0..<iris_inputs.length() {
    let pred = @bio.nn_predict(iris_net, iris_inputs[i])
    let target_class = if iris_targets[i][0] == 1.0 {
      0
    } else if iris_targets[i][1] == 1.0 {
      1
    } else {
      2
    }
    let match_str = if pred == target_class {
      correct = correct + 1
      "✓"
    } else {
      "✗"
    }
    println(
      "  " +
        match_str +
        " [" +
        iris_inputs[i][0].to_string() +
        ", " +
        iris_inputs[i][1].to_string() +
        "] => pred=" +
        class_names[pred] +
        " actual=" +
        class_names[target_class],
    )
  }
  println("Accuracy: " + correct.to_string() + "/" + iris_inputs.length().to_string())
  println("")

  // ------------------------------------------------------------------
  // 3. Network structure inspection
  // ------------------------------------------------------------------
  println("=== 3. Network Structure ===")
  let hidden = iris_net.hidden_layer()
  let output = iris_net.output_layer()
  println("Input layer: 2 units")
  println("Hidden layer: " + hidden.n_out().to_string() + " units, " + hidden.n_in().to_string() + " inputs")
  println("Output layer: " + output.n_out().to_string() + " units, " + output.n_in().to_string() + " inputs")
  println("Total weights: " + (hidden.n_out() * hidden.n_in() + output.n_out() * output.n_in()).to_string())
  println("")

  println("=== Demo complete ===")
}
