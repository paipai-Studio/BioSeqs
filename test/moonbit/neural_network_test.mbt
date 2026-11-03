///|
/// Tests for the NeuralNetwork module (Biopython Bio.NeuralNetwork).
///
/// Covers layer construction, network construction, forward pass, training,
/// prediction, sample data validation, and edge cases.

// ===========================================================================
// 1. Layer construction
// ===========================================================================

///|
test "nn_layer_creation" {
  let layer = @src.NnLayer::new(3, 2)
  assert_eq(layer.n_in(), 3)
  assert_eq(layer.n_out(), 2)
}

///|
test "nn_layer_weights_shape" {
  let layer = @src.NnLayer::new(4, 3)
  let w = layer.weights()
  assert_eq(w.length(), 3) // n_out rows
  for row in w {
    assert_eq(row.length(), 4) // n_in cols
  }
}

///|
test "nn_layer_biases_length" {
  let layer = @src.NnLayer::new(5, 2)
  let b = layer.biases()
  assert_eq(b.length(), 2)
}

///|
test "nn_layer_biases_zero" {
  let layer = @src.NnLayer::new(3, 2)
  let b = layer.biases()
  for v in b {
    assert_true(v == 0.0)
  }
}

// ===========================================================================
// 2. Network construction
// ===========================================================================

///|
test "nn_network_creation" {
  let net = @src.NeuralNetwork::new(2, 4, 1)
  let h = net.hidden_layer()
  let o = net.output_layer()
  assert_eq(h.n_in(), 2)
  assert_eq(h.n_out(), 4)
  assert_eq(o.n_in(), 4)
  assert_eq(o.n_out(), 1)
}

///|
test "nn_network_set_learning_rate" {
  let net = @src.NeuralNetwork::new(2, 4, 1)
  let net2 = net.set_learning_rate(0.5)
  assert_true(net2.learning_rate() == 0.5)
}

///|
test "nn_network_set_momentum" {
  let net = @src.NeuralNetwork::new(2, 4, 1)
  let net2 = net.set_momentum(0.9)
  assert_true(net2.momentum() == 0.9)
}

// ===========================================================================
// 3. Forward pass
// ===========================================================================

///|
test "nn_forward_returns_output" {
  let net = @src.NeuralNetwork::new(2, 4, 1)
  let output = @src.nn_forward(net, [1.0, 0.0])
  assert_eq(output.length(), 1)
}

///|
test "nn_forward_multi_output" {
  let net = @src.NeuralNetwork::new(2, 4, 3)
  let output = @src.nn_forward(net, [0.5, 0.5])
  assert_eq(output.length(), 3)
}

///|
test "nn_forward_output_in_range" {
  let net = @src.NeuralNetwork::new(2, 4, 1)
  let output = @src.nn_forward(net, [1.0, 1.0])
  for v in output {
    assert_true(v >= 0.0)
    assert_true(v <= 1.0)
  }
}

// ===========================================================================
// 4. Prediction
// ===========================================================================

///|
test "nn_predict_returns_index" {
  let net = @src.NeuralNetwork::new(2, 4, 3)
  let pred = @src.nn_predict(net, [0.5, 0.5])
  assert_true(pred >= 0)
  assert_true(pred < 3)
}

///|
test "nn_predict_binary" {
  let net = @src.NeuralNetwork::new(2, 4, 1)
  let pred = @src.nn_predict(net, [1.0, 0.0])
  assert_true(pred == 0)
}

// ===========================================================================
// 5. Training
// ===========================================================================

///|
test "nn_train_returns_error" {
  let (inputs, targets) = @src.nn_sample_xor_data()
  let net = @src.NeuralNetwork::new(2, 4, 1)
  let error = @src.nn_train(net, inputs, targets, 100)
  assert_true(error >= 0.0)
}

///|
test "nn_train_xor_reduces_error" {
  let (inputs, targets) = @src.nn_sample_xor_data()
  let net = @src.NeuralNetwork::new(2, 4, 1)
  // Train for a few epochs and record error.
  let error1 = @src.nn_train(net, inputs, targets, 10)
  // Continue training on the same network.
  let error2 = @src.nn_train(net, inputs, targets, 490)
  // More training should reduce error.
  assert_true(error2 <= error1)
}

///|
test "nn_train_example_returns_error" {
  let (inputs, targets) = @src.nn_sample_xor_data()
  let net = @src.NeuralNetwork::new(2, 4, 1)
  let error = @src.nn_train_example(net, inputs[0], targets[0])
  assert_true(error >= 0.0)
}

// ===========================================================================
// 6. Sample data
// ===========================================================================

///|
test "nn_sample_xor_data" {
  let (inputs, targets) = @src.nn_sample_xor_data()
  assert_eq(inputs.length(), 4)
  assert_eq(targets.length(), 4)
  for inp in inputs {
    assert_eq(inp.length(), 2)
  }
  for tgt in targets {
    assert_eq(tgt.length(), 1)
  }
}

///|
test "nn_sample_iris_data" {
  let (inputs, targets) = @src.nn_sample_iris_data()
  assert_eq(inputs.length(), 15)
  assert_eq(targets.length(), 15)
  for inp in inputs {
    assert_eq(inp.length(), 2)
  }
  for tgt in targets {
    assert_eq(tgt.length(), 3)
  }
}

///|
test "nn_sample_iris_targets_one_hot" {
  let (_, targets) = @src.nn_sample_iris_data()
  for tgt in targets {
    let mut sum = 0.0
    for v in tgt {
      sum = sum + v
    }
    assert_true(sum == 1.0)
  }
}

// ===========================================================================
// 7. Iris classification
// ===========================================================================

///|
test "nn_iris_train_and_predict" {
  let (inputs, targets) = @src.nn_sample_iris_data()
  let net = @src.NeuralNetwork::new(2, 6, 3)
  let _ = @src.nn_train(net, inputs, targets, 500)
  // Predict on training data.
  let mut correct = 0
  for i in 0..<inputs.length() {
    let pred = @src.nn_predict(net, inputs[i])
    // Check if predicted class matches target (argmax).
    let target_class = if targets[i][0] == 1.0 {
      0
    } else if targets[i][1] == 1.0 {
      1
    } else {
      2
    }
    if pred == target_class {
      correct = correct + 1
    }
  }
  // Should get at least 8 out of 15 correct.
  assert_true(correct >= 8)
}

// ===========================================================================
// 8. Edge cases
// ===========================================================================

///|
test "nn_train_single_example" {
  let net = @src.NeuralNetwork::new(2, 3, 1)
  let error = @src.nn_train_example(net, [1.0, 0.0], [1.0])
  assert_true(error >= 0.0)
}

///|
test "nn_forward_zero_input" {
  let net = @src.NeuralNetwork::new(2, 4, 1)
  let output = @src.nn_forward(net, [0.0, 0.0])
  assert_eq(output.length(), 1)
  assert_true(output[0] >= 0.0)
  assert_true(output[0] <= 1.0)
}

///|
test "nn_predict_3_class" {
  let net = @src.NeuralNetwork::new(2, 5, 3)
  let pred = @src.nn_predict(net, [1.0, 0.5])
  assert_true(pred >= 0)
  assert_true(pred < 3)
}
