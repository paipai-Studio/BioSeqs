///|
/// Bio.NeuralNetwork - Feedforward neural network with backpropagation
/// training, inspired by Biopython Bio.NeuralNetwork.BackPropagation.
/// Reference: Biopython NeuralNetwork module; Rumelhart et al. (1986)
/// Nature 323:533-536.

// ===========================================================================
// Layer struct
// ===========================================================================

///|
/// A layer in the neural network: weights[n_out][n_in], biases[n_out],
/// and the last activation values for backprop.
pub struct NnLayer {
  n_in : Int
  n_out : Int
  weights : Array[Array[Double]] // [n_out][n_in]
  biases : Array[Double] // [n_out]
  last_input : Array[Double] // cached input for backprop
  last_output : Array[Double] // cached output (activation) for backprop
  weight_grads : Array[Array[Double]] // gradient accumulator
  bias_grads : Array[Double] // gradient accumulator
}

///|
/// Construct a new layer with given input and output sizes. Weights are
/// initialized with small random values; biases to zero.
pub fn NnLayer::new(n_in : Int, n_out : Int) -> NnLayer {
  let weights : Array[Array[Double]] = Array::new()
  let weight_grads : Array[Array[Double]] = Array::new()
  for _ in 0..<n_out {
    let row = Array::make(n_in, 0.0)
    let grad_row = Array::make(n_in, 0.0)
    // Small random initialization in [-0.5, 0.5].
    for j in 0..<n_in {
      row[j] = nn_random_weight()
    }
    weights.push(row)
    weight_grads.push(grad_row)
  }
  let biases = Array::make(n_out, 0.0)
  let bias_grads = Array::make(n_out, 0.0)
  NnLayer::{
    n_in,
    n_out,
    weights,
    biases,
    last_input: Array::make(n_in, 0.0),
    last_output: Array::make(n_out, 0.0),
    weight_grads,
    bias_grads,
  }
}

///|
/// Accessor for the number of input units.
pub fn NnLayer::n_in(self : NnLayer) -> Int {
  self.n_in
}

///|
/// Accessor for the number of output units.
pub fn NnLayer::n_out(self : NnLayer) -> Int {
  self.n_out
}

///|
/// Accessor for the weight matrix [n_out][n_in].
pub fn NnLayer::weights(self : NnLayer) -> Array[Array[Double]] {
  self.weights
}

///|
/// Accessor for the bias vector [n_out].
pub fn NnLayer::biases(self : NnLayer) -> Array[Double] {
  self.biases
}

// ===========================================================================
// Network struct
// ===========================================================================

///|
/// A feedforward neural network with one hidden layer.
pub struct NeuralNetwork {
  input_layer : NnLayer
  hidden_layer : NnLayer
  output_layer : NnLayer
  learning_rate : Double
  momentum : Double
  // Previous weight deltas for momentum.
  prev_hidden_deltas : Array[Array[Double]]
  prev_output_deltas : Array[Array[Double]]
  prev_hidden_bias_deltas : Array[Double]
  prev_output_bias_deltas : Array[Double]
}

///|
/// Construct a new neural network with given architecture.
/// n_input: number of input features.
/// n_hidden: number of hidden units.
/// n_output: number of output classes.
pub fn NeuralNetwork::new(
  n_input : Int,
  n_hidden : Int,
  n_output : Int
) -> NeuralNetwork {
  let hidden = NnLayer::new(n_input, n_hidden)
  let output = NnLayer::new(n_hidden, n_output)
  let input = NnLayer::new(n_input, n_input) // pass-through (identity)
  let prev_hidden_deltas : Array[Array[Double]] = Array::new()
  for _ in 0..<n_hidden {
    prev_hidden_deltas.push(Array::make(n_input, 0.0))
  }
  let prev_output_deltas : Array[Array[Double]] = Array::new()
  for _ in 0..<n_output {
    prev_output_deltas.push(Array::make(n_hidden, 0.0))
  }
  NeuralNetwork::{
    input_layer: input,
    hidden_layer: hidden,
    output_layer: output,
    learning_rate: 0.1,
    momentum: 0.0,
    prev_hidden_deltas,
    prev_output_deltas,
    prev_hidden_bias_deltas: Array::make(n_hidden, 0.0),
    prev_output_bias_deltas: Array::make(n_output, 0.0),
  }
}

///|
/// Set the learning rate.
pub fn NeuralNetwork::set_learning_rate(
  self : NeuralNetwork,
  lr : Double
) -> NeuralNetwork {
  { ..self, learning_rate: lr }
}

///|
/// Set the momentum coefficient.
pub fn NeuralNetwork::set_momentum(
  self : NeuralNetwork,
  m : Double
) -> NeuralNetwork {
  { ..self, momentum: m }
}

///|
/// Accessor for the hidden layer.
pub fn NeuralNetwork::hidden_layer(self : NeuralNetwork) -> NnLayer {
  self.hidden_layer
}

///|
/// Accessor for the output layer.
pub fn NeuralNetwork::output_layer(self : NeuralNetwork) -> NnLayer {
  self.output_layer
}

///|
/// Accessor for the learning rate.
pub fn NeuralNetwork::learning_rate(self : NeuralNetwork) -> Double {
  self.learning_rate
}

///|
/// Accessor for the momentum coefficient.
pub fn NeuralNetwork::momentum(self : NeuralNetwork) -> Double {
  self.momentum
}

// ===========================================================================
// Internal helpers
// ===========================================================================

///|
/// Simple pseudo-random number generator state (LCG).
/// Uses a fixed seed for reproducibility.
let nn_rng_counter : Ref[Int] = Ref(42)

///|
/// Generate a pseudo-random weight in [-0.5, 0.5].
/// Uses an LCG with modulo to stay within Int range.
fn nn_random_weight() -> Double {
  // LCG: state = (state * 1103515245 + 12345) % 2147483648
  // Use bit masking to emulate modulo 2^31.
  nn_rng_counter.val = (nn_rng_counter.val * 1103515245 + 12345) & 0x7fffffff
  let val = (nn_rng_counter.val % 10000).to_double() / 10000.0
  val - 0.5
}

///|
/// Sigmoid activation function.
fn nn_sigmoid(x : Double) -> Double {
  if x >= 40.0 {
    return 1.0
  }
  if x <= -40.0 {
    return 0.0000000001
  }
  1.0 / (1.0 + @math.exp(-x))
}

///|
/// Derivative of sigmoid: sigmoid(x) * (1 - sigmoid(x)).
/// Given the activation output `a` (where a = sigmoid(x)), returns a*(1-a).
fn nn_sigmoid_deriv(a : Double) -> Double {
  a * (1.0 - a)
}

///|
/// Forward pass through a single layer. Stores input/output for backprop.
fn nn_layer_forward(layer : NnLayer, input : Array[Double]) -> Array[Double] {
  let output = Array::make(layer.n_out, 0.0)
  for i in 0..<layer.n_out {
    let mut s = layer.biases[i]
    for j in 0..<layer.n_in {
      s = s + layer.weights[i][j] * input[j]
    }
    output[i] = nn_sigmoid(s)
  }
  // Cache for backprop.
  for j in 0..<layer.n_in {
    layer.last_input[j] = input[j]
  }
  for i in 0..<layer.n_out {
    layer.last_output[i] = output[i]
  }
  output
}

// ===========================================================================
// Forward / backward passes
// ===========================================================================

///|
/// Forward pass through the entire network. Returns the output activations.
pub fn nn_forward(
  network : NeuralNetwork,
  input : Array[Double]
) -> Array[Double] {
  let hidden_out = nn_layer_forward(network.hidden_layer, input)
  nn_layer_forward(network.output_layer, hidden_out)
}

///|
/// Predict the class index (argmax of output).
pub fn nn_predict(network : NeuralNetwork, input : Array[Double]) -> Int {
  let output = nn_forward(network, input)
  let mut best = 0
  let mut best_val = output[0]
  for i in 1..<output.length() {
    if output[i] > best_val {
      best_val = output[i]
      best = i
    }
  }
  best
}

///|
/// Train the network on a single example using backpropagation.
/// `input` is the feature vector, `target` is the target output vector
/// (same length as n_output). Returns the squared error for this example.
pub fn nn_train_example(
  network : NeuralNetwork,
  input : Array[Double],
  target : Array[Double]
) -> Double {
  // Forward pass.
  let hidden_out = nn_layer_forward(network.hidden_layer, input)
  let output = nn_layer_forward(network.output_layer, hidden_out)
  // Compute output layer deltas: delta = (output - target) * sigmoid'(output)
  let output_deltas = Array::make(network.output_layer.n_out, 0.0)
  let mut error_sum = 0.0
  for i in 0..<network.output_layer.n_out {
    let err = output[i] - target[i]
    output_deltas[i] = err * nn_sigmoid_deriv(output[i])
    error_sum = error_sum + err * err
  }
  // Compute hidden layer deltas: delta = (sum_j w_ji * delta_j) * sigmoid'(hidden)
  let hidden_deltas = Array::make(network.hidden_layer.n_out, 0.0)
  for i in 0..<network.hidden_layer.n_out {
    let mut s = 0.0
    for j in 0..<network.output_layer.n_out {
      s = s + network.output_layer.weights[j][i] * output_deltas[j]
    }
    hidden_deltas[i] = s * nn_sigmoid_deriv(hidden_out[i])
  }
  // Update output layer weights and biases.
  let lr = network.learning_rate
  let mom = network.momentum
  for i in 0..<network.output_layer.n_out {
    for j in 0..<network.output_layer.n_in {
      let grad = output_deltas[i] * network.output_layer.last_input[j]
      let delta = lr * grad + mom * network.prev_output_deltas[i][j]
      network.output_layer.weights[i][j] = network.output_layer.weights[i][j] - delta
      network.prev_output_deltas[i][j] = delta
    }
    let bias_delta = lr * output_deltas[i] + mom * network.prev_output_bias_deltas[i]
    network.output_layer.biases[i] = network.output_layer.biases[i] - bias_delta
    network.prev_output_bias_deltas[i] = bias_delta
  }
  // Update hidden layer weights and biases.
  for i in 0..<network.hidden_layer.n_out {
    for j in 0..<network.hidden_layer.n_in {
      let grad = hidden_deltas[i] * network.hidden_layer.last_input[j]
      let delta = lr * grad + mom * network.prev_hidden_deltas[i][j]
      network.hidden_layer.weights[i][j] = network.hidden_layer.weights[i][j] - delta
      network.prev_hidden_deltas[i][j] = delta
    }
    let bias_delta = lr * hidden_deltas[i] + mom * network.prev_hidden_bias_deltas[i]
    network.hidden_layer.biases[i] = network.hidden_layer.biases[i] - bias_delta
    network.prev_hidden_bias_deltas[i] = bias_delta
  }
  error_sum / output.length().to_double()
}

///|
/// Train the network on a dataset for a number of epochs.
/// `inputs` is an array of feature vectors, `targets` is an array of target
/// output vectors. Returns the final mean squared error.
pub fn nn_train(
  network : NeuralNetwork,
  inputs : Array[Array[Double]],
  targets : Array[Array[Double]],
  n_epochs : Int
) -> Double {
  let n = inputs.length()
  let mut final_error = 0.0
  for _ in 0..<n_epochs {
    let mut epoch_error = 0.0
    for i in 0..<n {
      epoch_error = epoch_error + nn_train_example(network, inputs[i], targets[i])
    }
    final_error = epoch_error / n.to_double()
  }
  final_error
}

// ===========================================================================
// Sample data: XOR problem
// ===========================================================================

///|
/// Sample training data for the XOR problem.
/// Returns (inputs, targets) where inputs are 2D and targets are 1D.
pub fn nn_sample_xor_data() -> (Array[Array[Double]], Array[Array[Double]]) {
  let inputs = [
    [0.0, 0.0],
    [0.0, 1.0],
    [1.0, 0.0],
    [1.0, 1.0],
  ]
  let targets = [[0.0], [1.0], [1.0], [0.0]]
  (inputs, targets)
}

///|
/// Sample training data for the iris-like classification.
/// 3 classes, 2 features. Returns (inputs, targets) where targets are
/// one-hot encoded (3 outputs).
pub fn nn_sample_iris_data() -> (Array[Array[Double]], Array[Array[Double]]) {
  let inputs = [
    [5.1, 3.5],
    [4.9, 3.0],
    [4.7, 3.2],
    [4.6, 3.1],
    [5.0, 3.6],
    [7.0, 3.2],
    [6.4, 3.2],
    [6.9, 3.1],
    [5.5, 2.3],
    [6.5, 2.8],
    [6.3, 3.3],
    [5.8, 2.7],
    [6.1, 2.9],
    [5.1, 2.5],
    [6.3, 2.5],
  ]
  // Class 0: setosa (rows 0-4), Class 1: versicolor (rows 5-9), Class 2: virginica (rows 10-14)
  let targets = [
    [1.0, 0.0, 0.0],
    [1.0, 0.0, 0.0],
    [1.0, 0.0, 0.0],
    [1.0, 0.0, 0.0],
    [1.0, 0.0, 0.0],
    [0.0, 1.0, 0.0],
    [0.0, 1.0, 0.0],
    [0.0, 1.0, 0.0],
    [0.0, 1.0, 0.0],
    [0.0, 1.0, 0.0],
    [0.0, 0.0, 1.0],
    [0.0, 0.0, 1.0],
    [0.0, 0.0, 1.0],
    [0.0, 0.0, 1.0],
    [0.0, 0.0, 1.0],
  ]
  (inputs, targets)
}
