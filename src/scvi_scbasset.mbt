///|
/// SCBASSET: Chromatin accessibility prediction model using CNN on DNA sequences.
pub struct ScbassetConfig {
  n_filters : Int32
  n_conv_layers : Int32
  kernel_size : Int32
  n_dense_layers : Int32
  n_dense_hidden : Int32
  dropout_rate : Double
  seq_length : Int32
}

///|
pub fn ScbassetConfig::default() -> ScbassetConfig {
  {
    n_filters: 64,
    n_conv_layers: 2,
    kernel_size: 7,
    n_dense_layers: 1,
    n_dense_hidden: 128,
    dropout_rate: 0.3,
    seq_length: 1000,
  }
}

///|
pub struct ScbassetModel {
  config : ScbassetConfig
  conv_weights : Array[Array[Array[Double]]]
  dense_weights : Array[Array[Double]]
  output_bias : Array[Double]
  mut is_trained : Bool
}

///|
pub fn scbasset_init(config : ScbassetConfig) -> ScbassetModel {
  let conv = Array(capacity=config.n_conv_layers)
  for layer in 0..<config.n_conv_layers {
    let _n_in = if layer == 0 { 4 } else { config.n_filters }
    let filters = Array(capacity=config.n_filters)
    for _f in 0..<config.n_filters {
      let kernel = Array(capacity=config.kernel_size)
      for _k in 0..<config.kernel_size {
        kernel.push(0.01)
      }
      filters.push(kernel)
    }
    conv.push(filters)
  }
  let dense = Array(capacity=config.n_dense_layers)
  for _layer in 0..<config.n_dense_layers {
    let n_in = config.n_dense_hidden
    let row = Array(capacity=n_in)
    for _j in 0..<n_in {
      row.push(0.01)
    }
    dense.push(row)
  }
  {
    config,
    conv_weights: conv,
    dense_weights: dense,
    output_bias: [],
    is_trained: false,
  }
}

///|
pub fn scbasset_one_hot_encode(seq : String) -> Array[Array[Double]] {
  let encoded = Array(capacity=seq.length())
  for i in 0..<seq.length() {
    let c = seq[i]
    let row = [0.0, 0.0, 0.0, 0.0]
    let idx = match c.to_int() {
      65 | 97 => 0 // A/a
      67 | 99 => 1 // C/c
      71 | 103 => 2 // G/g
      84 | 116 => 3 // T/t
      _ => -1
    }
    if idx >= 0 {
      row[idx] = 1.0
    }
    encoded.push(row)
  }
  encoded
}

///|
pub fn scbasset_conv1d(
  input : Array[Array[Double]],
  weights : Array[Array[Double]],
  kernel_size : Int32,
) -> Array[Array[Double]] {
  let n_filters = weights.length()
  let seq_len = input.length()
  let out_len = seq_len - kernel_size + 1
  let output = Array(capacity=out_len)
  for pos in 0..<out_len {
    let row = Array(capacity=n_filters)
    for f in 0..<n_filters {
      let mut sum = 0.0
      for k in 0..<kernel_size {
        if k < weights[f].length() && pos + k < input.length() {
          for c in 0..<input[pos + k].length() {
            sum = sum + input[pos + k][c] * weights[f][k]
          }
        }
      }
      row.push(sum)
    }
    output.push(row)
  }
  output
}

///|
pub fn scbasset_predict_accessibility(
  model : ScbassetModel,
  seq : String,
) -> Double {
  let encoded = scbasset_one_hot_encode(seq)
  let mut x = encoded
  for layer in 0..<model.config.n_conv_layers {
    x = scbasset_conv1d(x, model.conv_weights[layer], model.config.kernel_size)
    for pos in 0..<x.length() {
      for f in 0..<x[pos].length() {
        if x[pos][f] < 0.0 {
          x[pos][f] = 0.0
        }
      }
    }
  }
  let mut sum = 0.0
  for pos in 0..<x.length() {
    for f in 0..<x[pos].length() {
      sum = sum + x[pos][f]
    }
  }
  sum = sum / x.length().to_double()
  for layer in 0..<model.dense_weights.length() {
    let mut dense_sum = 0.0
    for j in 0..<model.dense_weights[layer].length() {
      dense_sum = dense_sum + sum * model.dense_weights[layer][j]
    }
    sum = dense_sum
    if sum < 0.0 {
      sum = 0.0
    }
  }
  1.0 / (1.0 + @math.exp(-sum))
}

///|
pub fn scbasset_train(
  model : ScbassetModel,
  seqs : Array[String],
  labels : Array[Double],
  n_epochs : Int32,
) -> ScbassetModel {
  let m = model
  for _epoch in 0..<n_epochs {
    for i in 0..<seqs.length() {
      let pred = scbasset_predict_accessibility(m, seqs[i])
      let grad = pred - labels[i]
      for layer in 0..<m.conv_weights.length() {
        for f in 0..<m.conv_weights[layer].length() {
          for k in 0..<m.conv_weights[layer][f].length() {
            m.conv_weights[layer][f][k] = m.conv_weights[layer][f][k] -
              0.001 * grad * 0.01
          }
        }
      }
    }
  }
  m.is_trained = true
  m
}

///|
pub fn scbasset_interpret_sequence(
  _mdl : ScbassetModel,
  seq : String,
) -> Array[Double] {
  let encoded = scbasset_one_hot_encode(seq)
  let scores = Array(capacity=seq.length())
  for i in 0..<seq.length() {
    let mut score = 0.0
    for c in 0..<4 {
      score = score + encoded[i][c] * encoded[i][c]
    }
    scores.push(score)
  }
  scores
}
