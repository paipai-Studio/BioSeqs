///|
/// NN building blocks: FCLayers, Encoder, DecoderSCVI, LinearDecoderSCVI,
/// MultiEncoder, MultiDecoder, Embedding, one_hot_encode.
pub struct FCLayers {
  n_input : Int32
  n_output : Int32
  n_layers : Int32
  use_batch_norm : Bool
  use_layer_norm : Bool
  dropout_rate : Double
  weights : Array[Array[Array[Double]]]
  biases : Array[Array[Double]]
}

///|
pub fn fclayers_init(
  n_input : Int32,
  n_output : Int32,
  n_layers : Int32,
  use_batch_norm : Bool,
  dropout_rate : Double,
) -> FCLayers {
  let weights = Array(capacity=n_layers)
  let biases = Array(capacity=n_layers)
  for layer in 0..<n_layers {
    let in_dim = if layer == 0 { n_input } else { n_output }
    let w = Array(capacity=n_output)
    for _i in 0..<n_output {
      let row = Array(capacity=in_dim)
      for _j in 0..<in_dim {
        row.push(0.01)
      }
      w.push(row)
    }
    weights.push(w)
    let b = Array(capacity=n_output)
    for _i in 0..<n_output {
      b.push(0.0)
    }
    biases.push(b)
  }
  {
    n_input,
    n_output,
    n_layers,
    use_batch_norm,
    use_layer_norm: false,
    dropout_rate,
    weights,
    biases,
  }
}

///|
pub fn fclayers_forward(layers : FCLayers, x : Array[Double]) -> Array[Double] {
  let mut current = x
  for layer in 0..<layers.n_layers {
    let next = Array(capacity=layers.n_output)
    for i in 0..<layers.n_output {
      let mut sum = layers.biases[layer][i]
      for j in 0..<current.length() {
        if j < layers.weights[layer][i].length() {
          sum = sum + current[j] * layers.weights[layer][i][j]
        }
      }
      if layers.use_batch_norm {
        sum = sum
      }
      if layers.dropout_rate > 0.0 && sum > 0.0 {
        sum = sum * (1.0 - layers.dropout_rate)
      }
      next.push(sum)
    }
    current = next
  }
  current
}

///|
pub struct EncoderNN {
  n_input : Int32
  n_latent : Int32
  fc : FCLayers
  mean_weights : Array[Double]
  logvar_weights : Array[Double]
}

///|
pub fn encoder_nn_init(
  n_input : Int32,
  n_latent : Int32,
  n_hidden : Int32,
  n_layers : Int32,
) -> EncoderNN {
  let fc = fclayers_init(n_input, n_hidden, n_layers, true, 0.1)
  let mw = Array(capacity=n_latent)
  let lvw = Array(capacity=n_latent)
  for _i in 0..<n_latent {
    mw.push(0.01)
    lvw.push(0.01)
  }
  { n_input, n_latent, fc, mean_weights: mw, logvar_weights: lvw, }
}

///|
pub fn encoder_nn_forward(
  encoder : EncoderNN,
  x : Array[Double],
) -> (Array[Double], Array[Double]) {
  let h = fclayers_forward(encoder.fc, x)
  let mean = Array(capacity=encoder.n_latent)
  let logvar = Array(capacity=encoder.n_latent)
  for i in 0..<encoder.n_latent {
    let mut m_sum = 0.0
    let mut v_sum = 0.0
    for j in 0..<h.length() {
      m_sum = m_sum +
        h[j] * encoder.mean_weights[i % encoder.mean_weights.length()]
      v_sum = v_sum +
        h[j] * encoder.logvar_weights[i % encoder.logvar_weights.length()]
    }
    mean.push(m_sum)
    logvar.push(v_sum)
  }
  (mean, logvar)
}

///|
pub struct DecoderSCVI {
  n_latent : Int32
  n_output : Int32
  fc : FCLayers
  px_scale : Array[Double]
  px_dropout : Array[Double]
}

///|
pub fn decoder_scvi_init(
  n_latent : Int32,
  n_output : Int32,
  n_hidden : Int32,
  n_layers : Int32,
) -> DecoderSCVI {
  let fc = fclayers_init(n_latent, n_hidden, n_layers, true, 0.0)
  let scale = Array(capacity=n_output)
  let dropout = Array(capacity=n_output)
  for _i in 0..<n_output {
    scale.push(0.01)
    dropout.push(0.0)
  }
  { n_latent, n_output, fc, px_scale: scale, px_dropout: dropout, }
}

///|
pub fn decoder_scvi_forward(
  decoder : DecoderSCVI,
  z : Array[Double],
) -> (Array[Double], Array[Double]) {
  let h = fclayers_forward(decoder.fc, z)
  let rate = Array(capacity=decoder.n_output)
  let dropout = Array(capacity=decoder.n_output)
  for i in 0..<decoder.n_output {
    let mut sum = 0.0
    for j in 0..<h.length() {
      sum = sum + h[j] * decoder.px_scale[i % decoder.px_scale.length()]
    }
    rate.push(@math.exp(sum))
    dropout.push(
      1.0 /
      (1.0 + @math.exp(-decoder.px_dropout[i % decoder.px_dropout.length()])),
    )
  }
  (rate, dropout)
}

///|
pub struct LinearDecoderSCVI {
  n_latent : Int32
  n_output : Int32
  weights : Array[Array[Double]]
  bias : Array[Double]
}

///|
pub fn linear_decoder_scvi_init(
  n_latent : Int32,
  n_output : Int32,
) -> LinearDecoderSCVI {
  let w = Array(capacity=n_output)
  for _i in 0..<n_output {
    let row = Array(capacity=n_latent)
    for _j in 0..<n_latent {
      row.push(0.01)
    }
    w.push(row)
  }
  let b = Array(capacity=n_output)
  for _i in 0..<n_output {
    b.push(0.0)
  }
  { n_latent, n_output, weights: w, bias: b, }
}

///|
pub fn linear_decoder_scvi_forward(
  decoder : LinearDecoderSCVI,
  z : Array[Double],
) -> Array[Double] {
  let result = Array(capacity=decoder.n_output)
  for i in 0..<decoder.n_output {
    let mut sum = decoder.bias[i]
    for j in 0..<z.length() {
      if j < decoder.weights[i].length() {
        sum = sum + z[j] * decoder.weights[i][j]
      }
    }
    result.push(@math.exp(sum))
  }
  result
}

///|
pub fn linear_decoder_scvi_get_loadings(
  decoder : LinearDecoderSCVI,
) -> Array[Array[Double]] {
  decoder.weights
}

///|
pub struct MultiEncoderNN {
  n_input_rna : Int32
  n_input_atac : Int32
  n_latent : Int32
  rna_encoder : EncoderNN
  atac_encoder : EncoderNN
  joint_weights : Array[Array[Double]]
}

///|
pub fn multi_encoder_init(
  n_input_rna : Int32,
  n_input_atac : Int32,
  n_latent : Int32,
  n_hidden : Int32,
) -> MultiEncoderNN {
  let rna_enc = encoder_nn_init(n_input_rna, n_latent, n_hidden, 2)
  let atac_enc = encoder_nn_init(n_input_atac, n_latent, n_hidden, 2)
  let jw = Array(capacity=n_latent)
  for _i in 0..<n_latent {
    let row = Array(capacity=n_latent * 2)
    for _j in 0..<(n_latent * 2) {
      row.push(0.01)
    }
    jw.push(row)
  }
  {
    n_input_rna,
    n_input_atac,
    n_latent,
    rna_encoder: rna_enc,
    atac_encoder: atac_enc,
    joint_weights: jw,
  }
}

///|
pub fn multi_encoder_forward(
  encoder : MultiEncoderNN,
  rna : Array[Double],
  atac : Array[Double],
) -> Array[Double] {
  let (rna_mean, _) = encoder_nn_forward(encoder.rna_encoder, rna)
  let (atac_mean, _) = encoder_nn_forward(encoder.atac_encoder, atac)
  let combined = Array(capacity=rna_mean.length() + atac_mean.length())
  for i in 0..<rna_mean.length() {
    combined.push(rna_mean[i])
  }
  for i in 0..<atac_mean.length() {
    combined.push(atac_mean[i])
  }
  let z = Array(capacity=encoder.n_latent)
  for i in 0..<encoder.n_latent {
    let mut sum = 0.0
    for j in 0..<combined.length() {
      if j < encoder.joint_weights[i].length() {
        sum = sum + combined[j] * encoder.joint_weights[i][j]
      }
    }
    z.push(sum)
  }
  z
}

///|
pub struct EmbeddingLayer {
  n_categories : Int32
  n_dim : Int32
  weights : Array[Array[Double]]
}

///|
pub fn embedding_layer_init(
  n_categories : Int32,
  n_dim : Int32,
) -> EmbeddingLayer {
  let w = Array(capacity=n_categories)
  for _i in 0..<n_categories {
    let row = Array(capacity=n_dim)
    for _j in 0..<n_dim {
      row.push(0.01)
    }
    w.push(row)
  }
  { n_categories, n_dim, weights: w, }
}

///|
pub fn embedding_layer_forward(
  emb : EmbeddingLayer,
  indices : Array[Int32],
) -> Array[Array[Double]] {
  let result = Array(capacity=indices.length())
  for i in 0..<indices.length() {
    let idx = indices[i]
    if idx >= 0 && idx < emb.weights.length() {
      result.push(emb.weights[idx])
    } else {
      let zeros = Array(capacity=emb.n_dim)
      for _j in 0..<emb.n_dim {
        zeros.push(0.0)
      }
      result.push(zeros)
    }
  }
  result
}

///|
pub fn one_hot_encode(
  indices : Array[Int32],
  n_classes : Int32,
) -> Array[Array[Double]] {
  let result = Array(capacity=indices.length())
  for i in 0..<indices.length() {
    let row = Array(capacity=n_classes)
    for j in 0..<n_classes {
      row.push(if j == indices[i] { 1.0 } else { 0.0 })
    }
    result.push(row)
  }
  result
}

///|
pub fn _scvi_nn_dummy() -> Unit {
  let _ = FCLayers::{
    n_input: 0,
    n_output: 0,
    n_layers: 0,
    use_batch_norm: false,
    use_layer_norm: false,
    dropout_rate: 0.0,
    weights: [],
    biases: [],
  }
  let _ = EncoderNN::{
    n_input: 0,
    n_latent: 0,
    fc: fclayers_init(1, 1, 1, false, 0.0),
    mean_weights: [],
    logvar_weights: [],
  }
  let _ = DecoderSCVI::{
    n_latent: 0,
    n_output: 0,
    fc: fclayers_init(1, 1, 1, false, 0.0),
    px_scale: [],
    px_dropout: [],
  }
  let _ = LinearDecoderSCVI::{
    n_latent: 0,
    n_output: 0,
    weights: [],
    bias: [],
  }
  let _ = MultiEncoderNN::{
    n_input_rna: 0,
    n_input_atac: 0,
    n_latent: 0,
    rna_encoder: encoder_nn_init(1, 1, 1, 1),
    atac_encoder: encoder_nn_init(1, 1, 1, 1),
    joint_weights: [],
  }
  let _ = EmbeddingLayer::{ n_categories: 0, n_dim: 0, weights: [], }
}
