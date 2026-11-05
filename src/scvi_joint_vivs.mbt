///|
/// JointEmbeddingSCVI: Multi-modal joint embedding model for data integration.
/// VIVS: Variational inference for spatial data.
pub struct JointEmbeddingConfig {
  n_latent : Int32
  n_hidden : Int32
  n_modalities : Int32
  modality_weights : Array[Double]
}

///|
pub fn JointEmbeddingConfig::default() -> JointEmbeddingConfig {
  {
    n_latent: 10,
    n_hidden: 128,
    n_modalities: 2,
    modality_weights: [0.5, 0.5],
  }
}

///|
pub struct JointEmbeddingModel {
  config : JointEmbeddingConfig
  modality_dims : Array[Int32]
  encoders : Array[Array[Array[Double]]]
  joint_decoder : Array[Array[Double]]
  mut latent : Array[Array[Double]]
  mut is_trained : Bool
}

///|
pub fn joint_embedding_init(
  modality_dims : Array[Int32],
  config : JointEmbeddingConfig,
) -> JointEmbeddingModel {
  let encs = Array(capacity=config.n_modalities)
  for m in 0..<config.n_modalities {
    let dim = if m < modality_dims.length() { modality_dims[m] } else { 0 }
    let enc = Array(capacity=config.n_latent)
    for _i in 0..<config.n_latent {
      let row = Array(capacity=dim)
      for _j in 0..<dim {
        row.push(0.01)
      }
      enc.push(row)
    }
    encs.push(enc)
  }
  let total_dim = {
    let mut s = 0
    for i in 0..<modality_dims.length() {
      s = s + modality_dims[i]
    }
    s
  }
  let dec = Array(capacity=total_dim)
  for _i in 0..<total_dim {
    let row = Array(capacity=config.n_latent)
    for _j in 0..<config.n_latent {
      row.push(0.01)
    }
    dec.push(row)
  }
  {
    config,
    modality_dims,
    encoders: encs,
    joint_decoder: dec,
    latent: [],
    is_trained: false,
  }
}

///|
pub fn joint_embedding_encode(
  model : JointEmbeddingModel,
  modality_data : Array[Array[Double]],
) -> Array[Double] {
  let z = Array(capacity=model.config.n_latent)
  for k in 0..<model.config.n_latent {
    let mut sum = 0.0
    for m in 0..<model.config.n_modalities {
      if m < modality_data.length() && m < model.encoders.length() {
        let weight = if m < model.config.modality_weights.length() {
          model.config.modality_weights[m]
        } else {
          1.0
        }
        for j in 0..<modality_data[m].length() {
          if j < model.encoders[m][k].length() {
            sum = sum + weight * modality_data[m][j] * model.encoders[m][k][j]
          }
        }
      }
    }
    z.push(sum)
  }
  z
}

///|
pub fn joint_embedding_decode(
  model : JointEmbeddingModel,
  z : Array[Double],
) -> Array[Double] {
  let x_hat = Array(capacity=model.joint_decoder.length())
  for g in 0..<model.joint_decoder.length() {
    let mut sum = 0.0
    for j in 0..<z.length() {
      if j < model.joint_decoder[g].length() {
        sum = sum + z[j] * model.joint_decoder[g][j]
      }
    }
    x_hat.push(sum)
  }
  x_hat
}

///|
pub fn joint_embedding_train(
  model : JointEmbeddingModel,
  data : Array[Array[Array[Double]]],
  n_epochs : Int32,
) -> JointEmbeddingModel {
  let m = model
  let all_z = Array(capacity=data.length())
  for _epoch in 0..<n_epochs {
    for i in 0..<data.length() {
      let z = joint_embedding_encode(m, data[i])
      let x_hat = joint_embedding_decode(m, z)
      let mut idx = 0
      for mod in 0..<m.config.n_modalities {
        if mod < data[i].length() {
          for g in 0..<data[i][mod].length() {
            if idx < x_hat.length() {
              let grad = data[i][mod][g] - x_hat[idx]
              for j in 0..<m.config.n_latent {
                if j < m.joint_decoder[idx].length() {
                  m.joint_decoder[idx][j] = m.joint_decoder[idx][j] +
                    0.001 * grad * z[j]
                }
              }
            }
            idx = idx + 1
          }
        }
      }
    }
  }
  for i in 0..<data.length() {
    all_z.push(joint_embedding_encode(m, data[i]))
  }
  m.latent = all_z
  m.is_trained = true
  m
}

// VIVS: Variational inference for spatial data

///|
pub struct VivsConfig {
  n_latent : Int32
  n_hidden : Int32
  n_neighbors : Int32
  spatial_weight : Double
}

///|
pub fn VivsConfig::default() -> VivsConfig {
  { n_latent: 10, n_hidden: 128, n_neighbors: 6, spatial_weight: 0.1, }
}

///|
pub struct VivsModel {
  config : VivsConfig
  n_genes : Int32
  n_spots : Int32
  encoder : Array[Array[Double]]
  decoder : Array[Array[Double]]
  spatial_neighbors : Array[Array[Int32]]
  mut latent : Array[Array[Double]]
  mut is_trained : Bool
}

///|
pub fn vivs_init(
  n_genes : Int32,
  n_spots : Int32,
  spatial_neighbors : Array[Array[Int32]],
  config : VivsConfig,
) -> VivsModel {
  let make_w = (r : Int32, c : Int32) => {
    let w = Array(capacity=r)
    for _i in 0..<r {
      let row = Array(capacity=c)
      for _j in 0..<c {
        row.push(0.01)
      }
      w.push(row)
    }
    w
  }
  {
    config,
    n_genes,
    n_spots,
    spatial_neighbors,
    encoder: make_w(config.n_latent, n_genes),
    decoder: make_w(n_genes, config.n_latent),
    latent: [],
    is_trained: false,
  }
}

///|
pub fn vivs_encode(model : VivsModel, x : Array[Double]) -> Array[Double] {
  let z = Array(capacity=model.config.n_latent)
  for i in 0..<model.config.n_latent {
    let mut sum = 0.0
    for j in 0..<x.length() {
      if j < model.encoder[i].length() {
        sum = sum + x[j] * model.encoder[i][j]
      }
    }
    z.push(sum)
  }
  z
}

///|
pub fn vivspatial_loss(model : VivsModel, spot_idx : Int32) -> Double {
  let mut loss = 0.0
  if spot_idx < 0 || spot_idx >= model.spatial_neighbors.length() {
    return loss
  }
  if spot_idx >= model.latent.length() {
    return loss
  }
  let z = model.latent[spot_idx]
  for n in 0..<model.spatial_neighbors[spot_idx].length() {
    let neighbor = model.spatial_neighbors[spot_idx][n]
    if neighbor >= 0 && neighbor < model.latent.length() {
      for d in 0..<z.length() {
        if d < model.latent[neighbor].length() {
          let diff = z[d] - model.latent[neighbor][d]
          loss = loss + diff * diff
        }
      }
    }
  }
  loss * model.config.spatial_weight
}

///|
pub fn vivs_train(
  model : VivsModel,
  data : Array[Array[Double]],
  n_epochs : Int32,
) -> VivsModel {
  let m = model
  let all_z = Array(capacity=data.length())
  for _epoch in 0..<n_epochs {
    for i in 0..<data.length() {
      let z = vivs_encode(m, data[i])
      let x_hat = Array(capacity=m.n_genes)
      for g in 0..<m.n_genes {
        let mut sum = 0.0
        for j in 0..<z.length() {
          if j < m.decoder[g].length() {
            sum = sum + z[j] * m.decoder[g][j]
          }
        }
        x_hat.push(@math.exp(sum))
      }
      for g in 0..<m.n_genes {
        if g < data[i].length() && x_hat[g] > 0.0 {
          let grad = data[i][g] / x_hat[g] - 1.0
          for j in 0..<m.config.n_latent {
            if j < m.decoder[g].length() {
              m.decoder[g][j] = m.decoder[g][j] + 0.001 * grad * z[j]
            }
          }
        }
      }
    }
  }
  for i in 0..<data.length() {
    all_z.push(vivs_encode(m, data[i]))
  }
  m.latent = all_z
  m.is_trained = true
  m
}
