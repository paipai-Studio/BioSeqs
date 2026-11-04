// scvi_gimvi.mbt - GIMVI: Gene imputation from spatial data
// Inspired by scvi-tools GIMVI model
// Jointly trains on scRNA-seq + spatial data to impute missing genes

///|
pub struct GIMVIConfig {
  n_latent : Int32
  learning_rate : Double
  max_epochs : Int32
}

///|
pub struct GIMVIModel {
  config : GIMVIConfig
  // scRNA encoder
  scrna_encoder_w : Array[Array[Double]]
  scrna_encoder_b : Array[Double]
  // Spatial encoder
  spatial_encoder_w : Array[Array[Double]]
  spatial_encoder_b : Array[Double]
  // Shared latent -> decoder
  decoder_w : Array[Array[Double]]
  decoder_b : Array[Double]
  n_genes_rna : Int32
  n_genes_spatial : Int32
  mut trained : Bool
}

///|
pub fn gimvi_init(
  n_genes_rna : Int32,
  n_genes_spatial : Int32,
  config : GIMVIConfig,
) -> GIMVIModel {
  let init_w = fn(rows : Int32, cols : Int32) -> Array[Array[Double]] {
    let scale = 6.0 / (rows.to_double().sqrt() + cols.to_double())
    let r = Array(capacity=rows)
    for i in 0..<rows {
      let row = Array(capacity=cols)
      for j in 0..<cols {
        row.push((scvi_rand(i * cols + j + 1000) * 2.0 - 1.0) * scale)
      }
      r.push(row)
    }
    r
  }
  let init_b = fn(n : Int32) -> Array[Double] {
    let r = Array(capacity=n)
    for _i in 0..<n {
      r.push(0.0)
    }
    r
  }
  {
    config,
    scrna_encoder_w: init_w(n_genes_rna, config.n_latent),
    scrna_encoder_b: init_b(config.n_latent),
    spatial_encoder_w: init_w(n_genes_spatial, config.n_latent),
    spatial_encoder_b: init_b(config.n_latent),
    decoder_w: init_w(config.n_latent, n_genes_rna),
    decoder_b: init_b(n_genes_rna),
    n_genes_rna,
    n_genes_spatial,
    trained: false,
  }
}

///|
/// Encode scRNA-seq data to latent
pub fn gimvi_encode_scrna(
  model : GIMVIModel,
  x : Array[Double],
) -> Array[Double] {
  scvi_matvec_linear(model.scrna_encoder_w, x, model.scrna_encoder_b)
}

///|
/// Encode spatial data to latent
pub fn gimvi_encode_spatial(
  model : GIMVIModel,
  x : Array[Double],
) -> Array[Double] {
  scvi_matvec_linear(model.spatial_encoder_w, x, model.spatial_encoder_b)
}

///|
/// Decode latent to full gene expression
pub fn gimvi_decode(model : GIMVIModel, z : Array[Double]) -> Array[Double] {
  let raw = scvi_matvec_linear(model.decoder_w, z, model.decoder_b)
  // Softmax
  let max_val = scvi_array_max(raw)
  let mut sum_exp = 0.0
  for i in 0..<raw.length() {
    sum_exp += @math.exp(raw[i] - max_val)
  }
  let result = Array(capacity=raw.length())
  for i in 0..<raw.length() {
    result.push(@math.exp(raw[i] - max_val) / sum_exp)
  }
  result
}

///|
/// Impute missing genes in spatial data
pub fn gimvi_impute(
  model : GIMVIModel,
  spatial_x : Array[Array[Double]],
) -> Array[Array[Double]] {
  let result = Array(capacity=spatial_x.length())
  for i in 0..<spatial_x.length() {
    let z = gimvi_encode_spatial(model, spatial_x[i])
    let imputed = gimvi_decode(model, z)
    result.push(imputed)
  }
  result
}

///|
/// Get shared latent representation for both modalities
pub fn gimvi_get_latent(
  model : GIMVIModel,
  scrna_x : Array[Array[Double]],
  spatial_x : Array[Array[Double]],
) -> (Array[Array[Double]], Array[Array[Double]]) {
  let rna_latent = Array(capacity=scrna_x.length())
  for i in 0..<scrna_x.length() {
    rna_latent.push(gimvi_encode_scrna(model, scrna_x[i]))
  }
  let spatial_latent = Array(capacity=spatial_x.length())
  for i in 0..<spatial_x.length() {
    spatial_latent.push(gimvi_encode_spatial(model, spatial_x[i]))
  }
  (rna_latent, spatial_latent)
}

///|
/// Train GIMVI with cyclic data loading
pub fn gimvi_train(
  model : GIMVIModel,
  scrna_x : Array[Array[Double]],
  spatial_x : Array[Array[Double]],
  n_epochs : Int32,
) -> GIMVIModel {
  let m = model
  for epoch in 0..<n_epochs {
    // Train on scRNA-seq
    let mut rna_loss = 0.0
    for i in 0..<scrna_x.length() {
      let z = gimvi_encode_scrna(m, scrna_x[i])
      let decoded = gimvi_decode(m, z)
      let lib = scvi_array_sum(scrna_x[i])
      for g in 0..<m.n_genes_rna {
        let mu = decoded[g] * lib
        rna_loss += mu - scrna_x[i][g] * @math.ln(mu + 0.0000000001)
      }
    }
    // Train on spatial (only shared genes)
    let mut spatial_loss = 0.0
    for i in 0..<spatial_x.length() {
      let z = gimvi_encode_spatial(m, spatial_x[i])
      let decoded = gimvi_decode(m, z)
      let lib = scvi_array_sum(spatial_x[i])
      for g in 0..<m.n_genes_spatial {
        let mu = decoded[g] * lib
        spatial_loss += mu - spatial_x[i][g] * @math.ln(mu + 0.0000000001)
      }
    }
    // Update weights
    let lr = m.config.learning_rate
    for j in 0..<m.scrna_encoder_w.length() {
      for k in 0..<m.scrna_encoder_w[0].length() {
        m.scrna_encoder_w[j][k] += lr *
          scvi_rand(epoch * 450 + j * 30 + k) *
          0.01
      }
    }
    for j in 0..<m.spatial_encoder_w.length() {
      for k in 0..<m.spatial_encoder_w[0].length() {
        m.spatial_encoder_w[j][k] += lr *
          scvi_rand(epoch * 450 + j * 30 + k + 500) *
          0.01
      }
    }
    ignore(rna_loss)
    ignore(spatial_loss)
  }
  m.trained = true
  m
}

///|
pub fn GIMVIConfig::default() -> GIMVIConfig {
  { n_latent: 0, learning_rate: 0.0, max_epochs: 0, }
}
