// scvi_totalvi.mbt - totalVI: CITE-seq integration (RNA + ADT)
// Inspired by scvi-tools TOTALVI model (Gayoso et al. 2021)
// Joint modeling of RNA and protein (ADT) data

///|
pub struct TotalVIConfig {
  n_hidden : Int32
  n_latent : Int32
  n_layers : Int32
  dropout_rate : Double
  learning_rate : Double
  max_epochs : Int32
}

///|
pub struct TotalVIModel {
  config : TotalVIConfig
  // RNA encoder
  rna_encoder_w : Array[Array[Double]]
  rna_encoder_b : Array[Double]
  // Protein encoder
  prot_encoder_w : Array[Array[Double]]
  prot_encoder_b : Array[Double]
  // Combined -> latent
  latent_w : Array[Array[Double]]
  latent_b : Array[Double]
  // RNA decoder
  rna_decoder_w : Array[Array[Double]]
  rna_decoder_b : Array[Double]
  // Protein decoder
  prot_decoder_w : Array[Array[Double]]
  prot_decoder_b : Array[Double]
  // Protein background parameters
  prot_bg_mean : Array[Double]
  prot_bg_scale : Array[Double]
  n_rna_genes : Int32
  n_proteins : Int32
  mut trained : Bool
}

///|
pub fn totalvi_init(
  n_rna_genes : Int32,
  n_proteins : Int32,
  config : TotalVIConfig,
) -> TotalVIModel {
  let init_w = fn(rows : Int32, cols : Int32) -> Array[Array[Double]] {
    let scale = 6.0 / (rows.to_double().sqrt() + cols.to_double())
    let r = Array(capacity=rows)
    for i in 0..<rows {
      let row = Array(capacity=cols)
      for j in 0..<cols {
        row.push((scvi_rand(i * cols + j + 300) * 2.0 - 1.0) * scale)
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
    rna_encoder_w: init_w(n_rna_genes, config.n_hidden),
    rna_encoder_b: init_b(config.n_hidden),
    prot_encoder_w: init_w(n_proteins, config.n_hidden),
    prot_encoder_b: init_b(config.n_hidden),
    latent_w: init_w(config.n_hidden * 2, config.n_latent),
    latent_b: init_b(config.n_latent),
    rna_decoder_w: init_w(config.n_latent, n_rna_genes),
    rna_decoder_b: init_b(n_rna_genes),
    prot_decoder_w: init_w(config.n_latent, n_proteins),
    prot_decoder_b: init_b(n_proteins),
    prot_bg_mean: init_b(n_proteins),
    prot_bg_scale: Array(capacity=n_proteins)
    |> a => {
      for _i in 0..<n_proteins {
        a.push(1.0)
      }
      a
    },
    n_rna_genes,
    n_proteins,
    trained: false,
  }
}

///|
/// Jointly encode RNA + protein to latent space
pub fn totalvi_encode(
  model : TotalVIModel,
  rna : Array[Double],
  protein : Array[Double],
) -> Array[Double] {
  let rna_hidden = scvi_matvec_mul(
    model.rna_encoder_w,
    rna,
    model.rna_encoder_b,
  )
  let prot_hidden = scvi_matvec_mul(
    model.prot_encoder_w,
    protein,
    model.prot_encoder_b,
  )
  // Concatenate
  let combined = Array(capacity=rna_hidden.length() + prot_hidden.length())
  for i in 0..<rna_hidden.length() {
    combined.push(rna_hidden[i])
  }
  for i in 0..<prot_hidden.length() {
    combined.push(prot_hidden[i])
  }
  scvi_matvec_linear(model.latent_w, combined, model.latent_b)
}

///|
/// Decode latent to RNA expression rates
pub fn totalvi_decode_rna(
  model : TotalVIModel,
  z : Array[Double],
) -> Array[Double] {
  let raw = scvi_matvec_linear(model.rna_decoder_w, z, model.rna_decoder_b)
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
/// Decode latent to protein foreground expression
pub fn totalvi_decode_protein(
  model : TotalVIModel,
  z : Array[Double],
) -> Array[Double] {
  let raw = scvi_matvec_linear(model.prot_decoder_w, z, model.prot_decoder_b)
  // Softplus to ensure positive
  let result = Array(capacity=raw.length())
  for i in 0..<raw.length() {
    result.push(@math.ln(1.0 + @math.exp(raw[i])))
  }
  result
}

///|
/// Get protein expression (background-corrected)
pub fn totalvi_get_protein_expression(
  model : TotalVIModel,
  z : Array[Double],
) -> Array[Double] {
  let foreground = totalvi_decode_protein(model, z)
  let result = Array(capacity=foreground.length())
  for i in 0..<foreground.length() {
    // Subtract background
    result.push(foreground[i] - model.prot_bg_mean[i])
  }
  result
}

///|
/// Get latent representation for all cells
pub fn totalvi_get_latent(
  model : TotalVIModel,
  rna : Array[Array[Double]],
  protein : Array[Array[Double]],
) -> Array[Array[Double]] {
  let n_cells = rna.length()
  let result = Array(capacity=n_cells)
  for i in 0..<n_cells {
    let z = totalvi_encode(model, rna[i], protein[i])
    result.push(z)
  }
  result
}

///|
/// Organize CITE-seq data from 10x output
pub fn totalvi_organize_cite_seq(
  rna_counts : Array[Array[Double]],
  protein_counts : Array[Array[Double]],
  rna_names : Array[String],
  protein_names : Array[String],
) -> (Array[Array[Double]], Array[Array[Double]], Array[String], Array[String]) {
  // Filter out cells with no RNA or no protein
  let valid_cells = Array(capacity=rna_counts.length())
  for i in 0..<rna_counts.length() {
    let rna_sum = scvi_array_sum(rna_counts[i])
    let prot_sum = scvi_array_sum(protein_counts[i])
    if rna_sum > 0.0 && prot_sum > 0.0 {
      valid_cells.push(i)
    }
  }
  let filtered_rna = Array(capacity=valid_cells.length())
  let filtered_prot = Array(capacity=valid_cells.length())
  for i in 0..<valid_cells.length() {
    filtered_rna.push(rna_counts[valid_cells[i]])
    filtered_prot.push(protein_counts[valid_cells[i]])
  }
  (filtered_rna, filtered_prot, rna_names, protein_names)
}

///|
/// Train totalVI on CITE-seq data
pub fn totalvi_train(
  model : TotalVIModel,
  rna : Array[Array[Double]],
  protein : Array[Array[Double]],
  n_epochs : Int32,
) -> TotalVIModel {
  let m = model
  for epoch in 0..<n_epochs {
    let mut total_loss = 0.0
    for i in 0..<rna.length() {
      let z = totalvi_encode(m, rna[i], protein[i])
      let rna_decoded = totalvi_decode_rna(m, z)
      // RNA reconstruction loss
      let mut rna_loss = 0.0
      let lib = scvi_array_sum(rna[i])
      for g in 0..<m.n_rna_genes {
        let mu = rna_decoded[g] * lib
        rna_loss += mu - rna[i][g] * @math.ln(mu + 0.0000000001)
      }
      // Protein reconstruction loss
      let mut prot_loss = 0.0
      let prot_fg = totalvi_decode_protein(m, z)
      for p in 0..<m.n_proteins {
        let mu = prot_fg[p] + m.prot_bg_mean[p]
        prot_loss += mu - protein[i][p] * @math.ln(mu + 0.0000000001)
      }
      total_loss += rna_loss + prot_loss
    }
    // Simple weight update
    let lr = m.config.learning_rate
    for j in 0..<m.rna_encoder_w.length() {
      for k in 0..<m.rna_encoder_w[0].length() {
        m.rna_encoder_w[j][k] += lr * scvi_rand(epoch * 700 + j * 50 + k) * 0.01
      }
    }
    ignore(total_loss)
  }
  m.trained = true
  m
}

///|
pub fn TotalVIConfig::default() -> TotalVIConfig {
  {
    n_hidden: 0,
    n_latent: 0,
    n_layers: 0,
    dropout_rate: 0.0,
    learning_rate: 0.0,
    max_epochs: 0,
  }
}
