///|
/// CondSCVI — conditional single-cell VAE (Lopez et al., the
/// `scvi.model.CondSCVI` reference side of DestVI). Trains a cell-type
/// conditional NB-VAE on the single-cell reference: the one-hot cell type
/// embedding is injected into both encoder and decoder. The learned per-cell
/// type expression profiles (`gene × cell_type`) feed `destvi_init`, closing
/// the cells2spatial deconvolution loop with the existing DestVI model.

///|
pub struct CondScviConfig {
  n_hidden : Int32
  n_latent : Int32
  learning_rate : Double
  max_epochs : Int32
}

///|
pub struct CondScviModel {
  config : CondScviConfig
  // encoder: expression branch + cell-type embedding branch
  encoder_w1 : Array[Array[Double]] // n_hidden × n_genes
  encoder_b1 : Array[Double]
  enc_ct_w : Array[Array[Double]] // n_hidden × n_cell_types
  enc_ct_b : Array[Double]
  encoder_w_mu : Array[Array[Double]] // n_latent × n_hidden
  encoder_b_mu : Array[Double]
  encoder_w_lv : Array[Array[Double]]
  encoder_b_lv : Array[Double]
  // decoder: latent branch + cell-type embedding branch
  decoder_w1 : Array[Array[Double]] // n_hidden × n_latent
  decoder_b1 : Array[Double]
  dec_ct_w : Array[Array[Double]]
  dec_ct_b : Array[Double]
  decoder_w_out : Array[Array[Double]] // n_genes × n_hidden
  decoder_b_out : Array[Double]
  log_theta : Array[Double] // gene-wise inverse dispersion (fixed)
  // per-cell-type reference profiles (gene × cell_type), filled after training
  mut profiles : Array[Array[Double]]
  n_cell_types : Int32
  n_genes : Int32
  mut trained : Bool
  mut final_loss : Double
}

///|
pub fn CondScviConfig::default() -> CondScviConfig {
  CondScviConfig::{
    n_hidden: 32,
    n_latent: 5,
    learning_rate: 0.01,
    max_epochs: 60,
  }
}

///|
fn condscvi_xavier(
  rows : Int32,
  cols : Int32,
  base_seed : Int32,
) -> Array[Array[Double]] {
  let scale = 6.0 / (rows.to_double().sqrt() + cols.to_double())
  let result = Array(capacity=rows)
  for i = 0; i < rows; i = i + 1 {
    let row : Array[Double] = Array(capacity=cols)
    for j = 0; j < cols; j = j + 1 {
      row.push(scvi_rand(base_seed + i * cols + j) * 2.0 * scale - scale)
    }
    result.push(row)
  }
  result
}

///|
fn condscvi_zeros(n : Int32) -> Array[Double] {
  Array::make(n, 0.0)
}

///|
/// Initialize an untrained CondSCVI model (deterministic Xavier init).
pub fn scvi_condscvi_init_model(
  n_genes : Int32,
  n_cell_types : Int32,
  config~ : CondScviConfig = CondScviConfig::default(),
) -> CondScviModel {
  let h = config.n_hidden
  let z = config.n_latent
  CondScviModel::{
    config,
    encoder_w1: condscvi_xavier(h, n_genes, 42),
    encoder_b1: condscvi_zeros(h),
    enc_ct_w: condscvi_xavier(h, n_cell_types, 1000),
    enc_ct_b: condscvi_zeros(h),
    encoder_w_mu: condscvi_xavier(z, h, 2000),
    encoder_b_mu: condscvi_zeros(z),
    encoder_w_lv: condscvi_xavier(z, h, 3000),
    encoder_b_lv: condscvi_zeros(z),
    decoder_w1: condscvi_xavier(h, z, 4000),
    decoder_b1: condscvi_zeros(h),
    dec_ct_w: condscvi_xavier(h, n_cell_types, 5000),
    dec_ct_b: condscvi_zeros(h),
    decoder_w_out: condscvi_xavier(n_genes, h, 6000),
    decoder_b_out: condscvi_zeros(n_genes),
    log_theta: condscvi_zeros(n_genes),
    profiles: Array(capacity=0),
    n_cell_types,
    n_genes,
    trained: false,
    final_loss: 0.0,
  }
}

///|
/// One-hot cell-type vector.
fn condscvi_onehot(ct : Int32, k : Int32) -> Array[Double] {
  let v = Array::make(k, 0.0)
  let i = if ct < 0 { 0 } else if ct >= k { k - 1 } else { ct }
  v[i] = 1.0
  v
}

///|
/// Encoder input transform (mirrors scvi-tools CondSCVI): library-normalized
/// log1p expression. Bounds activations and keeps training stable.
fn condscvi_encoder_input(x : Array[Double]) -> Array[Double] {
  let lib = scvi_array_sum(x)
  let lib_size = if lib <= 0.0 { 1.0 } else { lib }
  let out = Array(capacity=x.length())
  for g = 0; g < x.length(); g = g + 1 {
    out.push(@math.ln(1.0 + x[g] / lib_size))
  }
  out
}

///|
/// Non-tape encoder forward: (mu, log_var) for one cell given its cell type.
pub fn scvi_condscvi_encode(
  model : CondScviModel,
  x : Array[Double],
  ct : Int32,
) -> (Array[Double], Array[Double]) {
  let xin = condscvi_encoder_input(x)
  let h_expr = scvi_matvec_mul(model.encoder_w1, xin, model.encoder_b1)
  let h_ct = scvi_matvec_mul(
    model.enc_ct_w,
    condscvi_onehot(ct, model.n_cell_types),
    model.enc_ct_b,
  )
  let act : Array[Double] = Array(capacity=h_expr.length())
  for i = 0; i < h_expr.length(); i = i + 1 {
    act.push(scvi_relu(h_expr[i] + h_ct[i]))
  }
  let mu = scvi_matvec_mul(model.encoder_w_mu, act, model.encoder_b_mu)
  let log_var = scvi_matvec_mul(model.encoder_w_lv, act, model.encoder_b_lv)
  (mu, log_var)
}

///|
/// Non-tape decoder forward: normalized rate vector for a latent code and
/// cell type.
pub fn scvi_condscvi_decode(
  model : CondScviModel,
  z : Array[Double],
  ct : Int32,
) -> Array[Double] {
  let d_expr = scvi_matvec_mul(model.decoder_w1, z, model.decoder_b1)
  let d_ct = scvi_matvec_mul(
    model.dec_ct_w,
    condscvi_onehot(ct, model.n_cell_types),
    model.dec_ct_b,
  )
  let act : Array[Double] = Array(capacity=d_expr.length())
  for i = 0; i < d_expr.length(); i = i + 1 {
    act.push(scvi_relu(d_expr[i] + d_ct[i]))
  }
  let raw = scvi_matvec_linear(model.decoder_w_out, act, model.decoder_b_out)
  let max_val = scvi_array_max(raw)
  let mut sum_exp = 0.0
  for i = 0; i < raw.length(); i = i + 1 {
    sum_exp += @math.exp(raw[i] - max_val)
  }
  let rates = Array(capacity=raw.length())
  for i = 0; i < raw.length(); i = i + 1 {
    rates.push(@math.exp(raw[i] - max_val) / sum_exp)
  }
  rates
}

///|
/// Latent representation (posterior mean) of every cell.
pub fn scvi_condscvi_latent(
  model : CondScviModel,
  x : Array[Array[Double]],
  cell_types : Array[Int32],
) -> Array[Array[Double]] {
  let out : Array[Array[Double]] = Array(capacity=x.length())
  for i = 0; i < x.length(); i = i + 1 {
    let ct = if i < cell_types.length() { cell_types[i] } else { 0 }
    let (mu, _lv) = scvi_condscvi_encode(model, x[i], ct)
    out.push(mu)
  }
  out
}

///|
/// Denoised expression per cell: rate vector scaled by the cell library size.
pub fn scvi_condscvi_normalized_expression(
  model : CondScviModel,
  x : Array[Array[Double]],
  cell_types : Array[Int32],
) -> Array[Array[Double]] {
  let out : Array[Array[Double]] = Array(capacity=x.length())
  for i = 0; i < x.length(); i = i + 1 {
    let ct = if i < cell_types.length() { cell_types[i] } else { 0 }
    let (mu, _lv) = scvi_condscvi_encode(model, x[i], ct)
    let rates = scvi_condscvi_decode(model, mu, ct)
    let lib = scvi_array_sum(x[i])
    let row : Array[Double] = Array(capacity=rates.length())
    for g = 0; g < rates.length(); g = g + 1 {
      row.push(rates[g] * lib)
    }
    out.push(row)
  }
  out
}

///|
/// Per-cell-type reference profiles (gene × cell_type): decode the type-mean
/// latent at the type-mean library size. This replaces the naive
/// `destvi_build_reference` mean-count profiles with a denoised, learned
/// signature suitable for `destvi_init`.
pub fn scvi_condscvi_per_cell_type_expression(
  model : CondScviModel,
  x : Array[Array[Double]],
  cell_types : Array[Int32],
) -> Array[Array[Double]] {
  // gene × cell_type
  let profiles : Array[Array[Double]] = Array(capacity=model.n_genes)
  for _g in 0..<model.n_genes {
    profiles.push(Array::make(model.n_cell_types, 0.0))
  }
  for ct = 0; ct < model.n_cell_types; ct = ct + 1 {
    let mut lib_sum = 0.0
    let z_sum = Array::make(model.config.n_latent, 0.0)
    let mut n_ct = 0
    for i = 0; i < x.length(); i = i + 1 {
      let c = if i < cell_types.length() { cell_types[i] } else { 0 }
      if c == ct {
        let (mu, _lv) = scvi_condscvi_encode(model, x[i], ct)
        for z = 0; z < z_sum.length() && z < mu.length(); z = z + 1 {
          z_sum[z] = z_sum[z] + mu[z]
        }
        lib_sum = lib_sum + scvi_array_sum(x[i])
        n_ct = n_ct + 1
      }
    }
    if n_ct > 0 {
      let z_mean : Array[Double] = Array(capacity=z_sum.length())
      for z = 0; z < z_sum.length(); z = z + 1 {
        z_mean.push(z_sum[z] / n_ct.to_double())
      }
      let rates = scvi_condscvi_decode(model, z_mean, ct)
      let lib_mean = lib_sum / n_ct.to_double()
      for g = 0; g < model.n_genes && g < rates.length(); g = g + 1 {
        profiles[g][ct] = rates[g] * lib_mean
      }
    }
  }
  profiles
}

///|
/// Train the conditional VAE with the autograd tape + Adam (mirrors
/// `scvi_train`): NB likelihood on softmax rates scaled by the cell library,
/// KL against the standard normal prior, deterministic reparameterization.
pub fn scvi_condscvi_train(
  model : CondScviModel,
  x : Array[Array[Double]],
  cell_types : Array[Int32],
) -> CondScviModel {
  if x.length() == 0 {
    model.profiles = Array(capacity=0)
    model.trained = true
    return model
  }
  let m = model
  let n_cells = x.length().to_double()
  let theta = Array(capacity=m.n_genes)
  for g in 0..<m.n_genes {
    theta.push(@math.exp(m.log_theta[g]))
  }
  // Adam template: parameter registration order defines the optimizer state.
  let template = ScviAgTape::new()
  let _register = fn(t : ScviAgTape, mm : CondScviModel) -> Unit {
    ignore(t.matrix_param(mm.encoder_w1))
    ignore(t.vector_param(mm.encoder_b1))
    ignore(t.matrix_param(mm.enc_ct_w))
    ignore(t.vector_param(mm.enc_ct_b))
    ignore(t.matrix_param(mm.encoder_w_mu))
    ignore(t.vector_param(mm.encoder_b_mu))
    ignore(t.matrix_param(mm.encoder_w_lv))
    ignore(t.vector_param(mm.encoder_b_lv))
    ignore(t.matrix_param(mm.decoder_w1))
    ignore(t.vector_param(mm.decoder_b1))
    ignore(t.matrix_param(mm.dec_ct_w))
    ignore(t.vector_param(mm.dec_ct_b))
    ignore(t.matrix_param(mm.decoder_w_out))
    ignore(t.vector_param(mm.decoder_b_out))
  }
  _register(template, m)
  let adam = scvi_adam_init(template, m.config.learning_rate)
  for _epoch in 0..<m.config.max_epochs {
    let t = ScviAgTape::new()
    let w1 = t.matrix_param(m.encoder_w1)
    let b1 = t.vector_param(m.encoder_b1)
    let w_ct = t.matrix_param(m.enc_ct_w)
    let b_ct = t.vector_param(m.enc_ct_b)
    let w_mu = t.matrix_param(m.encoder_w_mu)
    let b_mu = t.vector_param(m.encoder_b_mu)
    let w_lv = t.matrix_param(m.encoder_w_lv)
    let b_lv = t.vector_param(m.encoder_b_lv)
    let wd1 = t.matrix_param(m.decoder_w1)
    let bd1 = t.vector_param(m.decoder_b1)
    let w_ctd = t.matrix_param(m.dec_ct_w)
    let b_ctd = t.vector_param(m.dec_ct_b)
    let w_out = t.matrix_param(m.decoder_w_out)
    let b_out = t.vector_param(m.decoder_b_out)
    let mut total = t.leaf([0.0])
    for i in 0..<x.length() {
      let xi = x[i]
      let lib = scvi_array_sum(xi)
      let lib_size = if lib <= 0.0 { 1.0 } else { lib }
      let ct = if i < cell_types.length() { cell_types[i] } else { 0 }
      let oh = t.leaf(condscvi_onehot(ct, m.n_cell_types))
      let xin = t.leaf(condscvi_encoder_input(xi))
      // encoder: expression branch + cell-type branch
      let h = t.dense(xin, w1, b1)
      let h_ct = t.dense(oh, w_ct, b_ct)
      let act = t.relu(t.add(h, h_ct))
      let mu = t.dense(act, w_mu, b_mu)
      let log_var = t.dense(act, w_lv, b_lv)
      // deterministic reparameterization
      let d = t.dense(mu, wd1, bd1)
      let d_ct = t.dense(oh, w_ctd, b_ctd)
      let dact = t.relu(t.add(d, d_ct))
      let raw = t.dense(dact, w_out, b_out)
      let rates = t.softmax(raw)
      let mean_counts = t.scale(rates, lib_size)
      let nll = t.nb_nll(mean_counts, xi, theta)
      let kl = t.kl_divergence(mu, log_var)
      total = t.add(total, t.add(nll, kl))
    }
    let loss = t.scale(total, 1.0 / n_cells)
    t.backward(loss)
    ignore(t.clip_grad_norm(5.0))
    scvi_adam_step(t, adam)
    m.final_loss = t.value(loss)[0]
  }
  m.profiles = scvi_condscvi_per_cell_type_expression(m, x, cell_types)
  m.trained = true
  m
}

///|
/// cells2spatial bridge: use the learned CondSCVI reference profiles to run
/// the existing DestVI deconvolution and return spot × cell_type proportions.
pub fn scvi_cells2spatial(
  model : CondScviModel,
  spatial_x : Array[Array[Double]],
) -> Array[Array[Double]] {
  let n_spots = spatial_x.length()
  let dest = destvi_init(model.profiles, n_spots, DestVIConfig::{
    n_latent: model.config.n_latent,
    n_cell_types: model.n_cell_types,
    learning_rate: 0.01,
    max_epochs: 50,
  })
  destvi_deconvolve(dest, spatial_x)
}
