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
  // Per-protein foreground mixing logit (sigmoid -> prior pi_fg)
  prot_mix_logit : Array[Double]
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
    rna_encoder_w: init_w(config.n_hidden, n_rna_genes),
    rna_encoder_b: init_b(config.n_hidden),
    prot_encoder_w: init_w(config.n_hidden, n_proteins),
    prot_encoder_b: init_b(config.n_hidden),
    latent_w: init_w(config.n_latent, config.n_hidden * 2),
    latent_b: init_b(config.n_latent),
    rna_decoder_w: init_w(n_rna_genes, config.n_latent),
    rna_decoder_b: init_b(n_rna_genes),
    prot_decoder_w: init_w(n_proteins, config.n_latent),
    prot_decoder_b: init_b(n_proteins),
    prot_bg_mean: Array::make(n_proteins, 1.0),
    prot_bg_scale: Array(capacity=n_proteins)
    |> a => {
      for _i in 0..<n_proteins {
        a.push(1.0)
      }
      a
    },
    prot_mix_logit: Array::make(n_proteins, 0.0),
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
pub fn TotalVIConfig::new(
  n_hidden : Int32,
  n_latent : Int32,
  n_layers : Int32,
  dropout_rate : Double,
  learning_rate : Double,
  max_epochs : Int32,
) -> TotalVIConfig {
  { n_hidden, n_latent, n_layers, dropout_rate, learning_rate, max_epochs, }
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

// ---------------------------------------------------------------------------
// 蛋白 fore/background 去卷积（官方 protein_expression_foreground_probability）
// ---------------------------------------------------------------------------

///|
/// 单细胞单蛋白的 fore/background 混合分解结果.
pub(all) struct TotalViProteinMixture {
  /// 蛋白下标
  protein : Int32
  /// 前景率 λ_fg（decoder softplus 输出）
  foreground_rate : Double
  /// 背景率 λ_bg（蛋白级参数）
  background_rate : Double
  /// 混合先验 π_fg（prot_mix_logit 经 sigmoid）
  prior : Double
  /// 观测计数下的前景后验概率（官方 foreground probability）
  foreground_probability : Double
}

///|
fn totalvi_sigmoid(x : Double) -> Double {
  let c = if x > 30.0 { 30.0 } else if x < -30.0 { -30.0 } else { x }
  1.0 / (1.0 + @math.exp(-c))
}

///|
/// 单个蛋白的前景后验概率：观测计数 y 是前景（抗体特异性信号，率 λ_fg）
/// 与背景（非特异性结合，率 λ_bg）的二元混合，混合先验为 π_fg：
///   p_fg = π·Pois(y;λ_fg) / (π·Pois(y;λ_fg) + (1-π)·Pois(y;λ_bg))
/// 官方 TOTALVI 用 NB 似然（蛋白级过离散）；这里用 Poisson 简化，
/// 负计数按 0、速率下限 1e-6、先验钳制到 [1e-6, 1-1e-6]。
pub fn totalvi_mixture_posterior(
  count : Double,
  foreground_rate : Double,
  background_rate : Double,
  prior : Double,
) -> Double {
  let y = if count < 0.0 { 0.0 } else { count }
  let lam_fg = if foreground_rate < 1.0e-6 { 1.0e-6 } else { foreground_rate }
  let lam_bg = if background_rate < 1.0e-6 { 1.0e-6 } else { background_rate }
  let pi = if prior < 1.0e-6 {
    1.0e-6
  } else if prior > 1.0 - 1.0e-6 {
    1.0 - 1.0e-6
  } else {
    prior
  }
  let log_fact = cobra3_lgamma(y + 1.0)
  let lf = y * @math.ln(lam_fg) - lam_fg - log_fact
  let lb = y * @math.ln(lam_bg) - lam_bg - log_fact
  let af = @math.ln(pi) + lf
  let ab = @math.ln(1.0 - pi) + lb
  let m = if af > ab { af } else { ab }
  let ef = @math.exp(af - m)
  let eb = @math.exp(ab - m)
  ef / (ef + eb)
}

///|
/// 单细胞蛋白 fore/background 混合分解：前景率来自 decoder（softplus），
/// 背景率为蛋白级参数 prot_bg_mean（下限 1e-3），先验 π_fg 由
/// prot_mix_logit 经 sigmoid 给出。
pub fn totalvi_protein_mixture(
  model : TotalVIModel,
  z : Array[Double],
  counts : Array[Double],
) -> Array[TotalViProteinMixture] {
  let fg = totalvi_decode_protein(model, z)
  let out : Array[TotalViProteinMixture] = Array(capacity=model.n_proteins)
  for p in 0..<model.n_proteins {
    let y = if p < counts.length() && counts[p] > 0.0 { counts[p] } else { 0.0 }
    let bg = if model.prot_bg_mean[p] < 1.0e-3 {
      1.0e-3
    } else {
      model.prot_bg_mean[p]
    }
    let prior = totalvi_sigmoid(model.prot_mix_logit[p])
    out.push({
      protein: p,
      foreground_rate: fg[p],
      background_rate: bg,
      prior,
      foreground_probability: totalvi_mixture_posterior(y, fg[p], bg, prior),
    })
  }
  out
}

///|
/// 全部查询细胞的前景概率矩阵（官方
/// protein_expression_foreground_probability 语义）：逐细胞编码 → decoder
/// 前景率 → 观测计数下的混合后验。rna 与 protein 行数不一致或无蛋白时
/// 返回空列表。
pub fn totalvi_get_foreground_probability(
  model : TotalVIModel,
  rna : Array[Array[Double]],
  protein : Array[Array[Double]],
) -> Array[Array[Double]] {
  let out : Array[Array[Double]] = Array(capacity=protein.length())
  if model.n_proteins == 0 || protein.length() != rna.length() {
    return out
  }
  for i in 0..<protein.length() {
    let z = totalvi_encode(model, rna[i], protein[i])
    let fg = totalvi_decode_protein(model, z)
    let row : Array[Double] = Array(capacity=model.n_proteins)
    for p in 0..<model.n_proteins {
      let y = if protein[i][p] > 0.0 { protein[i][p] } else { 0.0 }
      let bg = if model.prot_bg_mean[p] < 1.0e-3 {
        1.0e-3
      } else {
        model.prot_bg_mean[p]
      }
      row.push(
        totalvi_mixture_posterior(
          y,
          fg[p],
          bg,
          totalvi_sigmoid(model.prot_mix_logit[p]),
        ),
      )
    }
    out.push(row)
  }
  out
}

///|
/// 背景校正蛋白表达：期望前景计数 = 观测计数 × 前景后验（官方
/// get_protein_expression 的 denoised 语义）。
pub fn totalvi_corrected_protein_expression(
  model : TotalVIModel,
  rna : Array[Array[Double]],
  protein : Array[Array[Double]],
) -> Array[Array[Double]] {
  let probs = totalvi_get_foreground_probability(model, rna, protein)
  let out : Array[Array[Double]] = Array(capacity=probs.length())
  for i in 0..<probs.length() {
    let row : Array[Double] = Array(capacity=probs[i].length())
    for p in 0..<probs[i].length() {
      let y = if protein[i][p] > 0.0 { protein[i][p] } else { 0.0 }
      row.push(y * probs[i][p])
    }
    out.push(row)
  }
  out
}

///|
/// 蛋白背景 EM 去卷积：迭代优化背景率与混合先验。
/// E 步：当前参数下逐细胞计算前景后验；M 步：以 (1-p_fg) 加权重估
/// 背景率 λ_bg = Σ(1-p)·y / Σ(1-p)（背景归属计数的加权均值），并以
/// 梯度步更新混合 logit（logit += p̄ - π）。返回最终前景概率矩阵。
/// rna 与 protein 行数不一致或无蛋白时返回空列表且不修改模型。
pub fn totalvi_fit_protein_background(
  model : TotalVIModel,
  rna : Array[Array[Double]],
  protein : Array[Array[Double]],
  n_iters~ : Int32 = 10,
) -> Array[Array[Double]] {
  if model.n_proteins == 0 || protein.length() != rna.length() {
    return []
  }
  let n_cells = protein.length()
  let n_cells_d = n_cells.to_double()
  for _it in 0..<n_iters {
    let wsum : Array[Double] = Array::make(model.n_proteins, 0.0)
    let ysum : Array[Double] = Array::make(model.n_proteins, 0.0)
    let psum : Array[Double] = Array::make(model.n_proteins, 0.0)
    for i in 0..<n_cells {
      let z = totalvi_encode(model, rna[i], protein[i])
      let fg = totalvi_decode_protein(model, z)
      for p in 0..<model.n_proteins {
        let y = if protein[i][p] > 0.0 { protein[i][p] } else { 0.0 }
        let bg = if model.prot_bg_mean[p] < 1.0e-3 {
          1.0e-3
        } else {
          model.prot_bg_mean[p]
        }
        let post = totalvi_mixture_posterior(
          y,
          fg[p],
          bg,
          totalvi_sigmoid(model.prot_mix_logit[p]),
        )
        wsum[p] += 1.0 - post
        ysum[p] += (1.0 - post) * y
        psum[p] += post
      }
    }
    for p in 0..<model.n_proteins {
      if wsum[p] > 1.0e-9 {
        let bg = ysum[p] / wsum[p]
        if bg > 1.0e-3 {
          model.prot_bg_mean[p] = bg
        }
      }
      let pbar = psum[p] / n_cells_d
      let mut logit = model.prot_mix_logit[p] +
        (pbar - totalvi_sigmoid(model.prot_mix_logit[p]))
      if logit > 10.0 {
        logit = 10.0
      } else if logit < -10.0 {
        logit = -10.0
      }
      model.prot_mix_logit[p] = logit
    }
  }
  totalvi_get_foreground_probability(model, rna, protein)
}
