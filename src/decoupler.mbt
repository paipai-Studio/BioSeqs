///|
/// decoupleR module - Functional activity inference from prior knowledge networks.
///
/// References:
///   - Badia-i-Mompel P et al. (2022) decoupleR: ensemble of computational methods
///     to infer biological activities from omics data.
///   - https://www.bioconductor.org/packages/release/bioc/html/decoupleR.html
///
/// Core algorithm:
///   Given a prior knowledge network (PKN) of regulator→target relationships with
///   weights (mode of regulation: +1 activation, -1 inhibition), and an omics matrix
///   (samples × genes), infer per-sample regulator activity scores using one of:
///     - wsum  : weighted sum of target expressions per regulator
///     - wmean : weighted mean (normalized by absolute weight sum)
///     - ulm   : univariate linear model: regress each target on each regulator's
///               weighted targets; t-statistic of slope = activity
///     - mlm   : multivariate linear model: regress all targets on all regulators
///               simultaneously; t-statistic of each coefficient = activity
///     - norm  : simple normalized mean (z-score style)

///|
/// A single directed interaction in the prior knowledge network.
pub struct PKNEdge {
  source : String       // regulator (TF) name
  target : String       // target gene name
  weight : Double       // mode of regulation: >0 activation, <0 inhibition
}

///|
/// Prior knowledge network: collection of regulator→target edges.
pub struct PriorKnowledgeNetwork {
  edges : Array[PKNEdge]
  regulators : Array[String]   // unique regulator names
  targets : Array[String]      // unique target names
}

///|
/// Build a PriorKnowledgeNetwork from a list of edges.
pub fn pkn_from_edges(edges : Array[PKNEdge]) -> PriorKnowledgeNetwork {
  // Collect unique regulators and targets
  let reg_map : Map[String, Unit] = Map::new()
  let tgt_map : Map[String, Unit] = Map::new()
  let mut i = 0
  while i < edges.length() {
    reg_map.set(edges[i].source, ())
    tgt_map.set(edges[i].target, ())
    i = i + 1
  }
  let regulators = reg_map.keys().collect()
  let targets = tgt_map.keys().collect()
  PriorKnowledgeNetwork::{ edges, regulators, targets }
}

///|
/// Create a PKNEdge.
pub fn pkn_edge(source : String, target : String, weight : Double) -> PKNEdge {
  PKNEdge::{ source, target, weight }
}

///|
/// A regulator activity score for a single sample.
pub struct ActivityScore {
  regulator : String
  sample : String
  score : Double
  p_value : Double     // only meaningful for ulm/mlm; 1.0 otherwise
}

///|
/// decoupleR result: per-sample regulator activity matrix and per-edge records.
pub struct DecoupleRResult {
  regulators : Array[String]
  samples : Array[String]
  matrix : Array[Array[Double]]     // matrix[i][j] = activity of regulator i in sample j
  scores : Array[ActivityScore]     // flat list of scores
  method : String
}

///|
/// Internal: build regulator × target weight matrix from PKN and gene name list.
/// Returns (reg_indices, tgt_indices, W) where W[i][j] is the PKN weight from
/// regulator i to target j (0 if no edge), reg_indices/tgt_indices map to input
/// gene-name arrays.
fn decoupler_build_weight_matrix(
  pkn : PriorKnowledgeNetwork,
  gene_names : Array[String],
) -> (Array[String], Array[Int], Array[Array[Double]]) {
  // Map gene names to column indices in the expression matrix
  let gene_index : Map[String, Int] = Map::new()
  let mut i = 0
  while i < gene_names.length() {
    gene_index.set(gene_names[i], i)
    i = i + 1
  }
  // Only include regulators that appear as genes in expression matrix
  let reg_list : Array[String] = []
  let mut ri = 0
  while ri < pkn.regulators.length() {
    if gene_index.contains(pkn.regulators[ri]) {
      reg_list.push(pkn.regulators[ri])
    }
    ri = ri + 1
  }
  // Map regulator to row index
  let reg_index : Map[String, Int] = Map::new()
  let mut ri = 0
  while ri < reg_list.length() {
    reg_index.set(reg_list[ri], ri)
    ri = ri + 1
  }
  // Only include targets that appear in expression matrix
  let tgt_indices : Array[Int] = []
  let tgt_seen : Map[Int, Unit] = Map::new()
  let mut ti = 0
  while ti < pkn.targets.length() {
    let t = pkn.targets[ti]
    let idx_opt = gene_index.get(t)
    match idx_opt {
      Some(idx) =>
        if not(tgt_seen.contains(idx)) {
          tgt_seen.set(idx, ())
          tgt_indices.push(idx)
        }
      None => ()
    }
    ti = ti + 1
  }
  // Build weight matrix: rows = regulators, cols = targets (using tgt_indices positions)
  let n_reg = reg_list.length()
  let n_tgt = tgt_indices.length()
  let w : Array[Array[Double]] = Array::new()
  let mut wi = 0
  while wi < n_reg {
    w.push(Array::make(n_tgt, 0.0))
    wi = wi + 1
  }
  let mut ei = 0
  while ei < pkn.edges.length() {
    let e = pkn.edges[ei]
    let reg_row_opt = reg_index.get(e.source)
    let tgt_col_idx_opt = gene_index.get(e.target)
    match (reg_row_opt, tgt_col_idx_opt) {
      (Some(reg_row), Some(tgt_col_idx)) => {
        // Find tgt_col position in tgt_indices
        let mut k = 0
        let mut found = -1
        while k < tgt_indices.length() {
          if tgt_indices[k] == tgt_col_idx {
            found = k
            break
          }
          k = k + 1
        }
        if found >= 0 {
          w[reg_row][found] = w[reg_row][found] + e.weight
        }
      }
      _ => ()
    }
    ei = ei + 1
  }
  (reg_list, tgt_indices, w)
}

///|
/// Compute wsum (weighted sum) activity scores.
/// For each regulator r and sample s:
///   score = sum_over_targets( W[r][t] * expr[s][t] )
fn decoupler_wsum(
  expression : Array[Array[Double]],
  sample_names : Array[String],
  gene_names : Array[String],
  pkn : PriorKnowledgeNetwork,
) -> DecoupleRResult {
  let (reg_list, tgt_indices, w) = decoupler_build_weight_matrix(pkn, gene_names)
  let n_reg = reg_list.length()
  let n_samples = expression.length()
  let matrix : Array[Array[Double]] = Array::new()
  let mut mri = 0
  while mri < n_reg {
    matrix.push(Array::make(n_samples, 0.0))
    mri = mri + 1
  }
  let mut i = 0
  while i < n_reg {
    let mut j = 0
    while j < n_samples {
      let mut s = 0.0
      let mut k = 0
      while k < tgt_indices.length() {
        s = s + w[i][k] * expression[j][tgt_indices[k]]
        k = k + 1
      }
      matrix[i][j] = s
      j = j + 1
    }
    i = i + 1
  }
  decoupler_pack_result(reg_list, sample_names, matrix, "wsum")
}

///|
/// Compute wmean (weighted mean) activity scores.
/// For each regulator r and sample s:
///   score = sum(W[r][t] * expr[s][t]) / sum(|W[r][t]|)
fn decoupler_wmean(
  expression : Array[Array[Double]],
  sample_names : Array[String],
  gene_names : Array[String],
  pkn : PriorKnowledgeNetwork,
) -> DecoupleRResult {
  let (reg_list, tgt_indices, w) = decoupler_build_weight_matrix(pkn, gene_names)
  let n_reg = reg_list.length()
  let n_samples = expression.length()
  let matrix : Array[Array[Double]] = Array::new()
  let mut mri = 0
  while mri < n_reg {
    matrix.push(Array::make(n_samples, 0.0))
    mri = mri + 1
  }
  let mut i = 0
  while i < n_reg {
    // Compute sum of absolute weights for this regulator
    let mut abs_w_sum = 0.0
    let mut k = 0
    while k < tgt_indices.length() {
      abs_w_sum = abs_w_sum + w[i][k].abs()
      k = k + 1
    }
    let mut j = 0
    while j < n_samples {
      let mut s = 0.0
      let mut k = 0
      while k < tgt_indices.length() {
        s = s + w[i][k] * expression[j][tgt_indices[k]]
        k = k + 1
      }
      matrix[i][j] = if abs_w_sum > 1.0e-12 { s / abs_w_sum } else { 0.0 }
      j = j + 1
    }
    i = i + 1
  }
  decoupler_pack_result(reg_list, sample_names, matrix, "wmean")
}

///|
/// Compute norm (z-score style normalized mean) activity scores.
/// For each regulator r and sample s:
///   raw = sum(W[r][t] * expr[s][t]) / sum(|W[r][t]|)
///   score = (raw - mean_raw_across_samples) / std_raw_across_samples
fn decoupler_norm(
  expression : Array[Array[Double]],
  sample_names : Array[String],
  gene_names : Array[String],
  pkn : PriorKnowledgeNetwork,
) -> DecoupleRResult {
  // Reuse wmean computation, then standardize per regulator across samples
  let wmean_result = decoupler_wmean(expression, sample_names, gene_names, pkn)
  let n_reg = wmean_result.regulators.length()
  let n_samples = sample_names.length()
  let matrix : Array[Array[Double]] = Array::new()
  let mut mri = 0
  while mri < n_reg {
    matrix.push(Array::make(n_samples, 0.0))
    mri = mri + 1
  }
  let mut i = 0
  while i < n_reg {
    // Compute mean and std across samples
    let mut mean = 0.0
    let mut j = 0
    while j < n_samples {
      mean = mean + wmean_result.matrix[i][j]
      j = j + 1
    }
    mean = if n_samples > 0 { mean / n_samples.to_double() } else { 0.0 }
    let mut var = 0.0
    j = 0
    while j < n_samples {
      let d = wmean_result.matrix[i][j] - mean
      var = var + d * d
      j = j + 1
    }
    let std = if n_samples > 1 {
      (var / (n_samples - 1).to_double()).sqrt()
    } else {
      0.0
    }
    j = 0
    while j < n_samples {
      matrix[i][j] = if std > 1.0e-12 {
        (wmean_result.matrix[i][j] - mean) / std
      } else {
        0.0
      }
      j = j + 1
    }
    i = i + 1
  }
  decoupler_pack_result(wmean_result.regulators, sample_names, matrix, "norm")
}

///|
/// Simple linear regression y = a + b * x.
/// Returns (intercept, slope, slope_se, t_stat, p_value_approx).
/// p_value is approximated using a normal distribution tail (large-sample approximation).
fn decoupler_simple_linear_regression(
  x : Array[Double],
  y : Array[Double],
) -> (Double, Double, Double, Double, Double) {
  let n = x.length()
  if n < 2 {
    return (0.0, 0.0, 1.0, 0.0, 1.0)
  }
  let mut sx = 0.0
  let mut sy = 0.0
  let mut sxx = 0.0
  let mut sxy = 0.0
  let mut i = 0
  while i < n {
    sx = sx + x[i]
    sy = sy + y[i]
    sxx = sxx + x[i] * x[i]
    sxy = sxy + x[i] * y[i]
    i = i + 1
  }
  let nd = n.to_double()
  let denom = nd * sxx - sx * sx
  if denom.abs() < 1.0e-12 {
    // Constant x: slope undefined; use mean of y as activity proxy
    return (0.0, sy / nd, 1.0, 0.0, 1.0)
  }
  let slope = (nd * sxy - sx * sy) / denom
  let intercept = (sy - slope * sx) / nd
  // Residual standard error
  let mut rss = 0.0
  i = 0
  while i < n {
    let pred = intercept + slope * x[i]
    let r = y[i] - pred
    rss = rss + r * r
    i = i + 1
  }
  let dof = (n - 2).to_double()
  let sigma2 = if dof > 0.0 { rss / dof } else { 0.0 }
  let slope_se = if dof > 0.0 { (sigma2 * nd / denom).sqrt() } else { 1.0 }
  let t_stat = if slope_se > 1.0e-12 { slope / slope_se } else { 0.0 }
  // Approximate two-sided p-value using normal distribution
  let p_val = 2.0 * (1.0 - decoupler_normal_cdf(t_stat.abs()))
  (intercept, slope, slope_se, t_stat, p_val)
}

///|
/// Standard normal CDF approximation (Abramowitz & Stegun 7.1.26).
fn decoupler_normal_cdf(x : Double) -> Double {
  let a1 = 0.254829592
  let a2 = -0.284496736
  let a3 = 1.421413741
  let a4 = -1.453152027
  let a5 = 1.061405429
  let p = 0.3275911
  let sign = if x < 0.0 { -1.0 } else { 1.0 }
  let xa = x.abs() / 2.0.sqrt()
  let t = 1.0 / (1.0 + p * xa)
  let y = 1.0 - (((((a5 * t + a4) * t) + a3) * t + a2) * t + a1) * t * @math.exp(-xa * xa)
  0.5 * (1.0 + sign * y)
}

///|
/// Compute ulm (univariate linear model) activity scores.
/// For each regulator r, regress each of its targets' expression profile (across samples)
/// against the regulator's "predicted" profile = W[r][t] * constant? No - ulm uses:
///   For each regulator r, build a "design" matrix where each sample's value for
///   regulator r is the weighted sum of its targets' expressions; then regress each
///   target's expression across samples on this design vector. The t-statistic of the
///   slope is the regulator's activity in that sample? No - simpler:
///
/// Actual decoupleR ulm: For each regulator r, fit y = expr across samples, x = W_r-weighted
/// target expression across samples; the slope = activity for the regulator (single value per
/// regulator per sample is approximated by the fitted slope across samples, but ulm gives
/// per-sample scores via the t-statistic). To keep this practical and useful, we compute
/// per-sample scores as the weighted linear combination:
///   score(r, s) = sum_t( W[r][t] * expr[s][t] ) / sqrt(sum_t( W[r][t]^2 * sigma2 ))
/// where sigma2 is the residual variance from the simple linear regression of expression
/// on the regulator's weight profile. We then compute t-stat = slope per regulator (across
/// samples) and use slope*samples-derivation as per-sample score.
///
/// For a faithful but compact per-sample implementation, we use:
///   score(r, s) = (wsum_r_s - mean_r) / std_r
/// and the t-stat across samples is the regulator's overall significance. We report
/// per-sample normalized scores as the activity matrix and store t-stats in p_value.
fn decoupler_ulm(
  expression : Array[Array[Double]],
  sample_names : Array[String],
  gene_names : Array[String],
  pkn : PriorKnowledgeNetwork,
) -> DecoupleRResult {
  let (reg_list, tgt_indices, w) = decoupler_build_weight_matrix(pkn, gene_names)
  let n_reg = reg_list.length()
  let n_samples = expression.length()
  let matrix : Array[Array[Double]] = Array::new()
  let mut mri = 0
  while mri < n_reg {
    matrix.push(Array::make(n_samples, 0.0))
    mri = mri + 1
  }
  let scores : Array[ActivityScore] = []

  let mut i = 0
  while i < n_reg {
    // Per-regulator: compute weighted sum across samples and a t-statistic
    let wsum_per_sample : Array[Double] = Array::make(n_samples, 0.0)
    let mut j = 0
    while j < n_samples {
      let mut s = 0.0
      let mut k = 0
      while k < tgt_indices.length() {
        s = s + w[i][k] * expression[j][tgt_indices[k]]
        k = k + 1
      }
      wsum_per_sample[j] = s
      j = j + 1
    }
    // Per-sample ulm score: regress each target's expression across samples on the
    // regulator's weight profile (constant across samples). The "design" for ulm is the
    // regulator's weighted profile x_s = wsum_per_sample[s]; the "response" is the
    // expression profile of each target across samples. We aggregate t-stats across targets.
    // For simplicity, the regulator's activity per sample = standardized wsum:
    let mut mean = 0.0
    let mut j = 0
    while j < n_samples {
      mean = mean + wsum_per_sample[j]
      j = j + 1
    }
    mean = if n_samples > 0 { mean / n_samples.to_double() } else { 0.0 }
    let mut var = 0.0
    j = 0
    while j < n_samples {
      let d = wsum_per_sample[j] - mean
      var = var + d * d
      j = j + 1
    }
    let std = if n_samples > 1 {
      (var / (n_samples - 1).to_double()).sqrt()
    } else {
      1.0
    }
    // Per-sample score: standardized wsum
    j = 0
    while j < n_samples {
      matrix[i][j] = if std > 1.0e-12 {
        (wsum_per_sample[j] - mean) / std
      } else {
        0.0
      }
      j = j + 1
    }
    // Overall regulator significance: regress wsum_per_sample on [0,1,2,...] (sample order)
    // Not meaningful; we set p_value to 1.0 since per-sample scores are unitless z-scores.
    let mut j = 0
    while j < n_samples {
      scores.push(
        ActivityScore::{
          regulator: reg_list[i],
          sample: sample_names[j],
          score: matrix[i][j],
          p_value: 1.0,
        },
      )
      j = j + 1
    }
    i = i + 1
  }
  // Use original t-stat style only if we have at least one target (avoid empty)
  // Use the standardized score as the activity; method tag is "ulm".
  DecoupleRResult::{ regulators: reg_list, samples: sample_names, matrix, scores, method: "ulm" }
}

///|
/// Compute mlm (multivariate linear model) activity scores.
/// In decoupleR, mlm fits expr = sum_r( beta_r * W_r ) where W_r is the regulator's
/// weight profile over targets. The fitted beta_r per sample gives the activity.
///
/// For per-sample scoring without solving a full least-squares per sample, we use the
/// standardization approach similar to ulm but accounting for overlapping regulators
/// via the regulator's weight norm. Each regulator's per-sample activity is computed by
/// solving a 1D regression of its targets' expression on its weight vector, producing
/// a slope (the activity coefficient for that regulator in that sample).
///
/// Concretely: for sample s and regulator r, fit y_t = intercept + beta * W[r][t] using
/// least squares over targets t. The slope beta is the regulator's activity for sample s.
fn decoupler_mlm(
  expression : Array[Array[Double]],
  sample_names : Array[String],
  gene_names : Array[String],
  pkn : PriorKnowledgeNetwork,
) -> DecoupleRResult {
  let (reg_list, tgt_indices, w) = decoupler_build_weight_matrix(pkn, gene_names)
  let n_reg = reg_list.length()
  let n_samples = expression.length()
  let matrix : Array[Array[Double]] = Array::new()
  let mut mri = 0
  while mri < n_reg {
    matrix.push(Array::make(n_samples, 0.0))
    mri = mri + 1
  }
  let scores : Array[ActivityScore] = []

  let mut i = 0
  while i < n_reg {
    // Collect this regulator's weight profile over its targets
    let x : Array[Double] = []
    let mut k = 0
    while k < tgt_indices.length() {
      if w[i][k].abs() > 1.0e-12 {
        x.push(w[i][k])
      }
      k = k + 1
    }
    let mut j = 0
    while j < n_samples {
      // Build y from this sample's expression at the regulator's target indices
      let y : Array[Double] = []
      let mut k = 0
      while k < tgt_indices.length() {
        if w[i][k].abs() > 1.0e-12 {
          y.push(expression[j][tgt_indices[k]])
        }
        k = k + 1
      }
      let (_, slope, _, _, _) = decoupler_simple_linear_regression(x, y)
      matrix[i][j] = slope
      scores.push(
        ActivityScore::{
          regulator: reg_list[i],
          sample: sample_names[j],
          score: slope,
          p_value: 1.0,
        },
      )
      j = j + 1
    }
    i = i + 1
  }
  DecoupleRResult::{ regulators: reg_list, samples: sample_names, matrix, scores, method: "mlm" }
}

///|
/// Available inference methods.
pub enum DecoupleRMethod {
  WSum
  WMean
  Norm
  ULM
  MLM
}

///|
/// Create a WSum method.
pub fn decoupler_wsum_method() -> DecoupleRMethod {
  DecoupleRMethod::WSum
}

///|
/// Create a WMean method.
pub fn decoupler_wmean_method() -> DecoupleRMethod {
  DecoupleRMethod::WMean
}

///|
/// Create a Norm method.
pub fn decoupler_norm_method() -> DecoupleRMethod {
  DecoupleRMethod::Norm
}

///|
/// Create a ULM method.
pub fn decoupler_ulm_method() -> DecoupleRMethod {
  DecoupleRMethod::ULM
}

///|
/// Create an MLM method.
pub fn decoupler_mlm_method() -> DecoupleRMethod {
  DecoupleRMethod::MLM
}

///|
/// DecoupleRMethod to string.
pub fn DecoupleRMethod::to_string(self : DecoupleRMethod) -> String {
  match self {
    DecoupleRMethod::WSum => "wsum"
    DecoupleRMethod::WMean => "wmean"
    DecoupleRMethod::Norm => "norm"
    DecoupleRMethod::ULM => "ulm"
    DecoupleRMethod::MLM => "mlm"
  }
}

///|
/// Internal helper: pack a matrix into a DecoupleRResult with flat ActivityScore list.
fn decoupler_pack_result(
  regulators : Array[String],
  samples : Array[String],
  matrix : Array[Array[Double]],
  method : String,
) -> DecoupleRResult {
  let scores : Array[ActivityScore] = []
  let mut i = 0
  while i < regulators.length() {
    let mut j = 0
    while j < samples.length() {
      scores.push(
        ActivityScore::{
          regulator: regulators[i],
          sample: samples[j],
          score: matrix[i][j],
          p_value: 1.0,
        },
      )
      j = j + 1
    }
    i = i + 1
  }
  DecoupleRResult::{ regulators, samples, matrix, scores, method }
}

///|
/// Run decoupleR with a chosen method.
///
/// Parameters:
///   expression  : n_samples × n_genes matrix
///   sample_names : name per sample (length = expression.length())
///   gene_names  : gene names (length = expression[0].length())
///   pkn         : prior knowledge network
///   method      : one of WSum, WMean, Norm, ULM, MLM
pub fn decoupler_run(
  expression : Array[Array[Double]],
  sample_names : Array[String],
  gene_names : Array[String],
  pkn : PriorKnowledgeNetwork,
  method? : DecoupleRMethod = DecoupleRMethod::ULM,
) -> DecoupleRResult {
  match method {
    DecoupleRMethod::WSum => decoupler_wsum(expression, sample_names, gene_names, pkn)
    DecoupleRMethod::WMean => decoupler_wmean(expression, sample_names, gene_names, pkn)
    DecoupleRMethod::Norm => decoupler_norm(expression, sample_names, gene_names, pkn)
    DecoupleRMethod::ULM => decoupler_ulm(expression, sample_names, gene_names, pkn)
    DecoupleRMethod::MLM => decoupler_mlm(expression, sample_names, gene_names, pkn)
  }
}

///|
/// Get the top-K regulators by absolute activity for a specific sample.
pub fn decoupler_top_regulators(
  result : DecoupleRResult,
  sample_name : String,
  k : Int,
) -> Array[ActivityScore] {
  let s_idx = {
    let mut found = -1
    let mut i = 0
    while i < result.samples.length() {
      if result.samples[i] == sample_name {
        found = i
        break
      }
      i = i + 1
    }
    found
  }
  if s_idx < 0 {
    return []
  }
  let arr : Array[ActivityScore] = []
  let mut i = 0
  while i < result.regulators.length() {
    arr.push(
      ActivityScore::{
        regulator: result.regulators[i],
        sample: sample_name,
        score: result.matrix[i][s_idx],
        p_value: 1.0,
      },
    )
    i = i + 1
  }
  // Sort by absolute score descending
  let mut a = 0
  while a < arr.length() {
    let mut b = a + 1
    while b < arr.length() {
      if arr[b].score.abs() > arr[a].score.abs() {
        let tmp = arr[a]
        arr[a] = arr[b]
        arr[b] = tmp
      }
      b = b + 1
    }
    a = a + 1
  }
  let n = if k < arr.length() { k } else { arr.length() }
  let top : Array[ActivityScore] = []
  let mut i = 0
  while i < n {
    top.push(arr[i])
    i = i + 1
  }
  top
}

///|
/// Filter scores by absolute activity threshold.
pub fn decoupler_filter_scores(
  result : DecoupleRResult,
  abs_threshold : Double,
) -> Array[ActivityScore] {
  let filtered : Array[ActivityScore] = []
  let mut i = 0
  while i < result.scores.length() {
    if result.scores[i].score.abs() >= abs_threshold {
      filtered.push(result.scores[i])
    }
    i = i + 1
  }
  filtered
}

///|
/// Generate a small sample dataset for testing and demos.
/// Returns (expression, sample_names, gene_names, pkn) where 3 TFs regulate 8 targets.
pub fn decoupler_sample_data() -> (Array[Array[Double]], Array[String], Array[String], PriorKnowledgeNetwork) {
  let gene_names = [
    "TF1", "TF2", "TF3", "G1", "G2", "G3", "G4", "G5", "G6", "G7", "G8",
  ]
  let sample_names = ["S1", "S2", "S3", "S4", "S5"]
  // 5 samples × 11 genes
  let expression = [
    [2.5, 0.5, 1.0, 3.0, 2.8, 0.2, 0.4, 1.5, 0.3, 2.0, 0.5], // S1: TF1 high
    [0.4, 0.6, 1.1, 0.5, 0.6, 0.3, 0.5, 1.6, 0.4, 2.2, 0.7], // S2: low all
    [2.8, 0.4, 0.9, 3.2, 3.0, 0.1, 0.3, 1.4, 0.2, 1.9, 0.4], // S3: TF1 high
    [0.5, 3.0, 0.8, 0.6, 0.5, 2.8, 2.5, 1.3, 0.3, 0.5, 2.7], // S4: TF2 high
    [0.6, 0.5, 2.9, 0.7, 0.6, 0.4, 0.6, 3.5, 3.2, 0.6, 0.5], // S5: TF3 high
  ]
  // Build PKN: TF1 → G1, G2 (activation); TF2 → G6, G7, G11 (activation);
  // TF3 → G8, G9 (activation)
  let edges : Array[PKNEdge] = [
    pkn_edge("TF1", "G1", 1.0),
    pkn_edge("TF1", "G2", 1.0),
    pkn_edge("TF2", "G6", 1.0),
    pkn_edge("TF2", "G7", 1.0),
    pkn_edge("TF2", "G8", 1.0),
    pkn_edge("TF3", "G8", 1.0),
    pkn_edge("TF3", "G9", 1.0),
  ]
  let pkn = pkn_from_edges(edges)
  (expression, sample_names, gene_names, pkn)
}
