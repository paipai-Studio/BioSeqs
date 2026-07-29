///|
/// MOFA2 - Multi-Omics Factor Analysis v2.
/// Inspired by the MOFA2 R package (Velten et al., 2021).
/// Implements sparse Bayesian factor analysis for multi-omics data integration
/// using an Expectation-Maximization (EM) algorithm with automatic relevance
/// determination (ARD) to infer the optimal number of latent factors.

///|
/// A single omics view.
pub struct MofaView {
  pub view_name : String
  pub data : Array[Array[Double]]
  pub feature_names : Array[String]
  pub sample_names : Array[String]
  pub view_type : String
}

///|
/// Get view name.
pub fn MofaView::get_view_name(self : MofaView) -> String {
  self.view_name
}

///|
/// Get data matrix.
pub fn MofaView::get_data(self : MofaView) -> Array[Array[Double]] {
  self.data
}

///|
/// Get feature names.
pub fn MofaView::get_feature_names(self : MofaView) -> Array[String] {
  self.feature_names
}

///|
/// Get number of features.
pub fn MofaView::n_features(self : MofaView) -> Int {
  self.feature_names.length()
}

///|
/// Get view type.
pub fn MofaView::get_view_type(self : MofaView) -> String {
  self.view_type
}

///|
/// Create a MOFA view.
pub fn MofaView::new(
  view_name : String,
  data : Array[Array[Double]],
  feature_names : Array[String],
  sample_names : Array[String],
  view_type? : String = "gaussian",
) -> MofaView {
  MofaView::{ view_name, data, feature_names, sample_names, view_type }
}

///|
/// MOFA2 hyperparameters.
pub struct MofaParams {
  n_factors : Int
  max_iterations : Int
  tolerance : Double
  ard_threshold : Double
  seed : Int
}

///|
/// Create default MOFA2 parameters.
pub fn MofaParams::new() -> MofaParams {
  MofaParams::{
    n_factors: 10,
    max_iterations: 100,
    tolerance: 1.0e-6,
    ard_threshold: 0.01,
    seed: 42,
  }
}

///|
/// Create MOFA2 parameters with custom values.
pub fn MofaParams::create(
  n_factors? : Int = 10,
  max_iterations? : Int = 100,
  tolerance? : Double = 1.0e-6,
  ard_threshold? : Double = 0.01,
  seed? : Int = 42,
) -> MofaParams {
  MofaParams::{
    n_factors: n_factors,
    max_iterations: max_iterations,
    tolerance: tolerance,
    ard_threshold: ard_threshold,
    seed: seed,
  }
}

///|
/// Result of MOFA2 factor analysis.
pub struct MofaResult {
  pub factors : Array[Array[Double]]
  pub loadings : Array[Array[Array[Double]]]
  pub factor_vars : Array[Double]
  pub view_vars : Array[Array[Double]]
  sigma : Array[Array[Double]]
  pub n_iterations : Int
  pub converged : Bool
  pub active_factors : Int
}

///|
/// Get factor matrix (samples x factors).
pub fn MofaResult::get_factors(self : MofaResult) -> Array[Array[Double]] {
  self.factors
}

///|
/// Get loading matrices (views x features x factors).
pub fn MofaResult::get_loadings(self : MofaResult) -> Array[Array[Array[Double]]] {
  self.loadings
}

///|
/// Get factor variance explained.
pub fn MofaResult::get_factor_vars(self : MofaResult) -> Array[Double] {
  self.factor_vars
}

///|
/// Get per-view variance explained.
pub fn MofaResult::get_view_vars(self : MofaResult) -> Array[Array[Double]] {
  self.view_vars
}

///|
/// Get number of iterations.
pub fn MofaResult::get_n_iterations(self : MofaResult) -> Int {
  self.n_iterations
}

///|
/// Check if algorithm converged.
pub fn MofaResult::is_converged(self : MofaResult) -> Bool {
  self.converged
}

///|
/// Get number of active factors.
pub fn MofaResult::get_active_factors(self : MofaResult) -> Int {
  self.active_factors
}

///|
/// Run MOFA2 factor analysis on multiple views.
pub fn mofa_run(views : Array[MofaView], params : MofaParams) -> MofaResult {
  let n_views = views.length()
  if n_views == 0 {
    return MofaResult::{
      factors: [],
      loadings: [],
      factor_vars: [],
      view_vars: [],
      sigma: [],
      n_iterations: 0,
      converged: false,
      active_factors: 0,
    }
  }
  
  let n_samples = views[0].data.length()
  let n_factors = params.n_factors
  
  // Center each view
  let centered_views : Array[Array[Array[Double]]] = Array::new()
  {
    let mut v = 0
    while v < n_views {
      let view_data = views[v].data
      let n_features = if view_data.length() > 0 { view_data[0].length() } else { 0 }
      let centered : Array[Array[Double]] = Array::new()
      
      // Compute feature means
      let means : Array[Double] = Array::make(n_features, 0.0)
      let mut s = 0
      while s < n_samples {
        let mut f = 0
        while f < n_features {
          means[f] = means[f] + view_data[s][f]
          f = f + 1
        }
        s = s + 1
      }
      {
        let mut f = 0
        while f < n_features {
          means[f] = means[f] / n_samples.to_double()
          f = f + 1
        }
      }
      
      // Center data
      s = 0
      while s < n_samples {
        let row : Array[Double] = Array::make(n_features, 0.0)
        let mut f = 0
        while f < n_features {
          row[f] = view_data[s][f] - means[f]
          f = f + 1
        }
        centered.push(row)
        s = s + 1
      }
      centered_views.push(centered)
      v = v + 1
    }
  }
  
  // Initialize factors Z (samples x factors)
  let mut z = mofa_init_factors(n_samples, n_factors, params.seed)
  
  // Initialize loadings W (views x features x factors) with small random values
  let loadings : Array[Array[Array[Double]]] = Array::new()
  {
    let mut v = 0
    while v < n_views {
      let d = if centered_views[v].length() > 0 {
        centered_views[v][0].length()
      } else {
        0
      }
      loadings.push(mofa_init_loadings(d, n_factors, params.seed + (v + 1) * 7919))
      v = v + 1
    }
  }
  
  // Initialize sigma (noise variance per view, feature)
  let sigma : Array[Array[Double]] = Array::new()
  {
    let mut v = 0
    while v < n_views {
      let d = if centered_views[v].length() > 0 {
        centered_views[v][0].length()
      } else {
        0
      }
      let view_sigma : Array[Double] = Array::make(d, 1.0)
      sigma.push(view_sigma)
      v = v + 1
    }
  }
  
  // ARD precision parameters
  let mut alpha : Array[Double] = Array::make(n_factors, 1.0)
  let mu = 0.0
  let tau = 1.0
  
  // EM iterations
  let mut iter = 0
  let mut converged = false
  let mut prev_elbo = -1.0e30
  
  while iter < params.max_iterations {
    // E-step: Update factors Z
    z = mofa_update_factors(centered_views, loadings, sigma, z, n_factors, n_views, n_samples)
    
    // M-step: Update loadings W
    {
      let mut v = 0
      while v < n_views {
        let n_features = if centered_views[v].length() > 0 {
          centered_views[v][0].length()
        } else {
          0
        }
        loadings[v] = mofa_update_loadings(
          centered_views[v], z, loadings[v], sigma[v],
          n_features, n_factors, n_samples
        )
        v = v + 1
      }
    }
    
    // M-step: Update sigma
    {
      let mut v = 0
      while v < n_views {
        let n_features = if centered_views[v].length() > 0 {
          centered_views[v][0].length()
        } else {
          0
        }
        sigma[v] = mofa_update_sigma(
          centered_views[v], z, loadings[v], sigma[v],
          n_features, n_factors, n_samples
        )
        v = v + 1
      }
    }
    
    // M-step: Update ARD alpha
    alpha = mofa_update_alpha(loadings, n_factors, n_views, mu, tau)
    
    // Compute ELBO (approximate)
    let elbo = mofa_compute_elbo(centered_views, z, loadings, sigma, alpha, n_factors, n_views, n_samples)
    let delta = (elbo - prev_elbo).abs()
    if delta < params.tolerance {
      converged = true
      break
    }
    prev_elbo = elbo
    iter = iter + 1
  }
  
  // Compute variance explained
  let (factor_vars, view_vars) = mofa_compute_variance_explained(
    centered_views, z, loadings, n_factors, n_views, n_samples
  )
  
  // Count active factors
  let active = mofa_count_active_factors(factor_vars, params.ard_threshold)
  
  MofaResult::{
    factors: z,
    loadings,
    factor_vars,
    view_vars,
    sigma,
    n_iterations: iter,
    converged,
    active_factors: active,
  }
}

///|
/// Initialize factor matrix with small random values.
fn mofa_init_factors(n_samples : Int, n_factors : Int, seed : Int) -> Array[Array[Double]] {
  let z : Array[Array[Double]] = Array::new()
  let mut s = seed
  let mut i = 0
  while i < n_samples {
    let row : Array[Double] = Array::new()
    let mut k = 0
    while k < n_factors {
      s = mofa_lcg(s)
      let r = (s.to_double() / 2147483647.0 - 0.5) * 0.5
      row.push(r)
      k = k + 1
    }
    z.push(row)
    i = i + 1
  }
  z
}

///|
/// Initialize loading matrix with small random values.
fn mofa_init_loadings(n_features : Int, n_factors : Int, seed : Int) -> Array[Array[Double]] {
  let w : Array[Array[Double]] = Array::new()
  let mut s = seed
  let mut f = 0
  while f < n_features {
    let row : Array[Double] = Array::new()
    let mut k = 0
    while k < n_factors {
      s = mofa_lcg(s)
      let r = (s.to_double() / 2147483647.0 - 0.5) * 0.5
      row.push(r)
      k = k + 1
    }
    w.push(row)
    f = f + 1
  }
  w
}

///|
/// Linear congruential generator.
fn mofa_lcg(seed : Int) -> Int {
  let state = seed * 1103515245 + 12345
  state & 0x7FFFFFFF
}

///|
/// E-step: Update factors Z.
fn mofa_update_factors(
  views : Array[Array[Array[Double]]],
  loadings : Array[Array[Array[Double]]],
  sigma : Array[Array[Double]],
  z : Array[Array[Double]],
  n_factors : Int,
  n_views : Int,
  n_samples : Int,
) -> Array[Array[Double]] {
  let new_z : Array[Array[Double]] = Array::new()
  let mut i = 0
  while i < n_samples {
    let row : Array[Double] = Array::make(n_factors, 0.0)
    let mut k = 0
    while k < n_factors {
      // Accumulate information from all views
      let mut precision = 0.0
      let mut mean_accum = 0.0
      let mut v = 0
      while v < n_views {
        let n_features = if views[v].length() > 0 { views[v][0].length() } else { 0 }
        let mut f = 0
        while f < n_features {
          let w_fk = loadings[v][f][k]
          let s2 = if sigma[v][f] > 1.0e-10 { sigma[v][f] } else { 1.0e-10 }
          precision = precision + w_fk * w_fk / s2
          let residual = views[v][i][f] - mofa_prediction_single(loadings, z, i, f, k, v)
          mean_accum = mean_accum + w_fk * residual / s2
          f = f + 1
        }
        v = v + 1
      }
      // Regularize toward zero
      precision = precision + 0.01
      row[k] = mean_accum / precision
      k = k + 1
    }
    new_z.push(row)
    i = i + 1
  }
  new_z
}

fn mofa_prediction_single(
  loadings : Array[Array[Array[Double]]],
  z : Array[Array[Double]],
  sample_idx : Int,
  feature_idx : Int,
  excluded_factor : Int,
  view_idx : Int,
) -> Double {
  let n_factors = z[0].length()
  let mut pred = 0.0
  let mut k = 0
  while k < n_factors {
    if k != excluded_factor {
      pred = pred + loadings[view_idx][feature_idx][k] * z[sample_idx][k]
    }
    k = k + 1
  }
  pred
}

///|
/// M-step: Update loadings W.
fn mofa_update_loadings(
  view_data : Array[Array[Double]],
  z : Array[Array[Double]],
  w : Array[Array[Double]],
  sigma_view : Array[Double],
  n_features : Int,
  n_factors : Int,
  n_samples : Int,
) -> Array[Array[Double]] {
  let new_w : Array[Array[Double]] = Array::new()
  let mut f = 0
  while f < n_features {
    let row : Array[Double] = Array::make(n_factors, 0.0)
    let mut k = 0
    while k < n_factors {
      let s2 = if sigma_view[f] > 1.0e-10 { sigma_view[f] } else { 1.0e-10 }
      let mut numerator = 0.0
      let mut denominator = 0.0
      let mut i = 0
      while i < n_samples {
        let residual = view_data[i][f] - mofa_prediction_loadings(w, z, i, f, k)
        numerator = numerator + z[i][k] * residual
        denominator = denominator + z[i][k] * z[i][k]
        i = i + 1
      }
      denominator = denominator + 0.01 * s2
      row[k] = if denominator > 1.0e-10 { numerator / denominator } else { w[f][k] }
      k = k + 1
    }
    new_w.push(row)
    f = f + 1
  }
  new_w
}

fn mofa_prediction_loadings(
  w : Array[Array[Double]],
  z : Array[Array[Double]],
  sample_idx : Int,
  feature_idx : Int,
  excluded_factor : Int,
) -> Double {
  let n_factors = w[0].length()
  let mut pred = 0.0
  let mut k = 0
  while k < n_factors {
    if k != excluded_factor {
      pred = pred + w[feature_idx][k] * z[sample_idx][k]
    }
    k = k + 1
  }
  pred
}

///|
/// M-step: Update sigma (noise variance).
fn mofa_update_sigma(
  view_data : Array[Array[Double]],
  z : Array[Array[Double]],
  w : Array[Array[Double]],
  sigma_view : Array[Double],
  n_features : Int,
  n_factors : Int,
  n_samples : Int,
) -> Array[Double] {
  let new_sigma : Array[Double] = Array::make(n_features, 0.0)
  let mut f = 0
  while f < n_features {
    let mut ss = 0.0
    let mut i = 0
    while i < n_samples {
      let mut pred = 0.0
      let mut k = 0
      while k < n_factors {
        pred = pred + w[f][k] * z[i][k]
        k = k + 1
      }
      let residual = view_data[i][f] - pred
      ss = ss + residual * residual
      i = i + 1
    }
    new_sigma[f] = if n_samples > 0 { ss / n_samples.to_double() } else { sigma_view[f] }
    if new_sigma[f] < 1.0e-6 {
      new_sigma[f] = 1.0e-6
    }
    f = f + 1
  }
  new_sigma
}

///|
/// M-step: Update ARD precision parameters.
fn mofa_update_alpha(
  loadings : Array[Array[Array[Double]]],
  n_factors : Int,
  n_views : Int,
  mu : Double,
  tau : Double,
) -> Array[Double] {
  let new_alpha : Array[Double] = Array::make(n_factors, 0.0)
  let mut k = 0
  while k < n_factors {
    let mut sum_sq = 0.0
    let mut total_features = 0
    let mut v = 0
    while v < n_views {
      let n_features = loadings[v].length()
      let mut f = 0
      while f < n_features {
        sum_sq = sum_sq + loadings[v][f][k] * loadings[v][f][k]
        f = f + 1
      }
      total_features = total_features + n_features
      v = v + 1
    }
    let expected = if total_features > 0 { sum_sq / total_features.to_double() } else { 0.0 }
    let prec = if expected > 1.0e-10 { 1.0 / expected } else { 1.0e6 }
    new_alpha[k] = prec
    k = k + 1
  }
  new_alpha
}

///|
/// Compute approximate ELBO.
fn mofa_compute_elbo(
  views : Array[Array[Array[Double]]],
  z : Array[Array[Double]],
  loadings : Array[Array[Array[Double]]],
  sigma : Array[Array[Double]],
  alpha : Array[Double],
  n_factors : Int,
  n_views : Int,
  n_samples : Int,
) -> Double {
  let mut elbo = 0.0
  let mut v = 0
  while v < n_views {
    let n_features = if views[v].length() > 0 { views[v][0].length() } else { 0 }
    let mut i = 0
    while i < n_samples {
      let mut f = 0
      while f < n_features {
        let mut pred = 0.0
        let mut k = 0
        while k < n_factors {
          pred = pred + loadings[v][f][k] * z[i][k]
          k = k + 1
        }
        let s2 = if sigma[v][f] > 1.0e-10 { sigma[v][f] } else { 1.0e-10 }
        let residual = views[v][i][f] - pred
        elbo = elbo - 0.5 * @math.ln(2.0 * 3.14159265 * s2) - residual * residual / (2.0 * s2)
        f = f + 1
      }
      i = i + 1
    }
    v = v + 1
  }
  elbo
}

///|
/// Compute variance explained by each factor per view.
fn mofa_compute_variance_explained(
  views : Array[Array[Array[Double]]],
  z : Array[Array[Double]],
  loadings : Array[Array[Array[Double]]],
  n_factors : Int,
  n_views : Int,
  n_samples : Int,
) -> (Array[Double], Array[Array[Double]]) {
  // Factor variance: variance of each factor across samples
  let factor_vars : Array[Double] = Array::make(n_factors, 0.0)
  {
    let mut k = 0
    while k < n_factors {
      let mut sum = 0.0
      let mut sum_sq = 0.0
      let mut i = 0
      while i < n_samples {
        sum = sum + z[i][k]
        sum_sq = sum_sq + z[i][k] * z[i][k]
        i = i + 1
      }
      let mean = sum / n_samples.to_double()
      let var = if n_samples > 1 { (sum_sq / n_samples.to_double()) - mean * mean } else { 0.0 }
      factor_vars[k] = if var > 0.0 { var } else { 0.0 }
      k = k + 1
    }
  }
  
  // View variance explained by each factor
  let view_vars : Array[Array[Double]] = Array::make(n_views, Array::make(n_factors, 0.0))
  {
    let mut v = 0
    while v < n_views {
      let n_features = if views[v].length() > 0 { views[v][0].length() } else { 0 }
      let mut k = 0
      while k < n_factors {
        let mut total_var = 0.0
        let mut f = 0
        while f < n_features {
          let loading = loadings[v][f][k]
          total_var = total_var + loading * loading * factor_vars[k]
          f = f + 1
        }
        view_vars[v][k] = total_var
        k = k + 1
      }
      v = v + 1
    }
  }
  
  // Normalize view_vars to percentages
  {
    let mut v = 0
    while v < n_views {
      let total = mofa_sum_double(view_vars[v])
      if total > 0.0 {
        let mut k = 0
        while k < n_factors {
          view_vars[v][k] = view_vars[v][k] / total * 100.0
          k = k + 1
        }
      }
      v = v + 1
    }
  }
  
  (factor_vars, view_vars)
}

///|
/// Count active factors based on variance threshold.
fn mofa_count_active_factors(factor_vars : Array[Double], threshold : Double) -> Int {
  let mut count = 0
  for v in factor_vars {
    if v > threshold {
      count = count + 1
    }
  }
  count
}

fn mofa_sum_double(arr : Array[Double]) -> Double {
  let mut s = 0.0
  for v in arr {
    s = s + v
  }
  s
}

///|
/// Create example data for MOFA2 testing.
pub fn mofa_create_example(
  n_samples : Int,
  n_factors : Int,
) -> Array[MofaView] {
  let mut rng = 42
  let views : Array[MofaView] = Array::new()
  
  // View 1: RNA-seq (high-dimensional)
  {
    let n_features = 50
    let data : Array[Array[Double]] = Array::new()
    let mut s = 0
    while s < n_samples {
      let row : Array[Double] = Array::make(n_features, 0.0)
      let mut f = 0
      while f < n_features {
        let mut val = 0.0
        let mut k = 0
        while k < n_factors {
          rng = mofa_lcg(rng)
          let loading = (rng.to_double() / 2147483647.0 - 0.5) * 2.0
          rng = mofa_lcg(rng)
          let factor_val = (rng.to_double() / 2147483647.0 - 0.5) * 2.0
          val = val + loading * factor_val
          k = k + 1
        }
        rng = mofa_lcg(rng)
        let noise = (rng.to_double() / 2147483647.0 - 0.5) * 0.5
        row[f] = val + noise
        f = f + 1
      }
      data.push(row)
      s = s + 1
    }
    let feature_names : Array[String] = Array::new()
    let mut f = 0
    while f < n_features {
      feature_names.push("gene_" + f.to_string())
      f = f + 1
    }
    views.push(MofaView::new("rna", data, feature_names, [], view_type="gaussian"))
  }
  
  // View 2: ATAC-seq (lower-dimensional, shares some factors)
  {
    let n_features = 30
    let data : Array[Array[Double]] = Array::new()
    let mut s = 0
    while s < n_samples {
      let row : Array[Double] = Array::make(n_features, 0.0)
      let mut f = 0
      while f < n_features {
        let mut val = 0.0
        let mut k = 0
        while k < n_factors {
          rng = mofa_lcg(rng)
          let loading = (rng.to_double() / 2147483647.0 - 0.5) * 1.5
          rng = mofa_lcg(rng)
          let factor_val = (rng.to_double() / 2147483647.0 - 0.5) * 2.0
          val = val + loading * factor_val
          k = k + 1
        }
        rng = mofa_lcg(rng)
        let noise = (rng.to_double() / 2147483647.0 - 0.5) * 0.3
        row[f] = val + noise
        f = f + 1
      }
      data.push(row)
      s = s + 1
    }
    let feature_names : Array[String] = Array::new()
    let mut f = 0
    while f < n_features {
      feature_names.push("peak_" + f.to_string())
      f = f + 1
    }
    views.push(MofaView::new("atac", data, feature_names, [], view_type="gaussian"))
  }
  
  // Set sample names
  {
    let mut v = 0
    while v < views.length() {
      let sample_names : Array[String] = Array::new()
      let mut s = 0
      while s < n_samples {
        sample_names.push("sample_" + s.to_string())
        s = s + 1
      }
      views[v] = MofaView::new(views[v].view_name, views[v].data, views[v].feature_names, sample_names, view_type=views[v].view_type)
      v = v + 1
    }
  }
  
  views
}
