///|
/// ProDy-style ENM parameter fitting.
///
/// Fit ANM/GNM parameters (cutoff, gamma) to experimental or MD data:
///   - Optimal gamma scaling: linear least-squares
///   - Cutoff optimization: grid search with B-factor correlation
///   - Joint (cutoff, gamma) optimization
///   - Mode matching: overlap between PCA and ANM modes
///   - Fit quality metrics: Pearson r, RMS error
///   - Cross-validation: leave-one-out

///|
/// Result of an ENM parameter fit.
pub struct ENMFitResult {
  cutoff : Double
  gamma : Double
  correlation : Double
  rms_error : Double
  predicted_bfactors : Array[Double]
}

///|
/// Fit the optimal gamma scaling factor.
///
/// Given predicted MSF (at gamma=1) and observed B-factors,
/// gamma_opt = sum(B_obs * B_pred) / sum(B_pred^2)
/// This minimizes sum((gamma * B_pred - B_obs)^2).
pub fn fit_gamma(
  predicted_msf : Array[Double],
  observed_bfactors : Array[Double],
) -> Double {
  let n = if predicted_msf.length() < observed_bfactors.length() {
    predicted_msf.length()
  } else {
    observed_bfactors.length()
  }
  if n == 0 {
    return 1.0
  }
  let mut num = 0.0
  let mut den = 0.0
  for i = 0; i < n; i = i + 1 {
    num = num + observed_bfactors[i] * predicted_msf[i]
    den = den + predicted_msf[i] * predicted_msf[i]
  }
  if den < 0.00000000000000000001 {
    return 1.0
  }
  num / den
}

///|
/// Fit quality: (Pearson correlation, RMS error).
pub fn fit_quality(
  predicted : Array[Double],
  observed : Array[Double],
) -> (Double, Double) {
  let n = if predicted.length() < observed.length() {
    predicted.length()
  } else {
    observed.length()
  }
  if n == 0 {
    return (0.0, 0.0)
  }
  let corr = stat_pearson_correlation(predicted, observed)
  let mut sum_sq = 0.0
  for i = 0; i < n; i = i + 1 {
    let diff = predicted[i] - observed[i]
    sum_sq = sum_sq + diff * diff
  }
  let rms = (sum_sq / n.to_double()).sqrt()
  (corr, rms)
}

///|
/// Scale predicted B-factors by gamma.
pub fn scale_bfactors(
  predicted : Array[Double],
  gamma : Double,
) -> Array[Double] {
  let result : Array[Double] = Array(capacity=0)
  for v in predicted {
    result.push(v * gamma)
  }
  result
}

///|
/// Optimize cutoff by grid search using ANM B-factor correlation.
///
/// For each cutoff, build ANM, predict B-factors, fit gamma,
/// and compute Pearson correlation with observed B-factors.
/// Returns the best result.
pub fn fit_cutoff(
  chain : Chain,
  observed_bfactors : Array[Double],
  cutoffs : Array[Double],
) -> ENMFitResult? {
  if cutoffs.length() == 0 {
    return None
  }
  let mut best : ENMFitResult? = None
  let mut best_corr = -2.0
  for ci = 0; ci < cutoffs.length(); ci = ci + 1 {
    let cutoff = cutoffs[ci]
    match ANM::new(chain, cutoff~) {
      Some(anm) => {
        let pred = ANM::predict_bfactors(anm)
        let gamma = fit_gamma(pred, observed_bfactors)
        let scaled = scale_bfactors(pred, gamma)
        let (corr, rms) = fit_quality(scaled, observed_bfactors)
        if corr > best_corr {
          best_corr = corr
          best = Some(ENMFitResult::{
            cutoff,
            gamma,
            correlation: corr,
            rms_error: rms,
            predicted_bfactors: scaled,
          })
        }
      }
      None => ()
    }
  }
  best
}

///|
/// Optimize both cutoff and gamma jointly.
///
/// For each cutoff, fit gamma analytically (least-squares),
/// then pick the (cutoff, gamma) with best correlation.
/// This is equivalent to fit_cutoff but returns all candidates.
pub fn fit_cutoff_gamma(
  chain : Chain,
  observed_bfactors : Array[Double],
  cutoffs : Array[Double],
) -> Array[ENMFitResult] {
  let results : Array[ENMFitResult] = Array(capacity=0)
  for ci = 0; ci < cutoffs.length(); ci = ci + 1 {
    let cutoff = cutoffs[ci]
    match ANM::new(chain, cutoff~) {
      Some(anm) => {
        let pred = ANM::predict_bfactors(anm)
        let gamma = fit_gamma(pred, observed_bfactors)
        let scaled = scale_bfactors(pred, gamma)
        let (corr, rms) = fit_quality(scaled, observed_bfactors)
        results.push(ENMFitResult::{
          cutoff,
          gamma,
          correlation: corr,
          rms_error: rms,
          predicted_bfactors: scaled,
        })
      }
      None => ()
    }
  }
  results
}

///|
/// Select the best fit from a list of results.
pub fn best_fit(results : Array[ENMFitResult]) -> ENMFitResult? {
  if results.length() == 0 {
    return None
  }
  let mut best_idx = 0
  let mut best_corr = results[0].correlation
  for i = 1; i < results.length(); i = i + 1 {
    if results[i].correlation > best_corr {
      best_corr = results[i].correlation
      best_idx = i
    }
  }
  Some(results[best_idx])
}

///|
/// Fit using GNM instead of ANM.
pub fn fit_cutoff_gnm(
  chain : Chain,
  observed_bfactors : Array[Double],
  cutoffs : Array[Double],
) -> ENMFitResult? {
  if cutoffs.length() == 0 {
    return None
  }
  let mut best : ENMFitResult? = None
  let mut best_corr = -2.0
  for ci = 0; ci < cutoffs.length(); ci = ci + 1 {
    let cutoff = cutoffs[ci]
    match GNM::new(chain, cutoff~) {
      Some(gnm) => {
        let pred = GNM::predict_bfactors(gnm)
        let gamma = fit_gamma(pred, observed_bfactors)
        let scaled = scale_bfactors(pred, gamma)
        let (corr, rms) = fit_quality(scaled, observed_bfactors)
        if corr > best_corr {
          best_corr = corr
          best = Some(ENMFitResult::{
            cutoff,
            gamma,
            correlation: corr,
            rms_error: rms,
            predicted_bfactors: scaled,
          })
        }
      }
      None => ()
    }
  }
  best
}

// ===== Mode matching =====

///|
/// Compute overlap matrix between two sets of modes.
///
/// overlap[i][j] = |<mode_a[i], mode_b[j]>|
/// Each mode should be normalized.
pub fn match_modes(
  modes_a : Array[Array[Double]],
  modes_b : Array[Array[Double]],
) -> Array[Array[Double]] {
  let na = modes_a.length()
  let nb = modes_b.length()
  let overlap : Array[Array[Double]] = Array::makei(na, { _ => Array::make(nb, 0.0) })
  for i = 0; i < na; i = i + 1 {
    for j = 0; j < nb; j = j + 1 {
      let mut dot = 0.0
      let len = if modes_a[i].length() < modes_b[j].length() {
        modes_a[i].length()
      } else {
        modes_b[j].length()
      }
      for k = 0; k < len; k = k + 1 {
        dot = dot + modes_a[i][k] * modes_b[j][k]
      }
      overlap[i][j] = dot.abs()
    }
  }
  overlap
}

///|
/// Effective overlap: sqrt(sum of squared max overlaps).
/// Measures how well mode set A is represented by mode set B.
pub fn enm_effective_overlap(overlap_matrix : Array[Array[Double]]) -> Double {
  let na = overlap_matrix.length()
  if na == 0 {
    return 0.0
  }
  let mut sum_sq = 0.0
  for i = 0; i < na; i = i + 1 {
    let mut max_val = 0.0
    for j = 0; j < overlap_matrix[i].length(); j = j + 1 {
      if overlap_matrix[i][j] > max_val {
        max_val = overlap_matrix[i][j]
      }
    }
    sum_sq = sum_sq + max_val * max_val
  }
  (sum_sq / na.to_double()).sqrt()
}

///|
/// Mode matching score: average of max overlaps.
pub fn mode_matching_score(overlap_matrix : Array[Array[Double]]) -> Double {
  let na = overlap_matrix.length()
  if na == 0 {
    return 0.0
  }
  let mut total = 0.0
  for i = 0; i < na; i = i + 1 {
    let mut max_val = 0.0
    for j = 0; j < overlap_matrix[i].length(); j = j + 1 {
      if overlap_matrix[i][j] > max_val {
        max_val = overlap_matrix[i][j]
      }
    }
    total = total + max_val
  }
  total / na.to_double()
}

///|
/// Cumulative overlap: for each mode in A, the norm of its
/// projection onto the subspace spanned by modes in B.
pub fn enm_cumulative_overlap(
  modes_a : Array[Array[Double]],
  modes_b : Array[Array[Double]],
) -> Array[Double] {
  let result : Array[Double] = Array(capacity=0)
  for a in modes_a {
    // Project mode_a onto span(modes_b)
    let mut proj_sq = 0.0
    for b in modes_b {
      let mut dot = 0.0
      let len = if a.length() < b.length() {
        a.length()
      } else {
        b.length()
      }
      for k = 0; k < len; k = k + 1 {
        dot = dot + a[k] * b[k]
      }
      proj_sq = proj_sq + dot * dot
    }
    result.push(proj_sq.sqrt())
  }
  result
}

// ===== Cross-validation =====

///|
/// Leave-one-out cross-validation of B-factor prediction.
///
/// For each residue i, fit gamma using all residues except i,
/// then predict B-factor(i) and compare to observed.
/// Returns residuals (predicted - observed) for each residue.
pub fn loo_cross_validate(
  predicted_msf : Array[Double],
  observed_bfactors : Array[Double],
) -> Array[Double] {
  let n = if predicted_msf.length() < observed_bfactors.length() {
    predicted_msf.length()
  } else {
    observed_bfactors.length()
  }
  let residuals : Array[Double] = Array::make(n, 0.0)
  for i = 0; i < n; i = i + 1 {
    // Fit gamma without residue i
    let mut num = 0.0
    let mut den = 0.0
    for j = 0; j < n; j = j + 1 {
      if j == i {
        continue
      }
      num = num + observed_bfactors[j] * predicted_msf[j]
      den = den + predicted_msf[j] * predicted_msf[j]
    }
    let gamma = if den < 0.00000000000000000001 {
      1.0
    } else {
      num / den
    }
    // Predict B-factor(i)
    let predicted = gamma * predicted_msf[i]
    residuals[i] = predicted - observed_bfactors[i]
  }
  residuals
}

///|
/// k-fold cross-validation.
/// Splits residues into k folds, fits on k-1 folds, predicts the held-out fold.
/// Returns (correlation, rms_error) over all predictions.
pub fn kfold_cross_validate(
  predicted_msf : Array[Double],
  observed_bfactors : Array[Double],
  k : Int,
) -> (Double, Double) {
  let n = if predicted_msf.length() < observed_bfactors.length() {
    predicted_msf.length()
  } else {
    observed_bfactors.length()
  }
  if n == 0 || k <= 0 {
    return (0.0, 0.0)
  }
  let fold_size = n / k
  let all_pred : Array[Double] = Array::make(n, 0.0)
  for fold = 0; fold < k; fold = fold + 1 {
    let start = fold * fold_size
    let end = if fold == k - 1 {
      n
    } else {
      (fold + 1) * fold_size
    }
    // Fit gamma without this fold
    let mut num = 0.0
    let mut den = 0.0
    for j = 0; j < n; j = j + 1 {
      if j >= start && j < end {
        continue
      }
      num = num + observed_bfactors[j] * predicted_msf[j]
      den = den + predicted_msf[j] * predicted_msf[j]
    }
    let gamma = if den < 0.00000000000000000001 {
      1.0
    } else {
      num / den
    }
    // Predict held-out fold
    for j = start; j < end; j = j + 1 {
      all_pred[j] = gamma * predicted_msf[j]
    }
  }
  fit_quality(all_pred, observed_bfactors)
}

// ===== Accessors =====

///|
pub fn ENMFitResult::get_cutoff(self : ENMFitResult) -> Double {
  self.cutoff
}

///|
pub fn ENMFitResult::get_gamma(self : ENMFitResult) -> Double {
  self.gamma
}

///|
pub fn ENMFitResult::get_correlation(self : ENMFitResult) -> Double {
  self.correlation
}

///|
pub fn ENMFitResult::get_rms_error(self : ENMFitResult) -> Double {
  self.rms_error
}

///|
pub fn ENMFitResult::get_predicted_bfactors(self : ENMFitResult) -> Array[Double] {
  self.predicted_bfactors
}

// ===== Advanced fitting =====

///|
/// Fit beta exponent for betaANM.
///
/// Given observed B-factors, search over beta values to find
/// the one that gives best B-factor correlation.
pub fn fit_beta(
  chain : Chain,
  observed_bfactors : Array[Double],
  betas : Array[Double],
) -> (Double, Double, Array[Double])? {
  // Returns (best_beta, best_corr, predicted_bfactors)
  if betas.length() == 0 {
    return None
  }
  let mut best_beta = betas[0]
  let mut best_corr = -2.0
  let mut best_pred : Array[Double] = Array(capacity=0)
  for bi = 0; bi < betas.length(); bi = bi + 1 {
    let beta = betas[bi]
    match BetaANM::new(chain, beta~) {
      Some(anm) => {
        let pred = BetaANM::predict_bfactors(anm)
        let gamma = fit_gamma(pred, observed_bfactors)
        let scaled = scale_bfactors(pred, gamma)
        let (corr, _) = fit_quality(scaled, observed_bfactors)
        if corr > best_corr {
          best_corr = corr
          best_beta = beta
          best_pred = scaled
        }
      }
      None => ()
    }
  }
  Some((best_beta, best_corr, best_pred))
}

///|
/// Compare ANM vs GNM B-factor prediction quality.
pub fn compare_anm_gnm(
  chain : Chain,
  observed_bfactors : Array[Double],
  cutoff : Double,
) -> (Double, Double) {
  // Returns (anm_corr, gnm_corr)
  let anm_corr = match ANM::new(chain, cutoff~) {
    Some(anm) => {
      let pred = ANM::predict_bfactors(anm)
      let gamma = fit_gamma(pred, observed_bfactors)
      let scaled = scale_bfactors(pred, gamma)
      let (corr, _) = fit_quality(scaled, observed_bfactors)
      corr
    }
    None => 0.0
  }
  let gnm_corr = match GNM::new(chain, cutoff~) {
    Some(gnm) => {
      let pred = GNM::predict_bfactors(gnm)
      let gamma = fit_gamma(pred, observed_bfactors)
      let scaled = scale_bfactors(pred, gamma)
      let (corr, _) = fit_quality(scaled, observed_bfactors)
      corr
    }
    None => 0.0
  }
  (anm_corr, gnm_corr)
}

///|
/// Compute fitting profile: correlation vs cutoff.
/// Returns array of (cutoff, correlation) pairs.
pub fn correlation_profile(
  chain : Chain,
  observed_bfactors : Array[Double],
  cutoffs : Array[Double],
) -> Array[(Double, Double)] {
  let result : Array[(Double, Double)] = Array(capacity=0)
  for ci = 0; ci < cutoffs.length(); ci = ci + 1 {
    let cutoff = cutoffs[ci]
    match ANM::new(chain, cutoff~) {
      Some(anm) => {
        let pred = ANM::predict_bfactors(anm)
        let gamma = fit_gamma(pred, observed_bfactors)
        let scaled = scale_bfactors(pred, gamma)
        let (corr, _) = fit_quality(scaled, observed_bfactors)
        result.push((cutoff, corr))
      }
      None => ()
    }
  }
  result
}
