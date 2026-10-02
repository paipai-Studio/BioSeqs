///|
/// Signature dynamics analysis.
/// Compute dynamics signatures for proteins and compare them.
/// A signature captures the essential dynamics character: slow mode profiles,
/// mean-square fluctuations, participation ratios.
/// Based on ProDy's signature dynamics approach.

///|
/// Dynamics signature for a protein.
pub(all) struct DynamicsSignature {
  n_residues : Int
  slowest_mode_profile : Array[Double] // |eigenvector| of slowest mode
  msf_profile : Array[Double] // mean-square fluctuation per residue
  participation_profile : Array[Double] // participation ratio per mode
  collectivity_profile : Array[Double] // collectivity per mode
  mean_fluctuation : Double
  max_fluctuation : Double
  min_fluctuation : Double
}

///|
/// Compute dynamics signature from eigenvalues and eigenvectors.
/// eigenvectors are sorted by ascending eigenvalue (slowest first).
pub fn dynamics_signature(
  eigenvalues : Array[Double],
  eigenvectors : Array[Array[Double]],
  n_modes? : Int = 5,
) -> DynamicsSignature {
  let n_modes_eff = if n_modes < eigenvalues.length() {
    n_modes
  } else {
    eigenvalues.length()
  }
  let n_res = if eigenvectors.length() > 0 {
    eigenvectors[0].length()
  } else {
    0
  }
  // Slowest mode profile (absolute values of first non-trivial eigenvector)
  let slowest_profile = Array::make(n_res, 0.0)
  if n_modes_eff > 0 && n_res > 0 {
    let vec = eigenvectors[0]
    for i = 0; i < n_res; i = i + 1 {
      slowest_profile[i] = vec[i].abs()
    }
  }
  // MSF profile: sum_k (1/lambda_k) * v_k^2
  let msf = Array::make(n_res, 0.0)
  for k = 0; k < n_modes_eff; k = k + 1 {
    let lam = eigenvalues[k]
    if lam < 1.0e-10 {
      continue
    }
    let vec = eigenvectors[k]
    for i = 0; i < n_res; i = i + 1 {
      msf[i] = msf[i] + vec[i] * vec[i] / lam
    }
  }
  // Participation ratio per mode
  let participation = Array::make(n_modes_eff, 0.0)
  for k = 0; k < n_modes_eff; k = k + 1 {
    let vec = eigenvectors[k]
    let mut sum_sq = 0.0
    let mut sum_4th = 0.0
    for i = 0; i < n_res; i = i + 1 {
      let s = vec[i] * vec[i]
      sum_sq = sum_sq + s
      sum_4th = sum_4th + s * s
    }
    if sum_4th > 1.0e-20 {
      participation[k] = (sum_sq * sum_sq) / (n_res.to_double() * sum_4th)
    }
  }
  // Collectivity per mode (using exponential form)
  let collectivity = Array::make(n_modes_eff, 0.0)
  for k = 0; k < n_modes_eff; k = k + 1 {
    let vec = eigenvectors[k]
    let mut sum_sq = 0.0
    let mut sum_4th = 0.0
    for i = 0; i < n_res; i = i + 1 {
      let s = vec[i] * vec[i]
      sum_sq = sum_sq + s
      sum_4th = sum_4th + s * s
    }
    if sum_4th > 1.0e-20 && sum_sq > 1.0e-20 {
      let pr = (sum_sq * sum_sq) / (n_res.to_double() * sum_4th)
      collectivity[k] = pr
    }
  }
  // Statistics
  let mut mean_fluc = 0.0
  let mut max_fluc = 0.0
  let mut min_fluc = if n_res > 0 { msf[0] } else { 0.0 }
  for i = 0; i < n_res; i = i + 1 {
    mean_fluc = mean_fluc + msf[i]
    if msf[i] > max_fluc {
      max_fluc = msf[i]
    }
    if msf[i] < min_fluc {
      min_fluc = msf[i]
    }
  }
  if n_res > 0 {
    mean_fluc = mean_fluc / n_res.to_double()
  }
  DynamicsSignature::{
    n_residues: n_res,
    slowest_mode_profile: slowest_profile,
    msf_profile: msf,
    participation_profile: participation,
    collectivity_profile: collectivity,
    mean_fluctuation: mean_fluc,
    max_fluctuation: max_fluc,
    min_fluctuation: min_fluc,
  }
}

///|
/// Compute rank of values (for rank correlation).
/// Returns ranks (1-based, ties get average rank).
pub fn rank_values(values : Array[Double]) -> Array[Double] {
  let n = values.length()
  let indices = Array::makei(n, { i => i })
  for sort_i = 1; sort_i < n; sort_i = sort_i + 1 {
    let key = indices[sort_i]
    let key_val = values[key]
    let mut sort_j = sort_i - 1
    while sort_j >= 0 && values[indices[sort_j]] > key_val {
      indices[sort_j + 1] = indices[sort_j]
      sort_j = sort_j - 1
    }
    indices[sort_j + 1] = key
  }
  let ranks = Array::make(n, 0.0)
  let mut i = 0
  while i < n {
    // Handle ties
    let mut j = i + 1
    while j < n && values[indices[j]] == values[indices[i]] {
      j = j + 1
    }
    let avg_rank = (i.to_double() + j.to_double() + 1.0) / 2.0
    for kk = i; kk < j; kk = kk + 1 {
      ranks[indices[kk]] = avg_rank
    }
    i = j
  }
  ranks
}

///|
/// Spearman rank correlation between two arrays.
pub fn rank_correlation(
  values1 : Array[Double],
  values2 : Array[Double],
) -> Double {
  let n = if values1.length() < values2.length() {
    values1.length()
  } else {
    values2.length()
  }
  if n < 2 {
    return 0.0
  }
  let ranks1 = rank_values(values1)
  let ranks2 = rank_values(values2)
  // Pearson correlation on ranks
  let mut mean1 = 0.0
  let mut mean2 = 0.0
  for i = 0; i < n; i = i + 1 {
    mean1 = mean1 + ranks1[i]
    mean2 = mean2 + ranks2[i]
  }
  mean1 = mean1 / n.to_double()
  mean2 = mean2 / n.to_double()
  let mut cov = 0.0
  let mut var1 = 0.0
  let mut var2 = 0.0
  for i = 0; i < n; i = i + 1 {
    let d1 = ranks1[i] - mean1
    let d2 = ranks2[i] - mean2
    cov = cov + d1 * d2
    var1 = var1 + d1 * d1
    var2 = var2 + d2 * d2
  }
  if var1 < 1.0e-20 || var2 < 1.0e-20 {
    return 0.0
  }
  cov / (var1 * var2).sqrt()
}

///|
/// Correlation between two dynamics signatures.
/// Compares MSF profiles using Pearson correlation.
pub fn signature_correlation(
  sig1 : DynamicsSignature,
  sig2 : DynamicsSignature,
) -> Double {
  let n = if sig1.msf_profile.length() < sig2.msf_profile.length() {
    sig1.msf_profile.length()
  } else {
    sig2.msf_profile.length()
  }
  if n == 0 {
    return 0.0
  }
  let mut mean1 = 0.0
  let mut mean2 = 0.0
  for i = 0; i < n; i = i + 1 {
    mean1 = mean1 + sig1.msf_profile[i]
    mean2 = mean2 + sig2.msf_profile[i]
  }
  mean1 = mean1 / n.to_double()
  mean2 = mean2 / n.to_double()
  let mut cov = 0.0
  let mut var1 = 0.0
  let mut var2 = 0.0
  for i = 0; i < n; i = i + 1 {
    let d1 = sig1.msf_profile[i] - mean1
    let d2 = sig2.msf_profile[i] - mean2
    cov = cov + d1 * d2
    var1 = var1 + d1 * d1
    var2 = var2 + d2 * d2
  }
  if var1 < 1.0e-20 || var2 < 1.0e-20 {
    return 0.0
  }
  cov / (var1 * var2).sqrt()
}

///|
/// Distance between two dynamics signatures.
/// Combines MSF profile distance and participation profile distance.
pub fn signature_distance(
  sig1 : DynamicsSignature,
  sig2 : DynamicsSignature,
) -> Double {
  // MSF profile distance (normalized RMSD)
  let n = if sig1.msf_profile.length() < sig2.msf_profile.length() {
    sig1.msf_profile.length()
  } else {
    sig2.msf_profile.length()
  }
  let mut msf_dist = 0.0
  if n > 0 {
    let mut sum_sq = 0.0
    for i = 0; i < n; i = i + 1 {
      let diff = sig1.msf_profile[i] - sig2.msf_profile[i]
      sum_sq = sum_sq + diff * diff
    }
    msf_dist = (sum_sq / n.to_double()).sqrt()
    // Normalize by max fluctuation
    let max_fluc = if sig1.max_fluctuation > sig2.max_fluctuation {
      sig1.max_fluctuation
    } else {
      sig2.max_fluctuation
    }
    if max_fluc > 1.0e-20 {
      msf_dist = msf_dist / max_fluc
    }
  }
  // Participation profile distance
  let np = if sig1.participation_profile.length() < sig2.participation_profile.length() {
    sig1.participation_profile.length()
  } else {
    sig2.participation_profile.length()
  }
  let mut part_dist = 0.0
  if np > 0 {
    let mut sum_sq = 0.0
    for i = 0; i < np; i = i + 1 {
      let diff = sig1.participation_profile[i] - sig2.participation_profile[i]
      sum_sq = sum_sq + diff * diff
    }
    part_dist = (sum_sq / np.to_double()).sqrt()
  }
  // Combined distance
  (msf_dist * msf_dist + part_dist * part_dist).sqrt()
}

///|
/// Rank-based signature comparison.
/// Uses Spearman rank correlation on MSF profiles.
pub fn rank_signature(
  sig1 : DynamicsSignature,
  sig2 : DynamicsSignature,
) -> Double {
  rank_correlation(sig1.msf_profile, sig2.msf_profile)
}
