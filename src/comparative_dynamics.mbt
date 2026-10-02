///|
/// Comparative dynamics analysis.
/// Compare the dynamics of two proteins or two models of the same protein.
/// Based on ProDy's comparative analysis tools.

///|
/// Result of comparing two sets of modes.
pub struct ComparativeResult {
  overlap_mat : Array[Array[Double]] // n1 x n2 overlap matrix
  cumulative_overlaps : Array[Double] // cumulative overlap for each mode in set 1
  effective_overlap : Double // average effective overlap
  rmsd_fluctuations : Double // RMSD between fluctuation profiles
  correlation_fluctuations : Double // Pearson correlation of fluctuations
}

///|
/// Compare two sets of modes.
/// modes1, modes2 are eigenvector arrays.
/// fluctuations1, fluctuations2 are per-residue fluctuation profiles.
pub fn compare_modes(
  modes1 : Array[Array[Double]],
  modes2 : Array[Array[Double]],
  fluctuations1 : Array[Double],
  fluctuations2 : Array[Double],
  k : Int, // number of modes from set 2 to use for cumulative overlap
) -> ComparativeResult {
  let om = overlap_matrix(modes1, modes2)
  // Cumulative overlap for each mode in set 1
  let cum_overlaps = Array::make(modes1.length(), 0.0)
  for i = 0; i < modes1.length(); i = i + 1 {
    cum_overlaps[i] = cumulative_overlap(modes1[i], modes2, k)
  }
  // Average effective overlap
  let mut sum_co = 0.0
  for co in cum_overlaps {
    sum_co = sum_co + co
  }
  let eff_ov = if modes1.length() > 0 {
    sum_co / modes1.length().to_double()
  } else {
    0.0
  }
  // RMSD and correlation of fluctuations
  let (rmsd_fluc, corr_fluc) = compare_fluctuations(fluctuations1, fluctuations2)
  ComparativeResult::{
    overlap_mat: om,
    cumulative_overlaps: cum_overlaps,
    effective_overlap: eff_ov,
    rmsd_fluctuations: rmsd_fluc,
    correlation_fluctuations: corr_fluc,
  }
}

///|
/// Compare two fluctuation profiles.
/// Returns (RMSD, Pearson correlation).
pub fn compare_fluctuations(
  fluc1 : Array[Double],
  fluc2 : Array[Double],
) -> (Double, Double) {
  let n = if fluc1.length() < fluc2.length() {
    fluc1.length()
  } else {
    fluc2.length()
  }
  if n == 0 {
    return (0.0, 0.0)
  }
  // RMSD
  let mut sum_sq = 0.0
  for i = 0; i < n; i = i + 1 {
    let diff = fluc1[i] - fluc2[i]
    sum_sq = sum_sq + diff * diff
  }
  let rmsd = (sum_sq / n.to_double()).sqrt()
  // Pearson correlation
  let mut mean1 = 0.0
  let mut mean2 = 0.0
  for i = 0; i < n; i = i + 1 {
    mean1 = mean1 + fluc1[i]
    mean2 = mean2 + fluc2[i]
  }
  mean1 = mean1 / n.to_double()
  mean2 = mean2 / n.to_double()
  let mut cov = 0.0
  let mut var1 = 0.0
  let mut var2 = 0.0
  for i = 0; i < n; i = i + 1 {
    let d1 = fluc1[i] - mean1
    let d2 = fluc2[i] - mean2
    cov = cov + d1 * d2
    var1 = var1 + d1 * d1
    var2 = var2 + d2 * d2
  }
  let corr = if var1 < 1.0e-20 || var2 < 1.0e-20 {
    0.0
  } else {
    cov / (var1 * var2).sqrt()
  }
  (rmsd, corr)
}

///|
/// Dynamics similarity score.
/// Combines overlap and fluctuation correlation into a single score [0, 1].
pub fn dynamics_similarity(
  modes1 : Array[Array[Double]],
  modes2 : Array[Array[Double]],
  fluctuations1 : Array[Double],
  fluctuations2 : Array[Double],
  k : Int,
) -> Double {
  let result = compare_modes(modes1, modes2, fluctuations1, fluctuations2, k)
  // Similarity = 0.5 * (effective_overlap + max(0, correlation))
  let corr_pos = if result.correlation_fluctuations > 0.0 {
    result.correlation_fluctuations
  } else {
    0.0
  }
  0.5 * (result.effective_overlap + corr_pos)
}

///|
/// Compare B-factor predictions from two models.
/// Returns (correlation, RMSD).
pub fn compare_bfactors(
  bfactors1 : Array[Double],
  bfactors2 : Array[Double],
) -> (Double, Double) {
  let n = if bfactors1.length() < bfactors2.length() {
    bfactors1.length()
  } else {
    bfactors2.length()
  }
  if n == 0 {
    return (0.0, 0.0)
  }
  // Normalize both to [0, 1]
  let mut max1 = 0.0
  let mut max2 = 0.0
  for i = 0; i < n; i = i + 1 {
    if bfactors1[i] > max1 {
      max1 = bfactors1[i]
    }
    if bfactors2[i] > max2 {
      max2 = bfactors2[i]
    }
  }
  let norm1 = Array::make(n, 0.0)
  let norm2 = Array::make(n, 0.0)
  for i = 0; i < n; i = i + 1 {
    norm1[i] = if max1 > 1.0e-20 { bfactors1[i] / max1 } else { 0.0 }
    norm2[i] = if max2 > 1.0e-20 { bfactors2[i] / max2 } else { 0.0 }
  }
  let (rmsd, corr) = compare_fluctuations(norm1, norm2)
  (corr, rmsd)
}

///|
/// Compute the degree of collectivity for a set of modes.
/// Returns average collectivity across the first k modes.
pub fn average_collectivity(
  eigenvectors : Array[Array[Double]],
  k : Int,
) -> Double {
  let max_k = if k < eigenvectors.length() { k } else {
    eigenvectors.length()
  }
  if max_k == 0 {
    return 0.0
  }
  let mut total = 0.0
  for m = 0; m < max_k; m = m + 1 {
    let vec = eigenvectors[m]
    let n = vec.length()
    let mut sum_sq = 0.0
    let mut sum_4th = 0.0
    for i = 0; i < n; i = i + 1 {
      let s = vec[i] * vec[i]
      sum_sq = sum_sq + s
      sum_4th = sum_4th + s * s
    }
    if sum_4th > 1.0e-20 {
      let pr = (sum_sq * sum_sq) / (n.to_double() * sum_4th)
      total = total + pr
    }
  }
  total / max_k.to_double()
}
