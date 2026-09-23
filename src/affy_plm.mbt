///|
/// Probe-level model fitting and quality control, ported from
/// Bioconductor's `affyPLM` package.
///
/// The reference package fits a robust probe-level model (PLM) to
/// Affymetrix probe intensities using iteratively reweighted least
/// squares (IRLS) with an M-estimator (Huber or Tukey biweight),
/// then derives NUSE (Normalized Unscaled Standard Errors) and RLE
/// (Relative Log Expression) quality statistics.  This port keeps
/// the same pipeline.

///|
/// Type of robust M-estimator for IRLS weighting.
pub(all) enum PlmEstimator {
  Huber(Double)
  Tukey(Double)
} derive(Eq)
pub extend PlmEstimator with Eq::{equal, not_equal}

///|
/// Default Huber estimator with k = 1.345.
pub fn plm_default_estimator() -> PlmEstimator {
  PlmEstimator::Huber(1.345)
}

///|
/// A probeset: one row of probe intensities per sample (probes x
/// samples matrix), plus a name.
pub struct PlmProbeset {
  name : String
  probes : Array[Array[Double]]
} derive(Eq)
pub extend PlmProbeset with Eq::{equal, not_equal}

///|
pub fn PlmProbeset::new(
  name~ : String,
  probes~ : Array[Array[Double]],
) -> PlmProbeset {
  PlmProbeset::{ name, probes }
}

///|
/// Result of fitting a PLM to one probeset.
pub struct PlmFit {
  probeset_name : String
  chip_effects : Array[Double]
  probe_effects : Array[Double]
  residuals : Array[Array[Double]]
  weights : Array[Array[Double]]
  n_iterations : Int
  sigma : Double
} derive(Eq)
pub extend PlmFit with Eq::{equal, not_equal}

///|
/// Absolute value.
fn plm_abs(x : Double) -> Double {
  if x < 0.0 {
    0.0 - x
  } else {
    x
  }
}

///|
/// sqrt via Newton's method.
fn plm_sqrt(x : Double) -> Double {
  if x <= 0.0 {
    return 0.0
  }
  let mut guess = x
  for _ in 0..<50 {
    guess = (guess + x / guess) / 2.0
  }
  guess
}

///|
/// Median absolute deviation (MAD) of a flat array.
fn plm_mad(arr : Array[Double]) -> Double {
  if arr.length() < 2 {
    return 1.0
  }
  let sorted = arr.copy()
  sorted.sort_by(fn(a, b) {
    if a < b {
      -1
    } else if a > b {
      1
    } else {
      0
    }
  })
  let n = sorted.length()
  let med = if n % 2 == 1 {
    sorted[n / 2]
  } else {
    (sorted[n / 2 - 1] + sorted[n / 2]) / 2.0
  }
  let devs : Array[Double] = Array(capacity=0)
  for v in arr {
    devs.push(plm_abs(v - med))
  }
  devs.sort_by(fn(a, b) {
    if a < b {
      -1
    } else if a > b {
      1
    } else {
      0
    }
  })
  let dn = devs.length()
  if dn % 2 == 1 {
    devs[dn / 2]
  } else {
    (devs[dn / 2 - 1] + devs[dn / 2]) / 2.0
  }
}

///|
/// Huber weight function: w(r) = 1 if |r| <= k, k/|r| otherwise.
fn plm_huber_weight(r : Double, k : Double) -> Double {
  let ar = plm_abs(r)
  if ar <= k {
    1.0
  } else if ar == 0.0 {
    1.0
  } else {
    k / ar
  }
}

///|
/// Tukey biweight weight function.
fn plm_tukey_weight(r : Double, k : Double) -> Double {
  let ar = plm_abs(r)
  if ar >= k {
    0.0
  } else {
    let u = ar / k
    let t = 1.0 - u * u
    t * t
  }
}

///|
/// Compute weight for a residual given the estimator and scale.
fn plm_weight(
  residual : Double,
  scale : Double,
  estimator : PlmEstimator,
) -> Double {
  if scale <= 0.0 {
    return 1.0
  }
  let r = residual / scale
  match estimator {
    PlmEstimator::Huber(k) => plm_huber_weight(r, k)
    PlmEstimator::Tukey(k) => plm_tukey_weight(r, k)
  }
}

///|
/// Fit a probe-level model using IRLS.
///
/// Model: log2(intensity[probe, sample]) = chip_effect[sample] +
/// probe_effect[probe] + error.  The fit alternates between weighted
/// least squares and weight recomputation until convergence or
/// `max_iter`.
pub fn fit_plm(
  probeset : PlmProbeset,
  estimator~ : PlmEstimator = plm_default_estimator(),
  max_iter~ : Int = 20,
  tol~ : Double = 1.0e-4,
) -> PlmFit {
  let n_probes = probeset.probes.length()
  let n_samples = if n_probes == 0 { 0 } else { probeset.probes[0].length() }
  if n_probes == 0 || n_samples == 0 {
    return PlmFit::{
      probeset_name: probeset.name,
      chip_effects: [],
      probe_effects: [],
      residuals: [],
      weights: [],
      n_iterations: 0,
      sigma: 0.0,
    }
  }
  let log_data : Array[Array[Double]] = Array(capacity=0)
  for i in 0..<n_probes {
    let row : Array[Double] = Array(capacity=0)
    for j in 0..<n_samples {
      let v = probeset.probes[i][j]
      row.push(plm_log2_safe(v))
    }
    log_data.push(row)
  }
  let chip = Array::make(n_samples, 0.0)
  let probe = Array::make(n_probes, 0.0)
  let mut grand_mean = 0.0
  for i in 0..<n_probes {
    for j in 0..<n_samples {
      grand_mean += log_data[i][j]
    }
  }
  grand_mean = grand_mean / (n_probes * n_samples).to_double()
  for j in 0..<n_samples {
    chip[j] = grand_mean
  }
  let weights : Array[Array[Double]] = Array(capacity=0)
  for i in 0..<n_probes {
    let row : Array[Double] = Array::make(n_samples, 1.0)
    weights.push(row)
  }
  let mut sigma = 1.0
  let mut iter = 0
  let residuals : Array[Array[Double]] = Array(capacity=0)
  for i in 0..<n_probes {
    let row : Array[Double] = Array::make(n_samples, 0.0)
    residuals.push(row)
  }
  while iter < max_iter {
    for i in 0..<n_probes {
      for j in 0..<n_samples {
        residuals[i][j] = log_data[i][j] - chip[j] - probe[i]
      }
    }
    let flat_resid : Array[Double] = Array(capacity=0)
    for i in 0..<n_probes {
      for j in 0..<n_samples {
        flat_resid.push(residuals[i][j])
      }
    }
    let new_sigma = plm_mad(flat_resid) * 1.4826
    if new_sigma > 1.0e-10 {
      sigma = new_sigma
    }
    for i in 0..<n_probes {
      for j in 0..<n_samples {
        weights[i][j] = plm_weight(residuals[i][j], sigma, estimator)
      }
    }
    for j in 0..<n_samples {
      let mut wsum = 0.0
      let mut wval = 0.0
      for i in 0..<n_probes {
        let w = weights[i][j]
        wsum += w
        wval += w * (log_data[i][j] - probe[i])
      }
      if wsum > 1.0e-10 {
        chip[j] = wval / wsum
      }
    }
    for i in 0..<n_probes {
      let mut wsum = 0.0
      let mut wval = 0.0
      for j in 0..<n_samples {
        let w = weights[i][j]
        wsum += w
        wval += w * (log_data[i][j] - chip[j])
      }
      if wsum > 1.0e-10 {
        probe[i] = wval / wsum
      }
    }
    let mut mean_probe = 0.0
    for i in 0..<n_probes {
      mean_probe += probe[i]
    }
    mean_probe = mean_probe / n_probes.to_double()
    for i in 0..<n_probes {
      probe[i] = probe[i] - mean_probe
    }
    for j in 0..<n_samples {
      chip[j] = chip[j] + mean_probe
    }
    iter += 1
    let mut max_change = 0.0
    for i in 0..<n_probes {
      for j in 0..<n_samples {
        let new_resid = log_data[i][j] - chip[j] - probe[i]
        let change = plm_abs(new_resid - residuals[i][j])
        if change > max_change {
          max_change = change
        }
      }
    }
    if max_change < tol {
      break
    }
  }
  PlmFit::{
    probeset_name: probeset.name,
    chip_effects: chip,
    probe_effects: probe,
    residuals,
    weights,
    n_iterations: iter,
    sigma,
  }
}

///|
/// log2 with a floor to avoid -inf for zero/negative values.
fn plm_log2_safe(x : Double) -> Double {
  if x <= 0.0 {
    return -20.0
  }
  let mut lo = -20.0
  let mut hi = 60.0
  for _ in 0..<60 {
    let mid = (lo + hi) / 2.0
    let pow = plm_pow2(mid)
    if pow < x {
      lo = mid
    } else {
      hi = mid
    }
  }
  (lo + hi) / 2.0
}

///|
fn plm_pow2(x : Double) -> Double {
  let ln2 = 0.6931471805599453
  plm_exp(x * ln2)
}

///|
fn plm_exp(x : Double) -> Double {
  if x > 40.0 {
    return plm_exp(40.0) * plm_exp(x - 40.0)
  }
  if x < -40.0 {
    return 0.0
  }
  let mut sum = 1.0
  let mut term = 1.0
  for i in 1..<30 {
    term = term * x / i.to_double()
    sum += term
    let abs_term = if term < 0.0 { 0.0 - term } else { term }
    if abs_term < 1.0e-12 {
      break
    }
  }
  sum
}

///|
/// NUSE (Normalized Unscaled Standard Errors) per sample.
///
/// For each probeset, compute the median weight per sample, then
/// normalize across samples so the median is 1.0.  Values > 1
/// indicate poorer probe behavior for that sample.
pub fn nuse(fits : Array[PlmFit]) -> Array[Double] {
  if fits.is_empty() {
    return []
  }
  let n_samples = fits[0].chip_effects.length()
  if n_samples == 0 {
    return []
  }
  let col_weights : Array[Array[Double]] = Array(capacity=0)
  for _ in 0..<n_samples {
    col_weights.push(Array(capacity=0))
  }
  for fit in fits {
    for j in 0..<n_samples {
      let mut wsum = 0.0
      let mut count = 0
      for i in 0..<fit.weights.length() {
        if j < fit.weights[i].length() {
          wsum += fit.weights[i][j]
          count += 1
        }
      }
      if count > 0 {
        col_weights[j].push(wsum / count.to_double())
      }
    }
  }
  let nuse_vals : Array[Double] = Array(capacity=0)
  for j in 0..<n_samples {
    let med = beclear_median_ref(col_weights[j])
    nuse_vals.push(med)
  }
  let overall = beclear_median_ref(nuse_vals)
  if overall <= 0.0 {
    return nuse_vals
  }
  let normalized : Array[Double] = Array(capacity=0)
  for v in nuse_vals {
    normalized.push(v / overall)
  }
  normalized
}

///|
/// RLE (Relative Log Expression) per sample.
///
/// For each probeset, compute the median chip effect, subtract it
/// from each sample's chip effect, then take the median across
/// probesets per sample.  Values near 0 indicate good quality.
pub fn rle(fits : Array[PlmFit]) -> Array[Double] {
  if fits.is_empty() {
    return []
  }
  let n_samples = fits[0].chip_effects.length()
  if n_samples == 0 {
    return []
  }
  let deviations : Array[Array[Double]] = Array(capacity=0)
  for _ in 0..<n_samples {
    deviations.push(Array(capacity=0))
  }
  for fit in fits {
    let chip = fit.chip_effects
    let med = beclear_median_ref(chip)
    for j in 0..<n_samples {
      deviations[j].push(chip[j] - med)
    }
  }
  let rle_vals : Array[Double] = Array(capacity=0)
  for j in 0..<n_samples {
    rle_vals.push(beclear_median_ref(deviations[j]))
  }
  rle_vals
}

///|
/// Median helper (reference shared with beclear).
fn beclear_median_ref(arr : Array[Double]) -> Double {
  if arr.is_empty() {
    return 0.0
  }
  let sorted = arr.copy()
  sorted.sort_by(fn(a, b) {
    if a < b {
      -1
    } else if a > b {
      1
    } else {
      0
    }
  })
  let n = sorted.length()
  if n % 2 == 1 {
    sorted[n / 2]
  } else {
    (sorted[n / 2 - 1] + sorted[n / 2]) / 2.0
  }
}

///|
/// Combined QC report: NUSE and RLE per sample.
pub struct PlmQCReport {
  nuse : Array[Double]
  rle : Array[Double]
  n_probesets : Int
} derive(Eq)
pub extend PlmQCReport with Eq::{equal, not_equal}

///|
pub fn affy_plm_qc(fits : Array[PlmFit]) -> PlmQCReport {
  PlmQCReport::{
    nuse: nuse(fits),
    rle: rle(fits),
    n_probesets: fits.length(),
  }
}

///|
/// A sample probeset for testing: 4 probes x 3 samples.
pub fn affy_plm_sample_probeset() -> PlmProbeset {
  let probes = [
    [100.0, 110.0, 95.0],
    [200.0, 210.0, 190.0],
    [50.0, 55.0, 48.0],
    [300.0, 310.0, 290.0],
  ]
  PlmProbeset::new(name="ps1", probes~)
}