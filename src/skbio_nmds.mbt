/// Non-metric Multidimensional Scaling (NMDS).
///
/// scikit-bio equivalent: `skbio.ordination.nmds()`
///
/// NMDS is an ordination method that preserves rank order of distances
/// rather than actual distances. Widely used in ecology when dissimilarity
/// measures are non-Euclidean (e.g., Bray-Curtis, Jaccard).
///
/// Algorithm: Kruskal's classic NMDS with monotone regression.

///|
/// NMDS result
pub struct SkbioNmdsResult {
  samples : Array[Array[Double]]
  stress : Double
  dimensions : Int32
  iterations : Int32
  converged : Bool
}

///|
pub fn SkbioNmdsResult::samples(self : SkbioNmdsResult) -> Array[Array[Double]] {
  self.samples
}

///|
pub fn SkbioNmdsResult::stress(self : SkbioNmdsResult) -> Double {
  self.stress
}

///|
pub fn SkbioNmdsResult::dimensions(self : SkbioNmdsResult) -> Int32 {
  self.dimensions
}

///|
pub fn SkbioNmdsResult::iterations(self : SkbioNmdsResult) -> Int32 {
  self.iterations
}

///|
pub fn SkbioNmdsResult::converged(self : SkbioNmdsResult) -> Bool {
  self.converged
}

///|
/// Compute Euclidean distance between two points
fn nmds_euclid(a : Array[Double], b : Array[Double]) -> Double {
  let mut sum = 0.0
  for i = 0; i < a.length(); i = i + 1 {
    let d = a[i] - b[i]
    sum = sum + d * d
  }
  sum.sqrt()
}

///|
/// Kruskal's monotone regression: find fitted distances that are
/// monotone with respect to the rank order of original distances.
fn nmds_monotone_regression(
  orig : Array[Double],
  config : Array[Double],
) -> Array[Double] {
  let n = orig.length()
  let indices : Array[Int32] = Array::make(n, 0)
  for i = 0; i < n; i = i + 1 {
    indices[i] = i
  }
  indices.sort_by(fn(a, b) {
    if orig[a] < orig[b] {
      -1
    } else if orig[a] > orig[b] {
      1
    } else {
      0
    }
  })
  // Pool-adjacent-violators algorithm
  let fitted = config.copy()
  let mut i : Int32 = 0
  while i < n {
    let mut j = i
    let mut sum = fitted[indices[i]]
    let mut count = 1.0
    while j + 1 < n && fitted[indices[j + 1]] < fitted[indices[j]] {
      j = j + 1
      sum = sum + fitted[indices[j]]
      count = count + 1.0
    }
    let avg = sum / count
    for k = i; k <= j; k = k + 1 {
      fitted[indices[k]] = avg
    }
    i = j + 1
  }
  fitted
}

///|
/// Perform NMDS on a distance matrix.
///
/// dm: n×n dissimilarity matrix (symmetric, zero diagonal)
/// dimensions: number of dimensions (typically 2 or 3)
/// max_iter: maximum iterations (default 300)
/// eps: convergence threshold for stress change (default 1e-6)
pub fn skbio_nmds(
  dm : Array[Array[Double]],
  dimensions~ : Int32 = 2,
  max_iter~ : Int32 = 300,
  eps~ : Double = 0.000001,
) -> SkbioNmdsResult {
  let n = dm.length()
  if n <= 1 {
    return {
      samples: [],
      stress: 0.0,
      dimensions,
      iterations: 0,
      converged: true,
    }
  }
  let dim = dimensions
  // Initialize configuration with pseudo-random values
  let mut seed = 42.0
  let samples : Array[Array[Double]] = Array::make(n, [])
  for i = 0; i < n; i = i + 1 {
    let coords : Array[Double] = Array::make(dim, 0.0)
    for j = 0; j < dim; j = j + 1 {
      seed = seed * 1103515245.0 + 12345.0
      seed = seed - (seed / 2147483647.0).floor() * 2147483647.0
      coords[j] = (seed / 2147483647.0 - 0.5) * 2.0
    }
    samples[i] = coords
  }
  // Collect original distances (upper triangle)
  let orig_dists : Array[Double] = Array(capacity=0)
  let orig_pairs : Array[(Int32, Int32)] = Array(capacity=0)
  for i = 0; i < n; i = i + 1 {
    for j = i + 1; j < n; j = j + 1 {
      orig_dists.push(dm[i][j])
      orig_pairs.push((i, j))
    }
  }
  let m = orig_dists.length()
  let mut prev_stress = 1000000000000.0
  let mut iter : Int32 = 0
  let mut converged = false
  let mut final_stress = 0.0
  while iter < max_iter {
    // Compute configuration distances
    let config_dists : Array[Double] = Array::make(m, 0.0)
    for k = 0; k < m; k = k + 1 {
      config_dists[k] = nmds_euclid(
        samples[orig_pairs[k].0],
        samples[orig_pairs[k].1],
      )
    }
    // Monotone regression
    let fitted = nmds_monotone_regression(orig_dists, config_dists)
    // Compute stress
    let mut numerator = 0.0
    let mut denominator = 0.0
    for k = 0; k < m; k = k + 1 {
      let diff = config_dists[k] - fitted[k]
      numerator = numerator + diff * diff
      denominator = denominator + fitted[k] * fitted[k]
    }
    let stress = if denominator > 0.0 {
      (numerator / denominator).sqrt()
    } else {
      0.0
    }
    final_stress = stress
    // Check convergence
    let mut stress_diff = prev_stress - stress
    if stress_diff < 0.0 {
      // abs
      stress_diff = -stress_diff
    }
    if stress_diff < eps {
      converged = true
      break
    }
    prev_stress = stress
    // Gradient descent step
    let lr = 0.01
    for k = 0; k < m; k = k + 1 {
      let i = orig_pairs[k].0
      let j = orig_pairs[k].1
      let d = config_dists[k]
      if d > 0.0 {
        let factor = lr * (d - fitted[k]) / d
        for c = 0; c < dim; c = c + 1 {
          let diff = samples[i][c] - samples[j][c]
          samples[i][c] = samples[i][c] - factor * diff
          samples[j][c] = samples[j][c] + factor * diff
        }
      }
    }
    iter = iter + 1
  }
  { samples, stress: final_stress, dimensions, iterations: iter, converged, }
}
