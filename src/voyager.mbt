///|
/// Voyager — Bioconductor spatial transcriptomics univariate, bivariate and
/// multivariate spatial autocorrelation statistics, inspired by the Voyager
/// Bioconductor package (Pardo, Spakov, Wolf et al.).
///
/// Implements:
///   - Spatial neighbour/weight construction (kNN, distance band, inverse
///     distance) with row-standardized / binary / globally-standardized /
///     Caussinus-Mestre style coding.
///   - Global Moran's I and Geary's c with the Cliff–Ord randomization
///     expectation and variance and a normal-approximation p-value.
///   - Local Moran's I (LISA) with quadrant classification and permutation
///     inference, local_res Geary's c, and local_res Getis–Ord Gi / Gi* with the
///     exact Ord–Getis randomization z-score.
///   - Bivariate Lee's L (global and local_res) and multivariate local_res Geary.
///   - Empirical variogram with spherical / exponential / Gaussian model
///     fitting by bounded grid search + linear least squares.
///   - Spatial Moran correlogram over distance bins.
///   - Deterministic permutation tests (seeded splitmix64 PRNG) and
///     Benjamini–Hochberg FDR.
///   - Immutable SpatialExperiment write-back of local_res results.

// ===========================================================================
// Error type and validation helpers
// ===========================================================================

///|
pub suberror VoyagerError {
  VoyagerError(String)
}

///|
fn voy_fail(message : String) -> Unit raise VoyagerError {
  raise VoyagerError::VoyagerError(message)
}

///|
fn voy_is_finite(value : Double) -> Bool {
  value == value && value.abs() <= 1.0e300
}

///|
fn voy_copy_double_matrix(
  values : Array[Array[Double]],
) -> Array[Array[Double]] {
  let output : Array[Array[Double]] = []
  for row in values {
    output.push(row.copy())
  }
  output
}

///|
fn voy_copy_string_map(values : Map[String, String]) -> Map[String, String] {
  let output : Map[String, String] = Map([], capacity=values.keys().length())
  for key in values.keys() {
    output[key] = values[key]
  }
  output
}

///|
fn voy_default_names(prefix : String, count : Int) -> Array[String] {
  let output : Array[String] = []
  for index in 0..<count {
    output.push(prefix + (index + 1).to_string())
  }
  output
}

///|
fn voy_validate_vector(
  values : Array[Double],
  label : String,
) -> Unit raise VoyagerError {
  if values.length() < 3 {
    voy_fail("Voyager " + label + " requires at least 3 observations")
  }
  for value in values {
    if !voy_is_finite(value) {
      voy_fail("Voyager " + label + " must be finite")
    }
  }
}

///|
fn voy_validate_weights(
  weights : VoyagerWeights,
  label : String,
) -> Unit raise VoyagerError {
  if weights.n < 3 {
    voy_fail("Voyager " + label + " weights require at least 3 spots")
  }
  if weights.neighbors.length() != weights.n {
    voy_fail("Voyager " + label + " weights neighbours must match spot count")
  }
  if weights.weights.length() != weights.n {
    voy_fail("Voyager " + label + " weights rows must match spot count")
  }
  let mut has_link = false
  for index in 0..<weights.n {
    if weights.neighbors[index].length() != weights.weights[index].length() {
      voy_fail(
        "Voyager " +
        label +
        " neighbour/weight length mismatch at spot " +
        index.to_string(),
      )
    }
    for value in weights.weights[index] {
      if !voy_is_finite(value) {
        voy_fail(
          "Voyager " +
          label +
          " weights must be finite at spot " +
          index.to_string(),
        )
      }
      if value < 0.0 {
        voy_fail(
          "Voyager " +
          label +
          " weights must be non-negative at spot " +
          index.to_string(),
        )
      }
    }
    if weights.neighbors[index].length() > 0 {
      has_link = true
    }
  }
  if !has_link {
    voy_fail("Voyager " + label + " weights must contain at least one link")
  }
}

///|
fn voy_validate_equal_length(
  a : Array[Double],
  b : Array[Double],
  label : String,
) -> Unit raise VoyagerError {
  if a.length() != b.length() {
    voy_fail("Voyager " + label + " vectors must have equal length")
  }
}

// ===========================================================================
// Deterministic PRNG (splitmix64) for reproducible permutations
// ===========================================================================

///|
priv struct VoyagerRng {
  mut state : Int64
}

///|
fn VoyagerRng::new(seed : Int) -> VoyagerRng {
  let init = if seed < 0 {
    seed.to_int64() + 9223372036854775807L
  } else {
    seed.to_int64()
  }
  VoyagerRng::{ state: init + 0x9E3779B97F4A7C15L, }
}

///|
fn VoyagerRng::next_uint64(self : VoyagerRng) -> Int64 {
  let mut z = self.state
  z = z + 0x9E3779B97F4A7C15L
  self.state = z
  z = z ^ (z >> 30)
  z = z * 0xBF58476D1CE4E5B9L
  z = z ^ (z >> 27)
  z = z * 0x94D049BB133111EBL
  z = z ^ (z >> 31)
  z
}

///|
fn VoyagerRng::next_int_range(self : VoyagerRng, upper_exclusive : Int) -> Int {
  if upper_exclusive <= 0 {
    return 0
  }
  let value = self.next_uint64()
  if upper_exclusive <= 2147483647 {
    (value & 2147483647L).to_int() % upper_exclusive
  } else {
    value.to_int() % upper_exclusive
  }
}

///|
/// Fisher–Yates shuffle (in place on a copy) using the deterministic PRNG.
fn voy_permute(values : Array[Double], rng : VoyagerRng) -> Array[Double] {
  let output = values.copy()
  let n = output.length()
  let mut i = n - 1
  while i > 0 {
    let j = rng.next_int_range(i + 1)
    let tmp = output[i]
    output[i] = output[j]
    output[j] = tmp
    i = i - 1
  }
  output
}

// ===========================================================================
// Benjamini–Hochberg FDR
// ===========================================================================

///|
fn voy_bh_adjust(p_values : Array[Double]) -> Array[Double] {
  let count = p_values.length()
  let output : Array[Double] = Array::make(count, 1.0)
  if count == 0 {
    return output
  }
  let indexed : Array[(Double, Int)] = []
  for index in 0..<count {
    indexed.push((p_values[index].max(0.0).min(1.0), index))
  }
  indexed.sort_by(fn(left : (Double, Int), right : (Double, Int)) -> Int {
    if left.0 < right.0 {
      -1
    } else if left.0 > right.0 {
      1
    } else {
      left.1 - right.1
    }
  })
  let mut running = 1.0
  let mut position = count - 1
  while position >= 0 {
    let candidate = (indexed[position].0 *
    count.to_double() /
    (position + 1).to_double()).min(1.0)
    running = running.min(candidate)
    output[indexed[position].1] = running
    position = position - 1
  }
  output
}

// ===========================================================================
// Normal CDF (Abramowitz–Stegun) for analytical p-values
// ===========================================================================

///|
fn voy_normal_cdf(x : Double) -> Double {
  let abs_x = if x < 0.0 { -x } else { x }
  let t = 1.0 / (1.0 + 0.2316419 * abs_x)
  let d = 0.3989422804014327 * @math.exp(-x * x / 2.0)
  let p = d *
    t *
    (
      0.31938153 +
      t *
      (-0.356563782 + t * (1.781477937 + t * (-1.821255978 + t * 1.330274429)))
    )
  if x > 0.0 {
    1.0 - p
  } else {
    p
  }
}

///|
fn voy_two_sided_p(z : Double) -> Double {
  let abs_z = if z < 0.0 { -z } else { z }
  2.0 * (1.0 - voy_normal_cdf(abs_z))
}

// ===========================================================================
// Distance and weights
// ===========================================================================

///|
fn voy_euclidean(a : Array[Double], b : Array[Double]) -> Double {
  let mut sum = 0.0
  let limit = if a.length() < b.length() { a.length() } else { b.length() }
  for index in 0..<limit {
    let diff = a[index] - b[index]
    sum = sum + diff * diff
  }
  sum.sqrt()
}

///|
fn voy_pairwise_distances(
  coords : Array[Array[Double]],
) -> Array[Array[Double]] {
  let n = coords.length()
  // Build a fresh inner array per row — `Array::make(n, Array::make(n, 0.0))`
  // would alias a single inner array across all rows.
  let output : Array[Array[Double]] = []
  for _ in 0..<n {
    output.push(Array::make(n, 0.0))
  }
  for i in 0..<n {
    for j in (i + 1)..<n {
      let d = voy_euclidean(coords[i], coords[j])
      output[i][j] = d
      output[j][i] = d
    }
  }
  output
}

///|
/// Weight-coding style for the spatial weights matrix.
///   "W" — row-standardized (each row sums to 1)
///   "B" — binary (w_ij = 1 for each link)
///   "C" — globally standardized (w_ij = 1 / S0, total sum = 1)
///   "S" — Caussinus–Mestre (w_ij = 1 / sqrt(k_i * k_j))
pub struct VoyagerWeights {
  style : String
  n : Int
  neighbors : Array[Array[Int]]
  weights : Array[Array[Double]]
} derive(Debug)
pub extend VoyagerWeights with @moonbitlang/core/debug.Debug::{to_repr}

///|
pub fn VoyagerWeights::new(
  neighbors : Array[Array[Int]],
  weights : Array[Array[Double]],
  style? : String = "W",
) -> VoyagerWeights raise VoyagerError {
  let n = neighbors.length()
  if weights.length() != n {
    voy_fail("Voyager weights neighbour/weight row counts must match")
  }
  for index in 0..<n {
    if neighbors[index].length() != weights[index].length() {
      voy_fail(
        "Voyager weights neighbour/weight length mismatch at spot " +
        index.to_string(),
      )
    }
    for value in weights[index] {
      if !voy_is_finite(value) || value < 0.0 {
        voy_fail(
          "Voyager weights entry must be finite and non-negative at spot " +
          index.to_string(),
        )
      }
    }
  }
  VoyagerWeights::{ style, n, neighbors, weights, }
}

///|
/// Total sum of all weights, S0 = sum_i sum_j w_ij.
pub fn voyager_weights_s0(weights : VoyagerWeights) -> Double {
  let mut total = 0.0
  for index in 0..<weights.n {
    for value in weights.weights[index] {
      total = total + value
    }
  }
  total
}

///|
/// S1 = 0.5 * sum_i sum_j (w_ij + w_ji)^2.
pub fn voyager_weights_s1(weights : VoyagerWeights) -> Double {
  let n = weights.n
  let mut total = 0.0
  for i in 0..<n {
    for k in 0..<weights.neighbors[i].length() {
      let j = weights.neighbors[i][k]
      let wij = weights.weights[i][k]
      let wji = voy_weights_lookup(weights, j, i)
      let s = wij + wji
      total = total + s * s
    }
  }
  total / 2.0
}

///|
/// S2 = sum_i (sum_j w_ij + sum_j w_ji)^2.
pub fn voyager_weights_s2(weights : VoyagerWeights) -> Double {
  let n = weights.n
  let row_sums : Array[Double] = Array::make(n, 0.0)
  let col_sums : Array[Double] = Array::make(n, 0.0)
  for i in 0..<n {
    for k in 0..<weights.neighbors[i].length() {
      let j = weights.neighbors[i][k]
      let wij = weights.weights[i][k]
      row_sums[i] = row_sums[i] + wij
      col_sums[j] = col_sums[j] + wij
    }
  }
  let mut total = 0.0
  for i in 0..<n {
    let s = row_sums[i] + col_sums[i]
    total = total + s * s
  }
  total
}

///|
fn voy_weights_lookup(weights : VoyagerWeights, i : Int, j : Int) -> Double {
  if i < 0 || i >= weights.n {
    return 0.0
  }
  let neighbors = weights.neighbors[i]
  let values = weights.weights[i]
  for k in 0..<neighbors.length() {
    if neighbors[k] == j {
      return values[k]
    }
  }
  0.0
}

///|
fn voy_apply_style(
  neighbors : Array[Array[Int]],
  raw : Array[Array[Double]],
  style : String,
) -> Array[Array[Double]] {
  let n = neighbors.length()
  match style {
    "B" => {
      let output : Array[Array[Double]] = []
      for index in 0..<n {
        let row : Array[Double] = []
        for _ in 0..<neighbors[index].length() {
          row.push(1.0)
        }
        output.push(row)
      }
      output
    }
    "C" => {
      let mut s0 = 0.0
      for index in 0..<n {
        for value in raw[index] {
          s0 = s0 + value
        }
      }
      let scale = if s0 > 0.0 { 1.0 / s0 } else { 0.0 }
      let output : Array[Array[Double]] = []
      for index in 0..<n {
        let row : Array[Double] = []
        for value in raw[index] {
          row.push(value * scale)
        }
        output.push(row)
      }
      output
    }
    "S" => {
      let k : Array[Double] = Array::make(n, 0.0)
      for index in 0..<n {
        for value in raw[index] {
          k[index] = k[index] + value
        }
      }
      let output : Array[Array[Double]] = []
      for index in 0..<n {
        let row : Array[Double] = []
        for nb in 0..<neighbors[index].length() {
          let j = neighbors[index][nb]
          let denom = (k[index] * k[j]).max(1.0e-300)
          row.push(1.0 / denom.sqrt())
        }
        output.push(row)
      }
      output
    }
    _ => {
      // "W" — row-standardized
      let output : Array[Array[Double]] = []
      for index in 0..<n {
        let mut row_sum = 0.0
        for value in raw[index] {
          row_sum = row_sum + value
        }
        let row : Array[Double] = []
        if row_sum > 0.0 {
          for value in raw[index] {
            row.push(value / row_sum)
          }
        } else {
          for _ in raw[index] {
            row.push(0.0)
          }
        }
        output.push(row)
      }
      output
    }
  }
}

///|
/// Build k-nearest-neighbour spatial weights.
/// `coords` is one row per spot. `k` is the number of neighbours (excluding self).
pub fn voyager_weights_knn(
  coords : Array[Array[Double]],
  k : Int,
  style? : String = "W",
) -> VoyagerWeights raise VoyagerError {
  if coords.length() < 3 {
    voy_fail("Voyager kNN weights require at least 3 spots")
  }
  if k < 1 {
    voy_fail("Voyager kNN weights require k >= 1")
  }
  let n = coords.length()
  for coordinate in coords {
    if coordinate.length() == 0 {
      voy_fail("Voyager kNN coordinates must not be empty")
    }
    for value in coordinate {
      if !voy_is_finite(value) {
        voy_fail("Voyager kNN coordinates must be finite")
      }
    }
  }
  let effective_k = k.min(n - 1)
  let distances = voy_pairwise_distances(coords)
  let neighbors : Array[Array[Int]] = []
  let raw : Array[Array[Double]] = []
  for i in 0..<n {
    let pairs : Array[(Double, Int)] = []
    for j in 0..<n {
      if j != i {
        pairs.push((distances[i][j], j))
      }
    }
    pairs.sort_by(fn(a : (Double, Int), b : (Double, Int)) -> Int {
      if a.0 < b.0 {
        -1
      } else if a.0 > b.0 {
        1
      } else {
        a.1 - b.1
      }
    })
    let nb : Array[Int] = []
    let w : Array[Double] = []
    for position in 0..<effective_k {
      nb.push(pairs[position].1)
      w.push(1.0)
    }
    neighbors.push(nb)
    raw.push(w)
  }
  let styled = voy_apply_style(neighbors, raw, style)
  VoyagerWeights::{ style, n, neighbors, weights: styled, }
}

///|
/// Build distance-band (fixed-radius) spatial weights. Spots within
/// `bandwidth` distance (exclusive of self) are linked.
pub fn voyager_weights_distance_band(
  coords : Array[Array[Double]],
  bandwidth : Double,
  style? : String = "W",
) -> VoyagerWeights raise VoyagerError {
  if coords.length() < 3 {
    voy_fail("Voyager distance-band weights require at least 3 spots")
  }
  if !voy_is_finite(bandwidth) || bandwidth <= 0.0 {
    voy_fail("Voyager distance-band bandwidth must be positive and finite")
  }
  let n = coords.length()
  for coordinate in coords {
    if coordinate.length() == 0 {
      voy_fail("Voyager distance-band coordinates must not be empty")
    }
    for value in coordinate {
      if !voy_is_finite(value) {
        voy_fail("Voyager distance-band coordinates must be finite")
      }
    }
  }
  let neighbors : Array[Array[Int]] = []
  let raw : Array[Array[Double]] = []
  for i in 0..<n {
    let nb : Array[Int] = []
    let w : Array[Double] = []
    for j in 0..<n {
      if i != j {
        let d = voy_euclidean(coords[i], coords[j])
        if d <= bandwidth {
          nb.push(j)
          w.push(1.0)
        }
      }
    }
    neighbors.push(nb)
    raw.push(w)
  }
  let styled = voy_apply_style(neighbors, raw, style)
  VoyagerWeights::{ style, n, neighbors, weights: styled, }
}

///|
/// Build inverse-distance spatial weights w_ij = 1 / d(i,j)^power,
/// constrained to a maximum `bandwidth` (0 = unlimited).
pub fn voyager_weights_inverse_distance(
  coords : Array[Array[Double]],
  power : Double,
  bandwidth? : Double = 0.0,
  style? : String = "W",
) -> VoyagerWeights raise VoyagerError {
  if coords.length() < 3 {
    voy_fail("Voyager inverse-distance weights require at least 3 spots")
  }
  if !voy_is_finite(power) || power <= 0.0 {
    voy_fail("Voyager inverse-distance power must be positive and finite")
  }
  let n = coords.length()
  for coordinate in coords {
    if coordinate.length() == 0 {
      voy_fail("Voyager inverse-distance coordinates must not be empty")
    }
    for value in coordinate {
      if !voy_is_finite(value) {
        voy_fail("Voyager inverse-distance coordinates must be finite")
      }
    }
  }
  let limited = bandwidth > 0.0 && voy_is_finite(bandwidth)
  let neighbors : Array[Array[Int]] = []
  let raw : Array[Array[Double]] = []
  for i in 0..<n {
    let nb : Array[Int] = []
    let w : Array[Double] = []
    for j in 0..<n {
      if i != j {
        let d = voy_euclidean(coords[i], coords[j])
        if d > 0.0 && (!limited || d <= bandwidth) {
          nb.push(j)
          w.push(1.0 / @math.pow(d, power))
        }
      }
    }
    neighbors.push(nb)
    raw.push(w)
  }
  let styled = voy_apply_style(neighbors, raw, style)
  VoyagerWeights::{ style, n, neighbors, weights: styled, }
}

// ===========================================================================
// Global Moran's I
// ===========================================================================

///|
pub struct VoyagerGlobalMoranI {
  feature : String
  estimate : Double
  expectation : Double
  variance : Double
  z_score : Double
  p_value : Double
  n : Int
  s0 : Double
  s1 : Double
  s2 : Double
} derive(Debug)
pub extend VoyagerGlobalMoranI with @moonbitlang/core/debug.Debug::{to_repr}

///|
/// Internal Moran's I numerator and denominator for a given (permuted) vector.
fn voy_moran_components(
  values : Array[Double],
  weights : VoyagerWeights,
  mean : Double,
) -> (Double, Double) {
  let n = values.length()
  let centered : Array[Double] = []
  for index in 0..<n {
    centered.push(values[index] - mean)
  }
  let mut numerator = 0.0
  for i in 0..<n {
    let neighbors = weights.neighbors[i]
    let w_values = weights.weights[i]
    let mut lag = 0.0
    for k in 0..<neighbors.length() {
      lag = lag + w_values[k] * centered[neighbors[k]]
    }
    numerator = numerator + centered[i] * lag
  }
  let mut denominator = 0.0
  for index in 0..<n {
    denominator = denominator + centered[index] * centered[index]
  }
  (numerator, denominator)
}

///|
/// Global Moran's I for a single feature vector with the analytical
/// Cliff–Ord randomization expectation and variance.
pub fn voyager_global_morans_i(
  values : Array[Double],
  weights : VoyagerWeights,
  feature? : String = "feature",
) -> VoyagerGlobalMoranI raise VoyagerError {
  voy_validate_vector(values, "Moran's I")
  voy_validate_weights(weights, "Moran's I")
  if values.length() != weights.n {
    voy_fail("Voyager Moran's I values must match weights length")
  }
  let n = values.length()
  let mut mean = 0.0
  for value in values {
    mean = mean + value
  }
  mean = mean / n.to_double()
  let s0 = voyager_weights_s0(weights)
  if s0 <= 0.0 {
    voy_fail("Voyager Moran's I requires non-zero total weight S0")
  }
  let s1 = voyager_weights_s1(weights)
  let s2 = voyager_weights_s2(weights)
  let (numerator, denominator) = voy_moran_components(values, weights, mean)
  if denominator <= 0.0 {
    voy_fail("Voyager Moran's I requires non-zero variance in values")
  }
  let estimate = n.to_double() / s0 * (numerator / denominator)
  let expectation = -1.0 / (n.to_double() - 1.0)
  // Randomization variance (Cliff & Ord 1981)
  let mut m2 = 0.0
  let mut m4 = 0.0
  for value in values {
    let dev = value - mean
    m2 = m2 + dev * dev
    m4 = m4 + dev * dev * dev * dev
  }
  let b2 = n.to_double() * m4 / (m2 * m2)
  let nd = n.to_double()
  let n1 = nd - 1.0
  let n2 = nd - 2.0
  let n3 = nd - 3.0
  let a = nd * ((nd * nd - 3.0 * nd + 3.0) * s1 - nd * s2 + 3.0 * s0 * s0)
  let b = b2 * ((nd * nd - nd) * s1 - 2.0 * nd * s2 + 6.0 * s0 * s0)
  let variance = (a - b) / (n1 * n2 * n3 * s0 * s0) - expectation * expectation
  let safe_variance = if variance > 0.0 { variance } else { 0.0 }
  let z_score = if safe_variance > 0.0 {
    (estimate - expectation) / safe_variance.sqrt()
  } else {
    0.0
  }
  let p_value = voy_two_sided_p(z_score)
  VoyagerGlobalMoranI::{
    feature,
    estimate,
    expectation,
    variance: safe_variance,
    z_score,
    p_value,
    n,
    s0,
    s1,
    s2,
  }
}

// ===========================================================================
// Global Geary's c
// ===========================================================================

///|
pub struct VoyagerGlobalGearyC {
  feature : String
  estimate : Double
  expectation : Double
  variance : Double
  z_score : Double
  p_value : Double
  n : Int
  s0 : Double
  s1 : Double
  s2 : Double
} derive(Debug)
pub extend VoyagerGlobalGearyC with @moonbitlang/core/debug.Debug::{to_repr}

///|
/// Global Geary's c with the Cliff–Ord randomization expectation and variance.
pub fn voyager_global_gearys_c(
  values : Array[Double],
  weights : VoyagerWeights,
  feature? : String = "feature",
) -> VoyagerGlobalGearyC raise VoyagerError {
  voy_validate_vector(values, "Geary's c")
  voy_validate_weights(weights, "Geary's c")
  if values.length() != weights.n {
    voy_fail("Voyager Geary's c values must match weights length")
  }
  let n = values.length()
  let mut mean = 0.0
  for value in values {
    mean = mean + value
  }
  mean = mean / n.to_double()
  let s0 = voyager_weights_s0(weights)
  if s0 <= 0.0 {
    voy_fail("Voyager Geary's c requires non-zero total weight S0")
  }
  let s1 = voyager_weights_s1(weights)
  let s2 = voyager_weights_s2(weights)
  let mut cross = 0.0
  let mut ssd = 0.0
  for i in 0..<n {
    let dev_i = values[i] - mean
    ssd = ssd + dev_i * dev_i
    let neighbors = weights.neighbors[i]
    let w_values = weights.weights[i]
    for k in 0..<neighbors.length() {
      let diff = values[i] - values[neighbors[k]]
      cross = cross + w_values[k] * diff * diff
    }
  }
  if ssd <= 0.0 {
    voy_fail("Voyager Geary's c requires non-zero variance in values")
  }
  let nd = n.to_double()
  let estimate = n1_value(n) * cross / (2.0 * s0 * ssd)
  let expectation = 1.0
  let mut m2 = 0.0
  let mut m4 = 0.0
  for value in values {
    let dev = value - mean
    m2 = m2 + dev * dev
    m4 = m4 + dev * dev * dev * dev
  }
  let b2 = nd * m4 / (m2 * m2)
  let n1 = nd - 1.0
  let n2 = nd - 2.0
  let n3 = nd - 3.0
  // Cliff–Ord (1981) randomization variance for Geary's c, as implemented in
  // spdep::geary.test.
  let term1 = n1 * s1 * (nd * nd - 3.0 * nd + 3.0 - n1 * b2)
  let term2 = 0.25 *
    n1 *
    (nd * nd + 3.0 * nd - 6.0 - (nd * nd - nd + 2.0) * b2) *
    s2
  let term3 = s0 * s0 * (nd * nd - 3.0 - n1 * n1 * b2)
  let variance = (term1 - term2 + term3) / (s0 * s0 * nd * n2 * n3)
  let safe_variance = if variance > 0.0 { variance } else { 0.0 }
  let z_score = if safe_variance > 0.0 {
    (estimate - expectation) / safe_variance.sqrt()
  } else {
    0.0
  }
  let p_value = voy_two_sided_p(z_score)
  VoyagerGlobalGearyC::{
    feature,
    estimate,
    expectation,
    variance: safe_variance,
    z_score,
    p_value,
    n,
    s0,
    s1,
    s2,
  }
}

///|
fn n1_value(n : Int) -> Double {
  (n - 1).to_double()
}

// ===========================================================================
// Permutation helpers shared by local_res statistics
// ===========================================================================

///|
fn voy_perm_p_value_two_sided(
  observed : Double,
  expectation : Double,
  permutations : Array[Double],
) -> (Double, Double, Double) {
  // Returns (p_value, perm_mean, perm_sd).
  let count = permutations.length()
  if count == 0 {
    return (1.0, observed, 0.0)
  }
  let mut sum = 0.0
  let mut sum_sq = 0.0
  for value in permutations {
    sum = sum + value
    sum_sq = sum_sq + value * value
  }
  let mean = sum / count.to_double()
  let var_sample = (sum_sq - sum * sum / count.to_double()) / count.to_double()
  let sd = if var_sample > 0.0 { var_sample.sqrt() } else { 0.0 }
  let obs_dev = (observed - expectation).abs()
  let mut extreme = 0
  for value in permutations {
    if (value - expectation).abs() >= obs_dev {
      extreme = extreme + 1
    }
  }
  let p = (extreme.to_double() + 1.0) / (count.to_double() + 1.0)
  (p, mean, sd)
}

// ===========================================================================
// Local Moran's I (LISA)
// ===========================================================================

///|
pub struct VoyagerLocalMoranI {
  feature : String
  local_i : Array[Double]
  expectations : Array[Double]
  variances : Array[Double]
  z_scores : Array[Double]
  p_values : Array[Double]
  fdr : Array[Double]
  quadrants : Array[String]
  perm_mean : Array[Double]
  perm_sd : Array[Double]
  n : Int
  permutations : Int
  style : String
} derive(Debug)
pub extend VoyagerLocalMoranI with @moonbitlang/core/debug.Debug::{to_repr}

///|
/// Quadrant classification from the Moran scatter plot. The spatial lag is the
/// row-standardized weighted sum of the neighbours' standardized values.
fn voy_moran_quadrant(
  z_i : Double,
  lag_i : Double,
  z_significant : Bool,
) -> String {
  if !z_significant {
    return "not significant"
  }
  if z_i > 0.0 && lag_i > 0.0 {
    "HH"
  } else if z_i < 0.0 && lag_i < 0.0 {
    "LL"
  } else if z_i > 0.0 && lag_i < 0.0 {
    "HL"
  } else {
    "LH"
  }
}

///|
fn voy_local_morans_point(
  values : Array[Double],
  weights : VoyagerWeights,
  mean : Double,
  m2 : Double,
) -> Array[Double] {
  let n = values.length()
  let centered : Array[Double] = []
  for index in 0..<n {
    centered.push(values[index] - mean)
  }
  let output : Array[Double] = Array::make(n, 0.0)
  if m2 <= 0.0 {
    return output
  }
  for i in 0..<n {
    let neighbors = weights.neighbors[i]
    let w_values = weights.weights[i]
    let mut lag = 0.0
    for k in 0..<neighbors.length() {
      lag = lag + w_values[k] * centered[neighbors[k]]
    }
    output[i] = centered[i] / m2 * lag
  }
  output
}

///|
/// Local Moran's I (Anselin 1995) for a single feature. Permutation inference
/// permutes the residual values across spots `permutations` times and computes
/// a pseudo p-value and per-location z-score. `fdr_threshold` controls BH
/// significance reporting via the `quadrants` field.
pub fn voyager_local_morans_i(
  values : Array[Double],
  weights : VoyagerWeights,
  permutations? : Int = 0,
  seed? : Int = 20240501,
  fdr_threshold? : Double = 0.05,
  feature? : String = "feature",
) -> VoyagerLocalMoranI raise VoyagerError {
  voy_validate_vector(values, "local_res Moran's I")
  voy_validate_weights(weights, "local_res Moran's I")
  if values.length() != weights.n {
    voy_fail("Voyager local_res Moran's I values must match weights length")
  }
  if permutations < 0 {
    voy_fail("Voyager local_res Moran's I permutations must be non-negative")
  }
  if !voy_is_finite(fdr_threshold) ||
    fdr_threshold <= 0.0 ||
    fdr_threshold > 1.0 {
    voy_fail("Voyager local_res Moran's I FDR threshold must be in (0, 1]")
  }
  let n = values.length()
  let mut mean = 0.0
  for value in values {
    mean = mean + value
  }
  mean = mean / n.to_double()
  let mut m2 = 0.0
  for value in values {
    let dev = value - mean
    m2 = m2 + dev * dev
  }
  m2 = m2 / n.to_double()
  if m2 <= 0.0 {
    voy_fail("Voyager local_res Moran's I requires non-zero variance in values")
  }
  let observed = voy_local_morans_point(values, weights, mean, m2)
  // Standardized values and spatial lag for quadrant classification.
  let z_scores_std : Array[Double] = []
  for value in values {
    z_scores_std.push((value - mean) / m2.sqrt())
  }
  let lag_std : Array[Double] = Array::make(n, 0.0)
  for i in 0..<n {
    let neighbors = weights.neighbors[i]
    let w_values = weights.weights[i]
    let mut s = 0.0
    for k in 0..<neighbors.length() {
      s = s + w_values[k] * z_scores_std[neighbors[k]]
    }
    lag_std[i] = s
  }
  let expectations : Array[Double] = Array::make(n, 0.0)
  let variances : Array[Double] = Array::make(n, 0.0)
  let z_scores : Array[Double] = Array::make(n, 0.0)
  let p_values : Array[Double] = Array::make(n, 1.0)
  let perm_mean : Array[Double] = Array::make(n, observed[0])
  let perm_sd : Array[Double] = Array::make(n, 0.0)
  for i in 0..<n {
    perm_mean[i] = observed[i]
  }
  let row_w : Array[Double] = Array::make(n, 0.0)
  let row_w_sq : Array[Double] = Array::make(n, 0.0)
  for i in 0..<n {
    let mut ws = 0.0
    let mut wsq = 0.0
    for value in weights.weights[i] {
      ws = ws + value
      wsq = wsq + value * value
    }
    row_w[i] = ws
    row_w_sq[i] = wsq
  }
  if permutations > 0 {
    let rng = VoyagerRng::new(seed)
    // perm_local[k][i] = local_res Moran's I at spot i under permutation k.
    let perm_local : Array[Array[Double]] = []
    for _ in 0..<permutations {
      let permuted = voy_permute(values, rng)
      let permuted_local = voy_local_morans_point(permuted, weights, mean, m2)
      perm_local.push(permuted_local)
    }
    for i in 0..<n {
      let exp_i = -row_w[i] / (n.to_double() - 1.0)
      expectations[i] = exp_i
      let permut_values : Array[Double] = []
      for k in 0..<permutations {
        permut_values.push(perm_local[k][i])
      }
      let (p, pmean, psd) = voy_perm_p_value_two_sided(
        observed[i],
        exp_i,
        permut_values,
      )
      p_values[i] = p
      perm_mean[i] = pmean
      perm_sd[i] = psd
      variances[i] = psd * psd
      z_scores[i] = if psd > 0.0 { (observed[i] - pmean) / psd } else { 0.0 }
    }
  } else {
    // Analytical approximation: use permutation-free variance under the
    // randomization assumption for reference only.
    for i in 0..<n {
      let exp_i = -row_w[i] / (n.to_double() - 1.0)
      expectations[i] = exp_i
      let nd = n.to_double()
      let n1 = nd - 1.0
      let var_i = (row_w_sq[i] * nd - row_w[i] * row_w[i]) /
        (n1 * (nd - 2.0).max(1.0))
      let safe_var = if var_i > 0.0 { var_i } else { 0.0 }
      variances[i] = safe_var
      z_scores[i] = if safe_var > 0.0 {
        (observed[i] - exp_i) / safe_var.sqrt()
      } else {
        0.0
      }
      p_values[i] = voy_two_sided_p(z_scores[i])
      perm_mean[i] = exp_i
      perm_sd[i] = safe_var.sqrt()
    }
  }
  let fdr = voy_bh_adjust(p_values)
  let quadrants : Array[String] = []
  for i in 0..<n {
    quadrants.push(
      voy_moran_quadrant(z_scores_std[i], lag_std[i], fdr[i] <= fdr_threshold),
    )
  }
  VoyagerLocalMoranI::{
    feature,
    local_i: observed,
    expectations,
    variances,
    z_scores,
    p_values,
    fdr,
    quadrants,
    perm_mean,
    perm_sd,
    n,
    permutations,
    style: weights.style,
  }
}

// ===========================================================================
// Local Geary's c
// ===========================================================================

///|
pub struct VoyagerLocalGearyC {
  feature : String
  local_c : Array[Double]
  expectations : Array[Double]
  z_scores : Array[Double]
  p_values : Array[Double]
  fdr : Array[Double]
  classifications : Array[String]
  n : Int
  permutations : Int
} derive(Debug)
pub extend VoyagerLocalGearyC with @moonbitlang/core/debug.Debug::{to_repr}

///|
fn voy_local_gearys_point(
  values : Array[Double],
  weights : VoyagerWeights,
  m2 : Double,
) -> Array[Double] {
  let n = values.length()
  let output : Array[Double] = Array::make(n, 0.0)
  if m2 <= 0.0 {
    return output
  }
  for i in 0..<n {
    let neighbors = weights.neighbors[i]
    let w_values = weights.weights[i]
    let mut acc = 0.0
    for k in 0..<neighbors.length() {
      let diff = values[i] - values[neighbors[k]]
      acc = acc + w_values[k] * diff * diff
    }
    output[i] = (n - 1).to_double() / m2 * acc
  }
  output
}

///|
/// Local Geary's c (Anselin 1995) with permutation inference. The
/// `classifications` field flags "hotspot-like" (low c, similar neighbours)
/// and "outlier-like" (high c, dissimilar neighbours) spots after BH-FDR.
pub fn voyager_local_gearys_c(
  values : Array[Double],
  weights : VoyagerWeights,
  permutations? : Int = 0,
  seed? : Int = 20240501,
  fdr_threshold? : Double = 0.05,
  feature? : String = "feature",
) -> VoyagerLocalGearyC raise VoyagerError {
  voy_validate_vector(values, "local_res Geary's c")
  voy_validate_weights(weights, "local_res Geary's c")
  if values.length() != weights.n {
    voy_fail("Voyager local_res Geary's c values must match weights length")
  }
  if permutations < 0 {
    voy_fail("Voyager local_res Geary's c permutations must be non-negative")
  }
  if !voy_is_finite(fdr_threshold) ||
    fdr_threshold <= 0.0 ||
    fdr_threshold > 1.0 {
    voy_fail("Voyager local_res Geary's c FDR threshold must be in (0, 1]")
  }
  let n = values.length()
  let mut mean = 0.0
  for value in values {
    mean = mean + value
  }
  mean = mean / n.to_double()
  let mut m2 = 0.0
  for value in values {
    let dev = value - mean
    m2 = m2 + dev * dev
  }
  m2 = m2 / n.to_double()
  if m2 <= 0.0 {
    voy_fail("Voyager local_res Geary's c requires non-zero variance in values")
  }
  let observed = voy_local_gearys_point(values, weights, m2)
  let expectations : Array[Double] = Array::make(n, 1.0)
  let z_scores : Array[Double] = Array::make(n, 0.0)
  let p_values : Array[Double] = Array::make(n, 1.0)
  if permutations > 0 {
    let rng = VoyagerRng::new(seed)
    let perm_local : Array[Array[Double]] = []
    for _ in 0..<permutations {
      let permuted = voy_permute(values, rng)
      perm_local.push(voy_local_gearys_point(permuted, weights, m2))
    }
    for i in 0..<n {
      let permut_values : Array[Double] = []
      for k in 0..<permutations {
        permut_values.push(perm_local[k][i])
      }
      let (p, pmean, psd) = voy_perm_p_value_two_sided(
        observed[i],
        1.0,
        permut_values,
      )
      p_values[i] = p
      z_scores[i] = if psd > 0.0 { (observed[i] - pmean) / psd } else { 0.0 }
    }
  } else {
    for i in 0..<n {
      // Without permutations there is no closed-form randomization variance
      // commonly used for local_res Geary; report only the point estimate and a
      // non-significant classification.
      p_values[i] = 1.0
      z_scores[i] = 0.0
    }
  }
  let fdr = voy_bh_adjust(p_values)
  let classifications : Array[String] = []
  for i in 0..<n {
    if fdr[i] <= fdr_threshold && observed[i] < 1.0 {
      classifications.push("similar")
    } else if fdr[i] <= fdr_threshold && observed[i] > 1.0 {
      classifications.push("dissimilar")
    } else {
      classifications.push("not significant")
    }
  }
  VoyagerLocalGearyC::{
    feature,
    local_c: observed,
    expectations,
    z_scores,
    p_values,
    fdr,
    classifications,
    n,
    permutations,
  }
}

// ===========================================================================
// Local Getis–Ord Gi / Gi*
// ===========================================================================

///|
pub struct VoyagerLocalGetisOrd {
  feature : String
  statistic : Array[Double]
  expectations : Array[Double]
  z_scores : Array[Double]
  p_values : Array[Double]
  fdr : Array[Double]
  classifications : Array[String]
  star : Bool
  n : Int
} derive(Debug)
pub extend VoyagerLocalGetisOrd with @moonbitlang/core/debug.Debug::{to_repr}

///|
/// Local Getis–Ord Gi (excludes self) or Gi* (includes self with w_ii = 1).
/// Significance uses the exact Ord–Getis randomization z-score for the
/// weighted sum.
pub fn voyager_local_getis_ord(
  values : Array[Double],
  weights : VoyagerWeights,
  star? : Bool = true,
  fdr_threshold? : Double = 0.05,
  feature? : String = "feature",
) -> VoyagerLocalGetisOrd raise VoyagerError {
  voy_validate_vector(values, "Getis–Ord")
  voy_validate_weights(weights, "Getis–Ord")
  if values.length() != weights.n {
    voy_fail("Voyager Getis–Ord values must match weights length")
  }
  if !voy_is_finite(fdr_threshold) ||
    fdr_threshold <= 0.0 ||
    fdr_threshold > 1.0 {
    voy_fail("Voyager Getis–Ord FDR threshold must be in (0, 1]")
  }
  let n = values.length()
  let mut total = 0.0
  let mut total_sq = 0.0
  for value in values {
    if value < 0.0 {
      voy_fail("Voyager Getis–Ord requires non-negative values")
    }
    total = total + value
    total_sq = total_sq + value * value
  }
  let mean = total / n.to_double()
  // Sample variance with denominator (n-1).
  let sample_var = (total_sq - total * total / n.to_double()) /
    (n.to_double() - 1.0)
  let safe_var = if sample_var > 0.0 { sample_var } else { 0.0 }
  let statistic : Array[Double] = Array::make(n, 0.0)
  let expectations : Array[Double] = Array::make(n, 0.0)
  let z_scores : Array[Double] = Array::make(n, 0.0)
  let p_values : Array[Double] = Array::make(n, 1.0)
  for i in 0..<n {
    let neighbors = weights.neighbors[i]
    let w_values = weights.weights[i]
    let mut weighted_sum = 0.0
    let mut w_i = 0.0
    let mut s1_i = 0.0
    for k in 0..<neighbors.length() {
      let j = neighbors[k]
      weighted_sum = weighted_sum + w_values[k] * values[j]
      w_i = w_i + w_values[k]
      s1_i = s1_i + w_values[k] * w_values[k]
    }
    let mut numerator_sum = weighted_sum
    let mut total_w = w_i
    let mut total_w_sq = s1_i
    if star {
      numerator_sum = numerator_sum + values[i]
      total_w = total_w + 1.0
      total_w_sq = total_w_sq + 1.0
    }
    let y = if star { total } else { total - values[i] }
    let n_used = if star { n.to_double() } else { (n - 1).to_double() }
    let stat = if y > 0.0 { numerator_sum / y } else { 0.0 }
    statistic[i] = stat
    expectations[i] = total_w / n_used
    let denom_root = if safe_var > 0.0 {
      let var_weighted = (n_used * total_w_sq - total_w * total_w) /
        (n_used - 1.0).max(1.0)
      if var_weighted > 0.0 {
        safe_var.sqrt() * var_weighted.sqrt()
      } else {
        0.0
      }
    } else {
      0.0
    }
    let z = if denom_root > 0.0 {
      (numerator_sum - total_w * mean) / denom_root
    } else {
      0.0
    }
    z_scores[i] = z
    p_values[i] = voy_two_sided_p(z)
  }
  let fdr = voy_bh_adjust(p_values)
  let classifications : Array[String] = []
  for i in 0..<n {
    if fdr[i] <= fdr_threshold && z_scores[i] > 0.0 {
      classifications.push("hotspot")
    } else if fdr[i] <= fdr_threshold && z_scores[i] < 0.0 {
      classifications.push("coldspot")
    } else {
      classifications.push("not significant")
    }
  }
  VoyagerLocalGetisOrd::{
    feature,
    statistic,
    expectations,
    z_scores,
    p_values,
    fdr,
    classifications,
    star,
    n,
  }
}

// ===========================================================================
// Lee's L (bivariate spatial association)
// ===========================================================================

///|
pub struct VoyagerLeesL {
  feature_x : String
  feature_y : String
  global_l : Double
  local_l : Array[Double]
  z_scores : Array[Double]
  p_values : Array[Double]
  fdr : Array[Double]
  n : Int
  s0 : Double
} derive(Debug)
pub extend VoyagerLeesL with @moonbitlang/core/debug.Debug::{to_repr}

///|
/// Lee's L (Lee 2001) bivariate spatial association. Returns the global L
/// (range approximately [-1, 1] for row-standardized weights) and per-location
/// local_res L values with permutation-based inference.
pub fn voyager_lees_l(
  values_x : Array[Double],
  values_y : Array[Double],
  weights : VoyagerWeights,
  permutations? : Int = 0,
  seed? : Int = 20240501,
  fdr_threshold? : Double = 0.05,
  feature_x? : String = "feature_x",
  feature_y? : String = "feature_y",
) -> VoyagerLeesL raise VoyagerError {
  voy_validate_vector(values_x, "Lee's L x")
  voy_validate_vector(values_y, "Lee's L y")
  voy_validate_weights(weights, "Lee's L")
  voy_validate_equal_length(values_x, values_y, "Lee's L")
  if values_x.length() != weights.n {
    voy_fail("Voyager Lee's L values must match weights length")
  }
  if !voy_is_finite(fdr_threshold) ||
    fdr_threshold <= 0.0 ||
    fdr_threshold > 1.0 {
    voy_fail("Voyager Lee's L FDR threshold must be in (0, 1]")
  }
  let n = values_x.length()
  let mut mean_x = 0.0
  let mut mean_y = 0.0
  for index in 0..<n {
    mean_x = mean_x + values_x[index]
    mean_y = mean_y + values_y[index]
  }
  mean_x = mean_x / n.to_double()
  mean_y = mean_y / n.to_double()
  let cx : Array[Double] = []
  let cy : Array[Double] = []
  let mut ss_x = 0.0
  let mut ss_y = 0.0
  for index in 0..<n {
    let dx = values_x[index] - mean_x
    let dy = values_y[index] - mean_y
    cx.push(dx)
    cy.push(dy)
    ss_x = ss_x + dx * dx
    ss_y = ss_y + dy * dy
  }
  if ss_x <= 0.0 || ss_y <= 0.0 {
    voy_fail("Voyager Lee's L requires non-zero variance in both features")
  }
  let s0 = voyager_weights_s0(weights)
  if s0 <= 0.0 {
    voy_fail("Voyager Lee's L requires non-zero total weight S0")
  }
  let lag_x : Array[Double] = Array::make(n, 0.0)
  let lag_y : Array[Double] = Array::make(n, 0.0)
  for i in 0..<n {
    let neighbors = weights.neighbors[i]
    let w_values = weights.weights[i]
    let mut lx = 0.0
    let mut ly = 0.0
    for k in 0..<neighbors.length() {
      lx = lx + w_values[k] * cx[neighbors[k]]
      ly = ly + w_values[k] * cy[neighbors[k]]
    }
    lag_x[i] = lx
    lag_y[i] = ly
  }
  let mut cross = 0.0
  let local_l : Array[Double] = Array::make(n, 0.0)
  for i in 0..<n {
    cross = cross + lag_x[i] * lag_y[i]
  }
  let denom = ss_x.sqrt() * ss_y.sqrt()
  let global_l = n.to_double() / (s0 * s0) * cross / denom
  let scale_local = n.to_double() / (s0 * s0) / denom
  for i in 0..<n {
    local_l[i] = scale_local * lag_x[i] * lag_y[i]
  }
  let z_scores : Array[Double] = Array::make(n, 0.0)
  let p_values : Array[Double] = Array::make(n, 1.0)
  if permutations > 0 {
    let rng = VoyagerRng::new(seed)
    let perm_local : Array[Array[Double]] = []
    for _ in 0..<permutations {
      let perm_y = voy_permute(cy, rng)
      let lag_py : Array[Double] = Array::make(n, 0.0)
      for i in 0..<n {
        let neighbors = weights.neighbors[i]
        let w_values = weights.weights[i]
        let mut s = 0.0
        for k in 0..<neighbors.length() {
          s = s + w_values[k] * perm_y[neighbors[k]]
        }
        lag_py[i] = s
      }
      let row : Array[Double] = Array::make(n, 0.0)
      for i in 0..<n {
        row[i] = scale_local * lag_x[i] * lag_py[i]
      }
      perm_local.push(row)
    }
    for i in 0..<n {
      let permut_values : Array[Double] = []
      for k in 0..<permutations {
        permut_values.push(perm_local[k][i])
      }
      let mut sum = 0.0
      let mut sum_sq = 0.0
      for value in permut_values {
        sum = sum + value
        sum_sq = sum_sq + value * value
      }
      let pmean = sum / permut_values.length().to_double()
      let pvar = (sum_sq - sum * sum / permut_values.length().to_double()) /
        permut_values.length().to_double()
      let psd = if pvar > 0.0 { pvar.sqrt() } else { 0.0 }
      let mut extreme = 0
      for value in permut_values {
        if (value - pmean).abs() >= (local_l[i] - pmean).abs() {
          extreme = extreme + 1
        }
      }
      p_values[i] = (extreme.to_double() + 1.0) /
        (permut_values.length().to_double() + 1.0)
      z_scores[i] = if psd > 0.0 { (local_l[i] - pmean) / psd } else { 0.0 }
    }
  }
  let fdr = voy_bh_adjust(p_values)
  VoyagerLeesL::{
    feature_x,
    feature_y,
    global_l,
    local_l,
    z_scores,
    p_values,
    fdr,
    n,
    s0,
  }
}

// ===========================================================================
// Multivariate local_res Geary
// ===========================================================================

///|
pub struct VoyagerMultivariateLocalGeary {
  feature_set : Array[String]
  local_statistic : Array[Double]
  expectations : Array[Double]
  z_scores : Array[Double]
  p_values : Array[Double]
  fdr : Array[Double]
  classifications : Array[String]
  n : Int
  n_features : Int
  permutations : Int
} derive(Debug)
pub extend VoyagerMultivariateLocalGeary with @moonbitlang/core/debug.Debug::{to_repr}

///|
fn voy_multivariate_local_geary_point(
  matrix : Array[Array[Double]],
  weights : VoyagerWeights,
  means : Array[Double],
) -> Array[Double] {
  // matrix[feature][spot]. Returns per-spot sum over features of weighted
  // squared feature-space distances to neighbours.
  let n_features = matrix.length()
  let n = if n_features > 0 { matrix[0].length() } else { 0 }
  let output : Array[Double] = Array::make(n, 0.0)
  for i in 0..<n {
    let neighbors = weights.neighbors[i]
    let w_values = weights.weights[i]
    let mut acc = 0.0
    for k in 0..<neighbors.length() {
      let j = neighbors[k]
      let mut dist_sq = 0.0
      for feature in 0..<n_features {
        let diff = matrix[feature][i] -
          means[feature] -
          (matrix[feature][j] - means[feature])
        dist_sq = dist_sq + diff * diff
      }
      acc = acc + w_values[k] * dist_sq
    }
    output[i] = acc
  }
  output
}

///|
/// Multivariate local_res Geary (co-located). Combines multiple features into a
/// single per-spot statistic measuring the average weighted squared distance
/// in feature space to the spot's neighbours.
pub fn voyager_multivariate_local_geary(
  matrix : Array[Array[Double]],
  weights : VoyagerWeights,
  permutations? : Int = 0,
  seed? : Int = 20240501,
  fdr_threshold? : Double = 0.05,
  feature_names? : Array[String] = [],
) -> VoyagerMultivariateLocalGeary raise VoyagerError {
  if matrix.length() == 0 {
    voy_fail(
      "Voyager multivariate local_res Geary requires at least one feature",
    )
  }
  let n_features = matrix.length()
  let n = matrix[0].length()
  if n < 3 {
    voy_fail("Voyager multivariate local_res Geary requires at least 3 spots")
  }
  for feature in 0..<n_features {
    if matrix[feature].length() != n {
      voy_fail(
        "Voyager multivariate local_res Geary features must share spot count",
      )
    }
    for value in matrix[feature] {
      if !voy_is_finite(value) {
        voy_fail("Voyager multivariate local_res Geary values must be finite")
      }
    }
  }
  voy_validate_weights(weights, "multivariate local_res Geary")
  if n != weights.n {
    voy_fail(
      "Voyager multivariate local_res Geary matrix must match weights length",
    )
  }
  let means : Array[Double] = Array::make(n_features, 0.0)
  for feature in 0..<n_features {
    let mut s = 0.0
    for value in matrix[feature] {
      s = s + value
    }
    means[feature] = s / n.to_double()
  }
  let observed = voy_multivariate_local_geary_point(matrix, weights, means)
  let expectations : Array[Double] = Array::make(n, 0.0)
  let z_scores : Array[Double] = Array::make(n, 0.0)
  let p_values : Array[Double] = Array::make(n, 1.0)
  if permutations > 0 {
    let rng = VoyagerRng::new(seed)
    let perm_local : Array[Array[Double]] = []
    for _ in 0..<permutations {
      let permuted : Array[Array[Double]] = []
      for feature in 0..<n_features {
        permuted.push(voy_permute(matrix[feature], rng))
      }
      perm_local.push(
        voy_multivariate_local_geary_point(permuted, weights, means),
      )
    }
    for i in 0..<n {
      let permut_values : Array[Double] = []
      for k in 0..<permutations {
        permut_values.push(perm_local[k][i])
      }
      let mut sum = 0.0
      let mut sum_sq = 0.0
      for value in permut_values {
        sum = sum + value
        sum_sq = sum_sq + value * value
      }
      let pmean = sum / permut_values.length().to_double()
      let pvar = (sum_sq - sum * sum / permut_values.length().to_double()) /
        permut_values.length().to_double()
      let psd = if pvar > 0.0 { pvar.sqrt() } else { 0.0 }
      let mut extreme = 0
      for value in permut_values {
        if (value - pmean).abs() >= (observed[i] - pmean).abs() {
          extreme = extreme + 1
        }
      }
      p_values[i] = (extreme.to_double() + 1.0) /
        (permut_values.length().to_double() + 1.0)
      z_scores[i] = if psd > 0.0 { (observed[i] - pmean) / psd } else { 0.0 }
    }
  }
  let fdr = voy_bh_adjust(p_values)
  let names = if feature_names.length() == n_features {
    feature_names
  } else {
    voy_default_names("feature", n_features)
  }
  let classifications : Array[String] = []
  for i in 0..<n {
    if fdr[i] <= fdr_threshold && z_scores[i] < 0.0 {
      classifications.push("similar")
    } else if fdr[i] <= fdr_threshold && z_scores[i] > 0.0 {
      classifications.push("dissimilar")
    } else {
      classifications.push("not significant")
    }
  }
  VoyagerMultivariateLocalGeary::{
    feature_set: names,
    local_statistic: observed,
    expectations,
    z_scores,
    p_values,
    fdr,
    classifications,
    n,
    n_features,
    permutations,
  }
}

// ===========================================================================
// Empirical variogram and model fitting
// ===========================================================================

///|
pub struct VoyagerVariogramPoint {
  lag : Double
  gamma : Double
  npairs : Int
} derive(Debug)
pub extend VoyagerVariogramPoint with @moonbitlang/core/debug.Debug::{to_repr}

///|
pub struct VoyagerVariogramModel {
  model_type : String
  nugget : Double
  sill : Double
  range : Double
  fitted_sse : Double
  n_points : Int
} derive(Debug)
pub extend VoyagerVariogramModel with @moonbitlang/core/debug.Debug::{to_repr}

///|
/// Empirical (sample) variogram. Pairs of spots are binned into `n_lags`
/// equidistant bins from 0 to `max_lag` (default: half the maximum pairwise
/// distance). The semivariance at each bin is half the mean squared
/// difference of values for the pairs assigned to that bin.
pub fn voyager_empirical_variogram(
  coords : Array[Array[Double]],
  values : Array[Double],
  n_lags? : Int = 15,
  max_lag? : Double = 0.0,
) -> Array[VoyagerVariogramPoint] raise VoyagerError {
  if coords.length() < 3 {
    voy_fail("Voyager variogram requires at least 3 spots")
  }
  if values.length() != coords.length() {
    voy_fail("Voyager variogram values must match coordinates length")
  }
  if n_lags < 1 {
    voy_fail("Voyager variogram n_lags must be positive")
  }
  for value in values {
    if !voy_is_finite(value) {
      voy_fail("Voyager variogram values must be finite")
    }
  }
  let n = coords.length()
  let distances = voy_pairwise_distances(coords)
  let mut max_dist = 0.0
  for i in 0..<n {
    for j in (i + 1)..<n {
      if distances[i][j] > max_dist {
        max_dist = distances[i][j]
      }
    }
  }
  if max_dist <= 0.0 {
    voy_fail("Voyager variogram requires non-coincident coordinates")
  }
  let upper = if max_lag > 0.0 && voy_is_finite(max_lag) {
    max_lag
  } else {
    max_dist / 2.0
  }
  let lag_width = upper / n_lags.to_double()
  if lag_width <= 0.0 {
    voy_fail("Voyager variogram lag width must be positive")
  }
  let sums : Array[Double] = Array::make(n_lags, 0.0)
  let counts : Array[Int] = Array::make(n_lags, 0)
  for i in 0..<n {
    for j in (i + 1)..<n {
      let d = distances[i][j]
      if d > 0.0 && d <= upper {
        let bin = (d / lag_width).to_int().min(n_lags - 1)
        let diff = values[i] - values[j]
        sums[bin] = sums[bin] + diff * diff
        counts[bin] = counts[bin] + 1
      }
    }
  }
  let output : Array[VoyagerVariogramPoint] = []
  for bin in 0..<n_lags {
    if counts[bin] > 0 {
      let lag = (bin.to_double() + 0.5) * lag_width
      let gamma = sums[bin] / (2.0 * counts[bin].to_double())
      output.push(VoyagerVariogramPoint::{ lag, gamma, npairs: counts[bin], })
    }
  }
  if output.length() == 0 {
    voy_fail("Voyager variogram produced no populated bins")
  }
  output
}

///|
/// Theoretical variogram shape function (excluding nugget). Returns the
/// bounded partial sill component for distance h.
fn voy_variogram_shape(
  model_type : String,
  h : Double,
  range : Double,
) -> Double {
  if range <= 0.0 {
    return 0.0
  }
  let ratio = h / range
  match model_type {
    "spherical" =>
      if h >= range {
        1.0
      } else {
        1.5 * ratio - 0.5 * ratio * ratio * ratio
      }
    "exponential" => 1.0 - @math.exp(-3.0 * ratio)
    "gaussian" => 1.0 - @math.exp(-3.0 * ratio * ratio)
    _ => 0.0
  }
}

///|
/// Fit a theoretical variogram model (spherical, exponential or gaussian) to
/// an empirical variogram by minimizing the sum of squared residuals. A
/// bounded grid search over `range` is combined with closed-form linear
/// least squares for (nugget, partial sill) at each candidate range.
pub fn voyager_fit_variogram(
  empirical : Array[VoyagerVariogramPoint],
  model_type? : String = "spherical",
  nugget? : Double = 0.0,
  fix_nugget? : Bool = false,
  range_grid? : Int = 40,
) -> VoyagerVariogramModel raise VoyagerError {
  if empirical.length() < 2 {
    voy_fail("Voyager variogram fit requires at least 2 empirical points")
  }
  let model : String = match model_type {
    "spherical" => "spherical"
    "exponential" => "exponential"
    "gaussian" => "gaussian"
    _ =>
      raise VoyagerError::VoyagerError(
        "Voyager variogram model must be spherical/exponential/gaussian",
      )
  }
  if !voy_is_finite(nugget) || nugget < 0.0 {
    voy_fail("Voyager variogram nugget must be non-negative and finite")
  }
  if range_grid < 4 {
    voy_fail("Voyager variogram range grid must have at least 4 points")
  }
  let mut max_lag = 0.0
  for point in empirical {
    if point.lag > max_lag {
      max_lag = point.lag
    }
  }
  if max_lag <= 0.0 {
    voy_fail("Voyager variogram empirical lags must be positive")
  }
  // Precompute shape values for each (point, candidate range) — but to keep
  // the implementation simple we iterate ranges and solve the linear LSQ.
  let mut best_sse = 1.0e300
  let mut best_nugget = if fix_nugget { nugget } else { 0.0 }
  let mut best_sill = 0.0
  let mut best_range = max_lag / 3.0
  let grid = range_grid
  let mut gi = 1
  while gi <= grid {
    let candidate_range = max_lag * gi.to_double() / grid.to_double()
    if candidate_range <= 0.0 {
      gi = gi + 1
      continue
    }
    // Linear least squares for y = nugget + partial_sill * shape(h).
    // If nugget is fixed, this is a one-parameter least squares.
    let mut sum_x = 0.0
    let mut sum_y = 0.0
    let mut sum_xx = 0.0
    let mut sum_xy = 0.0
    let mut count = 0.0
    for point in empirical {
      let shape = voy_variogram_shape(model, point.lag, candidate_range)
      let y = point.gamma
      if fix_nugget {
        // y - nugget = partial_sill * shape  -> one-parameter LSQ.
        sum_x = sum_x + shape
        sum_y = sum_y + (y - nugget)
        sum_xx = sum_xx + shape * shape
        sum_xy = sum_xy + shape * (y - nugget)
      } else {
        sum_x = sum_x + shape
        sum_y = sum_y + y
        sum_xx = sum_xx + shape * shape
        sum_xy = sum_xy + shape * y
      }
      count = count + 1.0
    }
    let (fit_nugget, fit_partial) = if fix_nugget {
      let denom = sum_xx
      let partial = if denom > 0.0 { sum_xy / denom } else { 0.0 }
      (nugget, partial)
    } else {
      // Two-parameter LSQ: y = a + b * shape, a = nugget, b = partial_sill.
      let denom = count * sum_xx - sum_x * sum_x
      if denom.abs() < 1.0e-300 {
        (0.0, 0.0)
      } else {
        let b = (count * sum_xy - sum_x * sum_y) / denom
        let a = (sum_y - b * sum_x) / count
        (a, b)
      }
    }
    // Enforce non-negative nugget and partial sill (physical constraint).
    let safe_nugget = if fit_nugget < 0.0 { 0.0 } else { fit_nugget }
    let safe_partial = if fit_partial < 0.0 { 0.0 } else { fit_partial }
    let mut sse = 0.0
    for point in empirical {
      let shape = voy_variogram_shape(model, point.lag, candidate_range)
      let predicted = safe_nugget + safe_partial * shape
      let residual = point.gamma - predicted
      sse = sse + residual * residual
    }
    if sse < best_sse {
      best_sse = sse
      best_nugget = safe_nugget
      best_sill = safe_partial
      best_range = candidate_range
    }
    gi = gi + 1
  }
  // Refine around the best range with a finer grid.
  let lo = (best_range * 0.5).max(max_lag / grid.to_double())
  let hi = (best_range * 1.5).min(max_lag)
  let fine = 20
  let mut fi = 0
  while fi <= fine {
    let candidate_range = lo + (hi - lo) * fi.to_double() / fine.to_double()
    if candidate_range <= 0.0 {
      fi = fi + 1
      continue
    }
    let mut sum_x = 0.0
    let mut sum_y = 0.0
    let mut sum_xx = 0.0
    let mut sum_xy = 0.0
    let mut count = 0.0
    for point in empirical {
      let shape = voy_variogram_shape(model, point.lag, candidate_range)
      let y = point.gamma
      if fix_nugget {
        sum_x = sum_x + shape
        sum_y = sum_y + (y - nugget)
        sum_xx = sum_xx + shape * shape
        sum_xy = sum_xy + shape * (y - nugget)
      } else {
        sum_x = sum_x + shape
        sum_y = sum_y + y
        sum_xx = sum_xx + shape * shape
        sum_xy = sum_xy + shape * y
      }
      count = count + 1.0
    }
    let (fit_nugget, fit_partial) = if fix_nugget {
      let denom = sum_xx
      let partial = if denom > 0.0 { sum_xy / denom } else { 0.0 }
      (nugget, partial)
    } else {
      let denom = count * sum_xx - sum_x * sum_x
      if denom.abs() < 1.0e-300 {
        (0.0, 0.0)
      } else {
        let b = (count * sum_xy - sum_x * sum_y) / denom
        let a = (sum_y - b * sum_x) / count
        (a, b)
      }
    }
    let safe_nugget = if fit_nugget < 0.0 { 0.0 } else { fit_nugget }
    let safe_partial = if fit_partial < 0.0 { 0.0 } else { fit_partial }
    let mut sse = 0.0
    for point in empirical {
      let shape = voy_variogram_shape(model, point.lag, candidate_range)
      let predicted = safe_nugget + safe_partial * shape
      let residual = point.gamma - predicted
      sse = sse + residual * residual
    }
    if sse < best_sse {
      best_sse = sse
      best_nugget = safe_nugget
      best_sill = safe_partial
      best_range = candidate_range
    }
    fi = fi + 1
  }
  VoyagerVariogramModel::{
    model_type: model,
    nugget: best_nugget,
    sill: best_sill,
    range: best_range,
    fitted_sse: best_sse,
    n_points: empirical.length(),
  }
}

///|
/// Evaluate a fitted variogram model at distance `h`.
pub fn voyager_variogram_predict(
  model : VoyagerVariogramModel,
  h : Double,
) -> Double {
  model.nugget +
  model.sill * voy_variogram_shape(model.model_type, h, model.range)
}

// ===========================================================================
// Spatial Moran correlogram
// ===========================================================================

///|
pub struct VoyagerCorrelogramPoint {
  lag : Double
  bin_min : Double
  bin_max : Double
  morans_i : Double
  expectation : Double
  npairs : Int
} derive(Debug)
pub extend VoyagerCorrelogramPoint with @moonbitlang/core/debug.Debug::{to_repr}

///|
/// Moran correlogram. Builds binary distance-band weights at each of `n_lags`
/// distance bins and computes Moran's I with its analytical expectation.
pub fn voyager_correlogram(
  coords : Array[Array[Double]],
  values : Array[Double],
  n_lags? : Int = 10,
  max_lag? : Double = 0.0,
  style? : String = "W",
) -> Array[VoyagerCorrelogramPoint] raise VoyagerError {
  if coords.length() < 3 {
    voy_fail("Voyager correlogram requires at least 3 spots")
  }
  if values.length() != coords.length() {
    voy_fail("Voyager correlogram values must match coordinates length")
  }
  if n_lags < 1 {
    voy_fail("Voyager correlogram n_lags must be positive")
  }
  for value in values {
    if !voy_is_finite(value) {
      voy_fail("Voyager correlogram values must be finite")
    }
  }
  let n = coords.length()
  let distances = voy_pairwise_distances(coords)
  let mut max_dist = 0.0
  for i in 0..<n {
    for j in (i + 1)..<n {
      if distances[i][j] > max_dist {
        max_dist = distances[i][j]
      }
    }
  }
  if max_dist <= 0.0 {
    voy_fail("Voyager correlogram requires non-coincident coordinates")
  }
  let upper = if max_lag > 0.0 && voy_is_finite(max_lag) {
    max_lag
  } else {
    max_dist / 2.0
  }
  let lag_width = upper / n_lags.to_double()
  let output : Array[VoyagerCorrelogramPoint] = []
  let mut mean = 0.0
  for value in values {
    mean = mean + value
  }
  mean = mean / n.to_double()
  let centered : Array[Double] = []
  let mut ssd = 0.0
  for value in values {
    let dev = value - mean
    centered.push(dev)
    ssd = ssd + dev * dev
  }
  let expectation = -1.0 / (n.to_double() - 1.0)
  for bin in 0..<n_lags {
    let bin_min = bin.to_double() * lag_width
    let bin_max = (bin + 1).to_double() * lag_width
    let neighbors : Array[Array[Int]] = []
    let raw : Array[Array[Double]] = []
    let mut npairs = 0
    for i in 0..<n {
      let nb : Array[Int] = []
      let w : Array[Double] = []
      for j in 0..<n {
        if i != j {
          let d = distances[i][j]
          if d > bin_min && d <= bin_max {
            nb.push(j)
            w.push(1.0)
            npairs = npairs + 1
          }
        }
      }
      neighbors.push(nb)
      raw.push(w)
    }
    // Skip distance bins with no observed pairs — they carry no information.
    if npairs == 0 {
      continue
    }
    let styled = voy_apply_style(neighbors, raw, style)
    let weights = VoyagerWeights::{ style, n, neighbors, weights: styled, }
    if ssd <= 0.0 {
      output.push(VoyagerCorrelogramPoint::{
        lag: (bin_min + bin_max) / 2.0,
        bin_min,
        bin_max,
        morans_i: 0.0,
        expectation,
        npairs: npairs / 2,
      })
      continue
    }
    let s0 = voyager_weights_s0(weights)
    if s0 <= 0.0 {
      output.push(VoyagerCorrelogramPoint::{
        lag: (bin_min + bin_max) / 2.0,
        bin_min,
        bin_max,
        morans_i: 0.0,
        expectation,
        npairs: npairs / 2,
      })
      continue
    }
    let mut numerator = 0.0
    for i in 0..<n {
      let nb = weights.neighbors[i]
      let wv = weights.weights[i]
      let mut lag = 0.0
      for k in 0..<nb.length() {
        lag = lag + wv[k] * centered[nb[k]]
      }
      numerator = numerator + centered[i] * lag
    }
    let morans_i = n.to_double() / s0 * (numerator / ssd)
    output.push(VoyagerCorrelogramPoint::{
      lag: (bin_min + bin_max) / 2.0,
      bin_min,
      bin_max,
      morans_i,
      expectation,
      npairs: npairs / 2,
    })
  }
  output
}

// ===========================================================================
// SpatialExperiment integration (immutable write-back)
// ===========================================================================

///|
pub struct VoyagerUnivariateSFEOutput {
  experiment : SpatialExperiment
  weights : VoyagerWeights
  results : Array[VoyagerUnivariateResult]
}

///|
pub struct VoyagerUnivariateResult {
  feature : String
  stat_method : String
  global_estimate : Double
  global_p_value : Double
  local_estimate : Array[Double]
  local_fdr : Array[Double]
  local_quadrants : Array[String]
}

///|
fn voy_copy_spatial_experiment(
  experiment : SpatialExperiment,
) -> SpatialExperiment {
  let output = SpatialExperiment::new()
  for assay_name in experiment.assay.keys() {
    output.assay[assay_name] = voy_copy_double_matrix(
      experiment.assay[assay_name],
    )
  }
  for row in experiment.row_data {
    output.row_data.push(voy_copy_string_map(row))
  }
  for column in experiment.col_data {
    output.col_data.push(voy_copy_string_map(column))
  }
  for coordinate in experiment.spatial_coords {
    output.spatial_coords.push(coordinate)
  }
  for image in experiment.images {
    output.images.push(image)
  }
  for key in experiment.metadata.keys() {
    output.metadata[key] = experiment.metadata[key]
  }
  output
}

///|
fn voy_experiment_gene_names(
  experiment : SpatialExperiment,
  genes : Int,
) -> Array[String] {
  let names : Array[String] = []
  for gene in 0..<genes {
    let generated = "gene_" + (gene + 1).to_string()
    if experiment.row_data.length() == 0 || gene >= experiment.row_data.length() {
      names.push(generated)
    } else {
      let row = experiment.row_data[gene]
      let name = match row.get("gene_name") {
        Some(value) =>
          if value.length() > 0 {
            value
          } else {
            match row.get("gene_id") {
              Some(identifier) =>
                if identifier.length() > 0 {
                  identifier
                } else {
                  generated
                }
              None => generated
            }
          }
        None =>
          match row.get("gene_id") {
            Some(identifier) =>
              if identifier.length() > 0 {
                identifier
              } else {
                generated
              }
            None => generated
          }
      }
      names.push(name)
    }
  }
  names
}

///|
fn voy_experiment_coordinates(
  experiment : SpatialExperiment,
  spots : Int,
) -> Array[Array[Double]] raise VoyagerError {
  if experiment.spatial_coords.length() != spots {
    voy_fail(
      "Voyager SpatialExperiment requires one spatial coordinate per spot",
    )
  }
  let mut use_z = false
  if spots > 0 {
    let first = experiment.spatial_coords[0].z
    for coordinate in experiment.spatial_coords {
      if (coordinate.z - first).abs() > 1.0e-14 {
        use_z = true
      }
    }
  }
  let coordinates : Array[Array[Double]] = []
  for coordinate in experiment.spatial_coords {
    if use_z {
      coordinates.push([coordinate.x, coordinate.y, coordinate.z])
    } else {
      coordinates.push([coordinate.x, coordinate.y])
    }
  }
  coordinates
}

///|
/// Run univariate spatial autocorrelation over one or more features of a
/// SpatialExperiment and return a *new* SpatialExperiment with the local_res
/// results written into `col_data` and global results into `row_data`/
/// `metadata`. The input experiment is left unchanged.
pub fn voyager_run_univariate_sfe(
  experiment : SpatialExperiment,
  feature_indices : Array[Int],
  assay_name? : String = "logcounts",
  weights? : VoyagerWeights? = None,
  k? : Int = 4,
  bandwidth? : Double = 0.0,
  stat_method? : String = "moran",
  permutations? : Int = 0,
  seed? : Int = 20240501,
  fdr_threshold? : Double = 0.05,
  output_prefix? : String = "voyager",
) -> VoyagerUnivariateSFEOutput raise VoyagerError {
  if experiment.col_data.length() == 0 {
    voy_fail("Voyager SpatialExperiment col_data must not be empty")
  }
  let expression : Array[Array[Double]] = match
    experiment.assay.get(assay_name) {
    Some(value) => value
    None =>
      raise VoyagerError::VoyagerError(
        "Voyager SpatialExperiment assay not found: " + assay_name,
      )
  }
  if expression.length() == 0 {
    voy_fail("Voyager SpatialExperiment assay must not be empty")
  }
  let spots = expression[0].length()
  if spots != experiment.col_data.length() {
    voy_fail("Voyager SpatialExperiment assay columns must match col_data rows")
  }
  if spots != experiment.spatial_coords.length() {
    voy_fail("Voyager SpatialExperiment coordinates must match col_data rows")
  }
  if experiment.row_data.length() != 0 &&
    experiment.row_data.length() != expression.length() {
    voy_fail(
      "Voyager SpatialExperiment rowData must match assay rows or be empty",
    )
  }
  let write_row_data = experiment.row_data.length() == expression.length()
  if feature_indices.length() == 0 {
    voy_fail("Voyager SpatialExperiment feature selection must not be empty")
  }
  for index in feature_indices {
    if index < 0 || index >= expression.length() {
      voy_fail(
        "Voyager SpatialExperiment feature index out of range: " +
        index.to_string(),
      )
    }
  }
  let coords = voy_experiment_coordinates(experiment, spots)
  let weights_resolved = match weights {
    Some(value) => {
      voy_validate_weights(value, "SpatialExperiment")
      if value.n != spots {
        voy_fail("Voyager SpatialExperiment weights must match spot count")
      }
      value
    }
    None =>
      if bandwidth > 0.0 && voy_is_finite(bandwidth) {
        voyager_weights_distance_band(coords, bandwidth)
      } else {
        voyager_weights_knn(coords, k)
      }
  }
  let copy = voy_copy_spatial_experiment(experiment)
  let gene_names = voy_experiment_gene_names(experiment, expression.length())
  let results : Array[VoyagerUnivariateResult] = []
  let normalized_method : String = match stat_method {
    "moran" => "moran"
    "geary" => "geary"
    "getis" => "getis"
    "getis_ord" => "getis"
    "gi" => "getis"
    _ =>
      raise VoyagerError::VoyagerError(
        "Voyager SpatialExperiment stat_method must be moran/geary/getis",
      )
  }
  for feature_index in feature_indices {
    let vector : Array[Double] = []
    for spot in 0..<spots {
      vector.push(expression[feature_index][spot])
    }
    let feature_name = gene_names[feature_index]
    let suffix = "." + feature_name
    match normalized_method {
      "moran" => {
        let global = voyager_global_morans_i(
          vector,
          weights_resolved,
          feature=feature_name,
        )
        let local_res = voyager_local_morans_i(
          vector,
          weights_resolved,
          permutations~,
          seed~,
          fdr_threshold~,
          feature=feature_name,
        )
        if write_row_data {
          copy.row_data[feature_index][output_prefix + ".moran.I" + suffix] = global.estimate.to_string()
          copy.row_data[feature_index][output_prefix + ".moran.p" + suffix] = global.p_value.to_string()
        }
        let local_str : Array[String] = []
        let fdr_str : Array[String] = []
        let quad_str : Array[String] = []
        for spot in 0..<spots {
          local_str.push(local_res.local_i[spot].to_string())
          fdr_str.push(local_res.fdr[spot].to_string())
          quad_str.push(local_res.quadrants[spot])
        }
        for spot in 0..<spots {
          copy.col_data[spot][output_prefix + ".moran.local" + suffix] = local_str[spot]
          copy.col_data[spot][output_prefix + ".moran.fdr" + suffix] = fdr_str[spot]
          copy.col_data[spot][output_prefix + ".moran.quadrant" + suffix] = quad_str[spot]
        }
        results.push(VoyagerUnivariateResult::{
          feature: feature_name,
          stat_method: "moran",
          global_estimate: global.estimate,
          global_p_value: global.p_value,
          local_estimate: local_res.local_i,
          local_fdr: local_res.fdr,
          local_quadrants: local_res.quadrants,
        })
      }
      "geary" => {
        let global = voyager_global_gearys_c(
          vector,
          weights_resolved,
          feature=feature_name,
        )
        let local_res = voyager_local_gearys_c(
          vector,
          weights_resolved,
          permutations~,
          seed~,
          fdr_threshold~,
          feature=feature_name,
        )
        if write_row_data {
          copy.row_data[feature_index][output_prefix + ".geary.C" + suffix] = global.estimate.to_string()
          copy.row_data[feature_index][output_prefix + ".geary.p" + suffix] = global.p_value.to_string()
        }
        let local_str : Array[String] = []
        let fdr_str : Array[String] = []
        let class_str : Array[String] = []
        for spot in 0..<spots {
          local_str.push(local_res.local_c[spot].to_string())
          fdr_str.push(local_res.fdr[spot].to_string())
          class_str.push(local_res.classifications[spot])
        }
        for spot in 0..<spots {
          copy.col_data[spot][output_prefix + ".geary.local" + suffix] = local_str[spot]
          copy.col_data[spot][output_prefix + ".geary.fdr" + suffix] = fdr_str[spot]
          copy.col_data[spot][output_prefix + ".geary.class" + suffix] = class_str[spot]
        }
        results.push(VoyagerUnivariateResult::{
          feature: feature_name,
          stat_method: "geary",
          global_estimate: global.estimate,
          global_p_value: global.p_value,
          local_estimate: local_res.local_c,
          local_fdr: local_res.fdr,
          local_quadrants: local_res.classifications,
        })
      }
      _ => {
        let local_res = voyager_local_getis_ord(
          vector,
          weights_resolved,
          star=true,
          fdr_threshold~,
          feature=feature_name,
        )
        let local_str : Array[String] = []
        let fdr_str : Array[String] = []
        let class_str : Array[String] = []
        let z_str : Array[String] = []
        for spot in 0..<spots {
          local_str.push(local_res.statistic[spot].to_string())
          fdr_str.push(local_res.fdr[spot].to_string())
          class_str.push(local_res.classifications[spot])
          z_str.push(local_res.z_scores[spot].to_string())
        }
        for spot in 0..<spots {
          copy.col_data[spot][output_prefix + ".getis.gi_star" + suffix] = local_str[spot]
          copy.col_data[spot][output_prefix + ".getis.fdr" + suffix] = fdr_str[spot]
          copy.col_data[spot][output_prefix + ".getis.class" + suffix] = class_str[spot]
          copy.col_data[spot][output_prefix + ".getis.z" + suffix] = z_str[spot]
        }
        results.push(VoyagerUnivariateResult::{
          feature: feature_name,
          stat_method: "getis",
          global_estimate: 0.0,
          global_p_value: 1.0,
          local_estimate: local_res.statistic,
          local_fdr: local_res.fdr,
          local_quadrants: local_res.classifications,
        })
      }
    }
  }
  copy.metadata[output_prefix + ".method"] = normalized_method
  copy.metadata[output_prefix + ".n_features"] = feature_indices
    .length()
    .to_string()
  copy.metadata[output_prefix + ".weights_style"] = weights_resolved.style
  VoyagerUnivariateSFEOutput::{
    experiment: copy,
    weights: weights_resolved,
    results,
  }
}

///|
/// Example SpatialExperiment with a smooth spatial gradient plus a hotspot,
/// suitable for exercising the univariate and bivariate workflows offline.
pub fn voyager_example_spatial_experiment() -> SpatialExperiment {
  let se = SpatialExperiment::new()
  let _ = se_set_metadata(se, "platform", "Visium")
  let _ = se_set_metadata(se, "sample_id", "voyager_example")
  // 6x6 grid of spots.
  let spots = 36
  let genes = 4
  let gene_names = ["gene_gradient", "gene_hotspot", "gene_noise", "gene_anti"]
  for gene in 0..<genes {
    let row : Map[String, String] = Map([], capacity=4)
    row["gene_id"] = "ENSG" + (gene + 1).to_string()
    row["gene_name"] = gene_names[gene]
    se.row_data.push(row)
  }
  // Build assay[feature][spot].
  let assay : Array[Array[Double]] = []
  for _ in 0..<genes {
    assay.push(Array::make(spots, 0.0))
  }
  let mut spot = 0
  for row_index in 0..<6 {
    for col_index in 0..<6 {
      let x = col_index.to_double()
      let y = row_index.to_double()
      se.col_data.push(se_make_col("spot_" + spot.to_string(), "tissue"))
      se.spatial_coords.push(SpatialCoord::new_2d(x, y))
      // Smooth gradient increasing with x.
      assay[0][spot] = x + y * 0.5
      // Hotspot centred near (3, 3).
      let dx = x - 3.0
      let dy = y - 3.0
      let dist_sq = dx * dx + dy * dy
      assay[1][spot] = 20.0 * @math.exp(-dist_sq / 4.0)
      // Pseudo-noise (deterministic so tests are stable).
      assay[2][spot] = ((spot * 37 + 11) % 13).to_double()
      // Anti-gradient (negatively correlated with feature 0).
      assay[3][spot] = 10.0 - (x + y * 0.5)
      spot = spot + 1
    }
  }
  let _ = se_add_assay(se, "logcounts", assay)
  se
}

///|
fn se_make_col(spot_id : String, tissue : String) -> Map[String, String] {
  let column : Map[String, String] = Map([], capacity=2)
  column["spot_id"] = spot_id
  column["tissue"] = tissue
  column
}
