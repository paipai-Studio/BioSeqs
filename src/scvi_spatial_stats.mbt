// scvi_spatial_stats.mbt - Spatial autocorrelation statistics
//
// Complements scVI spatial models and the existing spatial_nhood / spatialdecon
// modules with Moran's I, Geary's C, permutation testing and neighbour label
// enrichment (squidpy-style).

///|
fn scvi_spatial_distance(a : Array[Double], b : Array[Double]) -> Double {
  let mut sq = 0.0
  let n = if a.length() < b.length() { a.length() } else { b.length() }
  for i in 0..<n {
    let d = a[i] - b[i]
    sq += d * d
  }
  @math.pow(sq, 0.5)
}

///|
/// k-nearest-neighbour spatial weight matrix with inverse-distance weights,
/// row-normalised and excluding self-loops.
pub fn scvi_spatial_knn_weights(
  coords : Array[Array[Double]],
  k : Int32,
) -> Array[Array[Double]] {
  let n = coords.length()
  let out = Array(capacity=n)
  for _i in 0..<n {
    let row = Array::make(n, 0.0)
    out.push(row)
  }
  if n == 0 {
    return out
  }
  let kk = if k < 1 { 1 } else if k > n - 1 { n - 1 } else { k }
  for i in 0..<n {
    // Select k nearest neighbours by distance.
    let dist = Array(capacity=n)
    for j in 0..<n {
      if i == j {
        dist.push(1.0e30)
      } else {
        dist.push(scvi_spatial_distance(coords[i], coords[j]))
      }
    }
    for _pick in 0..<kk {
      let mut best = -1
      let mut best_d = 1.0e30
      for j in 0..<n {
        if dist[j] < best_d {
          best_d = dist[j]
          best = j
        }
      }
      if best >= 0 {
        let w = if best_d <= 0.0 { 1.0 } else { 1.0 / best_d }
        out[i][best] = w
        dist[best] = 1.0e30
      }
    }
    let mut s = 0.0
    for j in 0..<n {
      s += out[i][j]
    }
    if s > 0.0 {
      for j in 0..<n {
        out[i][j] = out[i][j] / s
      }
    }
  }
  out
}

///|
fn scvi_spatial_row_sums(weights : Array[Array[Double]]) -> Double {
  let mut s = 0.0
  for i in 0..<weights.length() {
    for j in 0..<weights[i].length() {
      s += weights[i][j]
    }
  }
  s
}

///|
fn scvi_spatial_mean(values : Array[Double]) -> Double {
  if values.length() == 0 {
    return 0.0
  }
  let mut s = 0.0
  for i in 0..<values.length() {
    s += values[i]
  }
  s / values.length().to_double()
}

///|
/// Moran's I spatial autocorrelation.
pub fn scvi_morans_i(
  values : Array[Double],
  weights : Array[Array[Double]],
) -> Double {
  let n = values.length()
  if n < 2 {
    return 0.0
  }
  let mean = scvi_spatial_mean(values)
  let mut num = 0.0
  let mut den = 0.0
  for i in 0..<n {
    let di = values[i] - mean
    den += di * di
    let wi = weights[i]
    for j in 0..<n {
      if j < wi.length() {
        num += wi[j] * di * (values[j] - mean)
      }
    }
  }
  let s0 = scvi_spatial_row_sums(weights)
  if den <= 0.0 || s0 <= 0.0 {
    return 0.0
  }
  n.to_double() / s0 * (num / den)
}

///|
/// Geary's C spatial autocorrelation.
pub fn scvi_gearys_c(
  values : Array[Double],
  weights : Array[Array[Double]],
) -> Double {
  let n = values.length()
  if n < 2 {
    return 0.0
  }
  let mean = scvi_spatial_mean(values)
  let mut num = 0.0
  let mut den = 0.0
  for i in 0..<n {
    let di = values[i] - mean
    den += di * di
    let wi = weights[i]
    for j in 0..<n {
      if j < wi.length() {
        let d = values[i] - values[j]
        num += wi[j] * d * d
      }
    }
  }
  let s0 = scvi_spatial_row_sums(weights)
  if den <= 0.0 || s0 <= 0.0 {
    return 0.0
  }
  (n - 1).to_double() / (2.0 * s0) * (num / den)
}

///|
/// Two-sided permutation test for Moran's I.
pub fn scvi_spatial_permutation_test(
  values : Array[Double],
  coords : Array[Array[Double]],
  k : Int32,
  n_perm : Int32,
  seed : Int32,
) -> Double {
  let weights = scvi_spatial_knn_weights(coords, k)
  let observed = scvi_morans_i(values, weights)
  let n = values.length()
  if n == 0 {
    return 1.0
  }
  let work = Array(capacity=n)
  for i in 0..<n {
    work.push(values[i])
  }
  let mut extreme = 0
  let perms = if n_perm < 1 { 1 } else { n_perm }
  for p in 0..<perms {
    // Fisher-Yates permutation driven by the deterministic RNG.
    for i in 0..<n {
      let j = (scvi_dist_uniform(seed + p * 104729 + i * 131) *
      (n - i).to_double()).to_int()
      let jj = if j < 0 { 0 } else if j >= n - i { n - i - 1 } else { j }
      let a = work[i]
      work[i] = work[i + jj]
      work[i + jj] = a
    }
    let perm_i = scvi_morans_i(work, weights)
    if perm_i.abs() >= observed.abs() {
      extreme = extreme + 1
    }
  }
  (extreme.to_double() + 1.0) / (perms.to_double() + 1.0)
}

///|
/// Moran's I and Geary's C in one call.
pub fn scvi_spatial_autocorrelation(
  values : Array[Double],
  coords : Array[Array[Double]],
  k : Int32,
) -> (Double, Double) {
  let weights = scvi_spatial_knn_weights(coords, k)
  (scvi_morans_i(values, weights), scvi_gearys_c(values, weights))
}

///|
/// Fraction of spatial neighbours that share the focal cell's label. Values
/// above the global label fraction indicate spatial enrichment.
pub fn scvi_spatial_neighbor_enrichment(
  labels : Array[Int32],
  coords : Array[Array[Double]],
  k : Int32,
) -> Double {
  let n = labels.length()
  if n == 0 {
    return 0.0
  }
  let weights = scvi_spatial_knn_weights(coords, k)
  let mut same = 0.0
  let mut total = 0.0
  for i in 0..<n {
    for j in 0..<n {
      if j < weights[i].length() && weights[i][j] > 0.0 && i != j {
        total += 1.0
        if labels[i] == labels[j] {
          same += 1.0
        }
      }
    }
  }
  if total <= 0.0 {
    return 0.0
  }
  same / total
}

///|
/// Global label fraction used as the null expectation for enrichment.
pub fn scvi_spatial_label_fraction(
  labels : Array[Int32],
  label : Int32,
) -> Double {
  if labels.length() == 0 {
    return 0.0
  }
  let mut count = 0
  for i in 0..<labels.length() {
    if labels[i] == label {
      count = count + 1
    }
  }
  count.to_double() / labels.length().to_double()
}
