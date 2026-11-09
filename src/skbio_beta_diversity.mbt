/// skbio_beta_diversity.mbt — Non-phylogenetic beta diversity metrics
/// Item 5: bray_curtis, jaccard, canberra, euclidean, manhattan, chebyshev, gower

///|
/// Bray-Curtis dissimilarity between two abundance vectors.
pub fn skbio_bray_curtis(a : Array[Double], b : Array[Double]) -> Double {
  let n = if a.length() < b.length() { a.length() } else { b.length() }
  let mut numerator = 0.0
  let mut sum_a = 0.0
  let mut sum_b = 0.0
  for i = 0; i < n; i = i + 1 {
    let diff = a[i] - b[i]
    numerator = numerator + (if diff < 0.0 { -diff } else { diff })
    sum_a = sum_a + a[i]
    sum_b = sum_b + b[i]
  }
  if sum_a + sum_b == 0.0 {
    return 0.0
  }
  numerator / (sum_a + sum_b)
}

///|
/// Jaccard dissimilarity between two presence/absence vectors.
pub fn skbio_jaccard(a : Array[Double], b : Array[Double]) -> Double {
  let n = if a.length() < b.length() { a.length() } else { b.length() }
  let mut intersection = 0
  let mut union = 0
  for i = 0; i < n; i = i + 1 {
    let pa = a[i] > 0.0
    let pb = b[i] > 0.0
    if pa && pb {
      intersection = intersection + 1
      union = union + 1
    } else if pa || pb {
      union = union + 1
    }
  }
  if union == 0 {
    return 0.0
  }
  1.0 - intersection.to_double() / union.to_double()
}

///|
/// Canberra distance between two vectors.
pub fn skbio_canberra(a : Array[Double], b : Array[Double]) -> Double {
  let n = if a.length() < b.length() { a.length() } else { b.length() }
  let mut total = 0.0
  for i = 0; i < n; i = i + 1 {
    let diff = a[i] - b[i]
    let abs_diff = if diff < 0.0 { -diff } else { diff }
    let abs_a = if a[i] < 0.0 { -a[i] } else { a[i] }
    let abs_b = if b[i] < 0.0 { -b[i] } else { b[i] }
    let denom = abs_a + abs_b
    if denom > 0.0 {
      total = total + abs_diff / denom
    }
  }
  total
}

///|
/// Euclidean distance between two vectors.
pub fn skbio_euclidean(a : Array[Double], b : Array[Double]) -> Double {
  let n = if a.length() < b.length() { a.length() } else { b.length() }
  let mut sum_sq = 0.0
  for i = 0; i < n; i = i + 1 {
    let diff = a[i] - b[i]
    sum_sq = sum_sq + diff * diff
  }
  @math.pow(sum_sq, 0.5)
}

///|
/// Manhattan (L1) distance between two vectors.
pub fn skbio_manhattan(a : Array[Double], b : Array[Double]) -> Double {
  let n = if a.length() < b.length() { a.length() } else { b.length() }
  let mut total = 0.0
  for i = 0; i < n; i = i + 1 {
    let diff = a[i] - b[i]
    total = total + (if diff < 0.0 { -diff } else { diff })
  }
  total
}

///|
/// Chebyshev (L-inf) distance: max absolute difference.
pub fn skbio_chebyshev(a : Array[Double], b : Array[Double]) -> Double {
  let n = if a.length() < b.length() { a.length() } else { b.length() }
  let mut max_diff = 0.0
  for i = 0; i < n; i = i + 1 {
    let diff = a[i] - b[i]
    let abs_diff = if diff < 0.0 { -diff } else { diff }
    if abs_diff > max_diff {
      max_diff = abs_diff
    }
  }
  max_diff
}

///|
/// Gower distance between two vectors (generalized, each feature range-normalized).
/// ranges: max-min for each feature (0 means feature is constant, skip).
pub fn skbio_gower(
  a : Array[Double],
  b : Array[Double],
  ranges : Array[Double],
) -> Double {
  let n = if a.length() < b.length() { a.length() } else { b.length() }
  let mut total = 0.0
  let mut count = 0.0
  for i = 0; i < n; i = i + 1 {
    if ranges[i] > 0.0 {
      let diff = a[i] - b[i]
      let abs_diff = if diff < 0.0 { -diff } else { diff }
      total = total + abs_diff / ranges[i]
      count = count + 1.0
    }
  }
  if count == 0.0 {
    return 0.0
  }
  total / count
}

///|
/// Beta diversity metric type.
pub enum SkbioBetaMetric {
  SkbioBrayCurtis
  SkbioJaccard
  SkbioCanberra
  SkbioEuclidean
  SkbioManhattan
  SkbioChebyshev
}

///|
/// Compute full beta diversity distance matrix between all pairs of samples.
/// samples: each row is a sample, each column is a feature/OTU.
pub fn skbio_beta_diversity(
  samples : Array[Array[Double]],
  metric~ : SkbioBetaMetric = SkbioBrayCurtis,
) -> Array[Array[Double]] {
  let n = samples.length()
  let matrix : Array[Array[Double]] = Array(capacity=n)
  for i = 0; i < n; i = i + 1 {
    let row : Array[Double] = Array(capacity=n)
    for j = 0; j < n; j = j + 1 {
      if i == j {
        row.push(0.0)
      } else {
        let d = match metric {
          SkbioBrayCurtis => skbio_bray_curtis(samples[i], samples[j])
          SkbioJaccard => skbio_jaccard(samples[i], samples[j])
          SkbioCanberra => skbio_canberra(samples[i], samples[j])
          SkbioEuclidean => skbio_euclidean(samples[i], samples[j])
          SkbioManhattan => skbio_manhattan(samples[i], samples[j])
          SkbioChebyshev => skbio_chebyshev(samples[i], samples[j])
        }
        row.push(d)
      }
    }
    matrix.push(row)
  }
  matrix
}

///|
pub fn skbio_beta_dummy() -> Unit {
  ignore(SkbioJaccard)
  ignore(SkbioCanberra)
  ignore(SkbioEuclidean)
  ignore(SkbioManhattan)
  ignore(SkbioChebyshev)
}
