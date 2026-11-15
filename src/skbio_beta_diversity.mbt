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

// ============================================================================
// Whittaker β 多样性（群落层面标量）
// ============================================================================

///|
/// Whittaker's β diversity (原始形式): β_W = γ / ᾱ − 1。
///
/// 其中：
/// - αᵢ = 样本 i 的观察丰富度（count of 列>0 的 feature）
/// - γ  = 跨所有样本的总观察丰富度（count of 至少一个样本>0 的列）
/// - ᾱ = mean(αᵢ)
///
/// 输入 `samples` 是 OTU table（样本 × feature）。返回单个 Double 值。
/// 当 ᾱ = 0（全部样本全零）时返回 0。
pub fn skbio_whittaker_beta(samples : Array[Array[Double]]) -> Double {
  let n = samples.length()
  if n == 0 { return 0.0 }
  let d = if n > 0 { samples[0].length() } else { 0 }
  let mut gamma_count = 0
  let mut alpha_sum = 0
  for f = 0; f < d; f = f + 1 {
    let mut any_present = false
    for i = 0; i < n; i = i + 1 {
      if samples[i][f] > 0.0 {
        any_present = true
        break
      }
    }
    if any_present {
      gamma_count = gamma_count + 1
    }
  }
  for i = 0; i < n; i = i + 1 {
    let mut present = 0
    for f = 0; f < d; f = f + 1 {
      if samples[i][f] > 0.0 { present = present + 1 }
    }
    alpha_sum = alpha_sum + present
  }
  let mean_alpha = alpha_sum.to_double() / n.to_double()
  if mean_alpha <= 0.0 {
    0.0
  } else {
    gamma_count.to_double() / mean_alpha - 1.0
  }
}

// ============================================================================
// subsample_counts（有放回/无放回稀释重采样）
// ============================================================================

///|
/// 有放回 multinomial 重采样：把 `counts` 视为 multinomial 概率，抽 `depth`
/// 次放回，返回新的 count 数组（长度与 counts 相同）。
pub fn skbio_subsample_counts_with_replacement(
  counts : Array[Double],
  depth : Int32,
  seed : Int32,
) -> Array[Double] {
  let d = counts.length()
  let mut total = 0.0
  for c in counts { total = total + c }
  if total <= 0.0 || depth <= 0 {
    return Array::make(d, 0.0)
  }
  // 归一化累积概率
  let cumul : Array[Double] = Array(capacity=d)
  let mut cums = 0.0
  for c in counts {
    cums = cums + c / total
    cumul.push(cums)
  }
  let out = Array::make(d, 0.0)
  let mut rng = seed
  let mut i = 0
  while i < depth {
    rng = (rng * 1103515245 + 12345) & 0x7FFFFFFF
    let u = rng.to_double() / 2147483647.0
    // binary search
    let mut lo = 0
    let mut hi = d - 1
    while lo < hi {
      let mid = (lo + hi) / 2
      if cumul[mid] >= u {
        hi = mid
      } else {
        lo = mid + 1
      }
    }
    out[lo] = out[lo] + 1.0
    i = i + 1
  }
  out
}

///|
/// 无放回重采样（稀疏曲线 / rarefaction 标准做法）。
///
/// 把 counts 展开为一个扁平的 feature index 数组，做 Fisher–Yates shuffle，
/// 取前 `depth` 个样本，然后把它们重新计数回到 feature 维度。返回新的
/// count 数组。当 depth > total 时返回全零。
pub fn skbio_subsample_counts(
  counts : Array[Double],
  depth : Int32,
  seed : Int32,
) -> Array[Double] {
  let d = counts.length()
  let d_i = (depth.to_double()).floor().to_int()
  let mut total : Int = 0
  for c in counts { total = total + c.floor().to_int() }
  if total <= 0 || d_i <= 0 || d_i > total {
    return Array::make(d, 0.0)
  }
  let flat : Array[Int] = Array(capacity=total)
  for f = 0; f < d; f = f + 1 {
    let nf = counts[f].floor().to_int()
    let mut k = 0
    while k < nf {
      flat.push(f)
      k = k + 1
    }
  }
  let mut rng = seed
  let mut i = total - 1
  while i >= (total - d_i) {
    rng = (rng * 1103515245 + 12345) & 0x7FFFFFFF
    let j = rng % (i + 1)
    let t = flat[i]; flat[i] = flat[j]; flat[j] = t
    i = i - 1
  }
  let out = Array::make(d, 0.0)
  let mut k = 0
  while k < d_i {
    let f = flat[k]
    out[f] = out[f] + 1.0
    k = k + 1
  }
  out
}

///|
/// isubsample — 对多个样本表（OTU table，每行一个样本）分别做无放回重采样。
/// depth 可以是单值（所有样本同深度）或 Array（每个样本一个深度）。
pub fn skbio_isubsample(
  samples : Array[Array[Double]],
  depths : Array[Int32],
  seed : Int32,
) -> Array[Array[Double]] {
  let n = samples.length()
  let lc = depths.length()
  let out : Array[Array[Double]] = Array(capacity=n)
  for i = 0; i < n; i = i + 1 {
    let dp = if lc == 1 { depths[0] } else {
      if i < lc { depths[i] } else { 0 }
    }
    out.push(skbio_subsample_counts(samples[i], dp, seed + i))
  }
  out
}
