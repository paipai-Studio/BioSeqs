///|
/// scikit-bio alpha diversity metrics.
///
/// Port of the formula-based alpha diversity metrics from
/// `skbio.diversity.alpha` (scikit-bio 0.7.3). Each metric takes a vector of
/// taxon counts; zero counts are dropped first (absent taxa), matching
/// scikit-bio's `_validate_alpha` wrapper. Metrics that scikit-bio marks
/// `empty=np.nan` return NaN for an all-zero community; count metrics return 0.
///
/// Included: sobs, singles, doubles, osd, dominance, simpson_d, inv_simpson,
/// enspie, berger_parker_d, simpson_e, heip_e, goods_coverage, margalef,
/// menhinick, mcintosh_d/e, brillouin_d, strong, renyi, tsallis, hill,
/// **shannon**, **chao1**, **ace** (plus **ace_default** with threshold=10).
///
/// Not yet included: `finite=True` variants, `gini_index`,
/// `kempton_taylor_q`, `esty_ci`/`lladser_ci`/`lladser_pe` (trie-based),
/// `fisher_alpha`.

///|
fn skbio_alpha_nan() -> Double {
  0.0 / 0.0
}

///|
fn skbio_alpha_sum(a : Array[Double]) -> Double {
  let mut s = 0.0
  for v in a {
    s = s + v
  }
  s
}

///|
/// Drop zero counts (absent taxa), matching scikit-bio's wrapper.
fn skbio_alpha_nonzero(counts : Array[Double]) -> Array[Double] {
  let out : Array[Double] = []
  for c in counts {
    if c != 0.0 {
      out.push(c)
    }
  }
  out
}

///|
/// Shannon entropy of a proportion vector (scikit-bio's `_entropy`).
fn skbio_alpha_entropy(probs : Array[Double]) -> Double {
  let mut s = 0.0
  for p in probs {
    if p > 0.0 {
      s = s - p * @math.ln(p)
    }
  }
  s
}

///|
/// `ln(n!)` via the factorial recurrence (exact for integer counts).
fn skbio_alpha_log_factorial(n : Int32) -> Double {
  let mut s = 0.0
  let mut k = 1
  while k <= n {
    s = s + @math.ln(k.to_double())
    k = k + 1
  }
  s
}

///|
/// Sort a Double array in descending order (for Strong's index).
fn skbio_alpha_sort_desc(a : Array[Double]) -> Array[Double] {
  let n = a.length()
  let b = Array::make(n, 0.0)
  for i = 0; i < n; i = i + 1 {
    b[i] = a[i]
  }
  for i = 1; i < n; i = i + 1 {
    let key = b[i]
    let mut j = i - 1
    while j >= 0 && b[j] < key {
      b[j + 1] = b[j]
      j = j - 1
    }
    b[j + 1] = key
  }
  b
}

// ---------------------------------------------------------------------------
// Count / richness metrics
// ---------------------------------------------------------------------------

///|
/// Observed species richness (number of nonzero taxa), a.k.a. `sobs` /
/// `observed_otus`.
pub fn skbio_alpha_sobs(counts : Array[Double]) -> Int32 {
  skbio_alpha_nonzero(counts).length()
}

///|
/// Number of singletons (taxa observed exactly once).
pub fn skbio_alpha_singles(counts : Array[Double]) -> Int32 {
  let c = skbio_alpha_nonzero(counts)
  let mut n = 0
  for v in c {
    if v == 1.0 {
      n = n + 1
    }
  }
  n
}

///|
/// Number of doubletons (taxa observed exactly twice).
pub fn skbio_alpha_doubles(counts : Array[Double]) -> Int32 {
  let c = skbio_alpha_nonzero(counts)
  let mut n = 0
  for v in c {
    if v == 2.0 {
      n = n + 1
    }
  }
  n
}

///|
/// Observed taxa, singletons, and doubletons (`osd`).
pub fn skbio_alpha_osd(counts : Array[Double]) -> (Int32, Int32, Int32) {
  (
    skbio_alpha_sobs(counts),
    skbio_alpha_singles(counts),
    skbio_alpha_doubles(counts),
  )
}

// ---------------------------------------------------------------------------
// Dominance / evenness metrics
// ---------------------------------------------------------------------------

///|
/// Simpson's dominance index (concentration), `sum(p_i^2)`.
pub fn skbio_alpha_dominance(counts : Array[Double]) -> Double {
  let c = skbio_alpha_nonzero(counts)
  let n = skbio_alpha_sum(c)
  if c.length() == 0 {
    return skbio_alpha_nan()
  }
  let mut d = 0.0
  for v in c {
    let p = v / n
    d = d + p * p
  }
  d
}

///|
/// Alias of `skbio_alpha_dominance` (scikit-bio's `simpson_d`).
pub fn skbio_alpha_simpson_d(counts : Array[Double]) -> Double {
  skbio_alpha_dominance(counts)
}

///|
/// Inverse Simpson index (`1 / D`), a.k.a. `inv_simpson` / `enspie`.
pub fn skbio_alpha_inv_simpson(counts : Array[Double]) -> Double {
  1.0 / skbio_alpha_dominance(counts)
}

///|
/// Alias of `skbio_alpha_inv_simpson` (scikit-bio's `enspie`).
pub fn skbio_alpha_enspie(counts : Array[Double]) -> Double {
  skbio_alpha_inv_simpson(counts)
}

///|
/// Berger-Parker dominance index (`max(p_i)`).
pub fn skbio_alpha_berger_parker_d(counts : Array[Double]) -> Double {
  let c = skbio_alpha_nonzero(counts)
  let n = skbio_alpha_sum(c)
  if c.length() == 0 {
    return skbio_alpha_nan()
  }
  let mut mx = 0.0
  for v in c {
    if v > mx {
      mx = v
    }
  }
  mx / n
}

///|
/// Simpson's evenness index (`1 / (D * S)`).
pub fn skbio_alpha_simpson_e(counts : Array[Double]) -> Double {
  let c = skbio_alpha_nonzero(counts)
  if c.length() == 0 {
    return skbio_alpha_nan()
  }
  1.0 / (c.length().to_double() * skbio_alpha_dominance(c))
}

///|
/// Heip's evenness measure (`(e^H - 1) / (S - 1)`; 1.0 when S == 1).
pub fn skbio_alpha_heip_e(counts : Array[Double]) -> Double {
  let c = skbio_alpha_nonzero(counts)
  let s = c.length()
  if s == 0 {
    return skbio_alpha_nan()
  }
  if s == 1 {
    return 1.0
  }
  let n = skbio_alpha_sum(c)
  let probs = skbio_alpha_probs(c, n)
  let h = skbio_alpha_entropy(probs)
  (@math.exp(h) - 1.0) / (s - 1).to_double()
}

///|
/// Good's coverage estimator (`1 - singletons / N`).
pub fn skbio_alpha_goods_coverage(counts : Array[Double]) -> Double {
  let c = skbio_alpha_nonzero(counts)
  let n = skbio_alpha_sum(c)
  if c.length() == 0 {
    return skbio_alpha_nan()
  }
  1.0 - skbio_alpha_singles(c).to_double() / n
}

///|
/// Margalef's richness index (`(S - 1) / ln(N)`; NaN when N == 1).
pub fn skbio_alpha_margalef(counts : Array[Double]) -> Double {
  let c = skbio_alpha_nonzero(counts)
  let n = skbio_alpha_sum(c)
  if n == 1.0 || c.length() == 0 {
    return skbio_alpha_nan()
  }
  (c.length() - 1).to_double() / @math.ln(n)
}

///|
/// Menhinick's richness index (`S / sqrt(N)`).
pub fn skbio_alpha_menhinick(counts : Array[Double]) -> Double {
  let c = skbio_alpha_nonzero(counts)
  if c.length() == 0 {
    return skbio_alpha_nan()
  }
  c.length().to_double() / skbio_alpha_sum(c).sqrt()
}

///|
/// McIntosh dominance index (`(N - U) / (N - sqrt(N))`, `U = sqrt(sum n_i^2)`).
pub fn skbio_alpha_mcintosh_d(counts : Array[Double]) -> Double {
  let c = skbio_alpha_nonzero(counts)
  let n = skbio_alpha_sum(c)
  if n == 1.0 || c.length() == 0 {
    return skbio_alpha_nan()
  }
  let mut sq = 0.0
  for v in c {
    sq = sq + v * v
  }
  let u = sq.sqrt()
  (n - u) / (n - n.sqrt())
}

///|
/// McIntosh evenness measure (`U / sqrt((N - S + 1)^2 + S - 1)`).
pub fn skbio_alpha_mcintosh_e(counts : Array[Double]) -> Double {
  let c = skbio_alpha_nonzero(counts)
  let s = c.length()
  let n = skbio_alpha_sum(c)
  if s == 0 {
    return skbio_alpha_nan()
  }
  let mut sq = 0.0
  for v in c {
    sq = sq + v * v
  }
  let numerator = sq.sqrt()
  let denom = ((n - s.to_double() + 1.0) * (n - s.to_double() + 1.0) +
  s.to_double() -
  1.0).sqrt()
  numerator / denom
}

///|
/// Brillouin's diversity index (`(ln N! - sum ln n_i!) / N`).
pub fn skbio_alpha_brillouin_d(counts : Array[Double]) -> Double {
  let c = skbio_alpha_nonzero(counts)
  let n = skbio_alpha_sum(c)
  if c.length() == 0 {
    return skbio_alpha_nan()
  }
  let mut sum_lfact = 0.0
  for v in c {
    sum_lfact = sum_lfact + skbio_alpha_log_factorial(v.to_int())
  }
  (skbio_alpha_log_factorial(n.to_int()) - sum_lfact) / n
}

///|
/// Strong's dominance index (max over cumulative ranked abundance minus i/S).
pub fn skbio_alpha_strong(counts : Array[Double]) -> Double {
  let c = skbio_alpha_nonzero(counts)
  let s = c.length()
  if s == 0 {
    return skbio_alpha_nan()
  }
  let n = skbio_alpha_sum(c)
  let sorted = skbio_alpha_sort_desc(c)
  let mut cum = 0.0
  let mut maxv = -1.0e308
  for i = 0; i < s; i = i + 1 {
    cum = cum + sorted[i]
    let v = cum / n - (i + 1).to_double() / s.to_double()
    if v > maxv {
      maxv = v
    }
  }
  maxv
}

// ---------------------------------------------------------------------------
// Generalized entropy
// ---------------------------------------------------------------------------

///|
fn skbio_alpha_probs(c : Array[Double], n : Double) -> Array[Double] {
  let out : Array[Double] = []
  for v in c {
    out.push(v / n)
  }
  out
}

///|
fn skbio_alpha_renyi_impl(c : Array[Double], order : Double) -> Double {
  let s = c.length()
  if s == 0 {
    return skbio_alpha_nan()
  }
  if s == 1 {
    return 0.0
  }
  let n = skbio_alpha_sum(c)
  let probs = skbio_alpha_probs(c, n)
  if order == 1.0 {
    skbio_alpha_entropy(probs)
  } else if order == 1.0 / 0.0 {
    let mut mx = 0.0
    for p in probs {
      if p > mx {
        mx = p
      }
    }
    -@math.ln(mx)
  } else {
    let mut sum = 0.0
    for p in probs {
      sum = sum + @math.pow(p, order)
    }
    @math.ln(sum) / (1.0 - order)
  }
}

///|
/// Renyi entropy of order `order`.
pub fn skbio_alpha_renyi(counts : Array[Double], order : Double) -> Double {
  skbio_alpha_renyi_impl(skbio_alpha_nonzero(counts), order)
}

///|
/// Tsallis entropy of order `order`.
pub fn skbio_alpha_tsallis(counts : Array[Double], order : Double) -> Double {
  let c = skbio_alpha_nonzero(counts)
  let s = c.length()
  if s == 0 {
    return skbio_alpha_nan()
  }
  if s == 1 {
    return 0.0
  }
  let n = skbio_alpha_sum(c)
  let probs = skbio_alpha_probs(c, n)
  if order == 1.0 {
    skbio_alpha_entropy(probs)
  } else if order == 1.0 / 0.0 {
    0.0
  } else {
    let mut sum = 0.0
    for p in probs {
      sum = sum + @math.pow(p, order)
    }
    (1.0 - sum) / (order - 1.0)
  }
}

///|
/// Hill number (effective number of species) of order `order`.
pub fn skbio_alpha_hill(counts : Array[Double], order : Double) -> Double {
  let c = skbio_alpha_nonzero(counts)
  let s = c.length()
  if s == 0 {
    return skbio_alpha_nan()
  }
  let n = skbio_alpha_sum(c)
  let probs = skbio_alpha_probs(c, n)
  if order == 1.0 {
    @math.exp(skbio_alpha_entropy(probs))
  } else if order == 1.0 / 0.0 {
    let mut mx = 0.0
    for p in probs {
      if p > mx {
        mx = p
      }
    }
    1.0 / mx
  } else {
    let mut sum = 0.0
    for p in probs {
      sum = sum + @math.pow(p, order)
    }
    @math.pow(sum, 1.0 / (1.0 - order))
  }
}

// ---------------------------------------------------------------------------
// Shannon / Chao1 / ACE（最常用的 alpha 多样性）
// ---------------------------------------------------------------------------

///|
/// Shannon-Wiener diversity index (`skbio.diversity.alpha.shannon`, base e).
/// H = -Σ pᵢ ln(pᵢ).  Empty input returns NaN (scikit-bio's convention).
pub fn skbio_alpha_shannon(counts : Array[Double]) -> Double {
  let c = skbio_alpha_nonzero(counts)
  if c.length() == 0 {
    return skbio_alpha_nan()
  }
  let n = skbio_alpha_sum(c)
  let probs = skbio_alpha_probs(c, n)
  skbio_alpha_entropy(probs)
}

///|
/// Chao1 richness estimator (`skbio.diversity.alpha.chao1`):
/// S_obs + f₁(f₁-1) / (2 f₂), where f₁ = singletons, f₂ = doubletons.
/// When f₂ = 0, falls back to S_obs + f₁(f₁-1)/2.
pub fn skbio_alpha_chao1(counts : Array[Double]) -> Double {
  let c = skbio_alpha_nonzero(counts)
  let s_obs = c.length().to_double()
  let f1 = skbio_alpha_singles(c).to_double()
  let f2 = skbio_alpha_doubles(c).to_double()
  if f2 > 0.0 {
    s_obs + f1 * (f1 - 1.0) / (2.0 * f2)
  } else if f1 > 0.0 {
    s_obs + f1 * (f1 - 1.0) / 2.0
  } else {
    s_obs
  }
}

///|
/// ACE richness estimator (`skbio.diversity.alpha.ace`). Counts taxa with
/// abundance ≤ `threshold` as "rare" (default 10), others as "abundant".
/// C_ACE = 1 - f₁ / N_rare;  S_ACE = S_abund + S_rare / C_ACE.
pub fn skbio_alpha_ace(counts : Array[Double], threshold : Int32) -> Double {
  let c = skbio_alpha_nonzero(counts)
  let t = threshold.to_double()
  let mut s_abund = 0.0
  let mut s_rare = 0.0
  let mut n_rare = 0.0
  let mut f1 = 0.0
  for v in c {
    if v <= t {
      s_rare = s_rare + 1.0
      n_rare = n_rare + v
      if v == 1.0 {
        f1 = f1 + 1.0
      }
    } else {
      s_abund = s_abund + 1.0
    }
  }
  let c_ace = if n_rare > 0.0 { 1.0 - f1 / n_rare } else { 0.0 }
  if c_ace > 0.0 {
    s_abund + s_rare / c_ace
  } else {
    s_abund + s_rare
  }
}

///|
/// Default ACE with threshold=10 (matches scikit-bio).
pub fn skbio_alpha_ace_default(counts : Array[Double]) -> Double {
  skbio_alpha_ace(counts, 10)
}

// ---------------------------------------------------------------------------
// pyramid / vandenbogaert（额外的覆盖度补充指标）
// ---------------------------------------------------------------------------

///|
/// Pyramid coverage (`skbio.diversity.alpha.pyramid`):
/// pyramid = S_obs / (S_obs + f₁)，where f₁ = singletons (count==1).
/// 衡量"我们已知的丰富度占已知+未知的比例"，接近 1 时说明样本已饱和。
pub fn skbio_alpha_pyramid(counts : Array[Double]) -> Double {
  let c = skbio_alpha_nonzero(counts)
  let s_obs = c.length().to_double()
  let f1 = skbio_alpha_singles(c).to_double()
  if s_obs + f1 > 0.0 {
    s_obs / (s_obs + f1)
  } else {
    0.0 / 0.0
  }
}

///|
/// Vandenbogaert's coverage (`skbio.diversity.alpha.vandenbogaert`):
/// C_vdb = 1 - f₁² / (f₂ + f₁²)，用 singletons 和 doubletons 的比例估计
/// 覆盖度。比 goods_coverage 更偏向稀有物种信息。
pub fn skbio_alpha_vandenbogaert(counts : Array[Double]) -> Double {
  let c = skbio_alpha_nonzero(counts)
  let f1 = skbio_alpha_singles(c).to_double()
  let f2 = skbio_alpha_doubles(c).to_double()
  let denom = f2 + f1 * f1
  if denom > 0.0 {
    1.0 - f1 * f1 / denom
  } else if f1 > 0.0 {
    0.0 // f₂=0 且有 singleton → 最差覆盖
  } else {
    1.0 // 没有 singleton → 完美覆盖（无未知物种）
  }
}

// ============================================================================
// subsample_power（测序饱和度功效曲线）
// ============================================================================

///|
/// Sampling-depth power curve (`skbio.stats.power.subsample_power`, with the
/// `observed_otus` metric): expected observed species at increasing depths,
/// via without-replacement rarefaction. Returns `(depth, expected_species)`.
pub fn skbio_subsample_power(
  counts : Array[Double],
  minimum : Int32,
  maximum : Int32,
  interval : Int32,
) -> Array[(Int32, Double)] {
  let out : Array[(Int32, Double)] = Array(capacity=0)
  let mut d = minimum
  while d <= maximum {
    out.push((d, vegan_rarefy(counts, d)))
    d = d + interval
  }
  out
}
