///|
/// scikit-bio alpha diversity metrics.
///
/// Port of the formula-based alpha diversity metrics from
/// `skbio.diversity.alpha` (scikit-bio 0.7.3) that were not already covered by
/// BioSeqs' Bioconductor (vegan) surface. Each metric takes a vector of taxon
/// counts; zero counts are dropped first (absent taxa), matching scikit-bio's
/// `_validate_alpha` wrapper. Metrics that scikit-bio marks `empty=np.nan`
/// return NaN for an all-zero community; count metrics return 0.
///
/// Not yet included: `finite=True` variants, `gini_index`,
/// `kempton_taylor_q`, `esty_ci`/`lladser_ci`/`lladser_pe` (trie-based), and
/// `ace`/`chao1`/`fisher_alpha` (already present in BioSeqs).

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
fn skbio_alpha_log_factorial(n : Int) -> Double {
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
pub fn skbio_alpha_sobs(counts : Array[Double]) -> Int {
  skbio_alpha_nonzero(counts).length()
}

///|
/// Number of singletons (taxa observed exactly once).
pub fn skbio_alpha_singles(counts : Array[Double]) -> Int {
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
pub fn skbio_alpha_doubles(counts : Array[Double]) -> Int {
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
pub fn skbio_alpha_osd(counts : Array[Double]) -> (Int, Int, Int) {
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
  let denom = ((n - s.to_double() + 1.0) * (n - s.to_double() + 1.0) + s.to_double() - 1.0).sqrt()
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