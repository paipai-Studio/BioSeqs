///|
/// RankProd - rank-product differential expression, ported from the
/// Bioconductor RankProd package (Hong, Breitling, McEntee, Wittner,
/// Chholiya & Benner; exact p-values by Eisinga, Breitling et al. 2017).
///
/// The rank product/sum is computed from every (or a random selection of)
/// cross-class sample pairing: per pair, features are ranked by log fold
/// change; the feature statistic is the geometric mean (rank product) or
/// arithmetic mean (rank sum) of its ranks. Analytical, permutation-free
/// p-values are evaluated from the exact rank-product bounds distribution
/// (``rankprodbounds``) or the dice-sum distribution for rank sums.
/// Multiple-testing results are reported as predicted false positives
/// (``pfp``), the rank-product analogue of the FDR.
///
/// The null resampling draws within-feature Fisher-Yates permutations with a
/// portable 32-bit Numerical-Recipes LCG, so identical ``seed`` values give
/// bit-identical results across R and MoonBit (see the frozen-value generator
/// /tmp/gen_rankprod_expected.R, which injects the same sampler into the real
/// package closure).
///
/// Matrices follow the package convention: features (genes) are rows and
/// samples are columns. Class labels use 0/1; a single class triggers the
/// paired one-sample analysis.

///|
/// Error raised by the rank-product functions.
pub suberror RankProdError {
  RankProdError(String)
}

///|
/// One side of a ``rankprod_top_genes`` selection (class1 < class2 or
/// class1 > class2). Indices are 0-based feature rows. Numeric columns mirror
/// the package ``topGene`` table, rounded to four significant digits as the
/// R function does via ``prettyNum``.
pub struct RankProdTopTable {
  indices : Array[Int]
  statistic : Array[Double]
  fold_change : Array[Double]
  pfp : Array[Double]
  p_value : Array[Double]
}

///|
/// Selected feature rows (0-based).
pub fn RankProdTopTable::indices(self : RankProdTopTable) -> Array[Int] {
  self.indices
}

///|
/// Rank product / rank sum statistic of the selected rows.
pub fn RankProdTopTable::statistic(self : RankProdTopTable) -> Array[Double] {
  self.statistic
}

///|
/// Fold change class1/class2 (``logbase^AveFC`` for logged input).
pub fn RankProdTopTable::fold_change(self : RankProdTopTable) -> Array[Double] {
  self.fold_change
}

///|
/// Predicted false-positive percentage of the selected rows.
pub fn RankProdTopTable::pfp(self : RankProdTopTable) -> Array[Double] {
  self.pfp
}

///|
/// Exact p-values of the selected rows.
pub fn RankProdTopTable::p_value(self : RankProdTopTable) -> Array[Double] {
  self.p_value
}

///|
/// Result of ``rank_products``. Arrays are indexed by feature row; features
/// that fail the ``min_valid_pairs`` threshold carry NaN values.
pub struct RankProdResult {
  statistic_down : Array[Double]
  statistic_up : Array[Double]
  p_value_down : Array[Double]
  p_value_up : Array[Double]
  pfp_down : Array[Double]
  pfp_up : Array[Double]
  rank_down : Array[Double]
  rank_up : Array[Double]
  average_fold_change : Array[Double]
  /// Geometric-mean (RP) / arithmetic-mean (RS) products of every pairing,
  /// rows in [observed, null draws ...], down direction (package
  /// RandomPairs_ranks).
  products_down : Array[Array[Double]]
  groups : Array[Int]
  n_features : Int
  nrep : Int
  n_valid : Int
  random_pairs : Int
  calculate_product : Bool
  logged : Bool
}

///|
/// Rank-product statistic (geometric mean of ranks) of every feature for the
/// class1 < class2 direction.
pub fn RankProdResult::statistic_down(self : RankProdResult) -> Array[Double] {
  self.statistic_down
}

///|
/// Rank-product statistic of the class1 > class2 (up-regulated) direction.
pub fn RankProdResult::statistic_up(self : RankProdResult) -> Array[Double] {
  self.statistic_up
}

///|
/// Exact analytical p-values, class1 < class2.
pub fn RankProdResult::p_value_down(self : RankProdResult) -> Array[Double] {
  self.p_value_down
}

///|
/// Exact analytical p-values, class1 > class2.
pub fn RankProdResult::p_value_up(self : RankProdResult) -> Array[Double] {
  self.p_value_up
}

///|
/// Predicted false-positive percentage (rank-product FDR analogue), down.
pub fn RankProdResult::pfp_down(self : RankProdResult) -> Array[Double] {
  self.pfp_down
}

///|
/// Predicted false-positive percentage, up.
pub fn RankProdResult::pfp_up(self : RankProdResult) -> Array[Double] {
  self.pfp_up
}

///|
/// R-style average ranks of the statistic (ties averaged), down direction.
pub fn RankProdResult::rank_down(self : RankProdResult) -> Array[Double] {
  self.rank_down
}

///|
/// R-style average ranks of the statistic (ties averaged), up direction.
pub fn RankProdResult::rank_up(self : RankProdResult) -> Array[Double] {
  self.rank_up
}

///|
/// Average fold change class1/class2 (log-unit difference when ``logged``,
/// raw ratio otherwise).
pub fn RankProdResult::average_fold_change(self : RankProdResult) -> Array[Double] {
  self.average_fold_change
}

///|
/// Per-pair products of ranks (random_pairs x n_features), down direction.
pub fn RankProdResult::products_down(self : RankProdResult) -> Array[Array[Double]] {
  self.products_down
}

///|
/// Remapped class labels per sample (0/1, or all 1 in the paired case).
pub fn RankProdResult::groups(self : RankProdResult) -> Array[Int] {
  self.groups
}

///|
/// Number of feature rows.
pub fn RankProdResult::n_features(self : RankProdResult) -> Int {
  self.n_features
}

///|
/// Number of within-pair comparisons the statistic averages over (size of the
/// smaller class, or all samples in the one-class paired case).
pub fn RankProdResult::nrep(self : RankProdResult) -> Int {
  self.nrep
}

///|
/// Number of features retained after the ``min_valid_pairs`` exclusion.
pub fn RankProdResult::n_valid(self : RankProdResult) -> Int {
  self.n_valid
}

///|
/// Number of (observed + null) pairings used.
pub fn RankProdResult::random_pairs(self : RankProdResult) -> Int {
  self.random_pairs
}

///|
/// True when the result is a rank product, false for a rank sum.
pub fn RankProdResult::is_rank_product(self : RankProdResult) -> Bool {
  self.calculate_product
}

// ---------------------------------------------------------------------------
// Portable RNG: Numerical Recipes 32-bit LCG + Fisher-Yates. The R frozen
// value generator uses the byte-identical construction.
// ---------------------------------------------------------------------------

struct RpLcg {
  mut state : Int64
}

fn rp_lcg_init(seed : Int) -> RpLcg {
  let m = 4294967296L // 2^32
  let s0 = seed.to_int64()
  let s = if s0 < 0 {
    (s0 % m + m) % m
  } else {
    s0 % m
  }
  { state: s }
}

fn rp_lcg_next(rng : RpLcg) -> Int64 {
  // 1664525 * state stays below 2^53, so the arithmetic is exact in double on
  // the R side as well.
  rng.state = (1664525L * rng.state + 1013904223L) & 4294967295L
  rng.state
}

/// Fisher-Yates permutation of ``values`` (already sorted ascending).
fn rp_lcg_shuffle(values : Array[Double], rng : RpLcg) -> Array[Double] {
  let n = values.length()
  if n <= 1 {
    return values.copy()
  }
  let out = values.copy()
  let mut i = n - 1
  while i >= 1 {
    let j = (rp_lcg_next(rng) % (i + 1).to_int64()).to_int()
    let tmp = out[i]
    out[i] = out[j]
    out[j] = tmp
    i -= 1
  }
  out
}

// ---------------------------------------------------------------------------
// Small numerical helpers
// ---------------------------------------------------------------------------

fn rp_nan() -> Double {
  0.0 / 0.0
}

fn rp_is_nan(x : Double) -> Bool {
  x.is_nan()
}

/// Stable ascending index order: ties keep the original index order; NaNs
/// sort last. Mirrors R's stable order() / sort.list() ties.
fn rp_stable_order(values : Array[Double]) -> Array[Int] {
  let idx : Array[Int] = []
  for i in 0..<values.length() {
    if !values[i].is_nan() {
      idx.push(i)
    }
  }
  // sort_by is documented unstable, but the index tiebreak makes every
  // comparison key unique, so the result is deterministic either way.
  idx.sort_by(fn(a, b) {
    let va = values[a]
    let vb = values[b]
    if va < vb {
      -1
    } else if va > vb {
      1
    } else {
      a - b
    }
  })
  for i in 0..<values.length() {
    if values[i].is_nan() {
      idx.push(i)
    }
  }
  idx
}

/// The reference implementation composes the per-pairing geometric mean as
/// ``exp(mean(log(ranks)))``; its libm rounds different summation orders 1-2
/// ULP apart (and a different libm differs again), even when the underlying
/// integer rank products are mathematically identical. R's ``order()`` then
/// sees bit-equal p-values/RPs and resolves them by ascending feature index.
/// Values within this numerical-noise band are therefore treated as tied.
let rp_tie_abs : Double = 1.0e-12
let rp_tie_rel : Double = 1.0e-12

/// Near-equality for cross-language tie detection: absolute tolerance for
/// p-values in [0, 1], relative tolerance for rank statistics >= 1.
fn rp_near(va : Double, vb : Double) -> Bool {
  let d = (va - vb).abs()
  if d <= rp_tie_abs {
    return true
  }
  let scale = if va.abs() > vb.abs() { va.abs() } else { vb.abs() }
  d <= rp_tie_rel * scale
}

/// Stable ascending order with ``rp_near`` ties resolved by index, mirroring
/// R ``order()`` on values that are equal up to libm noise.
fn rp_stable_order_tol(values : Array[Double]) -> Array[Int] {
  let idx : Array[Int] = []
  for i in 0..<values.length() {
    if !values[i].is_nan() {
      idx.push(i)
    }
  }
  idx.sort_by(fn(a, b) {
    let va = values[a]
    let vb = values[b]
    if rp_near(va, vb) {
      a - b
    } else if va < vb {
      -1
    } else {
      1
    }
  })
  idx
}

/// Position ranks 1..m of the finite values (smallest = 1), stable on ties;
/// NaN entries receive NaN. This is the package `output[order(v)] <- 1:m`
/// construction (NOT average ranks).
fn rp_position_ranks(column : Array[Double]) -> Array[Double] {
  let n = column.length()
  let order = rp_stable_order(column)
  let ranks : Array[Double] = Array::make(n, rp_nan())
  let mut pos = 1
  for i in order {
    if rp_is_nan(column[i]) {
      continue // missing positions stay NaN
    }
    ranks[i] = pos.to_double()
    pos += 1
  }
  ranks
}

/// R's default rank() with ties averaged (ties.method = "average"); NaNs are
/// left as NaN. With ``tol=true``, values within ``rp_near`` noise are treated
/// as tied, matching R's exact-double ties across libm implementations.
fn rp_average_ranks(values : Array[Double], tol~ : Bool = false) -> Array[Double] {
  let n = values.length()
  let ranks : Array[Double] = Array::make(n, rp_nan())
  let order = if tol {
    rp_stable_order_tol(values)
  } else {
    rp_stable_order(values)
  }
  let mut k = 0
  let total = order.length()
  while k < total {
    let i = order[k]
    if rp_is_nan(values[i]) {
      k += 1
      continue
    }
    let same = fn(u : Double, v : Double) -> Bool {
      if tol {
        rp_near(u, v)
      } else {
        u == v
      }
    }
    let mut j = k
    while j + 1 < total && !rp_is_nan(values[order[j + 1]]) &&
      same(values[order[j + 1]], values[i]) {
      j += 1
    }
    // positions k+1 .. j+1 (1-based) are tied
    let avg = ((k + 1) + (j + 1)).to_double() / 2.0
    let mut t = k
    while t <= j {
      ranks[order[t]] = avg
      t += 1
    }
    k = j + 1
  }
  ranks
}

/// Median of finite values, NaN removed (R median(x, na.rm = TRUE)); NaN when
/// no finite value exists.
fn rp_median_na(values : Array[Double]) -> Double {
  let finite : Array[Double] = values.filter(fn(v) { !rp_is_nan(v) })
  let n = finite.length()
  if n == 0 {
    return rp_nan()
  }
  finite.sort()
  if n % 2 == 1 {
    return finite[n / 2]
  }
  (finite[n / 2 - 1] + finite[n / 2]) / 2.0
}

/// Arithmetic mean of finite values.
fn rp_mean_na(values : Array[Double]) -> Double {
  let mut acc = 0.0
  let mut count = 0
  for v in values {
    if !rp_is_nan(v) {
      acc += v
      count += 1
    }
  }
  if count == 0 {
    return rp_nan()
  }
  acc / count.to_double()
}

/// Stable ascending sort of finite values (R sort() drops NaNs).
fn rp_sorted(values : Array[Double]) -> Array[Double] {
  let finite : Array[Double] = values.filter(fn(v) { !rp_is_nan(v) })
  finite.sort()
  finite
}

/// log10 of a raw ratio; non-positive or missing operands yield NaN.
fn rp_log10_ratio(v1 : Double, v2 : Double) -> Double {
  if rp_is_nan(v1) || rp_is_nan(v2) || v2 == 0.0 {
    return rp_nan()
  }
  let r = v1 / v2
  if r <= 0.0 {
    return rp_nan()
  }
  @math.log10(r)
}

/// log10 of a single raw value (paired case).
fn rp_log10_value(v : Double) -> Double {
  if rp_is_nan(v) || v <= 0.0 {
    return rp_nan()
  }
  @math.log10(v)
}

// ---------------------------------------------------------------------------
// Exact rank-product p-values: rankprodbounds (Eisinga et al. 2017)
// ---------------------------------------------------------------------------

struct RpParam {
  a : Array[Double]
  b : Array[Double]
  c : Array[Double]
  d : Double
  e : Double
}

/// Combine coefficient rows sharing the same (a, b) pair by summing c,
/// keeping first-occurrence order (package makeunique).
fn rp_make_unique(p : RpParam) -> RpParam {
  let ua : Array[Double] = []
  let ub : Array[Double] = []
  let uc : Array[Double] = []
  for i in 0..<p.a.length() {
    let mut found = -1
    for j in 0..<ua.length() {
      if ua[j] == p.a[i] && ub[j] == p.b[i] {
        found = j
        break
      }
    }
    if found < 0 {
      ua.push(p.a[i])
      ub.push(p.b[i])
      uc.push(p.c[i])
    } else {
      uc[found] += p.c[i]
    }
  }
  { a: ua, b: ub, c: uc, d: p.d, e: p.e }
}

/// Memoized coefficient recursion (package updateparam). ``delta`` is 0
/// ("lower") or 1 ("upper").
fn rp_update_param(
  cache : Array[Array[Option[RpParam]]],
  n : Double,
  k : Int,
  j : Int,
  delta : Double,
) -> RpParam {
  match cache[k][j] {
    Some(p) => p
    None => {
      let logn = @math.ln(n)
      let p = if j == 0 {
        { a: [], b: [], c: [], d: 0.0, e: @math.pow(n, k.to_double()) }
      } else if k == j {
        let p00 = rp_update_param(cache, n, k - 1, j - 1, delta)
        let newa0 = p00.a.map(fn(x) { x + 1.0 })
        let newc0 = p00.a.mapi(fn(i, a) { p00.c[i] / (a + 1.0) })
        let e = (1.0 - delta) * (1.0 - p00.e)
        let d = delta * p00.d + p00.e
        let q = {
          a: [1.0] + p00.a + newa0,
          b: [0.0] + p00.b + p00.b.copy(),
          c: [p00.d] + p00.c.map(fn(x) { delta * x }) + newc0,
          d,
          e,
        }
        rp_make_unique(q)
      } else {
        let p00 = rp_update_param(cache, n, k - 1, j - 1, delta)
        let p01 = rp_update_param(cache, n, k - 1, j, delta)
        let newa0 = p00.a.map(fn(x) { x + 1.0 })
        let newc0 = p00.a.mapi(fn(i, a) { p00.c[i] / (a + 1.0) })
        let newa1 = p01.a.map(fn(x) { x + 1.0 })
        let newc1 = p01.a.mapi(fn(i, a) { -p01.c[i] / (a + 1.0) })
        let lognnkj = ((k - j).to_double()) * logn
        let e = n * p01.e + (delta - 1.0) * (p00.e - p01.e)
        // lognminb aligned with [newa0..., newa1...]
        let lognminb : Array[Double] = []
        for b in p00.b {
          lognminb.push(-b * logn)
        }
        for b in p01.b {
          lognminb.push((1.0 - b) * logn)
        }
        let newa = newa0 + newa1
        let newc = newc0 + newc1
        let mut correction = 0.0
        for i in 0..<newa.length() {
          correction += newc[i] * @math.pow(lognminb[i], newa[i])
        }
        let d = delta * p00.d + (1.0 - delta) * p01.d / n +
          (p00.e - p01.e) / @math.exp(lognnkj) - correction
        let q = {
          a: [1.0, 1.0] + p00.a + p01.a + newa,
          b: [0.0, 1.0] + p00.b + p01.b + p00.b.copy() + p01.b.copy(),
          c: [p00.d, -p01.d] + p00.c.map(fn(x) { delta * x }) +
            p01.c.map(fn(x) { (1.0 - delta) * x / n }) + newc,
          d,
          e,
        }
        rp_make_unique(q)
      }
      cache[k][j] = Some(p)
      p
    }
  }
}

/// Exact CDF of the rank product at ``rho`` (the product of ranks) for ``n``
/// features and ``k`` paired comparisons under bound ``delta`` (0 = lower,
/// 1 = upper).
fn rp_rank_product_bound(
  rho : Double,
  n : Int,
  k : Int,
  delta : Double,
) -> Double raise RankProdError {
  if rho < 1.0 || rho > @math.pow(n.to_double(), k.to_double()) + 1.0 {
    raise RankProdError(
      "rank_product_pvalue: rank product out of bounds [1, n^k]",
    )
  }
  let nd = n.to_double()
  let logn = @math.ln(nd)
  let cache : Array[Array[Option[RpParam]]] = Array::make(k + 1, [])
  for kk in 0..<(k + 1) {
    cache[kk] = Array::make(kk + 1, None)
  }
  // allj = ceil(-log(rho)/logn + k)
  let jraw = -@math.ln(rho) / logn + k.to_double()
  let j = @math.ceil(jraw).to_int()
  let g = if j == 0 {
    let p = rp_update_param(cache, nd, k, 0, delta)
    p.e
  } else {
    let p = rp_update_param(cache, nd, k, j, delta)
    let mut total = p.e + p.d * rho
    for m in 0..<p.a.length() {
      let base = @math.ln(rho) -
        logn * (k.to_double() - j.to_double() + p.b[m])
      total += p.c[m] * rho * @math.pow(base, p.a[m])
    }
    total
  }
  g / @math.pow(nd, k.to_double())
}

///|
/// Exact rank-product p-value (geometric mean of the upper and lower bounds,
/// package ``rankprodbounds(..., Delta = "geometric")``). ``rank`` is the
/// geometric-mean rank statistic, ``n`` the number of features and ``k`` the
/// number of paired comparisons.
pub fn rank_product_pvalue(
  rank : Double,
  n : Int,
  k : Int,
) -> Double raise RankProdError {
  let (upper, lower, _) = rank_product_pvalue_bounds(rank, n, k)
  (upper * lower).sqrt()
}

///|
/// Exact rank-product p-value bounds. Returns ``(upper, lower, geometric)``
/// for the geometric-mean rank ``rank``; the geometric value is the package
/// default ``Delta = "geometric"`` p-value.
pub fn rank_product_pvalue_bounds(
  rank : Double,
  n : Int,
  k : Int,
) -> (Double, Double, Double) raise RankProdError {
  let rho = @math.pow(rank, k.to_double())
  let upper = rp_rank_product_bound(rho, n, k, 1.0)
  let lower = rp_rank_product_bound(rho, n, k, 0.0)
  let geom = (upper * lower).sqrt()
  (upper, lower, geom)
}

// ---------------------------------------------------------------------------
// Exact rank-sum p-values: dice.sum.cdf/pdf (sum of n dice, each 1..s)
// ---------------------------------------------------------------------------

/// Stirling approximation to lgamma for positive reals; matches base R well
/// enough for the non-integer ``choose`` fallback (all package call sites use
/// integer arguments).
fn rp_lgamma(x_in : Double) -> Double {
  let mut x = x_in
  if x <= 0.0 {
    return 0.0
  }
  let mut shift = 0.0
  while x < 7.0 {
    shift -= @math.ln(x)
    x += 1.0
  }
  let ln_2pi = 0.9189385332046727
  let inv_12x = 1.0 / (12.0 * x)
  let inv_360x3 = 1.0 / (360.0 * x * x * x)
  shift + ln_2pi + (x - 0.5) * @math.ln(x) - x + inv_12x - inv_360x3
}

/// Binomial coefficient C(n, k) matching R choose(): exact multiplicative
/// formula for integer n, gamma ratio otherwise.
fn rp_choose(n_in : Double, k : Int) -> Double {
  if k < 0 {
    return 0.0
  }
  let is_int = n_in == n_in.to_int64().to_double() &&
    n_in.abs() < 9.007199254740992e15
  if is_int {
    let ni = n_in.to_int()
    if ni >= 0 {
      if k > ni {
        return 0.0
      }
      let mut result = 1.0
      for i in 1..<(k + 1) {
        result = result * (ni - k + i).to_double() / i.to_double()
      }
      return result
    }
  }
  @math.exp(
    rp_lgamma(n_in + 1.0) - rp_lgamma(k.to_double() + 1.0) -
      rp_lgamma(n_in - k.to_double() + 1.0),
  )
}

///|
/// Exact CDF of the sum of ``n`` fair dice labelled 1..``s`` at total ``p``
/// (package ``dice.sum.cdf``): P(sum <= p), evaluated by inclusion-exclusion
/// in ordinary double precision.
pub fn rank_sum_pvalue(p : Double, n : Int, s : Int) -> Double {
  let a = (n * (s + 1) / 2).to_double() // trunc(n*(s+1)/2), positive domain
  let snd = s.to_double()
  let nnd = n.to_double()
  let denom = @math.pow(snd, nnd)
  if p <= a {
    let mm = @math.floor((p - nnd) / snd).to_int()
    let mut total = 0.0
    for i in 0..<(mm + 1) {
      let term = rp_choose(p - snd * i.to_double(), n)
      let sign = if i % 2 == 0 { 1.0 } else { -1.0 }
      total += sign * rp_choose(nnd, i) * term
    }
    total / denom
  } else {
    let mm = @math.floor((snd * nnd - p) / snd).to_int()
    let mut total = 0.0
    for i in 0..<(mm + 1) {
      let top = nnd * (snd + 1.0) - p - snd * i.to_double() - 1.0
      let term = rp_choose(top, n)
      let sign = if i % 2 == 0 { 1.0 } else { -1.0 }
      total += sign * rp_choose(nnd, i) * term
    }
    1.0 - total / denom
  }
}

///|
/// Exact probability mass of the sum of ``n`` dice 1..``s`` at total ``p``
/// (package ``dice.sum.pdf``).
pub fn rank_sum_pmf(p : Double, n : Int, s : Int) -> Double {
  let snd = s.to_double()
  let nnd = n.to_double()
  let denom = @math.pow(snd, nnd)
  let mm = @math.floor((p - nnd) / snd).to_int()
  let mut total = 0.0
  for i in 0..<(mm + 1) {
    let term = rp_choose(p - snd * i.to_double() - 1.0, n - 1)
    let sign = if i % 2 == 0 { 1.0 } else { -1.0 }
    total += sign * rp_choose(nnd, i) * term
  }
  total / denom
}

// ---------------------------------------------------------------------------
// Main entry point
// ---------------------------------------------------------------------------

///|
/// Rank product / rank sum analysis, mirroring Bioconductor
/// ``RankProducts``.
///
/// - ``data`` is feature-major (rows = genes, columns = samples).
/// - ``cl`` class labels per sample: two distinct values remapped to 0/1;
///   a single class selects the paired one-sample test.
/// - ``logged``: are the data already on a log scale? Otherwise raw ratios
///   are formed and rank-scored on the log10 scale (zeros are replaced by the
///   smallest positive value, as in the package).
/// - ``calculate_product``: true = rank product (geometric mean ranks);
///   false = rank sum (arithmetic mean ranks).
/// - ``random_pairs``: number of observed + null pairings. When negative the
///   package default ``npairs^2`` (rounded up to odd) is used; the one-class
///   case always uses one.
/// - ``seed`` drives the portable LCG null sampler.
/// - ``na_rm``: impute missing values with the pooled per-feature median.
/// - ``min_valid_pairs``: features with fewer finite values than this
///   threshold (default ``floor(nrep/2)``) in either class are excluded.
pub fn rank_products(
  data~ : Array[Array[Double]],
  cl~ : Array[Int],
  logged~ : Bool = true,
  calculate_product~ : Bool = true,
  random_pairs~ : Int = -1,
  seed~ : Int = 42,
  na_rm~ : Bool = true,
  min_valid_pairs~ : Int = -1,
) -> RankProdResult raise RankProdError {
  let ng = data.length()
  if ng == 0 {
    raise RankProdError("rank_products: empty data matrix")
  }
  let ns = cl.length()
  if ns == 0 {
    raise RankProdError("rank_products: empty class labels")
  }
  for row in data {
    if row.length() != ns {
      raise RankProdError(
        "rank_products: number of class labels must match data columns",
      )
    }
  }

  // Class labels -> unique values, remapped like xyCall.
  let levels : Array[Int] = []
  for c in cl {
    if !levels.contains(c) {
      levels.push(c)
    }
  }
  levels.sort()
  if levels.length() > 2 {
    raise RankProdError("rank_products: more than two class labels")
  }
  let groups = Array::make(ns, 1)
  let x_cols : Array[Int] = []
  let y_cols : Array[Int] = []
  let two_class = levels.length() == 2
  if two_class {
    let lo = levels[0]
    for i in 0..<ns {
      if cl[i] == lo {
        groups[i] = 0
        x_cols.push(i)
      } else {
        groups[i] = 1
        y_cols.push(i)
      }
    }
  } else {
    for i in 0..<ns {
      groups[i] = 1
      x_cols.push(i)
    }
  }
  let n1 = x_cols.length()
  let n2 = if two_class {
    y_cols.length()
  } else {
    0
  }
  let npairs = if n1 < n2 {
    n1
  } else {
    n2
  }
  let nrep = if npairs == 0 {
    ns
  } else {
    npairs
  }

  // RandomPairs bookkeeping (0 -> 1, even -> +1).
  let mut rpairs = if random_pairs < 0 {
    npairs * npairs
  } else {
    random_pairs
  }
  if rpairs == 0 {
    rpairs = 1
  }
  if rpairs % 2 == 0 {
    rpairs += 1
  }
  if !two_class {
    rpairs = 1
  }

  // Work on a mutable copy; mzmatch: +/- Inf and NaN both become NaN.
  let work : Array[Array[Double]] = data.map(fn(row) {
    row.map(fn(v) {
      if v.is_nan() || v.is_inf() {
        rp_nan()
      } else {
        v
      }
    })
  })

  // detection.limit = smallest non-zero finite value.
  let mut detection_limit = rp_nan()
  for row in work {
    for v in row {
      if !v.is_nan() && v != 0.0 {
        if detection_limit.is_nan() || v < detection_limit {
          detection_limit = v
        }
      }
    }
  }

  // Unlogged raw data: zero entries are replaced before group splitting.
  if !logged && !detection_limit.is_nan() {
    for row in work {
      for i in 0..<row.length() {
        if row[i] == 0.0 {
          row[i] = detection_limit
        }
      }
    }
  }

  let gr1 : Array[Array[Double]] = [] // feature-major values per class
  let gr2 : Array[Array[Double]] = []
  for g in 0..<ng {
    let r1 : Array[Double] = []
    for c in x_cols {
      r1.push(work[g][c])
    }
    gr1.push(r1)
    if two_class {
      let r2 : Array[Double] = []
      for c in y_cols {
        r2.push(work[g][c])
      }
      gr2.push(r2)
    }
  }

  // badrows / detection-limit replacement, only with active resampling.
  let min_valid = if min_valid_pairs < 0 {
    nrep / 2
  } else {
    min_valid_pairs
  }
  let badrows : Array[Int] = []
  if rpairs > 1 {
    for g in 0..<ng {
      let na1 = gr1[g].filter(fn(v) { rp_is_nan(v) }).length()
      let na2 = if two_class {
        gr2[g].filter(fn(v) { rp_is_nan(v) }).length()
      } else {
        0
      }
      let valid1 = n1 - na1
      let valid2 = if two_class {
        n2 - na2
      } else {
        n1
      }
      let minv = if valid1 < valid2 {
        valid1
      } else {
        valid2
      }
      if min_valid > minv {
        badrows.push(g)
      }
      // Package replaces an all-NA group with the detection limit only when
      // the OTHER group still carries finite values.
      if two_class {
        let all1 = na1 == n1
        let all2 = na2 == n2
        if all1 && !all2 {
          for i in 0..<gr1[g].length() {
            gr1[g][i] = detection_limit
          }
        }
        if all2 && !all1 {
          for i in 0..<gr2[g].length() {
            gr2[g][i] = detection_limit
          }
        }
      }
    }
  }

  // NAreplace: pooled per-feature median; gr1 filled first, then gr2.
  // Filling with the running median leaves the pooled median unchanged, so a
  // single median per row is equivalent to R's sequential cell replacement.
  if na_rm {
    for g in 0..<ng {
      let pool1 : Array[Double] = []
      for v in gr1[g] {
        if !rp_is_nan(v) {
          pool1.push(v)
        }
      }
      if two_class {
        for v in gr2[g] {
          if !rp_is_nan(v) {
            pool1.push(v)
          }
        }
      }
      let med1 = if pool1.length() > 0 {
        rp_median_na(pool1)
      } else {
        rp_nan()
      }
      for i in 0..<gr1[g].length() {
        if rp_is_nan(gr1[g][i]) {
          gr1[g][i] = med1
        }
      }
      if two_class {
        let pool2 : Array[Double] = []
        for v in gr1[g] {
          pool2.push(v)
        }
        for v in gr2[g] {
          if !rp_is_nan(v) {
            pool2.push(v)
          }
        }
        let med2 = rp_median_na(pool2)
        for i in 0..<gr2[g].length() {
          if rp_is_nan(gr2[g][i]) {
            gr2[g][i] = med2
          }
        }
      }
    }
  }

  // Average fold change from the (possibly imputed) group values.
  let ave_fc : Array[Double] = Array::make(ng, rp_nan())
  for g in 0..<ng {
    let m1 = rp_mean_na(gr1[g])
    if two_class {
      let m2 = rp_mean_na(gr2[g])
      ave_fc[g] = if logged {
        m1 - m2
      } else {
        m1 / m2
      }
    } else {
      ave_fc[g] = m1
    }
    if badrows.contains(g) {
      ave_fc[g] = rp_nan()
    }
  }

  // --- Observed pairing (package ratios.list[[1]]) -------------------------
  let observed : Array[Array[Double]] = Array::make(ng, [])
  for g in 0..<ng {
    let row : Array[Double] = []
    if two_class {
      for p in 0..<nrep {
        let v1 = gr1[g][p]
        let v2 = gr2[g][p]
        row.push(if logged {
          v1 - v2
        } else {
          rp_log10_ratio(v1, v2)
        })
      }
    } else {
      for p in 0..<nrep {
        let v1 = gr1[g][p]
        row.push(if logged {
          v1
        } else {
          rp_log10_value(v1)
        })
      }
    }
    observed[g] = row
  }

  // --- Null pairings. R loops feature-outer, pairing-inner and consumes the
  // RNG in that exact order; the draw for group1 precedes group2 within each
  // pairing. The LCG is reset once for the whole analysis.
  let rng = rp_lcg_init(seed)
  // null_lfc[rp - 1] holds pairing rp (rp = 1..rpairs-1), feature-major.
  let null_lfc : Array[Array[Array[Double]]] = []
  if rpairs > 1 {
    for _ in 1..<rpairs {
      null_lfc.push(Array::make(ng, []))
    }
    for g in 0..<ng {
      let s1 = rp_sorted(gr1[g])
      let s2 = if two_class {
        rp_sorted(gr2[g])
      } else {
        []
      }
      for rp in 1..<rpairs {
        let p1 = rp_lcg_shuffle(s1, rng)
        let row : Array[Double] = []
        if two_class {
          let p2 = rp_lcg_shuffle(s2, rng)
          for p in 0..<nrep {
            let v1 = p1[p]
            let v2 = p2[p]
            row.push(if logged {
              v1 - v2
            } else {
              rp_log10_ratio(v1, v2)
            })
          }
        } else {
          for p in 0..<nrep {
            let v = p1[p]
            row.push(if logged {
              v
            } else {
              rp_log10_value(v)
            })
          }
        }
        null_lfc[rp - 1][g] = row
      }
    }
  }

  // Rank each pair column, then aggregate per feature; up ranks come from the
  // negated fold-change matrix (package inputmatrix.Class2 <- inputmatrix*(-1)).
  let prod_down : Array[Array[Double]] = []
  let stat_matrix_down : Array[Array[Double]] = [] // rpairs x ng
  let stat_matrix_up : Array[Array[Double]] = []
  for rp in 0..<rpairs {
    let lfc = if rp == 0 {
      observed
    } else {
      null_lfc[rp - 1]
    }
    for g in badrows {
      for c in 0..<lfc[g].length() {
        lfc[g][c] = rp_nan()
      }
    }
    let down_cols : Array[Array[Double]] = Array::make(nrep, [])
    let up_cols : Array[Array[Double]] = Array::make(nrep, [])
    for c in 0..<nrep {
      let column : Array[Double] = Array::make(ng, rp_nan())
      for g in 0..<ng {
        column[g] = lfc[g][c]
      }
      down_cols[c] = rp_position_ranks(column)
      // Both the paired one-sample and two-class paths negate the SAME
      // input matrix for the reverse direction (package line
      // inputmatrix.Class2 <- inputmatrix * (-1)); the earlier 1/gr1
      // assignment to ratios.list.Class2 is dead code in RankProducts.
      for g in 0..<ng {
        if !column[g].is_nan() {
          column[g] = -column[g]
        }
      }
      up_cols[c] = rp_position_ranks(column)
    }
    let row_down : Array[Double] = []
    let row_up : Array[Double] = []
    let stat_of = fn(cols : Array[Array[Double]], g : Int) -> Double {
      let mut acc = 0.0
      let mut count = 0
      for col in cols {
        let r = col[g]
        if !rp_is_nan(r) {
          if calculate_product {
            acc += @math.ln(r)
          } else {
            acc += r
          }
          count += 1
        }
      }
      if count == 0 {
        rp_nan()
      } else if calculate_product {
        @math.exp(acc / count.to_double())
      } else {
        acc / count.to_double()
      }
    }
    for g in 0..<ng {
      row_down.push(stat_of(down_cols, g))
      row_up.push(stat_of(up_cols, g))
    }
    stat_matrix_down.push(row_down)
    stat_matrix_up.push(row_up)
    prod_down.push(row_down)
  }

  // Median across pairings; cap at ng.
  let ng_d = ng.to_double()
  let statistic_down : Array[Double] = []
  let statistic_up : Array[Double] = []
  for g in 0..<ng {
    let col_down : Array[Double] = stat_matrix_down.map(fn(row) { row[g] })
    let col_up : Array[Double] = stat_matrix_up.map(fn(row) { row[g] })
    let mut sd = rp_median_na(col_down)
    let mut su = rp_median_na(col_up)
    if !rp_is_nan(sd) && sd > ng_d {
      sd = ng_d
    }
    if !rp_is_nan(su) && su > ng_d {
      su = ng_d
    }
    statistic_down.push(sd)
    statistic_up.push(su)
  }

  let n_valid = ng - badrows.length()

  // Analytical p-values.
  let pval_down : Array[Double] = Array::make(ng, rp_nan())
  let pval_up : Array[Double] = Array::make(ng, rp_nan())
  for g in 0..<ng {
    if !rp_is_nan(statistic_down[g]) {
      pval_down[g] = if calculate_product {
        rank_product_pvalue(statistic_down[g], n_valid, nrep)
      } else {
        rank_sum_pvalue(statistic_down[g] * nrep.to_double(), nrep, n_valid)
      }
      pval_up[g] = if calculate_product {
        rank_product_pvalue(statistic_up[g], n_valid, nrep)
      } else {
        rank_sum_pvalue(statistic_up[g] * nrep.to_double(), nrep, n_valid)
      }
    }
  }

  // Stable position ranks of the p-values -> pfp = p * nvalid / position.
  // Ties within libm noise keep ascending feature index, like R order().
  let pfp_from = fn(pv : Array[Double]) -> Array[Double] {
    let order = rp_stable_order_tol(pv)
    let positions : Array[Int] = Array::make(ng, -1)
    let mut pos = 1
    for i in order {
      if rp_is_nan(pv[i]) {
        continue
      }
      positions[i] = pos
      pos += 1
    }
    let out : Array[Double] = Array::make(ng, rp_nan())
    for g in 0..<ng {
      if positions[g] > 0 {
        out[g] = pv[g] * n_valid.to_double() / positions[g].to_double()
      }
    }
    out
  }
  let pfp_down = pfp_from(pval_down)
  let pfp_up = pfp_from(pval_up)

  // R average ranks of the statistic (bad rows NaN); tolerance ties absorb
  // libm-level differences in the geometric-mean products.
  let rank_down = rp_average_ranks(statistic_down, tol=true)
  let rank_up = rp_average_ranks(statistic_up, tol=true)
  for g in badrows {
    rank_down[g] = rp_nan()
    rank_up[g] = rp_nan()
  }

  {
    statistic_down,
    statistic_up,
    p_value_down: pval_down,
    p_value_up: pval_up,
    pfp_down,
    pfp_up,
    rank_down,
    rank_up,
    average_fold_change: ave_fc,
    products_down: prod_down,
    groups,
    n_features: ng,
    nrep,
    n_valid,
    random_pairs: rpairs,
    calculate_product,
    logged,
  }
}

// ---------------------------------------------------------------------------
// Orirank: all k1*k2 pairwise ranks (average-tie ranks, package computeOrirank)
// ---------------------------------------------------------------------------

///|
/// All-pairs original ranks. Returns ``(down, up)`` matrices feature-major
/// over the ``k1 * k2`` cross-class pairings; ranks use R's average-tie
/// convention. Unlogged inputs are ranked as raw ratios (as in the package),
/// not log ratios.
pub fn rankprod_ori_rank(
  gr1~ : Array[Array[Double]],
  gr2~ : Array[Array[Double]],
  logged~ : Bool = true,
) -> (Array[Array[Double]], Array[Array[Double]]) raise RankProdError {
  let ng = gr1.length()
  if ng == 0 || gr2.length() != ng {
    raise RankProdError("rankprod_ori_rank: feature rows must match")
  }
  let k1 = gr1[0].length()
  let k2 = gr2[0].length()
  if k1 == 0 || k2 == 0 {
    raise RankProdError("rankprod_ori_rank: empty class columns")
  }
  let build = fn(
    a : Array[Array[Double]],
    ka : Int,
    b : Array[Array[Double]],
    kb : Int,
  ) -> Array[Array[Double]] {
    let np = ka * kb
    let ranks : Array[Array[Double]] = []
    for _ in 0..<ng {
      ranks.push(Array::make(np, 0.0))
    }
    let mut c = 0
    for i in 0..<ka {
      for j in 0..<kb {
        let col = Array::make(ng, 0.0).mapi(fn(g, _) {
          if logged {
            a[g][i] - b[g][j]
          } else {
            a[g][i] / b[g][j]
          }
        })
        let rc = rp_average_ranks(col)
        for g in 0..<ng {
          ranks[g][c] = rc[g]
        }
        c += 1
      }
    }
    ranks
  }
  let down = build(gr1, k1, gr2, k2)
  let up = build(gr2, k2, gr1, k1)
  (down, up)
}

// ---------------------------------------------------------------------------
// topGene
// ---------------------------------------------------------------------------

/// R ``signif(x, 4)``: round to four significant digits; mirrors the
/// ``as.numeric(prettyNum(..., digits = 4))`` columns of the package topGene
/// tables.
fn rp_signif4(x : Double) -> Double {
  if x.is_nan() || x.is_inf() || x == 0.0 {
    return x
  }
  let ax = x.abs()
  let mag = @math.floor(@math.log10(ax)).to_int()
  let factor = @math.pow(10.0, (4 - 1 - mag).to_double())
  let scaled = x * factor
  let rounded = if scaled >= 0.0 {
    @math.floor(scaled + 0.5)
  } else {
    -@math.floor(-scaled + 0.5)
  }
  rounded / factor
}

/// Indices of finite statistics in ascending order, index tiebreak (R
/// ``sort(..., index.return = TRUE)`` drops NaNs). Statistics within libm
/// noise are treated as tied and keep ascending index order.
fn rp_top_order(stat : Array[Double]) -> Array[Int] {
  rp_stable_order_tol(stat)
}

/// Build one top table for a direction.
fn rp_build_top(
  stat : Array[Double],
  pv : Array[Double],
  pf : Array[Double],
  avefc : Array[Double],
  logged : Bool,
  logbase : Double,
  cutoff : Double?,
  num_gene : Int?,
  use_pfp : Bool,
) -> RankProdTopTable {
  let order = rp_top_order(stat)
  let selected : Array[Int] = []
  match cutoff {
    Some(thr) => {
      let mut num_top = 0
      for t in 0..<order.length() {
        let g = order[t]
        let m = if use_pfp {
          pf[g]
        } else {
          pv[g]
        }
        if m < thr {
          num_top = t + 1
        }
      }
      for t in 0..<num_top {
        selected.push(order[t])
      }
    }
    None => {
      match num_gene {
        Some(k) => {
          let kk = if k < 0 {
            0
          } else if k > order.length() {
            order.length()
          } else {
            k
          }
          for t in 0..<kk {
            selected.push(order[t])
          }
        }
        None => ()
      }
    }
  }
  {
    indices: selected.copy(),
    statistic: selected.map(fn(g) { rp_signif4(stat[g]) }),
    fold_change: selected.map(fn(g) {
      let a = avefc[g]
      let raw = if logged {
        @math.pow(logbase, a)
      } else {
        a
      }
      rp_signif4(raw)
    }),
    pfp: selected.map(fn(g) { rp_signif4(pf[g]) }),
    p_value: selected.map(fn(g) { rp_signif4(pv[g]) }),
  }
}

///|
/// Significant-feature selection mirroring package ``topGene``. Exactly one
/// of ``cutoff`` (with ``method`` = "pfp" or "pval") or ``num_gene`` must be
/// given. Returns ``(table_down, table_up)`` for class1 < class2 and
/// class1 > class2. Fold changes are ``logbase^AveFC`` for logged input and
/// the raw ratio otherwise; table columns are rounded to four significant
/// digits as in the package.
pub fn rankprod_top_genes(
  result~ : RankProdResult,
  cutoff~ : Double? = None,
  num_gene~ : Int? = None,
  method~ : String = "pfp",
  logbase~ : Double = 2.0,
) -> (RankProdTopTable, RankProdTopTable) raise RankProdError {
  if cutoff is None && num_gene is None {
    raise RankProdError(
      "rankprod_top_genes: provide either cutoff or num_gene",
    )
  }
  if !(cutoff is None) && !(num_gene is None) {
    raise RankProdError(
      "rankprod_top_genes: cutoff and num_gene are mutually exclusive",
    )
  }
  let use_pfp = if method == "pfp" {
    true
  } else if method == "pval" {
    false
  } else {
    raise RankProdError(
      "rankprod_top_genes: method must be \"pfp\" or \"pval\"",
    )
  }
  let down = rp_build_top(
    result.statistic_down,
    result.p_value_down,
    result.pfp_down,
    result.average_fold_change,
    result.logged,
    logbase,
    cutoff,
    num_gene,
    use_pfp,
  )
  let up = rp_build_top(
    result.statistic_up,
    result.p_value_up,
    result.pfp_up,
    result.average_fold_change,
    result.logged,
    logbase,
    cutoff,
    num_gene,
    use_pfp,
  )
  (down, up)
}
