///|
/// Bioconductor limma gene-set enrichment tests: competitive tests allowing
/// for inter-gene correlation.
///
/// Ports (limma 3.68.5, pure R semantics):
/// - ``rankSumTestWithCorrelation``: Wilcoxon/Mann-Whitney rank-sum test with
///   a variance inflation factor derived from a common correlation among the
///   test statistics, plus average-rank tie correction and continuity
///   correction.
/// - ``cameraPR.default``: pre-ranked competitive gene-set test, either from
///   ranks (rank-sum) or from a two-sample t statistic with VIF.
/// - ``camera.default``: competitive test from a fitted linear model. The
///   moderated t statistics are computed with limma's scaled-F prior
///   (``squeezeVar``/``fitFDist``) and transformed by the Hill approximation
///   (``zscoreT``); inter-gene correlation can be estimated from the
///   standardized residuals.
///
/// All p-values use full-precision regularized incomplete beta (Student t
/// CDF) and Chebyshev erfc (standard normal CDF) implementations, matching R
/// to ~1e-12.
pub suberror GeneSetTestError {
  GeneSetTestError(String)
}

///|
/// A named set of genes identified by their zero-based positions into the
/// full statistic/expression vector.
pub struct CameraGeneSet {
  name : String
  indices : Array[Int]
}

///|
/// Create a gene set. Indices are copied; range/uniqueness is validated when
/// the set is used.
pub fn CameraGeneSet::create(
  name : String,
  indices : Array[Int],
) -> CameraGeneSet raise GeneSetTestError {
  if indices.length() == 0 {
    raise GeneSetTestError("gene set '" + name + "' is empty")
  }
  CameraGeneSet::{ name, indices: indices.copy(), }
}

///|
pub fn CameraGeneSet::name(self : CameraGeneSet) -> String {
  self.name
}

///|
/// Return a copy of the gene indices.
pub fn CameraGeneSet::indices(self : CameraGeneSet) -> Array[Int] {
  self.indices.copy()
}

///|
/// Result of a competitive gene-set test for one set.
pub struct CameraResult {
  name : String
  n_genes : Int
  correlation : Double?
  direction : String
  p_value : Double
  fdr : Double?
}

///|
pub fn CameraResult::name(self : CameraResult) -> String {
  self.name
}

///|
pub fn CameraResult::n_genes(self : CameraResult) -> Int {
  self.n_genes
}

///|
/// Per-set inter-gene correlation; ``None`` when a fixed correlation was used
/// (limma omits the Correlation column) or for one-gene sets under
/// estimation (limma records NA).
pub fn CameraResult::correlation(self : CameraResult) -> Double? {
  self.correlation
}

///|
/// ``"Up"`` or ``"Down"`` for directional tests; empty for non-directional.
pub fn CameraResult::direction(self : CameraResult) -> String {
  self.direction
}

///|
pub fn CameraResult::p_value(self : CameraResult) -> Double {
  self.p_value
}

///|
/// BH-adjusted p-value; ``None`` when only one set was tested.
pub fn CameraResult::fdr(self : CameraResult) -> Double? {
  self.fdr
}

// ---------------------------------------------------------------------------
// Numerical primitives
// ---------------------------------------------------------------------------

///|
/// Lanczos approximation for ln(gamma(x)) (Numerical Recipes form with the
/// 1/x normalization; the package-level ``lngamma`` in limma.mbt omits that
/// normalization and is biased by -ln(x), so it cannot be reused here).
fn gst_lngamma(x_in : Double) -> Double {
  let cof = [
    76.18009172947146, -86.50532032941677, 24.01409824083091, -1.231739572450155,
    0.1208650973866179e-2, -0.5395239384953e-5,
  ]
  let mut y = x_in
  let mut tmp = x_in + 5.5
  tmp -= (x_in + 0.5) * @math.ln(tmp)
  let mut ser = 1.000000000190015
  for c in cof {
    y += 1.0
    ser += c / y
  }
  -tmp + @math.ln(2.5066282746310005 * ser / x_in)
}

///|
/// Incomplete beta continued fraction (Numerical Recipes betacf), high
/// precision; replaces the looser beta routine elsewhere in the package.
fn gst_betacf(a : Double, b : Double, x : Double) -> Double {
  let fpmin = 1.0e-300
  let qab = a + b
  let qap = a + 1.0
  let qam = a - 1.0
  let mut c = 1.0
  let mut d = 1.0 - qab * x / qap
  if d.abs() < fpmin {
    d = fpmin
  }
  d = 1.0 / d
  let mut h = d
  let mut m = 1
  while m <= 300 {
    let m2 = 2 * m
    let mut aa = m.to_double() *
      (b - m.to_double()) *
      x /
      ((qam + m2.to_double()) * (a + m2.to_double()))
    d = 1.0 + aa * d
    if d.abs() < fpmin {
      d = fpmin
    }
    c = 1.0 + aa / c
    if c.abs() < fpmin {
      c = fpmin
    }
    d = 1.0 / d
    h *= d * c
    aa = -(a + m.to_double()) *
      (qab + m.to_double()) *
      x /
      ((a + m2.to_double()) * (qap + m2.to_double()))
    d = 1.0 + aa * d
    if d.abs() < fpmin {
      d = fpmin
    }
    c = 1.0 + aa / c
    if c.abs() < fpmin {
      c = fpmin
    }
    d = 1.0 / d
    let del = d * c
    h *= del
    if (del - 1.0).abs() < 1.0e-14 {
      m = 301
    } else {
      m += 1
    }
  }
  h
}

///|
/// Regularized incomplete beta function I_x(a, b).
fn gst_pbeta(x : Double, a : Double, b : Double) -> Double {
  if x <= 0.0 {
    return 0.0
  }
  if x >= 1.0 {
    return 1.0
  }
  let bt = @math.exp(
    gst_lngamma(a + b) -
    gst_lngamma(a) -
    gst_lngamma(b) +
    a * @math.ln(x) +
    b * @math.ln(1.0 - x),
  )
  if x < (a + 1.0) / (a + b + 2.0) {
    bt * gst_betacf(a, b, x) / a
  } else {
    1.0 - bt * gst_betacf(b, a, 1.0 - x) / b
  }
}

///|
/// Lower-tail Student t CDF P(T <= t).
fn gst_pt_lower(t : Double, df : Double) -> Double {
  let x = df / (df + t * t)
  if t > 0.0 {
    1.0 - 0.5 * gst_pbeta(x, df / 2.0, 0.5)
  } else {
    0.5 * gst_pbeta(x, df / 2.0, 0.5)
  }
}

///|
/// Coefficients for Numerical Recipes erfc Chebyshev approximation.
let gst_erfc_cof : Array[Double] = [
  -1.3026537197817094, 6.4196979235649026e-1, 1.9476473204185836e-2, -9.561514786808631e-3,
  -9.46595344482036e-4, 3.66839497852761e-4, 4.2523324806907e-5, -2.0278578112534e-5,
  -1.624290004647e-6, 1.303655835580e-6, 1.5626441722e-8, -8.5238095915e-8, 6.529054439e-9,
  5.059343495e-9, -9.91364156e-10, -2.27365122e-10, 9.6467911e-11, 2.394038e-12,
  -6.886027e-12, 8.94487e-13, 3.13092e-13, -1.12708e-13, 3.81e-14, 7.106e-15, -1.523e-15,
  -9.4e-17, 1.21e-16, -2.8e-18,
]

///|
/// erfc(z) for z >= 0 (Chebyshev, ~1e-16).
fn gst_erfc_cheb(z : Double) -> Double {
  let mut d = 0.0
  let mut dd = 0.0
  let t = 2.0 / (2.0 + z)
  let ty = 4.0 * t - 2.0
  let mut j = gst_erfc_cof.length() - 1
  while j >= 1 {
    let tmp = d
    d = ty * d - dd + gst_erfc_cof[j]
    dd = tmp
    j -= 1
  }
  t * @math.exp(-z * z + 0.5 * (gst_erfc_cof[0] + ty * d) - dd)
}

///|
/// Complementary error function for any real argument.
fn gst_erfc(x : Double) -> Double {
  if x >= 0.0 {
    gst_erfc_cheb(x)
  } else {
    2.0 - gst_erfc_cheb(-x)
  }
}

///|
/// Standard normal CDF.
fn gst_pnorm(x : Double) -> Double {
  0.5 * gst_erfc(-x / 1.4142135623730951)
}

///|
/// Digamma (log-gamma derivative) via recurrence + asymptotic series.
fn gst_digamma(x_in : Double) -> Double {
  let mut x = x_in
  let mut result = 0.0
  while x < 10.0 {
    result -= 1.0 / x
    x += 1.0
  }
  let inv = 1.0 / x
  let inv2 = inv * inv
  result += @math.ln(x) -
    0.5 * inv -
    inv2 *
    (
      1.0 / 12.0 -
      inv2 *
      (
        1.0 / 120.0 -
        inv2 *
        (
          1.0 / 252.0 -
          inv2 * (1.0 / 240.0 - inv2 * (5.0 / 660.0 - inv2 * 691.0 / 32760.0))
        )
      )
    )
  result
}

///|
/// Trigamma (polygamma order 1) via recurrence + asymptotic series.
fn gst_trigamma(x_in : Double) -> Double {
  let mut x = x_in
  let mut result = 0.0
  while x < 10.0 {
    result += 1.0 / (x * x)
    x += 1.0
  }
  let inv = 1.0 / x
  let inv2 = inv * inv
  result += inv +
    0.5 * inv2 +
    inv2 *
    inv *
    (
      1.0 / 6.0 +
      inv2 *
      (
        -1.0 / 30.0 +
        inv2 *
        (
          1.0 / 42.0 -
          inv2 *
          (
            1.0 / 30.0 +
            inv2 * (5.0 / 66.0 - inv2 * (691.0 / 2730.0 - inv2 * 7.0 / 6.0))
          )
        )
      )
    )
  result
}

///|
/// Polygamma order 2 via recurrence + asymptotic series.
fn gst_polygamma2(x_in : Double) -> Double {
  let mut x = x_in
  let mut result = 0.0
  // psi2(x+1) = psi2(x) + 2/x^3, so shifting up accumulates -2/x^3.
  while x < 10.0 {
    result -= 2.0 / (x * x * x)
    x += 1.0
  }
  let inv = 1.0 / x
  let inv2 = inv * inv
  result += -inv2 -
    inv2 * inv -
    0.5 * inv2 * inv2 +
    inv2 *
    inv2 *
    inv2 *
    (
      1.0 / 6.0 +
      inv2 * (-1.0 / 6.0 + inv2 * 0.3 - inv2 * (5.0 / 42.0 - inv2 * 5.0 / 30.0))
    )
  result
}

///|
/// Solve trigamma(y) = x for y (limma trigammaInverse, scalar path).
fn gst_trigamma_inverse(x : Double) -> Double {
  if x > 1.0e7 {
    return 1.0 / x.sqrt()
  }
  if x < 1.0e-6 {
    return 1.0 / x
  }
  let mut y = 0.5 + 1.0 / x
  let mut iter = 0
  while iter < 50 {
    let tri = gst_trigamma(y)
    let dif = tri * (1.0 - tri / x) / gst_polygamma2(y)
    y += dif
    iter += 1
    if (-dif / y).abs() < 1.0e-8 {
      iter = 50
    }
  }
  y
}

///|
/// limma fitFDist for a common numerator df without covariate: moment
/// estimation of a scaled F-distribution. Returns (scale s20, df2 prior).
fn gst_fit_f_dist(x : Array[Double], df1 : Double) -> (Double, Double) {
  let n = x.length()
  // Median (R median with average of middle pair).
  let sorted = x.copy()
  sorted.sort()
  let mid = n / 2
  let m_median = if n % 2 == 1 {
    sorted[mid]
  } else {
    0.5 * (sorted[mid - 1] + sorted[mid])
  }
  let m = if m_median == 0.0 { 1.0 } else { m_median }
  let floor = 1.0e-5 * m
  let mut sum_e = 0.0
  let mut sum_x = 0.0
  let zs = Array::make(n, 0.0)
  let mut i = 0
  while i < n {
    let xi = if x[i] < floor { floor } else { x[i] }
    sum_x += xi
    let e = @math.ln(xi) + @math.ln(gst_digamma(df1 / 2.0))
    zs[i] = e
    sum_e += e
    i += 1
  }
  let emean = sum_e / n.to_double()
  let mut evar = 0.0
  i = 0
  while i < n {
    let d = zs[i] - emean
    evar += d * d
    i += 1
  }
  evar = evar / (n - 1).to_double() - gst_trigamma(df1 / 2.0)
  if evar > 0.0 {
    let df2 = 2.0 * gst_trigamma_inverse(evar)
    let s20 = @math.exp(emean - @math.ln(gst_digamma(df2 / 2.0)))
    (s20, df2)
  } else {
    // df.prior = Inf; pooled variance scale
    (sum_x / n.to_double(), 1.0e300)
  }
}

///|
/// limma squeezeVar for common residual df without covariate/robust fitting.
/// Returns (df.prior, var.prior, var.post[]).
fn gst_squeeze_var(
  variance : Array[Double],
  df : Double,
) -> (Double, Double, Array[Double]) {
  let n = variance.length()
  if n < 3 {
    return (0.0, 0.0, variance.copy())
  }
  let (var_prior, df_prior) = gst_fit_f_dist(variance, df)
  let var_post = Array::make(n, 0.0)
  let mut i = 0
  while i < n {
    if df_prior > 1.0e100 {
      var_post[i] = var_prior
    } else {
      var_post[i] = (df * variance[i] + df_prior * var_prior) / (df + df_prior)
    }
    i += 1
  }
  (df_prior, var_prior, var_post)
}

///|
/// Hill's (1970) approximation of the inverse normal CDF applied to a
/// Student-t variate (limma zscoreT, method="hill").
fn gst_zscore_hill(x : Double, df : Double) -> Double {
  let a = df - 0.5
  let b = 48.0 * a * a
  let z0 = a * @math.ln(1.0 + x / df * x)
  let inner = (
      (((-0.4 * z0 - 3.3) * z0 - 24.0) * z0 - 85.5) /
      (0.8 * z0 * z0 + 100.0 + b) +
      z0 +
      3.0
    ) /
    b +
    1.0
  let z = inner * z0.sqrt()
  if x < 0.0 {
    -z
  } else {
    z
  }
}

// ---------------------------------------------------------------------------
// Ranks and validation helpers
// ---------------------------------------------------------------------------

///|
/// Average ranks (R rank default, ties.method="average"), ranks start at 1.
fn gst_average_ranks(values : Array[Double]) -> Array[Double] {
  let n = values.length()
  let order = Array::make(n, 0)
  let mut i = 0
  while i < n {
    order[i] = i
    i += 1
  }
  // Stable insertion sort on values.
  let mut a = 1
  while a < n {
    let key = order[a]
    let mut b = a - 1
    while b >= 0 && values[order[b]] > values[key] {
      order[b + 1] = order[b]
      b -= 1
    }
    order[b + 1] = key
    a += 1
  }
  let ranks = Array::make(n, 0.0)
  let mut start = 0
  while start < n {
    let mut finish = start + 1
    while finish < n && values[order[finish]] == values[order[start]] {
      finish += 1
    }
    let avg = (start + 1 + finish).to_double() / 2.0
    let mut k = start
    while k < finish {
      ranks[order[k]] = avg
      k += 1
    }
    start = finish
  }
  ranks
}

///|
fn gst_check_finite(
  values : Array[Double],
  what : String,
) -> Unit raise GeneSetTestError {
  for value in values {
    if value.is_nan() || value.abs() > 1.7976931348623157e308 {
      raise GeneSetTestError(what + " contains NA/NaN or infinite values")
    }
  }
}

///|
/// Validate gene-set indices against a universe of ``n_genes`` positions.
fn gst_validated_indices(
  set : CameraGeneSet,
  n_genes : Int,
) -> Array[Int] raise GeneSetTestError {
  let indices = set.indices
  let m = indices.length()
  if m == 0 {
    raise GeneSetTestError("gene set '" + set.name + "' is empty")
  }
  if m >= n_genes {
    raise GeneSetTestError(
      "gene set '" +
      set.name +
      "' must contain fewer genes than the universe (" +
      m.to_string() +
      " >= " +
      n_genes.to_string() +
      ")",
    )
  }
  for i in indices {
    if i < 0 || i >= n_genes {
      raise GeneSetTestError(
        "gene set '" +
        set.name +
        "' contains index " +
        i.to_string() +
        " outside [0, " +
        (n_genes - 1).to_string() +
        "]",
      )
    }
  }
  for i in 0..<m {
    let mut j = i + 1
    while j < m {
      if indices[i] == indices[j] {
        raise GeneSetTestError(
          "gene set '" +
          set.name +
          "' contains duplicated index " +
          indices[i].to_string(),
        )
      }
      j += 1
    }
  }
  indices
}

///|
/// BH adjustment (R p.adjust method="BH"); thin wrapper over the package
/// implementation.
fn gst_bh(p_values : Array[Double]) -> Array[Double] {
  adjust_p_values(p_values, "BH")
}

///|
/// Stable sort of set positions by p-value (R order, ties keep input order).
fn gst_order_by_p(p_values : Array[Double]) -> Array[Int] {
  let n = p_values.length()
  let order = Array::make(n, 0)
  let mut i = 0
  while i < n {
    order[i] = i
    i += 1
  }
  let mut a = 1
  while a < n {
    let key = order[a]
    let mut b = a - 1
    while b >= 0 && p_values[order[b]] > p_values[key] {
      order[b + 1] = order[b]
      b -= 1
    }
    order[b + 1] = key
    a += 1
  }
  order
}

// ---------------------------------------------------------------------------
// rankSumTestWithCorrelation
// ---------------------------------------------------------------------------

///|
/// Wilcoxon rank-sum test allowing for a common correlation among the
/// statistics.
///
/// Mirrors limma ``rankSumTestWithCorrelation``. Returns ``(down, up)``:
/// the p-values that the set is less than / greater than the rest. Note
/// that, as in R, the lower output tail uses the upper normal/t tail and
/// vice versa (R ranks run opposite to Mann-Whitney ranks).
///
/// - ``index``: zero-based positions of genes in the set.
/// - ``statistics``: gene-level statistics (ranked if called from ranks).
/// - ``correlation``: common pairwise correlation of statistics.
/// - ``df``: degrees of freedom; ``None`` denotes infinite (standard normal).
pub fn rank_sum_test_with_correlation(
  index : Array[Int],
  statistics : Array[Double],
  correlation? : Double = 0.0,
  df? : Double? = Option::None,
) -> (Double, Double) raise GeneSetTestError {
  let n = statistics.length()
  let n1 = index.length()
  if n1 == 0 {
    raise GeneSetTestError("index is empty")
  }
  if n1 >= n {
    raise GeneSetTestError("index must contain fewer genes than statistics")
  }
  gst_check_finite(statistics, "statistics")
  if correlation.is_nan() || correlation.abs() >= 1.0 {
    raise GeneSetTestError("correlation must be strictly between -1 and 1")
  }
  let n2 = n - n1
  let ranks = gst_average_ranks(statistics)
  let mut rank_sum = 0.0
  for i in index {
    if i < 0 || i >= n {
      raise GeneSetTestError("index " + i.to_string() + " out of range")
    }
    rank_sum += ranks[i]
  }
  let n1d = n1.to_double()
  let n2d = n2.to_double()
  let nd = n.to_double()
  let u = n1d * n2d + n1d * (n1d + 1.0) / 2.0 - rank_sum
  let mu = n1d * n2d / 2.0
  let mut sigma2 = if correlation == 0.0 || n1 == 1 {
    n1d * n2d * (nd + 1.0) / 12.0
  } else {
    let c1 = n1d * n2d
    let c2 = n1d * n2d * (n2d - 1.0)
    let c3 = n1d * (n1d - 1.0) * n2d * (n2d - 1.0)
    let c4 = n1d * (n1d - 1.0) * n2d
    (
      @math.asin(1.0) * c1 +
      @math.asin(0.5) * c2 +
      @math.asin(correlation / 2.0) * c3 +
      @math.asin((correlation + 1.0) / 2.0) * c4
    ) /
    2.0 /
    3.141592653589793
  }

  // Tie correction via NTIES (R rankSumTestWithCorrelation).
  let sorted = statistics.copy()
  sorted.sort()
  let mut tie_sum = 0.0
  let mut start = 0
  while start < n {
    let mut finish = start + 1
    while finish < n && sorted[finish] == sorted[start] {
      finish += 1
    }
    let t = (finish - start).to_double()
    if t > 1.0 {
      tie_sum += t * (t + 1.0) * (t - 1.0)
    }
    start = finish
  }
  if tie_sum > 0.0 {
    let adjustment = tie_sum / (nd * (nd + 1.0) * (nd - 1.0))
    sigma2 *= 1.0 - adjustment
  }

  let sigma = sigma2.sqrt()
  let z_lower = (u + 0.5 - mu) / sigma
  let z_upper = (u - 0.5 - mu) / sigma

  let down = match df {
    Option::None => 1.0 - gst_pnorm(z_upper)
    Some(d) => 1.0 - gst_pt_lower(z_upper, d)
  }
  let up = match df {
    Option::None => gst_pnorm(z_lower)
    Some(d) => gst_pt_lower(z_lower, d)
  }
  (down, up)
}

// ---------------------------------------------------------------------------
// cameraPR
// ---------------------------------------------------------------------------

///|
/// Pre-ranked competitive gene-set test (limma ``cameraPR.default``).
///
/// - ``statistic``: per-gene statistics, at least 3 genes required.
/// - ``sets``: named gene sets (zero-based indices).
/// - ``fixed_inter_gene_cor``: common VIF correlation used for every set.
/// - ``per_set_inter_gene_cor``: when given, one correlation per set;
///   overrides the fixed value and adds a Correlation column.
/// - ``use_ranks``: rank-sum test (default: two-sample t formulation).
/// - ``directional``: when false, statistics are converted to absolute
///   values and ranks are forced (as in R).
/// - ``sort``: sort results by p-value.
pub fn camera_pr(
  statistic : Array[Double],
  sets : Array[CameraGeneSet],
  fixed_inter_gene_cor? : Double = 0.01,
  per_set_inter_gene_cor? : Array[Double]? = Option::None,
  use_ranks? : Bool = false,
  directional? : Bool = true,
  sort? : Bool = true,
) -> Array[CameraResult] raise GeneSetTestError {
  let g = statistic.length()
  if g < 3 {
    raise GeneSetTestError("Too few genes in dataset: need at least 3")
  }
  gst_check_finite(statistic, "statistic")
  let nsets = sets.length()
  if nsets == 0 {
    raise GeneSetTestError("at least one gene set is required")
  }

  let (cors, fixed_cor) = match per_set_inter_gene_cor {
    Some(c) => {
      if c.length() != nsets {
        raise GeneSetTestError(
          "Length of inter_gene_cor doesn't match number of sets",
        )
      }
      (c.copy(), false)
    }
    Option::None => {
      if fixed_inter_gene_cor.is_nan() || fixed_inter_gene_cor.abs() >= 1.0 {
        raise GeneSetTestError(
          "inter_gene_cor must be strictly between -1 and 1",
        )
      }
      (Array::make(nsets, fixed_inter_gene_cor), true)
    }
  }
  for c in cors {
    if c.is_nan() || c.abs() >= 1.0 {
      raise GeneSetTestError(
        "NA inter.gene.cor not allowed and values must be strictly between -1 and 1",
      )
    }
  }

  // Copy statistic because non-directional tests may replace it by abs.
  let stats = statistic.copy()
  let ranks_mode = use_ranks || !directional
  if !directional {
    let mut min_stat = stats[0]
    for v in stats {
      if v < min_stat {
        min_stat = v
      }
    }
    if min_stat < 0.0 {
      let mut i = 0
      while i < g {
        stats[i] = stats[i].abs()
        i += 1
      }
    }
  }

  let df_camera = if ranks_mode {
    Option::None
  } else {
    Some((g - 2).to_double())
  }

  let mut sum = 0.0
  for v in stats {
    sum += v
  }
  let mean_stat = sum / g.to_double()
  let mut ss = 0.0
  for v in stats {
    let d = v - mean_stat
    ss += d * d
  }
  let var_stat = ss / (g - 1).to_double()

  let n_genes = Array::make(nsets, 0)
  let down = Array::make(nsets, 0.0)
  let up = Array::make(nsets, 0.0)
  let gd = g.to_double()

  let mut i = 0
  while i < nsets {
    let iset = gst_validated_indices(sets[i], g)
    let m = iset.length()
    n_genes[i] = m
    if ranks_mode {
      let (d, u) = rank_sum_test_with_correlation(
        iset,
        stats,
        correlation=cors[i],
        df=Option::None,
      )
      down[i] = d
      up[i] = u
    } else {
      let md = m.to_double()
      let m2 = g - m
      let m2d = m2.to_double()
      let vif = 1.0 + (md - 1.0) * cors[i]
      let mut mean_in = 0.0
      for gi in iset {
        mean_in += stats[gi]
      }
      mean_in /= md
      let delta = gd / m2d * (mean_in - mean_stat)
      let var_pooled = ((gd - 1.0) * var_stat - delta * delta * md * m2d / gd) /
        (gd - 2.0)
      let t = delta / (var_pooled * (vif / md + 1.0 / m2d)).sqrt()
      down[i] = gst_pt_lower(t, df_camera.unwrap())
      up[i] = 1.0 - down[i]
    }
    i += 1
  }

  gst_assemble_results(
    sets, n_genes, down, up, cors, fixed_cor, directional, sort,
  )
}

///|
/// Shared camera/cameraPR result assembly: direction, FDR, sorting.
fn gst_assemble_results(
  sets : Array[CameraGeneSet],
  n_genes : Array[Int],
  down : Array[Double],
  up : Array[Double],
  cors : Array[Double],
  fixed_cor : Bool,
  directional : Bool,
  sort : Bool,
) -> Array[CameraResult] {
  let nsets = sets.length()
  let p_values = Array::make(nsets, 0.0)
  let directions = Array::make(nsets, "")
  let mut i = 0
  while i < nsets {
    if directional {
      let two_sided = 2.0 * (if down[i] < up[i] { down[i] } else { up[i] })
      p_values[i] = two_sided
      directions[i] = if down[i] < up[i] { "Down" } else { "Up" }
    } else {
      p_values[i] = up[i]
    }
    i += 1
  }
  let fdrs = if nsets > 1 { gst_bh(p_values) } else { Array(capacity=0) }
  let order = if sort && nsets > 1 {
    gst_order_by_p(p_values)
  } else {
    let o = Array::make(nsets, 0)
    let mut k = 0
    while k < nsets {
      o[k] = k
      k += 1
    }
    o
  }
  let results = Array::make(nsets, CameraResult::{
    name: "",
    n_genes: 0,
    correlation: Option::None,
    direction: "",
    p_value: 0.0,
    fdr: Option::None,
  })
  let mut pos = 0
  for k in order {
    let correlation = if fixed_cor { Option::None } else { Some(cors[k]) }
    let fdr = if nsets > 1 { Some(fdrs[k]) } else { Option::None }
    results[pos] = CameraResult::{
      name: sets[k].name,
      n_genes: n_genes[k],
      correlation,
      direction: directions[k],
      p_value: p_values[k],
      fdr,
    }
    pos += 1
  }
  results
}

///|
/// Orthonormal basis of the nullspace of an n-by-p design matrix via
/// modified Gram-Schmidt on [design columns, identity columns]. The first
/// ``p`` accepted vectors span the design space; the remaining ``n - p``
/// span the residual space. Stored column-major: component ``i`` of column
/// ``k`` is at ``q[k * n + i]``. Only the residual columns are returned.
fn gst_qr_nullspace(
  design : Matrix,
  n : Int,
  p : Int,
) -> Array[Double] raise GeneSetTestError {
  let q = Array::make(n * n, 0.0)
  let mut ncols = 0

  let add_vector = fn(w_in : Array[Double]) -> Bool {
    let w = w_in.copy()
    let mut j = 0
    while j < ncols {
      let mut dot = 0.0
      let mut i = 0
      while i < n {
        dot += w[i] * q[j * n + i]
        i += 1
      }
      let mut i = 0
      while i < n {
        w[i] -= dot * q[j * n + i]
        i += 1
      }
      j += 1
    }
    let mut norm = 0.0
    for v in w {
      norm += v * v
    }
    norm = norm.sqrt()
    if norm < 1.0e-12 {
      return false
    }
    let mut i = 0
    while i < n {
      q[ncols * n + i] = w[i] / norm
      i += 1
    }
    ncols += 1
    true
  }

  // Design columns first.
  let mut c = 0
  while c < p {
    let col = Array::make(n, 0.0)
    let mut i = 0
    while i < n {
      col[i] = design.get(i, c)
      i += 1
    }
    if !add_vector(col) {
      raise GeneSetTestError("design matrix is not full column rank")
    }
    c += 1
  }
  // Identity vectors complete the basis.
  let mut e = 0
  while e < n {
    let col = Array::make(n, 0.0)
    col[e] = 1.0
    add_vector(col) |> ignore
    e += 1
  }
  if ncols != n {
    raise GeneSetTestError("failed to build residual-space basis")
  }

  let dfree = n - p
  let q2 = Array::make(dfree * n, 0.0)
  let mut k = 0
  while k < dfree {
    let mut i = 0
    while i < n {
      q2[k * n + i] = q[(p + k) * n + i]
      i += 1
    }
    k += 1
  }
  q2
}

// ---------------------------------------------------------------------------
// camera.default
// ---------------------------------------------------------------------------

///|
/// Competitive gene-set test directly from an expression matrix (limma
/// ``camera.default``).
///
/// A linear model is fit for every gene; moderated t statistics for the
/// selected contrast use limma's scaled-F variance squeezing, and are
/// transformed to standard-normal scale by the Hill approximation for the
/// (default) parametric test. Inter-gene correlation is either fixed or
/// estimated from the standardized residuals of each set.
///
/// - ``y``: genes (rows) by samples (cols) expression matrix.
/// - ``sets``: named gene sets (zero-based row indices).
/// - ``design``: samples by coefficients design matrix, full rank.
/// - ``contrast``: zero-based design column tested; defaults to the last.
/// - ``inter_gene_cor``: fixed correlation; ``None`` estimates per set.
/// - ``allow_negative_cor``: keep negative rank correlations (default false,
///   matching limma).
pub fn camera_default(
  y : Matrix,
  sets : Array[CameraGeneSet],
  design : Matrix,
  contrast? : Int? = Option::None,
  inter_gene_cor? : Double? = Some(0.01),
  allow_negative_cor? : Bool = false,
  use_ranks? : Bool = false,
  directional? : Bool = true,
  sort? : Bool = true,
) -> Array[CameraResult] raise GeneSetTestError {
  let g = y.rows
  let n = y.cols
  let p = design.cols
  if design.rows != n {
    raise GeneSetTestError("design must have one row per sample")
  }
  if p < 1 || n <= p {
    raise GeneSetTestError(
      "design must be full column rank with residual degrees of freedom",
    )
  }
  if g < 3 {
    raise GeneSetTestError("Too few genes in dataset: need at least 3")
  }
  let nsets = sets.length()
  if nsets == 0 {
    raise GeneSetTestError("at least one gene set is required")
  }
  let j_contrast = match contrast {
    Option::None => p - 1
    Some(j) => {
      if j < 0 || j >= p {
        raise GeneSetTestError(
          "contrast column " + j.to_string() + " out of range",
        )
      }
      j
    }
  }
  let fixed_cor_value = match inter_gene_cor {
    Option::None => Option::None
    Some(c) => {
      if c.is_nan() || c.abs() >= 1.0 {
        raise GeneSetTestError(
          "inter_gene_cor must be strictly between -1 and 1",
        )
      }
      Some(c)
    }
  }

  // Reorder design: contrast column last.
  let reordered_data = Array::make(n * p, 0.0)
  let col_order = Array::make(p, 0)
  let mut pos = 0
  let mut k = 0
  while k < p {
    if k != j_contrast {
      col_order[pos] = k
      pos += 1
    }
    k += 1
  }
  col_order[p - 1] = j_contrast
  let mut r = 0
  while r < n {
    let mut cc = 0
    while cc < p {
      reordered_data[r * p + cc] = design.get(r, col_order[cc])
      cc += 1
    }
    r += 1
  }
  let reordered = Matrix::new(n, p, reordered_data)
  let fit = lm_fit(y, reordered)
  let df_residual = fit.df_residual.to_double()

  let sigma2 = Array::make(g, 0.0)
  let mut gi = 0
  while gi < g {
    let s = fit.sigma[gi]
    sigma2[gi] = s * s
    gi += 1
  }
  let (df_prior, _, var_post) = gst_squeeze_var(sigma2, df_residual)
  let mut df_total = df_residual + df_prior
  let cap = g.to_double() * df_residual
  if df_total > cap {
    df_total = cap
  }
  // R rounding: signif(df.total, 10).
  df_total = (df_total * 1.0e10).floor() / 1.0e10

  // Moderated t for the contrast (last coefficient).
  // lm_fit stores the full SE (sigma * unscaled factor); the unscaled
  // factor is constant across genes for an unweighted design, so recover
  // it from the first gene.
  let mut unse = -1.0
  gi = 0
  while gi < g {
    if fit.sigma[gi] > 0.0 {
      unse = fit.stdev_unscaled.get(gi, p - 1) / fit.sigma[gi]
      gi = g
    }
    gi += 1
  }
  let modt = Array::make(g, 0.0)
  gi = 0
  while gi < g {
    let beta = fit.coefficients.get(gi, p - 1)
    let denom = if fit.sigma[gi] == 0.0 {
      (1.0e-8 * var_post[gi]).sqrt()
    } else {
      var_post[gi].sqrt() * unse
    }
    modt[gi] = beta / denom
    gi += 1
  }

  let ranks_mode = use_ranks || !directional
  let stats = modt.copy()
  if !directional {
    let mut ii = 0
    while ii < g {
      stats[ii] = stats[ii].abs()
      ii += 1
    }
  }

  // Statistics used for the parametric test: Hill normal-equivalents.
  let zstats = Array::make(g, 0.0)
  gi = 0
  while gi < g {
    zstats[gi] = gst_zscore_hill(modt[gi], df_total)
    gi += 1
  }
  let test_stats = if ranks_mode { stats } else { zstats }

  let mut sum = 0.0
  for v in test_stats {
    sum += v
  }
  let mean_stat = sum / g.to_double()
  let mut ssq = 0.0
  for v in test_stats {
    let d = v - mean_stat
    ssq += d * d
  }
  let var_stat = ssq / (g - 1).to_double()

  // Orthonormal basis of the residual space (R uses qr.qty effects; any
  // orthonormal basis spanning the nullspace gives identical VIF values).
  let dfree = n - p
  let residuals_needed = fixed_cor_value is None
  let q2 = if residuals_needed {
    gst_qr_nullspace(reordered, n, p)
  } else {
    Array::make(0, 0.0)
  }
  // Standardized residual coordinates U[gene, df] (limma camera).
  let resid_u = Array::make(g * dfree, 0.0)
  if residuals_needed {
    gi = 0
    while gi < g {
      let denom = sigma2[gi].max(1.0e-8).sqrt()
      let residual = Array::make(n, 0.0)
      let mut jj = 0
      while jj < n {
        residual[jj] = (y.get(gi, jj) - fit.fit.get(gi, jj)) / denom
        jj += 1
      }
      let mut kk = 0
      while kk < dfree {
        let mut coord = 0.0
        let mut jj = 0
        while jj < n {
          coord += q2[kk * n + jj] * residual[jj]
          jj += 1
        }
        resid_u[gi * dfree + kk] = coord
        kk += 1
      }
      gi += 1
    }
  }

  let n_genes = Array::make(nsets, 0)
  let down = Array::make(nsets, 0.0)
  let up = Array::make(nsets, 0.0)
  let cors = Array::make(nsets, 0.0)
  let has_cor = Array::make(nsets, false)
  let gd = g.to_double()
  // limma: fixed correlation + ranks -> df = Inf; fixed correlation +
  // parametric -> G-2; estimated correlation -> min(df.residual, G-2)
  // for both ranks and parametric paths.
  let gm2 = (g - 2).to_double()
  let df_camera : Double? = match fixed_cor_value {
    Some(_) => if ranks_mode { Option::None } else { Some(gm2) }
    Option::None => Some(if df_residual < gm2 { df_residual } else { gm2 })
  }

  let mut i = 0
  while i < nsets {
    let iset = gst_validated_indices(sets[i], g)
    let m = iset.length()
    n_genes[i] = m

    // Resolve correlation / VIF for this set.
    let mut correlation : Double? = Option::None
    match fixed_cor_value {
      Option::None => {
        let md = m.to_double()
        let mut col_sum_sq = 0.0
        let mut jj = 0
        while jj < dfree {
          let mut cm = 0.0
          for gene in iset {
            cm += resid_u[gene * dfree + jj]
          }
          cm /= md
          col_sum_sq += cm * cm
          jj += 1
        }
        let vif = md * col_sum_sq / df_residual
        if m > 1 {
          let c_est = (vif - 1.0) / (md - 1.0)
          correlation = Some(c_est)
          cors[i] = c_est
          has_cor[i] = true
        } else {
          correlation = Option::None
        }
      }
      Some(c) => {
        correlation = Some(c)
        cors[i] = c
        has_cor[i] = true
      }
    }

    if ranks_mode {
      // One-gene sets ignore correlation inside the rank-sum routine
      // (limma n1==1 branch), so 0.0 is identical there.
      let pass_cor = match correlation {
        Option::None => 0.0
        Some(c) => if !allow_negative_cor && c < 0.0 { 0.0 } else { c }
      }
      let (d, u) = rank_sum_test_with_correlation(
        iset,
        stats,
        correlation=pass_cor,
        df=df_camera,
      )
      down[i] = d
      up[i] = u
    } else {
      let c_use = correlation.unwrap_or(0.0)
      let mut vif = 1.0 + (m.to_double() - 1.0) * c_use
      if !allow_negative_cor && vif < 1.0 {
        vif = 1.0
      }
      let m2 = g - m
      let m2d = m2.to_double()
      let mut mean_in = 0.0
      for gene in iset {
        mean_in += test_stats[gene]
      }
      mean_in /= m.to_double()
      let delta = gd / m2d * (mean_in - mean_stat)
      let var_pooled = (
          (gd - 1.0) * var_stat - delta * delta * m.to_double() * m2d / gd
        ) /
        (gd - 2.0)
      let t = delta / (var_pooled * (vif / m.to_double() + 1.0 / m2d)).sqrt()
      down[i] = gst_pt_lower(t, df_camera.unwrap())
      up[i] = 1.0 - down[i]
    }
    i += 1
  }

  // Assemble with per-set correlation metadata.
  let p_values = Array::make(nsets, 0.0)
  let directions = Array::make(nsets, "")
  let mut i = 0
  while i < nsets {
    if directional {
      p_values[i] = 2.0 * (if down[i] < up[i] { down[i] } else { up[i] })
      directions[i] = if down[i] < up[i] { "Down" } else { "Up" }
    } else {
      p_values[i] = up[i]
    }
    i += 1
  }
  let fdrs = if nsets > 1 { gst_bh(p_values) } else { Array(capacity=0) }
  let order = if sort && nsets > 1 {
    gst_order_by_p(p_values)
  } else {
    let o = Array::make(nsets, 0)
    let mut k = 0
    while k < nsets {
      o[k] = k
      k += 1
    }
    o
  }
  let results = Array::make(nsets, CameraResult::{
    name: "",
    n_genes: 0,
    correlation: Option::None,
    direction: "",
    p_value: 0.0,
    fdr: Option::None,
  })
  let mut pos = 0
  for k in order {
    let corr = if fixed_cor_value is None {
      if has_cor[k] {
        Some(cors[k])
      } else {
        Option::None
      }
    } else {
      Option::None
    }
    results[pos] = CameraResult::{
      name: sets[k].name,
      n_genes: n_genes[k],
      correlation: corr,
      direction: directions[k],
      p_value: p_values[k],
      fdr: if nsets > 1 {
        Some(fdrs[k])
      } else {
        Option::None
      },
    }
    pos += 1
  }
  results
}
