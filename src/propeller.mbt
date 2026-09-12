///|
/// propeller - testing differences in cell type proportions, ported from
/// Bioconductor speckle (Belinda Phipson) with the limma empirical Bayes
/// machinery (lmFit / contrasts.fit / eBayes, including robust variance
/// shrinkage via fitFDistRobustly) reimplemented in MoonBit.
///
/// Pipeline mirroring the R package:
/// 1. ``get_transformed_props`` counts cell types per biological replicate and
///    applies the logit (default) or arcsin-square-root transformation.
/// 2. ``propeller_ttest`` compares two groups with moderated t-tests;
///    ``propeller_anova`` compares three or more groups with moderated F-tests.
/// 3. ``propeller`` is the one-call wrapper that builds the group design,
///    dispatches t-test/ANOVA and annotates baseline cell proportions.
///
/// Helper utilities ``estimate_beta_param``, ``estimate_beta_params_from_counts``
/// and ``propeller_norm_counts`` (speckle ``normCounts``) are included as well.
///
/// Matrices follow the package [Matrix] convention: cell types are rows and
/// samples are columns.

///|
/// Error raised by the propeller functions.
pub suberror PropellerError {
  PropellerError(String)
}

///|
/// Counts, raw proportions and transformed proportions for cell types (rows)
/// across samples (columns). Mirrors the list returned by speckle
/// ``getTransformedProps``.
pub struct TransformedProps {
  counts : Matrix
  proportions : Matrix
  transformed : Matrix
  cluster_names : Array[String]
  sample_names : Array[String]
}

///|
/// Moderated t-test model, exposed for advanced use with precomputed
/// matrices (``propeller_ttest_model``). Arrays are indexed by cell type.
pub struct PropellerTModel {
  means : Matrix
  prop_ratio : Array[Double]
  statistic : Array[Double]
  p_value : Array[Double]
  fdr : Array[Double]
  df_prior : Array[Double]
  s2_prior : Double
  df_total : Array[Double]
  s2_post : Array[Double]
}

///|
/// One row of ``propeller_ttest`` output.
pub struct PropellerTRow {
  cluster : String
  means : Array[Double]
  prop_ratio : Double
  statistic : Double
  p_value : Double
  fdr : Double
}

///|
/// Moderated ANOVA model, exposed for advanced use with precomputed matrices
/// (``propeller_anova_model``). Arrays are indexed by cell type.
pub struct PropellerAnovaModel {
  means : Matrix
  statistic : Array[Double]
  p_value : Array[Double]
  fdr : Array[Double]
  df_prior : Array[Double]
  s2_prior : Double
  df_total : Array[Double]
  s2_post : Array[Double]
}

///|
/// One row of ``propeller_anova`` output.
pub struct PropellerAnovaRow {
  cluster : String
  means : Array[Double]
  statistic : Double
  p_value : Double
  fdr : Double
}

///|
/// One row of the high-level ``propeller`` wrapper output.
pub struct PropellerRow {
  cluster : String
  kind : String
  baseline_prop : Double
  means : Array[Double]
  prop_ratio : Double?
  statistic : Double
  p_value : Double
  fdr : Double
}

///|
/// Method-of-moments Beta parameters for a proportion vector.
pub struct BetaMom {
  alpha : Double
  beta : Double
}

///|
/// Beta-binomial parameters estimated from a counts matrix.
pub struct BetaCountParams {
  n : Double
  alpha : Array[Double]
  beta : Array[Double]
  pi : Array[Double]
  dispersion : Array[Double]
  variance : Array[Double]
}

// ---------------------------------------------------------------------------
// Sorted unique levels and contingency tables (R table() semantics)
// ---------------------------------------------------------------------------

///|
/// Lexicographic comparison of level names, matching R ``table()`` /
/// ``sort(character)`` in the C locale. MoonBit's native ``String`` order is
/// length-first, so the default ``sort`` cannot be used for cluster, sample
/// or group names.
fn spk_str_cmp(a : String, b : String) -> Int {
  let n = if a.length() < b.length() { a.length() } else { b.length() }
  let mut cmp = 0
  let mut i = 0
  while i < n && cmp == 0 {
    let ca = a.unsafe_get(i).to_int()
    let cb = b.unsafe_get(i).to_int()
    if ca < cb {
      cmp = -1
    } else if ca > cb {
      cmp = 1
    }
    i += 1
  }
  if cmp != 0 {
    cmp
  } else if a.length() < b.length() {
    -1
  } else if a.length() > b.length() {
    1
  } else {
    0
  }
}

///|
fn spk_sorted_unique(values : Array[String]) -> Array[String] {
  let out = []
  for v in values {
    if !out.contains(v) {
      out.push(v)
    }
  }
  out.sort_by(spk_str_cmp)
  out
}

///|
/// Binary search for ``key`` in a lexicographically sorted array; returns
/// the index or -1.
fn spk_find_sorted(names : Array[String], key : String) -> Int {
  let mut lo = 0
  let mut hi = names.length()
  while lo < hi {
    let mid = lo + (hi - lo) / 2
    let c = spk_str_cmp(names[mid], key)
    if c < 0 {
      lo = mid + 1
    } else if c > 0 {
      hi = mid
    } else {
      return mid
    }
  }
  -1
}

///|
/// Counts, raw proportions and variance-stabilising transformed proportions
/// from cell-level ``clusters`` and ``samples`` vectors. ``transform`` is
/// ``"logit"`` (default, as in speckle) or ``"asin"`` (arcsin square root).
/// Rows of every output matrix correspond to clusters, columns to samples;
/// levels are ordered lexicographically, matching R ``table()``.
pub fn get_transformed_props(
  clusters : Array[String],
  samples : Array[String],
  transform? : String = "logit",
) -> TransformedProps raise PropellerError {
  let n = clusters.length()
  if n == 0 || samples.length() == 0 {
    raise PropellerError::PropellerError(
      "clusters and samples must be non-empty",
    )
  }
  if n != samples.length() {
    raise PropellerError::PropellerError(
      "clusters and samples must have the same length",
    )
  }
  if transform != "logit" && transform != "asin" {
    raise PropellerError::PropellerError(
      "transform must be \"logit\" or \"asin\"",
    )
  }
  let cluster_names = spk_sorted_unique(clusters)
  let sample_names = spk_sorted_unique(samples)
  let k = cluster_names.length()
  let s = sample_names.length()
  // table(sample, clusters) accumulated directly in transposed orientation.
  let counts = Matrix::zeros(k, s)
  for cell in 0..<n {
    let ic = spk_find_sorted(cluster_names, clusters[cell])
    let js = spk_find_sorted(sample_names, samples[cell])
    if ic >= 0 && js >= 0 {
      counts.set(ic, js, counts.get(ic, js) + 1.0)
    }
  }
  let proportions = Matrix::zeros(k, s)
  let transformed = Matrix::zeros(k, s)
  for j in 0..<s {
    let mut lib = 0.0
    for i in 0..<k {
      lib += counts.get(i, j)
    }
    if lib <= 0.0 {
      raise PropellerError::PropellerError(
        "sample " + sample_names[j] + " contains no cells",
      )
    }
    let lib_pseudo = lib + 0.5 * k.to_double()
    for i in 0..<k {
      let p = counts.get(i, j) / lib
      proportions.set(i, j, p)
      if transform == "asin" {
        transformed.set(i, j, @math.asin(p.sqrt()))
      } else {
        let pp = (counts.get(i, j) + 0.5) / lib_pseudo
        transformed.set(i, j, @math.ln(pp / (1.0 - pp)))
      }
    }
  }
  TransformedProps::{
    counts,
    proportions,
    transformed,
    cluster_names,
    sample_names,
  }
}

// ---------------------------------------------------------------------------
// Ordinary least squares (limma lmFit core for full-rank designs)
// ---------------------------------------------------------------------------

///|
priv struct SpkLM {
  beta : Matrix // g x p coefficients
  sigma2 : Array[Double] // residual variance per response row
  cov : Matrix // (X'X)^-1
  df_resid : Double
}

///|
/// OLS fit of every row of ``y`` (g x n) on design ``x`` (n x p):
/// beta = (X'X)^-1 X'y, sigma^2 = RSS/(n-p), unscaled SE = sqrt(diag(X'X)^-1).
fn spk_lm_fit(y : Matrix, x : Matrix) -> SpkLM raise PropellerError {
  let g = y.rows
  let n = x.rows
  let p = x.cols
  if y.cols != n || p == 0 || n <= p {
    raise PropellerError::PropellerError(
      "design matrix is incompatible with the data or not full rank",
    )
  }
  let xt = matrix_transpose(x)
  let xtx = mat_mult(xt, x)
  let cov = mat_inverse(xtx)
  if cov.rows != p {
    raise PropellerError::PropellerError("singular design matrix")
  }
  let map = mat_mult(cov, xt) // p x n
  let beta = Matrix::zeros(g, p)
  let sigma2 = Array::make(g, 0.0)
  for i in 0..<g {
    // beta row = map (p x n) times y column-vector.
    let mut yy = 0.0
    for j in 0..<p {
      let mut b = 0.0
      for a in 0..<n {
        b += map.get(j, a) * y.get(i, a)
      }
      beta.set(i, j, b)
    }
    for a in 0..<n {
      let yv = y.get(i, a)
      yy += yv * yv
    }
    // RSS = y'y - beta' X'X beta
    let mut cross = 0.0
    for j in 0..<p {
      let mut row = 0.0
      for l in 0..<p {
        row += xtx.get(j, l) * beta.get(i, l)
      }
      cross += beta.get(i, j) * row
    }
    let rss = yy - cross
    sigma2[i] = rss / (n - p).to_double()
  }
  SpkLM::{ beta, sigma2, cov, df_resid: (n - p).to_double(), }
}

///|
/// Columns of ``x`` selected by 0-based indices.
fn spk_select_columns(x : Matrix, cols : Array[Int]) -> Matrix {
  let out = Matrix::zeros(x.rows, cols.length())
  for i in 0..<x.rows {
    for j in 0..<cols.length() {
      out.set(i, j, x.get(i, cols[j]))
    }
  }
  out
}

// ---------------------------------------------------------------------------
// Empirical Bayes variance squeezing
// ---------------------------------------------------------------------------

///|
/// log(x) - digamma(x) (statmod logmdigamma).
fn spk_logmdigamma(x : Double) -> Double {
  @math.ln(x) - gst_digamma(x)
}

///|
/// limma fitFDist with a common numerator df, no covariate.
/// Returns (scale s20, df2); df2 is represented by a very large finite value
/// when R would return Inf.
fn spk_fit_f_dist(x : Array[Double], df1 : Double) -> (Double, Double) {
  let xs = []
  for v in x {
    if !v.is_nan() && v.abs() <= 1.7976931348623157e308 && v > -1.0e-15 {
      xs.push(if v < 0.0 { 0.0 } else { v })
    }
  }
  let n = xs.length()
  if n <= 1 {
    return (if n == 1 { xs[0] } else { 0.0 / 0.0 }, 0.0)
  }
  let sorted = xs.copy()
  sorted.sort()
  let med = if n % 2 == 1 {
    sorted[n / 2]
  } else {
    0.5 * (sorted[n / 2 - 1] + sorted[n / 2])
  }
  let m = if med == 0.0 { 1.0 } else { med }
  let floor = 1.0e-5 * m
  let nd = n.to_double()
  let mut sum_x = 0.0
  let es = Array::make(n, 0.0)
  let mut sum_e = 0.0
  for i in 0..<n {
    let xi = if xs[i] < floor { floor } else { xs[i] }
    sum_x += xi
    let e = @math.ln(xi) + spk_logmdigamma(df1 / 2.0)
    es[i] = e
    sum_e += e
  }
  let emean = sum_e / nd
  let mut ss = 0.0
  for e in es {
    let d = e - emean
    ss += d * d
  }
  let evar = ss / (n - 1).to_double() - gst_trigamma(df1 / 2.0)
  if evar > 0.0 {
    let df2 = 2.0 * gst_trigamma_inverse(evar)
    let s20 = @math.exp(emean - spk_logmdigamma(df2 / 2.0))
    (s20, df2)
  } else {
    (sum_x / nd, 1.0e300)
  }
}

///|
/// R mean(x, trim = t): equal counts removed from each tail.
fn spk_trimmed_mean(z : Array[Double], trim : Double) -> Double {
  let sorted = z.copy()
  sorted.sort()
  let n = sorted.length()
  let lo = (n.to_double() * trim).floor().to_int()
  let hi = n - lo // exclusive upper bound; drop ``lo`` values per tail
  let mut sum = 0.0
  for i in lo..<hi {
    sum += sorted[i]
  }
  sum / (hi - lo).to_double()
}

///|
/// R quantile type 7 (default).
fn spk_quantile7(sorted_asc : Array[Double], p : Double) -> Double {
  let n = sorted_asc.length()
  if n == 1 {
    return sorted_asc[0]
  }
  let h = (n - 1).to_double() * p + 1.0
  let hi = h.ceil().to_int()
  let lo = h.floor().to_int()
  if hi == lo {
    sorted_asc[lo - 1]
  } else {
    sorted_asc[lo - 1] +
    (h - lo.to_double()) * (sorted_asc[hi - 1] - sorted_asc[lo - 1])
  }
}

///|
/// Stable ascending order (R order).
fn spk_order(v : Array[Double]) -> Array[Int] {
  let n = v.length()
  let order = Array::make(n, 0)
  for i in 0..<n {
    order[i] = i
  }
  let mut a = 1
  while a < n {
    let key = order[a]
    let mut b = a - 1
    while b >= 0 && v[order[b]] > v[key] {
      order[b + 1] = order[b]
      b -= 1
    }
    order[b + 1] = key
    a += 1
  }
  order
}

///|
/// 128-node Gauss-Legendre quadrature on [0,1] (statmod gauss.quad.prob
/// uniform), computed by Newton iteration on Legendre polynomials.
fn spk_gauss_legendre_128() -> (Array[Double], Array[Double]) {
  let n = 128
  let nodes = Array::make(n, 0.0)
  let weights = Array::make(n, 0.0)
  let m = (n + 1) / 2
  for i in 1..<=m {
    let mut z = @math.cos(
      @math.PI * (i.to_double() - 0.25) / (n.to_double() + 0.5),
    )
    let mut pp = 0.0
    let mut iter = 0
    while iter < 100 {
      let mut p1 = 1.0
      let mut p2 = 0.0
      for j in 1..<=n {
        let p3 = p2
        p2 = p1
        let jd = j.to_double()
        p1 = ((2.0 * jd - 1.0) * z * p2 - (jd - 1.0) * p3) / jd
      }
      pp = n.to_double() * (z * p1 - p2) / (z * z - 1.0)
      let z1 = z
      z = z1 - p1 / pp
      if (z - z1).abs() < 3.0e-15 {
        iter = 100
      } else {
        iter += 1
      }
    }
    nodes[i - 1] = -z
    nodes[n - i] = z
    let w = 2.0 / ((1.0 - z * z) * pp * pp)
    weights[i - 1] = w
    weights[n - i] = w
  }
  // Map [-1,1] to [0,1].
  for i in 0..<n {
    nodes[i] = 0.5 * (nodes[i] + 1.0)
    weights[i] = 0.5 * weights[i]
  }
  (nodes, weights)
}

///|
/// log1p and expm1 with stable series for small arguments.
fn spk_log1p(x : Double) -> Double {
  if x.abs() < 1.0e-4 {
    // log(1+x) = x - x^2/2 + x^3/3 - ...
    let mut term = x
    let mut sum = 0.0
    let mut k = 1
    while k <= 20 {
      if k % 2 == 1 {
        sum += term
      } else {
        sum -= term
      }
      term *= x
      k += 1
    }
    sum
  } else {
    @math.ln(1.0 + x)
  }
}

///|
fn spk_expm1(x : Double) -> Double {
  if x.abs() < 1.0e-4 {
    // exp(x)-1 = x + x^2/2! + ...
    let mut term = x
    let mut sum = x
    let mut k = 2
    while k <= 20 {
      term *= x / k.to_double()
      sum += term
      k += 1
    }
    sum
  } else {
    @math.exp(x) - 1.0
  }
}

///|
/// Upper-tail chi-square probability P(chi^2_df > q), via the regularized
/// upper incomplete gamma function Q(df/2, q/2) (Numerical Recipes).
fn spk_pchisq_upper(q : Double, df : Double) -> Double {
  if q <= 0.0 {
    return 1.0
  }
  let a = df / 2.0
  let x = q / 2.0
  if x < a + 1.0 {
    1.0 - spk_gamma_p_series(a, x)
  } else {
    spk_gamma_q_cf(a, x)
  }
}

///|
fn spk_gamma_p_series(a : Double, x : Double) -> Double {
  if x <= 0.0 {
    return 0.0
  }
  let gln = gst_lngamma(a)
  if x < a + 1.0 {
    let mut ap = a
    let mut sum = 1.0 / a
    let mut del = sum
    for _ in 0..<200 {
      ap += 1.0
      del *= x / ap
      sum += del
      if del.abs() < 1.0e-15 * sum.abs() {
        break
      }
    }
    sum * @math.exp(-x + a * @math.ln(x) - gln)
  } else {
    1.0 - spk_gamma_q_cf(a, x)
  }
}

///|
fn spk_gamma_q_cf(a : Double, x : Double) -> Double {
  let gln = gst_lngamma(a)
  let fpmin = 1.0e-300
  let mut b = x + 1.0 - a
  let mut c = 1.0 / fpmin
  let mut d = 1.0 / b
  let mut h = d
  for i in 1..<=200 {
    let an = -i.to_double() * (i.to_double() - a)
    b += 2.0
    d = an * d + b
    if d.abs() < fpmin {
      d = fpmin
    }
    c = b + an / c
    if c.abs() < fpmin {
      c = fpmin
    }
    d = 1.0 / d
    let del = d * c
    h *= del
    if (del - 1.0).abs() < 1.0e-15 {
      break
    }
  }
  @math.exp(-x + a * @math.ln(x) - gln) * h
}

///|
/// F distribution lower-tail CDF P(F <= f).
fn spk_pf_lower(f : Double, df1 : Double, df2 : Double) -> Double {
  if f <= 0.0 {
    return 0.0
  }
  if df2 > 1.0e100 {
    // chi-square limit: df1*F ~ chi^2(df1)
    return 1.0 - spk_pchisq_upper(df1 * f, df1)
  }
  let x = df1 * f / (df1 * f + df2)
  gst_pbeta(x, df1 / 2.0, df2 / 2.0)
}

///|
/// F distribution upper-tail CDF P(F > f).
fn spk_pf_upper(f : Double, df1 : Double, df2 : Double) -> Double {
  if f <= 0.0 {
    return 1.0
  }
  if df2 > 1.0e100 {
    return spk_pchisq_upper(df1 * f, df1)
  }
  // P(F>f) = I_{df2/(df2+df1 f)}(df2/2, df1/2)
  let x = df2 / (df2 + df1 * f)
  gst_pbeta(x, df2 / 2.0, df1 / 2.0)
}

///|
/// Log of the F upper-tail probability (used by the finite-df2 robust path).
fn spk_log_pf_upper(f : Double, df1 : Double, df2 : Double) -> Double {
  let tail = spk_pf_upper(f, df1, df2)
  if tail > 1.0e-300 {
    @math.ln(tail)
  } else {
    // Stable evaluation of the regularized beta complement.
    let a = df2 / 2.0
    let b = df1 / 2.0
    let x = df2 / (df2 + df1 * f) // near 1
    let y = 1.0 - x
    let bt = gst_lngamma(a + b) -
      gst_lngamma(a) -
      gst_lngamma(b) +
      a * @math.ln(x) +
      b * @math.ln(y)
    let small = @math.exp(bt) * gst_betacf(a, b, x) / a
    spk_log1p(-small)
  }
}

///|
/// F quantile (lower tail) by bisection on the log scale.
fn spk_qf_lower(p : Double, df1 : Double, df2 : Double) -> Double {
  if p <= 0.0 {
    return 0.0
  }
  if p >= 1.0 {
    return 1.0e300
  }
  if df2 > 1.0e100 {
    return spk_qchisq_lower(p, df1) / df1
  }
  let mut lo = -700.0
  let mut hi = 0.0
  while spk_pf_lower(@math.exp(hi), df1, df2) < p {
    hi += 2.0
  }
  for _ in 0..<100 {
    let mid = 0.5 * (lo + hi)
    if spk_pf_lower(@math.exp(mid), df1, df2) < p {
      lo = mid
    } else {
      hi = mid
    }
  }
  @math.exp(0.5 * (lo + hi))
}

///|
/// Chi-square quantile (lower tail) by bisection.
fn spk_qchisq_lower(p : Double, df : Double) -> Double {
  let mut lo = 0.0
  let mut hi = 1.0
  while spk_gamma_p_series(df / 2.0, hi / 2.0) < p {
    hi *= 2.0
  }
  for _ in 0..<100 {
    let mid = 0.5 * (lo + hi)
    if spk_gamma_p_series(df / 2.0, mid / 2.0) < p {
      lo = mid
    } else {
      hi = mid
    }
  }
  0.5 * (lo + hi)
}

///|
/// F density, supporting infinite denominator df (chi-square limit).
fn spk_df(f : Double, df1 : Double, df2 : Double) -> Double {
  if df2 > 1.0e100 {
    // density for F = chi2(df1)/df1
    let y = df1 * f
    let a = df1 / 2.0
    return df1 *
      @math.exp(
        -y / 2.0 + (a - 1.0) * @math.ln(y) - a * @math.ln(2.0) - gst_lngamma(a),
      )
  }
  let a = df1 / 2.0
  let b = df2 / 2.0
  let lbeta = gst_lngamma(a) + gst_lngamma(b) - gst_lngamma(a + b)
  @math.exp(
    a * @math.ln(df1 / df2) +
    (a - 1.0) * @math.ln(f) -
    (a + b) * @math.ln(1.0 + df1 * f / df2) -
    lbeta,
  )
}

///|
/// Winsorized log(F) moments used by fitFDistRobustly.
/// Returns (mean, variance).
fn spk_winsorized_moments(
  df1 : Double,
  df2 : Double,
  nodes : Array[Double],
  weights : Array[Double],
) -> (Double, Double) {
  let wtp = [0.05, 0.10]
  let fq1 = spk_qf_lower(0.05, df1, df2)
  let fq2 = spk_qf_lower(0.90, df1, df2)
  let zq1 = @math.ln(fq1)
  let zq2 = @math.ln(fq2)
  let q1 = fq1 / (1.0 + fq1)
  let q2 = fq2 / (1.0 + fq2)
  let q21 = q2 - q1
  let mut sum_fz = 0.0
  let mut sum_fz2 = 0.0
  let mut sum_f = 0.0
  for a in 0..<nodes.length() {
    let node = q1 + q21 * nodes[a]
    let fnod = node / (1.0 - node)
    let znod = @math.ln(fnod)
    let dens = spk_df(fnod, df1, df2) / (1.0 - node) / (1.0 - node)
    let fw = weights[a] * dens
    sum_f += fw
    sum_fz += fw * znod
    sum_fz2 += fw * znod * znod
  }
  // m = q21 * sum(w f z) + sum(zq * tail.p)
  let mm = q21 * sum_fz + zq1 * wtp[0] + zq2 * wtp[1]
  let vv = q21 * (sum_fz2 - 2.0 * mm * sum_fz + mm * mm * sum_f) +
    (zq1 - mm) * (zq1 - mm) * wtp[0] +
    (zq2 - mm) * (zq2 - mm) * wtp[1]
  (mm, vv)
}

///|
/// limma fitFDistRobustly for a common numerator df, no covariate.
/// Returns (scale, df2, df2.shrunk per observation).
fn spk_fit_f_dist_robustly(
  x_in : Array[Double],
  df1 : Double,
) -> (Double, Double, Array[Double]) raise PropellerError {
  let n0 = x_in.length()
  let x = x_in.copy()
  if n0 == 2 {
    let (s, d) = spk_fit_f_dist(x, df1)
    return (s, d, Array::make(n0, d))
  }
  let n = n0
  let nd = n.to_double()
  // Offset very small variances away from zero.
  let sorted_x = x.copy()
  sorted_x.sort()
  let med = if n % 2 == 1 {
    sorted_x[n / 2]
  } else {
    0.5 * (sorted_x[n / 2 - 1] + sorted_x[n / 2])
  }
  if med <= 0.0 {
    raise PropellerError::PropellerError(
      "Variances are mostly <= 0 in robust eBayes",
    )
  }
  let tiny = med * 1.0e-12
  for i in 0..<n {
    if x[i] < tiny {
      x[i] = tiny
    }
  }
  let (non_scale, non_df2) = spk_fit_f_dist(x, df1)
  // Not enough observations for tail winsorization: identical to non-robust.
  if 0.10 < 1.0 / nd {
    return (non_scale, non_df2, Array::make(n, non_df2))
  }
  let z = Array::make(n, 0.0)
  for i in 0..<n {
    z[i] = @math.ln(x[i])
  }
  let ztrend = spk_trimmed_mean(z, 0.1)
  let zr = Array::make(n, 0.0)
  for i in 0..<n {
    zr[i] = z[i] - ztrend
  }
  let sorted_r = zr.copy()
  sorted_r.sort()
  let zlo = spk_quantile7(sorted_r, 0.05)
  let zhi = spk_quantile7(sorted_r, 0.90)
  let zw = Array::make(n, 0.0)
  let mut zwmean = 0.0
  for i in 0..<n {
    let v = if zr[i] < zlo { zlo } else if zr[i] > zhi { zhi } else { zr[i] }
    zw[i] = v
    zwmean += v
  }
  zwmean /= nd
  let mut zwvar = 0.0
  for v in zw {
    let d = v - zwmean
    zwvar += d * d
  }
  zwvar = zwvar / nd * nd / (nd - 1.0)
  let (gnodes, gweights) = spk_gauss_legendre_128()
  let (m_inf, v_inf) = spk_winsorized_moments(df1, 1.0e300, gnodes, gweights)
  let funval_inf = @math.ln(zwvar / v_inf)
  if funval_inf <= 0.0 {
    let ztc = ztrend + zwmean - m_inf
    let s20 = @math.exp(ztc)
    let tailp = Array::make(n, 0.0)
    let fstat = Array::make(n, 0.0)
    for i in 0..<n {
      fstat[i] = @math.exp(z[i] - ztc)
      tailp[i] = spk_pchisq_upper(fstat[i] * df1, df1)
    }
    let ranks = gst_average_ranks(fstat)
    let df_pooled = nd * df1
    let shrunk = Array::make(n, 1.0e300)
    for i in 0..<n {
      let emp = (nd - ranks[i] + 0.5) / nd
      let pno = tailp[i] / emp
      if pno < 1.0 {
        shrunk[i] = pno * df_pooled
      }
    }
    let order = spk_order(tailp)
    let mut running = 0.0
    for oi in order {
      if shrunk[oi] > running {
        running = shrunk[oi]
      }
      shrunk[oi] = running
    }
    return (s20, 1.0e300, shrunk)
  }
  if non_df2 > 1.0e100 {
    return (non_scale, non_df2, Array::make(n, non_df2))
  }
  // Estimate finite df2 by matching the winsorized variance.
  let linkfun = fn(d : Double) -> Double { d / (1.0 + d) }
  let linkinv = fn(u : Double) -> Double { u / (1.0 - u) }
  let rbx = linkfun(non_df2)
  let (_, v_low) = spk_winsorized_moments(df1, non_df2, gnodes, gweights)
  let funval_low = @math.ln(zwvar / v_low)
  let mut df2 = non_df2
  if funval_low < 0.0 {
    // uniroot on (rbx, 1), f.lower < 0, f.upper = funvalInf > 0.
    let mut lo = rbx
    let mut hi = 1.0 - 1.0e-12
    for _ in 0..<100 {
      let mid = 0.5 * (lo + hi)
      let (_, vv) = spk_winsorized_moments(df1, linkinv(mid), gnodes, gweights)
      if @math.ln(zwvar / vv) < 0.0 {
        lo = mid
      } else {
        hi = mid
      }
      if hi - lo < 1.0e-13 {
        break
      }
    }
    df2 = linkinv(0.5 * (lo + hi))
  }
  let (mm, _) = spk_winsorized_moments(df1, df2, gnodes, gweights)
  let ztc = ztrend + zwmean - mm
  let s20 = @math.exp(ztc)
  let fstat = Array::make(n, 0.0)
  let log_tail = Array::make(n, 0.0)
  let mut any_outlier = false
  let mut min_log_tail = 0.0
  for i in 0..<n {
    fstat[i] = @math.exp(z[i] - ztc)
    log_tail[i] = spk_log_pf_upper(fstat[i], df1, df2)
    if log_tail[i] < min_log_tail {
      min_log_tail = log_tail[i]
    }
  }
  let ranks = gst_average_ranks(fstat)
  let log_emp = Array::make(n, 0.0)
  let log_pno = Array::make(n, 0.0)
  for i in 0..<n {
    log_emp[i] = @math.ln((nd - ranks[i] + 0.5) / nd)
    log_pno[i] = if log_tail[i] - log_emp[i] < 0.0 {
      log_tail[i] - log_emp[i]
    } else {
      0.0
    }
    if log_pno[i] < 0.0 {
      any_outlier = true
    }
  }
  let shrunk = Array::make(n, df2)
  if any_outlier {
    let mut max_f = 0.0
    for f in fstat {
      if f > max_f {
        max_f = f
      }
    }
    let mut df2_out = if min_log_tail < -1.0e300 {
      0.0
    } else {
      @math.ln(0.5) / min_log_tail * df2
    }
    if df2_out > 0.0 {
      // One refinement iteration (R fitFDistRobustly).
      let new_log_tail = spk_log_pf_upper(max_f, df1, df2_out)
      df2_out = @math.ln(0.5) / new_log_tail * df2_out
    }
    for i in 0..<n {
      let pno = @math.exp(log_pno[i])
      let pout = -spk_expm1(log_pno[i])
      shrunk[i] = pno * df2 + pout * df2_out
    }
    // Enforce monotonicity in tail probability.
    let order = spk_order(log_tail)
    let ordered = Array::make(n, 0.0)
    for k in 0..<n {
      ordered[k] = shrunk[order[k]]
    }
    let mut cum = 0.0
    let cmean = Array::make(n, 0.0)
    for k in 0..<n {
      cum += ordered[k]
      cmean[k] = cum / (k + 1).to_double()
    }
    let mut imin = 0
    for k in 1..<n {
      if cmean[k] < cmean[imin] {
        imin = k
      }
    }
    for k in 0..<=imin {
      ordered[k] = cmean[imin]
    }
    let mut running = 0.0
    for k in 0..<n {
      if ordered[k] > running {
        running = ordered[k]
      }
      shrunk[order[k]] = running
    }
  }
  (s20, df2, shrunk)
}

///|
/// Posterior variances and total df from residual variances.
/// Returns (df.prior per row, s2.prior, s2.post, df.total per row).
fn spk_squeeze(
  variance : Array[Double],
  df : Double,
  robust : Bool,
) -> (Array[Double], Double, Array[Double], Array[Double]) raise PropellerError {
  let n = variance.length()
  let pooled = n.to_double() * df
  if n < 3 {
    // limma: no shrinkage.
    let dft = Array::make(n, df)
    return (Array::make(n, 0.0), 0.0, variance.copy(), dft)
  }
  let mut df_prior = Array::make(n, 0.0)
  let mut s2_prior = 0.0
  if robust {
    let (s, d2, shrunk) = spk_fit_f_dist_robustly(variance, df)
    s2_prior = s
    df_prior = shrunk
    let _ = d2
  } else {
    let (s, d2) = spk_fit_f_dist(variance, df)
    s2_prior = s
    for i in 0..<n {
      df_prior[i] = d2
    }
  }
  let s2_post = Array::make(n, 0.0)
  let df_total = Array::make(n, 0.0)
  for i in 0..<n {
    if df_prior[i] > 1.0e100 {
      s2_post[i] = s2_prior
      df_total[i] = pooled
    } else {
      s2_post[i] = (df * variance[i] + df_prior[i] * s2_prior) /
        (df + df_prior[i])
      let tot = df + df_prior[i]
      df_total[i] = if tot < pooled { tot } else { pooled }
    }
  }
  (df_prior, s2_prior, s2_post, df_total)
}

// ---------------------------------------------------------------------------
// Symmetric eigendecomposition (cyclic Jacobi) for the moderated F test
// ---------------------------------------------------------------------------

///|
/// Eigenvalues (descending) and orthonormal eigenvectors (columns) of a real
/// symmetric matrix, via cyclic Jacobi rotations.
fn spk_eigen_symmetric(a : Matrix) -> (Array[Double], Matrix) {
  let r = a.rows
  let v = Matrix::identity(r)
  let work = Matrix::zeros(r, r)
  for i in 0..<r {
    for j in 0..<r {
      work.set(i, j, a.get(i, j))
    }
  }
  for _sweep in 0..<100 {
    let mut off = 0.0
    for p in 0..<r {
      for q in (p + 1)..<r {
        off += work.get(p, q) * work.get(p, q)
      }
    }
    if off < 1.0e-30 {
      break
    }
    for p in 0..<r {
      for q in (p + 1)..<r {
        let apq = work.get(p, q)
        if apq.abs() > 1.0e-300 {
          let app = work.get(p, p)
          let aqq = work.get(q, q)
          let tau = (aqq - app) / (2.0 * apq)
          let t = if tau >= 0.0 {
            1.0 / (tau + (1.0 + tau * tau).sqrt())
          } else {
            -1.0 / (-tau + (1.0 + tau * tau).sqrt())
          }
          let c = 1.0 / (1.0 + t * t).sqrt()
          let s = t * c
          // Rotate rows/columns p,q.
          for k in 0..<r {
            let akp = work.get(k, p)
            let akq = work.get(k, q)
            work.set(k, p, c * akp - s * akq)
            work.set(k, q, s * akp + c * akq)
          }
          for k in 0..<r {
            let apk = work.get(p, k)
            let aqk = work.get(q, k)
            work.set(p, k, c * apk - s * aqk)
            work.set(q, k, s * apk + c * aqk)
          }
          for k in 0..<r {
            let vkp = v.get(k, p)
            let vkq = v.get(k, q)
            v.set(k, p, c * vkp - s * vkq)
            v.set(k, q, s * vkp + c * vkq)
          }
        }
      }
    }
  }
  let values = Array::make(r, 0.0)
  for i in 0..<r {
    values[i] = work.get(i, i)
  }
  // Sort eigenpairs descending.
  let mut a = 0
  while a < r {
    let mut b = a + 1
    while b < r {
      if values[b] > values[a] {
        let tv = values[a]
        values[a] = values[b]
        values[b] = tv
        for k in 0..<r {
          let te = v.get(k, a)
          v.set(k, a, v.get(k, b))
          v.set(k, b, te)
        }
      }
      b += 1
    }
    a += 1
  }
  (values, v)
}

// ---------------------------------------------------------------------------
// Moderated t-tests and F-tests
// ---------------------------------------------------------------------------

///|
/// Core moderated t-test on precomputed matrices. ``transformed`` and
/// ``proportions`` have cell types on rows and samples on columns; the design
/// has samples on rows. ``contrasts`` selects the two groups (typically
/// ``[1.0, -1.0]``) and may contain zeros for confounding covariates.
pub fn propeller_ttest_model(
  transformed : Matrix,
  proportions : Matrix,
  design : Matrix,
  contrasts : Array[Double],
  robust? : Bool = true,
) -> PropellerTModel raise PropellerError {
  let g = transformed.rows
  let p = design.cols
  if transformed.cols != design.rows {
    raise PropellerError::PropellerError(
      "transformed props and design have mismatched sample counts",
    )
  }
  if contrasts.length() != p {
    raise PropellerError::PropellerError(
      "contrasts length must equal the number of design columns",
    )
  }
  let use_robust = robust && g > 2
  let fit = spk_lm_fit(transformed, design)
  let (df_prior, s2_prior, s2_post, df_total) = spk_squeeze(
    fit.sigma2,
    fit.df_resid,
    use_robust,
  )
  // contrasts.fit: beta_c = C' beta, unscaled variance C' (X'X)^-1 C
  let mut var_c = 0.0
  for j in 0..<p {
    for l in 0..<p {
      var_c += contrasts[j] * contrasts[l] * fit.cov.get(j, l)
    }
  }
  let statistic = Array::make(g, 0.0)
  let p_value = Array::make(g, 0.0)
  for i in 0..<g {
    let mut bc = 0.0
    for j in 0..<p {
      bc += contrasts[j] * fit.beta.get(i, j)
    }
    let se = (var_c * s2_post[i]).sqrt()
    let t = bc / se
    statistic[i] = t
    p_value[i] = 2.0 * gst_pt_lower(-t.abs(), df_total[i])
  }
  // Group means and prop ratios: drop zero-contrast columns when confounded.
  let eff_cols = []
  let eff_cont = []
  for j in 0..<p {
    if contrasts[j] != 0.0 {
      eff_cols.push(j)
      eff_cont.push(contrasts[j])
    }
  }
  let means = if eff_cols.length() == p {
    spk_lm_fit(proportions, design).beta
  } else {
    spk_lm_fit(proportions, spk_select_columns(design, eff_cols)).beta
  }
  let prop_ratio = Array::make(g, 1.0)
  for i in 0..<g {
    let mut rr = 1.0
    for j in 0..<eff_cols.length() {
      rr *= @math.pow(means.get(i, j), eff_cont[j])
    }
    prop_ratio[i] = rr
  }
  let fdr = gst_bh(p_value)
  PropellerTModel::{
    means,
    prop_ratio,
    statistic,
    p_value,
    fdr,
    df_prior,
    s2_prior,
    df_total,
    s2_post,
  }
}

///|
/// Moderated t-test on the output of ``get_transformed_props``.
pub fn propeller_ttest(
  pl : TransformedProps,
  design : Matrix,
  contrasts : Array[Double],
  robust? : Bool = true,
  sort? : Bool = true,
) -> Array[PropellerTRow] raise PropellerError {
  let model = propeller_ttest_model(
    pl.transformed,
    pl.proportions,
    design,
    contrasts,
    robust~,
  )
  let g = pl.transformed.rows
  let idx = Array::make(g, 0)
  for i in 0..<g {
    idx[i] = i
  }
  if sort {
    let mut a = 1
    while a < g {
      let key = idx[a]
      let mut b = a - 1
      while b >= 0 && model.p_value[idx[b]] > model.p_value[key] {
        idx[b + 1] = idx[b]
        b -= 1
      }
      idx[b + 1] = key
      a += 1
    }
  }
  let out = []
  for i in idx {
    let ng = model.means.cols
    let means = Array::make(ng, 0.0)
    for j in 0..<ng {
      means[j] = model.means.get(i, j)
    }
    out.push(PropellerTRow::{
      cluster: pl.cluster_names[i],
      means,
      prop_ratio: model.prop_ratio[i],
      statistic: model.statistic[i],
      p_value: model.p_value[i],
      fdr: model.fdr[i],
    })
  }
  out
}

///|
/// Core moderated ANOVA F-test on precomputed matrices. ``coef`` gives the
/// 0-based design columns corresponding to the groups (default: all
/// columns). The first group column is replaced by an intercept internally,
/// exactly as speckle does.
pub fn propeller_anova_model(
  transformed : Matrix,
  proportions : Matrix,
  design : Matrix,
  coef? : Array[Int]? = None,
  robust? : Bool = true,
) -> PropellerAnovaModel raise PropellerError {
  let g = transformed.rows
  let n = design.rows
  let p = design.cols
  if transformed.cols != n {
    raise PropellerError::PropellerError(
      "transformed props and design have mismatched sample counts",
    )
  }
  let groups = match coef {
    Some(c) => c
    None => {
      let all = []
      for j in 0..<p {
        all.push(j)
      }
      all
    }
  }
  if groups.length() < 2 {
    raise PropellerError::PropellerError(
      "ANOVA requires at least two group columns",
    )
  }
  for c in groups {
    if c < 0 || c >= p {
      raise PropellerError::PropellerError("coef index out of range")
    }
  }
  // Mean proportions on the raw group columns.
  let group_design = spk_select_columns(design, groups)
  let means = spk_lm_fit(proportions, group_design).beta
  // Intercept-form design for the transformed model.
  let design2 = Matrix::zeros(n, p)
  for i in 0..<n {
    for j in 0..<p {
      design2.set(i, j, design.get(i, j))
    }
    design2.set(i, groups[0], 1.0)
  }
  let fit = spk_lm_fit(transformed, design2)
  let use_robust = robust && g > 2
  let (df_prior, s2_prior, s2_post, df_total) = spk_squeeze(
    fit.sigma2,
    fit.df_resid,
    use_robust,
  )
  // Tested columns: all group columns except the (now intercept) first.
  let tested = []
  for k in 1..<groups.length() {
    tested.push(groups[k])
  }
  let r = tested.length()
  // t statistics for the tested coefficients.
  let tmat = Matrix::zeros(g, r)
  for i in 0..<g {
    for a in 0..<r {
      let j = tested[a]
      let se = (fit.cov.get(j, j) * s2_post[i]).sqrt()
      tmat.set(i, a, fit.beta.get(i, j) / se)
    }
  }
  let statistic = Array::make(g, 0.0)
  let p_value = Array::make(g, 0.0)
  // classifyTestsF attributes df2 = df.prior + df.residual WITHOUT the pooled
  // pmin cap used for moderated t (limma contrasts/decidetests + ebayes).
  let f_df2 = Array::make(g, 0.0)
  for i in 0..<g {
    f_df2[i] = if df_prior[i] > 1.0e200 {
      1.0e300
    } else {
      fit.df_resid + df_prior[i]
    }
  }
  if r == 1 {
    for i in 0..<g {
      let f = tmat.get(i, 0) * tmat.get(i, 0)
      statistic[i] = f
      p_value[i] = spk_pf_upper(f, 1.0, f_df2[i])
    }
  } else {
    // Correlation matrix of the tested coefficients.
    let covsub = Matrix::zeros(r, r)
    for a in 0..<r {
      for b in 0..<r {
        covsub.set(a, b, fit.cov.get(tested[a], tested[b]))
      }
    }
    let cormat = Matrix::zeros(r, r)
    for a in 0..<r {
      for b in 0..<r {
        let va = covsub.get(a, a)
        let vb = covsub.get(b, b)
        cormat.set(a, b, covsub.get(a, b) / (va * vb).sqrt())
      }
    }
    let (evalues, evectors) = spk_eigen_symmetric(cormat)
    let mut rr = 0
    for a in 0..<r {
      if evalues[a] / evalues[0] > 1.0e-8 {
        rr += 1
      }
    }
    // Q = V[,1:rr] * (1/sqrt(lambda[1:rr])) / sqrt(rr)
    let qmat = Matrix::zeros(r, rr)
    for a in 0..<r {
      for b in 0..<rr {
        qmat.set(
          a,
          b,
          evectors.get(a, b) / evalues[b].sqrt() / rr.to_double().sqrt(),
        )
      }
    }
    for i in 0..<g {
      let mut f = 0.0
      for b in 0..<rr {
        let mut dot = 0.0
        for a in 0..<r {
          dot += tmat.get(i, a) * qmat.get(a, b)
        }
        f += dot * dot
      }
      statistic[i] = f
      p_value[i] = spk_pf_upper(f, rr.to_double(), f_df2[i])
    }
  }
  let fdr = gst_bh(p_value)
  PropellerAnovaModel::{
    means,
    statistic,
    p_value,
    fdr,
    df_prior,
    s2_prior,
    df_total,
    s2_post,
  }
}

///|
/// Moderated ANOVA F-test on the output of ``get_transformed_props``.
pub fn propeller_anova(
  pl : TransformedProps,
  design : Matrix,
  coef? : Array[Int]? = None,
  robust? : Bool = true,
  sort? : Bool = true,
) -> Array[PropellerAnovaRow] raise PropellerError {
  let model = propeller_anova_model(
    pl.transformed,
    pl.proportions,
    design,
    coef~,
    robust~,
  )
  let g = pl.transformed.rows
  let idx = Array::make(g, 0)
  for i in 0..<g {
    idx[i] = i
  }
  if sort {
    let mut a = 1
    while a < g {
      let key = idx[a]
      let mut b = a - 1
      while b >= 0 && model.p_value[idx[b]] > model.p_value[key] {
        idx[b + 1] = idx[b]
        b -= 1
      }
      idx[b + 1] = key
      a += 1
    }
  }
  let out = []
  for i in idx {
    let ng = model.means.cols
    let means = Array::make(ng, 0.0)
    for j in 0..<ng {
      means[j] = model.means.get(i, j)
    }
    out.push(PropellerAnovaRow::{
      cluster: pl.cluster_names[i],
      means,
      statistic: model.statistic[i],
      p_value: model.p_value[i],
      fdr: model.fdr[i],
    })
  }
  out
}

// ---------------------------------------------------------------------------
// High-level wrapper
// ---------------------------------------------------------------------------

///|
/// One-call differential proportions analysis, mirroring speckle
/// ``propeller``. Cells carry ``clusters`` (cell type), ``samples``
/// (biological replicate) and ``groups`` (condition). Samples are assigned
/// to groups from the cell-level table; two groups trigger moderated
/// t-tests, three or more groups moderated ANOVA. Results are sorted by
/// p-value and annotated with the global (baseline) proportion of each cell
/// type. ``transform`` defaults to ``"logit"``.
pub fn propeller(
  clusters : Array[String],
  samples : Array[String],
  groups : Array[String],
  transform? : String = "logit",
  robust? : Bool = true,
) -> Array[PropellerRow] raise PropellerError {
  let n = clusters.length()
  if n == 0 || samples.length() != n || groups.length() != n {
    raise PropellerError::PropellerError(
      "clusters, samples and groups must be non-empty and equal length",
    )
  }
  let pl = get_transformed_props(clusters, samples, transform~)
  let sample_names = spk_sorted_unique(samples)
  let group_names = spk_sorted_unique(groups)
  let ns = sample_names.length()
  let ng = group_names.length()
  if ng < 2 {
    raise PropellerError::PropellerError("at least two groups are required")
  }
  let coll = Array::make(ns * ng, 0.0)
  for cell in 0..<n {
    let i = spk_find_sorted(sample_names, samples[cell])
    let j = spk_find_sorted(group_names, groups[cell])
    if i >= 0 && j >= 0 {
      coll[i * ng + j] = 1.0
    }
  }
  let design = Matrix::new(ns, ng, coll)
  // Baseline proportions: table(clusters) / total cells.
  let k = pl.cluster_names.length()
  let baseline = Array::make(k, 0.0)
  for c in clusters {
    let i = spk_find_sorted(pl.cluster_names, c)
    if i >= 0 {
      baseline[i] += 1.0
    }
  }
  for i in 0..<k {
    baseline[i] /= n.to_double()
  }
  let rows = []
  if ng == 2 {
    let model = propeller_ttest_model(
      pl.transformed,
      pl.proportions,
      design,
      [1.0, -1.0],
      robust~,
    )
    let order = spk_order(model.p_value)
    for i in order {
      let means = Array::make(ng, 0.0)
      for j in 0..<ng {
        means[j] = model.means.get(i, j)
      }
      rows.push(PropellerRow::{
        cluster: pl.cluster_names[i],
        kind: "t",
        baseline_prop: baseline[i],
        means,
        prop_ratio: Some(model.prop_ratio[i]),
        statistic: model.statistic[i],
        p_value: model.p_value[i],
        fdr: model.fdr[i],
      })
    }
  } else {
    let model = propeller_anova_model(
      pl.transformed,
      pl.proportions,
      design,
      robust~,
    )
    let order = spk_order(model.p_value)
    for i in order {
      let means = Array::make(ng, 0.0)
      for j in 0..<ng {
        means[j] = model.means.get(i, j)
      }
      rows.push(PropellerRow::{
        cluster: pl.cluster_names[i],
        kind: "F",
        baseline_prop: baseline[i],
        means,
        prop_ratio: None,
        statistic: model.statistic[i],
        p_value: model.p_value[i],
        fdr: model.fdr[i],
      })
    }
  }
  rows
}

// ---------------------------------------------------------------------------
// Beta distribution helpers and library-size normalisation
// ---------------------------------------------------------------------------

///|
/// Method-of-moments estimates of Beta(alpha, beta) parameters for a vector
/// of proportions (speckle ``estimateBetaParam``). The sample variance uses
/// the n-1 denominator, matching R ``var()``.
pub fn estimate_beta_param(x : Array[Double]) -> BetaMom raise PropellerError {
  let n = x.length()
  if n < 2 {
    raise PropellerError::PropellerError(
      "at least two proportions are required",
    )
  }
  let mut mu = 0.0
  for v in x {
    mu += v
  }
  mu /= n.to_double()
  let mut vv = 0.0
  for v in x {
    let d = v - mu
    vv += d * d
  }
  vv /= (n - 1).to_double()
  let k = (1.0 - mu) / vv - 1.0 / mu
  BetaMom::{ alpha: k * mu * mu, beta: k * mu * (1.0 - mu), }
}

///|
/// Normalise a counts matrix (cell types x samples) to the median library
/// size (speckle ``normCounts``). With ``log=true``, scaled prior counts are
/// added before the log2 transformation.
pub fn propeller_norm_counts(
  counts : Matrix,
  log? : Bool = false,
  prior_count? : Double = 0.5,
  lib_size? : Array[Double]? = None,
) -> Matrix {
  let k = counts.rows
  let s = counts.cols
  let lib = Array::make(s, 0.0)
  match lib_size {
    Some(ls) if ls.length() == s =>
      for j in 0..<s {
        lib[j] = ls[j]
      }
    _ =>
      for j in 0..<s {
        let mut total = 0.0
        for i in 0..<k {
          total += counts.get(i, j)
        }
        lib[j] = total
      }
  }
  let sorted = lib.copy()
  sorted.sort()
  let m = if s % 2 == 1 {
    sorted[s / 2]
  } else {
    0.5 * (sorted[s / 2 - 1] + sorted[s / 2])
  }
  let out = Matrix::zeros(k, s)
  if log {
    let mut mean_lib = 0.0
    for v in lib {
      mean_lib += v
    }
    mean_lib /= s.to_double()
    for j in 0..<s {
      let pcs = lib[j] / mean_lib * prior_count
      let ls2 = lib[j] + 2.0 * pcs
      for i in 0..<k {
        out.set(i, j, @math.log2((counts.get(i, j) + pcs) / ls2 * m))
      }
    }
  } else {
    for j in 0..<s {
      for i in 0..<k {
        out.set(i, j, counts.get(i, j) / lib[j] * m)
      }
    }
  }
  out
}

///|
/// Estimate Beta-binomial parameters for every cell type from a counts
/// matrix (speckle ``estimateBetaParamsFromCounts``).
pub fn estimate_beta_params_from_counts(counts : Matrix) -> BetaCountParams {
  let k = counts.rows
  let s = counts.cols
  let nc = propeller_norm_counts(counts)
  let sd = s.to_double()
  let m1 = Array::make(k, 0.0)
  let m2 = Array::make(k, 0.0)
  let mut n = 0.0
  for j in 0..<s {
    let mut col_sum = 0.0
    for i in 0..<k {
      col_sum += nc.get(i, j)
    }
    n += col_sum
  }
  n /= sd
  for i in 0..<k {
    let mut sum = 0.0
    let mut sq = 0.0
    for j in 0..<s {
      let v = nc.get(i, j)
      sum += v
      sq += v * v
    }
    m1[i] = sum / sd
    m2[i] = sq / sd
  }
  let alpha = Array::make(k, 0.0)
  let beta = Array::make(k, 0.0)
  let pi = Array::make(k, 0.0)
  let disp = Array::make(k, 0.0)
  let variance = Array::make(k, 0.0)
  for i in 0..<k {
    let denom = n * (m2[i] / m1[i] - m1[i] - 1.0) + m1[i]
    let a = (n * m1[i] - m2[i]) / denom
    let b = (n - m1[i]) * (n - m2[i] / m1[i]) / denom
    alpha[i] = a
    beta[i] = b
    disp[i] = 1.0 / (a + b)
    pi[i] = a / (a + b)
    variance[i] = n *
      pi[i] *
      (1.0 - pi[i]) *
      (n * disp[i] + 1.0) /
      (1.0 + disp[i])
  }
  BetaCountParams::{ n, alpha, beta, pi, dispersion: disp, variance, }
}
