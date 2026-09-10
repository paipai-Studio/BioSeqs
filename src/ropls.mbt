///|
/// ropls - PCA, PLS(-DA) and OPLS(-DA) multivariate modeling.
///
/// A MoonBit port of the core algorithms of the Bioconductor ropls package
/// (Thevenot et al., 2015): NIPALS PCA, NIPALS PLS1/PLS2, orthogonal PLS
/// (Trygg and Wold, 2002) with VIP (including VIP4,p / VIP4,o for OPLS),
/// k-fold cross-validation (R2X / R2Y / Q2), and response permutation
/// testing (Szymanska et al., 2012).
///
/// All matrices are row-major `Array[Array[Double]]` with samples in rows
/// and variables in columns.

// ============================================================================
// Data structures
// ============================================================================

///|
/// Model response: continuous matrix (samples x responses) or class labels.
pub(all) enum RoplsY {
  Continuous(Array[Array[Double]])
  Classes(Array[String])
}

///|
/// Result of NIPALS PCA.
pub struct RoplsPcaResult {
  n : Int
  p : Int
  ncomp : Int
  scale : String
  /// Scores T (samples x components)
  scores : Array[Array[Double]]
  /// Loadings P (variables x components, unit norm)
  loadings : Array[Array[Double]]
  /// Variances associated with each component (sample variance of scores)
  variance : Array[Double]
  /// Per-component R2X
  r2x : Array[Double]
  /// Cumulative R2X
  r2x_cum : Array[Double]
  /// NIPALS iterations per component
  iters : Array[Int]
  x_mean : Array[Double]
  x_sd : Array[Double]
  /// Indices of removed near-zero variance variables
  dropped : Array[Int]
}

///|
/// A fitted (O)PLS(-DA) model.
pub struct RoplsModel {
  /// "PLS", "PLS-DA", "OPLS" or "OPLS-DA"
  kind : String
  n : Int
  p : Int
  q : Int
  /// Number of predictive components
  pre : Int
  /// Number of orthogonal components
  ort : Int
  scale : String
  crossval : Int
  /// Predictive scores T (n x pre)
  scores : Array[Array[Double]]
  /// X loadings P (p x pre)
  x_loadings : Array[Array[Double]]
  /// X weights W (p x pre, unit norm)
  weights : Array[Array[Double]]
  /// Weight star W* (p x pre)
  weight_star : Array[Array[Double]]
  y_scores : Array[Array[Double]]
  /// Y loadings C (q x pre)
  y_loadings : Array[Array[Double]]
  /// Regression coefficients B (p x q)
  coefficients : Array[Array[Double]]
  /// VIP values (one per X variable)
  vip : Array[Double]
  /// Orthogonal scores To (n x ort)
  ortho_scores : Array[Array[Double]]
  /// Orthogonal X loadings Po (p x ort)
  ortho_loadings : Array[Array[Double]]
  /// Orthogonal X weights Wo (p x ort, unit norm)
  ortho_weights : Array[Array[Double]]
  /// Orthogonal Y loadings Co (q x ort)
  ortho_y_loadings : Array[Array[Double]]
  /// Orthogonal VIP (one per X variable, OPLS only)
  ortho_vip : Array[Double]
  /// Per-row R2X: PLS has `pre` rows; OPLS has p1, o1..oort rows plus a
  /// final "sum" row (NaN per-component value, cumulative totals in *_cum)
  r2x : Array[Double]
  r2x_cum : Array[Double]
  r2y : Array[Double]
  r2y_cum : Array[Double]
  q2 : Array[Double]
  q2_cum : Array[Double]
  /// "R1" (significant), "NS" (Q2 below threshold) or "N4" (R2Y below 1%)
  significance : Array[String]
  total_r2x : Double
  total_r2y : Double
  total_q2 : Double
  rmsee : Double
  /// Fitted response in the original Y scale (n x q)
  y_pred : Array[Array[Double]]
  /// Fitted class labels (DA models only, otherwise empty)
  fitted_classes : Array[String]
  /// Sorted class levels (DA models only, otherwise empty)
  classes : Array[String]
  x_mean : Array[Double]
  x_sd : Array[Double]
  y_mean : Array[Double]
  y_sd : Array[Double]
  /// Indices of removed near-zero variance variables
  dropped : Array[Int]
}

///|
/// Result of response permutation testing.
pub struct RoplsPermResult {
  nperm : Int
  /// Observed model cumulative R2Y
  obs_r2y : Double
  /// Observed model cumulative Q2
  obs_q2 : Double
  /// Cumulative R2Y of each permuted model
  perm_r2y : Array[Double]
  /// Cumulative Q2 of each permuted model
  perm_q2 : Array[Double]
  /// Similarity with the original response (label agreement or correlation)
  sim : Array[Double]
  p_r2y : Double
  p_q2 : Double
}

// ============================================================================
// Small matrix / vector helpers (rop_ prefix: single src package)
// ============================================================================

///|
fn rop_dot(a : Array[Double], b : Array[Double]) -> Double {
  let mut s = 0.0
  for i = 0; i < a.length(); i = i + 1 {
    s = s + a[i] * b[i]
  }
  s
}

///|
fn rop_norm(v : Array[Double]) -> Double {
  rop_dot(v, v).sqrt()
}

///|
fn rop_copy_mat(m : Array[Array[Double]]) -> Array[Array[Double]] {
  m.map(fn(row) { row.copy() })
}

///|
fn[T] rop_take(a : Array[T], k : Int) -> Array[T] {
  let out : Array[T] = []
  for i in 0..<k {
    out.push(a[i])
  }
  out
}

///|
fn rop_mat_vec(m : Array[Array[Double]], v : Array[Double]) -> Array[Double] {
  let n = m.length()
  let p = v.length()
  let out = Array::make(n, 0.0)
  for i = 0; i < n; i = i + 1 {
    let mut s = 0.0
    for j = 0; j < p; j = j + 1 {
      s = s + m[i][j] * v[j]
    }
    out[i] = s
  }
  out
}

///|
fn rop_mat_t_vec(m : Array[Array[Double]], v : Array[Double]) -> Array[Double] {
  let n = m.length()
  let p = if n == 0 { 0 } else { m[0].length() }
  let out = Array::make(p, 0.0)
  for i = 0; i < n; i = i + 1 {
    for j = 0; j < p; j = j + 1 {
      out[j] = out[j] + m[i][j] * v[i]
    }
  }
  out
}

///|
/// m -= s * l' (in place)
fn rop_deflate(
  m : Array[Array[Double]],
  s : Array[Double],
  l : Array[Double],
) -> Unit {
  for i = 0; i < m.length(); i = i + 1 {
    for j = 0; j < l.length(); j = j + 1 {
      m[i][j] = m[i][j] - s[i] * l[j]
    }
  }
}

///|
/// Squared norm of the outer product s l', i.e. (s's)(l'l).
fn rop_outer_ss(s : Array[Double], l : Array[Double]) -> Double {
  rop_dot(s, s) * rop_dot(l, l)
}

///|
/// Invert a small square matrix by Gauss-Jordan elimination.
fn rop_invert(a : Array[Array[Double]]) -> Array[Array[Double]] raise {
  let k = a.length()
  let m : Array[Array[Double]] = Array::make(k, [])
  for i = 0; i < k; i = i + 1 {
    m[i] = Array::make(2 * k, 0.0)
    for j = 0; j < k; j = j + 1 {
      m[i][j] = a[i][j]
    }
    m[i][k + i] = 1.0
  }
  for col = 0; col < k; col = col + 1 {
    let mut piv = col
    let mut best = m[col][col].abs()
    for r = col + 1; r < k; r = r + 1 {
      if m[r][col].abs() > best {
        best = m[r][col].abs()
        piv = r
      }
    }
    if best < 1.0e-14 {
      raise SeqError("ropls: singular matrix in weight star computation")
    }
    if piv != col {
      let tmp = m[col]
      m[col] = m[piv]
      m[piv] = tmp
    }
    let pivval = m[col][col]
    for j = 0; j < 2 * k; j = j + 1 {
      m[col][j] = m[col][j] / pivval
    }
    for r = 0; r < k; r = r + 1 {
      if r != col {
        let f = m[r][col]
        if f != 0.0 {
          for j = 0; j < 2 * k; j = j + 1 {
            m[r][j] = m[r][j] - f * m[col][j]
          }
        }
      }
    }
  }
  let out : Array[Array[Double]] = Array::make(k, [])
  for i = 0; i < k; i = i + 1 {
    out[i] = Array::make(k, 0.0)
    for j = 0; j < k; j = j + 1 {
      out[i][j] = m[i][k + j]
    }
  }
  out
}

///|
/// Eigenvectors and eigenvalues (descending) of a small symmetric matrix
/// using cyclic Jacobi rotations.
fn rop_jacobi_eigen(
  aa : Array[Array[Double]],
) -> (Array[Double], Array[Array[Double]]) {
  let k = aa.length()
  let a = rop_copy_mat(aa)
  let v : Array[Array[Double]] = Array::make(k, [])
  for i = 0; i < k; i = i + 1 {
    v[i] = Array::make(k, 0.0)
    v[i][i] = 1.0
  }
  for _sweep in 0..<100 {
    let mut off = 0.0
    for p = 0; p < k; p = p + 1 {
      for q = p + 1; q < k; q = q + 1 {
        off = off + a[p][q] * a[p][q]
      }
    }
    if off < 1.0e-30 {
      break
    }
    for p = 0; p < k; p = p + 1 {
      for q = p + 1; q < k; q = q + 1 {
        if a[p][q].abs() > 1.0e-15 {
          let tau = (a[q][q] - a[p][p]) / (2.0 * a[p][q])
          let t = if tau >= 0.0 {
            1.0 / (tau + (1.0 + tau * tau).sqrt())
          } else {
            -1.0 / (-tau + (1.0 + tau * tau).sqrt())
          }
          let c = 1.0 / (1.0 + t * t).sqrt()
          let s = t * c
          for i = 0; i < k; i = i + 1 {
            let aip = a[i][p]
            let aiq = a[i][q]
            a[i][p] = c * aip - s * aiq
            a[i][q] = s * aip + c * aiq
          }
          for i = 0; i < k; i = i + 1 {
            let api = a[p][i]
            let aqi = a[q][i]
            a[p][i] = c * api - s * aqi
            a[q][i] = s * api + c * aqi
          }
          for i = 0; i < k; i = i + 1 {
            let vip = v[i][p]
            let viq = v[i][q]
            v[i][p] = c * vip - s * viq
            v[i][q] = s * vip + c * viq
          }
        }
      }
    }
  }
  let d = Array::make(k, 0.0)
  for i = 0; i < k; i = i + 1 {
    d[i] = a[i][i]
  }
  let idx : Array[Int] = []
  for i in 0..<k {
    idx.push(i)
  }
  idx.sort_by(fn(a : Int, b : Int) -> Int { d[b].compare(d[a]) })
  let dsorted = idx.map(fn(i) { d[i] })
  let vsorted : Array[Array[Double]] = Array::make(k, [])
  for i = 0; i < k; i = i + 1 {
    vsorted[i] = Array::make(k, 0.0)
    for j = 0; j < k; j = j + 1 {
      vsorted[i][j] = v[i][idx[j]]
    }
  }
  (dsorted, vsorted)
}

///|
fn rop_mean(v : Array[Double]) -> Double {
  let mut s = 0.0
  for x in v {
    s = s + x
  }
  s / v.length().to_double()
}

///|
/// Sample variance (denominator n-1).
fn rop_var(v : Array[Double]) -> Double {
  let n = v.length()
  if n < 2 {
    return 0.0
  }
  let m = rop_mean(v)
  let mut s = 0.0
  for x in v {
    s = s + (x - m) * (x - m)
  }
  s / (n - 1).to_double()
}

///|
fn rop_cor(a : Array[Double], b : Array[Double]) -> Double {
  let ma = rop_mean(a)
  let mb = rop_mean(b)
  let mut sa = 0.0
  let mut sb = 0.0
  let mut sc = 0.0
  for i = 0; i < a.length(); i = i + 1 {
    let da = a[i] - ma
    let db = b[i] - mb
    sa = sa + da * da
    sb = sb + db * db
    sc = sc + da * db
  }
  let den = (sa * sb).sqrt()
  if den == 0.0 {
    return 0.0
  }
  sc / den
}

///|
fn rop_check_scale(scale : String) -> Unit raise {
  if scale != "none" &&
    scale != "center" &&
    scale != "pareto" &&
    scale != "standard" {
    raise SeqError("ropls: scale must be one of none, center, pareto, standard")
  }
}

///|
fn rop_check_rect(x : Array[Array[Double]]) -> (Int, Int) raise {
  let n = x.length()
  if n == 0 {
    raise SeqError("ropls: empty X matrix")
  }
  let p = x[0].length()
  if p == 0 {
    raise SeqError("ropls: X matrix has no variables")
  }
  for row in x {
    if row.length() != p {
      raise SeqError("ropls: ragged X matrix")
    }
  }
  (n, p)
}

///|
/// Column means and scaling divisors (sample sd / sqrt sd / 1).
fn rop_scale_params(
  x : Array[Array[Double]],
  scale : String,
) -> (Array[Double], Array[Double]) {
  let n = x.length()
  let p = x[0].length()
  let mean = Array::make(p, 0.0)
  let sd = Array::make(p, 1.0)
  if scale != "none" {
    for j = 0; j < p; j = j + 1 {
      let mut s = 0.0
      for i = 0; i < n; i = i + 1 {
        s = s + x[i][j]
      }
      mean[j] = s / n.to_double()
    }
  }
  if scale == "pareto" || scale == "standard" {
    for j = 0; j < p; j = j + 1 {
      let col = Array::make(n, 0.0)
      for i = 0; i < n; i = i + 1 {
        col[i] = x[i][j]
      }
      let v = rop_var(col)
      let sdv = v.sqrt()
      sd[j] = if scale == "pareto" { sdv.sqrt() } else { sdv }
    }
  }
  (mean, sd)
}

///|
fn rop_apply_scale(
  x : Array[Array[Double]],
  mean : Array[Double],
  sd : Array[Double],
) -> Array[Array[Double]] {
  let n = x.length()
  let p = x[0].length()
  let out : Array[Array[Double]] = Array::make(n, [])
  for i = 0; i < n; i = i + 1 {
    out[i] = Array::make(p, 0.0)
    for j = 0; j < p; j = j + 1 {
      out[i][j] = (x[i][j] - mean[j]) / sd[j]
    }
  }
  out
}

///|
fn rop_sum_sq(m : Array[Array[Double]]) -> Double {
  let mut s = 0.0
  for row in m {
    for v in row {
      s = s + v * v
    }
  }
  s
}

///|
/// Indices of variables with variance below the ropls epsilon threshold.
fn rop_zero_var(x : Array[Array[Double]]) -> Array[Int] {
  let n = x.length()
  let p = x[0].length()
  let eps = 2.220446049250313e-16
  let dropped : Array[Int] = []
  for j = 0; j < p; j = j + 1 {
    let col = Array::make(n, 0.0)
    for i = 0; i < n; i = i + 1 {
      col[i] = x[i][j]
    }
    if rop_var(col) < eps {
      dropped.push(j)
    }
  }
  dropped
}

///|
fn rop_select_columns(
  x : Array[Array[Double]],
  keep : Array[Int],
) -> Array[Array[Double]] {
  let out : Array[Array[Double]] = Array::make(x.length(), [])
  for i = 0; i < x.length(); i = i + 1 {
    out[i] = Array::make(keep.length(), 0.0)
    for j = 0; j < keep.length(); j = j + 1 {
      out[i][j] = x[i][keep[j]]
    }
  }
  out
}

// ============================================================================
// NIPALS PCA
// ============================================================================

///|
/// NIPALS PCA as performed by ropls.
///
/// - `ncomp` number of components; `0` selects automatically (first count
///   with cumulative R2X above 50%, capped at `min(10, n, p)`).
/// - `scale` one of "standard" (default, mean 0 / variance 1), "pareto",
///   "center" or "none".
pub fn ropls_pca(
  x : Array[Array[Double]],
  ncomp? : Int = 0,
  scale? : String = "standard",
) -> RoplsPcaResult raise {
  rop_check_scale(scale)
  let (n0, p0) = rop_check_rect(x)
  if n0 < 2 {
    raise SeqError("ropls: PCA requires at least 2 samples")
  }
  let dropped = rop_zero_var(x)
  let pkept = p0 - dropped.length()
  if pkept == 0 {
    raise SeqError("ropls: all X variables have near-zero variance")
  }
  let keep : Array[Int] = []
  let drop_set = dropped.copy()
  for j = 0; j < p0; j = j + 1 {
    if not(drop_set.contains(j)) {
      keep.push(j)
    }
  }
  let xf = rop_select_columns(x, keep)
  let (mean0, sd0) = rop_scale_params(xf, scale)
  // Re-order means/sd back onto the original variable indexing is not needed:
  // results only expose kept variables, matching ropls post-filter matrices.
  let xm = rop_apply_scale(xf, mean0, sd0)
  let n = n0
  let p = pkept
  let ssx = rop_sum_sq(xm)
  let cap = if n < p {
    if n < 10 {
      n
    } else {
      10
    }
  } else if p < 10 {
    p
  } else {
    10
  }
  let autol = ncomp == 0
  let maxc = if autol { cap } else { ncomp }
  if maxc < 1 || maxc > cap {
    raise SeqError(
      "ropls: ncomp must be between 1 and min(10, n, p) (or 0 for automatic)",
    )
  }
  let tmat : Array[Array[Double]] = Array::make(n, [])
  let pmat : Array[Array[Double]] = Array::make(p, [])
  for i = 0; i < n; i = i + 1 {
    tmat[i] = Array::make(maxc, 0.0)
  }
  for j = 0; j < p; j = j + 1 {
    pmat[j] = Array::make(maxc, 0.0)
  }
  let varv : Array[Double] = []
  let r2v : Array[Double] = []
  let cumv : Array[Double] = []
  let iters : Array[Int] = []
  let mut selected = maxc
  for h in 0..<maxc {
    let mut told = xm.map(fn(row) { row[0] })
    let mut pold = Array::make(p, 0.0)
    let mut ite = 1
    let mut pnew : Array[Double] = []
    let mut tnew : Array[Double] = []
    while true {
      pnew = rop_mat_t_vec(xm, told)
      let tt = rop_dot(told, told)
      if tt > 0.0 {
        for j = 0; j < p; j = j + 1 {
          pnew[j] = pnew[j] / tt
        }
      }
      let pn = rop_norm(pnew)
      if pn > 0.0 {
        for j = 0; j < p; j = j + 1 {
          pnew[j] = pnew[j] / pn
        }
      }
      tnew = rop_mat_vec(xm, pnew)
      let mut diff = 0.0
      for j = 0; j < p; j = j + 1 {
        diff = diff + (pnew[j] - pold[j]) * (pnew[j] - pold[j])
      }
      if diff.sqrt() < 1.0e-6 || ite > 100 {
        break
      }
      told = tnew
      pold = pnew
      ite = ite + 1
    }
    let itef = ite
    for i = 0; i < n; i = i + 1 {
      tmat[i][h] = tnew[i]
    }
    for j = 0; j < p; j = j + 1 {
      pmat[j][h] = pnew[j]
    }
    varv.push(rop_dot(tnew, tnew) / (n - 1).to_double())
    let r2 = rop_outer_ss(tnew, pnew) / ssx
    r2v.push(r2)
    cumv.push(
      r2v.fold(init=0.0, fn(acc : Double, v : Double) -> Double { acc + v }),
    )
    iters.push(itef)
    rop_deflate(xm, tnew, pnew)
    if autol && cumv[h] > 0.5 {
      selected = h + 1
      break
    }
  }
  if autol && cumv[cumv.length() - 1] <= 0.5 {
    selected = maxc
  }
  let ts : Array[Array[Double]] = Array::make(n, [])
  let ps : Array[Array[Double]] = Array::make(p, [])
  for i = 0; i < n; i = i + 1 {
    ts[i] = rop_take(tmat[i], selected)
  }
  for j = 0; j < p; j = j + 1 {
    ps[j] = rop_take(pmat[j], selected)
  }
  RoplsPcaResult::{
    n,
    p,
    ncomp: selected,
    scale,
    scores: ts,
    loadings: ps,
    variance: rop_take(varv, selected),
    r2x: rop_take(r2v, selected),
    r2x_cum: rop_take(cumv, selected),
    iters: rop_take(iters, selected),
    x_mean: mean0,
    x_sd: sd0,
    dropped,
  }
}

// ============================================================================
// NIPALS (O)PLS machinery
// ============================================================================

///|
struct RoplsPrepared {
  n : Int
  p : Int
  q : Int
  xs : Array[Array[Double]]
  ys : Array[Array[Double]]
  xmean : Array[Double]
  xsd : Array[Double]
  ymean : Array[Double]
  ysd : Array[Double]
  ssx : Double
  ssy : Double
  dropped : Array[Int]
  classes : Array[String]
  da : Bool
}

///|
/// Lexicographic byte order (R's character ordering). MoonBit's built-in
/// String comparison is shortlex (length first), which differs from R.
fn rop_cmp_str(a : String, b : String) -> Int {
  let ca = a.to_array()
  let cb = b.to_array()
  let k = if ca.length() < cb.length() { ca.length() } else { cb.length() }
  for i in 0..<k {
    let oa = ca[i].to_int()
    let ob = cb[i].to_int()
    if oa < ob {
      return -1
    }
    if oa > ob {
      return 1
    }
  }
  ca.length().compare(cb.length())
}

///|
fn rop_dummy(
  labels : Array[String],
) -> (Array[String], Array[Array[Double]]) raise {
  let levels : Array[String] = []
  for lab in labels {
    if not(levels.contains(lab)) {
      levels.push(lab)
    }
  }
  levels.sort_by(fn(a : String, b : String) -> Int { rop_cmp_str(a, b) })
  if levels.length() < 2 {
    raise SeqError("ropls: classification requires at least 2 classes")
  }
  let n = labels.length()
  if levels.length() == 2 {
    let d : Array[Array[Double]] = Array::make(n, [])
    for i = 0; i < n; i = i + 1 {
      d[i] = [if labels[i] == levels[1] { 1.0 } else { 0.0 }]
    }
    (levels, d)
  } else {
    let k = levels.length()
    let d : Array[Array[Double]] = Array::make(n, [])
    for i = 0; i < n; i = i + 1 {
      d[i] = Array::make(k, 0.0)
      let mut found = false
      for j = 0; j < k; j = j + 1 {
        if labels[i] == levels[j] {
          d[i][j] = 1.0
          found = true
        }
      }
      if not(found) {
        raise SeqError("ropls: internal dummy coding error")
      }
    }
    (levels, d)
  }
}

///|
fn rop_prepare(
  xraw : Array[Array[Double]],
  resp : RoplsY,
  scale : String,
) -> RoplsPrepared raise {
  rop_check_scale(scale)
  let (n0, p0) = rop_check_rect(xraw)
  if n0 < 2 {
    raise SeqError("ropls: at least 2 samples required")
  }
  let (labels0, yraw, da) = match resp {
    RoplsY::Continuous(y) => {
      if y.length() != n0 {
        raise SeqError("ropls: X and Y sample counts differ")
      }
      if y.length() > 0 && y[0].length() == 0 {
        raise SeqError("ropls: empty Y response")
      }
      for row in y {
        if row.length() != y[0].length() {
          raise SeqError("ropls: ragged Y matrix")
        }
      }
      ([], y, false)
    }
    RoplsY::Classes(lab) => {
      if lab.length() != n0 {
        raise SeqError("ropls: X rows and labels differ")
      }
      let (lev, d) = rop_dummy(lab)
      (lev, d, true)
    }
  }
  let q0 = yraw[0].length()
  let dropped = rop_zero_var(xraw)
  let pkept = p0 - dropped.length()
  if pkept == 0 {
    raise SeqError("ropls: all X variables have near-zero variance")
  }
  let keep : Array[Int] = []
  for j = 0; j < p0; j = j + 1 {
    if not(dropped.contains(j)) {
      keep.push(j)
    }
  }
  let xf = rop_select_columns(xraw, keep)
  let (xmean, xsd) = rop_scale_params(xf, scale)
  let xs = rop_apply_scale(xf, xmean, xsd)
  // DA responses are always centered and scaled (as ropls does for
  // character responses); continuous responses follow the scale argument.
  let (ymean, ysd) = if da {
    rop_scale_params(yraw, "standard")
  } else {
    rop_scale_params(yraw, scale)
  }
  let ys = rop_apply_scale(yraw, ymean, ysd)
  let ssx = rop_sum_sq(xs)
  let ssy = rop_sum_sq(ys)
  RoplsPrepared::{
    n: n0,
    p: pkept,
    q: q0,
    xs,
    ys,
    xmean,
    xsd,
    ymean,
    ysd,
    ssx,
    ssy,
    dropped,
    classes: labels0,
    da,
  }
}

///|
/// One NIPALS PLS component on (deflated) X, Y.
/// Returns (w, t, c, u, p, iters). For a single response u is returned as
/// the initial y column, matching ropls.
fn rop_pls_component(
  x : Array[Array[Double]],
  y : Array[Array[Double]],
  tol : Double,
) -> (
  Array[Double],
  Array[Double],
  Array[Double],
  Array[Double],
  Array[Double],
  Int,
) {
  let n = x.length()
  let p = x[0].length()
  let q = y[0].length()
  let mut u = Array::make(n, 0.0)
  for i = 0; i < n; i = i + 1 {
    u[i] = y[i][0]
  }
  let mut told = Array::make(n, 0.0)
  let mut ite = 1
  let mut w : Array[Double] = []
  let mut t : Array[Double] = []
  let mut c : Array[Double] = []
  let mut maxg = 0
  while maxg < 100000 {
    maxg = maxg + 1
    let uu = rop_dot(u, u)
    w = rop_mat_t_vec(x, u)
    if uu > 0.0 {
      for j = 0; j < p; j = j + 1 {
        w[j] = w[j] / uu
      }
    }
    let wn = rop_norm(w)
    if wn > 0.0 {
      for j = 0; j < p; j = j + 1 {
        w[j] = w[j] / wn
      }
    }
    t = rop_mat_vec(x, w)
    let tt = rop_dot(t, t)
    c = rop_mat_t_vec(y, t)
    if tt > 0.0 {
      for k = 0; k < q; k = k + 1 {
        c[k] = c[k] / tt
      }
    }
    if q == 1 {
      break
    }
    let cc = rop_dot(c, c)
    let unew = rop_mat_vec(y, c)
    if cc > 0.0 {
      for i = 0; i < n; i = i + 1 {
        unew[i] = unew[i] / cc
      }
    }
    let mut d = 0.0
    for i = 0; i < n; i = i + 1 {
      let r = if t[i] != 0.0 { (told[i] - t[i]) / t[i] } else { 0.0 }
      d = d + r * r
    }
    // R breaks here and keeps the u that produced the converged t; only
    // iterate u forward when another round is required.
    if d.sqrt() < tol {
      break
    }
    u = unew
    told = t
    ite = ite + 1
  }
  let tt = rop_dot(t, t)
  let pv = rop_mat_t_vec(x, t)
  if tt > 0.0 {
    for j = 0; j < p; j = j + 1 {
      pv[j] = pv[j] / tt
    }
  }
  (w, t, c, u, pv, ite)
}

///|
/// Cross-validated PRESS for one component on the current (deflated) data.
fn rop_pls_press(
  xn : Array[Array[Double]],
  yn : Array[Array[Double]],
  folds : Array[Array[Int]],
) -> Double {
  let p = xn[0].length()
  let mut press = 0.0
  for fold in folds {
    let tr : Array[Int] = []
    let testset = fold
    for i = 0; i < xn.length(); i = i + 1 {
      if not(testset.contains(i)) {
        tr.push(i)
      }
    }
    let xtr : Array[Array[Double]] = tr.map(fn(i) { xn[i].copy() })
    let ytr : Array[Array[Double]] = tr.map(fn(i) { yn[i].copy() })
    let (w, _, c, _, _, _) = rop_pls_component(xtr, ytr, 1.0e-6)
    for i in testset {
      let mut xw = 0.0
      for j = 0; j < p; j = j + 1 {
        xw = xw + xn[i][j] * w[j]
      }
      for k = 0; k < c.length(); k = k + 1 {
        let e = yn[i][k] - xw * c[k]
        press = press + e * e
      }
    }
  }
  press
}

///|
fn rop_make_folds(n : Int, k : Int) -> Array[Array[Int]] raise {
  if k < 2 || k > n {
    raise SeqError(
      "ropls: crossval folds must be between 2 and the sample count",
    )
  }
  let folds : Array[Array[Int]] = Array::make(k, [])
  for i = 0; i < k; i = i + 1 {
    folds[i] = []
  }
  for i = 0; i < n; i = i + 1 {
    folds[i % k].push(i)
  }
  folds
}

///|
fn rop_signif(r2y_a : Double, q2_a : Double, threshold : Double) -> String {
  if r2y_a < 0.01 {
    "N4"
  } else if q2_a < threshold {
    "NS"
  } else {
    "R1"
  }
}

///|
/// Column-oriented storage of component matrices (p x h / n x h).
struct RoplsPlsFit {
  w : Array[Array[Double]]
  t : Array[Array[Double]]
  p : Array[Array[Double]]
  c : Array[Array[Double]]
  u : Array[Array[Double]]
  r2x : Array[Double]
  r2y : Array[Double]
  q2 : Array[Double]
  sig : Array[String]
  ncomp : Int
  iters : Array[Int]
}

///|
fn rop_fit_pls(
  pre0 : RoplsPrepared,
  ncomp : Int,
  crossval : Int,
  autol : Bool,
) -> RoplsPlsFit raise {
  let n = pre0.n
  let p = pre0.p
  let q = pre0.q
  let cap = if n < p {
    if n < 10 {
      n
    } else {
      10
    }
  } else if p < 10 {
    p
  } else {
    10
  }
  let maxc = if autol { cap } else { ncomp }
  if !autol && (ncomp < 1 || ncomp > cap) {
    raise SeqError("ropls: ncomp must be between 1 and min(10, n, p)")
  }
  let folds = rop_make_folds(n, crossval)
  let xn = rop_copy_mat(pre0.xs)
  let yn = rop_copy_mat(pre0.ys)
  let w : Array[Array[Double]] = Array::make(p, [])
  let pm : Array[Array[Double]] = Array::make(p, [])
  let t : Array[Array[Double]] = Array::make(n, [])
  let u : Array[Array[Double]] = Array::make(n, [])
  let c : Array[Array[Double]] = Array::make(q, [])
  for j = 0; j < p; j = j + 1 {
    w[j] = Array::make(maxc, 0.0)
    pm[j] = Array::make(maxc, 0.0)
  }
  for i = 0; i < n; i = i + 1 {
    t[i] = Array::make(maxc, 0.0)
    u[i] = Array::make(maxc, 0.0)
  }
  for k = 0; k < q; k = k + 1 {
    c[k] = Array::make(maxc, 0.0)
  }
  let r2x : Array[Double] = []
  let r2y : Array[Double] = []
  let q2v : Array[Double] = []
  let sig : Array[String] = []
  let iters : Array[Int] = []
  let threshold = if n > 100 { 0.0 } else { 0.05 }
  let mut rss = pre0.ssy
  let mut h = 0
  let mut selected = maxc
  let mut stopped = false
  while h < maxc {
    let (wv, tv, cv, uv, pv, ite) = rop_pls_component(xn, yn, 1.0e-6)
    let press = rop_pls_press(xn, yn, folds)
    let q2_a = 1.0 - press / rss
    let r2x_a = rop_outer_ss(tv, pv) / pre0.ssx
    let r2y_a = rop_outer_ss(tv, cv) / pre0.ssy
    let s = rop_signif(r2y_a, q2_a, threshold)
    for i = 0; i < n; i = i + 1 {
      t[i][h] = tv[i]
      u[i][h] = uv[i]
    }
    for j = 0; j < p; j = j + 1 {
      w[j][h] = wv[j]
      pm[j][h] = pv[j]
    }
    for k = 0; k < q; k = k + 1 {
      c[k][h] = cv[k]
    }
    r2x.push(r2x_a)
    r2y.push(r2y_a)
    q2v.push(q2_a)
    sig.push(s)
    iters.push(ite)
    if autol && s != "R1" {
      selected = h
      stopped = true
      break
    }
    // RSS of the current component fit, denominator for the next Q2
    rss = 0.0
    for i = 0; i < n; i = i + 1 {
      for k = 0; k < q; k = k + 1 {
        let e = yn[i][k] - tv[i] * cv[k]
        rss = rss + e * e
      }
    }
    rop_deflate(xn, tv, pv)
    rop_deflate(yn, tv, cv)
    h = h + 1
  }
  if autol && selected == 0 {
    raise SeqError(
      "ropls: no model built: the first predictive component is not significant",
    )
  }
  if autol && !stopped {
    selected = maxc
  }
  let trim = if autol { selected } else { ncomp }
  let ws : Array[Array[Double]] = Array::make(p, [])
  let ps : Array[Array[Double]] = Array::make(p, [])
  let ts : Array[Array[Double]] = Array::make(n, [])
  let us : Array[Array[Double]] = Array::make(n, [])
  let cs : Array[Array[Double]] = Array::make(q, [])
  for j = 0; j < p; j = j + 1 {
    ws[j] = rop_take(w[j], trim)
    ps[j] = rop_take(pm[j], trim)
  }
  for i = 0; i < n; i = i + 1 {
    ts[i] = rop_take(t[i], trim)
    us[i] = rop_take(u[i], trim)
  }
  for k = 0; k < q; k = k + 1 {
    cs[k] = rop_take(c[k], trim)
  }
  RoplsPlsFit::{
    w: ws,
    t: ts,
    p: ps,
    c: cs,
    u: us,
    r2x: rop_take(r2x, trim),
    r2y: rop_take(r2y, trim),
    q2: rop_take(q2v, trim),
    sig: rop_take(sig, trim),
    ncomp: trim,
    iters: rop_take(iters, trim),
  }
}

///|
/// Weight star: r = w for 1 component, r = w * solve(p' w) otherwise.
fn rop_weight_star(
  w : Array[Array[Double]],
  pm : Array[Array[Double]],
) -> Array[Array[Double]] raise {
  let p = w.length()
  let a = if p == 0 { 0 } else { w[0].length() }
  if a == 1 {
    let out : Array[Array[Double]] = Array::make(p, [])
    for j = 0; j < p; j = j + 1 {
      out[j] = [w[j][0]]
    }
    return out
  }
  let pw : Array[Array[Double]] = Array::make(a, [])
  for h1 = 0; h1 < a; h1 = h1 + 1 {
    pw[h1] = Array::make(a, 0.0)
    for h2 = 0; h2 < a; h2 = h2 + 1 {
      let mut s = 0.0
      for j = 0; j < p; j = j + 1 {
        s = s + pm[j][h1] * w[j][h2]
      }
      pw[h1][h2] = s
    }
  }
  let inv = rop_invert(pw)
  let out : Array[Array[Double]] = Array::make(p, [])
  for j = 0; j < p; j = j + 1 {
    out[j] = Array::make(a, 0.0)
    for h2 = 0; h2 < a; h2 = h2 + 1 {
      let mut s = 0.0
      for h1 = 0; h1 < a; h1 = h1 + 1 {
        s = s + w[j][h1] * inv[h1][h2]
      }
      out[j][h2] = s
    }
  }
  out
}

///|
/// b = r c' (p x q)
fn rop_coefficients(
  r : Array[Array[Double]],
  c : Array[Array[Double]],
) -> Array[Array[Double]] {
  let p = r.length()
  let a = if p == 0 { 0 } else { r[0].length() }
  let q = c.length()
  let b : Array[Array[Double]] = Array::make(p, [])
  for j = 0; j < p; j = j + 1 {
    b[j] = Array::make(q, 0.0)
    for k = 0; k < q; k = k + 1 {
      let mut s = 0.0
      for h = 0; h < a; h = h + 1 {
        s = s + r[j][h] * c[k][h]
      }
      b[j][k] = s
    }
  }
  b
}

///|
/// Fitted response back in the original Y scale.
fn rop_descale_y(
  ysfit : Array[Array[Double]],
  ymean : Array[Double],
  ysd : Array[Double],
) -> Array[Array[Double]] {
  let n = ysfit.length()
  let q = if n == 0 { 0 } else { ysfit[0].length() }
  let out : Array[Array[Double]] = Array::make(n, [])
  for i = 0; i < n; i = i + 1 {
    out[i] = Array::make(q, 0.0)
    for k = 0; k < q; k = k + 1 {
      out[i][k] = ysfit[i][k] * ysd[k] + ymean[k]
    }
  }
  out
}

///|
/// yfit = t c'
fn rop_tc(
  t : Array[Array[Double]],
  c : Array[Array[Double]],
) -> Array[Array[Double]] {
  let n = t.length()
  let q = c.length()
  let out : Array[Array[Double]] = Array::make(n, [])
  for i = 0; i < n; i = i + 1 {
    out[i] = Array::make(q, 0.0)
    for k = 0; k < q; k = k + 1 {
      let mut s = 0.0
      for h = 0; h < t[i].length(); h = h + 1 {
        s = s + t[i][h] * c[k][h]
      }
      out[i][k] = s
    }
  }
  out
}

///|
fn rop_rmse(a : Array[Array[Double]], b : Array[Array[Double]]) -> Double {
  let mut s = 0.0
  let mut cnt = 0
  for i = 0; i < a.length(); i = i + 1 {
    for k = 0; k < a[i].length(); k = k + 1 {
      let e = a[i][k] - b[i][k]
      s = s + e * e
      cnt = cnt + 1
    }
  }
  (s / cnt.to_double()).sqrt()
}

///|
fn rop_pls_vip(
  w : Array[Array[Double]],
  t : Array[Array[Double]],
  c : Array[Array[Double]],
) -> Array[Double] {
  let p = w.length()
  let a = if p == 0 { 0 } else { w[0].length() }
  let ssy : Array[Double] = []
  let mut total = 0.0
  for h = 0; h < a; h = h + 1 {
    let mut ss = 0.0
    for i = 0; i < t.length(); i = i + 1 {
      for k = 0; k < c.length(); k = k + 1 {
        let v = t[i][h] * c[k][h]
        ss = ss + v * v
      }
    }
    ssy.push(ss)
    total = total + ss
  }
  let vip = Array::make(p, 0.0)
  if total > 0.0 {
    for j = 0; j < p; j = j + 1 {
      let mut num = 0.0
      for h = 0; h < a; h = h + 1 {
        num = num + w[j][h] * w[j][h] * ssy[h]
      }
      vip[j] = (p.to_double() * num / total).sqrt()
    }
  }
  vip
}

// ============================================================================
// OPLS
// ============================================================================

///|
/// One predictive NIPALS pass on a training partition for OPLS.
/// Returns (w, t, c, u, p).
fn rop_opls_predictive(
  x : Array[Array[Double]],
  y : Array[Array[Double]],
) -> (Array[Double], Array[Double], Array[Double], Array[Double], Array[Double]) {
  let n = x.length()
  let p = x[0].length()
  let q = y[0].length()
  let mut u = Array::make(n, 0.0)
  for i = 0; i < n; i = i + 1 {
    u[i] = y[i][0]
  }
  let mut w : Array[Double] = []
  let mut t : Array[Double] = []
  let mut c : Array[Double] = []
  let mut maxg = 0
  while maxg < 100000 {
    maxg = maxg + 1
    let uu = rop_dot(u, u)
    w = rop_mat_t_vec(x, u)
    if uu > 0.0 {
      for j = 0; j < p; j = j + 1 {
        w[j] = w[j] / uu
      }
    }
    let wn = rop_norm(w)
    if wn > 0.0 {
      for j = 0; j < p; j = j + 1 {
        w[j] = w[j] / wn
      }
    }
    t = rop_mat_vec(x, w)
    let tt = rop_dot(t, t)
    c = rop_mat_t_vec(y, t)
    if tt > 0.0 {
      for k = 0; k < q; k = k + 1 {
        c[k] = c[k] / tt
      }
    }
    let cc = rop_dot(c, c)
    let unew = rop_mat_vec(y, c)
    if cc > 0.0 {
      for i = 0; i < n; i = i + 1 {
        unew[i] = unew[i] / cc
      }
    }
    let mut d = 0.0
    for i = 0; i < n; i = i + 1 {
      let r = if unew[i] != 0.0 { (u[i] - unew[i]) / unew[i] } else { 0.0 }
      d = d + r * r
    }
    u = unew
    if q == 1 || d.sqrt() < 1.0e-10 {
      break
    }
  }
  let tt = rop_dot(t, t)
  let pv = rop_mat_t_vec(x, t)
  if tt > 0.0 {
    for j = 0; j < p; j = j + 1 {
      pv[j] = pv[j] / tt
    }
  }
  (w, t, c, u, pv)
}

///|
/// tw = ww V, where columns of V correspond to non-vanishing singular values
/// of ww (p x q). Used to orthogonalize the OPLS weight for multi-response Y.
fn rop_opls_tw(
  x : Array[Array[Double]],
  y : Array[Array[Double]],
) -> Array[Array[Double]] {
  let p = x[0].length()
  let q = y[0].length()
  let eps = 2.220446049250313e-16
  // ww[:,j] = x' y_j / (y_j' y_j)
  let ww : Array[Array[Double]] = Array::make(p, [])
  for j = 0; j < p; j = j + 1 {
    ww[j] = Array::make(q, 0.0)
  }
  for jj = 0; jj < q; jj = jj + 1 {
    let mut yy = 0.0
    let col = Array::make(p, 0.0)
    for i = 0; i < x.length(); i = i + 1 {
      yy = yy + y[i][jj] * y[i][jj]
      for j = 0; j < p; j = j + 1 {
        col[j] = col[j] + x[i][j] * y[i][jj]
      }
    }
    if yy > 0.0 {
      for j = 0; j < p; j = j + 1 {
        ww[j][jj] = col[j] / yy
      }
    }
  }
  // Gram matrix ww' ww (q x q), eigenvalues are the squared singular values
  let gram : Array[Array[Double]] = Array::make(q, [])
  for k = 0; k < q; k = k + 1 {
    gram[k] = Array::make(q, 0.0)
    for l = 0; l < q; l = l + 1 {
      let mut s = 0.0
      for j = 0; j < p; j = j + 1 {
        s = s + ww[j][k] * ww[j][l]
      }
      gram[k][l] = s
    }
  }
  let (d, v) = rop_jacobi_eigen(gram)
  let mut total = 0.0
  for x in d {
    if x > 0.0 {
      total = total + x
    }
  }
  let keep : Array[Int] = []
  for k = 0; k < q; k = k + 1 {
    if d[k] > eps * total {
      keep.push(k)
    }
  }
  let tw : Array[Array[Double]] = Array::make(p, [])
  for j = 0; j < p; j = j + 1 {
    tw[j] = Array::make(keep.length(), 0.0)
    for h = 0; h < keep.length(); h = h + 1 {
      let mut s = 0.0
      for k = 0; k < q; k = k + 1 {
        s = s + ww[j][k] * v[k][keep[h]]
      }
      tw[j][h] = s
    }
  }
  tw
}

///|
struct RoplsOplsFit {
  w : Array[Double]
  t : Array[Double]
  pm : Array[Double]
  c : Array[Double]
  u : Array[Double]
  wo : Array[Array[Double]]
  to : Array[Array[Double]]
  po : Array[Array[Double]]
  co : Array[Array[Double]]
  r2x : Array[Double]
  r2x_cum : Array[Double]
  r2y : Array[Double]
  r2y_cum : Array[Double]
  q2 : Array[Double]
  q2_cum : Array[Double]
  sig : Array[String]
  ort : Int
}

///|
fn rop_fit_opls(
  pre0 : RoplsPrepared,
  ortho_req : Int,
  crossval : Int,
  autol : Bool,
) -> RoplsOplsFit raise {
  let n = pre0.n
  let p = pre0.p
  let q = pre0.q
  let cap = if n < p {
    if n < 10 {
      n
    } else {
      10
    }
  } else if p < 10 {
    p
  } else {
    10
  }
  let ort_max = if autol { cap - 1 } else { ortho_req }
  if !autol && (ortho_req < 1 || 1 + ortho_req > cap) {
    raise SeqError(
      "ropls: ortho must be >= 1 with predI + orthoI <= min(10, n, p)",
    )
  }
  if ort_max < 1 {
    raise SeqError(
      "ropls: OPLS requires room for at least one orthogonal component",
    )
  }
  let folds = rop_make_folds(n, crossval)
  let k = folds.length()
  // Training partitions (k folds + full data as the last item)
  let tra_x : Array[Array[Array[Double]]] = []
  let tra_y : Array[Array[Array[Double]]] = []
  let tes_x : Array[Array[Array[Double]]] = []
  let tes_y : Array[Array[Array[Double]]] = []
  for fold in folds {
    let tr : Array[Int] = []
    for i = 0; i < n; i = i + 1 {
      if not(fold.contains(i)) {
        tr.push(i)
      }
    }
    tra_x.push(tr.map(fn(i) { pre0.xs[i].copy() }))
    tra_y.push(tr.map(fn(i) { pre0.ys[i].copy() }))
    tes_x.push(fold.map(fn(i) { pre0.xs[i].copy() }))
    tes_y.push(fold.map(fn(i) { pre0.ys[i].copy() }))
  }
  tra_x.push(rop_copy_mat(pre0.xs))
  tra_y.push(rop_copy_mat(pre0.ys))
  // Orthogonal component storage (up to ort_max), plus predictive vectors.
  let wo : Array[Array[Double]] = Array::make(p, [])
  let po : Array[Array[Double]] = Array::make(p, [])
  for j = 0; j < p; j = j + 1 {
    wo[j] = Array::make(ort_max, 0.0)
    po[j] = Array::make(ort_max, 0.0)
  }
  let to : Array[Array[Double]] = Array::make(n, [])
  for i = 0; i < n; i = i + 1 {
    to[i] = Array::make(ort_max, 0.0)
  }
  let co : Array[Array[Double]] = Array::make(q, [])
  for kk = 0; kk < q; kk = kk + 1 {
    co[kk] = Array::make(ort_max, 0.0)
  }
  // Per-row metrics, rows: p1, o1..o(ort_max)
  let nrows = ort_max + 1
  let r2x_v = Array::make(nrows, 0.0)
  let r2x_cumv = Array::make(nrows, 0.0)
  let r2y_v = Array::make(nrows, 0.0)
  let r2y_cumv = Array::make(nrows, 0.0)
  let q2_v = Array::make(nrows, 0.0)
  let q2_cumv = Array::make(nrows, 0.0)
  let sigv : Array[String] = Array::make(nrows, "")
  let mut wfinal = Array::make(p, 0.0)
  let mut tfinal = Array::make(n, 0.0)
  let mut pfinal = Array::make(p, 0.0)
  let mut cfinal = Array::make(q, 0.0)
  let mut ufinal = Array::make(n, 0.0)
  let prk = Array::make(k, 0.0)
  let mut no_final = ort_max + 1
  let mut broke = false
  for noN in 1..<(ort_max + 2) {
    for cvN in 1..<(k + 2) {
      let xtr = tra_x[cvN - 1]
      let ytr = tra_y[cvN - 1]
      // For multiple responses, build the SVD-derived tw basis
      let tw = if q > 1 { rop_opls_tw(xtr, ytr) } else { [] }
      let (wv, tv, cv, uv, pv) = rop_opls_predictive(xtr, ytr)
      // Orthogonal weight
      let wov = pv.copy()
      if q > 1 {
        for h = 0; h < tw[0].length(); h = h + 1 {
          let mut num = 0.0
          let mut den = 0.0
          for j = 0; j < p; j = j + 1 {
            num = num + tw[j][h] * wov[j]
            den = den + tw[j][h] * tw[j][h]
          }
          for j = 0; j < p; j = j + 1 {
            wov[j] = wov[j] - num / den * tw[j][h]
          }
        }
      } else {
        let num = rop_dot(wv, pv)
        let den = rop_dot(wv, wv)
        for j = 0; j < p; j = j + 1 {
          wov[j] = wov[j] - num / den * wv[j]
        }
      }
      let won = rop_norm(wov)
      if won > 0.0 {
        for j = 0; j < p; j = j + 1 {
          wov[j] = wov[j] / won
        }
      }
      let tov = rop_mat_vec(xtr, wov)
      let tto = rop_dot(tov, tov)
      let cov = rop_mat_t_vec(ytr, tov)
      if tto > 0.0 {
        for kk = 0; kk < q; kk = kk + 1 {
          cov[kk] = cov[kk] / tto
        }
      }
      let pov = rop_mat_t_vec(xtr, tov)
      if tto > 0.0 {
        for j = 0; j < p; j = j + 1 {
          pov[j] = pov[j] / tto
        }
      }
      if cvN <= k {
        // Cross-validation: PRESS with the predictive component on the
        // progressively deflated held-out partition
        let xte = tes_x[cvN - 1]
        let yte = tes_y[cvN - 1]
        let mut foldpress = 0.0
        for i = 0; i < xte.length(); i = i + 1 {
          let mut xw = 0.0
          for j = 0; j < p; j = j + 1 {
            xw = xw + xte[i][j] * wv[j]
          }
          for kk = 0; kk < q; kk = kk + 1 {
            let e = yte[i][kk] - xw * cv[kk]
            foldpress = foldpress + e * e
          }
        }
        prk[cvN - 1] = foldpress
        let tote = rop_mat_vec(xte, wov)
        rop_deflate(tes_x[cvN - 1], tote, pov)
        if cvN == k {
          let mut press = 0.0
          for v in prk {
            press = press + v
          }
          let q2cum = 1.0 - press / pre0.ssy
          if noN == 1 {
            q2_v[0] = q2cum
            q2_cumv[0] = q2cum
          } else {
            q2_cumv[noN - 1] = q2cum - q2_v[0]
            let mut prev = 0.0
            for rr in 0..<(noN - 1) {
              prev = prev + q2_v[rr]
            }
            q2_v[noN - 1] = q2cum - prev
          }
        }
      } else {
        // Full data pass
        let r2y_cum_stage = rop_outer_ss(tv, cv) / pre0.ssy
        if noN == 1 {
          r2y_v[0] = r2y_cum_stage
          r2y_cumv[0] = r2y_cum_stage
        } else {
          r2y_cumv[noN - 1] = r2y_cum_stage - r2y_v[0]
          let mut prev = 0.0
          for rr in 0..<(noN - 1) {
            prev = prev + r2y_v[rr]
          }
          r2y_v[noN - 1] = r2y_cum_stage - prev
        }
        if noN <= ort_max {
          r2x_v[noN] = rop_outer_ss(tov, pov) / pre0.ssx
          for j = 0; j < p; j = j + 1 {
            wo[j][noN - 1] = wov[j]
            po[j][noN - 1] = pov[j]
          }
          for i = 0; i < n; i = i + 1 {
            to[i][noN - 1] = tov[i]
          }
          for kk = 0; kk < q; kk = kk + 1 {
            co[kk][noN - 1] = cov[kk]
          }
        }
        let s = rop_signif(r2y_v[noN - 1], q2_v[noN - 1], 0.01)
        sigv[noN - 1] = s
        if autol && s != "R1" {
          no_final = noN
          broke = true
          break
        } else {
          wfinal = wv.copy()
          tfinal = tv.copy()
          pfinal = pv.copy()
          cfinal = cv.copy()
          ufinal = uv.copy()
        }
      }
      if noN < ort_max + 1 {
        rop_deflate(tra_x[cvN - 1], tov, pov)
      }
    }
    if broke {
      break
    }
  }
  // Resolve the number of orthogonal components
  let ort = if autol {
    if sigv[0] != "R1" || (nrows > 1 && sigv[1] != "R1") {
      raise SeqError(
        "ropls: no OPLS model built: predictive or first orthogonal component not significant",
      )
    }
    let mut all_r1 = true
    for rr in 0..<nrows {
      if sigv[rr] != "R1" {
        all_r1 = false
      }
    }
    if all_r1 {
      ort_max
    } else {
      no_final - 2
    }
  } else {
    ortho_req
  }
  if ort < 0 {
    raise SeqError("ropls: no significant orthogonal OPLS component")
  }
  // Predictive R2X (computed after the last filtering stage in ropls)
  r2x_v[0] = rop_outer_ss(tfinal, pfinal) / pre0.ssx
  r2x_cumv[0] = r2x_v[0]
  for rr in 1..<(ort + 1) {
    r2x_cumv[rr] = r2x_cumv[rr - 1] + r2x_v[rr]
  }
  let wos : Array[Array[Double]] = Array::make(p, [])
  let pos : Array[Array[Double]] = Array::make(p, [])
  let tos : Array[Array[Double]] = Array::make(n, [])
  let cos : Array[Array[Double]] = Array::make(q, [])
  for j = 0; j < p; j = j + 1 {
    wos[j] = rop_take(wo[j], ort)
    pos[j] = rop_take(po[j], ort)
  }
  for i = 0; i < n; i = i + 1 {
    tos[i] = rop_take(to[i], ort)
  }
  for kk = 0; kk < q; kk = kk + 1 {
    cos[kk] = rop_take(co[kk], ort)
  }
  RoplsOplsFit::{
    w: wfinal,
    t: tfinal,
    pm: pfinal,
    c: cfinal,
    u: ufinal,
    wo: wos,
    to: tos,
    po: pos,
    co: cos,
    r2x: rop_take(r2x_v, ort + 1),
    r2x_cum: rop_take(r2x_cumv, ort + 1),
    r2y: rop_take(r2y_v, ort + 1),
    r2y_cum: rop_take(r2y_cumv, ort + 1),
    q2: rop_take(q2_v, ort + 1),
    q2_cum: rop_take(q2_cumv, ort + 1),
    sig: rop_take(sigv, ort + 1),
    ort,
  }
}

///|
/// OPLS VIP4,p and VIP4,o (Galindo-Prieto et al.).
fn rop_opls_vip(fit : RoplsOplsFit) -> (Array[Double], Array[Double]) {
  let p = fit.wo.length()
  let sxp = rop_outer_ss(fit.t, fit.pm)
  let sxo = {
    let mut s = 0.0
    for h = 0; h < fit.ort; h = h + 1 {
      let mut tt = 0.0
      let mut pp = 0.0
      for i = 0; i < fit.t.length(); i = i + 1 {
        tt = tt + fit.to[i][h] * fit.to[i][h]
      }
      for j = 0; j < p; j = j + 1 {
        pp = pp + fit.po[j][h] * fit.po[j][h]
      }
      s = s + tt * pp
    }
    s
  }
  let ssx = sxp + sxo
  let syp = rop_outer_ss(fit.t, fit.c)
  let syo = {
    let mut s = 0.0
    for h = 0; h < fit.ort; h = h + 1 {
      let mut tt = 0.0
      let mut cc = 0.0
      for i = 0; i < fit.t.length(); i = i + 1 {
        tt = tt + fit.to[i][h] * fit.to[i][h]
      }
      for k = 0; k < fit.co.length(); k = k + 1 {
        cc = cc + fit.co[k][h] * fit.co[k][h]
      }
      s = s + tt * cc
    }
    s
  }
  let ssy = syp + syo
  // Predictive VIP: single predictive component, normalized p loadings
  let mut pnorm = 0.0
  for j = 0; j < p; j = j + 1 {
    pnorm = pnorm + fit.pm[j] * fit.pm[j]
  }
  let vip = Array::make(p, 0.0)
  let ovip = Array::make(p, 0.0)
  if ssx > 0.0 && ssy > 0.0 && pnorm > 0.0 {
    let kp = p.to_double() / (sxp / ssx + syp / ssy)
    for j = 0; j < p; j = j + 1 {
      let pn2 = fit.pm[j] * fit.pm[j] / pnorm
      vip[j] = (kp * (pn2 * sxp / ssx + pn2 * syp / ssy)).sqrt()
    }
    let ko = p.to_double() / (sxo / ssx + syo / ssy)
    for j = 0; j < p; j = j + 1 {
      let mut numx = 0.0
      let mut numy = 0.0
      for h = 0; h < fit.ort; h = h + 1 {
        let mut pn = 0.0
        for jj = 0; jj < p; jj = jj + 1 {
          pn = pn + fit.po[jj][h] * fit.po[jj][h]
        }
        let pnj = fit.po[j][h]
        let pn2 = if pn > 0.0 { pnj * pnj / pn } else { 0.0 }
        let mut tt = 0.0
        let mut pp = 0.0
        for i = 0; i < fit.t.length(); i = i + 1 {
          tt = tt + fit.to[i][h] * fit.to[i][h]
        }
        for jj = 0; jj < p; jj = jj + 1 {
          pp = pp + fit.po[jj][h] * fit.po[jj][h]
        }
        let mut cc = 0.0
        for k = 0; k < fit.co.length(); k = k + 1 {
          cc = cc + fit.co[k][h] * fit.co[k][h]
        }
        let sxh = tt * pp
        let syh = tt * cc
        numx = numx + pn2 * sxh
        numy = numy + pn2 * syh
      }
      ovip[j] = (ko * (numx / ssx + numy / ssy)).sqrt()
    }
  }
  (vip, ovip)
}

// ============================================================================
// Public fitting entry
// ============================================================================

///|
fn rop_classify(
  raw : Array[Array[Double]],
  classes : Array[String],
) -> Array[String] {
  let k = classes.length()
  let out : Array[String] = []
  for row in raw {
    if k == 2 {
      out.push(if row[0] > 0.5 { classes[1] } else { classes[0] })
    } else {
      let mut best = 0
      for j = 1; j < k; j = j + 1 {
        if row[j] > row[best] {
          best = j
        }
      }
      out.push(classes[best])
    }
  }
  out
}

///|
/// Fit a PLS(-DA) (`ortho = 0`) or OPLS(-DA) (`ortho > 0`) model.
///
/// - `y` is either `RoplsY::Continuous(Y)` with Y a samples x responses
///   numeric matrix, or `RoplsY::Classes(labels)`.
/// - `ncomp` predictive components (ignored for OPLS which uses 1); `0`
///   enables automatic component selection for PLS.
/// - `ortho` number of Y-orthogonal components; `-1` selects automatically.
/// - `scale` one of "standard" (default), "pareto", "center", "none".
/// - `crossval` number of cross-validation folds (default 7).
pub fn ropls_fit(
  x : Array[Array[Double]],
  y~ : RoplsY,
  ncomp? : Int = 0,
  ortho? : Int = 0,
  scale? : String = "standard",
  crossval? : Int = 7,
) -> RoplsModel raise {
  let pre0 = rop_prepare(x, y, scale)
  let da = pre0.da
  if ortho != 0 {
    if pre0.q > 1 {
      raise SeqError("OPLS regression only available for a single 'y' response")
    }
    if da && pre0.classes.length() > 2 {
      raise SeqError(
        "OPLS-DA only available for binary classification (use PLS-DA for multiple classes)",
      )
    }
  }
  if ortho == 0 {
    let autol = ncomp == 0
    let fit = rop_fit_pls(pre0, if autol { 0 } else { ncomp }, crossval, autol)
    let rstar = rop_weight_star(fit.w, fit.p)
    let b = rop_coefficients(rstar, fit.c)
    let yscaled = rop_tc(fit.t, fit.c)
    let ypred = rop_descale_y(yscaled, pre0.ymean, pre0.ysd)
    let r2x_cum : Array[Double] = []
    let r2y_cum : Array[Double] = []
    let q2_cum : Array[Double] = []
    let mut accx = 0.0
    let mut accy = 0.0
    let mut prod = 1.0
    for h = 0; h < fit.ncomp; h = h + 1 {
      accx = accx + fit.r2x[h]
      accy = accy + fit.r2y[h]
      prod = prod * (1.0 - fit.q2[h])
      r2x_cum.push(accx)
      r2y_cum.push(accy)
      q2_cum.push(1.0 - prod)
    }
    let vip = rop_pls_vip(fit.w, fit.t, fit.c)
    let rmse = rop_rmse(ypred, rop_original_y(y))
    let rmsee = (rmse *
    rmse *
    pre0.n.to_double() /
    (pre0.n.to_double() - (1.0 + fit.ncomp.to_double()))).sqrt()
    let fitted = if da { rop_classify(ypred, pre0.classes) } else { [] }
    let kind = if da { "PLS-DA" } else { "PLS" }
    RoplsModel::{
      kind,
      n: pre0.n,
      p: pre0.p,
      q: pre0.q,
      pre: fit.ncomp,
      ort: 0,
      scale,
      crossval,
      scores: fit.t,
      x_loadings: fit.p,
      weights: fit.w,
      weight_star: rstar,
      y_scores: fit.u,
      y_loadings: fit.c,
      coefficients: b,
      vip,
      ortho_scores: [],
      ortho_loadings: [],
      ortho_weights: [],
      ortho_y_loadings: [],
      ortho_vip: [],
      r2x: fit.r2x,
      r2x_cum,
      r2y: fit.r2y,
      r2y_cum,
      q2: fit.q2,
      q2_cum,
      significance: fit.sig,
      total_r2x: accx,
      total_r2y: accy,
      total_q2: q2_cum[q2_cum.length() - 1],
      rmsee,
      y_pred: ypred,
      fitted_classes: fitted,
      classes: pre0.classes,
      x_mean: pre0.xmean,
      x_sd: pre0.xsd,
      y_mean: pre0.ymean,
      y_sd: pre0.ysd,
      dropped: pre0.dropped,
    }
  } else {
    let autol = ortho < 0
    let fit = rop_fit_opls(pre0, if autol { 0 } else { ortho }, crossval, autol)
    // Single predictive component matrices
    let p = pre0.p
    let q = pre0.q
    let n = pre0.n
    let wm : Array[Array[Double]] = Array::make(p, [])
    let pm : Array[Array[Double]] = Array::make(p, [])
    let tm : Array[Array[Double]] = Array::make(n, [])
    let um : Array[Array[Double]] = Array::make(n, [])
    let cm : Array[Array[Double]] = Array::make(q, [])
    for j = 0; j < p; j = j + 1 {
      wm[j] = [fit.w[j]]
      pm[j] = [fit.pm[j]]
    }
    for i = 0; i < n; i = i + 1 {
      tm[i] = [fit.t[i]]
      um[i] = [fit.u[i]]
    }
    for k = 0; k < q; k = k + 1 {
      cm[k] = [fit.c[k]]
    }
    let b = rop_coefficients(wm, cm)
    let yscaled = rop_tc(tm, cm)
    let ypred = rop_descale_y(yscaled, pre0.ymean, pre0.ysd)
    let (vip, ovip) = rop_opls_vip(fit)
    let rmse = rop_rmse(ypred, rop_original_y(y))
    let rmsee = (rmse *
    rmse *
    n.to_double() /
    (n.to_double() - (1.0 + 1.0 + fit.ort.to_double()))).sqrt()
    let mut total_r2y = 0.0
    let mut total_q2 = 0.0
    for h = 0; h < fit.r2y.length(); h = h + 1 {
      total_r2y = total_r2y + fit.r2y[h]
      total_q2 = total_q2 + fit.q2[h]
    }
    // R's modelDF appends a "sum" row: cumulative totals only, with
    // per-component columns set to NA and Signif. = "NA".
    let total_r2x = fit.r2x_cum[fit.r2x_cum.length() - 1]
    let r2x_tbl = fit.r2x.copy()
    let r2x_cum_tbl = fit.r2x_cum.copy()
    let r2y_tbl = fit.r2y.copy()
    let r2y_cum_tbl = fit.r2y_cum.copy()
    let q2_tbl = fit.q2.copy()
    let q2_cum_tbl = fit.q2_cum.copy()
    let sig_tbl = fit.sig.copy()
    let na = 0.0 / 0.0
    r2x_tbl.push(na)
    r2x_cum_tbl.push(total_r2x)
    r2y_tbl.push(na)
    r2y_cum_tbl.push(total_r2y)
    q2_tbl.push(na)
    q2_cum_tbl.push(total_q2)
    sig_tbl.push("NA")
    let fitted = if da { rop_classify(ypred, pre0.classes) } else { [] }
    let kind = if da { "OPLS-DA" } else { "OPLS" }
    RoplsModel::{
      kind,
      n,
      p,
      q,
      pre: 1,
      ort: fit.ort,
      scale,
      crossval,
      scores: tm,
      x_loadings: pm,
      weights: wm,
      weight_star: wm,
      y_scores: um,
      y_loadings: cm,
      coefficients: b,
      vip,
      ortho_scores: fit.to,
      ortho_loadings: fit.po,
      ortho_weights: fit.wo,
      ortho_y_loadings: fit.co,
      ortho_vip: ovip,
      r2x: r2x_tbl,
      r2x_cum: r2x_cum_tbl,
      r2y: r2y_tbl,
      r2y_cum: r2y_cum_tbl,
      q2: q2_tbl,
      q2_cum: q2_cum_tbl,
      significance: sig_tbl,
      total_r2x,
      total_r2y,
      total_q2,
      rmsee,
      y_pred: ypred,
      fitted_classes: fitted,
      classes: pre0.classes,
      x_mean: pre0.xmean,
      x_sd: pre0.xsd,
      y_mean: pre0.ymean,
      y_sd: pre0.ysd,
      dropped: pre0.dropped,
    }
  }
}

///|
/// Recover the numeric Y matrix in its original scale (dummy matrix for DA).
fn rop_original_y(resp : RoplsY) -> Array[Array[Double]] raise {
  match resp {
    RoplsY::Continuous(y) => y
    RoplsY::Classes(labels) => {
      let (_, d) = rop_dummy(labels)
      d
    }
  }
}

// ============================================================================
// Prediction
// ============================================================================

///|
fn rop_scale_new(
  xnew : Array[Array[Double]],
  model : RoplsModel,
) -> Array[Array[Double]] raise {
  let pfull = model.p + model.dropped.length()
  if xnew.length() == 0 || xnew[0].length() != pfull {
    raise SeqError(
      "ropls: newdata must have the same number of variables as the training data",
    )
  }
  let mut mat = xnew
  if model.dropped.length() > 0 {
    let keep : Array[Int] = []
    for j = 0; j < pfull; j = j + 1 {
      if not(model.dropped.contains(j)) {
        keep.push(j)
      }
    }
    mat = rop_select_columns(xnew, keep)
  }
  rop_apply_scale(mat, model.x_mean, model.x_sd)
}

///|
/// Predict responses for new X data, returned in the original Y scale
/// (dummy 0/1 scale for DA models).
pub fn ropls_predict(
  model : RoplsModel,
  xnew : Array[Array[Double]],
) -> Array[Array[Double]] raise {
  let xs = rop_scale_new(xnew, model)
  for h = 0; h < model.ort; h = h + 1 {
    let to = rop_mat_vec(xs, model.ortho_weights.map(fn(row) { row[h] }))
    let po = model.ortho_loadings.map(fn(row) { row[h] })
    rop_deflate(xs, to, po)
  }
  let ys = rop_coefficients_predict(xs, model.coefficients)
  rop_descale_y(ys, model.y_mean, model.y_sd)
}

///|
fn rop_coefficients_predict(
  xs : Array[Array[Double]],
  b : Array[Array[Double]],
) -> Array[Array[Double]] {
  let n = xs.length()
  let q = if b.length() == 0 { 0 } else { b[0].length() }
  let out : Array[Array[Double]] = Array::make(n, [])
  for i = 0; i < n; i = i + 1 {
    out[i] = Array::make(q, 0.0)
    for j = 0; j < b.length(); j = j + 1 {
      for k = 0; k < q; k = k + 1 {
        out[i][k] = out[i][k] + xs[i][j] * b[j][k]
      }
    }
  }
  out
}

///|
/// Predict class labels for new X data with a DA model.
pub fn ropls_predict_classes(
  model : RoplsModel,
  xnew : Array[Array[Double]],
) -> Array[String] raise {
  if model.classes.length() == 0 {
    raise SeqError("ropls: predict_classes requires a DA model")
  }
  let raw = ropls_predict(model, xnew)
  rop_classify(raw, model.classes)
}

// ============================================================================
// Permutation testing
// ============================================================================

///|
/// Response permutation test (Szymanska et al., 2012).
///
/// `perms` contains permutations as arrays of original row indices; the
/// response of sample i is replaced by `y[perms[i]]`. The same fixed numbers
/// of predictive / orthogonal components are used for every refit. Returns
/// permuted R2Y/Q2 values plus the (1 + count) / nperm p-values.
pub fn ropls_permutation(
  x : Array[Array[Double]],
  y~ : RoplsY,
  perms~ : Array[Array[Int]],
  ncomp? : Int = 1,
  ortho? : Int = 0,
  scale? : String = "standard",
  crossval? : Int = 7,
) -> RoplsPermResult raise {
  if perms.length() == 0 {
    raise SeqError("ropls: permutation test requires at least one permutation")
  }
  let observed = ropls_fit(x, y~, ncomp~, ortho~, scale~, crossval~)
  let n = observed.n
  let perm_r2y : Array[Double] = []
  let perm_q2 : Array[Double] = []
  let sims : Array[Double] = []
  let original_numeric = rop_original_y(y)
  for perm in perms {
    if perm.length() != n {
      raise SeqError("ropls: each permutation must cover all samples")
    }
    for idx in perm {
      if idx < 0 || idx >= n {
        raise SeqError("ropls: permutation index out of range")
      }
    }
    let perm_resp = match y {
      RoplsY::Continuous(_) => {
        let yp : Array[Array[Double]] = perm.map(fn(i) {
          original_numeric[i].copy()
        })
        RoplsY::Continuous(yp)
      }
      RoplsY::Classes(labels) => {
        let lp : Array[String] = perm.map(fn(i) { labels[i] })
        RoplsY::Classes(lp)
      }
    }
    let m = ropls_fit(x, y=perm_resp, ncomp~, ortho~, scale~, crossval~)
    perm_r2y.push(m.total_r2y)
    perm_q2.push(m.total_q2)
    // Similarity: label agreement for DA, correlation for a continuous
    // single response (matching ropls .similarityF)
    let sim = match y {
      RoplsY::Classes(labels) => {
        let mut agree = 0
        for i = 0; i < n; i = i + 1 {
          if labels[perm[i]] == labels[i] {
            agree = agree + 1
          }
        }
        agree.to_double() / n.to_double()
      }
      RoplsY::Continuous(_) =>
        if observed.q == 1 {
          let a = Array::make(n, 0.0)
          let b = Array::make(n, 0.0)
          for i = 0; i < n; i = i + 1 {
            a[i] = original_numeric[i][0]
            b[i] = original_numeric[perm[i]][0]
          }
          rop_cor(a, b)
        } else {
          let mut s = 0.0
          for k = 0; k < observed.q; k = k + 1 {
            let a = Array::make(n, 0.0)
            let b = Array::make(n, 0.0)
            for i = 0; i < n; i = i + 1 {
              a[i] = original_numeric[i][k]
              b[i] = original_numeric[perm[i]][k]
            }
            s = s + rop_cor(a, b)
          }
          s / observed.q.to_double()
        }
    }
    sims.push(sim)
  }
  let mut ge_r2 = 0
  let mut ge_q2 = 0
  for v in perm_r2y {
    if v >= observed.total_r2y {
      ge_r2 = ge_r2 + 1
    }
  }
  for v in perm_q2 {
    if v >= observed.total_q2 {
      ge_q2 = ge_q2 + 1
    }
  }
  let denom = perms.length().to_double()
  RoplsPermResult::{
    nperm: perms.length(),
    obs_r2y: observed.total_r2y,
    obs_q2: observed.total_q2,
    perm_r2y,
    perm_q2,
    sim: sims,
    p_r2y: (1.0 + ge_r2.to_double()) / denom,
    p_q2: (1.0 + ge_q2.to_double()) / denom,
  }
}
