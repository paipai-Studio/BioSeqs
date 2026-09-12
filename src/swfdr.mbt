///|
/// Bio.Swfdr module - science-wise false discovery rate estimation.
///
/// Port of the Bioconductor `swfdr` package (Bourgon 2018):
///   - `lm_pi0`: covariate-conditioned pi0 estimation via per-lambda
///     logistic/linear regression, followed over lambda by either a
///     `splines::bs` penalized unit spline or an exact replica of R's
///     `smooth.spline` (Cox-de Boor B-splines, sgram.f penalty matrix,
///     sbart.c Brent search in spar space).
///   - `lm_qvalue`: covariate-conditioned q-values (BH / pFDR).
///   - `calculateSwfdr`: EM estimation of pi0 and the Beta(alpha, beta)
///     distribution of alternative p-values, accounting for rounded and
///     truncated p-values.
///
/// Zero external dependencies; all spline, linear-algebra and incomplete
/// beta machinery is implemented locally and validated against R 4.3.3.

// ===========================================================================
// Public types and errors
// ===========================================================================

///|
/// Error raised by the swfdr functions, mirroring the stop() messages of the
/// Bioconductor package (input validation of p, lambda, df and X).
pub suberror SwfdrError {
  SwfdrError(String)
}

///|
/// Result of `swfdr_lm_pi0`, mirroring the R `lm_pi0` object.
///
/// Fields:
///   - lambda:     sorted unique lambda thresholds
///   - x_names:    covariate names (empty, names are not tracked)
///   - pi0_lambda: n_tests x n_lambda matrix of per-lambda pi0 estimates
///   - pi0:        smoothed pi0(x) estimate, one entry per test
///   - pi0_smooth: n_tests x n_lambda fitted smoother values
///                 (empty for smoothing="unit.spline")
///   - smoothing:  the smoothing method used
///   - type:       the regression type used
pub struct SwfdrPi0 {
  lambda : Array[Double]
  x_names : Array[String]
  pi0_lambda : Array[Array[Double]]
  pi0 : Array[Double]
  pi0_smooth : Array[Array[Double]]
  smoothing : String
  type_ : String
}

///|
pub fn SwfdrPi0::lambda(self : SwfdrPi0) -> Array[Double] {
  self.lambda
}

///|
pub fn SwfdrPi0::pi0_lambda(self : SwfdrPi0) -> Array[Array[Double]] {
  self.pi0_lambda
}

///|
pub fn SwfdrPi0::pi0(self : SwfdrPi0) -> Array[Double] {
  self.pi0
}

///|
pub fn SwfdrPi0::pi0_smooth(self : SwfdrPi0) -> Array[Array[Double]] {
  self.pi0_smooth
}

///|
pub fn SwfdrPi0::smoothing(self : SwfdrPi0) -> String {
  self.smoothing
}

///|
pub fn SwfdrPi0::regression_type(self : SwfdrPi0) -> String {
  self.type_
}

///|
/// Result of `swfdr_lm_qvalue`, mirroring the R `lm_qvalue` object.
pub struct SwfdrQvalue {
  pvalues : Array[Double]
  qvalues : Array[Double]
  pi0_result : SwfdrPi0
  pfdr : Bool
}

///|
pub fn SwfdrQvalue::pvalues(self : SwfdrQvalue) -> Array[Double] {
  self.pvalues
}

///|
pub fn SwfdrQvalue::qvalues(self : SwfdrQvalue) -> Array[Double] {
  self.qvalues
}

///|
pub fn SwfdrQvalue::pi0(self : SwfdrQvalue) -> SwfdrPi0 {
  self.pi0_result
}

///|
pub fn SwfdrQvalue::pfdr(self : SwfdrQvalue) -> Bool {
  self.pfdr
}

///|
/// Result of `swfdr_calculate`, mirroring the R `calculateSwfdr` object.
///
/// Fields:
///   - pi0:   estimated prior null probability (science-wise FDR)
///   - alpha: Beta(alpha, beta) alternative-distribution shape 1
///   - beta:  Beta(alpha, beta) alternative-distribution shape 2
///   - z:     posterior null indicator for non-rounded p-values (None where
///            the p-value is rounded)
///   - n0:    expected number of rounded null p-values in each of the six
///            bins (0, .005], (.005, .015], ..., (.045, .05]
///   - n:     observed number of rounded p-values in each of the six bins
pub struct SwfdrScience {
  pi0 : Double
  alpha : Double
  beta : Double
  z : Array[Double?]
  n0 : Array[Double]
  n : Array[Int]
}

///|
pub fn SwfdrScience::pi0(self : SwfdrScience) -> Double {
  self.pi0
}

///|
pub fn SwfdrScience::alpha(self : SwfdrScience) -> Double {
  self.alpha
}

///|
pub fn SwfdrScience::beta(self : SwfdrScience) -> Double {
  self.beta
}

///|
pub fn SwfdrScience::z(self : SwfdrScience) -> Array[Double?] {
  self.z
}

///|
pub fn SwfdrScience::n0(self : SwfdrScience) -> Array[Double] {
  self.n0
}

///|
pub fn SwfdrScience::n(self : SwfdrScience) -> Array[Int] {
  self.n
}

// ===========================================================================
// Dense linear algebra
// ===========================================================================

///|
/// Allocate an r x c zero matrix with independent rows.
fn _sw_mat(r : Int, c : Int) -> Array[Array[Double]] {
  let m : Array[Array[Double]] = []
  for _ in 0..<r {
    m.push(Array::make(c, 0.0))
  }
  m
}

///|
/// Cross product A^T A.
fn _sw_ata(a : Array[Array[Double]]) -> Array[Array[Double]] {
  let n = a.length()
  let p = if n == 0 { 0 } else { a[0].length() }
  let m = _sw_mat(p, p)
  let mut i = 0
  while i < n {
    let row = a[i]
    let mut j = 0
    while j < p {
      let aij = row[j]
      let mut k = 0
      while k < p {
        m[j][k] = m[j][k] + aij * row[k]
        k = k + 1
      }
      j = j + 1
    }
    i = i + 1
  }
  m
}

///|
/// Product A^T B where A is n x p and B n x q.
fn _sw_atb(
  a : Array[Array[Double]],
  b : Array[Array[Double]],
) -> Array[Array[Double]] {
  let n = a.length()
  let p = if n == 0 { 0 } else { a[0].length() }
  let q = if n == 0 { 0 } else { b[0].length() }
  let m = _sw_mat(p, q)
  let mut i = 0
  while i < n {
    let ra = a[i]
    let rb = b[i]
    let mut j = 0
    while j < p {
      let aij = ra[j]
      let mut k = 0
      while k < q {
        m[j][k] = m[j][k] + aij * rb[k]
        k = k + 1
      }
      j = j + 1
    }
    i = i + 1
  }
  m
}

///|
/// Matrix times matrix.
fn _sw_matmul(
  a : Array[Array[Double]],
  b : Array[Array[Double]],
) -> Array[Array[Double]] {
  let n = a.length()
  let p = if n == 0 { 0 } else { a[0].length() }
  let q = if p == 0 { 0 } else { b[0].length() }
  let m = _sw_mat(n, q)
  let mut i = 0
  while i < n {
    let ra = a[i]
    let mut k = 0
    while k < p {
      let aik = ra[k]
      if aik != 0.0 {
        let rb = b[k]
        let mut j = 0
        while j < q {
          m[i][j] = m[i][j] + aik * rb[j]
          j = j + 1
        }
      }
      k = k + 1
    }
    i = i + 1
  }
  m
}

///|
/// Dense Cholesky decomposition M = L L^T of a symmetric positive definite
/// matrix. Returns None when the matrix is not numerically positive
/// definite (mirroring R chol() failure in the spline search).
fn _sw_cholesky(m0 : Array[Array[Double]]) -> Array[Array[Double]]? {
  let n = m0.length()
  let l = _sw_mat(n, n)
  let mut i = 0
  while i < n {
    let mut j = 0
    while j <= i {
      let mut s = m0[i][j]
      let mut k = 0
      while k < j {
        s = s - l[i][k] * l[j][k]
        k = k + 1
      }
      if i == j {
        if s <= 1.0e-300 {
          return None
        }
        l[i][i] = s.sqrt()
      } else {
        l[i][j] = s / l[j][j]
      }
      j = j + 1
    }
    i = i + 1
  }
  Some(l)
}

///|
/// Solve M x = rhs for a symmetric positive definite M (Cholesky path).
fn _sw_chol_solve(
  m : Array[Array[Double]],
  rhs : Array[Double],
) -> Array[Double]? {
  match _sw_cholesky(m) {
    None => None
    Some(l) => {
      let n = l.length()
      let y = Array::make(n, 0.0)
      let mut i = 0
      while i < n {
        let mut s = rhs[i]
        let mut k = 0
        while k < i {
          s = s - l[i][k] * y[k]
          k = k + 1
        }
        y[i] = s / l[i][i]
        i = i + 1
      }
      let x = Array::make(n, 0.0)
      i = n - 1
      while i >= 0 {
        let mut s = y[i]
        let mut k = i + 1
        while k < n {
          s = s - l[k][i] * x[k]
          k = k + 1
        }
        x[i] = s / l[i][i]
        i = i - 1
      }
      Some(x)
    }
  }
}

///|
/// Inverse of a symmetric positive definite matrix (sizes here are at most
/// ~22 x 22), or None on failure.
fn _sw_chol_inverse(m : Array[Array[Double]]) -> Array[Array[Double]]? {
  let n = m.length()
  let inv = _sw_mat(n, n)
  let mut j = 0
  while j < n {
    let e = Array::make(n, 0.0)
    e[j] = 1.0
    match _sw_chol_solve(m, e) {
      None => return None
      Some(col) => {
        let mut i = 0
        while i < n {
          inv[i][j] = col[i]
          i = i + 1
        }
      }
    }
    j = j + 1
  }
  Some(inv)
}

// ===========================================================================
// Cox-de Boor B-spline basis (dense replica of splineDesign / sgram path)
// ===========================================================================

///|
/// R findInterval(x, knots, rightmost.closed=TRUE, all.inside=TRUE):
/// index i (0-based) with knots[i] <= x < knots[i+1], clamped to
/// [0, length(knots)-2].
fn _sw_find_interval(knots : Array[Double], x : Double) -> Int {
  let m = knots.length()
  let mut lo = 0
  let mut hi = m - 1
  while hi - lo > 1 {
    let mid = (lo + hi) / 2
    if knots[mid] <= x {
      lo = mid
    } else {
      hi = mid
    }
  }
  if lo > m - 2 {
    m - 2
  } else {
    lo
  }
}

///|
/// Build the Cox-de Boor recursion table up to the given order for one point.
/// Element k (0-based) has length length(knots) - (k + 1).
fn _sw_cxb_table(
  knots : Array[Double],
  x : Double,
  ord : Int,
) -> Array[Array[Double]] {
  let m = knots.length()
  let levels : Array[Array[Double]] = []
  let n1 = Array::make(m - 1, 0.0)
  let mut i = _sw_find_interval(knots, x)
  // Right multiple-boundary point: use the last nonzero span (R's
  // 1-based i = m - ord translates to 0-based m - ord - 1).
  if x == knots[m - 1] {
    i = m - ord - 1
  }
  n1[i] = 1.0
  levels.push(n1)
  let mut k = 2
  while k <= ord {
    let nk = m - k
    let v = Array::make(nk, 0.0)
    let prev = levels[k - 2]
    let mut j = 0
    while j < nk {
      let d1 = knots[j + k - 1] - knots[j]
      let d2 = knots[j + k] - knots[j + 1]
      let a = if d1 != 0.0 { (x - knots[j]) / d1 } else { 0.0 }
      let b = if d2 != 0.0 { (knots[j + k] - x) / d2 } else { 0.0 }
      v[j] = a * prev[j] + b * prev[j + 1]
      j = j + 1
    }
    levels.push(v)
    k = k + 1
  }
  levels
}

///|
/// Order-`ord` B-spline basis values (length length(knots) - ord) at x.
fn _sw_cxb_basis(knots : Array[Double], x : Double, ord : Int) -> Array[Double] {
  _sw_cxb_table(knots, x, ord)[ord - 1].copy()
}

///|
/// Second derivative of the order-4 basis at x (dense replica of bspvd via
/// the de Boor triangle). Each differentiation rebuilds a fresh table from
/// the previous full table, as required by the recursion.
fn _sw_cxb_deriv2(knots : Array[Double], x : Double) -> Array[Double] {
  let ord = 4
  let m = knots.length()
  let mut d = _sw_cxb_table(knots, x, ord)
  let mut s = 1
  while s <= 2 {
    let dn : Array[Array[Double]] = [[], [], [], []]
    let mut k = s + 1
    while k <= ord {
      let nk = m - k
      let v = Array::make(nk, 0.0)
      let prev = d[k - 2]
      let mut i = 0
      while i < nk {
        let d1 = knots[i + k - 1] - knots[i]
        let d2 = knots[i + k] - knots[i + 1]
        let a = if d1 != 0.0 { prev[i] / d1 } else { 0.0 }
        let b = if d2 != 0.0 { prev[i + 1] / d2 } else { 0.0 }
        v[i] = (k - 1).to_double() * (a - b)
        i = i + 1
      }
      dn[k - 1] = v
      k = k + 1
    }
    d = dn
    s = s + 1
  }
  d[3].copy()
}

///|
/// Dense replica of R's sgram.f: the penalty matrix
/// S = integral B''(t)^T B''(t) dt, integrating the piecewise-linear second
/// derivatives exactly over each nonzero knot span (the literal 0.3330
/// coefficient from the Fortran source is preserved).
fn _sw_sgram(knots : Array[Double]) -> Array[Array[Double]] {
  let nb = knots.length() - 4
  let smat = _sw_mat(nb, nb)
  let mut i = 0
  while i < nb {
    let t1 = knots[i]
    let t2 = knots[i + 1]
    let wpt = t2 - t1
    if wpt != 0.0 {
      let mid = (t1 + t2) / 2.0
      let brow = _sw_cxb_basis(knots, mid, 4)
      let idx : Array[Int] = []
      let mut j = 0
      while j < nb {
        if brow[j].abs() > 1.0e-13 {
          idx.push(j)
        }
        j = j + 1
      }
      // at most the 4 cubic basis functions are nonzero on a span
      let ii : Array[Int] = if idx.length() > 4 {
        let tail4 : Array[Int] = []
        let mut t = idx.length() - 4
        while t < idx.length() {
          tail4.push(idx[t])
          t = t + 1
        }
        tail4
      } else {
        idx
      }
      let dl = _sw_cxb_deriv2(knots, t1)
      let du = _sw_cxb_deriv2(knots, t2)
      let ni = ii.length()
      let mut a = 0
      while a < ni {
        let la = dl[ii[a]]
        let dda = du[ii[a]] - la
        let mut b = 0
        while b < ni {
          let lb = dl[ii[b]]
          let ddb = du[ii[b]] - lb
          smat[ii[a]][ii[b]] = smat[ii[a]][ii[b]] +
            wpt * (la * lb + 0.5 * (dda * lb + ddb * la) + 0.3330 * dda * ddb)
          b = b + 1
        }
        a = a + 1
      }
    }
    i = i + 1
  }
  smat
}

// ===========================================================================
// stats::smooth.spline replica (all.knots, n < 50, unit weights)
// ===========================================================================

///|
/// Fit a cubic smoothing spline to x, y with target degrees of freedom,
/// returning the fitted values at x. Exact replica of R
/// `smooth.spline(x, y, df=df, tol=1e-7, keep.data=FALSE)` for n < 50:
/// triple boundary knots, sbart.c Brent minimization over spar in
/// [-1.5, 1.5], and the documented Fortran quirk of returning the fit at
/// the *last evaluated trial point* rather than at the optimum.
fn _sw_smooth_spline(
  x : Array[Double],
  y : Array[Double],
  df : Double,
) -> Array[Double] {
  let n = x.length()
  if n < 4 || n != y.length() {
    return y.copy()
  }
  let x0v = x[0]
  let span = x[n - 1] - x[0]
  let ux = Array::make(n, 0.0)
  let mut i = 0
  while i < n {
    ux[i] = (x[i] - x0v) / span
    i = i + 1
  }
  // triple (not quadruple) boundary knots: c(rep(0,3), ux, rep(1,3))
  let knot : Array[Double] = []
  knot.push(0.0)
  knot.push(0.0)
  knot.push(0.0)
  for v in ux {
    knot.push(v)
  }
  knot.push(1.0)
  knot.push(1.0)
  knot.push(1.0)
  let nk = knot.length() - 4
  let bmat = _sw_mat(n, nk)
  i = 0
  while i < n {
    let row = _sw_cxb_basis(knot, ux[i], 4)
    let mut j = 0
    while j < nk {
      bmat[i][j] = row[j]
      j = j + 1
    }
    i = i + 1
  }
  let xx = _sw_ata(bmat)
  let sg = _sw_sgram(knot)
  // R: sum(diag(XX)[3:(nk-3)]) / sum(diag(SG)[3:(nk-3)]), i.e. 0-based 2..nk-4
  let mut num = 0.0
  let mut den = 0.0
  i = 2
  while i <= nk - 4 {
    num = num + xx[i][i]
    den = den + sg[i][i]
    i = i + 1
  }
  let ratio = num / den
  let lamof = fn(spar : Double, ratio : Double) -> Double {
    ratio * @math.exp((6.0 * spar - 2.0) * @math.ln(16.0))
  }
  // trace of the hat matrix B (XX + lambda SG)^{-1} B' = tr(M^{-1} XX)
  let traceof = fn(
    lambda : Double,
    xx : Array[Array[Double]],
    sg : Array[Array[Double]],
  ) -> Double? {
    let mm = _sw_mat(xx.length(), xx.length())
    let mut a = 0
    while a < xx.length() {
      let mut b = 0
      while b < xx.length() {
        mm[a][b] = xx[a][b] + lambda * sg[a][b]
        b = b + 1
      }
      a = a + 1
    }
    let inv = match _sw_chol_inverse(mm) {
      None => {
        // Near lambda = 0 two eigen-directions collapse; R's base solve still
        // returns a finite trace there. Stabilize with a proportional jitter
        // that is invisible once lambda reaches the df-matching scale.
        let mut dmax = 0.0
        let mut a2 = 0
        while a2 < xx.length() {
          if mm[a2][a2].abs() > dmax {
            dmax = mm[a2][a2].abs()
          }
          a2 = a2 + 1
        }
        let jit = 1.0e-14 * dmax
        let mut a3 = 0
        while a3 < xx.length() {
          mm[a3][a3] = mm[a3][a3] + jit
          a3 = a3 + 1
        }
        match _sw_chol_inverse(mm) {
          None => return None
          Some(iv) => iv
        }
      }
      Some(inv) => inv
    }
    let mut t = 0.0
    let mut j = 0
    while j < xx.length() {
      let mut k = 0
      while k < xx.length() {
        t = t + inv[j][k] * xx[k][j]
        k = k + 1
      }
      j = j + 1
    }
    Some(t)
  }
  // Brent search (sbart.c) over spar
  let gold = 0.381966011250105151795413165634
  let eps = 2.0e-8
  let tol0 = 1.0e-4
  let maxit = 500
  let mut a = -1.5
  let mut b = 1.5
  let mut xx0 = a + gold * (b - a)
  let mut v = xx0
  let mut w = xx0
  let mut e = 0.0
  let mut d = 0.0
  let mut last_lambda = lamof(xx0, ratio)
  let comp = fn(
    s : Double,
    ratio : Double,
    df : Double,
    xx : Array[Array[Double]],
    sg : Array[Array[Double]],
  ) -> Double? {
    let lambda = ratio * @math.exp((6.0 * s - 2.0) * @math.ln(16.0))
    match traceof(lambda, xx, sg) {
      None => None
      Some(tr) => Some(3.0 + (tr - df) * (tr - df))
    }
  }
  let mut ier = false
  let mut fx = match comp(xx0, ratio, df, xx, sg) {
    Some(val) => val
    None => {
      ier = true
      0.0
    }
  }
  let mut fv = fx
  let mut fw = fx
  let mut iter = 0
  while iter < maxit && !ier {
    let xm = (a + b) * 0.5
    let tol1 = eps * xx0.abs() + tol0 / 3.0
    let tol2 = 2.0 * tol1
    if (xx0 - xm).abs() <= tol2 - (b - a) * 0.5 {
      break
    }
    let mut use_golden = true
    if e.abs() > tol1 {
      let r = (xx0 - w) * (fx - fv)
      let q1 = (xx0 - v) * (fx - fw)
      let mut p = (xx0 - v) * q1 - (xx0 - w) * r
      let mut q = (q1 - r) * 2.0
      if q > 0.0 {
        p = -p
      }
      q = q.abs()
      let et = e
      e = d
      if p.abs() < (0.5 * q * et).abs() &&
        p > q * (a - xx0) &&
        p < q * (b - xx0) {
        d = p / q
        use_golden = false
        let uc = xx0 + d
        if uc - a < tol2 || b - uc < tol2 {
          d = if xm - xx0 >= 0.0 { tol1 } else { -tol1 }
        }
      }
    }
    if use_golden {
      e = if xx0 >= xm { a - xx0 } else { b - xx0 }
      d = gold * e
    }
    let step = if d.abs() >= tol1 { d } else if d >= 0.0 { tol1 } else { -tol1 }
    let u = xx0 + step
    let fu = match comp(u, ratio, df, xx, sg) {
      Some(val) => {
        last_lambda = lamof(u, ratio)
        val
      }
      None => {
        ier = true
        0.0
      }
    }
    if ier {
      break
    }
    if fu <= fx {
      if u >= xx0 {
        a = xx0
      } else {
        b = xx0
      }
      v = w
      fv = fw
      w = xx0
      fw = fx
      xx0 = u
      fx = fu
    } else {
      if u < xx0 {
        a = u
      } else {
        b = u
      }
      if fu <= fw || w == xx0 {
        v = w
        fv = fw
        w = u
        fw = fu
      } else if fu <= fv || v == xx0 || v == w {
        v = u
        fv = fu
      }
    }
    iter = iter + 1
  }
  if ier {
    // R falls back to the mean of y when the search hits a singular system
    // at small spar (df = 1 interpolation path).
    let mut mean = 0.0
    for v in y {
      mean = mean + v
    }
    return Array::make(n, mean / n.to_double())
  }
  // Fortran returns the fit at the last evaluated trial point.
  let mm = _sw_mat(nk, nk)
  i = 0
  while i < nk {
    let mut j = 0
    while j < nk {
      mm[i][j] = xx[i][j] + last_lambda * sg[i][j]
      j = j + 1
    }
    i = i + 1
  }
  let rhs = Array::make(nk, 0.0)
  i = 0
  while i < n {
    let mut j = 0
    while j < nk {
      rhs[j] = rhs[j] + bmat[i][j] * y[i]
      j = j + 1
    }
    i = i + 1
  }
  let mut coef = _sw_chol_solve(mm, rhs)
  if coef is None {
    // The two polynomial null-directions of B and SG leave M rank-deficient
    // at every lambda; R's Fortran solves this directly. A proportional
    // diagonal jitter fixes the nullspace without affecting the fit.
    let mut dmax = 0.0
    let mut jj = 0
    while jj < nk {
      if mm[jj][jj].abs() > dmax {
        dmax = mm[jj][jj].abs()
      }
      jj = jj + 1
    }
    jj = 0
    while jj < nk {
      mm[jj][jj] = mm[jj][jj] + 1.0e-14 * dmax
      jj = jj + 1
    }
    coef = _sw_chol_solve(mm, rhs)
  }
  let coef = match coef {
    Some(c) => c
    None => return y.copy()
  }
  let fit = Array::make(n, 0.0)
  i = 0
  while i < n {
    let mut s = 0.0
    let mut j = 0
    while j < nk {
      s = s + bmat[i][j] * coef[j]
      j = j + 1
    }
    fit[i] = s
    i = i + 1
  }
  fit
}

// ===========================================================================
// splines::bs penalized "unit spline" smoother (spline_pi0.R)
// ===========================================================================

///|
/// Construct the bs basis used by the unit spline smoother, mirroring
/// `bs(x, knots=x, Boundary.knots=c(0,1))`: an order-4 spline design on knots
/// c(rep(0,4), x, rep(1,4)) with the first column removed. Returns an
/// n x (n+3) matrix; `x` must be strictly increasing inside (0, 1).
pub fn swfdr_bs_unit(x : Array[Double]) -> Array[Array[Double]] {
  let n = x.length()
  let ak : Array[Double] = []
  ak.push(0.0)
  ak.push(0.0)
  ak.push(0.0)
  ak.push(0.0)
  for v in x {
    ak.push(v)
  }
  ak.push(1.0)
  ak.push(1.0)
  ak.push(1.0)
  ak.push(1.0)
  let nfull = ak.length() - 4 // n + 4
  let bmat = _sw_mat(n, nfull - 1)
  let mut i = 0
  while i < n {
    let row = _sw_cxb_basis(ak, x[i], 4)
    let mut j = 1
    while j < nfull {
      bmat[i][j - 1] = row[j]
      j = j + 1
    }
    i = i + 1
  }
  bmat
}

///|
/// Trace of the inverse of M + lambda*omega (used by optimize.multiplier).
fn _sw_pen_trace(
  xmat : Array[Array[Double]],
  omega : Array[Array[Double]],
  lambda : Double,
) -> Double? {
  let n = xmat.length()
  let mm = _sw_mat(n, n)
  let mut i = 0
  while i < n {
    let mut j = 0
    while j < n {
      mm[i][j] = xmat[i][j] + lambda * omega[i][j]
      j = j + 1
    }
    i = i + 1
  }
  match _sw_chol_inverse(mm) {
    None => None
    Some(inv) => {
      let mut t = 0.0
      let mut k = 0
      while k < n {
        t = t + inv[k][k]
        k = k + 1
      }
      Some(t)
    }
  }
}

///|
/// Recursive search (ported to a loop) from spline_pi0.R: find the Lagrange
/// multiplier making trace((X + lambda*omega)^{-1}) equal to `target`.
fn _sw_optimize_multiplier(
  xmat : Array[Array[Double]],
  omega : Array[Array[Double]],
  target : Double,
) -> Double {
  let tol = 1.0e-4
  let mut lo = 1.0e-6
  let mut mid = 1.0
  let mut hi = 2.0
  let mut vlo : Double? = None
  let mut vmi : Double? = None
  let mut vhi : Double? = None
  let mut iters = 0
  while iters < 10000 {
    iters = iters + 1
    let width = hi - lo
    if width / mid < tol {
      return mid
    }
    match vlo {
      None => vlo = _sw_pen_trace(xmat, omega, lo)
      Some(_) => ()
    }
    match vmi {
      None => vmi = _sw_pen_trace(xmat, omega, mid)
      Some(_) => ()
    }
    match vhi {
      None => vhi = _sw_pen_trace(xmat, omega, hi)
      Some(_) => ()
    }
    let slo = vlo.unwrap()
    let smi = vmi.unwrap()
    let shi = vhi.unwrap()
    if (slo - shi).abs() < tol {
      return (lo + mid + hi) / 3.0
    }
    if target < shi {
      let w = width
      let nlo = hi
      let nmi = hi + 2.0 * w
      let nhi = hi + 4.0 * w
      lo = nlo
      mid = nmi
      hi = nhi
      vlo = Some(smi)
      vmi = None
      vhi = None
    } else if target < smi {
      let nlo = mid
      let nmi = (mid + hi) / 2.0
      let nhi = hi
      lo = nlo
      mid = nmi
      hi = nhi
      vlo = Some(smi)
      vmi = None
      vhi = Some(shi)
    } else {
      let nlo = lo
      let nmi = (lo + mid) / 2.0
      let nhi = mid
      lo = nlo
      mid = nmi
      hi = nhi
      vlo = Some(slo)
      vmi = None
      vhi = Some(smi)
    }
  }
  mid
}

///|
/// Build the n x n unit-spline transformation matrix
/// B (B'B + multiplier * Omega)^{-1} B' with the multiplier chosen so that
/// the trace of the smoother equals 2 + `smooth_df` degrees of freedom.
pub fn swfdr_unit_spline_matrix(
  x : Array[Double],
  smooth_df : Int,
) -> Array[Array[Double]] {
  let bmat = swfdr_bs_unit(x)
  let nbs = bmat[0].length()
  let xtx = _sw_ata(bmat)
  // omega.sqrt rows i=3..nbs (1-based): (1,-2,1) at columns (i-2,i-1,i)
  let dmat = _sw_mat(nbs, nbs)
  let mut i = 3
  while i <= nbs {
    dmat[i - 1][i - 3] = 1.0
    dmat[i - 1][i - 2] = -2.0
    dmat[i - 1][i - 1] = 1.0
    i = i + 1
  }
  let omega = _sw_ata(dmat)
  let multiplier = _sw_optimize_multiplier(
    xtx,
    omega,
    (2 + smooth_df).to_double(),
  )
  let mm = _sw_mat(nbs, nbs)
  let mut a = 0
  while a < nbs {
    let mut b = 0
    while b < nbs {
      mm[a][b] = xtx[a][b] + multiplier * omega[a][b]
      b = b + 1
    }
    a = a + 1
  }
  let inv = match _sw_chol_inverse(mm) {
    Some(v) => v
    None =>
      // defensive: should not happen since mm is positive definite
      return _sw_mat(x.length(), x.length())
  }
  // transform = B %*% inv %*% t(B)
  let bt : Array[Array[Double]] = []
  let mut j = 0
  while j < nbs {
    let col : Array[Double] = []
    let mut r = 0
    while r < bmat.length() {
      col.push(bmat[r][j])
      r = r + 1
    }
    bt.push(col)
    j = j + 1
  }
  let tmp = _sw_matmul(inv, bt)
  _sw_matmul(bmat, tmp)
}

// ===========================================================================
// Logistic IRLS (glm binomial) and linear regression (lsfit)
// ===========================================================================

///|
/// Weighted least squares for design A (no implicit intercept) and response
/// rhs: coefficients = (A'A)^{-1} A'rhs via Cholesky, with a small ridge
/// fallback in the rank-deficient case.
fn _sw_wls_coef(a : Array[Array[Double]], rhs : Array[Double]) -> Array[Double] {
  let p = a[0].length()
  let ata = _sw_ata(a)
  let atr = Array::make(p, 0.0)
  let mut i = 0
  while i < a.length() {
    let mut j = 0
    while j < p {
      atr[j] = atr[j] + a[i][j] * rhs[i]
      j = j + 1
    }
    i = i + 1
  }
  let mut ridge = 0.0
  let mut attempt = 0
  while attempt < 6 {
    let mm = _sw_mat(p, p)
    let mut r = 0
    while r < p {
      let mut c = 0
      while c < p {
        mm[r][c] = ata[r][c]
        c = c + 1
      }
      mm[r][r] = mm[r][r] + ridge
      r = r + 1
    }
    match _sw_chol_solve(mm, atr) {
      Some(coef) => return coef
      None => ridge = if ridge == 0.0 { 1.0e-10 } else { ridge * 100.0 }
    }
    attempt = attempt + 1
  }
  Array::make(p, 0.0)
}

///|
/// R binomial()$linkinv (C_logit_linkinv, NOT the unclamped plogis):
/// the inverse logit, clamped to [.Machine$double.eps, 1 - double.eps]
/// once the linear predictor leaves the ~[-36.04, 36.04] range.
fn _sw_logistic_mu(e : Double) -> Double {
  let eps = 2.220446049250313e-16
  let raw = if e >= 0.0 {
    1.0 / (1.0 + @math.exp(-e))
  } else {
    let ex = @math.exp(e)
    ex / (1.0 + ex)
  }
  if raw < eps {
    eps
  } else if raw > 1.0 - eps {
    1.0 - eps
  } else {
    raw
  }
}

///|
fn _sw_fit_logistic(
  y : Array[Double],
  zmat : Array[Array[Double]],
) -> Array[Double] {
  let n = y.length()
  let p = if n == 0 { 0 } else { zmat[0].length() }
  let eta = Array::make(n, 0.0)
  let mut i = 0
  while i < n {
    let mu0 = (y[i] + 0.5) / 2.0
    eta[i] = @math.ln(mu0 / (1.0 - mu0))
    i = i + 1
  }
  let mut mu = Array::make(n, 0.0)
  let mut dev : Double = 1.0 / 0.0
  let mut it = 0
  while it < 25 {
    let eta_new = Array::make(n, 0.0)
    let mu_new = Array::make(n, 0.0)
    i = 0
    while i < n {
      mu[i] = _sw_logistic_mu(eta[i])
      i = i + 1
    }
    let wz = _sw_mat(n, p)
    let wr = Array::make(n, 0.0)
    let mut dev_pre = 0.0
    i = 0
    while i < n {
      let mui = mu[i]
      let mueta = mui * (1.0 - mui)
      if it == 0 {
        // R seeds devold with the deviance at the starting mu=(y+.5)/2.
        dev_pre = dev_pre -
          2.0 * (y[i] * @math.ln(mui) + (1.0 - y[i]) * @math.ln(1.0 - mui))
      }
      // glm.fit drops observations with zero derivative weight; they simply
      // contribute nothing to the weighted least squares fit.
      if mueta > 0.0 {
        let work = eta[i] + (y[i] - mui) / mueta
        let sw = mueta.sqrt()
        let mut j = 0
        while j < p {
          wz[i][j] = zmat[i][j] * sw
          j = j + 1
        }
        wr[i] = work * sw
      }
      i = i + 1
    }
    if it == 0 {
      dev = dev_pre
    }
    let coef = _sw_wls_coef(wz, wr)
    let mut dev_post = 0.0
    i = 0
    while i < n {
      let mut s = 0.0
      let mut j = 0
      while j < p {
        s = s + zmat[i][j] * coef[j]
        j = j + 1
      }
      eta_new[i] = s
      let mnew = _sw_logistic_mu(s)
      mu_new[i] = mnew
      dev_post = dev_post -
        2.0 * (y[i] * @math.ln(mnew) + (1.0 - y[i]) * @math.ln(1.0 - mnew))
      i = i + 1
    }
    let conv = (dev - dev_post).abs() / (dev_post.abs() + 0.1) < 1.0e-8
    mu = mu_new
    dev = dev_post
    if conv {
      break
    }
    // eta = eta_new
    let mut k = 0
    while k < n {
      eta[k] = eta_new[k]
      k = k + 1
    }
    it = it + 1
  }
  mu
}

///|
/// Linear regression fitted values matching lsfit(X, y) (implicit intercept):
/// fitted = c0 + X c.
fn _sw_fit_linear(
  xmat : Array[Array[Double]],
  y : Array[Double],
) -> Array[Double] {
  let n = y.length()
  let p = if n == 0 { 0 } else { xmat[0].length() }
  let z = _sw_mat(n, p + 1)
  let mut i = 0
  while i < n {
    z[i][0] = 1.0
    let mut j = 0
    while j < p {
      z[i][j + 1] = xmat[i][j]
      j = j + 1
    }
    i = i + 1
  }
  let coef = _sw_wls_coef(z, y)
  let fit = Array::make(n, 0.0)
  i = 0
  while i < n {
    let mut s = coef[0]
    let mut j = 0
    while j < p {
      s = s + xmat[i][j] * coef[j + 1]
      j = j + 1
    }
    fit[i] = s
    i = i + 1
  }
  fit
}

// ===========================================================================
// Special functions: lgamma, beta density and regularized incomplete beta
// ===========================================================================

///|
/// log Gamma(x) via the Lanczos approximation (Spousteau coefficients).
fn _sw_lgamma(value : Double) -> Double {
  let coefficients = [
    0.9999999999998099, 676.5203681218851, -1259.1392167224028, 771.3234287776531,
    -176.6150291621406, 12.507343278686905, -0.13857109526572012, 9.984369578019572e-6,
    1.5056327351493116e-7,
  ]
  if value < 0.5 {
    let pi = 3.141592653589793
    return @math.ln(pi) -
      @math.ln(@math.sin(pi * value).abs()) -
      _sw_lgamma(1.0 - value)
  }
  let shifted = value - 1.0
  let mut series = coefficients[0]
  let mut index = 1
  while index < coefficients.length() {
    series = series + coefficients[index] / (shifted + index.to_double())
    index = index + 1
  }
  let t = shifted + 7.5
  0.9189385332046727 + (shifted + 0.5) * @math.ln(t) - t + @math.ln(series)
}

///|
/// log B(a,b) = lgamma(a) + lgamma(b) - lgamma(a+b).
fn _sw_lbeta(a : Double, b : Double) -> Double {
  _sw_lgamma(a) + _sw_lgamma(b) - _sw_lgamma(a + b)
}

///|
/// Beta density dbeta(x, a, b).
fn _sw_dbeta(x : Double, a : Double, b : Double) -> Double {
  if x <= 0.0 || x >= 1.0 {
    if (x == 0.0 && a == 1.0) || (x == 1.0 && b == 1.0) {
      // boundary density may be finite only for shape 1; R otherwise gives
      // 0 or Inf, neither of which occurs in the swfdr data path
      @math.exp(-_sw_lbeta(a, b))
    } else {
      0.0
    }
  } else {
    @math.exp(
      (a - 1.0) * @math.ln(x) + (b - 1.0) * @math.ln(1.0 - x) - _sw_lbeta(a, b),
    )
  }
}

///|
/// Modified Lentz continued fraction for the regularized incomplete beta
/// (Numerical Recipes betacf).
fn _sw_beta_cf(a : Double, b : Double, x : Double) -> Double {
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
  while m <= 200 {
    let m2 = 2.0 * m.to_double()
    let mut aa = m.to_double() *
      (b - m.to_double()) *
      x /
      ((qam + m2) * (a + m2))
    d = 1.0 + aa * d
    if d.abs() < fpmin {
      d = fpmin
    }
    c = 1.0 + aa / c
    if c.abs() < fpmin {
      c = fpmin
    }
    d = 1.0 / d
    h = h * d * c
    aa = -(a + m.to_double()) *
      (qab + m.to_double()) *
      x /
      ((a + m2) * (qap + m2))
    d = 1.0 + aa * d
    if d.abs() < fpmin {
      d = fpmin
    }
    c = 1.0 + aa / c
    if c.abs() < fpmin {
      c = fpmin
    }
    d = 1.0 / d
    let delta = d * c
    h = h * delta
    // The first few convergents can satisfy |delta-1| ~ 0 identically
    // (e.g. B == 1 makes the m = 1 odd numerator zero, so the following
    // even step gives d*c == 1 algebraically). Only trust the criterion
    // once enough terms have been accumulated.
    if m >= 6 && (delta - 1.0).abs() < 3.0e-12 {
      break
    }
    m = m + 1
  }
  h
}

///|
/// Regularized incomplete beta I_x(a,b) = pbeta(x, a, b).
fn _sw_pbeta(x : Double, a : Double, b : Double) -> Double {
  if x <= 0.0 {
    return 0.0
  }
  if x >= 1.0 {
    return 1.0
  }
  // bt already carries the 1/a prefactor of the direct branch; the
  // complement identity I_x(a,b) = 1 - I_{1-x}(b,a) divides by b
  // instead, hence the a/b rescale there.
  let lnb = _sw_lgamma(a + b) -
    _sw_lgamma(a) -
    _sw_lgamma(b) +
    a * @math.ln(x) +
    b * @math.ln(1.0 - x) -
    @math.ln(a)
  let bt = @math.exp(lnb)
  if x < (a + 1.0) / (a + b + 2.0) {
    bt * _sw_beta_cf(a, b, x)
  } else {
    1.0 - bt * (a / b) * _sw_beta_cf(b, a, 1.0 - x)
  }
}

// ===========================================================================
// lm_pi0
// ===========================================================================

///|
/// Default lambda grid seq(0.05, 0.95, 0.05) (19 values).
fn _sw_default_lambda() -> Array[Double] {
  let v : Array[Double] = []
  let mut i = 1
  while i <= 19 {
    v.push(i.to_double() * 0.05)
    i = i + 1
  }
  v
}

///|
/// Validate p, lambda, smooth.df and the covariate matrix, mirroring the R
/// check_* functions. Returns the (possibly all-ones) covariate matrix with
/// one row per p-value.
fn _sw_check_inputs(
  p : Array[Double],
  lambda : Array[Double],
  xopt : Array[Array[Double]]?,
  smooth_df : Int,
) -> Array[Array[Double]] raise SwfdrError {
  if p.length() == 0 {
    raise SwfdrError("p must be numeric in range [0, 1]\n")
  }
  let mut i = 0
  // R: range(p, na.rm=TRUE) is checked before the NA test
  let mut have_finite = false
  let mut pmin = 0.0
  let mut pmax = 0.0
  while i < p.length() {
    let v = p[i]
    if v == v && !v.is_nan() {
      if !have_finite {
        pmin = v
        pmax = v
        have_finite = true
      } else {
        if v < pmin {
          pmin = v
        }
        if v > pmax {
          pmax = v
        }
      }
    }
    i = i + 1
  }
  if have_finite && (pmin < 0.0 || pmax > 1.0) {
    raise SwfdrError("p must be numeric in range [0, 1]\n")
  }
  i = 0
  while i < p.length() {
    if p[i] != p[i] {
      raise SwfdrError("p must not have missing values\n")
    }
    i = i + 1
  }
  if lambda.length() < 4 {
    raise SwfdrError("lambda must be of length >=4\n")
  }
  let mut lmin = lambda[0]
  let mut lmax = lambda[0]
  i = 1
  while i < lambda.length() {
    if lambda[i] < lmin {
      lmin = lambda[i]
    }
    if lambda[i] > lmax {
      lmax = lambda[i]
    }
    i = i + 1
  }
  if lmin < 0.0 || lmax >= 1.0 {
    raise SwfdrError("lambda values must all be in range [0, 1)\n")
  }
  if smooth_df <= 1 || smooth_df > lambda.length() {
    raise SwfdrError("df must be in range 1 < df < length(lambda)\n")
  }
  let xmat = match xopt {
    Some(x) => x
    None => {
      // X missing: warning in R, modeling then has no effect
      let ones : Array[Array[Double]] = []
      for _ in 0..<p.length() {
        ones.push([1.0])
      }
      ones
    }
  }
  if xmat.length() != p.length() {
    raise SwfdrError("incompatible X and p - different lengths\n")
  }
  i = 0
  while i < xmat.length() {
    let mut j = 0
    while j < xmat[i].length() {
      let v = xmat[i][j]
      if v.is_nan() || v.abs() > 1.0e300 {
        raise SwfdrError("X must not contain missing or non-finite values\n")
      }
      j = j + 1
    }
    i = i + 1
  }
  // pmax < max(lambda) only produces a warning in R; computation proceeds.
  xmat
}

///|
/// Sort and de-duplicate lambda (R sort(unique(lambda))).
fn _sw_sort_unique(x : Array[Double]) -> Array[Double] {
  let v = x.copy()
  v.sort()
  let out : Array[Double] = []
  let mut i = 0
  while i < v.length() {
    if out.length() == 0 || (out[out.length() - 1] - v[i]).abs() > 1.0e-14 {
      out.push(v[i])
    }
    i = i + 1
  }
  out
}

///|
/// Build a covariate matrix (one column per supplied vector).
pub fn swfdr_design_matrix(
  columns : Array[Array[Double]],
) -> Array[Array[Double]]? {
  if columns.length() == 0 {
    return None
  }
  let n = columns[0].length()
  let m = _sw_mat(n, columns.length())
  let mut j = 0
  while j < columns.length() {
    if columns[j].length() != n {
      return None
    }
    let mut i = 0
    while i < n {
      m[i][j] = columns[j][i]
      i = i + 1
    }
    j = j + 1
  }
  Some(m)
}

///|
/// Estimate pi0(x) from p-values and covariates, porting the Bioconductor
/// `lm_pi0` function.
///
/// At each lambda, y = (p >= lambda) is regressed on cbind(1, X) using
/// logistic (glm binomial) or linear (lsfit) regression, and divided by
/// (1 - lambda). The resulting pi0(lambda) curve is then smoothed over
/// lambda separately for every test, either with a `splines::bs` penalized
/// unit spline or an R-compatible smooth.spline.
///
/// Parameters:
///   - p:          p-values in [0, 1]
///   - lambda:     threshold grid (default seq(.05, .95, .05))
///   - X:          covariate matrix with one row per test (None = no
///                 covariates; a single all-ones column is used)
///   - type:       "logistic" (default) or "linear"
///   - smooth_df:  degrees of freedom for the smoother (default 3)
///   - threshold:  truncate estimates into [0, 1] (default true)
///   - smoothing:  "unit.spline" (default) or "smooth.spline"
pub fn swfdr_lm_pi0(
  p : Array[Double],
  lambda? : Array[Double]? = None,
  x? : Array[Array[Double]]? = None,
  type_? : String = "logistic",
  smooth_df? : Int = 3,
  threshold? : Bool = true,
  smoothing? : String = "unit.spline",
) -> SwfdrPi0 raise SwfdrError {
  let lams0 = match lambda {
    Some(v) => v
    None => _sw_default_lambda()
  }
  let lams = _sw_sort_unique(lams0)
  let df = smooth_df.to_double().round().to_int()
  let xmat = _sw_check_inputs(p, lams, x, df)
  if type_ != "logistic" && type_ != "linear" {
    raise SwfdrError("type must be \"logistic\" or \"linear\"\n")
  }
  if smoothing != "unit.spline" && smoothing != "smooth.spline" {
    raise SwfdrError("smoothing must be \"unit.spline\" or \"smooth.spline\"\n")
  }
  let n = p.length()
  let nl = lams.length()
  // design with intercept
  let p1 = if n == 0 { 0 } else { xmat[0].length() }
  let zmat = _sw_mat(n, p1 + 1)
  let mut i = 0
  while i < n {
    zmat[i][0] = 1.0
    let mut j = 0
    while j < p1 {
      zmat[i][j + 1] = xmat[i][j]
      j = j + 1
    }
    i = i + 1
  }
  let pi0_lambda = _sw_mat(n, nl)
  let mut li = 0
  while li < nl {
    let y = Array::make(n, 0.0)
    let mut k = 0
    while k < n {
      if p[k] >= lams[li] {
        y[k] = 1.0
      }
      k = k + 1
    }
    let fitted = if type_ == "logistic" {
      _sw_fit_logistic(y, zmat)
    } else {
      _sw_fit_linear(xmat, y)
    }
    let denom = 1.0 - lams[li]
    k = 0
    while k < n {
      pi0_lambda[k][li] = fitted[k] / denom
      k = k + 1
    }
    li = li + 1
  }
  if threshold {
    let mut r = 0
    while r < n {
      let mut c = 0
      while c < nl {
        if pi0_lambda[r][c] < 0.0 {
          pi0_lambda[r][c] = 0.0
        } else if pi0_lambda[r][c] > 1.0 {
          pi0_lambda[r][c] = 1.0
        }
        c = c + 1
      }
      r = r + 1
    }
  }
  let mut pi0_smooth : Array[Array[Double]] = []
  let pi0 = if smoothing == "smooth.spline" {
    pi0_smooth = _sw_mat(n, nl)
    let est = Array::make(n, 0.0)
    let mut r = 0
    while r < n {
      let fit = _sw_smooth_spline(lams, pi0_lambda[r], df.to_double())
      let mut c = 0
      while c < nl {
        let mut v = fit[c]
        if threshold {
          if v < 0.0 {
            v = 0.0
          } else if v > 1.0 {
            v = 1.0
          }
        }
        pi0_smooth[r][c] = v
        c = c + 1
      }
      est[r] = fit[nl - 1]
      r = r + 1
    }
    est
  } else {
    let tmat = swfdr_unit_spline_matrix(lams, df)
    let lastrow = tmat[nl - 1]
    let est = Array::make(n, 0.0)
    let mut r = 0
    while r < n {
      let mut s = 0.0
      let mut c = 0
      while c < nl {
        s = s + lastrow[c] * pi0_lambda[r][c]
        c = c + 1
      }
      est[r] = s
      r = r + 1
    }
    est
  }
  if threshold {
    let mut r = 0
    while r < n {
      if pi0[r] < 0.0 {
        pi0[r] = 0.0
      } else if pi0[r] > 1.0 {
        pi0[r] = 1.0
      }
      r = r + 1
    }
  }
  SwfdrPi0::{
    lambda: lams,
    x_names: [],
    pi0_lambda,
    pi0,
    pi0_smooth,
    smoothing,
    type_,
  }
}

// ===========================================================================
// lm_qvalue
// ===========================================================================

///|
/// Compute covariate-conditioned q-values, porting the Bioconductor
/// `lm_qvalue` function.
///
/// In decreasing-p order the running-minimum term p*n/rank (ordinary FDR) or
/// p*n/(rank*(1-(1-p)^n)) (pFDR finite-sample correction) is multiplied by
/// the covariate-specific pi0; the result is mapped back to the original
/// order. No final truncation at 1 is applied after the pi0 multiplication.
///
/// Parameters:
///   - p:          p-values in [0, 1]
///   - X:          covariate matrix (None if pi0 is supplied)
///   - pfdr:       use the pFDR correction (default false)
///   - pi0:        a precomputed SwfdrPi0 (estimated if omitted)
///   - smoothing:  smoothing method when pi0 is estimated
///   - lambda:     lambda grid when pi0 is estimated
///   - smooth_df:  smoother degrees of freedom
///   - threshold:  truncate pi0 estimates into [0, 1]
pub fn swfdr_lm_qvalue(
  p : Array[Double],
  x? : Array[Array[Double]]? = None,
  pfdr? : Bool = false,
  pi0? : SwfdrPi0? = None,
  smoothing? : String = "unit.spline",
  lambda? : Array[Double]? = None,
  smooth_df? : Int = 3,
  threshold? : Bool = true,
) -> SwfdrQvalue raise SwfdrError {
  let pi0_result = match pi0 {
    Some(v) => v
    None => swfdr_lm_pi0(p, lambda~, x~, smoothing~, smooth_df~, threshold~)
  }
  let n = p.length()
  // stable order: decreasing p, increasing original index for ties
  let order : Array[Int] = []
  let mut i = 0
  while i < n {
    order.push(i)
    i = i + 1
  }
  order.sort_by(fn(a : Int, b : Int) -> Int {
    if p[a] > p[b] {
      -1
    } else if p[a] < p[b] {
      1
    } else if a < b {
      -1
    } else if a > b {
      1
    } else {
      0
    }
  })
  let nd = n.to_double()
  let q = Array::make(n, 0.0)
  let mut running = 1.0e300
  let mut pos = 0
  while pos < n {
    let idx = order[pos]
    let rank = nd - pos.to_double() // i = n, n-1, ..., 1
    let mut term = p[idx] * nd / rank
    if pfdr {
      let onemp = 1.0 - @math.pow(1.0 - p[idx], nd)
      term = p[idx] * nd / (rank * onemp)
    }
    if term < running {
      running = term
    }
    let capped = if running > 1.0 { 1.0 } else { running }
    q[idx] = capped * pi0_result.pi0[idx]
    pos = pos + 1
  }
  SwfdrQvalue::{ pvalues: p.copy(), qvalues: q, pi0_result, pfdr, }
}

// ===========================================================================
// calculateSwfdr: EM with rounded/truncated p-values
// ===========================================================================

///|
/// Bin the rounded p-values into (0,.005], (.005,.015], ..., (.045,.05].
/// Values exactly 0 are moved to 1e-10 (first bin); values exactly 0.05 are
/// moved below 0 (dropped), reproducing the R source quirks.
fn _sw_rounded_counts(pp : Array[Double], rr : Array[Int]) -> Array[Int] {
  let edges = [0.0, 0.005, 0.015, 0.025, 0.035, 0.045, 0.05]
  let counts = Array::make(6, 0)
  let mut i = 0
  while i < pp.length() {
    if rr[i] > 0 {
      let mut v = pp[i]
      if v == 0.0 {
        v = 1.0e-10
      }
      if v == 0.05 {
        v = 0.05 - 1.0e10 - 10.0
      }
      if v > 0.0 && v <= 0.05 {
        let mut b = 0
        while b < 6 {
          if v > edges[b] && v <= edges[b + 1] {
            counts[b] = counts[b] + 1
            b = 6
          }
          b = b + 1
        }
      }
    }
    i = i + 1
  }
  counts
}

///|
/// Negative log-likelihood minimized in the M-step (R ll(a,b)).
fn _sw_em_negll(
  pp : Array[Double],
  tt : Array[Int],
  rr : Array[Int],
  z : Array[Double?],
  n1 : Array[Double],
  a : Double,
  b : Double,
) -> Double {
  let cut = _sw_pbeta(0.05, a, b)
  let mut tmp1 = 0.0
  let mut i = 0
  while i < pp.length() {
    let mut t = 0.0
    if rr[i] == 0 {
      if tt[i] == 0 {
        let d = _sw_dbeta(pp[i], a, b)
        t = @math.ln(if d / cut > 0.0 { d / cut } else { 1.0e-300 })
      } else {
        let pb = _sw_pbeta(pp[i], a, b)
        t = @math.ln(if pb / cut > 0.0 { pb / cut } else { 1.0e-300 })
      }
      match z[i] {
        Some(zv) => tmp1 = tmp1 - (1.0 - zv) * t
        None => ()
      }
    }
    i = i + 1
  }
  // bin alternative probabilities, bin1 .. bin6
  let upper = [0.05, 0.045, 0.035, 0.025, 0.015, 0.005]
  let lower = [0.045, 0.035, 0.025, 0.015, 0.005, 0.0]
  let prob = Array::make(6, 0.0)
  let mut k = 0
  while k < 6 {
    prob[5 - k] = (_sw_pbeta(upper[k], a, b) - _sw_pbeta(lower[k], a, b)) / cut
    k = k + 1
  }
  let mut tmp2 = 0.0
  k = 0
  while k < 6 {
    let pv = if prob[k] > 0.0 { prob[k] } else { 1.0e-300 }
    tmp2 = tmp2 - n1[k] * @math.ln(pv)
    k = k + 1
  }
  tmp1 + tmp2
}

///|
/// M-step bounded optimization, exactly mirroring
/// `stats4::mle(ll, start = list(a = 0.05, b = 100), lower = c(0.001, 1),
/// upper = c(1, 500), method = "L-BFGS-B")`: R's L-BFGS-B v2.3 with
/// factr = 1e7 and bounded central-difference gradients (ndeps = 1e-3).
/// The early-stopping trajectory of the algorithm matters for the EM
/// iterates, so the full L-BFGS-B machinery is ported in swfdr_lbfgs.mbt.
fn _sw_em_mstep(
  pp : Array[Double],
  tt : Array[Int],
  rr : Array[Int],
  z : Array[Double?],
  n1 : Array[Double],
) -> (Double, Double) {
  let x0 = [0.05, 100.0]
  let lo = [0.001, 1.0]
  let up = [1.0, 500.0]
  let nbd = [2, 2]
  let obj = fn(x : Array[Double]) -> Double {
    _sw_em_negll(pp, tt, rr, z, n1, x[0], x[1])
  }
  let (x, _f, _nfgv) = _sw_lbfgsb_minimize(x0, lo, up, nbd, obj)
  (x[0], x[1])
}

///|
/// Science-wise FDR estimation, porting the Bioconductor `calculateSwfdr`
/// function.
///
/// An EM algorithm fits pi0 (the prior null probability) and the
/// Beta(alpha, beta) distribution of alternative p-values, accounting for
/// rounded p-values (binned between 0 and 0.05) and truncated p-values
/// (reported only as upper bounds).
///
/// Parameters:
///   - p_values:           p-values
///   - truncated:          0/1 array, 1 = p-value is truncated (an upper bound)
///   - rounded:            0/1 array, 1 = p-value is rounded
///   - pi0:                initial pi0 (default 0.5)
///   - alpha:              initial beta shape 1 (default 1)
///   - beta:               initial beta shape 2 (default 50)
///   - num_em_iterations:  EM iterations (default 100)
pub fn swfdr_calculate(
  p_values : Array[Double],
  truncated : Array[Int],
  rounded : Array[Int],
  pi0? : Double = 0.5,
  alpha? : Double = 1.0,
  beta? : Double = 50.0,
  num_em_iterations? : Int = 100,
) -> SwfdrScience raise SwfdrError {
  let pp = p_values
  let tt = truncated
  let rr = rounded
  let m = pp.length()
  if tt.length() != m || rr.length() != m {
    raise SwfdrError("truncated and rounded must match p-values in length\n")
  }
  let counts = _sw_rounded_counts(pp, rr)
  let prob0 = [0.1, 0.2, 0.2, 0.2, 0.2, 0.1]
  let mut pi0v = pi0
  let mut alpha_v = alpha
  let mut beta_v = beta
  let z : Array[Double?] = Array::make(m, None)
  let n0 = Array::make(6, 0.0)
  let n1 = Array::make(6, 0.0)
  let mut it = 0
  while it < num_em_iterations {
    // ----- E-step -----
    let cut = _sw_pbeta(0.05, alpha_v, beta_v)
    let upper = [0.05, 0.045, 0.035, 0.025, 0.015, 0.005]
    let lower = [0.045, 0.035, 0.025, 0.015, 0.005, 0.0]
    let prob1 = Array::make(6, 0.0)
    let mut k = 0
    while k < 6 {
      prob1[5 - k] = (
          _sw_pbeta(upper[k], alpha_v, beta_v) -
          _sw_pbeta(lower[k], alpha_v, beta_v)
        ) /
        cut
      k = k + 1
    }
    let mut sum_n0 = 0.0
    let mut sum_n = 0.0
    k = 0
    while k < 6 {
      let pij0 = pi0v * prob0[k] / (prob0[k] * pi0v + prob1[k] * (1.0 - pi0v))
      n0[k] = counts[k].to_double() * pij0
      n1[k] = counts[k].to_double() - n0[k]
      sum_n0 = sum_n0 + n0[k]
      sum_n = sum_n + counts[k].to_double()
      k = k + 1
    }
    // z only for non-rounded p-values
    let mut sum_z = 0.0
    let mut n_nonrounded = 0
    let mut i = 0
    while i < m {
      if rr[i] == 0 {
        n_nonrounded = n_nonrounded + 1
        let zv = if tt[i] == 0 {
          let alt = _sw_dbeta(pp[i], alpha_v, beta_v) / cut
          pi0v * 20.0 / (pi0v * 20.0 + (1.0 - pi0v) * alt)
        } else {
          let alt = _sw_pbeta(pp[i], alpha_v, beta_v) / cut
          pi0v * 20.0 * pp[i] / (pi0v * 20.0 * pp[i] + (1.0 - pi0v) * alt)
        }
        z[i] = Some(zv)
        sum_z = sum_z + zv
      }
      i = i + 1
    }
    // ----- M-step -----
    pi0v = (sum_n0 + sum_z) / (sum_n + n_nonrounded.to_double())
    let (na, nb) = _sw_em_mstep(pp, tt, rr, z, n1)
    alpha_v = na
    beta_v = nb
    it = it + 1
  }
  SwfdrScience::{ pi0: pi0v, alpha: alpha_v, beta: beta_v, z, n0, n: counts, }
}
