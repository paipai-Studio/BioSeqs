///|
/// L-BFGS-B bounded minimizer.
///
/// Faithful 0-based MoonBit port of R's bundled L-BFGS-B v2.3
/// (`src/appl/lbfgsb.c`, Zhu/Byrd/Lu-Chen/Nocedal) together with:
///   - the Linpack helpers `dpofa` (Cholesky, with R's added 1e-14 positive
///     definiteness tolerance) and `dtrsl` (job 1: upper back-substitution,
///     job 11: upper transpose forward-substitution);
///   - the MINPACK-2 More'-Thuente line search (`dcsrch`/`dcstep`);
///   - the driver loop in `src/appl/optim.c` (`lbfgsb()`);
///   - R's bounded central-difference gradient (`optim_stats.c`, ndeps 1e-3).
///
/// The R `optim(..., method = "L-BFGS-B")` defaults are reproduced:
/// lmm = 5, factr = 1e7, pgtol = 0, maxit = 100.  The swfdr EM M-step
/// relies on the exact early-stopping trajectory of this algorithm, so a
/// generic coordinate optimizer is not an acceptable substitute.

// ===========================================================================
// BLAS / Linpack helpers (column-major, 0-based, unit increments)
// ===========================================================================

fn _lb_ddot(n : Int, x : Array[Double], x0 : Int, y : Array[Double], y0 : Int) -> Double {
  let mut s = 0.0
  let mut k = 0
  while k < n {
    s += x[x0 + k] * y[y0 + k]
    k += 1
  }
  s
}

fn _lb_dcopy(n : Int, x : Array[Double], x0 : Int, y : Array[Double], y0 : Int) -> Unit {
  let mut k = 0
  while k < n {
    y[y0 + k] = x[x0 + k]
    k += 1
  }
}

fn _lb_daxpy(
  n : Int,
  alpha : Double,
  x : Array[Double],
  x0 : Int,
  y : Array[Double],
  y0 : Int,
) -> Unit {
  let mut k = 0
  while k < n {
    y[y0 + k] += alpha * x[x0 + k]
    k += 1
  }
}

fn _lb_dscal(n : Int, alpha : Double, x : Array[Double], x0 : Int) -> Unit {
  let mut k = 0
  while k < n {
    x[x0 + k] *= alpha
    k += 1
  }
}

///|
/// Linpack dpofa: factor a symmetric positive definite matrix into an upper
/// triangular factor r with a = r'r, overwriting the upper triangle.
/// Returns 0 on success or the 1-based failing column.
fn _lb_dpofa(a : Array[Double], a0 : Int, lda : Int, n : Int) -> Int {
  let mut j = 0
  while j < n {
    let mut s = 0.0
    if j > 0 {
      let mut k = 0
      while k < j {
        let mut t =
          a[a0 + k + j * lda] - _lb_ddot(k, a, a0 + k * lda, a, a0 + j * lda)
        t /= a[a0 + k + k * lda]
        a[a0 + k + j * lda] = t
        s += t * t
        k += 1
      }
    }
    s = a[a0 + j + j * lda] - s
    if s <= 0.0 {
      return j + 1
    }
    a[a0 + j + j * lda] = s.sqrt()
    j += 1
  }
  0
}

///|
/// dtrsl job 11: solve trans(t) x = b for upper triangular t (forward
/// substitution, reading the upper triangle). In place on b.
/// Returns 0 on success or the 1-based index of a zero diagonal.
fn _lb_dtrsl_11(
  t : Array[Double],
  t0 : Int,
  lda : Int,
  n : Int,
  b : Array[Double],
  b0 : Int,
) -> Int {
  let mut j = 0
  while j < n {
    if t[t0 + j + j * lda] == 0.0 {
      return j + 1
    }
    j += 1
  }
  b[b0] /= t[t0]
  let mut j = 1
  while j < n {
    b[b0 + j] = (b[b0 + j] - _lb_ddot(j, t, t0 + j * lda, b, b0)) / t[t0 + j + j * lda]
    j += 1
  }
  0
}

///|
/// dtrsl job 01: solve t x = b for upper triangular t (back substitution,
/// reading the upper triangle). In place on b.
fn _lb_dtrsl_01(
  t : Array[Double],
  t0 : Int,
  lda : Int,
  n : Int,
  b : Array[Double],
  b0 : Int,
) -> Int {
  let mut j = 0
  while j < n {
    if t[t0 + j + j * lda] == 0.0 {
      return j + 1
    }
    j += 1
  }
  b[b0 + n - 1] /= t[t0 + n - 1 + (n - 1) * lda]
  let mut jj = 2
  while jj <= n {
    let j = n - jj
    _lb_daxpy(j + 1, -b[b0 + j + 1], t, t0 + (j + 1) * lda, b, b0)
    b[b0 + j] /= t[t0 + j + j * lda]
    jj += 1
  }
  0
}

// ===========================================================================
// More'-Thuente line search (MINPACK-2 dcsrch / dcstep)
// ===========================================================================

struct LbInterval {
  mut stx : Double
  mut fx : Double
  mut dx : Double
  mut sty : Double
  mut fy : Double
  mut dy : Double
  mut stp : Double
  mut brackt : Bool
  stpmin : Double
  stpmax : Double
}

///|
/// dcstep: safeguarded cubic/secant step update.
fn _lb_dcstep(v : LbInterval, fp : Double, dp : Double) -> Unit {
  let sgnd = dp * (v.dx / v.dx.abs())
  let mut stpf : Double = 0
  if fp > v.fx {
    // First case: higher function value.
    let theta = (v.fx - fp) * 3.0 / (v.stp - v.stx) + v.dx + dp
    let s = theta.abs().max(v.dx.abs()).max(dp.abs())
    let mut gamm = s * ((theta / s) * (theta / s) - v.dx / s * (dp / s)).sqrt()
    if v.stp < v.stx {
      gamm = -gamm
    }
    let p = gamm - v.dx + theta
    let q = gamm - v.dx + gamm + dp
    let r__ = p / q
    let stpc = v.stx + r__ * (v.stp - v.stx)
    let stpq =
      v.stx + v.dx / ((v.fx - fp) / (v.stp - v.stx) + v.dx) / 2.0 * (v.stp - v.stx)
    if (stpc - v.stx).abs() < (stpq - v.stx).abs() {
      stpf = stpc
    } else {
      stpf = stpc + (stpq - stpc) / 2.0
    }
    v.brackt = true
  } else if sgnd < 0.0 {
    // Second case: lower value, opposite-sign derivatives.
    let theta = (v.fx - fp) * 3.0 / (v.stp - v.stx) + v.dx + dp
    let s = theta.abs().max(v.dx.abs()).max(dp.abs())
    let mut gamm = s * ((theta / s) * (theta / s) - v.dx / s * (dp / s)).sqrt()
    if v.stp > v.stx {
      gamm = -gamm
    }
    let p = gamm - dp + theta
    let q = gamm - dp + gamm + v.dx
    let r__ = p / q
    let stpc = v.stp + r__ * (v.stx - v.stp)
    let stpq = v.stp + dp / (dp - v.dx) * (v.stx - v.stp)
    if (stpc - v.stp).abs() > (stpq - v.stp).abs() {
      stpf = stpc
    } else {
      stpf = stpq
    }
    v.brackt = true
  } else if dp.abs() < v.dx.abs() {
    // Third case: lower value, same-sign derivatives, |g| decreased.
    let theta = (v.fx - fp) * 3.0 / (v.stp - v.stx) + v.dx + dp
    let s = theta.abs().max(v.dx.abs()).max(dp.abs())
    let mut d1 = theta / s
    d1 = d1 * d1 - v.dx / s * (dp / s)
    let mut gamm = if d1 < 0.0 { 0.0 } else { s * d1.sqrt() }
    if v.stp > v.stx {
      gamm = -gamm
    }
    let p = gamm - dp + theta
    let q = gamm + (v.dx - dp) + gamm
    let r__ = p / q
    let stpc = if r__ < 0.0 && gamm != 0.0 {
      v.stp + r__ * (v.stx - v.stp)
    } else if v.stp > v.stx {
      v.stpmax
    } else {
      v.stpmin
    }
    let stpq = v.stp + dp / (dp - v.dx) * (v.stx - v.stp)
    if v.brackt {
      if (stpc - v.stp).abs() < (stpq - v.stp).abs() {
        stpf = stpc
      } else {
        stpf = stpq
      }
      let bound = v.stp + (v.sty - v.stp) * 0.66
      if v.stp > v.stx {
        stpf = bound.min(stpf)
      } else {
        stpf = bound.max(stpf)
      }
    } else {
      if (stpc - v.stp).abs() > (stpq - v.stp).abs() {
        stpf = stpc
      } else {
        stpf = stpq
      }
      stpf = v.stpmax.min(stpf)
      stpf = v.stpmin.max(stpf)
    }
  } else {
    // Fourth case.
    if v.brackt {
      let theta = (fp - v.fy) * 3.0 / (v.sty - v.stp) + v.dy + dp
      let s = theta.abs().max(v.dy.abs()).max(dp.abs())
      let mut gamm = s * ((theta / s) * (theta / s) - v.dy / s * (dp / s)).sqrt()
      if v.stp > v.sty {
        gamm = -gamm
      }
      let p = gamm - dp + theta
      let q = gamm - dp + gamm + v.dy
      let r__ = p / q
      stpf = v.stp + r__ * (v.sty - v.stp)
    } else if v.stp > v.stx {
      stpf = v.stpmax
    } else {
      stpf = v.stpmin
    }
  }
  // Update the bracketing interval.
  if fp > v.fx {
    v.sty = v.stp
    v.fy = fp
    v.dy = dp
  } else {
    if sgnd < 0.0 {
      v.sty = v.stx
      v.fy = v.fx
      v.dy = v.dx
    }
    v.stx = v.stp
    v.fx = fp
    v.dx = dp
  }
  v.stp = stpf
}

///|
/// dcsrch state. status: 1 = FG (evaluate f,g at stp), 0 = CONVERGENCE,
/// 2 = WARNING.
struct LbDcsrch {
  mut stage : Int
  mut brackt : Bool
  mut finit : Double
  mut ginit : Double
  mut gtest : Double
  mut stx : Double
  mut fx : Double
  mut gx : Double
  mut sty : Double
  mut fy : Double
  mut gy : Double
  mut stmin : Double
  mut stmax : Double
  mut width : Double
  mut width1 : Double
  mut status : Int
  ftol : Double
  gtol : Double
  xtol : Double
  stpmin : Double
  stpmax : Double
}

fn _lb_dcsrch_start(stpmax : Double) -> LbDcsrch {
  {
    stage: 1,
    brackt: false,
    finit: 0.0,
    ginit: 0.0,
    gtest: 0.0,
    stx: 0.0,
    fx: 0.0,
    gx: 0.0,
    sty: 0.0,
    fy: 0.0,
    gy: 0.0,
    stmin: 0.0,
    stmax: 0.0,
    width: stpmax,
    width1: stpmax / 0.5,
    status: 1,
    ftol: 0.001,
    gtol: 0.9,
    xtol: 0.1,
    stpmin: 0.0,
    stpmax,
  }
}

///|
/// One dcsrch call. The first invocation uses task START (f,g at step 0,
/// proposed stp); later calls pass the values evaluated at the returned
/// stp. The new stp is returned; status is stored in self.status.
fn _lb_dcsrch_apply(
  self : LbDcsrch,
  start : Bool,
  f : Double,
  g : Double,
  stp_in : Double,
) -> Double {
  if start {
    let stpv = stp_in
    if stpv < self.stpmin || stpv > self.stpmax || g >= 0.0 {
      self.status = 2
      return stpv
    }
    self.brackt = false
    self.stage = 1
    self.finit = f
    self.ginit = g
    self.gtest = self.ftol * g
    self.width = self.stpmax - self.stpmin
    self.width1 = self.width / 0.5
    self.stx = 0.0
    self.fx = self.finit
    self.gx = self.ginit
    self.sty = 0.0
    self.fy = self.finit
    self.gy = self.ginit
    self.stmin = 0.0
    self.stmax = stpv + stpv * 4.0
    self.status = 1
    return stpv
  }
  let mut stpv = stp_in
  let ftest = self.finit + stpv * self.gtest
  if self.stage == 1 && f <= ftest && g >= 0.0 {
    self.stage = 2
  }
  let mut warn = false
  if self.brackt && (stpv <= self.stmin || stpv >= self.stmax) {
    warn = true
  }
  if self.brackt && self.stmax - self.stmin <= self.xtol * self.stmax {
    warn = true
  }
  if stpv == self.stpmax && f <= ftest && g <= self.gtest {
    warn = true
  }
  if stpv == self.stpmin && (f > ftest || g >= self.gtest) {
    warn = true
  }
  if f <= ftest && g.abs() <= self.gtol * (-self.ginit) {
    self.status = 0
    return(stpv)
  }
  if warn {
    self.status = 2
    return(stpv)
  }
  if self.stage == 1 && f <= self.fx && f > ftest {
    let fm = f - stpv * self.gtest
    let fxm = self.fx - self.stx * self.gtest
    let fym = self.fy - self.sty * self.gtest
    let gm = g - self.gtest
    let gxm = self.gx - self.gtest
    let gym = self.gy - self.gtest
    let iv = {
      stx: self.stx,
      fx: fxm,
      dx: gxm,
      sty: self.sty,
      fy: fym,
      dy: gym,
      stp: stpv,
      brackt: self.brackt,
      stpmin: self.stmin,
      stpmax: self.stmax,
    }
    _lb_dcstep(iv, fm, gm)
    self.stx = iv.stx
    self.fx = iv.fx + iv.stx * self.gtest
    self.fy = iv.fy + iv.sty * self.gtest
    self.gx = iv.dx + self.gtest
    self.gy = iv.dy + self.gtest
    stpv = iv.stp
    self.brackt = iv.brackt
  } else {
    let iv = {
      stx: self.stx,
      fx: self.fx,
      dx: self.gx,
      sty: self.sty,
      fy: self.fy,
      dy: self.gy,
      stp: stpv,
      brackt: self.brackt,
      stpmin: self.stmin,
      stpmax: self.stmax,
    }
    _lb_dcstep(iv, f, g)
    self.stx = iv.stx
    self.fx = iv.fx
    self.gx = iv.dx
    self.sty = iv.sty
    self.fy = iv.fy
    self.gy = iv.dy
    stpv = iv.stp
    self.brackt = iv.brackt
  }
  if self.brackt {
    if (self.sty - self.stx).abs() >= self.width1 * 0.66 {
      stpv = self.stx + (self.sty - self.stx) * 0.5
    }
    self.width1 = self.width
    self.width = (self.sty - self.stx).abs()
  }
  if self.brackt {
    self.stmin = self.stx.min(self.sty)
    self.stmax = self.stx.max(self.sty)
  } else {
    self.stmin = stpv + (stpv - self.stx) * 1.1
    self.stmax = stpv + (stpv - self.stx) * 4.0
  }
  if stpv < self.stpmin {
    stpv = self.stpmin
  }
  if stpv > self.stpmax {
    stpv = self.stpmax
  }
  let narrow =
    self.brackt &&
    (stpv <= self.stmin || stpv >= self.stmax ||
      self.stmax - self.stmin <= self.xtol * self.stmax)
  if narrow {
    stpv = self.stx
  }
  self.status = 1
  stpv
}

// ===========================================================================
// L-BFGS-B main solver
// ===========================================================================

struct Lbfgsb {
  n : Int
  m : Int
  factr : Double
  pgtol : Double
  ndeps : Double
  // static mainlb state
  mut prjctd : Bool
  mut cnstnd : Bool
  mut boxed : Bool
  mut updatd : Bool
  mut nintol : Int
  mut iback : Int
  mut nskip : Int
  mut head : Int // 1-based column pointer
  mut col : Int
  mut itail : Int
  mut iter : Int
  mut iupdat : Int
  mut nint : Int
  mut nfgv : Int
  mut info : Int
  mut ifun : Int
  mut iword : Int
  mut nfree : Int
  mut nact : Int
  mut ileave : Int // 1-based Fortran value (n+1 ..)
  mut nenter : Int
  mut theta : Double
  mut fold : Double
  mut tol : Double
  mut dnorm : Double
  epsmch : Double
  mut gd : Double
  mut stpmx : Double
  mut sbgnrm : Double
  mut stp : Double
  mut gdold : Double
  mut dtd : Double
  // working storage
  ws : Array[Double]
  wy : Array[Double]
  sy : Array[Double]
  ss : Array[Double]
  wt : Array[Double]
  wn : Array[Double]
  snd : Array[Double]
  z : Array[Double]
  r : Array[Double]
  d : Array[Double]
  t : Array[Double]
  wa : Array[Double]
  indx : Array[Int]
  iwhere : Array[Int]
  indx2 : Array[Int]
}

fn _lbfgsb_new(n : Int, m : Int, factr : Double, pgtol : Double, ndeps : Double) -> Lbfgsb {
  {
    n,
    m,
    factr,
    pgtol,
    ndeps,
    prjctd: false,
    cnstnd: false,
    boxed: true,
    updatd: false,
    nintol: 0,
    iback: 0,
    nskip: 0,
    head: 1,
    col: 0,
    itail: 0,
    iter: 0,
    iupdat: 0,
    nint: 0,
    nfgv: 0,
    info: 0,
    ifun: 0,
    iword: 0,
    nfree: n,
    nact: 0,
    ileave: 0,
    nenter: 0,
    theta: 1.0,
    fold: 0.0,
    tol: factr * 2.220446049250313e-16,
    dnorm: 0.0,
    epsmch: 2.220446049250313e-16,
    gd: 0.0,
    stpmx: 0.0,
    sbgnrm: 0.0,
    stp: 0.0,
    gdold: 0.0,
    dtd: 0.0,
    ws: Array::make(n * m, 0.0),
    wy: Array::make(n * m, 0.0),
    sy: Array::make(m * m, 0.0),
    ss: Array::make(m * m, 0.0),
    wt: Array::make(m * m, 0.0),
    wn: Array::make(4 * m * m, 0.0),
    snd: Array::make(4 * m * m, 0.0),
    z: Array::make(n, 0.0),
    r: Array::make(n, 0.0),
    d: Array::make(n, 0.0),
    t: Array::make(n, 0.0),
    wa: Array::make(8 * m, 0.0),
    indx: Array::make(n, 0),
    iwhere: Array::make(n, 0),
    indx2: Array::make(n, 0),
  }
}

///|
/// active: initialize iwhere and project initial x into the feasible box.
fn _lbfgsb_active(
  self : Lbfgsb,
  l : Array[Double],
  u : Array[Double],
  nbd : Array[Int],
  x : Array[Double],
) -> Unit {
  let n = self.n
  self.prjctd = false
  self.cnstnd = false
  self.boxed = true
  let mut i = 0
  while i < n {
    if nbd[i] > 0 {
      if nbd[i] <= 2 && x[i] <= l[i] {
        if x[i] < l[i] {
          self.prjctd = true
          x[i] = l[i]
        }
      } else if nbd[i] >= 2 && x[i] >= u[i] {
        if x[i] > u[i] {
          self.prjctd = true
          x[i] = u[i]
        }
      }
    }
    i += 1
  }
  i = 0
  while i < n {
    if nbd[i] != 2 {
      self.boxed = false
    }
    if nbd[i] == 0 {
      self.iwhere[i] = -1
    } else {
      self.cnstnd = true
      if nbd[i] == 2 && u[i] - l[i] <= 0.0 {
        self.iwhere[i] = 3
      } else {
        self.iwhere[i] = 0
      }
    }
    i += 1
  }
}

///|
/// projgr: infinity norm of the projected gradient.
fn _lbfgsb_projgr(
  self : Lbfgsb,
  l : Array[Double],
  u : Array[Double],
  nbd : Array[Int],
  x : Array[Double],
  g : Array[Double],
) -> Unit {
  let n = self.n
  let mut sbg = 0.0
  let mut i = 0
  while i < n {
    let mut gi = g[i]
    if nbd[i] != 0 {
      if gi < 0.0 {
        if nbd[i] >= 2 {
          let d1 = x[i] - u[i]
          if gi < d1 {
            gi = d1
          }
        }
      } else if nbd[i] <= 2 {
        let d1 = x[i] - l[i]
        if gi > d1 {
          gi = d1
        }
      }
    }
    if sbg < gi.abs() {
      sbg = gi.abs()
    }
    i += 1
  }
  self.sbgnrm = sbg
}

///|
/// freev: entering/leaving variables and the free set at the GCP.
fn _lbfgsb_freev(self : Lbfgsb, wrk : (Bool) -> Unit) -> Unit {
  let n = self.n
  let mut nenter = 0
  let mut ileave = n + 1 // 1-based Fortran value
  if self.iter > 0 && self.cnstnd {
    let mut i = 0
    while i < self.nfree {
      let k = self.indx[i]
      if self.iwhere[k] > 0 {
        ileave -= 1
        self.indx2[ileave - 1] = k
      }
      i += 1
    }
    let mut i = self.nfree
    while i < n {
      let k = self.indx[i]
      if self.iwhere[k] <= 0 {
        nenter += 1
        self.indx2[nenter - 1] = k
      }
      i += 1
    }
  }
  wrk(ileave < n + 1 || nenter > 0 || self.updatd)
  let mut nfree = 0
  let mut iact = n
  let mut i = 0
  while i < n {
    if self.iwhere[i] <= 0 {
      self.indx[nfree] = i
      nfree += 1
    } else {
      iact -= 1
      self.indx[iact] = i
    }
    i += 1
  }
  self.nenter = nenter
  self.ileave = ileave
  self.nfree = nfree
}

///|
/// hpsolb: heapsort helper (nleft active breakpoints).
fn _lb_hpsolb(n : Int, t : Array[Double], iorder : Array[Int], iheap : Int) -> Unit {
  if n <= 0 {
    return
  }
  if iheap == 0 {
    let mut k = 1
    while k < n {
      let ddum = t[k]
      let indxin = iorder[k]
      let mut i = k
      let mut looping = true
      while looping {
        if i > 0 {
          let j = (i + 1) / 2 - 1
          if ddum < t[j] {
            t[i] = t[j]
            iorder[i] = iorder[j]
            i = j
          } else {
            looping = false
          }
        } else {
          looping = false
        }
      }
      t[i] = ddum
      iorder[i] = indxin
      k += 1
    }
  }
  if n > 1 {
    let mut i = 0
    let out = t[0]
    let indxou = iorder[0]
    let ddum = t[n - 1]
    let indxin = iorder[n - 1]
    let mut looping = true
    while looping {
      let j = 2 * i + 1
      if j <= n - 2 {
        let mut jm = j
        if t[jm + 1] < t[jm] {
          jm += 1
        }
        if t[jm] < ddum {
          t[i] = t[jm]
          iorder[i] = iorder[jm]
          i = jm
        } else {
          looping = false
        }
      } else {
        looping = false
      }
    }
    t[i] = ddum
    iorder[i] = indxin
    t[n - 1] = out
    iorder[n - 1] = indxou
  }
}

///|
/// bmv: product of the 2m x 2m middle matrix with vector v, returned in p.
fn _lbfgsb_bmv(
  self : Lbfgsb,
  col : Int,
  v : Array[Double],
  v0 : Int,
  p : Array[Double],
  p0 : Int,
) -> Int {
  let m = self.m
  let sy = self.sy
  let wt = self.wt
  if col == 0 {
    return 0
  }
  p[p0 + col] = v[v0 + col]
  let mut i = 1
  while i < col {
    let i2 = col + i
    let mut sum = 0.0
    let mut k = 0
    while k < i {
      sum += sy[i + k * m] * v[v0 + k] / sy[k + k * m]
      k += 1
    }
    p[p0 + i2] = v[v0 + i2] + sum
    i += 1
  }
  let mut info = _lb_dtrsl_11(wt, 0, m, col, p, p0 + col)
  if info != 0 {
    return info
  }
  let mut i = 0
  while i < col {
    p[p0 + i] = v[v0 + i] / sy[i + i * m].sqrt()
    i += 1
  }
  info = _lb_dtrsl_01(wt, 0, m, col, p, p0 + col)
  if info != 0 {
    return info
  }
  let mut i = 0
  while i < col {
    p[p0 + i] = -p[p0 + i] / sy[i + i * m].sqrt()
    i += 1
  }
  i = 0
  while i < col {
    let mut sum = 0.0
    let mut k = i + 1
    while k < col {
      sum += sy[k + i * m] * p[p0 + col + k] / sy[i + i * m]
      k += 1
    }
    p[p0 + i] += sum
    i += 1
  }
  0
}

///|
/// cauchy: generalized Cauchy point. p,c,wbp,v are 2m working vectors into
/// self.wa at the given offsets.
fn _lbfgsb_cauchy(
  self : Lbfgsb,
  x : Array[Double],
  l : Array[Double],
  u : Array[Double],
  nbd : Array[Int],
  g : Array[Double],
  iorder : Array[Int],
  t : Array[Double],
  xcp : Array[Double],
  p0 : Int,
  c0 : Int,
  wbp0 : Int,
  v0 : Int,
) -> Int {
  let n = self.n
  let m = self.m
  let col = self.col
  let theta = self.theta
  let head = self.head
  let wa = self.wa
  let wy = self.wy
  let ws = self.ws
  if self.sbgnrm <= 0.0 {
    _lb_dcopy(n, x, 0, xcp, 0)
    return 0
  }
  let mut bnded = true
  let mut nfree = n + 1
  let mut nbreak = 0
  let mut ibkmin = 0 // 1-based position
  let mut bkmin = 0.0
  let col2 = col * 2
  let mut f1 = 0.0
  let mut j = 0
  while j < col2 {
    wa[p0 + j] = 0.0
    j += 1
  }
  let mut i = 0
  while i < n {
    let neggi = -g[i]
    let mut tl = 0.0
    let mut tu = 0.0
    if self.iwhere[i] != 3 && self.iwhere[i] != -1 {
      if nbd[i] <= 2 {
        tl = x[i] - l[i]
      }
      if nbd[i] >= 2 {
        tu = u[i] - x[i]
      }
      let xlower = nbd[i] <= 2 && tl <= 0.0
      let xupper = nbd[i] >= 2 && tu <= 0.0
      self.iwhere[i] = 0
      if xlower {
        if neggi <= 0.0 {
          self.iwhere[i] = 1
        }
      } else if xupper {
        if neggi >= 0.0 {
          self.iwhere[i] = 2
        }
      } else if neggi.abs() <= 0.0 {
        self.iwhere[i] = -3
      }
    }
    let mut pointr = head
    if self.iwhere[i] != 0 && self.iwhere[i] != -1 {
      self.d[i] = 0.0
    } else {
      self.d[i] = neggi
      f1 -= neggi * neggi
      let mut jj = 0
      while jj < col {
        wa[p0 + jj] += wy[i + (pointr - 1) * n] * neggi
        wa[p0 + col + jj] += ws[i + (pointr - 1) * n] * neggi
        pointr = pointr % m + 1
        jj += 1
      }
      if nbd[i] <= 2 && nbd[i] != 0 && neggi < 0.0 {
        iorder[nbreak] = i
        t[nbreak] = tl / (-neggi)
        if nbreak == 0 || t[nbreak] < bkmin {
          bkmin = t[nbreak]
          ibkmin = nbreak + 1
        }
        nbreak += 1
      } else if nbd[i] >= 2 && neggi > 0.0 {
        iorder[nbreak] = i
        t[nbreak] = tu / neggi
        if nbreak == 0 || t[nbreak] < bkmin {
          bkmin = t[nbreak]
          ibkmin = nbreak + 1
        }
        nbreak += 1
      } else {
        nfree -= 1
        iorder[nfree] = i
        if neggi.abs() > 0.0 {
          bnded = false
        }
      }
    }
    i += 1
  }
  if theta != 1.0 {
    _lb_dscal(col, theta, wa, p0 + col)
  }
  _lb_dcopy(n, x, 0, xcp, 0)
  if nbreak == 0 && nfree == n + 1 {
    return 0
  }
  j = 0
  while j < col2 {
    wa[c0 + j] = 0.0
    j += 1
  }
  let mut f2 = -theta * f1
  let f2_org = f2
  if col > 0 {
    let info = self._lbfgsb_bmv(col, wa, p0, wa, v0)
    if info != 0 {
      return info
    }
    f2 -= _lb_ddot(col2, wa, v0, wa, p0)
  }
  let mut dtm = -f1 / f2
  let mut tsum = 0.0
  self.nint = 1
  if nbreak == 0 {
    // L888
    if dtm <= 0.0 {
      dtm = 0.0
    }
    tsum += dtm
    _lb_daxpy(n, tsum, self.d, 0, xcp, 0)
    if col > 0 {
      _lb_daxpy(col2, dtm, wa, p0, wa, c0)
    }
    return 0
  }
  let mut nleft = nbreak
  let mut cit = 1
  let mut tj = 0.0
  let mut searching = true
  while searching {
    let tj0 = tj
    let mut ibp : Int = 0
    if cit == 1 {
      tj = bkmin
      ibp = iorder[ibkmin - 1]
    } else {
      if cit == 2 {
        if ibkmin != nbreak {
          t[ibkmin - 1] = t[nbreak - 1]
          iorder[ibkmin - 1] = iorder[nbreak - 1]
        }
      }
      _lb_hpsolb(nleft, t, iorder, cit - 2)
      tj = t[nleft - 1]
      ibp = iorder[nleft - 1]
    }
    let dt = tj - tj0
    if dtm < dt {
      searching = false
    } else {
      tsum += dt
      nleft -= 1
      cit += 1
      let dibp = self.d[ibp]
      self.d[ibp] = 0.0
      let zibp = if dibp > 0.0 {
        xcp[ibp] = u[ibp]
        self.iwhere[ibp] = 2
        u[ibp] - x[ibp]
      } else {
        xcp[ibp] = l[ibp]
        self.iwhere[ibp] = 1
        l[ibp] - x[ibp]
      }
      if nleft == 0 && nbreak == n {
        dtm = dt
        tsum += dtm
        if col > 0 {
          _lb_daxpy(col2, dtm, wa, p0, wa, c0)
        }
        return 0
      }
      self.nint += 1
      let dibp2 = dibp * dibp
      f1 += dt * f2 + dibp2 - theta * dibp * zibp
      f2 -= theta * dibp2
      if col > 0 {
        _lb_daxpy(col2, dt, wa, p0, wa, c0)
        let mut pointr = head
        let mut jj = 0
        while jj < col {
          wa[wbp0 + jj] = wy[ibp + (pointr - 1) * n]
          wa[wbp0 + col + jj] = theta * ws[ibp + (pointr - 1) * n]
          pointr = pointr % m + 1
          jj += 1
        }
        let info = self._lbfgsb_bmv(col, wa, wbp0, wa, v0)
        if info != 0 {
          return info
        }
        let wmc = _lb_ddot(col2, wa, c0, wa, v0)
        let wmp = _lb_ddot(col2, wa, p0, wa, v0)
        let wmw = _lb_ddot(col2, wa, wbp0, wa, v0)
        _lb_daxpy(col2, -dibp, wa, wbp0, wa, p0)
        f1 += dibp * wmc
        f2 += 2.0 * dibp * wmp - dibp2 * wmw
      }
      let floor = self.epsmch * f2_org
      if f2 < floor {
        f2 = floor
      }
      if nleft > 0 {
        dtm = -f1 / f2
      } else if bnded {
        f1 = 0.0
        f2 = 0.0
        dtm = 0.0
        searching = false
      } else {
        dtm = -f1 / f2
        searching = false
      }
    }
  }
  // L888: minimizer inside the final segment.
  if dtm <= 0.0 {
    dtm = 0.0
  }
  tsum += dtm
  _lb_daxpy(n, tsum, self.d, 0, xcp, 0)
  if col > 0 {
    _lb_daxpy(col2, dtm, wa, p0, wa, c0)
  }
  0
}

///|
/// formt: T = theta*SS + L D^-1 L', Cholesky factorized in place.
fn _lbfgsb_formt(self : Lbfgsb, col : Int) -> Int {
  let m = self.m
  let wt = self.wt
  let sy = self.sy
  let ss = self.ss
  let theta = self.theta
  let mut j = 0
  while j < col {
    wt[j * m] = theta * ss[j * m]
    j += 1
  }
  let mut i = 1
  while i < col {
    let mut j = i
    while j < col {
      let mut ddum = 0.0
      let mut k = 0
      while k < i {
        ddum += sy[i + k * m] * sy[j + k * m] / sy[k + k * m]
        k += 1
      }
      wt[i + j * m] = ddum + theta * ss[i + j * m]
      j += 1
    }
    i += 1
  }
  if _lb_dpofa(wt, 0, m, col) != 0 {
    return -3
  }
  0
}

///|
/// matupd: update WS/WY and form SS, SY (lower triangle).
fn _lbfgsb_matupd(
  self : Lbfgsb,
  rr : Double,
  dr : Double,
) -> Unit {
  let n = self.n
  let m = self.m
  let mut col = self.col
  let mut itail = self.itail
  let mut head = self.head
  if self.iupdat <= m {
    col = self.iupdat
    itail = (head + self.iupdat - 2) % m + 1
  } else {
    itail = itail % m + 1
    head = head % m + 1
  }
  self.col = col
  self.itail = itail
  self.head = head
  let c0 = (itail - 1) * n
  _lb_dcopy(n, self.d, 0, self.ws, c0)
  _lb_dcopy(n, self.r, 0, self.wy, c0)
  self.theta = rr / dr
  if self.iupdat > m {
    let mut j = 0
    while j < col - 1 {
      // shift column j+1 of SS into column j (j entries)
      let mut k = 0
      while k <= j {
        self.ss[k + j * m] = self.ss[k + 1 + (j + 1) * m]
        k += 1
      }
      // shift column j+1 of SY into column j (col-j-1 entries)
      k = 0
      while k < col - j - 1 {
        self.sy[j + k + j * m] = self.sy[j + 1 + k + (j + 1) * m]
        k += 1
      }
      j += 1
    }
  }
  let mut pointr = head
  let mut j = 0
  while j < col - 1 {
    self.sy[col - 1 + j * m] = _lb_ddot(n, self.d, 0, self.wy, (pointr - 1) * n)
    self.ss[j + (col - 1) * m] = _lb_ddot(n, self.ws, (pointr - 1) * n, self.d, 0)
    pointr = pointr % m + 1
    j += 1
  }
  if self.stp == 1.0 {
    self.ss[col - 1 + (col - 1) * m] = self.dtd
  } else {
    self.ss[col - 1 + (col - 1) * m] = self.stp * self.stp * self.dtd
  }
  self.sy[col - 1 + (col - 1) * m] = dr
}

///|
/// formk: LEL^T factorization of the indefinite 2col x 2col matrix.
fn _lbfgsb_formk(self : Lbfgsb) -> Int {
  let n = self.n
  let m = self.m
  let m2 = 2 * m
  let col = self.col
  let wn = self.wn
  let wn1 = self.snd
  let ws = self.ws
  let wy = self.wy
  let sy = self.sy
  let theta = self.theta
  let nsub = self.nfree
  let indx = self.indx
  let indx2 = self.indx2
  let nenter = self.nenter
  let ileave0 = self.ileave - 1 // 0-based start of leaving indices
  let mut upcl = col
  if self.updatd {
    if self.iupdat > m {
      // jy is 1-based (1..m-1), kb = jy-1 is 0-based.
      let mut jy = 1
      while jy <= m - 1 {
        let len = m - jy
        let kb = jy - 1
        // wn1 block (1,1): (jy,jy) <- (jy+1,jy+1), 1-based
        let mut k = 0
        while k < len {
          wn1[kb + k + kb * m2] = wn1[kb + 1 + k + (kb + 1) * m2]
          k += 1
        }
        // wn1 block (2,2): (m+jy,m+jy) <- (m+jy+1,m+jy+1), 1-based
        let dr = m + kb // 0-based destination diagonal row
        k = 0
        while k < len {
          wn1[dr + k + dr * m2] = wn1[dr + 1 + k + (dr + 1) * m2]
          k += 1
        }
        // wn1 block (2,1): (m+1..2m-1, jy) <- (m+2..2m, jy+1), 1-based
        k = 0
        while k < m - 1 {
          wn1[m + k + kb * m2] = wn1[m + 1 + k + (kb + 1) * m2]
          k += 1
        }
        jy += 1
      }
    }
    let pbegin = 0
    let pend = nsub - 1
    let dbegin = nsub
    let dend = n - 1
    let iy = col - 1
    let isv = m + col - 1
    let mut ipntr = self.head + col - 1
    if ipntr > m {
      ipntr -= m
    }
    let mut jpntr = self.head
    let mut jy = 1
    while jy <= col {
      let js = m + jy - 1
      let mut temp1 = 0.0
      let mut temp2 = 0.0
      let mut temp3 = 0.0
      let mut k = pbegin
      while k <= pend {
        let k1 = indx[k]
        temp1 += wy[k1 + (ipntr - 1) * n] * wy[k1 + (jpntr - 1) * n]
        k += 1
      }
      k = dbegin
      while k <= dend {
        let k1 = indx[k]
        temp2 += ws[k1 + (ipntr - 1) * n] * ws[k1 + (jpntr - 1) * n]
        temp3 += ws[k1 + (ipntr - 1) * n] * wy[k1 + (jpntr - 1) * n]
        k += 1
      }
      wn1[iy + (jy - 1) * m2] = temp1
      wn1[isv + js * m2] = temp2
      wn1[isv + (jy - 1) * m2] = temp3
      jpntr = jpntr % m + 1
      jy += 1
    }
    // new column in block (2,1)
    let jy = col - 1
    let mut jpntrb = self.head + col - 1
    if jpntrb > m {
      jpntrb -= m
    }
    let mut ipntrb = self.head
    let mut ii = 0
    while ii < col {
      let isb = m + ii
      let mut temp3 = 0.0
      let mut k = pbegin
      while k <= pend {
        let k1 = indx[k]
        temp3 += ws[k1 + (ipntrb - 1) * n] * wy[k1 + (jpntrb - 1) * n]
        k += 1
      }
      ipntrb = ipntrb % m + 1
      wn1[isb + jy * m2] = temp3
      ii += 1
    }
    upcl = col - 1
  }
  // modify old blocks for changed free variables
  let mut ipntr = self.head
  let mut iy = 0
  while iy < upcl {
    let isv = m + iy
    let mut jpntr = self.head
    let mut jy = 0
    while jy <= iy {
      let js = m + jy
      let mut temp1 = 0.0
      let mut temp2 = 0.0
      let mut temp3 = 0.0
      let mut temp4 = 0.0
      let mut k = 0
      while k < nenter {
        let k1 = indx2[k]
        temp1 += wy[k1 + (ipntr - 1) * n] * wy[k1 + (jpntr - 1) * n]
        temp2 += ws[k1 + (ipntr - 1) * n] * ws[k1 + (jpntr - 1) * n]
        k += 1
      }
      k = ileave0
      while k < n {
        let k1 = indx2[k]
        temp3 += wy[k1 + (ipntr - 1) * n] * wy[k1 + (jpntr - 1) * n]
        temp4 += ws[k1 + (ipntr - 1) * n] * ws[k1 + (jpntr - 1) * n]
        k += 1
      }
      wn1[iy + jy * m2] = wn1[iy + jy * m2] + temp1 - temp3
      wn1[isv + js * m2] = wn1[isv + js * m2] - temp2 + temp4
      jpntr = jpntr % m + 1
      jy += 1
    }
    ipntr = ipntr % m + 1
    iy += 1
  }
  ipntr = self.head
  let mut is1 = m
  while is1 < m + upcl {
    let mut jpntr = self.head
    let mut jy = 0
    while jy < upcl {
      let mut temp1 = 0.0
      let mut temp3 = 0.0
      let mut k = 0
      while k < nenter {
        let k1 = indx2[k]
        temp1 += ws[k1 + (ipntr - 1) * n] * wy[k1 + (jpntr - 1) * n]
        k += 1
      }
      k = ileave0
      while k < n {
        let k1 = indx2[k]
        temp3 += ws[k1 + (ipntr - 1) * n] * wy[k1 + (jpntr - 1) * n]
        k += 1
      }
      if is1 + 1 <= jy + 1 + m {
        wn1[is1 + jy * m2] += temp1 - temp3
      } else {
        wn1[is1 + jy * m2] += -temp1 + temp3
      }
      jpntr = jpntr % m + 1
      jy += 1
    }
    ipntr = ipntr % m + 1
    is1 += 1
  }
  // form the upper triangle of wn
  let mut iiy = 0
  while iiy < col {
    let isv = col + iiy
    let is1b = m + iiy
    let mut jy = 0
    while jy <= iiy {
      let js = col + jy
      let js1b = m + jy
      wn[jy + iiy * m2] = wn1[iiy + jy * m2] / theta
      wn[js + isv * m2] = wn1[is1b + js1b * m2] * theta
      jy += 1
    }
    jy = 0
    while jy < iiy {
      wn[jy + isv * m2] = -wn1[is1b + jy * m2]
      jy += 1
    }
    jy = iiy
    while jy < col {
      wn[jy + isv * m2] = wn1[is1b + jy * m2]
      jy += 1
    }
    wn[iiy + iiy * m2] += sy[iiy + iiy * m]
    iiy += 1
  }
  // first Cholesky factor (1,1) block
  if _lb_dpofa(wn, 0, m2, col) != 0 {
    return -1
  }
  let mut js0 = col
  while js0 < 2 * col {
    if _lb_dtrsl_11(wn, 0, m2, col, wn, js0 * m2) != 0 {
      return -1
    }
    js0 += 1
  }
  // update (2,2) block
  let mut is0 = col
  while is0 < 2 * col {
    let mut js = is0
    while js < 2 * col {
      wn[is0 + js * m2] += _lb_ddot(col, wn, is0 * m2, wn, js * m2)
      js += 1
    }
    is0 += 1
  }
  // second Cholesky factor
  if _lb_dpofa(wn, col + col * m2, m2, col) != 0 {
    return -2
  }
  0
}

///|
/// cmprlb: r = -Z'B(xcp-xk) - Z'g.
fn _lbfgsb_cmprlb(
  self : Lbfgsb,
  x : Array[Double],
  g : Array[Double],
) -> Int {
  let n = self.n
  let m = self.m
  let col = self.col
  let nfree = self.nfree
  let indx = self.indx
  let wa = self.wa
  if !self.cnstnd && col > 0 {
    let mut i = 0
    while i < n {
      self.r[i] = -g[i]
      i += 1
    }
  } else {
    let mut i = 0
    while i < nfree {
      let k = indx[i]
      self.r[i] = -self.theta * (self.z[k] - x[k]) - g[k]
      i += 1
    }
    let info = self._lbfgsb_bmv(col, wa, 2 * m, wa, 0)
    if info != 0 {
      return -8
    }
    let mut pointr = self.head
    let mut j = 0
    while j < col {
      let a1 = wa[j]
      let a2 = self.theta * wa[col + j]
      let mut i = 0
      while i < nfree {
        let k = indx[i]
        self.r[i] += self.wy[k + (pointr - 1) * n] * a1 +
          self.ws[k + (pointr - 1) * n] * a2
        i += 1
      }
      pointr = pointr % m + 1
      j += 1
    }
  }
  0
}

///|
/// subsm: direct subspace minimization using the factorized matrix K.
/// On entry x is the Cauchy point and d (self.r) the reduced gradient;
/// x is updated to the subspace minimizer and d to the Newton direction.
fn _lbfgsb_subsm(
  self : Lbfgsb,
  l : Array[Double],
  u : Array[Double],
  nbd : Array[Int],
) -> Int {
  let n = self.n
  let m = self.m
  let ns = self.nfree
  let col = self.col
  if ns <= 0 {
    return 0
  }
  let indx = self.indx
  let wa = self.wa
  let wn = self.wn
  let m2 = 2 * m
  let col2 = 2 * col
  // wv = W'Zd
  let mut pointr = self.head
  let mut i = 0
  while i < col {
    let mut temp1 = 0.0
    let mut temp2 = 0.0
    let mut j = 0
    while j < ns {
      let k = indx[j]
      temp1 += self.wy[k + (pointr - 1) * n] * self.r[j]
      temp2 += self.ws[k + (pointr - 1) * n] * self.r[j]
      j += 1
    }
    wa[i] = temp1
    wa[col + i] = self.theta * temp2
    pointr = pointr % m + 1
    i += 1
  }
  // wv := K^-1 wv
  let mut info = _lb_dtrsl_11(wn, 0, m2, col2, wa, 0)
  if info != 0 {
    return info
  }
  i = 0
  while i < col {
    wa[i] = -wa[i]
    i += 1
  }
  info = _lb_dtrsl_01(wn, 0, m2, col2, wa, 0)
  if info != 0 {
    return info
  }
  // d = d/theta + Z'W wv / theta
  pointr = self.head
  let mut jy = 0
  while jy < col {
    let js = col + jy
    let mut ii = 0
    while ii < ns {
      let k = indx[ii]
      self.r[ii] += self.wy[k + (pointr - 1) * n] * wa[jy] / self.theta +
        self.ws[k + (pointr - 1) * n] * wa[js]
      ii += 1
    }
    pointr = pointr % m + 1
    jy += 1
  }
  i = 0
  while i < ns {
    self.r[i] /= self.theta
    i += 1
  }
  // backtrack to the feasible region
  let mut alpha = 1.0
  let mut temp1 = 1.0
  let mut ibd = 0
  i = 0
  while i < ns {
    let k = indx[i]
    let dk = self.r[i]
    if nbd[k] != 0 {
      if dk < 0.0 && nbd[k] <= 2 {
        let temp2 = l[k] - self.z[k]
        if temp2 >= 0.0 {
          temp1 = 0.0
        } else if dk * alpha < temp2 {
          temp1 = temp2 / dk
        }
      } else if dk > 0.0 && nbd[k] >= 2 {
        let temp2 = u[k] - self.z[k]
        if temp2 <= 0.0 {
          temp1 = 0.0
        } else if dk * alpha > temp2 {
          temp1 = temp2 / dk
        }
      }
      if temp1 < alpha {
        alpha = temp1
        ibd = i
      }
    }
    i += 1
  }
  if alpha < 1.0 {
    let dk = self.r[ibd]
    let k = indx[ibd]
    if dk > 0.0 {
      self.z[k] = u[k]
      self.r[ibd] = 0.0
    } else if dk < 0.0 {
      self.z[k] = l[k]
      self.r[ibd] = 0.0
    }
  }
  i = 0
  while i < ns {
    self.z[indx[i]] += alpha * self.r[i]
    i += 1
  }
  self.iword = if alpha < 1.0 { 1 } else { 0 }
  0
}

///|
/// R's bounded central-difference gradient (optim_stats.c fmingr).
fn _lbfgsb_numgrad(
  self : Lbfgsb,
  x : Array[Double],
  l : Array[Double],
  u : Array[Double],
  obj : (Array[Double]) -> Double,
) -> Array[Double] {
  let n = self.n
  let g = Array::make(n, 0.0)
  let mut i = 0
  while i < n {
    let xi = x[i]
    let eps0 = self.ndeps
    let mut epsused = eps0
    let mut tmp = xi + eps0
    if tmp > u[i] {
      tmp = u[i]
      epsused = tmp - xi
    }
    x[i] = tmp
    let val1 = obj(x)
    let mut eps = eps0
    tmp = xi - eps0
    if tmp < l[i] {
      tmp = l[i]
      eps = xi - tmp
    }
    x[i] = tmp
    let val2 = obj(x)
    x[i] = xi
    g[i] = (val1 - val2) / (epsused + eps)
    i += 1
  }
  g
}

///|
/// Reset the L-BFGS memory after a detected factorization failure.
fn _lbfgsb_refresh(self : Lbfgsb) -> Unit {
  self.info = 0
  self.col = 0
  self.head = 1
  self.theta = 1.0
  self.iupdat = 0
  self.updatd = false
}

///|
/// Run setulb/mainlb plus the optim.c driver. Returns (f, fail):
/// fail 0 = converged, 1 = maxit exceeded, 51 = warning, 52 = error.
fn _lbfgsb_run(
  self : Lbfgsb,
  x : Array[Double],
  l : Array[Double],
  u : Array[Double],
  nbd : Array[Int],
  obj : (Array[Double]) -> Double,
  maxit : Int,
) -> (Double, Int) {
  let n = self.n
  // START: initialize scalars, project x, classify variables.
  self.col = 0
  self.head = 1
  self.theta = 1.0
  self.iupdat = 0
  self.updatd = false
  self.iback = 0
  self.itail = 0
  self.ifun = 0
  self.iword = 0
  self.nact = 0
  self.ileave = 0
  self.nenter = 0
  self.iter = 0
  self.nfgv = 0
  self.nint = 0
  self.nintol = 0
  self.nskip = 0
  self.nfree = n
  self.tol = self.factr * self.epsmch
  self.info = 0
  self._lbfgsb_active(l, u, nbd, x)
  // FG_START: evaluate f0 and g0.
  let mut f = obj(x)
  let mut g = self._lbfgsb_numgrad(x, l, u, obj)
  self.nfgv = 1
  // L111: projected gradient test.
  self._lbfgsb_projgr(l, u, nbd, x, g)
  if self.sbgnrm <= self.pgtol {
    return (f, 0)
  }
  let mut driver_iter = 0
  let mut running = true
  let mut fail = 0
  while running {
    self.iword = -1
    let mut wrk = false
    if !self.cnstnd && self.col > 0 {
      _lb_dcopy(n, x, 0, self.z, 0)
      wrk = self.updatd
      self.nint = 0
    } else {
      // GCP
      self.info = self._lbfgsb_cauchy(
        x,
        l,
        u,
        nbd,
        g,
        self.indx2,
        self.t,
        self.z,
        0,
        2 * self.m,
        4 * self.m,
        6 * self.m,
      )
      if self.info != 0 {
        self._lbfgsb_refresh()
        continue
      }
      self.nintol += self.nint
      self._lbfgsb_freev(fn(v : Bool) {
        wrk = v
      })
      self.nact = n - self.nfree
    }
    // L333: subspace minimization.
    if self.nfree != 0 && self.col != 0 {
      if wrk {
        self.info = self._lbfgsb_formk()
      }
      if self.info != 0 {
        self._lbfgsb_refresh()
        continue
      }
      self.info = self._lbfgsb_cmprlb(x, g)
      if self.info == 0 {
        self.info = self._lbfgsb_subsm(l, u, nbd)
      }
      if self.info != 0 {
        self._lbfgsb_refresh()
        continue
      }
    }
    // L555: search direction d = z - x.
    let mut i = 0
    while i < n {
      self.d[i] = self.z[i] - x[i]
      i += 1
    }
    // lnsrlb setup.
    self.dtd = _lb_ddot(n, self.d, 0, self.d, 0)
    self.dnorm = self.dtd.sqrt()
    self.stpmx = 1.0e10
    if self.cnstnd {
      if self.iter == 0 {
        self.stpmx = 1.0
      } else {
        let mut ii = 0
        while ii < n {
          let a1 = self.d[ii]
          if nbd[ii] != 0 {
            if a1 < 0.0 && nbd[ii] <= 2 {
              let a2 = l[ii] - x[ii]
              if a2 >= 0.0 {
                self.stpmx = 0.0
              } else if a1 * self.stpmx < a2 {
                self.stpmx = a2 / a1
              }
            } else if a1 > 0.0 && nbd[ii] >= 2 {
              let a2 = u[ii] - x[ii]
              if a2 <= 0.0 {
                self.stpmx = 0.0
              } else if a1 * self.stpmx > a2 {
                self.stpmx = a2 / a1
              }
            }
          }
          ii += 1
        }
      }
    }
    self.stp = if self.iter == 0 && !self.boxed {
      (1.0 / self.dnorm).min(self.stpmx)
    } else {
      1.0
    }
    _lb_dcopy(n, x, 0, self.t, 0)
    _lb_dcopy(n, g, 0, self.r, 0)
    self.fold = f
    self.ifun = 0
    self.iback = 0
    let dc = _lb_dcsrch_start(self.stpmx)
    let mut first_call = true
    let mut ls_running = true
    // L556: line search with inline FG evaluations.
    while ls_running {
      self.gd = _lb_ddot(n, g, 0, self.d, 0)
      if self.ifun == 0 {
        self.gdold = self.gd
        if self.gd >= 0.0 {
          self.info = -4
          ls_running = false
        }
      }
      if ls_running {
        let new_stp = dc._lb_dcsrch_apply(first_call, f, self.gd, self.stp)
        first_call = false
        self.stp = new_stp
        if dc.status == 1 {
          // FG_LNSRCH: move x and evaluate.
          self.ifun += 1
          self.nfgv += 1
          self.iback = self.ifun - 1
          if self.stp == 1.0 {
            _lb_dcopy(n, self.z, 0, x, 0)
          } else {
            let mut ii = 0
            while ii < n {
              x[ii] = self.stp * self.d[ii] + self.t[ii]
              ii += 1
            }
          }
          f = obj(x)
          g = self._lbfgsb_numgrad(x, l, u, obj)
        } else {
          // CONVERGENCE or WARNING from the line search.
          ls_running = false
        }
      }
    }
    if self.info != 0 || self.iback >= 20 {
      // restore the previous iterate
      _lb_dcopy(n, self.t, 0, x, 0)
      _lb_dcopy(n, self.r, 0, g, 0)
      f = self.fold
      if self.col == 0 {
        if self.info == 0 {
          self.info = -9
          self.nfgv -= 1
          self.ifun -= 1
          self.iback -= 1
        }
        self.iter += 1
        fail = 52
        running = false
      } else {
        if self.info == 0 {
          self.nfgv -= 1
        }
        self._lbfgsb_refresh()
        // RESTART_FROM_LNSRCH: continue the outer loop.
      }
    } else {
      // NEW_X
      self.iter += 1
      self._lbfgsb_projgr(l, u, nbd, x, g)
      driver_iter += 1
      if driver_iter > maxit {
        fail = 1
        running = false
      }
    }
    if running {
      // L777: termination tests and L-BFGS update.
      if self.sbgnrm <= self.pgtol {
        fail = 0
        running = false
      } else {
        let ddum = self.fold.abs().max(f.abs()).max(1.0)
        if self.fold - f <= self.tol * ddum {
          if self.iback >= 10 {
            self.info = -5
          }
          fail = 0
          running = false
        }
      }
    }
    if running {
      let mut i = 0
      while i < n {
        self.r[i] = g[i] - self.r[i]
        i += 1
      }
      let rr = _lb_ddot(n, self.r, 0, self.r, 0)
      let (dr, ddum) = if self.stp == 1.0 {
        (self.gd - self.gdold, -self.gdold)
      } else {
        let drv = (self.gd - self.gdold) * self.stp
        _lb_dscal(n, self.stp, self.d, 0)
        (drv, -self.gdold * self.stp)
      }
      if dr <= self.epsmch * ddum {
        // skip the L-BFGS update
        self.nskip += 1
        self.updatd = false
      } else {
        self.updatd = true
        self.iupdat += 1
        self._lbfgsb_matupd(rr, dr)
        self.info = self._lbfgsb_formt(self.col)
        if self.info != 0 {
          self._lbfgsb_refresh()
        }
      }
    }
  }
  (f, fail)
}

///|
/// Bounded L-BFGS-B minimization using R `optim` defaults (lmm = 5,
/// factr = 1e7, pgtol = 0, maxit = 100, bounded central differences with
/// ndeps = 1e-3). Returns (x, f, number of function/gradient pairs).
fn _sw_lbfgsb_minimize(
  x0 : Array[Double],
  l : Array[Double],
  u : Array[Double],
  nbd : Array[Int],
  obj : (Array[Double]) -> Double,
) -> (Array[Double], Double, Int) {
  let n = x0.length()
  let x = x0.copy()
  let solver = _lbfgsb_new(n, 5, 1.0e7, 0.0, 1.0e-3)
  let (f, _fail) = solver._lbfgsb_run(x, l, u, nbd, obj, 100)
  (x, f, solver.nfgv)
}
