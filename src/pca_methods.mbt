///|
/// pcaMethods-style probabilistic / Bayesian PCA with missing data.
///
/// Port of Bioconductor's pcaMethods PPCA (Tipping & Bishop 1999) and BPCA
/// (Oba et al. 2003). Both estimate principal axes and scores from a data
/// matrix that may contain NaN (missing) entries; BPCA additionally provides
/// a Bayesian treatment with automatic relevance determination. The reference
/// equations are transcribed from pcaMethods' R implementation (`ppca.R`,
/// `bpca.R` + `BPCA_initmodel.R` + `BPCA_dostep.R`).

///|
pub struct PcaMethodsResult {
  scores : Array[Array[Double]]
  loadings : Array[Array[Double]]
  r2cum : Array[Double]
  tau : Double
  sigma2 : Double
  imputed : Array[Array[Double]]
}

// ---------------------------------------------------------------------------
// Small dense-matrix helpers (rows × columns)
// ---------------------------------------------------------------------------

///|
fn pm_nan() -> Double {
  0.0 / 0.0
}

///|
fn pm_mul(a : Array[Array[Double]], b : Array[Array[Double]]) -> Array[Array[Double]] {
  let m = a.length()
  let n = if m > 0 { a[0].length() } else { 0 }
  let p = if b.length() > 0 { b[0].length() } else { 0 }
  let c : Array[Array[Double]] = Array(capacity=0)
  for i = 0; i < m; i = i + 1 {
    let row = Array::make(p, 0.0)
    for k = 0; k < n; k = k + 1 {
      let aik = a[i][k]
      for j = 0; j < p; j = j + 1 {
        row[j] = row[j] + aik * b[k][j]
      }
    }
    c.push(row)
  }
  c
}

///|
fn pm_transpose(a : Array[Array[Double]]) -> Array[Array[Double]] {
  let m = a.length()
  let n = if m > 0 { a[0].length() } else { 0 }
  let t : Array[Array[Double]] = Array(capacity=0)
  for j = 0; j < n; j = j + 1 {
    let row = Array::make(m, 0.0)
    for i = 0; i < m; i = i + 1 {
      row[i] = a[i][j]
    }
    t.push(row)
  }
  t
}

///|
fn pm_identity(n : Int) -> Array[Array[Double]] {
  let m : Array[Array[Double]] = Array(capacity=0)
  for i = 0; i < n; i = i + 1 {
    let row = Array::make(n, 0.0)
    row[i] = 1.0
    m.push(row)
  }
  m
}

///|
fn pm_inverse(a : Array[Array[Double]]) -> Array[Array[Double]] {
  linalg_inverse(a)
}

///|
fn pm_trace(a : Array[Array[Double]]) -> Double {
  let n = a.length()
  let mut s = 0.0
  for i = 0; i < n; i = i + 1 {
    s = s + a[i][i]
  }
  s
}

///|
/// Sum of squared entries.
fn pm_sum_sq(a : Array[Array[Double]]) -> Double {
  let mut s = 0.0
  for i = 0; i < a.length(); i = i + 1 {
    for j = 0; j < a[i].length(); j = j + 1 {
      s = s + a[i][j] * a[i][j]
    }
  }
  s
}

///|
/// Sum of the elementwise product of two same-shaped matrices.
fn pm_sum_ew(a : Array[Array[Double]], b : Array[Array[Double]]) -> Double {
  let mut s = 0.0
  for i = 0; i < a.length(); i = i + 1 {
    for j = 0; j < a[i].length(); j = j + 1 {
      s = s + a[i][j] * b[i][j]
    }
  }
  s
}

///|
/// `diag(a^T a)`: squared 2-norms of the columns of `a`.
fn pm_col_norms_sq(a : Array[Array[Double]]) -> Array[Double] {
  let k = if a.length() > 0 { a[0].length() } else { 0 }
  let out = Array::make(k, 0.0)
  for j = 0; j < k; j = j + 1 {
    let mut s = 0.0
    for i = 0; i < a.length(); i = i + 1 {
      s = s + a[i][j] * a[i][j]
    }
    out[j] = s
  }
  out
}

///|
/// Unbiased covariance matrix (D×D) of an N×D matrix (R's `cov`).
fn pm_covariance(x : Array[Array[Double]]) -> Array[Array[Double]] {
  let n = x.length()
  let d = if n > 0 { x[0].length() } else { 0 }
  let means = Array::make(d, 0.0)
  for j = 0; j < d; j = j + 1 {
    let mut s = 0.0
    for i = 0; i < n; i = i + 1 {
      s = s + x[i][j]
    }
    means[j] = s / n.to_double()
  }
  let cov : Array[Array[Double]] = Array(capacity=0)
  for j = 0; j < d; j = j + 1 {
    cov.push(Array::make(d, 0.0))
  }
  for i = 0; i < n; i = i + 1 {
    for j = 0; j < d; j = j + 1 {
      let cj = x[i][j] - means[j]
      for l = 0; l < d; l = l + 1 {
        cov[j][l] = cov[j][l] + cj * (x[i][l] - means[l])
      }
    }
  }
  for j = 0; j < d; j = j + 1 {
    for l = 0; l < d; l = l + 1 {
      cov[j][l] = cov[j][l] / (n - 1).to_double()
    }
  }
  cov
}

///|
/// Orthonormal basis for the range of `a`, via its left singular vectors
/// (matching pcaMethods' `orth`, which uses SVD).
fn pm_orth(a : Array[Array[Double]]) -> Array[Array[Double]] {
  let d = a.length()
  let k = if d > 0 { a[0].length() } else { 0 }
  let ata = pm_mul(pm_transpose(a), a)
  let (evals, v) = pm_top_eigen(ata, k)
  let q : Array[Array[Double]] = Array(capacity=0)
  for i = 0; i < d; i = i + 1 {
    q.push(Array::make(k, 0.0))
  }
  for j = 0; j < k; j = j + 1 {
    let av = Array::make(d, 0.0)
    for i = 0; i < d; i = i + 1 {
      let mut sum = 0.0
      for l = 0; l < k; l = l + 1 {
        sum = sum + a[i][l] * v[l][j]
      }
      av[i] = sum
    }
    let sigma = evals[j].sqrt()
    for i = 0; i < d; i = i + 1 {
      q[i][j] = if sigma > 1.0e-15 { av[i] / sigma } else { 0.0 }
    }
  }
  q
}

///|
/// Correct Jacobi eigenvalue decomposition for a symmetric matrix. Columns of
/// the returned matrix are the eigenvectors.
fn pm_eigh_sym(a : Array[Array[Double]]) -> (Array[Double], Array[Array[Double]]) {
  linalg_jacobi(a, a.length(), 100 * a.length() * a.length())
}

///|
/// Top-`k` eigenvalues (descending) and eigenvectors (columns) of a symmetric
/// matrix.
fn pm_top_eigen(
  a : Array[Array[Double]],
  k : Int,
) -> (Array[Double], Array[Array[Double]]) {
  let n = a.length()
  let (evals, evecs) = pm_eigh_sym(a)
  let indices = Array::make(n, 0)
  for i = 0; i < n; i = i + 1 {
    indices[i] = i
  }
  for i = 0; i < n - 1; i = i + 1 {
    for j = i + 1; j < n; j = j + 1 {
      if evals[indices[j]] > evals[indices[i]] {
        let t = indices[i]
        indices[i] = indices[j]
        indices[j] = t
      }
    }
  }
  let top_evals = Array::make(k, 0.0)
  let top_evecs : Array[Array[Double]] = Array(capacity=0)
  for i = 0; i < n; i = i + 1 {
    let row = Array::make(k, 0.0)
    for j = 0; j < k; j = j + 1 {
      row[j] = evecs[i][indices[j]]
    }
    top_evecs.push(row)
  }
  for j = 0; j < k; j = j + 1 {
    top_evals[j] = evals[indices[j]]
  }
  (top_evals, top_evecs)
}

///|
/// Cumulative R2 for each successive component, over non-missing entries.
fn pm_r2cum(
  matrix : Array[Array[Double]],
  scores : Array[Array[Double]],
  loadings : Array[Array[Double]],
) -> Array[Double] {
  let n = matrix.length()
  let d = if n > 0 { matrix[0].length() } else { 0 }
  let k = if loadings.length() > 0 { loadings[0].length() } else { 0 }
  let mut tss = 0.0
  for i = 0; i < n; i = i + 1 {
    for j = 0; j < d; j = j + 1 {
      let v = matrix[i][j]
      if !v.is_nan() {
        tss = tss + v * v
      }
    }
  }
  let out = Array::make(k, 0.0)
  for c = 1; c <= k; c = c + 1 {
    // recon = scores[:,0:c] @ loadings[:,0:c].T
    let mut ss = 0.0
    for i = 0; i < n; i = i + 1 {
      for j = 0; j < d; j = j + 1 {
        let v = matrix[i][j]
        if !v.is_nan() {
          let mut r = 0.0
          for l = 0; l < c; l = l + 1 {
            r = r + scores[i][l] * loadings[j][l]
          }
          let diff = v - r
          ss = ss + diff * diff
        }
      }
    }
    out[c - 1] = 1.0 - ss / tss
  }
  out
}

// ---------------------------------------------------------------------------
// PPCA (Tipping-Bishop 1999), transcribed from pcaMethods ppca.R
// ---------------------------------------------------------------------------

///|
/// Probabilistic PCA with missing-value EM.
pub fn ppca(
  matrix~ : Array[Array[Double]],
  n_pcs~ : Int,
  max_iterations~ : Int,
  threshold~ : Double,
) -> PcaMethodsResult {
  let n = matrix.length()
  let d = if n > 0 { matrix[0].length() } else { 0 }
  let k = n_pcs
  // working copy M with NA -> 0
  let m : Array[Array[Double]] = Array(capacity=0)
  let mut missing = 0
  for i = 0; i < n; i = i + 1 {
    let row = Array::make(d, 0.0)
    for j = 0; j < d; j = j + 1 {
      let v = matrix[i][j]
      if v.is_nan() {
        missing = missing + 1
        row[j] = 0.0
      } else {
        row[j] = v
      }
    }
    m.push(row)
  }
  // deterministic init: C = first k rows transposed (D×k)
  let c0 = pm_transpose(Array::make(0, [])) // placeholder
  let c_init : Array[Array[Double]] = Array(capacity=0)
  for j = 0; j < d; j = j + 1 {
    let row = Array::make(k, 0.0)
    for l = 0; l < k; l = l + 1 {
      row[l] = m[l][j]
    }
    c_init.push(row)
  }
  let mut c = c_init
  let mut ctc = pm_mul(pm_transpose(c), c)
  let mut x = pm_mul(pm_mul(m, c), pm_inverse(ctc))
  let recon = pm_mul(x, pm_transpose(c))
  for i = 0; i < n; i = i + 1 {
    for j = 0; j < d; j = j + 1 {
      if matrix[i][j].is_nan() {
        recon[i][j] = 0.0
      }
    }
  }
  // ss
  let mut ss = 0.0
  for i = 0; i < n; i = i + 1 {
    for j = 0; j < d; j = j + 1 {
      let diff = recon[i][j] - m[i][j]
      ss = ss + diff * diff
    }
  }
  ss = ss / (n * d - missing).to_double()

  let mut count = 1
  let mut old = 1.0e308
  while count > 0 {
    let sx = pm_inverse(pm_plus_identity_scaled(pm_scale(ctc, 1.0 / ss), 1.0))
    let ss_old = ss
    if missing > 0 {
      let proj = pm_mul(x, pm_transpose(c))
      for i = 0; i < n; i = i + 1 {
        for j = 0; j < d; j = j + 1 {
          if matrix[i][j].is_nan() {
            m[i][j] = proj[i][j]
          }
        }
      }
    }
    x = pm_scale(pm_mul(pm_mul(m, c), sx), 1.0 / ss)
    let sum_xtx = pm_mul(pm_transpose(x), x)
    // c = (m^T x) * inv(sum_xtx + N*sx)
    let n_sx = pm_scale(sx, n.to_double())
    c = pm_mul(pm_mul(pm_transpose(m), x), pm_inverse(pm_add(sum_xtx, n_sx)))
    ctc = pm_mul(pm_transpose(c), c)
    // ss = (sum((c x^T - m^T)^2) + N*sum(ctc * sx) + missing*ss_old)/(N*D)
    let diff2 = pm_sub(pm_mul(c, pm_transpose(x)), pm_transpose(m))
    ss = pm_sum_sq(diff2) + n.to_double() * pm_sum_ew(ctc, sx) + missing.to_double() * ss_old
    ss = ss / (n * d).to_double()
    // objective (for convergence)
    let log_det = pm_log_det(sx)
    let objective = n.to_double() * (d.to_double() * @math.ln(ss) + pm_trace(sx) - log_det) +
      pm_trace(sum_xtx) - missing.to_double() * @math.ln(ss_old)
    let rel_ch = (1.0 - objective / old).abs()
    old = objective
    count = count + 1
    if rel_ch < threshold && count > 5 {
      count = 0
    } else if count > max_iterations {
      count = 0
    }
  }
  c = pm_orth(c)
  // rotate via eigen of cov(m c)
  let mc = pm_mul(m, c)
  let cov_mc = pm_covariance(mc)
  let (_, vecs) = pm_top_eigen(cov_mc, k)
  c = pm_mul(c, vecs)
  x = pm_mul(m, c)
  PcaMethodsResult::{
    scores: x,
    loadings: c,
    r2cum: pm_r2cum(matrix, x, c),
    tau: 0.0,
    sigma2: ss,
    imputed: m,
  }
}

// ---------------------------------------------------------------------------
// BPCA (Oba et al. 2003), transcribed from pcaMethods bpca.R
// ---------------------------------------------------------------------------

///|
/// Bayesian PCA with automatic relevance determination and missing-value imputation.
pub fn bpca(
  matrix~ : Array[Array[Double]],
  n_pcs~ : Int,
  max_steps~ : Int,
  threshold~ : Double,
) -> PcaMethodsResult {
  let rows = matrix.length()
  let cols = if rows > 0 { matrix[0].length() } else { 0 }
  let k = n_pcs
  // nan mask, yest (NA -> 0)
  let nan : Array[Array[Bool]] = Array(capacity=0)
  let yest : Array[Array[Double]] = Array(capacity=0)
  for i = 0; i < rows; i = i + 1 {
    let nr = Array::make(cols, false)
    let yr = Array::make(cols, 0.0)
    for j = 0; j < cols; j = j + 1 {
      let v = matrix[i][j]
      if v.is_nan() {
        nr[j] = true
        yr[j] = 0.0
      } else {
        yr[j] = v
      }
    }
    nan.push(nr)
    yest.push(yr)
  }
  // row_miss / row_nomiss
  let row_miss : Array[Int] = []
  let row_nomiss : Array[Int] = []
  for i = 0; i < rows; i = i + 1 {
    let mut has = false
    for j = 0; j < cols; j = j + 1 {
      if nan[i][j] {
        has = true
      }
    }
    if has {
      row_miss.push(i)
    } else {
      row_nomiss.push(i)
    }
  }
  // column means over observed
  let mean = Array::make(cols, 0.0)
  for j = 0; j < cols; j = j + 1 {
    let mut s = 0.0
    let mut cnt = 0
    for i = 0; i < rows; i = i + 1 {
      if !nan[i][j] {
        s = s + matrix[i][j]
        cnt = cnt + 1
      }
    }
    mean[j] = if cnt > 0 { s / cnt.to_double() } else { 0.0 }
  }
  // covariance of yest
  let covy = pm_covariance(yest)
  let (evals, evecs) = pm_top_eigen(covy, k)
  // PA = U sqrt(S)
  let mut pa : Array[Array[Double]] = Array(capacity=0)
  for j = 0; j < cols; j = j + 1 {
    let row = Array::make(k, 0.0)
    for l = 0; l < k; l = l + 1 {
      row[l] = evecs[j][l] * evals[l].sqrt()
    }
    pa.push(row)
  }
  let covy_trace = pm_trace(covy)
  let mut tau = 1.0 / (covy_trace - pm_sum(evals))
  tau = if tau > 1.0e10 { 1.0e10 } else if tau < 1.0e-10 { 1.0e-10 } else { tau }
  let galpha0 = 1.0e-10
  let balpha0 = 1.0
  let alpha = Array::make(k, 0.0)
  let pa_norms = pm_col_norms_sq(pa)
  for l = 0; l < k; l = l + 1 {
    alpha[l] = (2.0 * galpha0 + cols.to_double()) / (tau * pa_norms[l] + 2.0 * galpha0 / balpha0)
  }
  let gmu0 = 0.001
  let btau0 = 1.0
  let gtau0 = 1.0e-10
  let sigw = pm_identity(k)
  let mut tauold = 1000.0

  for step = 1; step <= max_steps + 1; step = step + 1 {
    let rx = pm_identity(k)
    // rx = I + tau * PA^T PA + SigW
    let patpa = pm_mul(pm_transpose(pa), pa)
    for i = 0; i < k; i = i + 1 {
      for j = 0; j < k; j = j + 1 {
        rx[i][j] = rx[i][j] + tau * patpa[i][j] + sigw[i][j]
      }
    }
    let rxinv = pm_inverse(rx)
    let tacc : Array[Array[Double]] = Array(capacity=0)
    for _j = 0; _j < cols; _j = _j + 1 {
      tacc.push(Array::make(k, 0.0))
    }
    let mut trs = 0.0
    // complete rows
    if row_nomiss.length() > 0 {
      // dy = y[complete] - mean  (n_c × d)
      let dy : Array[Array[Double]] = Array(capacity=0)
      for idx in row_nomiss {
        let row = Array::make(cols, 0.0)
        for j = 0; j < cols; j = j + 1 {
          row[j] = matrix[idx][j] - mean[j]
        }
        dy.push(row)
      }
      // x = tau * Rxinv * PA^T * dy^T  (k × n_c)
      let x = pm_scale(pm_mul(pm_mul(rxinv, pm_transpose(pa)), pm_transpose(dy)), tau)
      // T = dy^T x^T  (d × k)
      let t = pm_mul(pm_transpose(dy), pm_transpose(x))
      for j = 0; j < cols; j = j + 1 {
        for l = 0; l < k; l = l + 1 {
          tacc[j][l] = tacc[j][l] + t[j][l]
        }
      }
      trs = trs + pm_sum_sq(dy)
    }
    // incomplete rows
    for idx in row_miss {
      // dyo = observed entries of (row_i - mean), row vector 1×d_o
      let dyo_vals : Array[Double] = []
      let obs : Array[Int] = []
      for j = 0; j < cols; j = j + 1 {
        if !nan[idx][j] {
          obs.push(j)
          dyo_vals.push(matrix[idx][j] - mean[j])
        }
      }
      let d_o = obs.length()
      // Wm = PA[missing rows] (d_m × k), Wo = PA[observed rows] (d_o × k)
      let wm : Array[Array[Double]] = Array(capacity=0)
      let wo : Array[Array[Double]] = Array(capacity=0)
      for j = 0; j < cols; j = j + 1 {
        if nan[idx][j] {
          wm.push(pa[j])
        } else {
          wo.push(pa[j])
        }
      }
      // rxinv2 = inv(rx - tau * Wm^T Wm)
      let wmtwm = pm_mul(pm_transpose(wm), wm)
      let rx2 : Array[Array[Double]] = Array(capacity=0)
      for i = 0; i < k; i = i + 1 {
        let row = Array::make(k, 0.0)
        for j = 0; j < k; j = j + 1 {
          row[j] = rx[i][j] - tau * wmtwm[i][j]
        }
        rx2.push(row)
      }
      let rxinv2 = pm_inverse(rx2)
      // ex = tau * Wo^T dyo^T (k×1)
      let dyo_mat : Array[Array[Double]] = Array(capacity=0)
      dyo_mat.push(dyo_vals)  // 1×d_o
      let ex = pm_scale(pm_mul(pm_transpose(wo), pm_transpose(dyo_mat)), tau)
      let xx = pm_mul(rxinv2, ex)  // k×1
      // dym = Wm xx (d_m×1)
      let dym = pm_mul(wm, xx)     // d_m×1
      // full dy row (1×d): observed + imputed
      let dyrow = Array::make(cols, 0.0)
      let mut mcnt = 0
      for j = 0; j < cols; j = j + 1 {
        if nan[idx][j] {
          dyrow[j] = dym[mcnt][0]
          mcnt = mcnt + 1
        } else {
          dyrow[j] = matrix[idx][j] - mean[j]
        }
      }
      // yest[idx] = dy + mean
      for j = 0; j < cols; j = j + 1 {
        yest[idx][j] = dyrow[j] + mean[j]
      }
      // T += dy^T xx^T  (d×k)
      for j = 0; j < cols; j = j + 1 {
        for l = 0; l < k; l = l + 1 {
          tacc[j][l] = tacc[j][l] + dyrow[j] * xx[l][0]
        }
      }
      // T[missing,] += Wm Rxinv2  (d_m × k)
      let wmr = pm_mul(wm, rxinv2)
      mcnt = 0
      for j = 0; j < cols; j = j + 1 {
        if nan[idx][j] {
          for l = 0; l < k; l = l + 1 {
            tacc[j][l] = tacc[j][l] + wmr[mcnt][l]
          }
          mcnt = mcnt + 1
        }
      }
      // trS += sum(dy^2) + sum(missing)/tau + trace(Wm rxinv2 Wm^T)
      let mut n_miss = 0
      for j = 0; j < cols; j = j + 1 {
        if nan[idx][j] {
          n_miss = n_miss + 1
        }
      }
      let mut dysq = 0.0
      for j = 0; j < cols; j = j + 1 {
        dysq = dysq + dyrow[j] * dyrow[j]
      }
      let wmrxinv2wmt = pm_mul(pm_mul(wm, rxinv2), pm_transpose(wm))
      trs = trs + dysq + n_miss.to_double() / tau + pm_trace(wmrxinv2wmt)
    }
    // T /= rows; trS /= rows
    for j = 0; j < cols; j = j + 1 {
      for l = 0; l < k; l = l + 1 {
        tacc[j][l] = tacc[j][l] / rows.to_double()
      }
    }
    trs = trs / rows.to_double()
    // M-step
    let rxinv_full = pm_inverse(rx)
    // Dw = Rxinv + tau * T^T PA Rxinv + diag(alpha)/rows
    let tpat = pm_transpose(tacc)
    // hoisted: the k×k term T^T · PA · Rxinv is constant across (i, j)
    let tpat_pa_rxinv = pm_mul(tpat, pm_mul(pa, rxinv_full))
    let dw : Array[Array[Double]] = Array(capacity=0)
    for i = 0; i < k; i = i + 1 {
      let row = Array::make(k, 0.0)
      for j = 0; j < k; j = j + 1 {
        row[j] = rxinv_full[i][j] + tau * tpat_pa_rxinv[i][j] + (if i == j { alpha[i] } else { 0.0 }) / rows.to_double()
      }
      dw.push(row)
    }
    let dwinv = pm_inverse(dw)
    pa = pm_mul(tacc, dwinv)
    // tau = (cols + 2*gtau0/rows) / (trS - trace(T^T PA) + (mean·mean^T*gmu0 + 2*gtau0/btau0)/rows)
    let tpapa = pm_mul(pm_transpose(tacc), pa)
    let mut mean_sq = 0.0
    for j = 0; j < cols; j = j + 1 {
      mean_sq = mean_sq + mean[j] * mean[j]
    }
    let num = cols.to_double() + 2.0 * gtau0 / rows.to_double()
    let den = trs - pm_trace(tpapa) + (mean_sq * gmu0 + 2.0 * gtau0 / btau0) / rows.to_double()
    tau = num / den
    // SigW = Dwinv * (cols/rows)
    for i = 0; i < k; i = i + 1 {
      for j = 0; j < k; j = j + 1 {
        sigw[i][j] = dwinv[i][j] * cols.to_double() / rows.to_double()
      }
    }
    // alpha = (2*galpha0 + cols) / (tau*diag(PA^T PA) + diag(SigW) + 2*galpha0/balpha0)
    let pa_norms2 = pm_col_norms_sq(pa)
    for l = 0; l < k; l = l + 1 {
      alpha[l] = (2.0 * galpha0 + cols.to_double()) / (tau * pa_norms2[l] + sigw[l][l] + 2.0 * galpha0 / balpha0)
    }
    if step % 10 == 0 {
      let dtau = (@math.log10(tau) - @math.log10(tauold)).abs()
      if dtau < threshold {
        break
      }
      tauold = tau
    }
  }
  PcaMethodsResult::{
    scores: pm_mul(yest, pa),
    loadings: pa,
    r2cum: pm_r2cum(matrix, pm_mul(yest, pa), pa),
    tau,
    sigma2: 0.0,
    imputed: yest,
  }
}

// ---------------------------------------------------------------------------
// Small helpers for matrix arithmetic inline with `ppca`/`bpca`
// ---------------------------------------------------------------------------

///|
fn pm_plus_identity_scaled(a : Array[Array[Double]], s : Double) -> Array[Array[Double]] {
  let n = a.length()
  let out : Array[Array[Double]] = Array(capacity=0)
  for i = 0; i < n; i = i + 1 {
    let row = Array::make(n, 0.0)
    for j = 0; j < n; j = j + 1 {
      row[j] = a[i][j]
    }
    row[i] = row[i] + s
    out.push(row)
  }
  out
}

///|
fn pm_scale(a : Array[Array[Double]], s : Double) -> Array[Array[Double]] {
  let out : Array[Array[Double]] = Array(capacity=0)
  for i = 0; i < a.length(); i = i + 1 {
    let row = Array::make(a[i].length(), 0.0)
    for j = 0; j < a[i].length(); j = j + 1 {
      row[j] = a[i][j] * s
    }
    out.push(row)
  }
  out
}

///|
fn pm_add(a : Array[Array[Double]], b : Array[Array[Double]]) -> Array[Array[Double]] {
  let out : Array[Array[Double]] = Array(capacity=0)
  for i = 0; i < a.length(); i = i + 1 {
    let row = Array::make(a[i].length(), 0.0)
    for j = 0; j < a[i].length(); j = j + 1 {
      row[j] = a[i][j] + b[i][j]
    }
    out.push(row)
  }
  out
}

///|
fn pm_sub(a : Array[Array[Double]], b : Array[Array[Double]]) -> Array[Array[Double]] {
  let out : Array[Array[Double]] = Array(capacity=0)
  for i = 0; i < a.length(); i = i + 1 {
    let row = Array::make(a[i].length(), 0.0)
    for j = 0; j < a[i].length(); j = j + 1 {
      row[j] = a[i][j] - b[i][j]
    }
    out.push(row)
  }
  out
}

///|
fn pm_sum(a : Array[Double]) -> Double {
  let mut s = 0.0
  for v in a {
    s = s + v
  }
  s
}

///|
/// log-determinant of a positive-definite-ish matrix via LU (Gauss, partial
/// pivoting). Returns ln(det).
fn pm_log_det(a : Array[Array[Double]]) -> Double {
  let n = a.length()
  let m : Array[Array[Double]] = Array(capacity=0)
  for i = 0; i < n; i = i + 1 {
    let row = Array::make(n, 0.0)
    for j = 0; j < n; j = j + 1 {
      row[j] = a[i][j]
    }
    m.push(row)
  }
  let mut lndet = 0.0
  for col = 0; col < n; col = col + 1 {
    let mut pivot = col
    let mut best = m[col][col].abs()
    for r = col + 1; r < n; r = r + 1 {
      if m[r][col].abs() > best {
        best = m[r][col].abs()
        pivot = r
      }
    }
    if m[col][col] != m[pivot][col] {
      let tmp = m[col]
      m[col] = m[pivot]
      m[pivot] = tmp
    }
    let pv = m[col][col]
    lndet = lndet + @math.ln(pv.abs())
    for r = col + 1; r < n; r = r + 1 {
      let factor = m[r][col] / pv
      for j = col; j < n; j = j + 1 {
        m[r][j] = m[r][j] - factor * m[col][j]
      }
    }
  }
  lndet
}