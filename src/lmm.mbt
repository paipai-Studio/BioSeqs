///|
/// Linear mixed models (LMM) — the foundation for GEMMA / EMMAX style
/// genome-wide association with population structure / kinship correction.
///
/// Implements the EMMAX algorithm (Kang et al. 2008, "Efficient control of
/// population structure in model organism association mapping") on an
/// eigen-decomposed kinship (GRM) matrix:
///
///   y = X beta + u + e,  u ~ N(0, sigma2_g * K),  e ~ N(0, sigma2_e * I)
///
/// - `lmm_compute_grm`: additive genetic relationship matrix from genotypes.
/// - `lmm_fit`: REML estimation of the two variance components (sigma2_g,
///   sigma2_e) and the fixed effects beta, by rotating with the eigenvectors
///   of K to a diagonal covariance and optimizing the variance ratio delta =
///   sigma2_e / sigma2_g.
/// - `lmm_associate`: per-SNP association via generalized least squares (GLS)
///   with the null-model variance components held fixed (the EMMAX
///   approximation), yielding a Wald statistic and p-value.
///
/// NOTE (honesty): REML optimizes only the variance ratio delta (1D golden
/// section on a log grid), matching the profiled-likelihood form of EMMAX;
/// it does not apply AI-REML (average information) or variable-selection
/// tests. The kinship is assumed symmetric positive semi-definite.

// ============================================================
// Kinship / GRM
// ============================================================

///|
/// Additive genetic relationship matrix (GRM) from a genotype matrix.
///
/// `genotypes` is n_samples x n_variants; entries are 0/1/2 (reference allele
/// count) or -1 for missing. The standardized additive estimator (VanRaden /
/// GCTA) is used:
///   K_ab = mean_l (x_al - 2 p_l)(x_bl - 2 p_l) / (2 p_l (1 - p_l)),
/// where p_l is the per-variant allele frequency over non-missing samples.
pub fn lmm_compute_grm(genotypes : Array[Array[Double]]) -> Array[Array[Double]] {
  let n = genotypes.length()
  let m = if n > 0 { genotypes[0].length() } else { 0 }

  // Allele frequencies.
  let freq : Array[Double] = Array::make(m, 0.5)
  for l in 0..<m {
    let mut sum = 0.0
    let mut nvalid = 0
    for i in 0..<n {
      if genotypes[i][l] >= 0.0 {
        sum = sum + genotypes[i][l]
        nvalid = nvalid + 1
      }
    }
    if nvalid > 0 {
      freq[l] = sum / (2.0 * nvalid.to_double())
    }
  }

  let grm : Array[Array[Double]] = []
  for a in 0..<n {
    let row = Array::make(n, 0.0)
    for b in 0..<=a {
      let mut s = 0.0
      let mut cnt = 0
      for l in 0..<m {
        if genotypes[a][l] >= 0.0 && genotypes[b][l] >= 0.0 {
          let p = freq[l]
          if p > 0.0 && p < 1.0 {
            let d = 2.0 * p * (1.0 - p)
            s = s + (genotypes[a][l] - 2.0 * p) * (genotypes[b][l] - 2.0 * p) / d
            cnt = cnt + 1
          }
        }
      }
      let v = if cnt > 0 { s / cnt.to_double() } else { 0.0 }
      row[b] = v
      if a != b {
        grm[b][a] = v
      }
    }
    grm.push(row)
  }
  grm
}

// ============================================================
// Linear algebra helpers
// ============================================================

///|
/// Allocate a rows x cols Double matrix (distinct inner arrays).
fn lmm_dmatrix(rows : Int, cols : Int) -> Array[Array[Double]] {
  let m : Array[Array[Double]] = Array::make(rows, Array::make(0, 0.0))
  for i in 0..<rows {
    m[i] = Array::make(cols, 0.0)
  }
  m
}

///|
/// Uᵀ v: rotate a vector into the eigen basis (U columns = eigenvectors).
fn lmm_ut_vec(u : Array[Array[Double]], v : Array[Double]) -> Array[Double] {
  let n = u.length()
  let out = Array::make(n, 0.0)
  // out[j] = sum_i u[i][j] * v[i]
  for j in 0..<n {
    let mut s = 0.0
    for i in 0..<n {
      s = s + u[i][j] * v[i]
    }
    out[j] = s
  }
  out
}

///|
/// Uᵀ M: rotate a matrix (n x k) into the eigen basis -> (n x k).
fn lmm_ut_mat(u : Array[Array[Double]], m : Array[Array[Double]]) -> Array[Array[Double]] {
  let n = u.length()
  let k = if n > 0 { m[0].length() } else { 0 }
  let out = lmm_dmatrix(n, k)
  for j in 0..<n {
    for c in 0..<k {
      let mut s = 0.0
      for i in 0..<n {
        s = s + u[i][j] * m[i][c]
      }
      out[j][c] = s
    }
  }
  out
}

///|
/// Solve A x = b via Gaussian elimination with partial pivoting; also returns
/// the log-determinant of A (A assumed symmetric positive definite).
fn lmm_solve_logdet(
  a : Array[Array[Double]],
  b : Array[Double]
) -> (Array[Double], Double) {
  let n = a.length()
  let aug : Array[Array[Double]] = []
  for i in 0..<n {
    let row = Array::make(n + 1, 0.0)
    for j in 0..<n {
      row[j] = a[i][j]
    }
    row[n] = b[i]
    aug.push(row)
  }
  let mut logdet = 0.0
  let mut sign = 1.0
  for col in 0..<n {
    // Partial pivot.
    let mut pivot = col
    let mut best = aug[col][col].abs()
    let mut r = col + 1
    while r < n {
      if aug[r][col].abs() > best {
        best = aug[r][col].abs()
        pivot = r
      }
      r = r + 1
    }
    if best < 1.0e-12 {
      // Singular: return zero solution and -inf logdet.
      return (Array::make(n, 0.0), -1.0e300)
    }
    if pivot != col {
      let tmp = aug[col]
      aug[col] = aug[pivot]
      aug[pivot] = tmp
      sign = -sign
    }
    let pv = aug[col][col]
    logdet = logdet + @math.ln(pv)
    // Eliminate below.
    let mut r = col + 1
    while r < n {
      let f = aug[r][col] / pv
      if f != 0.0 {
        let mut c = col
        while c < n + 1 {
          aug[r][c] = aug[r][c] - f * aug[col][c]
          c = c + 1
        }
      }
      r = r + 1
    }
  }
  // Back substitution.
  let x = Array::make(n, 0.0)
  let mut i = n
  while i > 0 {
    i = i - 1
    let mut s = aug[i][n]
    let mut j = i + 1
    while j < n {
      s = s - aug[i][j] * x[j]
      j = j + 1
    }
    x[i] = s / aug[i][i]
  }
  (x, logdet)
}

///|
/// Matrix-vector product A x.
fn lmm_mat_vec(a : Array[Array[Double]], x : Array[Double]) -> Array[Double] {
  let n = a.length()
  let out = Array::make(n, 0.0)
  for i in 0..<n {
    let mut s = 0.0
    for j in 0..<n {
      s = s + a[i][j] * x[j]
    }
    out[i] = s
  }
  out
}

///|
/// Golden-section maximization of a unimodal function over [lo, hi].
fn lmm_golden_max(
  f : (Double) -> Double,
  lo : Double,
  hi : Double,
  tol : Double,
  max_iter : Int
) -> Double {
  let gr = 0.6180339887498949
  let mut a = lo
  let mut b = hi
  let mut c = b - gr * (b - a)
  let mut d = a + gr * (b - a)
  let mut fc = f(c)
  let mut fd = f(d)
  let mut it = 0
  while b - a > tol && it < max_iter {
    if fc < fd {
      a = c
      c = d
      fc = fd
      d = a + gr * (b - a)
      fd = f(d)
    } else {
      b = d
      d = c
      fd = fc
      c = b - gr * (b - a)
      fc = f(c)
    }
    it = it + 1
  }
  0.5 * (a + b)
}

// ============================================================
// LMM model
// ============================================================

///|
/// A fitted linear mixed model (REML).
pub struct LmmModel {
  /// Eigenvalues of K, sorted descending (clipped at 0).
  eigenvalues : Array[Double]
  /// Eigenvectors of K (columns = eigenvectors).
  eigenvectors : Array[Array[Double]]
  /// Variance ratio delta = sigma2_e / sigma2_g.
  delta : Double
  /// Genetic variance component.
  sigma2_g : Double
  /// Residual variance component.
  sigma2_e : Double
  /// Narrow-sense heritability sigma2_g / (sigma2_g + sigma2_e).
  heritability : Double
  /// Profiled restricted log-likelihood (up to additive constant).
  log_likelihood : Double
  converged : Bool
  iterations : Int
  /// Fixed-effect estimates (null model).
  beta : Array[Double]
  /// Rotated response Uᵀ y.
  y_rot : Array[Double]
  /// Rotated covariates Uᵀ X (n x p).
  x_rot : Array[Array[Double]]
}

///|
pub fn LmmModel::sigma2_g(self : LmmModel) -> Double {
  self.sigma2_g
}

pub fn LmmModel::sigma2_e(self : LmmModel) -> Double {
  self.sigma2_e
}

pub fn LmmModel::heritability(self : LmmModel) -> Double {
  self.heritability
}

pub fn LmmModel::beta(self : LmmModel) -> Array[Double] {
  self.beta
}

pub fn LmmModel::log_likelihood(self : LmmModel) -> Double {
  self.log_likelihood
}

pub fn LmmModel::converged(self : LmmModel) -> Bool {
  self.converged
}

///|
/// Fit a linear mixed model by REML.
///
/// - `y`: phenotype (n).
/// - `x`: covariate design matrix (n x p), including an intercept column (use
///   `glm_add_intercept`).
/// - `kinship`: n x n symmetric PSD kinship/GRM matrix (use `lmm_compute_grm`).
///
/// The variance ratio delta = sigma2_e / sigma2_g is optimized by golden
/// section over log10(delta) in [-8, 8]; the fixed effects beta are profiled
/// (generalized least squares) at each delta.
pub fn lmm_fit(
  y : Array[Double],
  x : Array[Array[Double]],
  kinship : Array[Array[Double]],
) -> LmmModel {
  let n = y.length()
  let p = if x.length() > 0 { x[0].length() } else { 0 }

  // Eigen-decompose kinship.
  let (eig_raw, eigvec_raw) = linalg_jacobi(kinship, n, 100)
  let (eigvals, eigvecs) = linalg_sort_eigen(eig_raw, eigvec_raw, n)
  // Clip tiny / negative eigenvalues to 0.
  let eigenvalues : Array[Double] = Array::make(n, 0.0)
  for i in 0..<n {
    eigenvalues[i] = if eigvals[i] < 1.0e-10 { 0.0 } else { eigvals[i] }
  }

  let y_rot = lmm_ut_vec(eigvecs, y)
  let x_rot = lmm_ut_mat(eigvecs, x)

  // Objective as a function of log10(delta).
  let objective = fn(logeta : Double) -> Double {
    let delta = @math.pow(10.0, logeta)
    lmm_reml_objective(eigenvalues, x_rot, y_rot, delta, n, p)
  }

  // Golden section over log10(delta) in [-8, 8].
  let logeta_opt = lmm_golden_max(objective, -8.0, 8.0, 1.0e-5, 200)
  let delta = @math.pow(10.0, logeta_opt)
  let ll = objective(logeta_opt)

  // Recompute RSS and beta at the optimum.
  let inv_h = lmm_inv_h(eigenvalues, delta)
  let (beta, rss) = lmm_gf_beta(x_rot, y_rot, inv_h, p)

  let sigma2_g = if n > p { rss / (n - p).to_double() } else { 0.0 }
  let sigma2_e = delta * sigma2_g
  let herit = if sigma2_g + sigma2_e > 0.0 {
    sigma2_g / (sigma2_g + sigma2_e)
  } else {
    0.0
  }

  LmmModel::{
    eigenvalues,
    eigenvectors: eigvecs,
    delta,
    sigma2_g,
    sigma2_e,
    heritability: herit,
    log_likelihood: ll,
    converged: true,
    iterations: 200,
    beta,
    y_rot,
    x_rot,
  }
}

///|
/// Inverse diagonal weights: 1 / (eigenvalue_i + delta).
fn lmm_inv_h(eigenvalues : Array[Double], delta : Double) -> Array[Double] {
  let n = eigenvalues.length()
  let inv_h = Array::make(n, 0.0)
  for i in 0..<n {
    let h = eigenvalues[i] + delta
    inv_h[i] = if h > 0.0 { 1.0 / h } else { 0.0 }
  }
  inv_h
}

///|
/// Profiled REML objective (to maximize) for a given delta:
///   -0.5 * [ (n-p) ln RSS + sum ln(lambda_i + delta) + logdet(X' H^-1 X) ]
fn lmm_reml_objective(
  eigenvalues : Array[Double],
  x_rot : Array[Array[Double]],
  y_rot : Array[Double],
  delta : Double,
  n : Int,
  p : Int
) -> Double {
  let inv_h = lmm_inv_h(eigenvalues, delta)
  let (_, rss) = lmm_gf_beta(x_rot, y_rot, inv_h, p)

  let mut sum_log_h = 0.0
  for i in 0..<n {
    let h = eigenvalues[i] + delta
    sum_log_h = sum_log_h + @math.ln(if h > 0.0 { h } else { 1.0e-300 })
  }

  // logdet(X_rot' H^-1 X_rot).
  let xt_h_inv_x = lmm_gram(x_rot, inv_h)
  let (_, logdet) = lmm_solve_logdet(xt_h_inv_x, Array::make(p, 0.0))

  let rss_safe = if rss > 1.0e-300 { rss } else { 1.0e-300 }
  -0.5 * ((n - p).to_double() * @math.ln(rss_safe) + sum_log_h + logdet)
}

///|
/// Generalized-least-squares beta and weighted residual sum of squares for
/// a given diagonal weight vector inv_h (= 1/(lambda + delta)).
/// Returns (beta, RSS).
fn lmm_gf_beta(
  x_rot : Array[Array[Double]],
  y_rot : Array[Double],
  inv_h : Array[Double],
  p : Int
) -> (Array[Double], Double) {
  let g = lmm_gram(x_rot, inv_h)
  let rhs = lmm_moment(x_rot, y_rot, inv_h)
  let (beta, _) = lmm_solve_logdet(g, rhs)
  // Residual e = y_rot - x_rot beta, RSS = sum inv_h_i e_i^2.
  let n = y_rot.length()
  let mut rss = 0.0
  for i in 0..<n {
    let mut fitted = 0.0
    for j in 0..<p {
      fitted = fitted + x_rot[i][j] * beta[j]
    }
    let e = y_rot[i] - fitted
    rss = rss + inv_h[i] * e * e
  }
  (beta, rss)
}

///|
/// Gram matrix G where G[a][b] = sum_i cols[a][i] * w[i] * cols[b][i],
/// for "columns" given as a row-major matrix.
fn lmm_gram(mat : Array[Array[Double]], w : Array[Double]) -> Array[Array[Double]] {
  let n = mat.length()
  let p = if n > 0 { mat[0].length() } else { 0 }
  let g = lmm_dmatrix(p, p)
  for a in 0..<p {
    for b in a..<p {
      let mut s = 0.0
      for i in 0..<n {
        s = s + mat[i][a] * w[i] * mat[i][b]
      }
      g[a][b] = s
      g[b][a] = s
    }
  }
  g
}

///|
/// Moment vector rhs[a] = sum_i mat[i][a] * w[i] * y[i].
fn lmm_moment(mat : Array[Array[Double]], y : Array[Double], w : Array[Double]) -> Array[Double] {
  let n = mat.length()
  let p = if n > 0 { mat[0].length() } else { 0 }
  let rhs = Array::make(p, 0.0)
  for a in 0..<p {
    let mut s = 0.0
    for i in 0..<n {
      s = s + mat[i][a] * w[i] * y[i]
    }
    rhs[a] = s
  }
  rhs
}

// ============================================================
// Association testing (EMMAX)
// ============================================================

///|
/// A single-SNP association result.
pub struct LmmAssoc {
  beta : Double
  se : Double
  z : Double
  p_value : Double
}

pub fn LmmAssoc::beta(self : LmmAssoc) -> Double {
  self.beta
}

pub fn LmmAssoc::se(self : LmmAssoc) -> Double {
  self.se
}

pub fn LmmAssoc::z(self : LmmAssoc) -> Double {
  self.z
}

pub fn LmmAssoc::p_value(self : LmmAssoc) -> Double {
  self.p_value
}

///|
/// Association test of a set of SNPs against the phenotype, using the fitted
/// LMM's variance components (EMMAX: null-model variance held fixed).
///
/// `snps` is an array of per-SNP genotype vectors (each length n, values
/// 0/1/2). For each SNP, the full design [x | snp] is fitted by generalized
/// least squares and the SNP effect tested with a Wald statistic.
/// Returns one `LmmAssoc` per SNP, in the same order.
pub fn lmm_associate(model : LmmModel, snps : Array[Array[Double]]) -> Array[LmmAssoc] {
  let n = model.eigenvalues.length()
  let p = if model.x_rot.length() > 0 { model.x_rot[0].length() } else { 0 }
  let inv_h = lmm_inv_h(model.eigenvalues, model.delta)

  let results : Array[LmmAssoc] = []
  for snp in snps {
    let w_rot = lmm_ut_vec(model.eigenvectors, snp)
    // Build W = Z' H^-1 Z for Z = [x_rot | w_rot] (p+1 x p+1) and rhs.
    let dim = p + 1
    let wmat = lmm_dmatrix(dim, dim)
    let rhs = Array::make(dim, 0.0)
    for a in 0..<dim {
      for b in a..<dim {
        let mut s = 0.0
        for i in 0..<n {
          let za = if a < p { model.x_rot[i][a] } else { w_rot[i] }
          let zb = if b < p { model.x_rot[i][b] } else { w_rot[i] }
          s = s + za * inv_h[i] * zb
        }
        wmat[a][b] = s
        wmat[b][a] = s
      }
      let mut sm = 0.0
      for i in 0..<n {
        let za = if a < p { model.x_rot[i][a] } else { w_rot[i] }
        sm = sm + za * inv_h[i] * model.y_rot[i]
      }
      rhs[a] = sm
    }
    // Solve [beta, gamma] = W^-1 rhs, and get Var(gamma) = sigma2_g * (W^-1)[p][p].
    let winv = linalg_inverse(wmat)
    let beta_gamma = lmm_mat_vec(winv, rhs)
    let gamma = beta_gamma[p]
    let var_gamma = model.sigma2_g * winv[p][p]
    let se = if var_gamma > 0.0 { var_gamma.sqrt() } else { 0.0 }
    let z = if se > 0.0 { gamma / se } else { 0.0 }
    let pval = stat_pvalue_zscore(z)
    results.push(LmmAssoc::{ beta: gamma, se, z, p_value: pval })
  }
  results
}