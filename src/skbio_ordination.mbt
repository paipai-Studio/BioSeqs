///|
/// scikit-bio ordination + distance-based beta-diversity helpers.
///
/// Port of `skbio.stats.ordination.center_distance_matrix` / `ca` and
/// `skbio.stats.distance.mantel`. `center_distance_matrix` performs the Gower
/// double-centering used to turn a (possibly non-Euclidean) distance matrix into
/// one suitable for PCoA; `ca` performs Correspondence Analysis (chi-square
/// standardization + SVD); `mantel` tests the correlation between two distance
/// matrices with a permutation-based p-value.

///|
pub struct CaResult {
  eigvals : Array[Double]
  sample_scores : Array[Array[Double]]
  feature_scores : Array[Array[Double]]
  proportion_explained : Array[Double]
}

///|
pub struct MantelResult {
  corr_coeff : Double
  p_value : Double
  n : Int
}

// ---------------------------------------------------------------------------
// Small matrix helpers
// ---------------------------------------------------------------------------

///|
fn skb_transpose(a : Array[Array[Double]]) -> Array[Array[Double]] {
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
fn skb_mul(a : Array[Array[Double]], b : Array[Array[Double]]) -> Array[Array[Double]] {
  let m = a.length()
  let k = if m > 0 { a[0].length() } else { 0 }
  let n = if b.length() > 0 { b[0].length() } else { 0 }
  let c : Array[Array[Double]] = Array(capacity=0)
  for i = 0; i < m; i = i + 1 {
    let row = Array::make(n, 0.0)
    for kk = 0; kk < k; kk = kk + 1 {
      let aik = a[i][kk]
      for j = 0; j < n; j = j + 1 {
        row[j] = row[j] + aik * b[kk][j]
      }
    }
    c.push(row)
  }
  c
}

///|
fn skb_sum(a : Array[Array[Double]]) -> Double {
  let mut s = 0.0
  for i = 0; i < a.length(); i = i + 1 {
    for j = 0; j < a[i].length(); j = j + 1 {
      s = s + a[i][j]
    }
  }
  s
}

///|
fn skb_row_sums(a : Array[Array[Double]]) -> Array[Double] {
  let out : Array[Double] = Array(capacity=0)
  for i = 0; i < a.length(); i = i + 1 {
    let mut s = 0.0
    for j = 0; j < a[i].length(); j = j + 1 {
      s = s + a[i][j]
    }
    out.push(s)
  }
  out
}

///|
fn skb_col_sums(a : Array[Array[Double]]) -> Array[Double] {
  let n = if a.length() > 0 { a[0].length() } else { 0 }
  let out = Array::make(n, 0.0)
  for i = 0; i < a.length(); i = i + 1 {
    for j = 0; j < n; j = j + 1 {
      out[j] = out[j] + a[i][j]
    }
  }
  out
}

// ---------------------------------------------------------------------------
// center_distance_matrix
// ---------------------------------------------------------------------------

///|
/// Gower double-centering of a distance matrix (Numerical Ecology eq. 9.20/9.21).
pub fn skbio_center_distance_matrix(
  matrix : Array[Array[Double]],
) -> Array[Array[Double]] {
  let n = matrix.length()
  // E = -D^2 / 2
  let e : Array[Array[Double]] = Array(capacity=0)
  for i = 0; i < n; i = i + 1 {
    let row = Array::make(n, 0.0)
    for j = 0; j < n; j = j + 1 {
      row[j] = -(matrix[i][j] * matrix[i][j]) / 2.0
    }
    e.push(row)
  }
  let row_means = skb_row_sums(e).map(fn(x) { x / n.to_double() })
  let col_means = skb_col_sums(e).map(fn(x) { x / n.to_double() })
  let matrix_mean = skb_sum(e) / (n * n).to_double()
  let f : Array[Array[Double]] = Array(capacity=0)
  for i = 0; i < n; i = i + 1 {
    let row = Array::make(n, 0.0)
    for j = 0; j < n; j = j + 1 {
      row[j] = e[i][j] - row_means[i] - col_means[j] + matrix_mean
    }
    f.push(row)
  }
  f
}

// ---------------------------------------------------------------------------
// Correspondence Analysis (ca)
// ---------------------------------------------------------------------------

///|
/// Sort eigenvalues descending and reorder eigenvector columns accordingly.
fn skb_sort_eigen(
  evals : Array[Double],
  evecs : Array[Array[Double]],
) -> (Array[Double], Array[Array[Double]]) {
  let n = evals.length()
  let idx = Array::make(n, 0)
  for i = 0; i < n; i = i + 1 {
    idx[i] = i
  }
  for i = 0; i < n - 1; i = i + 1 {
    for j = i + 1; j < n; j = j + 1 {
      if evals[idx[j]] > evals[idx[i]] {
        let t = idx[i]
        idx[i] = idx[j]
        idx[j] = t
      }
    }
  }
  let sorted_evals = Array::make(n, 0.0)
  for i = 0; i < n; i = i + 1 {
    sorted_evals[i] = evals[idx[i]]
  }
  let sorted_evecs : Array[Array[Double]] = Array(capacity=0)
  for i = 0; i < n; i = i + 1 {
    let row = Array::make(n, 0.0)
    for j = 0; j < n; j = j + 1 {
      row[j] = evecs[i][idx[j]]
    }
    sorted_evecs.push(row)
  }
  (sorted_evals, sorted_evecs)
}

///|
/// SVD of Q_bar via eigen-decomposition of Q_bar Q_bar^T.
/// Returns (U_hat: n×r, W: r, U: m×r) such that Q_bar = U_hat diag(W) U^T.
fn skb_svd_qbar(
  qbar : Array[Array[Double]],
) -> (Array[Array[Double]], Array[Double], Array[Array[Double]]) {
  let n = qbar.length()
  let m = if n > 0 { qbar[0].length() } else { 0 }
  let qbar_t = skb_transpose(qbar)
  let s = skb_mul(qbar, qbar_t)  // n × n
  let (evals, evecs) = eigh_symmetric(s, n)
  let (sorted_evals, sorted_evecs) = skb_sort_eigen(evals, evecs)

  // determine rank (drop numerically-zero singular values)
  let mut rank = n
  if n > 0 {
    let maxval = sorted_evals[0]
    let mut r = 0
    while r < n && sorted_evals[r] > maxval * 1.0e-13 {
      r = r + 1
    }
    rank = if r == 0 { 1 } else { r }
  }

  let w = Array::make(rank, 0.0)
  for k = 0; k < rank; k = k + 1 {
    w[k] = sorted_evals[k].sqrt()
  }

  // U_hat = sorted_evecs[:, 0:rank]  (n × rank)
  let u_hat : Array[Array[Double]] = Array(capacity=0)
  for i = 0; i < n; i = i + 1 {
    let row = Array::make(rank, 0.0)
    for k = 0; k < rank; k = k + 1 {
      row[k] = sorted_evecs[i][k]
    }
    u_hat.push(row)
  }

  // U = Q_bar^T U_hat / W  (m × rank)
  let u : Array[Array[Double]] = Array(capacity=0)
  for j = 0; j < m; j = j + 1 {
    let row = Array::make(rank, 0.0)
    u.push(row)
  }
  for k = 0; k < rank; k = k + 1 {
    // col = Q_bar^T @ U_hat[:,k]
    let qtu : Array[Double] = Array(capacity=0)
    for j = 0; j < m; j = j + 1 {
      let mut s2 = 0.0
      for i = 0; i < n; i = i + 1 {
        s2 = s2 + qbar_t[j][i] * u_hat[i][k]
      }
      qtu.push(s2)
    }
    let wk = w[k]
    for j = 0; j < m; j = j + 1 {
      u[j][k] = if wk > 1.0e-15 { qtu[j] / wk } else { 0.0 }
    }
  }

  (u_hat, w, u)
}

///|
/// Correspondence Analysis of a contingency table.
pub fn skbio_ca(x : Array[Array[Double]], scaling : Int) -> CaResult {
  let n = x.length()
  let m = if n > 0 { x[0].length() } else { 0 }
  let grand_total = skb_sum(x)
  // Q = x / grand_total
  let q : Array[Array[Double]] = Array(capacity=0)
  for i = 0; i < n; i = i + 1 {
    let row = Array::make(m, 0.0)
    for j = 0; j < m; j = j + 1 {
      row[j] = x[i][j] / grand_total
    }
    q.push(row)
  }
  let row_marginals = skb_row_sums(q)
  let col_marginals = skb_col_sums(q)
  // Q_bar = (Q - expected) / sqrt(expected)
  let qbar : Array[Array[Double]] = Array(capacity=0)
  for i = 0; i < n; i = i + 1 {
    let row = Array::make(m, 0.0)
    for j = 0; j < m; j = j + 1 {
      let expected = row_marginals[i] * col_marginals[j]
      row[j] = (q[i][j] - expected) / expected.sqrt()
    }
    qbar.push(row)
  }
  let (u_hat, w, u) = skb_svd_qbar(qbar)
  let rank = w.length()

  // V = col_marginals^-1/2 * U ; V_hat = row_marginals^-1/2 * U_hat
  let v : Array[Array[Double]] = Array(capacity=0)
  for j = 0; j < m; j = j + 1 {
    let row = Array::make(rank, 0.0)
    let denom = col_marginals[j].sqrt()
    for k = 0; k < rank; k = k + 1 {
      row[k] = if denom > 1.0e-15 { u[j][k] / denom } else { 0.0 }
    }
    v.push(row)
  }
  let v_hat : Array[Array[Double]] = Array(capacity=0)
  for i = 0; i < n; i = i + 1 {
    let row = Array::make(rank, 0.0)
    let denom = row_marginals[i].sqrt()
    for k = 0; k < rank; k = k + 1 {
      row[k] = if denom > 1.0e-15 { u_hat[i][k] / denom } else { 0.0 }
    }
    v_hat.push(row)
  }
  // F = V_hat * W ; F_hat = V * W
  let f : Array[Array[Double]] = Array(capacity=0)
  for i = 0; i < n; i = i + 1 {
    let row = Array::make(rank, 0.0)
    for k = 0; k < rank; k = k + 1 {
      row[k] = v_hat[i][k] * w[k]
    }
    f.push(row)
  }
  let f_hat : Array[Array[Double]] = Array(capacity=0)
  for j = 0; j < m; j = j + 1 {
    let row = Array::make(rank, 0.0)
    for k = 0; k < rank; k = k + 1 {
      row[k] = v[j][k] * w[k]
    }
    f_hat.push(row)
  }

  let eigvals = Array::make(rank, 0.0)
  for k = 0; k < rank; k = k + 1 {
    eigvals[k] = w[k] * w[k]
  }
  let eig_sum = skb_sum([eigvals])
  let proportion_explained = eigvals.map(fn(v) { v / eig_sum })

  let sample_scores = if scaling == 1 { f } else { v_hat }
  let feature_scores = if scaling == 1 { v } else { f_hat }

  CaResult::{
    eigvals,
    sample_scores,
    feature_scores,
    proportion_explained,
  }
}

// ---------------------------------------------------------------------------
// Mantel test
// ---------------------------------------------------------------------------

///|
/// Condensed (upper-triangle) form of a square distance matrix.
fn skb_condensed(matrix : Array[Array[Double]]) -> Array[Double] {
  let n = matrix.length()
  let out : Array[Double] = Array(capacity=0)
  for i = 0; i < n; i = i + 1 {
    for j = i + 1; j < n; j = j + 1 {
      out.push(matrix[i][j])
    }
  }
  out
}

///|
fn skb_mean(a : Array[Double]) -> Double {
  let mut s = 0.0
  for v in a {
    s = s + v
  }
  s / a.length().to_double()
}

///|
/// Pearson correlation coefficient of two vectors.
fn skb_pearson(a : Array[Double], b : Array[Double]) -> Double {
  let ma = skb_mean(a)
  let mb = skb_mean(b)
  let mut dot = 0.0
  let mut na = 0.0
  let mut nb = 0.0
  for i = 0; i < a.length(); i = i + 1 {
    let da = a[i] - ma
    let db = b[i] - mb
    dot = dot + da * db
    na = na + da * da
    nb = nb + db * db
  }
  if na == 0.0 || nb == 0.0 {
    return 0.0 / 0.0
  }
  let r = dot / (na.sqrt() * nb.sqrt())
  if r > 1.0 {
    1.0
  } else if r < -1.0 {
    -1.0
  } else {
    r
  }
}

///|
/// Deterministic LCG random source for the permutation test.
fn skb_lcg(seed : Int) -> Int {
  (seed * 1103515245 + 12345) & 0x7FFFFFFF
}

///|
/// Reorder a square distance matrix by permuting rows/columns together.
fn skb_permute_matrix(
  matrix : Array[Array[Double]],
  order : Array[Int],
) -> Array[Array[Double]] {
  let n = matrix.length()
  let out : Array[Array[Double]] = Array(capacity=0)
  for i = 0; i < n; i = i + 1 {
    let row = Array::make(n, 0.0)
    for j = 0; j < n; j = j + 1 {
      row[j] = matrix[order[i]][order[j]]
    }
    out.push(row)
  }
  out
}

///|
/// Mantel test: correlation coefficient between two distance matrices plus a
/// permutation p-value (two-sided). `permutations` is the number of
/// permutations; `seed` makes the permutation reproducible.
pub fn skbio_mantel(
  x : Array[Array[Double]],
  y : Array[Array[Double]],
  permutations : Int,
  seed : Int,
) -> MantelResult {
  let n = x.length()
  let y_flat = skb_condensed(y)
  let x_flat = skb_condensed(x)
  let corr = skb_pearson(x_flat, y_flat)

  // permutation test (two-sided)
  let mut count_better = 0
  let mut rng = seed
  let base_order = Array::make(n, 0)
  for i = 0; i < n; i = i + 1 {
    base_order[i] = i
  }
  for _p = 0; _p < permutations; _p = _p + 1 {
    // Fisher-Yates shuffle with LCG
    let order = Array::make(n, 0)
    for i = 0; i < n; i = i + 1 {
      order[i] = base_order[i]
    }
    let mut i = n - 1
    while i > 0 {
      rng = skb_lcg(rng)
      let j = rng % (i + 1)
      let t = order[i]
      order[i] = order[j]
      order[j] = t
      i = i - 1
    }
    let xp = skb_permute_matrix(x, order)
    let xp_flat = skb_condensed(xp)
    let rp = skb_pearson(xp_flat, y_flat)
    if rp.abs() >= corr.abs() {
      count_better = count_better + 1
    }
  }
  let p_value = if permutations == 0 {
    0.0 / 0.0
  } else {
    (count_better + 1).to_double() / (permutations + 1).to_double()
  }
  MantelResult::{ corr_coeff: corr, p_value, n }
}
///|
pub struct RdaResult {
  eigvals : Array[Double]
  proportion_explained : Array[Double]
  sample_scores : Array[Array[Double]]
  feature_scores : Array[Array[Double]]
}

///|
/// Column-centering (subtract each column's mean).
fn skb_center_columns(a : Array[Array[Double]]) -> Array[Array[Double]] {
  let n = a.length()
  let p = if n > 0 { a[0].length() } else { 0 }
  let means = skb_col_sums(a).map(fn(x) { x / n.to_double() })
  let out : Array[Array[Double]] = Array(capacity=0)
  for i = 0; i < n; i = i + 1 {
    let row = Array::make(p, 0.0)
    for j = 0; j < p; j = j + 1 {
      row[j] = a[i][j] - means[j]
    }
    out.push(row)
  }
  out
}

///|
/// Redundancy Analysis (RDA). `y` is the response table (samples × features),
/// `x` the explanatory table (samples × constraints). `scaling` is 1 or 2.
pub fn skbio_rda(
  y : Array[Array[Double]],
  x : Array[Array[Double]],
  scaling : Int,
) -> RdaResult {
  let n = y.length()
  let p = if n > 0 { y[0].length() } else { 0 }
  let y_c = skb_center_columns(y)
  let x_c = skb_center_columns(x)
  // B = (X'X)^-1 X'Y ; Y_hat = X B
  let xt = skb_transpose(x_c)
  let xtx_inv = linalg_inverse(skb_mul(xt, x_c))
  let b = skb_mul(xtx_inv, skb_mul(xt, y_c))
  let y_hat = skb_mul(x_c, b)
  // SVD of Y_hat -> constrained axes
  let (_l, s, u) = skb_svd_qbar(y_hat)
  let rank = s.length()
  let f = skb_mul(y_c, u)
  let z = skb_mul(y_hat, u)
  // residual
  let y_res : Array[Array[Double]] = Array(capacity=0)
  for i = 0; i < n; i = i + 1 {
    let row = Array::make(p, 0.0)
    for j = 0; j < p; j = j + 1 {
      row[j] = y_c[i][j] - y_hat[i][j]
    }
    y_res.push(row)
  }
  let (_l2, s_res, u_res) = skb_svd_qbar(y_res)
  let rank_res = s_res.length()
  let f_res = skb_mul(y_res, u_res)
  // concatenated eigenvalues (constrained block, then residual block)
  let eigvals : Array[Double] = Array(capacity=0)
  for v in s {
    eigvals.push(v)
  }
  for v in s_res {
    eigvals.push(v)
  }
  let mut sq_sum = 0.0
  let mut eig_sum = 0.0
  for v in eigvals {
    sq_sum = sq_sum + v * v
    eig_sum = eig_sum + v
  }
  // RDA's proportion explained is eigenvalue / sum(eigenvalues)
  let proportion_explained = eigvals.map(fn(v) { v / eig_sum })
  let cst = @math.pow(sq_sum, 0.25)
  let k_total = eigvals.length()
  let sf = Array::make(k_total, 0.0)
  for i = 0; i < k_total; i = i + 1 {
    sf[i] = if scaling == 1 { cst } else { eigvals[i] / cst }
  }
  // feature_scores = hstack(U, U_res) * sf ; sample_scores = hstack(F, F_res) / sf
  let feature_scores : Array[Array[Double]] = Array(capacity=0)
  for i = 0; i < p; i = i + 1 {
    let row : Array[Double] = Array(capacity=0)
    for k = 0; k < rank; k = k + 1 {
      row.push(u[i][k] * sf[k])
    }
    for k = 0; k < rank_res; k = k + 1 {
      row.push(u_res[i][k] * sf[rank + k])
    }
    feature_scores.push(row)
  }
  let sample_scores : Array[Array[Double]] = Array(capacity=0)
  for i = 0; i < n; i = i + 1 {
    let row : Array[Double] = Array(capacity=0)
    for k = 0; k < rank; k = k + 1 {
      row.push(f[i][k] / sf[k])
    }
    for k = 0; k < rank_res; k = k + 1 {
      row.push(f_res[i][k] / sf[rank + k])
    }
    sample_scores.push(row)
  }
  RdaResult::{ eigvals, proportion_explained, sample_scores, feature_scores }
}

///|
pub struct CcaResult {
  eigvals : Array[Double]
  proportion_explained : Array[Double]
  sample_scores : Array[Array[Double]]
  feature_scores : Array[Array[Double]]
}

///|
/// Weighted centering + weighted standard deviation scaling (column-wise).
fn skb_scale_weighted(
  a : Array[Array[Double]],
  weights : Array[Double],
) -> Array[Array[Double]] {
  let n = a.length()
  let p = if n > 0 { a[0].length() } else { 0 }
  let mut w_sum = 0.0
  for w in weights {
    w_sum = w_sum + w
  }
  let means = Array::make(p, 0.0)
  for j = 0; j < p; j = j + 1 {
    let mut s = 0.0
    for i = 0; i < n; i = i + 1 {
      s = s + weights[i] * a[i][j]
    }
    means[j] = s / w_sum
  }
  let stds = Array::make(p, 0.0)
  for j = 0; j < p; j = j + 1 {
    let mut s = 0.0
    for i = 0; i < n; i = i + 1 {
      let d = a[i][j] - means[j]
      s = s + weights[i] * d * d
    }
    stds[j] = (s / w_sum).sqrt()
  }
  let out : Array[Array[Double]] = Array(capacity=0)
  for i = 0; i < n; i = i + 1 {
    let row = Array::make(p, 0.0)
    for j = 0; j < p; j = j + 1 {
      let std = if stds[j] == 0.0 { 1.0 } else { stds[j] }
      row[j] = (a[i][j] - means[j]) / std
    }
    out.push(row)
  }
  out
}

///|
/// Canonical Correspondence Analysis (CCA).
pub fn skbio_cca(
  y : Array[Array[Double]],
  x : Array[Array[Double]],
  scaling : Int,
) -> CcaResult {
  let n = y.length()
  let p = if n > 0 { y[0].length() } else { 0 }
  let m = if x.length() > 0 { x[0].length() } else { 0 }
  // chi-square standardization of Y
  let grand_total = skb_sum(y)
  let q : Array[Array[Double]] = Array(capacity=0)
  for i = 0; i < n; i = i + 1 {
    let row = Array::make(p, 0.0)
    for j = 0; j < p; j = j + 1 {
      row[j] = y[i][j] / grand_total
    }
    q.push(row)
  }
  let row_marginals = skb_row_sums(q)
  let column_marginals = skb_col_sums(q)
  let q_bar : Array[Array[Double]] = Array(capacity=0)
  for i = 0; i < n; i = i + 1 {
    let row = Array::make(p, 0.0)
    for j = 0; j < p; j = j + 1 {
      let expected = row_marginals[i] * column_marginals[j]
      row[j] = (q[i][j] - expected) / expected.sqrt()
    }
    q_bar.push(row)
  }
  // weighted scale of X
  let x_scaled = skb_scale_weighted(x, row_marginals)
  let x_weighted : Array[Array[Double]] = Array(capacity=0)
  for i = 0; i < n; i = i + 1 {
    let row = Array::make(m, 0.0)
    for j = 0; j < m; j = j + 1 {
      row[j] = row_marginals[i].sqrt() * x_scaled[i][j]
    }
    x_weighted.push(row)
  }
  // B = (Xw'Xw)^-1 Xw' Q_bar ; Y_hat = Xw B
  let xt = skb_transpose(x_weighted)
  let xtx_inv = linalg_inverse(skb_mul(xt, x_weighted))
  let b = skb_mul(xtx_inv, skb_mul(xt, q_bar))
  let y_hat = skb_mul(x_weighted, b)
  let y_res : Array[Array[Double]] = Array(capacity=0)
  for i = 0; i < n; i = i + 1 {
    let row = Array::make(p, 0.0)
    for j = 0; j < p; j = j + 1 {
      row[j] = q_bar[i][j] - y_hat[i][j]
    }
    y_res.push(row)
  }
  // SVD of Y_hat
  let (_l, s, u) = skb_svd_qbar(y_hat)
  let rank = s.length()
  let qu = skb_mul(q_bar, u)
  let u_hat : Array[Array[Double]] = Array(capacity=0)
  for i = 0; i < n; i = i + 1 {
    let row = Array::make(rank, 0.0)
    for k = 0; k < rank; k = k + 1 {
      row[k] = qu[i][k] / s[k]
    }
    u_hat.push(row)
  }
  // SVD of Y_res
  let (_l2, s_res, u_res) = skb_svd_qbar(y_res)
  let rank_res = s_res.length()
  let yru = skb_mul(y_res, u_res)
  let u_hat_res : Array[Array[Double]] = Array(capacity=0)
  for i = 0; i < n; i = i + 1 {
    let row = Array::make(rank_res, 0.0)
    for k = 0; k < rank_res; k = k + 1 {
      row[k] = yru[i][k] / s_res[k]
    }
    u_hat_res.push(row)
  }
  // eigenvalues = squared singular values, concatenated
  let eigvals : Array[Double] = Array(capacity=0)
  for v in s {
    eigvals.push(v * v)
  }
  for v in s_res {
    eigvals.push(v * v)
  }
  let mut eig_sum = 0.0
  for v in eigvals {
    eig_sum = eig_sum + v
  }
  let proportion_explained = eigvals.map(fn(v) { v / eig_sum })
  // V = col_marg^-1/2 * U ; V_hat = row_marg^-1/2 * U_hat ; F = V_hat * s ; F_hat = V * s
  let v : Array[Array[Double]] = Array(capacity=0)
  for j = 0; j < p; j = j + 1 {
    let row = Array::make(rank, 0.0)
    for k = 0; k < rank; k = k + 1 {
      row[k] = u[j][k] / column_marginals[j].sqrt()
    }
    v.push(row)
  }
  let v_hat : Array[Array[Double]] = Array(capacity=0)
  for i = 0; i < n; i = i + 1 {
    let row = Array::make(rank, 0.0)
    for k = 0; k < rank; k = k + 1 {
      row[k] = u_hat[i][k] / row_marginals[i].sqrt()
    }
    v_hat.push(row)
  }
  let f : Array[Array[Double]] = Array(capacity=0)
  for i = 0; i < n; i = i + 1 {
    let row = Array::make(rank, 0.0)
    for k = 0; k < rank; k = k + 1 {
      row[k] = v_hat[i][k] * s[k]
    }
    f.push(row)
  }
  let f_hat : Array[Array[Double]] = Array(capacity=0)
  for j = 0; j < p; j = j + 1 {
    let row = Array::make(rank, 0.0)
    for k = 0; k < rank; k = k + 1 {
      row[k] = v[j][k] * s[k]
    }
    f_hat.push(row)
  }
  let v_res : Array[Array[Double]] = Array(capacity=0)
  for j = 0; j < p; j = j + 1 {
    let row = Array::make(rank_res, 0.0)
    for k = 0; k < rank_res; k = k + 1 {
      row[k] = u_res[j][k] / column_marginals[j].sqrt()
    }
    v_res.push(row)
  }
  let v_hat_res : Array[Array[Double]] = Array(capacity=0)
  for i = 0; i < n; i = i + 1 {
    let row = Array::make(rank_res, 0.0)
    for k = 0; k < rank_res; k = k + 1 {
      row[k] = u_hat_res[i][k] / row_marginals[i].sqrt()
    }
    v_hat_res.push(row)
  }
  let f_res : Array[Array[Double]] = Array(capacity=0)
  for i = 0; i < n; i = i + 1 {
    let row = Array::make(rank_res, 0.0)
    for k = 0; k < rank_res; k = k + 1 {
      row[k] = v_hat_res[i][k] * s_res[k]
    }
    f_res.push(row)
  }
  let f_hat_res : Array[Array[Double]] = Array(capacity=0)
  for j = 0; j < p; j = j + 1 {
    let row = Array::make(rank_res, 0.0)
    for k = 0; k < rank_res; k = k + 1 {
      row[k] = v_res[j][k] * s_res[k]
    }
    f_hat_res.push(row)
  }
  // concat feature/sample scores per scaling
  let feature_scores : Array[Array[Double]] = Array(capacity=0)
  let sample_scores : Array[Array[Double]] = Array(capacity=0)
  for j = 0; j < p; j = j + 1 {
    let row : Array[Double] = Array(capacity=0)
    let src_f = if scaling == 1 { v } else { f_hat }
    let src_r = if scaling == 1 { v_res } else { f_hat_res }
    for k = 0; k < rank; k = k + 1 {
      row.push(src_f[j][k])
    }
    for k = 0; k < rank_res; k = k + 1 {
      row.push(src_r[j][k])
    }
    feature_scores.push(row)
  }
  for i = 0; i < n; i = i + 1 {
    let row : Array[Double] = Array(capacity=0)
    let src_f = if scaling == 1 { f } else { v_hat }
    let src_r = if scaling == 1 { f_res } else { v_hat_res }
    for k = 0; k < rank; k = k + 1 {
      row.push(src_f[i][k])
    }
    for k = 0; k < rank_res; k = k + 1 {
      row.push(src_r[i][k])
    }
    sample_scores.push(row)
  }
  CcaResult::{ eigvals, proportion_explained, sample_scores, feature_scores }
}

///|
/// e_matrix: `E = -D^2 / 2` (element-wise), the first half of Gower centering.
pub fn skbio_e_matrix(dm : Array[Array[Double]]) -> Array[Array[Double]] {
  let n = dm.length()
  let out : Array[Array[Double]] = Array(capacity=0)
  for i = 0; i < n; i = i + 1 {
    let row = Array::make(n, 0.0)
    for j = 0; j < n; j = j + 1 {
      row[j] = dm[i][j] * dm[i][j] / -2.0
    }
    out.push(row)
  }
  out
}

///|
/// f_matrix: Gower double centering `E - row_mean - col_mean + grand_mean`.
pub fn skbio_f_matrix(e : Array[Array[Double]]) -> Array[Array[Double]] {
  let n = e.length()
  let row_means = skb_row_sums(e).map(fn(x) { x / n.to_double() })
  let col_means = skb_col_sums(e).map(fn(x) { x / n.to_double() })
  let grand_mean = skb_sum(e) / (n * n).to_double()
  let out : Array[Array[Double]] = Array(capacity=0)
  for i = 0; i < n; i = i + 1 {
    let row = Array::make(n, 0.0)
    for j = 0; j < n; j = j + 1 {
      row[j] = e[i][j] - row_means[i] - col_means[j] + grand_mean
    }
    out.push(row)
  }
  out
}

///|
/// Column-wise mean and standard deviation (unweighted, axis=0).
pub fn skbio_mean_and_std(a : Array[Array[Double]], ddof : Int) -> (Array[Double], Array[Double]) {
  let n = a.length()
  let p = if n > 0 { a[0].length() } else { 0 }
  let means = Array::make(p, 0.0)
  let stds = Array::make(p, 0.0)
  for j = 0; j < p; j = j + 1 {
    let mut s = 0.0
    for i = 0; i < n; i = i + 1 {
      s = s + a[i][j]
    }
    means[j] = s / n.to_double()
    let mut v = 0.0
    for i = 0; i < n; i = i + 1 {
      let d = a[i][j] - means[j]
      v = v + d * d
    }
    stds[j] = (v / (n - ddof).to_double()).sqrt()
  }
  (means, stds)
}

///|
/// Column-wise centering plus (optional) scaling to unit variance.
pub fn skbio_scale(
  a : Array[Array[Double]],
  with_mean : Bool,
  with_std : Bool,
  ddof : Int,
) -> Array[Array[Double]] {
  let n = a.length()
  let p = if n > 0 { a[0].length() } else { 0 }
  let (means, stds) = skbio_mean_and_std(a, ddof)
  let scaled_stds = Array::make(p, 0.0)
  for j = 0; j < p; j = j + 1 {
    scaled_stds[j] = if stds[j] == 0.0 { 1.0 } else { stds[j] }
  }
  let out : Array[Array[Double]] = Array(capacity=0)
  for i = 0; i < n; i = i + 1 {
    let row = Array::make(p, 0.0)
    for j = 0; j < p; j = j + 1 {
      let mut v = a[i][j]
      if with_mean {
        v = v - means[j]
      }
      if with_std {
        v = v / scaled_stds[j]
      }
      row[j] = v
    }
    out.push(row)
  }
  out
}

///|
/// Pearson correlation matrix between the columns of `x` and `y`.
pub fn skbio_corr(x : Array[Array[Double]], y : Array[Array[Double]]) -> Array[Array[Double]] {
  let n = x.length()
  let xs = skbio_scale(x, true, true, 0)
  let ys = skbio_scale(y, true, true, 0)
  let prod = skb_mul(skb_transpose(xs), ys)
  let out : Array[Array[Double]] = Array(capacity=0)
  for i = 0; i < prod.length(); i = i + 1 {
    let row = Array::make(prod[0].length(), 0.0)
    for j = 0; j < prod[0].length(); j = j + 1 {
      row[j] = prod[i][j] / n.to_double()
    }
    out.push(row)
  }
  out
}

///|
/// Numerical rank: count of singular values above the tolerance
/// `S.max() * max(m, n) * eps`.
pub fn skbio_svd_rank(s : Array[Double], m : Int, n : Int) -> Int {
  let mut smax = 0.0
  for v in s {
    if v > smax {
      smax = v
    }
  }
  let mx = if m > n { m } else { n }
  let tol = smax * mx.to_double() * 2.220446049250313e-16
  let mut rank = 0
  for v in s {
    if v > tol {
      rank = rank + 1
    }
  }
  rank
}

///|
pub struct PcoaResult {
  eigvals : Array[Double]
  proportion_explained : Array[Double]
  samples : Array[Array[Double]]
}

///|
/// Principal Coordinate Analysis (PCoA): Gower-centered distance matrix →
/// eigendecomposition → descending eigenvalues, proportion explained, and
/// sample coordinates (eigenvectors scaled by sqrt(eigenvalue)).
pub fn skbio_pcoa(dm : Array[Array[Double]]) -> PcoaResult {
  let n = dm.length()
  let b = skbio_center_distance_matrix(dm)
  let (raw_evals, raw_evecs) = eigh_symmetric(b, n)
  let (sorted_evals, sorted_evecs) = linalg_sort_eigen(raw_evals, raw_evecs, n)
  let evals = Array::make(n, 0.0)
  for j = 0; j < n; j = j + 1 {
    evals[j] = if sorted_evals[j] < 0.0 { 0.0 } else { sorted_evals[j] }
  }
  let mut sum = 0.0
  for v in evals {
    sum = sum + v
  }
  let proportion_explained = Array::make(n, 0.0)
  for j = 0; j < n; j = j + 1 {
    proportion_explained[j] = evals[j] / sum
  }
  let samples : Array[Array[Double]] = Array(capacity=0)
  for i = 0; i < n; i = i + 1 {
    let row = Array::make(n, 0.0)
    for j = 0; j < n; j = j + 1 {
      row[j] = sorted_evecs[i][j] * evals[j].sqrt()
    }
    samples.push(row)
  }
  PcoaResult::{ eigvals: evals, proportion_explained, samples }
}

///|
/// PCoA biplot: project feature matrix `y` (samples × features) onto the PCoA
/// axes. Returns the feature loadings (features × axes).
pub fn skbio_pcoa_biplot(
  samples : Array[Array[Double]],
  eigvals : Array[Double],
  y : Array[Array[Double]],
) -> Array[Array[Double]] {
  let n = samples.length()
  let d = samples[0].length()
  let p = y[0].length()
  let scaled = skbio_scale(samples, true, true, 1)
  // spc = (1/(N-1)) * y.T @ scaled  (p × d)
  let yt = skb_transpose(y)
  let spc = skb_mul(yt, scaled)
  let inv_sqrt : Array[Double] = Array(capacity=0)
  for j = 0; j < d; j = j + 1 {
    inv_sqrt.push(if eigvals[j] > 0.0 { 1.0 / eigvals[j].sqrt() } else { 0.0 })
  }
  let rn = (n - 1).to_double().sqrt()
  let out : Array[Array[Double]] = Array(capacity=0)
  for i = 0; i < p; i = i + 1 {
    let row = Array::make(d, 0.0)
    for j = 0; j < d; j = j + 1 {
      // Uproj = sqrt(N-1) * spc * diag(1/sqrt(eigval))  where spc = 1/(N-1)·Yᵀ·scale(coord)
      // => Yᵀ·scale(coord) · 1/sqrt(eigval) / sqrt(N-1)
      row[j] = spc[i][j] * inv_sqrt[j] / rn
    }
    out.push(row)
  }
  out
}
