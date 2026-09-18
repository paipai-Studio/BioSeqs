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
  let t : Array[Array[Double]] = Array::new()
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
  let c : Array[Array[Double]] = Array::new()
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
  let out : Array[Double] = Array::new()
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
  let e : Array[Array[Double]] = Array::new()
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
  let f : Array[Array[Double]] = Array::new()
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
  let sorted_evecs : Array[Array[Double]] = Array::new()
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
  let u_hat : Array[Array[Double]] = Array::new()
  for i = 0; i < n; i = i + 1 {
    let row = Array::make(rank, 0.0)
    for k = 0; k < rank; k = k + 1 {
      row[k] = sorted_evecs[i][k]
    }
    u_hat.push(row)
  }

  // U = Q_bar^T U_hat / W  (m × rank)
  let u : Array[Array[Double]] = Array::new()
  for j = 0; j < m; j = j + 1 {
    let row = Array::make(rank, 0.0)
    u.push(row)
  }
  for k = 0; k < rank; k = k + 1 {
    // col = Q_bar^T @ U_hat[:,k]
    let qtu : Array[Double] = Array::new()
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
  let q : Array[Array[Double]] = Array::new()
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
  let qbar : Array[Array[Double]] = Array::new()
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
  let v : Array[Array[Double]] = Array::new()
  for j = 0; j < m; j = j + 1 {
    let row = Array::make(rank, 0.0)
    let denom = col_marginals[j].sqrt()
    for k = 0; k < rank; k = k + 1 {
      row[k] = if denom > 1.0e-15 { u[j][k] / denom } else { 0.0 }
    }
    v.push(row)
  }
  let v_hat : Array[Array[Double]] = Array::new()
  for i = 0; i < n; i = i + 1 {
    let row = Array::make(rank, 0.0)
    let denom = row_marginals[i].sqrt()
    for k = 0; k < rank; k = k + 1 {
      row[k] = if denom > 1.0e-15 { u_hat[i][k] / denom } else { 0.0 }
    }
    v_hat.push(row)
  }
  // F = V_hat * W ; F_hat = V * W
  let f : Array[Array[Double]] = Array::new()
  for i = 0; i < n; i = i + 1 {
    let row = Array::make(rank, 0.0)
    for k = 0; k < rank; k = k + 1 {
      row[k] = v_hat[i][k] * w[k]
    }
    f.push(row)
  }
  let f_hat : Array[Array[Double]] = Array::new()
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
  let out : Array[Double] = Array::new()
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
  let out : Array[Array[Double]] = Array::new()
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