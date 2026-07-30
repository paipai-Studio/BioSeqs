///|
/// Missing value imputation utilities inspired by the Bioconductor {impute} package.
///
/// References:
///   - https://www.bioconductor.org/packages/release/bioc/html/impute.html
///   - Hastie T, Tibshirani R, Narasimhan B, Chu G. (2021) impute: impute: Imputation for
///     microarray data. R package version 1.66.0
///
/// Implemented imputers (all operate on a 2D numeric matrix with NA gaps):
///   * impute_na_by_mean / impute_na_by_median / impute_na_by_min / impute_na_by_zero
///   * impute_by_row_mean / impute_by_col_mean (margin-aware)
///   * impute_by_knn         - KNN with Euclidean distance and inverse-distance weighting
///   * impute_locf           - Last-Observation-Carried-Forward (time-course)
///   * impute_nocb           - Next-Observation-Carried-Back (time-course)
///
/// Helper functions:
///   * make_na_mask          - Mark NA positions for pre-processing checks
///   * impute_na_stats       - Summary of missingness per dataset
///   * apply_value_clamp     - Post-impute min/max clamping to preserve physiological range

// ---- utility: MoonBit does not expose a Double::INFINITY constant in core ----
fn impute_infinity() -> Double { 1.7976931348623157e+308 }

fn impute_is_na(x : Double) -> Bool { x.is_nan() }

// ---- struct definitions -----------------------------------------------------

pub struct ImputeNAStats {
  total_rows : Int
  total_cols : Int
  total_na : Int
  na_fraction : Double
  rows_with_na : Int
  cols_with_na : Int
  row_na_counts : Array[Int]
  col_na_counts : Array[Int]
}

pub struct KNNImputeParam {
  pub mut k : Int
  pub by_row : Bool
  pub eps : Double
  pub max_value : Option[Double]
  pub min_value : Option[Double]
}

pub fn KNNImputeParam::new() -> KNNImputeParam {
  KNNImputeParam::{ k: 10, by_row: true, eps: 0.000001, max_value: None, min_value: None, }
}

pub fn KNNImputeParam::with_k(k : Int) -> KNNImputeParam {
  let p = KNNImputeParam::new()
  if k >= 1 { p.k = k }
  p
}

// ---- generic helpers --------------------------------------------------------

fn impute_copy_matrix(mat : Array[Array[Double]]) -> Array[Array[Double]] {
  let n = mat.length()
  let out : Array[Array[Double]] = Array::new()
  for i = 0; i < n; i = i + 1 {
    let row = mat[i]
    let m = row.length()
    let new_row : Array[Double] = Array::new()
    for j = 0; j < m; j = j + 1 {
      new_row.push(row[j])
    }
    out.push(new_row)
  }
  out
}

fn impute_transpose(mat : Array[Array[Double]]) -> Array[Array[Double]] {
  if mat.length() == 0 { return Array::new() }
  let rows = mat.length()
  let cols = mat[0].length()
  let out : Array[Array[Double]] = Array::new()
  for j = 0; j < cols; j = j + 1 {
    let col : Array[Double] = Array::new()
    for i = 0; i < rows; i = i + 1 {
      col.push(mat[i][j])
    }
    out.push(col)
  }
  out
}

/// Mark which positions are missing (NA).  True => NA.
pub fn make_na_mask(mat : Array[Array[Double]]) -> Array[Array[Bool]] {
  let out : Array[Array[Bool]] = Array::new()
  for row in mat {
    let r : Array[Bool] = Array::new()
    for v in row {
      r.push(impute_is_na(v))
    }
    out.push(r)
  }
  out
}

/// Compute summary of missing data.
pub fn impute_na_stats(mat : Array[Array[Double]]) -> ImputeNAStats {
  let rows = mat.length()
  let cols = if rows == 0 { 0 } else { mat[0].length() }
  let mut total = 0
  let row_counts : Array[Int] = Array::new()
  let col_counts : Array[Int] = Array::make(cols, 0)
  let mut r_na = 0
  for i = 0; i < rows; i = i + 1 {
    let mut rc = 0
    for j = 0; j < cols; j = j + 1 {
      if impute_is_na(mat[i][j]) {
        total = total + 1
        rc = rc + 1
        col_counts[j] = col_counts[j] + 1
      }
    }
    row_counts.push(rc)
    if rc > 0 { r_na = r_na + 1 }
  }
  let mut c_na = 0
  for j = 0; j < cols; j = j + 1 {
    if col_counts[j] > 0 { c_na = c_na + 1 }
  }
  let denom = if rows * cols == 0 {
    1.0
  } else {
    (rows * cols).to_double()
  }
  ImputeNAStats::{
    total_rows: rows,
    total_cols: cols,
    total_na: total,
    na_fraction: total.to_double() / denom,
    rows_with_na: r_na,
    cols_with_na: c_na,
    row_na_counts: row_counts,
    col_na_counts: col_counts,
  }
}

pub fn impute_na_summary(stats : ImputeNAStats) -> String {
  let buf = StringBuilder::new()
  buf.write_string("ImputeNAStats [\{stats.total_rows} x \{stats.total_cols}]\n")
  let pct = stats.na_fraction * 100.0
  buf.write_string("  total NAs: \{stats.total_na} (\{pct}%)\n")
  buf.write_string("  rows with NA: \{stats.rows_with_na} / \{stats.total_rows}\n")
  buf.write_string("  cols with NA: \{stats.cols_with_na} / \{stats.total_cols}\n")
  buf.to_string()
}

fn impute_mean_clean(arr : Array[Double]) -> Double {
  let mut s = 0.0
  let mut c = 0
  for x in arr {
    if !impute_is_na(x) {
      s = s + x
      c = c + 1
    }
  }
  if c == 0 { 0.0 } else { s / c.to_double() }
}

fn impute_median_clean(arr : Array[Double]) -> Double {
  let clean : Array[Double] = Array::new()
  for x in arr {
    if !impute_is_na(x) { clean.push(x) }
  }
  clean.sort_by(fn(a : Double, b : Double) -> Int {
    if a < b { -1 } else if a > b { 1 } else { 0 }
  })
  let n = clean.length()
  if n == 0 { return 0.0 }
  let mid = n / 2
  if n % 2 == 0 {
    (clean[mid - 1] + clean[mid]) / 2.0
  } else {
    clean[mid]
  }
}

fn impute_min_clean(arr : Array[Double]) -> Double {
  let mut m = impute_infinity()
  for x in arr {
    if !impute_is_na(x) && x < m { m = x }
  }
  if impute_is_na(m) || m == impute_infinity() { 0.0 } else { m }
}

fn clamp_value(v : Double, p : KNNImputeParam) -> Double {
  let mut out = v
  match p.max_value {
    Some(hi) => if out > hi { out = hi }
    None => ()
  }
  match p.min_value {
    Some(lo) => if out < lo { out = lo }
    None => ()
  }
  out
}

// ---- simple column-wise imputers -------------------------------------------

pub fn impute_na_by_mean(mat : Array[Array[Double]]) -> Array[Array[Double]] {
  impute_by_col_mean(mat)
}

pub fn impute_na_by_median(mat : Array[Array[Double]]) -> Array[Array[Double]] {
  impute_by_col_median(mat)
}

pub fn impute_na_by_min(mat : Array[Array[Double]]) -> Array[Array[Double]] {
  impute_by_col_min(mat)
}

pub fn impute_na_by_zero(mat : Array[Array[Double]]) -> Array[Array[Double]] {
  let out = impute_copy_matrix(mat)
  for i = 0; i < out.length(); i = i + 1 {
    for j = 0; j < out[i].length(); j = j + 1 {
      if impute_is_na(out[i][j]) { out[i][j] = 0.0 }
    }
  }
  out
}

pub fn impute_by_row_mean(mat : Array[Array[Double]]) -> Array[Array[Double]] {
  let out = impute_copy_matrix(mat)
  for i = 0; i < out.length(); i = i + 1 {
    let mu = impute_mean_clean(out[i])
    for j = 0; j < out[i].length(); j = j + 1 {
      if impute_is_na(out[i][j]) { out[i][j] = mu }
    }
  }
  out
}

pub fn impute_by_col_mean(mat : Array[Array[Double]]) -> Array[Array[Double]] {
  if mat.length() == 0 { return Array::new() }
  let cols = mat[0].length()
  let means : Array[Double] = Array::new()
  for j = 0; j < cols; j = j + 1 {
    let c : Array[Double] = Array::new()
    for i = 0; i < mat.length(); i = i + 1 { c.push(mat[i][j]) }
    means.push(impute_mean_clean(c))
  }
  let out = impute_copy_matrix(mat)
  for i = 0; i < out.length(); i = i + 1 {
    for j = 0; j < cols; j = j + 1 {
      if impute_is_na(out[i][j]) { out[i][j] = means[j] }
    }
  }
  out
}

pub fn impute_by_row_median(mat : Array[Array[Double]]) -> Array[Array[Double]] {
  let out = impute_copy_matrix(mat)
  for i = 0; i < out.length(); i = i + 1 {
    let med = impute_median_clean(out[i])
    for j = 0; j < out[i].length(); j = j + 1 {
      if impute_is_na(out[i][j]) { out[i][j] = med }
    }
  }
  out
}

pub fn impute_by_col_median(mat : Array[Array[Double]]) -> Array[Array[Double]] {
  if mat.length() == 0 { return Array::new() }
  let cols = mat[0].length()
  let meds : Array[Double] = Array::new()
  for j = 0; j < cols; j = j + 1 {
    let c : Array[Double] = Array::new()
    for i = 0; i < mat.length(); i = i + 1 { c.push(mat[i][j]) }
    meds.push(impute_median_clean(c))
  }
  let out = impute_copy_matrix(mat)
  for i = 0; i < out.length(); i = i + 1 {
    for j = 0; j < cols; j = j + 1 {
      if impute_is_na(out[i][j]) { out[i][j] = meds[j] }
    }
  }
  out
}

pub fn impute_by_row_min(mat : Array[Array[Double]]) -> Array[Array[Double]] {
  let out = impute_copy_matrix(mat)
  for i = 0; i < out.length(); i = i + 1 {
    let mn = impute_min_clean(out[i])
    for j = 0; j < out[i].length(); j = j + 1 {
      if impute_is_na(out[i][j]) { out[i][j] = mn }
    }
  }
  out
}

pub fn impute_by_col_min(mat : Array[Array[Double]]) -> Array[Array[Double]] {
  if mat.length() == 0 { return Array::new() }
  let cols = mat[0].length()
  let mins : Array[Double] = Array::new()
  for j = 0; j < cols; j = j + 1 {
    let c : Array[Double] = Array::new()
    for i = 0; i < mat.length(); i = i + 1 { c.push(mat[i][j]) }
    mins.push(impute_min_clean(c))
  }
  let out = impute_copy_matrix(mat)
  for i = 0; i < out.length(); i = i + 1 {
    for j = 0; j < cols; j = j + 1 {
      if impute_is_na(out[i][j]) { out[i][j] = mins[j] }
    }
  }
  out
}

pub fn apply_value_clamp(
  mat : Array[Array[Double]],
  min_v : Double,
  max_v : Double,
) -> Array[Array[Double]] {
  let out = impute_copy_matrix(mat)
  for i = 0; i < out.length(); i = i + 1 {
    for j = 0; j < out[i].length(); j = j + 1 {
      if impute_is_na(out[i][j]) { continue }
      if out[i][j] < min_v { out[i][j] = min_v }
      if out[i][j] > max_v { out[i][j] = max_v }
    }
  }
  out
}

// ---- LOCF / NOCB ------------------------------------------------------------

pub fn impute_locf(mat : Array[Array[Double]], by_row? : Bool = true) -> Array[Array[Double]] {
  let base = if by_row {
    impute_copy_matrix(mat)
  } else {
    impute_transpose(mat)
  }
  for i = 0; i < base.length(); i = i + 1 {
    let mut last = 0.0 / 0.0 // NaN
    for j = 0; j < base[i].length(); j = j + 1 {
      if impute_is_na(base[i][j]) {
        if !impute_is_na(last) { base[i][j] = last }
      } else {
        last = base[i][j]
      }
    }
  }
  if by_row { base } else { impute_transpose(base) }
}

pub fn impute_nocb(mat : Array[Array[Double]], by_row? : Bool = true) -> Array[Array[Double]] {
  let base = if by_row {
    impute_copy_matrix(mat)
  } else {
    impute_transpose(mat)
  }
  for i = 0; i < base.length(); i = i + 1 {
    let row = base[i]
    let n = row.length()
    let mut nxt = 0.0 / 0.0
    for j0 = 0; j0 < n; j0 = j0 + 1 {
      let j = n - 1 - j0
      if impute_is_na(row[j]) {
        if !impute_is_na(nxt) { row[j] = nxt }
      } else {
        nxt = row[j]
      }
    }
  }
  if by_row { base } else { impute_transpose(base) }
}

// ---- KNN imputation ---------------------------------------------------------

// Euclidean distance on coordinates where BOTH vectors have non-NA values.
// Returns infinity when no valid overlap exists.
fn impute_pairwise_dist_ignore_na(
  a : Array[Double],
  b : Array[Double],
  eps : Double,
) -> Double {
  let n = a.length()
  let mut s = 0.0
  let mut k = 0
  for j = 0; j < n; j = j + 1 {
    if !impute_is_na(a[j]) && !impute_is_na(b[j]) {
      let d = a[j] - b[j]
      s = s + d * d
      k = k + 1
    }
  }
  if k == 0 {
    1.0 / 0.0
  } else {
    (s / k.to_double() + eps).sqrt()
  }
}

fn impute_compute_all_pair_distances(
  m : Array[Array[Double]],
  eps : Double,
) -> Array[Array[Double]] {
  let n = m.length()
  let out : Array[Array[Double]] = Array::new()
  for i = 0; i < n; i = i + 1 {
    let row : Array[Double] = Array::new()
    for j = 0; j < n; j = j + 1 {
      let d = if i == j {
        0.0
      } else {
        impute_pairwise_dist_ignore_na(m[i], m[j], eps)
      }
      row.push(d)
    }
    out.push(row)
  }
  out
}

fn impute_choose_k_neighbors(
  dist_row : Array[Double],
  k : Int,
  self_i : Int,
) -> Array[(Int, Double)] {
  let pairs : Array[(Int, Double)] = Array::new()
  let n = dist_row.length()
  for i = 0; i < n; i = i + 1 {
    if i == self_i { continue }
    pairs.push((i, dist_row[i]))
  }
  pairs.sort_by(fn(p : (Int, Double), q : (Int, Double)) -> Int {
    if p.1 < q.1 { -1 } else if p.1 > q.1 { 1 } else { 0 }
  })
  let m = if k < pairs.length() { k } else { pairs.length() }
  let out : Array[(Int, Double)] = Array::new()
  for i = 0; i < m; i = i + 1 {
    out.push(pairs[i])
  }
  out
}

fn impute_weighted_average(
  neighbors : Array[(Int, Double)],
  m : Array[Array[Double]],
  j : Int,
  eps : Double,
) -> Double {
  let mut num = 0.0
  let mut den = 0.0
  let n = neighbors.length()
  for idx = 0; idx < n; idx = idx + 1 {
    let ni = neighbors[idx].0
    let d = neighbors[idx].1
    let v = m[ni][j]
    if impute_is_na(v) { continue }
    let w = 1.0 / (d + eps)
    num = num + w * v
    den = den + w
  }
  if den == 0.0 { 0.0 } else { num / den }
}

pub fn impute_by_knn(
  mat : Array[Array[Double]],
  param : KNNImputeParam,
) -> Array[Array[Double]] {
  let (m, axis_len) = if param.by_row {
    (impute_copy_matrix(mat), mat.length())
  } else {
    (impute_transpose(mat), mat[0].length())
  }
  let dist = impute_compute_all_pair_distances(m, param.eps)
  for i = 0; i < axis_len; i = i + 1 {
    let neighbors = impute_choose_k_neighbors(dist[i], param.k, i)
    for j = 0; j < m[0].length(); j = j + 1 {
      if impute_is_na(m[i][j]) {
        let raw = impute_weighted_average(neighbors, m, j, param.eps)
        m[i][j] = clamp_value(raw, param)
      }
    }
  }
  if param.by_row { m } else { impute_transpose(m) }
}

pub fn impute_by_knn_simple(mat : Array[Array[Double]], k : Int) -> Array[Array[Double]] {
  let p = KNNImputeParam::new()
  p.k = k
  impute_by_knn(mat, p)
}
