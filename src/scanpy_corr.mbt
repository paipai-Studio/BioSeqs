/// Scanpy second-round: Correlation analysis
/// Inspired by scanpy.tl.corr2, correlation_matrix, partial_correlation

///|
/// Compute pairwise correlation between genes or between cells.
/// Returns a correlation matrix.
pub fn scanpy_correlation_matrix(
  x : Array[Array[Double]],
  method_~ : String = "pearson",
) -> Array[Array[Double]] {
  let n_rows = x.length()
  let n_cols = if n_rows > 0 { x[0].length() } else { 0 }
  // Compute correlation between columns (genes)
  let corr = Array(capacity=n_cols)
  for g1 in 0..<n_cols {
    let row = Array(capacity=n_cols)
    for g2 in 0..<n_cols {
      if g1 == g2 {
        row.push(1.0)
      } else if g1 < g2 {
        let r = compute_corr(x, g1, g2, method_)
        row.push(r)
      } else {
        row.push(0.0) // will be filled from symmetry
      }
    }
    corr.push(row)
  }
  // Fill lower triangle from symmetry
  for g1 in 0..<n_cols {
    for g2 in 0..<g1 {
      corr[g1][g2] = corr[g2][g1]
    }
  }
  corr
}

///|
/// corr2 — correlate a gene against an obs variable.
/// Returns correlation coefficient and p-value for each gene.
pub fn scanpy_corr2(
  x : Array[Array[Double]],
  obs_var : Array[Double],
) -> Array[(Double, Double)] {
  let n_cells = x.length()
  let n_genes = if n_cells > 0 { x[0].length() } else { 0 }
  let result = Array(capacity=n_genes)
  for g in 0..<n_genes {
    let gene_vals = Array(capacity=n_cells)
    for i in 0..<n_cells {
      gene_vals.push(x[i][g])
    }
    let r = pearson_r(gene_vals, obs_var)
    // Approximate p-value from t distribution
    let n_f = n_cells.to_double()
    let t = r * @math.pow((n_f - 2.0) / (1.0 - r * r + 0.0001), 0.5)
    let pval = normal_pvalue_two_sided_corr(t)
    result.push((r, pval))
  }
  result
}

///|
/// Partial correlation: correlation between two variables
/// controlling for a set of confounding variables.
/// Uses the precision matrix (inverse of correlation matrix) approach.
pub fn scanpy_partial_correlation(
  corr_matrix : Array[Array[Double]],
) -> Array[Array[Double]] {
  let n = corr_matrix.length()
  // Compute inverse of correlation matrix (precision matrix)
  // Using Gauss-Jordan elimination
  // Augment with identity
  let aug = Array(capacity=n)
  for i in 0..<n {
    let row = Array(capacity=2 * n)
    for j in 0..<n {
      row.push(corr_matrix[i][j])
    }
    for j in 0..<n {
      row.push(if i == j { 1.0 } else { 0.0 })
    }
    aug.push(row)
  }
  // Forward elimination
  for col in 0..<n {
    // Find pivot
    let mut max_row = col
    let mut max_val = 0.0
    for row in col..<n {
      let val = if aug[row][col] < 0.0 { -aug[row][col] } else { aug[row][col] }
      if val > max_val {
        max_val = val
        max_row = row
      }
    }
    if max_row != col {
      let tmp = aug[col]
      aug[col] = aug[max_row]
      aug[max_row] = tmp
    }
    let pivot = aug[col][col]
    if pivot == 0.0 {
      continue
    }
    for j in 0..<(2 * n) {
      aug[col][j] = aug[col][j] / pivot
    }
    for row in 0..<n {
      if row != col {
        let factor = aug[row][col]
        for j in 0..<(2 * n) {
          aug[row][j] = aug[row][j] - factor * aug[col][j]
        }
      }
    }
  }
  // Extract inverse (right half of augmented matrix)
  let inv = Array(capacity=n)
  for i in 0..<n {
    let row = Array(capacity=n)
    for j in 0..<n {
      row.push(aug[i][n + j])
    }
    inv.push(row)
  }
  // Partial correlation: -P_ij / sqrt(P_ii * P_jj)
  let pcorr = Array(capacity=n)
  for i in 0..<n {
    let row = Array(capacity=n)
    for j in 0..<n {
      if i == j {
        row.push(1.0)
      } else {
        let denom = @math.pow(inv[i][i] * inv[j][j] + 0.0001, 0.5)
        row.push(-inv[i][j] / denom)
      }
    }
    pcorr.push(row)
  }
  pcorr
}

// Helper: compute correlation between two columns

///|
fn compute_corr(
  x : Array[Array[Double]],
  g1 : Int32,
  g2 : Int32,
  method_ : String,
) -> Double {
  let n = x.length()
  let v1 = Array(capacity=n)
  let v2 = Array(capacity=n)
  for i in 0..<n {
    v1.push(x[i][g1])
    v2.push(x[i][g2])
  }
  match method_ {
    "pearson" => pearson_r(v1, v2)
    "spearman" => spearman_r(v1, v2)
    _ => pearson_r(v1, v2)
  }
}

// Helper: Pearson correlation

///|
fn pearson_r(a : Array[Double], b : Array[Double]) -> Double {
  let n = a.length()
  if n < 2 {
    return 0.0
  }
  let mut sa = 0.0
  let mut sb = 0.0
  let mut sab = 0.0
  let mut sa2 = 0.0
  let mut sb2 = 0.0
  for i in 0..<n {
    sa = sa + a[i]
    sb = sb + b[i]
    sab = sab + a[i] * b[i]
    sa2 = sa2 + a[i] * a[i]
    sb2 = sb2 + b[i] * b[i]
  }
  let nf = n.to_double()
  let num = sab - sa * sb / nf
  let den = @math.pow((sa2 - sa * sa / nf) * (sb2 - sb * sb / nf), 0.5)
  if den == 0.0 {
    0.0
  } else {
    num / den
  }
}

// Helper: Spearman rank correlation

///|
fn spearman_r(a : Array[Double], b : Array[Double]) -> Double {
  let n = a.length()
  if n < 2 {
    return 0.0
  }
  // Rank a and b
  let ra = rank_array(a)
  let rb = rank_array(b)
  pearson_r(ra, rb)
}

// Helper: rank an array (average ranks for ties)

///|
fn rank_array(a : Array[Double]) -> Array[Double] {
  let n = a.length()
  let indexed = Array(capacity=n)
  for i in 0..<n {
    indexed.push((i, a[i]))
  }
  indexed.sort_by(fn(x, y) {
    if x.1 < y.1 {
      -1
    } else if x.1 > y.1 {
      1
    } else {
      0
    }
  })
  let ranks = Array(capacity=n)
  for _i in 0..<n {
    ranks.push(0.0)
  }
  let mut i = 0
  while i < n {
    let mut j = i + 1
    while j < n && indexed[j].1 == indexed[i].1 {
      j = j + 1
    }
    let avg_rank = (i + j + 1).to_double() / 2.0
    for k in i..<j {
      ranks[indexed[k].0] = avg_rank
    }
    i = j
  }
  ranks
}

// Helper: two-sided p-value from normal approximation

///|
fn normal_pvalue_two_sided_corr(z : Double) -> Double {
  let abs_z = if z < 0.0 { -z } else { z }
  let t = 1.0 / (1.0 + 0.3275911 * abs_z / @math.pow(2.0, 0.5))
  let z2 = abs_z / @math.pow(2.0, 0.5)
  let erfc_val = t *
    (
      0.254829592 +
      t *
      (-0.284496736 + t * (1.421413741 + t * (-1.453152027 + t * 1.061405429)))
    ) *
    @math.pow(2.718281828459045, -z2 * z2)
  if erfc_val > 1.0 {
    1.0
  } else if erfc_val < 0.0 {
    0.0
  } else {
    erfc_val
  }
}
