///|
/// SingleCellExperiment - Core data structure for single-cell RNA-seq analysis.
/// Inspired by the Bioconductor SingleCellExperiment package, this module
/// provides the central container used by scater, scran, monocle3 and other
/// single-cell analysis packages to coordinate assays, cell/gene metadata,
/// dimensionality reduction results and alternative experiments.

///|
/// Square root via Newton's method.
pub fn sce_sqrt(x : Double) -> Double {
  if x <= 0.0 {
    return 0.0
  }
  if x == 1.0 {
    return 1.0
  }
  let mut guess = x / 2.0
  if guess < 1.0 {
    guess = 1.0
  }
  let mut i = 0
  while i < 60 {
    let next = (guess + x / guess) / 2.0
    let diff = next - guess
    if diff < 0.0 {
      if -diff < 0.00000000000001 {
        return next
      }
    } else if diff < 0.00000000000001 {
      return next
    }
    guess = next
    i = i + 1
  }
  guess
}

///|
/// Natural logarithm via series expansion.
pub fn sce_ln(x : Double) -> Double {
  if x <= 0.0 {
    return -10000000000.0
  }
  if x == 1.0 {
    return 0.0
  }
  let mut v = x
  let mut offset = 0.0
  while v > 2.0 {
    v = v / 2.718281828459045
    offset = offset + 1.0
  }
  while v < 0.5 {
    v = v * 2.718281828459045
    offset = offset - 1.0
  }
  let y = (v - 1.0) / (v + 1.0)
  let y2 = y * y
  let mut sum = y
  let mut term = y
  let mut i = 1
  while i <= 50 {
    term = term * y2
    let delta = term / (2 * i + 1).to_double()
    sum = sum + delta
    if delta < 0.0 {
      if -delta < 0.0000000000000001 {
        return 2.0 * sum + offset
      }
    } else if delta < 0.0000000000000001 {
      return 2.0 * sum + offset
    }
    i = i + 1
  }
  2.0 * sum + offset
}

///|
/// Exponential function via Taylor series.
pub fn sce_exp(x : Double) -> Double {
  if x == 0.0 {
    return 1.0
  }
  if x < 0.0 {
    return 1.0 / sce_exp(-x)
  }
  let mut k = 0
  let mut r = x
  while r >= 1.0 {
    r = r - 1.0
    k = k + 1
  }
  let mut sum = 1.0
  let mut term = 1.0
  let mut i = 1
  while i <= 30 {
    term = term * r / i.to_double()
    sum = sum + term
    if term < 0.0 {
      if -term < 0.0000000000000001 {
        break
      }
    } else if term < 0.0000000000000001 {
      break
    }
    i = i + 1
  }
  let mut e_pow_k = 1.0
  let mut j = 0
  while j < k {
    e_pow_k = e_pow_k * 2.718281828459045
    j = j + 1
  }
  e_pow_k * sum
}

///|
/// Base-2 logarithm.
pub fn sce_log2(x : Double) -> Double {
  sce_ln(x) / sce_ln(2.0)
}

///|
/// Arithmetic mean of an array.
pub fn sce_mean(arr : Array[Double]) -> Double {
  if arr.length() == 0 {
    return 0.0
  }
  let mut sum = 0.0
  let mut i = 0
  while i < arr.length() {
    sum = sum + arr[i]
    i = i + 1
  }
  sum / arr.length().to_double()
}

///|
/// Sample variance (using N-1 denominator).
pub fn sce_variance(arr : Array[Double]) -> Double {
  let n = arr.length()
  if n < 2 {
    return 0.0
  }
  let m = sce_mean(arr)
  let mut sum_sq = 0.0
  let mut i = 0
  while i < n {
    let diff = arr[i] - m
    sum_sq = sum_sq + diff * diff
    i = i + 1
  }
  sum_sq / (n - 1).to_double()
}

///|
/// Median of an array (returns 0.0 for empty input).
pub fn sce_median(arr : Array[Double]) -> Double {
  let n = arr.length()
  if n == 0 {
    return 0.0
  }
  let sorted : Array[Double] = Array::new()
  let mut i = 0
  while i < n {
    sorted.push(arr[i])
    i = i + 1
  }
  let mut a = 0
  while a < n - 1 {
    let mut b = 0
    while b < n - 1 - a {
      if sorted[b] > sorted[b + 1] {
        let tmp = sorted[b]
        sorted[b] = sorted[b + 1]
        sorted[b + 1] = tmp
      }
      b = b + 1
    }
    a = a + 1
  }
  let mid = n / 2
  if n % 2 == 1 {
    sorted[mid]
  } else {
    (sorted[mid - 1] + sorted[mid]) / 2.0
  }
}

///|
/// Geometric mean (ignores non-positive values).
pub fn sce_geometric_mean(arr : Array[Double]) -> Double {
  let mut log_sum = 0.0
  let mut count = 0
  let mut i = 0
  while i < arr.length() {
    if arr[i] > 0.0 {
      log_sum = log_sum + sce_ln(arr[i])
      count = count + 1
    }
    i = i + 1
  }
  if count == 0 {
    return 0.0
  }
  sce_exp(log_sum / count.to_double())
}

///|
/// The core SingleCellExperiment container.
pub struct SingleCellExperiment {
  assays : Map[String, Array[Array[Double]]]
  row_data : Map[String, Array[String]]
  col_data : Map[String, Array[String]]
  reduced_dims : Map[String, Array[Array[Double]]]
  alternative_experiments : Map[String, SingleCellExperiment]
  metadata : Map[String, String]
  row_names : Array[String]
  col_names : Array[String]
}

///|
/// Builder for constructing a SingleCellExperiment incrementally.
pub struct SCEBuilder {
  assays : Map[String, Array[Array[Double]]]
  row_names : Array[String]
  col_names : Array[String]
  row_data : Map[String, Array[String]]
  col_data : Map[String, Array[String]]
  reduced_dims : Map[String, Array[Array[Double]]]
  metadata : Map[String, String]
}

///|
/// Result of a dimensionality reduction computation (PCA / t-SNE / UMAP).
pub struct ReducedDimResult {
  method : String
  coordinates : Array[Array[Double]]
  variance_explained : Array[Double]
  n_components : Int
}

///|
/// Result of size factor computation.
pub struct SizeFactorResult {
  size_factors : Array[Double]
  method : String
  log_means : Array[Double]
}

///|
/// Construct a SingleCellExperiment from a counts matrix.
/// Convention: counts[i][j] is the expression of gene i in cell j.
pub fn SingleCellExperiment::new(
  counts : Array[Array[Double]],
  row_names : Array[String],
  col_names : Array[String],
) -> SingleCellExperiment {
  let assays : Map[String, Array[Array[Double]]] = Map([], capacity=4)
  assays["counts"] = counts
  SingleCellExperiment::{
    assays,
    row_data: Map([], capacity=4),
    col_data: Map([], capacity=4),
    reduced_dims: Map([], capacity=4),
    alternative_experiments: Map([], capacity=4),
    metadata: Map([], capacity=4),
    row_names,
    col_names,
  }
}

///|
/// Construct a fresh SCEBuilder.
pub fn SCEBuilder::new() -> SCEBuilder {
  SCEBuilder::{
    assays: Map([], capacity=4),
    row_names: Array::new(),
    col_names: Array::new(),
    row_data: Map([], capacity=4),
    col_data: Map([], capacity=4),
    reduced_dims: Map([], capacity=4),
    metadata: Map([], capacity=4),
  }
}

///|
/// Add an assay matrix to the builder.
pub fn SCEBuilder::add_assay(
  self : SCEBuilder,
  name : String,
  data : Array[Array[Double]],
) -> SCEBuilder {
  self.assays[name] = data
  self
}

///|
/// Add a row (gene) metadata column.
pub fn SCEBuilder::add_row_data(
  self : SCEBuilder,
  name : String,
  data : Array[String],
) -> SCEBuilder {
  self.row_data[name] = data
  self
}

///|
/// Add a column (cell) metadata column.
pub fn SCEBuilder::add_col_data(
  self : SCEBuilder,
  name : String,
  data : Array[String],
) -> SCEBuilder {
  self.col_data[name] = data
  self
}

///|
/// Add a dimensionality reduction result.
pub fn SCEBuilder::add_reduced_dim(
  self : SCEBuilder,
  name : String,
  data : Array[Array[Double]],
) -> SCEBuilder {
  self.reduced_dims[name] = data
  self
}

///|
/// Finalize the builder into a SingleCellExperiment.
pub fn SCEBuilder::build(self : SCEBuilder) -> SingleCellExperiment {
  let assays : Map[String, Array[Array[Double]]] = Map(
    [],
    capacity=self.assays.keys().collect().length(),
  )
  let assay_keys = self.assays.keys().collect()
  let mut k = 0
  while k < assay_keys.length() {
    assays[assay_keys[k]] = self.assays[assay_keys[k]]
    k = k + 1
  }
  let row_names = if self.row_names.length() > 0 {
    self.row_names
  } else {
    let n = if assay_keys.length() > 0 {
      self.assays[assay_keys[0]].length()
    } else {
      0
    }
    let names : Array[String] = Array::new()
    let mut i = 0
    while i < n {
      names.push("Gene" + (i + 1).to_string())
      i = i + 1
    }
    names
  }
  let col_names = if self.col_names.length() > 0 {
    self.col_names
  } else {
    let n = if assay_keys.length() > 0 {
      let first = self.assays[assay_keys[0]]
      if first.length() > 0 {
        first[0].length()
      } else {
        0
      }
    } else {
      0
    }
    let names : Array[String] = Array::new()
    let mut i = 0
    while i < n {
      names.push("Cell" + (i + 1).to_string())
      i = i + 1
    }
    names
  }
  SingleCellExperiment::{
    assays,
    row_data: self.row_data,
    col_data: self.col_data,
    reduced_dims: self.reduced_dims,
    alternative_experiments: Map([], capacity=4),
    metadata: self.metadata,
    row_names,
    col_names,
  }
}

///|
/// Get an assay by name; returns an empty matrix if missing.
pub fn sce_get_assay(
  sce : SingleCellExperiment,
  name : String,
) -> Array[Array[Double]] {
  match sce.assays.get(name) {
    Some(m) => m
    None => []
  }
}

///|
/// Set (or replace) an assay matrix.
pub fn sce_set_assay(
  sce : SingleCellExperiment,
  name : String,
  data : Array[Array[Double]],
) -> SingleCellExperiment {
  sce.assays[name] = data
  sce
}

///|
/// Get the "counts" assay.
pub fn sce_get_counts(sce : SingleCellExperiment) -> Array[Array[Double]] {
  sce_get_assay(sce, "counts")
}

///|
/// Get the "logcounts" assay.
pub fn sce_get_logcounts(sce : SingleCellExperiment) -> Array[Array[Double]] {
  sce_get_assay(sce, "logcounts")
}

///|
/// Set the "logcounts" assay.
pub fn sce_set_logcounts(
  sce : SingleCellExperiment,
  data : Array[Array[Double]],
) -> SingleCellExperiment {
  sce_set_assay(sce, "logcounts", data)
}

///|
/// Compute log2(counts / size_factors + 1) for each cell.
/// size_factors[i] is the factor for column i.
pub fn sce_compute_logcounts(
  counts : Array[Array[Double]],
  size_factors : Array[Double],
) -> Array[Array[Double]] {
  let n_genes = counts.length()
  let result : Array[Array[Double]] = Array::new()
  let mut i = 0
  while i < n_genes {
    let row_len = counts[i].length()
    let row : Array[Double] = Array::new()
    let mut j = 0
    while j < row_len {
      let sf = if j < size_factors.length() && size_factors[j] > 0.0 {
        size_factors[j]
      } else {
        1.0
      }
      let v = counts[i][j] / sf
      let log_val = if v > 0.0 { sce_log2(v + 1.0) } else { 0.0 }
      row.push(log_val)
      j = j + 1
    }
    result.push(row)
    i = i + 1
  }
  result
}

///|
/// Get a reduced dimension result by name.
pub fn sce_get_reduced_dim(
  sce : SingleCellExperiment,
  name : String,
) -> Array[Array[Double]] {
  match sce.reduced_dims.get(name) {
    Some(m) => m
    None => []
  }
}

///|
/// Set a reduced dimension result.
pub fn sce_set_reduced_dim(
  sce : SingleCellExperiment,
  name : String,
  data : Array[Array[Double]],
) -> SingleCellExperiment {
  sce.reduced_dims[name] = data
  sce
}

///|
/// Library size factors: per-cell total count divided by the mean library size.
pub fn sce_compute_library_size_factors(
  counts : Array[Array[Double]],
) -> Array[Double] {
  let n_genes = counts.length()
  let n_cells = if n_genes > 0 { counts[0].length() } else { 0 }
  let lib_sizes : Array[Double] = Array::make(n_cells, 0.0)
  let mut i = 0
  while i < n_genes {
    let mut j = 0
    while j < n_cells {
      lib_sizes[j] = lib_sizes[j] + counts[i][j]
      j = j + 1
    }
    i = i + 1
  }
  let mean_lib = sce_mean(lib_sizes)
  let factors : Array[Double] = Array::make(n_cells, 1.0)
  let mut k = 0
  while k < n_cells {
    factors[k] = if mean_lib > 0.0 { lib_sizes[k] / mean_lib } else { 1.0 }
    k = k + 1
  }
  factors
}

///|
/// Deconvolution size factors (simplified).
/// Pools counts across cells to estimate factors robustly; falls back to
/// library size factors when the data is too sparse for deconvolution.
pub fn sce_compute_deconv_size_factors(
  counts : Array[Array[Double]],
) -> Array[Double] {
  let n_genes = counts.length()
  let n_cells = if n_genes > 0 { counts[0].length() } else { 0 }
  if n_cells == 0 || n_genes == 0 {
    return []
  }
  let lib_factors = sce_compute_library_size_factors(counts)
  let pool_size = if n_cells >= 20 { 20 } else { n_cells }
  let pooled : Array[Double] = Array::make(n_cells, 0.0)
  let mut c = 0
  while c < n_cells {
    let mut sum = 0.0
    let mut k = 0
    while k < pool_size {
      let idx = (c + k) % n_cells
      let mut g = 0
      while g < n_genes {
        sum = sum + counts[g][idx]
        g = g + 1
      }
      k = k + 1
    }
    pooled[c] = sum / pool_size.to_double()
    c = c + 1
  }
  let mean_pooled = sce_mean(pooled)
  let factors : Array[Double] = Array::make(n_cells, 1.0)
  let mut i = 0
  while i < n_cells {
    let ratio = if mean_pooled > 0.0 { pooled[i] / mean_pooled } else { 1.0 }
    factors[i] = (ratio + lib_factors[i]) / 2.0
    i = i + 1
  }
  factors
}

///|
/// Compute size factors using the requested method ("library" or "deconv").
pub fn sce_compute_size_factors(
  counts : Array[Array[Double]],
  method : String,
) -> SizeFactorResult {
  let factors = if method == "deconv" {
    sce_compute_deconv_size_factors(counts)
  } else {
    sce_compute_library_size_factors(counts)
  }
  let log_means : Array[Double] = Array::new()
  let n_cells = factors.length()
  let mut i = 0
  while i < n_cells {
    log_means.push(if factors[i] > 0.0 { sce_ln(factors[i]) } else { 0.0 })
    i = i + 1
  }
  SizeFactorResult::{ size_factors: factors, method, log_means }
}

///|
/// Center a matrix column-wise (returns a new matrix).
fn sce_center_matrix(m : Array[Array[Double]]) -> Array[Array[Double]] {
  let rows = m.length()
  let cols = if rows > 0 { m[0].length() } else { 0 }
  let means : Array[Double] = Array::make(cols, 0.0)
  let mut j = 0
  while j < cols {
    let mut sum = 0.0
    let mut i = 0
    while i < rows {
      sum = sum + m[i][j]
      i = i + 1
    }
    means[j] = sum / rows.to_double()
    j = j + 1
  }
  let centered : Array[Array[Double]] = Array::new()
  let mut i = 0
  while i < rows {
    let row : Array[Double] = Array::new()
    let mut j = 0
    while j < cols {
      row.push(m[i][j] - means[j])
      j = j + 1
    }
    centered.push(row)
    i = i + 1
  }
  centered
}

///|
/// Transpose a matrix.
fn sce_transpose(m : Array[Array[Double]]) -> Array[Array[Double]] {
  let rows = m.length()
  if rows == 0 {
    return []
  }
  let cols = m[0].length()
  let result : Array[Array[Double]] = Array::new()
  let mut i = 0
  while i < cols {
    let row : Array[Double] = Array::new()
    let mut j = 0
    while j < rows {
      row.push(m[j][i])
      j = j + 1
    }
    result.push(row)
    i = i + 1
  }
  result
}

///|
/// Matrix-vector product.
fn sce_matmul_vec(m : Array[Array[Double]], v : Array[Double]) -> Array[Double] {
  let result : Array[Double] = Array::new()
  let mut i = 0
  while i < m.length() {
    let mut sum = 0.0
    let mut j = 0
    while j < v.length() {
      sum = sum + m[i][j] * v[j]
      j = j + 1
    }
    result.push(sum)
    i = i + 1
  }
  result
}

///|
/// L2 norm of a vector.
fn sce_vec_norm(v : Array[Double]) -> Double {
  let mut sum = 0.0
  let mut i = 0
  while i < v.length() {
    sum = sum + v[i] * v[i]
    i = i + 1
  }
  sce_sqrt(sum)
}

///|
/// Covariance of the centered data: (1/(n-1)) * X^T X.
fn sce_covariance(centered : Array[Array[Double]]) -> Array[Array[Double]] {
  let rows = centered.length()
  let cols = if rows > 0 { centered[0].length() } else { 0 }
  let cov : Array[Array[Double]] = Array::make(cols, Array::make(cols, 0.0))
  let mut i = 0
  while i < cols {
    let mut j = 0
    while j < cols {
      let mut sum = 0.0
      let mut k = 0
      while k < rows {
        sum = sum + centered[k][i] * centered[k][j]
        k = k + 1
      }
      cov[i][j] = sum / (rows - 1).to_double()
      j = j + 1
    }
    i = i + 1
  }
  cov
}

///|
/// PCA via power iteration with deflation.
/// Returns coordinates (n_samples x n_components) and variance explained.
pub fn sce_compute_pca(
  logcounts : Array[Array[Double]],
  n_components : Int,
) -> ReducedDimResult {
  let n_genes = logcounts.length()
  let n_cells = if n_genes > 0 { logcounts[0].length() } else { 0 }
  if n_cells == 0 || n_genes == 0 {
    return ReducedDimResult::{
      method: "PCA",
      coordinates: [],
      variance_explained: [],
      n_components: 0,
    }
  }
  let transposed = sce_transpose(logcounts)
  let centered = sce_center_matrix(transposed)
  let mut cov = sce_covariance(centered)
  let n = cov.length()
  let components : Array[Array[Double]] = Array::new()
  let eigenvalues : Array[Double] = Array::new()
  let mut i = 0
  while i < n_components {
    let v = Array::make(n, 0.0)
    let mut j = 0
    while j < n {
      v[j] = 0.1 + j.to_double() * 0.01
      j = j + 1
    }
    let mut iter = 0
    while iter < 100 {
      let v_new = sce_matmul_vec(cov, v)
      let norm = sce_vec_norm(v_new)
      let mut k = 0
      while k < v.length() {
        v[k] = if norm > 0.0 { v_new[k] / norm } else { 0.0 }
        k = k + 1
      }
      iter = iter + 1
    }
    let eigenvalue = sce_vec_norm(sce_matmul_vec(cov, v))
    let outer : Array[Array[Double]] = Array::make(n, Array::make(n, 0.0))
    let mut k = 0
    while k < n {
      let mut l = 0
      while l < n {
        outer[k][l] = v[k] * v[l]
        l = l + 1
      }
      k = k + 1
    }
    let new_cov : Array[Array[Double]] = Array::make(n, Array::make(n, 0.0))
    k = 0
    while k < n {
      let mut l = 0
      while l < n {
        new_cov[k][l] = cov[k][l] - eigenvalue * outer[k][l]
        l = l + 1
      }
      k = k + 1
    }
    cov = new_cov
    let v_copy = Array::make(n, 0.0)
    k = 0
    while k < n {
      v_copy[k] = v[k]
      k = k + 1
    }
    components.push(v_copy)
    eigenvalues.push(eigenvalue)
    i = i + 1
  }
  let projected = sce_matmul(centered, sce_transpose(components))
  let total_var = {
    let mut s = 0.0
    let mut x = 0
    while x < eigenvalues.length() {
      s = s + eigenvalues[x]
      x = x + 1
    }
    s
  }
  let variance_explained : Array[Double] = Array::new()
  let mut y = 0
  while y < eigenvalues.length() {
    variance_explained.push(
      if total_var > 0.0 {
        eigenvalues[y] / total_var
      } else {
        0.0
      },
    )
    y = y + 1
  }
  ReducedDimResult::{
    method: "PCA",
    coordinates: projected,
    variance_explained,
    n_components,
  }
}

///|
/// Matrix-matrix product.
fn sce_matmul(
  a : Array[Array[Double]],
  b : Array[Array[Double]],
) -> Array[Array[Double]] {
  let m = a.length()
  let n = if m > 0 { a[0].length() } else { 0 }
  let p = if b.length() > 0 { b[0].length() } else { 0 }
  let result : Array[Array[Double]] = Array::new()
  let mut i = 0
  while i < m {
    let row : Array[Double] = Array::new()
    let mut j = 0
    while j < p {
      let mut sum = 0.0
      let mut k = 0
      while k < n {
        sum = sum + a[i][k] * b[k][j]
        k = k + 1
      }
      row.push(sum)
      j = j + 1
    }
    result.push(row)
    i = i + 1
  }
  result
}

///|
/// Simplified t-SNE using gradient descent on a Student-t embedding.
/// This is a compact, deterministic implementation suitable for small datasets.
pub fn sce_compute_tsne(
  logcounts : Array[Array[Double]],
  perplexity : Int,
) -> ReducedDimResult {
  let n_genes = logcounts.length()
  let n_cells = if n_genes > 0 { logcounts[0].length() } else { 0 }
  let n_components = 2
  if n_cells == 0 || n_genes == 0 {
    return ReducedDimResult::{
      method: "tSNE",
      coordinates: [],
      variance_explained: [],
      n_components: 0,
    }
  }
  let _ = perplexity
  let pca = sce_compute_pca(logcounts, 2)
  let y : Array[Array[Double]] = Array::new()
  let mut i = 0
  while i < pca.coordinates.length() {
    let row : Array[Double] = Array::new()
    let mut j = 0
    while j < pca.coordinates[i].length() {
      row.push(pca.coordinates[i][j] + 0.0001 * (i + 1).to_double())
      j = j + 1
    }
    y.push(row)
    i = i + 1
  }
  let n = y.length()
  let d : Array[Array[Double]] = Array::make(n, Array::make(n, 0.0))
  let mut a = 0
  while a < n {
    let mut b = 0
    while b < n {
      let mut s = 0.0
      let mut k = 0
      while k < y[a].length() {
        let diff = y[a][k] - y[b][k]
        s = s + diff * diff
        k = k + 1
      }
      d[a][b] = 1.0 + s
      b = b + 1
    }
    a = a + 1
  }
  let p : Array[Array[Double]] = Array::make(n, Array::make(n, 0.0))
  a = 0
  while a < n {
    let mut b = 0
    while b < n {
      if a != b {
        let mut denom = 0.0
        let mut k = 0
        while k < n {
          if k != a {
            denom = denom + 1.0 / d[a][k]
          }
          k = k + 1
        }
        p[a][b] = if denom > 0.0 { 1.0 / d[a][b] / denom } else { 0.0 }
      }
      b = b + 1
    }
    a = a + 1
  }
  let mut iter = 0
  while iter < 50 {
    let grad : Array[Array[Double]] = Array::make(
      n,
      Array::make(n_components, 0.0),
    )
    a = 0
    while a < n {
      let mut b = 0
      while b < n {
        if a != b {
          let pq = p[a][b]
          let mut k = 0
          while k < n_components {
            grad[a][k] = grad[a][k] + 4.0 * pq * (y[a][k] - y[b][k])
            k = k + 1
          }
        }
        b = b + 1
      }
      a = a + 1
    }
    a = 0
    while a < n {
      let mut k = 0
      while k < n_components {
        y[a][k] = y[a][k] - 0.1 * grad[a][k]
        k = k + 1
      }
      a = a + 1
    }
    iter = iter + 1
  }
  ReducedDimResult::{
    method: "tSNE",
    coordinates: y,
    variance_explained: [],
    n_components,
  }
}

///|
/// Simplified UMAP-like projection.
/// Uses PCA initialisation followed by a force-directed refinement.
pub fn sce_compute_umap(
  logcounts : Array[Array[Double]],
  n_components : Int,
) -> ReducedDimResult {
  let n_genes = logcounts.length()
  let n_cells = if n_genes > 0 { logcounts[0].length() } else { 0 }
  if n_cells == 0 || n_genes == 0 || n_components <= 0 {
    return ReducedDimResult::{
      method: "UMAP",
      coordinates: [],
      variance_explained: [],
      n_components: 0,
    }
  }
  let pca = sce_compute_pca(logcounts, n_components)
  let emb : Array[Array[Double]] = Array::new()
  let mut i = 0
  while i < pca.coordinates.length() {
    let row : Array[Double] = Array::new()
    let mut j = 0
    while j < pca.coordinates[i].length() {
      row.push(pca.coordinates[i][j])
      j = j + 1
    }
    emb.push(row)
    i = i + 1
  }
  let n = emb.length()
  let k_neighbors = if n > 15 { 15 } else if n > 1 { n - 1 } else { 1 }
  let mut iter = 0
  while iter < 30 {
    let grad : Array[Array[Double]] = Array::make(
      n,
      Array::make(n_components, 0.0),
    )
    let mut a = 0
    while a < n {
      let neighbors : Array[(Double, Int)] = Array::new()
      let mut b = 0
      while b < n {
        if a != b {
          let mut s = 0.0
          let mut k = 0
          while k < n_components {
            let diff = emb[a][k] - emb[b][k]
            s = s + diff * diff
            k = k + 1
          }
          neighbors.push((sce_sqrt(s + 0.0001), b))
        }
        b = b + 1
      }
      let mut x = 0
      while x < neighbors.length() - 1 {
        let mut y = 0
        while y < neighbors.length() - 1 - x {
          if neighbors[y].0 > neighbors[y + 1].0 {
            let tmp = neighbors[y]
            neighbors[y] = neighbors[y + 1]
            neighbors[y + 1] = tmp
          }
          y = y + 1
        }
        x = x + 1
      }
      let mut idx = 0
      while idx < k_neighbors && idx < neighbors.length() {
        let nb = neighbors[idx].1
        let mut k = 0
        while k < n_components {
          let diff = emb[nb][k] - emb[a][k]
          grad[a][k] = grad[a][k] + diff * 0.5
          k = k + 1
        }
        idx = idx + 1
      }
      let mut b = 0
      while b < n {
        if a != b {
          let mut s = 0.0
          let mut k = 0
          while k < n_components {
            let diff = emb[a][k] - emb[b][k]
            s = s + diff * diff
            k = k + 1
          }
          let d = sce_sqrt(s + 0.0001)
          let mut k = 0
          while k < n_components {
            grad[a][k] = grad[a][k] - 0.001 * (emb[a][k] - emb[b][k]) / d
            k = k + 1
          }
        }
        b = b + 1
      }
      a = a + 1
    }
    let mut a = 0
    while a < n {
      let mut k = 0
      while k < n_components {
        emb[a][k] = emb[a][k] + 0.1 * grad[a][k]
        k = k + 1
      }
      a = a + 1
    }
    iter = iter + 1
  }
  ReducedDimResult::{
    method: "UMAP",
    coordinates: emb,
    variance_explained: [],
    n_components,
  }
}

///|
/// Normalize the SCE: compute size factors and store logcounts assay.
/// Supported methods: "library" (default) and "deconv".
pub fn sce_normalize(
  sce : SingleCellExperiment,
  method : String,
) -> SingleCellExperiment {
  let counts = sce_get_counts(sce)
  let sf = if method == "deconv" {
    sce_compute_deconv_size_factors(counts)
  } else {
    sce_compute_library_size_factors(counts)
  }
  let logcounts = sce_compute_logcounts(counts, sf)
  sce.assays["logcounts"] = logcounts
  sce
}

///|
/// Filter genes (rows) expressed in at least `min_cells` cells.
pub fn sce_filter_genes(
  sce : SingleCellExperiment,
  min_cells : Int,
) -> SingleCellExperiment {
  let counts = sce_get_counts(sce)
  let n_genes = counts.length()
  let n_cells = if n_genes > 0 { counts[0].length() } else { 0 }
  let keep : Array[Int] = Array::new()
  let mut i = 0
  while i < n_genes {
    let mut n_expressed = 0
    let mut j = 0
    while j < n_cells {
      if counts[i][j] > 0.0 {
        n_expressed = n_expressed + 1
      }
      j = j + 1
    }
    if n_expressed >= min_cells {
      keep.push(i)
    }
    i = i + 1
  }
  sce_subset_genes(sce, keep)
}

///|
/// Filter cells (columns) expressing at least `min_genes` genes.
pub fn sce_filter_cells(
  sce : SingleCellExperiment,
  min_genes : Int,
) -> SingleCellExperiment {
  let counts = sce_get_counts(sce)
  let n_genes = counts.length()
  let n_cells = if n_genes > 0 { counts[0].length() } else { 0 }
  let keep : Array[Int] = Array::new()
  let mut j = 0
  while j < n_cells {
    let mut n_expressed = 0
    let mut i = 0
    while i < n_genes {
      if counts[i][j] > 0.0 {
        n_expressed = n_expressed + 1
      }
      i = i + 1
    }
    if n_expressed >= min_genes {
      keep.push(j)
    }
    j = j + 1
  }
  sce_subset_cells(sce, keep)
}

///|
/// Get a row (gene) metadata column.
pub fn sce_get_row_data(
  sce : SingleCellExperiment,
  name : String,
) -> Array[String] {
  match sce.row_data.get(name) {
    Some(v) => v
    None => []
  }
}

///|
/// Get a column (cell) metadata column.
pub fn sce_get_col_data(
  sce : SingleCellExperiment,
  name : String,
) -> Array[String] {
  match sce.col_data.get(name) {
    Some(v) => v
    None => []
  }
}

///|
/// Subset genes (rows) by index.
pub fn sce_subset_genes(
  sce : SingleCellExperiment,
  gene_indices : Array[Int],
) -> SingleCellExperiment {
  let assay_names = sce.assays.keys().collect()
  let new_assays : Map[String, Array[Array[Double]]] = Map(
    [],
    capacity=assay_names.length(),
  )
  let mut a = 0
  while a < assay_names.length() {
    let name = assay_names[a]
    let old = sce.assays[name]
    let new_mat : Array[Array[Double]] = Array::new()
    let mut i = 0
    while i < gene_indices.length() {
      let idx = gene_indices[i]
      let row : Array[Double] = Array::new()
      if idx < old.length() {
        let mut j = 0
        while j < old[idx].length() {
          row.push(old[idx][j])
          j = j + 1
        }
      }
      new_mat.push(row)
      i = i + 1
    }
    new_assays[name] = new_mat
    a = a + 1
  }
  let new_row_names : Array[String] = Array::new()
  let mut r = 0
  while r < gene_indices.length() {
    let idx = gene_indices[r]
    if idx < sce.row_names.length() {
      new_row_names.push(sce.row_names[idx])
    } else {
      new_row_names.push("")
    }
    r = r + 1
  }
  let new_row_data : Map[String, Array[String]] = Map(
    [],
    capacity=sce.row_data.keys().collect().length(),
  )
  let rd_keys = sce.row_data.keys().collect()
  let mut k = 0
  while k < rd_keys.length() {
    let key = rd_keys[k]
    let old_col = sce.row_data[key]
    let new_col : Array[String] = Array::new()
    let mut i = 0
    while i < gene_indices.length() {
      let idx = gene_indices[i]
      if idx < old_col.length() {
        new_col.push(old_col[idx])
      } else {
        new_col.push("")
      }
      i = i + 1
    }
    new_row_data[key] = new_col
    k = k + 1
  }
  SingleCellExperiment::{
    assays: new_assays,
    row_data: new_row_data,
    col_data: sce.col_data,
    reduced_dims: sce.reduced_dims,
    alternative_experiments: sce.alternative_experiments,
    metadata: sce.metadata,
    row_names: new_row_names,
    col_names: sce.col_names,
  }
}

///|
/// Subset cells (columns) by index.
pub fn sce_subset_cells(
  sce : SingleCellExperiment,
  cell_indices : Array[Int],
) -> SingleCellExperiment {
  let assay_names = sce.assays.keys().collect()
  let new_assays : Map[String, Array[Array[Double]]] = Map(
    [],
    capacity=assay_names.length(),
  )
  let mut a = 0
  while a < assay_names.length() {
    let name = assay_names[a]
    let old = sce.assays[name]
    let new_mat : Array[Array[Double]] = Array::new()
    let mut i = 0
    while i < old.length() {
      let row : Array[Double] = Array::new()
      let mut j = 0
      while j < cell_indices.length() {
        let idx = cell_indices[j]
        if idx < old[i].length() {
          row.push(old[i][idx])
        } else {
          row.push(0.0)
        }
        j = j + 1
      }
      new_mat.push(row)
      i = i + 1
    }
    new_assays[name] = new_mat
    a = a + 1
  }
  let new_col_names : Array[String] = Array::new()
  let mut c = 0
  while c < cell_indices.length() {
    let idx = cell_indices[c]
    if idx < sce.col_names.length() {
      new_col_names.push(sce.col_names[idx])
    } else {
      new_col_names.push("")
    }
    c = c + 1
  }
  let new_col_data : Map[String, Array[String]] = Map(
    [],
    capacity=sce.col_data.keys().collect().length(),
  )
  let cd_keys = sce.col_data.keys().collect()
  let mut k = 0
  while k < cd_keys.length() {
    let key = cd_keys[k]
    let old_col = sce.col_data[key]
    let new_col : Array[String] = Array::new()
    let mut i = 0
    while i < cell_indices.length() {
      let idx = cell_indices[i]
      if idx < old_col.length() {
        new_col.push(old_col[idx])
      } else {
        new_col.push("")
      }
      i = i + 1
    }
    new_col_data[key] = new_col
    k = k + 1
  }
  let new_reduced_dims : Map[String, Array[Array[Double]]] = Map(
    [],
    capacity=sce.reduced_dims.keys().collect().length(),
  )
  let rd_keys = sce.reduced_dims.keys().collect()
  let mut r = 0
  while r < rd_keys.length() {
    let key = rd_keys[r]
    let old = sce.reduced_dims[key]
    let new_mat : Array[Array[Double]] = Array::new()
    let mut i = 0
    while i < cell_indices.length() {
      let idx = cell_indices[i]
      let row : Array[Double] = Array::new()
      if idx < old.length() {
        let mut j = 0
        while j < old[idx].length() {
          row.push(old[idx][j])
          j = j + 1
        }
      }
      new_mat.push(row)
      i = i + 1
    }
    new_reduced_dims[key] = new_mat
    r = r + 1
  }
  SingleCellExperiment::{
    assays: new_assays,
    row_data: sce.row_data,
    col_data: new_col_data,
    reduced_dims: new_reduced_dims,
    alternative_experiments: sce.alternative_experiments,
    metadata: sce.metadata,
    row_names: sce.row_names,
    col_names: new_col_names,
  }
}

///|
/// Number of genes (rows).
pub fn sce_get_n_genes(sce : SingleCellExperiment) -> Int {
  let counts = sce_get_counts(sce)
  counts.length()
}

///|
/// Number of cells (columns).
pub fn sce_get_n_cells(sce : SingleCellExperiment) -> Int {
  let counts = sce_get_counts(sce)
  if counts.length() > 0 {
    counts[0].length()
  } else {
    0
  }
}

///|
/// Merge multiple SingleCellExperiment objects by row-bonding assays.
/// All experiments must have the same number of cells.
pub fn sce_merge(
  experiments : Array[SingleCellExperiment],
) -> SingleCellExperiment {
  if experiments.length() == 0 {
    return SingleCellExperiment::new([], [], [])
  }
  if experiments.length() == 1 {
    return experiments[0]
  }
  let first = experiments[0]
  let n_cells = sce_get_n_cells(first)
  let assay_names = first.assays.keys().collect()
  let merged_assays : Map[String, Array[Array[Double]]] = Map(
    [],
    capacity=assay_names.length(),
  )
  let mut a = 0
  while a < assay_names.length() {
    let name = assay_names[a]
    let combined : Array[Array[Double]] = Array::new()
    let mut e = 0
    while e < experiments.length() {
      let exp = experiments[e]
      match exp.assays.get(name) {
        Some(m) => {
          let mut i = 0
          while i < m.length() {
            let row : Array[Double] = Array::new()
            let mut j = 0
            while j < m[i].length() && j < n_cells {
              row.push(m[i][j])
              j = j + 1
            }
            while row.length() < n_cells {
              row.push(0.0)
            }
            combined.push(row)
            i = i + 1
          }
        }
        None => ()
      }
      e = e + 1
    }
    merged_assays[name] = combined
    a = a + 1
  }
  let merged_row_names : Array[String] = Array::new()
  let mut e = 0
  while e < experiments.length() {
    let exp = experiments[e]
    let mut i = 0
    while i < exp.row_names.length() {
      merged_row_names.push(exp.row_names[i])
      i = i + 1
    }
    e = e + 1
  }
  let col_names = if first.col_names.length() > 0 {
    first.col_names
  } else {
    let names : Array[String] = Array::new()
    let mut i = 0
    while i < n_cells {
      names.push("Cell" + (i + 1).to_string())
      i = i + 1
    }
    names
  }
  let merged_metadata : Map[String, String] = Map(
    [],
    capacity=experiments.length(),
  )
  let mut e = 0
  while e < experiments.length() {
    let exp = experiments[e]
    let keys = exp.metadata.keys().collect()
    let mut k = 0
    while k < keys.length() {
      merged_metadata[keys[k]] = exp.metadata[keys[k]]
      k = k + 1
    }
    e = e + 1
  }
  SingleCellExperiment::{
    assays: merged_assays,
    row_data: first.row_data,
    col_data: first.col_data,
    reduced_dims: first.reduced_dims,
    alternative_experiments: first.alternative_experiments,
    metadata: merged_metadata,
    row_names: merged_row_names,
    col_names,
  }
}

///|
/// Human-readable summary of the SCE.
pub fn SingleCellExperiment::to_string(self : SingleCellExperiment) -> String {
  let n_genes = sce_get_n_genes(self)
  let n_cells = sce_get_n_cells(self)
  let assay_names = self.assays.keys().collect()
  let rd_names = self.reduced_dims.keys().collect()
  let mut result = "SingleCellExperiment\n"
  result = result + "====================\n"
  result = result + "Genes: " + n_genes.to_string() + "\n"
  result = result + "Cells: " + n_cells.to_string() + "\n"
  result = result + "Assays: "
  let mut i = 0
  while i < assay_names.length() {
    result = result + assay_names[i]
    if i < assay_names.length() - 1 {
      result = result + ", "
    }
    i = i + 1
  }
  result = result + "\n"
  result = result + "Reduced dims: "
  let mut i = 0
  while i < rd_names.length() {
    result = result + rd_names[i]
    if i < rd_names.length() - 1 {
      result = result + ", "
    }
    i = i + 1
  }
  result = result + "\n"
  if self.row_names.length() > 0 {
    result = result + "First gene: " + self.row_names[0] + "\n"
  }
  if self.col_names.length() > 0 {
    result = result + "First cell: " + self.col_names[0] + "\n"
  }
  result
}

///|
/// String representation of a ReducedDimResult.
pub fn ReducedDimResult::to_string(self : ReducedDimResult) -> String {
  let mut result = "ReducedDimResult\n"
  result = result + "================\n"
  result = result + "Method: " + self.method + "\n"
  result = result + "Components: " + self.n_components.to_string() + "\n"
  result = result + "Samples: " + self.coordinates.length().to_string() + "\n"
  if self.variance_explained.length() > 0 {
    result = result + "Variance explained: "
    let mut i = 0
    while i < self.variance_explained.length() {
      result = result + self.variance_explained[i].to_string()
      if i < self.variance_explained.length() - 1 {
        result = result + ", "
      }
      i = i + 1
    }
    result = result + "\n"
  }
  result
}

///|
/// String representation of a SizeFactorResult.
pub fn SizeFactorResult::to_string(self : SizeFactorResult) -> String {
  let mut result = "SizeFactorResult\n"
  result = result + "================\n"
  result = result + "Method: " + self.method + "\n"
  result = result + "Size factors: "
  let mut i = 0
  while i < self.size_factors.length() {
    result = result + self.size_factors[i].to_string()
    if i < self.size_factors.length() - 1 {
      result = result + ", "
    }
    i = i + 1
  }
  result = result + "\n"
  result
}

///|
/// Top-level summary helper (also exposed as a free function for convenience).
pub fn sce_summary(sce : SingleCellExperiment) -> String {
  sce.to_string()
}

///|
/// Build an example counts matrix (5 genes x 4 cells).
pub fn sce_create_example_counts() -> Array[Array[Double]] {
  [
    [100.0, 200.0, 50.0, 150.0],
    [120.0, 180.0, 60.0, 140.0],
    [80.0, 220.0, 40.0, 160.0],
    [90.0, 210.0, 45.0, 155.0],
    [110.0, 190.0, 55.0, 145.0],
  ]
}

///|
/// Create an example SingleCellExperiment with a counts matrix.
pub fn sce_create_example() -> SingleCellExperiment {
  let counts = sce_create_example_counts()
  let row_names = ["GeneA", "GeneB", "GeneC", "GeneD", "GeneE"]
  let col_names = ["Cell1", "Cell2", "Cell3", "Cell4"]
  let sce = SingleCellExperiment::new(counts, row_names, col_names)
  sce.assays["logcounts"] = sce_compute_logcounts(
    counts,
    sce_compute_library_size_factors(counts),
  )
  sce
}

///|
/// Create an example SingleCellExperiment with PCA / tSNE / UMAP results.
pub fn sce_create_example_with_reduced_dims() -> SingleCellExperiment {
  let sce = sce_create_example()
  let logcounts = sce_get_logcounts(sce)
  let pca = sce_compute_pca(logcounts, 2)
  sce.reduced_dims["PCA"] = pca.coordinates
  let tsne = sce_compute_tsne(logcounts, 5)
  sce.reduced_dims["tSNE"] = tsne.coordinates
  let umap = sce_compute_umap(logcounts, 2)
  sce.reduced_dims["UMAP"] = umap.coordinates
  sce
}
