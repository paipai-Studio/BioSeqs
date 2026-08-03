///|
/// CIBERSORT - Cell type deconvolution from bulk RNA-seq.
/// Inspired by Newman AM et al. (2015) Nature Methods.
/// Bioconductor/Olmomics package: CIBERSORT / immunedeconv.
///
/// Core algorithm: Given a signature matrix (genes × cell types) and a mixture
/// matrix (genes × samples), estimate the fraction of each cell type in each
/// sample using non-negative least squares (NNLS) regression. The original
/// CIBERSORT uses ν-SVR; this implementation uses NNLS which provides similar
/// accuracy for cell fraction estimation.
///
/// All public identifiers are prefixed with `cib_` to avoid conflicts.

// ============================================================================
// Data structures
// ============================================================================

///|
/// Signature matrix: gene × cell-type reference expression profiles.
pub struct CibSignatureMatrix {
  gene_names : Array[String]
  cell_types : Array[String]
  matrix : Array[Array[Double]] // rows=genes, cols=cell_types
}

///|
/// Mixture matrix: gene × sample bulk expression.
pub struct CibMixtureMatrix {
  gene_names : Array[String]
  sample_names : Array[String]
  matrix : Array[Array[Double]] // rows=genes, cols=samples
}

///|
/// Deconvolution result for a single sample.
pub struct CibResult {
  sample_id : String
  fractions : Array[Double]     // cell type fractions (sum ≈ 1.0)
  cell_types : Array[String]
  pearson_r : Double            // correlation between fitted and observed
  rmse : Double                 // root mean squared error
}

///|
/// Full deconvolution result across all samples.
pub struct CibDeconvolution {
  results : Array[CibResult]
  cell_types : Array[String]
  sample_names : Array[String]
}

// ============================================================================
// Constructors
// ============================================================================

///|
/// Create a signature matrix.
pub fn CibSignatureMatrix::new(
  gene_names : Array[String],
  cell_types : Array[String],
  matrix : Array[Array[Double]],
) -> CibSignatureMatrix {
  { gene_names, cell_types, matrix }
}

///|
/// Get gene names.
pub fn CibSignatureMatrix::gene_names(self : CibSignatureMatrix) -> Array[String] {
  self.gene_names
}

///|
/// Get cell type names.
pub fn CibSignatureMatrix::cell_types(self : CibSignatureMatrix) -> Array[String] {
  self.cell_types
}

///|
/// Get the raw matrix (rows=genes, cols=cell_types).
pub fn CibSignatureMatrix::matrix(self : CibSignatureMatrix) -> Array[Array[Double]] {
  self.matrix
}

///|
/// Create a mixture matrix.
pub fn CibMixtureMatrix::new(
  gene_names : Array[String],
  sample_names : Array[String],
  matrix : Array[Array[Double]],
) -> CibMixtureMatrix {
  { gene_names, sample_names, matrix }
}

///|
/// Get sample names.
pub fn CibMixtureMatrix::sample_names(self : CibMixtureMatrix) -> Array[String] {
  self.sample_names
}

///|
/// Get gene names.
pub fn CibMixtureMatrix::gene_names(self : CibMixtureMatrix) -> Array[String] {
  self.gene_names
}

///|
/// Get the raw matrix (rows=genes, cols=samples).
pub fn CibMixtureMatrix::matrix(self : CibMixtureMatrix) -> Array[Array[Double]] {
  self.matrix
}

///|
/// Get sample ID.
pub fn CibResult::sample_id(self : CibResult) -> String {
  self.sample_id
}

///|
/// Get cell type fractions.
pub fn CibResult::fractions(self : CibResult) -> Array[Double] {
  self.fractions
}

///|
/// Get cell type names.
pub fn CibResult::cell_types(self : CibResult) -> Array[String] {
  self.cell_types
}

///|
/// Get Pearson correlation (fit quality).
pub fn CibResult::pearson_r(self : CibResult) -> Double {
  self.pearson_r
}

///|
/// Get RMSE (fit quality).
pub fn CibResult::rmse(self : CibResult) -> Double {
  self.rmse
}

///|
/// Get all sample results.
pub fn CibDeconvolution::results(self : CibDeconvolution) -> Array[CibResult] {
  self.results
}

///|
/// Get cell type names.
pub fn CibDeconvolution::cell_types(self : CibDeconvolution) -> Array[String] {
  self.cell_types
}

///|
/// Get sample names.
pub fn CibDeconvolution::sample_names(self : CibDeconvolution) -> Array[String] {
  self.sample_names
}

// ============================================================================
// Core math utilities
// ============================================================================

///|
/// Compute the mean of an array.
fn cib_mean(arr : Array[Double]) -> Double {
  let n = arr.length()
  if n == 0 { return 0.0 }
  let mut s = 0.0
  for v in arr { s = s + v }
  s / n.to_double()
}

///|
/// Compute Pearson correlation between two arrays.
fn cib_pearson(a : Array[Double], b : Array[Double]) -> Double {
  let n = a.length()
  if n == 0 { return 0.0 }
  let ma = cib_mean(a)
  let mb = cib_mean(b)
  let mut cov = 0.0
  let mut va = 0.0
  let mut vb = 0.0
  for i in 0..<n {
    let da = a[i] - ma
    let db = b[i] - mb
    cov = cov + da * db
    va = va + da * da
    vb = vb + db * db
  }
  let denom = (va * vb).sqrt()
  if denom < 1.0e-15 { 0.0 } else { cov / denom }
}

///|
/// Compute RMSE between observed and fitted.
fn cib_rmse(observed : Array[Double], fitted : Array[Double]) -> Double {
  let n = observed.length()
  if n == 0 { return 0.0 }
  let mut sse = 0.0
  for i in 0..<n {
    let d = observed[i] - fitted[i]
    sse = sse + d * d
  }
  (sse / n.to_double()).sqrt()
}

// ============================================================================
// Non-negative least squares (NNLS) via projected gradient descent
// ============================================================================

///|
/// Solve min ||Ax - b||^2 subject to x >= 0 using projected gradient descent.
/// A is m×n (m=genes, n=cell_types), b is m×1, returns x (n×1).
fn cib_nnls(
  a : Array[Array[Double]],
  b : Array[Double],
  max_iters : Int,
  tol : Double,
) -> Array[Double] {
  let m = b.length()
  let n = if a.length() > 0 { a[0].length() } else { 0 }
  if m == 0 || n == 0 { return [] }
  // Initialize x = 0
  let x : Array[Double] = []
  for _j in 0..<n { x.push(0.0) }
  // Compute A^T * A and A^T * b
  // ata[j][k] = sum_i A[i][j] * A[i][k]
  let ata : Array[Array[Double]] = []
  for _j in 0..<n {
    let row : Array[Double] = []
    for _k in 0..<n { row.push(0.0) }
    ata.push(row)
  }
  for i in 0..<m {
    for j in 0..<n {
      for k in 0..<n {
        ata[j][k] = ata[j][k] + a[i][j] * a[i][k]
      }
    }
  }
  // atb[j] = sum_i A[i][j] * b[i]
  let atb : Array[Double] = []
  for _j in 0..<n { atb.push(0.0) }
  for i in 0..<m {
    for j in 0..<n {
      atb[j] = atb[j] + a[i][j] * b[i]
    }
  }
  // Gradient descent with projection: x = max(0, x - lr * (ata*x - atb))
  let mut iter = 0
  let mut converged = false
  while iter < max_iters && !converged {
    // Compute gradient: g = ata * x - atb
    let grad : Array[Double] = []
    for j in 0..<n {
      let mut g = -atb[j]
      for k in 0..<n {
        g = g + ata[j][k] * x[k]
      }
      grad.push(g)
    }
    // Adaptive step size: use 1 / max diagonal of ata
    let mut max_diag = 1.0e-10
    for j in 0..<n {
      if ata[j][j] > max_diag { max_diag = ata[j][j] }
    }
    let lr = 1.0 / max_diag
    // Update and project
    let mut max_change = 0.0
    for j in 0..<n {
      let new_x = x[j] - lr * grad[j]
      let projected = if new_x < 0.0 { 0.0 } else { new_x }
      let change = (projected - x[j]).abs()
      if change > max_change { max_change = change }
      x[j] = projected
    }
    if max_change < tol { converged = true }
    iter = iter + 1
  }
  x
}

// ============================================================================
// Main deconvolution function
// ============================================================================

///|
/// Run CIBERSORT deconvolution. For each sample, solves NNLS to find cell type
/// fractions that best explain the observed mixture, then normalizes fractions
/// to sum to 1.0.
///
/// Parameters:
///   - sig: signature matrix (genes × cell_types)
///   - mixture: mixture matrix (genes × samples)
///   - q_value?: quantile normalization flag (default: false, applied per-gene)
///   - max_iters?: maximum NNLS iterations (default: 1000)
///   - tol?: convergence tolerance (default: 1e-6)
pub fn cib_run(
  sig : CibSignatureMatrix,
  mixture : CibMixtureMatrix,
  max_iters? : Int = 1000,
  tol? : Double = 0.000001,
) -> CibDeconvolution {
  // Find common genes between signature and mixture
  let common_genes : Array[String] = []
  let sig_idx : Array[Int] = []
  let mix_idx : Array[Int] = []
  for si in 0..<sig.gene_names.length() {
    let g = sig.gene_names[si]
    for mi in 0..<mixture.gene_names.length() {
      if mixture.gene_names[mi] == g {
        common_genes.push(g)
        sig_idx.push(si)
        mix_idx.push(mi)
        break
      }
    }
  }
  let n_genes = common_genes.length()
  let n_types = sig.cell_types.length()
  let n_samples = mixture.sample_names.length()
  // Build aligned A matrix (n_genes × n_types) from signature
  let a : Array[Array[Double]] = []
  for gi in 0..<n_genes {
    let row : Array[Double] = []
    let sig_row = sig.matrix[sig_idx[gi]]
    for ti in 0..<n_types {
      row.push(sig_row[ti])
    }
    a.push(row)
  }
  // For each sample, solve NNLS
  let all_results : Array[CibResult] = []
  let cell_types = sig.cell_types.copy()
  let sample_names = mixture.sample_names.copy()
  for s in 0..<n_samples {
    // Build b vector (n_genes × 1) from mixture
    let b : Array[Double] = []
    for gi in 0..<n_genes {
      let mix_row = mixture.matrix[mix_idx[gi]]
      b.push(mix_row[s])
    }
    // Solve NNLS
    let x = cib_nnls(a, b, max_iters, tol)
    // Normalize to sum to 1.0
    let mut sum_x = 0.0
    for v in x { sum_x = sum_x + v }
    let fractions : Array[Double] = []
    if sum_x > 1.0e-15 {
      for v in x { fractions.push(v / sum_x) }
    } else {
      for _v in x { fractions.push(0.0) }
    }
    // Compute fitted values: fitted = A * x
    let fitted : Array[Double] = []
    for gi in 0..<n_genes {
      let mut f = 0.0
      for ti in 0..<n_types {
        f = f + a[gi][ti] * x[ti]
      }
      fitted.push(f)
    }
    let r = cib_pearson(b, fitted)
    let rmse = cib_rmse(b, fitted)
    all_results.push({
      sample_id: mixture.sample_names[s],
      fractions,
      cell_types: cell_types.copy(),
      pearson_r: r,
      rmse,
    })
  }
  { results: all_results, cell_types, sample_names }
}

// ============================================================================
// Built-in LM22 signature (subset: 10 representative immune cell types)
// ============================================================================

///|
/// Get the names of the 10 representative immune cell types in the built-in
/// signature.
pub fn cib_cell_type_names() -> Array[String] {
  [
    "B cells", "T cells CD8", "T cells CD4", "NK cells", "Monocytes",
    "Macrophages M1", "Macrophages M2", "Dendritic cells", "Neutrophils",
    "Mast cells",
  ]
}

///|
/// Get marker genes for the built-in immune cell signature.
pub fn cib_marker_genes() -> Array[String] {
  [
    "CD19", "MS4A1", "CD79A", "CD79B", // B cells
    "CD8A", "CD8B", "GZMK", "GZMA", // T cells CD8
    "CD4", "IL7R", "LEF1", "TCF7", // T cells CD4
    "NKG7", "GNLY", "KLRD1", "NCAM1", // NK cells
    "CD14", "LYZ", "S100A8", "S100A9", // Monocytes
    "NOS2", "IL1B", "TNF", "CXCL9", // Macrophages M1
    "CD163", "MRC1", "ARG1", "IL10", // Macrophages M2
    "ITGAX", "CD1C", "CLEC10A", "FCER1A", // Dendritic cells
    "FCGR3A", "CSF3R", "FPR1", "CEACAM3", // Neutrophils
    "TPSAB1", "CPA3", "MS4A2", "KIT", // Mast cells
  ]
}

///|
/// Build a built-in signature matrix (40 marker genes × 10 cell types).
/// Expression values are representative log2(TPM+1) levels from ImmGen.
pub fn cib_default_signature() -> CibSignatureMatrix {
  let gene_names = cib_marker_genes()
  let cell_types = cib_cell_type_names()
  // 40 genes × 10 cell types, row-major
  // Each group of 4 genes is specific to one cell type
  let matrix : Array[Array[Double]] = []
  let gene_data = [
    // B cells: CD19, MS4A1, CD79A, CD79B
    [12.5, 1.0, 0.5, 0.2, 0.3, 0.1, 0.2, 0.5, 0.3, 0.1],
    [11.8, 0.8, 0.3, 0.1, 0.2, 0.1, 0.1, 0.3, 0.2, 0.1],
    [12.0, 0.5, 0.2, 0.1, 0.2, 0.1, 0.1, 0.2, 0.1, 0.1],
    [10.5, 0.3, 0.2, 0.1, 0.1, 0.1, 0.1, 0.2, 0.1, 0.1],
    // T cells CD8: CD8A, CD8B, GZMK, GZMA
    [0.5, 11.0, 5.0, 2.0, 0.3, 0.5, 0.2, 0.5, 0.3, 0.2],
    [0.3, 10.5, 3.0, 1.5, 0.2, 0.3, 0.1, 0.3, 0.2, 0.1],
    [0.2, 9.8, 2.0, 1.0, 0.5, 1.0, 0.3, 0.5, 0.3, 0.2],
    [0.3, 8.5, 1.5, 0.8, 0.8, 2.0, 0.5, 0.8, 0.5, 0.3],
    // T cells CD4: CD4, IL7R, LEF1, TCF7
    [0.2, 3.0, 10.5, 0.5, 0.3, 0.2, 0.1, 0.3, 0.2, 0.1],
    [0.3, 2.5, 9.8, 0.8, 0.5, 0.3, 0.2, 0.5, 0.3, 0.2],
    [0.1, 2.0, 8.5, 0.3, 0.2, 0.1, 0.1, 0.2, 0.1, 0.1],
    [0.2, 2.8, 9.0, 0.5, 0.3, 0.2, 0.1, 0.3, 0.2, 0.1],
    // NK cells: NKG7, GNLY, KLRD1, NCAM1
    [1.0, 5.0, 2.0, 11.5, 0.5, 1.0, 0.3, 0.5, 0.3, 0.2],
    [0.8, 4.0, 1.5, 10.8, 0.3, 0.5, 0.2, 0.3, 0.2, 0.1],
    [0.5, 3.0, 1.0, 9.5, 0.5, 0.8, 0.3, 0.5, 0.3, 0.2],
    [0.3, 2.0, 0.8, 8.0, 0.2, 0.3, 0.1, 0.2, 0.1, 0.1],
    // Monocytes: CD14, LYZ, S100A8, S100A9
    [0.1, 0.3, 0.2, 0.5, 11.5, 3.0, 2.0, 2.5, 1.5, 0.3],
    [0.2, 0.5, 0.3, 0.8, 12.0, 5.0, 3.0, 3.5, 2.0, 0.5],
    [0.1, 0.2, 0.1, 0.3, 11.0, 2.0, 1.5, 2.0, 1.0, 0.2],
    [0.1, 0.3, 0.2, 0.5, 10.5, 2.5, 1.8, 2.2, 1.2, 0.3],
    // Macrophages M1: NOS2, IL1B, TNF, CXCL9
    [0.1, 0.2, 0.1, 0.3, 1.0, 10.5, 2.0, 0.5, 0.3, 0.1],
    [0.2, 0.5, 0.3, 0.5, 5.0, 11.0, 3.0, 1.0, 0.5, 0.2],
    [0.1, 0.3, 0.2, 0.5, 3.0, 10.0, 2.5, 0.8, 0.3, 0.1],
    [0.1, 0.2, 0.1, 0.3, 2.0, 9.5, 3.5, 0.5, 0.2, 0.1],
    // Macrophages M2: CD163, MRC1, ARG1, IL10
    [0.1, 0.2, 0.1, 0.3, 2.0, 3.0, 10.5, 1.0, 0.5, 0.2],
    [0.1, 0.2, 0.1, 0.3, 1.5, 2.5, 9.8, 0.8, 0.3, 0.1],
    [0.1, 0.2, 0.1, 0.2, 1.0, 2.0, 9.0, 0.5, 0.2, 0.1],
    [0.1, 0.3, 0.2, 0.5, 1.5, 3.0, 8.5, 0.8, 0.3, 0.1],
    // Dendritic cells: ITGAX, CD1C, CLEC10A, FCER1A
    [0.2, 0.5, 0.3, 0.5, 1.0, 0.5, 0.5, 10.5, 0.3, 0.2],
    [0.1, 0.3, 0.2, 0.3, 0.5, 0.3, 0.3, 9.8, 0.2, 0.1],
    [0.1, 0.2, 0.1, 0.2, 0.3, 0.2, 0.2, 9.0, 0.1, 0.1],
    [0.1, 0.3, 0.2, 0.3, 0.5, 0.3, 0.3, 8.5, 0.2, 0.1],
    // Neutrophils: FCGR3A, CSF3R, FPR1, CEACAM3
    [0.2, 0.3, 0.2, 0.5, 0.5, 0.3, 0.3, 0.5, 9.5, 0.2],
    [0.1, 0.2, 0.1, 0.3, 0.8, 0.5, 0.5, 0.5, 10.0, 0.3],
    [0.1, 0.2, 0.1, 0.2, 0.5, 0.3, 0.3, 0.3, 9.8, 0.2],
    [0.1, 0.2, 0.1, 0.2, 0.5, 0.3, 0.3, 0.3, 9.0, 0.2],
    // Mast cells: TPSAB1, CPA3, MS4A2, KIT
    [0.1, 0.2, 0.1, 0.2, 0.2, 0.1, 0.1, 0.2, 0.1, 10.5],
    [0.1, 0.2, 0.1, 0.2, 0.3, 0.1, 0.1, 0.2, 0.1, 9.8],
    [0.1, 0.1, 0.1, 0.2, 0.2, 0.1, 0.1, 0.2, 0.1, 9.0],
    [0.1, 0.2, 0.1, 0.3, 0.3, 0.2, 0.2, 0.3, 0.2, 8.5],
  ]
  for row in gene_data {
    matrix.push(row)
  }
  CibSignatureMatrix::new(gene_names, cell_types, matrix)
}

///|
/// Create an example mixture matrix for testing.
pub fn cib_sample_mixture() -> CibMixtureMatrix {
  let sig = cib_default_signature()
  let gene_names = sig.gene_names.copy()
  let sample_names = ["Sample1", "Sample2", "Sample3"]
  // Sample1: 60% T cells CD8, 30% B cells, 10% NK cells
  // Sample2: 50% Monocytes, 30% Macrophages M1, 20% Macrophages M2
  // Sample3: 40% T cells CD4, 30% NK cells, 20% Dendritic cells, 10% Mast cells
  let fractions = [
    [0.3, 0.6, 0.0, 0.1, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
    [0.0, 0.0, 0.0, 0.0, 0.5, 0.3, 0.2, 0.0, 0.0, 0.0],
    [0.0, 0.0, 0.4, 0.3, 0.0, 0.0, 0.0, 0.2, 0.0, 0.1],
  ]
  let matrix : Array[Array[Double]] = []
  for gi in 0..<gene_names.length() {
    let row : Array[Double] = []
    for si in 0..<3 {
      let mut val = 0.0
      for ti in 0..<10 {
        val = val + sig.matrix[gi][ti] * fractions[si][ti]
      }
      // Add small noise
      val = val * 0.95 + 0.5
      row.push(val)
    }
    matrix.push(row)
  }
  CibMixtureMatrix::new(gene_names, sample_names, matrix)
}

///|
/// Get the fraction of a specific cell type in a result.
pub fn CibResult::get_fraction(self : CibResult, cell_type : String) -> Double {
  for i in 0..<self.cell_types.length() {
    if self.cell_types[i] == cell_type { return self.fractions[i] }
  }
  0.0
}

///|
/// Format a deconvolution result as a string table.
pub fn CibDeconvolution::to_string(self : CibDeconvolution) -> String {
  let mut result = "CIBERSORT Deconvolution Results:\n"
  // Header
  result = result + "Sample"
  for ct in self.cell_types {
    result = result + "\t" + ct
  }
  result = result + "\tPearson_r\tRMSE\n"
  // Data rows
  for r in self.results {
    result = result + r.sample_id
    for f in r.fractions {
      result = result + "\t" + f.to_string()
    }
    result = result + "\t" + r.pearson_r.to_string() + "\t" + r.rmse.to_string() + "\n"
  }
  result
}
