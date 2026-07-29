///|
/// PROGENy - Pathway activity inference from gene expression data.
/// Inspired by Bioconductor PROGENy package.
/// Uses L2-regularized linear regression (Ridge) to estimate pathway activity.

///|
/// Matrix helper: transpose a matrix (rows become columns).
pub fn progeny_transpose(matrix : Array[Array[Double]]) -> Array[Array[Double]] {
  if matrix.length() == 0 {
    return []
  }
  let n_rows = matrix.length()
  let n_cols = matrix[0].length()
  let result : Array[Array[Double]] = []
  let mut j = 0
  while j < n_cols {
    let col : Array[Double] = []
    let mut i = 0
    while i < n_rows {
      col.push(matrix[i][j])
      i = i + 1
    }
    result.push(col)
    j = j + 1
  }
  result
}

///|
/// Matrix-vector multiplication: A * x.
pub fn progeny_mat_vec(
  matrix : Array[Array[Double]],
  vector : Array[Double],
) -> Array[Double] {
  let result : Array[Double] = []
  let mut i = 0
  while i < matrix.length() {
    let row = matrix[i]
    let mut sum = 0.0
    let mut j = 0
    while j < row.length() {
      sum = sum + row[j] * vector[j]
      j = j + 1
    }
    result.push(sum)
    i = i + 1
  }
  result
}

///|
/// Matrix multiplication: A * B.
pub fn progeny_mat_mul(
  a : Array[Array[Double]],
  b : Array[Array[Double]],
) -> Array[Array[Double]] {
  if a.length() == 0 || b.length() == 0 {
    return []
  }
  let n_rows_a = a.length()
  let n_cols_a = a[0].length()
  let n_cols_b = b[0].length()
  let b_t = progeny_transpose(b)
  let result : Array[Array[Double]] = []
  let mut i = 0
  while i < n_rows_a {
    let row_a = a[i]
    let row_result : Array[Double] = []
    let mut j = 0
    while j < n_cols_b {
      let col_b = b_t[j]
      let mut sum = 0.0
      let mut k = 0
      while k < n_cols_a {
        sum = sum + row_a[k] * col_b[k]
        k = k + 1
      }
      row_result.push(sum)
      j = j + 1
    }
    result.push(row_result)
    i = i + 1
  }
  result
}

///|
/// Create identity matrix of size n.
pub fn progeny_identity(n : Int, scale : Double) -> Array[Array[Double]] {
  let result : Array[Array[Double]] = []
  let mut i = 0
  while i < n {
    let row : Array[Double] = []
    let mut j = 0
    while j < n {
      if i == j {
        row.push(scale)
      } else {
        row.push(0.0)
      }
      j = j + 1
    }
    result.push(row)
    i = i + 1
  }
  result
}

///|
/// Add two matrices element-wise.
pub fn progeny_mat_add(
  a : Array[Array[Double]],
  b : Array[Array[Double]],
) -> Array[Array[Double]] {
  let result : Array[Array[Double]] = []
  let mut i = 0
  while i < a.length() {
    let row_a = a[i]
    let row_b = b[i]
    let row_result : Array[Double] = []
    let mut j = 0
    while j < row_a.length() {
      row_result.push(row_a[j] + row_b[j])
      j = j + 1
    }
    result.push(row_result)
    i = i + 1
  }
  result
}

///|
/// Gauss-Jordan elimination for matrix inversion.
pub fn progeny_invert(matrix : Array[Array[Double]]) -> Array[Array[Double]] {
  let n = matrix.length()
  if n == 0 {
    return []
  }
  // Augmented matrix [A | I]
  let aug : Array[Array[Double]] = []
  let mut i = 0
  while i < n {
    let row : Array[Double] = []
    let mut j = 0
    while j < n {
      row.push(matrix[i][j])
      j = j + 1
    }
    let mut j2 = 0
    while j2 < n {
      if i == j2 {
        row.push(1.0)
      } else {
        row.push(0.0)
      }
      j2 = j2 + 1
    }
    aug.push(row)
    i = i + 1
  }

  // Forward elimination with partial pivoting
  let mut col = 0
  while col < n {
    // Find pivot
    let mut max_val = aug[col][col]
    let mut max_row = col
    let mut row = col + 1
    while row < n {
      let val = aug[row][col]
      if val > max_val {
        max_val = val
        max_row = row
      }
      row = row + 1
    }
    // Swap rows
    let tmp = aug[col]
    aug[col] = aug[max_row]
    aug[max_row] = tmp

    // Scale pivot row
    let pivot = aug[col][col]
    if pivot == 0.0 {
      // Singular matrix, add small regularization
      let mut j = 0
      while j < 2 * n {
        aug[col][j] = aug[col][j] + 0.0001
        j = j + 1
      }
    }
    let pivot2 = aug[col][col]
    let mut j = 0
    while j < 2 * n {
      aug[col][j] = aug[col][j] / pivot2
      j = j + 1
    }

    // Eliminate column
    let mut r = 0
    while r < n {
      if r != col {
        let factor = aug[r][col]
        let mut j2 = 0
        while j2 < 2 * n {
          aug[r][j2] = aug[r][j2] - factor * aug[col][j2]
          j2 = j2 + 1
        }
      }
      r = r + 1
    }
    col = col + 1
  }

  // Extract inverse from augmented matrix
  let inverse : Array[Array[Double]] = []
  let mut i2 = 0
  while i2 < n {
    let row : Array[Double] = []
    let mut j = n
    while j < 2 * n {
      row.push(aug[i2][j])
      j = j + 1
    }
    inverse.push(row)
    i2 = i2 + 1
  }
  inverse
}

///|
pub struct PathwayGeneSet {
  pathway_name : String
  genes : Array[String]
  weights : Array[Double]
}

///|
pub fn PathwayGeneSet::new(
  pathway_name : String,
  genes : Array[String],
  weights : Array[Double],
) -> PathwayGeneSet {
  PathwayGeneSet::{ pathway_name, genes, weights }
}

///|
pub struct ProgenyData {
  expression : Array[Array[Double]]
  sample_names : Array[String]
  gene_names : Array[String]
  pathways : Array[PathwayGeneSet]
  lambda : Double
}

///|
pub fn ProgenyData::new(
  expression : Array[Array[Double]],
  sample_names : Array[String],
  gene_names : Array[String],
  pathways : Array[PathwayGeneSet],
  lambda : Double,
) -> ProgenyData {
  ProgenyData::{ expression, sample_names, gene_names, pathways, lambda }
}

///|
pub struct ProgenyResult {
  sample_names : Array[String]
  pathway_names : Array[String]
  activity_matrix : Array[Array[Double]]
}

///|
pub fn ProgenyResult::new(
  sample_names : Array[String],
  pathway_names : Array[String],
  activity_matrix : Array[Array[Double]],
) -> ProgenyResult {
  ProgenyResult::{ sample_names, pathway_names, activity_matrix }
}

///|
/// Build the pathway weight matrix (genes × pathways).
pub fn build_pathway_matrix(
  gene_names : Array[String],
  pathways : Array[PathwayGeneSet],
) -> Array[Array[Double]] {
  let n_genes = gene_names.length()
  let n_pathways = pathways.length()
  let matrix : Array[Array[Double]] = []
  let mut p = 0
  while p < n_pathways {
    let pathway = pathways[p]
    let row : Array[Double] = []
    let mut g = 0
    while g < n_genes {
      let gene_name = gene_names[g]
      let mut weight = 0.0
      let mut idx = 0
      while idx < pathway.genes.length() {
        if pathway.genes[idx] == gene_name {
          weight = pathway.weights[idx]
          break
        }
        idx = idx + 1
      }
      row.push(weight)
      g = g + 1
    }
    matrix.push(row)
    p = p + 1
  }
  // Transpose to get genes × pathways
  progeny_transpose(matrix)
}

///|
/// Run PROGENy to infer pathway activity from gene expression.
/// Uses L2-regularized linear regression (Ridge regression).
/// Solution: β = (X^T X + λI)^(-1) X^T y
pub fn run_progeny(data : ProgenyData) -> ProgenyResult {
  let n_samples = data.sample_names.length()
  let n_genes = data.gene_names.length()
  let n_pathways = data.pathways.length()

  // Build pathway weight matrix X (genes × pathways)
  let x = build_pathway_matrix(data.gene_names, data.pathways)

  // Compute X^T X (pathways × pathways)
  let x_t = progeny_transpose(x)
  let x_t_x = progeny_mat_mul(x_t, x)

  // Add regularization: X^T X + λI
  let reg = progeny_identity(n_pathways, data.lambda)
  let x_t_x_reg = progeny_mat_add(x_t_x, reg)

  // Compute (X^T X + λI)^(-1)
  let x_t_x_inv = progeny_invert(x_t_x_reg)

  // Compute X^T y for each sample
  // X^T is (pathways × genes), expression is (samples × genes)
  // Transpose expression to (genes × samples) for multiplication
  let x_t_y_base = progeny_mat_mul(x_t, progeny_transpose(data.expression))
  // x_t_y_base is (pathways × samples)

  // Compute pathway activity: β = (X^T X + λI)^(-1) * X^T y
  let activity_matrix = progeny_mat_mul(x_t_x_inv, x_t_y_base)

  // activity_matrix is (pathways × samples), transpose to (samples × pathways)
  let result_matrix = progeny_transpose(activity_matrix)

  // Extract pathway names
  let pathway_names : Array[String] = []
  let mut p = 0
  while p < n_pathways {
    pathway_names.push(data.pathways[p].pathway_name)
    p = p + 1
  }

  ProgenyResult::{
    sample_names: data.sample_names.copy(),
    pathway_names,
    activity_matrix: result_matrix,
  }
}

///|
/// Get pathway activity for a specific pathway.
pub fn get_pathway_activity(
  result : ProgenyResult,
  pathway_name : String,
) -> Array[Double] {
  let mut idx = 0
  while idx < result.pathway_names.length() {
    if result.pathway_names[idx] == pathway_name {
      let activity : Array[Double] = []
      let mut i = 0
      while i < result.activity_matrix.length() {
        activity.push(result.activity_matrix[i][idx])
        i = i + 1
      }
      return activity
    }
    idx = idx + 1
  }
  []
}

///|
/// Get sample activity profile for a specific sample.
pub fn get_sample_profile(
  result : ProgenyResult,
  sample_name : String,
) -> Array[Double] {
  let mut idx = 0
  while idx < result.sample_names.length() {
    if result.sample_names[idx] == sample_name {
      return result.activity_matrix[idx].copy()
    }
    idx = idx + 1
  }
  []
}
