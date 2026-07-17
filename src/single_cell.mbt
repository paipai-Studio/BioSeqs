///|
/// Single Cell Analysis module, inspired by Bioconductor SingleCellExperiment.
/// Provides functions for single-cell RNA-seq data analysis.

///|
pub struct SingleCellExperiment {
  assays : Array[Array[Array[Double]]]
  col_data : Array[Array[String]]
  row_data : Array[Array[String]]
  reduced_dimensions : Array[(String, Array[Array[Double]])]
}

///|
pub fn SingleCellExperiment::new(
  assays : Array[Array[Array[Double]]],
  col_data : Array[Array[String]],
  row_data : Array[Array[String]],
) -> SingleCellExperiment {
  SingleCellExperiment::{
    assays,
    col_data,
    row_data,
    reduced_dimensions: Array::new(),
  }
}

///|
pub struct QCMetrics {
  n_umi : Array[Double]
  n_genes : Array[Int]
  mito_percent : Array[Double]
}

///|
pub fn QCMetrics::new(
  n_umi : Array[Double],
  n_genes : Array[Int],
  mito_percent : Array[Double],
) -> QCMetrics {
  QCMetrics::{ n_umi, n_genes, mito_percent }
}

///|
pub fn calculate_qc_metrics(
  counts : Array[Array[Double]],
  gene_names : Array[String],
) -> QCMetrics {
  let n_cells = counts.length()
  let n_genes = if n_cells > 0 { counts[0].length() } else { 0 }

  let n_umi : Array[Double] = Array::new()
  let n_gene_arr : Array[Int] = Array::new()
  let mito_percent_arr : Array[Double] = Array::new()

  let mut i = 0
  while i < n_cells {
    let mut total = 0.0
    let mut detected = 0
    let mut mito_total = 0.0

    let mut j = 0
    while j < n_genes {
      let count = counts[i][j]
      total = total + count
      if count > 0.0 {
        detected = detected + 1
      }
      if j < gene_names.length() {
        let gene_name = gene_names[j]
        let gene_upper = gene_name.to_upper()
        if gene_upper.starts_with("MT-") || gene_upper.starts_with("MT") {
          mito_total = mito_total + count
        }
      }
      j = j + 1
    }

    n_umi.push(total)
    n_gene_arr.push(detected)

    let mito_pct = if total > 0.0 { mito_total / total * 100.0 } else { 0.0 }
    mito_percent_arr.push(mito_pct)

    i = i + 1
  }

  QCMetrics::new(n_umi, n_gene_arr, mito_percent_arr)
}

///|
pub fn filter_cells(
  counts : Array[Array[Double]],
  qc : QCMetrics,
  min_n_umi : Int,
  max_n_umi : Int,
  min_n_genes : Int,
  max_mito_percent : Double,
) -> Array[Array[Double]] {
  let filtered : Array[Array[Double]] = Array::new()

  let mut i = 0
  while i < counts.length() {
    if qc.n_umi[i] >= min_n_umi.to_double() &&
      qc.n_umi[i] <= max_n_umi.to_double() &&
      qc.n_genes[i] >= min_n_genes &&
      qc.mito_percent[i] <= max_mito_percent {
      filtered.push(counts[i])
    }
    i = i + 1
  }

  filtered
}

///|
pub fn log_normalize(
  counts : Array[Array[Double]],
  size_factors : Array[Double]?,
) -> Array[Array[Double]] {
  let n_cells = counts.length()
  let n_genes = if n_cells > 0 { counts[0].length() } else { 0 }

  let sf = match size_factors {
    None => {
      let sf_arr : Array[Double] = Array::new()
      let mut i = 0
      while i < n_cells {
        let mut sum = 0.0
        let mut j = 0
        while j < n_genes {
          sum = sum + counts[i][j]
          j = j + 1
        }
        sf_arr.push(sum / 10000.0)
        i = i + 1
      }
      sf_arr
    }
    Some(sf_val) => sf_val
  }

  let normalized : Array[Array[Double]] = Array::new()

  let mut i = 0
  while i < n_cells {
    let row : Array[Double] = Array::new()
    let mut j = 0
    while j < n_genes {
      let norm = if sf[i] > 0.0 { counts[i][j] / sf[i] } else { 0.0 }
      let log_val = if norm > 0.0 { @math.ln(norm) + 1.0 } else { 0.0 }
      row.push(log_val)
      j = j + 1
    }
    normalized.push(row)
    i = i + 1
  }

  normalized
}

///|
fn matmul(
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
fn matmul_vec(m : Array[Array[Double]], v : Array[Double]) -> Array[Double] {
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
fn transpose(m : Array[Array[Double]]) -> Array[Array[Double]] {
  let rows = m.length()
  let cols = if rows > 0 { m[0].length() } else { 0 }

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
fn vec_norm(v : Array[Double]) -> Double {
  let mut sum = 0.0
  let mut i = 0
  while i < v.length() {
    sum = sum + v[i] * v[i]
    i = i + 1
  }
  sum.sqrt()
}

///|
fn center_matrix(m : Array[Array[Double]]) -> Array[Array[Double]] {
  let rows = m.length()
  let cols = if rows > 0 { m[0].length() } else { 0 }

  let means : Array[Double] = Array::new()
  let mut j = 0
  while j < cols {
    let mut sum = 0.0
    let mut i = 0
    while i < rows {
      sum = sum + m[i][j]
      i = i + 1
    }
    means.push(sum / rows.to_double())
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
pub fn pca_power_iteration(
  data : Array[Array[Double]],
  n_components : Int,
) -> (Array[Array[Double]], Array[Double]) {
  let centered = center_matrix(data)
  let mut cov = matmul(transpose(centered), centered)

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
      let v_new = matmul_vec(cov, v)
      let norm = vec_norm(v_new)

      let mut k = 0
      while k < v.length() {
        v[k] = if norm > 0.0 { v_new[k] / norm } else { 0.0 }
        k = k + 1
      }

      iter = iter + 1
    }

    let eigenvalue = vec_norm(matmul_vec(cov, v))

    let outer : Array[Array[Double]] = Array::new()
    let mut k = 0
    while k < n {
      let row : Array[Double] = Array::new()
      let mut l = 0
      while l < n {
        row.push(v[k] * v[l])
        l = l + 1
      }
      outer.push(row)
      k = k + 1
    }

    let new_cov : Array[Array[Double]] = Array::new()
    k = 0
    while k < n {
      let row : Array[Double] = Array::new()
      let mut l = 0
      while l < n {
        row.push(cov[k][l] - eigenvalue * outer[k][l])
        l = l + 1
      }
      new_cov.push(row)
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

  let transformed = matmul(centered, transpose(components))

  (transformed, eigenvalues)
}

///|
pub fn find_variable_genes(
  normalized : Array[Array[Double]],
  n_top_genes : Int,
) -> Array[Int] {
  let n_genes = if normalized.length() > 0 { normalized[0].length() } else { 0 }
  let n_cells = normalized.length()

  let means : Array[Double] = Array::new()
  let vars : Array[Double] = Array::new()

  let mut i = 0
  while i < n_genes {
    let mut sum = 0.0
    let mut j = 0
    while j < n_cells {
      sum = sum + normalized[j][i]
      j = j + 1
    }
    let mean = sum / n_cells.to_double()
    means.push(mean)

    let mut var_sum = 0.0
    j = 0
    while j < n_cells {
      let diff = normalized[j][i] - mean
      var_sum = var_sum + diff * diff
      j = j + 1
    }
    vars.push(var_sum / (n_cells - 1).to_double())

    i = i + 1
  }

  let gene_vars : Array[(Double, Int)] = Array::new()
  i = 0
  while i < n_genes {
    gene_vars.push((vars[i], i))
    i = i + 1
  }

  let mut j = 0
  while j < gene_vars.length() - 1 {
    let mut k = j + 1
    while k < gene_vars.length() {
      if gene_vars[j].0 < gene_vars[k].0 {
        let tmp = gene_vars[j]
        gene_vars[j] = gene_vars[k]
        gene_vars[k] = tmp
      }
      k = k + 1
    }
    j = j + 1
  }

  let top_indices : Array[Int] = Array::new()
  i = 0
  while i < n_top_genes && i < gene_vars.length() {
    top_indices.push(gene_vars[i].1)
    i = i + 1
  }

  top_indices
}

///|
pub fn select_genes_by_indices(
  data : Array[Array[Double]],
  indices : Array[Int],
) -> Array[Array[Double]] {
  let n_cells = data.length()
  let n_selected = indices.length()

  let selected : Array[Array[Double]] = Array::new()
  let mut i = 0
  while i < n_cells {
    let row : Array[Double] = Array::new()
    let mut j = 0
    while j < n_selected {
      let idx = indices[j]
      if idx < data[i].length() {
        row.push(data[i][idx])
      } else {
        row.push(0.0)
      }
      j = j + 1
    }
    selected.push(row)
    i = i + 1
  }

  selected
}
