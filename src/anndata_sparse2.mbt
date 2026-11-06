/// anndata_sparse2.mbt — Sparse matrix format conversions and operations
/// Items 25-27: CSR/CSC/COO explicit conversion, eliminate_zeros, sort_indices

///|
/// CSR representation
pub struct AnnDataCSR {
  data : Array[Double]
  indices : Array[Int32]
  indptr : Array[Int32]
  n_rows : Int32
  n_cols : Int32
}

///|
/// CSC representation
pub struct AnnDataCSC {
  data : Array[Double]
  indices : Array[Int32]
  indptr : Array[Int32]
  n_rows : Int32
  n_cols : Int32
}

///|
/// COO representation
pub struct AnnDataCOO {
  data : Array[Double]
  row : Array[Int32]
  col : Array[Int32]
  n_rows : Int32
  n_cols : Int32
}

///|
// === COO construction from dense matrix ===

///|
/// Convert dense matrix to COO format.
pub fn anndata_dense_to_coo(matrix : Array[Array[Double]]) -> AnnDataCOO {
  let data : Array[Double] = Array(capacity=0)
  let row : Array[Int32] = Array(capacity=0)
  let col : Array[Int32] = Array(capacity=0)
  let n_rows = matrix.length()
  let n_cols = if n_rows > 0 { matrix[0].length() } else { 0 }
  for i = 0; i < n_rows; i = i + 1 {
    for j = 0; j < n_cols; j = j + 1 {
      if matrix[i][j] != 0.0 {
        data.push(matrix[i][j])
        row.push(i)
        col.push(j)
      }
    }
  }
  AnnDataCOO::{ data, row, col, n_rows, n_cols, }
}

///|
// === COO → CSR ===

///|
/// Convert COO to CSR format.
pub fn anndata_coo_to_csr(coo : AnnDataCOO) -> AnnDataCSR {
  let n_rows = coo.n_rows
  let nnz = coo.data.length()
  // Count per row
  let counts : Array[Int32] = Array(capacity=n_rows + 1)
  for _i = 0; _i < n_rows + 1; _i = _i + 1 {
    counts.push(0)
  }
  for r in coo.row {
    counts[r + 1] = counts[r + 1] + 1
  }
  // Cumulative sum → indptr
  let indptr : Array[Int32] = Array(capacity=n_rows + 1)
  indptr.push(0)
  let mut cum = 0
  for i = 0; i < n_rows; i = i + 1 {
    cum = cum + counts[i + 1]
    indptr.push(cum)
  }
  // Fill data and indices
  let data : Array[Double] = Array(capacity=nnz)
  let indices : Array[Int32] = Array(capacity=nnz)
  for _i = 0; _i < nnz; _i = _i + 1 {
    data.push(0.0)
    indices.push(0)
  }
  // Position tracker
  let pos : Array[Int32] = Array(capacity=n_rows)
  for i = 0; i < n_rows; i = i + 1 {
    pos.push(indptr[i])
  }
  for k = 0; k < nnz; k = k + 1 {
    let r = coo.row[k]
    let p = pos[r]
    data[p] = coo.data[k]
    indices[p] = coo.col[k]
    pos[r] = pos[r] + 1
  }
  AnnDataCSR::{ data, indices, indptr, n_rows: coo.n_rows, n_cols: coo.n_cols, }
}

///|
// === CSR → CSC ===

///|
/// Convert CSR to CSC format (transpose of CSR of transposed matrix).
pub fn anndata_csr_to_csc(csr : AnnDataCSR) -> AnnDataCSC {
  let n_rows = csr.n_rows
  let n_cols = csr.n_cols
  let nnz = csr.data.length()
  // Count per column
  let counts : Array[Int32] = Array(capacity=n_cols + 1)
  for _i = 0; _i < n_cols + 1; _i = _i + 1 {
    counts.push(0)
  }
  for idx in csr.indices {
    counts[idx + 1] = counts[idx + 1] + 1
  }
  // Cumulative sum → indptr
  let indptr : Array[Int32] = Array(capacity=n_cols + 1)
  indptr.push(0)
  let mut cum = 0
  for j = 0; j < n_cols; j = j + 1 {
    cum = cum + counts[j + 1]
    indptr.push(cum)
  }
  // Fill
  let data : Array[Double] = Array(capacity=nnz)
  let indices : Array[Int32] = Array(capacity=nnz)
  for _i = 0; _i < nnz; _i = _i + 1 {
    data.push(0.0)
    indices.push(0)
  }
  let pos : Array[Int32] = Array(capacity=n_cols)
  for j = 0; j < n_cols; j = j + 1 {
    pos.push(indptr[j])
  }
  for r = 0; r < n_rows; r = r + 1 {
    let start = csr.indptr[r]
    let end = csr.indptr[r + 1]
    for k = start; k < end; k = k + 1 {
      let c = csr.indices[k]
      let p = pos[c]
      data[p] = csr.data[k]
      indices[p] = r
      pos[c] = pos[c] + 1
    }
  }
  AnnDataCSC::{ data, indices, indptr, n_rows: csr.n_rows, n_cols: csr.n_cols, }
}

///|
// === CSR → dense ===

///|
/// Convert CSR back to dense matrix.
pub fn anndata_csr_to_dense(csr : AnnDataCSR) -> Array[Array[Double]] {
  let n_rows = csr.n_rows
  let n_cols = csr.n_cols
  let result : Array[Array[Double]] = Array(capacity=n_rows)
  for i = 0; i < n_rows; i = i + 1 {
    let row : Array[Double] = Array(capacity=n_cols)
    for _j = 0; _j < n_cols; _j = _j + 1 {
      row.push(0.0)
    }
    let start = csr.indptr[i]
    let end = csr.indptr[i + 1]
    for k = start; k < end; k = k + 1 {
      row[csr.indices[k]] = csr.data[k]
    }
    result.push(row)
  }
  result
}

///|
// === eliminate_zeros ===

///|
/// Remove explicit zeros from COO format.
pub fn anndata_coo_eliminate_zeros(coo : AnnDataCOO) -> AnnDataCOO {
  let data : Array[Double] = Array(capacity=0)
  let row : Array[Int32] = Array(capacity=0)
  let col : Array[Int32] = Array(capacity=0)
  for k = 0; k < coo.data.length(); k = k + 1 {
    if coo.data[k] != 0.0 {
      data.push(coo.data[k])
      row.push(coo.row[k])
      col.push(coo.col[k])
    }
  }
  AnnDataCOO::{ data, row, col, n_rows: coo.n_rows, n_cols: coo.n_cols, }
}

///|
/// Remove explicit zeros from CSR format.
pub fn anndata_csr_eliminate_zeros(csr : AnnDataCSR) -> AnnDataCSR {
  let n_rows = csr.n_rows
  let data : Array[Double] = Array(capacity=0)
  let indices : Array[Int32] = Array(capacity=0)
  let indptr : Array[Int32] = Array(capacity=n_rows + 1)
  indptr.push(0)
  for i = 0; i < n_rows; i = i + 1 {
    let start = csr.indptr[i]
    let end = csr.indptr[i + 1]
    for k = start; k < end; k = k + 1 {
      if csr.data[k] != 0.0 {
        data.push(csr.data[k])
        indices.push(csr.indices[k])
      }
    }
    indptr.push(data.length())
  }
  AnnDataCSR::{ data, indices, indptr, n_rows: csr.n_rows, n_cols: csr.n_cols, }
}

///|
// === sort_indices ===

///|
/// Sort CSR indices within each row (column indices in ascending order).
pub fn anndata_csr_sort_indices(csr : AnnDataCSR) -> AnnDataCSR {
  let n_rows = csr.n_rows
  let data : Array[Double] = Array(capacity=csr.data.length())
  let indices : Array[Int32] = Array(capacity=csr.indices.length())
  let indptr : Array[Int32] = Array(capacity=n_rows + 1)
  for i = 0; i < n_rows + 1; i = i + 1 {
    indptr.push(csr.indptr[i])
  }
  for i = 0; i < n_rows; i = i + 1 {
    let start = csr.indptr[i]
    let end = csr.indptr[i + 1]
    // Collect (col, data) pairs and sort by col
    let pairs : Array[(Int32, Double)] = Array(capacity=end - start)
    for k = start; k < end; k = k + 1 {
      pairs.push((csr.indices[k], csr.data[k]))
    }
    pairs.sort_by(fn(a, b) {
      if a.0 < b.0 {
        -1
      } else if a.0 > b.0 {
        1
      } else {
        0
      }
    })
    for p in pairs {
      indices.push(p.0)
      data.push(p.1)
    }
  }
  AnnDataCSR::{ data, indices, indptr, n_rows: csr.n_rows, n_cols: csr.n_cols, }
}

///|
/// Count non-zero elements in CSR.
pub fn anndata_csr_nnz(csr : AnnDataCSR) -> Int32 {
  csr.data.length()
}

///|
/// Density of CSR matrix.
pub fn anndata_csr_density(csr : AnnDataCSR) -> Double {
  let total = csr.n_rows * csr.n_cols
  if total == 0 {
    0.0
  } else {
    csr.data.length().to_double() / total.to_double()
  }
}

///|
pub fn anndata_sparse2_dummy() -> Unit {
  ignore(AnnDataCSR::{
    data: Array(capacity=0),
    indices: Array(capacity=0),
    indptr: Array(capacity=0),
    n_rows: 0,
    n_cols: 0,
  })
  ignore(AnnDataCSC::{
    data: Array(capacity=0),
    indices: Array(capacity=0),
    indptr: Array(capacity=0),
    n_rows: 0,
    n_cols: 0,
  })
  ignore(AnnDataCOO::{
    data: Array(capacity=0),
    row: Array(capacity=0),
    col: Array(capacity=0),
    n_rows: 0,
    n_cols: 0,
  })
}
