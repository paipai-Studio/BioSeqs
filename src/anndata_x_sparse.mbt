/// anndata_x_sparse.mbt — first-class CSR/CSC access to the AnnData X matrix.
///
/// `AnnData.x_sparse` stores a CSR-like SparseArray; this module bridges it
/// with the explicit `AnnDataCSR` / `AnnDataCSC` formats (anndata_sparse2),
/// mirroring `scipy.sparse` interop in Python anndata: construct an AnnData
/// straight from CSR/CSC, request X as CSR/CSC, slice CSR along rows/cols
/// while staying sparse, and swap X inplace.

///|
/// Assemble a CSR matrix from its raw parts (data / indices / indptr /
/// shape). This is the public constructor for `AnnDataCSR`; validation
/// happens later in `anndata_csr_to_sparse`, so malformed parts are
/// representable (and detected on conversion).
pub fn anndata_csr_make(
  data : Array[Double],
  indices : Array[Int32],
  indptr : Array[Int32],
  n_rows : Int32,
  n_cols : Int32,
) -> AnnDataCSR {
  AnnDataCSR::{ data, indices, indptr, n_rows, n_cols, }
}

///|
/// Convert a dense matrix to CSR format.
pub fn anndata_dense_to_csr(matrix : Array[Array[Double]]) -> AnnDataCSR {
  anndata_coo_to_csr(anndata_dense_to_coo(matrix))
}

///|
/// Convert a dense matrix to CSC format.
pub fn anndata_dense_to_csc(matrix : Array[Array[Double]]) -> AnnDataCSC {
  anndata_csr_to_csc(anndata_dense_to_csr(matrix))
}

///|
/// Convert CSR to CSC (inverse of `anndata_csr_to_csc` via transpose rule:
/// CSC of A is the CSR of A transposed; here computed directly).
pub fn anndata_csc_to_csr(csc : AnnDataCSC) -> AnnDataCSR {
  let n_rows = csc.n_rows
  let n_cols = csc.n_cols
  let nnz = csc.data.length()
  let counts : Array[Int32] = Array(capacity=n_rows + 1)
  for _i = 0; _i < n_rows + 1; _i = _i + 1 {
    counts.push(0)
  }
  for idx in csc.indices {
    counts[idx + 1] = counts[idx + 1] + 1
  }
  let indptr : Array[Int32] = Array(capacity=n_rows + 1)
  indptr.push(0)
  let mut cum = 0
  for i = 0; i < n_rows; i = i + 1 {
    cum = cum + counts[i + 1]
    indptr.push(cum)
  }
  let data : Array[Double] = Array(capacity=nnz)
  let indices : Array[Int32] = Array(capacity=nnz)
  for _i = 0; _i < nnz; _i = _i + 1 {
    data.push(0.0)
    indices.push(0)
  }
  let pos : Array[Int32] = Array(capacity=n_rows)
  for i = 0; i < n_rows; i = i + 1 {
    pos.push(indptr[i])
  }
  for c = 0; c < n_cols; c = c + 1 {
    let start = csc.indptr[c]
    let end_ = csc.indptr[c + 1]
    for k = start; k < end_; k = k + 1 {
      let r = csc.indices[k]
      let p = pos[r]
      data[p] = csc.data[k]
      indices[p] = c
      pos[r] = pos[r] + 1
    }
  }
  AnnDataCSR::{ data, indices, indptr, n_rows: csc.n_rows, n_cols: csc.n_cols, }
}

///|
/// Convert CSC back to a dense matrix.
pub fn anndata_csc_to_dense(csc : AnnDataCSC) -> Array[Array[Double]] {
  let n_rows = csc.n_rows
  let n_cols = csc.n_cols
  let out : Array[Array[Double]] = Array(capacity=n_rows)
  for _i = 0; _i < n_rows; _i = _i + 1 {
    let row : Array[Double] = Array(capacity=n_cols)
    for _j = 0; _j < n_cols; _j = _j + 1 {
      row.push(0.0)
    }
    out.push(row)
  }
  for j = 0; j < n_cols; j = j + 1 {
    let start = csc.indptr[j]
    let end_ = csc.indptr[j + 1]
    for k = start; k < end_; k = k + 1 {
      out[csc.indices[k]][j] = csc.data[k]
    }
  }
  out
}

///|
/// Build the SparseArray backing for a CSR matrix (sparsity preserved, no
/// dense round-trip). Returns `None` when the CSR is malformed (indptr
/// length mismatch, data/indices length mismatch, out-of-range indices).
pub fn anndata_csr_to_sparse(csr : AnnDataCSR) -> SparseArray? {
  if csr.indptr.length() != csr.n_rows + 1 {
    return None
  }
  if csr.indices.length() != csr.data.length() {
    return None
  }
  for idx in csr.indices {
    if idx < 0 || idx >= csr.n_cols {
      return None
    }
  }
  let coords : Array[Array[Int32]] = Array(capacity=csr.data.length())
  for r = 0; r < csr.n_rows; r = r + 1 {
    let start = csr.indptr[r]
    let end_ = csr.indptr[r + 1]
    if start < 0 || end_ < start || end_ > csr.data.length() {
      return None
    }
    for k = start; k < end_; k = k + 1 {
      let c : Array[Int32] = Array(capacity=2)
      c.push(r)
      c.push(csr.indices[k])
      coords.push(c)
    }
  }
  let arr = SparseArray::from_coo([csr.n_rows, csr.n_cols], coords, csr.data) catch {
    _ => return None
  }
  Some(arr)
}

///|
/// Construct an AnnData directly from a CSR matrix (X stays sparse).
/// Returns `None` when the CSR is malformed.
pub fn AnnData::new_csr(
  csr : AnnDataCSR,
  obs_names : Array[String],
  var_names : Array[String],
) -> AnnData? {
  match anndata_csr_to_sparse(csr) {
    Some(s) => Some(AnnData::new_sparse(s, obs_names, var_names))
    None => None
  }
}

///|
/// Construct an AnnData directly from a CSC matrix (X stays sparse).
pub fn AnnData::new_csc(
  csc : AnnDataCSC,
  obs_names : Array[String],
  var_names : Array[String],
) -> AnnData? {
  match anndata_csr_to_sparse(anndata_csc_to_csr(csc)) {
    Some(s) => Some(AnnData::new_sparse(s, obs_names, var_names))
    None => None
  }
}

///|
/// Request X as an explicit CSR matrix (sparse storage is converted without
/// a dense round-trip; dense storage is converted via COO).
pub fn anndata_x_as_csr(adata : AnnData) -> AnnDataCSR {
  match adata.x_sparse {
    Some(s) => {
      let dims = s.dim()
      let coords = s.nzcoordinates()
      let row : Array[Int32] = Array(capacity=coords.length())
      let col : Array[Int32] = Array(capacity=coords.length())
      for c in coords {
        row.push(c[0])
        col.push(c[1])
      }
      let coo = AnnDataCOO::{
        data: s.nzvalues(),
        row,
        col,
        n_rows: dims[0],
        n_cols: dims[1],
      }
      anndata_coo_to_csr(coo)
    }
    None => anndata_dense_to_csr(adata.x)
  }
}

///|
/// Request X as an explicit CSC matrix.
pub fn anndata_x_as_csc(adata : AnnData) -> AnnDataCSC {
  anndata_csr_to_csc(anndata_x_as_csr(adata))
}

///|
/// Replace X inplace with a CSR-backed sparse storage. Returns `self`;
/// when the CSR is malformed the original X is kept unchanged.
pub fn AnnData::set_x_csr(self : AnnData, csr : AnnDataCSR) -> AnnData {
  match anndata_csr_to_sparse(csr) {
    Some(s) => {
      self.x = []
      self.x_sparse = Some(s)
      self
    }
    None => self
  }
}

///|
/// Replace X inplace with a CSC-backed sparse storage. Returns `self`.
pub fn AnnData::set_x_csc(self : AnnData, csc : AnnDataCSC) -> AnnData {
  self.set_x_csr(anndata_csc_to_csr(csc))
}

///|
/// Storage layout of X: "sparse" or "dense".
pub fn AnnData::x_layout(self : AnnData) -> String {
  match self.x_sparse {
    Some(_) => "sparse"
    None => "dense"
  }
}

///|
/// Get a single CSR row as a dense vector.
pub fn anndata_csr_get_row(csr : AnnDataCSR, i : Int32) -> Array[Double] {
  let row : Array[Double] = Array(capacity=csr.n_cols)
  for _j = 0; _j < csr.n_cols; _j = _j + 1 {
    row.push(0.0)
  }
  if i >= 0 && i < csr.n_rows {
    let start = csr.indptr[i]
    let end_ = csr.indptr[i + 1]
    for k = start; k < end_; k = k + 1 {
      row[csr.indices[k]] = csr.data[k]
    }
  }
  row
}

///|
/// Slice CSR rows while staying sparse (`csr[row_idx]` in scipy).
pub fn anndata_csr_slice_rows(
  csr : AnnDataCSR,
  idx : Array[Int32],
) -> AnnDataCSR {
  let data : Array[Double] = Array(capacity=0)
  let indices : Array[Int32] = Array(capacity=0)
  let indptr : Array[Int32] = Array(capacity=idx.length() + 1)
  indptr.push(0)
  for i in idx {
    if i >= 0 && i < csr.n_rows {
      let start = csr.indptr[i]
      let end_ = csr.indptr[i + 1]
      for k = start; k < end_; k = k + 1 {
        data.push(csr.data[k])
        indices.push(csr.indices[k])
      }
    }
    indptr.push(data.length())
  }
  AnnDataCSR::{
    data,
    indices,
    indptr,
    n_rows: idx.length(),
    n_cols: csr.n_cols,
  }
}

///|
/// Slice CSR columns while staying sparse (`csr[:, col_idx]` in scipy).
pub fn anndata_csr_slice_cols(
  csr : AnnDataCSR,
  idx : Array[Int32],
) -> AnnDataCSR {
  // Map old column -> new position
  let remap : Map[Int32, Int32] = Map([], capacity=idx.length())
  for pos = 0; pos < idx.length(); pos = pos + 1 {
    remap[idx[pos]] = pos
  }
  let data : Array[Double] = Array(capacity=0)
  let indices : Array[Int32] = Array(capacity=0)
  let indptr : Array[Int32] = Array(capacity=csr.n_rows + 1)
  indptr.push(0)
  for i = 0; i < csr.n_rows; i = i + 1 {
    let start = csr.indptr[i]
    let end_ = csr.indptr[i + 1]
    // collect (new_col, value) pairs for this row, keep ascending col order
    let pairs : Array[(Int32, Double)] = Array(capacity=0)
    for k = start; k < end_; k = k + 1 {
      match remap.get(csr.indices[k]) {
        Some(new_col) => pairs.push((new_col, csr.data[k]))
        None => ()
      }
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
    indptr.push(data.length())
  }
  AnnDataCSR::{
    data,
    indices,
    indptr,
    n_rows: csr.n_rows,
    n_cols: idx.length(),
  }
}
