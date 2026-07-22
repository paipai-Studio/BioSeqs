///|
/// Matrix - Sparse and dense matrix operations for BioSeqs.
/// Inspired by Bioconductor Matrix package.
/// Provides matrix operations for genomic data analysis.

///|
pub enum BiocMatrixType {
  BiocDenseMatrix
  BiocSparseMatrixCSC
  BiocSparseMatrixCSR
  BiocDiagonalMatrix
}

///|
pub struct BiocMatrix {
  matrix_type : BiocMatrixType
  nrow : Int
  ncol : Int
  data : Array[Double]
  row_indices : Array[Int]
  col_indices : Array[Int]
  ptr : Array[Int]
}

///|
pub fn BiocMatrix::new(matrix_type : BiocMatrixType, nrow : Int, ncol : Int) -> BiocMatrix {
  BiocMatrix::{
    matrix_type,
    nrow,
    ncol,
    data: Array::new(),
    row_indices: Array::new(),
    col_indices: Array::new(),
    ptr: Array::new()
  }
}

///|
pub fn BiocMatrix::dense_from_array(data : Array[Array[Double]]) -> BiocMatrix {
  let nrow = data.length()
  let ncol = if nrow > 0 { data[0].length() } else { 0 }
  
  let flat_data : Array[Double] = Array::new()
  let mut i = 0
  while i < nrow {
    let mut j = 0
    while j < ncol {
      flat_data.push(data[i][j])
      j = j + 1
    }
    i = i + 1
  }
  
  BiocMatrix::{
    matrix_type: BiocDenseMatrix,
    nrow,
    ncol,
    data: flat_data,
    row_indices: Array::new(),
    col_indices: Array::new(),
    ptr: Array::new()
  }
}

///|
pub fn BiocMatrix::csc_from_triplets(row_idx : Array[Int], col_idx : Array[Int], values : Array[Double], nrow : Int, ncol : Int) -> BiocMatrix {
  let sorted_indices = sort_triplets(row_idx, col_idx, values)
  let sorted_rows = sorted_indices.0
  let sorted_cols = sorted_indices.1
  let sorted_values = sorted_indices.2
  
  let ptr : Array[Int] = Array::new()
  ptr.push(0)
  
  let mut current_col = 0
  let mut count = 0
  let mut i = 0
  
  while i < sorted_cols.length() {
    if sorted_cols[i] == current_col {
      count = count + 1
    } else {
      ptr.push(ptr[ptr.length() - 1] + count)
      current_col = sorted_cols[i]
      count = 1
    }
    i = i + 1
  }
  ptr.push(ptr[ptr.length() - 1] + count)
  
  BiocMatrix::{
    matrix_type: BiocSparseMatrixCSC,
    nrow,
    ncol,
    data: sorted_values,
    row_indices: sorted_rows,
    col_indices: sorted_cols,
    ptr: ptr
  }
}

///|
pub fn BiocMatrix::csr_from_triplets(row_idx : Array[Int], col_idx : Array[Int], values : Array[Double], nrow : Int, ncol : Int) -> BiocMatrix {
  let sorted_indices = sort_triplets_by_row(row_idx, col_idx, values)
  let sorted_rows = sorted_indices.0
  let sorted_cols = sorted_indices.1
  let sorted_values = sorted_indices.2
  
  let ptr : Array[Int] = Array::new()
  ptr.push(0)
  
  let mut current_row = 0
  let mut count = 0
  let mut i = 0
  
  while i < sorted_rows.length() {
    if sorted_rows[i] == current_row {
      count = count + 1
    } else {
      ptr.push(ptr[ptr.length() - 1] + count)
      current_row = sorted_rows[i]
      count = 1
    }
    i = i + 1
  }
  ptr.push(ptr[ptr.length() - 1] + count)
  
  BiocMatrix::{
    matrix_type: BiocSparseMatrixCSR,
    nrow,
    ncol,
    data: sorted_values,
    row_indices: sorted_rows,
    col_indices: sorted_cols,
    ptr: ptr
  }
}

///|
fn sort_triplets(row_idx : Array[Int], col_idx : Array[Int], values : Array[Double]) -> (Array[Int], Array[Int], Array[Double]) {
  let n = row_idx.length()
  let indices : Array[Int] = Array::new()
  let mut i = 0
  while i < n {
    indices.push(i)
    i = i + 1
  }
  
  let mut sorted = false
  while !sorted {
    sorted = true
    i = 0
    while i < n - 1 {
      let j = i + 1
      if col_idx[indices[i]] > col_idx[indices[j]] || 
         (col_idx[indices[i]] == col_idx[indices[j]] && row_idx[indices[i]] > row_idx[indices[j]]) {
        let temp = indices[i]
        indices[i] = indices[j]
        indices[j] = temp
        sorted = false
      }
      i = i + 1
    }
  }
  
  let sorted_rows : Array[Int] = Array::new()
  let sorted_cols : Array[Int] = Array::new()
  let sorted_values : Array[Double] = Array::new()
  
  i = 0
  while i < n {
    sorted_rows.push(row_idx[indices[i]])
    sorted_cols.push(col_idx[indices[i]])
    sorted_values.push(values[indices[i]])
    i = i + 1
  }
  
  (sorted_rows, sorted_cols, sorted_values)
}

///|
fn sort_triplets_by_row(row_idx : Array[Int], col_idx : Array[Int], values : Array[Double]) -> (Array[Int], Array[Int], Array[Double]) {
  let n = row_idx.length()
  let indices : Array[Int] = Array::new()
  let mut i = 0
  while i < n {
    indices.push(i)
    i = i + 1
  }
  
  let mut sorted = false
  while !sorted {
    sorted = true
    i = 0
    while i < n - 1 {
      let j = i + 1
      if row_idx[indices[i]] > row_idx[indices[j]] || 
         (row_idx[indices[i]] == row_idx[indices[j]] && col_idx[indices[i]] > col_idx[indices[j]]) {
        let temp = indices[i]
        indices[i] = indices[j]
        indices[j] = temp
        sorted = false
      }
      i = i + 1
    }
  }
  
  let sorted_rows : Array[Int] = Array::new()
  let sorted_cols : Array[Int] = Array::new()
  let sorted_values : Array[Double] = Array::new()
  
  i = 0
  while i < n {
    sorted_rows.push(row_idx[indices[i]])
    sorted_cols.push(col_idx[indices[i]])
    sorted_values.push(values[indices[i]])
    i = i + 1
  }
  
  (sorted_rows, sorted_cols, sorted_values)
}

///|
pub fn BiocMatrix::get(self : BiocMatrix, i : Int, j : Int) -> Double {
  if i < 0 || i >= self.nrow || j < 0 || j >= self.ncol {
    return 0.0
  }
  
  match self.matrix_type {
    BiocDenseMatrix => {
      let idx = i * self.ncol + j
      if idx < self.data.length() {
        self.data[idx]
      } else {
        0.0
      }
    }
    BiocSparseMatrixCSC => {
      let start = self.ptr[j]
      let end = self.ptr[j + 1]
      let mut k = start
      while k < end {
        if self.row_indices[k] == i {
          return self.data[k]
        }
        k = k + 1
      }
      0.0
    }
    BiocSparseMatrixCSR => {
      let start = self.ptr[i]
      let end = self.ptr[i + 1]
      let mut k = start
      while k < end {
        if self.col_indices[k] == j {
          return self.data[k]
        }
        k = k + 1
      }
      0.0
    }
    BiocDiagonalMatrix => {
      if i == j && i < self.data.length() {
        self.data[i]
      } else {
        0.0
      }
    }
  }
}

///|
pub fn BiocMatrix::set(self : BiocMatrix, i : Int, j : Int, value : Double) -> BiocMatrix {
  if i < 0 || i >= self.nrow || j < 0 || j >= self.ncol {
    return self
  }
  
  match self.matrix_type {
    BiocDenseMatrix => {
      let idx = i * self.ncol + j
      let new_data = self.data.copy()
      if idx < new_data.length() {
        new_data[idx] = value
      }
      BiocMatrix::{
        matrix_type: BiocDenseMatrix,
        nrow: self.nrow,
        ncol: self.ncol,
        data: new_data,
        row_indices: self.row_indices,
        col_indices: self.col_indices,
        ptr: self.ptr
      }
    }
    _ => self
  }
}

///|
pub fn BiocMatrix::transpose(self : BiocMatrix) -> BiocMatrix {
  match self.matrix_type {
    BiocDenseMatrix => {
      let result_data : Array[Array[Double]] = Array::new()
      let mut j = 0
      while j < self.ncol {
        let row : Array[Double] = Array::new()
        let mut i = 0
        while i < self.nrow {
          row.push(self.get(i, j))
          i = i + 1
        }
        result_data.push(row)
        j = j + 1
      }
      BiocMatrix::dense_from_array(result_data)
    }
    BiocSparseMatrixCSC => {
      BiocMatrix::csr_from_triplets(self.col_indices, self.row_indices, self.data, self.ncol, self.nrow)
    }
    BiocSparseMatrixCSR => {
      BiocMatrix::csc_from_triplets(self.col_indices, self.row_indices, self.data, self.ncol, self.nrow)
    }
    BiocDiagonalMatrix => self
  }
}

///|
pub fn BiocMatrix::add(self : BiocMatrix, other : BiocMatrix) -> BiocMatrix? {
  if self.nrow != other.nrow || self.ncol != other.ncol {
    return None
  }
  
  let result_data : Array[Array[Double]] = Array::new()
  let mut i = 0
  while i < self.nrow {
    let row : Array[Double] = Array::new()
    let mut j = 0
    while j < self.ncol {
      row.push(self.get(i, j) + other.get(i, j))
      j = j + 1
    }
    result_data.push(row)
    i = i + 1
  }
  
  Some(BiocMatrix::dense_from_array(result_data))
}

///|
pub fn BiocMatrix::multiply(self : BiocMatrix, other : BiocMatrix) -> BiocMatrix? {
  if self.ncol != other.nrow {
    return None
  }
  
  let result_data : Array[Array[Double]] = Array::new()
  let mut i = 0
  while i < self.nrow {
    let row : Array[Double] = Array::new()
    let mut j = 0
    while j < other.ncol {
      let mut sum = 0.0
      let mut k = 0
      while k < self.ncol {
        sum = sum + self.get(i, k) * other.get(k, j)
        k = k + 1
      }
      row.push(sum)
      j = j + 1
    }
    result_data.push(row)
    i = i + 1
  }
  
  Some(BiocMatrix::dense_from_array(result_data))
}

///|
pub fn BiocMatrix::row_sums(self : BiocMatrix) -> Array[Double] {
  let result : Array[Double] = Array::new()
  let mut i = 0
  while i < self.nrow {
    let mut sum = 0.0
    let mut j = 0
    while j < self.ncol {
      sum = sum + self.get(i, j)
      j = j + 1
    }
    result.push(sum)
    i = i + 1
  }
  result
}

///|
pub fn BiocMatrix::col_sums(self : BiocMatrix) -> Array[Double] {
  let result : Array[Double] = Array::new()
  let mut j = 0
  while j < self.ncol {
    let mut sum = 0.0
    let mut i = 0
    while i < self.nrow {
      sum = sum + self.get(i, j)
      i = i + 1
    }
    result.push(sum)
    j = j + 1
  }
  result
}

///|
pub fn BiocMatrix::row_means(self : BiocMatrix) -> Array[Double] {
  let sums = self.row_sums()
  let result : Array[Double] = Array::new()
  let mut i = 0
  while i < sums.length() {
    result.push(sums[i] / self.ncol.to_double())
    i = i + 1
  }
  result
}

///|
pub fn BiocMatrix::col_means(self : BiocMatrix) -> Array[Double] {
  let sums = self.col_sums()
  let result : Array[Double] = Array::new()
  let mut i = 0
  while i < sums.length() {
    result.push(sums[i] / self.nrow.to_double())
    i = i + 1
  }
  result
}

///|
pub fn BiocMatrix::norm(self : BiocMatrix, p : Double) -> Double {
  let mut sum = 0.0
  
  match self.matrix_type {
    BiocDenseMatrix => {
      let mut i = 0
      while i < self.data.length() {
        sum = sum + @math.pow(@double.abs(self.data[i]), p)
        i = i + 1
      }
    }
    BiocSparseMatrixCSC | BiocSparseMatrixCSR => {
      let mut i = 0
      while i < self.data.length() {
        sum = sum + @math.pow(@double.abs(self.data[i]), p)
        i = i + 1
      }
    }
    BiocDiagonalMatrix => {
      let mut i = 0
      while i < self.data.length() {
        sum = sum + @math.pow(@double.abs(self.data[i]), p)
        i = i + 1
      }
    }
  }
  
  @math.pow(sum, 1.0 / p)
}

///|
pub fn BiocMatrix::nrow(self : BiocMatrix) -> Int {
  self.nrow
}

///|
pub fn BiocMatrix::ncol(self : BiocMatrix) -> Int {
  self.ncol
}

///|
pub fn BiocMatrix::dim(self : BiocMatrix) -> (Int, Int) {
  (self.nrow, self.ncol)
}

///|
pub fn BiocMatrix::nnz(self : BiocMatrix) -> Int {
  self.data.length()
}

///|
pub fn create_example_dense_matrix() -> BiocMatrix {
  BiocMatrix::dense_from_array([
    [1.0, 2.0, 3.0],
    [4.0, 5.0, 6.0],
    [7.0, 8.0, 9.0]
  ])
}

///|
pub fn create_example_sparse_matrix() -> BiocMatrix {
  BiocMatrix::csc_from_triplets(
    [0, 1, 2, 0, 1, 2],
    [0, 0, 0, 1, 1, 1],
    [1.0, 4.0, 7.0, 2.0, 5.0, 8.0],
    3,
    2
  )
}
