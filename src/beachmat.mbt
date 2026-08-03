///|
/// beachmat - Bioconductor API for matrix access.
///
/// Provides a C++/API-like interface for accessing matrix data from
/// SummarizedExperiment objects, supporting block processing,
/// random access, and different input matrix types.
///
/// Key concepts:
///   - BmatParam: configuration for matrix access (block size, storage mode)
///   - Bmat: interface for accessing a matrix (column/row-wise)
///   - Block processing: iterate over blocks of rows/columns
///   - Linear iterator: iterate element by element
///
/// All public identifiers are prefixed with `bmat_` to avoid conflicts.

// ============================================================================
// Data structures
// ============================================================================

///|
/// Storage modes for matrix data.
enum BmatStorageMode {
  Dense
  Sparse
  Delayed
}

///|
/// Parameters for matrix access.
pub struct BmatParam {
  row_block_size : Int
  col_block_size : Int
  cache_data : Bool
  storage_mode : BmatStorageMode
}

///|
/// Create default parameters.
pub fn BmatParam::new() -> BmatParam {
  { row_block_size: 0, col_block_size: 0, cache_data: true, storage_mode: BmatStorageMode::Dense }
}

///|
/// Create parameters with specific block sizes.
pub fn BmatParam::with_blocks(row_block : Int, col_block : Int) -> BmatParam {
  { row_block_size: row_block, col_block_size: col_block, cache_data: true, storage_mode: BmatStorageMode::Dense }
}

///|
/// Create parameters for column-wise access.
pub fn BmatParam::column_param(block_size : Int) -> BmatParam {
  { row_block_size: 0, col_block_size: block_size, cache_data: true, storage_mode: BmatStorageMode::Dense }
}

///|
/// Create parameters for row-wise access.
pub fn BmatParam::row_param(block_size : Int) -> BmatParam {
  { row_block_size: block_size, col_block_size: 0, cache_data: true, storage_mode: BmatStorageMode::Dense }
}

///|
/// Get row block size.
pub fn BmatParam::row_block_size(self : BmatParam) -> Int {
  self.row_block_size
}

///|
/// Get column block size.
pub fn BmatParam::col_block_size(self : BmatParam) -> Int {
  self.col_block_size
}

// ============================================================================
// Bmat - Matrix access interface
// ============================================================================

///|
/// A matrix accessor that provides block-wise and random access.
pub struct Bmat {
  dim_rows : Int
  dim_cols : Int
  data : Array[Double]
  param : BmatParam
}

///|
/// Create a Bmat from a 2D array (row-major).
pub fn Bmat::new(data : Array[Array[Double]], param : BmatParam) -> Bmat {
  let dim_rows = data.length()
  if dim_rows == 0 {
    return { dim_rows: 0, dim_cols: 0, data: Array::make(0, 0.0), param }
  }
  let dim_cols = data[0].length()
  let flat : Array[Double] = []
  for i in 0..<dim_rows {
    for j in 0..<dim_cols {
      flat.push(data[i][j])
    }
  }
  { dim_rows, dim_cols, data: flat, param }
}

///|
/// Create a zero-filled Bmat.
pub fn Bmat::zeros(dim_rows : Int, dim_cols : Int) -> Bmat {
  let data = Array::make(dim_rows * dim_cols, 0.0)
  { dim_rows, dim_cols, data, param: BmatParam::new() }
}

///|
/// Get number of rows.
pub fn Bmat::dim_rows(self : Bmat) -> Int {
  self.dim_rows
}

///|
/// Get number of columns.
pub fn Bmat::dim_cols(self : Bmat) -> Int {
  self.dim_cols
}

///|
/// Get element at (row, col).
pub fn Bmat::get(self : Bmat, row : Int, col : Int) -> Double {
  self.data[row * self.dim_cols + col]
}

///|
/// Set element at (row, col).
pub fn Bmat::set(self : Bmat, row : Int, col : Int, val : Double) -> Unit {
  self.data[row * self.dim_cols + col] = val
}

///|
/// Get a column as an array.
pub fn Bmat::get_col(self : Bmat, col : Int) -> Array[Double] {
  let result : Array[Double] = []
  for i in 0..<self.dim_rows {
    result.push(self.data[i * self.dim_cols + col])
  }
  result
}

///|
/// Get a row as an array.
pub fn Bmat::get_row(self : Bmat, row : Int) -> Array[Double] {
  let result : Array[Double] = []
  for j in 0..<self.dim_cols {
    result.push(self.data[row * self.dim_cols + j])
  }
  result
}

// ============================================================================
// Block iteration
// ============================================================================

///|
/// A block of data from a Bmat.
pub struct BmatBlock {
  row_start : Int
  row_end : Int
  col_start : Int
  col_end : Int
  block_data : Array[Array[Double]]
}

///|
/// Get start row.
pub fn BmatBlock::row_start(self : BmatBlock) -> Int {
  self.row_start
}

///|
/// Get end row.
pub fn BmatBlock::row_end(self : BmatBlock) -> Int {
  self.row_end
}

///|
/// Get start col.
pub fn BmatBlock::col_start(self : BmatBlock) -> Int {
  self.col_start
}

///|
/// Get end col.
pub fn BmatBlock::col_end(self : BmatBlock) -> Int {
  self.col_end
}

///|
/// Get block data as 2D array.
pub fn BmatBlock::get_data(self : BmatBlock) -> Array[Array[Double]] {
  self.block_data
}

///|
/// Number of rows in the block.
pub fn BmatBlock::block_nrow(self : BmatBlock) -> Int {
  self.row_end - self.row_start
}

///|
/// Number of columns in the block.
pub fn BmatBlock::block_ncol(self : BmatBlock) -> Int {
  self.col_end - self.col_start
}

///|
/// Iterate over column blocks and apply a function to each block.
pub fn bmat_apply_col_blocks(
  bmat : Bmat,
  col_block_size : Int,
  apply_fn : (BmatBlock) -> Unit,
) -> Unit {
  let n_cols = bmat.dim_cols
  let mut col = 0
  while col < n_cols {
    let end = if col + col_block_size < n_cols { col + col_block_size } else { n_cols }
    let block_data : Array[Array[Double]] = []
    let mut i = 0
    while i < bmat.dim_rows {
      let row_data : Array[Double] = []
      let mut j = col
      while j < end {
        row_data.push(bmat.get(i, j))
        j = j + 1
      }
      block_data.push(row_data)
      i = i + 1
    }
    let block = BmatBlock::{ row_start: 0, row_end: bmat.dim_rows, col_start: col, col_end: end, block_data }
    apply_fn(block)
    col = end
  }
}

///|
/// Iterate over row blocks and apply a function to each block.
pub fn bmat_apply_row_blocks(
  bmat : Bmat,
  row_block_size : Int,
  apply_fn : (BmatBlock) -> Unit,
) -> Unit {
  let n_rows = bmat.dim_rows
  let mut row = 0
  while row < n_rows {
    let end = if row + row_block_size < n_rows { row + row_block_size } else { n_rows }
    let block_data : Array[Array[Double]] = []
    let mut i = row
    while i < end {
      let row_data : Array[Double] = []
      let mut j = 0
      while j < bmat.dim_cols {
        row_data.push(bmat.get(i, j))
        j = j + 1
      }
      block_data.push(row_data)
      i = i + 1
    }
    let block = BmatBlock::{ row_start: row, row_end: end, col_start: 0, col_end: bmat.dim_cols, block_data }
    apply_fn(block)
    row = end
  }
}

///|
/// Iterate over the full matrix and apply a function.
pub fn bmat_foreach(
  bmat : Bmat,
  iter_fn : (Int, Int, Double) -> Unit,
) -> Unit {
  let mut i = 0
  while i < bmat.dim_rows {
    let mut j = 0
    while j < bmat.dim_cols {
      iter_fn(i, j, bmat.get(i, j))
      j = j + 1
    }
    i = i + 1
  }
}

// ============================================================================
// Linear iterator
// ============================================================================

///|
/// A linear iterator over all elements in a Bmat.
pub struct BmatIterator {
  data : Array[Double]
  total_cols : Int
  total : Int
  mut pos : Int
}

///|
/// Create a new iterator.
pub fn BmatIterator::new(bmat : Bmat) -> BmatIterator {
  { data: bmat.data, total_cols: bmat.dim_cols, total: bmat.dim_rows * bmat.dim_cols, pos: 0 }
}

///|
/// Check if the iterator has more elements.
pub fn BmatIterator::has_next(self : BmatIterator) -> Bool {
  self.pos < self.total
}

///|
/// Get the next value and advance.
pub fn BmatIterator::next(self : BmatIterator) -> Double? {
  if self.pos >= self.total {
    return None
  }
  let val = self.data[self.pos]
  self.pos = self.pos + 1
  Some(val)
}

///|
/// Get current row index.
pub fn BmatIterator::cur_row(self : BmatIterator) -> Int {
  if self.total_cols == 0 { 0 } else { self.pos / self.total_cols }
}

///|
/// Get current column index.
pub fn BmatIterator::cur_col(self : BmatIterator) -> Int {
  if self.total_cols == 0 { 0 } else { self.pos % self.total_cols }
}

// ============================================================================
// Utility functions
// ============================================================================

///|
/// Create a Bmat from a simple 2D array.
pub fn bmat_from_array(data : Array[Array[Double]]) -> Bmat {
  Bmat::new(data, BmatParam::new())
}

///|
/// Convert a Bmat to a 2D array.
pub fn bmat_to_array(bmat : Bmat) -> Array[Array[Double]] {
  let result : Array[Array[Double]] = []
  for i in 0..<bmat.dim_rows {
    let row : Array[Double] = []
    for j in 0..<bmat.dim_cols {
      row.push(bmat.get(i, j))
    }
    result.push(row)
  }
  result
}

///|
/// Subset a Bmat to a specific range.
pub fn bmat_subset(
  bmat : Bmat,
  row_start : Int,
  row_end : Int,
  col_start : Int,
  col_end : Int,
) -> Bmat {
  let new_data : Array[Array[Double]] = []
  let mut i = row_start
  while i < row_end {
    let row : Array[Double] = []
    let mut j = col_start
    while j < col_end {
      row.push(bmat.get(i, j))
      j = j + 1
    }
    new_data.push(row)
    i = i + 1
  }
  Bmat::new(new_data, BmatParam::new())
}

///|
/// Transpose a Bmat.
pub fn bmat_transpose(bmat : Bmat) -> Bmat {
  let new_data : Array[Array[Double]] = []
  let mut j = 0
  while j < bmat.dim_cols {
    let row : Array[Double] = []
    let mut i = 0
    while i < bmat.dim_rows {
      row.push(bmat.get(i, j))
      i = i + 1
    }
    new_data.push(row)
    j = j + 1
  }
  Bmat::new(new_data, BmatParam::new())
}

///|
/// Bind columns (same row count).
pub fn bmat_bind_cols(left : Bmat, right : Bmat) -> Bmat {
  let n_rows = left.dim_rows
  let ncol_left = left.dim_cols
  let ncol_right = right.dim_cols
  let new_data : Array[Array[Double]] = []
  let mut i = 0
  while i < n_rows {
    let row : Array[Double] = []
    let mut j = 0
    while j < ncol_left {
      row.push(left.get(i, j))
      j = j + 1
    }
    let mut k = 0
    while k < ncol_right {
      row.push(right.get(i, k))
      k = k + 1
    }
    new_data.push(row)
    i = i + 1
  }
  Bmat::new(new_data, BmatParam::new())
}

///|
/// Bind rows (same col count).
pub fn bmat_bind_rows(top : Bmat, bottom : Bmat) -> Bmat {
  let n_cols = top.dim_cols
  let new_data : Array[Array[Double]] = []
  let mut i = 0
  while i < top.dim_rows {
    let row : Array[Double] = []
    let mut j = 0
    while j < n_cols {
      row.push(top.get(i, j))
      j = j + 1
    }
    new_data.push(row)
    i = i + 1
  }
  let mut k = 0
  while k < bottom.dim_rows {
    let row : Array[Double] = []
    let mut j = 0
    while j < n_cols {
      row.push(bottom.get(k, j))
      j = j + 1
    }
    new_data.push(row)
    k = k + 1
  }
  Bmat::new(new_data, BmatParam::new())
}

///|
/// Apply a function to each element.
pub fn bmat_apply_elementwise(
  bmat : Bmat,
  map_fn : (Double) -> Double,
) -> Bmat {
  let new_data : Array[Array[Double]] = []
  let mut i = 0
  while i < bmat.dim_rows {
    let row : Array[Double] = []
    let mut j = 0
    while j < bmat.dim_cols {
      row.push(map_fn(bmat.get(i, j)))
      j = j + 1
    }
    new_data.push(row)
    i = i + 1
  }
  Bmat::new(new_data, BmatParam::new())
}

///|
/// Pretty-print a Bmat.
pub fn Bmat::to_string(self : Bmat) -> String {
  let mut s = "Bmat (" + self.dim_rows.to_string() + " x " + self.dim_cols.to_string() + "):\n"
  let mut i = 0
  while i < self.dim_rows {
    let mut j = 0
    while j < self.dim_cols {
      let val = self.get(i, j)
      s = s + val.to_string()
      if j < self.dim_cols - 1 {
        s = s + "\t"
      }
      j = j + 1
    }
    if i < self.dim_rows - 1 {
      s = s + "\n"
    }
    i = i + 1
  }
  s
}
