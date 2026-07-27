///|
/// DelayedMatrixStats - Efficient row/column statistics on DelayedArray matrices.
/// Inspired by Bioconductor's DelayedMatrixStats package.
/// Provides memory-efficient statistical computations for large genomic matrices.

///|
pub struct DelayedMatrix {
  data : Array[Array[Double]]
  nrow : Int
  ncol : Int
  row_names : Array[String]
  col_names : Array[String]
}

///|
pub struct MatrixStatsResult {
  values : Array[Double]
  stat_type : String
  dimension : String
  row_names : Array[String]
  col_names : Array[String]
}

///|
pub fn MatrixStatsResult::new(
  values : Array[Double],
  stat_type : String,
  dimension : String,
  row_names : Array[String],
  col_names : Array[String],
) -> MatrixStatsResult {
  MatrixStatsResult::{ values, stat_type, dimension, row_names, col_names }
}

///|
pub fn DelayedMatrix::nrow(self : DelayedMatrix) -> Int {
  self.nrow
}

///|
pub fn DelayedMatrix::ncol(self : DelayedMatrix) -> Int {
  self.ncol
}

///|
pub fn DelayedMatrix::dim(self : DelayedMatrix) -> Array[Int] {
  [self.nrow, self.ncol]
}

///|
pub fn DelayedMatrix::get(self : DelayedMatrix, i : Int, j : Int) -> Double {
  self.data[i][j]
}

///|
pub fn create_delayed_matrix(
  data : Array[Array[Double]],
  row_names : Array[String],
  col_names : Array[String],
) -> DelayedMatrix {
  let nrow = data.length()
  let ncol = if nrow > 0 { data[0].length() } else { 0 }
  DelayedMatrix::{ data, nrow, ncol, row_names, col_names }
}

///|
pub fn delayed_matrix_subset(
  matrix : DelayedMatrix,
  rows : Array[Int],
  cols : Array[Int],
) -> DelayedMatrix {
  let new_data : Array[Array[Double]] = Array::new()
  let new_row_names : Array[String] = Array::new()
  let mut i = 0
  while i < rows.length() {
    let row_data : Array[Double] = Array::new()
    let mut j = 0
    while j < cols.length() {
      row_data.push(matrix.data[rows[i]][cols[j]])
      j = j + 1
    }
    new_data.push(row_data)
    if matrix.row_names.length() > 0 {
      new_row_names.push(matrix.row_names[rows[i]])
    }
    i = i + 1
  }
  let new_col_names : Array[String] = Array::new()
  if matrix.col_names.length() > 0 {
    let mut j = 0
    while j < cols.length() {
      new_col_names.push(matrix.col_names[cols[j]])
      j = j + 1
    }
  }
  DelayedMatrix::{
    data: new_data,
    nrow: rows.length(),
    ncol: cols.length(),
    row_names: new_row_names,
    col_names: new_col_names,
  }
}

///|
fn get_row_values(matrix : DelayedMatrix, i : Int) -> Array[Double] {
  matrix.data[i]
}

///|
fn get_col_values(matrix : DelayedMatrix, j : Int) -> Array[Double] {
  let result : Array[Double] = Array::new()
  let mut i = 0
  while i < matrix.nrow {
    result.push(matrix.data[i][j])
    i = i + 1
  }
  result
}

///|
fn compute_single_stat(values : Array[Double], stat : String) -> Double {
  let n = values.length()
  if n == 0 {
    return 0.0
  }

  let valid : Array[Double] = Array::new()
  let mut i = 0
  while i < n {
    if !values[i].is_nan() {
      valid.push(values[i])
    }
    i = i + 1
  }

  let nv = valid.length()
  if nv == 0 {
    return 0.0
  }

  match stat {
    "mean" => {
      let mut sum = 0.0
      let mut j = 0
      while j < nv {
        sum = sum + valid[j]
        j = j + 1
      }
      sum / nv.to_double()
    }
    "sum" => {
      let mut sum = 0.0
      let mut j = 0
      while j < nv {
        sum = sum + valid[j]
        j = j + 1
      }
      sum
    }
    "min" => {
      let mut min_val = valid[0]
      let mut j = 1
      while j < nv {
        if valid[j] < min_val {
          min_val = valid[j]
        }
        j = j + 1
      }
      min_val
    }
    "max" => {
      let mut max_val = valid[0]
      let mut j = 1
      while j < nv {
        if valid[j] > max_val {
          max_val = valid[j]
        }
        j = j + 1
      }
      max_val
    }
    "median" => {
      let sorted = valid.copy()
      sorted.sort()
      if nv % 2 == 0 {
        (sorted[nv / 2 - 1] + sorted[nv / 2]) / 2.0
      } else {
        sorted[nv / 2]
      }
    }
    "var" => {
      let mean = compute_single_stat(valid, "mean")
      let mut sum_sq = 0.0
      let mut j = 0
      while j < nv {
        let diff = valid[j] - mean
        sum_sq = sum_sq + diff * diff
        j = j + 1
      }
      if nv <= 1 {
        0.0
      } else {
        sum_sq / (nv.to_double() - 1.0)
      }
    }
    "sd" => {
      let variance = compute_single_stat(valid, "var")
      variance.sqrt()
    }
    "n" | "nna" | "nn" => {
      nv.to_double()
    }
    _ => 0.0
  }
}

///|
pub fn row_stats(matrix : DelayedMatrix, stat : String) -> Array[Double] {
  let result : Array[Double] = Array::new()
  let mut i = 0
  while i < matrix.nrow {
    let vals = get_row_values(matrix, i)
    result.push(compute_single_stat(vals, stat))
    i = i + 1
  }
  result
}

///|
pub fn col_stats(matrix : DelayedMatrix, stat : String) -> Array[Double] {
  let result : Array[Double] = Array::new()
  let mut j = 0
  while j < matrix.ncol {
    let vals = get_col_values(matrix, j)
    result.push(compute_single_stat(vals, stat))
    j = j + 1
  }
  result
}

///|
pub fn row_medians(matrix : DelayedMatrix) -> Array[Double] {
  row_stats(matrix, "median")
}

///|
pub fn col_medians(matrix : DelayedMatrix) -> Array[Double] {
  col_stats(matrix, "median")
}

///|
pub fn dms_row_means(matrix : DelayedMatrix) -> Array[Double] {
  row_stats(matrix, "mean")
}

///|
pub fn dms_col_means(matrix : DelayedMatrix) -> Array[Double] {
  col_stats(matrix, "mean")
}

///|
pub fn dms_row_vars(matrix : DelayedMatrix) -> Array[Double] {
  row_stats(matrix, "var")
}

///|
pub fn col_vars(matrix : DelayedMatrix) -> Array[Double] {
  col_stats(matrix, "var")
}

///|
pub fn row_allsums(matrix : DelayedMatrix) -> Array[Double] {
  row_stats(matrix, "sum")
}

///|
pub fn row_n(matrix : DelayedMatrix) -> Array[Int] {
  let result : Array[Int] = Array::new()
  let mut i = 0
  while i < matrix.nrow {
    let vals = get_row_values(matrix, i)
    let mut count = 0
    let mut j = 0
    while j < vals.length() {
      if !vals[j].is_nan() {
        count = count + 1
      }
      j = j + 1
    }
    result.push(count)
    i = i + 1
  }
  result
}

///|
pub fn col_n(matrix : DelayedMatrix) -> Array[Int] {
  let result : Array[Int] = Array::new()
  let mut j = 0
  while j < matrix.ncol {
    let vals = get_col_values(matrix, j)
    let mut count = 0
    let mut i = 0
    while i < vals.length() {
      if !vals[i].is_nan() {
        count = count + 1
      }
      i = i + 1
    }
    result.push(count)
    j = j + 1
  }
  result
}