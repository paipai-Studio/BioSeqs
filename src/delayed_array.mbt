///|
pub struct DelayedArray {
  data : Array[Array[Double]]
  dim : Array[Int]
  delayed_ops : Array[String]
  chunk_size : Int
}

///|
pub fn DelayedArray::new(data : Array[Array[Double]]) -> DelayedArray {
  let rows = data.length()
  let cols = if rows > 0 { data[0].length() } else { 0 }
  DelayedArray::{ data, dim: [rows, cols], delayed_ops: [], chunk_size: 100 }
}

///|
pub fn DelayedArray::dim(self : DelayedArray) -> Array[Int] {
  self.dim
}

///|
pub fn DelayedArray::nrow(self : DelayedArray) -> Int {
  self.dim[0]
}

///|
pub fn DelayedArray::ncol(self : DelayedArray) -> Int {
  self.dim[1]
}

///|
pub fn da_slice(
  darr : DelayedArray,
  row_idx : Array[Int],
  col_idx : Array[Int],
) -> DelayedArray {
  let new_data : Array[Array[Double]] = Array::new()
  let mut i = 0
  while i < row_idx.length() {
    let row_data : Array[Double] = Array::new()
    let mut j = 0
    while j < col_idx.length() {
      row_data.push(darr.data[row_idx[i]][col_idx[j]])
      j = j + 1
    }
    new_data.push(row_data)
    i = i + 1
  }
  let result = DelayedArray::new(new_data)
  result.delayed_ops.push("slice")
  result
}

///|
pub fn da_transpose(darr : DelayedArray) -> DelayedArray {
  let rows = darr.nrow()
  let cols = darr.ncol()
  let new_data : Array[Array[Double]] = Array::new()
  let mut i = 0
  while i < cols {
    let row_data : Array[Double] = Array::new()
    let mut j = 0
    while j < rows {
      row_data.push(darr.data[j][i])
      j = j + 1
    }
    new_data.push(row_data)
    i = i + 1
  }
  let result = DelayedArray::new(new_data)
  result.delayed_ops.push("transpose")
  result
}

///|
pub fn da_subset_rows(
  darr : DelayedArray,
  indices : Array[Int],
) -> DelayedArray {
  let new_data : Array[Array[Double]] = Array::new()
  let mut i = 0
  while i < indices.length() {
    new_data.push(darr.data[indices[i]])
    i = i + 1
  }
  let result = DelayedArray::new(new_data)
  result.delayed_ops.push("subset_rows")
  result
}

///|
pub fn da_subset_cols(
  darr : DelayedArray,
  indices : Array[Int],
) -> DelayedArray {
  let new_data : Array[Array[Double]] = Array::new()
  let mut i = 0
  while i < darr.nrow() {
    let row_data : Array[Double] = Array::new()
    let mut j = 0
    while j < indices.length() {
      row_data.push(darr.data[i][indices[j]])
      j = j + 1
    }
    new_data.push(row_data)
    i = i + 1
  }
  let result = DelayedArray::new(new_data)
  result.delayed_ops.push("subset_cols")
  result
}

///|
pub fn da_row_sums(darr : DelayedArray) -> Array[Double] {
  let result : Array[Double] = Array::new()
  let mut i = 0
  while i < darr.nrow() {
    let mut sum = 0.0
    let mut j = 0
    while j < darr.ncol() {
      sum = sum + darr.data[i][j]
      j = j + 1
    }
    result.push(sum)
    i = i + 1
  }
  result
}

///|
pub fn da_col_sums(darr : DelayedArray) -> Array[Double] {
  let result : Array[Double] = Array::new()
  let mut j = 0
  while j < darr.ncol() {
    let mut sum = 0.0
    let mut i = 0
    while i < darr.nrow() {
      sum = sum + darr.data[i][j]
      i = i + 1
    }
    result.push(sum)
    j = j + 1
  }
  result
}

///|
pub fn da_row_means(darr : DelayedArray) -> Array[Double] {
  let result : Array[Double] = Array::new()
  let mut i = 0
  while i < darr.nrow() {
    let mut sum = 0.0
    let mut j = 0
    while j < darr.ncol() {
      sum = sum + darr.data[i][j]
      j = j + 1
    }
    result.push(sum / darr.ncol().to_double())
    i = i + 1
  }
  result
}

///|
pub fn da_col_means(darr : DelayedArray) -> Array[Double] {
  let result : Array[Double] = Array::new()
  let mut j = 0
  while j < darr.ncol() {
    let mut sum = 0.0
    let mut i = 0
    while i < darr.nrow() {
      sum = sum + darr.data[i][j]
      i = i + 1
    }
    result.push(sum / darr.nrow().to_double())
    j = j + 1
  }
  result
}

///|
pub fn da_compute(darr : DelayedArray) -> Array[Array[Double]] {
  darr.data
}

///|
pub fn da_summary(darr : DelayedArray) -> String {
  let mut result = "DelayedArray Summary\n"
  result = result + "====================\n"
  result = result +
    "Dimensions: " +
    darr.nrow().to_string() +
    " x " +
    darr.ncol().to_string() +
    "\n"
  result = result + "Chunk size: " + darr.chunk_size.to_string() + "\n"
  result = result +
    "Delayed operations: " +
    darr.delayed_ops.length().to_string() +
    "\n"

  if darr.delayed_ops.length() > 0 {
    result = result + "Operations: "
    let mut i = 0
    while i < darr.delayed_ops.length() {
      result = result + darr.delayed_ops[i]
      if i < darr.delayed_ops.length() - 1 {
        result = result + ", "
      }
      i = i + 1
    }
    result = result + "\n"
  }

  let row_sums = da_row_sums(darr)
  let col_sums = da_col_sums(darr)

  let mut total_sum = 0.0
  let mut i = 0
  while i < row_sums.length() {
    total_sum = total_sum + row_sums[i]
    i = i + 1
  }

  result = result + "Total sum: " + total_sum.to_string() + "\n"
  result = result + "Row sums (first 5): "
  let mut j = 0
  while j < 5 && j < row_sums.length() {
    result = result + row_sums[j].to_string()
    if j < 4 && j < row_sums.length() - 1 {
      result = result + ", "
    }
    j = j + 1
  }
  result = result + "\n"

  result = result + "Column sums (first 5): "
  j = 0
  while j < 5 && j < col_sums.length() {
    result = result + col_sums[j].to_string()
    if j < 4 && j < col_sums.length() - 1 {
      result = result + ", "
    }
    j = j + 1
  }
  result = result + "\n"

  result
}

///|
pub fn da_create_example_data() -> DelayedArray {
  let data : Array[Array[Double]] = Array::new()
  let mut i = 0
  while i < 100 {
    let row : Array[Double] = Array::new()
    let mut j = 0
    while j < 50 {
      row.push((i * j).to_double() / 1000.0)
      j = j + 1
    }
    data.push(row)
    i = i + 1
  }
  DelayedArray::new(data)
}
