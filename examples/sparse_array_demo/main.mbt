///|
fn format_doubles(values : Array[Double]) -> String {
  let mut result = "["
  let mut index = 0
  while index < values.length() {
    if index > 0 {
      result = result + ", "
    }
    result = result + values[index].to_string()
    index = index + 1
  }
  result + "]"
}

///|
fn format_ints(values : Array[Int]) -> String {
  let mut result = "["
  let mut index = 0
  while index < values.length() {
    if index > 0 {
      result = result + ", "
    }
    result = result + values[index].to_string()
    index = index + 1
  }
  result + "]"
}

///|
fn format_matrix(values : Array[Array[Double]]) -> String {
  let mut result = "["
  let mut row = 0
  while row < values.length() {
    if row > 0 {
      result = result + ", "
    }
    result = result + format_doubles(values[row])
    row = row + 1
  }
  result + "]"
}

///|
fn main {
  println("=== Bioconductor SparseArray Demo ===\n")

  let tensor = @src.sparse_array_sample() catch {
    _ => abort("failed to create sparse tensor")
  }
  let tensor_summary = tensor.summary() catch {
    _ => abort("failed to summarize sparse tensor")
  }
  let tensor_value = tensor.get([0, 4, 1]) catch {
    _ => abort("failed to access sparse tensor")
  }
  println("1. Multidimensional sparse tensor")
  println("   " + tensor_summary)
  println("   value[0, 4, 1] = " + tensor_value.to_string())
  println(
    "   stored coordinates = " + tensor.nzcoordinates().length().to_string(),
  )

  let permuted = tensor.aperm([2, 0, 1]) catch {
    _ => abort("failed to permute sparse tensor")
  }
  let permuted_summary = permuted.summary() catch {
    _ => abort("failed to summarize permuted tensor")
  }
  let permuted_value = permuted.get([1, 0, 4]) catch {
    _ => abort("failed to access permuted tensor")
  }
  println("\n2. Dimension permutation")
  println("   " + permuted_summary)
  println("   old [0, 4, 1] -> new [1, 0, 4] = " + permuted_value.to_string())

  let sliced = tensor.slice([0, 0, 1], [4, 5, 2]) catch {
    _ => abort("failed to slice sparse tensor")
  }
  let sliced_summary = sliced.summary() catch {
    _ => abort("failed to summarize sparse slice")
  }
  println("\n3. Sparse slice without dense materialization")
  println("   " + sliced_summary)

  let counts = @src.SparseArray::from_matrix([
    [1.0, 0.0, 2.0],
    [0.0, 3.0, 0.0],
    [4.0, 0.0, 5.0],
  ]) catch {
    _ => abort("failed to create sparse count matrix")
  }
  let counts_summary = counts.summary() catch {
    _ => abort("failed to summarize sparse count matrix")
  }
  let row_sums = counts.row_sums() catch {
    _ => abort("failed to compute row sums")
  }
  let column_sums = counts.column_sums() catch {
    _ => abort("failed to compute column sums")
  }
  let row_nonzero_counts = counts.row_nonzero_counts() catch {
    _ => abort("failed to count row nonzero values")
  }
  println("\n4. Matrix summaries")
  println("   " + counts_summary)
  println("   row sums: " + format_doubles(row_sums))
  println("   column sums: " + format_doubles(column_sums))
  println("   row nonzero counts: " + format_ints(row_nonzero_counts))

  let design = @src.SparseArray::from_matrix([
    [1.0, 0.0],
    [0.0, 1.0],
    [1.0, 1.0],
  ]) catch {
    _ => abort("failed to create sparse design matrix")
  }
  let product = counts.matmul(design) catch {
    _ => abort("failed to multiply sparse matrices")
  }
  println("\n5. Sparse matrix multiplication")
  println("   result: " + format_matrix(product))

  let normalized = counts.scale(0.5) catch {
    _ => abort("failed to scale sparse matrix")
  }
  let combined = counts.add(normalized) catch {
    _ => abort("failed to add sparse matrices")
  }
  let combined_summary = combined.summary() catch {
    _ => abort("failed to summarize sparse arithmetic result")
  }
  let dense_combined = combined.to_dense_matrix() catch {
    _ => abort("failed to materialize sparse arithmetic result")
  }
  println("\n6. Sparse arithmetic")
  println("   counts + counts * 0.5: " + combined_summary)
  println("   dense view: " + format_matrix(dense_combined))
}
