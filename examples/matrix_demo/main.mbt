///|
/// Matrix Demo - Matrix operations demo.

pub fn main() {
  println("=== Matrix Demo ===")
  
  let dense = @src.create_example_dense_matrix()
  println("\nDense Matrix (3x3):")
  println("nrow: " + dense.nrow().to_string())
  println("ncol: " + dense.ncol().to_string())
  println("Element [0,0]: " + dense.get(0, 0).to_string())
  println("Element [2,2]: " + dense.get(2, 2).to_string())
  
  let row_sums = dense.row_sums()
  println("\nRow sums: " + array_to_string(row_sums))
  
  let col_means = dense.col_means()
  println("Column means: " + array_to_string(col_means))
  
  let t = dense.transpose()
  println("\nTransposed Matrix:")
  println("Element [0,1]: " + t.get(0, 1).to_string())
  
  let sparse = @src.create_example_sparse_matrix()
  println("\nSparse Matrix (CSC, 3x2):")
  println("nnz: " + sparse.nnz().to_string())
  println("Element [0,0]: " + sparse.get(0, 0).to_string())
  
  let mat1 = @src.BiocMatrix::dense_from_array([[1.0, 2.0], [3.0, 4.0]])
  let mat2 = @src.BiocMatrix::dense_from_array([[2.0, 0.0], [0.0, 2.0]])
  let product = mat1.multiply(mat2)
  
  if product is Some(_) {
    println("\nMatrix multiplication result:")
    println("Element [0,0]: " + product.unwrap().get(0, 0).to_string())
    println("Element [1,1]: " + product.unwrap().get(1, 1).to_string())
  }
  
  let norm = dense.norm(2.0)
  println("\nFrobenius norm: " + norm.to_string())
}

fn array_to_string(arr : Array[Double]) -> String {
  let mut result = "["
  let mut i = 0
  while i < arr.length() {
    if i > 0 {
      result = result + ", "
    }
    result = result + arr[i].to_string()
    i = i + 1
  }
  result = result + "]"
  result
}
