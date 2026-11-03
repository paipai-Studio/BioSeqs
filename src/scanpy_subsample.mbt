/// Scanpy-inspired: pp.subsample — random subsampling of cells or genes

///|
/// Randomly subsample cells (rows) from a matrix.
///
/// Parameters:
/// - `x`: cells × genes matrix
/// - `fraction`: fraction of cells to keep (0-1)
/// - `seed`: random seed for reproducibility
pub fn scanpy_subsample_cells(
  x : Array[Array[Double]],
  fraction : Double,
  seed~ : Int32 = 0,
) -> Array[Array[Double]] {
  let n = x.length()
  if n == 0 || fraction >= 1.0 {
    return x
  }
  let target = (n.to_double() * fraction).floor().to_int()
  if target <= 0 {
    return []
  }
  // Fisher-Yates partial shuffle to select `target` indices
  let indices = Array(capacity=n)
  for i in 0..<n {
    indices.push(i)
  }
  let mut rng = seed
  for i in 0..<target {
    rng = (rng * 1103515245 + 12345) & 0x7FFFFFFF
    let j = i + rng.reinterpret_as_uint().reinterpret_as_int() % (n - i)
    let j_s = j
    let tmp = indices[i]
    indices[i] = indices[j_s]
    indices[j_s] = tmp
  }
  let result = Array(capacity=target)
  for i in 0..<target {
    result.push(x[indices[i]])
  }
  result
}

///|
/// Randomly subsample genes (columns) from a matrix.
///
/// Parameters:
/// - `x`: cells × genes matrix
/// - `fraction`: fraction of genes to keep (0-1)
/// - `seed`: random seed
pub fn scanpy_subsample_genes(
  x : Array[Array[Double]],
  fraction : Double,
  seed~ : Int32 = 0,
) -> Array[Array[Double]] {
  let n_cells = x.length()
  if n_cells == 0 {
    return []
  }
  let n_genes = x[0].length()
  if fraction >= 1.0 {
    return x
  }
  let target = (n_genes.to_double() * fraction).floor().to_int()
  if target <= 0 {
    return []
  }
  let indices = Array(capacity=n_genes)
  for i in 0..<n_genes {
    indices.push(i)
  }
  let mut rng = seed
  for i in 0..<target {
    rng = (rng * 1103515245 + 12345) & 0x7FFFFFFF
    let j = i + rng.reinterpret_as_uint().reinterpret_as_int() % (n_genes - i)
    let j_s = j
    let tmp = indices[i]
    indices[i] = indices[j_s]
    indices[j_s] = tmp
  }
  let result = Array(capacity=n_cells)
  for c in 0..<n_cells {
    let row = Array(capacity=target)
    for i in 0..<target {
      row.push(x[c][indices[i]])
    }
    result.push(row)
  }
  result
}

///|
/// Subsample both cells and genes.
pub fn scanpy_subsample(
  x : Array[Array[Double]],
  cell_fraction : Double,
  gene_fraction : Double,
  seed~ : Int32 = 0,
) -> Array[Array[Double]] {
  let after_cells = scanpy_subsample_cells(x, cell_fraction, seed~)
  scanpy_subsample_genes(after_cells, gene_fraction, seed=seed + 1)
}
