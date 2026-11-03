/// Scanpy-inspired: pp.downsample_counts — downsample counts to fixed total per cell

///|
/// Downsample counts per cell to a fixed total.
/// Each cell's counts are multinomial-sampled to `total_counts`.
///
/// Parameters:
/// - `x`: cells × genes count matrix
/// - `total_counts`: target total counts per cell
/// - `seed`: random seed
pub fn scanpy_downsample_counts(
  x : Array[Array[Double]],
  total_counts : Int32,
  seed~ : Int32 = 0,
) -> Array[Array[Double]] {
  let n_cells = x.length()
  if n_cells == 0 {
    return []
  }
  let n_genes = x[0].length()
  let result = Array(capacity=n_cells)
  let mut rng = seed
  for i in 0..<n_cells {
    // Compute current total
    let mut current_total = 0.0
    for j in 0..<n_genes {
      current_total += x[i][j]
    }
    let row = Array(capacity=n_genes)
    if current_total <= 0.0 || total_counts <= 0 {
      for _ in 0..<n_genes {
        row.push(0.0)
      }
    } else {
      // Multinomial sampling: for each count, pick a gene proportionally
      let _remaining = if total_counts.to_double() < current_total {
        total_counts
      } else {
        current_total.floor().to_int()
      }
      // Simple approach: scale counts and round
      let scale = total_counts.to_double() / current_total
      let mut assigned = 0
      for j in 0..<n_genes {
        let scaled = x[i][j] * scale
        let floored = scaled.floor().to_int()
        row.push(floored.to_double())
        assigned += floored
      }
      // Distribute remaining counts randomly
      let deficit = total_counts - assigned
      for _ in 0..<deficit {
        rng = (rng * 1103515245 + 12345) & 0x7FFFFFFF
        let j = rng.reinterpret_as_uint() % n_genes.reinterpret_as_uint()
        let j_s = j
        row[j_s.reinterpret_as_int()] += 1.0
      }
    }
    result.push(row)
  }
  result
}

///|
/// Downsample to a fraction of the median library size.
pub fn scanpy_downsample_to_median(
  x : Array[Array[Double]],
  seed~ : Int32 = 0,
) -> Array[Array[Double]] {
  let n_cells = x.length()
  if n_cells == 0 {
    return []
  }
  let n_genes = x[0].length()
  // Compute library sizes
  let lib_sizes = Array(capacity=n_cells)
  for i in 0..<n_cells {
    let mut s = 0.0
    for j in 0..<n_genes {
      s += x[i][j]
    }
    lib_sizes.push(s)
  }
  // Compute median
  let sorted = Array::from_iter(lib_sizes.iter())
  sorted.sort_by(fn(a, b) { if a < b { -1 } else if a > b { 1 } else { 0 } })
  let median = sorted[sorted.length() / 2]
  scanpy_downsample_counts(x, median.floor().to_int(), seed~)
}
