/// Scanpy second-round: Spatial transcriptomics analysis
/// Inspired by scanpy.tl.spatial_neighbors, spatial_autocorr, sepal

///|
/// Spatial neighbor graph result.
pub struct SpatialNeighborsResult {
  n_neighbors : Int32
  distances : Array[Array[Double]] // pairwise distances
  connectivities : Array[Array[Double]] // adjacency matrix
}

///|
/// Build spatial neighbor graph from spatial coordinates.
/// Uses k-nearest-neighbors or radius-based approach.
pub fn scanpy_spatial_neighbors(
  coords : Array[Array[Double]],
  n_neighbors~ : Int32 = 6,
  radius~ : Double = 0.0,
) -> SpatialNeighborsResult {
  let n = coords.length()
  let distances = Array(capacity=n)
  let connectivities = Array(capacity=n)
  for i in 0..<n {
    let row_d = Array(capacity=n)
    let row_c = Array(capacity=n)
    for j in 0..<n {
      let mut dist = 0.0
      let ndim = coords[i].length()
      for d in 0..<ndim {
        let diff = coords[i][d] - coords[j][d]
        dist = dist + diff * diff
      }
      dist = @math.pow(dist, 0.5)
      row_d.push(dist)
      if radius > 0.0 {
        if dist <= radius && i != j {
          row_c.push(1.0)
        } else {
          row_c.push(0.0)
        }
      } else {
        row_c.push(0.0)
      }
    }
    // KNN: mark n_neighbors nearest as connected
    if radius <= 0.0 {
      let indices = Array(capacity=n)
      for j in 0..<n {
        indices.push(j)
      }
      indices.sort_by(fn(a, b) {
        if row_d[a] < row_d[b] {
          -1
        } else if row_d[a] > row_d[b] {
          1
        } else {
          0
        }
      })
      let k = if n_neighbors < n { n_neighbors } else { n }
      for kk in 0..<k {
        if kk < indices.length() && indices[kk] != i {
          row_c[indices[kk]] = 1.0
        }
      }
    }
    distances.push(row_d)
    connectivities.push(row_c)
  }
  SpatialNeighborsResult::{ n_neighbors, distances, connectivities, }
}

///|
/// Moran's I — global spatial autocorrelation.
/// I = (n / S0) * (sum_i sum_j w_ij * (x_i - mean) * (x_j - mean)) / sum_i (x_i - mean)^2
pub fn scanpy_morans_i(
  values : Array[Double],
  weights : Array[Array[Double]],
) -> Double {
  let n = values.length()
  if n < 2 {
    return 0.0
  }
  let mut mean = 0.0
  for i in 0..<n {
    mean = mean + values[i]
  }
  mean = mean / n.to_double()
  let mut s0 = 0.0
  let mut numerator = 0.0
  let mut denominator = 0.0
  for i in 0..<n {
    let di = values[i] - mean
    denominator = denominator + di * di
    for j in 0..<n {
      s0 = s0 + weights[i][j]
      numerator = numerator + weights[i][j] * di * (values[j] - mean)
    }
  }
  if s0 == 0.0 || denominator == 0.0 {
    return 0.0
  }
  n.to_double() / s0 * numerator / denominator
}

///|
/// Geary's C — local spatial autocorrelation.
/// C = ((n-1) / (2*S0)) * (sum_i sum_j w_ij * (x_i - x_j)^2) / sum_i (x_i - mean)^2
pub fn scanpy_gearys_c(
  values : Array[Double],
  weights : Array[Array[Double]],
) -> Double {
  let n = values.length()
  if n < 2 {
    return 0.0
  }
  let mut mean = 0.0
  for i in 0..<n {
    mean = mean + values[i]
  }
  mean = mean / n.to_double()
  let mut s0 = 0.0
  let mut numerator = 0.0
  let mut denominator = 0.0
  for i in 0..<n {
    denominator = denominator + (values[i] - mean) * (values[i] - mean)
    for j in 0..<n {
      s0 = s0 + weights[i][j]
      let diff = values[i] - values[j]
      numerator = numerator + weights[i][j] * diff * diff
    }
  }
  if s0 == 0.0 || denominator == 0.0 {
    return 0.0
  }
  (n - 1).to_double() / (2.0 * s0) * numerator / denominator
}

///|
/// Sepal — detect spatially variable genes.
/// Returns (gene_index, morans_i, gearys_c) sorted by Moran's I descending.
pub fn scanpy_sepal(
  x : Array[Array[Double]],
  weights : Array[Array[Double]],
) -> Array[(Int32, Double, Double)] {
  let n_genes = if x.length() > 0 { x[0].length() } else { 0 }
  let results = Array(capacity=n_genes)
  for g in 0..<n_genes {
    let values = Array(capacity=x.length())
    for i in 0..<x.length() {
      values.push(x[i][g])
    }
    let mi = scanpy_morans_i(values, weights)
    let gc = scanpy_gearys_c(values, weights)
    results.push((g, mi, gc))
  }
  results.sort_by(fn(a, b) {
    if a.1 > b.1 {
      -1
    } else if a.1 < b.1 {
      1
    } else {
      0
    }
  })
  results
}

///|
/// Generate spatial count plot data: bin expression values into spatial grid.
/// Returns (grid_x, grid_y, value) tuples for visualization.
pub fn scanpy_spatial_count_plot(
  coords : Array[Array[Double]],
  values : Array[Double],
  grid_size~ : Int32 = 50,
) -> Array[(Int32, Int32, Double)] {
  let n = coords.length()
  if n == 0 {
    return []
  }
  let mut min_x = coords[0][0]
  let mut max_x = coords[0][0]
  let mut min_y = coords[0][1]
  let mut max_y = coords[0][1]
  for i in 0..<n {
    if coords[i][0] < min_x {
      min_x = coords[i][0]
    }
    if coords[i][0] > max_x {
      max_x = coords[i][0]
    }
    if coords[i][1] < min_y {
      min_y = coords[i][1]
    }
    if coords[i][1] > max_y {
      max_y = coords[i][1]
    }
  }
  let range_x = max_x - min_x
  let range_y = max_y - min_y
  if range_x == 0.0 || range_y == 0.0 {
    return []
  }
  // Accumulate values into grid
  let grid = Array(capacity=grid_size * grid_size)
  for _i in 0..<(grid_size * grid_size) {
    grid.push(0.0)
  }
  let counts = Array(capacity=grid_size * grid_size)
  for _i in 0..<(grid_size * grid_size) {
    grid.push(0.0)
    counts.push(0)
  }
  ignore(grid)
  let grid2 = Array(capacity=grid_size * grid_size)
  for _i in 0..<(grid_size * grid_size) {
    grid2.push(0.0)
  }
  let counts2 = Array(capacity=grid_size * grid_size)
  for _i in 0..<(grid_size * grid_size) {
    counts2.push(0)
  }
  for i in 0..<n {
    let gx = ((coords[i][0] - min_x) / range_x * (grid_size - 1).to_double())
      .floor()
      .to_int()
    let gy = ((coords[i][1] - min_y) / range_y * (grid_size - 1).to_double())
      .floor()
      .to_int()
    let idx = gy * grid_size + gx
    if idx >= 0 && idx < grid2.length() {
      grid2[idx] = grid2[idx] + values[i]
      counts2[idx] = counts2[idx] + 1
    }
  }
  let result = Array(capacity=grid_size * grid_size)
  for gy in 0..<grid_size {
    for gx in 0..<grid_size {
      let idx = gy * grid_size + gx
      let v = if counts2[idx] > 0 {
        grid2[idx] / counts2[idx].to_double()
      } else {
        0.0
      }
      result.push((gx, gy, v))
    }
  }
  result
}
