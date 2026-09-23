///|
/// Portable core of Bioconductor FlowSOM 2.21.0.
/// Input matrices use cell x marker orientation.

///|
pub suberror FlowSomError {
  FlowSomError(String)
}

///|
pub enum FlowSomDistance {
  FlowSomManhattan
  FlowSomEuclidean
  FlowSomChebyshev
  FlowSomCosine
} derive(Eq, Debug)
pub extend FlowSomDistance with Eq::{equal, not_equal}
pub extend FlowSomDistance with @moonbitlang/core/debug.Debug::{to_repr}

///|
pub enum FlowSomInitialization {
  FlowSomRandom
  FlowSomKwsp
  FlowSomPca
} derive(Eq, Debug)
pub extend FlowSomInitialization with Eq::{equal, not_equal}
pub extend FlowSomInitialization with @moonbitlang/core/debug.Debug::{to_repr}

///|
pub fn flowsom_manhattan_distance() -> FlowSomDistance {
  FlowSomManhattan
}

///|
pub fn flowsom_euclidean_distance() -> FlowSomDistance {
  FlowSomEuclidean
}

///|
pub fn flowsom_chebyshev_distance() -> FlowSomDistance {
  FlowSomChebyshev
}

///|
pub fn flowsom_cosine_distance() -> FlowSomDistance {
  FlowSomCosine
}

///|
pub fn flowsom_random_initialization() -> FlowSomInitialization {
  FlowSomRandom
}

///|
pub fn flowsom_kwsp_initialization() -> FlowSomInitialization {
  FlowSomKwsp
}

///|
pub fn flowsom_pca_initialization() -> FlowSomInitialization {
  FlowSomPca
}

///|
pub struct FlowSomConfig {
  xdim : Int
  ydim : Int
  rlen : Int
  mst_runs : Int
  alpha_start : Double
  alpha_end : Double
  radius_start : Double
  radius_end : Double
  distance : FlowSomDistance
  initialization : FlowSomInitialization
  importance : Array[Double]
  metaclusters : Int
  meta_starts : Int
  outlier_mad : Double
  seed : Int
} derive(Debug)
pub extend FlowSomConfig with @moonbitlang/core/debug.Debug::{to_repr}

///|
pub fn FlowSomConfig::create(
  xdim? : Int = 10,
  ydim? : Int = 10,
  rlen? : Int = 10,
  mst_runs? : Int = 1,
  alpha_start? : Double = 0.05,
  alpha_end? : Double = 0.01,
  radius_start? : Double = -1.0,
  radius_end? : Double = 0.0,
  distance? : FlowSomDistance = FlowSomEuclidean,
  initialization? : FlowSomInitialization = FlowSomKwsp,
  importance? : Array[Double] = [],
  metaclusters? : Int = 0,
  meta_starts? : Int = 10,
  outlier_mad? : Double = 4.0,
  seed? : Int = 1,
) -> FlowSomConfig raise FlowSomError {
  if xdim <= 0 || ydim <= 0 {
    flowsom_fail("FlowSOM grid dimensions must be positive")
  }
  if rlen <= 0 || mst_runs <= 0 {
    flowsom_fail("FlowSOM rlen and mst_runs must be positive")
  }
  if !flowsom_is_finite(alpha_start) ||
    !flowsom_is_finite(alpha_end) ||
    alpha_start <= 0.0 ||
    alpha_end < 0.0 ||
    alpha_start < alpha_end {
    flowsom_fail(
      "FlowSOM learning rates must be finite, non-negative and decreasing",
    )
  }
  if !flowsom_is_finite(radius_start) ||
    !flowsom_is_finite(radius_end) ||
    radius_start < -1.0 ||
    radius_end < 0.0 ||
    (radius_start >= 0.0 && radius_start < radius_end) {
    flowsom_fail("FlowSOM radii must be finite, non-negative and decreasing")
  }
  if metaclusters < 0 || metaclusters > xdim * ydim {
    flowsom_fail("FlowSOM metaclusters must be between zero and node count")
  }
  if meta_starts <= 0 {
    flowsom_fail("FlowSOM meta-clustering starts must be positive")
  }
  if !flowsom_is_finite(outlier_mad) || outlier_mad < 0.0 {
    flowsom_fail("FlowSOM outlier MAD multiplier must be non-negative")
  }
  if seed <= 0 {
    flowsom_fail("FlowSOM random seed must be positive")
  }
  for value in importance {
    if !flowsom_is_finite(value) || value <= 0.0 {
      flowsom_fail("FlowSOM marker importance values must be positive")
    }
  }
  FlowSomConfig::{
    xdim,
    ydim,
    rlen,
    mst_runs,
    alpha_start,
    alpha_end,
    radius_start,
    radius_end,
    distance,
    initialization,
    importance: importance.copy(),
    metaclusters,
    meta_starts,
    outlier_mad,
    seed,
  }
}

///|
pub struct FlowSomMstEdge {
  from : Int
  to : Int
  weight : Double
} derive(Eq, Debug)
pub extend FlowSomMstEdge with Eq::{equal, not_equal}
pub extend FlowSomMstEdge with @moonbitlang/core/debug.Debug::{to_repr}

///|
pub struct FlowSomMapping {
  clusters : Array[Int]
  distances : Array[Double]
} derive(Debug)
pub extend FlowSomMapping with @moonbitlang/core/debug.Debug::{to_repr}

///|
pub struct FlowSomOutlierReport {
  median_distances : Array[Double?]
  mad_distances : Array[Double?]
  thresholds : Array[Double?]
  counts : Array[Int]
  maximum_distances : Array[Double?]
  per_cell : Array[Bool]
} derive(Debug)
pub extend FlowSomOutlierReport with @moonbitlang/core/debug.Debug::{to_repr}

///|
pub struct FlowSomModel {
  data : Array[Array[Double]]
  marker_names : Array[String]
  codes : Array[Array[Double]]
  grid : Array[(Int, Int)]
  mapping : FlowSomMapping
  mst_edges : Array[FlowSomMstEdge]
  topology_distances : Array[Array[Int]]
  node_counts : Array[Int]
  node_percentages : Array[Double]
  node_medians : Array[Array[Double?]]
  node_cvs : Array[Array[Double?]]
  node_sds : Array[Array[Double?]]
  node_mads : Array[Array[Double?]]
  meta_clusters : Array[Int]
  cell_meta_clusters : Array[Int]
  outliers : FlowSomOutlierReport
  quantization_error : Double
  topographic_error : Double
  iterations : Int
  config : FlowSomConfig
} derive(Debug)
pub extend FlowSomModel with @moonbitlang/core/debug.Debug::{to_repr}

///|
pub struct FlowSomMappedData {
  mapping : FlowSomMapping
  meta_clusters : Array[Int]
  outliers : Array[Bool]
} derive(Debug)
pub extend FlowSomMappedData with @moonbitlang/core/debug.Debug::{to_repr}

///|
pub struct FlowSomPurity {
  mean : Double
  worst : Double
  below_075 : Int
} derive(Debug)
pub extend FlowSomPurity with @moonbitlang/core/debug.Debug::{to_repr}

///|
pub struct FlowSomSceOutput {
  experiment : SingleCellExperiment
  model : FlowSomModel
}

///|
fn flowsom_fail(message : String) -> Unit raise FlowSomError {
  raise FlowSomError(message)
}

///|
fn flowsom_is_finite(value : Double) -> Bool {
  value == value && value.abs() <= 1.0e300
}

///|
fn flowsom_double_matrix(
  rows : Int,
  columns : Int,
  value : Double,
) -> Array[Array[Double]] {
  let output : Array[Array[Double]] = []
  for _ in 0..<rows {
    output.push(Array::make(columns, value))
  }
  output
}

///|
fn flowsom_int_matrix(
  rows : Int,
  columns : Int,
  value : Int,
) -> Array[Array[Int]] {
  let output : Array[Array[Int]] = []
  for _ in 0..<rows {
    output.push(Array::make(columns, value))
  }
  output
}

///|
fn flowsom_optional_matrix(rows : Int, columns : Int) -> Array[Array[Double?]] {
  let output : Array[Array[Double?]] = []
  for _ in 0..<rows {
    output.push(Array::make(columns, None))
  }
  output
}

///|
fn flowsom_copy_matrix(data : Array[Array[Double]]) -> Array[Array[Double]] {
  data.map(fn(row) { row.copy() })
}

///|
fn flowsom_matrix_shape(
  data : Array[Array[Double]],
) -> (Int, Int) raise FlowSomError {
  if data.length() == 0 {
    flowsom_fail("FlowSOM data must contain at least one cell")
  }
  let markers = data[0].length()
  if markers == 0 {
    flowsom_fail("FlowSOM data must contain at least one marker")
  }
  for row in data {
    if row.length() != markers {
      flowsom_fail("FlowSOM data must be rectangular")
    }
    for value in row {
      if !flowsom_is_finite(value) {
        flowsom_fail("FlowSOM data must contain only finite values")
      }
    }
  }
  (data.length(), markers)
}

///|
fn flowsom_lcg(state : Int) -> (Int, Double) {
  let normalized = (state - 1) % 2147483646 + 1
  let quotient = normalized / 44488
  let remainder = normalized % 44488
  let candidate = 48271 * remainder - 3399 * quotient
  let next = if candidate > 0 { candidate } else { candidate + 2147483647 }
  (next, next.to_double() / 2147483647.0)
}

///|
pub fn flowsom_distance(
  left : Array[Double],
  right : Array[Double],
  distance : FlowSomDistance,
) -> Double raise FlowSomError {
  if left.length() == 0 || left.length() != right.length() {
    flowsom_fail("FlowSOM distance vectors must be non-empty and equal-length")
  }
  for value in left {
    if !flowsom_is_finite(value) {
      flowsom_fail("FlowSOM distance vectors must be finite")
    }
  }
  for value in right {
    if !flowsom_is_finite(value) {
      flowsom_fail("FlowSOM distance vectors must be finite")
    }
  }
  flowsom_distance_unchecked(left, right, distance)
}

///|
fn flowsom_distance_unchecked(
  left : Array[Double],
  right : Array[Double],
  distance : FlowSomDistance,
) -> Double {
  match distance {
    FlowSomManhattan => {
      let mut total = 0.0
      for index in 0..<left.length() {
        total = total + (left[index] - right[index]).abs()
      }
      total
    }
    FlowSomEuclidean => {
      let mut total = 0.0
      for index in 0..<left.length() {
        let delta = left[index] - right[index]
        total = total + delta * delta
      }
      total.sqrt()
    }
    FlowSomChebyshev => {
      let mut maximum = 0.0
      for index in 0..<left.length() {
        maximum = maximum.max((left[index] - right[index]).abs())
      }
      maximum
    }
    FlowSomCosine => {
      let mut dot = 0.0
      let mut left_norm = 0.0
      let mut right_norm = 0.0
      for index in 0..<left.length() {
        dot = dot + left[index] * right[index]
        left_norm = left_norm + left[index] * left[index]
        right_norm = right_norm + right[index] * right[index]
      }
      if left_norm == 0.0 && right_norm == 0.0 {
        0.0
      } else if left_norm == 0.0 || right_norm == 0.0 {
        1.0
      } else {
        1.0 - dot / (left_norm * right_norm).sqrt()
      }
    }
  }
}

///|
fn flowsom_weight_data(
  data : Array[Array[Double]],
  importance : Array[Double],
) -> Array[Array[Double]] {
  if importance.length() == 0 {
    return flowsom_copy_matrix(data)
  }
  data.map(fn(row) { row.mapi(fn(index, value) { value * importance[index] }) })
}

///|
fn flowsom_grid(xdim : Int, ydim : Int) -> Array[(Int, Int)] {
  let grid : Array[(Int, Int)] = []
  for y in 0..<ydim {
    for x in 0..<xdim {
      grid.push((x, y))
    }
  }
  grid
}

///|
fn flowsom_grid_distances(grid : Array[(Int, Int)]) -> Array[Array[Int]] {
  let distances = flowsom_int_matrix(grid.length(), grid.length(), 0)
  for left in 0..<grid.length() {
    let (left_x, left_y) = grid[left]
    for right in 0..<grid.length() {
      let (right_x, right_y) = grid[right]
      distances[left][right] = (left_x - right_x)
        .abs()
        .max((left_y - right_y).abs())
    }
  }
  distances
}

///|
fn flowsom_default_radius(grid : Array[(Int, Int)]) -> Double {
  let distances : Array[Double] = []
  for left in 0..<grid.length() {
    let (left_x, left_y) = grid[left]
    for right in (left + 1)..<grid.length() {
      let (right_x, right_y) = grid[right]
      distances.push(
        (left_x - right_x).abs().max((left_y - right_y).abs()).to_double(),
      )
    }
  }
  if distances.length() == 0 {
    return 0.0
  }
  distances.sort_by(fn(left : Double, right : Double) -> Int {
    if left < right {
      -1
    } else if left > right {
      1
    } else {
      0
    }
  })
  let position = (0.67 * (distances.length() - 1).to_double()).to_int()
  distances[position]
}

///|
pub fn flowsom_initialize_kwsp(
  data : Array[Array[Double]],
  nodes : Int,
  seed? : Int = 1,
) -> Array[Array[Double]] raise FlowSomError {
  let (cells, _) = flowsom_matrix_shape(data)
  if nodes <= 0 || nodes > cells {
    flowsom_fail("FlowSOM KWSP nodes must be between one and cell count")
  }
  if seed <= 0 {
    flowsom_fail("FlowSOM random seed must be positive")
  }
  let (_, first_draw) = flowsom_lcg(seed)
  let first = (first_draw * cells.to_double()).to_int().min(cells - 1)
  let selected : Array[Int] = [first]
  let distances : Array[Double] = []
  for cell in 0..<cells {
    let distance = flowsom_distance_unchecked(
      data[cell],
      data[first],
      FlowSomEuclidean,
    )
    distances.push(distance * distance)
  }
  while selected.length() < nodes {
    let mut best = -1
    let mut best_distance = -1.0
    for cell in 0..<cells {
      if !selected.contains(cell) &&
        (
          distances[cell] > best_distance ||
          (distances[cell] == best_distance && (best < 0 || cell < best))
        ) {
        best = cell
        best_distance = distances[cell]
      }
    }
    selected.push(best)
    for cell in 0..<cells {
      let distance = flowsom_distance_unchecked(
        data[cell],
        data[best],
        FlowSomEuclidean,
      )
      distances[cell] = distances[cell].min(distance * distance)
    }
  }
  selected.map(fn(index) { data[index].copy() })
}

///|
fn flowsom_initialize_random(
  data : Array[Array[Double]],
  nodes : Int,
  seed : Int,
) -> Array[Array[Double]] raise FlowSomError {
  let (cells, _) = flowsom_matrix_shape(data)
  if nodes <= 0 || nodes > cells {
    flowsom_fail("FlowSOM random initialization requires nodes <= cell count")
  }
  let indices = Array::make(cells, 0)
  for index in 0..<cells {
    indices[index] = index
  }
  let mut state = seed
  for index in 0..<nodes {
    let draw = flowsom_lcg(state)
    state = draw.0
    let target = (index + (draw.1 * (cells - index).to_double()).to_int()).min(
      cells - 1,
    )
    let temporary = indices[index]
    indices[index] = indices[target]
    indices[target] = temporary
  }
  let output : Array[Array[Double]] = []
  for index in 0..<nodes {
    output.push(data[indices[index]].copy())
  }
  output
}

///|
fn flowsom_column_means(data : Array[Array[Double]]) -> Array[Double] {
  let means = Array::make(data[0].length(), 0.0)
  for row in data {
    for marker in 0..<row.length() {
      means[marker] = means[marker] + row[marker]
    }
  }
  for marker in 0..<means.length() {
    means[marker] = means[marker] / data.length().to_double()
  }
  means
}

///|
fn flowsom_covariance(
  data : Array[Array[Double]],
  means : Array[Double],
) -> Array[Array[Double]] {
  let markers = means.length()
  let covariance = flowsom_double_matrix(markers, markers, 0.0)
  let denominator = (data.length() - 1).max(1).to_double()
  for left in 0..<markers {
    for right in left..<markers {
      let mut total = 0.0
      for row in data {
        total = total + (row[left] - means[left]) * (row[right] - means[right])
      }
      let value = total / denominator
      covariance[left][right] = value
      covariance[right][left] = value
    }
  }
  covariance
}

///|
fn flowsom_vector_norm(values : Array[Double]) -> Double {
  let mut total = 0.0
  for value in values {
    total = total + value * value
  }
  total.sqrt()
}

///|
fn flowsom_power_component(
  covariance : Array[Array[Double]],
  excluded : Array[Double]?,
) -> (Array[Double], Double) {
  let size = covariance.length()
  let vector = Array::make(size, 0.0)
  for index in 0..<size {
    vector[index] = 1.0 / (index + 1).to_double()
  }
  match excluded {
    Some(previous) => {
      let mut projection = 0.0
      for index in 0..<size {
        projection = projection + vector[index] * previous[index]
      }
      for index in 0..<size {
        vector[index] = vector[index] - projection * previous[index]
      }
    }
    None => ()
  }
  let initial_norm = flowsom_vector_norm(vector)
  if initial_norm > 0.0 {
    for index in 0..<size {
      vector[index] = vector[index] / initial_norm
    }
  }
  for _ in 0..<100 {
    let next = Array::make(size, 0.0)
    for row in 0..<size {
      for column in 0..<size {
        next[row] = next[row] + covariance[row][column] * vector[column]
      }
    }
    match excluded {
      Some(previous) => {
        let mut projection = 0.0
        for index in 0..<size {
          projection = projection + next[index] * previous[index]
        }
        for index in 0..<size {
          next[index] = next[index] - projection * previous[index]
        }
      }
      None => ()
    }
    let norm = flowsom_vector_norm(next)
    if norm <= 1.0e-12 {
      break
    }
    let mut difference = 0.0
    for index in 0..<size {
      next[index] = next[index] / norm
      difference = difference + (next[index] - vector[index]).abs()
      vector[index] = next[index]
    }
    if difference <= 1.0e-10 {
      break
    }
  }
  let product = Array::make(size, 0.0)
  for row in 0..<size {
    for column in 0..<size {
      product[row] = product[row] + covariance[row][column] * vector[column]
    }
  }
  let mut eigenvalue = 0.0
  for index in 0..<size {
    eigenvalue = eigenvalue + vector[index] * product[index]
  }
  (vector, eigenvalue.max(0.0))
}

///|
pub fn flowsom_initialize_pca(
  data : Array[Array[Double]],
  xdim : Int,
  ydim : Int,
) -> Array[Array[Double]] raise FlowSomError {
  let (_, markers) = flowsom_matrix_shape(data)
  if xdim <= 0 || ydim <= 0 {
    flowsom_fail("FlowSOM PCA grid dimensions must be positive")
  }
  let means = flowsom_column_means(data)
  let covariance = flowsom_covariance(data, means)
  let first = flowsom_power_component(covariance, None)
  let second = if markers > 1 {
    flowsom_power_component(covariance, Some(first.0))
  } else {
    (Array::make(markers, 0.0), 0.0)
  }
  let first_scale = 5.0 * first.1.sqrt()
  let second_scale = 5.0 * second.1.sqrt()
  let codes : Array[Array[Double]] = []
  for y in 0..<ydim {
    let y_position = if ydim == 1 {
      0.0
    } else {
      2.0 * y.to_double() / (ydim - 1).to_double() - 1.0
    }
    for x in 0..<xdim {
      let x_position = if xdim == 1 {
        0.0
      } else {
        2.0 * x.to_double() / (xdim - 1).to_double() - 1.0
      }
      let code = Array::make(markers, 0.0)
      for marker in 0..<markers {
        code[marker] = means[marker] +
          first.0[marker] * first_scale * x_position +
          second.0[marker] * second_scale * y_position
      }
      codes.push(code)
    }
  }
  codes
}

///|
fn flowsom_mst_unchecked(codes : Array[Array[Double]]) -> Array[FlowSomMstEdge] {
  let nodes = codes.length()
  if nodes <= 1 {
    return []
  }
  let included = Array::make(nodes, false)
  let best = Array::make(nodes, 1.0e300)
  let parents = Array::make(nodes, -1)
  best[0] = 0.0
  let edges : Array[FlowSomMstEdge] = []
  for _ in 0..<nodes {
    let mut selected = -1
    let mut selected_distance = 1.0e301
    for node in 0..<nodes {
      if !included[node] &&
        (
          best[node] < selected_distance ||
          (best[node] == selected_distance && (selected < 0 || node < selected))
        ) {
        selected = node
        selected_distance = best[node]
      }
    }
    included[selected] = true
    if parents[selected] >= 0 {
      edges.push(FlowSomMstEdge::{
        from: parents[selected],
        to: selected,
        weight: best[selected],
      })
    }
    for candidate in 0..<nodes {
      if !included[candidate] {
        let distance = flowsom_distance_unchecked(
          codes[selected],
          codes[candidate],
          FlowSomEuclidean,
        )
        if distance < best[candidate] ||
          (
            distance == best[candidate] &&
            (parents[candidate] < 0 || selected < parents[candidate])
          ) {
          best[candidate] = distance
          parents[candidate] = selected
        }
      }
    }
  }
  edges
}

///|
pub fn flowsom_build_mst(
  codes : Array[Array[Double]],
) -> Array[FlowSomMstEdge] raise FlowSomError {
  ignore(flowsom_matrix_shape(codes))
  flowsom_mst_unchecked(codes)
}

///|
fn flowsom_topology_from_mst(
  nodes : Int,
  edges : Array[FlowSomMstEdge],
) -> Array[Array[Int]] {
  let infinity = nodes + 1
  let distances = flowsom_int_matrix(nodes, nodes, infinity)
  for node in 0..<nodes {
    distances[node][node] = 0
  }
  for edge in edges {
    distances[edge.from][edge.to] = 1
    distances[edge.to][edge.from] = 1
  }
  for middle in 0..<nodes {
    for left in 0..<nodes {
      for right in 0..<nodes {
        distances[left][right] = distances[left][right].min(
          distances[left][middle] + distances[middle][right],
        )
      }
    }
  }
  distances
}

///|
pub fn flowsom_mst_distances(
  codes : Array[Array[Double]],
) -> Array[Array[Int]] raise FlowSomError {
  ignore(flowsom_matrix_shape(codes))
  flowsom_topology_from_mst(codes.length(), flowsom_mst_unchecked(codes))
}

///|
fn flowsom_map_unchecked(
  codes : Array[Array[Double]],
  data : Array[Array[Double]],
  distance : FlowSomDistance,
) -> FlowSomMapping {
  let clusters : Array[Int] = []
  let distances : Array[Double] = []
  for row in data {
    let mut best = 0
    let mut best_distance = flowsom_distance_unchecked(row, codes[0], distance)
    for node in 1..<codes.length() {
      let current = flowsom_distance_unchecked(row, codes[node], distance)
      if current < best_distance {
        best = node
        best_distance = current
      }
    }
    clusters.push(best)
    distances.push(best_distance)
  }
  FlowSomMapping::{ clusters, distances, }
}

///|
pub fn flowsom_map_data(
  codes : Array[Array[Double]],
  data : Array[Array[Double]],
  distance? : FlowSomDistance = FlowSomEuclidean,
) -> FlowSomMapping raise FlowSomError {
  let (_, code_markers) = flowsom_matrix_shape(codes)
  let (_, data_markers) = flowsom_matrix_shape(data)
  if code_markers != data_markers {
    flowsom_fail("FlowSOM codes and data must have equal marker counts")
  }
  flowsom_map_unchecked(codes, data, distance)
}

///|
fn flowsom_training_stage(
  data : Array[Array[Double]],
  codes : Array[Array[Double]],
  topology : Array[Array[Int]],
  alpha_start : Double,
  alpha_end : Double,
  radius_start : Double,
  radius_end : Double,
  rlen : Int,
  distance : FlowSomDistance,
  initial_state : Int,
) -> (Int, Int) {
  let total_iterations = rlen * data.length()
  let radius_step = (radius_start - radius_end) / total_iterations.to_double()
  let mut radius = radius_start
  let mut state = initial_state
  let mut completed = 0
  let mut epoch_change = 1.0
  let mut stop = false
  for iteration in 0..<total_iterations {
    if iteration % data.length() == 0 {
      if epoch_change < 1.0 {
        stop = true
      }
      epoch_change = 0.0
    }
    if stop {
      break
    }
    let draw = flowsom_lcg(state)
    state = draw.0
    let selected = (draw.1 * data.length().to_double())
      .to_int()
      .min(data.length() - 1)
    let mut nearest = 0
    let mut nearest_distance = flowsom_distance_unchecked(
      data[selected],
      codes[0],
      distance,
    )
    for node in 1..<codes.length() {
      let current = flowsom_distance_unchecked(
        data[selected],
        codes[node],
        distance,
      )
      if current < nearest_distance {
        nearest = node
        nearest_distance = current
      }
    }
    if radius < 1.0 {
      radius = 0.5
    }
    let alpha = alpha_start -
      (alpha_start - alpha_end) *
      iteration.to_double() /
      total_iterations.to_double()
    for node in 0..<codes.length() {
      if topology[node][nearest].to_double() <= radius {
        for marker in 0..<codes[node].length() {
          let change = data[selected][marker] - codes[node][marker]
          epoch_change = epoch_change + change.abs()
          codes[node][marker] = codes[node][marker] + change * alpha
        }
      }
    }
    radius = radius - radius_step
    completed = completed + 1
  }
  (state, completed)
}

///|
fn flowsom_median(values : Array[Double]) -> Double {
  let sorted = values.copy()
  sorted.sort_by(fn(left : Double, right : Double) -> Int {
    if left < right {
      -1
    } else if left > right {
      1
    } else {
      0
    }
  })
  let middle = sorted.length() / 2
  if sorted.length() % 2 == 0 {
    (sorted[middle - 1] + sorted[middle]) / 2.0
  } else {
    sorted[middle]
  }
}

///|
fn flowsom_mad(values : Array[Double], center : Double) -> Double {
  let deviations = values.map(fn(value) { (value - center).abs() })
  1.4826 * flowsom_median(deviations)
}

///|
priv struct FlowSomNodeSummary {
  counts : Array[Int]
  percentages : Array[Double]
  medians : Array[Array[Double?]]
  cvs : Array[Array[Double?]]
  sds : Array[Array[Double?]]
  mads : Array[Array[Double?]]
}

///|
fn flowsom_node_summary(
  data : Array[Array[Double]],
  clusters : Array[Int],
  nodes : Int,
) -> FlowSomNodeSummary {
  let markers = data[0].length()
  let counts = Array::make(nodes, 0)
  for cluster in clusters {
    counts[cluster] = counts[cluster] + 1
  }
  let percentages = counts.map(fn(count) {
    count.to_double() / data.length().to_double()
  })
  let medians = flowsom_optional_matrix(nodes, markers)
  let cvs = flowsom_optional_matrix(nodes, markers)
  let sds = flowsom_optional_matrix(nodes, markers)
  let mads = flowsom_optional_matrix(nodes, markers)
  for node in 0..<nodes {
    if counts[node] > 0 {
      for marker in 0..<markers {
        let values : Array[Double] = []
        for cell in 0..<data.length() {
          if clusters[cell] == node {
            values.push(data[cell][marker])
          }
        }
        let median = flowsom_median(values)
        let mut mean = 0.0
        for value in values {
          mean = mean + value
        }
        mean = mean / values.length().to_double()
        let mut variance = 0.0
        for value in values {
          let delta = value - mean
          variance = variance + delta * delta
        }
        let sd = if values.length() > 1 {
          (variance / (values.length() - 1).to_double()).sqrt()
        } else {
          0.0
        }
        medians[node][marker] = Some(median)
        sds[node][marker] = Some(sd)
        mads[node][marker] = Some(flowsom_mad(values, median))
        if mean != 0.0 {
          cvs[node][marker] = Some(sd / mean)
        }
      }
    }
  }
  FlowSomNodeSummary::{ counts, percentages, medians, cvs, sds, mads, }
}

///|
fn flowsom_outlier_report(
  mapping : FlowSomMapping,
  nodes : Int,
  mad_allowed : Double,
) -> FlowSomOutlierReport {
  let medians : Array[Double?] = Array::make(nodes, None)
  let mads : Array[Double?] = Array::make(nodes, None)
  let thresholds : Array[Double?] = Array::make(nodes, None)
  let maximums : Array[Double?] = Array::make(nodes, None)
  let counts = Array::make(nodes, 0)
  for node in 0..<nodes {
    let values : Array[Double] = []
    for cell in 0..<mapping.clusters.length() {
      if mapping.clusters[cell] == node {
        values.push(mapping.distances[cell])
      }
    }
    if values.length() > 0 {
      let median = flowsom_median(values)
      let mad = flowsom_mad(values, median)
      let mut maximum = values[0]
      for value in values {
        maximum = maximum.max(value)
      }
      medians[node] = Some(median)
      mads[node] = Some(mad)
      thresholds[node] = Some(median + mad_allowed * mad)
      maximums[node] = Some(maximum)
    }
  }
  let per_cell = Array::make(mapping.clusters.length(), false)
  for cell in 0..<mapping.clusters.length() {
    let node = mapping.clusters[cell]
    match thresholds[node] {
      Some(threshold) =>
        if mapping.distances[cell] > threshold {
          per_cell[cell] = true
          counts[node] = counts[node] + 1
        }
      None => ()
    }
  }
  FlowSomOutlierReport::{
    median_distances: medians,
    mad_distances: mads,
    thresholds,
    counts,
    maximum_distances: maximums,
    per_cell,
  }
}

///|
pub fn flowsom_meta_cluster_codes(
  codes : Array[Array[Double]],
  clusters : Int,
  starts? : Int = 10,
  seed? : Int = 1,
) -> Array[Int] raise FlowSomError {
  let (nodes, _) = flowsom_matrix_shape(codes)
  if clusters <= 0 || clusters > nodes {
    flowsom_fail(
      "FlowSOM meta-cluster count must be between one and node count",
    )
  }
  if starts <= 0 || seed <= 0 {
    flowsom_fail("FlowSOM meta-clustering starts and seed must be positive")
  }
  if clusters == 1 {
    return Array::make(nodes, 0)
  }
  let config = BlusterKmeansConfig::create(clusters, starts~, seed~) catch {
    BlusterError(message) =>
      raise FlowSomError("FlowSOM meta-clustering setup failed: " + message)
  }
  let result = bluster_cluster_rows_kmeans(codes, config) catch {
    BlusterError(message) =>
      raise FlowSomError("FlowSOM meta-clustering failed: " + message)
  }
  result.clusters
}

///|
fn flowsom_line_residual(
  values : Array[Double],
  start : Int,
  end : Int,
) -> Double {
  let count = end - start
  if count <= 1 {
    return 0.0
  }
  let mut mean_x = 0.0
  let mut mean_y = 0.0
  for index in start..<end {
    mean_x = mean_x + (index + 1).to_double()
    mean_y = mean_y + values[index]
  }
  mean_x = mean_x / count.to_double()
  mean_y = mean_y / count.to_double()
  let mut numerator = 0.0
  let mut denominator = 0.0
  for index in start..<end {
    let centered_x = (index + 1).to_double() - mean_x
    numerator = numerator + centered_x * (values[index] - mean_y)
    denominator = denominator + centered_x * centered_x
  }
  let slope = if denominator > 0.0 { numerator / denominator } else { 0.0 }
  let intercept = mean_y - slope * mean_x
  let mut residual = 0.0
  for index in start..<end {
    residual = residual +
      (values[index] - (intercept + slope * (index + 1).to_double())).abs()
  }
  residual
}

///|
pub fn flowsom_determine_meta_clusters(
  codes : Array[Array[Double]],
  maximum : Int,
  starts? : Int = 10,
  seed? : Int = 1,
) -> Int raise FlowSomError {
  let (nodes, _) = flowsom_matrix_shape(codes)
  if maximum <= 0 || maximum > nodes {
    flowsom_fail("FlowSOM maximum meta-clusters must be within node count")
  }
  if starts <= 0 || seed <= 0 {
    flowsom_fail("FlowSOM meta-clustering starts and seed must be positive")
  }
  if maximum <= 2 {
    return maximum
  }
  let sse = Array::make(maximum, 0.0)
  for count in 1..<=maximum {
    if count == 1 {
      let means = flowsom_column_means(codes)
      let mut total = 0.0
      for code in codes {
        for marker in 0..<code.length() {
          let delta = code[marker] - means[marker]
          total = total + delta * delta
        }
      }
      sse[count - 1] = total
    } else {
      let config = BlusterKmeansConfig::create(
        count,
        starts~,
        seed=seed + count * 7919,
      ) catch {
        BlusterError(message) =>
          raise FlowSomError(
            "FlowSOM automatic meta-clustering setup failed: " + message,
          )
      }
      let result = bluster_cluster_rows_kmeans(codes, config) catch {
        BlusterError(message) =>
          raise FlowSomError(
            "FlowSOM automatic meta-clustering failed: " + message,
          )
      }
      sse[count - 1] = result.total_within_sum_squares
    }
  }
  let smoothed = sse.copy()
  for index in 1..<(maximum - 1) {
    smoothed[index] = 0.8 * sse[index] +
      0.1 * sse[index - 1] +
      0.1 * sse[index + 1]
  }
  let mut optimum = 2
  let mut best_residual = 1.0e300
  for split in 1..<(maximum - 1) {
    let residual = flowsom_line_residual(smoothed, 0, split) +
      flowsom_line_residual(smoothed, split, maximum)
    if residual < best_residual {
      best_residual = residual
      optimum = split + 1
    }
  }
  optimum
}

///|
fn flowsom_quality(
  data : Array[Array[Double]],
  codes : Array[Array[Double]],
  mapping : FlowSomMapping,
  topology : Array[Array[Int]],
  distance : FlowSomDistance,
) -> (Double, Double) {
  let mut quantization = 0.0
  let mut topographic_errors = 0
  for cell in 0..<data.length() {
    quantization = quantization + mapping.distances[cell]
    if codes.length() > 1 {
      let first = mapping.clusters[cell]
      let mut second = if first == 0 { 1 } else { 0 }
      let mut second_distance = flowsom_distance_unchecked(
        data[cell],
        codes[second],
        distance,
      )
      for node in 0..<codes.length() {
        if node != first {
          let current = flowsom_distance_unchecked(
            data[cell],
            codes[node],
            distance,
          )
          if current < second_distance {
            second = node
            second_distance = current
          }
        }
      }
      if topology[first][second] > 1 {
        topographic_errors = topographic_errors + 1
      }
    }
  }
  (
    quantization / data.length().to_double(),
    topographic_errors.to_double() / data.length().to_double(),
  )
}

///|
pub fn flowsom_train(
  data : Array[Array[Double]],
  config : FlowSomConfig,
  marker_names? : Array[String] = [],
) -> FlowSomModel raise FlowSomError {
  let (cells, markers) = flowsom_matrix_shape(data)
  let nodes = config.xdim * config.ydim
  if marker_names.length() != 0 && marker_names.length() != markers {
    flowsom_fail("FlowSOM marker names must match the data columns")
  }
  if config.importance.length() != 0 && config.importance.length() != markers {
    flowsom_fail("FlowSOM importance must contain one value per marker")
  }
  if config.initialization != FlowSomPca && nodes > cells {
    flowsom_fail(
      "FlowSOM random and KWSP initialization require at least one cell per node",
    )
  }
  let names = if marker_names.length() == 0 {
    let generated : Array[String] = []
    for marker in 0..<markers {
      generated.push("marker." + (marker + 1).to_string())
    }
    generated
  } else {
    let seen : Array[String] = []
    for name in marker_names {
      if name.length() == 0 || seen.contains(name) {
        flowsom_fail("FlowSOM marker names must be non-empty and unique")
      }
      seen.push(name)
    }
    marker_names.copy()
  }
  let weighted = flowsom_weight_data(data, config.importance)
  let grid = flowsom_grid(config.xdim, config.ydim)
  let codes = match config.initialization {
    FlowSomRandom => flowsom_initialize_random(weighted, nodes, config.seed)
    FlowSomKwsp => flowsom_initialize_kwsp(weighted, nodes, seed=config.seed)
    FlowSomPca => flowsom_initialize_pca(weighted, config.xdim, config.ydim)
  }
  let mut topology = flowsom_grid_distances(grid)
  let resolved_radius = if config.radius_start < 0.0 {
    flowsom_default_radius(grid)
  } else {
    config.radius_start
  }
  let mut state = config.seed
  let mut iterations = 0
  for stage in 0..<config.mst_runs {
    let start_fraction = stage.to_double() / config.mst_runs.to_double()
    let end_fraction = (stage + 1).to_double() / config.mst_runs.to_double()
    let stage_alpha_start = config.alpha_start +
      (config.alpha_end - config.alpha_start) * start_fraction
    let stage_alpha_end = config.alpha_start +
      (config.alpha_end - config.alpha_start) * end_fraction
    let stage_radius_start = resolved_radius +
      (config.radius_end - resolved_radius) * start_fraction
    let stage_radius_end = resolved_radius +
      (config.radius_end - resolved_radius) * end_fraction
    let trained = flowsom_training_stage(
      weighted,
      codes,
      topology,
      stage_alpha_start,
      stage_alpha_end,
      stage_radius_start,
      stage_radius_end,
      config.rlen,
      config.distance,
      state,
    )
    state = trained.0
    iterations = iterations + trained.1
    topology = flowsom_topology_from_mst(nodes, flowsom_mst_unchecked(codes))
  }
  let mst_edges = flowsom_mst_unchecked(codes)
  let mapping = flowsom_map_unchecked(codes, weighted, config.distance)
  let summary = flowsom_node_summary(data, mapping.clusters, nodes)
  let meta_clusters = if config.metaclusters > 0 {
    flowsom_meta_cluster_codes(
      codes,
      config.metaclusters,
      starts=config.meta_starts,
      seed=config.seed + 104729,
    )
  } else {
    []
  }
  let cell_meta_clusters = if meta_clusters.length() > 0 {
    mapping.clusters.map(fn(cluster) { meta_clusters[cluster] })
  } else {
    []
  }
  let outliers = flowsom_outlier_report(mapping, nodes, config.outlier_mad)
  let quality = flowsom_quality(
    weighted,
    codes,
    mapping,
    topology,
    config.distance,
  )
  FlowSomModel::{
    data: flowsom_copy_matrix(data),
    marker_names: names,
    codes: flowsom_copy_matrix(codes),
    grid,
    mapping,
    mst_edges,
    topology_distances: topology,
    node_counts: summary.counts,
    node_percentages: summary.percentages,
    node_medians: summary.medians,
    node_cvs: summary.cvs,
    node_sds: summary.sds,
    node_mads: summary.mads,
    meta_clusters,
    cell_meta_clusters,
    outliers,
    quantization_error: quality.0,
    topographic_error: quality.1,
    iterations,
    config,
  }
}

///|
pub fn flowsom_map_new(
  model : FlowSomModel,
  data : Array[Array[Double]],
) -> FlowSomMappedData raise FlowSomError {
  let (_, markers) = flowsom_matrix_shape(data)
  if markers != model.marker_names.length() {
    flowsom_fail("FlowSOM new data marker count does not match the model")
  }
  let weighted = flowsom_weight_data(data, model.config.importance)
  let mapping = flowsom_map_unchecked(
    model.codes,
    weighted,
    model.config.distance,
  )
  let meta_clusters = if model.meta_clusters.length() > 0 {
    mapping.clusters.map(fn(cluster) { model.meta_clusters[cluster] })
  } else {
    []
  }
  let outliers = Array::make(data.length(), false)
  for cell in 0..<data.length() {
    match model.outliers.thresholds[mapping.clusters[cell]] {
      Some(threshold) => outliers[cell] = mapping.distances[cell] > threshold
      None => ()
    }
  }
  FlowSomMappedData::{ mapping, meta_clusters, outliers, }
}

///|
pub fn flowsom_marker_outliers(
  model : FlowSomModel,
  data : Array[Array[Double]],
  mad_allowed? : Double = 4.0,
) -> Array[Array[Int]] raise FlowSomError {
  let (_, markers) = flowsom_matrix_shape(data)
  if markers != model.marker_names.length() {
    flowsom_fail("FlowSOM marker outlier data does not match the model")
  }
  if !flowsom_is_finite(mad_allowed) || mad_allowed < 0.0 {
    flowsom_fail("FlowSOM marker outlier MAD multiplier must be non-negative")
  }
  let mapped = flowsom_map_new(model, data)
  let output = flowsom_int_matrix(data.length(), markers, 0)
  for cell in 0..<data.length() {
    let node = mapped.mapping.clusters[cell]
    for marker in 0..<markers {
      match (model.node_medians[node][marker], model.node_mads[node][marker]) {
        (Some(median), Some(mad)) => {
          let upper = median + mad_allowed * mad
          let lower = median - mad_allowed * mad
          if data[cell][marker] > upper {
            output[cell][marker] = 1
          } else if data[cell][marker] < lower {
            output[cell][marker] = -1
          }
        }
        _ => ()
      }
    }
  }
  output
}

///|
pub fn flowsom_node_positive_percentages(
  model : FlowSomModel,
  cutoffs : Array[Double],
) -> Array[Array[Double?]] raise FlowSomError {
  if cutoffs.length() != model.marker_names.length() {
    flowsom_fail("FlowSOM positive cutoffs must match marker count")
  }
  for cutoff in cutoffs {
    if !flowsom_is_finite(cutoff) {
      flowsom_fail("FlowSOM positive cutoffs must be finite")
    }
  }
  let nodes = model.codes.length()
  let markers = model.marker_names.length()
  let output = flowsom_optional_matrix(nodes, markers)
  for node in 0..<nodes {
    if model.node_counts[node] > 0 {
      for marker in 0..<markers {
        let mut positive = 0
        for cell in 0..<model.data.length() {
          if model.mapping.clusters[cell] == node &&
            model.data[cell][marker] > cutoffs[marker] {
            positive = positive + 1
          }
        }
        output[node][marker] = Some(
          positive.to_double() / model.node_counts[node].to_double(),
        )
      }
    }
  }
  output
}

///|
pub fn flowsom_meta_counts(
  model : FlowSomModel,
) -> Array[Int] raise FlowSomError {
  if model.meta_clusters.length() == 0 {
    flowsom_fail("FlowSOM model does not contain meta-clusters")
  }
  let mut groups = 0
  for cluster in model.meta_clusters {
    groups = groups.max(cluster + 1)
  }
  let counts = Array::make(groups, 0)
  for cluster in model.cell_meta_clusters {
    counts[cluster] = counts[cluster] + 1
  }
  counts
}

///|
pub fn flowsom_meta_medians(
  model : FlowSomModel,
) -> Array[Array[Double?]] raise FlowSomError {
  if model.cell_meta_clusters.length() == 0 {
    flowsom_fail("FlowSOM model does not contain meta-clusters")
  }
  let mut groups = 0
  for cluster in model.meta_clusters {
    groups = groups.max(cluster + 1)
  }
  flowsom_node_summary(model.data, model.cell_meta_clusters, groups).medians
}

///|
fn flowsom_unique_labels(labels : Array[Int]) -> Array[Int] {
  let output : Array[Int] = []
  for label in labels {
    if !output.contains(label) {
      output.push(label)
    }
  }
  output.sort_by(fn(left : Int, right : Int) -> Int {
    if left < right {
      -1
    } else if left > right {
      1
    } else {
      0
    }
  })
  output
}

///|
fn flowsom_label_index(labels : Array[Int], value : Int) -> Int {
  for index in 0..<labels.length() {
    if labels[index] == value {
      return index
    }
  }
  -1
}

///|
pub fn flowsom_purity(
  real_clusters : Array[Int],
  predicted_clusters : Array[Int],
  weighted? : Bool = true,
) -> FlowSomPurity raise FlowSomError {
  if real_clusters.length() == 0 ||
    real_clusters.length() != predicted_clusters.length() {
    flowsom_fail(
      "FlowSOM purity cluster vectors must be non-empty and equal-length",
    )
  }
  let real = flowsom_unique_labels(real_clusters)
  let predicted = flowsom_unique_labels(predicted_clusters)
  let table = flowsom_int_matrix(predicted.length(), real.length(), 0)
  for index in 0..<real_clusters.length() {
    let row = flowsom_label_index(predicted, predicted_clusters[index])
    let column = flowsom_label_index(real, real_clusters[index])
    table[row][column] = table[row][column] + 1
  }
  let mut mean = 0.0
  let mut worst = 1.0
  let mut below = 0
  for row in 0..<predicted.length() {
    let mut total = 0
    let mut maximum = 0
    for column in 0..<real.length() {
      total = total + table[row][column]
      maximum = maximum.max(table[row][column])
    }
    let purity = maximum.to_double() / total.to_double()
    worst = worst.min(purity)
    if purity < 0.75 {
      below = below + 1
    }
    if weighted {
      mean = mean +
        purity * total.to_double() / real_clusters.length().to_double()
    } else {
      mean = mean + purity / predicted.length().to_double()
    }
  }
  FlowSomPurity::{ mean, worst, below_075: below, }
}

///|
pub fn flowsom_f_measure(
  real_clusters : Array[Int],
  predicted_clusters : Array[Int],
) -> Double raise FlowSomError {
  if real_clusters.length() == 0 ||
    real_clusters.length() != predicted_clusters.length() {
    flowsom_fail(
      "FlowSOM F-measure cluster vectors must be non-empty and equal-length",
    )
  }
  let real = flowsom_unique_labels(real_clusters)
  let predicted = flowsom_unique_labels(predicted_clusters)
  let table = flowsom_int_matrix(real.length(), predicted.length(), 0)
  for index in 0..<real_clusters.length() {
    let row = flowsom_label_index(real, real_clusters[index])
    let column = flowsom_label_index(predicted, predicted_clusters[index])
    table[row][column] = table[row][column] + 1
  }
  let row_sums = Array::make(real.length(), 0)
  let column_sums = Array::make(predicted.length(), 0)
  for row in 0..<real.length() {
    for column in 0..<predicted.length() {
      row_sums[row] = row_sums[row] + table[row][column]
      column_sums[column] = column_sums[column] + table[row][column]
    }
  }
  let mut score = 0.0
  for row in 0..<real.length() {
    let mut best = 0.0
    for column in 0..<predicted.length() {
      let precision = if column_sums[column] > 0 {
        table[row][column].to_double() / column_sums[column].to_double()
      } else {
        0.0
      }
      let recall = table[row][column].to_double() / row_sums[row].to_double()
      let value = if precision + recall > 0.0 {
        2.0 * precision * recall / (precision + recall)
      } else {
        0.0
      }
      best = best.max(value)
    }
    score = score +
      best * row_sums[row].to_double() / real_clusters.length().to_double()
  }
  score
}

///|
fn flowsom_config_with_metaclusters(
  config : FlowSomConfig,
  metaclusters : Int,
) -> FlowSomConfig {
  FlowSomConfig::{
    xdim: config.xdim,
    ydim: config.ydim,
    rlen: config.rlen,
    mst_runs: config.mst_runs,
    alpha_start: config.alpha_start,
    alpha_end: config.alpha_end,
    radius_start: config.radius_start,
    radius_end: config.radius_end,
    distance: config.distance,
    initialization: config.initialization,
    importance: config.importance.copy(),
    metaclusters,
    meta_starts: config.meta_starts,
    outlier_mad: config.outlier_mad,
    seed: config.seed,
  }
}

///|
pub fn flowsom_with_metaclusters(
  model : FlowSomModel,
  clusters : Int,
  starts? : Int = 10,
  seed? : Int = 1,
) -> FlowSomModel raise FlowSomError {
  let meta_clusters = flowsom_meta_cluster_codes(
    model.codes,
    clusters,
    starts~,
    seed~,
  )
  let cell_meta_clusters = model.mapping.clusters.map(fn(cluster) {
    meta_clusters[cluster]
  })
  FlowSomModel::{
    data: flowsom_copy_matrix(model.data),
    marker_names: model.marker_names.copy(),
    codes: flowsom_copy_matrix(model.codes),
    grid: model.grid.copy(),
    mapping: FlowSomMapping::{
      clusters: model.mapping.clusters.copy(),
      distances: model.mapping.distances.copy(),
    },
    mst_edges: model.mst_edges.copy(),
    topology_distances: model.topology_distances.map(fn(row) { row.copy() }),
    node_counts: model.node_counts.copy(),
    node_percentages: model.node_percentages.copy(),
    node_medians: model.node_medians.map(fn(row) { row.copy() }),
    node_cvs: model.node_cvs.map(fn(row) { row.copy() }),
    node_sds: model.node_sds.map(fn(row) { row.copy() }),
    node_mads: model.node_mads.map(fn(row) { row.copy() }),
    meta_clusters,
    cell_meta_clusters,
    outliers: FlowSomOutlierReport::{
      median_distances: model.outliers.median_distances.copy(),
      mad_distances: model.outliers.mad_distances.copy(),
      thresholds: model.outliers.thresholds.copy(),
      counts: model.outliers.counts.copy(),
      maximum_distances: model.outliers.maximum_distances.copy(),
      per_cell: model.outliers.per_cell.copy(),
    },
    quantization_error: model.quantization_error,
    topographic_error: model.topographic_error,
    iterations: model.iterations,
    config: flowsom_config_with_metaclusters(model.config, clusters),
  }
}

///|
pub fn flowsom_with_automatic_metaclusters(
  model : FlowSomModel,
  maximum : Int,
  starts? : Int = 10,
  seed? : Int = 1,
) -> FlowSomModel raise FlowSomError {
  let clusters = flowsom_determine_meta_clusters(
    model.codes,
    maximum,
    starts~,
    seed~,
  )
  flowsom_with_metaclusters(model, clusters, starts~, seed~)
}

///|
fn flowsom_resolve_marker_indices(
  marker_count : Int,
  marker_indices : Array[Int],
) -> Array[Int] raise FlowSomError {
  if marker_count <= 0 {
    flowsom_fail("FlowSOM input must contain at least one marker")
  }
  if marker_indices.length() == 0 {
    let all = Array::make(marker_count, 0)
    for marker in 0..<marker_count {
      all[marker] = marker
    }
    return all
  }
  let selected : Array[Int] = []
  for marker in marker_indices {
    if marker < 0 || marker >= marker_count {
      flowsom_fail("FlowSOM marker index is outside the input matrix")
    }
    if selected.contains(marker) {
      flowsom_fail("FlowSOM marker indices must be unique")
    }
    selected.push(marker)
  }
  selected
}

///|
pub fn flowsom_train_flow_frame(
  frame : FlowFrame,
  config : FlowSomConfig,
  marker_indices? : Array[Int] = [],
) -> FlowSomModel raise FlowSomError {
  let (_, markers) = flowsom_matrix_shape(frame.exprs)
  let selected = flowsom_resolve_marker_indices(markers, marker_indices)
  let available_names = flow_core_get_markers(frame)
  let names : Array[String] = []
  let data : Array[Array[Double]] = []
  for marker in selected {
    if marker < available_names.length() && available_names[marker].length() > 0 {
      names.push(available_names[marker])
    } else {
      names.push("marker." + (marker + 1).to_string())
    }
  }
  for source in frame.exprs {
    let row : Array[Double] = []
    for marker in selected {
      row.push(source[marker])
    }
    data.push(row)
  }
  flowsom_train(data, config, marker_names=names)
}

///|
fn flowsom_copy_sce(source : SingleCellExperiment) -> SingleCellExperiment {
  let assays : Map[String, Array[Array[Double]]] = Map([])
  for key in source.assays.keys() {
    assays[key] = flowsom_copy_matrix(source.assays[key])
  }
  let row_data : Map[String, Array[String]] = Map([])
  for key in source.row_data.keys() {
    row_data[key] = source.row_data[key].copy()
  }
  let col_data : Map[String, Array[String]] = Map([])
  for key in source.col_data.keys() {
    col_data[key] = source.col_data[key].copy()
  }
  let reduced_dims : Map[String, Array[Array[Double]]] = Map([])
  for key in source.reduced_dims.keys() {
    reduced_dims[key] = flowsom_copy_matrix(source.reduced_dims[key])
  }
  let alternatives : Map[String, SingleCellExperiment] = Map([])
  for key in source.alternative_experiments.keys() {
    alternatives[key] = flowsom_copy_sce(source.alternative_experiments[key])
  }
  let metadata : Map[String, String] = Map([])
  for key in source.metadata.keys() {
    metadata[key] = source.metadata[key]
  }
  SingleCellExperiment::{
    assays,
    row_data,
    col_data,
    reduced_dims,
    alternative_experiments: alternatives,
    metadata,
    row_names: source.row_names.copy(),
    col_names: source.col_names.copy(),
  }
}

///|
pub fn flowsom_cluster_sce(
  experiment : SingleCellExperiment,
  config : FlowSomConfig,
  assay_name? : String = "counts",
  marker_indices? : Array[Int] = [],
  cluster_column? : String = "FlowSOM.cluster",
  meta_column? : String = "FlowSOM.metacluster",
) -> FlowSomSceOutput raise FlowSomError {
  if assay_name.length() == 0 ||
    cluster_column.length() == 0 ||
    meta_column.length() == 0 {
    flowsom_fail("FlowSOM SCE assay and output names must be non-empty")
  }
  if cluster_column == meta_column {
    flowsom_fail("FlowSOM SCE cluster and meta-cluster columns must differ")
  }
  let assay = match experiment.assays.get(assay_name) {
    Some(value) => value
    None => {
      flowsom_fail("FlowSOM SCE assay '" + assay_name + "' was not found")
      []
    }
  }
  if assay.length() == 0 {
    flowsom_fail("FlowSOM SCE assay must contain at least one marker")
  }
  if experiment.row_names.length() != assay.length() {
    flowsom_fail("FlowSOM SCE assay rows must match row names")
  }
  let cells = experiment.col_names.length()
  if cells == 0 {
    flowsom_fail("FlowSOM SCE must contain at least one cell")
  }
  for row in assay {
    if row.length() != cells {
      flowsom_fail("FlowSOM SCE assay columns must match column names")
    }
  }
  let selected = flowsom_resolve_marker_indices(assay.length(), marker_indices)
  let names = selected.map(fn(marker) { experiment.row_names[marker] })
  let data : Array[Array[Double]] = []
  for cell in 0..<cells {
    let row : Array[Double] = []
    for marker in selected {
      row.push(assay[marker][cell])
    }
    data.push(row)
  }
  let model = flowsom_train(data, config, marker_names=names)
  let output = flowsom_copy_sce(experiment)
  output.col_data[cluster_column] = model.mapping.clusters.map(fn(cluster) {
    (cluster + 1).to_string()
  })
  if model.cell_meta_clusters.length() > 0 {
    output.col_data[meta_column] = model.cell_meta_clusters.map(fn(cluster) {
      (cluster + 1).to_string()
    })
  }
  output.metadata["FlowSOM.assay"] = assay_name
  output.metadata["FlowSOM.grid"] = config.xdim.to_string() +
    "x" +
    config.ydim.to_string()
  output.metadata["FlowSOM.quantization_error"] = model.quantization_error.to_string()
  output.metadata["FlowSOM.topographic_error"] = model.topographic_error.to_string()
  FlowSomSceOutput::{ experiment: output, model, }
}

///|
pub fn flowsom_summary(model : FlowSomModel) -> String {
  let mut occupied = 0
  for count in model.node_counts {
    if count > 0 {
      occupied = occupied + 1
    }
  }
  let meta = if model.meta_clusters.length() == 0 {
    "none"
  } else {
    let mut groups = 0
    for cluster in model.meta_clusters {
      groups = groups.max(cluster + 1)
    }
    groups.to_string()
  }
  "FlowSOM model\n" +
  "  cells: " +
  model.data.length().to_string() +
  "\n  markers: " +
  model.marker_names.length().to_string() +
  "\n  nodes: " +
  model.codes.length().to_string() +
  "\n  occupied nodes: " +
  occupied.to_string() +
  "\n  metaclusters: " +
  meta +
  "\n  quantization error: " +
  model.quantization_error.to_string() +
  "\n  topographic error: " +
  model.topographic_error.to_string()
}
