///|
/// Portable advanced core of Bioconductor slingshot 2.21.0.
/// Coordinates use cell x dimension orientation.

///|
pub suberror SlingshotAdvancedError {
  SlingshotAdvancedError(String)
}

///|
pub enum SlingshotClusterDistance {
  SlingshotCenterEuclidean
  SlingshotScaledDiagonal
  SlingshotScaledFull
} derive(Eq, Debug)

///|
pub enum SlingshotCurveExtension {
  SlingshotExtendNone
  SlingshotExtendLine
  SlingshotExtendPc1
} derive(Eq, Debug)

///|
pub fn slingshot_center_euclidean() -> SlingshotClusterDistance {
  SlingshotCenterEuclidean
}

///|
pub fn slingshot_scaled_diagonal() -> SlingshotClusterDistance {
  SlingshotScaledDiagonal
}

///|
pub fn slingshot_scaled_full() -> SlingshotClusterDistance {
  SlingshotScaledFull
}

///|
pub fn slingshot_extend_none() -> SlingshotCurveExtension {
  SlingshotExtendNone
}

///|
pub fn slingshot_extend_line() -> SlingshotCurveExtension {
  SlingshotExtendLine
}

///|
pub fn slingshot_extend_pc1() -> SlingshotCurveExtension {
  SlingshotExtendPc1
}

///|
pub struct SlingshotAdvancedConfig {
  start_clusters : Array[String]
  end_clusters : Array[String]
  distance : SlingshotClusterDistance
  use_median : Bool
  omega : Double
  automatic_omega : Bool
  omega_scale : Double
  shrink : Double
  extension : SlingshotCurveExtension
  reweight : Bool
  reassign : Bool
  max_iterations : Int
  tolerance : Double
  smoother_span : Double
  curve_points : Int
  covariance_ridge : Double
} derive(Debug)

///|
pub fn SlingshotAdvancedConfig::create(
  start_clusters? : Array[String] = [],
  end_clusters? : Array[String] = [],
  distance? : SlingshotClusterDistance = SlingshotScaledFull,
  use_median? : Bool = false,
  omega? : Double = -1.0,
  automatic_omega? : Bool = false,
  omega_scale? : Double = 1.5,
  shrink? : Double = 1.0,
  extension? : SlingshotCurveExtension = SlingshotExtendLine,
  reweight? : Bool = true,
  reassign? : Bool = true,
  max_iterations? : Int = 15,
  tolerance? : Double = 0.001,
  smoother_span? : Double = 0.3,
  curve_points? : Int = 100,
  covariance_ridge? : Double = 1.0e-6,
) -> SlingshotAdvancedConfig raise SlingshotAdvancedError {
  if !slingshot_adv_is_finite(omega) || omega < -1.0 {
    slingshot_adv_fail("slingshot omega must be -1 or non-negative")
  }
  if !slingshot_adv_is_finite(omega_scale) || omega_scale <= 0.0 {
    slingshot_adv_fail("slingshot omega scale must be positive")
  }
  if !slingshot_adv_is_finite(shrink) || shrink < 0.0 || shrink > 1.0 {
    slingshot_adv_fail("slingshot shrink must be between zero and one")
  }
  if max_iterations <= 0 {
    slingshot_adv_fail("slingshot iteration limit must be positive")
  }
  if !slingshot_adv_is_finite(tolerance) || tolerance <= 0.0 {
    slingshot_adv_fail("slingshot tolerance must be positive")
  }
  if !slingshot_adv_is_finite(smoother_span) ||
    smoother_span <= 0.0 ||
    smoother_span > 1.0 {
    slingshot_adv_fail("slingshot smoother span must be in (0, 1]")
  }
  if curve_points < 2 {
    slingshot_adv_fail("slingshot curve point count must be at least two")
  }
  if !slingshot_adv_is_finite(covariance_ridge) || covariance_ridge <= 0.0 {
    slingshot_adv_fail("slingshot covariance ridge must be positive")
  }
  for name in start_clusters {
    if name.length() == 0 {
      slingshot_adv_fail("slingshot start cluster names cannot be empty")
    }
  }
  for name in end_clusters {
    if name.length() == 0 {
      slingshot_adv_fail("slingshot end cluster names cannot be empty")
    }
  }
  SlingshotAdvancedConfig::{
    start_clusters: start_clusters.copy(),
    end_clusters: end_clusters.copy(),
    distance,
    use_median,
    omega,
    automatic_omega,
    omega_scale,
    shrink,
    extension,
    reweight,
    reassign,
    max_iterations,
    tolerance,
    smoother_span,
    curve_points,
    covariance_ridge,
  }
}

///|
pub struct SlingshotClusterSummary {
  name : String
  center : Array[Double]
  covariance : Array[Array[Double]]
  size : Double
} derive(Debug)

///|
pub struct SlingshotLineageEdge {
  from : Int
  to : Int
  distance : Double
} derive(Eq, Debug)

///|
pub struct SlingshotLineageModel {
  coordinates : Array[Array[Double]]
  cluster_names : Array[String]
  cluster_weights : Array[Array[Double]]
  clusters : Array[SlingshotClusterSummary]
  distance_matrix : Array[Array[Double]]
  edges : Array[SlingshotLineageEdge]
  components : Array[Int]
  roots : Array[Int]
  lineages : Array[Array[Int]]
  initial_weights : Array[Array[Double]]
  omega_threshold : Double
} derive(Debug)

///|
pub struct SlingshotAdvancedCurve {
  name : String
  lineage : Array[Int]
  points : Array[Array[Double]]
  pseudotime : Array[Double]
  distances : Array[Double]
  weights : Array[Double]
  length : Double
} derive(Debug)

///|
pub struct SlingshotAdvancedResult {
  lineage_model : SlingshotLineageModel
  curves : Array[SlingshotAdvancedCurve]
  pseudotime : Array[Array[Double?]]
  weights : Array[Array[Double]]
  average_pseudotime : Array[Double]
  branch_ids : Array[String]
  iterations : Int
  converged : Bool
  total_distance : Double
  config : SlingshotAdvancedConfig
} derive(Debug)

///|
pub struct SlingshotProjection {
  pseudotime : Array[Array[Double?]]
  weights : Array[Array[Double]]
  distances : Array[Array[Double]]
  average_pseudotime : Array[Double]
  branch_ids : Array[String]
} derive(Debug)

///|
pub struct SlingshotAdvancedSceOutput {
  experiment : SingleCellExperiment
  result : SlingshotAdvancedResult
}

///|
fn slingshot_adv_fail(message : String) -> Unit raise SlingshotAdvancedError {
  raise SlingshotAdvancedError(message)
}

///|
fn slingshot_adv_is_finite(value : Double) -> Bool {
  value == value && value.abs() <= 1.0e300
}

///|
fn slingshot_adv_double_matrix(
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
fn slingshot_adv_optional_matrix(
  rows : Int,
  columns : Int,
) -> Array[Array[Double?]] {
  let output : Array[Array[Double?]] = []
  for _ in 0..<rows {
    output.push(Array::make(columns, None))
  }
  output
}

///|
fn slingshot_adv_copy_matrix(
  matrix : Array[Array[Double]],
) -> Array[Array[Double]] {
  matrix.map(fn(row) { row.copy() })
}

///|
fn slingshot_adv_shape(
  matrix : Array[Array[Double]],
  name : String,
) -> (Int, Int) raise SlingshotAdvancedError {
  if matrix.length() == 0 {
    slingshot_adv_fail(name + " must contain at least one row")
  }
  let columns = matrix[0].length()
  if columns == 0 {
    slingshot_adv_fail(name + " must contain at least one column")
  }
  for row in matrix {
    if row.length() != columns {
      slingshot_adv_fail(name + " must be rectangular")
    }
    for value in row {
      if !slingshot_adv_is_finite(value) {
        slingshot_adv_fail(name + " must contain only finite values")
      }
    }
  }
  (matrix.length(), columns)
}

///|
fn slingshot_adv_unique_labels(labels : Array[String]) -> Array[String] {
  let output : Array[String] = []
  for label in labels {
    if label != "-1" && !output.contains(label) {
      output.push(label)
    }
  }
  output
}

///|
fn slingshot_adv_label_index(names : Array[String], name : String) -> Int {
  for index in 0..<names.length() {
    if names[index] == name {
      return index
    }
  }
  -1
}

///|
fn slingshot_adv_validate_cluster_names(
  cluster_names : Array[String],
) -> Unit raise SlingshotAdvancedError {
  if cluster_names.length() == 0 {
    slingshot_adv_fail("slingshot must contain at least one cluster")
  }
  for index in 0..<cluster_names.length() {
    if cluster_names[index].length() == 0 || cluster_names[index] == "-1" {
      slingshot_adv_fail("slingshot cluster names must be non-empty and not -1")
    }
    for other in 0..<index {
      if cluster_names[other] == cluster_names[index] {
        slingshot_adv_fail("slingshot cluster names must be unique")
      }
    }
  }
}

///|
fn slingshot_adv_normalize_weights(
  weights : Array[Array[Double]],
  cluster_names : Array[String],
  cells : Int,
) -> Array[Array[Double]] raise SlingshotAdvancedError {
  if weights.length() != cells {
    slingshot_adv_fail("slingshot cluster weights must have one row per cell")
  }
  let output = slingshot_adv_double_matrix(cells, cluster_names.length(), 0.0)
  let totals = Array::make(cluster_names.length(), 0.0)
  for cell in 0..<cells {
    if weights[cell].length() != cluster_names.length() {
      slingshot_adv_fail(
        "slingshot cluster weight columns must match cluster names",
      )
    }
    let mut total = 0.0
    for cluster in 0..<cluster_names.length() {
      let value = weights[cell][cluster]
      if !slingshot_adv_is_finite(value) || value < 0.0 {
        slingshot_adv_fail(
          "slingshot cluster weights must be finite and non-negative",
        )
      }
      total = total + value
    }
    if total > 0.0 {
      for cluster in 0..<cluster_names.length() {
        output[cell][cluster] = weights[cell][cluster] / total
        totals[cluster] = totals[cluster] + output[cell][cluster]
      }
    }
  }
  for cluster in 0..<cluster_names.length() {
    if totals[cluster] <= 0.0 {
      slingshot_adv_fail(
        "slingshot cluster '" + cluster_names[cluster] + "' has no cells",
      )
    }
  }
  output
}

///|
fn slingshot_adv_weighted_mean(
  coordinates : Array[Array[Double]],
  weights : Array[Array[Double]],
  cluster : Int,
) -> Array[Double] {
  let dimensions = coordinates[0].length()
  let center = Array::make(dimensions, 0.0)
  let mut total = 0.0
  for cell in 0..<coordinates.length() {
    let weight = weights[cell][cluster]
    total = total + weight
    for dimension in 0..<dimensions {
      center[dimension] = center[dimension] +
        weight * coordinates[cell][dimension]
    }
  }
  for dimension in 0..<dimensions {
    center[dimension] = center[dimension] / total
  }
  center
}

///|
fn slingshot_adv_weighted_median(
  coordinates : Array[Array[Double]],
  weights : Array[Array[Double]],
  cluster : Int,
) -> Array[Double] {
  let dimensions = coordinates[0].length()
  let center : Array[Double] = []
  for dimension in 0..<dimensions {
    let entries : Array[(Double, Double, Int)] = []
    let mut total = 0.0
    for cell in 0..<coordinates.length() {
      let weight = weights[cell][cluster]
      if weight > 0.0 {
        entries.push((coordinates[cell][dimension], weight, cell))
        total = total + weight
      }
    }
    entries.sort_by(fn(
      left : (Double, Double, Int),
      right : (Double, Double, Int),
    ) -> Int {
      if left.0 < right.0 {
        -1
      } else if left.0 > right.0 {
        1
      } else if left.2 < right.2 {
        -1
      } else if left.2 > right.2 {
        1
      } else {
        0
      }
    })
    let mut cumulative = 0.0
    let mut selected = entries[entries.length() - 1].0
    for entry in entries {
      cumulative = cumulative + entry.1
      if cumulative >= 0.5 * total {
        selected = entry.0
        break
      }
    }
    center.push(selected)
  }
  center
}

///|
fn slingshot_adv_covariance(
  coordinates : Array[Array[Double]],
  weights : Array[Array[Double]],
  cluster : Int,
  center : Array[Double],
  ridge : Double,
) -> Array[Array[Double]] {
  let dimensions = center.length()
  let output = slingshot_adv_double_matrix(dimensions, dimensions, 0.0)
  let mut total = 0.0
  let mut sum_squares = 0.0
  for cell in 0..<coordinates.length() {
    let weight = weights[cell][cluster]
    total = total + weight
    sum_squares = sum_squares + weight * weight
    for row in 0..<dimensions {
      let left = coordinates[cell][row] - center[row]
      for column in 0..<dimensions {
        output[row][column] = output[row][column] +
          weight * left * (coordinates[cell][column] - center[column])
      }
    }
  }
  let denominator = (total - sum_squares / total).max(1.0)
  for row in 0..<dimensions {
    for column in 0..<dimensions {
      output[row][column] = output[row][column] / denominator
    }
    output[row][row] = output[row][row] + ridge
  }
  output
}

///|
fn slingshot_adv_cluster_summaries(
  coordinates : Array[Array[Double]],
  weights : Array[Array[Double]],
  cluster_names : Array[String],
  config : SlingshotAdvancedConfig,
) -> Array[SlingshotClusterSummary] {
  let output : Array[SlingshotClusterSummary] = []
  for cluster in 0..<cluster_names.length() {
    let center = if config.use_median {
      slingshot_adv_weighted_median(coordinates, weights, cluster)
    } else {
      slingshot_adv_weighted_mean(coordinates, weights, cluster)
    }
    let covariance = slingshot_adv_covariance(
      coordinates,
      weights,
      cluster,
      center,
      config.covariance_ridge,
    )
    let mut size = 0.0
    for cell in 0..<coordinates.length() {
      size = size + weights[cell][cluster]
    }
    output.push(SlingshotClusterSummary::{
      name: cluster_names[cluster],
      center,
      covariance,
      size,
    })
  }
  output
}

///|
fn slingshot_adv_inverse(
  matrix : Array[Array[Double]],
) -> Array[Array[Double]]? {
  let size = matrix.length()
  if size == 0 {
    return None
  }
  let augmented = slingshot_adv_double_matrix(size, size * 2, 0.0)
  for row in 0..<size {
    if matrix[row].length() != size {
      return None
    }
    for column in 0..<size {
      augmented[row][column] = matrix[row][column]
    }
    augmented[row][size + row] = 1.0
  }
  for pivot in 0..<size {
    let mut best = pivot
    for row in (pivot + 1)..<size {
      if augmented[row][pivot].abs() > augmented[best][pivot].abs() {
        best = row
      }
    }
    if augmented[best][pivot].abs() < 1.0e-12 {
      return None
    }
    if best != pivot {
      let temporary = augmented[pivot]
      augmented[pivot] = augmented[best]
      augmented[best] = temporary
    }
    let scale = augmented[pivot][pivot]
    for column in 0..<(size * 2) {
      augmented[pivot][column] = augmented[pivot][column] / scale
    }
    for row in 0..<size {
      if row != pivot {
        let factor = augmented[row][pivot]
        for column in 0..<(size * 2) {
          augmented[row][column] = augmented[row][column] -
            factor * augmented[pivot][column]
        }
      }
    }
  }
  let inverse = slingshot_adv_double_matrix(size, size, 0.0)
  for row in 0..<size {
    for column in 0..<size {
      inverse[row][column] = augmented[row][size + column]
    }
  }
  Some(inverse)
}

///|
fn slingshot_adv_euclidean(
  left : Array[Double],
  right : Array[Double],
) -> Double {
  let mut squared = 0.0
  for dimension in 0..<left.length() {
    let difference = left[dimension] - right[dimension]
    squared = squared + difference * difference
  }
  squared.sqrt()
}

///|
fn slingshot_adv_cluster_distance(
  left : SlingshotClusterSummary,
  right : SlingshotClusterSummary,
  kind : SlingshotClusterDistance,
  ridge : Double,
) -> Double {
  let dimensions = left.center.length()
  if kind == SlingshotCenterEuclidean {
    return slingshot_adv_euclidean(left.center, right.center)
  }
  let delta = Array::make(dimensions, 0.0)
  for dimension in 0..<dimensions {
    delta[dimension] = left.center[dimension] - right.center[dimension]
  }
  if kind == SlingshotScaledDiagonal {
    let mut squared = 0.0
    for dimension in 0..<dimensions {
      let variance = 0.5 *
        (
          left.covariance[dimension][dimension] +
          right.covariance[dimension][dimension]
        ) +
        ridge
      squared = squared + delta[dimension] * delta[dimension] / variance
    }
    return squared.max(0.0).sqrt()
  }
  let pooled = slingshot_adv_double_matrix(dimensions, dimensions, 0.0)
  for row in 0..<dimensions {
    for column in 0..<dimensions {
      pooled[row][column] = 0.5 *
        (left.covariance[row][column] + right.covariance[row][column])
    }
    pooled[row][row] = pooled[row][row] + ridge
  }
  match slingshot_adv_inverse(pooled) {
    Some(inverse) => {
      let mut squared = 0.0
      for row in 0..<dimensions {
        for column in 0..<dimensions {
          squared = squared + delta[row] * inverse[row][column] * delta[column]
        }
      }
      squared.max(0.0).sqrt()
    }
    None =>
      slingshot_adv_cluster_distance(
        left,
        right,
        SlingshotScaledDiagonal,
        ridge,
      )
  }
}

///|
fn slingshot_adv_distance_matrix(
  clusters : Array[SlingshotClusterSummary],
  config : SlingshotAdvancedConfig,
) -> Array[Array[Double]] {
  let output = slingshot_adv_double_matrix(
    clusters.length(),
    clusters.length(),
    0.0,
  )
  for left in 0..<clusters.length() {
    for right in (left + 1)..<clusters.length() {
      let distance = slingshot_adv_cluster_distance(
        clusters[left],
        clusters[right],
        config.distance,
        config.covariance_ridge,
      )
      output[left][right] = distance
      output[right][left] = distance
    }
  }
  output
}

///|
fn slingshot_adv_find(parent : Array[Int], node : Int) -> Int {
  let mut current = node
  while parent[current] != current {
    current = parent[current]
  }
  current
}

///|
fn slingshot_adv_union(
  parent : Array[Int],
  ranks : Array[Int],
  left : Int,
  right : Int,
) -> Unit {
  let left_root = slingshot_adv_find(parent, left)
  let right_root = slingshot_adv_find(parent, right)
  if left_root == right_root {
    return
  }
  if ranks[left_root] < ranks[right_root] {
    parent[left_root] = right_root
  } else if ranks[left_root] > ranks[right_root] {
    parent[right_root] = left_root
  } else {
    parent[right_root] = left_root
    ranks[left_root] = ranks[left_root] + 1
  }
}

///|
fn slingshot_adv_sorted_edges(
  distances : Array[Array[Double]],
) -> Array[SlingshotLineageEdge] {
  let output : Array[SlingshotLineageEdge] = []
  for left in 0..<distances.length() {
    for right in (left + 1)..<distances.length() {
      output.push(SlingshotLineageEdge::{
        from: left,
        to: right,
        distance: distances[left][right],
      })
    }
  }
  output.sort_by(fn(
    left : SlingshotLineageEdge,
    right : SlingshotLineageEdge,
  ) -> Int {
    if left.distance < right.distance {
      -1
    } else if left.distance > right.distance {
      1
    } else if left.from < right.from {
      -1
    } else if left.from > right.from {
      1
    } else if left.to < right.to {
      -1
    } else if left.to > right.to {
      1
    } else {
      0
    }
  })
  output
}

///|
fn slingshot_adv_endpoint_mask(
  cluster_names : Array[String],
  config : SlingshotAdvancedConfig,
) -> Array[Bool] raise SlingshotAdvancedError {
  let output = Array::make(cluster_names.length(), false)
  for name in config.start_clusters {
    if slingshot_adv_label_index(cluster_names, name) < 0 {
      slingshot_adv_fail("slingshot start cluster '" + name + "' was not found")
    }
  }
  for name in config.end_clusters {
    let index = slingshot_adv_label_index(cluster_names, name)
    if index < 0 {
      slingshot_adv_fail("slingshot end cluster '" + name + "' was not found")
    }
    output[index] = true
  }
  output
}

///|
fn slingshot_adv_kruskal(
  distances : Array[Array[Double]],
  endpoints : Array[Bool],
  threshold : Double,
) -> Array[SlingshotLineageEdge] {
  let size = distances.length()
  let parent = Array::make(size, 0)
  let ranks = Array::make(size, 0)
  let degrees = Array::make(size, 0)
  for index in 0..<size {
    parent[index] = index
  }
  let output : Array[SlingshotLineageEdge] = []
  for edge in slingshot_adv_sorted_edges(distances) {
    if threshold >= 0.0 && edge.distance > threshold {
      continue
    }
    if slingshot_adv_find(parent, edge.from) ==
      slingshot_adv_find(parent, edge.to) {
      continue
    }
    if (endpoints[edge.from] && degrees[edge.from] >= 1) ||
      (endpoints[edge.to] && degrees[edge.to] >= 1) {
      continue
    }
    output.push(edge)
    degrees[edge.from] = degrees[edge.from] + 1
    degrees[edge.to] = degrees[edge.to] + 1
    slingshot_adv_union(parent, ranks, edge.from, edge.to)
  }
  output
}

///|
fn slingshot_adv_median(values : Array[Double]) -> Double {
  if values.length() == 0 {
    return 0.0
  }
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
    0.5 * (sorted[middle - 1] + sorted[middle])
  } else {
    sorted[middle]
  }
}

///|
fn slingshot_adv_resolve_omega(
  distances : Array[Array[Double]],
  endpoints : Array[Bool],
  config : SlingshotAdvancedConfig,
) -> Double {
  if config.automatic_omega {
    let unrestricted = slingshot_adv_kruskal(distances, endpoints, -1.0)
    let lengths = unrestricted.map(fn(edge) { edge.distance })
    return slingshot_adv_median(lengths) * config.omega_scale
  }
  config.omega
}

///|
fn slingshot_adv_components(
  nodes : Int,
  edges : Array[SlingshotLineageEdge],
) -> Array[Int] {
  let parent = Array::make(nodes, 0)
  let ranks = Array::make(nodes, 0)
  for node in 0..<nodes {
    parent[node] = node
  }
  for edge in edges {
    slingshot_adv_union(parent, ranks, edge.from, edge.to)
  }
  let root_to_component : Map[Int, Int] = Map([])
  let output = Array::make(nodes, 0)
  let mut next = 0
  for node in 0..<nodes {
    let root = slingshot_adv_find(parent, node)
    match root_to_component.get(root) {
      Some(component) => output[node] = component
      None => {
        root_to_component[root] = next
        output[node] = next
        next = next + 1
      }
    }
  }
  output
}

///|
fn slingshot_adv_adjacency(
  nodes : Int,
  edges : Array[SlingshotLineageEdge],
) -> Array[Array[Int]] {
  let output : Array[Array[Int]] = []
  for _ in 0..<nodes {
    output.push([])
  }
  for edge in edges {
    output[edge.from].push(edge.to)
    output[edge.to].push(edge.from)
  }
  for neighbors in output {
    neighbors.sort()
  }
  output
}

///|
fn slingshot_adv_path(
  adjacency : Array[Array[Int]],
  start : Int,
  end : Int,
) -> Array[Int] {
  let previous = Array::make(adjacency.length(), -1)
  let visited = Array::make(adjacency.length(), false)
  let queue : Array[Int] = [start]
  visited[start] = true
  let mut head = 0
  while head < queue.length() {
    let current = queue[head]
    head = head + 1
    if current == end {
      break
    }
    for neighbor in adjacency[current] {
      if !visited[neighbor] {
        visited[neighbor] = true
        previous[neighbor] = current
        queue.push(neighbor)
      }
    }
  }
  if !visited[end] {
    return []
  }
  let reverse : Array[Int] = []
  let mut current = end
  while current >= 0 {
    reverse.push(current)
    if current == start {
      break
    }
    current = previous[current]
  }
  let output : Array[Int] = []
  let mut index = reverse.length() - 1
  while index >= 0 {
    output.push(reverse[index])
    index = index - 1
  }
  output
}

///|
fn slingshot_adv_choose_root(
  nodes : Array[Int],
  adjacency : Array[Array[Int]],
  cluster_names : Array[String],
  starts : Array[String],
) -> Int {
  for name in starts {
    let index = slingshot_adv_label_index(cluster_names, name)
    if nodes.contains(index) {
      return index
    }
  }
  if nodes.length() == 1 {
    return nodes[0]
  }
  let leaves = nodes.filter(fn(node) { adjacency[node].length() <= 1 })
  let candidates = if leaves.length() > 0 { leaves } else { nodes }
  let mut best = candidates[0]
  let mut best_mean = -1.0
  for candidate in candidates {
    let mut total = 0.0
    let mut comparisons = 0
    for leaf in leaves {
      if leaf != candidate {
        total = total +
          slingshot_adv_path(adjacency, candidate, leaf).length().to_double()
        comparisons = comparisons + 1
      }
    }
    let mean = if comparisons > 0 {
      total / comparisons.to_double()
    } else {
      0.0
    }
    if mean > best_mean || (mean == best_mean && candidate < best) {
      best = candidate
      best_mean = mean
    }
  }
  best
}

///|
fn slingshot_adv_lineages(
  cluster_names : Array[String],
  edges : Array[SlingshotLineageEdge],
  components : Array[Int],
  config : SlingshotAdvancedConfig,
) -> (Array[Int], Array[Array[Int]]) {
  let adjacency = slingshot_adv_adjacency(cluster_names.length(), edges)
  let mut component_count = 0
  for component in components {
    component_count = component_count.max(component + 1)
  }
  let roots : Array[Int] = []
  let lineages : Array[Array[Int]] = []
  for component in 0..<component_count {
    let nodes : Array[Int] = []
    for node in 0..<components.length() {
      if components[node] == component {
        nodes.push(node)
      }
    }
    let root = slingshot_adv_choose_root(
      nodes,
      adjacency,
      cluster_names,
      config.start_clusters,
    )
    roots.push(root)
    if nodes.length() == 1 {
      lineages.push([root])
      continue
    }
    let leaves = nodes.filter(fn(node) {
      node != root && adjacency[node].length() <= 1
    })
    for leaf in leaves {
      let name = cluster_names[leaf]
      let is_other_start = config.start_clusters.contains(name) &&
        !config.end_clusters.contains(name)
      if !is_other_start {
        let path = slingshot_adv_path(adjacency, root, leaf)
        if path.length() > 0 {
          lineages.push(path)
        }
      }
    }
  }
  lineages.sort_by(fn(left : Array[Int], right : Array[Int]) -> Int {
    if left.length() > right.length() {
      -1
    } else if left.length() < right.length() {
      1
    } else {
      for index in 0..<left.length() {
        if left[index] < right[index] {
          return -1
        } else if left[index] > right[index] {
          return 1
        }
      }
      0
    }
  })
  (roots, lineages)
}

///|
fn slingshot_adv_initial_weights(
  cluster_weights : Array[Array[Double]],
  lineages : Array[Array[Int]],
) -> Array[Array[Double]] {
  let output = slingshot_adv_double_matrix(
    cluster_weights.length(),
    lineages.length(),
    0.0,
  )
  for cell in 0..<cluster_weights.length() {
    for lineage in 0..<lineages.length() {
      for cluster in lineages[lineage] {
        output[cell][lineage] = output[cell][lineage] +
          cluster_weights[cell][cluster]
      }
      output[cell][lineage] = output[cell][lineage].min(1.0)
    }
  }
  output
}

///|
pub fn slingshot_infer_weighted_lineages(
  coordinates : Array[Array[Double]],
  cluster_weights : Array[Array[Double]],
  cluster_names : Array[String],
  config : SlingshotAdvancedConfig,
) -> SlingshotLineageModel raise SlingshotAdvancedError {
  let (cells, _) = slingshot_adv_shape(coordinates, "slingshot coordinates")
  slingshot_adv_validate_cluster_names(cluster_names)
  let normalized = slingshot_adv_normalize_weights(
    cluster_weights, cluster_names, cells,
  )
  let clusters = slingshot_adv_cluster_summaries(
    coordinates, normalized, cluster_names, config,
  )
  let distances = slingshot_adv_distance_matrix(clusters, config)
  let endpoints = slingshot_adv_endpoint_mask(cluster_names, config)
  let omega_threshold = slingshot_adv_resolve_omega(
    distances, endpoints, config,
  )
  let edges = slingshot_adv_kruskal(distances, endpoints, omega_threshold)
  let components = slingshot_adv_components(cluster_names.length(), edges)
  let (roots, lineages) = slingshot_adv_lineages(
    cluster_names, edges, components, config,
  )
  let initial_weights = slingshot_adv_initial_weights(normalized, lineages)
  SlingshotLineageModel::{
    coordinates: slingshot_adv_copy_matrix(coordinates),
    cluster_names: cluster_names.copy(),
    cluster_weights: normalized,
    clusters,
    distance_matrix: distances,
    edges,
    components,
    roots,
    lineages,
    initial_weights,
    omega_threshold,
  }
}

///|
pub fn slingshot_infer_lineages(
  coordinates : Array[Array[Double]],
  cluster_labels : Array[String],
  config : SlingshotAdvancedConfig,
) -> SlingshotLineageModel raise SlingshotAdvancedError {
  let (cells, _) = slingshot_adv_shape(coordinates, "slingshot coordinates")
  if cluster_labels.length() != cells {
    slingshot_adv_fail("slingshot labels must have one value per cell")
  }
  let cluster_names = slingshot_adv_unique_labels(cluster_labels)
  slingshot_adv_validate_cluster_names(cluster_names)
  let weights = slingshot_adv_double_matrix(cells, cluster_names.length(), 0.0)
  for cell in 0..<cells {
    if cluster_labels[cell].length() == 0 {
      slingshot_adv_fail("slingshot cluster labels cannot be empty")
    }
    if cluster_labels[cell] != "-1" {
      let cluster = slingshot_adv_label_index(
        cluster_names,
        cluster_labels[cell],
      )
      weights[cell][cluster] = 1.0
    }
  }
  slingshot_infer_weighted_lineages(coordinates, weights, cluster_names, config)
}

///|
fn slingshot_adv_polyline_length(points : Array[Array[Double]]) -> Double {
  let mut total = 0.0
  for index in 1..<points.length() {
    total = total + slingshot_adv_euclidean(points[index - 1], points[index])
  }
  total
}

///|
fn slingshot_adv_interpolate(
  left : Array[Double],
  right : Array[Double],
  fraction : Double,
) -> Array[Double] {
  let output = Array::make(left.length(), 0.0)
  for dimension in 0..<left.length() {
    output[dimension] = left[dimension] +
      fraction * (right[dimension] - left[dimension])
  }
  output
}

///|
fn slingshot_adv_resample(
  points : Array[Array[Double]],
  count : Int,
) -> Array[Array[Double]] {
  if count == 1 {
    return [points[0].copy()]
  }
  if points.length() == 1 {
    let output : Array[Array[Double]] = []
    for _ in 0..<count {
      output.push(points[0].copy())
    }
    return output
  }
  let cumulative = Array::make(points.length(), 0.0)
  for index in 1..<points.length() {
    cumulative[index] = cumulative[index - 1] +
      slingshot_adv_euclidean(points[index - 1], points[index])
  }
  let total = cumulative[cumulative.length() - 1]
  if total <= 0.0 {
    let output : Array[Array[Double]] = []
    for _ in 0..<count {
      output.push(points[0].copy())
    }
    return output
  }
  let output : Array[Array[Double]] = []
  for target_index in 0..<count {
    if target_index == 0 {
      output.push(points[0].copy())
      continue
    }
    if target_index == count - 1 {
      output.push(points[points.length() - 1].copy())
      continue
    }
    let target = total * target_index.to_double() / (count - 1).to_double()
    let mut segment = 0
    while segment + 2 < cumulative.length() && cumulative[segment + 1] < target {
      segment = segment + 1
    }
    let width = cumulative[segment + 1] - cumulative[segment]
    let fraction = if width > 0.0 {
      (target - cumulative[segment]) / width
    } else {
      0.0
    }
    output.push(
      slingshot_adv_interpolate(points[segment], points[segment + 1], fraction),
    )
  }
  output
}

///|
fn slingshot_adv_first_pc(covariance : Array[Array[Double]]) -> Array[Double] {
  let dimensions = covariance.length()
  let mut vector = Array::make(dimensions, 1.0)
  let norm0 = dimensions.to_double().sqrt()
  for dimension in 0..<dimensions {
    vector[dimension] = vector[dimension] / norm0
  }
  for _ in 0..<100 {
    let next = Array::make(dimensions, 0.0)
    for row in 0..<dimensions {
      for column in 0..<dimensions {
        next[row] = next[row] + covariance[row][column] * vector[column]
      }
    }
    let mut norm = 0.0
    for value in next {
      norm = norm + value * value
    }
    norm = norm.sqrt()
    if norm <= 1.0e-12 {
      break
    }
    let mut change = 0.0
    for dimension in 0..<dimensions {
      next[dimension] = next[dimension] / norm
      change = change + (next[dimension] - vector[dimension]).abs()
    }
    vector = next
    if change < 1.0e-10 {
      break
    }
  }
  vector
}

///|
fn slingshot_adv_dot(left : Array[Double], right : Array[Double]) -> Double {
  let mut output = 0.0
  for index in 0..<left.length() {
    output = output + left[index] * right[index]
  }
  output
}

///|
fn slingshot_adv_endpoint_extension(
  model : SlingshotLineageModel,
  lineage : Array[Int],
  config : SlingshotAdvancedConfig,
) -> Array[Array[Double]] {
  let points : Array[Array[Double]] = []
  for cluster in lineage {
    points.push(model.clusters[cluster].center.copy())
  }
  if lineage.length() == 1 {
    let cluster = model.clusters[lineage[0]]
    let direction = slingshot_adv_first_pc(cluster.covariance)
    let scale = cluster.covariance[0][0].max(config.covariance_ridge).sqrt()
    let before = cluster.center.copy()
    let after = cluster.center.copy()
    for dimension in 0..<direction.length() {
      before[dimension] = before[dimension] - direction[dimension] * scale
      after[dimension] = after[dimension] + direction[dimension] * scale
    }
    return [before, cluster.center.copy(), after]
  }
  if config.extension == SlingshotExtendNone {
    return points
  }
  let first = points[0]
  let second = points[1]
  let penultimate = points[points.length() - 2]
  let last = points[points.length() - 1]
  let start_direction = Array::make(first.length(), 0.0)
  let end_direction = Array::make(first.length(), 0.0)
  if config.extension == SlingshotExtendLine {
    for dimension in 0..<first.length() {
      start_direction[dimension] = first[dimension] - second[dimension]
      end_direction[dimension] = last[dimension] - penultimate[dimension]
    }
  } else {
    let start_pc = slingshot_adv_first_pc(model.clusters[lineage[0]].covariance)
    let end_pc = slingshot_adv_first_pc(
      model.clusters[lineage[lineage.length() - 1]].covariance,
    )
    let start_leg = Array::make(first.length(), 0.0)
    let end_leg = Array::make(first.length(), 0.0)
    for dimension in 0..<first.length() {
      start_leg[dimension] = first[dimension] - second[dimension]
      end_leg[dimension] = last[dimension] - penultimate[dimension]
    }
    let start_sign = if slingshot_adv_dot(start_pc, start_leg) >= 0.0 {
      1.0
    } else {
      -1.0
    }
    let end_sign = if slingshot_adv_dot(end_pc, end_leg) >= 0.0 {
      1.0
    } else {
      -1.0
    }
    let start_scale = slingshot_adv_euclidean(first, second)
    let end_scale = slingshot_adv_euclidean(last, penultimate)
    for dimension in 0..<first.length() {
      start_direction[dimension] = start_sign *
        start_pc[dimension] *
        start_scale
      end_direction[dimension] = end_sign * end_pc[dimension] * end_scale
    }
  }
  let before = first.copy()
  let after = last.copy()
  for dimension in 0..<first.length() {
    before[dimension] = before[dimension] + start_direction[dimension]
    after[dimension] = after[dimension] + end_direction[dimension]
  }
  let output : Array[Array[Double]] = [before]
  for point in points {
    output.push(point)
  }
  output.push(after)
  output
}

///|
fn slingshot_adv_project_point(
  point : Array[Double],
  curve : Array[Array[Double]],
) -> (Double, Double) {
  if curve.length() == 1 {
    return (0.0, slingshot_adv_euclidean(point, curve[0]))
  }
  let cumulative = Array::make(curve.length(), 0.0)
  for index in 1..<curve.length() {
    cumulative[index] = cumulative[index - 1] +
      slingshot_adv_euclidean(curve[index - 1], curve[index])
  }
  let mut best_distance = 1.0e300
  let mut best_time = 0.0
  for segment in 0..<(curve.length() - 1) {
    let left = curve[segment]
    let right = curve[segment + 1]
    let mut length_squared = 0.0
    let mut numerator = 0.0
    for dimension in 0..<point.length() {
      let direction = right[dimension] - left[dimension]
      length_squared = length_squared + direction * direction
      numerator = numerator + (point[dimension] - left[dimension]) * direction
    }
    let fraction = if length_squared > 0.0 {
      (numerator / length_squared).max(0.0).min(1.0)
    } else {
      0.0
    }
    let projected = slingshot_adv_interpolate(left, right, fraction)
    let distance = slingshot_adv_euclidean(point, projected)
    let segment_length = length_squared.sqrt()
    let time = cumulative[segment] + fraction * segment_length
    if distance < best_distance ||
      (distance == best_distance && time < best_time) {
      best_distance = distance
      best_time = time
    }
  }
  (best_time, best_distance)
}

///|
fn slingshot_adv_project_all(
  coordinates : Array[Array[Double]],
  curve : Array[Array[Double]],
) -> (Array[Double], Array[Double]) {
  let pseudotime = Array::make(coordinates.length(), 0.0)
  let distances = Array::make(coordinates.length(), 0.0)
  for cell in 0..<coordinates.length() {
    let (time, distance) = slingshot_adv_project_point(coordinates[cell], curve)
    pseudotime[cell] = time
    distances[cell] = distance
  }
  (pseudotime, distances)
}

///|
fn slingshot_adv_local_linear(
  times : Array[Double],
  values : Array[Double],
  weights : Array[Double],
  target : Double,
  bandwidth : Double,
) -> Double {
  let mut s0 = 0.0
  let mut s1 = 0.0
  let mut s2 = 0.0
  let mut t0 = 0.0
  let mut t1 = 0.0
  for index in 0..<times.length() {
    if weights[index] > 0.0 {
      let centered = times[index] - target
      let scaled = centered / bandwidth
      let kernel = @math.exp(-0.5 * scaled * scaled) * weights[index]
      s0 = s0 + kernel
      s1 = s1 + kernel * centered
      s2 = s2 + kernel * centered * centered
      t0 = t0 + kernel * values[index]
      t1 = t1 + kernel * centered * values[index]
    }
  }
  let determinant = s0 * s2 - s1 * s1
  if determinant.abs() > 1.0e-12 {
    (t0 * s2 - t1 * s1) / determinant
  } else if s0 > 0.0 {
    t0 / s0
  } else {
    0.0
  }
}

///|
fn slingshot_adv_smooth_curve(
  coordinates : Array[Array[Double]],
  pseudotime : Array[Double],
  weights : Array[Double],
  points : Int,
  span : Double,
) -> Array[Array[Double]] {
  let dimensions = coordinates[0].length()
  let mut minimum = 1.0e300
  let mut maximum = -1.0e300
  for cell in 0..<coordinates.length() {
    if weights[cell] > 0.0 {
      minimum = minimum.min(pseudotime[cell])
      maximum = maximum.max(pseudotime[cell])
    }
  }
  if minimum > maximum {
    minimum = 0.0
    maximum = 1.0
  }
  let range = (maximum - minimum).max(1.0e-8)
  let bandwidth = (range * span).max(1.0e-8)
  let output : Array[Array[Double]] = []
  for index in 0..<points {
    let target = minimum + range * index.to_double() / (points - 1).to_double()
    let point = Array::make(dimensions, 0.0)
    for dimension in 0..<dimensions {
      let values = coordinates.map(fn(row) { row[dimension] })
      point[dimension] = slingshot_adv_local_linear(
        pseudotime, values, weights, target, bandwidth,
      )
    }
    output.push(point)
  }
  output
}

///|
fn slingshot_adv_distance_ranks(
  distances : Array[Array[Double]],
  weights : Array[Array[Double]],
) -> Array[Array[Double]] {
  let cells = distances.length()
  let lineages = if cells > 0 { distances[0].length() } else { 0 }
  let entries : Array[(Double, Int, Int, Double)] = []
  let mut total = 0.0
  for cell in 0..<cells {
    let mut row_total = 0.0
    for lineage in 0..<lineages {
      row_total = row_total + weights[cell][lineage]
    }
    if row_total > 0.0 {
      for lineage in 0..<lineages {
        if weights[cell][lineage] > 0.0 {
          let probability = weights[cell][lineage] / row_total
          entries.push((distances[cell][lineage], cell, lineage, probability))
          total = total + probability
        }
      }
    }
  }
  entries.sort_by(fn(
    left : (Double, Int, Int, Double),
    right : (Double, Int, Int, Double),
  ) -> Int {
    if left.0 < right.0 {
      -1
    } else if left.0 > right.0 {
      1
    } else if left.1 < right.1 {
      -1
    } else if left.1 > right.1 {
      1
    } else if left.2 < right.2 {
      -1
    } else if left.2 > right.2 {
      1
    } else {
      0
    }
  })
  let output = slingshot_adv_double_matrix(cells, lineages, 1.0)
  let mut cumulative = 0.0
  for entry in entries {
    cumulative = cumulative + entry.3
    output[entry.1][entry.2] = if total > 0.0 {
      cumulative / total
    } else {
      1.0
    }
  }
  output
}

///|
fn slingshot_adv_update_weights(
  distances : Array[Array[Double]],
  current : Array[Array[Double]],
  reweight : Bool,
  reassign : Bool,
) -> Array[Array[Double]] {
  if !reweight && !reassign {
    return slingshot_adv_copy_matrix(current)
  }
  let cells = distances.length()
  let lineages = distances[0].length()
  let ranks = slingshot_adv_distance_ranks(distances, current)
  let output = slingshot_adv_copy_matrix(current)
  if reweight {
    for cell in 0..<cells {
      let mut maximum = 0.0
      for lineage in 0..<lineages {
        if current[cell][lineage] > 0.0 {
          let candidate = 1.0 - ranks[cell][lineage] * ranks[cell][lineage]
          output[cell][lineage] = candidate.max(0.0)
          maximum = maximum.max(output[cell][lineage])
        } else {
          output[cell][lineage] = 0.0
        }
      }
      if maximum > 0.0 {
        for lineage in 0..<lineages {
          output[cell][lineage] = output[cell][lineage] / maximum
        }
      }
    }
  }
  if reassign {
    for cell in 0..<cells {
      let mut nearest = 0
      for lineage in 1..<lineages {
        if distances[cell][lineage] < distances[cell][nearest] {
          nearest = lineage
        }
      }
      for lineage in 0..<lineages {
        if ranks[cell][lineage] < 0.5 {
          output[cell][lineage] = 1.0
        } else if ranks[cell][lineage] > 0.9 && output[cell][lineage] < 0.1 {
          output[cell][lineage] = 0.0
        }
      }
      let mut total = 0.0
      for lineage in 0..<lineages {
        total = total + output[cell][lineage]
      }
      if total <= 0.0 {
        output[cell][nearest] = 1.0
      }
    }
  }
  output
}

///|
fn slingshot_adv_shared_prefix(left : Array[Int], right : Array[Int]) -> Int {
  let limit = left.length().min(right.length())
  let mut shared = 0
  while shared < limit && left[shared] == right[shared] {
    shared = shared + 1
  }
  shared
}

///|
fn slingshot_adv_shrink_curves(
  curves : Array[Array[Array[Double]]],
  lineages : Array[Array[Int]],
  shrink : Double,
) -> Array[Array[Array[Double]]] {
  if shrink <= 0.0 || curves.length() <= 1 {
    return curves.map(fn(curve) { slingshot_adv_copy_matrix(curve) })
  }
  let output = curves.map(fn(curve) { slingshot_adv_copy_matrix(curve) })
  for lineage in 0..<curves.length() {
    for point in 0..<curves[lineage].length() {
      let fraction = point.to_double() /
        (curves[lineage].length() - 1).to_double()
      let partners : Array[Int] = [lineage]
      let thresholds : Array[Double] = [1.0]
      for other in 0..<curves.length() {
        if other != lineage {
          let shared = slingshot_adv_shared_prefix(
            lineages[lineage],
            lineages[other],
          )
          if shared > 0 {
            let left_denominator = (lineages[lineage].length() - 1).max(1)
            let right_denominator = (lineages[other].length() - 1).max(1)
            let threshold = ((shared - 1).to_double() /
            left_denominator.to_double()).min(
              (shared - 1).to_double() / right_denominator.to_double(),
            )
            if fraction <= threshold.max(1.0e-12) {
              partners.push(other)
              thresholds.push(threshold.max(1.0e-12))
            }
          }
        }
      }
      if partners.length() > 1 {
        let average = Array::make(curves[lineage][point].length(), 0.0)
        for partner in partners {
          let partner_point = (fraction *
            (curves[partner].length() - 1).to_double())
            .round()
            .to_int()
            .max(0)
            .min(curves[partner].length() - 1)
          for dimension in 0..<average.length() {
            average[dimension] = average[dimension] +
              curves[partner][partner_point][dimension]
          }
        }
        for dimension in 0..<average.length() {
          average[dimension] = average[dimension] /
            partners.length().to_double()
        }
        let threshold = thresholds[1]
        let taper = 0.5 * (1.0 + @math.cos(@math.PI * fraction / threshold))
        let amount = shrink * taper.max(0.0).min(1.0)
        for dimension in 0..<average.length() {
          output[lineage][point][dimension] = amount * average[dimension] +
            (1.0 - amount) * curves[lineage][point][dimension]
        }
      }
    }
  }
  output
}

///|
fn slingshot_adv_curves_to_distance_matrix(
  curves : Array[SlingshotAdvancedCurve],
) -> Array[Array[Double]] {
  let cells = curves[0].distances.length()
  let output = slingshot_adv_double_matrix(cells, curves.length(), 0.0)
  for lineage in 0..<curves.length() {
    for cell in 0..<cells {
      output[cell][lineage] = curves[lineage].distances[cell]
    }
  }
  output
}

///|
fn slingshot_adv_total_distance(
  distances : Array[Array[Double]],
  weights : Array[Array[Double]],
) -> Double {
  let mut output = 0.0
  for cell in 0..<distances.length() {
    for lineage in 0..<distances[cell].length() {
      if weights[cell][lineage] > 0.0 {
        output = output + distances[cell][lineage]
      }
    }
  }
  output
}

///|
fn slingshot_adv_build_output(
  model : SlingshotLineageModel,
  curves : Array[SlingshotAdvancedCurve],
  weights : Array[Array[Double]],
  iterations : Int,
  converged : Bool,
  total_distance : Double,
  config : SlingshotAdvancedConfig,
) -> SlingshotAdvancedResult {
  let cells = model.coordinates.length()
  let pseudotime = slingshot_adv_optional_matrix(cells, curves.length())
  let average = Array::make(cells, 0.0)
  for cell in 0..<cells {
    let mut weight_total = 0.0
    for lineage in 0..<curves.length() {
      if weights[cell][lineage] > 0.0 {
        pseudotime[cell][lineage] = Some(curves[lineage].pseudotime[cell])
        average[cell] = average[cell] +
          weights[cell][lineage] * curves[lineage].pseudotime[cell]
        weight_total = weight_total + weights[cell][lineage]
      }
    }
    if weight_total > 0.0 {
      average[cell] = average[cell] / weight_total
    }
  }
  let branch_ids = slingshot_adv_branch_ids_from_weights(weights, -1.0)
  SlingshotAdvancedResult::{
    lineage_model: model,
    curves,
    pseudotime,
    weights,
    average_pseudotime: average,
    branch_ids,
    iterations,
    converged,
    total_distance,
    config,
  }
}

///|
pub fn slingshot_fit_simultaneous_curves(
  model : SlingshotLineageModel,
  config : SlingshotAdvancedConfig,
) -> SlingshotAdvancedResult raise SlingshotAdvancedError {
  if model.lineages.length() == 0 {
    slingshot_adv_fail("slingshot lineage model contains no lineages")
  }
  let raw_curves : Array[Array[Array[Double]]] = []
  for lineage in model.lineages {
    raw_curves.push(
      slingshot_adv_resample(
        slingshot_adv_endpoint_extension(model, lineage, config),
        config.curve_points,
      ),
    )
  }
  let curves : Array[SlingshotAdvancedCurve] = []
  for lineage in 0..<model.lineages.length() {
    let (pseudotime, distances) = slingshot_adv_project_all(
      model.coordinates,
      raw_curves[lineage],
    )
    curves.push(SlingshotAdvancedCurve::{
      name: "Lineage" + (lineage + 1).to_string(),
      lineage: model.lineages[lineage].copy(),
      points: raw_curves[lineage],
      pseudotime,
      distances,
      weights: model.initial_weights.map(fn(row) { row[lineage] }),
      length: slingshot_adv_polyline_length(raw_curves[lineage]),
    })
  }
  let mut weights = slingshot_adv_copy_matrix(model.initial_weights)
  let mut distance_matrix = slingshot_adv_curves_to_distance_matrix(curves)
  let mut total_distance = slingshot_adv_total_distance(
    distance_matrix, weights,
  )
  let mut converged = false
  let mut iterations = 0
  while iterations < config.max_iterations && !converged {
    iterations = iterations + 1
    let previous = total_distance
    weights = slingshot_adv_update_weights(
      distance_matrix,
      weights,
      config.reweight,
      config.reassign,
    )
    let smoothed : Array[Array[Array[Double]]] = []
    for lineage in 0..<curves.length() {
      smoothed.push(
        slingshot_adv_smooth_curve(
          model.coordinates,
          curves[lineage].pseudotime,
          weights.map(fn(row) { row[lineage] }),
          config.curve_points,
          config.smoother_span,
        ),
      )
    }
    let shrunk = slingshot_adv_shrink_curves(
      smoothed,
      model.lineages,
      config.shrink,
    )
    for lineage in 0..<curves.length() {
      let (pseudotime, distances) = slingshot_adv_project_all(
        model.coordinates,
        shrunk[lineage],
      )
      curves[lineage] = SlingshotAdvancedCurve::{
        name: curves[lineage].name,
        lineage: curves[lineage].lineage,
        points: shrunk[lineage],
        pseudotime,
        distances,
        weights: weights.map(fn(row) { row[lineage] }),
        length: slingshot_adv_polyline_length(shrunk[lineage]),
      }
    }
    distance_matrix = slingshot_adv_curves_to_distance_matrix(curves)
    total_distance = slingshot_adv_total_distance(distance_matrix, weights)
    let scale = previous.abs().max(1.0e-12)
    converged = (previous - total_distance).abs() <= config.tolerance * scale
  }
  weights = slingshot_adv_update_weights(
    distance_matrix,
    weights,
    config.reweight,
    config.reassign,
  )
  for lineage in 0..<curves.length() {
    curves[lineage] = SlingshotAdvancedCurve::{
      name: curves[lineage].name,
      lineage: curves[lineage].lineage,
      points: curves[lineage].points,
      pseudotime: curves[lineage].pseudotime,
      distances: curves[lineage].distances,
      weights: weights.map(fn(row) { row[lineage] }),
      length: curves[lineage].length,
    }
  }
  slingshot_adv_build_output(
    model, curves, weights, iterations, converged, total_distance, config,
  )
}

///|
pub fn slingshot_advanced(
  coordinates : Array[Array[Double]],
  cluster_labels : Array[String],
  config : SlingshotAdvancedConfig,
) -> SlingshotAdvancedResult raise SlingshotAdvancedError {
  slingshot_fit_simultaneous_curves(
    slingshot_infer_lineages(coordinates, cluster_labels, config),
    config,
  )
}

///|
pub fn slingshot_advanced_weighted(
  coordinates : Array[Array[Double]],
  cluster_weights : Array[Array[Double]],
  cluster_names : Array[String],
  config : SlingshotAdvancedConfig,
) -> SlingshotAdvancedResult raise SlingshotAdvancedError {
  slingshot_fit_simultaneous_curves(
    slingshot_infer_weighted_lineages(
      coordinates, cluster_weights, cluster_names, config,
    ),
    config,
  )
}

///|
fn slingshot_adv_branch_ids_from_weights(
  weights : Array[Array[Double]],
  threshold : Double,
) -> Array[String] {
  let lineages = if weights.length() > 0 { weights[0].length() } else { 0 }
  let cutoff = if threshold < 0.0 && lineages > 0 {
    1.0 / lineages.to_double()
  } else {
    threshold
  }
  let output : Array[String] = []
  for row in weights {
    let mut label = ""
    for lineage in 0..<row.length() {
      if row[lineage] >= cutoff {
        if label.length() > 0 {
          label = label + ","
        }
        label = label + (lineage + 1).to_string()
      }
    }
    output.push(label)
  }
  output
}

///|
pub fn slingshot_branch_ids(
  result : SlingshotAdvancedResult,
  threshold? : Double = -1.0,
) -> Array[String] raise SlingshotAdvancedError {
  if !slingshot_adv_is_finite(threshold) || threshold > 1.0 || threshold < -1.0 {
    slingshot_adv_fail("slingshot branch threshold must be in [0, 1]")
  }
  slingshot_adv_branch_ids_from_weights(result.weights, threshold)
}

///|
pub fn slingshot_curve_weight_probabilities(
  result : SlingshotAdvancedResult,
) -> Array[Array[Double]] {
  let output = slingshot_adv_copy_matrix(result.weights)
  for row in output {
    let mut total = 0.0
    for value in row {
      total = total + value
    }
    if total > 0.0 {
      for lineage in 0..<row.length() {
        row[lineage] = row[lineage] / total
      }
    }
  }
  output
}

///|
fn slingshot_adv_quantile(
  values : Array[Double],
  probability : Double,
) -> Double {
  if values.length() == 0 {
    return 1.0
  }
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
  let position = probability * (sorted.length() - 1).to_double()
  let lower = position.floor().to_int()
  let upper = position.ceil().to_int()
  let fraction = position - lower.to_double()
  sorted[lower] * (1.0 - fraction) + sorted[upper] * fraction
}

///|
pub fn slingshot_predict(
  result : SlingshotAdvancedResult,
  new_coordinates : Array[Array[Double]],
) -> SlingshotProjection raise SlingshotAdvancedError {
  let (cells, dimensions) = slingshot_adv_shape(
    new_coordinates, "slingshot prediction coordinates",
  )
  if dimensions != result.lineage_model.coordinates[0].length() {
    slingshot_adv_fail(
      "slingshot prediction dimensions must match the fitted model",
    )
  }
  let lineages = result.curves.length()
  let times = slingshot_adv_double_matrix(cells, lineages, 0.0)
  let distances = slingshot_adv_double_matrix(cells, lineages, 0.0)
  let weights = slingshot_adv_double_matrix(cells, lineages, 0.0)
  let thresholds = Array::make(lineages, 1.0)
  for lineage in 0..<lineages {
    let active : Array[Double] = []
    for cell in 0..<result.weights.length() {
      if result.weights[cell][lineage] > 0.0 {
        active.push(result.curves[lineage].distances[cell])
      }
    }
    thresholds[lineage] = slingshot_adv_quantile(active, 0.9).max(1.0e-8)
    for cell in 0..<cells {
      let (time, distance) = slingshot_adv_project_point(
        new_coordinates[cell],
        result.curves[lineage].points,
      )
      times[cell][lineage] = time
      distances[cell][lineage] = distance
      let scaled = distance / thresholds[lineage]
      weights[cell][lineage] = (1.0 - scaled * scaled).max(0.0)
    }
  }
  for cell in 0..<cells {
    let mut maximum = 0.0
    let mut nearest = 0
    for lineage in 0..<lineages {
      maximum = maximum.max(weights[cell][lineage])
      if distances[cell][lineage] < distances[cell][nearest] {
        nearest = lineage
      }
    }
    if maximum > 0.0 {
      for lineage in 0..<lineages {
        weights[cell][lineage] = weights[cell][lineage] / maximum
      }
    } else {
      weights[cell][nearest] = 1.0
    }
  }
  let pseudotime = slingshot_adv_optional_matrix(cells, lineages)
  let average = Array::make(cells, 0.0)
  for cell in 0..<cells {
    let mut total = 0.0
    for lineage in 0..<lineages {
      if weights[cell][lineage] > 0.0 {
        pseudotime[cell][lineage] = Some(times[cell][lineage])
        average[cell] = average[cell] +
          weights[cell][lineage] * times[cell][lineage]
        total = total + weights[cell][lineage]
      }
    }
    if total > 0.0 {
      average[cell] = average[cell] / total
    }
  }
  SlingshotProjection::{
    pseudotime,
    weights,
    distances,
    average_pseudotime: average,
    branch_ids: slingshot_adv_branch_ids_from_weights(weights, -1.0),
  }
}

///|
fn slingshot_adv_copy_sce(
  source : SingleCellExperiment,
) -> SingleCellExperiment {
  let assays : Map[String, Array[Array[Double]]] = Map([])
  for key in source.assays.keys() {
    assays[key] = slingshot_adv_copy_matrix(source.assays[key])
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
    reduced_dims[key] = slingshot_adv_copy_matrix(source.reduced_dims[key])
  }
  let alternatives : Map[String, SingleCellExperiment] = Map([])
  for key in source.alternative_experiments.keys() {
    alternatives[key] = slingshot_adv_copy_sce(
      source.alternative_experiments[key],
    )
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
pub fn slingshot_advanced_sce(
  experiment : SingleCellExperiment,
  config : SlingshotAdvancedConfig,
  reduced_dim_name? : String = "PCA",
  cluster_column? : String = "cluster",
  output_prefix? : String = "slingshot",
) -> SlingshotAdvancedSceOutput raise SlingshotAdvancedError {
  if reduced_dim_name.length() == 0 ||
    cluster_column.length() == 0 ||
    output_prefix.length() == 0 {
    slingshot_adv_fail("slingshot SCE names must be non-empty")
  }
  let coordinates = match experiment.reduced_dims.get(reduced_dim_name) {
    Some(value) => value
    None => {
      slingshot_adv_fail(
        "slingshot reduced dimension '" + reduced_dim_name + "' was not found",
      )
      []
    }
  }
  let labels = match experiment.col_data.get(cluster_column) {
    Some(value) => value
    None => {
      slingshot_adv_fail(
        "slingshot cluster column '" + cluster_column + "' was not found",
      )
      []
    }
  }
  if coordinates.length() != experiment.col_names.length() ||
    labels.length() != experiment.col_names.length() {
    slingshot_adv_fail(
      "slingshot SCE coordinates and labels must match the cell count",
    )
  }
  let result = slingshot_advanced(coordinates, labels, config)
  let output = slingshot_adv_copy_sce(experiment)
  output.col_data[output_prefix + ".branch"] = result.branch_ids.copy()
  output.col_data[output_prefix + ".pseudotime"] = result.average_pseudotime.map(fn(
      value,
    ) {
      value.to_string()
    },
  )
  let pseudotime_numeric = slingshot_adv_double_matrix(
    result.pseudotime.length(),
    result.curves.length(),
    -1.0,
  )
  for cell in 0..<result.pseudotime.length() {
    for lineage in 0..<result.curves.length() {
      match result.pseudotime[cell][lineage] {
        Some(value) => pseudotime_numeric[cell][lineage] = value
        None => ()
      }
    }
  }
  output.reduced_dims[output_prefix + ".pseudotime"] = pseudotime_numeric
  output.reduced_dims[output_prefix + ".weights"] = slingshot_adv_copy_matrix(
    result.weights,
  )
  output.metadata[output_prefix + ".reduced_dim"] = reduced_dim_name
  output.metadata[output_prefix + ".cluster_column"] = cluster_column
  output.metadata[output_prefix + ".lineages"] = result.curves
    .length()
    .to_string()
  output.metadata[output_prefix + ".converged"] = result.converged.to_string()
  SlingshotAdvancedSceOutput::{ experiment: output, result, }
}

///|
pub fn slingshot_advanced_summary(result : SlingshotAdvancedResult) -> String {
  let mut occupied = 0
  for row in result.lineage_model.cluster_weights {
    let mut total = 0.0
    for weight in row {
      total = total + weight
    }
    if total > 0.0 {
      occupied = occupied + 1
    }
  }
  "Slingshot 2.21.0 advanced result\n" +
  "cells: " +
  result.lineage_model.coordinates.length().to_string() +
  "\nclustered cells: " +
  occupied.to_string() +
  "\nclusters: " +
  result.lineage_model.cluster_names.length().to_string() +
  "\nlineages: " +
  result.curves.length().to_string() +
  "\nforest components: " +
  result.lineage_model.roots.length().to_string() +
  "\niterations: " +
  result.iterations.to_string() +
  "\nconverged: " +
  result.converged.to_string()
}
