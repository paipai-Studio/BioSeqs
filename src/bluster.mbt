///|
/// Portable clustering and diagnostics inspired by Bioconductor bluster 1.23.0.
/// All input matrices use observation x variable orientation.

///|
pub suberror BlusterError {
  BlusterError(String)
}

///|
pub enum BlusterDistance {
  BlusterEuclidean
  BlusterManhattan
  BlusterCosine
} derive(Eq, Debug)

///|
pub enum BlusterSnnWeight {
  BlusterRank
  BlusterNumber
  BlusterJaccard
} derive(Eq, Debug)

///|
pub enum BlusterLinkDenominator {
  BlusterLinkMinimum
  BlusterLinkMaximum
  BlusterLinkUnion
} derive(Eq, Debug)

///|
pub fn bluster_euclidean_distance() -> BlusterDistance {
  BlusterEuclidean
}

///|
pub fn bluster_manhattan_distance() -> BlusterDistance {
  BlusterManhattan
}

///|
pub fn bluster_cosine_distance() -> BlusterDistance {
  BlusterCosine
}

///|
pub fn bluster_rank_weight() -> BlusterSnnWeight {
  BlusterRank
}

///|
pub fn bluster_number_weight() -> BlusterSnnWeight {
  BlusterNumber
}

///|
pub fn bluster_jaccard_weight() -> BlusterSnnWeight {
  BlusterJaccard
}

///|
pub fn bluster_link_minimum() -> BlusterLinkDenominator {
  BlusterLinkMinimum
}

///|
pub fn bluster_link_maximum() -> BlusterLinkDenominator {
  BlusterLinkMaximum
}

///|
pub fn bluster_link_union() -> BlusterLinkDenominator {
  BlusterLinkUnion
}

///|
pub struct BlusterGraph {
  weights : Array[Array[Double]]
  directed : Bool
  k : Int
  graph_kind : String
} derive(Debug)

///|
pub fn BlusterGraph::node_count(self : BlusterGraph) -> Int {
  self.weights.length()
}

///|
pub fn BlusterGraph::edge_weight(
  self : BlusterGraph,
  left : Int,
  right : Int,
) -> Double {
  self.weights[left][right]
}

///|
pub fn BlusterGraph::edge_count(self : BlusterGraph) -> Int {
  let mut output = 0
  for left in 0..<self.weights.length() {
    let start = if self.directed { 0 } else { left + 1 }
    for right in start..<self.weights.length() {
      if left != right && self.weights[left][right] > 0.0 {
        output = output + 1
      }
    }
  }
  output
}

///|
pub struct BlusterKmeansConfig {
  centers : Int
  max_iterations : Int
  starts : Int
  tolerance : Double
  seed : Int
} derive(Eq, Debug)

///|
pub fn BlusterKmeansConfig::create(
  centers : Int,
  max_iterations? : Int = 100,
  starts? : Int = 10,
  tolerance? : Double = 1.0e-8,
  seed? : Int = 1,
) -> BlusterKmeansConfig raise BlusterError {
  if centers <= 0 {
    bluster_fail("bluster k-means centers must be positive")
  }
  if max_iterations <= 0 {
    bluster_fail("bluster k-means max_iterations must be positive")
  }
  if starts <= 0 {
    bluster_fail("bluster k-means starts must be positive")
  }
  if !bluster_is_finite(tolerance) || tolerance < 0.0 {
    bluster_fail("bluster k-means tolerance must be finite and non-negative")
  }
  if seed <= 0 {
    bluster_fail("bluster random seed must be positive")
  }
  BlusterKmeansConfig::{ centers, max_iterations, starts, tolerance, seed }
}

///|
pub struct BlusterKmeansResult {
  clusters : Array[Int]
  centers : Array[Array[Double]]
  within_sum_squares : Array[Double]
  total_within_sum_squares : Double
  iterations : Int
} derive(Debug)

///|
pub struct BlusterGraphConfig {
  k : Int
  shared : Bool
  snn_weight : BlusterSnnWeight
  distance : BlusterDistance
  resolution : Double
  max_iterations : Int
  seed : Int
} derive(Eq, Debug)

///|
pub fn BlusterGraphConfig::create(
  k? : Int = 10,
  shared? : Bool = true,
  snn_weight? : BlusterSnnWeight = BlusterRank,
  distance? : BlusterDistance = BlusterEuclidean,
  resolution? : Double = 1.0,
  max_iterations? : Int = 100,
  seed? : Int = 1,
) -> BlusterGraphConfig raise BlusterError {
  if k <= 0 {
    bluster_fail("bluster graph k must be positive")
  }
  if !bluster_is_finite(resolution) || resolution <= 0.0 {
    bluster_fail("bluster graph resolution must be finite and positive")
  }
  if max_iterations <= 0 {
    bluster_fail("bluster graph max_iterations must be positive")
  }
  if seed <= 0 {
    bluster_fail("bluster random seed must be positive")
  }
  BlusterGraphConfig::{
    k,
    shared,
    snn_weight,
    distance,
    resolution,
    max_iterations,
    seed,
  }
}

///|
pub struct BlusterGraphClusterResult {
  clusters : Array[Int]
  graph : BlusterGraph
  modularity : Double
  iterations : Int
} derive(Debug)

///|
pub struct BlusterTwoStepConfig {
  first_centers : Int
  first_starts : Int
  first_max_iterations : Int
  second_k : Int
  snn_weight : BlusterSnnWeight
  distance : BlusterDistance
  resolution : Double
  second_max_iterations : Int
  seed : Int
} derive(Eq, Debug)

///|
pub fn BlusterTwoStepConfig::create(
  first_centers? : Int = 0,
  first_starts? : Int = 5,
  first_max_iterations? : Int = 100,
  second_k? : Int = 10,
  snn_weight? : BlusterSnnWeight = BlusterRank,
  distance? : BlusterDistance = BlusterEuclidean,
  resolution? : Double = 1.0,
  second_max_iterations? : Int = 100,
  seed? : Int = 1,
) -> BlusterTwoStepConfig raise BlusterError {
  if first_centers < 0 {
    bluster_fail("bluster two-step first_centers cannot be negative")
  }
  if first_starts <= 0 ||
    first_max_iterations <= 0 ||
    second_k <= 0 ||
    second_max_iterations <= 0 {
    bluster_fail(
      "bluster two-step iteration, start and k values must be positive",
    )
  }
  if !bluster_is_finite(resolution) || resolution <= 0.0 {
    bluster_fail("bluster two-step resolution must be finite and positive")
  }
  if seed <= 0 {
    bluster_fail("bluster random seed must be positive")
  }
  BlusterTwoStepConfig::{
    first_centers,
    first_starts,
    first_max_iterations,
    second_k,
    snn_weight,
    distance,
    resolution,
    second_max_iterations,
    seed,
  }
}

///|
pub struct BlusterTwoStepResult {
  clusters : Array[Int]
  first : BlusterKmeansResult
  centroids : Array[Array[Double]]
  centroid_clusters : Array[Int]
  graph : BlusterGraph
  modularity : Double
} derive(Debug)

///|
pub struct BlusterRandResult {
  cluster_names : Array[Int]
  correct : Array[Array[Double?]]
  totals : Array[Array[Double?]]
  ratios : Array[Array[Double?]]
  index : Double?
  adjusted : Bool
} derive(Debug)

///|
pub struct BlusterSilhouetteResult {
  clusters : Array[Int]
  other_clusters : Array[Int]
  widths : Array[Double]
} derive(Debug)

///|
pub struct BlusterClusterMetric {
  cluster_names : Array[Int]
  values : Array[Double?]
} derive(Debug)

///|
pub struct BlusterPurityResult {
  purity : Array[Double]
  maximum_clusters : Array[Int]
  radius : Double
} derive(Debug)

///|
pub struct BlusterModularityResult {
  cluster_names : Array[Int]
  observed : Array[Array[Double?]]
  expected : Array[Array[Double?]]
  modularity : Array[Array[Double?]]
  ratios : Array[Array[Double?]]
  total_weight : Double
} derive(Debug)

///|
pub struct BlusterNestedResult {
  reference_names : Array[Int]
  alternative_names : Array[Int]
  proportions : Array[Array[Double]]
  alternative_mapping : Array[Int]
  alternative_maximum : Array[Double]
  reference_scores : Array[Double]
} derive(Debug)

///|
pub struct BlusterStabilityResult {
  cluster_names : Array[Int]
  ratios : Array[Array[Double?]]
  adjusted_indices : Array[Double?]
} derive(Debug)

///|
pub struct BlusterSceOutput {
  experiment : SingleCellExperiment
  clustering : BlusterGraphClusterResult
}

///|
fn bluster_fail(message : String) -> Unit raise BlusterError {
  raise BlusterError::BlusterError(message)
}

///|
fn bluster_is_finite(value : Double) -> Bool {
  value == value && value.abs() <= 1.0e300
}

///|
fn bluster_matrix_shape(
  data : Array[Array[Double]],
  allow_empty? : Bool = false,
) -> (Int, Int) raise BlusterError {
  let observations = data.length()
  if observations == 0 {
    if allow_empty {
      return (0, 0)
    }
    bluster_fail("bluster data must contain at least one observation")
  }
  let variables = data[0].length()
  if variables == 0 && !allow_empty {
    bluster_fail("bluster data must contain at least one variable")
  }
  for row in data {
    if row.length() != variables {
      bluster_fail("bluster data must be rectangular")
    }
    for value in row {
      if !bluster_is_finite(value) {
        bluster_fail("bluster data must contain only finite values")
      }
    }
  }
  (observations, variables)
}

///|
fn bluster_validate_labels(
  labels : Array[Int],
  expected : Int,
) -> Unit raise BlusterError {
  if labels.length() != expected {
    bluster_fail("bluster cluster labels must match the number of observations")
  }
}

///|
fn bluster_copy_matrix(data : Array[Array[Double]]) -> Array[Array[Double]] {
  let output : Array[Array[Double]] = []
  for row in data {
    output.push(row.copy())
  }
  output
}

///|
fn bluster_double_matrix(
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
fn bluster_int_matrix(
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
fn bluster_optional_double_matrix(
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
fn bluster_copy_sce(source : SingleCellExperiment) -> SingleCellExperiment {
  let assays : Map[String, Array[Array[Double]]] = Map([])
  for key in source.assays.keys() {
    assays[key] = bluster_copy_matrix(source.assays[key])
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
    reduced_dims[key] = bluster_copy_matrix(source.reduced_dims[key])
  }
  let alternatives : Map[String, SingleCellExperiment] = Map([])
  for key in source.alternative_experiments.keys() {
    alternatives[key] = bluster_copy_sce(source.alternative_experiments[key])
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
fn bluster_lcg(state : Int) -> (Int, Double) {
  let normalized = (state - 1) % 2147483646 + 1
  let quotient = normalized / 44488
  let remainder = normalized % 44488
  let candidate = 48271 * remainder - 3399 * quotient
  let next = if candidate > 0 { candidate } else { candidate + 2147483647 }
  (next, next.to_double() / 2147483647.0)
}

///|
fn bluster_distance_value(
  left : Array[Double],
  right : Array[Double],
  distance_kind : BlusterDistance,
) -> Double {
  match distance_kind {
    BlusterEuclidean => {
      let mut total = 0.0
      for variable in 0..<left.length() {
        let difference = left[variable] - right[variable]
        total = total + difference * difference
      }
      total.sqrt()
    }
    BlusterManhattan => {
      let mut total = 0.0
      for variable in 0..<left.length() {
        total = total + (left[variable] - right[variable]).abs()
      }
      total
    }
    BlusterCosine => {
      let mut product = 0.0
      let mut left_norm = 0.0
      let mut right_norm = 0.0
      for variable in 0..<left.length() {
        product = product + left[variable] * right[variable]
        left_norm = left_norm + left[variable] * left[variable]
        right_norm = right_norm + right[variable] * right[variable]
      }
      if left_norm <= 0.0 || right_norm <= 0.0 {
        if left_norm <= 0.0 && right_norm <= 0.0 {
          0.0
        } else {
          1.0
        }
      } else {
        let raw = 1.0 - product / (left_norm.sqrt() * right_norm.sqrt())
        if raw < 0.0 && raw > -1.0e-12 {
          0.0
        } else {
          raw
        }
      }
    }
  }
}

///|
fn bluster_squared_euclidean(
  left : Array[Double],
  right : Array[Double],
) -> Double {
  let mut output = 0.0
  for variable in 0..<left.length() {
    let difference = left[variable] - right[variable]
    output = output + difference * difference
  }
  output
}

///|
fn bluster_sorted_unique(labels : Array[Int]) -> Array[Int] {
  let output : Array[Int] = []
  for label in labels {
    let mut present = false
    for existing in output {
      if existing == label {
        present = true
        break
      }
    }
    if !present {
      let mut position = output.length()
      while position > 0 && output[position - 1] > label {
        position = position - 1
      }
      output.push(label)
      let mut index = output.length() - 1
      while index > position {
        output[index] = output[index - 1]
        index = index - 1
      }
      output[position] = label
    }
  }
  output
}

///|
fn bluster_label_position(labels : Array[Int], target : Int) -> Int {
  for index in 0..<labels.length() {
    if labels[index] == target {
      return index
    }
  }
  -1
}

///|
fn bluster_relabel(clusters : Array[Int]) -> Array[Int] {
  let names = bluster_sorted_unique(clusters)
  clusters.map(fn(label) { bluster_label_position(names, label) })
}

///|
fn bluster_neighbor_lists(
  data : Array[Array[Double]],
  k : Int,
  distance_kind : BlusterDistance,
) -> (Array[Array[Int]], Array[Array[Double]]) {
  let observations = data.length()
  let effective_k = k.min((observations - 1).max(0))
  let indices : Array[Array[Int]] = []
  let distances : Array[Array[Double]] = []
  for observation in 0..<observations {
    let candidates : Array[(Double, Int)] = []
    for other in 0..<observations {
      if other != observation {
        candidates.push(
          (
            bluster_distance_value(
              data[observation],
              data[other],
              distance_kind,
            ),
            other,
          ),
        )
      }
    }
    for current in 1..<candidates.length() {
      let value = candidates[current]
      let mut position = current
      while position > 0 &&
            (
              candidates[position - 1].0 > value.0 ||
              (
                candidates[position - 1].0 == value.0 &&
                candidates[position - 1].1 > value.1
              )
            ) {
        candidates[position] = candidates[position - 1]
        position = position - 1
      }
      candidates[position] = value
    }
    let local_indices : Array[Int] = []
    let local_distances : Array[Double] = []
    for position in 0..<effective_k {
      local_distances.push(candidates[position].0)
      local_indices.push(candidates[position].1)
    }
    indices.push(local_indices)
    distances.push(local_distances)
  }
  (indices, distances)
}

///|
fn bluster_neighbor_rank(
  observation : Int,
  target : Int,
  neighbors : Array[Int],
) -> Int {
  if observation == target {
    return 0
  }
  for position in 0..<neighbors.length() {
    if neighbors[position] == target {
      return position + 1
    }
  }
  -1
}

///|
pub fn bluster_make_snn_graph(
  data : Array[Array[Double]],
  k? : Int = 10,
  weighting? : BlusterSnnWeight = BlusterRank,
  distance? : BlusterDistance = BlusterEuclidean,
) -> BlusterGraph raise BlusterError {
  let (observations, _) = bluster_matrix_shape(data)
  if k <= 0 {
    bluster_fail("bluster SNN k must be positive")
  }
  let effective_k = k.min((observations - 1).max(0))
  let (neighbors, _) = bluster_neighbor_lists(data, effective_k, distance)
  let weights = bluster_double_matrix(observations, observations, 0.0)
  for left in 0..<observations {
    for right in (left + 1)..<observations {
      let mut shared = 0
      let mut rank_weight = 0.0
      for candidate in 0..<observations {
        let left_rank = bluster_neighbor_rank(left, candidate, neighbors[left])
        if left_rank >= 0 {
          let right_rank = bluster_neighbor_rank(
            right,
            candidate,
            neighbors[right],
          )
          if right_rank >= 0 {
            shared = shared + 1
            let score = effective_k.to_double() -
              (left_rank + right_rank).to_double() / 2.0
            if score > rank_weight {
              rank_weight = score
            }
          }
        }
      }
      if shared > 0 {
        let weight = match weighting {
          BlusterRank => rank_weight.max(1.0e-6)
          BlusterNumber => shared.to_double()
          BlusterJaccard => {
            let union = 2 * (effective_k + 1) - shared
            shared.to_double() / union.to_double()
          }
        }
        weights[left][right] = weight
        weights[right][left] = weight
      }
    }
  }
  BlusterGraph::{ weights, directed: false, k: effective_k, graph_kind: "SNN" }
}

///|
pub fn bluster_make_knn_graph(
  data : Array[Array[Double]],
  k? : Int = 10,
  directed? : Bool = false,
  distance? : BlusterDistance = BlusterEuclidean,
) -> BlusterGraph raise BlusterError {
  let (observations, _) = bluster_matrix_shape(data)
  if k <= 0 {
    bluster_fail("bluster KNN k must be positive")
  }
  let effective_k = k.min((observations - 1).max(0))
  let (neighbors, _) = bluster_neighbor_lists(data, effective_k, distance)
  let weights = bluster_double_matrix(observations, observations, 0.0)
  for observation in 0..<observations {
    for neighbor in neighbors[observation] {
      weights[observation][neighbor] = 1.0
      if !directed {
        weights[neighbor][observation] = 1.0
      }
    }
  }
  BlusterGraph::{ weights, directed, k: effective_k, graph_kind: "KNN" }
}

///|
fn bluster_initialize_centers(
  data : Array[Array[Double]],
  centers : Int,
  seed : Int,
) -> (Array[Array[Double]], Int) {
  let observations = data.length()
  let selected = Array::make(observations, false)
  let output : Array[Array[Double]] = []
  let mut state = seed
  let (next, first_uniform) = bluster_lcg(state)
  state = next
  let first = (first_uniform * observations.to_double())
    .to_int()
    .min(observations - 1)
  output.push(data[first].copy())
  selected[first] = true
  while output.length() < centers {
    let nearest = Array::make(observations, 0.0)
    let mut total = 0.0
    for observation in 0..<observations {
      let mut best = 1.0e300
      for center in output {
        let distance = bluster_squared_euclidean(data[observation], center)
        if distance < best {
          best = distance
        }
      }
      nearest[observation] = best
      total = total + best
    }
    let mut chosen = -1
    if total > 0.0 {
      let (following, uniform) = bluster_lcg(state)
      state = following
      let target = uniform * total
      let mut cumulative = 0.0
      for observation in 0..<observations {
        cumulative = cumulative + nearest[observation]
        if chosen < 0 && cumulative >= target && !selected[observation] {
          chosen = observation
        }
      }
    }
    if chosen < 0 {
      for observation in 0..<observations {
        if !selected[observation] {
          chosen = observation
          break
        }
      }
    }
    if chosen < 0 {
      chosen = output.length() % observations
    }
    output.push(data[chosen].copy())
    selected[chosen] = true
  }
  (output, state)
}

///|
fn bluster_assign_centers(
  data : Array[Array[Double]],
  centers : Array[Array[Double]],
) -> Array[Int] {
  let output = Array::make(data.length(), 0)
  for observation in 0..<data.length() {
    let mut best_center = 0
    let mut best_distance = bluster_squared_euclidean(
      data[observation],
      centers[0],
    )
    for center in 1..<centers.length() {
      let distance = bluster_squared_euclidean(
        data[observation],
        centers[center],
      )
      if distance < best_distance {
        best_distance = distance
        best_center = center
      }
    }
    output[observation] = best_center
  }
  output
}

///|
fn bluster_update_centers(
  data : Array[Array[Double]],
  clusters : Array[Int],
  previous : Array[Array[Double]],
) -> Array[Array[Double]] {
  let centers = previous.length()
  let variables = data[0].length()
  let output = bluster_double_matrix(centers, variables, 0.0)
  let sizes = Array::make(centers, 0)
  for observation in 0..<data.length() {
    let cluster = clusters[observation]
    sizes[cluster] = sizes[cluster] + 1
    for variable in 0..<variables {
      output[cluster][variable] = output[cluster][variable] +
        data[observation][variable]
    }
  }
  for center in 0..<centers {
    if sizes[center] > 0 {
      for variable in 0..<variables {
        output[center][variable] = output[center][variable] /
          sizes[center].to_double()
      }
    } else {
      output[center] = previous[center].copy()
    }
  }
  output
}

///|
fn bluster_kmeans_statistics(
  data : Array[Array[Double]],
  clusters : Array[Int],
  centers : Array[Array[Double]],
) -> (Array[Double], Double) {
  let within = Array::make(centers.length(), 0.0)
  let mut total = 0.0
  for observation in 0..<data.length() {
    let cluster = clusters[observation]
    let value = bluster_squared_euclidean(data[observation], centers[cluster])
    within[cluster] = within[cluster] + value
    total = total + value
  }
  (within, total)
}

///|
fn bluster_kmeans_once(
  data : Array[Array[Double]],
  config : BlusterKmeansConfig,
  seed : Int,
) -> BlusterKmeansResult {
  let (initial, _) = bluster_initialize_centers(data, config.centers, seed)
  let mut centers = initial
  let mut clusters = Array::make(data.length(), 0)
  let mut iterations = 0
  for iteration in 0..<config.max_iterations {
    clusters = bluster_assign_centers(data, centers)
    let updated = bluster_update_centers(data, clusters, centers)
    let mut maximum_shift = 0.0
    for center in 0..<centers.length() {
      let shift = bluster_squared_euclidean(centers[center], updated[center]).sqrt()
      if shift > maximum_shift {
        maximum_shift = shift
      }
    }
    centers = updated
    iterations = iteration + 1
    if maximum_shift <= config.tolerance {
      break
    }
  }
  clusters = bluster_assign_centers(data, centers)
  let (within, total) = bluster_kmeans_statistics(data, clusters, centers)
  BlusterKmeansResult::{
    clusters,
    centers,
    within_sum_squares: within,
    total_within_sum_squares: total,
    iterations,
  }
}

///|
pub fn bluster_cluster_rows_kmeans(
  data : Array[Array[Double]],
  config : BlusterKmeansConfig,
) -> BlusterKmeansResult raise BlusterError {
  let (observations, _) = bluster_matrix_shape(data)
  if config.centers > observations {
    bluster_fail("bluster k-means centers cannot exceed observations")
  }
  let mut best : BlusterKmeansResult? = None
  for start in 0..<config.starts {
    let raw_seed = (config.seed + start * 104729) % 2147483646 + 1
    let candidate = bluster_kmeans_once(data, config, raw_seed)
    match best {
      None => best = Some(candidate)
      Some(current) =>
        if candidate.total_within_sum_squares < current.total_within_sum_squares {
          best = Some(candidate)
        }
    }
  }
  match best {
    Some(output) => output
    None => {
      bluster_fail("bluster k-means did not produce a result")
      bluster_kmeans_once(data, config, config.seed)
    }
  }
}

///|
fn bluster_undirected_weight(
  graph : BlusterGraph,
  left : Int,
  right : Int,
) -> Double {
  if graph.directed {
    graph.weights[left][right] + graph.weights[right][left]
  } else {
    graph.weights[left][right]
  }
}

///|
fn bluster_graph_total_weight(graph : BlusterGraph) -> Double {
  let mut total = 0.0
  for left in 0..<graph.weights.length() {
    total = total + graph.weights[left][left]
    for right in (left + 1)..<graph.weights.length() {
      total = total + bluster_undirected_weight(graph, left, right)
    }
  }
  total
}

///|
fn bluster_global_modularity(
  graph : BlusterGraph,
  clusters : Array[Int],
  resolution : Double,
) -> Double {
  let observations = graph.weights.length()
  let total = bluster_graph_total_weight(graph)
  if total <= 0.0 {
    return 0.0
  }
  let strengths = Array::make(observations, 0.0)
  for left in 0..<observations {
    strengths[left] = strengths[left] + 2.0 * graph.weights[left][left]
    for right in 0..<observations {
      if left != right {
        strengths[left] = strengths[left] +
          bluster_undirected_weight(graph, left, right)
      }
    }
  }
  let double_total = 2.0 * total
  let mut output = 0.0
  for left in 0..<observations {
    for right in 0..<observations {
      if clusters[left] == clusters[right] {
        let observed = if left == right {
          2.0 * graph.weights[left][left]
        } else {
          bluster_undirected_weight(graph, left, right)
        }
        output = output +
          observed -
          resolution * strengths[left] * strengths[right] / double_total
      }
    }
  }
  output / double_total
}

///|
fn bluster_shuffled_indices(length : Int, seed : Int) -> (Array[Int], Int) {
  let output = Array::make(length, 0)
  for index in 0..<length {
    output[index] = index
  }
  let mut state = seed
  let mut position = length - 1
  while position > 0 {
    let (next, uniform) = bluster_lcg(state)
    state = next
    let selected = (uniform * (position + 1).to_double()).to_int().min(position)
    let temporary = output[position]
    output[position] = output[selected]
    output[selected] = temporary
    position = position - 1
  }
  (output, state)
}

///|
fn bluster_louvain(
  graph : BlusterGraph,
  resolution : Double,
  max_iterations : Int,
  seed : Int,
) -> (Array[Int], Double, Int) {
  let observations = graph.weights.length()
  let clusters = Array::make(observations, 0)
  for observation in 0..<observations {
    clusters[observation] = observation
  }
  let mut current_modularity = bluster_global_modularity(
    graph, clusters, resolution,
  )
  let mut state = seed
  let mut completed = 0
  for iteration in 0..<max_iterations {
    let (order, next_state) = bluster_shuffled_indices(observations, state)
    state = next_state
    let mut changed = false
    for observation in order {
      let original = clusters[observation]
      let candidates : Array[Int] = [original]
      for other in 0..<observations {
        if other != observation &&
          bluster_undirected_weight(graph, observation, other) > 0.0 {
          let candidate = clusters[other]
          let mut present = false
          for existing in candidates {
            if existing == candidate {
              present = true
              break
            }
          }
          if !present {
            candidates.push(candidate)
          }
        }
      }
      let mut best_cluster = original
      let mut best_modularity = current_modularity
      for candidate in candidates {
        if candidate != original {
          clusters[observation] = candidate
          let value = bluster_global_modularity(graph, clusters, resolution)
          if value > best_modularity + 1.0e-12 ||
            (
              (value - best_modularity).abs() <= 1.0e-12 &&
              candidate < best_cluster
            ) {
            best_modularity = value
            best_cluster = candidate
          }
        }
      }
      clusters[observation] = best_cluster
      if best_cluster != original {
        changed = true
        current_modularity = best_modularity
      }
    }
    completed = iteration + 1
    if !changed {
      break
    }
  }
  let relabeled = bluster_relabel(clusters)
  (
    relabeled,
    bluster_global_modularity(graph, relabeled, resolution),
    completed,
  )
}

///|
pub fn bluster_cluster_rows_graph(
  data : Array[Array[Double]],
  config : BlusterGraphConfig,
) -> BlusterGraphClusterResult raise BlusterError {
  ignore(bluster_matrix_shape(data))
  let graph = if config.shared {
    bluster_make_snn_graph(
      data,
      k=config.k,
      weighting=config.snn_weight,
      distance=config.distance,
    )
  } else {
    bluster_make_knn_graph(
      data,
      k=config.k,
      directed=false,
      distance=config.distance,
    )
  }
  let (clusters, modularity, iterations) = bluster_louvain(
    graph,
    config.resolution,
    config.max_iterations,
    config.seed,
  )
  BlusterGraphClusterResult::{ clusters, graph, modularity, iterations }
}

///|
pub fn bluster_cluster_rows_two_step(
  data : Array[Array[Double]],
  config : BlusterTwoStepConfig,
) -> BlusterTwoStepResult raise BlusterError {
  let (observations, _) = bluster_matrix_shape(data)
  let automatic = observations.to_double().sqrt().round().to_int().max(1)
  let first_centers = (if config.first_centers == 0 {
    automatic
  } else {
    config.first_centers
  }).min(observations)
  let first_config = BlusterKmeansConfig::create(
    first_centers,
    max_iterations=config.first_max_iterations,
    starts=config.first_starts,
    seed=config.seed,
  )
  let first = bluster_cluster_rows_kmeans(data, first_config)
  let centroids = bluster_copy_matrix(first.centers)
  let second_config = BlusterGraphConfig::create(
    k=config.second_k,
    shared=true,
    snn_weight=config.snn_weight,
    distance=config.distance,
    resolution=config.resolution,
    max_iterations=config.second_max_iterations,
    seed=config.seed,
  )
  let second = bluster_cluster_rows_graph(first.centers, second_config)
  let clusters = first.clusters.map(fn(cluster) { second.clusters[cluster] })
  BlusterTwoStepResult::{
    clusters: bluster_relabel(clusters),
    first,
    centroids,
    centroid_clusters: second.clusters,
    graph: second.graph,
    modularity: second.modularity,
  }
}

///|
fn bluster_choose_two(value : Int) -> Double {
  if value < 2 {
    0.0
  } else {
    value.to_double() * (value - 1).to_double() / 2.0
  }
}

///|
fn bluster_pairwise_rand_levels(
  reference : Array[Int],
  alternative : Array[Int],
  levels : Array[Int],
  adjusted : Bool,
) -> BlusterRandResult raise BlusterError {
  if reference.length() != alternative.length() {
    bluster_fail("bluster Rand inputs must have equal lengths")
  }
  let cluster_count = levels.length()
  let correct = bluster_optional_double_matrix(cluster_count, cluster_count)
  let totals = bluster_optional_double_matrix(cluster_count, cluster_count)
  let ratios = bluster_optional_double_matrix(cluster_count, cluster_count)
  let alternative_names = bluster_sorted_unique(alternative)
  let reference_sizes = Array::make(cluster_count, 0)
  let alternative_sizes = Array::make(alternative_names.length(), 0)
  let table = bluster_int_matrix(cluster_count, alternative_names.length(), 0)
  for observation in 0..<reference.length() {
    let reference_position = bluster_label_position(
      levels,
      reference[observation],
    )
    let alternative_position = bluster_label_position(
      alternative_names,
      alternative[observation],
    )
    if reference_position >= 0 {
      reference_sizes[reference_position] = reference_sizes[reference_position] +
        1
      if alternative_position >= 0 {
        table[reference_position][alternative_position] = table[reference_position][alternative_position] +
          1
      }
    }
    if alternative_position >= 0 {
      alternative_sizes[alternative_position] = alternative_sizes[alternative_position] +
        1
    }
  }
  let all_pairs = bluster_choose_two(reference.length())
  let mut same_alternative = 0.0
  for size in alternative_sizes {
    same_alternative = same_alternative + bluster_choose_two(size)
  }
  let same_probability = if all_pairs > 0.0 {
    same_alternative / all_pairs
  } else {
    0.0
  }
  let mut index_correct = 0.0
  let mut index_total = 0.0
  for right in 0..<cluster_count {
    for left in 0..<(right + 1) {
      let raw_total = if left == right {
        bluster_choose_two(reference_sizes[left])
      } else {
        reference_sizes[left].to_double() * reference_sizes[right].to_double()
      }
      let mut raw_correct = 0.0
      if left == right {
        for alternative_position in 0..<alternative_names.length() {
          raw_correct = raw_correct +
            bluster_choose_two(table[left][alternative_position])
        }
      } else {
        let mut same = 0.0
        for alternative_position in 0..<alternative_names.length() {
          same = same +
            table[left][alternative_position].to_double() *
            table[right][alternative_position].to_double()
        }
        raw_correct = raw_total - same
      }
      let expected = if !adjusted {
        0.0
      } else if left == right {
        raw_total * same_probability
      } else {
        raw_total * (1.0 - same_probability)
      }
      let adjusted_correct = raw_correct - expected
      let adjusted_total = raw_total - expected
      correct[left][right] = Some(adjusted_correct)
      totals[left][right] = Some(adjusted_total)
      if adjusted_total.abs() > 1.0e-15 {
        ratios[left][right] = Some(adjusted_correct / adjusted_total)
        index_correct = index_correct + adjusted_correct
        index_total = index_total + adjusted_total
      }
    }
  }
  let index = if index_total.abs() > 1.0e-15 {
    Some(index_correct / index_total)
  } else {
    None
  }
  BlusterRandResult::{
    cluster_names: levels.copy(),
    correct,
    totals,
    ratios,
    index,
    adjusted,
  }
}

///|
pub fn bluster_pairwise_rand(
  reference : Array[Int],
  alternative : Array[Int],
  adjusted? : Bool = true,
) -> BlusterRandResult raise BlusterError {
  bluster_pairwise_rand_levels(
    reference,
    alternative,
    bluster_sorted_unique(reference),
    adjusted,
  )
}

///|
pub fn bluster_compare_clusterings(
  clusterings : Array[Array[Int]],
  adjusted? : Bool = true,
) -> Array[Array[Double?]] raise BlusterError {
  if clusterings.length() == 0 {
    return []
  }
  let observations = clusterings[0].length()
  for clustering in clusterings {
    if clustering.length() != observations {
      bluster_fail("bluster clusterings must have equal lengths")
    }
  }
  let output = bluster_optional_double_matrix(
    clusterings.length(),
    clusterings.length(),
  )
  for left in 0..<clusterings.length() {
    output[left][left] = Some(1.0)
    for right in (left + 1)..<clusterings.length() {
      let result = bluster_pairwise_rand(
        clusterings[left],
        clusterings[right],
        adjusted~,
      )
      output[left][right] = result.index
      output[right][left] = result.index
    }
  }
  output
}

///|
pub fn bluster_approx_silhouette(
  data : Array[Array[Double]],
  clusters : Array[Int],
) -> BlusterSilhouetteResult raise BlusterError {
  let (observations, variables) = bluster_matrix_shape(data)
  bluster_validate_labels(clusters, observations)
  let names = bluster_sorted_unique(clusters)
  let centroids = bluster_double_matrix(names.length(), variables, 0.0)
  let sizes = Array::make(names.length(), 0)
  for observation in 0..<observations {
    let position = bluster_label_position(names, clusters[observation])
    sizes[position] = sizes[position] + 1
    for variable in 0..<variables {
      centroids[position][variable] = centroids[position][variable] +
        data[observation][variable]
    }
  }
  for cluster in 0..<names.length() {
    for variable in 0..<variables {
      centroids[cluster][variable] = centroids[cluster][variable] /
        sizes[cluster].to_double()
    }
  }
  let variances = Array::make(names.length(), 0.0)
  for observation in 0..<observations {
    let position = bluster_label_position(names, clusters[observation])
    variances[position] = variances[position] +
      bluster_squared_euclidean(data[observation], centroids[position])
  }
  for cluster in 0..<names.length() {
    variances[cluster] = variances[cluster] / sizes[cluster].to_double()
  }
  let other_clusters = Array::make(observations, 0)
  let widths = Array::make(observations, 0.0)
  for observation in 0..<observations {
    let self_position = bluster_label_position(names, clusters[observation])
    let self_distance = (bluster_squared_euclidean(
      data[observation],
      centroids[self_position],
    ) +
    variances[self_position]).sqrt()
    let mut other_position = self_position
    let mut other_distance = 1.0e300
    for candidate in 0..<names.length() {
      if candidate != self_position {
        let distance = (bluster_squared_euclidean(
          data[observation],
          centroids[candidate],
        ) +
        variances[candidate]).sqrt()
        if distance < other_distance {
          other_distance = distance
          other_position = candidate
        }
      }
    }
    if names.length() == 1 {
      other_clusters[observation] = clusters[observation]
      widths[observation] = 0.0
    } else {
      other_clusters[observation] = names[other_position]
      let denominator = self_distance.max(other_distance)
      widths[observation] = if denominator > 0.0 {
        (other_distance - self_distance) / denominator
      } else {
        0.0
      }
    }
  }
  BlusterSilhouetteResult::{ clusters: clusters.copy(), other_clusters, widths }
}

///|
pub fn bluster_cluster_rmsd(
  data : Array[Array[Double]],
  clusters : Array[Int],
  sum_squares? : Bool = false,
) -> BlusterClusterMetric raise BlusterError {
  let (observations, variables) = bluster_matrix_shape(data, allow_empty=true)
  bluster_validate_labels(clusters, observations)
  let names = bluster_sorted_unique(clusters)
  let values : Array[Double?] = Array::make(names.length(), None)
  for cluster in 0..<names.length() {
    let indices : Array[Int] = []
    for observation in 0..<observations {
      if clusters[observation] == names[cluster] {
        indices.push(observation)
      }
    }
    if indices.length() >= 2 {
      let means = Array::make(variables, 0.0)
      for observation in indices {
        for variable in 0..<variables {
          means[variable] = means[variable] + data[observation][variable]
        }
      }
      for variable in 0..<variables {
        means[variable] = means[variable] / indices.length().to_double()
      }
      let mut total = 0.0
      for observation in indices {
        total = total + bluster_squared_euclidean(data[observation], means)
      }
      values[cluster] = if sum_squares {
        Some(total)
      } else {
        Some((total / (indices.length() - 1).to_double()).sqrt())
      }
    }
  }
  BlusterClusterMetric::{ cluster_names: names, values }
}

///|
fn bluster_median(values : Array[Double]) -> Double {
  if values.length() == 0 {
    return 0.0
  }
  let sorted = values.copy()
  for current in 1..<sorted.length() {
    let value = sorted[current]
    let mut position = current
    while position > 0 && sorted[position - 1] > value {
      sorted[position] = sorted[position - 1]
      position = position - 1
    }
    sorted[position] = value
  }
  if sorted.length() % 2 == 1 {
    sorted[sorted.length() / 2]
  } else {
    (sorted[sorted.length() / 2 - 1] + sorted[sorted.length() / 2]) / 2.0
  }
}

///|
pub fn bluster_neighbor_purity(
  data : Array[Array[Double]],
  clusters : Array[Int],
  k? : Int = 50,
  balanced? : Bool = true,
  custom_weights? : Array[Double]? = None,
  distance? : BlusterDistance = BlusterEuclidean,
) -> BlusterPurityResult raise BlusterError {
  let (observations, _) = bluster_matrix_shape(data)
  bluster_validate_labels(clusters, observations)
  if k <= 0 {
    bluster_fail("bluster neighbor purity k must be positive")
  }
  let effective_k = k.min((observations - 1).max(0))
  let (_, neighbor_distances) = bluster_neighbor_lists(
    data, effective_k, distance,
  )
  let kth : Array[Double] = []
  if effective_k > 0 {
    for neighbor_row in neighbor_distances {
      kth.push(neighbor_row[effective_k - 1])
    }
  }
  let radius = bluster_median(kth)
  let names = bluster_sorted_unique(clusters)
  let sizes = Array::make(names.length(), 0)
  for label in clusters {
    let position = bluster_label_position(names, label)
    sizes[position] = sizes[position] + 1
  }
  let weights = Array::make(observations, 1.0)
  match custom_weights {
    Some(given) => {
      if given.length() != observations {
        bluster_fail("bluster custom purity weights must match observations")
      }
      for observation in 0..<observations {
        if !bluster_is_finite(given[observation]) || given[observation] < 0.0 {
          bluster_fail(
            "bluster custom purity weights must be finite and non-negative",
          )
        }
        weights[observation] = given[observation]
      }
    }
    None =>
      if balanced {
        for observation in 0..<observations {
          let position = bluster_label_position(names, clusters[observation])
          weights[observation] = 1.0 / sizes[position].to_double()
        }
      }
  }
  let purity = Array::make(observations, 0.0)
  let maximum_clusters = Array::make(observations, names[0])
  for observation in 0..<observations {
    let totals = Array::make(names.length(), 0.0)
    for other in 0..<observations {
      let local_distance = bluster_distance_value(
        data[observation],
        data[other],
        distance,
      )
      if local_distance <= radius + 1.0e-12 {
        let position = bluster_label_position(names, clusters[other])
        totals[position] = totals[position] + weights[other]
      }
    }
    let own = bluster_label_position(names, clusters[observation])
    let mut total = 0.0
    let mut maximum = 0
    for position in 0..<totals.length() {
      total = total + totals[position]
      if totals[position] > totals[maximum] {
        maximum = position
      }
    }
    purity[observation] = if total > 0.0 { totals[own] / total } else { 0.0 }
    maximum_clusters[observation] = names[maximum]
  }
  BlusterPurityResult::{ purity, maximum_clusters, radius }
}

///|
pub fn bluster_pairwise_modularity(
  graph : BlusterGraph,
  clusters : Array[Int],
) -> BlusterModularityResult raise BlusterError {
  let observations = graph.weights.length()
  bluster_validate_labels(clusters, observations)
  let names = bluster_sorted_unique(clusters)
  let cluster_count = names.length()
  let observed = bluster_optional_double_matrix(cluster_count, cluster_count)
  let expected = bluster_optional_double_matrix(cluster_count, cluster_count)
  let modularity = bluster_optional_double_matrix(cluster_count, cluster_count)
  let ratios = bluster_optional_double_matrix(cluster_count, cluster_count)
  let total_weight = bluster_graph_total_weight(graph)
  let strengths = Array::make(cluster_count, 0.0)
  for observation in 0..<observations {
    let cluster = bluster_label_position(names, clusters[observation])
    strengths[cluster] = strengths[cluster] +
      2.0 * graph.weights[observation][observation]
    for other in 0..<observations {
      if other != observation {
        strengths[cluster] = strengths[cluster] +
          bluster_undirected_weight(graph, observation, other)
      }
    }
  }
  for right in 0..<cluster_count {
    for left in 0..<(right + 1) {
      let mut actual = 0.0
      for first in 0..<observations {
        if clusters[first] == names[left] {
          if left == right {
            actual = actual + graph.weights[first][first]
            for second in (first + 1)..<observations {
              if clusters[second] == names[right] {
                actual = actual +
                  bluster_undirected_weight(graph, first, second)
              }
            }
          } else {
            for second in 0..<observations {
              if clusters[second] == names[right] {
                actual = actual +
                  bluster_undirected_weight(graph, first, second)
              }
            }
          }
        }
      }
      let null_weight = if total_weight > 0.0 {
        if left == right {
          strengths[left] * strengths[right] / (4.0 * total_weight)
        } else {
          strengths[left] * strengths[right] / (2.0 * total_weight)
        }
      } else {
        0.0
      }
      observed[left][right] = Some(actual)
      expected[left][right] = Some(null_weight)
      modularity[left][right] = if total_weight > 0.0 {
        Some((actual - null_weight) / total_weight)
      } else {
        None
      }
      ratios[left][right] = if null_weight > 0.0 {
        Some(actual / null_weight)
      } else {
        None
      }
    }
  }
  BlusterModularityResult::{
    cluster_names: names,
    observed,
    expected,
    modularity,
    ratios,
    total_weight,
  }
}

///|
pub fn bluster_merge_communities(
  graph : BlusterGraph,
  clusters : Array[Int],
  number : Int,
) -> Array[Int] raise BlusterError {
  bluster_validate_labels(clusters, graph.weights.length())
  if number <= 0 {
    bluster_fail("bluster merged community number must be positive")
  }
  let output = clusters.copy()
  while bluster_sorted_unique(output).length() > number {
    let names = bluster_sorted_unique(output)
    let mut best_labels = output.copy()
    let mut best_modularity = -1.0e300
    for left in 0..<names.length() {
      for right in (left + 1)..<names.length() {
        let candidate = output.copy()
        for observation in 0..<candidate.length() {
          if candidate[observation] == names[right] {
            candidate[observation] = names[left]
          }
        }
        let value = bluster_global_modularity(graph, candidate, 1.0)
        if value > best_modularity {
          best_modularity = value
          best_labels = candidate
        }
      }
    }
    for observation in 0..<output.length() {
      output[observation] = best_labels[observation]
    }
  }
  output
}

///|
pub fn bluster_nested_clusters(
  reference : Array[Int],
  alternative : Array[Int],
) -> BlusterNestedResult raise BlusterError {
  if reference.length() != alternative.length() {
    bluster_fail("bluster nested clusterings must have equal lengths")
  }
  let reference_names = bluster_sorted_unique(reference)
  let alternative_names = bluster_sorted_unique(alternative)
  let table = bluster_int_matrix(
    alternative_names.length(),
    reference_names.length(),
    0,
  )
  for observation in 0..<reference.length() {
    let alternative_position = bluster_label_position(
      alternative_names,
      alternative[observation],
    )
    let reference_position = bluster_label_position(
      reference_names,
      reference[observation],
    )
    table[alternative_position][reference_position] = table[alternative_position][reference_position] +
      1
  }
  let proportions = bluster_double_matrix(
    alternative_names.length(),
    reference_names.length(),
    0.0,
  )
  let alternative_mapping = Array::make(
    alternative_names.length(),
    if reference_names.length() > 0 {
      reference_names[0]
    } else {
      0
    },
  )
  let alternative_maximum = Array::make(alternative_names.length(), 0.0)
  for alternative_position in 0..<alternative_names.length() {
    let mut total = 0
    for count in table[alternative_position] {
      total = total + count
    }
    let mut maximum = 0
    for reference_position in 0..<reference_names.length() {
      proportions[alternative_position][reference_position] = if total > 0 {
        table[alternative_position][reference_position].to_double() /
        total.to_double()
      } else {
        0.0
      }
      if proportions[alternative_position][reference_position] >
        proportions[alternative_position][maximum] {
        maximum = reference_position
      }
    }
    if reference_names.length() > 0 {
      alternative_mapping[alternative_position] = reference_names[maximum]
      alternative_maximum[alternative_position] = proportions[alternative_position][maximum]
    }
  }
  let reference_scores = Array::make(reference_names.length(), 0.0)
  for reference_position in 0..<reference_names.length() {
    let mut numerator = 0.0
    let mut denominator = 0
    for alternative_position in 0..<alternative_names.length() {
      numerator = numerator +
        proportions[alternative_position][reference_position] *
        table[alternative_position][reference_position].to_double()
      denominator = denominator +
        table[alternative_position][reference_position]
    }
    reference_scores[reference_position] = if denominator > 0 {
      numerator / denominator.to_double()
    } else {
      0.0
    }
  }
  BlusterNestedResult::{
    reference_names,
    alternative_names,
    proportions,
    alternative_mapping,
    alternative_maximum,
    reference_scores,
  }
}

///|
pub fn bluster_link_clusters_matrix(
  left : Array[Int],
  right : Array[Int],
  denominator? : BlusterLinkDenominator = BlusterLinkUnion,
) -> Array[Array[Double]] raise BlusterError {
  if left.length() != right.length() {
    bluster_fail("bluster linked clusterings must have equal lengths")
  }
  let left_names = bluster_sorted_unique(left)
  let right_names = bluster_sorted_unique(right)
  let left_sizes = Array::make(left_names.length(), 0)
  let right_sizes = Array::make(right_names.length(), 0)
  let overlap = bluster_int_matrix(left_names.length(), right_names.length(), 0)
  for observation in 0..<left.length() {
    let left_position = bluster_label_position(left_names, left[observation])
    let right_position = bluster_label_position(right_names, right[observation])
    left_sizes[left_position] = left_sizes[left_position] + 1
    right_sizes[right_position] = right_sizes[right_position] + 1
    overlap[left_position][right_position] = overlap[left_position][right_position] +
      1
  }
  let output = bluster_double_matrix(
    left_names.length(),
    right_names.length(),
    0.0,
  )
  for left_position in 0..<left_names.length() {
    for right_position in 0..<right_names.length() {
      let common = overlap[left_position][right_position]
      let scale = match denominator {
        BlusterLinkMinimum =>
          left_sizes[left_position].min(right_sizes[right_position])
        BlusterLinkMaximum =>
          left_sizes[left_position].max(right_sizes[right_position])
        BlusterLinkUnion =>
          left_sizes[left_position] + right_sizes[right_position] - common
      }
      output[left_position][right_position] = if scale > 0 {
        common.to_double() / scale.to_double()
      } else {
        0.0
      }
    }
  }
  output
}

///|
fn bluster_average_optional(
  values : Array[Double?],
  use_mean : Bool,
) -> Double? {
  let retained : Array[Double] = []
  for value in values {
    match value {
      Some(actual) => if bluster_is_finite(actual) { retained.push(actual) }
      None => ()
    }
  }
  if retained.length() == 0 {
    None
  } else if use_mean {
    let mut total = 0.0
    for value in retained {
      total = total + value
    }
    Some(total / retained.length().to_double())
  } else {
    Some(bluster_median(retained))
  }
}

///|
pub fn bluster_bootstrap_kmeans_stability(
  data : Array[Array[Double]],
  config : BlusterKmeansConfig,
  iterations? : Int = 20,
  use_mean? : Bool = false,
  seed? : Int = 1,
) -> BlusterStabilityResult raise BlusterError {
  let (observations, _) = bluster_matrix_shape(data)
  if iterations <= 0 {
    bluster_fail("bluster bootstrap iterations must be positive")
  }
  if seed <= 0 {
    bluster_fail("bluster random seed must be positive")
  }
  let original = bluster_cluster_rows_kmeans(data, config)
  let names = bluster_sorted_unique(original.clusters)
  let collated : Array[Array[Array[Double?]]] = []
  let indices : Array[Double?] = []
  let mut state = seed
  for iteration in 0..<iterations {
    let chosen = Array::make(observations, 0)
    let resampled : Array[Array[Double]] = []
    for position in 0..<observations {
      let (next, uniform) = bluster_lcg(state)
      state = next
      let selected = (uniform * observations.to_double())
        .to_int()
        .min(observations - 1)
      chosen[position] = selected
      resampled.push(data[selected].copy())
    }
    let iteration_config = BlusterKmeansConfig::{
      ..config,
      seed: (config.seed + (iteration + 1) * 104729) % 2147483646 + 1,
    }
    let reclustered = bluster_cluster_rows_kmeans(resampled, iteration_config)
    let reference = chosen.map(fn(index) { original.clusters[index] })
    let comparison = bluster_pairwise_rand_levels(
      reference,
      reclustered.clusters,
      names,
      true,
    )
    collated.push(comparison.ratios)
    indices.push(comparison.index)
  }
  let ratios = bluster_optional_double_matrix(names.length(), names.length())
  for right in 0..<names.length() {
    for left in 0..<(right + 1) {
      let values : Array[Double?] = []
      for iteration in 0..<iterations {
        values.push(collated[iteration][left][right])
      }
      ratios[left][right] = bluster_average_optional(values, use_mean)
    }
  }
  BlusterStabilityResult::{
    cluster_names: names,
    ratios,
    adjusted_indices: indices,
  }
}

///|
pub fn bluster_cluster_sce(
  experiment : SingleCellExperiment,
  config : BlusterGraphConfig,
  reduced_dim_name? : String = "PCA",
  output_column? : String = "bluster.cluster",
) -> BlusterSceOutput raise BlusterError {
  if reduced_dim_name.length() == 0 || output_column.length() == 0 {
    bluster_fail("bluster SCE names cannot be empty")
  }
  let coordinates = match experiment.reduced_dims.get(reduced_dim_name) {
    Some(value) => value
    None => {
      bluster_fail(
        "bluster reduced dimension '" + reduced_dim_name + "' was not found",
      )
      []
    }
  }
  let (observations, _) = bluster_matrix_shape(coordinates)
  if observations != experiment.col_names.length() {
    bluster_fail("bluster reduced dimension rows must match SCE column names")
  }
  let clustering = bluster_cluster_rows_graph(coordinates, config)
  let output = bluster_copy_sce(experiment)
  output.col_data[output_column] = clustering.clusters.map(fn(cluster) {
    (cluster + 1).to_string()
  })
  output.metadata["bluster.reduced_dim"] = reduced_dim_name
  output.metadata["bluster.graph"] = clustering.graph.graph_kind
  output.metadata["bluster.modularity"] = clustering.modularity.to_string()
  BlusterSceOutput::{ experiment: output, clustering }
}
