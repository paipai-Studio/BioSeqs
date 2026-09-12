// Spatially aware clustering inspired by the Bioconductor Banksy package.
//
// Expression matrices use genes in rows and spots in columns. Neighborhood
// harmonics follow Banksy's H0, H1, ... convention:
//   H0 = weighted neighborhood mean
//   Hm = |sum_j weight_ij * centered(x_j) * exp(i * m * phi_ij)|

///|
pub suberror BanksyError {
  BanksyError(String)
}

///|
pub(all) enum BanksySpatialMode {
  BanksyKnnMedian
  BanksyKnnInverseDistance
  BanksyKnnInversePower
  BanksyKnnRank
  BanksyKnnUniform
  BanksyRadiusGaussian
} derive(Eq, Debug)

///|
pub fn banksy_knn_median() -> BanksySpatialMode {
  BanksyKnnMedian
}

///|
pub fn banksy_knn_inverse_distance() -> BanksySpatialMode {
  BanksyKnnInverseDistance
}

///|
pub fn banksy_knn_inverse_power() -> BanksySpatialMode {
  BanksyKnnInversePower
}

///|
pub fn banksy_knn_rank() -> BanksySpatialMode {
  BanksyKnnRank
}

///|
pub fn banksy_knn_uniform() -> BanksySpatialMode {
  BanksyKnnUniform
}

///|
pub fn banksy_radius_gaussian() -> BanksySpatialMode {
  BanksyRadiusGaussian
}

///|
pub struct BanksyComputeConfig {
  max_harmonic : Int
  k_geom : Array[Int]
  spatial_mode : BanksySpatialMode
  distance_power : Double
  sigma : Double
  alpha : Double
  k_spatial : Int
  sample_size : Int
  sample_renormalize : Bool
  seed : Int
  center_harmonics : Bool
} derive(Debug)

///|
pub fn BanksyComputeConfig::create(
  max_harmonic? : Int = 0,
  k_geom? : Array[Int] = [15],
  spatial_mode? : BanksySpatialMode = BanksyKnnMedian,
  distance_power? : Double = 2.0,
  sigma? : Double = 1.5,
  alpha? : Double = 0.05,
  k_spatial? : Int = 100,
  sample_size? : Int = 0,
  sample_renormalize? : Bool = true,
  seed? : Int = 1,
  center_harmonics? : Bool = true,
) -> BanksyComputeConfig raise BanksyError {
  if max_harmonic < 0 {
    raise BanksyError("BANKSY maximum harmonic must be non-negative")
  }
  if k_geom.length() != 1 && k_geom.length() != max_harmonic + 1 {
    raise BanksyError(
      "BANKSY k_geom must contain one value or one value per harmonic",
    )
  }
  for value in k_geom {
    if value < 1 {
      raise BanksyError("BANKSY k_geom values must be positive")
    }
  }
  if !banksy_is_finite(distance_power) || distance_power <= 0.0 {
    raise BanksyError("BANKSY distance power must be finite and positive")
  }
  if !banksy_is_finite(sigma) || sigma <= 0.0 {
    raise BanksyError("BANKSY Gaussian sigma must be finite and positive")
  }
  if !banksy_is_finite(alpha) || alpha <= 0.0 || alpha >= 1.0 {
    raise BanksyError("BANKSY radial alpha must be finite and in (0, 1)")
  }
  if k_spatial < 1 {
    raise BanksyError("BANKSY radial neighbor count must be positive")
  }
  if sample_size < 0 {
    raise BanksyError("BANKSY sample size must be non-negative")
  }
  if seed < 0 {
    raise BanksyError("BANKSY seed must be non-negative")
  }
  BanksyComputeConfig::{
    max_harmonic,
    k_geom: k_geom.copy(),
    spatial_mode,
    distance_power,
    sigma,
    alpha,
    k_spatial,
    sample_size,
    sample_renormalize,
    seed,
    center_harmonics,
  }
}

///|
pub fn BanksyComputeConfig::default() -> BanksyComputeConfig {
  BanksyComputeConfig::{
    max_harmonic: 0,
    k_geom: [15],
    spatial_mode: BanksyKnnMedian,
    distance_power: 2.0,
    sigma: 1.5,
    alpha: 0.05,
    k_spatial: 100,
    sample_size: 0,
    sample_renormalize: true,
    seed: 1,
    center_harmonics: true,
  }
}

///|
pub struct BanksyNeighbor {
  index : Int
  distance : Double
  weight : Double
  angle : Double
} derive(Eq, Debug)

///|
pub struct BanksyModel {
  expression : Array[Array[Double]]
  coordinates : Array[Array[Double]]
  gene_names : Array[String]
  spot_names : Array[String]
  harmonics : Array[Array[Array[Double]]]
  neighbors : Array[Array[Array[BanksyNeighbor]]]
  config : BanksyComputeConfig
} derive(Debug)

///|
pub struct BanksyMatrix {
  values : Array[Array[Double]]
  feature_names : Array[String]
  component_weights : Array[Double]
  lambda : Double
  max_harmonic : Int
  scaled : Bool
} derive(Debug)

///|
pub struct BanksyPcaResult {
  scores : Array[Array[Double]]
  loadings : Array[Array[Double]]
  eigenvalues : Array[Double]
  variance_explained : Array[Double]
  n_components : Int
} derive(Debug)

///|
pub struct BanksyClusterConfig {
  n_clusters : Int
  n_starts : Int
  max_iterations : Int
  tolerance : Double
  seed : Int
} derive(Eq, Debug)

///|
pub fn BanksyClusterConfig::create(
  n_clusters? : Int = 5,
  n_starts? : Int = 10,
  max_iterations? : Int = 100,
  tolerance? : Double = 1.0e-6,
  seed? : Int = 1,
) -> BanksyClusterConfig raise BanksyError {
  if n_clusters < 1 {
    raise BanksyError("BANKSY cluster count must be positive")
  }
  if n_starts < 1 {
    raise BanksyError("BANKSY k-means start count must be positive")
  }
  if max_iterations < 1 {
    raise BanksyError("BANKSY k-means iteration limit must be positive")
  }
  if !banksy_is_finite(tolerance) || tolerance <= 0.0 {
    raise BanksyError("BANKSY k-means tolerance must be finite and positive")
  }
  if seed < 0 {
    raise BanksyError("BANKSY k-means seed must be non-negative")
  }
  BanksyClusterConfig::{
    n_clusters,
    n_starts,
    max_iterations,
    tolerance,
    seed,
  }
}

///|
pub fn BanksyClusterConfig::default() -> BanksyClusterConfig {
  BanksyClusterConfig::{
    n_clusters: 5,
    n_starts: 10,
    max_iterations: 100,
    tolerance: 1.0e-6,
    seed: 1,
  }
}

///|
pub struct BanksyClusterResult {
  labels : Array[Int]
  centroids : Array[Array[Double]]
  cluster_sizes : Array[Int]
  inertia : Double
  iterations : Int
  converged : Bool
  silhouette : Double
} derive(Debug)

///|
pub struct BanksySmoothingResult {
  labels : Array[Int]
  iterations : Int
  changed_labels : Int
  converged : Bool
} derive(Eq, Debug)

///|
pub struct BanksyRun {
  lambda : Double
  matrix : BanksyMatrix
  pca : BanksyPcaResult
  clustering : BanksyClusterResult
  smoothing : BanksySmoothingResult?
  spatial_agreement : Double
} derive(Debug)

///|
pub struct BanksySpatialExperimentOutput {
  experiment : SpatialExperiment
  model : BanksyModel
  run : BanksyRun
}

///|
priv struct BanksyDistanceIndex {
  index : Int
  distance : Double
}

///|
priv struct BanksyKmeansFit {
  labels : Array[Int]
  centroids : Array[Array[Double]]
  inertia : Double
  iterations : Int
  converged : Bool
}

///|
fn banksy_is_finite(value : Double) -> Bool {
  value == value && value.abs() <= 1.0e300
}

///|
fn banksy_zero_matrix(rows : Int, columns : Int) -> Array[Array[Double]] {
  let output : Array[Array[Double]] = []
  for _ in 0..<rows {
    output.push(Array::make(columns, 0.0))
  }
  output
}

///|
fn banksy_copy_matrix(matrix : Array[Array[Double]]) -> Array[Array[Double]] {
  let output : Array[Array[Double]] = []
  for row in matrix {
    output.push(row.copy())
  }
  output
}

///|
fn banksy_validate_names(
  names : Array[String],
  expected : Int,
  prefix : String,
  label : String,
) -> Array[String] raise BanksyError {
  if names.length() == 0 {
    let generated : Array[String] = []
    for index in 0..<expected {
      generated.push(prefix + (index + 1).to_string())
    }
    return generated
  }
  if names.length() != expected {
    raise BanksyError(
      "BANKSY " + label + " names do not match matrix dimensions",
    )
  }
  let seen : Map[String, Bool] = Map([])
  let copied : Array[String] = []
  for name in names {
    if name.length() == 0 {
      raise BanksyError("BANKSY " + label + " names must not be empty")
    }
    if seen.contains(name) {
      raise BanksyError("BANKSY " + label + " names must be unique")
    }
    seen[name] = true
    copied.push(name)
  }
  copied
}

///|
fn banksy_validate_input(
  expression : Array[Array[Double]],
  coordinates : Array[Array[Double]],
  gene_names : Array[String],
  spot_names : Array[String],
) -> (Int, Int, Array[String], Array[String]) raise BanksyError {
  if expression.length() == 0 {
    raise BanksyError("BANKSY expression must contain at least one gene")
  }
  let spots = expression[0].length()
  if spots < 2 {
    raise BanksyError("BANKSY expression must contain at least two spots")
  }
  for row in expression {
    if row.length() != spots {
      raise BanksyError("BANKSY expression must be rectangular")
    }
    for value in row {
      if !banksy_is_finite(value) {
        raise BanksyError("BANKSY expression values must be finite")
      }
    }
  }
  if coordinates.length() != spots {
    raise BanksyError("BANKSY coordinates must contain one row per spot")
  }
  let dimensions = coordinates[0].length()
  if dimensions < 2 {
    raise BanksyError("BANKSY coordinates require at least two dimensions")
  }
  for coordinate in coordinates {
    if coordinate.length() != dimensions {
      raise BanksyError("BANKSY coordinates must be rectangular")
    }
    for value in coordinate {
      if !banksy_is_finite(value) {
        raise BanksyError("BANKSY coordinates must be finite")
      }
    }
  }
  let genes = expression.length()
  let checked_genes = banksy_validate_names(gene_names, genes, "gene_", "gene")
  let checked_spots = banksy_validate_names(spot_names, spots, "spot_", "spot")
  (genes, spots, checked_genes, checked_spots)
}

///|
fn banksy_distance(
  coordinates : Array[Array[Double]],
  first : Int,
  second : Int,
) -> Double {
  let mut squared = 0.0
  for dimension in 0..<coordinates[first].length() {
    let difference = coordinates[first][dimension] -
      coordinates[second][dimension]
    squared = squared + difference * difference
  }
  squared.sqrt()
}

///|
fn banksy_angle(
  coordinates : Array[Array[Double]],
  from : Int,
  to : Int,
) -> Double {
  let dx = coordinates[to][0] - coordinates[from][0]
  let dy = coordinates[to][1] - coordinates[from][1]
  let angle = @math.atan2(dy, dx)
  if angle < 0.0 {
    angle + 2.0 * 3.141592653589793
  } else {
    angle
  }
}

///|
fn banksy_sorted_distances(
  coordinates : Array[Array[Double]],
  source : Int,
) -> Array[BanksyDistanceIndex] {
  let distances : Array[BanksyDistanceIndex] = []
  for target in 0..<coordinates.length() {
    if target != source {
      distances.push(BanksyDistanceIndex::{
        index: target,
        distance: banksy_distance(coordinates, source, target),
      })
    }
  }
  distances.sort_by(fn(
    left : BanksyDistanceIndex,
    right : BanksyDistanceIndex,
  ) -> Int {
    if left.distance < right.distance {
      -1
    } else if left.distance > right.distance {
      1
    } else if left.index < right.index {
      -1
    } else if left.index > right.index {
      1
    } else {
      0
    }
  })
  distances
}

///|
fn banksy_median(values : Array[Double]) -> Double {
  if values.length() == 0 {
    return 0.0
  }
  let sorted = values.copy()
  sorted.sort()
  let middle = sorted.length() / 2
  if sorted.length() % 2 == 1 {
    sorted[middle]
  } else {
    (sorted[middle - 1] + sorted[middle]) / 2.0
  }
}

///|
fn banksy_sample_order(
  length : Int,
  sample_size : Int,
  seed : Int,
  source : Int,
) -> Array[Int] {
  let order : Array[Int] = []
  for index in 0..<length {
    order.push(index)
  }
  if sample_size <= 0 || sample_size >= length {
    return order
  }
  let scores : Array[Int] = Array::make(length, 0)
  for index in 0..<length {
    let a = (index + 1) * 1103515245
    let b = (source + 1) * 12345
    let c = (seed + 1) * 48271
    scores[index] = (a + b + c) & 0x7fffffff
  }
  order.sort_by(fn(left : Int, right : Int) -> Int {
    if scores[left] < scores[right] {
      -1
    } else if scores[left] > scores[right] {
      1
    } else if left < right {
      -1
    } else if left > right {
      1
    } else {
      0
    }
  })
  let selected : Array[Int] = []
  for index in 0..<sample_size {
    selected.push(order[index])
  }
  selected.sort()
  selected
}

///|
fn banksy_normalize_neighbors(
  neighbors : Array[BanksyNeighbor],
) -> Array[BanksyNeighbor] {
  let mut total = 0.0
  for neighbor in neighbors {
    total = total + neighbor.weight
  }
  if total <= 0.0 {
    return neighbors
  }
  let normalized : Array[BanksyNeighbor] = []
  for neighbor in neighbors {
    normalized.push(BanksyNeighbor::{
      ..neighbor,
      weight: neighbor.weight / total,
    })
  }
  normalized
}

///|
fn banksy_global_nearest_median(coordinates : Array[Array[Double]]) -> Double {
  let nearest : Array[Double] = []
  for source in 0..<coordinates.length() {
    let sorted = banksy_sorted_distances(coordinates, source)
    if sorted.length() > 0 {
      nearest.push(sorted[0].distance)
    }
  }
  banksy_median(nearest).max(1.0e-12)
}

///|
fn banksy_neighbor_graph(
  coordinates : Array[Array[Double]],
  config : BanksyComputeConfig,
  harmonic : Int,
) -> Array[Array[BanksyNeighbor]] raise BanksyError {
  let spots = coordinates.length()
  let k_geom = if config.k_geom.length() == 1 {
    config.k_geom[0]
  } else {
    config.k_geom[harmonic]
  }
  let requested = match config.spatial_mode {
    BanksyRadiusGaussian => config.k_spatial
    _ => k_geom
  }
  if requested >= spots {
    raise BanksyError(
      "BANKSY neighbor count must be smaller than the number of spots",
    )
  }
  if config.sample_size > requested {
    raise BanksyError(
      "BANKSY sample size must not exceed the available neighbor count",
    )
  }
  let nearest_median = if config.spatial_mode is BanksyRadiusGaussian {
    banksy_global_nearest_median(coordinates)
  } else {
    1.0
  }
  let radius = config.sigma *
    (-coordinates[0].length().to_double() * @math.ln(config.alpha)).sqrt() *
    nearest_median
  let graph : Array[Array[BanksyNeighbor]] = []
  for source in 0..<spots {
    let sorted = banksy_sorted_distances(coordinates, source)
    let raw : Array[BanksyNeighbor] = []
    let local_distances : Array[Double] = []
    for rank in 0..<requested {
      local_distances.push(sorted[rank].distance)
    }
    let local_median = banksy_median(local_distances).max(1.0e-12)
    for rank in 0..<requested {
      let candidate = sorted[rank]
      if config.spatial_mode is BanksyRadiusGaussian &&
        candidate.distance >= radius {
        continue
      }
      let safe_distance = candidate.distance.max(1.0e-12)
      let weight = match config.spatial_mode {
        BanksyKnnMedian =>
          @math.exp(
            -(candidate.distance * candidate.distance) /
            (local_median * local_median),
          )
        BanksyKnnInverseDistance => 1.0 / safe_distance
        BanksyKnnInversePower =>
          1.0 / @math.pow(safe_distance, config.distance_power)
        BanksyKnnRank => {
          let rank_value = (rank + 1).to_double()
          let width = requested.to_double() / 1.5
          @math.exp(-(rank_value * rank_value) / (2.0 * width * width))
        }
        BanksyKnnUniform => 1.0
        BanksyRadiusGaussian =>
          @math.exp(
            -0.5 *
            (candidate.distance / (nearest_median * config.sigma)) *
            (candidate.distance / (nearest_median * config.sigma)),
          )
      }
      raw.push(BanksyNeighbor::{
        index: candidate.index,
        distance: candidate.distance,
        weight,
        angle: banksy_angle(coordinates, source, candidate.index),
      })
    }
    if raw.length() == 0 {
      let candidate = sorted[0]
      raw.push(BanksyNeighbor::{
        index: candidate.index,
        distance: candidate.distance,
        weight: 1.0,
        angle: banksy_angle(coordinates, source, candidate.index),
      })
    }
    let normalized = banksy_normalize_neighbors(raw)
    if config.sample_size > 0 && config.sample_size < normalized.length() {
      let order = banksy_sample_order(
        normalized.length(),
        config.sample_size,
        config.seed + harmonic * 7919,
        source,
      )
      let sampled : Array[BanksyNeighbor] = []
      for index in order {
        sampled.push(normalized[index])
      }
      if config.sample_renormalize {
        graph.push(banksy_normalize_neighbors(sampled))
      } else {
        graph.push(sampled)
      }
    } else {
      graph.push(normalized)
    }
  }
  graph
}

///|
fn banksy_compute_harmonic(
  expression : Array[Array[Double]],
  graph : Array[Array[BanksyNeighbor]],
  harmonic : Int,
  center : Bool,
) -> Array[Array[Double]] {
  let genes = expression.length()
  let spots = graph.length()
  let output = banksy_zero_matrix(genes, spots)
  for gene in 0..<genes {
    for spot in 0..<spots {
      let neighborhood = graph[spot]
      let mut local_mean = 0.0
      if harmonic > 0 && center {
        for neighbor in neighborhood {
          local_mean = local_mean + expression[gene][neighbor.index]
        }
        local_mean = local_mean / neighborhood.length().to_double()
      }
      if harmonic == 0 {
        let mut value = 0.0
        for neighbor in neighborhood {
          value = value + neighbor.weight * expression[gene][neighbor.index]
        }
        output[gene][spot] = value
      } else {
        let mut real = 0.0
        let mut imaginary = 0.0
        for neighbor in neighborhood {
          let value = expression[gene][neighbor.index] - local_mean
          let phase = harmonic.to_double() * neighbor.angle
          real = real + neighbor.weight * value * @math.cos(phase)
          imaginary = imaginary + neighbor.weight * value * @math.sin(phase)
        }
        output[gene][spot] = (real * real + imaginary * imaginary).sqrt()
      }
    }
  }
  output
}

///|
pub fn banksy_compute_with_config(
  expression : Array[Array[Double]],
  coordinates : Array[Array[Double]],
  config : BanksyComputeConfig,
  gene_names? : Array[String] = [],
  spot_names? : Array[String] = [],
) -> BanksyModel raise BanksyError {
  let (_, _, checked_genes, checked_spots) = banksy_validate_input(
    expression, coordinates, gene_names, spot_names,
  )
  let harmonics : Array[Array[Array[Double]]] = []
  let graphs : Array[Array[Array[BanksyNeighbor]]] = []
  for harmonic in 0..<=config.max_harmonic {
    let graph = banksy_neighbor_graph(coordinates, config, harmonic)
    let values = banksy_compute_harmonic(
      expression,
      graph,
      harmonic,
      harmonic > 0 && config.center_harmonics,
    )
    graphs.push(graph)
    harmonics.push(values)
  }
  BanksyModel::{
    expression: banksy_copy_matrix(expression),
    coordinates: banksy_copy_matrix(coordinates),
    gene_names: checked_genes,
    spot_names: checked_spots,
    harmonics,
    neighbors: graphs,
    config,
  }
}

///|
pub fn banksy_compute(
  expression : Array[Array[Double]],
  coordinates : Array[Array[Double]],
  k_geom? : Int = 15,
  compute_agf? : Bool = false,
  gene_names? : Array[String] = [],
  spot_names? : Array[String] = [],
) -> BanksyModel raise BanksyError {
  let config = BanksyComputeConfig::create(
    max_harmonic=if compute_agf { 1 } else { 0 },
    k_geom=[k_geom],
  )
  banksy_compute_with_config(
    expression,
    coordinates,
    config,
    gene_names~,
    spot_names~,
  )
}

///|
pub fn BanksyModel::n_genes(self : BanksyModel) -> Int {
  self.expression.length()
}

///|
pub fn BanksyModel::n_spots(self : BanksyModel) -> Int {
  if self.expression.length() == 0 {
    0
  } else {
    self.expression[0].length()
  }
}

///|
pub fn BanksyModel::max_harmonic(self : BanksyModel) -> Int {
  self.harmonics.length() - 1
}

///|
pub fn BanksyModel::harmonic(
  self : BanksyModel,
  order : Int,
) -> Array[Array[Double]]? {
  if order < 0 || order >= self.harmonics.length() {
    None
  } else {
    Some(banksy_copy_matrix(self.harmonics[order]))
  }
}

///|
pub fn BanksyModel::summary(self : BanksyModel) -> String {
  "BANKSY model: " +
  self.n_genes().to_string() +
  " genes, " +
  self.n_spots().to_string() +
  " spots, H0-H" +
  self.max_harmonic().to_string()
}

///|
fn banksy_validate_groups(
  groups : Array[String],
  spots : Int,
) -> Array[String] raise BanksyError {
  if groups.length() == 0 {
    return []
  }
  if groups.length() != spots {
    raise BanksyError("BANKSY scaling groups must contain one label per spot")
  }
  for group in groups {
    if group.length() == 0 {
      raise BanksyError("BANKSY scaling group labels must not be empty")
    }
  }
  groups.copy()
}

///|
fn banksy_scale_component(
  matrix : Array[Array[Double]],
  groups : Array[String],
) -> Array[Array[Double]] {
  let output = banksy_copy_matrix(matrix)
  if matrix.length() == 0 {
    return output
  }
  let spots = matrix[0].length()
  if groups.length() == 0 {
    for row in 0..<matrix.length() {
      let mut mean = 0.0
      for value in matrix[row] {
        mean = mean + value
      }
      mean = mean / spots.to_double()
      let mut sum_square = 0.0
      for value in matrix[row] {
        let difference = value - mean
        sum_square = sum_square + difference * difference
      }
      let standard_deviation = if spots > 1 {
        (sum_square / (spots - 1).to_double()).sqrt()
      } else {
        0.0
      }
      for spot in 0..<spots {
        output[row][spot] = if standard_deviation > 1.0e-14 {
          (matrix[row][spot] - mean) / standard_deviation
        } else {
          0.0
        }
      }
    }
    return output
  }
  let unique : Array[String] = []
  let seen : Map[String, Bool] = Map([])
  for group in groups {
    if !seen.contains(group) {
      seen[group] = true
      unique.push(group)
    }
  }
  for row in 0..<matrix.length() {
    for group in unique {
      let indices : Array[Int] = []
      for spot in 0..<spots {
        if groups[spot] == group {
          indices.push(spot)
        }
      }
      let mut mean = 0.0
      for spot in indices {
        mean = mean + matrix[row][spot]
      }
      mean = mean / indices.length().to_double()
      let mut sum_square = 0.0
      for spot in indices {
        let difference = matrix[row][spot] - mean
        sum_square = sum_square + difference * difference
      }
      let standard_deviation = if indices.length() > 1 {
        (sum_square / (indices.length() - 1).to_double()).sqrt()
      } else {
        0.0
      }
      for spot in indices {
        output[row][spot] = if standard_deviation > 1.0e-14 {
          (matrix[row][spot] - mean) / standard_deviation
        } else {
          0.0
        }
      }
    }
  }
  output
}

///|
pub fn banksy_component_weights(
  lambda : Double,
  max_harmonic : Int,
) -> Array[Double] raise BanksyError {
  if !banksy_is_finite(lambda) || lambda < 0.0 || lambda > 1.0 {
    raise BanksyError("BANKSY lambda must be finite and in [0, 1]")
  }
  if max_harmonic < 0 {
    raise BanksyError("BANKSY maximum harmonic must be non-negative")
  }
  let raw : Array[Double] = []
  let mut total = 0.0
  for harmonic in 0..<=max_harmonic {
    let value = @math.pow(2.0, -harmonic.to_double())
    raw.push(value)
    total = total + value
  }
  let weights : Array[Double] = [(1.0 - lambda).sqrt()]
  for value in raw {
    weights.push((lambda * value / total).sqrt())
  }
  weights
}

///|
pub fn banksy_get_matrix(
  model : BanksyModel,
  lambda? : Double = 0.2,
  max_harmonic? : Int = -1,
  scale? : Bool = false,
  groups? : Array[String] = [],
) -> BanksyMatrix raise BanksyError {
  let maximum = if max_harmonic < 0 {
    model.max_harmonic()
  } else {
    max_harmonic
  }
  if maximum < 0 || maximum > model.max_harmonic() {
    raise BanksyError("BANKSY requested harmonic was not computed")
  }
  let checked_groups = banksy_validate_groups(groups, model.n_spots())
  let weights = banksy_component_weights(lambda, maximum)
  let values : Array[Array[Double]] = []
  let feature_names : Array[String] = []
  let own = if scale {
    banksy_scale_component(model.expression, checked_groups)
  } else {
    banksy_copy_matrix(model.expression)
  }
  for gene in 0..<model.n_genes() {
    let row : Array[Double] = []
    for value in own[gene] {
      row.push(value * weights[0])
    }
    values.push(row)
    feature_names.push(model.gene_names[gene])
  }
  for harmonic in 0..<=maximum {
    let component = if scale {
      banksy_scale_component(model.harmonics[harmonic], checked_groups)
    } else {
      banksy_copy_matrix(model.harmonics[harmonic])
    }
    for gene in 0..<model.n_genes() {
      let row : Array[Double] = []
      for value in component[gene] {
        row.push(value * weights[harmonic + 1])
      }
      values.push(row)
      feature_names.push(model.gene_names[gene] + "_H" + harmonic.to_string())
    }
  }
  BanksyMatrix::{
    values,
    feature_names,
    component_weights: weights,
    lambda,
    max_harmonic: maximum,
    scaled: scale,
  }
}

///|
pub fn BanksyMatrix::n_features(self : BanksyMatrix) -> Int {
  self.values.length()
}

///|
pub fn BanksyMatrix::n_spots(self : BanksyMatrix) -> Int {
  if self.values.length() == 0 {
    0
  } else {
    self.values[0].length()
  }
}

///|
pub fn BanksyMatrix::spot_matrix(self : BanksyMatrix) -> Array[Array[Double]] {
  let output = banksy_zero_matrix(self.n_spots(), self.n_features())
  for feature in 0..<self.n_features() {
    for spot in 0..<self.n_spots() {
      output[spot][feature] = self.values[feature][spot]
    }
  }
  output
}

///|
fn banksy_center_columns(data : Array[Array[Double]]) -> Array[Array[Double]] {
  if data.length() == 0 {
    return []
  }
  let output = banksy_copy_matrix(data)
  for column in 0..<data[0].length() {
    let mut mean = 0.0
    for row in 0..<data.length() {
      mean = mean + data[row][column]
    }
    mean = mean / data.length().to_double()
    for row in 0..<data.length() {
      output[row][column] = data[row][column] - mean
    }
  }
  output
}

///|
fn banksy_identity(size : Int) -> Array[Array[Double]] {
  let matrix = banksy_zero_matrix(size, size)
  for index in 0..<size {
    matrix[index][index] = 1.0
  }
  matrix
}

///|
fn banksy_symmetric_eigen(
  input : Array[Array[Double]],
) -> (Array[Double], Array[Array[Double]]) {
  let size = input.length()
  let matrix = banksy_copy_matrix(input)
  let vectors = banksy_identity(size)
  let limit = (size * size * 30).max(60)
  for _ in 0..<limit {
    let mut p = 0
    let mut q = if size > 1 { 1 } else { 0 }
    let mut maximum = 0.0
    for row in 0..<size {
      for column in (row + 1)..<size {
        let value = matrix[row][column].abs()
        if value > maximum {
          maximum = value
          p = row
          q = column
        }
      }
    }
    if maximum < 1.0e-12 || size < 2 {
      break
    }
    let app = matrix[p][p]
    let aqq = matrix[q][q]
    let apq = matrix[p][q]
    let tau = (aqq - app) / (2.0 * apq)
    let t = if tau >= 0.0 {
      1.0 / (tau + (1.0 + tau * tau).sqrt())
    } else {
      -1.0 / (-tau + (1.0 + tau * tau).sqrt())
    }
    let cosine = 1.0 / (1.0 + t * t).sqrt()
    let sine = t * cosine
    for index in 0..<size {
      if index != p && index != q {
        let aip = matrix[index][p]
        let aiq = matrix[index][q]
        matrix[index][p] = cosine * aip - sine * aiq
        matrix[p][index] = matrix[index][p]
        matrix[index][q] = sine * aip + cosine * aiq
        matrix[q][index] = matrix[index][q]
      }
    }
    matrix[p][p] = app - t * apq
    matrix[q][q] = aqq + t * apq
    matrix[p][q] = 0.0
    matrix[q][p] = 0.0
    for index in 0..<size {
      let vip = vectors[index][p]
      let viq = vectors[index][q]
      vectors[index][p] = cosine * vip - sine * viq
      vectors[index][q] = sine * vip + cosine * viq
    }
  }
  let values : Array[Double] = []
  for index in 0..<size {
    values.push(matrix[index][index])
  }
  let order : Array[Int] = []
  for index in 0..<size {
    order.push(index)
  }
  order.sort_by(fn(left : Int, right : Int) -> Int {
    if values[left] > values[right] {
      -1
    } else if values[left] < values[right] {
      1
    } else if left < right {
      -1
    } else if left > right {
      1
    } else {
      0
    }
  })
  let sorted_values : Array[Double] = []
  let sorted_vectors = banksy_zero_matrix(size, size)
  for column in 0..<size {
    let source = order[column]
    sorted_values.push(values[source])
    for row in 0..<size {
      sorted_vectors[row][column] = vectors[row][source]
    }
  }
  (sorted_values, sorted_vectors)
}

///|
pub fn banksy_run_pca(
  matrix : BanksyMatrix,
  n_components? : Int = 20,
) -> BanksyPcaResult raise BanksyError {
  let spots = matrix.n_spots()
  let features = matrix.n_features()
  if spots < 2 || features < 1 {
    raise BanksyError("BANKSY PCA requires at least two spots and one feature")
  }
  let maximum = (spots - 1).min(features)
  if n_components < 1 || n_components > maximum {
    raise BanksyError(
      "BANKSY PCA component count must be in [1, min(spots - 1, features)]",
    )
  }
  let centered = banksy_center_columns(matrix.spot_matrix())
  let gram = banksy_zero_matrix(spots, spots)
  for first in 0..<spots {
    for second in first..<spots {
      let mut value = 0.0
      for feature in 0..<features {
        value = value + centered[first][feature] * centered[second][feature]
      }
      value = value / (spots - 1).to_double()
      gram[first][second] = value
      gram[second][first] = value
    }
  }
  let (all_values, all_vectors) = banksy_symmetric_eigen(gram)
  let mut total_variance = 0.0
  for value in all_values {
    if value > 0.0 {
      total_variance = total_variance + value
    }
  }
  if total_variance <= 1.0e-14 {
    raise BanksyError("BANKSY PCA input has no non-constant variation")
  }
  let component_count = n_components
  let scores = banksy_zero_matrix(spots, component_count)
  let loadings = banksy_zero_matrix(features, component_count)
  let eigenvalues : Array[Double] = []
  let variance_explained : Array[Double] = []
  for component in 0..<component_count {
    let eigenvalue = all_values[component].max(0.0)
    let singular = (eigenvalue * (spots - 1).to_double()).sqrt()
    eigenvalues.push(eigenvalue)
    variance_explained.push(eigenvalue / total_variance)
    for spot in 0..<spots {
      scores[spot][component] = all_vectors[spot][component] * singular
    }
    if singular > 1.0e-14 {
      for feature in 0..<features {
        let mut value = 0.0
        for spot in 0..<spots {
          value = value + centered[spot][feature] * all_vectors[spot][component]
        }
        loadings[feature][component] = value / singular
      }
    }
    let mut largest = 0
    for feature in 1..<features {
      if loadings[feature][component].abs() > loadings[largest][component].abs() {
        largest = feature
      }
    }
    if loadings[largest][component] < 0.0 {
      for spot in 0..<spots {
        scores[spot][component] = -scores[spot][component]
      }
      for feature in 0..<features {
        loadings[feature][component] = -loadings[feature][component]
      }
    }
  }
  BanksyPcaResult::{
    scores,
    loadings,
    eigenvalues,
    variance_explained,
    n_components: component_count,
  }
}

///|
pub fn BanksyPcaResult::n_spots(self : BanksyPcaResult) -> Int {
  self.scores.length()
}

///|
pub fn BanksyPcaResult::summary(self : BanksyPcaResult) -> String {
  let mut cumulative = 0.0
  for value in self.variance_explained {
    cumulative = cumulative + value
  }
  "BANKSY PCA: " +
  self.n_components.to_string() +
  " components, explained variance=" +
  cumulative.to_string()
}

///|
fn banksy_squared_distance(
  left : Array[Double],
  right : Array[Double],
) -> Double {
  let mut value = 0.0
  for index in 0..<left.length() {
    let difference = left[index] - right[index]
    value = value + difference * difference
  }
  value
}

///|
fn banksy_initialize_centroids(
  data : Array[Array[Double]],
  clusters : Int,
  first_index : Int,
) -> Array[Array[Double]] {
  let centroids : Array[Array[Double]] = [data[first_index].copy()]
  let selected = Array::make(data.length(), false)
  selected[first_index] = true
  while centroids.length() < clusters {
    let mut best_index = -1
    let mut best_distance = -1.0
    for point in 0..<data.length() {
      if selected[point] {
        continue
      }
      let mut minimum = banksy_squared_distance(data[point], centroids[0])
      for center in 1..<centroids.length() {
        minimum = minimum.min(
          banksy_squared_distance(data[point], centroids[center]),
        )
      }
      if minimum > best_distance {
        best_distance = minimum
        best_index = point
      }
    }
    if best_index < 0 {
      break
    }
    selected[best_index] = true
    centroids.push(data[best_index].copy())
  }
  centroids
}

///|
fn banksy_assign_clusters(
  data : Array[Array[Double]],
  centroids : Array[Array[Double]],
) -> (Array[Int], Double) {
  let labels = Array::make(data.length(), 0)
  let mut inertia = 0.0
  for point in 0..<data.length() {
    let mut best = 0
    let mut best_distance = banksy_squared_distance(data[point], centroids[0])
    for cluster in 1..<centroids.length() {
      let distance = banksy_squared_distance(data[point], centroids[cluster])
      if distance < best_distance {
        best = cluster
        best_distance = distance
      }
    }
    labels[point] = best
    inertia = inertia + best_distance
  }
  (labels, inertia)
}

///|
fn banksy_update_centroids(
  data : Array[Array[Double]],
  labels : Array[Int],
  old_centroids : Array[Array[Double]],
) -> Array[Array[Double]] {
  let clusters = old_centroids.length()
  let dimensions = data[0].length()
  let centroids = banksy_zero_matrix(clusters, dimensions)
  let counts = Array::make(clusters, 0)
  for point in 0..<data.length() {
    let cluster = labels[point]
    counts[cluster] = counts[cluster] + 1
    for dimension in 0..<dimensions {
      centroids[cluster][dimension] = centroids[cluster][dimension] +
        data[point][dimension]
    }
  }
  for cluster in 0..<clusters {
    if counts[cluster] > 0 {
      for dimension in 0..<dimensions {
        centroids[cluster][dimension] = centroids[cluster][dimension] /
          counts[cluster].to_double()
      }
    } else {
      let mut farthest = 0
      let mut farthest_distance = -1.0
      for point in 0..<data.length() {
        let distance = banksy_squared_distance(
          data[point],
          old_centroids[labels[point]],
        )
        if distance > farthest_distance {
          farthest = point
          farthest_distance = distance
        }
      }
      centroids[cluster] = data[farthest].copy()
    }
  }
  centroids
}

///|
fn banksy_kmeans_once(
  data : Array[Array[Double]],
  config : BanksyClusterConfig,
  start : Int,
) -> BanksyKmeansFit {
  let first = (config.seed % data.length() + start) % data.length()
  let mut centroids = banksy_initialize_centroids(
    data,
    config.n_clusters,
    first,
  )
  let mut labels = Array::make(data.length(), -1)
  let mut iterations = 0
  let mut converged = false
  for iteration in 0..<config.max_iterations {
    let (next_labels, _) = banksy_assign_clusters(data, centroids)
    let next_centroids = banksy_update_centroids(data, next_labels, centroids)
    let mut shift = 0.0
    for cluster in 0..<centroids.length() {
      shift = shift.max(
        banksy_squared_distance(centroids[cluster], next_centroids[cluster]).sqrt(),
      )
    }
    let mut changed = false
    for point in 0..<labels.length() {
      if labels[point] != next_labels[point] {
        changed = true
      }
    }
    labels = next_labels
    centroids = next_centroids
    iterations = iteration + 1
    if !changed || shift < config.tolerance {
      converged = true
      break
    }
  }
  let (final_labels, inertia) = banksy_assign_clusters(data, centroids)
  BanksyKmeansFit::{
    labels: final_labels,
    centroids,
    inertia,
    iterations,
    converged,
  }
}

///|
pub fn banksy_silhouette(
  data : Array[Array[Double]],
  labels : Array[Int],
) -> Double raise BanksyError {
  if data.length() != labels.length() {
    raise BanksyError("BANKSY silhouette labels must match data rows")
  }
  if data.length() < 2 {
    return 0.0
  }
  let mut maximum_label = -1
  for label in labels {
    if label < 0 {
      raise BanksyError("BANKSY cluster labels must be non-negative")
    }
    maximum_label = maximum_label.max(label)
  }
  let clusters = maximum_label + 1
  let mut total = 0.0
  for point in 0..<data.length() {
    let own = labels[point]
    let mut own_distance = 0.0
    let mut own_count = 0
    let other_sums = Array::make(clusters, 0.0)
    let other_counts = Array::make(clusters, 0)
    for other in 0..<data.length() {
      if other == point {
        continue
      }
      let distance = banksy_squared_distance(data[point], data[other]).sqrt()
      if labels[other] == own {
        own_distance = own_distance + distance
        own_count = own_count + 1
      } else {
        other_sums[labels[other]] = other_sums[labels[other]] + distance
        other_counts[labels[other]] = other_counts[labels[other]] + 1
      }
    }
    if own_count == 0 {
      continue
    }
    let a = own_distance / own_count.to_double()
    let mut b = 1.0e300
    for cluster in 0..<clusters {
      if cluster != own && other_counts[cluster] > 0 {
        b = b.min(other_sums[cluster] / other_counts[cluster].to_double())
      }
    }
    if b < 1.0e299 {
      let denominator = a.max(b)
      if denominator > 0.0 {
        total = total + (b - a) / denominator
      }
    }
  }
  total / data.length().to_double()
}

///|
pub fn banksy_cluster(
  pca : BanksyPcaResult,
  config? : BanksyClusterConfig = BanksyClusterConfig::default(),
) -> BanksyClusterResult raise BanksyError {
  if pca.scores.length() == 0 {
    raise BanksyError("BANKSY clustering requires non-empty PCA scores")
  }
  if config.n_clusters > pca.scores.length() {
    raise BanksyError("BANKSY cluster count must not exceed the spot count")
  }
  let mut best : BanksyKmeansFit? = None
  for start in 0..<config.n_starts {
    let fit = banksy_kmeans_once(pca.scores, config, start)
    best = match best {
      None => Some(fit)
      Some(current) =>
        if fit.inertia < current.inertia {
          Some(fit)
        } else {
          Some(current)
        }
    }
  }
  let fit = match best {
    Some(value) => value
    None => raise BanksyError("BANKSY k-means did not produce a fit")
  }
  let sizes = Array::make(config.n_clusters, 0)
  for label in fit.labels {
    sizes[label] = sizes[label] + 1
  }
  let silhouette = banksy_silhouette(pca.scores, fit.labels)
  BanksyClusterResult::{
    labels: fit.labels,
    centroids: fit.centroids,
    cluster_sizes: sizes,
    inertia: fit.inertia,
    iterations: fit.iterations,
    converged: fit.converged,
    silhouette,
  }
}

///|
pub fn BanksyClusterResult::n_clusters(self : BanksyClusterResult) -> Int {
  self.cluster_sizes.length()
}

///|
pub fn BanksyClusterResult::summary(self : BanksyClusterResult) -> String {
  "BANKSY clustering: " +
  self.n_clusters().to_string() +
  " clusters, inertia=" +
  self.inertia.to_string() +
  ", silhouette=" +
  self.silhouette.to_string()
}

///|
fn banksy_plain_knn(
  coordinates : Array[Array[Double]],
  k : Int,
) -> Array[Array[Int]] raise BanksyError {
  if k < 1 || k >= coordinates.length() {
    raise BanksyError(
      "BANKSY smoothing neighbor count must be in [1, spots - 1]",
    )
  }
  let output : Array[Array[Int]] = []
  for source in 0..<coordinates.length() {
    let sorted = banksy_sorted_distances(coordinates, source)
    let neighbors : Array[Int] = []
    for index in 0..<k {
      neighbors.push(sorted[index].index)
    }
    output.push(neighbors)
  }
  output
}

///|
pub fn banksy_smooth_labels(
  labels : Array[Int],
  coordinates : Array[Array[Double]],
  k? : Int = 15,
  proportion_threshold? : Double = 0.5,
  max_iterations? : Int = 10,
) -> BanksySmoothingResult raise BanksyError {
  if labels.length() != coordinates.length() || labels.length() < 2 {
    raise BanksyError(
      "BANKSY smoothing labels and coordinates must contain the same spots",
    )
  }
  if !banksy_is_finite(proportion_threshold) ||
    proportion_threshold < 0.0 ||
    proportion_threshold > 1.0 {
    raise BanksyError("BANKSY smoothing proportion threshold must be in [0, 1]")
  }
  if max_iterations == 0 || max_iterations < -1 {
    raise BanksyError(
      "BANKSY smoothing iterations must be positive or -1 for convergence",
    )
  }
  let mut maximum_label = -1
  for label in labels {
    if label < 0 {
      raise BanksyError("BANKSY smoothing labels must be non-negative")
    }
    maximum_label = maximum_label.max(label)
  }
  let neighbors = banksy_plain_knn(coordinates, k)
  let raw = labels.copy()
  let mut current = labels.copy()
  let limit = if max_iterations == -1 { 10000 } else { max_iterations }
  let mut iterations = 0
  let mut converged = false
  for iteration in 0..<limit {
    let next = current.copy()
    let mut changes = 0
    for spot in 0..<current.length() {
      let counts = Array::make(maximum_label + 1, 0)
      counts[current[spot]] = counts[current[spot]] + 1
      for neighbor in neighbors[spot] {
        counts[current[neighbor]] = counts[current[neighbor]] + 1
      }
      let mut best = current[spot]
      for label in 0..<counts.length() {
        if counts[label] > counts[best] {
          best = label
        }
      }
      let proportion = counts[best].to_double() / (k + 1).to_double()
      if proportion > proportion_threshold {
        next[spot] = best
      }
      if next[spot] != current[spot] {
        changes = changes + 1
      }
    }
    current = next
    iterations = iteration + 1
    if changes == 0 {
      converged = true
      break
    }
  }
  let mut total_changed = 0
  for spot in 0..<raw.length() {
    if raw[spot] != current[spot] {
      total_changed = total_changed + 1
    }
  }
  BanksySmoothingResult::{
    labels: current,
    iterations,
    changed_labels: total_changed,
    converged,
  }
}

///|
pub fn banksy_spatial_agreement(
  labels : Array[Int],
  coordinates : Array[Array[Double]],
  k? : Int = 6,
) -> Double raise BanksyError {
  if labels.length() != coordinates.length() || labels.length() < 2 {
    raise BanksyError("BANKSY spatial agreement labels must match coordinates")
  }
  let neighbors = banksy_plain_knn(coordinates, k)
  let mut matches = 0
  let mut total = 0
  for spot in 0..<labels.length() {
    for neighbor in neighbors[spot] {
      if labels[spot] == labels[neighbor] {
        matches = matches + 1
      }
      total = total + 1
    }
  }
  if total == 0 {
    0.0
  } else {
    matches.to_double() / total.to_double()
  }
}

///|
pub fn banksy_adjusted_rand_index(
  left : Array[Int],
  right : Array[Int],
) -> Double raise BanksyError {
  if left.length() != right.length() || left.length() < 2 {
    raise BanksyError(
      "BANKSY ARI label vectors must have equal length of at least two",
    )
  }
  let left_counts : Map[Int, Int] = Map([])
  let right_counts : Map[Int, Int] = Map([])
  let cells : Map[String, Int] = Map([])
  for index in 0..<left.length() {
    if left[index] < 0 || right[index] < 0 {
      raise BanksyError("BANKSY ARI labels must be non-negative")
    }
    left_counts[left[index]] = match left_counts.get(left[index]) {
      Some(value) => value + 1
      None => 1
    }
    right_counts[right[index]] = match right_counts.get(right[index]) {
      Some(value) => value + 1
      None => 1
    }
    let key = left[index].to_string() + ":" + right[index].to_string()
    cells[key] = match cells.get(key) {
      Some(value) => value + 1
      None => 1
    }
  }
  let choose_two = fn(value : Int) -> Double {
    (value * (value - 1) / 2).to_double()
  }
  let mut cell_sum = 0.0
  for key in cells.keys() {
    cell_sum = cell_sum + choose_two(cells[key])
  }
  let mut left_sum = 0.0
  for key in left_counts.keys() {
    left_sum = left_sum + choose_two(left_counts[key])
  }
  let mut right_sum = 0.0
  for key in right_counts.keys() {
    right_sum = right_sum + choose_two(right_counts[key])
  }
  let total_pairs = choose_two(left.length())
  let expected = left_sum * right_sum / total_pairs
  let maximum = (left_sum + right_sum) / 2.0
  if (maximum - expected).abs() <= 1.0e-14 {
    if (cell_sum - expected).abs() <= 1.0e-14 {
      1.0
    } else {
      0.0
    }
  } else {
    (cell_sum - expected) / (maximum - expected)
  }
}

///|
pub fn banksy_run(
  model : BanksyModel,
  lambda? : Double = 0.2,
  n_components? : Int = 20,
  cluster_config? : BanksyClusterConfig = BanksyClusterConfig::default(),
  max_harmonic? : Int = -1,
  scale? : Bool = true,
  groups? : Array[String] = [],
  smooth? : Bool = false,
  smoothing_k? : Int = 15,
  smoothing_threshold? : Double = 0.5,
  smoothing_iterations? : Int = 10,
) -> BanksyRun raise BanksyError {
  let matrix = banksy_get_matrix(model, lambda~, max_harmonic~, scale~, groups~)
  let pca = banksy_run_pca(matrix, n_components~)
  let clustering = banksy_cluster(pca, config=cluster_config)
  let smoothing = if smooth {
    Some(
      banksy_smooth_labels(
        clustering.labels,
        model.coordinates,
        k=smoothing_k,
        proportion_threshold=smoothing_threshold,
        max_iterations=smoothing_iterations,
      ),
    )
  } else {
    None
  }
  let final_labels = match smoothing {
    Some(result) => result.labels
    None => clustering.labels
  }
  let agreement_k = smoothing_k.min(model.n_spots() - 1).max(1)
  let agreement = banksy_spatial_agreement(
    final_labels,
    model.coordinates,
    k=agreement_k,
  )
  BanksyRun::{
    lambda,
    matrix,
    pca,
    clustering,
    smoothing,
    spatial_agreement: agreement,
  }
}

///|
pub fn banksy_parameter_sweep(
  model : BanksyModel,
  lambdas : Array[Double],
  n_components? : Int = 20,
  cluster_config? : BanksyClusterConfig = BanksyClusterConfig::default(),
  max_harmonic? : Int = -1,
  scale? : Bool = true,
  groups? : Array[String] = [],
) -> Array[BanksyRun] raise BanksyError {
  if lambdas.length() == 0 {
    raise BanksyError("BANKSY parameter sweep requires at least one lambda")
  }
  let runs : Array[BanksyRun] = []
  for lambda in lambdas {
    runs.push(
      banksy_run(
        model,
        lambda~,
        n_components~,
        cluster_config~,
        max_harmonic~,
        scale~,
        groups~,
      ),
    )
  }
  runs
}

///|
pub fn BanksyRun::final_labels(self : BanksyRun) -> Array[Int] {
  match self.smoothing {
    Some(result) => result.labels.copy()
    None => self.clustering.labels.copy()
  }
}

///|
pub fn BanksyRun::summary(self : BanksyRun) -> String {
  "BANKSY run: lambda=" +
  self.lambda.to_string() +
  ", " +
  self.clustering.summary() +
  ", spatial agreement=" +
  self.spatial_agreement.to_string()
}

///|
fn banksy_copy_string_map(source : Map[String, String]) -> Map[String, String] {
  let copied : Map[String, String] = Map([])
  for key in source.keys() {
    copied[key] = source[key]
  }
  copied
}

///|
fn banksy_copy_spatial_experiment(
  experiment : SpatialExperiment,
) -> SpatialExperiment {
  let copied = SpatialExperiment::new()
  for key in experiment.assay.keys() {
    copied.assay[key] = banksy_copy_matrix(experiment.assay[key])
  }
  for row in experiment.row_data {
    copied.row_data.push(banksy_copy_string_map(row))
  }
  for column in experiment.col_data {
    copied.col_data.push(banksy_copy_string_map(column))
  }
  for coordinate in experiment.spatial_coords {
    copied.spatial_coords.push(coordinate)
  }
  for image in experiment.images {
    copied.images.push(
      SpatialImage::new(
        image.id,
        banksy_copy_matrix(image.data),
        image.scale_factor,
      ),
    )
  }
  for key in experiment.metadata.keys() {
    copied.metadata[key] = experiment.metadata[key]
  }
  copied
}

///|
fn banksy_experiment_gene_names(
  experiment : SpatialExperiment,
  genes : Int,
) -> Array[String] raise BanksyError {
  if experiment.row_data.length() != 0 && experiment.row_data.length() != genes {
    raise BanksyError("BANKSY SpatialExperiment rowData must match assay rows")
  }
  let names : Array[String] = []
  for gene in 0..<genes {
    let generated = "gene_" + (gene + 1).to_string()
    if experiment.row_data.length() == 0 {
      names.push(generated)
    } else {
      let row = experiment.row_data[gene]
      let name = match row.get("gene_name") {
        Some(value) =>
          if value.length() > 0 {
            value
          } else {
            match row.get("gene_id") {
              Some(identifier) =>
                if identifier.length() > 0 {
                  identifier
                } else {
                  generated
                }
              None => generated
            }
          }
        None =>
          match row.get("gene_id") {
            Some(identifier) =>
              if identifier.length() > 0 {
                identifier
              } else {
                generated
              }
            None => generated
          }
      }
      names.push(name)
    }
  }
  names
}

///|
fn banksy_experiment_spot_names(
  experiment : SpatialExperiment,
  spots : Int,
) -> Array[String] raise BanksyError {
  if experiment.col_data.length() != 0 && experiment.col_data.length() != spots {
    raise BanksyError(
      "BANKSY SpatialExperiment colData must match assay columns",
    )
  }
  let names : Array[String] = []
  for spot in 0..<spots {
    let generated = "spot_" + (spot + 1).to_string()
    if experiment.col_data.length() == 0 {
      names.push(generated)
    } else {
      let column = experiment.col_data[spot]
      let name = match column.get("spot_id") {
        Some(value) => if value.length() > 0 { value } else { generated }
        None => generated
      }
      names.push(name)
    }
  }
  names
}

///|
fn banksy_experiment_coordinates(
  experiment : SpatialExperiment,
  spots : Int,
) -> Array[Array[Double]] raise BanksyError {
  if experiment.spatial_coords.length() != spots {
    raise BanksyError(
      "BANKSY SpatialExperiment requires one spatial coordinate per spot",
    )
  }
  let mut use_z = false
  if spots > 0 {
    let first = experiment.spatial_coords[0].z
    for coordinate in experiment.spatial_coords {
      if (coordinate.z - first).abs() > 1.0e-14 {
        use_z = true
      }
    }
  }
  let coordinates : Array[Array[Double]] = []
  for coordinate in experiment.spatial_coords {
    if use_z {
      coordinates.push([coordinate.x, coordinate.y, coordinate.z])
    } else {
      coordinates.push([coordinate.x, coordinate.y])
    }
  }
  coordinates
}

///|
pub fn banksy_spatial_experiment(
  experiment : SpatialExperiment,
  assay_name? : String = "logcounts",
  config? : BanksyComputeConfig = BanksyComputeConfig::default(),
  lambda? : Double = 0.2,
  n_components? : Int = 20,
  cluster_config? : BanksyClusterConfig = BanksyClusterConfig::default(),
  max_harmonic? : Int = -1,
  scale? : Bool = true,
  group_key? : String = "",
  smooth? : Bool = true,
  smoothing_k? : Int = 15,
  smoothing_threshold? : Double = 0.5,
  smoothing_iterations? : Int = 10,
) -> BanksySpatialExperimentOutput raise BanksyError {
  let expression = match experiment.assay.get(assay_name) {
    Some(value) => value
    None =>
      raise BanksyError(
        "BANKSY SpatialExperiment assay not found: " + assay_name,
      )
  }
  if expression.length() == 0 {
    raise BanksyError("BANKSY SpatialExperiment assay must not be empty")
  }
  let spots = expression[0].length()
  let gene_names = banksy_experiment_gene_names(experiment, expression.length())
  let spot_names = banksy_experiment_spot_names(experiment, spots)
  let coordinates = banksy_experiment_coordinates(experiment, spots)
  let groups : Array[String] = []
  if group_key.length() > 0 {
    if experiment.col_data.length() != spots {
      raise BanksyError(
        "BANKSY grouping requires complete SpatialExperiment colData",
      )
    }
    for column in experiment.col_data {
      match column.get(group_key) {
        Some(value) =>
          if value.length() == 0 {
            raise BanksyError(
              "BANKSY SpatialExperiment grouping values must not be empty",
            )
          } else {
            groups.push(value)
          }
        None =>
          raise BanksyError(
            "BANKSY SpatialExperiment group field not found: " + group_key,
          )
      }
    }
  }
  let model = banksy_compute_with_config(
    expression,
    coordinates,
    config,
    gene_names~,
    spot_names~,
  )
  let run = banksy_run(
    model,
    lambda~,
    n_components~,
    cluster_config~,
    max_harmonic~,
    scale~,
    groups~,
    smooth~,
    smoothing_k~,
    smoothing_threshold~,
    smoothing_iterations~,
  )
  let enriched = banksy_copy_spatial_experiment(experiment)
  if enriched.row_data.length() == 0 {
    for gene_name in gene_names {
      enriched.row_data.push(Map([("gene_name", gene_name)]))
    }
  }
  if enriched.col_data.length() == 0 {
    for spot_name in spot_names {
      enriched.col_data.push(Map([("spot_id", spot_name)]))
    }
  }
  for harmonic in 0..<=model.max_harmonic() {
    enriched.assay["H" + harmonic.to_string()] = banksy_copy_matrix(
      model.harmonics[harmonic],
    )
  }
  let labels = run.final_labels()
  for spot in 0..<spots {
    enriched.col_data[spot]["banksy_cluster"] = (labels[spot] + 1).to_string()
    enriched.col_data[spot]["banksy_cluster_raw"] = (run.clustering.labels[spot] +
    1).to_string()
  }
  enriched.metadata["banksy_assay"] = assay_name
  enriched.metadata["banksy_lambda"] = lambda.to_string()
  enriched.metadata["banksy_max_harmonic"] = model.max_harmonic().to_string()
  enriched.metadata["banksy_clusters"] = run.clustering.n_clusters().to_string()
  enriched.metadata["banksy_spatial_agreement"] = run.spatial_agreement.to_string()
  BanksySpatialExperimentOutput::{ experiment: enriched, model, run, }
}

///|
pub fn banksy_example_data() -> (
  Array[Array[Double]],
  Array[Array[Double]],
  Array[String],
  Array[String],
) {
  let width = 6
  let height = 4
  let spots = width * height
  let coordinates : Array[Array[Double]] = []
  let spot_names : Array[String] = []
  let expression = banksy_zero_matrix(5, spots)
  for y in 0..<height {
    for x in 0..<width {
      let spot = y * width + x
      coordinates.push([x.to_double(), y.to_double()])
      spot_names.push("spot_" + (spot + 1).to_string())
      let left = if x < width / 2 { 8.0 } else { 1.0 }
      let right = if x >= width / 2 { 8.0 } else { 1.0 }
      let upper = if y < height / 2 { 6.0 } else { 2.0 }
      let boundary = if x == 2 || x == 3 { 7.0 } else { 1.0 }
      let variation = ((x * 3 + y * 5) % 4).to_double() * 0.15
      expression[0][spot] = left + variation
      expression[1][spot] = right + variation
      expression[2][spot] = upper + variation
      expression[3][spot] = boundary + variation
      expression[4][spot] = 3.0 + variation
    }
  }
  (
    expression,
    coordinates,
    ["left_marker", "right_marker", "upper_marker", "boundary", "baseline"],
    spot_names,
  )
}
