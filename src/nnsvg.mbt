// Nearest-neighbor Gaussian-process detection of spatially variable genes.
//
// This module follows the statistical model and result conventions of the
// Bioconductor nnSVG package:
//   y = X beta + w + epsilon
//   cov(w_i, w_j) = sigma_sq * exp(-distance(i, j) / length_scale)
//   epsilon ~ N(0, tau_sq I)
//
// The Gaussian-process likelihood is evaluated with a predecessor k-nearest
// neighbor (NNGP) factorization. Expression matrices use genes in rows and
// spatial spots in columns.

///|
pub suberror NnsvgError {
  NnsvgError(String)
}

///|
pub enum NnsvgOrdering {
  Ammd
  SumCoordinates
} derive(Eq, Debug)
pub extend NnsvgOrdering with Eq::{equal, not_equal}
pub extend NnsvgOrdering with @moonbitlang/core/debug.Debug::{to_repr}

///|
pub fn nnsvg_ammd_ordering() -> NnsvgOrdering {
  Ammd
}

///|
pub fn nnsvg_sum_coordinates_ordering() -> NnsvgOrdering {
  SumCoordinates
}

///|
pub struct NnsvgConfig {
  n_neighbors : Int
  ordering : NnsvgOrdering
  length_scale_min : Double
  length_scale_max : Double
  length_scale_grid : Int
  proportion_grid : Int
  refinement_steps : Int
  minimum_variance : Double
  fdr_threshold : Double
} derive(Eq, Debug)
pub extend NnsvgConfig with Eq::{equal, not_equal}
pub extend NnsvgConfig with @moonbitlang/core/debug.Debug::{to_repr}

///|
pub fn NnsvgConfig::create(
  n_neighbors? : Int = 10,
  ordering? : NnsvgOrdering = Ammd,
  length_scale_min? : Double = 0.02,
  length_scale_max? : Double = 2.0,
  length_scale_grid? : Int = 10,
  proportion_grid? : Int = 10,
  refinement_steps? : Int = 2,
  minimum_variance? : Double = 1.0e-10,
  fdr_threshold? : Double = 0.05,
) -> NnsvgConfig raise NnsvgError {
  if n_neighbors < 1 {
    raise NnsvgError("nnSVG neighbor count must be positive")
  }
  if !nnsvg_is_finite(length_scale_min) || length_scale_min <= 0.0 {
    raise NnsvgError("nnSVG minimum length scale must be finite and positive")
  }
  if !nnsvg_is_finite(length_scale_max) || length_scale_max <= length_scale_min {
    raise NnsvgError(
      "nnSVG maximum length scale must be finite and exceed the minimum",
    )
  }
  if length_scale_grid < 2 || proportion_grid < 2 {
    raise NnsvgError("nnSVG parameter grids must contain at least two points")
  }
  if refinement_steps < 0 {
    raise NnsvgError("nnSVG refinement step count must be non-negative")
  }
  if !nnsvg_is_finite(minimum_variance) || minimum_variance <= 0.0 {
    raise NnsvgError("nnSVG minimum variance must be finite and positive")
  }
  if !nnsvg_is_finite(fdr_threshold) ||
    fdr_threshold <= 0.0 ||
    fdr_threshold > 1.0 {
    raise NnsvgError("nnSVG FDR threshold must be finite and in (0, 1]")
  }
  NnsvgConfig::{
    n_neighbors,
    ordering,
    length_scale_min,
    length_scale_max,
    length_scale_grid,
    proportion_grid,
    refinement_steps,
    minimum_variance,
    fdr_threshold,
  }
}

///|
pub fn NnsvgConfig::default() -> NnsvgConfig {
  NnsvgConfig::{
    n_neighbors: 10,
    ordering: Ammd,
    length_scale_min: 0.02,
    length_scale_max: 2.0,
    length_scale_grid: 10,
    proportion_grid: 10,
    refinement_steps: 2,
    minimum_variance: 1.0e-10,
    fdr_threshold: 0.05,
  }
}

///|
pub struct NnsvgNeighborGraph {
  coordinates : Array[Array[Double]]
  order : Array[Int]
  inverse_order : Array[Int]
  neighbors : Array[Array[Int]]
  distances : Array[Array[Double]]
  n_neighbors : Int
  ordering : NnsvgOrdering
  range_scale : Double
} derive(Debug)
pub extend NnsvgNeighborGraph with @moonbitlang/core/debug.Debug::{to_repr}

///|
pub struct NnsvgGeneResult {
  gene_name : String
  sigma_sq : Double
  tau_sq : Double
  phi : Double
  length_scale : Double
  log_likelihood : Double
  mean : Double
  variance : Double
  spatial_covariance : Double
  proportion_spatial_variance : Double
  linear_log_likelihood : Double
  likelihood_ratio : Double
  rank : Double
  p_value : Double
  adjusted_p_value : Double
  beta : Array[Double]
  converged : Bool
} derive(Debug)
pub extend NnsvgGeneResult with @moonbitlang/core/debug.Debug::{to_repr}

///|
pub struct NnsvgResult {
  gene_names : Array[String]
  genes : Array[NnsvgGeneResult]
  graph : NnsvgNeighborGraph
  design_columns : Int
  config : NnsvgConfig
} derive(Debug)
pub extend NnsvgResult with @moonbitlang/core/debug.Debug::{to_repr}

///|
pub struct NnsvgFilterResult {
  counts : Array[Array[Double]]
  gene_names : Array[String]
  kept_indices : Array[Int]
  removed_low_expression : Array[Int]
  removed_mitochondrial : Array[Int]
} derive(Debug)
pub extend NnsvgFilterResult with @moonbitlang/core/debug.Debug::{to_repr}

///|
pub struct NnsvgSpatialFilterOutput {
  experiment : SpatialExperiment
  filter : NnsvgFilterResult
}

///|
pub struct NnsvgSpatialExperimentOutput {
  experiment : SpatialExperiment
  result : NnsvgResult
}

///|
priv struct NnsvgConditionalSystem {
  coefficients : Array[Array[Double]]
  conditional_variances : Array[Double]
}

///|
priv struct NnsvgProfileFit {
  beta : Array[Double]
  total_variance : Double
  proportion_spatial_variance : Double
  length_scale : Double
  log_likelihood : Double
}

///|
priv struct NnsvgRawGeneFit {
  gene_name : String
  sigma_sq : Double
  tau_sq : Double
  phi : Double
  length_scale : Double
  log_likelihood : Double
  mean : Double
  variance : Double
  spatial_covariance : Double
  proportion_spatial_variance : Double
  linear_log_likelihood : Double
  likelihood_ratio : Double
  p_value : Double
  beta : Array[Double]
  converged : Bool
}

///|
fn nnsvg_is_finite(value : Double) -> Bool {
  value == value && value.abs() <= 1.0e300
}

///|
fn nnsvg_zero_matrix(rows : Int, columns : Int) -> Array[Array[Double]] {
  let matrix : Array[Array[Double]] = []
  for _ in 0..<rows {
    matrix.push(Array::make(columns, 0.0))
  }
  matrix
}

///|
fn nnsvg_copy_matrix(matrix : Array[Array[Double]]) -> Array[Array[Double]] {
  let copy : Array[Array[Double]] = []
  for row in matrix {
    copy.push(row.copy())
  }
  copy
}

///|
fn nnsvg_distance(
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
fn nnsvg_validate_expression(
  expression : Array[Array[Double]],
) -> (Int, Int) raise NnsvgError {
  if expression.length() == 0 {
    raise NnsvgError("nnSVG expression must contain at least one gene")
  }
  let spots = expression[0].length()
  if spots < 3 {
    raise NnsvgError("nnSVG expression must contain at least three spots")
  }
  for row in expression {
    if row.length() != spots {
      raise NnsvgError("nnSVG expression must be rectangular")
    }
    for value in row {
      if !nnsvg_is_finite(value) {
        raise NnsvgError("nnSVG expression must contain finite values")
      }
    }
  }
  (expression.length(), spots)
}

///|
fn nnsvg_validate_coordinates(
  coordinates : Array[Array[Double]],
  expected_spots : Int,
) -> Int raise NnsvgError {
  if coordinates.length() != expected_spots {
    raise NnsvgError(
      "nnSVG coordinate rows must match the number of expression spots",
    )
  }
  if coordinates.length() == 0 || coordinates[0].length() < 2 {
    raise NnsvgError("nnSVG coordinates must contain at least two dimensions")
  }
  let dimensions = coordinates[0].length()
  for row in coordinates {
    if row.length() != dimensions {
      raise NnsvgError("nnSVG coordinates must be rectangular")
    }
    for value in row {
      if !nnsvg_is_finite(value) {
        raise NnsvgError("nnSVG coordinates must contain finite values")
      }
    }
  }
  dimensions
}

///|
fn nnsvg_prepare_gene_names(
  gene_names : Array[String],
  genes : Int,
) -> Array[String] raise NnsvgError {
  if gene_names.length() == 0 {
    let generated : Array[String] = []
    for index in 0..<genes {
      generated.push("gene_" + (index + 1).to_string())
    }
    return generated
  }
  if gene_names.length() != genes {
    raise NnsvgError(
      "nnSVG gene name count must match the expression row count",
    )
  }
  let prepared : Array[String] = []
  for index in 0..<genes {
    if gene_names[index].length() == 0 {
      prepared.push("gene_" + (index + 1).to_string())
    } else {
      prepared.push(gene_names[index])
    }
  }
  prepared
}

///|
fn nnsvg_cholesky(matrix : Array[Array[Double]]) -> Array[Array[Double]]? {
  let size = matrix.length()
  if size == 0 {
    return None
  }
  let lower = nnsvg_zero_matrix(size, size)
  for row in 0..<size {
    if matrix[row].length() != size {
      return None
    }
    for column in 0..<(row + 1) {
      let mut value = matrix[row][column]
      for inner in 0..<column {
        value = value - lower[row][inner] * lower[column][inner]
      }
      if row == column {
        if value <= 1.0e-12 || !nnsvg_is_finite(value) {
          return None
        }
        lower[row][column] = value.sqrt()
      } else {
        lower[row][column] = value / lower[column][column]
      }
    }
  }
  Some(lower)
}

///|
fn nnsvg_cholesky_solve(
  lower : Array[Array[Double]],
  right_hand_side : Array[Double],
) -> Array[Double] {
  let size = lower.length()
  let forward = Array::make(size, 0.0)
  for row in 0..<size {
    let mut value = right_hand_side[row]
    for column in 0..<row {
      value = value - lower[row][column] * forward[column]
    }
    forward[row] = value / lower[row][row]
  }
  let solution = Array::make(size, 0.0)
  let mut row = size - 1
  while row >= 0 {
    let mut value = forward[row]
    for column in (row + 1)..<size {
      value = value - lower[column][row] * solution[column]
    }
    solution[row] = value / lower[row][row]
    row = row - 1
  }
  solution
}

///|
fn nnsvg_prepare_design(
  design : Array[Array[Double]],
  spots : Int,
) -> Array[Array[Double]] raise NnsvgError {
  if design.length() == 0 {
    let intercept : Array[Array[Double]] = []
    for _ in 0..<spots {
      intercept.push([1.0])
    }
    return intercept
  }
  if design.length() != spots {
    raise NnsvgError("nnSVG design rows must match the number of spots")
  }
  if design[0].length() == 0 {
    raise NnsvgError("nnSVG design must contain at least one column")
  }
  let columns = design[0].length()
  if columns >= spots {
    raise NnsvgError("nnSVG design must have fewer columns than spatial spots")
  }
  let prepared : Array[Array[Double]] = []
  for row in design {
    if row.length() != columns {
      raise NnsvgError("nnSVG design must be rectangular")
    }
    let copied = row.copy()
    for value in copied {
      if !nnsvg_is_finite(value) {
        raise NnsvgError("nnSVG design must contain finite values")
      }
    }
    prepared.push(copied)
  }
  let information = nnsvg_zero_matrix(columns, columns)
  for row in prepared {
    for first in 0..<columns {
      for second in 0..<columns {
        information[first][second] = information[first][second] +
          row[first] * row[second]
      }
    }
  }
  match nnsvg_cholesky(information) {
    Some(_) => prepared
    None => raise NnsvgError("nnSVG design matrix must have full column rank")
  }
}

///|
fn nnsvg_scale_coordinates(
  coordinates : Array[Array[Double]],
) -> (Array[Array[Double]], Double) raise NnsvgError {
  let dimensions = coordinates[0].length()
  let minima = coordinates[0].copy()
  let maxima = coordinates[0].copy()
  for row in 1..<coordinates.length() {
    for dimension in 0..<dimensions {
      minima[dimension] = minima[dimension].min(coordinates[row][dimension])
      maxima[dimension] = maxima[dimension].max(coordinates[row][dimension])
    }
  }
  let mut range_scale = 0.0
  for dimension in 0..<dimensions {
    range_scale = range_scale.max(maxima[dimension] - minima[dimension])
  }
  if !nnsvg_is_finite(range_scale) || range_scale <= 1.0e-14 {
    raise NnsvgError(
      "nnSVG coordinates must vary in at least one spatial dimension",
    )
  }
  let scaled = nnsvg_zero_matrix(coordinates.length(), dimensions)
  for row in 0..<coordinates.length() {
    for dimension in 0..<dimensions {
      scaled[row][dimension] = (coordinates[row][dimension] - minima[dimension]) /
        range_scale
    }
  }
  (scaled, range_scale)
}

///|
fn nnsvg_sum_coordinate_order(coordinates : Array[Array[Double]]) -> Array[Int] {
  let order : Array[Int] = []
  for index in 0..<coordinates.length() {
    order.push(index)
  }
  order.sort_by(fn(left : Int, right : Int) -> Int {
    let mut left_sum = 0.0
    let mut right_sum = 0.0
    for dimension in 0..<coordinates[left].length() {
      left_sum = left_sum + coordinates[left][dimension]
      right_sum = right_sum + coordinates[right][dimension]
    }
    if left_sum < right_sum {
      -1
    } else if left_sum > right_sum {
      1
    } else {
      left - right
    }
  })
  order
}

///|
fn nnsvg_ammd_order(coordinates : Array[Array[Double]]) -> Array[Int] {
  let spots = coordinates.length()
  let order : Array[Int] = []
  let selected = Array::make(spots, false)
  let minimum_distances = Array::make(spots, 1.0e300)
  let sum_order = nnsvg_sum_coordinate_order(coordinates)
  let first = sum_order[0]
  order.push(first)
  selected[first] = true
  for index in 0..<spots {
    if index != first {
      minimum_distances[index] = nnsvg_distance(coordinates, index, first)
    }
  }
  while order.length() < spots {
    let mut best = -1
    let mut best_distance = -1.0
    for candidate in 0..<spots {
      if !selected[candidate] &&
        (
          minimum_distances[candidate] > best_distance + 1.0e-14 ||
          (
            (minimum_distances[candidate] - best_distance).abs() <= 1.0e-14 &&
            (best < 0 || candidate < best)
          )
        ) {
        best = candidate
        best_distance = minimum_distances[candidate]
      }
    }
    order.push(best)
    selected[best] = true
    for candidate in 0..<spots {
      if !selected[candidate] {
        minimum_distances[candidate] = minimum_distances[candidate].min(
          nnsvg_distance(coordinates, candidate, best),
        )
      }
    }
  }
  order
}

///|
pub fn nnsvg_build_neighbors(
  coordinates : Array[Array[Double]],
  n_neighbors? : Int = 10,
  ordering? : NnsvgOrdering = Ammd,
) -> NnsvgNeighborGraph raise NnsvgError {
  if n_neighbors < 1 {
    raise NnsvgError("nnSVG neighbor count must be positive")
  }
  if coordinates.length() < 3 {
    raise NnsvgError("nnSVG coordinates must contain at least three spots")
  }
  ignore(nnsvg_validate_coordinates(coordinates, coordinates.length()))
  let (scaled, range_scale) = nnsvg_scale_coordinates(coordinates)
  let order = match ordering {
    Ammd => nnsvg_ammd_order(scaled)
    SumCoordinates => nnsvg_sum_coordinate_order(scaled)
  }
  let spots = scaled.length()
  let inverse_order = Array::make(spots, 0)
  for position in 0..<spots {
    inverse_order[order[position]] = position
  }
  let neighbors : Array[Array[Int]] = []
  let distances : Array[Array[Double]] = []
  for _ in 0..<spots {
    neighbors.push([])
    distances.push([])
  }
  let retained_neighbors = n_neighbors.min(spots - 1)
  for position in 1..<spots {
    let original = order[position]
    let candidates : Array[(Double, Int)] = []
    for predecessor in 0..<position {
      let neighbor = order[predecessor]
      candidates.push((nnsvg_distance(scaled, original, neighbor), neighbor))
    }
    candidates.sort_by(fn(left : (Double, Int), right : (Double, Int)) -> Int {
      if left.0 < right.0 {
        -1
      } else if left.0 > right.0 {
        1
      } else {
        inverse_order[left.1] - inverse_order[right.1]
      }
    })
    let count = retained_neighbors.min(candidates.length())
    for index in 0..<count {
      distances[original].push(candidates[index].0)
      neighbors[original].push(candidates[index].1)
    }
  }
  NnsvgNeighborGraph::{
    coordinates: scaled,
    order,
    inverse_order,
    neighbors,
    distances,
    n_neighbors: retained_neighbors,
    ordering,
    range_scale,
  }
}

///|
fn nnsvg_conditional_system(
  graph : NnsvgNeighborGraph,
  proportion : Double,
  length_scale : Double,
  minimum_variance : Double,
) -> NnsvgConditionalSystem? {
  let spots = graph.coordinates.length()
  let coefficients : Array[Array[Double]] = []
  for _ in 0..<spots {
    coefficients.push([])
  }
  let conditional_variances = Array::make(spots, 1.0)
  for position in 0..<spots {
    let original = graph.order[position]
    let local_neighbors = graph.neighbors[original]
    if local_neighbors.length() > 0 {
      if proportion <= 0.0 {
        coefficients[original] = Array::make(local_neighbors.length(), 0.0)
      } else {
        let covariance = nnsvg_zero_matrix(
          local_neighbors.length(),
          local_neighbors.length(),
        )
        let cross_covariance = Array::make(local_neighbors.length(), 0.0)
        for first in 0..<local_neighbors.length() {
          covariance[first][first] = 1.0
          cross_covariance[first] = proportion *
            @math.exp(
              -nnsvg_distance(
                graph.coordinates,
                original,
                local_neighbors[first],
              ) /
              length_scale,
            )
          for second in 0..<first {
            let value = proportion *
              @math.exp(
                -nnsvg_distance(
                  graph.coordinates,
                  local_neighbors[first],
                  local_neighbors[second],
                ) /
                length_scale,
              )
            covariance[first][second] = value
            covariance[second][first] = value
          }
        }
        let lower = match nnsvg_cholesky(covariance) {
          Some(value) => value
          None => return None
        }
        let local_coefficients = nnsvg_cholesky_solve(lower, cross_covariance)
        let mut reduction = 0.0
        for index in 0..<local_neighbors.length() {
          reduction = reduction +
            cross_covariance[index] * local_coefficients[index]
        }
        let conditional_variance = 1.0 - reduction
        if !nnsvg_is_finite(conditional_variance) ||
          conditional_variance <= minimum_variance {
          return None
        }
        coefficients[original] = local_coefficients
        conditional_variances[original] = conditional_variance
      }
    }
  }
  Some(NnsvgConditionalSystem::{ coefficients, conditional_variances, })
}

///|
fn nnsvg_profile_likelihood(
  values : Array[Double],
  design : Array[Array[Double]],
  graph : NnsvgNeighborGraph,
  proportion : Double,
  length_scale : Double,
  minimum_variance : Double,
) -> NnsvgProfileFit? {
  let system = match
    nnsvg_conditional_system(graph, proportion, length_scale, minimum_variance) {
    Some(value) => value
    None => return None
  }
  let columns = design[0].length()
  let information = nnsvg_zero_matrix(columns, columns)
  let score = Array::make(columns, 0.0)
  let transformed_response = Array::make(values.length(), 0.0)
  let transformed_design = nnsvg_zero_matrix(values.length(), columns)
  let mut log_determinant = 0.0
  for position in 0..<graph.order.length() {
    let original = graph.order[position]
    let local_neighbors = graph.neighbors[original]
    let local_coefficients = system.coefficients[original]
    let mut response = values[original]
    for column in 0..<columns {
      transformed_design[original][column] = design[original][column]
    }
    for neighbor_index in 0..<local_neighbors.length() {
      let neighbor = local_neighbors[neighbor_index]
      let coefficient = local_coefficients[neighbor_index]
      response = response - coefficient * values[neighbor]
      for column in 0..<columns {
        transformed_design[original][column] = transformed_design[original][column] -
          coefficient * design[neighbor][column]
      }
    }
    transformed_response[original] = response
    let conditional_variance = system.conditional_variances[original]
    let weight = 1.0 / conditional_variance
    log_determinant = log_determinant + @math.ln(conditional_variance)
    for first in 0..<columns {
      score[first] = score[first] +
        transformed_design[original][first] * response * weight
      for second in 0..<columns {
        information[first][second] = information[first][second] +
          transformed_design[original][first] *
          transformed_design[original][second] *
          weight
      }
    }
  }
  let lower = match nnsvg_cholesky(information) {
    Some(value) => value
    None => return None
  }
  let beta = nnsvg_cholesky_solve(lower, score)
  let mut residual_sum = 0.0
  for position in 0..<graph.order.length() {
    let original = graph.order[position]
    let mut residual = transformed_response[original]
    for column in 0..<columns {
      residual = residual - transformed_design[original][column] * beta[column]
    }
    residual_sum = residual_sum +
      residual * residual / system.conditional_variances[original]
  }
  let total_variance = (residual_sum / values.length().to_double()).max(
    minimum_variance,
  )
  let log_likelihood = -0.5 *
    (
      values.length().to_double() *
      (1.8378770664093453 + 1.0 + @math.ln(total_variance)) +
      log_determinant
    )
  if !nnsvg_is_finite(log_likelihood) {
    return None
  }
  Some(NnsvgProfileFit::{
    beta,
    total_variance,
    proportion_spatial_variance: proportion,
    length_scale,
    log_likelihood,
  })
}

///|
fn nnsvg_update_best(
  best : NnsvgProfileFit?,
  candidate : NnsvgProfileFit?,
) -> NnsvgProfileFit? {
  match candidate {
    None => best
    Some(value) =>
      match best {
        None => Some(value)
        Some(current) =>
          if value.log_likelihood > current.log_likelihood {
            Some(value)
          } else {
            Some(current)
          }
      }
  }
}

///|
fn nnsvg_optimize_spatial(
  values : Array[Double],
  design : Array[Array[Double]],
  graph : NnsvgNeighborGraph,
  config : NnsvgConfig,
) -> (NnsvgProfileFit, NnsvgProfileFit, Bool) raise NnsvgError {
  let linear = match
    nnsvg_profile_likelihood(
      values,
      design,
      graph,
      0.0,
      config.length_scale_min,
      config.minimum_variance,
    ) {
    Some(value) => value
    None => raise NnsvgError("nnSVG could not fit the non-spatial model")
  }
  if linear.total_variance <= config.minimum_variance * 1.000001 {
    return (linear, linear, true)
  }
  let log_minimum = @math.ln(config.length_scale_min)
  let log_maximum = @math.ln(config.length_scale_max)
  let mut best : NnsvgProfileFit? = None
  for length_index in 0..<config.length_scale_grid {
    let length_fraction = length_index.to_double() /
      (config.length_scale_grid - 1).to_double()
    let length_scale = @math.exp(
      log_minimum + length_fraction * (log_maximum - log_minimum),
    )
    for proportion_index in 0..<config.proportion_grid {
      let proportion = 0.01 +
        0.98 *
        proportion_index.to_double() /
        (config.proportion_grid - 1).to_double()
      best = nnsvg_update_best(
        best,
        nnsvg_profile_likelihood(
          values,
          design,
          graph,
          proportion,
          length_scale,
          config.minimum_variance,
        ),
      )
    }
  }
  let initial = match best {
    Some(value) => value
    None => return (linear, linear, false)
  }
  let mut refined = initial
  let mut log_step = (log_maximum - log_minimum) /
    (config.length_scale_grid - 1).to_double()
  let mut proportion_step = 0.98 / (config.proportion_grid - 1).to_double()
  for _ in 0..<config.refinement_steps {
    log_step = log_step * 0.5
    proportion_step = proportion_step * 0.5
    let center_log_length = @math.ln(refined.length_scale)
    let mut local_best : NnsvgProfileFit? = Some(refined)
    for length_offset in 0..<3 {
      let candidate_log_length = (center_log_length +
        (length_offset - 1).to_double() * log_step)
        .max(log_minimum)
        .min(log_maximum)
      let candidate_length = @math.exp(candidate_log_length)
      for proportion_offset in 0..<3 {
        let candidate_proportion = (refined.proportion_spatial_variance +
          (proportion_offset - 1).to_double() * proportion_step)
          .max(0.001)
          .min(0.999)
        local_best = nnsvg_update_best(
          local_best,
          nnsvg_profile_likelihood(
            values,
            design,
            graph,
            candidate_proportion,
            candidate_length,
            config.minimum_variance,
          ),
        )
      }
    }
    match local_best {
      Some(value) => refined = value
      None => ()
    }
  }
  let selected = if refined.log_likelihood > linear.log_likelihood + 1.0e-10 {
    refined
  } else {
    linear
  }
  (selected, linear, true)
}

///|
fn nnsvg_mean(values : Array[Double]) -> Double {
  let mut total = 0.0
  for value in values {
    total = total + value
  }
  total / values.length().to_double()
}

///|
fn nnsvg_variance(values : Array[Double], mean : Double) -> Double {
  let mut total = 0.0
  for value in values {
    let difference = value - mean
    total = total + difference * difference
  }
  total / (values.length() - 1).to_double()
}

///|
fn nnsvg_bh_adjust(p_values : Array[Double]) -> Array[Double] {
  let count = p_values.length()
  let order : Array[Int] = []
  for index in 0..<count {
    order.push(index)
  }
  order.sort_by(fn(left : Int, right : Int) -> Int {
    if p_values[left] < p_values[right] {
      -1
    } else if p_values[left] > p_values[right] {
      1
    } else {
      left - right
    }
  })
  let adjusted = Array::make(count, 1.0)
  let mut running = 1.0
  let mut rank = count
  while rank > 0 {
    let original = order[rank - 1]
    let candidate = (p_values[original] * count.to_double() / rank.to_double()).min(
      1.0,
    )
    running = running.min(candidate)
    adjusted[original] = running.max(p_values[original])
    rank = rank - 1
  }
  adjusted
}

///|
fn nnsvg_likelihood_ranks(likelihood_ratios : Array[Double]) -> Array[Double] {
  let count = likelihood_ratios.length()
  let order : Array[Int] = []
  for index in 0..<count {
    order.push(index)
  }
  order.sort_by(fn(left : Int, right : Int) -> Int {
    if likelihood_ratios[left] > likelihood_ratios[right] {
      -1
    } else if likelihood_ratios[left] < likelihood_ratios[right] {
      1
    } else {
      left - right
    }
  })
  let ranks = Array::make(count, 0.0)
  let mut start = 0
  while start < count {
    let mut end = start + 1
    while end < count &&
          (likelihood_ratios[order[end]] - likelihood_ratios[order[start]]).abs() <=
          1.0e-12 {
      end = end + 1
    }
    let average_rank = ((start + 1).to_double() + end.to_double()) / 2.0
    for position in start..<end {
      ranks[order[position]] = average_rank
    }
    start = end
  }
  ranks
}

///|
pub fn nnsvg_with_config(
  expression : Array[Array[Double]],
  coordinates : Array[Array[Double]],
  config : NnsvgConfig,
  design? : Array[Array[Double]] = [],
  gene_names? : Array[String] = [],
) -> NnsvgResult raise NnsvgError {
  let (gene_count, spot_count) = nnsvg_validate_expression(expression)
  ignore(nnsvg_validate_coordinates(coordinates, spot_count))
  let prepared_names = nnsvg_prepare_gene_names(gene_names, gene_count)
  let prepared_design = nnsvg_prepare_design(design, spot_count)
  let graph = nnsvg_build_neighbors(
    coordinates,
    n_neighbors=config.n_neighbors,
    ordering=config.ordering,
  )
  let raw_results : Array[NnsvgRawGeneFit] = []
  let p_values = Array::make(gene_count, 1.0)
  let likelihood_ratios = Array::make(gene_count, 0.0)
  for gene in 0..<gene_count {
    let values = expression[gene]
    let mean = nnsvg_mean(values)
    let variance = nnsvg_variance(values, mean)
    let (spatial, linear, converged) = nnsvg_optimize_spatial(
      values, prepared_design, graph, config,
    )
    let likelihood_ratio = (2.0 *
    (spatial.log_likelihood - linear.log_likelihood)).max(0.0)
    let p_value = @math.exp(-0.5 * likelihood_ratio).max(0.0).min(1.0)
    let proportion = spatial.proportion_spatial_variance
    let sigma_sq = spatial.total_variance * proportion
    let tau_sq = spatial.total_variance * (1.0 - proportion)
    let length_scale = if proportion <= 0.0 {
      0.0
    } else {
      spatial.length_scale
    }
    let phi = if length_scale <= 0.0 { 0.0 } else { 1.0 / length_scale }
    p_values[gene] = p_value
    likelihood_ratios[gene] = likelihood_ratio
    raw_results.push(NnsvgRawGeneFit::{
      gene_name: prepared_names[gene],
      sigma_sq,
      tau_sq,
      phi,
      length_scale,
      log_likelihood: spatial.log_likelihood,
      mean,
      variance,
      spatial_covariance: sigma_sq,
      proportion_spatial_variance: proportion,
      linear_log_likelihood: linear.log_likelihood,
      likelihood_ratio,
      p_value,
      beta: spatial.beta,
      converged,
    })
  }
  let adjusted = nnsvg_bh_adjust(p_values)
  let ranks = nnsvg_likelihood_ranks(likelihood_ratios)
  let genes : Array[NnsvgGeneResult] = []
  for index in 0..<gene_count {
    let raw = raw_results[index]
    genes.push(NnsvgGeneResult::{
      gene_name: raw.gene_name,
      sigma_sq: raw.sigma_sq,
      tau_sq: raw.tau_sq,
      phi: raw.phi,
      length_scale: raw.length_scale,
      log_likelihood: raw.log_likelihood,
      mean: raw.mean,
      variance: raw.variance,
      spatial_covariance: raw.spatial_covariance,
      proportion_spatial_variance: raw.proportion_spatial_variance,
      linear_log_likelihood: raw.linear_log_likelihood,
      likelihood_ratio: raw.likelihood_ratio,
      rank: ranks[index],
      p_value: raw.p_value,
      adjusted_p_value: adjusted[index],
      beta: raw.beta,
      converged: raw.converged,
    })
  }
  NnsvgResult::{
    gene_names: prepared_names,
    genes,
    graph,
    design_columns: prepared_design[0].length(),
    config,
  }
}

///|
pub fn nnsvg(
  expression : Array[Array[Double]],
  coordinates : Array[Array[Double]],
  design? : Array[Array[Double]] = [],
  gene_names? : Array[String] = [],
  n_neighbors? : Int = 10,
) -> NnsvgResult raise NnsvgError {
  let config = NnsvgConfig::create(n_neighbors~)
  nnsvg_with_config(expression, coordinates, config, design~, gene_names~)
}

///|
pub fn NnsvgResult::n_genes(self : NnsvgResult) -> Int {
  self.genes.length()
}

///|
pub fn NnsvgResult::n_spots(self : NnsvgResult) -> Int {
  self.graph.coordinates.length()
}

///|
pub fn NnsvgResult::gene(
  self : NnsvgResult,
  gene_name : String,
) -> NnsvgGeneResult? {
  for gene in self.genes {
    if gene.gene_name == gene_name {
      return Some(gene)
    }
  }
  None
}

///|
pub fn NnsvgResult::significant(
  self : NnsvgResult,
  threshold? : Double = -1.0,
) -> Array[NnsvgGeneResult] {
  let selected_threshold = if threshold < 0.0 {
    self.config.fdr_threshold
  } else {
    threshold
  }
  let selected : Array[NnsvgGeneResult] = []
  for gene in self.genes {
    if gene.adjusted_p_value <= selected_threshold {
      selected.push(gene)
    }
  }
  selected.sort_by(fn(left : NnsvgGeneResult, right : NnsvgGeneResult) -> Int {
    if left.rank < right.rank {
      -1
    } else if left.rank > right.rank {
      1
    } else {
      0
    }
  })
  selected
}

///|
pub fn NnsvgResult::top(
  self : NnsvgResult,
  count : Int,
) -> Array[NnsvgGeneResult] {
  if count <= 0 {
    return []
  }
  let sorted = self.genes.copy()
  sorted.sort_by(fn(left : NnsvgGeneResult, right : NnsvgGeneResult) -> Int {
    if left.rank < right.rank {
      -1
    } else if left.rank > right.rank {
      1
    } else {
      0
    }
  })
  let output : Array[NnsvgGeneResult] = []
  for index in 0..<count.min(sorted.length()) {
    output.push(sorted[index])
  }
  output
}

///|
pub fn NnsvgResult::summary(self : NnsvgResult) -> String {
  "nnSVG(" +
  self.n_genes().to_string() +
  " genes x " +
  self.n_spots().to_string() +
  " spots, neighbors=" +
  self.graph.n_neighbors.to_string() +
  ", significant=" +
  self.significant().length().to_string() +
  ", FDR=" +
  self.config.fdr_threshold.to_string() +
  ")"
}

///|
pub fn nnsvg_filter_genes(
  counts : Array[Array[Double]],
  gene_names : Array[String],
  minimum_count? : Double = 3.0,
  minimum_spot_percentage? : Double = 0.5,
  filter_mitochondrial? : Bool = true,
) -> NnsvgFilterResult raise NnsvgError {
  if counts.length() == 0 || counts[0].length() == 0 {
    raise NnsvgError("nnSVG gene filtering requires a non-empty count matrix")
  }
  if gene_names.length() != counts.length() {
    raise NnsvgError(
      "nnSVG filtering gene names must match the count matrix rows",
    )
  }
  if !nnsvg_is_finite(minimum_count) || minimum_count < 0.0 {
    raise NnsvgError(
      "nnSVG filtering minimum count must be finite and non-negative",
    )
  }
  if !nnsvg_is_finite(minimum_spot_percentage) ||
    minimum_spot_percentage < 0.0 ||
    minimum_spot_percentage > 100.0 {
    raise NnsvgError(
      "nnSVG filtering spot percentage must be finite and in [0, 100]",
    )
  }
  let spots = counts[0].length()
  let required_spots = (minimum_spot_percentage / 100.0 * spots.to_double())
    .ceil()
    .to_int()
  let filtered_counts : Array[Array[Double]] = []
  let filtered_names : Array[String] = []
  let kept_indices : Array[Int] = []
  let removed_low_expression : Array[Int] = []
  let removed_mitochondrial : Array[Int] = []
  for gene in 0..<counts.length() {
    if counts[gene].length() != spots {
      raise NnsvgError("nnSVG filtering count matrix must be rectangular")
    }
    let mut detected = 0
    for value in counts[gene] {
      if !nnsvg_is_finite(value) || value < 0.0 {
        raise NnsvgError(
          "nnSVG filtering counts must contain finite non-negative values",
        )
      }
      if value >= minimum_count {
        detected = detected + 1
      }
    }
    let low_expression = detected < required_spots
    let mitochondrial = filter_mitochondrial &&
      (
        gene_names[gene].starts_with("MT-") ||
        gene_names[gene].starts_with("mt-")
      )
    if low_expression {
      removed_low_expression.push(gene)
    }
    if mitochondrial {
      removed_mitochondrial.push(gene)
    }
    if !low_expression && !mitochondrial {
      kept_indices.push(gene)
      filtered_counts.push(counts[gene].copy())
      filtered_names.push(gene_names[gene])
    }
  }
  NnsvgFilterResult::{
    counts: filtered_counts,
    gene_names: filtered_names,
    kept_indices,
    removed_low_expression,
    removed_mitochondrial,
  }
}

///|
fn nnsvg_copy_string_map(source : Map[String, String]) -> Map[String, String] {
  let copy : Map[String, String] = Map([])
  for key in source.keys() {
    copy[key] = source[key]
  }
  copy
}

///|
fn nnsvg_copy_spatial_experiment(
  experiment : SpatialExperiment,
) -> SpatialExperiment {
  let copy = SpatialExperiment::new()
  for key in experiment.assay.keys() {
    copy.assay[key] = nnsvg_copy_matrix(experiment.assay[key])
  }
  for row in experiment.row_data {
    copy.row_data.push(nnsvg_copy_string_map(row))
  }
  for column in experiment.col_data {
    copy.col_data.push(nnsvg_copy_string_map(column))
  }
  for coordinate in experiment.spatial_coords {
    copy.spatial_coords.push(coordinate)
  }
  for image in experiment.images {
    copy.images.push(image)
  }
  for key in experiment.metadata.keys() {
    copy.metadata[key] = experiment.metadata[key]
  }
  copy
}

///|
fn nnsvg_spatial_gene_names(
  experiment : SpatialExperiment,
  genes : Int,
) -> Array[String] raise NnsvgError {
  if experiment.row_data.length() != 0 && experiment.row_data.length() != genes {
    raise NnsvgError(
      "nnSVG SpatialExperiment rowData must match the assay rows",
    )
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
          if value.length() == 0 {
            match row.get("gene_id") {
              Some(identifier) =>
                if identifier.length() == 0 {
                  generated
                } else {
                  identifier
                }
              None => generated
            }
          } else {
            value
          }
        None =>
          match row.get("gene_id") {
            Some(identifier) =>
              if identifier.length() == 0 {
                generated
              } else {
                identifier
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
fn nnsvg_spatial_coordinates(
  experiment : SpatialExperiment,
) -> Array[Array[Double]] {
  let coordinates : Array[Array[Double]] = []
  let mut has_z_variation = false
  if experiment.spatial_coords.length() > 0 {
    let first_z = experiment.spatial_coords[0].z
    for coordinate in experiment.spatial_coords {
      if (coordinate.z - first_z).abs() > 1.0e-14 {
        has_z_variation = true
      }
    }
  }
  for coordinate in experiment.spatial_coords {
    if has_z_variation {
      coordinates.push([coordinate.x, coordinate.y, coordinate.z])
    } else {
      coordinates.push([coordinate.x, coordinate.y])
    }
  }
  coordinates
}

///|
pub fn nnsvg_filter_spatial_experiment(
  experiment : SpatialExperiment,
  count_assay? : String = "counts",
  minimum_count? : Double = 3.0,
  minimum_spot_percentage? : Double = 0.5,
  filter_mitochondrial? : Bool = true,
) -> NnsvgSpatialFilterOutput raise NnsvgError {
  let counts = match experiment.assay.get(count_assay) {
    Some(value) => value
    None =>
      raise NnsvgError(
        "nnSVG SpatialExperiment assay not found: " + count_assay,
      )
  }
  if counts.length() == 0 {
    raise NnsvgError("nnSVG SpatialExperiment count assay is empty")
  }
  let gene_names = nnsvg_spatial_gene_names(experiment, counts.length())
  let filter = nnsvg_filter_genes(
    counts,
    gene_names,
    minimum_count~,
    minimum_spot_percentage~,
    filter_mitochondrial~,
  )
  let filtered = SpatialExperiment::new()
  for assay_name in experiment.assay.keys() {
    let source = experiment.assay[assay_name]
    if source.length() != counts.length() {
      raise NnsvgError(
        "nnSVG SpatialExperiment assays must share the same gene rows",
      )
    }
    let subset : Array[Array[Double]] = []
    for index in filter.kept_indices {
      subset.push(source[index].copy())
    }
    filtered.assay[assay_name] = subset
  }
  for index in filter.kept_indices {
    if experiment.row_data.length() == 0 {
      filtered.row_data.push(Map([("gene_name", gene_names[index])]))
    } else {
      filtered.row_data.push(nnsvg_copy_string_map(experiment.row_data[index]))
    }
  }
  for column in experiment.col_data {
    filtered.col_data.push(nnsvg_copy_string_map(column))
  }
  for coordinate in experiment.spatial_coords {
    filtered.spatial_coords.push(coordinate)
  }
  for image in experiment.images {
    filtered.images.push(image)
  }
  for key in experiment.metadata.keys() {
    filtered.metadata[key] = experiment.metadata[key]
  }
  filtered.metadata["nnsvg_filter_kept"] = filter.kept_indices
    .length()
    .to_string()
  filtered.metadata["nnsvg_filter_removed"] = (counts.length() -
  filter.kept_indices.length()).to_string()
  NnsvgSpatialFilterOutput::{ experiment: filtered, filter, }
}

///|
pub fn nnsvg_spatial_experiment(
  experiment : SpatialExperiment,
  assay_name? : String = "logcounts",
  config? : NnsvgConfig = NnsvgConfig::default(),
  design? : Array[Array[Double]] = [],
) -> NnsvgSpatialExperimentOutput raise NnsvgError {
  let expression = match experiment.assay.get(assay_name) {
    Some(value) => value
    None =>
      raise NnsvgError("nnSVG SpatialExperiment assay not found: " + assay_name)
  }
  if expression.length() == 0 {
    raise NnsvgError("nnSVG SpatialExperiment assay is empty")
  }
  let gene_names = nnsvg_spatial_gene_names(experiment, expression.length())
  let coordinates = nnsvg_spatial_coordinates(experiment)
  let result = nnsvg_with_config(
    expression,
    coordinates,
    config,
    design~,
    gene_names~,
  )
  let enriched = nnsvg_copy_spatial_experiment(experiment)
  if enriched.row_data.length() == 0 {
    for gene_name in gene_names {
      enriched.row_data.push(Map([("gene_name", gene_name)]))
    }
  }
  for index in 0..<result.genes.length() {
    let gene = result.genes[index]
    let row = enriched.row_data[index]
    row["nnsvg_sigma_sq"] = gene.sigma_sq.to_string()
    row["nnsvg_tau_sq"] = gene.tau_sq.to_string()
    row["nnsvg_phi"] = gene.phi.to_string()
    row["nnsvg_length_scale"] = gene.length_scale.to_string()
    row["nnsvg_loglik"] = gene.log_likelihood.to_string()
    row["nnsvg_mean"] = gene.mean.to_string()
    row["nnsvg_var"] = gene.variance.to_string()
    row["nnsvg_spcov"] = gene.spatial_covariance.to_string()
    row["nnsvg_prop_sv"] = gene.proportion_spatial_variance.to_string()
    row["nnsvg_loglik_lm"] = gene.linear_log_likelihood.to_string()
    row["nnsvg_LR_stat"] = gene.likelihood_ratio.to_string()
    row["nnsvg_rank"] = gene.rank.to_string()
    row["nnsvg_pval"] = gene.p_value.to_string()
    row["nnsvg_padj"] = gene.adjusted_p_value.to_string()
  }
  enriched.metadata["nnsvg_assay"] = assay_name
  enriched.metadata["nnsvg_neighbors"] = result.graph.n_neighbors.to_string()
  enriched.metadata["nnsvg_genes"] = result.n_genes().to_string()
  enriched.metadata["nnsvg_significant"] = result
    .significant()
    .length()
    .to_string()
  NnsvgSpatialExperimentOutput::{ experiment: enriched, result, }
}

///|
pub fn nnsvg_example_data() -> (
  Array[Array[Double]],
  Array[Array[Double]],
  Array[String],
) {
  let coordinates = [
    [0.0, 0.0],
    [1.0, 0.0],
    [2.0, 0.0],
    [3.0, 0.0],
    [0.0, 1.0],
    [1.0, 1.0],
    [2.0, 1.0],
    [3.0, 1.0],
    [0.0, 2.0],
    [1.0, 2.0],
    [2.0, 2.0],
    [3.0, 2.0],
  ]
  let expression = [
    [0.0, 0.5, 1.0, 1.5, 0.5, 1.0, 1.5, 2.0, 1.0, 1.5, 2.0, 2.5],
    [2.5, 2.0, 1.5, 1.0, 2.0, 1.5, 1.0, 0.5, 1.5, 1.0, 0.5, 0.0],
    [1.0, 2.0, 1.0, 2.0, 2.0, 1.0, 2.0, 1.0, 1.0, 2.0, 1.0, 2.0],
    [3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0, 3.0],
  ]
  let gene_names = [
    "x_y_gradient", "reverse_gradient", "checkerboard", "constant",
  ]
  (expression, coordinates, gene_names)
}
