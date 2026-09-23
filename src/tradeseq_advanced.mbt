// Advanced trajectory differential expression inspired by tradeSeq 1.27.0.
//
// Count matrices use feature x cell orientation. Pseudotime and lineage
// weights use cell x lineage orientation, matching slingshot output.

///|
pub suberror TradeSeqAdvancedError {
  TradeSeqAdvancedError(String)
}

///|
pub struct TradeSeqAdvancedConfig {
  n_knots : Int
  smoothing_penalty : Double
  max_iterations : Int
  tolerance : Double
  minimum_dispersion : Double
  maximum_dispersion : Double
  ridge : Double
  fdr_threshold : Double
  test_points : Int
} derive(Debug)
pub extend TradeSeqAdvancedConfig with @moonbitlang/core/debug.Debug::{to_repr}

///|
pub fn TradeSeqAdvancedConfig::create(
  n_knots? : Int = 6,
  smoothing_penalty? : Double = 1.0,
  max_iterations? : Int = 60,
  tolerance? : Double = 1.0e-6,
  minimum_dispersion? : Double = 1.0e-4,
  maximum_dispersion? : Double = 100.0,
  ridge? : Double = 1.0e-8,
  fdr_threshold? : Double = 0.05,
  test_points? : Int = 12,
) -> TradeSeqAdvancedConfig raise TradeSeqAdvancedError {
  if n_knots < 3 {
    raise TradeSeqAdvancedError("n_knots must be at least 3")
  }
  if !tradeseq_adv_is_finite(smoothing_penalty) || smoothing_penalty < 0.0 {
    raise TradeSeqAdvancedError(
      "smoothing_penalty must be finite and non-negative",
    )
  }
  if max_iterations < 1 {
    raise TradeSeqAdvancedError("max_iterations must be positive")
  }
  if !tradeseq_adv_is_finite(tolerance) || tolerance <= 0.0 {
    raise TradeSeqAdvancedError("tolerance must be finite and positive")
  }
  if !tradeseq_adv_is_finite(minimum_dispersion) || minimum_dispersion <= 0.0 {
    raise TradeSeqAdvancedError(
      "minimum_dispersion must be finite and positive",
    )
  }
  if !tradeseq_adv_is_finite(maximum_dispersion) ||
    maximum_dispersion < minimum_dispersion {
    raise TradeSeqAdvancedError(
      "maximum_dispersion must be finite and at least minimum_dispersion",
    )
  }
  if !tradeseq_adv_is_finite(ridge) || ridge <= 0.0 {
    raise TradeSeqAdvancedError("ridge must be finite and positive")
  }
  if !tradeseq_adv_is_finite(fdr_threshold) ||
    fdr_threshold <= 0.0 ||
    fdr_threshold > 1.0 {
    raise TradeSeqAdvancedError("fdr_threshold must be in (0, 1]")
  }
  if test_points < 2 {
    raise TradeSeqAdvancedError("test_points must be at least 2")
  }
  TradeSeqAdvancedConfig::{
    n_knots,
    smoothing_penalty,
    max_iterations,
    tolerance,
    minimum_dispersion,
    maximum_dispersion,
    ridge,
    fdr_threshold,
    test_points,
  }
}

///|
pub fn TradeSeqAdvancedConfig::default() -> TradeSeqAdvancedConfig {
  TradeSeqAdvancedConfig::{
    n_knots: 6,
    smoothing_penalty: 1.0,
    max_iterations: 60,
    tolerance: 1.0e-6,
    minimum_dispersion: 1.0e-4,
    maximum_dispersion: 100.0,
    ridge: 1.0e-8,
    fdr_threshold: 0.05,
    test_points: 12,
  }
}

///|
pub struct TradeSeqGeneModel {
  gene_id : String
  coefficients : Array[Double]
  covariance : Array[Array[Double]]
  fitted_values : Array[Double]
  dispersion : Double
  log_likelihood : Double
  aic : Double
  effective_df : Double
  converged : Bool
  iterations : Int
} derive(Debug)
pub extend TradeSeqGeneModel with @moonbitlang/core/debug.Debug::{to_repr}

///|
pub struct TradeSeqAdvancedFit {
  gene_names : Array[String]
  lineage_names : Array[String]
  models : Array[TradeSeqGeneModel]
  pseudotime : Array[Array[Double]]
  cell_weights : Array[Array[Double]]
  offsets : Array[Double]
  lineage_minimum : Array[Double]
  lineage_maximum : Array[Double]
  config : TradeSeqAdvancedConfig
  n_genes : Int
  n_cells : Int
  n_lineages : Int
} derive(Debug)
pub extend TradeSeqAdvancedFit with @moonbitlang/core/debug.Debug::{to_repr}

///|
pub struct TradeSeqAdvancedTestResult {
  gene_id : String
  wald_statistic : Double
  degrees_freedom : Int
  p_value : Double
  adjusted_p_value : Double
  log2_fold_change : Double
  significant : Bool
} derive(Debug)
pub extend TradeSeqAdvancedTestResult with @moonbitlang/core/debug.Debug::{to_repr}

///|
pub struct TradeSeqAdvancedTestTable {
  test_name : String
  results : Array[TradeSeqAdvancedTestResult]
  fdr_threshold : Double
} derive(Debug)
pub extend TradeSeqAdvancedTestTable with @moonbitlang/core/debug.Debug::{to_repr}

///|
pub struct TradeSeqKnotEvaluation {
  candidates : Array[Int]
  gene_aic : Array[Array[Double]]
  mean_aic : Array[Double]
  selected_n_knots : Int
} derive(Debug)
pub extend TradeSeqKnotEvaluation with @moonbitlang/core/debug.Debug::{to_repr}

///|
pub struct TradeSeqSmoothPrediction {
  gene_id : String
  pseudotime : Array[Array[Double]]
  fitted : Array[Array[Double]]
  standard_errors : Array[Array[Double]]
} derive(Debug)
pub extend TradeSeqSmoothPrediction with @moonbitlang/core/debug.Debug::{to_repr}

///|
pub struct TradeSeqAdvancedSceOutput {
  experiment : SingleCellExperiment
  fit : TradeSeqAdvancedFit
}

///|
priv struct TradeSeqAdvObservation {
  cell : Int
  lineage : Int
  weight : Double
  basis : Array[Double]
}

///|
fn tradeseq_adv_is_finite(value : Double) -> Bool {
  value == value && value.abs() < 1.0e300
}

///|
fn tradeseq_adv_clamp(value : Double, lower : Double, upper : Double) -> Double {
  if value < lower {
    lower
  } else if value > upper {
    upper
  } else {
    value
  }
}

///|
fn tradeseq_adv_zero_matrix(rows : Int, columns : Int) -> Array[Array[Double]] {
  let output : Array[Array[Double]] = []
  for _ in 0..<rows {
    output.push(Array::make(columns, 0.0))
  }
  output
}

///|
fn tradeseq_adv_copy_matrix(
  matrix : Array[Array[Double]],
) -> Array[Array[Double]] {
  let output : Array[Array[Double]] = []
  for row in matrix {
    output.push(row.copy())
  }
  output
}

///|
fn tradeseq_adv_prepare_names(
  supplied : Array[String],
  length : Int,
  prefix : String,
) -> Array[String] raise TradeSeqAdvancedError {
  if supplied.length() != 0 && supplied.length() != length {
    raise TradeSeqAdvancedError(
      prefix + " names do not match the matrix dimension",
    )
  }
  let output : Array[String] = []
  for index in 0..<length {
    let name = if supplied.length() == 0 {
      prefix + (index + 1).to_string()
    } else {
      supplied[index]
    }
    if name.trim().length() == 0 {
      raise TradeSeqAdvancedError(prefix + " names cannot be empty")
    }
    if output.contains(name) {
      raise TradeSeqAdvancedError(prefix + " names must be unique: " + name)
    }
    output.push(name)
  }
  output
}

///|
fn tradeseq_adv_validate_counts(
  counts : Array[Array[Double]],
) -> (Int, Int) raise TradeSeqAdvancedError {
  if counts.length() == 0 {
    raise TradeSeqAdvancedError("counts must contain at least one gene")
  }
  let cells = counts[0].length()
  if cells == 0 {
    raise TradeSeqAdvancedError("counts must contain at least one cell")
  }
  for row in counts {
    if row.length() != cells {
      raise TradeSeqAdvancedError("counts must be rectangular")
    }
    for value in row {
      if !tradeseq_adv_is_finite(value) || value < 0.0 {
        raise TradeSeqAdvancedError(
          "counts must contain finite non-negative values",
        )
      }
      if (value - value.round()).abs() > 1.0e-8 {
        raise TradeSeqAdvancedError("counts must contain integer values")
      }
    }
  }
  (counts.length(), cells)
}

///|
fn tradeseq_adv_prepare_trajectory(
  pseudotime : Array[Array[Double]],
  cell_weights : Array[Array[Double]],
  cells : Int,
) -> (Array[Array[Double]], Array[Array[Double]], Array[Double], Array[Double]) raise TradeSeqAdvancedError {
  if pseudotime.length() != cells || cell_weights.length() != cells {
    raise TradeSeqAdvancedError(
      "pseudotime and cell_weights rows must match count cells",
    )
  }
  if cells == 0 || pseudotime[0].length() == 0 {
    raise TradeSeqAdvancedError("at least one lineage is required")
  }
  let lineages = pseudotime[0].length()
  if cell_weights[0].length() != lineages {
    raise TradeSeqAdvancedError(
      "pseudotime and cell_weights must have identical dimensions",
    )
  }
  let prepared_time = tradeseq_adv_copy_matrix(pseudotime)
  let prepared_weights = tradeseq_adv_zero_matrix(cells, lineages)
  for cell in 0..<cells {
    if pseudotime[cell].length() != lineages ||
      cell_weights[cell].length() != lineages {
      raise TradeSeqAdvancedError(
        "pseudotime and cell_weights must be rectangular",
      )
    }
    let mut total = 0.0
    for lineage in 0..<lineages {
      let weight = cell_weights[cell][lineage]
      if !tradeseq_adv_is_finite(weight) || weight < 0.0 {
        raise TradeSeqAdvancedError(
          "cell_weights must contain finite non-negative values",
        )
      }
      if weight > 0.0 && !tradeseq_adv_is_finite(pseudotime[cell][lineage]) {
        raise TradeSeqAdvancedError("active lineage pseudotime must be finite")
      }
      total = total + weight
    }
    if total <= 0.0 {
      raise TradeSeqAdvancedError(
        "every cell must have at least one positive lineage weight",
      )
    }
    for lineage in 0..<lineages {
      prepared_weights[cell][lineage] = cell_weights[cell][lineage] / total
    }
  }
  let minimum = Array::make(lineages, 1.0e300)
  let maximum = Array::make(lineages, -1.0e300)
  for cell in 0..<cells {
    for lineage in 0..<lineages {
      if prepared_weights[cell][lineage] > 0.0 {
        let value = prepared_time[cell][lineage]
        if value < minimum[lineage] {
          minimum[lineage] = value
        }
        if value > maximum[lineage] {
          maximum[lineage] = value
        }
      }
    }
  }
  for lineage in 0..<lineages {
    if minimum[lineage] > 1.0e299 {
      raise TradeSeqAdvancedError("every lineage must contain active cells")
    }
    if maximum[lineage] - minimum[lineage] <= 1.0e-12 {
      raise TradeSeqAdvancedError(
        "every lineage must span more than one pseudotime value",
      )
    }
  }
  (prepared_time, prepared_weights, minimum, maximum)
}

///|
fn tradeseq_adv_prepare_offsets(
  counts : Array[Array[Double]],
  supplied : Array[Double],
) -> Array[Double] raise TradeSeqAdvancedError {
  let cells = counts[0].length()
  if supplied.length() != 0 {
    if supplied.length() != cells {
      raise TradeSeqAdvancedError("offsets must match count cells")
    }
    for value in supplied {
      if !tradeseq_adv_is_finite(value) {
        raise TradeSeqAdvancedError("offsets must contain finite values")
      }
    }
    return supplied.copy()
  }
  let offsets = Array::make(cells, 0.0)
  let mut mean_log_library = 0.0
  for cell in 0..<cells {
    let mut library = 0.0
    for gene in 0..<counts.length() {
      library = library + counts[gene][cell]
    }
    offsets[cell] = @math.ln(library.max(1.0))
    mean_log_library = mean_log_library + offsets[cell]
  }
  mean_log_library = mean_log_library / cells.to_double()
  for cell in 0..<cells {
    offsets[cell] = offsets[cell] - mean_log_library
  }
  offsets
}

///|
fn tradeseq_adv_bspline_basis(x_value : Double, count : Int) -> Array[Double] {
  let x = tradeseq_adv_clamp(x_value, 0.0, 1.0)
  let degree = (count - 1).min(3)
  if x >= 1.0 {
    let output = Array::make(count, 0.0)
    output[count - 1] = 1.0
    return output
  }
  let knot_count = count + degree + 1
  let knots = Array::make(knot_count, 0.0)
  for index in 0..<knot_count {
    if index <= degree {
      knots[index] = 0.0
    } else if index >= count {
      knots[index] = 1.0
    } else {
      knots[index] = (index - degree).to_double() / (count - degree).to_double()
    }
  }
  let mut previous = Array::make(knot_count - 1, 0.0)
  for index in 0..<(knot_count - 1) {
    if x >= knots[index] && x < knots[index + 1] {
      previous[index] = 1.0
    }
  }
  let mut level = 1
  while level <= degree {
    let current = Array::make(knot_count - level - 1, 0.0)
    for index in 0..<current.length() {
      let left_denominator = knots[index + level] - knots[index]
      let right_denominator = knots[index + level + 1] - knots[index + 1]
      let left = if left_denominator > 0.0 {
        (x - knots[index]) / left_denominator * previous[index]
      } else {
        0.0
      }
      let right = if right_denominator > 0.0 {
        (knots[index + level + 1] - x) / right_denominator * previous[index + 1]
      } else {
        0.0
      }
      current[index] = left + right
    }
    previous = current
    level = level + 1
  }
  let output : Array[Double] = []
  for index in 0..<count {
    output.push(previous[index])
  }
  output
}

///|
fn tradeseq_adv_normalized_time(
  value : Double,
  minimum : Double,
  maximum : Double,
) -> Double {
  tradeseq_adv_clamp((value - minimum) / (maximum - minimum), 0.0, 1.0)
}

///|
fn tradeseq_adv_observations(
  pseudotime : Array[Array[Double]],
  weights : Array[Array[Double]],
  minimum : Array[Double],
  maximum : Array[Double],
  n_knots : Int,
) -> Array[TradeSeqAdvObservation] {
  let output : Array[TradeSeqAdvObservation] = []
  for cell in 0..<pseudotime.length() {
    for lineage in 0..<pseudotime[cell].length() {
      if weights[cell][lineage] > 1.0e-12 {
        let normalized = tradeseq_adv_normalized_time(
          pseudotime[cell][lineage],
          minimum[lineage],
          maximum[lineage],
        )
        output.push(TradeSeqAdvObservation::{
          cell,
          lineage,
          weight: weights[cell][lineage],
          basis: tradeseq_adv_bspline_basis(normalized, n_knots),
        })
      }
    }
  }
  output
}

///|
fn tradeseq_adv_add_penalty(
  matrix : Array[Array[Double]],
  lineages : Int,
  n_knots : Int,
  penalty : Double,
  ridge : Double,
) -> Unit {
  let parameters = lineages * n_knots
  for parameter in 0..<parameters {
    matrix[parameter][parameter] = matrix[parameter][parameter] + ridge
  }
  if penalty <= 0.0 || n_knots < 3 {
    return
  }
  for lineage in 0..<lineages {
    let start = lineage * n_knots
    for position in 0..<(n_knots - 2) {
      let indices = [
        start + position,
        start + position + 1,
        start + position + 2,
      ]
      let values = [1.0, -2.0, 1.0]
      for left in 0..<3 {
        for right in 0..<3 {
          matrix[indices[left]][indices[right]] = matrix[indices[left]][indices[right]] +
            penalty * values[left] * values[right]
        }
      }
    }
  }
}

///|
fn tradeseq_adv_inverse(matrix : Array[Array[Double]]) -> Array[Array[Double]]? {
  if matrix.length() == 0 || matrix[0].length() != matrix.length() {
    return None
  }
  let size = matrix.length()
  let augmented = tradeseq_adv_zero_matrix(size, size * 2)
  for row in 0..<size {
    if matrix[row].length() != size {
      return None
    }
    for column in 0..<size {
      augmented[row][column] = matrix[row][column]
    }
    augmented[row][size + row] = 1.0
  }
  for column in 0..<size {
    let mut pivot = column
    let mut pivot_value = augmented[pivot][column].abs()
    for row in (column + 1)..<size {
      let value = augmented[row][column].abs()
      if value > pivot_value {
        pivot = row
        pivot_value = value
      }
    }
    if pivot_value < 1.0e-14 || !tradeseq_adv_is_finite(pivot_value) {
      return None
    }
    if pivot != column {
      let temporary = augmented[column]
      augmented[column] = augmented[pivot]
      augmented[pivot] = temporary
    }
    let divisor = augmented[column][column]
    for current in 0..<(size * 2) {
      augmented[column][current] = augmented[column][current] / divisor
    }
    for row in 0..<size {
      if row != column {
        let factor = augmented[row][column]
        if factor.abs() > 0.0 {
          for current in 0..<(size * 2) {
            augmented[row][current] = augmented[row][current] -
              factor * augmented[column][current]
          }
        }
      }
    }
  }
  let output = tradeseq_adv_zero_matrix(size, size)
  for row in 0..<size {
    for column in 0..<size {
      output[row][column] = augmented[row][size + column]
    }
  }
  Some(output)
}

///|
fn tradeseq_adv_matrix_vector(
  matrix : Array[Array[Double]],
  vector : Array[Double],
) -> Array[Double] {
  let output = Array::make(matrix.length(), 0.0)
  for row in 0..<matrix.length() {
    let mut value = 0.0
    for column in 0..<vector.length() {
      value = value + matrix[row][column] * vector[column]
    }
    output[row] = value
  }
  output
}

///|
fn tradeseq_adv_log_gamma(value : Double) -> Double {
  let coefficients = [
    0.99999999999980993, 676.5203681218851, -1259.1392167224028, 771.32342877765313,
    -176.61502916214059, 12.507343278686905, -0.13857109526572012, 9.9843695780195716e-6,
    1.5056327351493116e-7,
  ]
  let shifted = value - 1.0
  let mut sum = coefficients[0]
  for index in 1..<coefficients.length() {
    sum = sum + coefficients[index] / (shifted + index.to_double())
  }
  let t = shifted + 7.5
  0.9189385332046727 + (shifted + 0.5) * @math.ln(t) - t + @math.ln(sum)
}

///|
fn tradeseq_adv_nb_log_likelihood(
  response : Array[Double],
  fitted : Array[Double],
  observations : Array[TradeSeqAdvObservation],
  dispersion : Double,
) -> Double {
  let size = 1.0 / dispersion
  let mut result = 0.0
  for index in 0..<observations.length() {
    let y = response[observations[index].cell]
    let mu = fitted[index].max(1.0e-12)
    let contribution = tradeseq_adv_log_gamma(y + size) -
      tradeseq_adv_log_gamma(size) -
      tradeseq_adv_log_gamma(y + 1.0) +
      size * @math.ln(size / (size + mu)) +
      y * @math.ln(mu / (size + mu))
    result = result + observations[index].weight * contribution
  }
  result
}

///|
fn tradeseq_adv_initial_coefficients(
  response : Array[Double],
  offsets : Array[Double],
  observations : Array[TradeSeqAdvObservation],
  lineages : Int,
  n_knots : Int,
) -> Array[Double] {
  let weighted_sum = Array::make(lineages, 0.0)
  let weight_sum = Array::make(lineages, 0.0)
  for observation in observations {
    weighted_sum[observation.lineage] = weighted_sum[observation.lineage] +
      observation.weight *
      response[observation.cell] *
      @math.exp(-offsets[observation.cell])
    weight_sum[observation.lineage] = weight_sum[observation.lineage] +
      observation.weight
  }
  let coefficients = Array::make(lineages * n_knots, 0.0)
  for lineage in 0..<lineages {
    let mean = if weight_sum[lineage] > 0.0 {
      weighted_sum[lineage] / weight_sum[lineage]
    } else {
      0.0
    }
    let value = @math.ln(mean.max(1.0e-8))
    for knot in 0..<n_knots {
      coefficients[lineage * n_knots + knot] = value
    }
  }
  coefficients
}

///|
fn tradeseq_adv_eta(
  coefficients : Array[Double],
  observation : TradeSeqAdvObservation,
  n_knots : Int,
  offset : Double,
) -> Double {
  let start = observation.lineage * n_knots
  let mut value = offset
  for knot in 0..<n_knots {
    value = value + observation.basis[knot] * coefficients[start + knot]
  }
  tradeseq_adv_clamp(value, -25.0, 25.0)
}

///|
fn tradeseq_adv_fit_gene(
  gene_id : String,
  response : Array[Double],
  offsets : Array[Double],
  observations : Array[TradeSeqAdvObservation],
  lineages : Int,
  config : TradeSeqAdvancedConfig,
) -> TradeSeqGeneModel {
  let n_knots = config.n_knots
  let parameters = lineages * n_knots
  let coefficients = tradeseq_adv_initial_coefficients(
    response, offsets, observations, lineages, n_knots,
  )
  let mut dispersion = 0.1
  let mut converged = false
  let mut iterations = 0
  let expanded_fitted = Array::make(observations.length(), 0.0)
  let mut final_information = tradeseq_adv_zero_matrix(parameters, parameters)
  let mut iteration = 0
  while iteration < config.max_iterations {
    let information = tradeseq_adv_zero_matrix(parameters, parameters)
    let right_hand_side = Array::make(parameters, 0.0)
    let mut dispersion_numerator = 0.0
    let mut dispersion_denominator = 0.0
    for index in 0..<observations.length() {
      let observation = observations[index]
      let eta = tradeseq_adv_eta(
        coefficients,
        observation,
        n_knots,
        offsets[observation.cell],
      )
      let mu = @math.exp(eta).max(1.0e-10)
      expanded_fitted[index] = mu
      let y = response[observation.cell]
      let residual = y - mu
      dispersion_numerator = dispersion_numerator +
        observation.weight * (residual * residual - mu)
      dispersion_denominator = dispersion_denominator +
        observation.weight * mu * mu
    }
    if iteration > 0 && dispersion_denominator > 0.0 {
      let estimate = tradeseq_adv_clamp(
        dispersion_numerator / dispersion_denominator,
        config.minimum_dispersion,
        config.maximum_dispersion,
      )
      dispersion = 0.7 * dispersion + 0.3 * estimate
    }
    for index in 0..<observations.length() {
      let observation = observations[index]
      let mu = expanded_fitted[index]
      let y = response[observation.cell]
      let eta = tradeseq_adv_eta(
        coefficients,
        observation,
        n_knots,
        offsets[observation.cell],
      )
      let irls_weight = observation.weight * mu / (1.0 + dispersion * mu)
      let working = eta + (y - mu) / mu - offsets[observation.cell]
      let start = observation.lineage * n_knots
      for left in 0..<n_knots {
        let left_index = start + left
        right_hand_side[left_index] = right_hand_side[left_index] +
          irls_weight * observation.basis[left] * working
        for right in 0..<n_knots {
          let right_index = start + right
          information[left_index][right_index] = information[left_index][right_index] +
            irls_weight * observation.basis[left] * observation.basis[right]
        }
      }
    }
    tradeseq_adv_add_penalty(
      information,
      lineages,
      n_knots,
      config.smoothing_penalty,
      config.ridge,
    )
    let inverse = match tradeseq_adv_inverse(information) {
      Some(value) => value
      None => {
        iterations = iteration + 1
        break
      }
    }
    let updated = tradeseq_adv_matrix_vector(inverse, right_hand_side)
    let mut maximum_change = 0.0
    for parameter in 0..<parameters {
      let change = (updated[parameter] - coefficients[parameter]).abs()
      if change > maximum_change {
        maximum_change = change
      }
      coefficients[parameter] = updated[parameter]
    }
    final_information = information
    iterations = iteration + 1
    let mut coefficient_scale = 0.0
    for value in coefficients {
      coefficient_scale = coefficient_scale.max(value.abs())
    }
    if maximum_change <= config.tolerance * (1.0 + coefficient_scale) {
      converged = true
      break
    }
    iteration = iteration + 1
  }
  for index in 0..<observations.length() {
    let observation = observations[index]
    expanded_fitted[index] = @math.exp(
      tradeseq_adv_eta(
        coefficients,
        observation,
        n_knots,
        offsets[observation.cell],
      ),
    )
  }
  let fitted_values = Array::make(response.length(), 0.0)
  for index in 0..<observations.length() {
    let observation = observations[index]
    fitted_values[observation.cell] = fitted_values[observation.cell] +
      observation.weight * expanded_fitted[index]
  }
  if final_information.length() == 0 {
    final_information = tradeseq_adv_zero_matrix(parameters, parameters)
    for parameter in 0..<parameters {
      final_information[parameter][parameter] = 1.0 / config.ridge
    }
  }
  let covariance = match tradeseq_adv_inverse(final_information) {
    Some(value) => value
    None => {
      let fallback = tradeseq_adv_zero_matrix(parameters, parameters)
      for parameter in 0..<parameters {
        fallback[parameter][parameter] = 1.0 / config.ridge
      }
      fallback
    }
  }
  let log_likelihood = tradeseq_adv_nb_log_likelihood(
    response, expanded_fitted, observations, dispersion,
  )
  TradeSeqGeneModel::{
    gene_id,
    coefficients,
    covariance,
    fitted_values,
    dispersion,
    log_likelihood,
    aic: -2.0 * log_likelihood + 2.0 * parameters.to_double(),
    effective_df: parameters.to_double(),
    converged,
    iterations,
  }
}

///|
pub fn tradeseq_fit_advanced(
  counts : Array[Array[Double]],
  pseudotime : Array[Array[Double]],
  cell_weights : Array[Array[Double]],
  gene_names? : Array[String] = [],
  lineage_names? : Array[String] = [],
  offsets? : Array[Double] = [],
  config? : TradeSeqAdvancedConfig = TradeSeqAdvancedConfig::default(),
) -> TradeSeqAdvancedFit raise TradeSeqAdvancedError {
  let (genes, cells) = tradeseq_adv_validate_counts(counts)
  let (prepared_time, prepared_weights, minimum, maximum) = tradeseq_adv_prepare_trajectory(
    pseudotime, cell_weights, cells,
  )
  let lineages = prepared_time[0].length()
  let prepared_gene_names = tradeseq_adv_prepare_names(
    gene_names, genes, "Gene",
  )
  let prepared_lineage_names = tradeseq_adv_prepare_names(
    lineage_names, lineages, "Lineage",
  )
  let prepared_offsets = tradeseq_adv_prepare_offsets(counts, offsets)
  let observations = tradeseq_adv_observations(
    prepared_time,
    prepared_weights,
    minimum,
    maximum,
    config.n_knots,
  )
  let models : Array[TradeSeqGeneModel] = []
  for gene in 0..<genes {
    models.push(
      tradeseq_adv_fit_gene(
        prepared_gene_names[gene],
        counts[gene],
        prepared_offsets,
        observations,
        lineages,
        config,
      ),
    )
  }
  TradeSeqAdvancedFit::{
    gene_names: prepared_gene_names,
    lineage_names: prepared_lineage_names,
    models,
    pseudotime: prepared_time,
    cell_weights: prepared_weights,
    offsets: prepared_offsets,
    lineage_minimum: minimum,
    lineage_maximum: maximum,
    config,
    n_genes: genes,
    n_cells: cells,
    n_lineages: lineages,
  }
}

///|
fn tradeseq_adv_contrast_curve(
  lineages : Int,
  n_knots : Int,
  lineage : Int,
  first : Double,
  second : Double,
) -> Array[Double] {
  let output = Array::make(lineages * n_knots, 0.0)
  let first_basis = tradeseq_adv_bspline_basis(first, n_knots)
  let second_basis = tradeseq_adv_bspline_basis(second, n_knots)
  for knot in 0..<n_knots {
    output[lineage * n_knots + knot] = first_basis[knot] - second_basis[knot]
  }
  output
}

///|
fn tradeseq_adv_contrast_lineages(
  lineages : Int,
  n_knots : Int,
  first_lineage : Int,
  second_lineage : Int,
  position : Double,
) -> Array[Double] {
  let output = Array::make(lineages * n_knots, 0.0)
  let basis = tradeseq_adv_bspline_basis(position, n_knots)
  for knot in 0..<n_knots {
    output[first_lineage * n_knots + knot] = basis[knot]
    output[second_lineage * n_knots + knot] = -basis[knot]
  }
  output
}

///|
fn tradeseq_adv_independent_rows(
  rows : Array[Array[Double]],
) -> Array[Array[Double]] {
  let selected : Array[Array[Double]] = []
  let orthonormal : Array[Array[Double]] = []
  for row in rows {
    let residual = row.copy()
    for direction in orthonormal {
      let mut projection = 0.0
      for column in 0..<residual.length() {
        projection = projection + residual[column] * direction[column]
      }
      for column in 0..<residual.length() {
        residual[column] = residual[column] - projection * direction[column]
      }
    }
    let mut norm_squared = 0.0
    for value in residual {
      norm_squared = norm_squared + value * value
    }
    if norm_squared > 1.0e-12 {
      let norm = norm_squared.sqrt()
      for column in 0..<residual.length() {
        residual[column] = residual[column] / norm
      }
      selected.push(row.copy())
      orthonormal.push(residual)
    }
  }
  selected
}

///|
fn tradeseq_adv_log_gamma_q(a : Double, x : Double) -> Double {
  if x <= 0.0 {
    return 1.0
  }
  if x < a + 1.0 {
    let mut term = 1.0 / a
    let mut sum = term
    let mut index = 1
    while index < 300 {
      term = term * x / (a + index.to_double())
      sum = sum + term
      if term.abs() < sum.abs() * 1.0e-12 {
        break
      }
      index = index + 1
    }
    return tradeseq_adv_clamp(
      1.0 - sum * @math.exp(-x + a * @math.ln(x) - tradeseq_adv_log_gamma(a)),
      0.0,
      1.0,
    )
  }
  let tiny = 1.0e-300
  let mut b = x + 1.0 - a
  let mut c = 1.0 / tiny
  let mut d = 1.0 / b
  let mut value = d
  let mut index = 1
  while index < 300 {
    let an = -index.to_double() * (index.to_double() - a)
    b = b + 2.0
    d = an * d + b
    if d.abs() < tiny {
      d = tiny
    }
    c = b + an / c
    if c.abs() < tiny {
      c = tiny
    }
    d = 1.0 / d
    let delta = d * c
    value = value * delta
    if (delta - 1.0).abs() < 1.0e-12 {
      break
    }
    index = index + 1
  }
  tradeseq_adv_clamp(
    value * @math.exp(-x + a * @math.ln(x) - tradeseq_adv_log_gamma(a)),
    0.0,
    1.0,
  )
}

///|
fn tradeseq_adv_chi_square_q(degrees : Int, statistic : Double) -> Double {
  if statistic <= 0.0 || degrees <= 0 {
    return 1.0
  }
  tradeseq_adv_log_gamma_q(0.5 * degrees.to_double(), 0.5 * statistic)
}

///|
fn tradeseq_adv_wald(
  model : TradeSeqGeneModel,
  raw_contrasts : Array[Array[Double]],
  l2fc : Double,
) -> (Double, Int, Double, Double) {
  let contrasts = tradeseq_adv_independent_rows(raw_contrasts)
  let degrees = contrasts.length()
  if degrees == 0 {
    return (0.0, 0, 1.0, 0.0)
  }
  let effects = Array::make(degrees, 0.0)
  let adjusted_effects = Array::make(degrees, 0.0)
  let threshold = l2fc.max(0.0) * @math.ln(2.0)
  let mut maximum_effect = 0.0
  for contrast in 0..<degrees {
    let mut effect = 0.0
    for parameter in 0..<model.coefficients.length() {
      effect = effect +
        contrasts[contrast][parameter] * model.coefficients[parameter]
    }
    effects[contrast] = effect
    adjusted_effects[contrast] = if effect > threshold {
      effect - threshold
    } else if effect < -threshold {
      effect + threshold
    } else {
      0.0
    }
    if effect.abs() > maximum_effect.abs() {
      maximum_effect = effect
    }
  }
  let variance = tradeseq_adv_zero_matrix(degrees, degrees)
  for left in 0..<degrees {
    for right in 0..<degrees {
      let mut value = 0.0
      for first in 0..<model.coefficients.length() {
        for second in 0..<model.coefficients.length() {
          value = value +
            contrasts[left][first] *
            model.covariance[first][second] *
            contrasts[right][second]
        }
      }
      variance[left][right] = value
    }
  }
  let mut maximum_diagonal = 0.0
  for index in 0..<degrees {
    maximum_diagonal = maximum_diagonal.max(variance[index][index].abs())
  }
  let regularization = maximum_diagonal * 1.0e-8 + 1.0e-10
  for index in 0..<degrees {
    variance[index][index] = variance[index][index] + regularization
  }
  let inverse = match tradeseq_adv_inverse(variance) {
    Some(value) => value
    None => return (0.0, degrees, 1.0, maximum_effect / @math.ln(2.0))
  }
  let projected = tradeseq_adv_matrix_vector(inverse, adjusted_effects)
  let mut statistic = 0.0
  for index in 0..<degrees {
    statistic = statistic + adjusted_effects[index] * projected[index]
  }
  statistic = statistic.max(0.0)
  (
    statistic,
    degrees,
    tradeseq_adv_chi_square_q(degrees, statistic),
    maximum_effect / @math.ln(2.0),
  )
}

///|
fn tradeseq_adv_bh(values : Array[Double]) -> Array[Double] {
  let indices : Array[Int] = []
  for index in 0..<values.length() {
    indices.push(index)
  }
  for left in 0..<indices.length() {
    let mut best = left
    for right in (left + 1)..<indices.length() {
      if values[indices[right]] < values[indices[best]] {
        best = right
      }
    }
    if best != left {
      let temporary = indices[left]
      indices[left] = indices[best]
      indices[best] = temporary
    }
  }
  let output = Array::make(values.length(), 1.0)
  let mut running = 1.0
  let mut position = indices.length() - 1
  while position >= 0 {
    let index = indices[position]
    let adjusted = (values[index] *
    indices.length().to_double() /
    (position + 1).to_double()).min(1.0)
    running = running.min(adjusted)
    output[index] = running
    position = position - 1
  }
  output
}

///|
fn tradeseq_adv_test_table(
  fit : TradeSeqAdvancedFit,
  test_name : String,
  contrasts : Array[Array[Double]],
  l2fc : Double,
) -> TradeSeqAdvancedTestTable {
  let statistics = Array::make(fit.n_genes, 0.0)
  let degrees = Array::make(fit.n_genes, 0)
  let p_values = Array::make(fit.n_genes, 1.0)
  let fold_changes = Array::make(fit.n_genes, 0.0)
  for gene in 0..<fit.n_genes {
    let (statistic, df, p_value, fold_change) = tradeseq_adv_wald(
      fit.models[gene],
      contrasts,
      l2fc,
    )
    statistics[gene] = statistic
    degrees[gene] = df
    p_values[gene] = p_value
    fold_changes[gene] = fold_change
  }
  let adjusted = tradeseq_adv_bh(p_values)
  let results : Array[TradeSeqAdvancedTestResult] = []
  for gene in 0..<fit.n_genes {
    results.push(TradeSeqAdvancedTestResult::{
      gene_id: fit.gene_names[gene],
      wald_statistic: statistics[gene],
      degrees_freedom: degrees[gene],
      p_value: p_values[gene],
      adjusted_p_value: adjusted[gene],
      log2_fold_change: fold_changes[gene],
      significant: adjusted[gene] <= fit.config.fdr_threshold,
    })
  }
  TradeSeqAdvancedTestTable::{
    test_name,
    results,
    fdr_threshold: fit.config.fdr_threshold,
  }
}

///|
pub fn tradeseq_association_test_advanced(
  fit : TradeSeqAdvancedFit,
  l2fc? : Double = 0.0,
) -> TradeSeqAdvancedTestTable raise TradeSeqAdvancedError {
  if !tradeseq_adv_is_finite(l2fc) || l2fc < 0.0 {
    raise TradeSeqAdvancedError("l2fc must be finite and non-negative")
  }
  let contrasts : Array[Array[Double]] = []
  for lineage in 0..<fit.n_lineages {
    for point in 1..<fit.config.test_points {
      contrasts.push(
        tradeseq_adv_contrast_curve(
          fit.n_lineages,
          fit.config.n_knots,
          lineage,
          point.to_double() / (fit.config.test_points - 1).to_double(),
          0.0,
        ),
      )
    }
  }
  tradeseq_adv_test_table(fit, "associationTest", contrasts, l2fc)
}

///|
pub fn tradeseq_start_vs_end_test_advanced(
  fit : TradeSeqAdvancedFit,
  l2fc? : Double = 0.0,
) -> TradeSeqAdvancedTestTable raise TradeSeqAdvancedError {
  if !tradeseq_adv_is_finite(l2fc) || l2fc < 0.0 {
    raise TradeSeqAdvancedError("l2fc must be finite and non-negative")
  }
  let contrasts : Array[Array[Double]] = []
  for lineage in 0..<fit.n_lineages {
    contrasts.push(
      tradeseq_adv_contrast_curve(
        fit.n_lineages,
        fit.config.n_knots,
        lineage,
        1.0,
        0.0,
      ),
    )
  }
  tradeseq_adv_test_table(fit, "startVsEndTest", contrasts, l2fc)
}

///|
pub fn tradeseq_diff_end_test_advanced(
  fit : TradeSeqAdvancedFit,
  l2fc? : Double = 0.0,
) -> TradeSeqAdvancedTestTable raise TradeSeqAdvancedError {
  if fit.n_lineages < 2 {
    raise TradeSeqAdvancedError("diffEndTest requires at least two lineages")
  }
  if !tradeseq_adv_is_finite(l2fc) || l2fc < 0.0 {
    raise TradeSeqAdvancedError("l2fc must be finite and non-negative")
  }
  let contrasts : Array[Array[Double]] = []
  for first in 0..<fit.n_lineages {
    for second in (first + 1)..<fit.n_lineages {
      contrasts.push(
        tradeseq_adv_contrast_lineages(
          fit.n_lineages,
          fit.config.n_knots,
          first,
          second,
          1.0,
        ),
      )
    }
  }
  tradeseq_adv_test_table(fit, "diffEndTest", contrasts, l2fc)
}

///|
pub fn tradeseq_pattern_test_advanced(
  fit : TradeSeqAdvancedFit,
  n_points? : Int = 0,
  l2fc? : Double = 0.0,
) -> TradeSeqAdvancedTestTable raise TradeSeqAdvancedError {
  if fit.n_lineages < 2 {
    raise TradeSeqAdvancedError("patternTest requires at least two lineages")
  }
  let points = if n_points == 0 { fit.config.test_points } else { n_points }
  if points < 2 {
    raise TradeSeqAdvancedError("n_points must be at least 2")
  }
  if !tradeseq_adv_is_finite(l2fc) || l2fc < 0.0 {
    raise TradeSeqAdvancedError("l2fc must be finite and non-negative")
  }
  let contrasts : Array[Array[Double]] = []
  for first in 0..<fit.n_lineages {
    for second in (first + 1)..<fit.n_lineages {
      for point in 0..<points {
        contrasts.push(
          tradeseq_adv_contrast_lineages(
            fit.n_lineages,
            fit.config.n_knots,
            first,
            second,
            point.to_double() / (points - 1).to_double(),
          ),
        )
      }
    }
  }
  tradeseq_adv_test_table(fit, "patternTest", contrasts, l2fc)
}

///|
pub fn tradeseq_early_de_test_advanced(
  fit : TradeSeqAdvancedFit,
  start_fraction : Double,
  end_fraction : Double,
  n_points? : Int = 0,
  l2fc? : Double = 0.0,
) -> TradeSeqAdvancedTestTable raise TradeSeqAdvancedError {
  if fit.n_lineages < 2 {
    raise TradeSeqAdvancedError("earlyDETest requires at least two lineages")
  }
  if !tradeseq_adv_is_finite(start_fraction) ||
    !tradeseq_adv_is_finite(end_fraction) ||
    start_fraction < 0.0 ||
    end_fraction > 1.0 ||
    end_fraction <= start_fraction {
    raise TradeSeqAdvancedError(
      "earlyDE interval must satisfy 0 <= start < end <= 1",
    )
  }
  let points = if n_points == 0 { fit.config.test_points } else { n_points }
  if points < 2 {
    raise TradeSeqAdvancedError("n_points must be at least 2")
  }
  if !tradeseq_adv_is_finite(l2fc) || l2fc < 0.0 {
    raise TradeSeqAdvancedError("l2fc must be finite and non-negative")
  }
  let contrasts : Array[Array[Double]] = []
  for first in 0..<fit.n_lineages {
    for second in (first + 1)..<fit.n_lineages {
      for point in 0..<points {
        let fraction = start_fraction +
          (end_fraction - start_fraction) *
          point.to_double() /
          (points - 1).to_double()
        contrasts.push(
          tradeseq_adv_contrast_lineages(
            fit.n_lineages,
            fit.config.n_knots,
            first,
            second,
            fraction,
          ),
        )
      }
    }
  }
  tradeseq_adv_test_table(fit, "earlyDETest", contrasts, l2fc)
}

///|
pub fn tradeseq_predict_smooth_advanced(
  fit : TradeSeqAdvancedFit,
  gene_id : String,
  n_points? : Int = 100,
) -> TradeSeqSmoothPrediction raise TradeSeqAdvancedError {
  if n_points < 2 {
    raise TradeSeqAdvancedError("n_points must be at least 2")
  }
  let mut gene = -1
  for index in 0..<fit.gene_names.length() {
    if fit.gene_names[index] == gene_id {
      gene = index
      break
    }
  }
  if gene < 0 {
    raise TradeSeqAdvancedError("unknown gene: " + gene_id)
  }
  let model = fit.models[gene]
  let pseudotime = tradeseq_adv_zero_matrix(fit.n_lineages, n_points)
  let fitted = tradeseq_adv_zero_matrix(fit.n_lineages, n_points)
  let standard_errors = tradeseq_adv_zero_matrix(fit.n_lineages, n_points)
  let parameters = fit.n_lineages * fit.config.n_knots
  for lineage in 0..<fit.n_lineages {
    for point in 0..<n_points {
      let fraction = point.to_double() / (n_points - 1).to_double()
      pseudotime[lineage][point] = fit.lineage_minimum[lineage] +
        fraction * (fit.lineage_maximum[lineage] - fit.lineage_minimum[lineage])
      let basis = tradeseq_adv_bspline_basis(fraction, fit.config.n_knots)
      let design = Array::make(parameters, 0.0)
      let mut eta = 0.0
      for knot in 0..<fit.config.n_knots {
        let parameter = lineage * fit.config.n_knots + knot
        design[parameter] = basis[knot]
        eta = eta + basis[knot] * model.coefficients[parameter]
      }
      eta = tradeseq_adv_clamp(eta, -25.0, 25.0)
      let mean = @math.exp(eta)
      let mut variance = 0.0
      for first in 0..<parameters {
        for second in 0..<parameters {
          variance = variance +
            design[first] * model.covariance[first][second] * design[second]
        }
      }
      fitted[lineage][point] = mean
      standard_errors[lineage][point] = mean * variance.max(0.0).sqrt()
    }
  }
  TradeSeqSmoothPrediction::{ gene_id, pseudotime, fitted, standard_errors, }
}

///|
pub fn tradeseq_evaluate_k_advanced(
  counts : Array[Array[Double]],
  pseudotime : Array[Array[Double]],
  cell_weights : Array[Array[Double]],
  candidates : Array[Int],
  gene_names? : Array[String] = [],
  lineage_names? : Array[String] = [],
  offsets? : Array[Double] = [],
  config? : TradeSeqAdvancedConfig = TradeSeqAdvancedConfig::default(),
) -> TradeSeqKnotEvaluation raise TradeSeqAdvancedError {
  if candidates.length() < 2 {
    raise TradeSeqAdvancedError("at least two knot candidates are required")
  }
  let unique : Array[Int] = []
  for candidate in candidates {
    if candidate < 3 {
      raise TradeSeqAdvancedError("knot candidates must be at least 3")
    }
    if unique.contains(candidate) {
      raise TradeSeqAdvancedError("knot candidates must be unique")
    }
    unique.push(candidate)
  }
  let genes = counts.length()
  let gene_aic = tradeseq_adv_zero_matrix(genes, candidates.length())
  let mean_aic = Array::make(candidates.length(), 0.0)
  for candidate_index in 0..<candidates.length() {
    let candidate_config = TradeSeqAdvancedConfig::create(
      n_knots=candidates[candidate_index],
      smoothing_penalty=config.smoothing_penalty,
      max_iterations=config.max_iterations,
      tolerance=config.tolerance,
      minimum_dispersion=config.minimum_dispersion,
      maximum_dispersion=config.maximum_dispersion,
      ridge=config.ridge,
      fdr_threshold=config.fdr_threshold,
      test_points=config.test_points,
    )
    let fit = tradeseq_fit_advanced(
      counts,
      pseudotime,
      cell_weights,
      gene_names~,
      lineage_names~,
      offsets~,
      config=candidate_config,
    )
    for gene in 0..<genes {
      gene_aic[gene][candidate_index] = fit.models[gene].aic
      mean_aic[candidate_index] = mean_aic[candidate_index] +
        fit.models[gene].aic
    }
    mean_aic[candidate_index] = mean_aic[candidate_index] / genes.to_double()
  }
  let mut best = 0
  for index in 1..<candidates.length() {
    if mean_aic[index] < mean_aic[best] {
      best = index
    }
  }
  TradeSeqKnotEvaluation::{
    candidates: candidates.copy(),
    gene_aic,
    mean_aic,
    selected_n_knots: candidates[best],
  }
}

///|
pub fn tradeseq_fit_from_slingshot(
  counts : Array[Array[Double]],
  gene_names : Array[String],
  slingshot : SlingshotAdvancedResult,
  offsets? : Array[Double] = [],
  config? : TradeSeqAdvancedConfig = TradeSeqAdvancedConfig::default(),
) -> TradeSeqAdvancedFit raise TradeSeqAdvancedError {
  let cells = slingshot.weights.length()
  if cells == 0 {
    raise TradeSeqAdvancedError("slingshot result contains no cells")
  }
  let lineages = slingshot.weights[0].length()
  let pseudotime = tradeseq_adv_zero_matrix(cells, lineages)
  let weights = tradeseq_adv_zero_matrix(cells, lineages)
  for cell in 0..<cells {
    if slingshot.weights[cell].length() != lineages ||
      slingshot.pseudotime[cell].length() != lineages {
      raise TradeSeqAdvancedError("slingshot matrices must be rectangular")
    }
    for lineage in 0..<lineages {
      match slingshot.pseudotime[cell][lineage] {
        Some(value) => {
          pseudotime[cell][lineage] = value
          weights[cell][lineage] = slingshot.weights[cell][lineage]
        }
        None => {
          pseudotime[cell][lineage] = 0.0
          weights[cell][lineage] = 0.0
        }
      }
    }
  }
  let lineage_names : Array[String] = []
  for lineage in 0..<lineages {
    lineage_names.push(slingshot.curves[lineage].name)
  }
  tradeseq_fit_advanced(
    counts,
    pseudotime,
    weights,
    gene_names~,
    lineage_names~,
    offsets~,
    config~,
  )
}

///|
fn tradeseq_adv_copy_string_map(
  source : Map[String, Array[String]],
) -> Map[String, Array[String]] {
  let output : Map[String, Array[String]] = Map([])
  for key in source.keys() {
    output[key] = source[key].copy()
  }
  output
}

///|
fn tradeseq_adv_copy_sce(source : SingleCellExperiment) -> SingleCellExperiment {
  let assays : Map[String, Array[Array[Double]]] = Map([])
  for key in source.assays.keys() {
    assays[key] = tradeseq_adv_copy_matrix(source.assays[key])
  }
  let reduced_dims : Map[String, Array[Array[Double]]] = Map([])
  for key in source.reduced_dims.keys() {
    reduced_dims[key] = tradeseq_adv_copy_matrix(source.reduced_dims[key])
  }
  let alternatives : Map[String, SingleCellExperiment] = Map([])
  for key in source.alternative_experiments.keys() {
    alternatives[key] = tradeseq_adv_copy_sce(
      source.alternative_experiments[key],
    )
  }
  let metadata : Map[String, String] = Map([])
  for key in source.metadata.keys() {
    metadata[key] = source.metadata[key]
  }
  SingleCellExperiment::{
    assays,
    row_data: tradeseq_adv_copy_string_map(source.row_data),
    col_data: tradeseq_adv_copy_string_map(source.col_data),
    reduced_dims,
    alternative_experiments: alternatives,
    metadata,
    row_names: source.row_names.copy(),
    col_names: source.col_names.copy(),
  }
}

///|
fn tradeseq_adv_test_column(
  table : TradeSeqAdvancedTestTable,
  adjusted : Bool,
) -> Array[String] {
  table.results.map(fn(result) {
    if adjusted {
      result.adjusted_p_value.to_string()
    } else {
      result.p_value.to_string()
    }
  })
}

///|
pub fn tradeseq_advanced_sce(
  experiment : SingleCellExperiment,
  assay_name? : String = "counts",
  pseudotime_name? : String = "slingshot.pseudotime",
  weights_name? : String = "slingshot.weights",
  output_prefix? : String = "tradeSeq",
  config? : TradeSeqAdvancedConfig = TradeSeqAdvancedConfig::default(),
) -> TradeSeqAdvancedSceOutput raise TradeSeqAdvancedError {
  if output_prefix.trim().length() == 0 {
    raise TradeSeqAdvancedError("output_prefix cannot be empty")
  }
  let counts = match experiment.assays.get(assay_name) {
    Some(value) => value
    None => raise TradeSeqAdvancedError("SCE assay not found: " + assay_name)
  }
  let pseudotime = match experiment.reduced_dims.get(pseudotime_name) {
    Some(value) => value
    None =>
      raise TradeSeqAdvancedError(
        "SCE pseudotime matrix not found: " + pseudotime_name,
      )
  }
  let weights = match experiment.reduced_dims.get(weights_name) {
    Some(value) => value
    None =>
      raise TradeSeqAdvancedError(
        "SCE lineage weight matrix not found: " + weights_name,
      )
  }
  let fit = tradeseq_fit_advanced(
    counts,
    pseudotime,
    weights,
    gene_names=experiment.row_names,
    config~,
  )
  let association = tradeseq_association_test_advanced(fit)
  let start_end = tradeseq_start_vs_end_test_advanced(fit)
  let output = tradeseq_adv_copy_sce(experiment)
  let fitted : Array[Array[Double]] = []
  for model in fit.models {
    fitted.push(model.fitted_values.copy())
  }
  output.assays[output_prefix + ".fitted"] = fitted
  output.row_data[output_prefix + ".association.pvalue"] = tradeseq_adv_test_column(
    association, false,
  )
  output.row_data[output_prefix + ".association.padj"] = tradeseq_adv_test_column(
    association, true,
  )
  output.row_data[output_prefix + ".startVsEnd.pvalue"] = tradeseq_adv_test_column(
    start_end, false,
  )
  output.row_data[output_prefix + ".startVsEnd.padj"] = tradeseq_adv_test_column(
    start_end, true,
  )
  output.row_data[output_prefix + ".dispersion"] = fit.models.map(fn(model) {
    model.dispersion.to_string()
  })
  output.row_data[output_prefix + ".converged"] = fit.models.map(fn(model) {
    model.converged.to_string()
  })
  output.metadata[output_prefix + ".nKnots"] = fit.config.n_knots.to_string()
  output.metadata[output_prefix + ".lineages"] = fit.lineage_names.join(",")
  TradeSeqAdvancedSceOutput::{ experiment: output, fit, }
}

///|
pub fn tradeseq_advanced_summary(fit : TradeSeqAdvancedFit) -> String {
  let mut converged = 0
  let mut mean_dispersion = 0.0
  for model in fit.models {
    if model.converged {
      converged = converged + 1
    }
    mean_dispersion = mean_dispersion + model.dispersion
  }
  mean_dispersion = mean_dispersion / fit.n_genes.to_double()
  "tradeSeq 1.27.0 advanced fit\n" +
  "genes: " +
  fit.n_genes.to_string() +
  "\n" +
  "cells: " +
  fit.n_cells.to_string() +
  "\n" +
  "lineages: " +
  fit.n_lineages.to_string() +
  "\n" +
  "knots: " +
  fit.config.n_knots.to_string() +
  "\n" +
  "converged genes: " +
  converged.to_string() +
  "/" +
  fit.n_genes.to_string() +
  "\n" +
  "mean dispersion: " +
  mean_dispersion.to_string()
}
