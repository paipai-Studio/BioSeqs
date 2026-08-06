// Advanced single-cell hurdle models inspired by Bioconductor MAST 1.39.0.
//
// Expression matrices use the Bioconductor feature x cell convention. Each
// feature is fitted with a Bayesian logistic model for detection and a
// Gaussian model for the positive expression values.

///|
pub suberror MastAdvancedError {
  MastAdvancedError(String)
}

///|
pub struct MastAdvancedConfig {
  max_iterations : Int
  tolerance : Double
  cauchy_prior_scale : Double
  ridge : Double
  empirical_bayes : Bool
  ebayes_use_full_model : Bool
  maximum_prior_df : Double
  minimum_detected : Int
  minimum_probability : Double
  fdr_threshold : Double
} derive(Debug)

///|
pub fn MastAdvancedConfig::create(
  max_iterations? : Int = 100,
  tolerance? : Double = 1.0e-8,
  cauchy_prior_scale? : Double = 2.5,
  ridge? : Double = 1.0e-8,
  empirical_bayes? : Bool = true,
  ebayes_use_full_model? : Bool = false,
  maximum_prior_df? : Double = 1000.0,
  minimum_detected? : Int = 2,
  minimum_probability? : Double = 1.0e-9,
  fdr_threshold? : Double = 0.05,
) -> MastAdvancedConfig raise MastAdvancedError {
  if max_iterations < 1 {
    raise MastAdvancedError("max_iterations must be positive")
  }
  if tolerance <= 0.0 || !mast_adv_is_finite(tolerance) {
    raise MastAdvancedError("tolerance must be finite and positive")
  }
  if cauchy_prior_scale <= 0.0 || !mast_adv_is_finite(cauchy_prior_scale) {
    raise MastAdvancedError("cauchy_prior_scale must be finite and positive")
  }
  if ridge < 0.0 || !mast_adv_is_finite(ridge) {
    raise MastAdvancedError("ridge must be finite and non-negative")
  }
  if maximum_prior_df <= 0.0 || !mast_adv_is_finite(maximum_prior_df) {
    raise MastAdvancedError("maximum_prior_df must be finite and positive")
  }
  if minimum_detected < 1 {
    raise MastAdvancedError("minimum_detected must be positive")
  }
  if minimum_probability <= 0.0 || minimum_probability >= 0.5 {
    raise MastAdvancedError("minimum_probability must be in (0, 0.5)")
  }
  if fdr_threshold <= 0.0 || fdr_threshold > 1.0 {
    raise MastAdvancedError("fdr_threshold must be in (0, 1]")
  }
  MastAdvancedConfig::{
    max_iterations,
    tolerance,
    cauchy_prior_scale,
    ridge,
    empirical_bayes,
    ebayes_use_full_model,
    maximum_prior_df,
    minimum_detected,
    minimum_probability,
    fdr_threshold,
  }
}

///|
pub fn MastAdvancedConfig::default() -> MastAdvancedConfig {
  MastAdvancedConfig::{
    max_iterations: 100,
    tolerance: 1.0e-8,
    cauchy_prior_scale: 2.5,
    ridge: 1.0e-8,
    empirical_bayes: true,
    ebayes_use_full_model: false,
    maximum_prior_df: 1000.0,
    minimum_detected: 2,
    minimum_probability: 1.0e-9,
    fdr_threshold: 0.05,
  }
}

///|
pub struct MastAdvancedData {
  expression : Array[Array[Double]]
  gene_names : Array[String]
  cell_names : Array[String]
  design : Array[Array[Double]]
  coefficient_names : Array[String]
  cdr : Array[Double]
  detection_limit : Double
  n_genes : Int
  n_cells : Int
  n_coefficients : Int
} derive(Debug)

///|
pub struct MastAdvancedResult {
  gene_names : Array[String]
  coefficient_names : Array[String]
  tested_coefficients : Array[Int]
  tested_names : Array[String]
  detection_rates : Array[Double]
  n_detected : Array[Int]
  discrete_coefficients : Array[Array[Double]]
  continuous_coefficients : Array[Array[Double]]
  discrete_standard_errors : Array[Array[Double]]
  continuous_standard_errors : Array[Array[Double]]
  discrete_statistics : Array[Double?]
  continuous_statistics : Array[Double?]
  hurdle_statistics : Array[Double?]
  discrete_p_values : Array[Double?]
  continuous_p_values : Array[Double?]
  hurdle_p_values : Array[Double?]
  discrete_fdr : Array[Double?]
  continuous_fdr : Array[Double?]
  hurdle_fdr : Array[Double?]
  marginal_log_fc : Array[Double?]
  marginal_log_fc_variance : Array[Double?]
  continuous_variance : Array[Double?]
  moderated_variance : Array[Double?]
  discrete_converged : Array[Bool]
  continuous_converged : Array[Bool]
  prior_variance : Double
  prior_df : Double
  empirical_bayes : Bool
  fdr_threshold : Double
} derive(Debug)

///|
pub struct MastAdvancedSCEOutput {
  experiment : SingleCellExperiment
  result : MastAdvancedResult
}

///|
priv struct MastAdvCoreFit {
  coefficients : Array[Double]
  covariance : Array[Array[Double]]
  log_likelihood : Double
  residual_df : Int
  raw_variance : Double
  moderated_variance : Double
  converged : Bool
  fitted : Bool
}

///|
fn mast_adv_is_finite(value : Double) -> Bool {
  value == value && value.abs() < 1.0e300
}

///|
fn mast_adv_clamp(value : Double, lower : Double, upper : Double) -> Double {
  if value < lower {
    lower
  } else if value > upper {
    upper
  } else {
    value
  }
}

///|
fn mast_adv_zero_matrix(rows : Int, columns : Int) -> Array[Array[Double]] {
  let output : Array[Array[Double]] = []
  for _ in 0..<rows {
    output.push(Array::make(columns, 0.0))
  }
  output
}

///|
fn mast_adv_copy_matrix(matrix : Array[Array[Double]]) -> Array[Array[Double]] {
  let output : Array[Array[Double]] = []
  for row in matrix {
    output.push(row.copy())
  }
  output
}

///|
fn mast_adv_copy_string_map(
  source : Map[String, Array[String]],
) -> Map[String, Array[String]] {
  let output : Map[String, Array[String]] = Map([])
  let keys = source.keys().collect()
  for key in keys {
    output[key] = source[key].copy()
  }
  output
}

///|
fn mast_adv_prepare_names(
  supplied : Array[String],
  length : Int,
  prefix : String,
) -> Array[String] raise MastAdvancedError {
  if supplied.length() != 0 && supplied.length() != length {
    raise MastAdvancedError(prefix + " names do not match the matrix dimension")
  }
  let output : Array[String] = []
  for index in 0..<length {
    let name = if supplied.length() == 0 {
      prefix + (index + 1).to_string()
    } else {
      supplied[index]
    }
    if name.trim().length() == 0 {
      raise MastAdvancedError(prefix + " names cannot be empty")
    }
    if output.contains(name) {
      raise MastAdvancedError(prefix + " names must be unique: " + name)
    }
    output.push(name)
  }
  output
}

///|
fn mast_adv_validate_expression(
  expression : Array[Array[Double]],
) -> (Int, Int) raise MastAdvancedError {
  if expression.length() == 0 {
    raise MastAdvancedError("expression must contain at least one gene")
  }
  let cells = expression[0].length()
  if cells == 0 {
    raise MastAdvancedError("expression must contain at least one cell")
  }
  for row in expression {
    if row.length() != cells {
      raise MastAdvancedError("expression must be rectangular")
    }
    for value in row {
      if !mast_adv_is_finite(value) || value < 0.0 {
        raise MastAdvancedError(
          "expression must contain finite non-negative values",
        )
      }
    }
  }
  (expression.length(), cells)
}

///|
fn mast_adv_matrix_rank(matrix : Array[Array[Double]]) -> Int {
  if matrix.length() == 0 || matrix[0].length() == 0 {
    return 0
  }
  let work = mast_adv_copy_matrix(matrix)
  let rows = work.length()
  let columns = work[0].length()
  let mut rank = 0
  let mut column = 0
  while column < columns && rank < rows {
    let mut pivot = rank
    let mut pivot_value = work[pivot][column].abs()
    for candidate in (rank + 1)..<rows {
      let value = work[candidate][column].abs()
      if value > pivot_value {
        pivot = candidate
        pivot_value = value
      }
    }
    if pivot_value > 1.0e-10 {
      if pivot != rank {
        let temporary = work[rank]
        work[rank] = work[pivot]
        work[pivot] = temporary
      }
      let divisor = work[rank][column]
      for current in column..<columns {
        work[rank][current] = work[rank][current] / divisor
      }
      for row in 0..<rows {
        if row != rank {
          let factor = work[row][column]
          for current in column..<columns {
            work[row][current] = work[row][current] -
              factor * work[rank][current]
          }
        }
      }
      rank = rank + 1
    }
    column = column + 1
  }
  rank
}

///|
fn mast_adv_calculate_cdr(
  expression : Array[Array[Double]],
  detection_limit : Double,
) -> Array[Double] {
  let genes = expression.length()
  let cells = expression[0].length()
  let output = Array::make(cells, 0.0)
  for cell in 0..<cells {
    let mut detected = 0
    for gene in 0..<genes {
      if expression[gene][cell] > detection_limit {
        detected = detected + 1
      }
    }
    output[cell] = detected.to_double() / genes.to_double()
  }
  output
}

///|
pub fn MastAdvancedData::create(
  expression : Array[Array[Double]],
  design : Array[Array[Double]],
  coefficient_names : Array[String],
  gene_names? : Array[String] = [],
  cell_names? : Array[String] = [],
  detection_limit? : Double = 0.0,
  include_cdr? : Bool = true,
  cdr_name? : String = "cngeneson",
) -> MastAdvancedData raise MastAdvancedError {
  let (genes, cells) = mast_adv_validate_expression(expression)
  if !mast_adv_is_finite(detection_limit) || detection_limit < 0.0 {
    raise MastAdvancedError("detection_limit must be finite and non-negative")
  }
  if design.length() != cells {
    raise MastAdvancedError("design rows must match expression cells")
  }
  if design.length() == 0 || design[0].length() == 0 {
    raise MastAdvancedError("design must contain at least one coefficient")
  }
  let supplied_coefficients = design[0].length()
  if coefficient_names.length() != supplied_coefficients {
    raise MastAdvancedError("coefficient_names must match the design columns")
  }
  let prepared_coefficient_names = mast_adv_prepare_names(
    coefficient_names, supplied_coefficients, "coefficient",
  )
  for row in design {
    if row.length() != supplied_coefficients {
      raise MastAdvancedError("design must be rectangular")
    }
    for value in row {
      if !mast_adv_is_finite(value) {
        raise MastAdvancedError("design must contain finite values")
      }
    }
  }
  let cdr = mast_adv_calculate_cdr(expression, detection_limit)
  let prepared_design : Array[Array[Double]] = []
  for cell in 0..<cells {
    let row = design[cell].copy()
    if include_cdr {
      row.push(cdr[cell])
    }
    prepared_design.push(row)
  }
  let final_coefficient_names = prepared_coefficient_names.copy()
  if include_cdr {
    if cdr_name.trim().length() == 0 {
      raise MastAdvancedError("cdr_name cannot be empty")
    }
    if final_coefficient_names.contains(cdr_name) {
      raise MastAdvancedError("cdr_name duplicates a coefficient name")
    }
    final_coefficient_names.push(cdr_name)
  }
  if mast_adv_matrix_rank(prepared_design) != final_coefficient_names.length() {
    raise MastAdvancedError("design matrix must have full column rank")
  }
  MastAdvancedData::{
    expression: mast_adv_copy_matrix(expression),
    gene_names: mast_adv_prepare_names(gene_names, genes, "Gene"),
    cell_names: mast_adv_prepare_names(cell_names, cells, "Cell"),
    design: prepared_design,
    coefficient_names: final_coefficient_names,
    cdr,
    detection_limit,
    n_genes: genes,
    n_cells: cells,
    n_coefficients: final_coefficient_names.length(),
  }
}

///|
pub fn mast_advanced_from_groups(
  expression : Array[Array[Double]],
  groups : Array[String],
  reference? : String = "",
  gene_names? : Array[String] = [],
  cell_names? : Array[String] = [],
  detection_limit? : Double = 0.0,
  include_cdr? : Bool = true,
) -> MastAdvancedData raise MastAdvancedError {
  let (_, cells) = mast_adv_validate_expression(expression)
  if groups.length() != cells {
    raise MastAdvancedError("groups must match expression cells")
  }
  let levels : Array[String] = []
  for group in groups {
    if group.trim().length() == 0 {
      raise MastAdvancedError("group labels cannot be empty")
    }
    if !levels.contains(group) {
      levels.push(group)
    }
  }
  if levels.length() < 2 {
    raise MastAdvancedError("at least two group levels are required")
  }
  let selected_reference = if reference == "" { levels[0] } else { reference }
  if !levels.contains(selected_reference) {
    raise MastAdvancedError("reference group is absent")
  }
  let encoded_levels : Array[String] = []
  for level in levels {
    if level != selected_reference {
      encoded_levels.push(level)
    }
  }
  let design : Array[Array[Double]] = []
  for group in groups {
    let row = Array::make(encoded_levels.length() + 1, 0.0)
    row[0] = 1.0
    for level in 0..<encoded_levels.length() {
      if group == encoded_levels[level] {
        row[level + 1] = 1.0
      }
    }
    design.push(row)
  }
  let coefficient_names : Array[String] = ["(Intercept)"]
  for level in encoded_levels {
    coefficient_names.push("group:" + level)
  }
  MastAdvancedData::create(
    expression,
    design,
    coefficient_names,
    gene_names~,
    cell_names~,
    detection_limit~,
    include_cdr~,
  )
}

///|
fn mast_adv_cholesky(matrix : Array[Array[Double]]) -> Array[Array[Double]]? {
  let size = matrix.length()
  if size == 0 {
    return None
  }
  let lower = mast_adv_zero_matrix(size, size)
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
        if value <= 1.0e-14 || !mast_adv_is_finite(value) {
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
fn mast_adv_cholesky_solve(
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
  let output = Array::make(size, 0.0)
  let mut row = size - 1
  while row >= 0 {
    let mut value = forward[row]
    for column in (row + 1)..<size {
      value = value - lower[column][row] * output[column]
    }
    output[row] = value / lower[row][row]
    row = row - 1
  }
  output
}

///|
fn mast_adv_cholesky_inverse(
  lower : Array[Array[Double]],
) -> Array[Array[Double]] {
  let size = lower.length()
  let inverse = mast_adv_zero_matrix(size, size)
  for column in 0..<size {
    let basis = Array::make(size, 0.0)
    basis[column] = 1.0
    let solution = mast_adv_cholesky_solve(lower, basis)
    for row in 0..<size {
      inverse[row][column] = solution[row]
    }
  }
  inverse
}

///|
fn mast_adv_regularized_cholesky(
  matrix : Array[Array[Double]],
  ridge : Double,
) -> Array[Array[Double]]? {
  let mut added = ridge
  for attempt in 0..<14 {
    let candidate = mast_adv_copy_matrix(matrix)
    for index in 0..<candidate.length() {
      candidate[index][index] = candidate[index][index] + added
    }
    match mast_adv_cholesky(candidate) {
      Some(lower) => return Some(lower)
      None =>
        added = if attempt == 0 && added == 0.0 {
          1.0e-10
        } else if added == 0.0 {
          1.0e-10
        } else {
          added * 10.0
        }
    }
  }
  None
}

///|
fn mast_adv_logistic(value : Double) -> Double {
  if value >= 0.0 {
    let inverse = @math.exp(-value)
    1.0 / (1.0 + inverse)
  } else {
    let exponential = @math.exp(value)
    exponential / (1.0 + exponential)
  }
}

///|
fn mast_adv_linear_predictor(
  design : Array[Array[Double]],
  coefficients : Array[Double],
) -> Array[Double] {
  let output = Array::make(design.length(), 0.0)
  for row in 0..<design.length() {
    let mut value = 0.0
    for column in 0..<coefficients.length() {
      value = value + design[row][column] * coefficients[column]
    }
    output[row] = value
  }
  output
}

///|
fn mast_adv_logistic_fit(
  design : Array[Array[Double]],
  detected : Array[Bool],
  config : MastAdvancedConfig,
) -> MastAdvCoreFit {
  let observations = design.length()
  let parameters = design[0].length()
  let coefficients = Array::make(parameters, 0.0)
  let mut positives = 0
  for value in detected {
    if value {
      positives = positives + 1
    }
  }
  if positives == 0 || positives == observations {
    return MastAdvCoreFit::{
      coefficients,
      covariance: mast_adv_zero_matrix(parameters, parameters),
      log_likelihood: 0.0,
      residual_df: 0,
      raw_variance: 0.0,
      moderated_variance: 0.0,
      converged: true,
      fitted: false,
    }
  }
  coefficients[0] = @math.ln(
    positives.to_double() / (observations - positives).to_double(),
  )
  let scale_squared = config.cauchy_prior_scale * config.cauchy_prior_scale
  let mut converged = false
  for _ in 0..<config.max_iterations {
    let eta = mast_adv_linear_predictor(design, coefficients)
    let information = mast_adv_zero_matrix(parameters, parameters)
    let right_hand_side = Array::make(parameters, 0.0)
    for observation in 0..<observations {
      let probability = mast_adv_clamp(
        mast_adv_logistic(eta[observation]),
        config.minimum_probability,
        1.0 - config.minimum_probability,
      )
      let weight = probability * (1.0 - probability)
      let response = if detected[observation] { 1.0 } else { 0.0 }
      let working = eta[observation] + (response - probability) / weight
      for left in 0..<parameters {
        right_hand_side[left] = right_hand_side[left] +
          design[observation][left] * weight * working
        for right in 0..<parameters {
          information[left][right] = information[left][right] +
            design[observation][left] * weight * design[observation][right]
        }
      }
    }
    for parameter in 1..<parameters {
      information[parameter][parameter] = information[parameter][parameter] +
        2.0 /
        (scale_squared + coefficients[parameter] * coefficients[parameter])
    }
    let lower = match mast_adv_regularized_cholesky(information, config.ridge) {
      Some(value) => value
      None => break
    }
    let updated = mast_adv_cholesky_solve(lower, right_hand_side)
    let mut maximum_change = 0.0
    for parameter in 0..<parameters {
      maximum_change = maximum_change.max(
        (updated[parameter] - coefficients[parameter]).abs(),
      )
      coefficients[parameter] = mast_adv_clamp(updated[parameter], -40.0, 40.0)
    }
    if maximum_change < config.tolerance {
      converged = true
      break
    }
  }
  let eta = mast_adv_linear_predictor(design, coefficients)
  let final_information = mast_adv_zero_matrix(parameters, parameters)
  let mut log_likelihood = 0.0
  for observation in 0..<observations {
    let probability = mast_adv_clamp(
      mast_adv_logistic(eta[observation]),
      config.minimum_probability,
      1.0 - config.minimum_probability,
    )
    let weight = probability * (1.0 - probability)
    let contribution = if detected[observation] {
      @math.ln(probability)
    } else {
      @math.ln(1.0 - probability)
    }
    log_likelihood = log_likelihood + contribution
    for left in 0..<parameters {
      for right in 0..<parameters {
        final_information[left][right] = final_information[left][right] +
          design[observation][left] * weight * design[observation][right]
      }
    }
  }
  for parameter in 1..<parameters {
    final_information[parameter][parameter] = final_information[parameter][parameter] +
      2.0 / (scale_squared + coefficients[parameter] * coefficients[parameter])
  }
  let covariance = match
    mast_adv_regularized_cholesky(final_information, config.ridge) {
    Some(lower) => mast_adv_cholesky_inverse(lower)
    None => mast_adv_zero_matrix(parameters, parameters)
  }
  MastAdvCoreFit::{
    coefficients,
    covariance,
    log_likelihood,
    residual_df: (positives.min(observations - positives) - parameters).max(0),
    raw_variance: 1.0,
    moderated_variance: 1.0,
    converged,
    fitted: converged,
  }
}

///|
fn mast_adv_positive_data(
  expression : Array[Double],
  design : Array[Array[Double]],
  detection_limit : Double,
) -> (Array[Double], Array[Array[Double]]) {
  let response : Array[Double] = []
  let positive_design : Array[Array[Double]] = []
  for cell in 0..<expression.length() {
    if expression[cell] > detection_limit {
      response.push(expression[cell])
      positive_design.push(design[cell].copy())
    }
  }
  (response, positive_design)
}

///|
fn mast_adv_continuous_fit(
  expression : Array[Double],
  design : Array[Array[Double]],
  detection_limit : Double,
  prior_variance : Double,
  prior_df : Double,
  config : MastAdvancedConfig,
) -> MastAdvCoreFit {
  let parameters = design[0].length()
  let (response, positive_design) = mast_adv_positive_data(
    expression, design, detection_limit,
  )
  let observations = response.length()
  if observations <= parameters ||
    observations < config.minimum_detected ||
    mast_adv_matrix_rank(positive_design) != parameters {
    return MastAdvCoreFit::{
      coefficients: Array::make(parameters, 0.0),
      covariance: mast_adv_zero_matrix(parameters, parameters),
      log_likelihood: 0.0,
      residual_df: 0,
      raw_variance: 0.0,
      moderated_variance: 0.0,
      converged: false,
      fitted: false,
    }
  }
  let information = mast_adv_zero_matrix(parameters, parameters)
  let right_hand_side = Array::make(parameters, 0.0)
  for observation in 0..<observations {
    for left in 0..<parameters {
      right_hand_side[left] = right_hand_side[left] +
        positive_design[observation][left] * response[observation]
      for right in 0..<parameters {
        information[left][right] = information[left][right] +
          positive_design[observation][left] *
          positive_design[observation][right]
      }
    }
  }
  let lower = match mast_adv_regularized_cholesky(information, config.ridge) {
    Some(value) => value
    None =>
      return MastAdvCoreFit::{
        coefficients: Array::make(parameters, 0.0),
        covariance: mast_adv_zero_matrix(parameters, parameters),
        log_likelihood: 0.0,
        residual_df: 0,
        raw_variance: 0.0,
        moderated_variance: 0.0,
        converged: false,
        fitted: false,
      }
  }
  let coefficients = mast_adv_cholesky_solve(lower, right_hand_side)
  let fitted = mast_adv_linear_predictor(positive_design, coefficients)
  let mut residual_sum_squares = 0.0
  for observation in 0..<observations {
    let residual = response[observation] - fitted[observation]
    residual_sum_squares = residual_sum_squares + residual * residual
  }
  let residual_df = observations - parameters
  let raw_variance = residual_sum_squares / residual_df.to_double()
  let moderated_variance = if prior_df > 0.0 {
    (residual_sum_squares + prior_variance * prior_df) /
    (residual_df.to_double() + prior_df)
  } else {
    raw_variance
  }
  let likelihood_variance = (if prior_df > 0.0 {
    (residual_sum_squares + prior_variance * prior_df) /
    (observations.to_double() + prior_df)
  } else {
    residual_sum_squares / observations.to_double()
  }).max(1.0e-15)
  let covariance = mast_adv_cholesky_inverse(lower)
  for row in 0..<parameters {
    for column in 0..<parameters {
      covariance[row][column] = covariance[row][column] * moderated_variance
    }
  }
  let log_likelihood = -0.5 *
    observations.to_double() *
    (@math.ln(2.0 * 3.141592653589793 * likelihood_variance) + 1.0)
  MastAdvCoreFit::{
    coefficients,
    covariance,
    log_likelihood,
    residual_df,
    raw_variance,
    moderated_variance,
    converged: true,
    fitted: true,
  }
}

///|
fn mast_adv_log_gamma(value : Double) -> Double {
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
fn mast_adv_ebayes_log_likelihood(
  residual_df : Array[Double],
  residual_ss : Array[Double],
  log_a : Double,
  log_b : Double,
) -> Double {
  let a = @math.exp(log_a)
  let b = @math.exp(log_b)
  let mut value = 0.0
  for index in 0..<residual_df.length() {
    let half_df = 0.5 * residual_df[index]
    value = value -
      (
        mast_adv_log_gamma(half_df) +
        mast_adv_log_gamma(a) -
        mast_adv_log_gamma(half_df + a)
      ) -
      half_df * @math.ln(b) -
      (half_df + a) * @math.ln(1.0 + residual_ss[index] / (2.0 * b))
  }
  value
}

///|
fn mast_adv_ebayes_sufficient(
  data : MastAdvancedData,
  use_full_model : Bool,
) -> (Array[Double], Array[Double]) {
  let residual_df : Array[Double] = []
  let residual_ss : Array[Double] = []
  for gene in 0..<data.n_genes {
    let (response, positive_design) = mast_adv_positive_data(
      data.expression[gene],
      data.design,
      data.detection_limit,
    )
    let observations = response.length()
    if use_full_model {
      let parameters = data.n_coefficients
      if observations > parameters &&
        mast_adv_matrix_rank(positive_design) == parameters {
        let information = mast_adv_zero_matrix(parameters, parameters)
        let right_hand_side = Array::make(parameters, 0.0)
        for observation in 0..<observations {
          for left in 0..<parameters {
            right_hand_side[left] = right_hand_side[left] +
              positive_design[observation][left] * response[observation]
            for right in 0..<parameters {
              information[left][right] = information[left][right] +
                positive_design[observation][left] *
                positive_design[observation][right]
            }
          }
        }
        match mast_adv_regularized_cholesky(information, 1.0e-10) {
          Some(lower) => {
            let coefficients = mast_adv_cholesky_solve(lower, right_hand_side)
            let fitted = mast_adv_linear_predictor(
              positive_design, coefficients,
            )
            let mut sum_squares = 0.0
            for observation in 0..<observations {
              let residual = response[observation] - fitted[observation]
              sum_squares = sum_squares + residual * residual
            }
            if sum_squares > 1.0e-15 {
              residual_df.push((observations - parameters).to_double())
              residual_ss.push(sum_squares)
            }
          }
          None => ()
        }
      }
    } else if observations > 1 {
      let mut mean = 0.0
      for value in response {
        mean = mean + value
      }
      mean = mean / observations.to_double()
      let mut sum_squares = 0.0
      for value in response {
        let residual = value - mean
        sum_squares = sum_squares + residual * residual
      }
      if sum_squares > 1.0e-15 {
        residual_df.push((observations - 1).to_double())
        residual_ss.push(sum_squares)
      }
    }
  }
  (residual_df, residual_ss)
}

///|
fn mast_adv_estimate_ebayes(
  data : MastAdvancedData,
  config : MastAdvancedConfig,
) -> (Double, Double) {
  if !config.empirical_bayes {
    return (0.0, 0.0)
  }
  let (residual_df, residual_ss) = mast_adv_ebayes_sufficient(
    data,
    config.ebayes_use_full_model,
  )
  if residual_df.length() < 2 {
    return (0.0, 0.0)
  }
  let mut mean_variance = 0.0
  for index in 0..<residual_df.length() {
    mean_variance = mean_variance + residual_ss[index] / residual_df[index]
  }
  mean_variance = (mean_variance / residual_df.length().to_double()).max(
    1.0e-12,
  )
  let mut log_a = @math.ln(2.0)
  let mut log_b = @math.ln(2.0 * mean_variance)
  let mut step_a = 1.0
  let mut step_b = 1.0
  let mut best = mast_adv_ebayes_log_likelihood(
    residual_df, residual_ss, log_a, log_b,
  )
  for _ in 0..<160 {
    let mut improved = false
    let candidates = [
      (step_a, 0.0),
      (-step_a, 0.0),
      (0.0, step_b),
      (0.0, -step_b),
      (step_a, step_b),
      (step_a, -step_b),
      (-step_a, step_b),
      (-step_a, -step_b),
    ]
    for candidate in candidates {
      let candidate_a = mast_adv_clamp(log_a + candidate.0, -6.0, 12.0)
      let candidate_b = mast_adv_clamp(log_b + candidate.1, -24.0, 24.0)
      let likelihood = mast_adv_ebayes_log_likelihood(
        residual_df, residual_ss, candidate_a, candidate_b,
      )
      if likelihood > best {
        best = likelihood
        log_a = candidate_a
        log_b = candidate_b
        improved = true
      }
    }
    if !improved {
      step_a = step_a * 0.7
      step_b = step_b * 0.7
    }
    if step_a < 1.0e-5 && step_b < 1.0e-5 {
      break
    }
  }
  let a = @math.exp(log_a)
  let b = @math.exp(log_b)
  ((b / a).max(1.0e-12), (2.0 * a).min(config.maximum_prior_df))
}

///|
fn mast_adv_select_columns(
  design : Array[Array[Double]],
  keep : Array[Int],
) -> Array[Array[Double]] {
  let output : Array[Array[Double]] = []
  for row in design {
    let selected : Array[Double] = []
    for column in keep {
      selected.push(row[column])
    }
    output.push(selected)
  }
  output
}

///|
fn mast_adv_chi_square_q(degrees : Double, statistic : Double) -> Double {
  if statistic <= 0.0 {
    return 1.0
  }
  mast_adv_upper_gamma_q(0.5 * degrees, 0.5 * statistic)
}

///|
fn mast_adv_upper_gamma_q(a : Double, x : Double) -> Double {
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
    return mast_adv_clamp(
      1.0 - sum * @math.exp(-x + a * @math.ln(x) - mast_adv_log_gamma(a)),
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
  mast_adv_clamp(
    value * @math.exp(-x + a * @math.ln(x) - mast_adv_log_gamma(a)),
    0.0,
    1.0,
  )
}

///|
fn mast_adv_bh(values : Array[Double?]) -> Array[Double?] {
  let indices : Array[Int] = []
  for index in 0..<values.length() {
    if values[index] is Some(_) {
      indices.push(index)
    }
  }
  for left in 0..<indices.length() {
    let mut best = left
    for right in (left + 1)..<indices.length() {
      let right_value = values[indices[right]].unwrap_or(1.0)
      let best_value = values[indices[best]].unwrap_or(1.0)
      if right_value < best_value {
        best = right
      }
    }
    if best != left {
      let temporary = indices[left]
      indices[left] = indices[best]
      indices[best] = temporary
    }
  }
  let output : Array[Double?] = Array::make(values.length(), None)
  let count = indices.length()
  let mut running = 1.0
  let mut position = count - 1
  while position >= 0 {
    let index = indices[position]
    let adjusted = (values[index].unwrap_or(1.0) *
    count.to_double() /
    (position + 1).to_double()).min(1.0)
    running = running.min(adjusted)
    output[index] = Some(running)
    position = position - 1
  }
  output
}

///|
fn mast_adv_standard_errors(covariance : Array[Array[Double]]) -> Array[Double] {
  let output = Array::make(covariance.length(), 0.0)
  for index in 0..<covariance.length() {
    output[index] = covariance[index][index].max(0.0).sqrt()
  }
  output
}

///|
fn mast_adv_marginal_effect(
  discrete : MastAdvCoreFit,
  continuous : MastAdvCoreFit,
  tested : Int,
) -> (Double, Double)? {
  if !discrete.fitted ||
    !continuous.fitted ||
    tested <= 0 ||
    tested >= discrete.coefficients.length() {
    return None
  }
  let continuous_baseline = continuous.coefficients[0]
  let continuous_comparison = continuous_baseline +
    continuous.coefficients[tested]
  let discrete_baseline = mast_adv_logistic(discrete.coefficients[0])
  let discrete_comparison = mast_adv_logistic(
    discrete.coefficients[0] + discrete.coefficients[tested],
  )
  let effect = continuous_comparison * discrete_comparison -
    continuous_baseline * discrete_baseline
  let continuous_gradient = Array::make(continuous.coefficients.length(), 0.0)
  continuous_gradient[0] = discrete_comparison - discrete_baseline
  continuous_gradient[tested] = discrete_comparison
  let discrete_gradient = Array::make(discrete.coefficients.length(), 0.0)
  discrete_gradient[0] = continuous_comparison *
    discrete_comparison *
    (1.0 - discrete_comparison) -
    continuous_baseline * discrete_baseline * (1.0 - discrete_baseline)
  discrete_gradient[tested] = continuous_comparison *
    discrete_comparison *
    (1.0 - discrete_comparison)
  let mut variance = 0.0
  for left in 0..<continuous_gradient.length() {
    for right in 0..<continuous_gradient.length() {
      variance = variance +
        continuous_gradient[left] *
        continuous.covariance[left][right] *
        continuous_gradient[right] +
        discrete_gradient[left] *
        discrete.covariance[left][right] *
        discrete_gradient[right]
    }
  }
  Some((effect, variance.max(0.0)))
}

///|
pub fn mast_advanced_zlm(
  data : MastAdvancedData,
  tested_coefficients : Array[Int],
  config? : MastAdvancedConfig = MastAdvancedConfig::default(),
) -> MastAdvancedResult raise MastAdvancedError {
  if tested_coefficients.length() == 0 {
    raise MastAdvancedError("tested_coefficients cannot be empty")
  }
  let tested = tested_coefficients.copy()
  for index in 0..<tested.length() {
    if tested[index] <= 0 || tested[index] >= data.n_coefficients {
      raise MastAdvancedError(
        "tested coefficients must be non-intercept design columns",
      )
    }
    for previous in 0..<index {
      if tested[previous] == tested[index] {
        raise MastAdvancedError("tested coefficients must be unique")
      }
    }
  }
  let kept : Array[Int] = []
  for coefficient in 0..<data.n_coefficients {
    if !tested.contains(coefficient) {
      kept.push(coefficient)
    }
  }
  if kept.length() == 0 {
    raise MastAdvancedError("the reduced model must retain a coefficient")
  }
  let reduced_design = mast_adv_select_columns(data.design, kept)
  if mast_adv_matrix_rank(reduced_design) != kept.length() {
    raise MastAdvancedError("the reduced design matrix is rank deficient")
  }
  let (prior_variance, prior_df) = mast_adv_estimate_ebayes(data, config)
  let detection_rates = Array::make(data.n_genes, 0.0)
  let n_detected = Array::make(data.n_genes, 0)
  let discrete_coefficients : Array[Array[Double]] = []
  let continuous_coefficients : Array[Array[Double]] = []
  let discrete_standard_errors : Array[Array[Double]] = []
  let continuous_standard_errors : Array[Array[Double]] = []
  let discrete_statistics : Array[Double?] = Array::make(data.n_genes, None)
  let continuous_statistics : Array[Double?] = Array::make(data.n_genes, None)
  let hurdle_statistics : Array[Double?] = Array::make(data.n_genes, None)
  let discrete_p_values : Array[Double?] = Array::make(data.n_genes, None)
  let continuous_p_values : Array[Double?] = Array::make(data.n_genes, None)
  let hurdle_p_values : Array[Double?] = Array::make(data.n_genes, None)
  let marginal_log_fc : Array[Double?] = Array::make(data.n_genes, None)
  let marginal_log_fc_variance : Array[Double?] = Array::make(
    data.n_genes,
    None,
  )
  let continuous_variance : Array[Double?] = Array::make(data.n_genes, None)
  let moderated_variance : Array[Double?] = Array::make(data.n_genes, None)
  let discrete_converged = Array::make(data.n_genes, false)
  let continuous_converged = Array::make(data.n_genes, false)
  for gene in 0..<data.n_genes {
    let detected = Array::make(data.n_cells, false)
    let mut detected_count = 0
    for cell in 0..<data.n_cells {
      detected[cell] = data.expression[gene][cell] > data.detection_limit
      if detected[cell] {
        detected_count = detected_count + 1
      }
    }
    n_detected[gene] = detected_count
    detection_rates[gene] = detected_count.to_double() /
      data.n_cells.to_double()
    let discrete_full = mast_adv_logistic_fit(data.design, detected, config)
    let discrete_reduced = mast_adv_logistic_fit(
      reduced_design, detected, config,
    )
    let continuous_full = mast_adv_continuous_fit(
      data.expression[gene],
      data.design,
      data.detection_limit,
      prior_variance,
      prior_df,
      config,
    )
    let continuous_reduced = mast_adv_continuous_fit(
      data.expression[gene],
      reduced_design,
      data.detection_limit,
      prior_variance,
      prior_df,
      config,
    )
    discrete_coefficients.push(discrete_full.coefficients)
    continuous_coefficients.push(continuous_full.coefficients)
    discrete_standard_errors.push(
      mast_adv_standard_errors(discrete_full.covariance),
    )
    continuous_standard_errors.push(
      mast_adv_standard_errors(continuous_full.covariance),
    )
    discrete_converged[gene] = discrete_full.converged
    continuous_converged[gene] = continuous_full.converged
    if continuous_full.fitted {
      continuous_variance[gene] = Some(continuous_full.raw_variance)
      moderated_variance[gene] = Some(continuous_full.moderated_variance)
    }
    let discrete_testable = discrete_full.fitted && discrete_reduced.fitted
    let continuous_testable = continuous_full.fitted &&
      continuous_reduced.fitted &&
      continuous_full.residual_df > 0 &&
      continuous_reduced.residual_df > 0
    let mut combined_statistic = 0.0
    let mut combined_degrees = 0
    if discrete_testable {
      let statistic = (2.0 *
      (discrete_full.log_likelihood - discrete_reduced.log_likelihood)).max(0.0)
      discrete_statistics[gene] = Some(statistic)
      discrete_p_values[gene] = Some(
        mast_adv_chi_square_q(tested.length().to_double(), statistic),
      )
      combined_statistic = combined_statistic + statistic
      combined_degrees = combined_degrees + tested.length()
    }
    if continuous_testable {
      let statistic = (2.0 *
      (continuous_full.log_likelihood - continuous_reduced.log_likelihood)).max(
        0.0,
      )
      continuous_statistics[gene] = Some(statistic)
      continuous_p_values[gene] = Some(
        mast_adv_chi_square_q(tested.length().to_double(), statistic),
      )
      combined_statistic = combined_statistic + statistic
      combined_degrees = combined_degrees + tested.length()
    }
    if combined_degrees > 0 {
      hurdle_statistics[gene] = Some(combined_statistic)
      hurdle_p_values[gene] = Some(
        mast_adv_chi_square_q(combined_degrees.to_double(), combined_statistic),
      )
    }
    match mast_adv_marginal_effect(discrete_full, continuous_full, tested[0]) {
      Some(effect) => {
        marginal_log_fc[gene] = Some(effect.0)
        marginal_log_fc_variance[gene] = Some(effect.1)
      }
      None => ()
    }
  }
  let tested_names : Array[String] = []
  for coefficient in tested {
    tested_names.push(data.coefficient_names[coefficient])
  }
  MastAdvancedResult::{
    gene_names: data.gene_names.copy(),
    coefficient_names: data.coefficient_names.copy(),
    tested_coefficients: tested,
    tested_names,
    detection_rates,
    n_detected,
    discrete_coefficients,
    continuous_coefficients,
    discrete_standard_errors,
    continuous_standard_errors,
    discrete_statistics,
    continuous_statistics,
    hurdle_statistics,
    discrete_p_values,
    continuous_p_values,
    hurdle_p_values,
    discrete_fdr: mast_adv_bh(discrete_p_values),
    continuous_fdr: mast_adv_bh(continuous_p_values),
    hurdle_fdr: mast_adv_bh(hurdle_p_values),
    marginal_log_fc,
    marginal_log_fc_variance,
    continuous_variance,
    moderated_variance,
    discrete_converged,
    continuous_converged,
    prior_variance,
    prior_df,
    empirical_bayes: config.empirical_bayes && prior_df > 0.0,
    fdr_threshold: config.fdr_threshold,
  }
}

///|
pub fn MastAdvancedResult::index_of(
  self : MastAdvancedResult,
  gene_name : String,
) -> Int? {
  for index in 0..<self.gene_names.length() {
    if self.gene_names[index] == gene_name {
      return Some(index)
    }
  }
  None
}

///|
pub fn MastAdvancedResult::n_tested(self : MastAdvancedResult) -> Int {
  let mut count = 0
  for value in self.hurdle_p_values {
    if value is Some(_) {
      count = count + 1
    }
  }
  count
}

///|
pub fn MastAdvancedResult::n_significant(
  self : MastAdvancedResult,
  threshold? : Double = -1.0,
) -> Int {
  let selected = if threshold < 0.0 { self.fdr_threshold } else { threshold }
  let mut count = 0
  for value in self.hurdle_fdr {
    match value {
      Some(adjusted) => if adjusted <= selected { count = count + 1 }
      None => ()
    }
  }
  count
}

///|
pub fn MastAdvancedResult::significant_genes(
  self : MastAdvancedResult,
  threshold? : Double = -1.0,
) -> Array[String] {
  let selected = if threshold < 0.0 { self.fdr_threshold } else { threshold }
  let output : Array[String] = []
  for index in 0..<self.gene_names.length() {
    match self.hurdle_fdr[index] {
      Some(adjusted) =>
        if adjusted <= selected {
          output.push(self.gene_names[index])
        }
      None => ()
    }
  }
  output
}

///|
pub fn MastAdvancedResult::top_genes(
  self : MastAdvancedResult,
  number : Int,
) -> Array[String] {
  if number <= 0 {
    return []
  }
  let indices : Array[Int] = []
  for index in 0..<self.gene_names.length() {
    if self.hurdle_p_values[index] is Some(_) {
      indices.push(index)
    }
  }
  for left in 0..<indices.length() {
    let mut best = left
    for right in (left + 1)..<indices.length() {
      if self.hurdle_p_values[indices[right]].unwrap_or(1.0) <
        self.hurdle_p_values[indices[best]].unwrap_or(1.0) {
        best = right
      }
    }
    if best != left {
      let temporary = indices[left]
      indices[left] = indices[best]
      indices[best] = temporary
    }
  }
  let output : Array[String] = []
  let limit = number.min(indices.length())
  for position in 0..<limit {
    output.push(self.gene_names[indices[position]])
  }
  output
}

///|
pub fn MastAdvancedResult::summary(self : MastAdvancedResult) -> String {
  "MAST advanced zlm: " +
  self.gene_names.length().to_string() +
  " genes, " +
  self.n_tested().to_string() +
  " tested, " +
  self.n_significant().to_string() +
  " significant at FDR " +
  self.fdr_threshold.to_string() +
  "; contrast=" +
  self.tested_names.join(",") +
  "; prior_df=" +
  self.prior_df.to_string()
}

///|
fn mast_adv_option_string(value : Double?) -> String {
  match value {
    Some(number) => number.to_string()
    None => "NA"
  }
}

///|
fn mast_adv_copy_sce(experiment : SingleCellExperiment) -> SingleCellExperiment {
  let assays : Map[String, Array[Array[Double]]] = Map([])
  let assay_names = experiment.assays.keys().collect()
  for name in assay_names {
    assays[name] = mast_adv_copy_matrix(experiment.assays[name])
  }
  let reduced_dims : Map[String, Array[Array[Double]]] = Map([])
  let reduced_names = experiment.reduced_dims.keys().collect()
  for name in reduced_names {
    reduced_dims[name] = mast_adv_copy_matrix(experiment.reduced_dims[name])
  }
  let alternatives : Map[String, SingleCellExperiment] = Map([])
  let alternative_names = experiment.alternative_experiments.keys().collect()
  for name in alternative_names {
    alternatives[name] = mast_adv_copy_sce(
      experiment.alternative_experiments[name],
    )
  }
  let metadata : Map[String, String] = Map([])
  let metadata_names = experiment.metadata.keys().collect()
  for name in metadata_names {
    metadata[name] = experiment.metadata[name]
  }
  SingleCellExperiment::{
    assays,
    row_data: mast_adv_copy_string_map(experiment.row_data),
    col_data: mast_adv_copy_string_map(experiment.col_data),
    reduced_dims,
    alternative_experiments: alternatives,
    metadata,
    row_names: experiment.row_names.copy(),
    col_names: experiment.col_names.copy(),
  }
}

///|
pub fn mast_advanced_zlm_sce(
  experiment : SingleCellExperiment,
  group_column : String,
  assay_name? : String = "logcounts",
  reference? : String = "",
  output_prefix? : String = "mast",
  detection_limit? : Double = 0.0,
  include_cdr? : Bool = true,
  config? : MastAdvancedConfig = MastAdvancedConfig::default(),
) -> MastAdvancedSCEOutput raise MastAdvancedError {
  if output_prefix.trim().length() == 0 {
    raise MastAdvancedError("output_prefix cannot be empty")
  }
  let expression = match experiment.assays.get(assay_name) {
    Some(value) => value
    None => raise MastAdvancedError("SCE assay not found: " + assay_name)
  }
  let groups = match experiment.col_data.get(group_column) {
    Some(value) => value
    None =>
      raise MastAdvancedError("SCE column data not found: " + group_column)
  }
  let data = mast_advanced_from_groups(
    expression,
    groups,
    reference~,
    gene_names=experiment.row_names,
    cell_names=experiment.col_names,
    detection_limit~,
    include_cdr~,
  )
  let tested : Array[Int] = []
  let upper = if include_cdr {
    data.n_coefficients - 1
  } else {
    data.n_coefficients
  }
  for coefficient in 1..<upper {
    tested.push(coefficient)
  }
  let result = mast_advanced_zlm(data, tested, config~)
  let enriched = mast_adv_copy_sce(experiment)
  let detection_rate : Array[String] = []
  let detected : Array[String] = []
  let discrete_p : Array[String] = []
  let continuous_p : Array[String] = []
  let hurdle_p : Array[String] = []
  let hurdle_fdr : Array[String] = []
  let marginal_log_fc : Array[String] = []
  let classification : Array[String] = []
  for gene in 0..<result.gene_names.length() {
    detection_rate.push(result.detection_rates[gene].to_string())
    detected.push(result.n_detected[gene].to_string())
    discrete_p.push(mast_adv_option_string(result.discrete_p_values[gene]))
    continuous_p.push(mast_adv_option_string(result.continuous_p_values[gene]))
    hurdle_p.push(mast_adv_option_string(result.hurdle_p_values[gene]))
    hurdle_fdr.push(mast_adv_option_string(result.hurdle_fdr[gene]))
    marginal_log_fc.push(mast_adv_option_string(result.marginal_log_fc[gene]))
    classification.push(
      match result.hurdle_fdr[gene] {
        Some(value) =>
          if value <= config.fdr_threshold {
            "significant"
          } else {
            "not_significant"
          }
        None => "not_tested"
      },
    )
  }
  enriched.row_data[output_prefix + ".detectionRate"] = detection_rate
  enriched.row_data[output_prefix + ".nDetected"] = detected
  enriched.row_data[output_prefix + ".discreteP"] = discrete_p
  enriched.row_data[output_prefix + ".continuousP"] = continuous_p
  enriched.row_data[output_prefix + ".hurdleP"] = hurdle_p
  enriched.row_data[output_prefix + ".hurdleFdr"] = hurdle_fdr
  enriched.row_data[output_prefix + ".logFC"] = marginal_log_fc
  enriched.row_data[output_prefix + ".class"] = classification
  enriched.col_data[output_prefix + ".cdr"] = data.cdr.map(fn(value) {
    value.to_string()
  })
  enriched.metadata[output_prefix + ".assay"] = assay_name
  enriched.metadata[output_prefix + ".group"] = group_column
  enriched.metadata[output_prefix + ".reference"] = if reference == "" {
    groups[0]
  } else {
    reference
  }
  enriched.metadata[output_prefix + ".contrast"] = result.tested_names.join(",")
  enriched.metadata[output_prefix + ".priorVariance"] = result.prior_variance.to_string()
  enriched.metadata[output_prefix + ".priorDf"] = result.prior_df.to_string()
  enriched.metadata[output_prefix + ".tested"] = result.n_tested().to_string()
  enriched.metadata[output_prefix + ".significant"] = result
    .n_significant()
    .to_string()
  MastAdvancedSCEOutput::{ experiment: enriched, result }
}

///|
pub fn mast_advanced_example() -> (MastAdvancedData, Array[Int]) raise MastAdvancedError {
  let genes = 8
  let cells = 48
  let expression = mast_adv_zero_matrix(genes, cells)
  let groups : Array[String] = []
  for cell in 0..<cells {
    let treatment = cell >= 24
    groups.push(if treatment { "stimulated" } else { "control" })
    for gene in 0..<genes {
      let detected = if gene == 0 {
        if treatment {
          cell % 8 != 0
        } else {
          cell % 3 == 0
        }
      } else if gene == 1 {
        cell % 4 != 0
      } else if gene == 2 {
        cell % 5 != 0
      } else {
        (cell + gene) % 6 != 0
      }
      if detected {
        let baseline = 2.0 + gene.to_double() * 0.15
        let condition_effect = if treatment && gene == 0 {
          1.8
        } else if treatment && gene == 1 {
          1.2
        } else {
          0.0
        }
        let noise = ((cell * 17 + gene * 11) % 9).to_double() * 0.04
        expression[gene][cell] = baseline + condition_effect + noise
      }
    }
  }
  let gene_names = [
    "both_components", "continuous_only", "stable_1", "stable_2", "stable_3", "stable_4",
    "stable_5", "stable_6",
  ]
  let data = mast_advanced_from_groups(
    expression,
    groups,
    reference="control",
    gene_names~,
  )
  (data, [1])
}
