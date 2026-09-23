///|
/// Gene-level variance decomposition and repeated-measures differential
/// expression inspired by Bioconductor variancePartition.

///|
pub suberror VariancePartitionError {
  VariancePartitionError(String)
}

///|
pub struct VariancePartitionFixedEffect {
  name : String
  column_names : Array[String]
  columns : Array[Array[Double]]
} derive(Debug)
pub extend VariancePartitionFixedEffect with @moonbitlang/core/debug.Debug::{to_repr}

///|
pub struct VariancePartitionRandomEffect {
  name : String
  levels : Array[String]
  level_indices : Array[Int]
} derive(Debug)
pub extend VariancePartitionRandomEffect with @moonbitlang/core/debug.Debug::{to_repr}

///|
pub struct VariancePartitionDesign {
  sample_names : Array[String]
  fixed_effects : Array[VariancePartitionFixedEffect]
  random_effects : Array[VariancePartitionRandomEffect]
  matrix : Array[Array[Double]]
  coefficient_names : Array[String]
  fixed_starts : Array[Int]
  fixed_widths : Array[Int]
  intercept : Bool
} derive(Debug)
pub extend VariancePartitionDesign with @moonbitlang/core/debug.Debug::{to_repr}

///|
pub struct VariancePartitionBLUP {
  effect_name : String
  levels : Array[String]
  estimates : Array[Double]
} derive(Debug)
pub extend VariancePartitionBLUP with @moonbitlang/core/debug.Debug::{to_repr}

///|
pub struct VariancePartitionFit {
  gene_name : String
  coefficient_names : Array[String]
  coefficients : Array[Double]
  coefficient_covariance : Array[Array[Double]]
  fixed_effect_names : Array[String]
  fixed_variances : Array[Double]
  random_effect_names : Array[String]
  random_variances : Array[Double]
  residual_variance : Double
  component_names : Array[String]
  variance_fractions : Array[Double]
  blups : Array[VariancePartitionBLUP]
  fitted_values : Array[Double]
  residuals : Array[Double]
  log_likelihood : Double
  converged : Bool
  iterations : Int
  reml : Bool
} derive(Debug)
pub extend VariancePartitionFit with @moonbitlang/core/debug.Debug::{to_repr}

///|
pub struct VariancePartitionResult {
  gene_names : Array[String]
  component_names : Array[String]
  fractions : Array[Array[Double]]
  median_fractions : Array[Double]
  fits : Array[VariancePartitionFit]
} derive(Debug)
pub extend VariancePartitionResult with @moonbitlang/core/debug.Debug::{to_repr}

///|
pub struct DreamContrast {
  name : String
  coefficients : Array[Double]
} derive(Debug)
pub extend DreamContrast with @moonbitlang/core/debug.Debug::{to_repr}

///|
pub struct DreamGeneResult {
  gene_name : String
  estimate : Double
  standard_error : Double
  statistic : Double
  degrees_of_freedom : Double
  p_value : Double
  adjusted_p_value : Double
  converged : Bool
  random_variances : Array[Double]
  residual_variance : Double
} derive(Debug)
pub extend DreamGeneResult with @moonbitlang/core/debug.Debug::{to_repr}

///|
pub struct DreamResult {
  contrast_name : String
  coefficient_names : Array[String]
  genes : Array[DreamGeneResult]
} derive(Debug)
pub extend DreamResult with @moonbitlang/core/debug.Debug::{to_repr}

///|
priv struct VariancePartitionEvaluation {
  objective : Double
  beta : Array[Double]
  beta_covariance : Array[Array[Double]]
  fixed_fitted : Array[Double]
  marginal_residual : Array[Double]
  covariance_cholesky : Array[Array[Double]]
}

///|
priv struct VariancePartitionCoreFit {
  evaluation : VariancePartitionEvaluation
  log_variances : Array[Double]
  variances : Array[Double]
  weights : Array[Double]
  converged : Bool
  iterations : Int
}

///|
fn vp_is_finite(value : Double) -> Bool {
  value == value && value.abs() < 1.0e300
}

///|
fn vp_find_string(values : Array[String], target : String) -> Int {
  for index in 0..<values.length() {
    if values[index] == target {
      return index
    }
  }
  -1
}

///|
fn vp_zero_matrix(rows : Int, columns : Int) -> Array[Array[Double]] {
  let output : Array[Array[Double]] = []
  for _ in 0..<rows {
    output.push(Array::make(columns, 0.0))
  }
  output
}

///|
fn vp_copy_matrix(matrix : Array[Array[Double]]) -> Array[Array[Double]] {
  let output : Array[Array[Double]] = []
  for row in matrix {
    output.push(row.copy())
  }
  output
}

///|
fn vp_dot(left : Array[Double], right : Array[Double]) -> Double {
  let mut total = 0.0
  for index in 0..<left.length() {
    total = total + left[index] * right[index]
  }
  total
}

///|
fn vp_mean(values : Array[Double]) -> Double {
  if values.length() == 0 {
    return 0.0
  }
  let mut total = 0.0
  for value in values {
    total = total + value
  }
  total / values.length().to_double()
}

///|
fn vp_sample_variance(values : Array[Double]) -> Double {
  if values.length() < 2 {
    return 0.0
  }
  let mean = vp_mean(values)
  let mut total = 0.0
  for value in values {
    let difference = value - mean
    total = total + difference * difference
  }
  total / (values.length() - 1).to_double()
}

///|
fn vp_matrix_vector(
  matrix : Array[Array[Double]],
  vector : Array[Double],
) -> Array[Double] {
  let output = Array::make(matrix.length(), 0.0)
  for row in 0..<matrix.length() {
    for column in 0..<vector.length() {
      output[row] = output[row] + matrix[row][column] * vector[column]
    }
  }
  output
}

///|
fn vp_cholesky(matrix : Array[Array[Double]]) -> Array[Array[Double]]? {
  let size = matrix.length()
  if size == 0 {
    return None
  }
  let lower = vp_zero_matrix(size, size)
  for row in 0..<size {
    if matrix[row].length() != size {
      return None
    }
    let column_limit = row + 1
    for column in 0..<column_limit {
      let mut value = matrix[row][column]
      for inner in 0..<column {
        value = value - lower[row][inner] * lower[column][inner]
      }
      if row == column {
        if value <= 1.0e-12 || !vp_is_finite(value) {
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
fn vp_cholesky_solve(
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
fn vp_cholesky_inverse(lower : Array[Array[Double]]) -> Array[Array[Double]] {
  let size = lower.length()
  let inverse = vp_zero_matrix(size, size)
  for column in 0..<size {
    let basis = Array::make(size, 0.0)
    basis[column] = 1.0
    let solution = vp_cholesky_solve(lower, basis)
    for row in 0..<size {
      inverse[row][column] = solution[row]
    }
  }
  inverse
}

///|
fn vp_log_determinant_from_cholesky(lower : Array[Array[Double]]) -> Double {
  let mut value = 0.0
  for index in 0..<lower.length() {
    value = value + 2.0 * @math.ln(lower[index][index])
  }
  value
}

///|
fn vp_general_inverse(matrix : Array[Array[Double]]) -> Array[Array[Double]]? {
  let size = matrix.length()
  if size == 0 {
    return None
  }
  let augmented = vp_zero_matrix(size, 2 * size)
  for row in 0..<size {
    if matrix[row].length() != size {
      return None
    }
    for column in 0..<size {
      augmented[row][column] = matrix[row][column]
    }
    augmented[row][size + row] = 1.0
  }
  for pivot_column in 0..<size {
    let mut pivot_row = pivot_column
    let mut pivot_value = augmented[pivot_row][pivot_column].abs()
    for candidate in (pivot_column + 1)..<size {
      let value = augmented[candidate][pivot_column].abs()
      if value > pivot_value {
        pivot_row = candidate
        pivot_value = value
      }
    }
    if pivot_value < 1.0e-12 {
      return None
    }
    if pivot_row != pivot_column {
      let temporary = augmented[pivot_column]
      augmented[pivot_column] = augmented[pivot_row]
      augmented[pivot_row] = temporary
    }
    let divisor = augmented[pivot_column][pivot_column]
    for column in 0..<(2 * size) {
      augmented[pivot_column][column] = augmented[pivot_column][column] /
        divisor
    }
    for row in 0..<size {
      if row != pivot_column {
        let factor = augmented[row][pivot_column]
        for column in 0..<(2 * size) {
          augmented[row][column] = augmented[row][column] -
            factor * augmented[pivot_column][column]
        }
      }
    }
  }
  let inverse = vp_zero_matrix(size, size)
  for row in 0..<size {
    for column in 0..<size {
      inverse[row][column] = augmented[row][size + column]
    }
  }
  Some(inverse)
}

///|
fn vp_validate_name(
  name : String,
  kind : String,
) -> Unit raise VariancePartitionError {
  if name == "" {
    raise VariancePartitionError(
      "variancePartition " + kind + " name must not be empty",
    )
  }
}

///|
pub fn vp_numeric_effect(
  name : String,
  values : Array[Double],
) -> VariancePartitionFixedEffect raise VariancePartitionError {
  vp_validate_name(name, "fixed effect")
  if values.length() == 0 {
    raise VariancePartitionError(
      "variancePartition numeric effect must contain observations",
    )
  }
  let columns : Array[Array[Double]] = []
  for value in values {
    if !vp_is_finite(value) {
      raise VariancePartitionError(
        "variancePartition numeric effects must be finite",
      )
    }
    columns.push([value])
  }
  VariancePartitionFixedEffect::{ name, column_names: [name], columns, }
}

///|
pub fn vp_categorical_effect(
  name : String,
  values : Array[String],
  reference? : String = "",
) -> VariancePartitionFixedEffect raise VariancePartitionError {
  vp_validate_name(name, "fixed effect")
  if values.length() == 0 {
    raise VariancePartitionError(
      "variancePartition categorical effect must contain observations",
    )
  }
  let levels : Array[String] = []
  for value in values {
    if value == "" {
      raise VariancePartitionError(
        "variancePartition categorical levels must not be empty",
      )
    }
    if !levels.contains(value) {
      levels.push(value)
    }
  }
  if levels.length() < 2 {
    raise VariancePartitionError(
      "variancePartition categorical effect requires at least two levels",
    )
  }
  let selected_reference = if reference == "" { levels[0] } else { reference }
  if !levels.contains(selected_reference) {
    raise VariancePartitionError(
      "variancePartition categorical reference level is absent",
    )
  }
  let encoded_levels : Array[String] = []
  for level in levels {
    if level != selected_reference {
      encoded_levels.push(level)
    }
  }
  let columns : Array[Array[Double]] = []
  for value in values {
    let row = Array::make(encoded_levels.length(), 0.0)
    for column in 0..<encoded_levels.length() {
      if value == encoded_levels[column] {
        row[column] = 1.0
      }
    }
    columns.push(row)
  }
  let column_names : Array[String] = []
  for level in encoded_levels {
    column_names.push(name + ":" + level)
  }
  VariancePartitionFixedEffect::{ name, column_names, columns, }
}

///|
pub fn vp_random_effect(
  name : String,
  values : Array[String],
) -> VariancePartitionRandomEffect raise VariancePartitionError {
  vp_validate_name(name, "random effect")
  if values.length() == 0 {
    raise VariancePartitionError(
      "variancePartition random effect must contain observations",
    )
  }
  let levels : Array[String] = []
  let indices : Array[Int] = []
  for value in values {
    if value == "" {
      raise VariancePartitionError(
        "variancePartition random-effect levels must not be empty",
      )
    }
    let existing = vp_find_string(levels, value)
    if existing >= 0 {
      indices.push(existing)
    } else {
      levels.push(value)
      indices.push(levels.length() - 1)
    }
  }
  if levels.length() < 2 {
    raise VariancePartitionError(
      "variancePartition random effect requires at least two levels",
    )
  }
  if levels.length() >= values.length() {
    raise VariancePartitionError(
      "variancePartition random effects require replicated levels",
    )
  }
  VariancePartitionRandomEffect::{ name, levels, level_indices: indices, }
}

///|
pub fn vp_design(
  sample_names : Array[String],
  fixed_effects : Array[VariancePartitionFixedEffect],
  random_effects : Array[VariancePartitionRandomEffect],
  intercept? : Bool = true,
) -> VariancePartitionDesign raise VariancePartitionError {
  let observations = sample_names.length()
  if observations < 3 {
    raise VariancePartitionError(
      "variancePartition design requires at least three samples",
    )
  }
  let seen_samples : Map[String, Bool] = Map([], capacity=observations)
  for sample in sample_names {
    if sample == "" || seen_samples.contains(sample) {
      raise VariancePartitionError(
        "variancePartition sample names must be non-empty and unique",
      )
    }
    seen_samples[sample] = true
  }
  let seen_effects : Map[String, Bool] = Map(
    [],
    capacity=fixed_effects.length() + random_effects.length(),
  )
  let coefficient_names : Array[String] = []
  if intercept {
    coefficient_names.push("(Intercept)")
  }
  let fixed_starts : Array[Int] = []
  let fixed_widths : Array[Int] = []
  let copied_fixed : Array[VariancePartitionFixedEffect] = []
  let mut coefficient_count = if intercept { 1 } else { 0 }
  for effect in fixed_effects {
    vp_validate_name(effect.name, "fixed effect")
    if seen_effects.contains(effect.name) {
      raise VariancePartitionError(
        "variancePartition effect names must be unique",
      )
    }
    seen_effects[effect.name] = true
    if effect.columns.length() != observations ||
      effect.column_names.length() == 0 {
      raise VariancePartitionError(
        "variancePartition fixed-effect dimensions do not match samples",
      )
    }
    let width = effect.column_names.length()
    for row in effect.columns {
      if row.length() != width {
        raise VariancePartitionError(
          "variancePartition fixed-effect matrix must be rectangular",
        )
      }
      for value in row {
        if !vp_is_finite(value) {
          raise VariancePartitionError(
            "variancePartition fixed-effect values must be finite",
          )
        }
      }
    }
    fixed_starts.push(coefficient_count)
    fixed_widths.push(width)
    coefficient_count = coefficient_count + width
    for column_name in effect.column_names {
      coefficient_names.push(column_name)
    }
    copied_fixed.push(VariancePartitionFixedEffect::{
      name: effect.name,
      column_names: effect.column_names.copy(),
      columns: vp_copy_matrix(effect.columns),
    })
  }
  let copied_random : Array[VariancePartitionRandomEffect] = []
  for effect in random_effects {
    vp_validate_name(effect.name, "random effect")
    if seen_effects.contains(effect.name) {
      raise VariancePartitionError(
        "variancePartition effect names must be unique",
      )
    }
    seen_effects[effect.name] = true
    if effect.level_indices.length() != observations {
      raise VariancePartitionError(
        "variancePartition random-effect dimensions do not match samples",
      )
    }
    if effect.levels.length() < 2 || effect.levels.length() >= observations {
      raise VariancePartitionError(
        "variancePartition random effects require replicated levels",
      )
    }
    for index in effect.level_indices {
      if index < 0 || index >= effect.levels.length() {
        raise VariancePartitionError(
          "variancePartition random-effect level index is invalid",
        )
      }
    }
    copied_random.push(VariancePartitionRandomEffect::{
      name: effect.name,
      levels: effect.levels.copy(),
      level_indices: effect.level_indices.copy(),
    })
  }
  if coefficient_count == 0 || coefficient_count >= observations {
    raise VariancePartitionError(
      "variancePartition fixed-effect design has no residual degrees of freedom",
    )
  }
  let matrix = vp_zero_matrix(observations, coefficient_count)
  for sample in 0..<observations {
    let mut column = 0
    if intercept {
      matrix[sample][0] = 1.0
      column = 1
    }
    for effect in copied_fixed {
      for value in effect.columns[sample] {
        matrix[sample][column] = value
        column = column + 1
      }
    }
  }
  let information = vp_zero_matrix(coefficient_count, coefficient_count)
  for sample in 0..<observations {
    for left in 0..<coefficient_count {
      for right in 0..<coefficient_count {
        information[left][right] = information[left][right] +
          matrix[sample][left] * matrix[sample][right]
      }
    }
  }
  if vp_cholesky(information) is None {
    raise VariancePartitionError(
      "variancePartition fixed-effect design is rank deficient",
    )
  }
  VariancePartitionDesign::{
    sample_names: sample_names.copy(),
    fixed_effects: copied_fixed,
    random_effects: copied_random,
    matrix,
    coefficient_names,
    fixed_starts,
    fixed_widths,
    intercept,
  }
}

///|
pub fn vp_coefficient_contrast(
  design : VariancePartitionDesign,
  name : String,
  coefficient_name : String,
) -> DreamContrast raise VariancePartitionError {
  vp_validate_name(name, "contrast")
  let index = vp_find_string(design.coefficient_names, coefficient_name)
  if index < 0 {
    raise VariancePartitionError(
      "variancePartition contrast coefficient is absent from the design",
    )
  }
  let coefficients = Array::make(design.coefficient_names.length(), 0.0)
  coefficients[index] = 1.0
  DreamContrast::{ name, coefficients, }
}

///|
pub fn vp_contrast(
  name : String,
  coefficients : Array[Double],
) -> DreamContrast raise VariancePartitionError {
  vp_validate_name(name, "contrast")
  if coefficients.length() == 0 {
    raise VariancePartitionError(
      "variancePartition contrast must contain coefficients",
    )
  }
  let mut magnitude = 0.0
  for coefficient in coefficients {
    if !vp_is_finite(coefficient) {
      raise VariancePartitionError(
        "variancePartition contrast coefficients must be finite",
      )
    }
    magnitude = magnitude + coefficient.abs()
  }
  if magnitude == 0.0 {
    raise VariancePartitionError(
      "variancePartition contrast must not be all zero",
    )
  }
  DreamContrast::{ name, coefficients: coefficients.copy(), }
}

///|
fn vp_scale_weights(
  weights : Array[Double],
  observations : Int,
) -> Array[Double] raise VariancePartitionError {
  if weights.length() == 0 {
    return Array::make(observations, 1.0)
  }
  if weights.length() != observations {
    raise VariancePartitionError(
      "variancePartition precision weights must match the samples",
    )
  }
  let mut total = 0.0
  for weight in weights {
    if !vp_is_finite(weight) || weight <= 0.0 {
      raise VariancePartitionError(
        "variancePartition precision weights must be finite and positive",
      )
    }
    total = total + weight
  }
  let mean = total / observations.to_double()
  let scaled = Array::make(observations, 0.0)
  for index in 0..<observations {
    scaled[index] = weights[index] / mean
  }
  scaled
}

///|
fn vp_covariance_matrix(
  design : VariancePartitionDesign,
  variances : Array[Double],
  weights : Array[Double],
) -> Array[Array[Double]] {
  let observations = design.sample_names.length()
  let covariance = vp_zero_matrix(observations, observations)
  for effect_index in 0..<design.random_effects.length() {
    let effect = design.random_effects[effect_index]
    let variance = variances[effect_index]
    for row in 0..<observations {
      for column in 0..<observations {
        if effect.level_indices[row] == effect.level_indices[column] {
          covariance[row][column] = covariance[row][column] + variance
        }
      }
    }
  }
  let residual_variance = variances[variances.length() - 1]
  for observation in 0..<observations {
    covariance[observation][observation] = covariance[observation][observation] +
      residual_variance / weights[observation]
  }
  covariance
}

///|
fn vp_evaluate_model(
  expression : Array[Double],
  design : VariancePartitionDesign,
  log_variances : Array[Double],
  weights : Array[Double],
  reml : Bool,
) -> VariancePartitionEvaluation? {
  let observations = expression.length()
  let parameters = design.coefficient_names.length()
  let variances = Array::make(log_variances.length(), 0.0)
  for index in 0..<log_variances.length() {
    variances[index] = @math.exp(log_variances[index])
  }
  let covariance = vp_covariance_matrix(design, variances, weights)
  let covariance_cholesky = match vp_cholesky(covariance) {
    Some(value) => value
    None => return None
  }
  let inverse_expression = vp_cholesky_solve(covariance_cholesky, expression)
  let inverse_design = vp_zero_matrix(observations, parameters)
  for column in 0..<parameters {
    let design_column = Array::make(observations, 0.0)
    for row in 0..<observations {
      design_column[row] = design.matrix[row][column]
    }
    let solved = vp_cholesky_solve(covariance_cholesky, design_column)
    for row in 0..<observations {
      inverse_design[row][column] = solved[row]
    }
  }
  let information = vp_zero_matrix(parameters, parameters)
  let right_hand_side = Array::make(parameters, 0.0)
  for left in 0..<parameters {
    for row in 0..<observations {
      right_hand_side[left] = right_hand_side[left] +
        design.matrix[row][left] * inverse_expression[row]
    }
    for right in 0..<parameters {
      for row in 0..<observations {
        information[left][right] = information[left][right] +
          design.matrix[row][left] * inverse_design[row][right]
      }
    }
  }
  let information_cholesky = match vp_cholesky(information) {
    Some(value) => value
    None => return None
  }
  let beta = vp_cholesky_solve(information_cholesky, right_hand_side)
  let beta_covariance = vp_cholesky_inverse(information_cholesky)
  let fixed_fitted = vp_matrix_vector(design.matrix, beta)
  let marginal_residual = Array::make(observations, 0.0)
  for index in 0..<observations {
    marginal_residual[index] = expression[index] - fixed_fitted[index]
  }
  let solved_residual = vp_cholesky_solve(
    covariance_cholesky, marginal_residual,
  )
  let quadratic = vp_dot(marginal_residual, solved_residual)
  let degrees = if reml { observations - parameters } else { observations }
  let mut objective = 0.5 *
    (
      degrees.to_double() * @math.ln(6.283185307179586) +
      vp_log_determinant_from_cholesky(covariance_cholesky) +
      quadratic
    )
  if reml {
    objective = objective +
      0.5 * vp_log_determinant_from_cholesky(information_cholesky)
  }
  if !vp_is_finite(objective) {
    return None
  }
  Some(VariancePartitionEvaluation::{
    objective,
    beta,
    beta_covariance,
    fixed_fitted,
    marginal_residual,
    covariance_cholesky,
  })
}

///|
fn vp_fit_core(
  expression : Array[Double],
  design : VariancePartitionDesign,
  weights : Array[Double],
  reml : Bool,
  max_iterations : Int,
  tolerance : Double,
) -> VariancePartitionCoreFit raise VariancePartitionError {
  let observations = design.sample_names.length()
  if expression.length() != observations {
    raise VariancePartitionError(
      "variancePartition expression values must match the samples",
    )
  }
  for value in expression {
    if !vp_is_finite(value) {
      raise VariancePartitionError(
        "variancePartition expression values must be finite",
      )
    }
  }
  if max_iterations <= 0 {
    raise VariancePartitionError(
      "variancePartition max_iterations must be positive",
    )
  }
  if !vp_is_finite(tolerance) || tolerance <= 0.0 || tolerance >= 1.0 {
    raise VariancePartitionError(
      "variancePartition tolerance must be between zero and one",
    )
  }
  let scaled_weights = vp_scale_weights(weights, observations)
  let component_count = design.random_effects.length() + 1
  let observed_variance = vp_sample_variance(expression).max(1.0e-8)
  let initial = Array::make(component_count, 0.0)
  if component_count == 1 {
    initial[0] = observed_variance
  } else {
    let random_total = observed_variance * 0.5
    for index in 0..<(component_count - 1) {
      initial[index] = random_total / (component_count - 1).to_double()
    }
    initial[component_count - 1] = observed_variance * 0.5
  }
  let log_variances = Array::make(component_count, 0.0)
  for index in 0..<component_count {
    log_variances[index] = @math.ln(initial[index].max(1.0e-10))
  }
  let lower_bound = @math.ln((observed_variance * 1.0e-8).max(1.0e-12))
  let upper_bound = @math.ln((observed_variance * 1.0e4).max(1.0e-4))
  let mut current = match
    vp_evaluate_model(expression, design, log_variances, scaled_weights, reml) {
    Some(value) => value
    None =>
      raise VariancePartitionError(
        "variancePartition initial mixed model is not positive definite",
      )
  }
  let mut step = 1.0
  let mut iterations = 0
  while iterations < max_iterations && step > tolerance {
    let mut improved = false
    for parameter in 0..<component_count {
      let original = log_variances[parameter]
      let mut best_parameter = original
      let mut best_evaluation = current
      let lower_candidate = (original - step).max(lower_bound)
      let upper_candidate = (original + step).min(upper_bound)
      if lower_candidate != original {
        log_variances[parameter] = lower_candidate
        match
          vp_evaluate_model(
            expression, design, log_variances, scaled_weights, reml,
          ) {
          Some(candidate) =>
            if candidate.objective + 1.0e-10 < best_evaluation.objective {
              best_parameter = lower_candidate
              best_evaluation = candidate
            }
          None => ()
        }
      }
      if upper_candidate != original {
        log_variances[parameter] = upper_candidate
        match
          vp_evaluate_model(
            expression, design, log_variances, scaled_weights, reml,
          ) {
          Some(candidate) =>
            if candidate.objective + 1.0e-10 < best_evaluation.objective {
              best_parameter = upper_candidate
              best_evaluation = candidate
            }
          None => ()
        }
      }
      log_variances[parameter] = best_parameter
      if best_parameter != original {
        current = best_evaluation
        improved = true
      }
    }
    if !improved {
      step = step * 0.5
    }
    iterations = iterations + 1
  }
  let variances = Array::make(component_count, 0.0)
  for index in 0..<component_count {
    variances[index] = @math.exp(log_variances[index])
  }
  VariancePartitionCoreFit::{
    evaluation: current,
    log_variances,
    variances,
    weights: scaled_weights,
    converged: step <= tolerance,
    iterations,
  }
}

///|
fn vp_blups(
  core : VariancePartitionCoreFit,
  design : VariancePartitionDesign,
) -> Array[VariancePartitionBLUP] {
  let solved_residual = vp_cholesky_solve(
    core.evaluation.covariance_cholesky,
    core.evaluation.marginal_residual,
  )
  let output : Array[VariancePartitionBLUP] = []
  for effect_index in 0..<design.random_effects.length() {
    let effect = design.random_effects[effect_index]
    let estimates = Array::make(effect.levels.length(), 0.0)
    for observation in 0..<effect.level_indices.length() {
      let level = effect.level_indices[observation]
      estimates[level] = estimates[level] +
        core.variances[effect_index] * solved_residual[observation]
    }
    output.push(VariancePartitionBLUP::{
      effect_name: effect.name,
      levels: effect.levels.copy(),
      estimates,
    })
  }
  output
}

///|
fn vp_public_fit(
  gene_name : String,
  expression : Array[Double],
  design : VariancePartitionDesign,
  core : VariancePartitionCoreFit,
  reml : Bool,
) -> VariancePartitionFit {
  let fixed_variances : Array[Double] = []
  let fixed_names : Array[String] = []
  for effect_index in 0..<design.fixed_effects.length() {
    let effect = design.fixed_effects[effect_index]
    let prediction = Array::make(expression.length(), 0.0)
    let start = design.fixed_starts[effect_index]
    let width = design.fixed_widths[effect_index]
    for observation in 0..<expression.length() {
      for column in 0..<width {
        prediction[observation] = prediction[observation] +
          effect.columns[observation][column] *
          core.evaluation.beta[start + column]
      }
    }
    fixed_names.push(effect.name)
    fixed_variances.push(vp_sample_variance(prediction))
  }
  let random_names : Array[String] = []
  let random_variances : Array[Double] = []
  for index in 0..<design.random_effects.length() {
    random_names.push(design.random_effects[index].name)
    random_variances.push(core.variances[index])
  }
  let residual_variance = core.variances[core.variances.length() - 1]
  let component_names : Array[String] = []
  let component_variances : Array[Double] = []
  for index in 0..<fixed_names.length() {
    component_names.push(fixed_names[index])
    component_variances.push(fixed_variances[index])
  }
  for index in 0..<random_names.length() {
    component_names.push(random_names[index])
    component_variances.push(random_variances[index])
  }
  component_names.push("Residuals")
  component_variances.push(residual_variance)
  let mut total_variance = 0.0
  for value in component_variances {
    total_variance = total_variance + value
  }
  let fractions = Array::make(component_variances.length(), 0.0)
  if total_variance <= 1.0e-20 {
    fractions[fractions.length() - 1] = 1.0
  } else {
    for index in 0..<fractions.length() {
      fractions[index] = component_variances[index] / total_variance
    }
  }
  let blups = vp_blups(core, design)
  let fitted_values = core.evaluation.fixed_fitted.copy()
  for effect_index in 0..<design.random_effects.length() {
    let effect = design.random_effects[effect_index]
    for observation in 0..<expression.length() {
      fitted_values[observation] = fitted_values[observation] +
        blups[effect_index].estimates[effect.level_indices[observation]]
    }
  }
  let residuals = Array::make(expression.length(), 0.0)
  for observation in 0..<expression.length() {
    residuals[observation] = expression[observation] -
      fitted_values[observation]
  }
  VariancePartitionFit::{
    gene_name,
    coefficient_names: design.coefficient_names.copy(),
    coefficients: core.evaluation.beta.copy(),
    coefficient_covariance: vp_copy_matrix(core.evaluation.beta_covariance),
    fixed_effect_names: fixed_names,
    fixed_variances,
    random_effect_names: random_names,
    random_variances,
    residual_variance,
    component_names,
    variance_fractions: fractions,
    blups,
    fitted_values,
    residuals,
    log_likelihood: -core.evaluation.objective,
    converged: core.converged,
    iterations: core.iterations,
    reml,
  }
}

///|
pub fn vp_fit_gene(
  expression : Array[Double],
  gene_name? : String = "Gene",
  design : VariancePartitionDesign,
  weights? : Array[Double] = [],
  reml? : Bool = true,
  max_iterations? : Int = 100,
  tolerance? : Double = 0.005,
) -> VariancePartitionFit raise VariancePartitionError {
  if gene_name == "" {
    raise VariancePartitionError(
      "variancePartition gene name must not be empty",
    )
  }
  let core = vp_fit_core(
    expression, design, weights, reml, max_iterations, tolerance,
  )
  vp_public_fit(gene_name, expression, design, core, reml)
}

///|
pub fn VariancePartitionFit::fraction(
  self : VariancePartitionFit,
  component : String,
) -> Double? {
  let index = vp_find_string(self.component_names, component)
  if index < 0 {
    None
  } else {
    Some(self.variance_fractions[index])
  }
}

///|
pub fn VariancePartitionFit::blup(
  self : VariancePartitionFit,
  effect : String,
  level : String,
) -> Double? {
  for values in self.blups {
    if values.effect_name == effect {
      let index = vp_find_string(values.levels, level)
      if index >= 0 {
        return Some(values.estimates[index])
      }
      return None
    }
  }
  None
}

///|
fn vp_validate_expression_matrix(
  expression : Array[Array[Double]],
  gene_names : Array[String],
  samples : Int,
  weights : Array[Array[Double]],
) -> Array[String] raise VariancePartitionError {
  if expression.length() == 0 {
    raise VariancePartitionError(
      "variancePartition expression matrix must contain genes",
    )
  }
  let names : Array[String] = []
  if gene_names.length() == 0 {
    for index in 0..<expression.length() {
      names.push("Gene" + (index + 1).to_string())
    }
  } else {
    if gene_names.length() != expression.length() {
      raise VariancePartitionError(
        "variancePartition gene names must match expression rows",
      )
    }
    let seen : Map[String, Bool] = Map([], capacity=gene_names.length())
    for name in gene_names {
      if name == "" || seen.contains(name) {
        raise VariancePartitionError(
          "variancePartition gene names must be non-empty and unique",
        )
      }
      seen[name] = true
      names.push(name)
    }
  }
  for row in expression {
    if row.length() != samples {
      raise VariancePartitionError(
        "variancePartition expression matrix must be rectangular and match the design",
      )
    }
  }
  if weights.length() > 0 {
    if weights.length() != expression.length() {
      raise VariancePartitionError(
        "variancePartition weights must match expression rows",
      )
    }
    for row in weights {
      if row.length() != samples {
        raise VariancePartitionError(
          "variancePartition weights must match expression columns",
        )
      }
    }
  }
  names
}

///|
fn vp_median(values : Array[Double]) -> Double {
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
  if sorted.length() % 2 == 1 {
    sorted[middle]
  } else {
    (sorted[middle - 1] + sorted[middle]) / 2.0
  }
}

///|
pub fn fit_extract_variance_partition(
  expression : Array[Array[Double]],
  gene_names? : Array[String] = [],
  design : VariancePartitionDesign,
  weights? : Array[Array[Double]] = [],
  reml? : Bool = false,
  max_iterations? : Int = 100,
  tolerance? : Double = 0.005,
) -> VariancePartitionResult raise VariancePartitionError {
  let names = vp_validate_expression_matrix(
    expression,
    gene_names,
    design.sample_names.length(),
    weights,
  )
  let fits : Array[VariancePartitionFit] = []
  let fractions : Array[Array[Double]] = []
  for gene in 0..<expression.length() {
    let gene_weights = if weights.length() == 0 { [] } else { weights[gene] }
    let fit = vp_fit_gene(
      expression[gene],
      gene_name=names[gene],
      design,
      weights=gene_weights,
      reml~,
      max_iterations~,
      tolerance~,
    )
    fits.push(fit)
    fractions.push(fit.variance_fractions.copy())
  }
  let component_names = fits[0].component_names.copy()
  let median_fractions = Array::make(component_names.length(), 0.0)
  for component in 0..<component_names.length() {
    let values : Array[Double] = []
    for row in fractions {
      values.push(row[component])
    }
    median_fractions[component] = vp_median(values)
  }
  VariancePartitionResult::{
    gene_names: names,
    component_names,
    fractions,
    median_fractions,
    fits,
  }
}

///|
pub fn fit_extract_variance_partition_se(
  experiment : SummarizedExperiment,
  assay_name : String,
  gene_names? : Array[String] = [],
  design : VariancePartitionDesign,
  weights_assay? : String = "",
  reml? : Bool = false,
  max_iterations? : Int = 100,
  tolerance? : Double = 0.005,
) -> VariancePartitionResult raise VariancePartitionError {
  let expression = match se_assay(experiment, assay_name) {
    Some(value) => value
    None =>
      raise VariancePartitionError(
        "variancePartition assay is absent from SummarizedExperiment",
      )
  }
  let weights = if weights_assay == "" {
    []
  } else {
    match se_assay(experiment, weights_assay) {
      Some(value) => value
      None =>
        raise VariancePartitionError(
          "variancePartition weights assay is absent from SummarizedExperiment",
        )
    }
  }
  fit_extract_variance_partition(
    expression,
    design,
    gene_names~,
    weights~,
    reml~,
    max_iterations~,
    tolerance~,
  )
}

///|
pub fn VariancePartitionResult::fraction(
  self : VariancePartitionResult,
  gene : String,
  component : String,
) -> Double? {
  let gene_index = vp_find_string(self.gene_names, gene)
  let component_index = vp_find_string(self.component_names, component)
  if gene_index < 0 || component_index < 0 {
    None
  } else {
    Some(self.fractions[gene_index][component_index])
  }
}

///|
pub fn VariancePartitionResult::summary(
  self : VariancePartitionResult,
) -> String {
  let mut output = "variancePartition Summary\n"
  output = output + "Genes: " + self.gene_names.length().to_string() + "\n"
  for index in 0..<self.component_names.length() {
    output = output +
      self.component_names[index] +
      " median fraction: " +
      self.median_fractions[index].to_string() +
      "\n"
  }
  output
}

///|
fn vp_contrast_variance(
  covariance : Array[Array[Double]],
  contrast : Array[Double],
) -> Double {
  let mut value = 0.0
  for row in 0..<contrast.length() {
    for column in 0..<contrast.length() {
      value = value + contrast[row] * covariance[row][column] * contrast[column]
    }
  }
  value
}

///|
fn vp_objective_at(
  expression : Array[Double],
  design : VariancePartitionDesign,
  log_variances : Array[Double],
  weights : Array[Double],
  reml : Bool,
) -> Double? {
  match vp_evaluate_model(expression, design, log_variances, weights, reml) {
    Some(value) => Some(value.objective)
    None => None
  }
}

///|
fn vp_satterthwaite_df(
  expression : Array[Double],
  design : VariancePartitionDesign,
  contrast : Array[Double],
  core : VariancePartitionCoreFit,
  reml : Bool,
) -> Double {
  let residual_df = (expression.length() - design.coefficient_names.length())
    .to_double()
    .max(1.0)
  let base_variance = vp_contrast_variance(
    core.evaluation.beta_covariance,
    contrast,
  )
  if base_variance <= 1.0e-20 {
    return residual_df
  }
  let parameters = core.log_variances.length()
  let step = 0.05
  let base_objective = core.evaluation.objective
  let hessian = vp_zero_matrix(parameters, parameters)
  let gradient = Array::make(parameters, 0.0)
  for left in 0..<parameters {
    let plus = core.log_variances.copy()
    let minus = core.log_variances.copy()
    plus[left] = plus[left] + step
    minus[left] = minus[left] - step
    let plus_evaluation = vp_evaluate_model(
      expression,
      design,
      plus,
      core.weights,
      reml,
    )
    let minus_evaluation = vp_evaluate_model(
      expression,
      design,
      minus,
      core.weights,
      reml,
    )
    match (plus_evaluation, minus_evaluation) {
      (Some(positive), Some(negative)) => {
        let positive_variance = vp_contrast_variance(
          positive.beta_covariance,
          contrast,
        )
        let negative_variance = vp_contrast_variance(
          negative.beta_covariance,
          contrast,
        )
        gradient[left] = (positive_variance - negative_variance) / (2.0 * step)
        hessian[left][left] = (
            positive.objective - 2.0 * base_objective + negative.objective
          ) /
          (step * step)
      }
      _ => return residual_df
    }
    for right in (left + 1)..<parameters {
      let plus_plus = core.log_variances.copy()
      let plus_minus = core.log_variances.copy()
      let minus_plus = core.log_variances.copy()
      let minus_minus = core.log_variances.copy()
      plus_plus[left] = plus_plus[left] + step
      plus_plus[right] = plus_plus[right] + step
      plus_minus[left] = plus_minus[left] + step
      plus_minus[right] = plus_minus[right] - step
      minus_plus[left] = minus_plus[left] - step
      minus_plus[right] = minus_plus[right] + step
      minus_minus[left] = minus_minus[left] - step
      minus_minus[right] = minus_minus[right] - step
      let fpp = vp_objective_at(
        expression,
        design,
        plus_plus,
        core.weights,
        reml,
      )
      let fpm = vp_objective_at(
        expression,
        design,
        plus_minus,
        core.weights,
        reml,
      )
      let fmp = vp_objective_at(
        expression,
        design,
        minus_plus,
        core.weights,
        reml,
      )
      let fmm = vp_objective_at(
        expression,
        design,
        minus_minus,
        core.weights,
        reml,
      )
      match (fpp, fpm, fmp, fmm) {
        (Some(a), Some(b), Some(c), Some(d)) => {
          let value = (a - b - c + d) / (4.0 * step * step)
          hessian[left][right] = value
          hessian[right][left] = value
        }
        _ => return residual_df
      }
    }
  }
  for index in 0..<parameters {
    hessian[index][index] = hessian[index][index] + 1.0e-8
  }
  let covariance_parameters = match vp_general_inverse(hessian) {
    Some(value) => value
    None => return residual_df
  }
  let mut variance_of_variance = 0.0
  for left in 0..<parameters {
    for right in 0..<parameters {
      variance_of_variance = variance_of_variance +
        gradient[left] * covariance_parameters[left][right] * gradient[right]
    }
  }
  if !vp_is_finite(variance_of_variance) || variance_of_variance <= 1.0e-20 {
    return residual_df
  }
  let estimated = 2.0 * base_variance * base_variance / variance_of_variance
  if !vp_is_finite(estimated) {
    residual_df
  } else {
    estimated.max(1.0).min(residual_df)
  }
}

///|
fn vp_log_gamma(value : Double) -> Double {
  let coefficients = [
    76.18009172947146, -86.50532032941677, 24.01409824083091, -1.231739572450155,
    0.001208650973866179, -0.000005395239384953,
  ]
  let mut shifted = value
  let temporary = value + 5.5
  let leading = (value + 0.5) * @math.ln(temporary) - temporary
  let mut series = 1.000000000190015
  for coefficient in coefficients {
    shifted = shifted + 1.0
    series = series + coefficient / shifted
  }
  leading + @math.ln(2.5066282746310005 * series)
}

///|
fn vp_beta_continued_fraction(x : Double, a : Double, b : Double) -> Double {
  let maximum_iterations = 200
  let epsilon = 3.0e-14
  let minimum = 1.0e-300
  let mut c = 1.0
  let mut d = 1.0 - (a + b) * x / (a + 1.0)
  if d.abs() < minimum {
    d = minimum
  }
  d = 1.0 / d
  let mut fraction = d
  let mut iteration = 1
  while iteration <= maximum_iterations {
    let doubled = (2 * iteration).to_double()
    let iteration_value = iteration.to_double()
    let mut numerator = iteration_value *
      (b - iteration_value) *
      x /
      ((a + doubled - 1.0) * (a + doubled))
    d = 1.0 + numerator * d
    if d.abs() < minimum {
      d = minimum
    }
    c = 1.0 + numerator / c
    if c.abs() < minimum {
      c = minimum
    }
    d = 1.0 / d
    fraction = fraction * d * c
    numerator = -(a + iteration_value) *
      (a + b + iteration_value) *
      x /
      ((a + doubled) * (a + doubled + 1.0))
    d = 1.0 + numerator * d
    if d.abs() < minimum {
      d = minimum
    }
    c = 1.0 + numerator / c
    if c.abs() < minimum {
      c = minimum
    }
    d = 1.0 / d
    let delta = d * c
    fraction = fraction * delta
    if (delta - 1.0).abs() <= epsilon {
      break
    }
    iteration = iteration + 1
  }
  fraction
}

///|
fn vp_regularized_beta(x : Double, a : Double, b : Double) -> Double {
  if x <= 0.0 {
    return 0.0
  }
  if x >= 1.0 {
    return 1.0
  }
  let log_front = vp_log_gamma(a + b) -
    vp_log_gamma(a) -
    vp_log_gamma(b) +
    a * @math.ln(x) +
    b * @math.ln(1.0 - x)
  let front = @math.exp(log_front)
  let result = if x < (a + 1.0) / (a + b + 2.0) {
    front * vp_beta_continued_fraction(x, a, b) / a
  } else {
    1.0 - front * vp_beta_continued_fraction(1.0 - x, b, a) / b
  }
  result.max(0.0).min(1.0)
}

///|
fn vp_student_two_sided_p(statistic : Double, degrees : Double) -> Double {
  if !vp_is_finite(statistic) || degrees <= 0.0 {
    return 1.0
  }
  let x = degrees / (degrees + statistic * statistic)
  vp_regularized_beta(x, degrees / 2.0, 0.5)
}

///|
fn vp_bh_adjust(results : Array[DreamGeneResult]) -> Array[DreamGeneResult] {
  let count = results.length()
  let order : Array[Int] = []
  for index in 0..<count {
    order.push(index)
  }
  order.sort_by(fn(left : Int, right : Int) -> Int {
    if results[left].p_value < results[right].p_value {
      -1
    } else if results[left].p_value > results[right].p_value {
      1
    } else if left < right {
      -1
    } else if left > right {
      1
    } else {
      0
    }
  })
  let adjusted = Array::make(count, 1.0)
  let mut running = 1.0
  let mut rank = count
  while rank > 0 {
    let original = order[rank - 1]
    let candidate = (results[original].p_value *
    count.to_double() /
    rank.to_double()).min(1.0)
    running = running.min(candidate)
    adjusted[original] = running.max(results[original].p_value)
    rank = rank - 1
  }
  let output : Array[DreamGeneResult] = []
  for index in 0..<count {
    let value = results[index]
    output.push(DreamGeneResult::{
      gene_name: value.gene_name,
      estimate: value.estimate,
      standard_error: value.standard_error,
      statistic: value.statistic,
      degrees_of_freedom: value.degrees_of_freedom,
      p_value: value.p_value,
      adjusted_p_value: adjusted[index],
      converged: value.converged,
      random_variances: value.random_variances.copy(),
      residual_variance: value.residual_variance,
    })
  }
  output
}

///|
pub fn dream(
  expression : Array[Array[Double]],
  gene_names? : Array[String] = [],
  design : VariancePartitionDesign,
  contrast : DreamContrast,
  weights? : Array[Array[Double]] = [],
  reml? : Bool = true,
  max_iterations? : Int = 100,
  tolerance? : Double = 0.005,
) -> DreamResult raise VariancePartitionError {
  let names = vp_validate_expression_matrix(
    expression,
    gene_names,
    design.sample_names.length(),
    weights,
  )
  if contrast.coefficients.length() != design.coefficient_names.length() {
    raise VariancePartitionError(
      "variancePartition contrast must match the fixed-effect design",
    )
  }
  let raw_results : Array[DreamGeneResult] = []
  for gene in 0..<expression.length() {
    let gene_weights = if weights.length() == 0 { [] } else { weights[gene] }
    let core = vp_fit_core(
      expression[gene],
      design,
      gene_weights,
      reml,
      max_iterations,
      tolerance,
    )
    let estimate = vp_dot(contrast.coefficients, core.evaluation.beta)
    let contrast_variance = vp_contrast_variance(
      core.evaluation.beta_covariance,
      contrast.coefficients,
    ).max(0.0)
    let standard_error = contrast_variance.sqrt()
    let statistic = if standard_error <= 1.0e-20 {
      0.0
    } else {
      estimate / standard_error
    }
    let degrees = vp_satterthwaite_df(
      expression[gene],
      design,
      contrast.coefficients,
      core,
      reml,
    )
    let p_value = vp_student_two_sided_p(statistic, degrees)
    let random_variances : Array[Double] = []
    for index in 0..<design.random_effects.length() {
      random_variances.push(core.variances[index])
    }
    raw_results.push(DreamGeneResult::{
      gene_name: names[gene],
      estimate,
      standard_error,
      statistic,
      degrees_of_freedom: degrees,
      p_value,
      adjusted_p_value: 1.0,
      converged: core.converged,
      random_variances,
      residual_variance: core.variances[core.variances.length() - 1],
    })
  }
  DreamResult::{
    contrast_name: contrast.name,
    coefficient_names: design.coefficient_names.copy(),
    genes: vp_bh_adjust(raw_results),
  }
}

///|
pub fn dream_se(
  experiment : SummarizedExperiment,
  assay_name : String,
  gene_names? : Array[String] = [],
  design : VariancePartitionDesign,
  contrast : DreamContrast,
  weights_assay? : String = "",
  reml? : Bool = true,
  max_iterations? : Int = 100,
  tolerance? : Double = 0.005,
) -> DreamResult raise VariancePartitionError {
  let expression = match se_assay(experiment, assay_name) {
    Some(value) => value
    None =>
      raise VariancePartitionError(
        "variancePartition assay is absent from SummarizedExperiment",
      )
  }
  let weights = if weights_assay == "" {
    []
  } else {
    match se_assay(experiment, weights_assay) {
      Some(value) => value
      None =>
        raise VariancePartitionError(
          "variancePartition weights assay is absent from SummarizedExperiment",
        )
    }
  }
  dream(
    expression,
    design,
    contrast,
    gene_names~,
    weights~,
    reml~,
    max_iterations~,
    tolerance~,
  )
}

///|
pub fn DreamResult::top_table(
  self : DreamResult,
  maximum_fdr? : Double = 1.0,
) -> Array[DreamGeneResult] raise VariancePartitionError {
  if !vp_is_finite(maximum_fdr) || maximum_fdr < 0.0 || maximum_fdr > 1.0 {
    raise VariancePartitionError(
      "variancePartition maximum FDR must be between zero and one",
    )
  }
  let selected : Array[DreamGeneResult] = []
  for gene in self.genes {
    if gene.adjusted_p_value <= maximum_fdr {
      selected.push(gene)
    }
  }
  selected.sort_by(fn(left : DreamGeneResult, right : DreamGeneResult) -> Int {
    if left.p_value < right.p_value {
      -1
    } else if left.p_value > right.p_value {
      1
    } else if left.gene_name < right.gene_name {
      -1
    } else if left.gene_name > right.gene_name {
      1
    } else {
      0
    }
  })
  selected
}

///|
pub fn DreamResult::summary(self : DreamResult) -> String {
  let mut significant = 0
  for gene in self.genes {
    if gene.adjusted_p_value <= 0.05 {
      significant = significant + 1
    }
  }
  "dream Differential Expression Summary\n" +
  "Contrast: " +
  self.contrast_name +
  "\nGenes: " +
  self.genes.length().to_string() +
  "\nFDR <= 0.05: " +
  significant.to_string() +
  "\n"
}
