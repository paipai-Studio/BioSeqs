// SpatialDecon-inspired mixed-cell deconvolution for spatial expression data.

///|
pub suberror SpatialDeconError {
  SpatialDeconError(String)
}

///|
pub struct SpatialDeconProfile {
  gene_names : Array[String]
  cell_types : Array[String]
  values : Array[Array[Double]]
} derive(Debug)

///|
pub struct SpatialDeconData {
  gene_names : Array[String]
  spot_ids : Array[String]
  values : Array[Array[Double]]
  background : Array[Array[Double]]
  weights : Array[Array[Double]]
} derive(Debug)

///|
pub struct SpatialDeconConfig {
  residual_threshold : Double
  lower_threshold : Double
  signal_floor : Double
  rescale_profile : Bool
  profile_quantile : Double
  profile_target : Double
  refit_outliers : Bool
  max_iterations : Int
  line_search_steps : Int
  tolerance : Double
  ridge : Double
} derive(Debug)

///|
pub struct SpatialDeconAbundance {
  cell_type : String
  spot_id : String
  abundance : Double
  proportion : Double
  standard_error : Double
  p_value : Double
} derive(Debug)

///|
pub struct SpatialDeconResult {
  gene_names : Array[String]
  spot_ids : Array[String]
  cell_types : Array[String]
  beta : Array[Array[Double]]
  standard_errors : Array[Array[Double]]
  t_statistics : Array[Array[Double]]
  p_values : Array[Array[Double]]
  proportions : Array[Array[Double]]
  cells_per_100 : Array[Array[Double]]
  cell_counts : Array[Array[Double]]
  fitted : Array[Array[Double]]
  residuals : Array[Array[Double]]
  outliers : Array[Array[Bool]]
  covariances : Array[Array[Array[Double]]]
  objectives : Array[Double]
  rmse : Array[Double]
  correlations : Array[Double]
  iterations : Array[Int]
  converged : Array[Bool]
  config : SpatialDeconConfig
} derive(Debug)

///|
pub struct SpatialDeconCellMerge {
  name : String
  members : Array[String]
} derive(Debug)

///|
pub struct SpatialDeconReverseResult {
  gene_names : Array[String]
  spot_ids : Array[String]
  cell_types : Array[String]
  coefficients : Array[Array[Double]]
  fitted : Array[Array[Double]]
  residuals : Array[Array[Double]]
  correlations : Array[Double]
  residual_sds : Array[Double]
  objectives : Array[Double]
  iterations : Array[Int]
  converged : Array[Bool]
} derive(Debug)

///|
pub struct SpatialDeconSpatialExperimentOutput {
  experiment : SpatialExperiment
  result : SpatialDeconResult
}

///|
priv struct SdAlignedData {
  gene_names : Array[String]
  values : Array[Array[Double]]
  background : Array[Array[Double]]
  weights : Array[Array[Double]]
  profile : Array[Array[Double]]
}

///|
priv struct SdDerivatives {
  objective : Double
  gradient : Array[Double]
  hessian : Array[Array[Double]]
  expected_hessian : Array[Array[Double]]
  fitted : Array[Double]
}

///|
priv struct SdSolve {
  beta : Array[Double]
  covariance : Array[Array[Double]]
  fitted : Array[Double]
  objective : Double
  iterations : Int
  converged : Bool
}

///|
fn sd_is_finite(value : Double) -> Bool {
  value == value && value.abs() <= 1.0e300
}

///|
fn sd_copy_matrix(matrix : Array[Array[Double]]) -> Array[Array[Double]] {
  let copied : Array[Array[Double]] = []
  for row in matrix {
    copied.push(row.copy())
  }
  copied
}

///|
fn sd_copy_bool_matrix(matrix : Array[Array[Bool]]) -> Array[Array[Bool]] {
  let copied : Array[Array[Bool]] = []
  for row in matrix {
    copied.push(row.copy())
  }
  copied
}

///|
fn sd_zero_matrix(rows : Int, columns : Int) -> Array[Array[Double]] {
  let matrix : Array[Array[Double]] = []
  for _ in 0..<rows {
    matrix.push(Array::make(columns, 0.0))
  }
  matrix
}

///|
fn sd_bool_matrix(rows : Int, columns : Int) -> Array[Array[Bool]] {
  let matrix : Array[Array[Bool]] = []
  for _ in 0..<rows {
    matrix.push(Array::make(columns, false))
  }
  matrix
}

///|
fn sd_contains_string(values : Array[String], target : String) -> Bool {
  for value in values {
    if value == target {
      return true
    }
  }
  false
}

///|
fn sd_string_index(values : Array[String], target : String) -> Int {
  for index in 0..<values.length() {
    if values[index] == target {
      return index
    }
  }
  -1
}

///|
fn sd_validate_names(
  values : Array[String],
  label : String,
) -> Unit raise SpatialDeconError {
  if values.length() == 0 {
    raise SpatialDeconError(label + " must not be empty")
  }
  for index in 0..<values.length() {
    if values[index].length() == 0 {
      raise SpatialDeconError(label + " must not contain empty names")
    }
    for previous in 0..<index {
      if values[previous] == values[index] {
        raise SpatialDeconError(
          label + " contains duplicate name: " + values[index],
        )
      }
    }
  }
}

///|
fn sd_validate_matrix(
  matrix : Array[Array[Double]],
  rows : Int,
  columns : Int,
  label : String,
  allow_zero : Bool,
) -> Unit raise SpatialDeconError {
  if matrix.length() != rows {
    raise SpatialDeconError(label + " row count does not match gene names")
  }
  for row in matrix {
    if row.length() != columns {
      raise SpatialDeconError(label + " must be a rectangular matrix")
    }
    for value in row {
      if !sd_is_finite(value) {
        raise SpatialDeconError(label + " must contain only finite values")
      }
      if value < 0.0 || (!allow_zero && value <= 0.0) {
        raise SpatialDeconError(
          label + " contains an invalid non-positive value",
        )
      }
    }
  }
}

///|
pub fn SpatialDeconProfile::create(
  gene_names : Array[String],
  cell_types : Array[String],
  values : Array[Array[Double]],
) -> SpatialDeconProfile raise SpatialDeconError {
  sd_validate_names(gene_names, "SpatialDecon profile gene names")
  sd_validate_names(cell_types, "SpatialDecon profile cell types")
  sd_validate_matrix(
    values,
    gene_names.length(),
    cell_types.length(),
    "SpatialDecon profile",
    true,
  )
  for cell_type in 0..<cell_types.length() {
    let mut total = 0.0
    for gene in 0..<gene_names.length() {
      total = total + values[gene][cell_type]
    }
    if total <= 0.0 {
      raise SpatialDeconError(
        "SpatialDecon profile cell type has no positive expression: " +
        cell_types[cell_type],
      )
    }
  }
  SpatialDeconProfile::{
    gene_names: gene_names.copy(),
    cell_types: cell_types.copy(),
    values: sd_copy_matrix(values),
  }
}

///|
pub fn SpatialDeconData::create(
  gene_names : Array[String],
  spot_ids : Array[String],
  values : Array[Array[Double]],
  background? : Array[Array[Double]] = [],
  weights? : Array[Array[Double]] = [],
) -> SpatialDeconData raise SpatialDeconError {
  sd_validate_names(gene_names, "SpatialDecon data gene names")
  sd_validate_names(spot_ids, "SpatialDecon spot identifiers")
  sd_validate_matrix(
    values,
    gene_names.length(),
    spot_ids.length(),
    "SpatialDecon expression",
    true,
  )
  let final_background = if background.length() == 0 {
    sd_zero_matrix(gene_names.length(), spot_ids.length())
  } else {
    sd_validate_matrix(
      background,
      gene_names.length(),
      spot_ids.length(),
      "SpatialDecon background",
      true,
    )
    sd_copy_matrix(background)
  }
  let final_weights = if weights.length() == 0 {
    let output = sd_zero_matrix(gene_names.length(), spot_ids.length())
    for gene in 0..<gene_names.length() {
      for spot in 0..<spot_ids.length() {
        output[gene][spot] = 1.0
      }
    }
    output
  } else {
    sd_validate_matrix(
      weights,
      gene_names.length(),
      spot_ids.length(),
      "SpatialDecon weights",
      false,
    )
    sd_copy_matrix(weights)
  }
  SpatialDeconData::{
    gene_names: gene_names.copy(),
    spot_ids: spot_ids.copy(),
    values: sd_copy_matrix(values),
    background: final_background,
    weights: final_weights,
  }
}

///|
pub fn SpatialDeconData::with_scalar_background(
  gene_names : Array[String],
  spot_ids : Array[String],
  values : Array[Array[Double]],
  background : Double,
  weights? : Array[Array[Double]] = [],
) -> SpatialDeconData raise SpatialDeconError {
  if !sd_is_finite(background) || background < 0.0 {
    raise SpatialDeconError(
      "SpatialDecon scalar background must be finite and non-negative",
    )
  }
  let matrix = sd_zero_matrix(gene_names.length(), spot_ids.length())
  for gene in 0..<gene_names.length() {
    for spot in 0..<spot_ids.length() {
      matrix[gene][spot] = background
    }
  }
  SpatialDeconData::create(
    gene_names,
    spot_ids,
    values,
    background=matrix,
    weights~,
  )
}

///|
pub fn SpatialDeconConfig::default() -> SpatialDeconConfig {
  SpatialDeconConfig::{
    residual_threshold: 3.0,
    lower_threshold: 0.5,
    signal_floor: 1.0e-4,
    rescale_profile: true,
    profile_quantile: 0.99,
    profile_target: 2.0,
    refit_outliers: true,
    max_iterations: 1000,
    line_search_steps: 30,
    tolerance: 1.0e-8,
    ridge: 1.0e-8,
  }
}

///|
pub fn SpatialDeconConfig::create(
  residual_threshold? : Double = 3.0,
  lower_threshold? : Double = 0.5,
  signal_floor? : Double = 1.0e-4,
  rescale_profile? : Bool = true,
  profile_quantile? : Double = 0.99,
  profile_target? : Double = 2.0,
  refit_outliers? : Bool = true,
  max_iterations? : Int = 1000,
  line_search_steps? : Int = 30,
  tolerance? : Double = 1.0e-8,
  ridge? : Double = 1.0e-8,
) -> SpatialDeconConfig raise SpatialDeconError {
  if !sd_is_finite(residual_threshold) || residual_threshold <= 0.0 {
    raise SpatialDeconError(
      "SpatialDecon residual threshold must be finite and positive",
    )
  }
  if !sd_is_finite(lower_threshold) || lower_threshold <= 0.0 {
    raise SpatialDeconError(
      "SpatialDecon lower threshold must be finite and positive",
    )
  }
  if !sd_is_finite(signal_floor) || signal_floor <= 0.0 {
    raise SpatialDeconError(
      "SpatialDecon signal floor must be finite and positive",
    )
  }
  if !sd_is_finite(profile_quantile) ||
    profile_quantile <= 0.0 ||
    profile_quantile > 1.0 {
    raise SpatialDeconError("SpatialDecon profile quantile must be in (0, 1]")
  }
  if !sd_is_finite(profile_target) || profile_target <= 0.0 {
    raise SpatialDeconError(
      "SpatialDecon profile target must be finite and positive",
    )
  }
  if max_iterations <= 0 {
    raise SpatialDeconError("SpatialDecon maximum iterations must be positive")
  }
  if line_search_steps <= 0 {
    raise SpatialDeconError("SpatialDecon line-search steps must be positive")
  }
  if !sd_is_finite(tolerance) || tolerance <= 0.0 || tolerance >= 1.0 {
    raise SpatialDeconError("SpatialDecon tolerance must be in (0, 1)")
  }
  if !sd_is_finite(ridge) || ridge <= 0.0 {
    raise SpatialDeconError("SpatialDecon ridge must be finite and positive")
  }
  SpatialDeconConfig::{
    residual_threshold,
    lower_threshold,
    signal_floor,
    rescale_profile,
    profile_quantile,
    profile_target,
    refit_outliers,
    max_iterations,
    line_search_steps,
    tolerance,
    ridge,
  }
}

///|
pub fn SpatialDeconCellMerge::create(
  name : String,
  members : Array[String],
) -> SpatialDeconCellMerge raise SpatialDeconError {
  if name.length() == 0 {
    raise SpatialDeconError(
      "SpatialDecon merged cell type name must not be empty",
    )
  }
  sd_validate_names(members, "SpatialDecon merge members")
  SpatialDeconCellMerge::{ name, members: members.copy(), }
}

///|
pub fn SpatialDeconProfile::n_genes(self : SpatialDeconProfile) -> Int {
  self.gene_names.length()
}

///|
pub fn SpatialDeconProfile::n_cell_types(self : SpatialDeconProfile) -> Int {
  self.cell_types.length()
}

///|
pub fn SpatialDeconData::n_genes(self : SpatialDeconData) -> Int {
  self.gene_names.length()
}

///|
pub fn SpatialDeconData::n_spots(self : SpatialDeconData) -> Int {
  self.spot_ids.length()
}

///|
fn sd_quantile(values : Array[Double], probability : Double) -> Double {
  if values.length() == 0 {
    return 0.0
  }
  let sorted = values.copy()
  for index in 1..<sorted.length() {
    let value = sorted[index]
    let mut position = index
    while position > 0 && sorted[position - 1] > value {
      sorted[position] = sorted[position - 1]
      position = position - 1
    }
    sorted[position] = value
  }
  if sorted.length() == 1 {
    return sorted[0]
  }
  let location = (sorted.length() - 1).to_double() * probability
  let lower = location.to_int()
  let fraction = location - lower.to_double()
  if lower + 1 >= sorted.length() {
    sorted[sorted.length() - 1]
  } else {
    sorted[lower] * (1.0 - fraction) + sorted[lower + 1] * fraction
  }
}

///|
fn sd_profile_scale(
  profile : Array[Array[Double]],
  config : SpatialDeconConfig,
) -> Array[Array[Double]] raise SpatialDeconError {
  let output = sd_copy_matrix(profile)
  if !config.rescale_profile {
    return output
  }
  let flattened : Array[Double] = []
  for row in profile {
    for value in row {
      flattened.push(value)
    }
  }
  let reference = sd_quantile(flattened, config.profile_quantile)
  if reference <= 0.0 {
    raise SpatialDeconError(
      "SpatialDecon profile quantile is zero and cannot be rescaled",
    )
  }
  let factor = config.profile_target / reference
  for gene in 0..<output.length() {
    for cell_type in 0..<output[gene].length() {
      output[gene][cell_type] = output[gene][cell_type] * factor
    }
  }
  output
}

///|
fn sd_align_data(
  data : SpatialDeconData,
  profile : SpatialDeconProfile,
  config : SpatialDeconConfig,
) -> SdAlignedData raise SpatialDeconError {
  let genes : Array[String] = []
  let values : Array[Array[Double]] = []
  let background : Array[Array[Double]] = []
  let weights : Array[Array[Double]] = []
  let profile_rows : Array[Array[Double]] = []
  for data_index in 0..<data.gene_names.length() {
    let profile_index = sd_string_index(
      profile.gene_names,
      data.gene_names[data_index],
    )
    if profile_index >= 0 {
      genes.push(data.gene_names[data_index])
      values.push(data.values[data_index].copy())
      background.push(data.background[data_index].copy())
      weights.push(data.weights[data_index].copy())
      profile_rows.push(profile.values[profile_index].copy())
    }
  }
  if genes.length() == 0 {
    raise SpatialDeconError(
      "SpatialDecon data and profile have no shared gene names",
    )
  }
  if genes.length() < profile.cell_types.length() {
    raise SpatialDeconError(
      "SpatialDecon requires at least as many shared genes as cell types",
    )
  }
  SdAlignedData::{
    gene_names: genes,
    values,
    background,
    weights,
    profile: sd_profile_scale(profile_rows, config),
  }
}

///|
fn sd_dot(left : Array[Double], right : Array[Double]) -> Double {
  let mut total = 0.0
  for index in 0..<left.length() {
    total = total + left[index] * right[index]
  }
  total
}

///|
fn sd_derivatives(
  y : Array[Double],
  background : Array[Double],
  weights : Array[Double],
  profile : Array[Array[Double]],
  included : Array[Bool],
  beta : Array[Double],
  epsilon : Double,
  signal_floor : Double,
) -> SdDerivatives {
  let cell_types = beta.length()
  let gradient = Array::make(cell_types, 0.0)
  let hessian = sd_zero_matrix(cell_types, cell_types)
  let expected_hessian = sd_zero_matrix(cell_types, cell_types)
  let fitted = Array::make(y.length(), 0.0)
  let mut objective = 0.0
  for gene in 0..<y.length() {
    let signal = sd_dot(profile[gene], beta)
    let mean = background[gene] +
      (if signal < signal_floor { signal_floor } else { signal })
    fitted[gene] = mean
    if included[gene] {
      let observed = if y[gene] < epsilon { epsilon } else { y[gene] }
      let residual = @math.ln(observed) - @math.ln(mean)
      objective = objective + weights[gene] * residual * residual
      let gradient_factor = -weights[gene] * residual / mean
      let hessian_factor = weights[gene] * (residual + 1.0) / (mean * mean)
      let expected_factor = weights[gene] / (mean * mean)
      for left in 0..<cell_types {
        gradient[left] = gradient[left] + profile[gene][left] * gradient_factor
        for right in 0..<cell_types {
          let cross = profile[gene][left] * profile[gene][right]
          hessian[left][right] = hessian[left][right] + cross * hessian_factor
          expected_hessian[left][right] = expected_hessian[left][right] +
            cross * expected_factor
        }
      }
    }
  }
  SdDerivatives::{ objective, gradient, hessian, expected_hessian, fitted, }
}

///|
fn sd_solve_linear(
  matrix : Array[Array[Double]],
  values : Array[Double],
) -> Array[Double]? {
  let size = values.length()
  if size == 0 || matrix.length() != size {
    return None
  }
  let a = sd_copy_matrix(matrix)
  let b = values.copy()
  for column in 0..<size {
    let mut pivot = column
    let mut pivot_value = a[column][column].abs()
    for row in (column + 1)..<size {
      if a[row][column].abs() > pivot_value {
        pivot = row
        pivot_value = a[row][column].abs()
      }
    }
    if pivot_value <= 1.0e-14 || !sd_is_finite(pivot_value) {
      return None
    }
    if pivot != column {
      let temporary_row = a[column]
      a[column] = a[pivot]
      a[pivot] = temporary_row
      let temporary_value = b[column]
      b[column] = b[pivot]
      b[pivot] = temporary_value
    }
    let diagonal = a[column][column]
    for row in (column + 1)..<size {
      let factor = a[row][column] / diagonal
      a[row][column] = 0.0
      for next in (column + 1)..<size {
        a[row][next] = a[row][next] - factor * a[column][next]
      }
      b[row] = b[row] - factor * b[column]
    }
  }
  let solution = Array::make(size, 0.0)
  let mut row = size - 1
  while row >= 0 {
    let mut value = b[row]
    for column in (row + 1)..<size {
      value = value - a[row][column] * solution[column]
    }
    solution[row] = value / a[row][row]
    if !sd_is_finite(solution[row]) {
      return None
    }
    row = row - 1
  }
  Some(solution)
}

///|
fn sd_inverse(matrix : Array[Array[Double]]) -> Array[Array[Double]]? {
  let size = matrix.length()
  if size == 0 {
    return None
  }
  let augmented : Array[Array[Double]] = []
  for row in 0..<size {
    if matrix[row].length() != size {
      return None
    }
    let output = Array::make(size * 2, 0.0)
    for column in 0..<size {
      output[column] = matrix[row][column]
      output[size + column] = if row == column { 1.0 } else { 0.0 }
    }
    augmented.push(output)
  }
  for column in 0..<size {
    let mut pivot = column
    let mut pivot_value = augmented[column][column].abs()
    for row in (column + 1)..<size {
      if augmented[row][column].abs() > pivot_value {
        pivot = row
        pivot_value = augmented[row][column].abs()
      }
    }
    if pivot_value <= 1.0e-14 || !sd_is_finite(pivot_value) {
      return None
    }
    if pivot != column {
      let temporary = augmented[column]
      augmented[column] = augmented[pivot]
      augmented[pivot] = temporary
    }
    let diagonal = augmented[column][column]
    for entry in 0..<(size * 2) {
      augmented[column][entry] = augmented[column][entry] / diagonal
    }
    for row in 0..<size {
      if row != column {
        let factor = augmented[row][column]
        for entry in 0..<(size * 2) {
          augmented[row][entry] = augmented[row][entry] -
            factor * augmented[column][entry]
        }
      }
    }
  }
  let inverse = sd_zero_matrix(size, size)
  for row in 0..<size {
    for column in 0..<size {
      inverse[row][column] = augmented[row][size + column]
      if !sd_is_finite(inverse[row][column]) {
        return None
      }
    }
  }
  Some(inverse)
}

///|
fn sd_covariance(
  derivatives : SdDerivatives,
  ridge : Double,
) -> Array[Array[Double]] {
  let size = derivatives.gradient.length()
  let observed = sd_copy_matrix(derivatives.hessian)
  for index in 0..<size {
    observed[index][index] = observed[index][index] + ridge
  }
  match sd_inverse(observed) {
    Some(inverse) => {
      let mut valid = true
      for index in 0..<size {
        if inverse[index][index] <= 0.0 || !sd_is_finite(inverse[index][index]) {
          valid = false
        }
      }
      if valid {
        return inverse
      }
    }
    None => ()
  }
  let expected = sd_copy_matrix(derivatives.expected_hessian)
  for index in 0..<size {
    expected[index][index] = expected[index][index] + ridge
  }
  match sd_inverse(expected) {
    Some(inverse) => inverse
    None => {
      let fallback = sd_zero_matrix(size, size)
      for index in 0..<size {
        fallback[index][index] = 1.0 /
          (if expected[index][index] > ridge {
            expected[index][index]
          } else {
            ridge
          })
      }
      fallback
    }
  }
}

///|
fn sd_projected_gradient_norm(
  beta : Array[Double],
  gradient : Array[Double],
  tolerance : Double,
) -> Double {
  let mut maximum = 0.0
  for index in 0..<beta.length() {
    let value = if beta[index] <= tolerance && gradient[index] > 0.0 {
      0.0
    } else {
      gradient[index].abs()
    }
    if value > maximum {
      maximum = value
    }
  }
  maximum
}

///|
fn sd_solve_one(
  y : Array[Double],
  background : Array[Double],
  weights : Array[Double],
  profile : Array[Array[Double]],
  included : Array[Bool],
  epsilon : Double,
  config : SpatialDeconConfig,
  initial : Array[Double],
) -> SdSolve {
  let cell_types = if profile.length() == 0 { 0 } else { profile[0].length() }
  let beta = if initial.length() == cell_types {
    initial.copy()
  } else {
    let mut mean_signal = 0.0
    let mut mean_profile = 0.0
    let mut points = 0
    for gene in 0..<y.length() {
      if included[gene] {
        let signal = y[gene] - background[gene]
        mean_signal = mean_signal + (if signal > 0.0 { signal } else { 0.0 })
        for cell_type in 0..<cell_types {
          mean_profile = mean_profile + profile[gene][cell_type]
        }
        points = points + 1
      }
    }
    let denominator = if points > 0 {
      mean_profile / points.to_double()
    } else {
      0.0
    }
    let start = if denominator > 0.0 && points > 0 {
      mean_signal / points.to_double() / (denominator * cell_types.to_double())
    } else {
      epsilon
    }
    Array::make(cell_types, if start > 0.0 { start } else { epsilon })
  }
  let mut iterations = 0
  let mut converged = false
  let mut derivatives = sd_derivatives(
    y,
    background,
    weights,
    profile,
    included,
    beta,
    epsilon,
    config.signal_floor,
  )
  while iterations < config.max_iterations && !converged {
    if sd_projected_gradient_norm(beta, derivatives.gradient, config.tolerance) <=
      config.tolerance {
      converged = true
      break
    }
    let system = sd_copy_matrix(derivatives.hessian)
    for index in 0..<cell_types {
      system[index][index] = system[index][index] + config.ridge
    }
    let negative_gradient = Array::make(cell_types, 0.0)
    for index in 0..<cell_types {
      negative_gradient[index] = -derivatives.gradient[index]
    }
    let direction = match sd_solve_linear(system, negative_gradient) {
      Some(candidate) => candidate
      None => negative_gradient.copy()
    }
    let mut directional_derivative = sd_dot(derivatives.gradient, direction)
    if !sd_is_finite(directional_derivative) || directional_derivative >= 0.0 {
      let mut scale = config.ridge
      for index in 0..<cell_types {
        if derivatives.expected_hessian[index][index] > scale {
          scale = derivatives.expected_hessian[index][index]
        }
      }
      for index in 0..<cell_types {
        direction[index] = -derivatives.gradient[index] / scale
      }
      directional_derivative = sd_dot(derivatives.gradient, direction)
    }
    let candidate = beta.copy()
    let mut step = 1.0
    let mut accepted = false
    let mut candidate_derivatives = derivatives
    for _ in 0..<config.line_search_steps {
      for index in 0..<cell_types {
        let proposed = beta[index] + step * direction[index]
        candidate[index] = if proposed > 0.0 { proposed } else { 0.0 }
      }
      candidate_derivatives = sd_derivatives(
        y,
        background,
        weights,
        profile,
        included,
        candidate,
        epsilon,
        config.signal_floor,
      )
      if candidate_derivatives.objective <
        derivatives.objective -
        config.tolerance * step * directional_derivative.abs() {
        accepted = true
        break
      }
      step = step * 0.5
    }
    if !accepted {
      break
    }
    let mut maximum_change = 0.0
    let mut maximum_value = 0.0
    for index in 0..<cell_types {
      let change = (candidate[index] - beta[index]).abs()
      if change > maximum_change {
        maximum_change = change
      }
      beta[index] = candidate[index]
      if beta[index].abs() > maximum_value {
        maximum_value = beta[index].abs()
      }
    }
    derivatives = candidate_derivatives
    iterations = iterations + 1
    if maximum_change <= config.tolerance * (1.0 + maximum_value) {
      converged = true
    }
  }
  let final_derivatives = sd_derivatives(
    y,
    background,
    weights,
    profile,
    included,
    beta,
    epsilon,
    config.signal_floor,
  )
  SdSolve::{
    beta,
    covariance: sd_covariance(final_derivatives, config.ridge),
    fitted: final_derivatives.fitted,
    objective: final_derivatives.objective,
    iterations,
    converged,
  }
}

///|
fn sd_erf(value : Double) -> Double {
  let sign = if value < 0.0 { -1.0 } else { 1.0 }
  let x = value.abs()
  let t = 1.0 / (1.0 + 0.3275911 * x)
  let polynomial = (
      (((1.061405429 * t - 1.453152027) * t + 1.421413741) * t - 0.284496736) *
      t +
      0.254829592
    ) *
    t
  sign * (1.0 - polynomial * @math.exp(-x * x))
}

///|
fn sd_normal_cdf(value : Double) -> Double {
  0.5 * (1.0 + sd_erf(value / 1.4142135623730951))
}

///|
fn sd_mean(values : Array[Double]) -> Double {
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
fn sd_correlation(left : Array[Double], right : Array[Double]) -> Double {
  if left.length() == 0 || left.length() != right.length() {
    return 0.0
  }
  let left_mean = sd_mean(left)
  let right_mean = sd_mean(right)
  let mut covariance = 0.0
  let mut left_variance = 0.0
  let mut right_variance = 0.0
  for index in 0..<left.length() {
    let left_delta = left[index] - left_mean
    let right_delta = right[index] - right_mean
    covariance = covariance + left_delta * right_delta
    left_variance = left_variance + left_delta * left_delta
    right_variance = right_variance + right_delta * right_delta
  }
  let denominator = (left_variance * right_variance).sqrt()
  if denominator <= 1.0e-14 {
    0.0
  } else {
    covariance / denominator
  }
}

///|
fn sd_sample_sd(values : Array[Double]) -> Double {
  if values.length() <= 1 {
    return 0.0
  }
  let mean = sd_mean(values)
  let mut total = 0.0
  for value in values {
    let delta = value - mean
    total = total + delta * delta
  }
  (total / (values.length() - 1).to_double()).sqrt()
}

///|
fn sd_proportions(beta : Array[Array[Double]]) -> Array[Array[Double]] {
  let cell_types = beta.length()
  let spots = if cell_types == 0 { 0 } else { beta[0].length() }
  let output = sd_zero_matrix(cell_types, spots)
  for spot in 0..<spots {
    let mut total = 0.0
    for cell_type in 0..<cell_types {
      total = total + beta[cell_type][spot]
    }
    if total > 0.0 {
      for cell_type in 0..<cell_types {
        output[cell_type][spot] = beta[cell_type][spot] / total
      }
    }
  }
  output
}

///|
fn sd_cells_per_100(beta : Array[Array[Double]]) -> Array[Array[Double]] {
  let cell_types = beta.length()
  let spots = if cell_types == 0 { 0 } else { beta[0].length() }
  let output = sd_zero_matrix(cell_types, spots)
  let mut maximum_total = 0.0
  for spot in 0..<spots {
    let mut total = 0.0
    for cell_type in 0..<cell_types {
      total = total + beta[cell_type][spot]
    }
    if total > maximum_total {
      maximum_total = total
    }
  }
  if maximum_total > 0.0 {
    for cell_type in 0..<cell_types {
      for spot in 0..<spots {
        output[cell_type][spot] = beta[cell_type][spot] / maximum_total * 100.0
      }
    }
  }
  output
}

///|
fn sd_scaled_cell_counts(
  cells_per_100 : Array[Array[Double]],
  nuclei_counts : Array[Double],
) -> Array[Array[Double]] {
  if nuclei_counts.length() == 0 {
    return []
  }
  let cell_types = cells_per_100.length()
  let spots = if cell_types == 0 { 0 } else { cells_per_100[0].length() }
  let output = sd_zero_matrix(cell_types, spots)
  for cell_type in 0..<cell_types {
    for spot in 0..<spots {
      output[cell_type][spot] = cells_per_100[cell_type][spot] /
        100.0 *
        nuclei_counts[spot]
    }
  }
  output
}

///|
pub fn spatial_decon(
  data : SpatialDeconData,
  profile : SpatialDeconProfile,
  nuclei_counts? : Array[Double] = [],
  config? : SpatialDeconConfig = SpatialDeconConfig::default(),
) -> SpatialDeconResult raise SpatialDeconError {
  if nuclei_counts.length() != 0 &&
    nuclei_counts.length() != data.spot_ids.length() {
    raise SpatialDeconError(
      "SpatialDecon nuclei counts must match the number of spots",
    )
  }
  for count in nuclei_counts {
    if !sd_is_finite(count) || count < 0.0 {
      raise SpatialDeconError(
        "SpatialDecon nuclei counts must be finite and non-negative",
      )
    }
  }
  let aligned = sd_align_data(data, profile, config)
  let genes = aligned.gene_names.length()
  let spots = data.spot_ids.length()
  let cell_types = profile.cell_types.length()
  let mut epsilon = 1.0e300
  for gene in 0..<genes {
    for spot in 0..<spots {
      let value = aligned.values[gene][spot]
      if value > 0.0 && value < epsilon {
        epsilon = value
      }
    }
  }
  if epsilon > 1.0e299 {
    raise SpatialDeconError(
      "SpatialDecon expression must contain at least one positive value",
    )
  }
  let outliers = sd_bool_matrix(genes, spots)
  let initial_solutions : Array[SdSolve] = []
  for spot in 0..<spots {
    let y = Array::make(genes, 0.0)
    let background = Array::make(genes, 0.0)
    let weights = Array::make(genes, 0.0)
    for gene in 0..<genes {
      y[gene] = aligned.values[gene][spot]
      background[gene] = aligned.background[gene][spot]
      weights[gene] = aligned.weights[gene][spot]
    }
    let included = Array::make(genes, true)
    let solution = sd_solve_one(
      y,
      background,
      weights,
      aligned.profile,
      included,
      epsilon,
      config,
      [],
    )
    initial_solutions.push(solution)
    if config.refit_outliers {
      let mut retained = 0
      for gene in 0..<genes {
        let observed = if y[gene] > config.lower_threshold {
          y[gene]
        } else {
          config.lower_threshold
        }
        let fitted = if solution.fitted[gene] > config.lower_threshold {
          solution.fitted[gene]
        } else {
          config.lower_threshold
        }
        let residual = @math.ln(observed) / @math.ln(2.0) -
          @math.ln(fitted) / @math.ln(2.0)
        outliers[gene][spot] = residual.abs() > config.residual_threshold
        if !outliers[gene][spot] {
          retained = retained + 1
        }
      }
      if retained <= cell_types {
        for gene in 0..<genes {
          outliers[gene][spot] = false
        }
      }
    }
  }
  let beta = sd_zero_matrix(cell_types, spots)
  let standard_errors = sd_zero_matrix(cell_types, spots)
  let t_statistics = sd_zero_matrix(cell_types, spots)
  let p_values = sd_zero_matrix(cell_types, spots)
  let fitted = sd_zero_matrix(genes, spots)
  let residuals = sd_zero_matrix(genes, spots)
  let covariances : Array[Array[Array[Double]]] = []
  let objectives = Array::make(spots, 0.0)
  let rmse = Array::make(spots, 0.0)
  let correlations = Array::make(spots, 0.0)
  let iterations = Array::make(spots, 0)
  let converged = Array::make(spots, false)
  for spot in 0..<spots {
    let y = Array::make(genes, 0.0)
    let background = Array::make(genes, 0.0)
    let weights = Array::make(genes, 0.0)
    let included = Array::make(genes, true)
    for gene in 0..<genes {
      y[gene] = aligned.values[gene][spot]
      background[gene] = aligned.background[gene][spot]
      weights[gene] = aligned.weights[gene][spot]
      included[gene] = !outliers[gene][spot]
    }
    let solution = if config.refit_outliers {
      sd_solve_one(
        y,
        background,
        weights,
        aligned.profile,
        included,
        epsilon,
        config,
        initial_solutions[spot].beta,
      )
    } else {
      initial_solutions[spot]
    }
    covariances.push(sd_copy_matrix(solution.covariance))
    objectives[spot] = solution.objective
    iterations[spot] = solution.iterations
    converged[spot] = solution.converged
    for cell_type in 0..<cell_types {
      beta[cell_type][spot] = if solution.beta[cell_type] > 0.0 {
        solution.beta[cell_type]
      } else {
        0.0
      }
      let variance = solution.covariance[cell_type][cell_type]
      let se = if variance > 0.0 { variance.sqrt() } else { 0.0 }
      standard_errors[cell_type][spot] = se
      let statistic = if se > 1.0e-14 {
        beta[cell_type][spot] / se
      } else {
        0.0
      }
      t_statistics[cell_type][spot] = statistic
      let probability = if statistic <= 0.0 {
        1.0
      } else {
        2.0 * (1.0 - sd_normal_cdf(statistic))
      }
      p_values[cell_type][spot] = if probability < 0.0 {
        0.0
      } else if probability > 1.0 {
        1.0
      } else {
        probability
      }
    }
    let observed_for_correlation = Array::make(genes, 0.0)
    let fitted_for_correlation = Array::make(genes, 0.0)
    let mut squared_error = 0.0
    let mut retained = 0
    for gene in 0..<genes {
      fitted[gene][spot] = solution.fitted[gene]
      observed_for_correlation[gene] = y[gene]
      fitted_for_correlation[gene] = solution.fitted[gene]
      if outliers[gene][spot] {
        residuals[gene][spot] = @double.not_a_number
      } else {
        let observed = if y[gene] > config.lower_threshold {
          y[gene]
        } else {
          config.lower_threshold
        }
        let expected = if solution.fitted[gene] > config.lower_threshold {
          solution.fitted[gene]
        } else {
          config.lower_threshold
        }
        let residual = @math.ln(observed) / @math.ln(2.0) -
          @math.ln(expected) / @math.ln(2.0)
        residuals[gene][spot] = residual
        squared_error = squared_error + residual * residual
        retained = retained + 1
      }
    }
    rmse[spot] = if retained > 0 {
      (squared_error / retained.to_double()).sqrt()
    } else {
      0.0
    }
    correlations[spot] = sd_correlation(
      observed_for_correlation, fitted_for_correlation,
    )
  }
  let proportions = sd_proportions(beta)
  let cells_per_100 = sd_cells_per_100(beta)
  SpatialDeconResult::{
    gene_names: aligned.gene_names.copy(),
    spot_ids: data.spot_ids.copy(),
    cell_types: profile.cell_types.copy(),
    beta,
    standard_errors,
    t_statistics,
    p_values,
    proportions,
    cells_per_100,
    cell_counts: sd_scaled_cell_counts(cells_per_100, nuclei_counts),
    fitted,
    residuals,
    outliers,
    covariances,
    objectives,
    rmse,
    correlations,
    iterations,
    converged,
    config,
  }
}

///|
pub fn SpatialDeconResult::n_genes(self : SpatialDeconResult) -> Int {
  self.gene_names.length()
}

///|
pub fn SpatialDeconResult::n_spots(self : SpatialDeconResult) -> Int {
  self.spot_ids.length()
}

///|
pub fn SpatialDeconResult::n_cell_types(self : SpatialDeconResult) -> Int {
  self.cell_types.length()
}

///|
pub fn SpatialDeconResult::abundance_for(
  self : SpatialDeconResult,
  cell_type : String,
  spot_id : String,
) -> Double? {
  let cell_index = sd_string_index(self.cell_types, cell_type)
  let spot_index = sd_string_index(self.spot_ids, spot_id)
  if cell_index < 0 || spot_index < 0 {
    None
  } else {
    Some(self.beta[cell_index][spot_index])
  }
}

///|
pub fn SpatialDeconResult::top_cell_types(
  self : SpatialDeconResult,
  spot_id : String,
  limit? : Int = 5,
) -> Array[SpatialDeconAbundance] raise SpatialDeconError {
  if limit <= 0 {
    raise SpatialDeconError("SpatialDecon ranking limit must be positive")
  }
  let spot = sd_string_index(self.spot_ids, spot_id)
  if spot < 0 {
    raise SpatialDeconError("SpatialDecon spot not found: " + spot_id)
  }
  let values : Array[SpatialDeconAbundance] = []
  for cell_type in 0..<self.cell_types.length() {
    values.push(SpatialDeconAbundance::{
      cell_type: self.cell_types[cell_type],
      spot_id,
      abundance: self.beta[cell_type][spot],
      proportion: self.proportions[cell_type][spot],
      standard_error: self.standard_errors[cell_type][spot],
      p_value: self.p_values[cell_type][spot],
    })
  }
  for index in 1..<values.length() {
    let value = values[index]
    let mut position = index
    while position > 0 &&
          (
            values[position - 1].abundance < value.abundance ||
            (
              values[position - 1].abundance == value.abundance &&
              values[position - 1].cell_type > value.cell_type
            )
          ) {
      values[position] = values[position - 1]
      position = position - 1
    }
    values[position] = value
  }
  let count = if limit < values.length() { limit } else { values.length() }
  let top = []
  for index in 0..<count {
    top.push(values[index])
  }
  top
}

///|
pub fn SpatialDeconResult::summary(self : SpatialDeconResult) -> String {
  let mut converged_count = 0
  let mut outlier_count = 0
  for value in self.converged {
    if value {
      converged_count = converged_count + 1
    }
  }
  for row in self.outliers {
    for value in row {
      if value {
        outlier_count = outlier_count + 1
      }
    }
  }
  "SpatialDeconResult(genes=" +
  self.n_genes().to_string() +
  ", spots=" +
  self.n_spots().to_string() +
  ", cell_types=" +
  self.n_cell_types().to_string() +
  ", converged=" +
  converged_count.to_string() +
  "/" +
  self.n_spots().to_string() +
  ", outliers=" +
  outlier_count.to_string() +
  ")"
}

///|
pub fn collapse_spatial_decon(
  result : SpatialDeconResult,
  merges : Array[SpatialDeconCellMerge],
) -> SpatialDeconResult raise SpatialDeconError {
  if merges.length() == 0 {
    raise SpatialDeconError("SpatialDecon cell-type merges must not be empty")
  }
  let merged_names : Array[String] = []
  let used_members : Array[String] = []
  for merge in merges {
    if sd_contains_string(merged_names, merge.name) {
      raise SpatialDeconError(
        "SpatialDecon duplicate merged cell type: " + merge.name,
      )
    }
    merged_names.push(merge.name)
    for source_type in merge.members {
      if !sd_contains_string(result.cell_types, source_type) {
        raise SpatialDeconError(
          "SpatialDecon merge member not found: " + source_type,
        )
      }
      if sd_contains_string(used_members, source_type) {
        raise SpatialDeconError(
          "SpatialDecon merge member appears more than once: " + source_type,
        )
      }
      used_members.push(source_type)
    }
  }
  let groups = merges.length()
  let spots = result.spot_ids.length()
  let beta = sd_zero_matrix(groups, spots)
  let covariances : Array[Array[Array[Double]]] = []
  for spot in 0..<spots {
    let covariance = sd_zero_matrix(groups, groups)
    for left in 0..<groups {
      for left_member in merges[left].members {
        let left_index = sd_string_index(result.cell_types, left_member)
        beta[left][spot] = beta[left][spot] + result.beta[left_index][spot]
        for right in 0..<groups {
          for right_member in merges[right].members {
            let right_index = sd_string_index(result.cell_types, right_member)
            covariance[left][right] = covariance[left][right] +
              result.covariances[spot][left_index][right_index]
          }
        }
      }
    }
    covariances.push(covariance)
  }
  let standard_errors = sd_zero_matrix(groups, spots)
  let t_statistics = sd_zero_matrix(groups, spots)
  let p_values = sd_zero_matrix(groups, spots)
  for group in 0..<groups {
    for spot in 0..<spots {
      let variance = covariances[spot][group][group]
      let se = if variance > 0.0 { variance.sqrt() } else { 0.0 }
      standard_errors[group][spot] = se
      let statistic = if se > 1.0e-14 { beta[group][spot] / se } else { 0.0 }
      t_statistics[group][spot] = statistic
      p_values[group][spot] = if statistic <= 0.0 {
        1.0
      } else {
        let probability = 2.0 * (1.0 - sd_normal_cdf(statistic))
        if probability < 0.0 {
          0.0
        } else {
          probability
        }
      }
    }
  }
  let proportions = sd_proportions(beta)
  let cells_per_100 = sd_cells_per_100(beta)
  let cell_counts = if result.cell_counts.length() == 0 {
    []
  } else {
    let output = sd_zero_matrix(groups, spots)
    for group in 0..<groups {
      for source_type in merges[group].members {
        let member_index = sd_string_index(result.cell_types, source_type)
        for spot in 0..<spots {
          output[group][spot] = output[group][spot] +
            result.cell_counts[member_index][spot]
        }
      }
    }
    output
  }
  SpatialDeconResult::{
    gene_names: result.gene_names.copy(),
    spot_ids: result.spot_ids.copy(),
    cell_types: merged_names,
    beta,
    standard_errors,
    t_statistics,
    p_values,
    proportions,
    cells_per_100,
    cell_counts,
    fitted: sd_copy_matrix(result.fitted),
    residuals: sd_copy_matrix(result.residuals),
    outliers: sd_copy_bool_matrix(result.outliers),
    covariances,
    objectives: result.objectives.copy(),
    rmse: result.rmse.copy(),
    correlations: result.correlations.copy(),
    iterations: result.iterations.copy(),
    converged: result.converged.copy(),
    config: result.config,
  }
}

///|
pub fn spatial_decon_background(
  gene_names : Array[String],
  spot_ids : Array[String],
  values : Array[Array[Double]],
  probe_pools : Array[String],
  negative_genes : Array[String],
) -> Array[Array[Double]] raise SpatialDeconError {
  sd_validate_names(gene_names, "SpatialDecon background gene names")
  sd_validate_names(spot_ids, "SpatialDecon background spot identifiers")
  sd_validate_names(negative_genes, "SpatialDecon negative-probe names")
  sd_validate_matrix(
    values,
    gene_names.length(),
    spot_ids.length(),
    "SpatialDecon background expression",
    true,
  )
  if probe_pools.length() != gene_names.length() {
    raise SpatialDeconError(
      "SpatialDecon probe pools must match the number of genes",
    )
  }
  for pool in probe_pools {
    if pool.length() == 0 {
      raise SpatialDeconError("SpatialDecon probe pool names must not be empty")
    }
  }
  let pools : Array[String] = []
  for pool in probe_pools {
    if !sd_contains_string(pools, pool) {
      pools.push(pool)
    }
  }
  let output = sd_zero_matrix(gene_names.length(), spot_ids.length())
  for pool in pools {
    let negatives : Array[Int] = []
    for gene in 0..<gene_names.length() {
      if probe_pools[gene] == pool &&
        sd_contains_string(negative_genes, gene_names[gene]) {
        negatives.push(gene)
      }
    }
    if negatives.length() == 0 {
      raise SpatialDeconError(
        "SpatialDecon probe pool has no specified negative probe: " + pool,
      )
    }
    let means = Array::make(spot_ids.length(), 0.0)
    for spot in 0..<spot_ids.length() {
      for gene in negatives {
        means[spot] = means[spot] + values[gene][spot]
      }
      means[spot] = means[spot] / negatives.length().to_double()
    }
    for gene in 0..<gene_names.length() {
      if probe_pools[gene] == pool {
        for spot in 0..<spot_ids.length() {
          output[gene][spot] = means[spot]
        }
      }
    }
  }
  output
}

///|
pub fn create_spatial_decon_profile(
  gene_names : Array[String],
  cell_ids : Array[String],
  cell_types : Array[String],
  counts : Array[Array[Double]],
  normalize? : Bool = false,
  scaling_factor? : Double = 5.0,
  min_cells? : Int = 15,
  min_genes? : Int = 100,
  gene_filter? : Array[String] = [],
) -> SpatialDeconProfile raise SpatialDeconError {
  sd_validate_names(gene_names, "SpatialDecon single-cell gene names")
  sd_validate_names(cell_ids, "SpatialDecon single-cell identifiers")
  if cell_types.length() != cell_ids.length() {
    raise SpatialDeconError(
      "SpatialDecon cell types must match single-cell identifiers",
    )
  }
  for cell_type in cell_types {
    if cell_type.length() == 0 {
      raise SpatialDeconError(
        "SpatialDecon single-cell types must not be empty",
      )
    }
  }
  sd_validate_matrix(
    counts,
    gene_names.length(),
    cell_ids.length(),
    "SpatialDecon single-cell counts",
    true,
  )
  if !sd_is_finite(scaling_factor) || scaling_factor <= 0.0 {
    raise SpatialDeconError(
      "SpatialDecon profile scaling factor must be finite and positive",
    )
  }
  if min_cells < 0 || min_genes < 0 {
    raise SpatialDeconError(
      "SpatialDecon profile filtering thresholds must be non-negative",
    )
  }
  if gene_filter.length() > 0 {
    sd_validate_names(gene_filter, "SpatialDecon profile gene filter")
  }
  let libraries = Array::make(cell_ids.length(), 0.0)
  let expressed = Array::make(cell_ids.length(), 0)
  for cell in 0..<cell_ids.length() {
    for gene in 0..<gene_names.length() {
      libraries[cell] = libraries[cell] + counts[gene][cell]
      if counts[gene][cell] > 0.0 {
        expressed[cell] = expressed[cell] + 1
      }
    }
  }
  let positive_libraries : Array[Double] = []
  for library in libraries {
    if library > 0.0 {
      positive_libraries.push(library)
    }
  }
  if positive_libraries.length() == 0 {
    raise SpatialDeconError(
      "SpatialDecon single-cell counts contain no positive library",
    )
  }
  let median_library = sd_quantile(positive_libraries, 0.5)
  let ordered_types : Array[String] = []
  for cell_type in cell_types {
    if !sd_contains_string(ordered_types, cell_type) {
      ordered_types.push(cell_type)
    }
  }
  let retained_types : Array[String] = []
  let retained_cells : Array[Array[Int]] = []
  for cell_type in ordered_types {
    let indices : Array[Int] = []
    for cell in 0..<cell_ids.length() {
      if cell_types[cell] == cell_type &&
        expressed[cell] > min_genes &&
        libraries[cell] > 0.0 {
        indices.push(cell)
      }
    }
    if indices.length() > min_cells {
      retained_types.push(cell_type)
      retained_cells.push(indices)
    }
  }
  if retained_types.length() == 0 {
    raise SpatialDeconError(
      "SpatialDecon profile filtering removed every cell type",
    )
  }
  let selected_genes : Array[Int] = []
  for gene in 0..<gene_names.length() {
    if gene_filter.length() == 0 ||
      sd_contains_string(gene_filter, gene_names[gene]) {
      let mut positive = false
      for group in retained_cells {
        for cell in group {
          if counts[gene][cell] > 0.0 {
            positive = true
          }
        }
      }
      if positive {
        selected_genes.push(gene)
      }
    }
  }
  if selected_genes.length() == 0 {
    raise SpatialDeconError(
      "SpatialDecon profile gene filtering removed every gene",
    )
  }
  let output_genes : Array[String] = []
  let values : Array[Array[Double]] = []
  for gene in selected_genes {
    output_genes.push(gene_names[gene])
    let row = Array::make(retained_types.length(), 0.0)
    for cell_type in 0..<retained_types.length() {
      for cell in retained_cells[cell_type] {
        let normalized = if normalize {
          counts[gene][cell] * median_library / libraries[cell]
        } else {
          counts[gene][cell]
        }
        row[cell_type] = row[cell_type] + normalized
      }
      row[cell_type] = row[cell_type] /
        retained_cells[cell_type].length().to_double() *
        scaling_factor
    }
    values.push(row)
  }
  SpatialDeconProfile::create(output_genes, retained_types, values)
}

///|
pub fn reverse_spatial_decon(
  data : SpatialDeconData,
  result : SpatialDeconResult,
  epsilon? : Double = 0.0,
) -> SpatialDeconReverseResult raise SpatialDeconError {
  if data.spot_ids != result.spot_ids {
    raise SpatialDeconError(
      "SpatialDecon reverse analysis requires matching spot identifiers",
    )
  }
  if epsilon < 0.0 || !sd_is_finite(epsilon) {
    raise SpatialDeconError(
      "SpatialDecon reverse epsilon must be finite and non-negative",
    )
  }
  let active_types : Array[String] = []
  let active_indices : Array[Int] = []
  for cell_type in 0..<result.cell_types.length() {
    if sd_sample_sd(result.beta[cell_type]) > 1.0e-14 {
      active_types.push(result.cell_types[cell_type])
      active_indices.push(cell_type)
    }
  }
  if active_types.length() == 0 {
    raise SpatialDeconError(
      "SpatialDecon reverse analysis requires a varying cell score",
    )
  }
  let positive : Array[Double] = []
  for row in data.values {
    for value in row {
      if value > 0.0 {
        positive.push(value)
      }
    }
  }
  if positive.length() == 0 {
    raise SpatialDeconError(
      "SpatialDecon reverse expression must contain positive values",
    )
  }
  let effective_epsilon = if epsilon > 0.0 {
    epsilon
  } else {
    sd_quantile(positive, 0.01)
  }
  let spots = data.spot_ids.length()
  let design : Array[Array[Double]] = []
  for spot in 0..<spots {
    let row : Array[Double] = [1.0]
    for index in active_indices {
      row.push(result.beta[index][spot])
    }
    design.push(row)
  }
  let genes = data.gene_names.length()
  let coefficients = sd_zero_matrix(genes, active_types.length() + 1)
  let fitted = sd_zero_matrix(genes, spots)
  let residuals = sd_zero_matrix(genes, spots)
  let correlations = Array::make(genes, 0.0)
  let residual_sds = Array::make(genes, 0.0)
  let objectives = Array::make(genes, 0.0)
  let iterations = Array::make(genes, 0)
  let converged = Array::make(genes, false)
  let background = Array::make(spots, 0.0)
  let weights = Array::make(spots, 1.0)
  let included = Array::make(spots, true)
  for gene in 0..<genes {
    let solution = sd_solve_one(
      data.values[gene],
      background,
      weights,
      design,
      included,
      effective_epsilon,
      result.config,
      [],
    )
    coefficients[gene] = solution.beta.copy()
    fitted[gene] = solution.fitted.copy()
    objectives[gene] = solution.objective
    iterations[gene] = solution.iterations
    converged[gene] = solution.converged
    let gene_residuals = Array::make(spots, 0.0)
    for spot in 0..<spots {
      let observed = if data.values[gene][spot] > effective_epsilon {
        data.values[gene][spot]
      } else {
        effective_epsilon
      }
      let expected = if fitted[gene][spot] > effective_epsilon {
        fitted[gene][spot]
      } else {
        effective_epsilon
      }
      gene_residuals[spot] = @math.ln(observed) / @math.ln(2.0) -
        @math.ln(expected) / @math.ln(2.0)
      residuals[gene][spot] = gene_residuals[spot]
    }
    correlations[gene] = sd_correlation(data.values[gene], fitted[gene])
    residual_sds[gene] = sd_sample_sd(gene_residuals)
  }
  SpatialDeconReverseResult::{
    gene_names: data.gene_names.copy(),
    spot_ids: data.spot_ids.copy(),
    cell_types: active_types,
    coefficients,
    fitted,
    residuals,
    correlations,
    residual_sds,
    objectives,
    iterations,
    converged,
  }
}

///|
fn sd_copy_string_map(source : Map[String, String]) -> Map[String, String] {
  let output : Map[String, String] = Map([], capacity=source.length())
  for key in source.keys() {
    output[key] = source[key]
  }
  output
}

///|
fn sd_copy_spatial_experiment(
  experiment : SpatialExperiment,
) -> SpatialExperiment {
  let output = SpatialExperiment::new()
  for key in experiment.assay.keys() {
    output.assay[key] = sd_copy_matrix(experiment.assay[key])
  }
  for row in experiment.row_data {
    output.row_data.push(sd_copy_string_map(row))
  }
  for column in experiment.col_data {
    output.col_data.push(sd_copy_string_map(column))
  }
  for coordinate in experiment.spatial_coords {
    output.spatial_coords.push(coordinate)
  }
  for image in experiment.images {
    output.images.push(
      SpatialImage::new(
        image.id,
        sd_copy_matrix(image.data),
        image.scale_factor,
      ),
    )
  }
  for key in experiment.metadata.keys() {
    output.metadata[key] = experiment.metadata[key]
  }
  output
}

///|
fn sd_experiment_gene_names(
  experiment : SpatialExperiment,
  rows : Int,
  key : String,
) -> Array[String] raise SpatialDeconError {
  if experiment.row_data.length() != 0 && experiment.row_data.length() != rows {
    raise SpatialDeconError(
      "SpatialDecon SpatialExperiment row data must match assay rows",
    )
  }
  let names : Array[String] = []
  for row in 0..<rows {
    let generated = "gene_" + (row + 1).to_string()
    if key.length() == 0 || experiment.row_data.length() == 0 {
      names.push(generated)
    } else {
      match experiment.row_data[row].get(key) {
        Some(value) =>
          if value.length() == 0 {
            raise SpatialDeconError(
              "SpatialDecon SpatialExperiment gene identifier is empty",
            )
          } else {
            names.push(value)
          }
        None =>
          raise SpatialDeconError(
            "SpatialDecon SpatialExperiment gene field not found: " + key,
          )
      }
    }
  }
  names
}

///|
fn sd_experiment_spot_ids(
  experiment : SpatialExperiment,
  columns : Int,
  key : String,
) -> Array[String] raise SpatialDeconError {
  if experiment.col_data.length() != 0 &&
    experiment.col_data.length() != columns {
    raise SpatialDeconError(
      "SpatialDecon SpatialExperiment column data must match assay columns",
    )
  }
  let names : Array[String] = []
  for column in 0..<columns {
    let generated = "spot_" + (column + 1).to_string()
    if key.length() == 0 || experiment.col_data.length() == 0 {
      names.push(generated)
    } else {
      match experiment.col_data[column].get(key) {
        Some(value) =>
          if value.length() == 0 {
            raise SpatialDeconError(
              "SpatialDecon SpatialExperiment spot identifier is empty",
            )
          } else {
            names.push(value)
          }
        None =>
          raise SpatialDeconError(
            "SpatialDecon SpatialExperiment spot field not found: " + key,
          )
      }
    }
  }
  names
}

///|
pub fn spatial_decon_spatial_experiment(
  experiment : SpatialExperiment,
  profile : SpatialDeconProfile,
  assay_name? : String = "counts",
  background_assay_name? : String = "",
  weight_assay_name? : String = "",
  gene_id_key? : String = "gene_id",
  spot_id_key? : String = "spot_id",
  abundance_prefix? : String = "SpatialDecon:",
  proportion_prefix? : String = "SpatialDeconProp:",
  scalar_background? : Double = 0.0,
  nuclei_counts? : Array[Double] = [],
  config? : SpatialDeconConfig = SpatialDeconConfig::default(),
) -> SpatialDeconSpatialExperimentOutput raise SpatialDeconError {
  if assay_name.length() == 0 ||
    abundance_prefix.length() == 0 ||
    proportion_prefix.length() == 0 {
    raise SpatialDeconError(
      "SpatialDecon SpatialExperiment names and prefixes must not be empty",
    )
  }
  if !sd_is_finite(scalar_background) || scalar_background < 0.0 {
    raise SpatialDeconError(
      "SpatialDecon SpatialExperiment background must be non-negative",
    )
  }
  let values = match experiment.assay.get(assay_name) {
    Some(matrix) => matrix
    None =>
      raise SpatialDeconError(
        "SpatialDecon SpatialExperiment assay not found: " + assay_name,
      )
  }
  if values.length() == 0 || values[0].length() == 0 {
    raise SpatialDeconError(
      "SpatialDecon SpatialExperiment assay must not be empty",
    )
  }
  let rows = values.length()
  let columns = values[0].length()
  let gene_names = sd_experiment_gene_names(experiment, rows, gene_id_key)
  let spot_ids = sd_experiment_spot_ids(experiment, columns, spot_id_key)
  let background = if background_assay_name.length() == 0 {
    let output = sd_zero_matrix(rows, columns)
    for row in 0..<rows {
      for column in 0..<columns {
        output[row][column] = scalar_background
      }
    }
    output
  } else {
    match experiment.assay.get(background_assay_name) {
      Some(matrix) => sd_copy_matrix(matrix)
      None =>
        raise SpatialDeconError(
          "SpatialDecon background assay not found: " + background_assay_name,
        )
    }
  }
  let weights = if weight_assay_name.length() == 0 {
    []
  } else {
    match experiment.assay.get(weight_assay_name) {
      Some(matrix) => sd_copy_matrix(matrix)
      None =>
        raise SpatialDeconError(
          "SpatialDecon weight assay not found: " + weight_assay_name,
        )
    }
  }
  let data = SpatialDeconData::create(
    gene_names,
    spot_ids,
    values,
    background~,
    weights~,
  )
  let result = spatial_decon(data, profile, nuclei_counts~, config~)
  let enriched = sd_copy_spatial_experiment(experiment)
  if enriched.col_data.length() == 0 {
    for _ in 0..<columns {
      enriched.col_data.push(Map([], capacity=result.cell_types.length() * 2))
    }
  }
  for spot in 0..<columns {
    for cell_type in 0..<result.cell_types.length() {
      enriched.col_data[spot][abundance_prefix + result.cell_types[cell_type]] = result.beta[cell_type][spot].to_string()
      enriched.col_data[spot][proportion_prefix + result.cell_types[cell_type]] = result.proportions[cell_type][spot].to_string()
    }
  }
  enriched.metadata["spatialdecon_genes"] = result.n_genes().to_string()
  enriched.metadata["spatialdecon_spots"] = result.n_spots().to_string()
  enriched.metadata["spatialdecon_cell_types"] = result
    .n_cell_types()
    .to_string()
  let mut converged_count = 0
  for value in result.converged {
    if value {
      converged_count = converged_count + 1
    }
  }
  enriched.metadata["spatialdecon_converged"] = converged_count.to_string()
  SpatialDeconSpatialExperimentOutput::{ experiment: enriched, result, }
}

///|
pub fn spatial_decon_example_data() -> (
  SpatialDeconData,
  SpatialDeconProfile,
  Array[Double],
) raise SpatialDeconError {
  let genes = ["T1", "T2", "T3", "B1", "B2", "B3", "M1", "M2", "M3"]
  let spots = ["region_1", "region_2", "region_3", "region_4", "region_5"]
  let cell_types = ["T_cell", "B_cell", "Myeloid"]
  let profile_values = [
    [8.0, 0.5, 0.3],
    [6.0, 0.6, 0.4],
    [5.0, 0.4, 0.5],
    [0.4, 8.0, 0.5],
    [0.5, 6.0, 0.4],
    [0.3, 5.0, 0.6],
    [0.5, 0.4, 8.0],
    [0.4, 0.6, 6.0],
    [0.6, 0.5, 5.0],
  ]
  let true_beta = [
    [4.0, 1.0, 2.0, 0.5, 3.0],
    [1.0, 4.0, 2.0, 1.0, 0.5],
    [0.5, 1.0, 2.0, 4.0, 3.0],
  ]
  let background = sd_zero_matrix(genes.length(), spots.length())
  let values = sd_zero_matrix(genes.length(), spots.length())
  for gene in 0..<genes.length() {
    for spot in 0..<spots.length() {
      background[gene][spot] = 0.5 + spot.to_double() * 0.05
      let mut signal = 0.0
      for cell_type in 0..<cell_types.length() {
        signal = signal +
          profile_values[gene][cell_type] * true_beta[cell_type][spot]
      }
      let adjustment = 1.0 + ((gene * 3 + spot * 2) % 5 - 2).to_double() * 0.01
      values[gene][spot] = background[gene][spot] + signal * adjustment
    }
  }
  let profile = SpatialDeconProfile::create(genes, cell_types, profile_values)
  let data = SpatialDeconData::create(genes, spots, values, background~)
  (data, profile, [100.0, 120.0, 110.0, 90.0, 105.0])
}
