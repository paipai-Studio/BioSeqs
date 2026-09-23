// Approximate posterior estimation for negative-binomial GLM coefficients.
//
// This module follows Bioconductor apeglm's nbinom workflow. Coefficients and
// prior scales use the natural-log scale; convenience output methods can
// convert effects to log2 scale.

///|
pub suberror ApeglmError {
  ApeglmError(String)
}

///|
pub struct ApeglmConfig {
  coefficient : Int
  threshold : Double
  interval_level : Double
  prior_scale : Double
  prior_df : Double
  prior_no_shrink_scale : Double
  adaptive_prior : Bool
  multiplier : Double
  max_iterations : Int
  tolerance : Double
  max_step : Double
  random_starts : Int
} derive(Eq, Debug)
pub extend ApeglmConfig with Eq::{equal, not_equal}
pub extend ApeglmConfig with @moonbitlang/core/debug.Debug::{to_repr}

///|
pub fn ApeglmConfig::create(
  coefficient? : Int = -1,
  threshold? : Double = 0.0,
  interval_level? : Double = 0.95,
  prior_scale? : Double = 1.0,
  prior_df? : Double = 1.0,
  prior_no_shrink_scale? : Double = 15.0,
  adaptive_prior? : Bool = true,
  multiplier? : Double = 1.0,
  max_iterations? : Int = 100,
  tolerance? : Double = 1.0e-7,
  max_step? : Double = 2.0,
  random_starts? : Int = 5,
) -> ApeglmConfig raise ApeglmError {
  if coefficient < -1 {
    raise ApeglmError("apeglm coefficient must be -1 or a zero-based index")
  }
  if !apeglm_is_finite(threshold) || threshold < 0.0 {
    raise ApeglmError("apeglm threshold must be finite and non-negative")
  }
  if !apeglm_is_finite(interval_level) ||
    interval_level <= 0.0 ||
    interval_level >= 1.0 {
    raise ApeglmError("apeglm interval level must be finite and in (0, 1)")
  }
  if !apeglm_is_finite(prior_scale) || prior_scale <= 0.0 {
    raise ApeglmError("apeglm prior scale must be finite and positive")
  }
  if !apeglm_is_finite(prior_df) || prior_df <= 0.0 {
    raise ApeglmError(
      "apeglm prior degrees of freedom must be finite and positive",
    )
  }
  if !apeglm_is_finite(prior_no_shrink_scale) || prior_no_shrink_scale <= 0.0 {
    raise ApeglmError(
      "apeglm no-shrink prior scale must be finite and positive",
    )
  }
  if !apeglm_is_finite(multiplier) || multiplier <= 0.0 {
    raise ApeglmError("apeglm prior multiplier must be finite and positive")
  }
  if max_iterations < 1 {
    raise ApeglmError("apeglm maximum iterations must be positive")
  }
  if !apeglm_is_finite(tolerance) || tolerance <= 0.0 {
    raise ApeglmError("apeglm tolerance must be finite and positive")
  }
  if !apeglm_is_finite(max_step) || max_step <= 0.0 {
    raise ApeglmError("apeglm maximum Newton step must be finite and positive")
  }
  if random_starts < 1 {
    raise ApeglmError("apeglm start count must be positive")
  }
  ApeglmConfig::{
    coefficient,
    threshold,
    interval_level,
    prior_scale,
    prior_df,
    prior_no_shrink_scale,
    adaptive_prior,
    multiplier,
    max_iterations,
    tolerance,
    max_step,
    random_starts,
  }
}

///|
pub fn ApeglmConfig::default() -> ApeglmConfig {
  ApeglmConfig::{
    coefficient: -1,
    threshold: 0.0,
    interval_level: 0.95,
    prior_scale: 1.0,
    prior_df: 1.0,
    prior_no_shrink_scale: 15.0,
    adaptive_prior: true,
    multiplier: 1.0,
    max_iterations: 100,
    tolerance: 1.0e-7,
    max_step: 2.0,
    random_starts: 5,
  }
}

///|
pub struct ApeglmPriorControl {
  no_shrink : Array[Int]
  prior_mean : Double
  prior_scale : Double
  prior_variance : Double
  prior_df : Double
  prior_no_shrink_scale : Double
  adaptive : Bool
} derive(Debug)
pub extend ApeglmPriorControl with @moonbitlang/core/debug.Debug::{to_repr}

///|
pub struct ApeglmGeneResult {
  gene_name : String
  base_mean : Double
  mle : Array[Double]
  mle_sd : Array[Double]
  map : Array[Double]
  posterior_sd : Array[Double]
  fsr : Double
  s_value : Double
  threshold_probability : Double
  interval_lower : Double
  interval_upper : Double
  log_posterior : Double
  converged : Bool
  iterations : Int
} derive(Debug)
pub extend ApeglmGeneResult with @moonbitlang/core/debug.Debug::{to_repr}

///|
pub struct ApeglmResult {
  gene_names : Array[String]
  coefficient_names : Array[String]
  coefficient : Int
  genes : Array[ApeglmGeneResult]
  prior_control : ApeglmPriorControl
  threshold : Double
  interval_level : Double
  config : ApeglmConfig
} derive(Debug)
pub extend ApeglmResult with @moonbitlang/core/debug.Debug::{to_repr}

///|
pub struct ApeglmSummary {
  gene_count : Int
  converged_count : Int
  selected_count : Int
  low_fsr_count : Int
  prior_scale : Double
  median_absolute_effect : Double
  coefficient_name : String
} derive(Eq, Debug)
pub extend ApeglmSummary with Eq::{equal, not_equal}
pub extend ApeglmSummary with @moonbitlang/core/debug.Debug::{to_repr}

///|
pub struct ApeglmSummarizedExperimentOutput {
  experiment : SummarizedExperiment
  result : ApeglmResult
}

///|
priv struct ApeglmEvaluation {
  value : Double
  gradient : Array[Double]
  information : Array[Array[Double]]
}

///|
priv struct ApeglmOptimization {
  beta : Array[Double]
  covariance : Array[Array[Double]]
  value : Double
  converged : Bool
  iterations : Int
}

///|
fn apeglm_is_finite(value : Double) -> Bool {
  value == value && value.abs() <= 1.0e300
}

///|
fn apeglm_clamp(value : Double, lower : Double, upper : Double) -> Double {
  value.max(lower).min(upper)
}

///|
fn apeglm_zero_matrix(rows : Int, columns : Int) -> Array[Array[Double]] {
  let output : Array[Array[Double]] = []
  for _ in 0..<rows {
    output.push(Array::make(columns, 0.0))
  }
  output
}

///|
fn apeglm_copy_matrix(matrix : Array[Array[Double]]) -> Array[Array[Double]] {
  let output : Array[Array[Double]] = []
  for row in matrix {
    output.push(row.copy())
  }
  output
}

///|
fn apeglm_mean(values : Array[Double]) -> Double {
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
fn apeglm_normal_cdf(value : Double) -> Double {
  let absolute = value.abs()
  let t = 1.0 / (1.0 + 0.2316419 * absolute)
  let density = 0.3989422804014327 * @math.exp(-0.5 * absolute * absolute)
  let tail = density *
    t *
    (
      0.319381530 +
      t *
      (-0.356563782 + t * (1.781477937 + t * (-1.821255978 + t * 1.330274429)))
    )
  if value >= 0.0 {
    apeglm_clamp(1.0 - tail, 0.0, 1.0)
  } else {
    apeglm_clamp(tail, 0.0, 1.0)
  }
}

///|
fn apeglm_normal_quantile(probability : Double) -> Double {
  let target = apeglm_clamp(probability, 1.0e-12, 1.0 - 1.0e-12)
  let mut lower = -8.0
  let mut upper = 8.0
  for _ in 0..<80 {
    let middle = 0.5 * (lower + upper)
    if apeglm_normal_cdf(middle) < target {
      lower = middle
    } else {
      upper = middle
    }
  }
  0.5 * (lower + upper)
}

///|
fn apeglm_cholesky(matrix : Array[Array[Double]]) -> Array[Array[Double]]? {
  let size = matrix.length()
  if size == 0 {
    return None
  }
  let lower = apeglm_zero_matrix(size, size)
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
        if value <= 1.0e-12 || !apeglm_is_finite(value) {
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
fn apeglm_regularized_cholesky(
  matrix : Array[Array[Double]],
) -> Array[Array[Double]]? {
  let mut ridge = 0.0
  for attempt in 0..<14 {
    let candidate = apeglm_copy_matrix(matrix)
    for index in 0..<candidate.length() {
      candidate[index][index] = candidate[index][index] + ridge
    }
    match apeglm_cholesky(candidate) {
      Some(lower) => return Some(lower)
      None => ridge = if attempt == 0 { 1.0e-10 } else { ridge * 10.0 }
    }
  }
  None
}

///|
fn apeglm_cholesky_solve(
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
fn apeglm_cholesky_inverse(
  lower : Array[Array[Double]],
) -> Array[Array[Double]] {
  let size = lower.length()
  let inverse = apeglm_zero_matrix(size, size)
  for column in 0..<size {
    let basis = Array::make(size, 0.0)
    basis[column] = 1.0
    let solution = apeglm_cholesky_solve(lower, basis)
    for row in 0..<size {
      inverse[row][column] = solution[row]
    }
  }
  inverse
}

///|
fn apeglm_solve_information(
  information : Array[Array[Double]],
  gradient : Array[Double],
) -> Array[Double] {
  match apeglm_regularized_cholesky(information) {
    Some(lower) => apeglm_cholesky_solve(lower, gradient)
    None => {
      let output = Array::make(gradient.length(), 0.0)
      for index in 0..<gradient.length() {
        output[index] = gradient[index] /
          (information[index][index].abs() + 1.0e-6)
      }
      output
    }
  }
}

///|
fn apeglm_covariance(
  information : Array[Array[Double]],
) -> Array[Array[Double]] {
  match apeglm_regularized_cholesky(information) {
    Some(lower) => apeglm_cholesky_inverse(lower)
    None => {
      let output = apeglm_zero_matrix(
        information.length(),
        information.length(),
      )
      for index in 0..<information.length() {
        output[index][index] = 1.0 / (information[index][index].abs() + 1.0e-6)
      }
      output
    }
  }
}

///|
fn apeglm_validate_matrix(
  matrix : Array[Array[Double]],
  rows : Int,
  columns : Int,
  label : String,
  non_negative : Bool,
) -> Unit raise ApeglmError {
  if matrix.length() != rows {
    raise ApeglmError(label + " row count does not match the expected shape")
  }
  for row in matrix {
    if row.length() != columns {
      raise ApeglmError(label + " must be rectangular with the expected shape")
    }
    for value in row {
      if !apeglm_is_finite(value) || (non_negative && value < 0.0) {
        raise ApeglmError(label + " contains an invalid value")
      }
    }
  }
}

///|
fn apeglm_validate_inputs(
  counts : Array[Array[Double]],
  design : Array[Array[Double]],
  dispersions : Array[Double],
  offsets : Array[Array[Double]],
  weights : Array[Array[Double]],
  gene_names : Array[String],
  coefficient_names : Array[String],
  config : ApeglmConfig,
) -> (Int, Int, Int, Int) raise ApeglmError {
  if counts.length() == 0 {
    raise ApeglmError("apeglm counts must contain at least one feature")
  }
  let samples = counts[0].length()
  if samples < 2 {
    raise ApeglmError("apeglm counts must contain at least two samples")
  }
  apeglm_validate_matrix(
    counts,
    counts.length(),
    samples,
    "apeglm counts",
    true,
  )
  for row in counts {
    for value in row {
      if (value - value.round()).abs() > 1.0e-8 {
        raise ApeglmError("apeglm requires integer-valued counts")
      }
    }
  }
  if design.length() != samples {
    raise ApeglmError("apeglm design rows must match the sample count")
  }
  if design.length() == 0 || design[0].length() < 2 {
    raise ApeglmError(
      "apeglm negative-binomial design requires an intercept and a coefficient",
    )
  }
  let coefficients = design[0].length()
  apeglm_validate_matrix(design, samples, coefficients, "apeglm design", false)
  let coefficient = if config.coefficient < 0 {
    coefficients - 1
  } else {
    config.coefficient
  }
  if coefficient <= 0 || coefficient >= coefficients {
    raise ApeglmError(
      "apeglm coefficient must select a non-intercept design column",
    )
  }
  if dispersions.length() != counts.length() {
    raise ApeglmError("apeglm dispersions must match the feature count")
  }
  for dispersion in dispersions {
    if !apeglm_is_finite(dispersion) || dispersion <= 0.0 {
      raise ApeglmError("apeglm dispersions must be finite and positive")
    }
  }
  if offsets.length() > 0 {
    apeglm_validate_matrix(
      offsets,
      counts.length(),
      samples,
      "apeglm offsets",
      false,
    )
  }
  if weights.length() > 0 {
    apeglm_validate_matrix(
      weights,
      counts.length(),
      samples,
      "apeglm weights",
      true,
    )
    for row in weights {
      let mut positive = false
      for value in row {
        if value > 0.0 {
          positive = true
        }
      }
      if !positive {
        raise ApeglmError(
          "apeglm each feature must have at least one positive weight",
        )
      }
    }
  }
  if gene_names.length() > 0 && gene_names.length() != counts.length() {
    raise ApeglmError("apeglm gene names must match the feature count")
  }
  if coefficient_names.length() > 0 &&
    coefficient_names.length() != coefficients {
    raise ApeglmError(
      "apeglm coefficient names must match the design column count",
    )
  }
  (counts.length(), samples, coefficients, coefficient)
}

///|
fn apeglm_names(
  names : Array[String],
  count : Int,
  prefix : String,
) -> Array[String] raise ApeglmError {
  if names.length() == 0 {
    let generated : Array[String] = []
    for index in 0..<count {
      generated.push(prefix + (index + 1).to_string())
    }
    return generated
  }
  let output = names.copy()
  for index in 0..<output.length() {
    if output[index].length() == 0 {
      raise ApeglmError("apeglm names cannot be empty")
    }
    for previous in 0..<index {
      if output[index] == output[previous] {
        raise ApeglmError("apeglm names must be unique: " + output[index])
      }
    }
  }
  output
}

///|
fn apeglm_offsets(
  offsets : Array[Array[Double]],
  genes : Int,
  samples : Int,
) -> Array[Array[Double]] {
  if offsets.length() == 0 {
    apeglm_zero_matrix(genes, samples)
  } else {
    apeglm_copy_matrix(offsets)
  }
}

///|
fn apeglm_weights(
  weights : Array[Array[Double]],
  genes : Int,
  samples : Int,
) -> Array[Array[Double]] {
  if weights.length() == 0 {
    let output : Array[Array[Double]] = []
    for _ in 0..<genes {
      output.push(Array::make(samples, 1.0))
    }
    output
  } else {
    apeglm_copy_matrix(weights)
  }
}

///|
pub fn apeglm_log_likelihood_nb(
  counts : Array[Double],
  design : Array[Array[Double]],
  beta : Array[Double],
  dispersion : Double,
  offsets? : Array[Double] = [],
  weights? : Array[Double] = [],
) -> Double raise ApeglmError {
  if counts.length() == 0 || design.length() != counts.length() {
    raise ApeglmError(
      "apeglm log-likelihood counts and design must have matching rows",
    )
  }
  let coefficients = beta.length()
  for row in design {
    if row.length() != coefficients {
      raise ApeglmError("apeglm log-likelihood design columns must match beta")
    }
  }
  if !apeglm_is_finite(dispersion) || dispersion <= 0.0 {
    raise ApeglmError("apeglm log-likelihood dispersion must be positive")
  }
  if offsets.length() > 0 && offsets.length() != counts.length() {
    raise ApeglmError("apeglm log-likelihood offsets have invalid length")
  }
  if weights.length() > 0 && weights.length() != counts.length() {
    raise ApeglmError("apeglm log-likelihood weights have invalid length")
  }
  let size = 1.0 / dispersion
  let mut result = 0.0
  for sample in 0..<counts.length() {
    let count = counts[sample]
    if !apeglm_is_finite(count) || count < 0.0 {
      raise ApeglmError(
        "apeglm log-likelihood counts must be finite and non-negative",
      )
    }
    let mut eta = if offsets.length() == 0 { 0.0 } else { offsets[sample] }
    for coefficient in 0..<coefficients {
      eta = eta + design[sample][coefficient] * beta[coefficient]
    }
    let bounded_eta = apeglm_clamp(eta, -40.0, 40.0)
    let mean = @math.exp(bounded_eta)
    let weight = if weights.length() == 0 { 1.0 } else { weights[sample] }
    result = result +
      weight * (count * bounded_eta - (count + size) * @math.ln(size + mean))
  }
  result
}

///|
fn apeglm_evaluate(
  counts : Array[Double],
  design : Array[Array[Double]],
  dispersion : Double,
  offsets : Array[Double],
  weights : Array[Double],
  beta : Array[Double],
  coefficient : Int,
  with_prior : Bool,
  prior_scale : Double,
  prior_df : Double,
  no_shrink_scale : Double,
) -> ApeglmEvaluation {
  let coefficients = beta.length()
  let gradient = Array::make(coefficients, 0.0)
  let information = apeglm_zero_matrix(coefficients, coefficients)
  let size = 1.0 / dispersion
  let mut value = 0.0
  for sample in 0..<counts.length() {
    let mut eta = offsets[sample]
    for index in 0..<coefficients {
      eta = eta + design[sample][index] * beta[index]
    }
    eta = apeglm_clamp(eta, -40.0, 40.0)
    let mean = @math.exp(eta)
    let denominator = size + mean
    let count = counts[sample]
    let weight = weights[sample]
    value = value +
      weight * (count * eta - (count + size) * @math.ln(denominator))
    let residual = weight * (count - (count + size) * mean / denominator)
    let curvature = weight *
      (count + size) *
      size *
      mean /
      (denominator * denominator)
    for left in 0..<coefficients {
      gradient[left] = gradient[left] + design[sample][left] * residual
      for right in 0..<coefficients {
        information[left][right] = information[left][right] +
          design[sample][left] * design[sample][right] * curvature
      }
    }
  }
  if with_prior {
    let normal_variance = no_shrink_scale * no_shrink_scale
    let t_variance_term = prior_df * prior_scale * prior_scale
    for index in 0..<coefficients {
      let estimate = beta[index]
      if index == coefficient {
        let denominator = t_variance_term + estimate * estimate
        value = value -
          0.5 *
          (prior_df + 1.0) *
          @math.ln(1.0 + estimate * estimate / t_variance_term)
        gradient[index] = gradient[index] -
          (prior_df + 1.0) * estimate / denominator
        information[index][index] = information[index][index] +
          (prior_df + 1.0) *
          (t_variance_term - estimate * estimate) /
          (denominator * denominator)
      } else {
        value = value - 0.5 * estimate * estimate / normal_variance
        gradient[index] = gradient[index] - estimate / normal_variance
        information[index][index] = information[index][index] +
          1.0 / normal_variance
      }
    }
  }
  ApeglmEvaluation::{ value, gradient, information, }
}

///|
fn apeglm_initial_beta(
  counts : Array[Double],
  design : Array[Array[Double]],
  offsets : Array[Double],
  weights : Array[Double],
) -> Array[Double] {
  let coefficients = design[0].length()
  let beta = Array::make(coefficients, 0.0)
  let mut count_total = 0.0
  let mut exposure_total = 0.0
  for sample in 0..<counts.length() {
    count_total = count_total + weights[sample] * counts[sample]
    exposure_total = exposure_total +
      weights[sample] * @math.exp(apeglm_clamp(offsets[sample], -40.0, 40.0))
  }
  beta[0] = @math.ln((count_total + 0.5) / (exposure_total + 0.5))
  beta
}

///|
fn apeglm_optimize(
  counts : Array[Double],
  design : Array[Array[Double]],
  dispersion : Double,
  offsets : Array[Double],
  weights : Array[Double],
  initial : Array[Double],
  coefficient : Int,
  with_prior : Bool,
  prior_scale : Double,
  config : ApeglmConfig,
) -> ApeglmOptimization {
  let beta = initial.copy()
  let mut evaluation = apeglm_evaluate(
    counts,
    design,
    dispersion,
    offsets,
    weights,
    beta,
    coefficient,
    with_prior,
    prior_scale,
    config.prior_df,
    config.prior_no_shrink_scale,
  )
  let mut converged = false
  let mut iterations = 0
  for iteration in 0..<config.max_iterations {
    iterations = iteration + 1
    let step = apeglm_solve_information(
      evaluation.information,
      evaluation.gradient,
    )
    let mut maximum = 0.0
    for value in step {
      maximum = maximum.max(value.abs())
    }
    let step_scale = if maximum > config.max_step {
      config.max_step / maximum
    } else {
      1.0
    }
    let mut line_scale = step_scale
    let mut accepted = false
    let mut accepted_change = 0.0
    let mut accepted_beta = beta.copy()
    let mut accepted_evaluation = evaluation
    for _ in 0..<28 {
      let candidate = beta.copy()
      let mut change = 0.0
      for index in 0..<candidate.length() {
        let delta = line_scale * step[index]
        candidate[index] = apeglm_clamp(candidate[index] + delta, -30.0, 30.0)
        change = change.max((candidate[index] - beta[index]).abs())
      }
      let next = apeglm_evaluate(
        counts,
        design,
        dispersion,
        offsets,
        weights,
        candidate,
        coefficient,
        with_prior,
        prior_scale,
        config.prior_df,
        config.prior_no_shrink_scale,
      )
      if apeglm_is_finite(next.value) &&
        next.value >= evaluation.value - 1.0e-10 {
        accepted = true
        accepted_change = change
        accepted_beta = candidate
        accepted_evaluation = next
        break
      }
      line_scale = line_scale * 0.5
    }
    if !accepted {
      let mut maximum_gradient = 0.0
      for value in evaluation.gradient {
        maximum_gradient = maximum_gradient.max(value.abs())
      }
      converged = maximum_gradient <= config.tolerance * 10.0
      break
    }
    for index in 0..<beta.length() {
      beta[index] = accepted_beta[index]
    }
    evaluation = accepted_evaluation
    let mut beta_scale = 1.0
    for value in beta {
      beta_scale = beta_scale.max(value.abs())
    }
    if accepted_change <= config.tolerance * beta_scale {
      converged = true
      break
    }
  }
  ApeglmOptimization::{
    beta,
    covariance: apeglm_covariance(evaluation.information),
    value: evaluation.value,
    converged,
    iterations,
  }
}

///|
fn apeglm_target_starts(
  mle : Array[Double],
  coefficient : Int,
  prior_scale : Double,
  count : Int,
) -> Array[Array[Double]] {
  let starts : Array[Array[Double]] = [mle.copy()]
  if starts.length() < count {
    let center = mle.copy()
    center[coefficient] = 0.0
    starts.push(center)
  }
  let locations = [
    2.0 * prior_scale,
    -2.0 * prior_scale,
    4.0 * prior_scale,
    -4.0 * prior_scale,
    prior_scale,
    -prior_scale,
  ]
  let mut index = 0
  while starts.length() < count && index < locations.length() {
    let candidate = mle.copy()
    candidate[coefficient] = locations[index]
    starts.push(candidate)
    index = index + 1
  }
  let mut shell = 8.0
  while starts.length() < count {
    let candidate = mle.copy()
    let sign = if starts.length() % 2 == 0 { 1.0 } else { -1.0 }
    candidate[coefficient] = apeglm_clamp(
      sign * shell * prior_scale,
      -30.0,
      30.0,
    )
    starts.push(candidate)
    if sign < 0.0 {
      shell = shell * 2.0
    }
  }
  starts
}

///|
fn apeglm_fit_map(
  counts : Array[Double],
  design : Array[Array[Double]],
  dispersion : Double,
  offsets : Array[Double],
  weights : Array[Double],
  mle : Array[Double],
  coefficient : Int,
  prior_scale : Double,
  config : ApeglmConfig,
) -> ApeglmOptimization {
  let starts = apeglm_target_starts(
    mle,
    coefficient,
    prior_scale,
    config.random_starts,
  )
  let mut best = apeglm_optimize(
    counts,
    design,
    dispersion,
    offsets,
    weights,
    starts[0],
    coefficient,
    true,
    prior_scale,
    config,
  )
  for index in 1..<starts.length() {
    let candidate = apeglm_optimize(
      counts,
      design,
      dispersion,
      offsets,
      weights,
      starts[index],
      coefficient,
      true,
      prior_scale,
      config,
    )
    if candidate.value > best.value {
      best = candidate
    }
  }
  best
}

///|
fn apeglm_standard_errors(covariance : Array[Array[Double]]) -> Array[Double] {
  let output = Array::make(covariance.length(), 0.0)
  for index in 0..<covariance.length() {
    output[index] = covariance[index][index].max(1.0e-16).sqrt()
  }
  output
}

///|
fn apeglm_prior_objective(
  variance : Double,
  estimates : Array[Double],
  standard_errors : Array[Double],
) -> Double {
  let mut numerator = 0.0
  let mut denominator = 0.0
  for index in 0..<estimates.length() {
    let sampling_variance = standard_errors[index] * standard_errors[index]
    let inverse_information = 1.0 /
      (2.0 * (variance + sampling_variance) * (variance + sampling_variance))
    numerator = numerator +
      (estimates[index] * estimates[index] - sampling_variance) *
      inverse_information
    denominator = denominator + inverse_information
  }
  if denominator <= 0.0 {
    -variance
  } else {
    numerator / denominator - variance
  }
}

///|
pub fn apeglm_prior_variance(
  estimates : Array[Double],
  standard_errors : Array[Double],
  minimum? : Double = 1.0e-6,
  maximum? : Double = 400.0,
) -> Double raise ApeglmError {
  if estimates.length() == 0 || estimates.length() != standard_errors.length() {
    raise ApeglmError(
      "apeglm prior adaptation requires matching non-empty MLE and SE arrays",
    )
  }
  if !apeglm_is_finite(minimum) ||
    !apeglm_is_finite(maximum) ||
    minimum <= 0.0 ||
    maximum <= minimum {
    raise ApeglmError("apeglm prior variance bounds are invalid")
  }
  for index in 0..<estimates.length() {
    if !apeglm_is_finite(estimates[index]) ||
      !apeglm_is_finite(standard_errors[index]) ||
      standard_errors[index] <= 0.0 {
      raise ApeglmError(
        "apeglm prior adaptation requires finite MLEs and positive SEs",
      )
    }
  }
  let lower_value = apeglm_prior_objective(minimum, estimates, standard_errors)
  if lower_value <= 0.0 {
    return minimum
  }
  let upper_value = apeglm_prior_objective(maximum, estimates, standard_errors)
  if upper_value >= 0.0 {
    return maximum
  }
  let mut lower = minimum
  let mut upper = maximum
  for _ in 0..<100 {
    let middle = 0.5 * (lower + upper)
    if apeglm_prior_objective(middle, estimates, standard_errors) > 0.0 {
      lower = middle
    } else {
      upper = middle
    }
  }
  0.5 * (lower + upper)
}

///|
pub fn apeglm_s_values(local_fsr : Array[Double]) -> Array[Double] {
  let count = local_fsr.length()
  let order : Array[Int] = []
  for index in 0..<count {
    order.push(index)
  }
  for left in 0..<count {
    for right in 0..<(count - left - 1) {
      let first = local_fsr[order[right]]
      let second = local_fsr[order[right + 1]]
      if first > second || (first == second && order[right] > order[right + 1]) {
        let temporary = order[right]
        order[right] = order[right + 1]
        order[right + 1] = temporary
      }
    }
  }
  let output = Array::make(count, 0.0)
  let mut cumulative = 0.0
  for rank in 0..<count {
    cumulative = cumulative + local_fsr[order[rank]]
    output[order[rank]] = cumulative / (rank + 1).to_double()
  }
  output
}

///|
fn apeglm_base_mean(counts : Array[Double], offsets : Array[Double]) -> Double {
  let normalized : Array[Double] = []
  for sample in 0..<counts.length() {
    normalized.push(
      counts[sample] / @math.exp(apeglm_clamp(offsets[sample], -40.0, 40.0)),
    )
  }
  apeglm_mean(normalized)
}

///|
fn apeglm_gene_with_s_value(
  gene : ApeglmGeneResult,
  s_value : Double,
) -> ApeglmGeneResult {
  ApeglmGeneResult::{
    gene_name: gene.gene_name,
    base_mean: gene.base_mean,
    mle: gene.mle,
    mle_sd: gene.mle_sd,
    map: gene.map,
    posterior_sd: gene.posterior_sd,
    fsr: gene.fsr,
    s_value,
    threshold_probability: gene.threshold_probability,
    interval_lower: gene.interval_lower,
    interval_upper: gene.interval_upper,
    log_posterior: gene.log_posterior,
    converged: gene.converged,
    iterations: gene.iterations,
  }
}

///|
pub fn apeglm_fit(
  counts : Array[Array[Double]],
  design : Array[Array[Double]],
  dispersions : Array[Double],
  config? : ApeglmConfig = ApeglmConfig::default(),
  offsets? : Array[Array[Double]] = [],
  weights? : Array[Array[Double]] = [],
  gene_names? : Array[String] = [],
  coefficient_names? : Array[String] = [],
) -> ApeglmResult raise ApeglmError {
  let (gene_count, sample_count, coefficient_count, coefficient) = apeglm_validate_inputs(
    counts, design, dispersions, offsets, weights, gene_names, coefficient_names,
    config,
  )
  let names = apeglm_names(gene_names, gene_count, "feature_")
  let coefficients = apeglm_names(
    coefficient_names, coefficient_count, "coefficient_",
  )
  let model_offsets = apeglm_offsets(offsets, gene_count, sample_count)
  let model_weights = apeglm_weights(weights, gene_count, sample_count)
  let mle_fits : Array[ApeglmOptimization] = []
  let target_mle : Array[Double] = []
  let target_se : Array[Double] = []
  for gene in 0..<gene_count {
    let initial = apeglm_initial_beta(
      counts[gene],
      design,
      model_offsets[gene],
      model_weights[gene],
    )
    let fit = apeglm_optimize(
      counts[gene],
      design,
      dispersions[gene],
      model_offsets[gene],
      model_weights[gene],
      initial,
      coefficient,
      false,
      config.prior_scale,
      config,
    )
    let standard_errors = apeglm_standard_errors(fit.covariance)
    mle_fits.push(fit)
    target_mle.push(fit.beta[coefficient])
    target_se.push(standard_errors[coefficient])
  }
  let prior_variance = if config.adaptive_prior {
    apeglm_prior_variance(target_mle, target_se)
  } else {
    config.prior_scale * config.prior_scale
  }
  let prior_scale = if config.adaptive_prior {
    (config.multiplier * prior_variance.sqrt()).min(1.0).max(0.001)
  } else {
    config.prior_scale
  }
  let raw_genes : Array[ApeglmGeneResult] = []
  let local_fsr : Array[Double] = []
  let alpha = 1.0 - config.interval_level
  let quantile = apeglm_normal_quantile(1.0 - alpha / 2.0)
  for gene in 0..<gene_count {
    let mle_fit = mle_fits[gene]
    let mle_sd = apeglm_standard_errors(mle_fit.covariance)
    let map_fit = apeglm_fit_map(
      counts[gene],
      design,
      dispersions[gene],
      model_offsets[gene],
      model_weights[gene],
      mle_fit.beta,
      coefficient,
      prior_scale,
      config,
    )
    let posterior_sd = apeglm_standard_errors(map_fit.covariance)
    let estimate = map_fit.beta[coefficient]
    let standard_error = posterior_sd[coefficient].max(1.0e-12)
    let fsr = apeglm_clamp(
      apeglm_normal_cdf(-estimate.abs() / standard_error),
      0.0,
      0.5,
    )
    let threshold_probability = if estimate >= 0.0 {
      apeglm_normal_cdf((config.threshold - estimate) / standard_error)
    } else {
      1.0 - apeglm_normal_cdf((-config.threshold - estimate) / standard_error)
    }
    local_fsr.push(fsr)
    raw_genes.push(ApeglmGeneResult::{
      gene_name: names[gene],
      base_mean: apeglm_base_mean(counts[gene], model_offsets[gene]),
      mle: mle_fit.beta.copy(),
      mle_sd,
      map: map_fit.beta.copy(),
      posterior_sd,
      fsr,
      s_value: 0.0,
      threshold_probability: apeglm_clamp(threshold_probability, 0.0, 1.0),
      interval_lower: estimate - quantile * standard_error,
      interval_upper: estimate + quantile * standard_error,
      log_posterior: map_fit.value,
      converged: map_fit.converged,
      iterations: map_fit.iterations,
    })
  }
  let s_values = apeglm_s_values(local_fsr)
  let genes : Array[ApeglmGeneResult] = []
  for index in 0..<gene_count {
    genes.push(apeglm_gene_with_s_value(raw_genes[index], s_values[index]))
  }
  let no_shrink : Array[Int] = []
  for index in 0..<coefficient_count {
    if index != coefficient {
      no_shrink.push(index)
    }
  }
  ApeglmResult::{
    gene_names: names,
    coefficient_names: coefficients,
    coefficient,
    genes,
    prior_control: ApeglmPriorControl::{
      no_shrink,
      prior_mean: 0.0,
      prior_scale,
      prior_variance,
      prior_df: config.prior_df,
      prior_no_shrink_scale: config.prior_no_shrink_scale,
      adaptive: config.adaptive_prior,
    },
    threshold: config.threshold,
    interval_level: config.interval_level,
    config,
  }
}

///|
pub fn ApeglmGeneResult::effect(
  self : ApeglmGeneResult,
  coefficient : Int,
  log2_scale? : Bool = false,
) -> Double {
  let value = self.map[coefficient]
  if log2_scale {
    value / @math.ln(2.0)
  } else {
    value
  }
}

///|
pub fn ApeglmGeneResult::standard_error(
  self : ApeglmGeneResult,
  coefficient : Int,
  log2_scale? : Bool = false,
) -> Double {
  let value = self.posterior_sd[coefficient]
  if log2_scale {
    value / @math.ln(2.0)
  } else {
    value
  }
}

///|
pub fn ApeglmResult::map_matrix(self : ApeglmResult) -> Array[Array[Double]] {
  let output : Array[Array[Double]] = []
  for gene in self.genes {
    output.push(gene.map.copy())
  }
  output
}

///|
pub fn ApeglmResult::sd_matrix(self : ApeglmResult) -> Array[Array[Double]] {
  let output : Array[Array[Double]] = []
  for gene in self.genes {
    output.push(gene.posterior_sd.copy())
  }
  output
}

///|
pub fn ApeglmResult::gene(
  self : ApeglmResult,
  index : Int,
) -> ApeglmGeneResult? {
  if index < 0 || index >= self.genes.length() {
    None
  } else {
    Some(self.genes[index])
  }
}

///|
pub fn ApeglmResult::find_gene(
  self : ApeglmResult,
  name : String,
) -> ApeglmGeneResult? {
  for gene in self.genes {
    if gene.gene_name == name {
      return Some(gene)
    }
  }
  None
}

///|
pub fn ApeglmResult::ranked(self : ApeglmResult) -> Array[ApeglmGeneResult] {
  let output = self.genes.copy()
  for left in 0..<output.length() {
    for right in 0..<(output.length() - left - 1) {
      let first = output[right]
      let second = output[right + 1]
      let first_effect = first.map[self.coefficient].abs()
      let second_effect = second.map[self.coefficient].abs()
      if first.s_value > second.s_value ||
        (first.s_value == second.s_value && first_effect < second_effect) {
        output[right] = second
        output[right + 1] = first
      }
    }
  }
  output
}

///|
pub fn ApeglmResult::select(
  self : ApeglmResult,
  maximum_s_value? : Double = 0.05,
  minimum_absolute_effect? : Double = 0.0,
  use_threshold_probability? : Bool = false,
) -> Array[ApeglmGeneResult] {
  let output : Array[ApeglmGeneResult] = []
  for gene in self.genes {
    let probability = if use_threshold_probability {
      gene.threshold_probability
    } else {
      gene.s_value
    }
    if probability <= maximum_s_value &&
      gene.map[self.coefficient].abs() >= minimum_absolute_effect {
      output.push(gene)
    }
  }
  output
}

///|
pub fn ApeglmResult::summary(
  self : ApeglmResult,
  maximum_s_value? : Double = 0.05,
) -> ApeglmSummary {
  let mut converged_count = 0
  let mut low_fsr_count = 0
  let effects : Array[Double] = []
  for gene in self.genes {
    if gene.converged {
      converged_count = converged_count + 1
    }
    if gene.fsr <= maximum_s_value {
      low_fsr_count = low_fsr_count + 1
    }
    effects.push(gene.map[self.coefficient].abs())
  }
  effects.sort()
  let median_absolute_effect = if effects.length() == 0 {
    0.0
  } else if effects.length() % 2 == 1 {
    effects[effects.length() / 2]
  } else {
    let middle = effects.length() / 2
    0.5 * (effects[middle - 1] + effects[middle])
  }
  ApeglmSummary::{
    gene_count: self.genes.length(),
    converged_count,
    selected_count: self.select(maximum_s_value~).length(),
    low_fsr_count,
    prior_scale: self.prior_control.prior_scale,
    median_absolute_effect,
    coefficient_name: self.coefficient_names[self.coefficient],
  }
}

///|
pub fn ApeglmResult::to_tsv(
  self : ApeglmResult,
  log2_scale? : Bool = true,
) -> String {
  let output = StringBuilder()
  let scale_label = if log2_scale { "log2" } else { "ln" }
  output.write_string(
    "gene\tbaseMean\t" +
    scale_label +
    "MLE\t" +
    scale_label +
    "MAP\tposteriorSD\tFSR\tsvalue\tFSOS\tintervalLower\tintervalUpper\tconverged\n",
  )
  let divisor = if log2_scale { @math.ln(2.0) } else { 1.0 }
  for gene in self.genes {
    output.write_string(gene.gene_name)
    output.write_string("\t")
    output.write_string(gene.base_mean.to_string())
    output.write_string("\t")
    output.write_string((gene.mle[self.coefficient] / divisor).to_string())
    output.write_string("\t")
    output.write_string((gene.map[self.coefficient] / divisor).to_string())
    output.write_string("\t")
    output.write_string(
      (gene.posterior_sd[self.coefficient] / divisor).to_string(),
    )
    output.write_string("\t")
    output.write_string(gene.fsr.to_string())
    output.write_string("\t")
    output.write_string(gene.s_value.to_string())
    output.write_string("\t")
    output.write_string(gene.threshold_probability.to_string())
    output.write_string("\t")
    output.write_string((gene.interval_lower / divisor).to_string())
    output.write_string("\t")
    output.write_string((gene.interval_upper / divisor).to_string())
    output.write_string("\t")
    output.write_string(gene.converged.to_string())
    output.write_string("\n")
  }
  output.to_string()
}

///|
pub fn apeglm_from_deseq2(
  dataset : DESeqDataSet,
  config? : ApeglmConfig = ApeglmConfig::default(),
  coefficient_names? : Array[String] = [],
) -> ApeglmResult raise ApeglmError {
  let counts : Array[Array[Double]] = []
  for row in dataset.counts {
    let converted : Array[Double] = []
    for value in row {
      converted.push(value.to_double())
    }
    counts.push(converted)
  }
  let samples = if counts.length() == 0 { 0 } else { counts[0].length() }
  if dataset.size_factors.length() != samples {
    raise ApeglmError("apeglm DESeq2 size factors must match the sample count")
  }
  let offsets : Array[Array[Double]] = []
  for _ in 0..<counts.length() {
    let row : Array[Double] = []
    for factor in dataset.size_factors {
      if !apeglm_is_finite(factor) || factor <= 0.0 {
        raise ApeglmError(
          "apeglm DESeq2 size factors must be finite and positive",
        )
      }
      row.push(@math.ln(factor))
    }
    offsets.push(row)
  }
  apeglm_fit(
    counts,
    dataset.design,
    dataset.dispersions,
    config~,
    offsets~,
    gene_names=dataset.row_names,
    coefficient_names~,
  )
}

///|
fn apeglm_copy_summarized_experiment(
  experiment : SummarizedExperiment,
) -> SummarizedExperiment {
  let assays : Map[String, Array[Array[Double]]] = Map([])
  for key in experiment.assays.keys() {
    assays[key] = apeglm_copy_matrix(experiment.assays[key])
  }
  let metadata : Map[String, String] = Map([])
  for key in experiment.metadata.keys() {
    metadata[key] = experiment.metadata[key]
  }
  SummarizedExperiment::{
    assays,
    row_ranges: experiment.row_ranges.copy(),
    col_data: experiment.col_data.copy(),
    metadata,
  }
}

///|
pub fn apeglm_summarized_experiment(
  experiment : SummarizedExperiment,
  design : Array[Array[Double]],
  dispersions : Array[Double],
  config? : ApeglmConfig = ApeglmConfig::default(),
  assay_name? : String = "counts",
  offsets? : Array[Array[Double]] = [],
  weights? : Array[Array[Double]] = [],
  gene_names? : Array[String] = [],
  coefficient_names? : Array[String] = [],
) -> ApeglmSummarizedExperimentOutput raise ApeglmError {
  let counts = match se_assay(experiment, assay_name) {
    Some(value) => value
    None =>
      raise ApeglmError(
        "apeglm SummarizedExperiment assay not found: " + assay_name,
      )
  }
  let result = apeglm_fit(
    counts,
    design,
    dispersions,
    config~,
    offsets~,
    weights~,
    gene_names~,
    coefficient_names~,
  )
  let enriched = apeglm_copy_summarized_experiment(experiment)
  let map : Array[Array[Double]] = []
  let standard_error : Array[Array[Double]] = []
  let fsr : Array[Array[Double]] = []
  let s_value : Array[Array[Double]] = []
  let threshold_probability : Array[Array[Double]] = []
  for gene in result.genes {
    map.push([gene.map[result.coefficient]])
    standard_error.push([gene.posterior_sd[result.coefficient]])
    fsr.push([gene.fsr])
    s_value.push([gene.s_value])
    threshold_probability.push([gene.threshold_probability])
  }
  enriched.assays["apeglm_map"] = map
  enriched.assays["apeglm_sd"] = standard_error
  enriched.assays["apeglm_fsr"] = fsr
  enriched.assays["apeglm_svalue"] = s_value
  enriched.assays["apeglm_fsos"] = threshold_probability
  enriched.metadata["apeglm_coefficient"] = result.coefficient_names[result.coefficient]
  enriched.metadata["apeglm_prior_scale"] = result.prior_control.prior_scale.to_string()
  enriched.metadata["apeglm_scale"] = "natural-log"
  ApeglmSummarizedExperimentOutput::{ experiment: enriched, result, }
}

///|
pub fn apeglm_example_data() -> (
  Array[Array[Double]],
  Array[Array[Double]],
  Array[Double],
  Array[String],
  Array[String],
) {
  let counts = [
    [48.0, 52.0, 45.0, 54.0, 198.0, 220.0, 205.0, 230.0],
    [160.0, 148.0, 171.0, 155.0, 39.0, 44.0, 35.0, 41.0],
    [72.0, 68.0, 75.0, 70.0, 79.0, 74.0, 77.0, 73.0],
    [0.0, 1.0, 0.0, 0.0, 4.0, 0.0, 5.0, 0.0],
    [12.0, 20.0, 15.0, 18.0, 27.0, 19.0, 31.0, 23.0],
    [310.0, 290.0, 325.0, 305.0, 640.0, 615.0, 670.0, 650.0],
    [9.0, 8.0, 11.0, 10.0, 8.0, 12.0, 9.0, 11.0],
    [3.0, 0.0, 2.0, 1.0, 0.0, 1.0, 0.0, 2.0],
  ]
  let design = [
    [1.0, 0.0, 0.0],
    [1.0, 1.0, 0.0],
    [1.0, 0.0, 0.0],
    [1.0, 1.0, 0.0],
    [1.0, 0.0, 1.0],
    [1.0, 1.0, 1.0],
    [1.0, 0.0, 1.0],
    [1.0, 1.0, 1.0],
  ]
  let dispersions = [0.08, 0.08, 0.1, 0.8, 0.25, 0.05, 0.2, 0.9]
  let genes = [
    "strong_up", "strong_down", "stable", "low_noisy", "moderate_up", "abundant_up",
    "null_low", "sparse_null",
  ]
  let coefficient_names = ["intercept", "batch", "condition"]
  (counts, design, dispersions, genes, coefficient_names)
}
