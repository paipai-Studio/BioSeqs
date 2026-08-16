// Dirichlet-multinomial mixture clustering and classification inspired by
// Bioconductor DirichletMultinomial. Count matrices use sample x taxon
// orientation, matching the upstream dmn() API.

///|
pub suberror DirichletMultinomialError {
  DirichletMultinomialError(String)
}

///|
pub(all) enum DmnCriterion {
  DmnLaplace
  DmnBic
  DmnAic
} derive(Eq, Debug)

///|
pub struct DmnConfig {
  components : Int
  max_iterations : Int
  optimizer_iterations : Int
  soft_kmeans_iterations : Int
  tolerance : Double
  optimizer_tolerance : Double
  soft_beta : Double
  seed : Int
  prior_shape : Double
  prior_rate : Double
  min_alpha : Double
} derive(Eq, Debug)

///|
pub fn DmnConfig::create(
  components? : Int = 1,
  max_iterations? : Int = 100,
  optimizer_iterations? : Int = 250,
  soft_kmeans_iterations? : Int = 200,
  tolerance? : Double = 1.0e-6,
  optimizer_tolerance? : Double = 1.0e-5,
  soft_beta? : Double = 50.0,
  seed? : Int = 42,
  prior_shape? : Double = 0.1,
  prior_rate? : Double = 0.1,
  min_alpha? : Double = 1.0e-8,
) -> DmnConfig raise DirichletMultinomialError {
  if components < 1 {
    raise DirichletMultinomialError(
      "DirichletMultinomial component count must be positive",
    )
  }
  if max_iterations < 1 ||
    optimizer_iterations < 1 ||
    soft_kmeans_iterations < 1 {
    raise DirichletMultinomialError(
      "DirichletMultinomial iteration limits must be positive",
    )
  }
  if !dmn_is_finite(tolerance) || tolerance <= 0.0 {
    raise DirichletMultinomialError(
      "DirichletMultinomial EM tolerance must be finite and positive",
    )
  }
  if !dmn_is_finite(optimizer_tolerance) || optimizer_tolerance <= 0.0 {
    raise DirichletMultinomialError(
      "DirichletMultinomial optimizer tolerance must be finite and positive",
    )
  }
  if !dmn_is_finite(soft_beta) || soft_beta <= 0.0 {
    raise DirichletMultinomialError(
      "DirichletMultinomial soft k-means beta must be finite and positive",
    )
  }
  if seed <= 0 {
    raise DirichletMultinomialError(
      "DirichletMultinomial random seed must be positive",
    )
  }
  if !dmn_is_finite(prior_shape) || prior_shape <= 0.0 {
    raise DirichletMultinomialError(
      "DirichletMultinomial Gamma prior shape must be finite and positive",
    )
  }
  if !dmn_is_finite(prior_rate) || prior_rate <= 0.0 {
    raise DirichletMultinomialError(
      "DirichletMultinomial Gamma prior rate must be finite and positive",
    )
  }
  if !dmn_is_finite(min_alpha) || min_alpha <= 0.0 || min_alpha >= 1.0 {
    raise DirichletMultinomialError(
      "DirichletMultinomial minimum alpha must be finite and in (0, 1)",
    )
  }
  DmnConfig::{
    components,
    max_iterations,
    optimizer_iterations,
    soft_kmeans_iterations,
    tolerance,
    optimizer_tolerance,
    soft_beta,
    seed,
    prior_shape,
    prior_rate,
    min_alpha,
  }
}

///|
pub fn DmnConfig::default() -> DmnConfig {
  DmnConfig::{
    components: 1,
    max_iterations: 100,
    optimizer_iterations: 250,
    soft_kmeans_iterations: 200,
    tolerance: 1.0e-6,
    optimizer_tolerance: 1.0e-5,
    soft_beta: 50.0,
    seed: 42,
    prior_shape: 0.1,
    prior_rate: 0.1,
    min_alpha: 1.0e-8,
  }
}

///|
pub fn DmnConfig::with_components(
  self : DmnConfig,
  components : Int,
  seed? : Int = self.seed,
) -> DmnConfig raise DirichletMultinomialError {
  DmnConfig::create(
    components~,
    max_iterations=self.max_iterations,
    optimizer_iterations=self.optimizer_iterations,
    soft_kmeans_iterations=self.soft_kmeans_iterations,
    tolerance=self.tolerance,
    optimizer_tolerance=self.optimizer_tolerance,
    soft_beta=self.soft_beta,
    seed~,
    prior_shape=self.prior_shape,
    prior_rate=self.prior_rate,
    min_alpha=self.min_alpha,
  )
}

///|
pub struct DmnGoodnessOfFit {
  negative_log_evidence : Double
  log_determinant : Double
  laplace : Double
  bic : Double
  aic : Double
} derive(Eq, Debug)

///|
pub struct DmnFit {
  estimates : Array[Array[Double]]
  lower : Array[Array[Double]]
  upper : Array[Array[Double]]
  weights : Array[Double]
  responsibilities : Array[Array[Double]]
  sample_names : Array[String]
  taxon_names : Array[String]
  goodness : DmnGoodnessOfFit
  log_likelihoods : Array[Double]
  converged : Bool
  iterations : Int
  config : DmnConfig
}

///|
pub struct DmnModelSelection {
  fits : Array[DmnFit]
  scores : Array[Double]
  best_index : Int
  criterion : DmnCriterion
}

///|
pub struct DmnGroupFit {
  group_names : Array[String]
  models : Array[DmnFit]
  priors : Array[Double]
  taxon_names : Array[String]
}

///|
pub struct DmnCrossValidation {
  probabilities : Array[Array[Double]]
  predictions : Array[String]
  truth : Array[String]
  group_names : Array[String]
  fold_ids : Array[Int]
  correct : Int
  accuracy : Double
}

///|
pub struct DmnRoc {
  thresholds : Array[Double]
  sensitivity : Array[Double]
  specificity : Array[Double]
  auc : Double
  best_threshold : Double
} derive(Eq, Debug)

///|
fn dmn_is_finite(value : Double) -> Bool {
  value == value && value.abs() <= 1.0e300
}

///|
fn dmn_make_matrix(
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
fn dmn_copy_matrix(values : Array[Array[Double]]) -> Array[Array[Double]] {
  let output : Array[Array[Double]] = []
  for row in values {
    output.push(row.copy())
  }
  output
}

///|
fn dmn_default_names(prefix : String, count : Int) -> Array[String] {
  let names : Array[String] = []
  for index in 0..<count {
    names.push(prefix + (index + 1).to_string())
  }
  names
}

///|
fn dmn_validate_names(
  names : Array[String],
  expected : Int,
  label : String,
) -> Unit raise DirichletMultinomialError {
  if names.length() != expected {
    raise DirichletMultinomialError(
      "DirichletMultinomial " + label + " name count does not match data",
    )
  }
  let seen : Map[String, Bool] = Map([])
  for name in names {
    let trimmed = name.trim().to_owned()
    if trimmed.length() == 0 {
      raise DirichletMultinomialError(
        "DirichletMultinomial " + label + " names must not be empty",
      )
    }
    if seen.contains(trimmed) {
      raise DirichletMultinomialError(
        "DirichletMultinomial " + label + " names must be unique",
      )
    }
    seen[trimmed] = true
  }
}

///|
fn dmn_validate_counts(
  counts : Array[Array[Int]],
) -> (Int, Int) raise DirichletMultinomialError {
  if counts.length() == 0 {
    raise DirichletMultinomialError(
      "DirichletMultinomial count matrix must contain samples",
    )
  }
  let taxa = counts[0].length()
  if taxa == 0 {
    raise DirichletMultinomialError(
      "DirichletMultinomial count matrix must contain taxa",
    )
  }
  for sample in 0..<counts.length() {
    if counts[sample].length() != taxa {
      raise DirichletMultinomialError(
        "DirichletMultinomial count matrix must be rectangular",
      )
    }
    let mut total = 0
    for count in counts[sample] {
      if count < 0 {
        raise DirichletMultinomialError(
          "DirichletMultinomial counts must be non-negative",
        )
      }
      total = total + count
    }
    if total == 0 {
      raise DirichletMultinomialError(
        "DirichletMultinomial every sample must have a positive library size",
      )
    }
  }
  (counts.length(), taxa)
}

///|
fn dmn_log_gamma(value : Double) -> Double {
  let coefficients = [
    0.9999999999998099, 676.5203681218851, -1259.1392167224028, 771.3234287776531,
    -176.6150291621406, 12.507343278686905, -0.13857109526572012, 9.984369578019572e-6,
    1.5056327351493116e-7,
  ]
  if value < 0.5 {
    let pi = 3.141592653589793
    return @math.ln(pi) -
      @math.ln(@math.sin(pi * value).abs()) -
      dmn_log_gamma(1.0 - value)
  }
  let shifted = value - 1.0
  let mut series = coefficients[0]
  for index in 1..<coefficients.length() {
    series = series + coefficients[index] / (shifted + index.to_double())
  }
  let t = shifted + 7.5
  0.9189385332046727 + (shifted + 0.5) * @math.ln(t) - t + @math.ln(series)
}

///|
fn dmn_digamma(value : Double) -> Double {
  let mut x = value
  let mut result = 0.0
  while x < 8.0 {
    result = result - 1.0 / x
    x = x + 1.0
  }
  let inverse = 1.0 / x
  let inverse2 = inverse * inverse
  result +
  @math.ln(x) -
  0.5 * inverse -
  inverse2 * (1.0 / 12.0 - inverse2 * (1.0 / 120.0 - inverse2 * (1.0 / 252.0)))
}

///|
fn dmn_trigamma(value : Double) -> Double {
  let mut x = value
  let mut result = 0.0
  while x < 8.0 {
    result = result + 1.0 / (x * x)
    x = x + 1.0
  }
  let inverse = 1.0 / x
  let inverse2 = inverse * inverse
  result +
  inverse +
  0.5 * inverse2 +
  inverse2 * inverse / 6.0 -
  inverse2 * inverse2 * inverse / 30.0 +
  inverse2 * inverse2 * inverse2 * inverse / 42.0
}

///|
fn dmn_log_sum_exp(values : Array[Double]) -> Double {
  let mut maximum = values[0]
  for value in values {
    if value > maximum {
      maximum = value
    }
  }
  let mut total = 0.0
  for value in values {
    total = total + @math.exp(value - maximum)
  }
  maximum + @math.ln(total)
}

///|
fn dmn_row_total(row : Array[Int]) -> Int {
  let mut total = 0
  for value in row {
    total = total + value
  }
  total
}

///|
fn dmn_component_log_kernel(
  counts : Array[Int],
  alpha : Array[Double],
) -> Double {
  let mut alpha_total = 0.0
  let mut count_total = 0
  let mut value = 0.0
  for taxon in 0..<alpha.length() {
    alpha_total = alpha_total + alpha[taxon]
    count_total = count_total + counts[taxon]
    value = value +
      dmn_log_gamma(alpha[taxon] + counts[taxon].to_double()) -
      dmn_log_gamma(alpha[taxon])
  }
  value +
  dmn_log_gamma(alpha_total) -
  dmn_log_gamma(alpha_total + count_total.to_double())
}

///|
pub fn dirichlet_multinomial_log_pmf(
  counts : Array[Int],
  alpha : Array[Double],
) -> Double raise DirichletMultinomialError {
  if counts.length() == 0 || counts.length() != alpha.length() {
    raise DirichletMultinomialError(
      "Dirichlet-multinomial count and alpha vectors must have equal positive length",
    )
  }
  let mut total = 0
  let mut coefficient = 0.0
  for index in 0..<counts.length() {
    if counts[index] < 0 {
      raise DirichletMultinomialError(
        "Dirichlet-multinomial counts must be non-negative",
      )
    }
    if !dmn_is_finite(alpha[index]) || alpha[index] <= 0.0 {
      raise DirichletMultinomialError(
        "Dirichlet-multinomial alpha values must be finite and positive",
      )
    }
    total = total + counts[index]
    coefficient = coefficient - dmn_log_gamma(counts[index].to_double() + 1.0)
  }
  coefficient = coefficient + dmn_log_gamma(total.to_double() + 1.0)
  coefficient + dmn_component_log_kernel(counts, alpha)
}

///|
pub fn dirichlet_multinomial_mean(
  total_count : Int,
  alpha : Array[Double],
) -> Array[Double] raise DirichletMultinomialError {
  if total_count < 0 || alpha.length() == 0 {
    raise DirichletMultinomialError(
      "Dirichlet-multinomial mean requires a non-negative count and positive alpha length",
    )
  }
  let mut concentration = 0.0
  for value in alpha {
    if !dmn_is_finite(value) || value <= 0.0 {
      raise DirichletMultinomialError(
        "Dirichlet-multinomial alpha values must be finite and positive",
      )
    }
    concentration = concentration + value
  }
  alpha.map(fn(value) { total_count.to_double() * value / concentration })
}

///|
pub fn dirichlet_multinomial_covariance(
  total_count : Int,
  alpha : Array[Double],
) -> Array[Array[Double]] raise DirichletMultinomialError {
  let mean = dirichlet_multinomial_mean(total_count, alpha)
  let taxa = alpha.length()
  let mut concentration = 0.0
  for value in alpha {
    concentration = concentration + value
  }
  let inflation = (total_count.to_double() + concentration) /
    (1.0 + concentration)
  let covariance = dmn_make_matrix(taxa, taxa, 0.0)
  for left in 0..<taxa {
    let probability_left = alpha[left] / concentration
    for right in 0..<taxa {
      let probability_right = alpha[right] / concentration
      covariance[left][right] = if left == right {
        mean[left] * (1.0 - probability_left) * inflation
      } else {
        -total_count.to_double() *
        probability_left *
        probability_right *
        inflation
      }
    }
  }
  covariance
}

///|
fn dmn_lcg_next(state : Int) -> Int {
  state * 48271 % 2147483647
}

///|
fn dmn_soft_kmeans(
  counts : Array[Array[Int]],
  config : DmnConfig,
) -> (Array[Array[Double]], Array[Array[Double]]) {
  let samples = counts.length()
  let taxa = counts[0].length()
  let components = config.components
  let proportions = dmn_make_matrix(samples, taxa, 0.0)
  for sample in 0..<samples {
    let total = dmn_row_total(counts[sample]).to_double()
    for taxon in 0..<taxa {
      proportions[sample][taxon] = counts[sample][taxon].to_double() / total
    }
  }

  let responsibilities = dmn_make_matrix(samples, components, 0.0)
  let mut state = config.seed % 2147483647
  for sample in 0..<samples {
    let component = if sample < components {
      sample
    } else {
      state = dmn_lcg_next(state)
      state % components
    }
    responsibilities[sample][component] = 1.0
  }

  let centers = dmn_make_matrix(components, taxa, 0.0)
  let mut iteration = 0
  let mut change = 1.0e300
  while iteration < config.soft_kmeans_iterations && change > config.tolerance {
    change = 0.0
    for component in 0..<components {
      let mut weight = 0.0
      for sample in 0..<samples {
        weight = weight + responsibilities[sample][component]
      }
      if weight <= 1.0e-12 {
        let selected = component % samples
        for taxon in 0..<taxa {
          centers[component][taxon] = proportions[selected][taxon]
        }
      } else {
        for taxon in 0..<taxa {
          let previous = centers[component][taxon]
          let mut value = 0.0
          for sample in 0..<samples {
            value = value +
              responsibilities[sample][component] * proportions[sample][taxon]
          }
          centers[component][taxon] = value / weight
          let difference = centers[component][taxon] - previous
          change = change + difference * difference
        }
      }
    }
    change = change.sqrt()

    for sample in 0..<samples {
      let scores : Array[Double] = []
      for component in 0..<components {
        let mut distance = 0.0
        for taxon in 0..<taxa {
          let difference = centers[component][taxon] -
            proportions[sample][taxon]
          distance = distance + difference * difference
        }
        scores.push(-config.soft_beta * distance.sqrt())
      }
      let normalizer = dmn_log_sum_exp(scores)
      for component in 0..<components {
        responsibilities[sample][component] = @math.exp(
          scores[component] - normalizer,
        )
      }
    }
    iteration = iteration + 1
  }
  (responsibilities, centers)
}

///|
fn dmn_alpha_objective_gradient(
  lambda : Array[Double],
  counts : Array[Array[Int]],
  membership : Array[Double],
  config : DmnConfig,
) -> (Double, Array[Double]) {
  let samples = counts.length()
  let taxa = lambda.length()
  let alpha : Array[Double] = []
  let mut concentration = 0.0
  let mut lambda_sum = 0.0
  for value in lambda {
    let bounded = value.max(-18.0).min(18.0)
    let parameter = @math.exp(bounded).max(config.min_alpha)
    alpha.push(parameter)
    concentration = concentration + parameter
    lambda_sum = lambda_sum + bounded
  }
  let mut weight = 0.0
  for value in membership {
    weight = weight + value
  }

  let mut objective = config.prior_rate * concentration -
    config.prior_shape * lambda_sum -
    weight * dmn_log_gamma(concentration)
  for parameter in alpha {
    objective = objective + weight * dmn_log_gamma(parameter)
  }
  for sample in 0..<samples {
    if membership[sample] > 0.0 {
      let total = dmn_row_total(counts[sample]).to_double()
      objective = objective +
        membership[sample] * dmn_log_gamma(concentration + total)
      for taxon in 0..<taxa {
        objective = objective -
          membership[sample] *
          dmn_log_gamma(alpha[taxon] + counts[sample][taxon].to_double())
      }
    }
  }

  let mut common = -weight * dmn_digamma(concentration)
  for sample in 0..<samples {
    common = common +
      membership[sample] *
      dmn_digamma(concentration + dmn_row_total(counts[sample]).to_double())
  }
  let gradient : Array[Double] = []
  for taxon in 0..<taxa {
    let mut value = weight * dmn_digamma(alpha[taxon]) + common
    for sample in 0..<samples {
      value = value -
        membership[sample] *
        dmn_digamma(alpha[taxon] + counts[sample][taxon].to_double())
    }
    gradient.push(
      alpha[taxon] * (config.prior_rate + value) - config.prior_shape,
    )
  }
  (objective, gradient)
}

///|
fn dmn_dot(left : Array[Double], right : Array[Double]) -> Double {
  let mut value = 0.0
  for index in 0..<left.length() {
    value = value + left[index] * right[index]
  }
  value
}

///|
fn dmn_identity(size : Int) -> Array[Array[Double]] {
  let output = dmn_make_matrix(size, size, 0.0)
  for index in 0..<size {
    output[index][index] = 1.0
  }
  output
}

///|
fn dmn_matrix_vector(
  matrix : Array[Array[Double]],
  vector : Array[Double],
) -> Array[Double] {
  let output = Array::make(matrix.length(), 0.0)
  for row in 0..<matrix.length() {
    output[row] = dmn_dot(matrix[row], vector)
  }
  output
}

///|
fn dmn_optimize_alpha(
  counts : Array[Array[Int]],
  membership : Array[Double],
  initial : Array[Double],
  config : DmnConfig,
) -> Array[Double] {
  let taxa = initial.length()
  let mut lambda = initial.map(fn(value) {
    @math.ln(value.max(config.min_alpha))
  })
  let mut inverse_hessian = dmn_identity(taxa)
  let (initial_objective, initial_gradient) = dmn_alpha_objective_gradient(
    lambda, counts, membership, config,
  )
  let mut objective = initial_objective
  let mut gradient = initial_gradient
  let mut iteration = 0

  while iteration < config.optimizer_iterations {
    let mut maximum_gradient = 0.0
    for value in gradient {
      maximum_gradient = maximum_gradient.max(value.abs())
    }
    if maximum_gradient <= config.optimizer_tolerance {
      break
    }

    let transformed = dmn_matrix_vector(inverse_hessian, gradient)
    let direction = transformed.map(fn(value) { -value })
    let mut directional_derivative = dmn_dot(gradient, direction)
    if !dmn_is_finite(directional_derivative) || directional_derivative >= 0.0 {
      let norm = dmn_dot(gradient, gradient).sqrt().max(1.0)
      for index in 0..<taxa {
        direction[index] = -gradient[index] / norm
      }
      directional_derivative = dmn_dot(gradient, direction)
      inverse_hessian = dmn_identity(taxa)
    }

    let mut step = 1.0
    let candidate = lambda.copy()
    let mut candidate_objective = objective
    let mut candidate_gradient = gradient.copy()
    let mut accepted = false
    for _ in 0..<40 {
      for index in 0..<taxa {
        candidate[index] = (lambda[index] + step * direction[index])
          .max(-18.0)
          .min(18.0)
      }
      let (next_objective, next_gradient) = dmn_alpha_objective_gradient(
        candidate, counts, membership, config,
      )
      if dmn_is_finite(next_objective) &&
        next_objective <= objective + 1.0e-4 * step * directional_derivative {
        candidate_objective = next_objective
        candidate_gradient = next_gradient
        accepted = true
        break
      }
      step = step * 0.5
    }
    if !accepted {
      break
    }

    let displacement : Array[Double] = []
    let gradient_change : Array[Double] = []
    for index in 0..<taxa {
      displacement.push(candidate[index] - lambda[index])
      gradient_change.push(candidate_gradient[index] - gradient[index])
    }
    let curvature = dmn_dot(gradient_change, displacement)
    if dmn_is_finite(curvature) && curvature > 1.0e-12 {
      let rho = 1.0 / curvature
      let hessian_gradient = dmn_matrix_vector(inverse_hessian, gradient_change)
      let quadratic = dmn_dot(gradient_change, hessian_gradient)
      for row in 0..<taxa {
        for column in 0..<taxa {
          inverse_hessian[row][column] = inverse_hessian[row][column] +
            (1.0 + quadratic * rho) *
            rho *
            displacement[row] *
            displacement[column] -
            rho *
            (
              displacement[row] * hessian_gradient[column] +
              hessian_gradient[row] * displacement[column]
            )
        }
      }
    } else {
      inverse_hessian = dmn_identity(taxa)
    }
    lambda = candidate.copy()
    objective = candidate_objective
    gradient = candidate_gradient
    iteration = iteration + 1
  }

  lambda.map(fn(value) {
    @math.exp(value.max(-18.0).min(18.0)).max(config.min_alpha)
  })
}

///|
fn dmn_expectation(
  counts : Array[Array[Int]],
  estimates : Array[Array[Double]],
  weights : Array[Double],
) -> (Array[Array[Double]], Double) {
  let samples = counts.length()
  let components = estimates.length()
  let responsibilities = dmn_make_matrix(samples, components, 0.0)
  let mut log_likelihood = 0.0
  for sample in 0..<samples {
    let scores : Array[Double] = []
    for component in 0..<components {
      scores.push(
        @math.ln(weights[component].max(1.0e-300)) +
        dmn_component_log_kernel(counts[sample], estimates[component]),
      )
    }
    let normalizer = dmn_log_sum_exp(scores)
    let mut coefficient = dmn_log_gamma(
      dmn_row_total(counts[sample]).to_double() + 1.0,
    )
    for count in counts[sample] {
      coefficient = coefficient - dmn_log_gamma(count.to_double() + 1.0)
    }
    log_likelihood = log_likelihood + normalizer + coefficient
    for component in 0..<components {
      responsibilities[sample][component] = @math.exp(
        scores[component] - normalizer,
      )
    }
  }
  (responsibilities, log_likelihood)
}

///|
fn dmn_negative_log_evidence(
  log_likelihood : Double,
  estimates : Array[Array[Double]],
  config : DmnConfig,
) -> Double {
  let components = estimates.length()
  let taxa = estimates[0].length()
  let mut alpha_total = 0.0
  let mut lambda_total = 0.0
  for component in estimates {
    for alpha in component {
      alpha_total = alpha_total + alpha
      lambda_total = lambda_total + @math.ln(alpha)
    }
  }
  -log_likelihood +
  (components * taxa).to_double() * dmn_log_gamma(config.prior_shape) -
  config.prior_shape *
  (components * taxa).to_double() *
  @math.ln(config.prior_rate) +
  config.prior_rate * alpha_total -
  config.prior_shape * lambda_total
}

///|
fn dmn_hessian(
  counts : Array[Array[Int]],
  membership : Array[Double],
  alpha : Array[Double],
  config : DmnConfig,
) -> Array[Array[Double]] {
  let taxa = alpha.length()
  let samples = counts.length()
  let hessian = dmn_make_matrix(taxa, taxa, 0.0)
  let mut concentration = 0.0
  let mut weight = 0.0
  for value in alpha {
    concentration = concentration + value
  }
  for value in membership {
    weight = weight + value
  }

  let mut common_digamma = -weight * dmn_digamma(concentration)
  let mut common_trigamma = -weight * dmn_trigamma(concentration)
  for sample in 0..<samples {
    let shifted = concentration + dmn_row_total(counts[sample]).to_double()
    common_digamma = common_digamma + membership[sample] * dmn_digamma(shifted)
    common_trigamma = common_trigamma +
      membership[sample] * dmn_trigamma(shifted)
  }

  let first : Array[Double] = []
  let second : Array[Double] = []
  for taxon in 0..<taxa {
    let mut first_value = common_digamma + weight * dmn_digamma(alpha[taxon])
    let mut second_value = common_trigamma + weight * dmn_trigamma(alpha[taxon])
    for sample in 0..<samples {
      let shifted = alpha[taxon] + counts[sample][taxon].to_double()
      first_value = first_value - membership[sample] * dmn_digamma(shifted)
      second_value = second_value - membership[sample] * dmn_trigamma(shifted)
    }
    first.push(first_value + config.prior_rate)
    second.push(second_value)
  }

  for row in 0..<taxa {
    for column in 0..<taxa {
      hessian[row][column] = if row == column {
        alpha[row] * first[row] + alpha[row] * alpha[row] * second[row]
      } else {
        alpha[row] * alpha[column] * common_trigamma
      }
    }
  }
  hessian
}

///|
fn dmn_log_abs_determinant(matrix : Array[Array[Double]]) -> Double {
  let size = matrix.length()
  let work = dmn_copy_matrix(matrix)
  for index in 0..<size {
    work[index][index] = work[index][index] + 1.0e-10
  }
  let mut result = 0.0
  for pivot in 0..<size {
    let mut selected = pivot
    let mut maximum = work[pivot][pivot].abs()
    for row in (pivot + 1)..<size {
      if work[row][pivot].abs() > maximum {
        maximum = work[row][pivot].abs()
        selected = row
      }
    }
    if maximum < 1.0e-14 {
      return 0.0
    }
    if selected != pivot {
      let swap = work[pivot]
      work[pivot] = work[selected]
      work[selected] = swap
    }
    let diagonal = work[pivot][pivot]
    result = result + @math.ln(diagonal.abs())
    for row in (pivot + 1)..<size {
      let factor = work[row][pivot] / diagonal
      for column in (pivot + 1)..<size {
        work[row][column] = work[row][column] - factor * work[pivot][column]
      }
    }
  }
  result
}

///|
fn dmn_inverse(matrix : Array[Array[Double]]) -> Array[Array[Double]]? {
  let size = matrix.length()
  let augmented = dmn_make_matrix(size, size * 2, 0.0)
  for row in 0..<size {
    for column in 0..<size {
      augmented[row][column] = matrix[row][column]
    }
    augmented[row][row] = augmented[row][row] + 1.0e-10
    augmented[row][size + row] = 1.0
  }
  for pivot in 0..<size {
    let mut selected = pivot
    let mut maximum = augmented[pivot][pivot].abs()
    for row in (pivot + 1)..<size {
      if augmented[row][pivot].abs() > maximum {
        maximum = augmented[row][pivot].abs()
        selected = row
      }
    }
    if maximum < 1.0e-14 {
      return None
    }
    if selected != pivot {
      let swap = augmented[pivot]
      augmented[pivot] = augmented[selected]
      augmented[selected] = swap
    }
    let diagonal = augmented[pivot][pivot]
    for column in 0..<(size * 2) {
      augmented[pivot][column] = augmented[pivot][column] / diagonal
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
  let inverse = dmn_make_matrix(size, size, 0.0)
  for row in 0..<size {
    for column in 0..<size {
      inverse[row][column] = augmented[row][size + column]
    }
  }
  Some(inverse)
}

///|
fn dmn_sort_components(
  estimates : Array[Array[Double]],
  weights : Array[Double],
  responsibilities : Array[Array[Double]],
) -> (Array[Array[Double]], Array[Double], Array[Array[Double]]) {
  let order : Array[Int] = []
  for component in 0..<weights.length() {
    order.push(component)
  }
  order.sort_by(fn(left : Int, right : Int) -> Int {
    if weights[left] > weights[right] {
      -1
    } else if weights[left] < weights[right] {
      1
    } else {
      left - right
    }
  })
  let sorted_estimates : Array[Array[Double]] = []
  let sorted_weights : Array[Double] = []
  for component in order {
    sorted_estimates.push(estimates[component].copy())
    sorted_weights.push(weights[component])
  }
  let sorted_responsibilities = dmn_make_matrix(
    responsibilities.length(),
    weights.length(),
    0.0,
  )
  for sample in 0..<responsibilities.length() {
    for component in 0..<order.length() {
      sorted_responsibilities[sample][component] = responsibilities[sample][order[component]]
    }
  }
  (sorted_estimates, sorted_weights, sorted_responsibilities)
}

///|
pub fn dirichlet_multinomial_fit(
  counts : Array[Array[Int]],
  sample_names? : Array[String] = [],
  taxon_names? : Array[String] = [],
  config? : DmnConfig = DmnConfig::default(),
) -> DmnFit raise DirichletMultinomialError {
  let (samples, taxa) = dmn_validate_counts(counts)
  if config.components > samples {
    raise DirichletMultinomialError(
      "DirichletMultinomial component count must not exceed sample count",
    )
  }
  let resolved_sample_names = if sample_names.length() == 0 {
    dmn_default_names("sample", samples)
  } else {
    dmn_validate_names(sample_names, samples, "sample")
    sample_names.copy()
  }
  let resolved_taxon_names = if taxon_names.length() == 0 {
    dmn_default_names("taxon", taxa)
  } else {
    dmn_validate_names(taxon_names, taxa, "taxon")
    taxon_names.copy()
  }

  let (initial_responsibilities, centers) = dmn_soft_kmeans(counts, config)
  let mut responsibilities = initial_responsibilities
  let estimates : Array[Array[Double]] = []
  for component in 0..<config.components {
    let membership = responsibilities.map(fn(row) { row[component] })
    estimates.push(
      dmn_optimize_alpha(counts, membership, centers[component], config),
    )
  }
  let weights = Array::make(config.components, 0.0)
  for component in 0..<config.components {
    for sample in 0..<samples {
      weights[component] = weights[component] +
        responsibilities[sample][component]
    }
    weights[component] = (weights[component] / samples.to_double()).max(1.0e-12)
  }
  let mut weight_total = 0.0
  for value in weights {
    weight_total = weight_total + value
  }
  for component in 0..<weights.length() {
    weights[component] = weights[component] / weight_total
  }

  let log_likelihoods : Array[Double] = []
  let mut previous = 1.0e300
  let mut converged = false
  let mut iteration = 0
  while iteration < config.max_iterations {
    let (next_responsibilities, _) = dmn_expectation(counts, estimates, weights)
    responsibilities = next_responsibilities
    for component in 0..<config.components {
      let membership = responsibilities.map(fn(row) { row[component] })
      estimates[component] = dmn_optimize_alpha(
        counts,
        membership,
        estimates[component],
        config,
      )
      let mut component_weight = 0.0
      for value in membership {
        component_weight = component_weight + value
      }
      weights[component] = (component_weight / samples.to_double()).max(1.0e-12)
    }
    weight_total = 0.0
    for value in weights {
      weight_total = weight_total + value
    }
    for component in 0..<weights.length() {
      weights[component] = weights[component] / weight_total
    }
    let (_, updated_log_likelihood) = dmn_expectation(
      counts, estimates, weights,
    )
    let negative_evidence = dmn_negative_log_evidence(
      updated_log_likelihood, estimates, config,
    )
    log_likelihoods.push(updated_log_likelihood)
    let change = (previous - negative_evidence).abs()
    if iteration > 0 && change <= config.tolerance * (1.0 + previous.abs()) {
      converged = true
      iteration = iteration + 1
      break
    }
    previous = negative_evidence
    iteration = iteration + 1
  }

  let (final_responsibilities, final_log_likelihood) = dmn_expectation(
    counts, estimates, weights,
  )
  let (sorted_estimates, sorted_weights, sorted_responsibilities) = dmn_sort_components(
    estimates, weights, final_responsibilities,
  )
  let lower = dmn_make_matrix(config.components, taxa, 0.0)
  let upper = dmn_make_matrix(config.components, taxa, 0.0)
  let mut log_determinant = 0.0
  for component in 0..<config.components {
    let membership = sorted_responsibilities.map(fn(row) { row[component] })
    let hessian = dmn_hessian(
      counts,
      membership,
      sorted_estimates[component],
      config,
    )
    log_determinant = log_determinant + dmn_log_abs_determinant(hessian)
    if component > 0 {
      log_determinant = log_determinant +
        @math.ln(samples.to_double()) -
        @math.ln(sorted_weights[component].max(1.0e-300))
    }
    match dmn_inverse(hessian) {
      Some(covariance) =>
        for taxon in 0..<taxa {
          let standard_error = covariance[taxon][taxon].abs().sqrt().min(100.0)
          let center = @math.ln(sorted_estimates[component][taxon])
          lower[component][taxon] = @math.exp(center - 1.96 * standard_error)
          upper[component][taxon] = @math.exp(center + 1.96 * standard_error)
        }
      None =>
        for taxon in 0..<taxa {
          lower[component][taxon] = sorted_estimates[component][taxon]
          upper[component][taxon] = sorted_estimates[component][taxon]
        }
    }
  }
  let negative_log_evidence = dmn_negative_log_evidence(
    final_log_likelihood, sorted_estimates, config,
  )
  let parameters = (config.components * taxa + config.components - 1).to_double()
  let goodness = DmnGoodnessOfFit::{
    negative_log_evidence,
    log_determinant,
    laplace: negative_log_evidence +
    0.5 * log_determinant -
    0.5 * parameters * @math.ln(2.0 * 3.141592653589793),
    bic: negative_log_evidence +
    0.5 * @math.ln(samples.to_double()) * parameters,
    aic: negative_log_evidence + parameters,
  }
  DmnFit::{
    estimates: sorted_estimates,
    lower,
    upper,
    weights: sorted_weights,
    responsibilities: sorted_responsibilities,
    sample_names: resolved_sample_names,
    taxon_names: resolved_taxon_names,
    goodness,
    log_likelihoods,
    converged,
    iterations: iteration,
    config,
  }
}

///|
pub fn DmnFit::component_count(self : DmnFit) -> Int {
  self.estimates.length()
}

///|
pub fn DmnFit::sample_count(self : DmnFit) -> Int {
  self.responsibilities.length()
}

///|
pub fn DmnFit::taxon_count(self : DmnFit) -> Int {
  self.taxon_names.length()
}

///|
pub fn DmnFit::mixture(self : DmnFit) -> Array[Array[Double]] {
  dmn_copy_matrix(self.responsibilities)
}

///|
pub fn DmnFit::mixture_weights(self : DmnFit) -> Array[Double] {
  self.weights.copy()
}

///|
pub fn DmnFit::concentrations(self : DmnFit) -> Array[Double] {
  self.estimates.map(fn(component) {
    let mut total = 0.0
    for value in component {
      total = total + value
    }
    total
  })
}

///|
pub fn DmnFit::component_proportions(self : DmnFit) -> Array[Array[Double]] {
  let output = dmn_copy_matrix(self.estimates)
  for component in 0..<output.length() {
    let mut total = 0.0
    for value in output[component] {
      total = total + value
    }
    for taxon in 0..<output[component].length() {
      output[component][taxon] = output[component][taxon] / total
    }
  }
  output
}

///|
pub fn DmnFit::fitted(
  self : DmnFit,
  scale? : Bool = false,
) -> Array[Array[Double]] {
  let taxa = self.taxon_count()
  let components = self.component_count()
  let output = dmn_make_matrix(taxa, components, 0.0)
  for component in 0..<components {
    let mut concentration = 0.0
    for value in self.estimates[component] {
      concentration = concentration + value
    }
    for taxon in 0..<taxa {
      output[taxon][component] = if scale {
        self.estimates[component][taxon] / concentration
      } else {
        self.estimates[component][taxon]
      }
    }
  }
  output
}

///|
pub fn DmnFit::assignments(self : DmnFit) -> Array[Int] {
  self.responsibilities.map(fn(row) {
    let mut selected = 0
    for component in 1..<row.length() {
      if row[component] > row[selected] {
        selected = component
      }
    }
    selected
  })
}

///|
pub fn DmnFit::negative_log_evidence(
  self : DmnFit,
  counts : Array[Array[Int]],
) -> Array[Array[Double]] raise DirichletMultinomialError {
  let (_, taxa) = dmn_validate_counts(counts)
  if taxa != self.taxon_count() {
    raise DirichletMultinomialError(
      "DirichletMultinomial prediction taxon count does not match the model",
    )
  }
  let output = dmn_make_matrix(counts.length(), self.component_count(), 0.0)
  for sample in 0..<counts.length() {
    for component in 0..<self.component_count() {
      output[sample][component] = -dmn_component_log_kernel(
        counts[sample],
        self.estimates[component],
      )
    }
  }
  output
}

///|
pub fn DmnFit::predict(
  self : DmnFit,
  counts : Array[Array[Int]],
) -> Array[Array[Double]] raise DirichletMultinomialError {
  let negative_evidence = self.negative_log_evidence(counts)
  let output = dmn_make_matrix(counts.length(), self.component_count(), 0.0)
  for sample in 0..<counts.length() {
    let scores : Array[Double] = []
    for component in 0..<self.component_count() {
      scores.push(
        @math.ln(self.weights[component].max(1.0e-300)) -
        negative_evidence[sample][component],
      )
    }
    let normalizer = dmn_log_sum_exp(scores)
    for component in 0..<self.component_count() {
      output[sample][component] = @math.exp(scores[component] - normalizer)
    }
  }
  output
}

///|
pub fn DmnFit::predict_assignments(
  self : DmnFit,
  counts : Array[Array[Int]],
) -> Array[Int] raise DirichletMultinomialError {
  self
  .predict(counts)
  .map(fn(row) {
    let mut selected = 0
    for component in 1..<row.length() {
      if row[component] > row[selected] {
        selected = component
      }
    }
    selected
  })
}

///|
pub fn DmnFit::summary(self : DmnFit) -> String {
  "DMN fit: k=" +
  self.component_count().to_string() +
  ", samples=" +
  self.sample_count().to_string() +
  ", taxa=" +
  self.taxon_count().to_string() +
  ", iterations=" +
  self.iterations.to_string() +
  ", converged=" +
  self.converged.to_string() +
  ", Laplace=" +
  self.goodness.laplace.to_string() +
  ", BIC=" +
  self.goodness.bic.to_string() +
  ", AIC=" +
  self.goodness.aic.to_string()
}

///|
fn dmn_criterion_score(fit : DmnFit, criterion : DmnCriterion) -> Double {
  match criterion {
    DmnLaplace => fit.goodness.laplace
    DmnBic => fit.goodness.bic
    DmnAic => fit.goodness.aic
  }
}

///|
pub fn dirichlet_multinomial_select(
  counts : Array[Array[Int]],
  max_components : Int,
  criterion? : DmnCriterion = DmnLaplace,
  sample_names? : Array[String] = [],
  taxon_names? : Array[String] = [],
  config? : DmnConfig = DmnConfig::default(),
) -> DmnModelSelection raise DirichletMultinomialError {
  let (samples, _) = dmn_validate_counts(counts)
  if max_components < 1 || max_components > samples {
    raise DirichletMultinomialError(
      "DirichletMultinomial model-selection range must be within sample count",
    )
  }
  let fits : Array[DmnFit] = []
  let scores : Array[Double] = []
  let mut best_index = 0
  for components in 1..<=max_components {
    let fit_config = config.with_components(
      components,
      seed=config.seed + components * 7919,
    )
    let fit = dirichlet_multinomial_fit(
      counts,
      sample_names~,
      taxon_names~,
      config=fit_config,
    )
    fits.push(fit)
    let score = dmn_criterion_score(fit, criterion)
    scores.push(score)
    if score < scores[best_index] {
      best_index = components - 1
    }
  }
  DmnModelSelection::{ fits, scores, best_index, criterion }
}

///|
pub fn DmnModelSelection::best_fit(self : DmnModelSelection) -> DmnFit {
  self.fits[self.best_index]
}

///|
pub fn DmnModelSelection::best_component_count(self : DmnModelSelection) -> Int {
  self.best_index + 1
}

///|
fn dmn_unique_groups(groups : Array[String]) -> Array[String] {
  let output : Array[String] = []
  let seen : Map[String, Bool] = Map([])
  for group in groups {
    if !seen.contains(group) {
      seen[group] = true
      output.push(group)
    }
  }
  output
}

///|
pub fn dirichlet_multinomial_group_fit(
  counts : Array[Array[Int]],
  groups : Array[String],
  components_by_group? : Array[Int] = [],
  max_components? : Int = 3,
  taxon_names? : Array[String] = [],
  config? : DmnConfig = DmnConfig::default(),
) -> DmnGroupFit raise DirichletMultinomialError {
  let (samples, taxa) = dmn_validate_counts(counts)
  if groups.length() != samples {
    raise DirichletMultinomialError(
      "DirichletMultinomial group label count does not match samples",
    )
  }
  for group in groups {
    if group.trim().length() == 0 {
      raise DirichletMultinomialError(
        "DirichletMultinomial group labels must not be empty",
      )
    }
  }
  let group_names = dmn_unique_groups(groups)
  if group_names.length() < 2 {
    raise DirichletMultinomialError(
      "DirichletMultinomial group fitting requires at least two groups",
    )
  }
  if components_by_group.length() > 0 &&
    components_by_group.length() != group_names.length() {
    raise DirichletMultinomialError(
      "DirichletMultinomial group component counts must match group count",
    )
  }
  if components_by_group.length() == 0 && max_components < 1 {
    raise DirichletMultinomialError(
      "DirichletMultinomial maximum group component count must be positive",
    )
  }
  let resolved_taxon_names = if taxon_names.length() == 0 {
    dmn_default_names("taxon", taxa)
  } else {
    dmn_validate_names(taxon_names, taxa, "taxon")
    taxon_names.copy()
  }
  let models : Array[DmnFit] = []
  let priors : Array[Double] = []
  for group_index in 0..<group_names.length() {
    let group_counts : Array[Array[Int]] = []
    let group_samples : Array[String] = []
    for sample in 0..<samples {
      if groups[sample] == group_names[group_index] {
        group_counts.push(counts[sample].copy())
        group_samples.push("sample" + (sample + 1).to_string())
      }
    }
    priors.push(group_counts.length().to_double() / samples.to_double())
    let model = if components_by_group.length() > 0 {
      if components_by_group[group_index] < 1 ||
        components_by_group[group_index] > group_counts.length() {
        raise DirichletMultinomialError(
          "DirichletMultinomial group component count exceeds group sample count",
        )
      }
      dirichlet_multinomial_fit(
        group_counts,
        sample_names=group_samples,
        taxon_names=resolved_taxon_names,
        config=config.with_components(
          components_by_group[group_index],
          seed=config.seed + group_index * 3571 + 1,
        ),
      )
    } else {
      let maximum = max_components.min(group_counts.length())
      dirichlet_multinomial_select(
        group_counts,
        maximum,
        sample_names=group_samples,
        taxon_names=resolved_taxon_names,
        config=config.with_components(
          1,
          seed=config.seed + group_index * 3571 + 1,
        ),
      ).best_fit()
    }
    models.push(model)
  }
  DmnGroupFit::{
    group_names,
    models,
    priors,
    taxon_names: resolved_taxon_names,
  }
}

///|
pub fn DmnGroupFit::predict(
  self : DmnGroupFit,
  counts : Array[Array[Int]],
) -> Array[Array[Double]] raise DirichletMultinomialError {
  let (_, taxa) = dmn_validate_counts(counts)
  if taxa != self.taxon_names.length() {
    raise DirichletMultinomialError(
      "DirichletMultinomial group prediction taxon count does not match models",
    )
  }
  let output = dmn_make_matrix(counts.length(), self.group_names.length(), 0.0)
  for sample in 0..<counts.length() {
    let group_scores : Array[Double] = []
    for group in 0..<self.group_names.length() {
      let component_scores : Array[Double] = []
      for component in 0..<self.models[group].component_count() {
        component_scores.push(
          @math.ln(self.models[group].weights[component].max(1.0e-300)) +
          dmn_component_log_kernel(
            counts[sample],
            self.models[group].estimates[component],
          ),
        )
      }
      group_scores.push(
        @math.ln(self.priors[group].max(1.0e-300)) +
        dmn_log_sum_exp(component_scores),
      )
    }
    let normalizer = dmn_log_sum_exp(group_scores)
    for group in 0..<self.group_names.length() {
      output[sample][group] = @math.exp(group_scores[group] - normalizer)
    }
  }
  output
}

///|
pub fn DmnGroupFit::predict_assignments(
  self : DmnGroupFit,
  counts : Array[Array[Int]],
) -> Array[String] raise DirichletMultinomialError {
  self
  .predict(counts)
  .map(fn(row) {
    let mut selected = 0
    for group in 1..<row.length() {
      if row[group] > row[selected] {
        selected = group
      }
    }
    self.group_names[selected]
  })
}

///|
pub fn dirichlet_multinomial_cross_validate(
  counts : Array[Array[Int]],
  groups : Array[String],
  folds : Int,
  components_by_group? : Array[Int] = [],
  taxon_names? : Array[String] = [],
  config? : DmnConfig = DmnConfig::default(),
) -> DmnCrossValidation raise DirichletMultinomialError {
  let (samples, _) = dmn_validate_counts(counts)
  if groups.length() != samples {
    raise DirichletMultinomialError(
      "DirichletMultinomial cross-validation group count does not match samples",
    )
  }
  let group_names = dmn_unique_groups(groups)
  if folds < 2 {
    raise DirichletMultinomialError(
      "DirichletMultinomial cross-validation requires at least two folds",
    )
  }
  let group_seen : Map[String, Int] = Map([])
  let group_sizes : Map[String, Int] = Map([])
  for group in groups {
    group_sizes[group] = match group_sizes.get(group) {
      Some(value) => value + 1
      None => 1
    }
  }
  for group in group_names {
    if group_sizes[group] < folds {
      raise DirichletMultinomialError(
        "DirichletMultinomial each group must contain at least one sample per fold",
      )
    }
  }
  let fold_ids = Array::make(samples, 0)
  for sample in 0..<samples {
    let seen = match group_seen.get(groups[sample]) {
      Some(value) => value
      None => 0
    }
    fold_ids[sample] = seen % folds
    group_seen[groups[sample]] = seen + 1
  }

  let probabilities = dmn_make_matrix(samples, group_names.length(), 0.0)
  let predictions = Array::make(samples, "")
  for fold in 0..<folds {
    let training_counts : Array[Array[Int]] = []
    let training_groups : Array[String] = []
    let testing_counts : Array[Array[Int]] = []
    let testing_indices : Array[Int] = []
    for sample in 0..<samples {
      if fold_ids[sample] == fold {
        testing_counts.push(counts[sample].copy())
        testing_indices.push(sample)
      } else {
        training_counts.push(counts[sample].copy())
        training_groups.push(groups[sample])
      }
    }
    let classifier = dirichlet_multinomial_group_fit(
      training_counts,
      training_groups,
      components_by_group~,
      taxon_names~,
      config=config.with_components(1, seed=config.seed + fold * 104729 + 1),
    )
    let fold_probabilities = classifier.predict(testing_counts)
    let fold_predictions = classifier.predict_assignments(testing_counts)
    for index in 0..<testing_indices.length() {
      let original = testing_indices[index]
      for group in 0..<group_names.length() {
        match classifier.group_names.search(group_names[group]) {
          Some(classifier_group) =>
            probabilities[original][group] = fold_probabilities[index][classifier_group]
          None => ()
        }
      }
      predictions[original] = fold_predictions[index]
    }
  }
  let mut correct = 0
  for sample in 0..<samples {
    if predictions[sample] == groups[sample] {
      correct = correct + 1
    }
  }
  DmnCrossValidation::{
    probabilities,
    predictions,
    truth: groups.copy(),
    group_names,
    fold_ids,
    correct,
    accuracy: correct.to_double() / samples.to_double(),
  }
}

///|
pub fn dirichlet_multinomial_roc(
  truth : Array[Bool],
  scores : Array[Double],
) -> DmnRoc raise DirichletMultinomialError {
  if truth.length() == 0 || truth.length() != scores.length() {
    raise DirichletMultinomialError(
      "DirichletMultinomial ROC truth and score vectors must have equal positive length",
    )
  }
  let mut positives = 0
  let mut negatives = 0
  let mut maximum = scores[0]
  let mut minimum = scores[0]
  for index in 0..<scores.length() {
    if !dmn_is_finite(scores[index]) {
      raise DirichletMultinomialError(
        "DirichletMultinomial ROC scores must be finite",
      )
    }
    if truth[index] {
      positives = positives + 1
    } else {
      negatives = negatives + 1
    }
    maximum = maximum.max(scores[index])
    minimum = minimum.min(scores[index])
  }
  if positives == 0 || negatives == 0 {
    raise DirichletMultinomialError(
      "DirichletMultinomial ROC requires both positive and negative samples",
    )
  }
  let sorted_scores = scores.copy()
  sorted_scores.sort_by(fn(left : Double, right : Double) -> Int {
    if left > right {
      -1
    } else if left < right {
      1
    } else {
      0
    }
  })
  let thresholds : Array[Double] = [maximum + 1.0e-12]
  for value in sorted_scores {
    if (value - thresholds[thresholds.length() - 1]).abs() > 1.0e-12 {
      thresholds.push(value)
    }
  }
  thresholds.push(minimum - 1.0e-12)
  let sensitivity : Array[Double] = []
  let specificity : Array[Double] = []
  let mut best_threshold = thresholds[0]
  let mut best_youden = -1.0e300
  for threshold in thresholds {
    let mut true_positive = 0
    let mut false_positive = 0
    for index in 0..<scores.length() {
      if scores[index] >= threshold {
        if truth[index] {
          true_positive = true_positive + 1
        } else {
          false_positive = false_positive + 1
        }
      }
    }
    let current_sensitivity = true_positive.to_double() / positives.to_double()
    let current_specificity = 1.0 -
      false_positive.to_double() / negatives.to_double()
    sensitivity.push(current_sensitivity)
    specificity.push(current_specificity)
    let youden = current_sensitivity + current_specificity - 1.0
    if youden > best_youden {
      best_youden = youden
      best_threshold = threshold
    }
  }
  let mut concordance = 0.0
  for positive in 0..<scores.length() {
    if truth[positive] {
      for negative in 0..<scores.length() {
        if !truth[negative] {
          if scores[positive] > scores[negative] {
            concordance = concordance + 1.0
          } else if scores[positive] == scores[negative] {
            concordance = concordance + 0.5
          }
        }
      }
    }
  }
  DmnRoc::{
    thresholds,
    sensitivity,
    specificity,
    auc: concordance / (positives * negatives).to_double(),
    best_threshold,
  }
}

///|
pub fn dirichlet_multinomial_fit_se(
  experiment : SummarizedExperiment,
  assay_name : String,
  sample_names? : Array[String] = [],
  taxon_names? : Array[String] = [],
  config? : DmnConfig = DmnConfig::default(),
) -> DmnFit raise DirichletMultinomialError {
  let assay = match se_assay(experiment, assay_name) {
    Some(value) => value
    None =>
      raise DirichletMultinomialError(
        "DirichletMultinomial assay was not found in SummarizedExperiment",
      )
  }
  if assay.length() == 0 || assay[0].length() == 0 {
    raise DirichletMultinomialError(
      "DirichletMultinomial assay must contain taxa and samples",
    )
  }
  let taxa = assay.length()
  let samples = assay[0].length()
  let counts : Array[Array[Int]] = []
  for sample in 0..<samples {
    let row : Array[Int] = []
    for taxon in 0..<taxa {
      if assay[taxon].length() != samples {
        raise DirichletMultinomialError(
          "DirichletMultinomial SummarizedExperiment assay must be rectangular",
        )
      }
      let value = assay[taxon][sample]
      if !dmn_is_finite(value) ||
        value < 0.0 ||
        (value - value.round()).abs() > 1.0e-8 {
        raise DirichletMultinomialError(
          "DirichletMultinomial assay must contain non-negative integer counts",
        )
      }
      row.push(value.round().to_int())
    }
    counts.push(row)
  }
  dirichlet_multinomial_fit(counts, sample_names~, taxon_names~, config~)
}
