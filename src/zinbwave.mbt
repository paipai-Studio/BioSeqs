// Zero-inflated negative binomial factor models for single-cell counts.
//
// This module follows the model used by Bioconductor zinbwave:
//   log(mu)   = X beta_mu + (V gamma_mu)^T + W alpha_mu + O_mu
//   logit(pi) = X beta_pi + (V gamma_pi)^T + W alpha_pi + O_pi
// where W contains inferred cell-level factors. Count matrices use the
// BioSeqs convention of genes in rows and cells in columns.

///|
pub suberror ZinbWaveError {
  ZinbWaveError(String)
}

///|
pub struct ZinbWaveConfig {
  n_factors : Int
  max_iterations : Int
  inner_iterations : Int
  tolerance : Double
  ridge : Double
  dispersion_shrinkage : Double
  minimum_mean : Double
  minimum_probability : Double
} derive(Eq, Debug)
pub extend ZinbWaveConfig with Eq::{equal, not_equal}
pub extend ZinbWaveConfig with @moonbitlang/core/debug.Debug::{to_repr}

///|
pub fn ZinbWaveConfig::create(
  n_factors? : Int = 2,
  max_iterations? : Int = 50,
  inner_iterations? : Int = 2,
  tolerance? : Double = 1.0e-5,
  ridge? : Double = 1.0e-3,
  dispersion_shrinkage? : Double = 5.0,
  minimum_mean? : Double = 1.0e-8,
  minimum_probability? : Double = 1.0e-6,
) -> ZinbWaveConfig raise ZinbWaveError {
  if n_factors < 0 {
    raise ZinbWaveError("zinbwave factor count must be non-negative")
  }
  if max_iterations < 1 || inner_iterations < 1 {
    raise ZinbWaveError("zinbwave iteration counts must be positive")
  }
  if !zinbwave_is_finite(tolerance) || tolerance <= 0.0 {
    raise ZinbWaveError("zinbwave tolerance must be finite and positive")
  }
  if !zinbwave_is_finite(ridge) || ridge < 0.0 {
    raise ZinbWaveError(
      "zinbwave ridge penalty must be finite and non-negative",
    )
  }
  if !zinbwave_is_finite(dispersion_shrinkage) || dispersion_shrinkage < 0.0 {
    raise ZinbWaveError(
      "zinbwave dispersion shrinkage must be finite and non-negative",
    )
  }
  if !zinbwave_is_finite(minimum_mean) || minimum_mean <= 0.0 {
    raise ZinbWaveError("zinbwave minimum mean must be finite and positive")
  }
  if !zinbwave_is_finite(minimum_probability) ||
    minimum_probability <= 0.0 ||
    minimum_probability >= 0.5 {
    raise ZinbWaveError(
      "zinbwave minimum probability must be finite and in (0, 0.5)",
    )
  }
  ZinbWaveConfig::{
    n_factors,
    max_iterations,
    inner_iterations,
    tolerance,
    ridge,
    dispersion_shrinkage,
    minimum_mean,
    minimum_probability,
  }
}

///|
pub fn ZinbWaveConfig::default() -> ZinbWaveConfig {
  ZinbWaveConfig::{
    n_factors: 2,
    max_iterations: 50,
    inner_iterations: 2,
    tolerance: 1.0e-5,
    ridge: 1.0e-3,
    dispersion_shrinkage: 5.0,
    minimum_mean: 1.0e-8,
    minimum_probability: 1.0e-6,
  }
}

///|
pub struct ZinbWaveDiagnostics {
  iterations : Int
  converged : Bool
  initial_log_likelihood : Double
  final_log_likelihood : Double
  log_likelihoods : Array[Double]
  maximum_parameter_change : Double
} derive(Eq, Debug)
pub extend ZinbWaveDiagnostics with Eq::{equal, not_equal}
pub extend ZinbWaveDiagnostics with @moonbitlang/core/debug.Debug::{to_repr}

///|
pub struct ZinbWaveGeneSummary {
  gene_name : String
  dispersion : Double
  inverse_dispersion : Double
  observed_zero_fraction : Double
  fitted_zero_fraction : Double
  mean_count : Double
  mean_fitted : Double
} derive(Eq, Debug)
pub extend ZinbWaveGeneSummary with Eq::{equal, not_equal}
pub extend ZinbWaveGeneSummary with @moonbitlang/core/debug.Debug::{to_repr}

///|
pub struct ZinbWaveModel {
  gene_names : Array[String]
  cell_names : Array[String]
  factors : Array[Array[Double]]
  cell_design : Array[Array[Double]]
  gene_design : Array[Array[Double]]
  beta_mu : Array[Array[Double]]
  beta_pi : Array[Array[Double]]
  gamma_mu : Array[Array[Double]]
  gamma_pi : Array[Array[Double]]
  mean_offsets : Array[Array[Double]]
  zero_offsets : Array[Array[Double]]
  dispersions : Array[Double]
  means : Array[Array[Double]]
  zero_probabilities : Array[Array[Double]]
  observational_weights : Array[Array[Double]]
  deviance_residuals : Array[Array[Double]]
  log_likelihood : Double
  diagnostics : ZinbWaveDiagnostics
  config : ZinbWaveConfig
}

///|
pub struct ZinbWaveSCEOutput {
  experiment : SingleCellExperiment
  model : ZinbWaveModel
}

///|
fn zinbwave_is_finite(value : Double) -> Bool {
  value == value && value.abs() <= 1.0e300
}

///|
fn zinbwave_clamp(value : Double, lower : Double, upper : Double) -> Double {
  value.max(lower).min(upper)
}

///|
fn zinbwave_sigmoid(value : Double, minimum : Double) -> Double {
  let probability = if value >= 0.0 {
    let inverse = @math.exp(-value.min(40.0))
    1.0 / (1.0 + inverse)
  } else {
    let exponential = @math.exp(value.max(-40.0))
    exponential / (1.0 + exponential)
  }
  zinbwave_clamp(probability, minimum, 1.0 - minimum)
}

///|
fn zinbwave_logit(probability : Double, minimum : Double) -> Double {
  let bounded = zinbwave_clamp(probability, minimum, 1.0 - minimum)
  @math.ln(bounded / (1.0 - bounded))
}

///|
fn zinbwave_zero_matrix(rows : Int, columns : Int) -> Array[Array[Double]] {
  let matrix : Array[Array[Double]] = []
  for _ in 0..<rows {
    matrix.push(Array::make(columns, 0.0))
  }
  matrix
}

///|
fn zinbwave_copy_matrix(matrix : Array[Array[Double]]) -> Array[Array[Double]] {
  let copy : Array[Array[Double]] = []
  for row in matrix {
    copy.push(row.copy())
  }
  copy
}

///|
fn zinbwave_validate_counts(
  counts : Array[Array[Double]],
) -> (Int, Int) raise ZinbWaveError {
  if counts.length() == 0 {
    raise ZinbWaveError("zinbwave counts must contain at least one gene")
  }
  let cells = counts[0].length()
  if cells < 2 {
    raise ZinbWaveError("zinbwave counts must contain at least two cells")
  }
  for row in counts {
    if row.length() != cells {
      raise ZinbWaveError("zinbwave counts must be rectangular")
    }
    for value in row {
      if !zinbwave_is_finite(value) || value < 0.0 {
        raise ZinbWaveError(
          "zinbwave counts must contain finite non-negative values",
        )
      }
      if (value - value.round()).abs() > 1.0e-8 {
        raise ZinbWaveError("zinbwave requires integer-valued counts")
      }
    }
  }
  (counts.length(), cells)
}

///|
fn zinbwave_validate_covariates(
  covariates : Array[Array[Double]],
  expected_rows : Int,
  label : String,
) -> Int raise ZinbWaveError {
  if covariates.length() == 0 {
    return 0
  }
  if covariates.length() != expected_rows {
    raise ZinbWaveError(label + " rows must match the associated dimension")
  }
  let columns = covariates[0].length()
  if columns == 0 {
    raise ZinbWaveError(label + " must contain at least one column")
  }
  for row in covariates {
    if row.length() != columns {
      raise ZinbWaveError(label + " must be rectangular")
    }
    for value in row {
      if !zinbwave_is_finite(value) {
        raise ZinbWaveError(label + " must contain finite values")
      }
    }
  }
  columns
}

///|
fn zinbwave_design_with_intercept(
  covariates : Array[Array[Double]],
  rows : Int,
) -> Array[Array[Double]] {
  let columns = if covariates.length() == 0 {
    0
  } else {
    covariates[0].length()
  }
  let design : Array[Array[Double]] = []
  for row in 0..<rows {
    let values = Array::make(columns + 1, 0.0)
    values[0] = 1.0
    for column in 0..<columns {
      values[column + 1] = covariates[row][column]
    }
    design.push(values)
  }
  design
}

///|
fn zinbwave_validate_offsets(
  offsets : Array[Array[Double]],
  genes : Int,
  cells : Int,
  label : String,
) -> Array[Array[Double]] raise ZinbWaveError {
  if offsets.length() == 0 {
    return []
  }
  if offsets.length() != genes {
    raise ZinbWaveError(label + " rows must match the gene count")
  }
  let result : Array[Array[Double]] = []
  for row in offsets {
    if row.length() != cells {
      raise ZinbWaveError(label + " columns must match the cell count")
    }
    let values : Array[Double] = []
    for value in row {
      if !zinbwave_is_finite(value) {
        raise ZinbWaveError(label + " must contain finite values")
      }
      values.push(value)
    }
    result.push(values)
  }
  result
}

///|
fn zinbwave_default_mean_offsets(
  counts : Array[Array[Double]],
) -> Array[Array[Double]] raise ZinbWaveError {
  let genes = counts.length()
  let cells = counts[0].length()
  let library_sizes = Array::make(cells, 0.0)
  for gene in 0..<genes {
    for cell in 0..<cells {
      library_sizes[cell] = library_sizes[cell] + counts[gene][cell]
    }
  }
  let mut mean_log = 0.0
  for size in library_sizes {
    if size <= 0.0 {
      raise ZinbWaveError(
        "zinbwave cannot estimate offsets for an empty cell library",
      )
    }
    mean_log = mean_log + @math.ln(size)
  }
  mean_log = mean_log / cells.to_double()
  let offsets = zinbwave_zero_matrix(genes, cells)
  for gene in 0..<genes {
    for cell in 0..<cells {
      offsets[gene][cell] = @math.ln(library_sizes[cell]) - mean_log
    }
  }
  offsets
}

///|
fn zinbwave_prepare_names(
  supplied : Array[String],
  length : Int,
  prefix : String,
) -> Array[String] raise ZinbWaveError {
  if supplied.length() != 0 && supplied.length() != length {
    raise ZinbWaveError(prefix + " names must match the matrix dimension")
  }
  let names : Array[String] = []
  for index in 0..<length {
    let name = if supplied.length() == 0 {
      prefix + (index + 1).to_string()
    } else {
      supplied[index]
    }
    if name.trim().length() == 0 {
      raise ZinbWaveError(prefix + " names cannot be empty")
    }
    for previous in names {
      if previous == name {
        raise ZinbWaveError(prefix + " names must be unique: " + name)
      }
    }
    names.push(name)
  }
  names
}

///|
fn zinbwave_solve(
  matrix : Array[Array[Double]],
  right_hand_side : Array[Double],
) -> Array[Double] {
  let size = right_hand_side.length()
  if size == 0 {
    return []
  }
  let augmented : Array[Array[Double]] = []
  for row in 0..<size {
    let values = Array::make(size + 1, 0.0)
    for column in 0..<size {
      values[column] = matrix[row][column]
    }
    values[size] = right_hand_side[row]
    augmented.push(values)
  }
  for column in 0..<size {
    let mut pivot_row = column
    let mut pivot_size = augmented[column][column].abs()
    for row in (column + 1)..<size {
      if augmented[row][column].abs() > pivot_size {
        pivot_row = row
        pivot_size = augmented[row][column].abs()
      }
    }
    if pivot_row != column {
      let temporary = augmented[column]
      augmented[column] = augmented[pivot_row]
      augmented[pivot_row] = temporary
    }
    if augmented[column][column].abs() < 1.0e-10 {
      augmented[column][column] = if augmented[column][column] < 0.0 {
        -1.0e-10
      } else {
        1.0e-10
      }
    }
    for row in (column + 1)..<size {
      let ratio = augmented[row][column] / augmented[column][column]
      for item in column..<(size + 1) {
        augmented[row][item] = augmented[row][item] -
          ratio * augmented[column][item]
      }
    }
  }
  let solution = Array::make(size, 0.0)
  let mut row = size - 1
  while row >= 0 {
    let mut value = augmented[row][size]
    for column in (row + 1)..<size {
      value = value - augmented[row][column] * solution[column]
    }
    solution[row] = value / augmented[row][row]
    row = row - 1
  }
  solution
}

///|
fn zinbwave_weighted_ridge(
  design : Array[Array[Double]],
  response : Array[Double],
  weights : Array[Double],
  ridge : Double,
  penalize_intercept : Bool,
) -> Array[Double] {
  let observations = design.length()
  let parameters = design[0].length()
  let information = zinbwave_zero_matrix(parameters, parameters)
  let score = Array::make(parameters, 0.0)
  for observation in 0..<observations {
    let weight = weights[observation].max(1.0e-12)
    for left in 0..<parameters {
      score[left] = score[left] +
        design[observation][left] * weight * response[observation]
      for right in 0..<parameters {
        information[left][right] = information[left][right] +
          design[observation][left] * weight * design[observation][right]
      }
    }
  }
  for parameter in 0..<parameters {
    information[parameter][parameter] = information[parameter][parameter] +
      1.0e-8
    if penalize_intercept || parameter > 0 {
      information[parameter][parameter] = information[parameter][parameter] +
        ridge
    }
  }
  zinbwave_solve(information, score)
}

///|
fn zinbwave_initialize_factors(
  counts : Array[Array[Double]],
  offsets : Array[Array[Double]],
  n_factors : Int,
) -> Array[Array[Double]] {
  let genes = counts.length()
  let cells = counts[0].length()
  let factors = zinbwave_zero_matrix(cells, 0)
  if n_factors == 0 {
    return factors
  }
  let transformed = zinbwave_zero_matrix(cells, genes)
  for gene in 0..<genes {
    let mut center = 0.0
    for cell in 0..<cells {
      transformed[cell][gene] = @math.ln(counts[gene][cell] + 1.0) -
        offsets[gene][cell]
      center = center + transformed[cell][gene]
    }
    center = center / cells.to_double()
    let mut sum_square = 0.0
    for cell in 0..<cells {
      transformed[cell][gene] = transformed[cell][gene] - center
      sum_square = sum_square +
        transformed[cell][gene] * transformed[cell][gene]
    }
    let scale = (sum_square / cells.to_double()).sqrt().max(1.0)
    for cell in 0..<cells {
      transformed[cell][gene] = transformed[cell][gene] / scale
    }
  }
  let gram = zinbwave_zero_matrix(cells, cells)
  for left in 0..<cells {
    for right in left..<cells {
      let mut value = 0.0
      for gene in 0..<genes {
        value = value + transformed[left][gene] * transformed[right][gene]
      }
      value = value / genes.to_double()
      gram[left][right] = value
      gram[right][left] = value
    }
  }
  for factor in 0..<n_factors {
    let vector = Array::make(cells, 0.0)
    for cell in 0..<cells {
      vector[cell] = ((cell + 1) * (factor + 2) % (cells + 3) + 1).to_double()
    }
    let mut center = stat_mean(vector)
    for cell in 0..<cells {
      vector[cell] = vector[cell] - center
    }
    for _ in 0..<80 {
      let next = Array::make(cells, 0.0)
      for left in 0..<cells {
        for right in 0..<cells {
          next[left] = next[left] + gram[left][right] * vector[right]
        }
      }
      for previous in 0..<factor {
        let mut projection = 0.0
        for cell in 0..<cells {
          projection = projection + next[cell] * factors[cell][previous]
        }
        projection = projection / cells.to_double()
        for cell in 0..<cells {
          next[cell] = next[cell] - projection * factors[cell][previous]
        }
      }
      center = stat_mean(next)
      let mut norm = 0.0
      for cell in 0..<cells {
        next[cell] = next[cell] - center
        norm = norm + next[cell] * next[cell]
      }
      norm = (norm / cells.to_double()).sqrt()
      if norm < 1.0e-10 {
        for cell in 0..<cells {
          next[cell] = if cell == factor % cells { 1.0 } else { 0.0 }
        }
        center = stat_mean(next)
        norm = 0.0
        for cell in 0..<cells {
          next[cell] = next[cell] - center
          norm = norm + next[cell] * next[cell]
        }
        norm = (norm / cells.to_double()).sqrt().max(1.0e-10)
      }
      for cell in 0..<cells {
        vector[cell] = next[cell] / norm
      }
    }
    for cell in 0..<cells {
      factors[cell].push(vector[cell])
    }
  }
  factors
}

///|
fn zinbwave_cell_latent_design(
  cell_design : Array[Array[Double]],
  factors : Array[Array[Double]],
) -> Array[Array[Double]] {
  let cells = cell_design.length()
  let known = cell_design[0].length()
  let latent = if factors.length() == 0 { 0 } else { factors[0].length() }
  let design : Array[Array[Double]] = []
  for cell in 0..<cells {
    let row = Array::make(known + latent, 0.0)
    for column in 0..<known {
      row[column] = cell_design[cell][column]
    }
    for factor in 0..<latent {
      row[known + factor] = factors[cell][factor]
    }
    design.push(row)
  }
  design
}

///|
fn zinbwave_initialize_coefficients(
  counts : Array[Array[Double]],
  mean_offsets : Array[Array[Double]],
  cell_latent_design : Array[Array[Double]],
  ridge : Double,
  minimum_probability : Double,
) -> (Array[Array[Double]], Array[Array[Double]]) {
  let genes = counts.length()
  let cells = counts[0].length()
  let beta_mu : Array[Array[Double]] = []
  let beta_pi : Array[Array[Double]] = []
  for gene in 0..<genes {
    let mean_response = Array::make(cells, 0.0)
    let mean_weights = Array::make(cells, 1.0)
    let zero_response = Array::make(cells, 0.0)
    for cell in 0..<cells {
      mean_response[cell] = @math.ln(counts[gene][cell] + 0.5) -
        mean_offsets[gene][cell]
      if counts[gene][cell] == 0.0 {
        mean_weights[cell] = 0.25
        zero_response[cell] = 1.0
      }
    }
    beta_mu.push(
      zinbwave_weighted_ridge(
        cell_latent_design, mean_response, mean_weights, ridge, false,
      ),
    )
    let zero_fraction = zinbwave_clamp(
      stat_mean(zero_response),
      minimum_probability,
      1.0 - minimum_probability,
    )
    let coefficients = Array::make(cell_latent_design[0].length(), 0.0)
    coefficients[0] = zinbwave_logit(zero_fraction, minimum_probability)
    beta_pi.push(coefficients)
  }
  (beta_mu, beta_pi)
}

///|
fn zinbwave_initial_dispersions(
  counts : Array[Array[Double]],
  mean_offsets : Array[Array[Double]],
) -> Array[Double] {
  let genes = counts.length()
  let cells = counts[0].length()
  let raw = Array::make(genes, 0.1)
  for gene in 0..<genes {
    let normalized = Array::make(cells, 0.0)
    for cell in 0..<cells {
      normalized[cell] = counts[gene][cell] /
        @math.exp(mean_offsets[gene][cell])
    }
    let mean = stat_mean(normalized)
    let variance = stat_variance(normalized)
    raw[gene] = if mean > 1.0e-8 {
      ((variance - mean) / (mean * mean)).max(1.0e-4).min(100.0)
    } else {
      1.0
    }
  }
  raw
}

///|
fn zinbwave_linear_predictors(
  beta_mu : Array[Array[Double]],
  beta_pi : Array[Array[Double]],
  gamma_mu : Array[Array[Double]],
  gamma_pi : Array[Array[Double]],
  cell_latent_design : Array[Array[Double]],
  gene_design : Array[Array[Double]],
  mean_offsets : Array[Array[Double]],
  zero_offsets : Array[Array[Double]],
  config : ZinbWaveConfig,
) -> (Array[Array[Double]], Array[Array[Double]]) {
  let genes = beta_mu.length()
  let cells = cell_latent_design.length()
  let means = zinbwave_zero_matrix(genes, cells)
  let probabilities = zinbwave_zero_matrix(genes, cells)
  for gene in 0..<genes {
    for cell in 0..<cells {
      let mut mean_linear = mean_offsets[gene][cell]
      let mut zero_linear = zero_offsets[gene][cell]
      for parameter in 0..<cell_latent_design[cell].length() {
        mean_linear = mean_linear +
          cell_latent_design[cell][parameter] * beta_mu[gene][parameter]
        zero_linear = zero_linear +
          cell_latent_design[cell][parameter] * beta_pi[gene][parameter]
      }
      for parameter in 0..<gene_design[gene].length() {
        mean_linear = mean_linear +
          gene_design[gene][parameter] * gamma_mu[cell][parameter]
        zero_linear = zero_linear +
          gene_design[gene][parameter] * gamma_pi[cell][parameter]
      }
      means[gene][cell] = @math.exp(zinbwave_clamp(mean_linear, -30.0, 30.0)).max(
        config.minimum_mean,
      )
      probabilities[gene][cell] = zinbwave_sigmoid(
        zero_linear,
        config.minimum_probability,
      )
    }
  }
  (means, probabilities)
}

///|
fn zinbwave_lgamma(value : Double) -> Double {
  let coefficients = [
    676.5203681218851, -1259.1392167224028, 771.3234287776531, -176.6150291621406,
    12.507343278686905, -0.13857109526572012, 0.000009984369578019572, 0.00000015056327351493116,
  ]
  if value < 0.5 {
    return @math.ln(3.141592653589793 / @math.sin(3.141592653589793 * value)) -
      zinbwave_lgamma(1.0 - value)
  }
  let shifted = value - 1.0
  let mut series = 0.9999999999998099
  for index in 0..<coefficients.length() {
    series = series + coefficients[index] / (shifted + index.to_double() + 1.0)
  }
  let base = shifted + 7.5
  0.9189385332046727 +
  (shifted + 0.5) * @math.ln(base) -
  base +
  @math.ln(series)
}

///|
fn zinbwave_nb_log_probability(
  count : Double,
  mean : Double,
  dispersion : Double,
) -> Double {
  let theta = 1.0 / dispersion.max(1.0e-8)
  zinbwave_lgamma(count + theta) -
  zinbwave_lgamma(theta) -
  zinbwave_lgamma(count + 1.0) +
  theta * (@math.ln(theta) - @math.ln(theta + mean)) +
  count * (@math.ln(mean) - @math.ln(theta + mean))
}

///|
fn zinbwave_responsibilities(
  counts : Array[Array[Double]],
  means : Array[Array[Double]],
  probabilities : Array[Array[Double]],
  dispersions : Array[Double],
) -> Array[Array[Double]] {
  let genes = counts.length()
  let cells = counts[0].length()
  let weights = zinbwave_zero_matrix(genes, cells)
  for gene in 0..<genes {
    let theta = 1.0 / dispersions[gene].max(1.0e-8)
    for cell in 0..<cells {
      if counts[gene][cell] > 0.0 {
        weights[gene][cell] = 1.0
      } else {
        let nb_zero = @math.exp(
          theta * (@math.ln(theta) - @math.ln(theta + means[gene][cell])),
        )
        let probability = probabilities[gene][cell]
        let mixture = probability + (1.0 - probability) * nb_zero
        weights[gene][cell] = ((1.0 - probability) *
        nb_zero /
        mixture.max(1.0e-300)).max(1.0e-12)
      }
    }
  }
  weights
}

///|
fn zinbwave_log_likelihood(
  counts : Array[Array[Double]],
  means : Array[Array[Double]],
  probabilities : Array[Array[Double]],
  dispersions : Array[Double],
) -> Double {
  let mut total = 0.0
  for gene in 0..<counts.length() {
    let theta = 1.0 / dispersions[gene].max(1.0e-8)
    for cell in 0..<counts[gene].length() {
      let probability = probabilities[gene][cell]
      if counts[gene][cell] == 0.0 {
        let nb_zero = @math.exp(
          theta * (@math.ln(theta) - @math.ln(theta + means[gene][cell])),
        )
        total = total +
          @math.ln((probability + (1.0 - probability) * nb_zero).max(1.0e-300))
      } else {
        total = total +
          @math.ln((1.0 - probability).max(1.0e-300)) +
          zinbwave_nb_log_probability(
            counts[gene][cell],
            means[gene][cell],
            dispersions[gene],
          )
      }
    }
  }
  total
}

///|
fn zinbwave_update_gene_coefficients(
  counts : Array[Array[Double]],
  means : Array[Array[Double]],
  probabilities : Array[Array[Double]],
  responsibilities : Array[Array[Double]],
  dispersions : Array[Double],
  cell_latent_design : Array[Array[Double]],
  gene_design : Array[Array[Double]],
  gamma_mu : Array[Array[Double]],
  gamma_pi : Array[Array[Double]],
  mean_offsets : Array[Array[Double]],
  zero_offsets : Array[Array[Double]],
  beta_mu : Array[Array[Double]],
  beta_pi : Array[Array[Double]],
  config : ZinbWaveConfig,
) -> Double {
  let genes = counts.length()
  let cells = counts[0].length()
  let mut maximum_change = 0.0
  for gene in 0..<genes {
    let mean_response = Array::make(cells, 0.0)
    let mean_weights = Array::make(cells, 0.0)
    let zero_response = Array::make(cells, 0.0)
    let zero_weights = Array::make(cells, 0.0)
    for cell in 0..<cells {
      let mean = means[gene][cell]
      let mut current_mean_linear = @math.ln(mean)
      let probability = probabilities[gene][cell]
      let current_zero_linear = zinbwave_logit(
        probability,
        config.minimum_probability,
      )
      let mut gene_mean_effect = mean_offsets[gene][cell]
      let mut gene_zero_effect = zero_offsets[gene][cell]
      for parameter in 0..<gene_design[gene].length() {
        gene_mean_effect = gene_mean_effect +
          gene_design[gene][parameter] * gamma_mu[cell][parameter]
        gene_zero_effect = gene_zero_effect +
          gene_design[gene][parameter] * gamma_pi[cell][parameter]
      }
      current_mean_linear = current_mean_linear +
        (counts[gene][cell] - mean) / mean
      mean_response[cell] = current_mean_linear - gene_mean_effect
      mean_weights[cell] = (responsibilities[gene][cell] *
      mean /
      (1.0 + dispersions[gene] * mean)).max(1.0e-10)
      let logistic_variance = (probability * (1.0 - probability)).max(1.0e-8)
      let structural_zero = 1.0 - responsibilities[gene][cell]
      zero_response[cell] = current_zero_linear +
        (structural_zero - probability) / logistic_variance -
        gene_zero_effect
      zero_weights[cell] = logistic_variance.max(1.0e-10)
    }
    let updated_mean = zinbwave_weighted_ridge(
      cell_latent_design,
      mean_response,
      mean_weights,
      config.ridge,
      false,
    )
    let updated_zero = zinbwave_weighted_ridge(
      cell_latent_design,
      zero_response,
      zero_weights,
      config.ridge,
      false,
    )
    for parameter in 0..<updated_mean.length() {
      maximum_change = maximum_change.max(
        (updated_mean[parameter] - beta_mu[gene][parameter]).abs(),
      )
      maximum_change = maximum_change.max(
        (updated_zero[parameter] - beta_pi[gene][parameter]).abs(),
      )
      beta_mu[gene][parameter] = updated_mean[parameter]
      beta_pi[gene][parameter] = updated_zero[parameter]
    }
  }
  maximum_change
}

///|
fn zinbwave_update_cell_gene_effects(
  counts : Array[Array[Double]],
  means : Array[Array[Double]],
  probabilities : Array[Array[Double]],
  responsibilities : Array[Array[Double]],
  dispersions : Array[Double],
  cell_latent_design : Array[Array[Double]],
  gene_design : Array[Array[Double]],
  beta_mu : Array[Array[Double]],
  beta_pi : Array[Array[Double]],
  mean_offsets : Array[Array[Double]],
  zero_offsets : Array[Array[Double]],
  gamma_mu : Array[Array[Double]],
  gamma_pi : Array[Array[Double]],
  config : ZinbWaveConfig,
) -> Double {
  let genes = counts.length()
  let cells = counts[0].length()
  let mut maximum_change = 0.0
  for cell in 0..<cells {
    let mean_response = Array::make(genes, 0.0)
    let mean_weights = Array::make(genes, 0.0)
    let zero_response = Array::make(genes, 0.0)
    let zero_weights = Array::make(genes, 0.0)
    for gene in 0..<genes {
      let mean = means[gene][cell]
      let probability = probabilities[gene][cell]
      let mut cell_mean_effect = mean_offsets[gene][cell]
      let mut cell_zero_effect = zero_offsets[gene][cell]
      for parameter in 0..<cell_latent_design[cell].length() {
        cell_mean_effect = cell_mean_effect +
          cell_latent_design[cell][parameter] * beta_mu[gene][parameter]
        cell_zero_effect = cell_zero_effect +
          cell_latent_design[cell][parameter] * beta_pi[gene][parameter]
      }
      mean_response[gene] = @math.ln(mean) +
        (counts[gene][cell] - mean) / mean -
        cell_mean_effect
      mean_weights[gene] = (responsibilities[gene][cell] *
      mean /
      (1.0 + dispersions[gene] * mean)).max(1.0e-10)
      let logistic_variance = (probability * (1.0 - probability)).max(1.0e-8)
      let structural_zero = 1.0 - responsibilities[gene][cell]
      zero_response[gene] = zinbwave_logit(
          probability,
          config.minimum_probability,
        ) +
        (structural_zero - probability) / logistic_variance -
        cell_zero_effect
      zero_weights[gene] = logistic_variance.max(1.0e-10)
    }
    let updated_mean = zinbwave_weighted_ridge(
      gene_design,
      mean_response,
      mean_weights,
      config.ridge,
      false,
    )
    let updated_zero = zinbwave_weighted_ridge(
      gene_design,
      zero_response,
      zero_weights,
      config.ridge,
      false,
    )
    for parameter in 0..<updated_mean.length() {
      maximum_change = maximum_change.max(
        (updated_mean[parameter] - gamma_mu[cell][parameter]).abs(),
      )
      maximum_change = maximum_change.max(
        (updated_zero[parameter] - gamma_pi[cell][parameter]).abs(),
      )
      gamma_mu[cell][parameter] = updated_mean[parameter]
      gamma_pi[cell][parameter] = updated_zero[parameter]
    }
  }
  maximum_change
}

///|
fn zinbwave_update_factors(
  counts : Array[Array[Double]],
  means : Array[Array[Double]],
  probabilities : Array[Array[Double]],
  responsibilities : Array[Array[Double]],
  dispersions : Array[Double],
  cell_design : Array[Array[Double]],
  gene_design : Array[Array[Double]],
  beta_mu : Array[Array[Double]],
  beta_pi : Array[Array[Double]],
  gamma_mu : Array[Array[Double]],
  gamma_pi : Array[Array[Double]],
  mean_offsets : Array[Array[Double]],
  zero_offsets : Array[Array[Double]],
  factors : Array[Array[Double]],
  config : ZinbWaveConfig,
) -> Double {
  let latent = config.n_factors
  if latent == 0 {
    return 0.0
  }
  let genes = counts.length()
  let cells = counts[0].length()
  let known = cell_design[0].length()
  let mut maximum_change = 0.0
  for cell in 0..<cells {
    let information = zinbwave_zero_matrix(latent, latent)
    let score = Array::make(latent, 0.0)
    for gene in 0..<genes {
      let mean = means[gene][cell]
      let probability = probabilities[gene][cell]
      let mut fixed_mean = mean_offsets[gene][cell]
      let mut fixed_zero = zero_offsets[gene][cell]
      for parameter in 0..<known {
        fixed_mean = fixed_mean +
          cell_design[cell][parameter] * beta_mu[gene][parameter]
        fixed_zero = fixed_zero +
          cell_design[cell][parameter] * beta_pi[gene][parameter]
      }
      for parameter in 0..<gene_design[gene].length() {
        fixed_mean = fixed_mean +
          gene_design[gene][parameter] * gamma_mu[cell][parameter]
        fixed_zero = fixed_zero +
          gene_design[gene][parameter] * gamma_pi[cell][parameter]
      }
      let mean_target = @math.ln(mean) +
        (counts[gene][cell] - mean) / mean -
        fixed_mean
      let mean_weight = (responsibilities[gene][cell] *
      mean /
      (1.0 + dispersions[gene] * mean)).max(1.0e-10)
      let logistic_variance = (probability * (1.0 - probability)).max(1.0e-8)
      let structural_zero = 1.0 - responsibilities[gene][cell]
      let zero_target = zinbwave_logit(probability, config.minimum_probability) +
        (structural_zero - probability) / logistic_variance -
        fixed_zero
      for left in 0..<latent {
        let mean_loading = beta_mu[gene][known + left]
        let zero_loading = beta_pi[gene][known + left]
        score[left] = score[left] +
          mean_loading * mean_weight * mean_target +
          zero_loading * logistic_variance * zero_target
        for right in 0..<latent {
          information[left][right] = information[left][right] +
            mean_loading * mean_weight * beta_mu[gene][known + right] +
            zero_loading * logistic_variance * beta_pi[gene][known + right]
        }
      }
    }
    for factor in 0..<latent {
      information[factor][factor] = information[factor][factor] +
        config.ridge +
        1.0e-8
    }
    let updated = zinbwave_solve(information, score)
    for factor in 0..<latent {
      maximum_change = maximum_change.max(
        (updated[factor] - factors[cell][factor]).abs(),
      )
      factors[cell][factor] = updated[factor]
    }
  }
  // Center and scale W while preserving the linear predictors in beta.
  for factor in 0..<latent {
    let column = Array::make(cells, 0.0)
    for cell in 0..<cells {
      column[cell] = factors[cell][factor]
    }
    let center = stat_mean(column)
    let mut sum_square = 0.0
    for cell in 0..<cells {
      factors[cell][factor] = factors[cell][factor] - center
      sum_square = sum_square + factors[cell][factor] * factors[cell][factor]
    }
    let scale = (sum_square / cells.to_double()).sqrt().max(1.0e-8)
    for cell in 0..<cells {
      factors[cell][factor] = factors[cell][factor] / scale
    }
    for gene in 0..<genes {
      let mean_loading = beta_mu[gene][known + factor]
      let zero_loading = beta_pi[gene][known + factor]
      beta_mu[gene][0] = beta_mu[gene][0] + center * mean_loading
      beta_pi[gene][0] = beta_pi[gene][0] + center * zero_loading
      beta_mu[gene][known + factor] = mean_loading * scale
      beta_pi[gene][known + factor] = zero_loading * scale
    }
  }
  maximum_change
}

///|
fn zinbwave_update_dispersions(
  counts : Array[Array[Double]],
  means : Array[Array[Double]],
  responsibilities : Array[Array[Double]],
  previous : Array[Double],
  shrinkage : Double,
) -> Array[Double] {
  let genes = counts.length()
  let cells = counts[0].length()
  let raw = Array::make(genes, 0.1)
  for gene in 0..<genes {
    let mut numerator = 0.0
    let mut denominator = 0.0
    for cell in 0..<cells {
      let mean = means[gene][cell]
      let difference = counts[gene][cell] - mean
      numerator = numerator +
        responsibilities[gene][cell] *
        (difference * difference - counts[gene][cell])
      denominator = denominator + responsibilities[gene][cell] * mean * mean
    }
    raw[gene] = if denominator > 1.0e-12 {
      (numerator / denominator).max(1.0e-4).min(100.0)
    } else {
      previous[gene]
    }
  }
  let common = stat_median(raw).max(1.0e-4)
  let information = cells.to_double()
  let updated = Array::make(genes, 0.1)
  for gene in 0..<genes {
    updated[gene] = zinbwave_clamp(
      (information * raw[gene] + shrinkage * common) / (information + shrinkage),
      1.0e-4,
      100.0,
    )
  }
  updated
}

///|
fn zinbwave_deviance_residuals(
  counts : Array[Array[Double]],
  means : Array[Array[Double]],
  weights : Array[Array[Double]],
  dispersions : Array[Double],
) -> Array[Array[Double]] {
  let genes = counts.length()
  let cells = counts[0].length()
  let residuals = zinbwave_zero_matrix(genes, cells)
  for gene in 0..<genes {
    let theta = 1.0 / dispersions[gene].max(1.0e-8)
    for cell in 0..<cells {
      let count = counts[gene][cell]
      let mean = means[gene][cell]
      let first = if count > 0.0 { count * @math.ln(count / mean) } else { 0.0 }
      let second = (count + theta) * @math.ln((count + theta) / (mean + theta))
      let deviance = (2.0 * (first - second)).max(0.0)
      let sign = if count >= mean { 1.0 } else { -1.0 }
      residuals[gene][cell] = sign * (weights[gene][cell] * deviance).sqrt()
    }
  }
  residuals
}

///|
pub fn zinbwave_fit(
  counts : Array[Array[Double]],
  config? : ZinbWaveConfig = ZinbWaveConfig::default(),
  cell_covariates? : Array[Array[Double]] = [],
  gene_covariates? : Array[Array[Double]] = [],
  mean_offsets? : Array[Array[Double]] = [],
  zero_offsets? : Array[Array[Double]] = [],
  gene_names? : Array[String] = [],
  cell_names? : Array[String] = [],
) -> ZinbWaveModel raise ZinbWaveError {
  let (genes, cells) = zinbwave_validate_counts(counts)
  ignore(
    zinbwave_validate_covariates(
      cell_covariates, cells, "zinbwave cell covariates",
    ),
  )
  ignore(
    zinbwave_validate_covariates(
      gene_covariates, genes, "zinbwave gene covariates",
    ),
  )
  if config.n_factors >= cells {
    raise ZinbWaveError(
      "zinbwave factor count must be smaller than the number of cells",
    )
  }
  if config.n_factors > genes {
    raise ZinbWaveError(
      "zinbwave factor count cannot exceed the number of genes",
    )
  }
  let prepared_gene_names = zinbwave_prepare_names(gene_names, genes, "Gene")
  let prepared_cell_names = zinbwave_prepare_names(cell_names, cells, "Cell")
  let prepared_mean_offsets = match
    zinbwave_validate_offsets(
      mean_offsets, genes, cells, "zinbwave mean offsets",
    ) {
    [] => zinbwave_default_mean_offsets(counts)
    provided => provided
  }
  let prepared_zero_offsets = match
    zinbwave_validate_offsets(
      zero_offsets, genes, cells, "zinbwave zero offsets",
    ) {
    [] => zinbwave_zero_matrix(genes, cells)
    provided => provided
  }
  let cell_design = zinbwave_design_with_intercept(cell_covariates, cells)
  let gene_design = zinbwave_design_with_intercept(gene_covariates, genes)
  let factors = zinbwave_initialize_factors(
    counts,
    prepared_mean_offsets,
    config.n_factors,
  )
  let mut cell_latent_design = zinbwave_cell_latent_design(cell_design, factors)
  let (beta_mu, beta_pi) = zinbwave_initialize_coefficients(
    counts,
    prepared_mean_offsets,
    cell_latent_design,
    config.ridge,
    config.minimum_probability,
  )
  let gamma_mu = zinbwave_zero_matrix(cells, gene_design[0].length())
  let gamma_pi = zinbwave_zero_matrix(cells, gene_design[0].length())
  let mut dispersions = zinbwave_initial_dispersions(
    counts, prepared_mean_offsets,
  )
  let (initial_means, initial_probabilities) = zinbwave_linear_predictors(
    beta_mu, beta_pi, gamma_mu, gamma_pi, cell_latent_design, gene_design, prepared_mean_offsets,
    prepared_zero_offsets, config,
  )
  let initial_log_likelihood = zinbwave_log_likelihood(
    counts, initial_means, initial_probabilities, dispersions,
  )
  let history : Array[Double] = [initial_log_likelihood]
  let mut means = initial_means
  let mut probabilities = initial_probabilities
  let mut responsibilities = zinbwave_responsibilities(
    counts, means, probabilities, dispersions,
  )
  let mut converged = false
  let mut completed_iterations = 0
  let mut maximum_parameter_change = 0.0
  for iteration in 0..<config.max_iterations {
    maximum_parameter_change = 0.0
    for _ in 0..<config.inner_iterations {
      maximum_parameter_change = maximum_parameter_change.max(
        zinbwave_update_gene_coefficients(
          counts, means, probabilities, responsibilities, dispersions, cell_latent_design,
          gene_design, gamma_mu, gamma_pi, prepared_mean_offsets, prepared_zero_offsets,
          beta_mu, beta_pi, config,
        ),
      )
      let after_gene_update = zinbwave_linear_predictors(
        beta_mu, beta_pi, gamma_mu, gamma_pi, cell_latent_design, gene_design, prepared_mean_offsets,
        prepared_zero_offsets, config,
      )
      means = after_gene_update.0
      probabilities = after_gene_update.1
      maximum_parameter_change = maximum_parameter_change.max(
        zinbwave_update_cell_gene_effects(
          counts, means, probabilities, responsibilities, dispersions, cell_latent_design,
          gene_design, beta_mu, beta_pi, prepared_mean_offsets, prepared_zero_offsets,
          gamma_mu, gamma_pi, config,
        ),
      )
      let after_covariate_update = zinbwave_linear_predictors(
        beta_mu, beta_pi, gamma_mu, gamma_pi, cell_latent_design, gene_design, prepared_mean_offsets,
        prepared_zero_offsets, config,
      )
      means = after_covariate_update.0
      probabilities = after_covariate_update.1
      maximum_parameter_change = maximum_parameter_change.max(
        zinbwave_update_factors(
          counts, means, probabilities, responsibilities, dispersions, cell_design,
          gene_design, beta_mu, beta_pi, gamma_mu, gamma_pi, prepared_mean_offsets,
          prepared_zero_offsets, factors, config,
        ),
      )
      cell_latent_design = zinbwave_cell_latent_design(cell_design, factors)
      let after_factor_update = zinbwave_linear_predictors(
        beta_mu, beta_pi, gamma_mu, gamma_pi, cell_latent_design, gene_design, prepared_mean_offsets,
        prepared_zero_offsets, config,
      )
      means = after_factor_update.0
      probabilities = after_factor_update.1
    }
    dispersions = zinbwave_update_dispersions(
      counts,
      means,
      responsibilities,
      dispersions,
      config.dispersion_shrinkage,
    )
    responsibilities = zinbwave_responsibilities(
      counts, means, probabilities, dispersions,
    )
    let log_likelihood = zinbwave_log_likelihood(
      counts, means, probabilities, dispersions,
    )
    history.push(log_likelihood)
    completed_iterations = iteration + 1
    let previous = history[history.length() - 2]
    let relative_change = (log_likelihood - previous).abs() /
      (1.0 + previous.abs())
    if relative_change < config.tolerance &&
      maximum_parameter_change < config.tolerance.sqrt() {
      converged = true
      break
    }
  }
  let residuals = zinbwave_deviance_residuals(
    counts, means, responsibilities, dispersions,
  )
  let final_log_likelihood = history[history.length() - 1]
  ZinbWaveModel::{
    gene_names: prepared_gene_names,
    cell_names: prepared_cell_names,
    factors,
    cell_design,
    gene_design,
    beta_mu,
    beta_pi,
    gamma_mu,
    gamma_pi,
    mean_offsets: prepared_mean_offsets,
    zero_offsets: prepared_zero_offsets,
    dispersions,
    means,
    zero_probabilities: probabilities,
    observational_weights: responsibilities,
    deviance_residuals: residuals,
    log_likelihood: final_log_likelihood,
    diagnostics: ZinbWaveDiagnostics::{
      iterations: completed_iterations,
      converged,
      initial_log_likelihood,
      final_log_likelihood,
      log_likelihoods: history,
      maximum_parameter_change,
    },
    config,
  }
}

///|
pub fn ZinbWaveModel::n_genes(self : ZinbWaveModel) -> Int {
  self.gene_names.length()
}

///|
pub fn ZinbWaveModel::n_cells(self : ZinbWaveModel) -> Int {
  self.cell_names.length()
}

///|
pub fn ZinbWaveModel::n_factors(self : ZinbWaveModel) -> Int {
  self.config.n_factors
}

///|
pub fn ZinbWaveModel::n_parameters(self : ZinbWaveModel) -> Int {
  let genes = self.n_genes()
  let cells = self.n_cells()
  genes * self.beta_mu[0].length() * 2 +
  cells * self.gamma_mu[0].length() * 2 +
  cells * self.n_factors() +
  genes
}

///|
pub fn ZinbWaveModel::aic(self : ZinbWaveModel) -> Double {
  2.0 * self.n_parameters().to_double() - 2.0 * self.log_likelihood
}

///|
pub fn ZinbWaveModel::bic(self : ZinbWaveModel) -> Double {
  let observations = (self.n_genes() * self.n_cells()).to_double()
  @math.ln(observations) * self.n_parameters().to_double() -
  2.0 * self.log_likelihood
}

///|
pub fn ZinbWaveModel::mean_zero_probability(self : ZinbWaveModel) -> Double {
  let mut total = 0.0
  let mut count = 0
  for row in self.zero_probabilities {
    for value in row {
      total = total + value
      count = count + 1
    }
  }
  if count == 0 {
    0.0
  } else {
    total / count.to_double()
  }
}

///|
pub fn ZinbWaveModel::expected_mean(
  self : ZinbWaveModel,
  gene : Int,
  cell : Int,
) -> Double? {
  if gene < 0 || gene >= self.n_genes() || cell < 0 || cell >= self.n_cells() {
    None
  } else {
    Some(self.means[gene][cell])
  }
}

///|
pub fn ZinbWaveModel::zero_probability_at(
  self : ZinbWaveModel,
  gene : Int,
  cell : Int,
) -> Double? {
  if gene < 0 || gene >= self.n_genes() || cell < 0 || cell >= self.n_cells() {
    None
  } else {
    Some(self.zero_probabilities[gene][cell])
  }
}

///|
pub fn ZinbWaveModel::cell_factors(
  self : ZinbWaveModel,
  cell : Int,
) -> Array[Double]? {
  if cell < 0 || cell >= self.n_cells() {
    None
  } else {
    Some(self.factors[cell].copy())
  }
}

///|
pub fn ZinbWaveModel::normalized_values(
  self : ZinbWaveModel,
) -> Array[Array[Double]] {
  let normalized = zinbwave_zero_matrix(self.n_genes(), self.n_cells())
  for gene in 0..<self.n_genes() {
    for cell in 0..<self.n_cells() {
      normalized[gene][cell] = self.means[gene][cell] /
        @math.exp(self.mean_offsets[gene][cell])
    }
  }
  normalized
}

///|
pub fn ZinbWaveModel::impute_zeros(
  self : ZinbWaveModel,
  counts : Array[Array[Double]],
) -> Array[Array[Double]] raise ZinbWaveError {
  let (genes, cells) = zinbwave_validate_counts(counts)
  if genes != self.n_genes() || cells != self.n_cells() {
    raise ZinbWaveError("zinbwave imputation counts do not match the model")
  }
  let imputed = zinbwave_copy_matrix(counts)
  for gene in 0..<genes {
    for cell in 0..<cells {
      if counts[gene][cell] == 0.0 {
        imputed[gene][cell] = (1.0 - self.observational_weights[gene][cell]) *
          self.means[gene][cell]
      }
    }
  }
  imputed
}

///|
pub fn ZinbWaveModel::gene_summary(
  self : ZinbWaveModel,
  counts : Array[Array[Double]],
  gene : Int,
) -> ZinbWaveGeneSummary? {
  if gene < 0 ||
    gene >= self.n_genes() ||
    counts.length() != self.n_genes() ||
    counts[gene].length() != self.n_cells() {
    return None
  }
  let mut observed_zeros = 0.0
  let mut fitted_zeros = 0.0
  let mut observed_mean = 0.0
  let mut fitted_mean = 0.0
  let theta = 1.0 / self.dispersions[gene]
  for cell in 0..<self.n_cells() {
    if counts[gene][cell] == 0.0 {
      observed_zeros = observed_zeros + 1.0
    }
    let nb_zero = @math.exp(
      theta * (@math.ln(theta) - @math.ln(theta + self.means[gene][cell])),
    )
    fitted_zeros = fitted_zeros +
      self.zero_probabilities[gene][cell] +
      (1.0 - self.zero_probabilities[gene][cell]) * nb_zero
    observed_mean = observed_mean + counts[gene][cell]
    fitted_mean = fitted_mean + self.means[gene][cell]
  }
  Some(ZinbWaveGeneSummary::{
    gene_name: self.gene_names[gene],
    dispersion: self.dispersions[gene],
    inverse_dispersion: 1.0 / self.dispersions[gene],
    observed_zero_fraction: observed_zeros / self.n_cells().to_double(),
    fitted_zero_fraction: fitted_zeros / self.n_cells().to_double(),
    mean_count: observed_mean / self.n_cells().to_double(),
    mean_fitted: fitted_mean / self.n_cells().to_double(),
  })
}

///|
pub fn ZinbWaveModel::summary(self : ZinbWaveModel) -> String {
  "ZINB-WaVE(" +
  self.n_genes().to_string() +
  " genes x " +
  self.n_cells().to_string() +
  " cells, K=" +
  self.n_factors().to_string() +
  ", logLik=" +
  self.log_likelihood.to_string() +
  ", mean pi=" +
  self.mean_zero_probability().to_string() +
  ", iterations=" +
  self.diagnostics.iterations.to_string() +
  ", converged=" +
  self.diagnostics.converged.to_string() +
  ")"
}

///|
fn zinbwave_copy_sce(sce : SingleCellExperiment) -> SingleCellExperiment {
  let assays : Map[String, Array[Array[Double]]] = Map([])
  for key in sce.assays.keys() {
    assays[key] = sce.assays[key]
  }
  let row_data : Map[String, Array[String]] = Map([])
  for key in sce.row_data.keys() {
    row_data[key] = sce.row_data[key]
  }
  let col_data : Map[String, Array[String]] = Map([])
  for key in sce.col_data.keys() {
    col_data[key] = sce.col_data[key]
  }
  let reduced_dims : Map[String, Array[Array[Double]]] = Map([])
  for key in sce.reduced_dims.keys() {
    reduced_dims[key] = sce.reduced_dims[key]
  }
  let alternative_experiments : Map[String, SingleCellExperiment] = Map([])
  for key in sce.alternative_experiments.keys() {
    alternative_experiments[key] = sce.alternative_experiments[key]
  }
  let metadata : Map[String, String] = Map([])
  for key in sce.metadata.keys() {
    metadata[key] = sce.metadata[key]
  }
  SingleCellExperiment::{
    assays,
    row_data,
    col_data,
    reduced_dims,
    alternative_experiments,
    metadata,
    row_names: sce.row_names.copy(),
    col_names: sce.col_names.copy(),
  }
}

///|
pub fn zinbwave_sce(
  sce : SingleCellExperiment,
  input_assay? : String = "counts",
  reduced_dim_name? : String = "zinbwave",
  config? : ZinbWaveConfig = ZinbWaveConfig::default(),
  cell_covariates? : Array[Array[Double]] = [],
  gene_covariates? : Array[Array[Double]] = [],
  mean_offsets? : Array[Array[Double]] = [],
  zero_offsets? : Array[Array[Double]] = [],
) -> ZinbWaveSCEOutput raise ZinbWaveError {
  let counts = sce_get_assay(sce, input_assay)
  if counts.length() == 0 {
    raise ZinbWaveError(
      "zinbwave SingleCellExperiment assay not found: " + input_assay,
    )
  }
  let model = zinbwave_fit(
    counts,
    config~,
    cell_covariates~,
    gene_covariates~,
    mean_offsets~,
    zero_offsets~,
    gene_names=sce.row_names,
    cell_names=sce.col_names,
  )
  let enriched = zinbwave_copy_sce(sce)
  enriched.assays["zinbwave_weights"] = zinbwave_copy_matrix(
    model.observational_weights,
  )
  enriched.assays["zinbwave_residuals"] = zinbwave_copy_matrix(
    model.deviance_residuals,
  )
  enriched.assays["zinbwave_normalized"] = model.normalized_values()
  enriched.assays["zinbwave_imputed"] = model.impute_zeros(counts)
  enriched.reduced_dims[reduced_dim_name] = zinbwave_copy_matrix(model.factors)
  enriched.metadata["zinbwave_factors"] = model.n_factors().to_string()
  enriched.metadata["zinbwave_iterations"] = model.diagnostics.iterations.to_string()
  enriched.metadata["zinbwave_converged"] = model.diagnostics.converged.to_string()
  enriched.metadata["zinbwave_log_likelihood"] = model.log_likelihood.to_string()
  ZinbWaveSCEOutput::{ experiment: enriched, model, }
}

///|
pub fn zinbwave_example_data() -> (
  Array[Array[Double]],
  Array[String],
  Array[String],
  Array[Array[Double]],
) {
  let counts = [
    [18.0, 16.0, 20.0, 15.0, 0.0, 1.0, 0.0, 0.0],
    [12.0, 10.0, 14.0, 11.0, 0.0, 0.0, 1.0, 0.0],
    [0.0, 0.0, 1.0, 0.0, 17.0, 19.0, 15.0, 18.0],
    [1.0, 0.0, 0.0, 0.0, 11.0, 13.0, 10.0, 12.0],
    [6.0, 5.0, 7.0, 6.0, 5.0, 7.0, 6.0, 5.0],
    [0.0, 2.0, 0.0, 1.0, 0.0, 1.0, 0.0, 2.0],
  ]
  let gene_names = [
    "A_marker_1", "A_marker_2", "B_marker_1", "B_marker_2", "housekeeping", "dropout_gene",
  ]
  let cell_names = ["A1", "A2", "A3", "A4", "B1", "B2", "B3", "B4"]
  let batch = [[0.0], [0.0], [1.0], [1.0], [0.0], [0.0], [1.0], [1.0]]
  (counts, gene_names, cell_names, batch)
}
