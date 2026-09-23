// Cellular Latent Dirichlet Allocation for single-cell count matrices.
//
// This module implements the celda_CG hierarchy from Bioconductor celda.
// Matrices use feature x cell orientation and labels are zero-based.

///|
pub suberror CeldaError {
  CeldaError(String)
}

///|
pub(all) enum CeldaAlgorithm {
  CeldaEm
  CeldaGibbs
} derive(Eq, Debug)
pub extend CeldaAlgorithm with Eq::{equal, not_equal}
pub extend CeldaAlgorithm with @moonbitlang/core/debug.Debug::{to_repr}

///|
pub(all) enum CeldaCGCriterion {
  CeldaBic
  CeldaPerplexity
  CeldaLogLikelihood
} derive(Eq, Debug)
pub extend CeldaCGCriterion with Eq::{equal, not_equal}
pub extend CeldaCGCriterion with @moonbitlang/core/debug.Debug::{to_repr}

///|
pub struct CeldaCGConfig {
  cell_populations : Int
  feature_modules : Int
  alpha : Double
  beta : Double
  delta : Double
  gamma : Double
  algorithm : CeldaAlgorithm
  max_iterations : Int
  stop_iterations : Int
  chains : Int
  tolerance : Double
  seed : Int
} derive(Eq, Debug)
pub extend CeldaCGConfig with Eq::{equal, not_equal}
pub extend CeldaCGConfig with @moonbitlang/core/debug.Debug::{to_repr}

///|
pub fn CeldaCGConfig::create(
  cell_populations : Int,
  feature_modules : Int,
  alpha? : Double = 1.0,
  beta? : Double = 1.0,
  delta? : Double = 1.0,
  gamma? : Double = 1.0,
  algorithm? : CeldaAlgorithm = CeldaEm,
  max_iterations? : Int = 200,
  stop_iterations? : Int = 10,
  chains? : Int = 3,
  tolerance? : Double = 1.0e-8,
  seed? : Int = 12345,
) -> CeldaCGConfig raise CeldaError {
  if cell_populations < 1 {
    raise CeldaError::CeldaError("celda_CG cell_populations must be positive")
  }
  if feature_modules < 1 {
    raise CeldaError::CeldaError("celda_CG feature_modules must be positive")
  }
  if !celda_is_finite(alpha) || alpha <= 0.0 {
    raise CeldaError::CeldaError("celda_CG alpha must be finite and positive")
  }
  if !celda_is_finite(beta) || beta <= 0.0 {
    raise CeldaError::CeldaError("celda_CG beta must be finite and positive")
  }
  if !celda_is_finite(delta) || delta <= 0.0 {
    raise CeldaError::CeldaError("celda_CG delta must be finite and positive")
  }
  if !celda_is_finite(gamma) || gamma <= 0.0 {
    raise CeldaError::CeldaError("celda_CG gamma must be finite and positive")
  }
  if max_iterations < 1 || stop_iterations < 1 {
    raise CeldaError::CeldaError("celda_CG iteration limits must be positive")
  }
  if chains < 1 {
    raise CeldaError::CeldaError("celda_CG chains must be positive")
  }
  if !celda_is_finite(tolerance) || tolerance <= 0.0 {
    raise CeldaError::CeldaError(
      "celda_CG tolerance must be finite and positive",
    )
  }
  if seed <= 0 {
    raise CeldaError::CeldaError("celda_CG seed must be positive")
  }
  CeldaCGConfig::{
    cell_populations,
    feature_modules,
    alpha,
    beta,
    delta,
    gamma,
    algorithm,
    max_iterations,
    stop_iterations,
    chains,
    tolerance,
    seed,
  }
}

///|
pub fn CeldaCGConfig::with_dimensions(
  self : CeldaCGConfig,
  cell_populations : Int,
  feature_modules : Int,
  seed? : Int = self.seed,
) -> CeldaCGConfig raise CeldaError {
  CeldaCGConfig::create(
    cell_populations,
    feature_modules,
    alpha=self.alpha,
    beta=self.beta,
    delta=self.delta,
    gamma=self.gamma,
    algorithm=self.algorithm,
    max_iterations=self.max_iterations,
    stop_iterations=self.stop_iterations,
    chains=self.chains,
    tolerance=self.tolerance,
    seed~,
  )
}

///|
pub struct CeldaCGDiagnostics {
  iterations : Int
  converged : Bool
  log_likelihoods : Array[Double]
  chain_scores : Array[Double]
  best_chain : Int
} derive(Eq, Debug)
pub extend CeldaCGDiagnostics with Eq::{equal, not_equal}
pub extend CeldaCGDiagnostics with @moonbitlang/core/debug.Debug::{to_repr}

///|
pub struct CeldaFeatureScore {
  feature_index : Int
  feature_name : String
  probability : Double
  total_count : Double
} derive(Eq, Debug)
pub extend CeldaFeatureScore with Eq::{equal, not_equal}
pub extend CeldaFeatureScore with @moonbitlang/core/debug.Debug::{to_repr}

///|
pub struct CeldaCGResult {
  cell_clusters : Array[Int]
  feature_modules : Array[Int]
  sample_indices : Array[Int]
  sample_names : Array[String]
  feature_names : Array[String]
  cell_names : Array[String]
  theta : Array[Array[Double]]
  phi : Array[Array[Double]]
  psi : Array[Array[Double]]
  eta : Array[Double]
  fitted_counts : Array[Array[Double]]
  feature_totals : Array[Double]
  cell_totals : Array[Double]
  log_likelihood : Double
  perplexity : Double
  aic : Double
  bic : Double
  diagnostics : CeldaCGDiagnostics
  config : CeldaCGConfig
}

///|
pub struct CeldaCGPrediction {
  probabilities : Array[Array[Double]]
  assignments : Array[Int]
} derive(Eq, Debug)
pub extend CeldaCGPrediction with Eq::{equal, not_equal}
pub extend CeldaCGPrediction with @moonbitlang/core/debug.Debug::{to_repr}

///|
pub struct CeldaCGGridEntry {
  cell_populations : Int
  feature_modules : Int
  score : Double
  result : CeldaCGResult
}

///|
pub struct CeldaCGGridSearch {
  entries : Array[CeldaCGGridEntry]
  best_index : Int
  criterion : CeldaCGCriterion
}

///|
pub struct CeldaCGSCEOutput {
  experiment : SingleCellExperiment
  result : CeldaCGResult
}

///|
priv struct CeldaCGSufficient {
  population_by_sample : Array[Array[Double]]
  module_by_population : Array[Array[Double]]
  feature_totals : Array[Double]
  module_totals : Array[Double]
  features_by_module : Array[Int]
  cell_totals : Array[Double]
}

///|
priv struct CeldaRng {
  mut state : Int
}

///|
fn CeldaRng::new(seed : Int) -> CeldaRng {
  let normalized = (seed - 1) % 2147483646 + 1
  CeldaRng::{ state: normalized, }
}

///|
fn CeldaRng::uniform(self : CeldaRng) -> Double {
  let quotient = self.state / 44488
  let remainder = self.state % 44488
  let candidate = 48271 * remainder - 3399 * quotient
  self.state = if candidate > 0 { candidate } else { candidate + 2147483647 }
  self.state.to_double() / 2147483647.0
}

///|
fn celda_is_finite(value : Double) -> Bool {
  value == value && value.abs() <= 1.0e300
}

///|
fn celda_log_gamma(value : Double) -> Double {
  let coefficients = [
    0.9999999999998099, 676.5203681218851, -1259.1392167224028, 771.3234287776531,
    -176.6150291621406, 12.507343278686905, -0.13857109526572012, 9.984369578019572e-6,
    1.5056327351493116e-7,
  ]
  if value < 0.5 {
    let pi = 3.141592653589793
    return @math.ln(pi) -
      @math.ln(@math.sin(pi * value).abs()) -
      celda_log_gamma(1.0 - value)
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
fn celda_zero_matrix(rows : Int, columns : Int) -> Array[Array[Double]] {
  let output : Array[Array[Double]] = []
  for _ in 0..<rows {
    output.push(Array::make(columns, 0.0))
  }
  output
}

///|
fn celda_copy_matrix(values : Array[Array[Double]]) -> Array[Array[Double]] {
  let output : Array[Array[Double]] = []
  for row in values {
    output.push(row.copy())
  }
  output
}

///|
fn celda_default_names(prefix : String, count : Int) -> Array[String] {
  let output : Array[String] = []
  for index in 0..<count {
    output.push(prefix + (index + 1).to_string())
  }
  output
}

///|
fn celda_prepare_names(
  supplied : Array[String],
  count : Int,
  prefix : String,
) -> Array[String] raise CeldaError {
  if supplied.length() == 0 {
    return celda_default_names(prefix, count)
  }
  if supplied.length() != count {
    raise CeldaError::CeldaError(
      "celda_CG " + prefix + " names must match the matrix dimension",
    )
  }
  let seen : Map[String, Bool] = Map([])
  let output : Array[String] = []
  for name in supplied {
    let trimmed = name.trim().to_owned()
    if trimmed.length() == 0 {
      raise CeldaError::CeldaError(
        "celda_CG " + prefix + " names cannot be empty",
      )
    }
    if seen.contains(trimmed) {
      raise CeldaError::CeldaError(
        "celda_CG " + prefix + " names must be unique: " + trimmed,
      )
    }
    seen[trimmed] = true
    output.push(trimmed)
  }
  output
}

///|
fn celda_validate_counts(
  counts : Array[Array[Double]],
) -> (Int, Int) raise CeldaError {
  if counts.length() == 0 {
    raise CeldaError::CeldaError(
      "celda_CG counts must contain at least one feature",
    )
  }
  let cells = counts[0].length()
  if cells == 0 {
    raise CeldaError::CeldaError(
      "celda_CG counts must contain at least one cell",
    )
  }
  let cell_totals = Array::make(cells, 0.0)
  for feature in 0..<counts.length() {
    if counts[feature].length() != cells {
      raise CeldaError::CeldaError(
        "celda_CG counts must be a rectangular feature-by-cell matrix",
      )
    }
    let mut feature_total = 0.0
    for cell in 0..<cells {
      let value = counts[feature][cell]
      if !celda_is_finite(value) || value < 0.0 || value.floor() != value {
        raise CeldaError::CeldaError(
          "celda_CG counts must be finite non-negative integers",
        )
      }
      feature_total = feature_total + value
      cell_totals[cell] = cell_totals[cell] + value
    }
    if feature_total <= 0.0 {
      raise CeldaError::CeldaError(
        "celda_CG each feature must have a positive total count",
      )
    }
  }
  for total in cell_totals {
    if total <= 0.0 {
      raise CeldaError::CeldaError(
        "celda_CG each cell must have a positive library size",
      )
    }
  }
  (counts.length(), cells)
}

///|
fn celda_validate_labels(
  labels : Array[Int],
  count : Int,
  groups : Int,
  label : String,
) -> Unit raise CeldaError {
  if labels.length() != count {
    raise CeldaError::CeldaError(
      "celda_CG " + label + " assignments must match the data dimension",
    )
  }
  let observed = Array::make(groups, 0)
  for value in labels {
    if value < 0 || value >= groups {
      raise CeldaError::CeldaError(
        "celda_CG " + label + " assignments are out of range",
      )
    }
    observed[value] = observed[value] + 1
  }
  for size in observed {
    if size == 0 {
      raise CeldaError::CeldaError(
        "celda_CG predefined " + label + " assignments must use every label",
      )
    }
  }
}

///|
fn celda_encode_samples(
  labels : Array[String],
  cells : Int,
) -> (Array[String], Array[Int]) raise CeldaError {
  let source = if labels.length() == 0 {
    Array::make(cells, "Sample1")
  } else {
    labels
  }
  if source.length() != cells {
    raise CeldaError::CeldaError(
      "celda_CG sample labels must match the number of cells",
    )
  }
  let names : Array[String] = []
  let indices : Array[Int] = []
  for raw in source {
    let label = raw.trim().to_owned()
    if label.length() == 0 {
      raise CeldaError::CeldaError("celda_CG sample labels cannot be empty")
    }
    let mut found = -1
    for index in 0..<names.length() {
      if names[index] == label {
        found = index
      }
    }
    if found < 0 {
      names.push(label)
      found = names.length() - 1
    }
    indices.push(found)
  }
  (names, indices)
}

///|
fn celda_group_sizes(labels : Array[Int], groups : Int) -> Array[Int] {
  let sizes = Array::make(groups, 0)
  for label in labels {
    sizes[label] = sizes[label] + 1
  }
  sizes
}

///|
fn celda_balanced_random_labels(
  count : Int,
  groups : Int,
  rng : CeldaRng,
) -> Array[Int] {
  let labels : Array[Int] = []
  for index in 0..<count {
    labels.push(index % groups)
  }
  let mut index = count - 1
  while index > 0 {
    let target = (rng.uniform() * (index + 1).to_double()).floor().to_int()
    let temporary = labels[index]
    labels[index] = labels[target]
    labels[target] = temporary
    index = index - 1
  }
  labels
}

///|
fn celda_squared_distance(
  left : Array[Double],
  right : Array[Double],
) -> Double {
  let mut total = 0.0
  for index in 0..<left.length() {
    let difference = left[index] - right[index]
    total = total + difference * difference
  }
  total
}

///|
fn celda_cell_profiles(counts : Array[Array[Double]]) -> Array[Array[Double]] {
  let features = counts.length()
  let cells = counts[0].length()
  let output = celda_zero_matrix(cells, features)
  for cell in 0..<cells {
    let mut total = 0.0
    for feature in 0..<features {
      total = total + counts[feature][cell]
    }
    for feature in 0..<features {
      output[cell][feature] = counts[feature][cell] / total
    }
  }
  output
}

///|
fn celda_feature_profiles(
  counts : Array[Array[Double]],
) -> Array[Array[Double]] {
  let features = counts.length()
  let cells = counts[0].length()
  let output = celda_zero_matrix(features, cells)
  for feature in 0..<features {
    let mut total = 0.0
    for cell in 0..<cells {
      total = total + counts[feature][cell]
    }
    for cell in 0..<cells {
      output[feature][cell] = counts[feature][cell] / total
    }
  }
  output
}

///|
fn celda_farthest_labels(
  profiles : Array[Array[Double]],
  groups : Int,
  offset : Int,
) -> Array[Int] {
  let count = profiles.length()
  let centers : Array[Int] = [offset % count]
  while centers.length() < groups {
    let mut best_index = -1
    let mut best_distance = -1.0
    for index in 0..<count {
      let mut already_center = false
      for center in centers {
        if center == index {
          already_center = true
        }
      }
      if !already_center {
        let mut nearest = 1.0e300
        for center in centers {
          let distance = celda_squared_distance(
            profiles[index],
            profiles[center],
          )
          if distance < nearest {
            nearest = distance
          }
        }
        if nearest > best_distance {
          best_distance = nearest
          best_index = index
        }
      }
    }
    centers.push(best_index)
  }
  let labels = Array::make(count, 0)
  for index in 0..<count {
    let mut best_group = 0
    let mut best_distance = celda_squared_distance(
      profiles[index],
      profiles[centers[0]],
    )
    for group in 1..<groups {
      let distance = celda_squared_distance(
        profiles[index],
        profiles[centers[group]],
      )
      if distance < best_distance {
        best_distance = distance
        best_group = group
      }
    }
    labels[index] = best_group
  }
  let sizes = celda_group_sizes(labels, groups)
  for group in 0..<groups {
    if sizes[group] == 0 {
      labels[centers[group]] = group
    }
  }
  labels
}

///|
fn celda_sufficient(
  counts : Array[Array[Double]],
  sample_indices : Array[Int],
  sample_count : Int,
  cell_clusters : Array[Int],
  feature_modules : Array[Int],
  populations : Int,
  modules : Int,
) -> CeldaCGSufficient {
  let features = counts.length()
  let cells = counts[0].length()
  let population_by_sample = celda_zero_matrix(populations, sample_count)
  let module_by_population = celda_zero_matrix(modules, populations)
  let feature_totals = Array::make(features, 0.0)
  let module_totals = Array::make(modules, 0.0)
  let features_by_module = Array::make(modules, 1)
  let cell_totals = Array::make(cells, 0.0)
  for cell in 0..<cells {
    population_by_sample[cell_clusters[cell]][sample_indices[cell]] = population_by_sample[cell_clusters[cell]][sample_indices[cell]] +
      1.0
  }
  for feature in 0..<features {
    let module_index = feature_modules[feature]
    features_by_module[module_index] = features_by_module[module_index] + 1
    for cell in 0..<cells {
      let value = counts[feature][cell]
      feature_totals[feature] = feature_totals[feature] + value
      module_totals[module_index] = module_totals[module_index] + value
      cell_totals[cell] = cell_totals[cell] + value
      module_by_population[module_index][cell_clusters[cell]] = module_by_population[module_index][cell_clusters[cell]] +
        value
    }
  }
  CeldaCGSufficient::{
    population_by_sample,
    module_by_population,
    feature_totals,
    module_totals,
    features_by_module,
    cell_totals,
  }
}

///|
fn celda_log_likelihood_from_sufficient(
  sufficient : CeldaCGSufficient,
  sample_count : Int,
  config : CeldaCGConfig,
) -> Double {
  let populations = config.cell_populations
  let modules = config.feature_modules
  let mut theta = sample_count.to_double() *
    celda_log_gamma(populations.to_double() * config.alpha) -
    (sample_count * populations).to_double() * celda_log_gamma(config.alpha)
  for sample in 0..<sample_count {
    let mut denominator = 0.0
    for population in 0..<populations {
      let value = sufficient.population_by_sample[population][sample] +
        config.alpha
      theta = theta + celda_log_gamma(value)
      denominator = denominator + value
    }
    theta = theta - celda_log_gamma(denominator)
  }
  let mut phi = populations.to_double() *
    celda_log_gamma(modules.to_double() * config.beta) -
    (populations * modules).to_double() * celda_log_gamma(config.beta)
  for population in 0..<populations {
    let mut denominator = 0.0
    for module_ in 0..<modules {
      let value = sufficient.module_by_population[module_][population] +
        config.beta
      phi = phi + celda_log_gamma(value)
      denominator = denominator + value
    }
    phi = phi - celda_log_gamma(denominator)
  }
  let mut total_features_with_pseudogenes = 0
  let mut psi = 0.0
  for module_ in 0..<modules {
    total_features_with_pseudogenes = total_features_with_pseudogenes +
      sufficient.features_by_module[module_]
    psi = psi +
      celda_log_gamma(
        sufficient.features_by_module[module_].to_double() * config.delta,
      ) -
      celda_log_gamma(
        sufficient.module_totals[module_] +
        sufficient.features_by_module[module_].to_double() * config.delta,
      )
  }
  for total in sufficient.feature_totals {
    psi = psi + celda_log_gamma(total + config.delta)
  }
  psi = psi -
    total_features_with_pseudogenes.to_double() * celda_log_gamma(config.delta)
  let mut eta = celda_log_gamma(modules.to_double() * config.gamma) -
    modules.to_double() * celda_log_gamma(config.gamma)
  let mut eta_denominator = 0.0
  for size in sufficient.features_by_module {
    eta = eta + celda_log_gamma(size.to_double() + config.gamma)
    eta_denominator = eta_denominator + size.to_double() + config.gamma
  }
  eta = eta - celda_log_gamma(eta_denominator)
  theta + phi + psi + eta
}

///|
fn celda_internal_log_likelihood(
  counts : Array[Array[Double]],
  sample_indices : Array[Int],
  sample_count : Int,
  cell_clusters : Array[Int],
  feature_modules : Array[Int],
  config : CeldaCGConfig,
) -> Double {
  let sufficient = celda_sufficient(
    counts,
    sample_indices,
    sample_count,
    cell_clusters,
    feature_modules,
    config.cell_populations,
    config.feature_modules,
  )
  celda_log_likelihood_from_sufficient(sufficient, sample_count, config)
}

///|
pub fn celda_cg_log_likelihood(
  counts : Array[Array[Double]],
  sample_labels : Array[String],
  cell_clusters : Array[Int],
  feature_modules : Array[Int],
  config : CeldaCGConfig,
) -> Double raise CeldaError {
  let (features, cells) = celda_validate_counts(counts)
  if config.cell_populations > cells {
    raise CeldaError::CeldaError(
      "celda_CG cell_populations cannot exceed the number of cells",
    )
  }
  if config.feature_modules > features {
    raise CeldaError::CeldaError(
      "celda_CG feature_modules cannot exceed the number of features",
    )
  }
  celda_validate_labels(cell_clusters, cells, config.cell_populations, "cell")
  celda_validate_labels(
    feature_modules,
    features,
    config.feature_modules,
    "feature",
  )
  let (sample_names, sample_indices) = celda_encode_samples(
    sample_labels, cells,
  )
  celda_internal_log_likelihood(
    counts,
    sample_indices,
    sample_names.length(),
    cell_clusters,
    feature_modules,
    config,
  )
}

///|
fn celda_softmax(scores : Array[Double]) -> Array[Double] {
  let mut maximum = scores[0]
  for score in scores {
    if score > maximum {
      maximum = score
    }
  }
  let probabilities = Array::make(scores.length(), 0.0)
  let mut total = 0.0
  for index in 0..<scores.length() {
    let value = @math.exp((scores[index] - maximum).max(-700.0))
    probabilities[index] = value
    total = total + value
  }
  for index in 0..<probabilities.length() {
    probabilities[index] = probabilities[index] / total
  }
  probabilities
}

///|
fn celda_choose_label(
  scores : Array[Double],
  current : Int,
  algorithm : CeldaAlgorithm,
  rng : CeldaRng,
) -> Int {
  match algorithm {
    CeldaEm => {
      let mut best = current
      let mut best_score = scores[current]
      for index in 0..<scores.length() {
        if scores[index] > best_score + 1.0e-12 {
          best = index
          best_score = scores[index]
        }
      }
      best
    }
    CeldaGibbs => {
      let probabilities = celda_softmax(scores)
      let target = rng.uniform()
      let mut cumulative = 0.0
      for index in 0..<probabilities.length() {
        cumulative = cumulative + probabilities[index]
        if target <= cumulative {
          return index
        }
      }
      probabilities.length() - 1
    }
  }
}

///|
fn celda_update_features(
  counts : Array[Array[Double]],
  sample_indices : Array[Int],
  sample_count : Int,
  cell_clusters : Array[Int],
  feature_modules : Array[Int],
  config : CeldaCGConfig,
  rng : CeldaRng,
) -> Unit {
  let sizes = celda_group_sizes(feature_modules, config.feature_modules)
  for feature in 0..<feature_modules.length() {
    let current = feature_modules[feature]
    if sizes[current] > 1 {
      let scores = Array::make(config.feature_modules, -1.0e300)
      for candidate in 0..<config.feature_modules {
        feature_modules[feature] = candidate
        scores[candidate] = celda_internal_log_likelihood(
          counts, sample_indices, sample_count, cell_clusters, feature_modules, config,
        )
      }
      feature_modules[feature] = current
      let selected = celda_choose_label(scores, current, config.algorithm, rng)
      if selected != current {
        sizes[current] = sizes[current] - 1
        sizes[selected] = sizes[selected] + 1
        feature_modules[feature] = selected
      }
    }
  }
}

///|
fn celda_update_cells(
  counts : Array[Array[Double]],
  sample_indices : Array[Int],
  sample_count : Int,
  cell_clusters : Array[Int],
  feature_modules : Array[Int],
  config : CeldaCGConfig,
  rng : CeldaRng,
) -> Unit {
  let sizes = celda_group_sizes(cell_clusters, config.cell_populations)
  for cell in 0..<cell_clusters.length() {
    let current = cell_clusters[cell]
    if sizes[current] > 1 {
      let scores = Array::make(config.cell_populations, -1.0e300)
      for candidate in 0..<config.cell_populations {
        cell_clusters[cell] = candidate
        scores[candidate] = celda_internal_log_likelihood(
          counts, sample_indices, sample_count, cell_clusters, feature_modules, config,
        )
      }
      cell_clusters[cell] = current
      let selected = celda_choose_label(scores, current, config.algorithm, rng)
      if selected != current {
        sizes[current] = sizes[current] - 1
        sizes[selected] = sizes[selected] + 1
        cell_clusters[cell] = selected
      }
    }
  }
}

///|
fn celda_reorder_labels(
  labels : Array[Int],
  weights : Array[Double],
  groups : Int,
) -> Array[Int] {
  let totals = Array::make(groups, 0.0)
  for index in 0..<labels.length() {
    totals[labels[index]] = totals[labels[index]] + weights[index]
  }
  let order : Array[Int] = []
  for group in 0..<groups {
    order.push(group)
  }
  for left in 0..<groups {
    let mut best = left
    for right in (left + 1)..<groups {
      if totals[order[right]] > totals[order[best]] {
        best = right
      }
    }
    if best != left {
      let temporary = order[left]
      order[left] = order[best]
      order[best] = temporary
    }
  }
  let mapping = Array::make(groups, 0)
  for index in 0..<groups {
    mapping[order[index]] = index
  }
  let output : Array[Int] = []
  for label in labels {
    output.push(mapping[label])
  }
  output
}

///|
fn celda_factorize(
  counts : Array[Array[Double]],
  sample_indices : Array[Int],
  sample_count : Int,
  cell_clusters : Array[Int],
  feature_modules : Array[Int],
  config : CeldaCGConfig,
) -> (
  Array[Array[Double]],
  Array[Array[Double]],
  Array[Array[Double]],
  Array[Double],
  Array[Array[Double]],
  Array[Double],
  Array[Double],
) {
  let sufficient = celda_sufficient(
    counts,
    sample_indices,
    sample_count,
    cell_clusters,
    feature_modules,
    config.cell_populations,
    config.feature_modules,
  )
  let theta = celda_zero_matrix(config.cell_populations, sample_count)
  for sample in 0..<sample_count {
    let mut denominator = 0.0
    for population in 0..<config.cell_populations {
      denominator = denominator +
        sufficient.population_by_sample[population][sample] +
        config.alpha
    }
    for population in 0..<config.cell_populations {
      theta[population][sample] = (
          sufficient.population_by_sample[population][sample] + config.alpha
        ) /
        denominator
    }
  }
  let phi = celda_zero_matrix(config.feature_modules, config.cell_populations)
  for population in 0..<config.cell_populations {
    let mut denominator = 0.0
    for module_ in 0..<config.feature_modules {
      denominator = denominator +
        sufficient.module_by_population[module_][population] +
        config.beta
    }
    for module_ in 0..<config.feature_modules {
      phi[module_][population] = (
          sufficient.module_by_population[module_][population] + config.beta
        ) /
        denominator
    }
  }
  let psi = celda_zero_matrix(counts.length(), config.feature_modules)
  for feature in 0..<counts.length() {
    let module_ = feature_modules[feature]
    let denominator = sufficient.module_totals[module_] +
      sufficient.features_by_module[module_].to_double() * config.delta
    psi[feature][module_] = (sufficient.feature_totals[feature] + config.delta) /
      denominator
  }
  let eta = Array::make(config.feature_modules, 0.0)
  let mut eta_denominator = 0.0
  for size in sufficient.features_by_module {
    eta_denominator = eta_denominator + size.to_double() + config.gamma
  }
  for module_ in 0..<config.feature_modules {
    eta[module_] = (
        sufficient.features_by_module[module_].to_double() + config.gamma
      ) /
      eta_denominator
  }
  let fitted = celda_zero_matrix(counts.length(), counts[0].length())
  for cell in 0..<counts[0].length() {
    let population = cell_clusters[cell]
    let raw = Array::make(counts.length(), 0.0)
    let mut denominator = 0.0
    for feature in 0..<counts.length() {
      let module_ = feature_modules[feature]
      raw[feature] = psi[feature][module_] * phi[module_][population]
      denominator = denominator + raw[feature]
    }
    for feature in 0..<counts.length() {
      fitted[feature][cell] = sufficient.cell_totals[cell] *
        raw[feature] /
        denominator
    }
  }
  (
    theta,
    phi,
    psi,
    eta,
    fitted,
    sufficient.feature_totals,
    sufficient.cell_totals,
  )
}

///|
fn celda_perplexity(
  counts : Array[Array[Double]],
  fitted : Array[Array[Double]],
  cell_totals : Array[Double],
) -> Double {
  let mut negative_log_likelihood = 0.0
  let mut total_counts = 0.0
  for feature in 0..<counts.length() {
    for cell in 0..<counts[feature].length() {
      let value = counts[feature][cell]
      if value > 0.0 {
        let probability = (fitted[feature][cell] / cell_totals[cell]).max(
          1.0e-300,
        )
        negative_log_likelihood = negative_log_likelihood -
          value * @math.ln(probability)
        total_counts = total_counts + value
      }
    }
  }
  @math.exp(negative_log_likelihood / total_counts)
}

///|
fn celda_information_criteria(
  log_likelihood : Double,
  features : Int,
  cells : Int,
  samples : Int,
  config : CeldaCGConfig,
) -> (Double, Double) {
  let parameter_count = samples * (config.cell_populations - 1) +
    config.cell_populations * (config.feature_modules - 1) +
    features +
    config.feature_modules -
    1
  let aic = -2.0 * log_likelihood + 2.0 * parameter_count.to_double()
  let observations = (features * cells).to_double().max(2.0)
  let bic = -2.0 * log_likelihood +
    parameter_count.to_double() * @math.ln(observations)
  (aic, bic)
}

///|
pub fn celda_cg_fit(
  counts : Array[Array[Double]],
  sample_labels : Array[String],
  config : CeldaCGConfig,
  feature_names? : Array[String] = [],
  cell_names? : Array[String] = [],
  initial_cell_clusters? : Array[Int] = [],
  initial_feature_modules? : Array[Int] = [],
) -> CeldaCGResult raise CeldaError {
  let (features, cells) = celda_validate_counts(counts)
  if config.cell_populations > cells {
    raise CeldaError::CeldaError(
      "celda_CG cell_populations cannot exceed the number of cells",
    )
  }
  if config.feature_modules > features {
    raise CeldaError::CeldaError(
      "celda_CG feature_modules cannot exceed the number of features",
    )
  }
  let prepared_features = celda_prepare_names(
    feature_names, features, "Feature",
  )
  let prepared_cells = celda_prepare_names(cell_names, cells, "Cell")
  let (sample_names, sample_indices) = celda_encode_samples(
    sample_labels, cells,
  )
  if initial_cell_clusters.length() > 0 {
    celda_validate_labels(
      initial_cell_clusters,
      cells,
      config.cell_populations,
      "cell",
    )
  }
  if initial_feature_modules.length() > 0 {
    celda_validate_labels(
      initial_feature_modules,
      features,
      config.feature_modules,
      "feature",
    )
  }
  let cell_profiles = celda_cell_profiles(counts)
  let feature_profiles = celda_feature_profiles(counts)
  let chain_scores : Array[Double] = []
  let mut global_score = -1.0e300
  let mut global_cells : Array[Int] = []
  let mut global_features : Array[Int] = []
  let mut global_history : Array[Double] = []
  let mut global_iterations = 0
  let mut global_converged = false
  let mut best_chain = 0
  for chain in 0..<config.chains {
    let chain_seed = (config.seed - 1 + chain * 104729) % 2147483646 + 1
    let rng = CeldaRng::new(chain_seed)
    let cell_clusters = if initial_cell_clusters.length() > 0 {
      initial_cell_clusters.copy()
    } else if chain == 0 {
      celda_farthest_labels(cell_profiles, config.cell_populations, chain_seed)
    } else {
      celda_balanced_random_labels(cells, config.cell_populations, rng)
    }
    let feature_modules = if initial_feature_modules.length() > 0 {
      initial_feature_modules.copy()
    } else if chain == 0 {
      celda_farthest_labels(
        feature_profiles,
        config.feature_modules,
        chain_seed + 17,
      )
    } else {
      celda_balanced_random_labels(features, config.feature_modules, rng)
    }
    let initial_score = celda_internal_log_likelihood(
      counts,
      sample_indices,
      sample_names.length(),
      cell_clusters,
      feature_modules,
      config,
    )
    let history : Array[Double] = [initial_score]
    let mut chain_best = initial_score
    let mut best_cells = cell_clusters.copy()
    let mut best_features = feature_modules.copy()
    let mut without_improvement = 0
    let mut iterations = 0
    while iterations < config.max_iterations &&
          without_improvement < config.stop_iterations {
      celda_update_features(
        counts,
        sample_indices,
        sample_names.length(),
        cell_clusters,
        feature_modules,
        config,
        rng,
      )
      celda_update_cells(
        counts,
        sample_indices,
        sample_names.length(),
        cell_clusters,
        feature_modules,
        config,
        rng,
      )
      let score = celda_internal_log_likelihood(
        counts,
        sample_indices,
        sample_names.length(),
        cell_clusters,
        feature_modules,
        config,
      )
      history.push(score)
      iterations = iterations + 1
      if score > chain_best + config.tolerance {
        chain_best = score
        best_cells = cell_clusters.copy()
        best_features = feature_modules.copy()
        without_improvement = 0
      } else {
        without_improvement = without_improvement + 1
      }
    }
    chain_scores.push(chain_best)
    if chain == 0 || chain_best > global_score {
      global_score = chain_best
      global_cells = best_cells
      global_features = best_features
      global_history = history
      global_iterations = iterations
      global_converged = without_improvement >= config.stop_iterations
      best_chain = chain
    }
  }
  let cell_weights = Array::make(cells, 1.0)
  let feature_weights = Array::make(features, 0.0)
  for feature in 0..<features {
    for cell in 0..<cells {
      feature_weights[feature] = feature_weights[feature] +
        counts[feature][cell]
    }
  }
  let reordered_cells = celda_reorder_labels(
    global_cells,
    cell_weights,
    config.cell_populations,
  )
  let reordered_features = celda_reorder_labels(
    global_features,
    feature_weights,
    config.feature_modules,
  )
  global_score = celda_internal_log_likelihood(
    counts,
    sample_indices,
    sample_names.length(),
    reordered_cells,
    reordered_features,
    config,
  )
  let (theta, phi, psi, eta, fitted, feature_totals, cell_totals) = celda_factorize(
    counts,
    sample_indices,
    sample_names.length(),
    reordered_cells,
    reordered_features,
    config,
  )
  let perplexity = celda_perplexity(counts, fitted, cell_totals)
  let (aic, bic) = celda_information_criteria(
    global_score,
    features,
    cells,
    sample_names.length(),
    config,
  )
  CeldaCGResult::{
    cell_clusters: reordered_cells,
    feature_modules: reordered_features,
    sample_indices,
    sample_names,
    feature_names: prepared_features,
    cell_names: prepared_cells,
    theta,
    phi,
    psi,
    eta,
    fitted_counts: fitted,
    feature_totals,
    cell_totals,
    log_likelihood: global_score,
    perplexity,
    aic,
    bic,
    diagnostics: CeldaCGDiagnostics::{
      iterations: global_iterations,
      converged: global_converged,
      log_likelihoods: global_history,
      chain_scores,
      best_chain,
    },
    config,
  }
}

///|
pub fn CeldaCGResult::n_features(self : CeldaCGResult) -> Int {
  self.feature_names.length()
}

///|
pub fn CeldaCGResult::n_cells(self : CeldaCGResult) -> Int {
  self.cell_names.length()
}

///|
pub fn CeldaCGResult::n_samples(self : CeldaCGResult) -> Int {
  self.sample_names.length()
}

///|
pub fn CeldaCGResult::population_sizes(self : CeldaCGResult) -> Array[Int] {
  celda_group_sizes(self.cell_clusters, self.config.cell_populations)
}

///|
pub fn CeldaCGResult::module_sizes(self : CeldaCGResult) -> Array[Int] {
  celda_group_sizes(self.feature_modules, self.config.feature_modules)
}

///|
pub fn CeldaCGResult::cell_population(
  self : CeldaCGResult,
  name : String,
) -> Int? {
  for index in 0..<self.cell_names.length() {
    if self.cell_names[index] == name {
      return Some(self.cell_clusters[index])
    }
  }
  None
}

///|
pub fn CeldaCGResult::feature_module(
  self : CeldaCGResult,
  name : String,
) -> Int? {
  for index in 0..<self.feature_names.length() {
    if self.feature_names[index] == name {
      return Some(self.feature_modules[index])
    }
  }
  None
}

///|
pub fn CeldaCGResult::cells_in_population(
  self : CeldaCGResult,
  population : Int,
) -> Array[String] raise CeldaError {
  if population < 0 || population >= self.config.cell_populations {
    raise CeldaError::CeldaError("celda_CG population index is out of range")
  }
  let output : Array[String] = []
  for index in 0..<self.cell_clusters.length() {
    if self.cell_clusters[index] == population {
      output.push(self.cell_names[index])
    }
  }
  output
}

///|
pub fn CeldaCGResult::features_in_module(
  self : CeldaCGResult,
  module_ : Int,
) -> Array[String] raise CeldaError {
  if module_ < 0 || module_ >= self.config.feature_modules {
    raise CeldaError::CeldaError("celda_CG module index is out of range")
  }
  let output : Array[String] = []
  for index in 0..<self.feature_modules.length() {
    if self.feature_modules[index] == module_ {
      output.push(self.feature_names[index])
    }
  }
  output
}

///|
pub fn CeldaCGResult::top_features(
  self : CeldaCGResult,
  module_ : Int,
  limit? : Int = 10,
) -> Array[CeldaFeatureScore] raise CeldaError {
  if module_ < 0 || module_ >= self.config.feature_modules {
    raise CeldaError::CeldaError("celda_CG module index is out of range")
  }
  if limit < 1 {
    raise CeldaError::CeldaError("celda_CG top feature limit must be positive")
  }
  let scores : Array[CeldaFeatureScore] = []
  for feature in 0..<self.feature_modules.length() {
    if self.feature_modules[feature] == module_ {
      scores.push(CeldaFeatureScore::{
        feature_index: feature,
        feature_name: self.feature_names[feature],
        probability: self.psi[feature][module_],
        total_count: self.feature_totals[feature],
      })
    }
  }
  for left in 0..<scores.length() {
    let mut best = left
    for right in (left + 1)..<scores.length() {
      if scores[right].probability > scores[best].probability {
        best = right
      }
    }
    if best != left {
      let temporary = scores[left]
      scores[left] = scores[best]
      scores[best] = temporary
    }
  }
  let output : Array[CeldaFeatureScore] = []
  let count = if limit < scores.length() { limit } else { scores.length() }
  for index in 0..<count {
    output.push(scores[index])
  }
  output
}

///|
pub fn CeldaCGResult::summary(self : CeldaCGResult) -> String {
  "celda_CG: " +
  self.n_features().to_string() +
  " features x " +
  self.n_cells().to_string() +
  " cells, K=" +
  self.config.cell_populations.to_string() +
  ", L=" +
  self.config.feature_modules.to_string() +
  ", logLik=" +
  self.log_likelihood.to_string() +
  ", perplexity=" +
  self.perplexity.to_string() +
  ", converged=" +
  self.diagnostics.converged.to_string()
}

///|
pub fn celda_cg_predict_cells(
  result : CeldaCGResult,
  counts : Array[Array[Double]],
) -> CeldaCGPrediction raise CeldaError {
  let (features, cells) = celda_validate_counts(counts)
  if features != result.n_features() {
    raise CeldaError::CeldaError(
      "celda_CG prediction features must match the fitted model",
    )
  }
  let probabilities = celda_zero_matrix(cells, result.config.cell_populations)
  let assignments = Array::make(cells, 0)
  for cell in 0..<cells {
    let scores = Array::make(result.config.cell_populations, 0.0)
    for population in 0..<result.config.cell_populations {
      let raw = Array::make(features, 0.0)
      let mut denominator = 0.0
      for feature in 0..<features {
        let module_ = result.feature_modules[feature]
        raw[feature] = result.psi[feature][module_] *
          result.phi[module_][population]
        denominator = denominator + raw[feature]
      }
      let mut score = 0.0
      for feature in 0..<features {
        if counts[feature][cell] > 0.0 {
          score = score +
            counts[feature][cell] *
            @math.ln((raw[feature] / denominator).max(1.0e-300))
        }
      }
      let mut prior = 0.0
      for sample in 0..<result.n_samples() {
        prior = prior + result.theta[population][sample]
      }
      prior = prior / result.n_samples().to_double()
      scores[population] = score + @math.ln(prior.max(1.0e-300))
    }
    let posterior = celda_softmax(scores)
    probabilities[cell] = posterior
    let mut best = 0
    for population in 1..<result.config.cell_populations {
      if posterior[population] > posterior[best] {
        best = population
      }
    }
    assignments[cell] = best
  }
  CeldaCGPrediction::{ probabilities, assignments, }
}

///|
pub fn celda_cg_grid_search(
  counts : Array[Array[Double]],
  sample_labels : Array[String],
  cell_populations : Array[Int],
  feature_modules : Array[Int],
  base_config : CeldaCGConfig,
  criterion? : CeldaCGCriterion = CeldaBic,
  feature_names? : Array[String] = [],
  cell_names? : Array[String] = [],
) -> CeldaCGGridSearch raise CeldaError {
  if cell_populations.length() == 0 || feature_modules.length() == 0 {
    raise CeldaError::CeldaError(
      "celda_CG grid search dimensions cannot be empty",
    )
  }
  let entries : Array[CeldaCGGridEntry] = []
  let seen : Map[String, Bool] = Map([])
  let mut best_index = 0
  for k in cell_populations {
    for l in feature_modules {
      let key = k.to_string() + ":" + l.to_string()
      if seen.contains(key) {
        raise CeldaError::CeldaError(
          "celda_CG grid search dimensions must be unique",
        )
      }
      seen[key] = true
      let seed = (base_config.seed - 1 + entries.length() * 7919) % 2147483646 +
        1
      let config = base_config.with_dimensions(k, l, seed~)
      let result = celda_cg_fit(
        counts,
        sample_labels,
        config,
        feature_names~,
        cell_names~,
      )
      let score = match criterion {
        CeldaBic => result.bic
        CeldaPerplexity => result.perplexity
        CeldaLogLikelihood => result.log_likelihood
      }
      entries.push(CeldaCGGridEntry::{
        cell_populations: k,
        feature_modules: l,
        score,
        result,
      })
      let current = entries.length() - 1
      let better = match criterion {
        CeldaLogLikelihood => entries[current].score > entries[best_index].score
        _ => entries[current].score < entries[best_index].score
      }
      if current == 0 || better {
        best_index = current
      }
    }
  }
  CeldaCGGridSearch::{ entries, best_index, criterion, }
}

///|
pub fn CeldaCGGridSearch::best_result(
  self : CeldaCGGridSearch,
) -> CeldaCGResult {
  self.entries[self.best_index].result
}

///|
fn celda_copy_sce(sce : SingleCellExperiment) -> SingleCellExperiment {
  let assays : Map[String, Array[Array[Double]]] = Map([])
  for key in sce.assays.keys() {
    assays[key] = celda_copy_matrix(sce.assays[key])
  }
  let row_data : Map[String, Array[String]] = Map([])
  for key in sce.row_data.keys() {
    row_data[key] = sce.row_data[key].copy()
  }
  let col_data : Map[String, Array[String]] = Map([])
  for key in sce.col_data.keys() {
    col_data[key] = sce.col_data[key].copy()
  }
  let reduced_dims : Map[String, Array[Array[Double]]] = Map([])
  for key in sce.reduced_dims.keys() {
    reduced_dims[key] = celda_copy_matrix(sce.reduced_dims[key])
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
pub fn celda_cg_sce(
  sce : SingleCellExperiment,
  config : CeldaCGConfig,
  sample_column? : String = "",
  input_assay? : String = "counts",
  fitted_assay? : String = "celda_fitted",
) -> CeldaCGSCEOutput raise CeldaError {
  let counts = sce_get_assay(sce, input_assay)
  if counts.length() == 0 {
    raise CeldaError::CeldaError(
      "celda_CG SingleCellExperiment assay not found: " + input_assay,
    )
  }
  let samples = if sample_column.length() == 0 {
    []
  } else {
    let values = sce_get_col_data(sce, sample_column)
    if values.length() == 0 {
      raise CeldaError::CeldaError(
        "celda_CG SingleCellExperiment column data not found: " + sample_column,
      )
    }
    values
  }
  let result = celda_cg_fit(
    counts,
    samples,
    config,
    feature_names=sce.row_names,
    cell_names=sce.col_names,
  )
  let enriched = celda_copy_sce(sce)
  enriched.assays[fitted_assay] = celda_copy_matrix(result.fitted_counts)
  let module_labels : Array[String] = []
  for value in result.feature_modules {
    module_labels.push((value + 1).to_string())
  }
  let population_labels : Array[String] = []
  let output_samples : Array[String] = []
  for cell in 0..<result.cell_clusters.length() {
    population_labels.push((result.cell_clusters[cell] + 1).to_string())
    output_samples.push(result.sample_names[result.sample_indices[cell]])
  }
  enriched.row_data["celda_feature_module"] = module_labels
  enriched.col_data["celda_cell_population"] = population_labels
  enriched.col_data["celda_sample_label"] = output_samples
  enriched.metadata["celda_model"] = "celda_CG"
  enriched.metadata["celda_K"] = config.cell_populations.to_string()
  enriched.metadata["celda_L"] = config.feature_modules.to_string()
  enriched.metadata["celda_log_likelihood"] = result.log_likelihood.to_string()
  enriched.metadata["celda_perplexity"] = result.perplexity.to_string()
  enriched.metadata["celda_converged"] = result.diagnostics.converged.to_string()
  enriched.metadata["celda_iterations"] = result.diagnostics.iterations.to_string()
  CeldaCGSCEOutput::{ experiment: enriched, result, }
}

///|
pub fn celda_cg_example_data() -> (
  Array[Array[Double]],
  Array[String],
  Array[String],
  Array[String],
) {
  let counts = [
    [90.0, 84.0, 88.0, 82.0, 4.0, 3.0, 5.0, 4.0, 3.0, 2.0, 4.0, 3.0],
    [72.0, 69.0, 75.0, 70.0, 3.0, 2.0, 4.0, 3.0, 2.0, 2.0, 3.0, 2.0],
    [55.0, 50.0, 58.0, 53.0, 2.0, 3.0, 2.0, 3.0, 2.0, 3.0, 2.0, 3.0],
    [3.0, 4.0, 3.0, 4.0, 92.0, 87.0, 95.0, 90.0, 4.0, 3.0, 5.0, 4.0],
    [2.0, 3.0, 2.0, 3.0, 74.0, 71.0, 78.0, 73.0, 3.0, 2.0, 4.0, 3.0],
    [3.0, 2.0, 3.0, 2.0, 58.0, 54.0, 61.0, 56.0, 2.0, 3.0, 2.0, 3.0],
    [4.0, 3.0, 5.0, 4.0, 3.0, 4.0, 3.0, 4.0, 94.0, 89.0, 97.0, 92.0],
    [3.0, 2.0, 4.0, 3.0, 2.0, 3.0, 2.0, 3.0, 76.0, 72.0, 80.0, 75.0],
    [2.0, 3.0, 2.0, 3.0, 3.0, 2.0, 3.0, 2.0, 59.0, 55.0, 63.0, 57.0],
    [22.0, 20.0, 21.0, 23.0, 21.0, 22.0, 20.0, 23.0, 22.0, 21.0, 23.0, 20.0],
    [17.0, 18.0, 16.0, 19.0, 18.0, 17.0, 19.0, 16.0, 17.0, 19.0, 16.0, 18.0],
    [12.0, 11.0, 13.0, 12.0, 11.0, 13.0, 12.0, 11.0, 13.0, 12.0, 11.0, 13.0],
  ]
  let features = [
    "T_marker_1", "T_marker_2", "T_marker_3", "B_marker_1", "B_marker_2", "B_marker_3",
    "Mono_marker_1", "Mono_marker_2", "Mono_marker_3", "RPL", "RPS", "ACTB",
  ]
  let cells = [
    "T1", "T2", "T3", "T4", "B1", "B2", "B3", "B4", "M1", "M2", "M3", "M4",
  ]
  let samples = [
    "Donor1", "Donor1", "Donor2", "Donor2", "Donor1", "Donor1", "Donor2", "Donor2",
    "Donor1", "Donor1", "Donor2", "Donor2",
  ]
  (counts, features, cells, samples)
}
