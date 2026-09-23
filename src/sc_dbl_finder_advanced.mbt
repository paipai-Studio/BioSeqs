// Advanced scDblFinder workflow inspired by Bioconductor scDblFinder 1.27.6.

///|
pub suberror ScDblFinderError {
  ScDblFinderError(String)
}

///|
pub struct ScDblFinderConfig {
  expected_doublet_rate : Double
  doublet_rate_per_1000 : Double
  rate_uncertainty : Double
  artificial_doublets : Int
  selected_features : Int
  dimensions : Int
  neighbors : Int
  iterations : Int
  classifier_steps : Int
  learning_rate : Double
  regularization : Double
  stringency : Double
  random_fraction : Double
  half_size_fraction : Double
  unidentifiable_threshold : Double
  cluster_count : Int
  seed : Int
} derive(Debug)
pub extend ScDblFinderConfig with @moonbitlang/core/debug.Debug::{to_repr}

///|
pub struct ScDblFinderArtificialData {
  counts : Array[Array[Double]]
  names : Array[String]
  parent_one : Array[Int]
  parent_two : Array[Int]
  origins : Array[String]
} derive(Debug)
pub extend ScDblFinderArtificialData with @moonbitlang/core/debug.Debug::{to_repr}

///|
pub struct ScDblFinderThreshold {
  threshold : Double
  cost : Double
  false_positive_rate : Double
  false_negative_rate : Double
  called_rate : Double
} derive(Debug)
pub extend ScDblFinderThreshold with @moonbitlang/core/debug.Debug::{to_repr}

///|
pub struct ScDblFinderCellResult {
  cell_name : String
  sample : String
  cluster : String
  score : Double
  is_doublet : Bool
  weighted_ratio : Double
  neighbor_ratio : Double
  cxds_score : Double
  nearest_class : String
  origin : String
  origin_ambiguity : Double
} derive(Debug)
pub extend ScDblFinderCellResult with @moonbitlang/core/debug.Debug::{to_repr}

///|
pub struct ScDblFinderResult {
  cell_names : Array[String]
  gene_names : Array[String]
  samples : Array[String]
  clusters : Array[String]
  known_doublets : Array[Bool]
  scores : Array[Double]
  calls : Array[Bool]
  weighted_ratios : Array[Double]
  neighbor_ratios : Array[Double]
  cxds_scores : Array[Double]
  distances_to_nearest : Array[Double]
  distances_to_doublet : Array[Double]
  distances_to_real : Array[Double]
  nearest_classes : Array[String]
  origins : Array[String]
  origin_ambiguities : Array[Double]
  library_sizes : Array[Double]
  detected_features : Array[Int]
  features_above_two : Array[Int]
  pca : Array[Array[Double]]
  selected_gene_names : Array[String]
  sample_ids : Array[String]
  thresholds : Array[ScDblFinderThreshold]
  expected_doublet_rates : Array[Double]
  artificial_doublets : Int
  classifier_iterations : Int
  config : ScDblFinderConfig
} derive(Debug)
pub extend ScDblFinderResult with @moonbitlang/core/debug.Debug::{to_repr}

///|
pub struct ScDblFinderEnrichment {
  combination : String
  observed : Int
  expected : Double
  log2_enrichment : Double
  p_value : Double
  adjusted_p_value : Double
} derive(Debug)
pub extend ScDblFinderEnrichment with @moonbitlang/core/debug.Debug::{to_repr}

///|
pub struct ScDblFinderSingleCellExperimentOutput {
  experiment : SingleCellExperiment
  result : ScDblFinderResult
}

///|
priv struct ScdfaKnnMetrics {
  weighted : Array[Double]
  ratios : Array[Double]
  nearest : Array[Double]
  nearest_doublet : Array[Double]
  nearest_real : Array[Double]
  nearest_class : Array[String]
  origins : Array[String]
  ambiguities : Array[Double]
}

///|
priv struct ScdfaCaptureFit {
  global_indices : Array[Int]
  scores : Array[Double]
  calls : Array[Bool]
  weighted_ratios : Array[Double]
  neighbor_ratios : Array[Double]
  cxds_scores : Array[Double]
  distances_to_nearest : Array[Double]
  distances_to_doublet : Array[Double]
  distances_to_real : Array[Double]
  nearest_classes : Array[String]
  origins : Array[String]
  ambiguities : Array[Double]
  library_sizes : Array[Double]
  detected_features : Array[Int]
  features_above_two : Array[Int]
  pca : Array[Array[Double]]
  selected_indices : Array[Int]
  threshold : ScDblFinderThreshold
  expected_rate : Double
  artificial_count : Int
}

///|
fn scdfa_fail(message : String) -> Unit raise ScDblFinderError {
  raise ScDblFinderError(message)
}

///|
fn scdfa_finite(value : Double) -> Bool {
  value == value && value.abs() <= 1.0e300
}

///|
fn scdfa_copy_matrix(values : Array[Array[Double]]) -> Array[Array[Double]] {
  let result = []
  for row in values {
    result.push(row.copy())
  }
  result
}

///|
fn scdfa_contains_string(values : Array[String], query : String) -> Bool {
  for value in values {
    if value == query {
      return true
    }
  }
  false
}

///|
fn scdfa_string_index(values : Array[String], query : String) -> Int {
  for index in 0..<values.length() {
    if values[index] == query {
      return index
    }
  }
  -1
}

///|
fn scdfa_unique_strings(values : Array[String]) -> Array[String] {
  let result = []
  for value in values {
    if !scdfa_contains_string(result, value) {
      result.push(value)
    }
  }
  result
}

///|
fn scdfa_zero_matrix(rows : Int, columns : Int) -> Array[Array[Double]] {
  let result = []
  for _ in 0..<rows {
    result.push(Array::make(columns, 0.0))
  }
  result
}

///|
fn scdfa_median(values : Array[Double]) -> Double {
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
  let middle = sorted.length() / 2
  if sorted.length() % 2 == 0 {
    (sorted[middle - 1] + sorted[middle]) / 2.0
  } else {
    sorted[middle]
  }
}

///|
fn scdfa_quantile(values : Array[Double], probability : Double) -> Double {
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
  let location = probability * (sorted.length() - 1).to_double()
  let lower = location.to_int()
  let upper = if lower + 1 < sorted.length() { lower + 1 } else { lower }
  let fraction = location - lower.to_double()
  sorted[lower] * (1.0 - fraction) + sorted[upper] * fraction
}

///|
fn scdfa_log2(value : Double) -> Double {
  @math.ln(value) / 0.6931471805599453
}

///|
fn scdfa_sigmoid(value : Double) -> Double {
  if value >= 35.0 {
    1.0
  } else if value <= -35.0 {
    0.0
  } else {
    1.0 / (1.0 + @math.exp(-value))
  }
}

///|
fn scdfa_validate_names(
  values : Array[String],
  expected : Int,
  label : String,
) -> Unit raise ScDblFinderError {
  if values.length() != expected {
    scdfa_fail(
      "scDblFinder " + label + " length must be " + expected.to_string(),
    )
  }
  let seen : Array[String] = []
  for value in values {
    if value == "" {
      scdfa_fail("scDblFinder " + label + " cannot contain empty names")
    }
    if scdfa_contains_string(seen, value) {
      scdfa_fail("scDblFinder " + label + " contains duplicate: " + value)
    }
    seen.push(value)
  }
}

///|
fn scdfa_validate_data(data : SingleCellData) -> Unit raise ScDblFinderError {
  let cells = data.counts.length()
  if cells < 3 {
    scdfa_fail("scDblFinder requires at least three cells")
  }
  if data.gene_names.length() < 2 {
    scdfa_fail("scDblFinder requires at least two genes")
  }
  scdfa_validate_names(data.cell_names, cells, "cell names")
  scdfa_validate_names(data.gene_names, data.gene_names.length(), "gene names")
  for cell in 0..<cells {
    if data.counts[cell].length() != data.gene_names.length() {
      scdfa_fail("scDblFinder count matrix must be rectangular (cell x gene)")
    }
    let mut total = 0.0
    for value in data.counts[cell] {
      if !scdfa_finite(value) || value < 0.0 {
        scdfa_fail("scDblFinder counts must contain finite non-negative values")
      }
      total = total + value
    }
    if total <= 0.0 {
      scdfa_fail(
        "scDblFinder cell " + data.cell_names[cell] + " has zero library size",
      )
    }
  }
}

///|
pub fn SingleCellData::create(
  counts : Array[Array[Double]],
  cell_names : Array[String],
  gene_names : Array[String],
) -> SingleCellData raise ScDblFinderError {
  let result = SingleCellData::new(
    scdfa_copy_matrix(counts),
    cell_names.copy(),
    gene_names.copy(),
  )
  scdfa_validate_data(result)
  result
}

///|
pub fn ScDblFinderConfig::create(
  expected_doublet_rate? : Double = -1.0,
  doublet_rate_per_1000? : Double = 0.008,
  rate_uncertainty? : Double = -1.0,
  artificial_doublets? : Int = 0,
  selected_features? : Int = 1000,
  dimensions? : Int = 20,
  neighbors? : Int = 15,
  iterations? : Int = 3,
  classifier_steps? : Int = 300,
  learning_rate? : Double = 0.15,
  regularization? : Double = 0.01,
  stringency? : Double = 0.5,
  random_fraction? : Double = 0.1,
  half_size_fraction? : Double = 0.25,
  unidentifiable_threshold? : Double = 0.1,
  cluster_count? : Int = 0,
  seed? : Int = 1,
) -> ScDblFinderConfig raise ScDblFinderError {
  if !scdfa_finite(expected_doublet_rate) ||
    expected_doublet_rate > 1.0 ||
    (expected_doublet_rate < 0.0 && expected_doublet_rate != -1.0) {
    scdfa_fail("scDblFinder expected_doublet_rate must be -1 or in [0, 1]")
  }
  if !scdfa_finite(doublet_rate_per_1000) ||
    doublet_rate_per_1000 < 0.0 ||
    doublet_rate_per_1000 > 1.0 {
    scdfa_fail("scDblFinder doublet_rate_per_1000 must be in [0, 1]")
  }
  if !scdfa_finite(rate_uncertainty) ||
    rate_uncertainty > 1.0 ||
    (rate_uncertainty < 0.0 && rate_uncertainty != -1.0) {
    scdfa_fail("scDblFinder rate_uncertainty must be -1 or in [0, 1]")
  }
  if artificial_doublets < 0 {
    scdfa_fail("scDblFinder artificial_doublets cannot be negative")
  }
  if selected_features < 2 || dimensions < 1 || neighbors < 1 {
    scdfa_fail(
      "scDblFinder selected_features, dimensions, and neighbors must be positive",
    )
  }
  if iterations < 1 || classifier_steps < 1 {
    scdfa_fail("scDblFinder iterations and classifier_steps must be positive")
  }
  if !scdfa_finite(learning_rate) || learning_rate <= 0.0 || learning_rate > 1.0 {
    scdfa_fail("scDblFinder learning_rate must be in (0, 1]")
  }
  if !scdfa_finite(regularization) || regularization < 0.0 {
    scdfa_fail("scDblFinder regularization cannot be negative")
  }
  if !scdfa_finite(stringency) || stringency <= 0.0 || stringency >= 1.0 {
    scdfa_fail("scDblFinder stringency must be in (0, 1)")
  }
  if !scdfa_finite(random_fraction) ||
    random_fraction < 0.0 ||
    random_fraction > 1.0 {
    scdfa_fail("scDblFinder random_fraction must be in [0, 1]")
  }
  if !scdfa_finite(half_size_fraction) ||
    half_size_fraction < 0.0 ||
    half_size_fraction > 1.0 {
    scdfa_fail("scDblFinder half_size_fraction must be in [0, 1]")
  }
  if !scdfa_finite(unidentifiable_threshold) ||
    unidentifiable_threshold < 0.0 ||
    unidentifiable_threshold > 1.0 {
    scdfa_fail("scDblFinder unidentifiable_threshold must be in [0, 1]")
  }
  if cluster_count < 0 || cluster_count == 1 {
    scdfa_fail("scDblFinder cluster_count must be 0 or at least 2")
  }
  ScDblFinderConfig::{
    expected_doublet_rate,
    doublet_rate_per_1000,
    rate_uncertainty,
    artificial_doublets,
    selected_features,
    dimensions,
    neighbors,
    iterations,
    classifier_steps,
    learning_rate,
    regularization,
    stringency,
    random_fraction,
    half_size_fraction,
    unidentifiable_threshold,
    cluster_count,
    seed,
  }
}

///|
pub fn ScDblFinderConfig::default() -> ScDblFinderConfig {
  ScDblFinderConfig::{
    expected_doublet_rate: -1.0,
    doublet_rate_per_1000: 0.008,
    rate_uncertainty: -1.0,
    artificial_doublets: 0,
    selected_features: 1000,
    dimensions: 20,
    neighbors: 15,
    iterations: 3,
    classifier_steps: 300,
    learning_rate: 0.15,
    regularization: 0.01,
    stringency: 0.5,
    random_fraction: 0.1,
    half_size_fraction: 0.25,
    unidentifiable_threshold: 0.1,
    cluster_count: 0,
    seed: 1,
  }
}

///|
pub fn scdf_expected_doublet_rate(
  cells : Int,
  rate_per_1000? : Double = 0.008,
) -> Double raise ScDblFinderError {
  if cells < 1 {
    scdfa_fail("scDblFinder cell count must be positive")
  }
  if !scdfa_finite(rate_per_1000) || rate_per_1000 < 0.0 || rate_per_1000 > 1.0 {
    scdfa_fail("scDblFinder rate_per_1000 must be in [0, 1]")
  }
  let result = rate_per_1000 * cells.to_double() / 1000.0
  if result > 1.0 {
    1.0
  } else {
    result
  }
}

///|
pub fn scdf_homotypic_proportion(
  clusters : Array[String],
) -> Double raise ScDblFinderError {
  if clusters.length() == 0 {
    scdfa_fail("scDblFinder clusters cannot be empty")
  }
  for cluster in clusters {
    if cluster == "" {
      scdfa_fail("scDblFinder clusters cannot contain empty labels")
    }
  }
  let labels = scdfa_unique_strings(clusters)
  let mut result = 0.0
  for label in labels {
    let mut count = 0
    for cluster in clusters {
      if cluster == label {
        count = count + 1
      }
    }
    let proportion = count.to_double() / clusters.length().to_double()
    result = result + proportion * proportion
  }
  result
}

///|
fn scdfa_library_sizes(counts : Array[Array[Double]]) -> Array[Double] {
  let result = []
  for row in counts {
    let mut total = 0.0
    for value in row {
      total = total + value
    }
    result.push(total)
  }
  result
}

///|
fn scdfa_detected_features(
  counts : Array[Array[Double]],
  threshold : Double,
) -> Array[Int] {
  let result = []
  for row in counts {
    let mut total = 0
    for value in row {
      if value > threshold {
        total = total + 1
      }
    }
    result.push(total)
  }
  result
}

///|
fn scdfa_selected_features(
  counts : Array[Array[Double]],
  number : Int,
) -> Array[Int] {
  let cells = counts.length()
  let genes = counts[0].length()
  let libraries = scdfa_library_sizes(counts)
  let positive = []
  for value in libraries {
    if value > 0.0 {
      positive.push(value)
    }
  }
  let median_library = scdfa_median(positive)
  let variances = Array::make(genes, 0.0)
  for gene in 0..<genes {
    let values = Array::make(cells, 0.0)
    let mut mean = 0.0
    for cell in 0..<cells {
      let normalized = counts[cell][gene] * median_library / libraries[cell]
      let value = scdfa_log2(normalized + 1.0)
      values[cell] = value
      mean = mean + value
    }
    mean = mean / cells.to_double()
    let mut variance = 0.0
    for value in values {
      let difference = value - mean
      variance = variance + difference * difference
    }
    variances[gene] = variance / (cells - 1).to_double()
  }
  let indices = []
  for gene in 0..<genes {
    indices.push(gene)
  }
  for index in 1..<indices.length() {
    let value = indices[index]
    let mut position = index
    while position > 0 &&
          (
            variances[indices[position - 1]] < variances[value] ||
            (
              variances[indices[position - 1]] == variances[value] &&
              indices[position - 1] > value
            )
          ) {
      indices[position] = indices[position - 1]
      position = position - 1
    }
    indices[position] = value
  }
  let limit = if number < genes { number } else { genes }
  let selected = []
  for index in 0..<limit {
    selected.push(indices[index])
  }
  selected
}

///|
fn scdfa_processed_coordinates(
  counts : Array[Array[Double]],
  real_cells : Int,
  selected : Array[Int],
  dimensions : Int,
) -> Array[Array[Double]] {
  let libraries = scdfa_library_sizes(counts)
  let real_libraries = []
  for cell in 0..<real_cells {
    real_libraries.push(libraries[cell])
  }
  let median_library = scdfa_median(real_libraries)
  let logcounts = scdfa_zero_matrix(selected.length(), counts.length())
  for feature in 0..<selected.length() {
    let gene = selected[feature]
    for cell in 0..<counts.length() {
      let normalized = counts[cell][gene] * median_library / libraries[cell]
      logcounts[feature][cell] = scdfa_log2(normalized + 1.0)
    }
  }
  let mut components = dimensions
  if components > selected.length() {
    components = selected.length()
  }
  if components >= counts.length() {
    components = counts.length() - 1
  }
  if components < 1 {
    components = 1
  }
  sce_compute_pca(logcounts, components).coordinates
}

///|
fn scdfa_origin_name(
  first : String,
  second : String,
  labels : Array[String],
) -> String {
  let left = scdfa_string_index(labels, first)
  let right = scdfa_string_index(labels, second)
  if left <= right {
    first + "+" + second
  } else {
    second + "+" + first
  }
}

///|
pub fn scdf_generate_artificial_doublets(
  data : SingleCellData,
  number : Int,
  clusters? : Array[String] = [],
  random_fraction? : Double = 0.1,
  half_size_fraction? : Double = 0.25,
  seed? : Int = 1,
) -> ScDblFinderArtificialData raise ScDblFinderError {
  scdfa_validate_data(data)
  if number < 1 {
    scdfa_fail("scDblFinder artificial doublet number must be positive")
  }
  if clusters.length() > 0 && clusters.length() != data.counts.length() {
    scdfa_fail("scDblFinder cluster length must match cells")
  }
  if !scdfa_finite(random_fraction) ||
    random_fraction < 0.0 ||
    random_fraction > 1.0 {
    scdfa_fail("scDblFinder random_fraction must be in [0, 1]")
  }
  if !scdfa_finite(half_size_fraction) ||
    half_size_fraction < 0.0 ||
    half_size_fraction > 1.0 {
    scdfa_fail("scDblFinder half_size_fraction must be in [0, 1]")
  }
  let libraries = scdfa_library_sizes(data.counts)
  let eligible : Array[Int] = []
  let lower = if libraries.length() >= 10 {
    scdfa_quantile(libraries, 0.05)
  } else {
    0.0
  }
  let upper = if libraries.length() >= 10 {
    scdfa_quantile(libraries, 0.95)
  } else {
    1.0e300
  }
  for cell in 0..<libraries.length() {
    if libraries[cell] > 0.0 &&
      libraries[cell] >= lower &&
      libraries[cell] <= upper {
      eligible.push(cell)
    }
  }
  if eligible.length() < 2 {
    scdfa_fail("scDblFinder has fewer than two eligible cells")
  }
  let random_candidates : Array[(Int, Int)] = []
  let cross_candidates : Array[(Int, Int)] = []
  for left in 0..<eligible.length() {
    for right in (left + 1)..<eligible.length() {
      let pair = (eligible[left], eligible[right])
      random_candidates.push(pair)
      if clusters.length() > 0 && clusters[pair.0] != clusters[pair.1] {
        cross_candidates.push(pair)
      }
    }
  }
  let use_clusters = clusters.length() > 0 && cross_candidates.length() > 0
  let labels = if use_clusters { scdfa_unique_strings(clusters) } else { [] }
  let random_number = if use_clusters {
    (number.to_double() * random_fraction).to_int()
  } else {
    number
  }
  let half_number = (number.to_double() * half_size_fraction).to_int()
  let counts = []
  let names = []
  let parent_one = []
  let parent_two = []
  let origins = []
  let offset = if seed < 0 { -seed } else { seed }
  for artificial in 0..<number {
    let random = artificial < random_number
    let candidates = if random { random_candidates } else { cross_candidates }
    let candidate_index = (
        offset + artificial * 37 + artificial * artificial * 11
      ) %
      candidates.length()
    let pair = candidates[candidate_index]
    let row = []
    let scale = if artificial < half_number { 0.5 } else { 1.0 }
    for gene in 0..<data.gene_names.length() {
      row.push((data.counts[pair.0][gene] + data.counts[pair.1][gene]) * scale)
    }
    counts.push(row)
    names.push("artificial_" + (artificial + 1).to_string())
    parent_one.push(pair.0)
    parent_two.push(pair.1)
    origins.push(
      if use_clusters {
        scdfa_origin_name(clusters[pair.0], clusters[pair.1], labels)
      } else {
        ""
      },
    )
  }
  ScDblFinderArtificialData::{ counts, names, parent_one, parent_two, origins, }
}

///|
fn scdfa_distance(first : Array[Double], second : Array[Double]) -> Double {
  let mut total = 0.0
  for index in 0..<first.length() {
    let difference = first[index] - second[index]
    total = total + difference * difference
  }
  total.sqrt()
}

///|
fn scdfa_insert_neighbor(
  distances : Array[Double],
  indices : Array[Int],
  distance : Double,
  index : Int,
) -> Unit {
  let mut position = distances.length()
  for current in 0..<distances.length() {
    if distance < distances[current] ||
      (distance == distances[current] && index < indices[current]) {
      position = current
      break
    }
  }
  if position < distances.length() {
    let mut cursor = distances.length() - 1
    while cursor > position {
      distances[cursor] = distances[cursor - 1]
      indices[cursor] = indices[cursor - 1]
      cursor = cursor - 1
    }
    distances[position] = distance
    indices[position] = index
  }
}

///|
fn scdfa_knn_metrics(
  coordinates : Array[Array[Double]],
  real_cells : Int,
  artificial_origins : Array[String],
  neighbors : Int,
) -> ScdfaKnnMetrics {
  let total = coordinates.length()
  let k = if neighbors < total - 1 { neighbors } else { total - 1 }
  let weighted = Array::make(total, 0.0)
  let ratios = Array::make(total, 0.0)
  let nearest = Array::make(total, 0.0)
  let nearest_doublet = Array::make(total, 0.0)
  let nearest_real = Array::make(total, 0.0)
  let nearest_class = Array::make(total, "")
  let origins = Array::make(total, "")
  let ambiguities = Array::make(total, 1.0)
  let origin_labels = scdfa_unique_strings(artificial_origins)
  for cell in 0..<total {
    let distances = Array::make(k, 1.0e300)
    let indices = Array::make(k, -1)
    let mut closest_artificial = 1.0e300
    let mut closest_real = 1.0e300
    for other in 0..<total {
      if other != cell {
        let distance = scdfa_distance(coordinates[cell], coordinates[other])
        scdfa_insert_neighbor(distances, indices, distance, other)
        if other < real_cells {
          if distance < closest_real {
            closest_real = distance
          }
        } else if distance < closest_artificial {
          closest_artificial = distance
        }
      }
    }
    nearest[cell] = distances[0]
    nearest_doublet[cell] = closest_artificial
    nearest_real[cell] = closest_real
    nearest_class[cell] = if indices[0] >= real_cells {
      "artificialDoublet"
    } else {
      "cell"
    }
    let votes = Array::make(origin_labels.length(), 0.0)
    let mut artificial_count = 0
    let mut weighted_artificial = 0.0
    let mut total_weight = 0.0
    let mut origin_total = 0.0
    for rank in 0..<k {
      let distance = if distances[rank] > 1.0e-12 {
        distances[rank]
      } else {
        1.0e-12
      }
      let rank_weight = (k - rank).to_double().sqrt() / distance
      total_weight = total_weight + rank_weight
      if indices[rank] >= real_cells {
        artificial_count = artificial_count + 1
        weighted_artificial = weighted_artificial + rank_weight
        let origin = artificial_origins[indices[rank] - real_cells]
        let origin_index = scdfa_string_index(origin_labels, origin)
        if origin_index >= 0 && origin != "" {
          votes[origin_index] = votes[origin_index] + rank_weight
          origin_total = origin_total + rank_weight
        }
      }
    }
    ratios[cell] = artificial_count.to_double() / k.to_double()
    weighted[cell] = if total_weight > 0.0 {
      weighted_artificial / total_weight
    } else {
      0.0
    }
    if origin_total > 0.0 {
      let mut best = 0
      for origin in 1..<votes.length() {
        if votes[origin] > votes[best] {
          best = origin
        }
      }
      origins[cell] = origin_labels[best]
      ambiguities[cell] = 1.0 - votes[best] / origin_total
    }
  }
  ScdfaKnnMetrics::{
    weighted,
    ratios,
    nearest,
    nearest_doublet,
    nearest_real,
    nearest_class,
    origins,
    ambiguities,
  }
}

///|
fn scdfa_cxds_scores(
  counts : Array[Array[Double]],
  real_cells : Int,
  selected : Array[Int],
) -> Array[Double] {
  let feature_count = if selected.length() < 64 {
    selected.length()
  } else {
    64
  }
  let probabilities = Array::make(feature_count, 0.0)
  for feature in 0..<feature_count {
    let gene = selected[feature]
    let mut expressed = 0
    for cell in 0..<real_cells {
      if counts[cell][gene] >= 1.0 {
        expressed = expressed + 1
      }
    }
    probabilities[feature] = expressed.to_double() / real_cells.to_double()
  }
  let pair_scores = scdfa_zero_matrix(feature_count, feature_count)
  for left in 0..<feature_count {
    for right in (left + 1)..<feature_count {
      let mut observed = 0
      for cell in 0..<real_cells {
        if counts[cell][selected[left]] >= 1.0 &&
          counts[cell][selected[right]] >= 1.0 {
          observed = observed + 1
        }
      }
      let observed_rate = (observed.to_double() + 0.5) /
        (real_cells.to_double() + 1.0)
      let expected_rate = probabilities[left] * probabilities[right]
      let surprise = if expected_rate > observed_rate {
        @math.ln((expected_rate + 1.0e-12) / (observed_rate + 1.0e-12))
      } else {
        0.0
      }
      pair_scores[left][right] = surprise
      pair_scores[right][left] = surprise
    }
  }
  let raw = Array::make(counts.length(), 0.0)
  for cell in 0..<counts.length() {
    let mut total = 0.0
    for left in 0..<feature_count {
      if counts[cell][selected[left]] >= 1.0 {
        for right in (left + 1)..<feature_count {
          if counts[cell][selected[right]] >= 1.0 {
            total = total + pair_scores[left][right]
          }
        }
      }
    }
    raw[cell] = total
  }
  let mut minimum = raw[0]
  let mut maximum = raw[0]
  for value in raw {
    if value < minimum {
      minimum = value
    }
    if value > maximum {
      maximum = value
    }
  }
  if maximum > minimum {
    for index in 0..<raw.length() {
      raw[index] = (raw[index] - minimum) / (maximum - minimum)
    }
  } else {
    for index in 0..<raw.length() {
      raw[index] = 0.0
    }
  }
  raw
}

///|
fn scdfa_classifier_features(
  counts : Array[Array[Double]],
  coordinates : Array[Array[Double]],
  metrics : ScdfaKnnMetrics,
  cxds : Array[Double],
) -> Array[Array[Double]] {
  let libraries = scdfa_library_sizes(counts)
  let detected = scdfa_detected_features(counts, 0.0)
  let above_two = scdfa_detected_features(counts, 2.0)
  let extra_pcs = if coordinates[0].length() < 5 {
    coordinates[0].length()
  } else {
    5
  }
  let features = []
  for cell in 0..<counts.length() {
    let proximity = metrics.nearest_real[cell] /
      (metrics.nearest_real[cell] + metrics.nearest_doublet[cell] + 1.0e-12)
    let row = [
      metrics.weighted[cell],
      metrics.ratios[cell],
      proximity,
      @math.ln(libraries[cell] + 1.0),
      detected[cell].to_double(),
      above_two[cell].to_double(),
      cxds[cell],
    ]
    for component in 0..<extra_pcs {
      row.push(coordinates[cell][component])
    }
    features.push(row)
  }
  features
}

///|
fn scdfa_standardize(features : Array[Array[Double]]) -> Array[Array[Double]] {
  let rows = features.length()
  let columns = features[0].length()
  let means = Array::make(columns, 0.0)
  let scales = Array::make(columns, 0.0)
  for column in 0..<columns {
    for row in 0..<rows {
      means[column] = means[column] + features[row][column]
    }
    means[column] = means[column] / rows.to_double()
    for row in 0..<rows {
      let difference = features[row][column] - means[column]
      scales[column] = scales[column] + difference * difference
    }
    scales[column] = (scales[column] / rows.to_double()).sqrt()
    if scales[column] < 1.0e-12 {
      scales[column] = 1.0
    }
  }
  let result = scdfa_zero_matrix(rows, columns)
  for row in 0..<rows {
    for column in 0..<columns {
      result[row][column] = (features[row][column] - means[column]) /
        scales[column]
    }
  }
  result
}

///|
fn scdfa_train_classifier(
  features : Array[Array[Double]],
  labels : Array[Bool],
  included : Array[Bool],
  steps : Int,
  learning_rate : Double,
  regularization : Double,
) -> Array[Double] {
  let columns = features[0].length()
  let coefficients = Array::make(columns, 0.0)
  let mut intercept = 0.0
  let mut positives = 0
  let mut negatives = 0
  for row in 0..<features.length() {
    if included[row] {
      if labels[row] {
        positives = positives + 1
      } else {
        negatives = negatives + 1
      }
    }
  }
  if positives == 0 || negatives == 0 {
    return Array::make(features.length(), 0.5)
  }
  for step in 0..<steps {
    let gradient = Array::make(columns, 0.0)
    let mut intercept_gradient = 0.0
    for row in 0..<features.length() {
      if included[row] {
        let mut linear = intercept
        for column in 0..<columns {
          linear = linear + coefficients[column] * features[row][column]
        }
        let probability = scdfa_sigmoid(linear)
        let target = if labels[row] { 1.0 } else { 0.0 }
        let class_weight = if labels[row] {
          0.5 / positives.to_double()
        } else {
          0.5 / negatives.to_double()
        }
        let error = (probability - target) * class_weight
        intercept_gradient = intercept_gradient + error
        for column in 0..<columns {
          gradient[column] = gradient[column] + error * features[row][column]
        }
      }
    }
    let step_size = learning_rate / (1.0 + step.to_double() * 0.002)
    intercept = intercept - step_size * intercept_gradient
    for column in 0..<columns {
      coefficients[column] = coefficients[column] -
        step_size * (gradient[column] + regularization * coefficients[column])
    }
  }
  let result = Array::make(features.length(), 0.0)
  for row in 0..<features.length() {
    let mut linear = intercept
    for column in 0..<columns {
      linear = linear + coefficients[column] * features[row][column]
    }
    result[row] = scdfa_sigmoid(linear)
  }
  result
}

///|
fn scdfa_threshold_cost(
  real_scores : Array[Double],
  artificial_scores : Array[Double],
  threshold : Double,
  expected_rate : Double,
  uncertainty : Double,
  stringency : Double,
) -> ScDblFinderThreshold {
  let mut called_real = 0
  for score in real_scores {
    if score >= threshold {
      called_real = called_real + 1
    }
  }
  let mut missed_artificial = 0
  for score in artificial_scores {
    if score < threshold {
      missed_artificial = missed_artificial + 1
    }
  }
  let called_rate = called_real.to_double() / real_scores.length().to_double()
  let lower = if expected_rate > uncertainty {
    expected_rate - uncertainty
  } else {
    0.0
  }
  let upper = if expected_rate + uncertainty < 1.0 {
    expected_rate + uncertainty
  } else {
    1.0
  }
  let deviation = if called_rate >= lower && called_rate <= upper {
    0.0
  } else if called_rate < lower {
    (lower - called_rate) / (lower + 1.0 / real_scores.length().to_double())
  } else {
    (called_rate - upper) / (upper + 1.0 / real_scores.length().to_double())
  }
  let false_positive_rate = called_rate
  let false_negative_rate = missed_artificial.to_double() /
    artificial_scores.length().to_double()
  ScDblFinderThreshold::{
    threshold,
    cost: deviation * deviation +
    2.0 * stringency * false_positive_rate +
    2.0 * (1.0 - stringency) * false_negative_rate,
    false_positive_rate,
    false_negative_rate,
    called_rate,
  }
}

///|
pub fn scdf_optimize_threshold(
  real_scores : Array[Double],
  artificial_scores : Array[Double],
  expected_rate : Double,
  uncertainty? : Double = -1.0,
  stringency? : Double = 0.5,
) -> ScDblFinderThreshold raise ScDblFinderError {
  if real_scores.length() == 0 || artificial_scores.length() == 0 {
    scdfa_fail("scDblFinder thresholding requires both score classes")
  }
  if !scdfa_finite(expected_rate) || expected_rate < 0.0 || expected_rate > 1.0 {
    scdfa_fail("scDblFinder expected rate must be in [0, 1]")
  }
  if !scdfa_finite(uncertainty) ||
    uncertainty > 1.0 ||
    (uncertainty < 0.0 && uncertainty != -1.0) {
    scdfa_fail("scDblFinder uncertainty must be -1 or in [0, 1]")
  }
  if !scdfa_finite(stringency) || stringency <= 0.0 || stringency >= 1.0 {
    scdfa_fail("scDblFinder stringency must be in (0, 1)")
  }
  let all_scores = []
  for score in real_scores {
    if !scdfa_finite(score) || score < 0.0 || score > 1.0 {
      scdfa_fail("scDblFinder scores must be finite and in [0, 1]")
    }
    all_scores.push(score)
  }
  for score in artificial_scores {
    if !scdfa_finite(score) || score < 0.0 || score > 1.0 {
      scdfa_fail("scDblFinder scores must be finite and in [0, 1]")
    }
    all_scores.push(score)
  }
  for index in 1..<all_scores.length() {
    let value = all_scores[index]
    let mut position = index
    while position > 0 && all_scores[position - 1] > value {
      all_scores[position] = all_scores[position - 1]
      position = position - 1
    }
    all_scores[position] = value
  }
  let effective_uncertainty = if uncertainty >= 0.0 {
    uncertainty
  } else {
    0.3 * expected_rate + 0.025
  }
  let mut best = scdfa_threshold_cost(
    real_scores, artificial_scores, 0.0, expected_rate, effective_uncertainty, stringency,
  )
  for index in 0..<all_scores.length() {
    let candidate = if index + 1 < all_scores.length() {
      (all_scores[index] + all_scores[index + 1]) / 2.0
    } else {
      all_scores[index]
    }
    let current = scdfa_threshold_cost(
      real_scores, artificial_scores, candidate, expected_rate, effective_uncertainty,
      stringency,
    )
    if current.cost < best.cost ||
      (current.cost == best.cost && current.threshold > best.threshold) {
      best = current
    }
  }
  best
}

///|
fn scdfa_top_indices(
  values : Array[Double],
  candidates : Array[Int],
  limit : Int,
  descending : Bool,
) -> Array[Int] {
  let ordered = candidates.copy()
  for index in 1..<ordered.length() {
    let value = ordered[index]
    let mut position = index
    while position > 0 &&
          (if descending {
            values[ordered[position - 1]] < values[value]
          } else {
            values[ordered[position - 1]] > values[value]
          }) {
      ordered[position] = ordered[position - 1]
      position = position - 1
    }
    ordered[position] = value
  }
  let result = []
  let count = if limit < ordered.length() { limit } else { ordered.length() }
  for index in 0..<count {
    result.push(ordered[index])
  }
  result
}

///|
fn scdfa_fit_capture(
  data : SingleCellData,
  global_indices : Array[Int],
  clusters : Array[String],
  known_doublets : Array[Bool],
  config : ScDblFinderConfig,
) -> ScdfaCaptureFit raise ScDblFinderError {
  let local_counts = []
  let local_names = []
  let local_clusters = []
  let local_known = []
  for index in global_indices {
    local_counts.push(data.counts[index].copy())
    local_names.push(data.cell_names[index])
    if clusters.length() > 0 {
      local_clusters.push(clusters[index])
    }
    if known_doublets.length() > 0 {
      local_known.push(known_doublets[index])
    } else {
      local_known.push(false)
    }
  }
  let local_data = SingleCellData::create(
    local_counts,
    local_names,
    data.gene_names,
  )
  let selected = scdfa_selected_features(
    local_data.counts,
    config.selected_features,
  )
  let artificial_number = if config.artificial_doublets > 0 {
    config.artificial_doublets
  } else {
    let proposed = local_data.counts.length() * 2
    if proposed < 150 {
      150
    } else if proposed > 1500 {
      1500
    } else {
      proposed
    }
  }
  let artificial = scdf_generate_artificial_doublets(
    local_data,
    artificial_number,
    clusters=local_clusters,
    random_fraction=config.random_fraction,
    half_size_fraction=config.half_size_fraction,
    seed=config.seed + global_indices[0],
  )
  let combined = scdfa_copy_matrix(local_data.counts)
  for row in artificial.counts {
    combined.push(row.copy())
  }
  let real_cells = local_data.counts.length()
  let coordinates = scdfa_processed_coordinates(
    combined,
    real_cells,
    selected,
    config.dimensions,
  )
  let metrics = scdfa_knn_metrics(
    coordinates,
    real_cells,
    artificial.origins,
    config.neighbors,
  )
  let cxds = scdfa_cxds_scores(combined, real_cells, selected)
  let raw_features = scdfa_classifier_features(
    combined, coordinates, metrics, cxds,
  )
  let features = scdfa_standardize(raw_features)
  let labels = Array::make(combined.length(), false)
  let included = Array::make(combined.length(), true)
  for cell in 0..<real_cells {
    labels[cell] = local_known[cell]
  }
  for cell in real_cells..<combined.length() {
    labels[cell] = true
  }
  let expected_rate = if config.expected_doublet_rate >= 0.0 {
    config.expected_doublet_rate
  } else {
    scdf_expected_doublet_rate(
      real_cells,
      rate_per_1000=config.doublet_rate_per_1000,
    )
  }
  let heterotypic_rate = if local_clusters.length() > 0 {
    expected_rate * (1.0 - scdf_homotypic_proportion(local_clusters))
  } else {
    expected_rate
  }
  let mut scores = Array::make(combined.length(), 0.0)
  let mut threshold = ScDblFinderThreshold::{
    threshold: 0.5,
    cost: 0.0,
    false_positive_rate: 0.0,
    false_negative_rate: 0.0,
    called_rate: 0.0,
  }
  for iteration in 0..<config.iterations {
    scores = scdfa_train_classifier(
      features,
      labels,
      included,
      config.classifier_steps,
      config.learning_rate,
      config.regularization,
    )
    let real_scores = []
    let artificial_scores = []
    for cell in 0..<real_cells {
      real_scores.push(scores[cell])
    }
    for cell in real_cells..<scores.length() {
      artificial_scores.push(scores[cell])
    }
    threshold = scdf_optimize_threshold(
      real_scores,
      artificial_scores,
      heterotypic_rate,
      uncertainty=config.rate_uncertainty,
      stringency=config.stringency,
    )
    if iteration + 1 < config.iterations {
      for cell in 0..<included.length() {
        included[cell] = true
      }
      let suspect_real = []
      for cell in 0..<real_cells {
        if scores[cell] >= threshold.threshold && !local_known[cell] {
          suspect_real.push(cell)
        }
      }
      let real_limit = (real_cells.to_double() * 0.2).to_int()
      let excluded_real = scdfa_top_indices(
        scores, suspect_real, real_limit, true,
      )
      for cell in excluded_real {
        included[cell] = false
      }
      let weak_artificial = []
      for cell in real_cells..<combined.length() {
        if scores[cell] < config.unidentifiable_threshold {
          weak_artificial.push(cell)
        }
      }
      let artificial_limit = (artificial_number.to_double() * 0.1).to_int()
      let excluded_artificial = scdfa_top_indices(
        scores, weak_artificial, artificial_limit, false,
      )
      for cell in excluded_artificial {
        included[cell] = false
      }
    }
  }
  let calls = Array::make(real_cells, false)
  for cell in 0..<real_cells {
    calls[cell] = scores[cell] >= threshold.threshold
  }
  let libraries = scdfa_library_sizes(local_data.counts)
  let detected = scdfa_detected_features(local_data.counts, 0.0)
  let above_two = scdfa_detected_features(local_data.counts, 2.0)
  let real_scores = []
  let real_weighted = []
  let real_ratios = []
  let real_cxds = []
  let real_nearest = []
  let real_doublet_distance = []
  let real_real_distance = []
  let real_nearest_class = []
  let real_origins = []
  let real_ambiguities = []
  let real_pca = []
  for cell in 0..<real_cells {
    real_scores.push(scores[cell])
    real_weighted.push(metrics.weighted[cell])
    real_ratios.push(metrics.ratios[cell])
    real_cxds.push(cxds[cell])
    real_nearest.push(metrics.nearest[cell])
    real_doublet_distance.push(metrics.nearest_doublet[cell])
    real_real_distance.push(metrics.nearest_real[cell])
    real_nearest_class.push(metrics.nearest_class[cell])
    real_origins.push(metrics.origins[cell])
    real_ambiguities.push(metrics.ambiguities[cell])
    real_pca.push(coordinates[cell].copy())
  }
  ScdfaCaptureFit::{
    global_indices: global_indices.copy(),
    scores: real_scores,
    calls,
    weighted_ratios: real_weighted,
    neighbor_ratios: real_ratios,
    cxds_scores: real_cxds,
    distances_to_nearest: real_nearest,
    distances_to_doublet: real_doublet_distance,
    distances_to_real: real_real_distance,
    nearest_classes: real_nearest_class,
    origins: real_origins,
    ambiguities: real_ambiguities,
    library_sizes: libraries,
    detected_features: detected,
    features_above_two: above_two,
    pca: real_pca,
    selected_indices: selected,
    threshold,
    expected_rate: heterotypic_rate,
    artificial_count: artificial_number,
  }
}

///|
fn scdfa_auto_clusters(
  data : SingleCellData,
  count : Int,
  config : ScDblFinderConfig,
) -> Array[String] {
  let selected = scdfa_selected_features(data.counts, config.selected_features)
  let coordinates = scdfa_processed_coordinates(
    data.counts,
    data.counts.length(),
    selected,
    config.dimensions,
  )
  let clusters = if count < data.counts.length() {
    count
  } else {
    data.counts.length()
  }
  let centroids = []
  centroids.push(coordinates[0].copy())
  while centroids.length() < clusters {
    let mut best_index = 0
    let mut best_distance = -1.0
    for cell in 0..<coordinates.length() {
      let mut nearest = scdfa_distance(coordinates[cell], centroids[0])
      for centroid in 1..<centroids.length() {
        let distance = scdfa_distance(coordinates[cell], centroids[centroid])
        if distance < nearest {
          nearest = distance
        }
      }
      if nearest > best_distance {
        best_distance = nearest
        best_index = cell
      }
    }
    centroids.push(coordinates[best_index].copy())
  }
  let assignments = Array::make(coordinates.length(), 0)
  for _ in 0..<50 {
    let mut changed = false
    for cell in 0..<coordinates.length() {
      let mut best = 0
      let mut best_distance = scdfa_distance(coordinates[cell], centroids[0])
      for cluster in 1..<clusters {
        let distance = scdfa_distance(coordinates[cell], centroids[cluster])
        if distance < best_distance {
          best = cluster
          best_distance = distance
        }
      }
      if assignments[cell] != best {
        assignments[cell] = best
        changed = true
      }
    }
    let sums = scdfa_zero_matrix(clusters, coordinates[0].length())
    let sizes = Array::make(clusters, 0)
    for cell in 0..<coordinates.length() {
      let cluster = assignments[cell]
      sizes[cluster] = sizes[cluster] + 1
      for component in 0..<coordinates[cell].length() {
        sums[cluster][component] = sums[cluster][component] +
          coordinates[cell][component]
      }
    }
    for cluster in 0..<clusters {
      if sizes[cluster] > 0 {
        for component in 0..<centroids[cluster].length() {
          centroids[cluster][component] = sums[cluster][component] /
            sizes[cluster].to_double()
        }
      }
    }
    if !changed {
      break
    }
  }
  let result = []
  for assignment in assignments {
    result.push("cluster_" + (assignment + 1).to_string())
  }
  result
}

///|
pub fn scdf_fast_cluster(
  data : SingleCellData,
  cluster_count : Int,
  dimensions? : Int = 20,
  selected_features? : Int = 1000,
) -> Array[String] raise ScDblFinderError {
  scdfa_validate_data(data)
  if cluster_count < 2 || cluster_count > data.counts.length() {
    scdfa_fail(
      "scDblFinder cluster_count must be between 2 and the number of cells",
    )
  }
  let config = ScDblFinderConfig::create(
    cluster_count~,
    dimensions~,
    selected_features~,
  )
  scdfa_auto_clusters(data, cluster_count, config)
}

///|
pub fn sc_dbl_finder(
  data : SingleCellData,
  clusters? : Array[String] = [],
  samples? : Array[String] = [],
  known_doublets? : Array[Bool] = [],
  config? : ScDblFinderConfig = ScDblFinderConfig::default(),
) -> ScDblFinderResult raise ScDblFinderError {
  scdfa_validate_data(data)
  let cells = data.counts.length()
  if clusters.length() > 0 && clusters.length() != cells {
    scdfa_fail("scDblFinder cluster length must match cells")
  }
  if samples.length() > 0 && samples.length() != cells {
    scdfa_fail("scDblFinder sample length must match cells")
  }
  if known_doublets.length() > 0 && known_doublets.length() != cells {
    scdfa_fail("scDblFinder known_doublets length must match cells")
  }
  for cluster in clusters {
    if cluster == "" {
      scdfa_fail("scDblFinder clusters cannot contain empty labels")
    }
  }
  let effective_samples = if samples.length() > 0 {
    samples.copy()
  } else {
    Array::make(cells, "all")
  }
  for sample in effective_samples {
    if sample == "" {
      scdfa_fail("scDblFinder samples cannot contain empty labels")
    }
  }
  let effective_clusters = if clusters.length() > 0 {
    clusters.copy()
  } else if config.cluster_count > 1 {
    scdfa_auto_clusters(data, config.cluster_count, config)
  } else {
    []
  }
  let effective_known = if known_doublets.length() > 0 {
    known_doublets.copy()
  } else {
    Array::make(cells, false)
  }
  let sample_ids = scdfa_unique_strings(effective_samples)
  let scores = Array::make(cells, 0.0)
  let calls = Array::make(cells, false)
  let weighted = Array::make(cells, 0.0)
  let ratios = Array::make(cells, 0.0)
  let cxds = Array::make(cells, 0.0)
  let nearest = Array::make(cells, 0.0)
  let nearest_doublet = Array::make(cells, 0.0)
  let nearest_real = Array::make(cells, 0.0)
  let nearest_classes = Array::make(cells, "")
  let origins = Array::make(cells, "")
  let ambiguities = Array::make(cells, 1.0)
  let libraries = Array::make(cells, 0.0)
  let detected = Array::make(cells, 0)
  let above_two = Array::make(cells, 0)
  let pca : Array[Array[Double]] = Array::make(cells, [])
  let thresholds = []
  let expected_rates = []
  let selected_indices : Array[Int] = []
  let mut artificial_total = 0
  for sample_id in sample_ids {
    let indices = []
    for cell in 0..<cells {
      if effective_samples[cell] == sample_id {
        indices.push(cell)
      }
    }
    if indices.length() < 3 {
      scdfa_fail(
        "scDblFinder sample " + sample_id + " has fewer than three cells",
      )
    }
    let fit = scdfa_fit_capture(
      data, indices, effective_clusters, effective_known, config,
    )
    thresholds.push(fit.threshold)
    expected_rates.push(fit.expected_rate)
    artificial_total = artificial_total + fit.artificial_count
    for index in fit.selected_indices {
      if !selected_indices.contains(index) {
        selected_indices.push(index)
      }
    }
    for local_index in 0..<fit.global_indices.length() {
      let global_index = fit.global_indices[local_index]
      scores[global_index] = fit.scores[local_index]
      calls[global_index] = fit.calls[local_index]
      weighted[global_index] = fit.weighted_ratios[local_index]
      ratios[global_index] = fit.neighbor_ratios[local_index]
      cxds[global_index] = fit.cxds_scores[local_index]
      nearest[global_index] = fit.distances_to_nearest[local_index]
      nearest_doublet[global_index] = fit.distances_to_doublet[local_index]
      nearest_real[global_index] = fit.distances_to_real[local_index]
      nearest_classes[global_index] = fit.nearest_classes[local_index]
      origins[global_index] = fit.origins[local_index]
      ambiguities[global_index] = fit.ambiguities[local_index]
      libraries[global_index] = fit.library_sizes[local_index]
      detected[global_index] = fit.detected_features[local_index]
      above_two[global_index] = fit.features_above_two[local_index]
      pca[global_index] = fit.pca[local_index].copy()
    }
  }
  let selected_gene_names = []
  for index in selected_indices {
    selected_gene_names.push(data.gene_names[index])
  }
  ScDblFinderResult::{
    cell_names: data.cell_names.copy(),
    gene_names: data.gene_names.copy(),
    samples: effective_samples,
    clusters: effective_clusters,
    known_doublets: effective_known,
    scores,
    calls,
    weighted_ratios: weighted,
    neighbor_ratios: ratios,
    cxds_scores: cxds,
    distances_to_nearest: nearest,
    distances_to_doublet: nearest_doublet,
    distances_to_real: nearest_real,
    nearest_classes,
    origins,
    origin_ambiguities: ambiguities,
    library_sizes: libraries,
    detected_features: detected,
    features_above_two: above_two,
    pca,
    selected_gene_names,
    sample_ids,
    thresholds,
    expected_doublet_rates: expected_rates,
    artificial_doublets: artificial_total,
    classifier_iterations: config.iterations,
    config,
  }
}

///|
pub fn ScDblFinderResult::n_cells(self : ScDblFinderResult) -> Int {
  self.cell_names.length()
}

///|
pub fn ScDblFinderResult::n_doublets(self : ScDblFinderResult) -> Int {
  let mut total = 0
  for call in self.calls {
    if call {
      total = total + 1
    }
  }
  total
}

///|
pub fn ScDblFinderResult::doublet_rate(self : ScDblFinderResult) -> Double {
  self.n_doublets().to_double() / self.n_cells().to_double()
}

///|
pub fn ScDblFinderResult::cell(
  self : ScDblFinderResult,
  cell_name : String,
) -> ScDblFinderCellResult? {
  let index = scdfa_string_index(self.cell_names, cell_name)
  if index < 0 {
    return None
  }
  Some(ScDblFinderCellResult::{
    cell_name,
    sample: self.samples[index],
    cluster: if self.clusters.length() > 0 {
      self.clusters[index]
    } else {
      ""
    },
    score: self.scores[index],
    is_doublet: self.calls[index],
    weighted_ratio: self.weighted_ratios[index],
    neighbor_ratio: self.neighbor_ratios[index],
    cxds_score: self.cxds_scores[index],
    nearest_class: self.nearest_classes[index],
    origin: self.origins[index],
    origin_ambiguity: self.origin_ambiguities[index],
  })
}

///|
pub fn ScDblFinderResult::top_doublets(
  self : ScDblFinderResult,
  limit? : Int = 10,
) -> Array[ScDblFinderCellResult] raise ScDblFinderError {
  if limit < 1 {
    scdfa_fail("scDblFinder top-doublet limit must be positive")
  }
  let candidates = []
  for index in 0..<self.n_cells() {
    candidates.push(index)
  }
  let ordered = scdfa_top_indices(self.scores, candidates, limit, true)
  let result = []
  for index in ordered {
    match self.cell(self.cell_names[index]) {
      Some(value) => result.push(value)
      None => ()
    }
  }
  result
}

///|
pub fn ScDblFinderResult::summary(self : ScDblFinderResult) -> String {
  "scDblFinder Result\n" +
  "cells=" +
  self.n_cells().to_string() +
  ", doublets=" +
  self.n_doublets().to_string() +
  ", rate=" +
  self.doublet_rate().to_string() +
  "\nsamples=" +
  self.sample_ids.length().to_string() +
  ", selected_genes=" +
  self.selected_gene_names.length().to_string() +
  ", artificial_doublets=" +
  self.artificial_doublets.to_string() +
  ", iterations=" +
  self.classifier_iterations.to_string()
}

///|
pub fn ScDblFinderResult::filter_singlets(
  self : ScDblFinderResult,
  data : SingleCellData,
) -> SingleCellData raise ScDblFinderError {
  scdfa_validate_data(data)
  if data.cell_names != self.cell_names || data.gene_names != self.gene_names {
    scdfa_fail("scDblFinder result and data identifiers do not match")
  }
  let counts = []
  let names = []
  for cell in 0..<self.n_cells() {
    if !self.calls[cell] {
      counts.push(data.counts[cell].copy())
      names.push(data.cell_names[cell])
    }
  }
  if counts.length() < 3 {
    scdfa_fail("scDblFinder filtering would leave fewer than three cells")
  }
  SingleCellData::create(counts, names, data.gene_names)
}

///|
fn scdfa_erf(value : Double) -> Double {
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
fn scdfa_normal_cdf(value : Double) -> Double {
  0.5 * (1.0 + scdfa_erf(value / 1.4142135623730951))
}

///|
fn scdfa_poisson_upper(observed : Int, expected : Double) -> Double {
  if observed <= 0 {
    return 1.0
  }
  if expected <= 0.0 {
    return 0.0
  }
  if expected > 100.0 {
    let z = (observed.to_double() - 0.5 - expected) / expected.sqrt()
    return 1.0 - scdfa_normal_cdf(z)
  }
  let mut term = @math.exp(-expected)
  let mut cumulative = term
  for count in 1..<observed {
    term = term * expected / count.to_double()
    cumulative = cumulative + term
  }
  let result = 1.0 - cumulative
  if result < 0.0 {
    0.0
  } else if result > 1.0 {
    1.0
  } else {
    result
  }
}

///|
fn scdfa_bh_adjust(values : Array[Double]) -> Array[Double] {
  let indices = []
  for index in 0..<values.length() {
    indices.push(index)
  }
  for index in 1..<indices.length() {
    let value = indices[index]
    let mut position = index
    while position > 0 && values[indices[position - 1]] > values[value] {
      indices[position] = indices[position - 1]
      position = position - 1
    }
    indices[position] = value
  }
  let adjusted = Array::make(values.length(), 1.0)
  let mut running = 1.0
  let mut rank = values.length()
  while rank > 0 {
    let index = indices[rank - 1]
    let candidate = values[index] *
      values.length().to_double() /
      rank.to_double()
    if candidate < running {
      running = candidate
    }
    adjusted[index] = if running > 1.0 { 1.0 } else { running }
    rank = rank - 1
  }
  adjusted
}

///|
pub fn scdf_pairwise_enrichment(
  result : ScDblFinderResult,
) -> Array[ScDblFinderEnrichment] raise ScDblFinderError {
  if result.clusters.length() == 0 {
    scdfa_fail("scDblFinder enrichment requires cluster-aware fitting")
  }
  let labels = scdfa_unique_strings(result.clusters)
  if labels.length() < 2 {
    scdfa_fail("scDblFinder enrichment requires at least two clusters")
  }
  let proportions = Array::make(labels.length(), 0.0)
  for cluster in result.clusters {
    let index = scdfa_string_index(labels, cluster)
    proportions[index] = proportions[index] + 1.0
  }
  for index in 0..<proportions.length() {
    proportions[index] = proportions[index] / result.n_cells().to_double()
  }
  let pair_names = []
  let pair_weights = []
  let observed = []
  for left in 0..<labels.length() {
    for right in (left + 1)..<labels.length() {
      pair_names.push(labels[left] + "+" + labels[right])
      pair_weights.push(2.0 * proportions[left] * proportions[right])
      observed.push(0)
    }
  }
  for cell in 0..<result.n_cells() {
    if result.calls[cell] {
      let index = scdfa_string_index(pair_names, result.origins[cell])
      if index >= 0 {
        observed[index] = observed[index] + 1
      }
    }
  }
  let heterotypic_probability = {
    let mut total = 0.0
    for weight in pair_weights {
      total = total + weight
    }
    total
  }
  let global_rate = if result.expected_doublet_rates.length() > 0 {
    let mut weighted = 0.0
    for sample in 0..<result.sample_ids.length() {
      let mut sample_cells = 0
      for value in result.samples {
        if value == result.sample_ids[sample] {
          sample_cells = sample_cells + 1
        }
      }
      weighted = weighted +
        result.expected_doublet_rates[sample] * sample_cells.to_double()
    }
    weighted / result.n_cells().to_double()
  } else {
    0.0
  }
  let p_values = []
  let expected_values = []
  for pair in 0..<pair_names.length() {
    let expected = if heterotypic_probability > 0.0 {
      global_rate *
      result.n_cells().to_double() *
      pair_weights[pair] /
      heterotypic_probability
    } else {
      0.0
    }
    expected_values.push(expected)
    p_values.push(scdfa_poisson_upper(observed[pair], expected))
  }
  let adjusted = scdfa_bh_adjust(p_values)
  let enrichment = []
  for pair in 0..<pair_names.length() {
    enrichment.push(ScDblFinderEnrichment::{
      combination: pair_names[pair],
      observed: observed[pair],
      expected: expected_values[pair],
      log2_enrichment: scdfa_log2(
        (observed[pair] + 1).to_double() / (expected_values[pair] + 1.0),
      ),
      p_value: p_values[pair],
      adjusted_p_value: adjusted[pair],
    })
  }
  for index in 1..<enrichment.length() {
    let value = enrichment[index]
    let mut position = index
    while position > 0 && enrichment[position - 1].p_value > value.p_value {
      enrichment[position] = enrichment[position - 1]
      position = position - 1
    }
    enrichment[position] = value
  }
  enrichment
}

///|
fn scdfa_copy_sce(source : SingleCellExperiment) -> SingleCellExperiment {
  let assays : Map[String, Array[Array[Double]]] = Map([])
  for key in source.assays.keys() {
    assays[key] = scdfa_copy_matrix(source.assays[key])
  }
  let row_data : Map[String, Array[String]] = Map([])
  for key in source.row_data.keys() {
    row_data[key] = source.row_data[key].copy()
  }
  let col_data : Map[String, Array[String]] = Map([])
  for key in source.col_data.keys() {
    col_data[key] = source.col_data[key].copy()
  }
  let reduced_dims : Map[String, Array[Array[Double]]] = Map([])
  for key in source.reduced_dims.keys() {
    reduced_dims[key] = scdfa_copy_matrix(source.reduced_dims[key])
  }
  let alternatives : Map[String, SingleCellExperiment] = Map([])
  for key in source.alternative_experiments.keys() {
    alternatives[key] = source.alternative_experiments[key]
  }
  let metadata : Map[String, String] = Map([])
  for key in source.metadata.keys() {
    metadata[key] = source.metadata[key]
  }
  SingleCellExperiment::{
    assays,
    row_data,
    col_data,
    reduced_dims,
    alternative_experiments: alternatives,
    metadata,
    row_names: source.row_names.copy(),
    col_names: source.col_names.copy(),
  }
}

///|
pub fn sc_dbl_finder_single_cell_experiment(
  experiment : SingleCellExperiment,
  assay_name? : String = "counts",
  cluster_column? : String = "",
  sample_column? : String = "",
  known_doublet_column? : String = "",
  output_prefix? : String = "scDblFinder",
  config? : ScDblFinderConfig = ScDblFinderConfig::default(),
) -> ScDblFinderSingleCellExperimentOutput raise ScDblFinderError {
  if assay_name == "" || output_prefix == "" {
    scdfa_fail("scDblFinder assay name and output prefix cannot be empty")
  }
  if !experiment.assays.contains(assay_name) {
    scdfa_fail("scDblFinder assay not found: " + assay_name)
  }
  let matrix = experiment.assays[assay_name]
  if matrix.length() == 0 || matrix[0].length() == 0 {
    scdfa_fail("scDblFinder assay cannot be empty")
  }
  let genes = matrix.length()
  let cells = matrix[0].length()
  if experiment.row_names.length() != genes ||
    experiment.col_names.length() != cells {
    scdfa_fail(
      "scDblFinder SingleCellExperiment names must match assay dimensions",
    )
  }
  let counts = scdfa_zero_matrix(cells, genes)
  for gene in 0..<genes {
    if matrix[gene].length() != cells {
      scdfa_fail("scDblFinder SingleCellExperiment assay must be rectangular")
    }
    for cell in 0..<cells {
      counts[cell][gene] = matrix[gene][cell]
    }
  }
  let data = SingleCellData::create(
    counts,
    experiment.col_names,
    experiment.row_names,
  )
  let clusters = if cluster_column == "" {
    []
  } else {
    if !experiment.col_data.contains(cluster_column) {
      scdfa_fail("scDblFinder cluster column not found: " + cluster_column)
    }
    experiment.col_data[cluster_column].copy()
  }
  let samples = if sample_column == "" {
    []
  } else {
    if !experiment.col_data.contains(sample_column) {
      scdfa_fail("scDblFinder sample column not found: " + sample_column)
    }
    experiment.col_data[sample_column].copy()
  }
  let known = if known_doublet_column == "" {
    []
  } else {
    if !experiment.col_data.contains(known_doublet_column) {
      scdfa_fail(
        "scDblFinder known-doublet column not found: " + known_doublet_column,
      )
    }
    let values = experiment.col_data[known_doublet_column]
    let parsed = []
    for value in values {
      if value == "true" || value == "doublet" || value == "1" {
        parsed.push(true)
      } else if value == "false" || value == "singlet" || value == "0" {
        parsed.push(false)
      } else {
        scdfa_fail(
          "scDblFinder known-doublet column contains invalid value: " + value,
        )
      }
    }
    parsed
  }
  let result = sc_dbl_finder(
    data,
    clusters~,
    samples~,
    known_doublets=known,
    config~,
  )
  let output = scdfa_copy_sce(experiment)
  let score_values = []
  let call_values = []
  let ratio_values = []
  let weighted_values = []
  let origin_values = []
  for cell in 0..<result.n_cells() {
    score_values.push(result.scores[cell].to_string())
    call_values.push(if result.calls[cell] { "doublet" } else { "singlet" })
    ratio_values.push(result.neighbor_ratios[cell].to_string())
    weighted_values.push(result.weighted_ratios[cell].to_string())
    origin_values.push(result.origins[cell])
  }
  output.col_data[output_prefix + ".score"] = score_values
  output.col_data[output_prefix + ".class"] = call_values
  output.col_data[output_prefix + ".ratio"] = ratio_values
  output.col_data[output_prefix + ".weighted"] = weighted_values
  output.col_data[output_prefix + ".mostLikelyOrigin"] = origin_values
  output.reduced_dims[output_prefix + ".PCA"] = scdfa_copy_matrix(result.pca)
  let selected = []
  for gene in output.row_names {
    selected.push(
      if scdfa_contains_string(result.selected_gene_names, gene) {
        "true"
      } else {
        "false"
      },
    )
  }
  output.row_data[output_prefix + ".selected"] = selected
  output.metadata[output_prefix + ".doublets"] = result.n_doublets().to_string()
  output.metadata[output_prefix + ".rate"] = result.doublet_rate().to_string()
  output.metadata[output_prefix + ".artificial"] = result.artificial_doublets.to_string()
  ScDblFinderSingleCellExperimentOutput::{ experiment: output, result, }
}

///|
pub fn scdf_create_advanced_example() -> (
  SingleCellData,
  Array[String],
  Array[String],
  Array[Bool],
) raise ScDblFinderError {
  let genes : Array[String] = []
  for gene in 0..<24 {
    genes.push("gene_" + (gene + 1).to_string())
  }
  let counts : Array[Array[Double]] = []
  let names : Array[String] = []
  let clusters : Array[String] = []
  let samples : Array[String] = []
  let truth : Array[Bool] = []
  for cluster in 0..<3 {
    for cell in 0..<12 {
      let row : Array[Double] = []
      for gene in 0..<24 {
        let block = gene / 8
        let base = if block == cluster { 14.0 } else { 1.0 }
        let variation = ((cell * 3 + gene * 5 + cluster) % 5).to_double() * 0.4
        row.push(base + variation)
      }
      counts.push(row)
      names.push("singlet_" + (cluster * 12 + cell + 1).to_string())
      clusters.push("type_" + (cluster + 1).to_string())
      samples.push(if cell % 2 == 0 { "capture_A" } else { "capture_B" })
      truth.push(false)
    }
  }
  let doublet_pairs = [(0, 1), (0, 2), (1, 2), (0, 1), (0, 2), (1, 2)]
  for doublet in 0..<doublet_pairs.length() {
    let pair = doublet_pairs[doublet]
    let row : Array[Double] = []
    for gene in 0..<24 {
      let first = if gene / 8 == pair.0 { 14.0 } else { 1.0 }
      let second = if gene / 8 == pair.1 { 14.0 } else { 1.0 }
      let variation = ((doublet * 7 + gene * 2) % 4).to_double() * 0.3
      row.push((first + second + variation) * 0.6)
    }
    counts.push(row)
    names.push("known_doublet_" + (doublet + 1).to_string())
    clusters.push("type_" + (pair.0 + 1).to_string())
    samples.push(if doublet % 2 == 0 { "capture_A" } else { "capture_B" })
    truth.push(true)
  }
  (SingleCellData::create(counts, names, genes), clusters, samples, truth)
}
