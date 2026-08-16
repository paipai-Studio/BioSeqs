// Portable implementation of the core Bioconductor SingleR 2.15.2 workflow.
//
// Expression matrices follow the Bioconductor gene x sample/cell convention.
// The compatibility API remains in single_r.mbt; this module provides strict
// training, iterative fine-tuning, pruning, cluster annotation, multi-reference
// recomputation, and SingleCellExperiment integration.

///|
pub suberror SingleRAdvancedError {
  SingleRAdvancedError(String)
}

///|
pub struct SingleRAdvancedConfig {
  quantile : Double
  fine_tune : Bool
  tune_threshold : Double
  prune : Bool
  nmads : Double
  minimum_delta_median : Double
  minimum_delta_next : Double
  marker_count : Int
  combine_marker_count : Int
} derive(Debug)

///|
pub fn SingleRAdvancedConfig::create(
  quantile? : Double = 0.8,
  fine_tune? : Bool = true,
  tune_threshold? : Double = 0.05,
  prune? : Bool = true,
  nmads? : Double = 3.0,
  minimum_delta_median? : Double = -1.0e300,
  minimum_delta_next? : Double = 0.0,
  marker_count? : Int = 0,
  combine_marker_count? : Int = 10,
) -> SingleRAdvancedConfig raise SingleRAdvancedError {
  if !single_r_adv_is_finite(quantile) || quantile < 0.0 || quantile > 1.0 {
    raise SingleRAdvancedError("quantile must be finite and in [0, 1]")
  }
  if !single_r_adv_is_finite(tune_threshold) || tune_threshold < 0.0 {
    raise SingleRAdvancedError("tune_threshold must be finite and non-negative")
  }
  if !single_r_adv_is_finite(nmads) || nmads < 0.0 {
    raise SingleRAdvancedError("nmads must be finite and non-negative")
  }
  if !single_r_adv_is_finite(minimum_delta_median) {
    raise SingleRAdvancedError("minimum_delta_median must be finite")
  }
  if !single_r_adv_is_finite(minimum_delta_next) {
    raise SingleRAdvancedError("minimum_delta_next must be finite")
  }
  if marker_count < 0 {
    raise SingleRAdvancedError("marker_count must be non-negative")
  }
  if combine_marker_count < 1 {
    raise SingleRAdvancedError("combine_marker_count must be positive")
  }
  SingleRAdvancedConfig::{
    quantile,
    fine_tune,
    tune_threshold,
    prune,
    nmads,
    minimum_delta_median,
    minimum_delta_next,
    marker_count,
    combine_marker_count,
  }
}

///|
pub fn SingleRAdvancedConfig::default() -> SingleRAdvancedConfig {
  SingleRAdvancedConfig::{
    quantile: 0.8,
    fine_tune: true,
    tune_threshold: 0.05,
    prune: true,
    nmads: 3.0,
    minimum_delta_median: -1.0e300,
    minimum_delta_next: 0.0,
    marker_count: 0,
    combine_marker_count: 10,
  }
}

///|
pub struct SingleRAdvancedReference {
  expression : Array[Array[Double]]
  labels : Array[String]
  gene_names : Array[String]
  sample_names : Array[String]
  label_names : Array[String]
  reference_name : String
  n_genes : Int
  n_samples : Int
} derive(Debug)

///|
pub struct SingleRAdvancedData {
  expression : Array[Array[Double]]
  gene_names : Array[String]
  cell_names : Array[String]
  n_genes : Int
  n_cells : Int
} derive(Debug)

///|
pub struct SingleRAdvancedTraining {
  reference : SingleRAdvancedReference
  test_gene_names : Array[String]
  common_gene_names : Array[String]
  reference_gene_indices : Array[Int]
  test_gene_indices : Array[Int]
  label_names : Array[String]
  sample_label_indices : Array[Int]
  pairwise_markers : Array[Array[Array[Int]]]
  common_markers : Array[Int]
  marker_count : Int
} derive(Debug)

///|
pub struct SingleRAdvancedResult {
  cell_names : Array[String]
  label_names : Array[String]
  scores : Array[Array[Double]]
  labels : Array[String]
  first_labels : Array[String]
  delta_next : Array[Double]
  delta_median : Array[Double]
  pruned_labels : Array[String?]
  fine_tuned : Array[Bool]
  pruning_thresholds : Array[Double?]
  common_gene_names : Array[String]
  marker_gene_names : Array[String]
  reference_name : String
  quantile : Double
} derive(Debug)

///|
pub struct SingleRAdvancedCombinedResult {
  cell_names : Array[String]
  labels : Array[String]
  pruned_labels : Array[String?]
  references : Array[String]
  selected_references : Array[String]
  scores : Array[Array[Double]]
  delta_next : Array[Double]
  per_reference : Array[SingleRAdvancedResult]
  marker_gene_names : Array[String]
} derive(Debug)

///|
pub struct SingleRAdvancedSCEOutput {
  experiment : SingleCellExperiment
  result : SingleRAdvancedResult
}

///|
fn single_r_adv_is_finite(value : Double) -> Bool {
  value == value && value.abs() <= 1.0e300
}

///|
fn single_r_adv_copy_matrix(
  matrix : Array[Array[Double]],
) -> Array[Array[Double]] {
  let output : Array[Array[Double]] = []
  for row in matrix {
    output.push(row.copy())
  }
  output
}

///|
fn single_r_adv_prepare_names(
  supplied : Array[String],
  count : Int,
  prefix : String,
) -> Array[String] raise SingleRAdvancedError {
  if supplied.length() != 0 && supplied.length() != count {
    raise SingleRAdvancedError(
      prefix + " names must match the matrix dimension",
    )
  }
  let output : Array[String] = []
  let seen : Map[String, Bool] = Map([])
  for index in 0..<count {
    let name = if supplied.length() == 0 {
      prefix + (index + 1).to_string()
    } else {
      supplied[index]
    }
    if name.trim().length() == 0 {
      raise SingleRAdvancedError(prefix + " names cannot be empty")
    }
    if seen.contains(name) {
      raise SingleRAdvancedError(prefix + " names must be unique: " + name)
    }
    seen[name] = true
    output.push(name)
  }
  output
}

///|
fn single_r_adv_validate_matrix(
  expression : Array[Array[Double]],
  column_name : String,
) -> (Int, Int) raise SingleRAdvancedError {
  if expression.length() == 0 {
    raise SingleRAdvancedError("expression must contain at least one gene")
  }
  let columns = expression[0].length()
  if columns == 0 {
    raise SingleRAdvancedError(
      "expression must contain at least one " + column_name,
    )
  }
  for row in expression {
    if row.length() != columns {
      raise SingleRAdvancedError("expression must be rectangular")
    }
    for value in row {
      if !single_r_adv_is_finite(value) {
        raise SingleRAdvancedError("expression values must be finite")
      }
    }
  }
  (expression.length(), columns)
}

///|
fn single_r_adv_unique_labels(labels : Array[String]) -> Array[String] {
  let output : Array[String] = []
  for label in labels {
    if !output.contains(label) {
      output.push(label)
    }
  }
  output
}

///|
pub fn SingleRAdvancedReference::create(
  expression : Array[Array[Double]],
  labels : Array[String],
  gene_names? : Array[String] = [],
  sample_names? : Array[String] = [],
  reference_name? : String = "reference",
) -> SingleRAdvancedReference raise SingleRAdvancedError {
  let (n_genes, n_samples) = single_r_adv_validate_matrix(
    expression, "reference sample",
  )
  if labels.length() != n_samples {
    raise SingleRAdvancedError(
      "labels must contain one value per reference sample",
    )
  }
  let copied_labels : Array[String] = []
  for label in labels {
    if label.trim().length() == 0 {
      raise SingleRAdvancedError("reference labels cannot be empty")
    }
    copied_labels.push(label)
  }
  if reference_name.trim().length() == 0 {
    raise SingleRAdvancedError("reference_name cannot be empty")
  }
  let prepared_genes = single_r_adv_prepare_names(gene_names, n_genes, "gene")
  let prepared_samples = single_r_adv_prepare_names(
    sample_names, n_samples, "reference sample",
  )
  SingleRAdvancedReference::{
    expression: single_r_adv_copy_matrix(expression),
    labels: copied_labels,
    gene_names: prepared_genes,
    sample_names: prepared_samples,
    label_names: single_r_adv_unique_labels(copied_labels),
    reference_name,
    n_genes,
    n_samples,
  }
}

///|
pub fn SingleRAdvancedData::create(
  expression : Array[Array[Double]],
  gene_names? : Array[String] = [],
  cell_names? : Array[String] = [],
) -> SingleRAdvancedData raise SingleRAdvancedError {
  let (n_genes, n_cells) = single_r_adv_validate_matrix(expression, "cell")
  SingleRAdvancedData::{
    expression: single_r_adv_copy_matrix(expression),
    gene_names: single_r_adv_prepare_names(gene_names, n_genes, "gene"),
    cell_names: single_r_adv_prepare_names(cell_names, n_cells, "cell"),
    n_genes,
    n_cells,
  }
}

///|
pub fn SingleRAdvancedReference::copy_expression(
  self : SingleRAdvancedReference,
) -> Array[Array[Double]] {
  single_r_adv_copy_matrix(self.expression)
}

///|
pub fn SingleRAdvancedData::copy_expression(
  self : SingleRAdvancedData,
) -> Array[Array[Double]] {
  single_r_adv_copy_matrix(self.expression)
}

///|
fn single_r_adv_label_index(labels : Array[String], query : String) -> Int {
  for index in 0..<labels.length() {
    if labels[index] == query {
      return index
    }
  }
  -1
}

///|
fn single_r_adv_sorted_values(values : Array[Double]) -> Array[Double] {
  let output = values.copy()
  for index in 1..<output.length() {
    let value = output[index]
    let mut position = index
    while position > 0 && output[position - 1] > value {
      output[position] = output[position - 1]
      position = position - 1
    }
    output[position] = value
  }
  output
}

///|
fn single_r_adv_quantile(
  values : Array[Double],
  probability : Double,
) -> Double {
  if values.length() == 0 {
    return 0.0
  }
  let sorted = single_r_adv_sorted_values(values)
  if sorted.length() == 1 {
    return sorted[0]
  }
  let position = probability * (sorted.length() - 1).to_double()
  let lower = position.to_int()
  let upper = if lower + 1 < sorted.length() { lower + 1 } else { lower }
  let fraction = position - lower.to_double()
  sorted[lower] * (1.0 - fraction) + sorted[upper] * fraction
}

///|
fn single_r_adv_median(values : Array[Double]) -> Double {
  single_r_adv_quantile(values, 0.5)
}

///|
fn single_r_adv_rank(values : Array[Double]) -> Array[Double] {
  let ordered : Array[(Double, Int)] = []
  for index in 0..<values.length() {
    ordered.push((values[index], index))
  }
  for index in 1..<ordered.length() {
    let value = ordered[index]
    let mut position = index
    while position > 0 &&
          (
            ordered[position - 1].0 > value.0 ||
            (
              ordered[position - 1].0 == value.0 &&
              ordered[position - 1].1 > value.1
            )
          ) {
      ordered[position] = ordered[position - 1]
      position = position - 1
    }
    ordered[position] = value
  }
  let ranks = Array::make(values.length(), 0.0)
  let mut start = 0
  while start < ordered.length() {
    let mut end = start + 1
    while end < ordered.length() && ordered[end].0 == ordered[start].0 {
      end = end + 1
    }
    let average = (start + 1 + end).to_double() / 2.0
    for position in start..<end {
      ranks[ordered[position].1] = average
    }
    start = end
  }
  ranks
}

///|
pub fn single_r_advanced_spearman(
  left : Array[Double],
  right : Array[Double],
) -> Double raise SingleRAdvancedError {
  if left.length() != right.length() {
    raise SingleRAdvancedError(
      "Spearman correlation vectors must have equal length",
    )
  }
  if left.length() < 2 {
    return 0.0
  }
  for value in left {
    if !single_r_adv_is_finite(value) {
      raise SingleRAdvancedError("Spearman correlation vectors must be finite")
    }
  }
  for value in right {
    if !single_r_adv_is_finite(value) {
      raise SingleRAdvancedError("Spearman correlation vectors must be finite")
    }
  }
  let left_rank = single_r_adv_rank(left)
  let right_rank = single_r_adv_rank(right)
  let mean = (left.length() + 1).to_double() / 2.0
  let mut cross = 0.0
  let mut left_square = 0.0
  let mut right_square = 0.0
  for index in 0..<left.length() {
    let left_delta = left_rank[index] - mean
    let right_delta = right_rank[index] - mean
    cross = cross + left_delta * right_delta
    left_square = left_square + left_delta * left_delta
    right_square = right_square + right_delta * right_delta
  }
  if left_square <= 0.0 || right_square <= 0.0 {
    0.0
  } else {
    cross / (left_square * right_square).sqrt()
  }
}

///|
fn single_r_adv_default_marker_count(label_count : Int) -> Int {
  if label_count <= 1 {
    return 500
  }
  let exponent = @math.ln(label_count.to_double()) / @math.ln(2.0)
  let value = 500.0 * @math.pow(2.0 / 3.0, exponent)
  let converted = value.to_int()
  if converted < 1 {
    1
  } else {
    converted
  }
}

///|
fn single_r_adv_build_gene_index(names : Array[String]) -> Map[String, Int] {
  let output : Map[String, Int] = Map([])
  for index in 0..<names.length() {
    output[names[index]] = index
  }
  output
}

///|
fn single_r_adv_common_genes(
  reference : SingleRAdvancedReference,
  test_gene_names : Array[String],
  restrict_genes : Array[String],
) -> (Array[String], Array[Int], Array[Int]) raise SingleRAdvancedError {
  let test_names = single_r_adv_prepare_names(
    test_gene_names,
    test_gene_names.length(),
    "test gene",
  )
  if test_names.length() == 0 {
    raise SingleRAdvancedError("test_gene_names cannot be empty")
  }
  let reference_index = single_r_adv_build_gene_index(reference.gene_names)
  let restriction : Map[String, Bool] = Map([])
  for name in restrict_genes {
    if name.trim().length() == 0 {
      raise SingleRAdvancedError("restrict_genes cannot contain empty names")
    }
    restriction[name] = true
  }
  let common : Array[String] = []
  let reference_indices : Array[Int] = []
  let test_indices : Array[Int] = []
  for test_index in 0..<test_names.length() {
    let name = test_names[test_index]
    if (restrict_genes.length() == 0 || restriction.contains(name)) &&
      reference_index.contains(name) {
      common.push(name)
      reference_indices.push(reference_index[name])
      test_indices.push(test_index)
    }
  }
  if common.length() < 2 {
    raise SingleRAdvancedError(
      "reference and test must share at least two eligible genes",
    )
  }
  (common, reference_indices, test_indices)
}

///|
fn single_r_adv_label_medians(
  reference : SingleRAdvancedReference,
  reference_gene_indices : Array[Int],
  label_names : Array[String],
  sample_label_indices : Array[Int],
) -> Array[Array[Double]] {
  let output : Array[Array[Double]] = []
  for label in 0..<label_names.length() {
    let medians : Array[Double] = []
    for common_gene in 0..<reference_gene_indices.length() {
      let values : Array[Double] = []
      let reference_gene = reference_gene_indices[common_gene]
      for sample in 0..<reference.n_samples {
        if sample_label_indices[sample] == label {
          values.push(reference.expression[reference_gene][sample])
        }
      }
      medians.push(single_r_adv_median(values))
    }
    output.push(medians)
  }
  output
}

///|
fn single_r_adv_sort_gene_effects(
  effects : Array[(Int, Double)],
) -> Array[(Int, Double)] {
  let output = effects.copy()
  for index in 1..<output.length() {
    let value = output[index]
    let mut position = index
    while position > 0 &&
          (
            output[position - 1].1 < value.1 ||
            (
              output[position - 1].1 == value.1 &&
              output[position - 1].0 > value.0
            )
          ) {
      output[position] = output[position - 1]
      position = position - 1
    }
    output[position] = value
  }
  output
}

///|
fn single_r_adv_select_markers(
  medians : Array[Array[Double]],
  marker_count : Int,
) -> Array[Array[Array[Int]]] {
  let label_count = medians.length()
  let gene_count = if label_count == 0 { 0 } else { medians[0].length() }
  let pairwise : Array[Array[Array[Int]]] = []
  for target in 0..<label_count {
    let target_markers : Array[Array[Int]] = []
    for other in 0..<label_count {
      let selected : Array[Int] = []
      if target != other {
        let effects : Array[(Int, Double)] = []
        for gene in 0..<gene_count {
          let effect = medians[target][gene] - medians[other][gene]
          if effect > 0.0 {
            effects.push((gene, effect))
          }
        }
        let sorted = single_r_adv_sort_gene_effects(effects)
        let limit = if marker_count < sorted.length() {
          marker_count
        } else {
          sorted.length()
        }
        for index in 0..<limit {
          selected.push(sorted[index].0)
        }
      }
      target_markers.push(selected)
    }
    pairwise.push(target_markers)
  }
  pairwise
}

///|
fn single_r_adv_union_markers(
  pairwise : Array[Array[Array[Int]]],
  gene_count : Int,
) -> Array[Int] {
  let selected = Array::make(gene_count, false)
  for target in pairwise {
    for comparison in target {
      for gene in comparison {
        selected[gene] = true
      }
    }
  }
  let output : Array[Int] = []
  for gene in 0..<gene_count {
    if selected[gene] {
      output.push(gene)
    }
  }
  if output.length() < 2 {
    for gene in 0..<gene_count {
      if !output.contains(gene) {
        output.push(gene)
      }
    }
  }
  output
}

///|
pub fn single_r_advanced_train(
  reference : SingleRAdvancedReference,
  test_gene_names : Array[String],
  restrict_genes? : Array[String] = [],
  marker_count? : Int = 0,
) -> SingleRAdvancedTraining raise SingleRAdvancedError {
  if marker_count < 0 {
    raise SingleRAdvancedError("marker_count must be non-negative")
  }
  let (common, reference_indices, test_indices) = single_r_adv_common_genes(
    reference, test_gene_names, restrict_genes,
  )
  let label_names = reference.label_names.copy()
  let sample_label_indices : Array[Int] = []
  for label in reference.labels {
    let index = single_r_adv_label_index(label_names, label)
    if index < 0 {
      raise SingleRAdvancedError("internal reference label mismatch")
    }
    sample_label_indices.push(index)
  }
  let selected_count = if marker_count == 0 {
    single_r_adv_default_marker_count(label_names.length())
  } else {
    marker_count
  }
  let medians = single_r_adv_label_medians(
    reference, reference_indices, label_names, sample_label_indices,
  )
  let pairwise = single_r_adv_select_markers(medians, selected_count)
  let common_markers = single_r_adv_union_markers(pairwise, common.length())
  SingleRAdvancedTraining::{
    reference,
    test_gene_names: test_gene_names.copy(),
    common_gene_names: common,
    reference_gene_indices: reference_indices,
    test_gene_indices: test_indices,
    label_names,
    sample_label_indices,
    pairwise_markers: pairwise,
    common_markers,
    marker_count: selected_count,
  }
}

///|
pub fn SingleRAdvancedTraining::marker_gene_names(
  self : SingleRAdvancedTraining,
) -> Array[String] {
  let output : Array[String] = []
  for index in self.common_markers {
    output.push(self.common_gene_names[index])
  }
  output
}

///|
pub fn SingleRAdvancedTraining::markers_between(
  self : SingleRAdvancedTraining,
  target : String,
  other : String,
) -> Array[String] raise SingleRAdvancedError {
  let target_index = single_r_adv_label_index(self.label_names, target)
  let other_index = single_r_adv_label_index(self.label_names, other)
  if target_index < 0 || other_index < 0 {
    raise SingleRAdvancedError("unknown label in marker query")
  }
  let output : Array[String] = []
  for index in self.pairwise_markers[target_index][other_index] {
    output.push(self.common_gene_names[index])
  }
  output
}

///|
fn single_r_adv_same_strings(
  left : Array[String],
  right : Array[String],
) -> Bool {
  if left.length() != right.length() {
    return false
  }
  for index in 0..<left.length() {
    if left[index] != right[index] {
      return false
    }
  }
  true
}

///|
fn single_r_adv_correlation(
  data : SingleRAdvancedData,
  cell : Int,
  training : SingleRAdvancedTraining,
  sample : Int,
  features : Array[Int],
) -> Double raise SingleRAdvancedError {
  let test_values : Array[Double] = []
  let reference_values : Array[Double] = []
  for feature in features {
    test_values.push(data.expression[training.test_gene_indices[feature]][cell])
    reference_values.push(
      training.reference.expression[training.reference_gene_indices[feature]][sample],
    )
  }
  single_r_advanced_spearman(test_values, reference_values)
}

///|
fn single_r_adv_score_labels(
  data : SingleRAdvancedData,
  cell : Int,
  training : SingleRAdvancedTraining,
  features : Array[Int],
  labels : Array[Int],
  quantile : Double,
) -> Array[Double] raise SingleRAdvancedError {
  let output : Array[Double] = []
  for label in labels {
    let correlations : Array[Double] = []
    for sample in 0..<training.reference.n_samples {
      if training.sample_label_indices[sample] == label {
        correlations.push(
          single_r_adv_correlation(data, cell, training, sample, features),
        )
      }
    }
    output.push(single_r_adv_quantile(correlations, quantile))
  }
  output
}

///|
fn single_r_adv_best_score(
  labels : Array[Int],
  scores : Array[Double],
) -> (Int, Double, Double) {
  if labels.length() == 0 {
    return (-1, 0.0, 0.0)
  }
  let mut best_position = 0
  let mut second = -2.0
  for position in 1..<scores.length() {
    if scores[position] > scores[best_position] {
      second = scores[best_position].max(second)
      best_position = position
    } else if scores[position] > second {
      second = scores[position]
    }
  }
  if scores.length() == 1 {
    second = scores[best_position]
  }
  (labels[best_position], scores[best_position], scores[best_position] - second)
}

///|
fn single_r_adv_candidates(
  labels : Array[Int],
  scores : Array[Double],
  threshold : Double,
) -> Array[Int] {
  if labels.length() == 0 {
    return []
  }
  let mut maximum = scores[0]
  for score in scores {
    if score > maximum {
      maximum = score
    }
  }
  let output : Array[Int] = []
  for position in 0..<labels.length() {
    if scores[position] >= maximum - threshold {
      output.push(labels[position])
    }
  }
  output
}

///|
fn single_r_adv_candidate_markers(
  training : SingleRAdvancedTraining,
  candidates : Array[Int],
) -> Array[Int] {
  let chosen = Array::make(training.common_gene_names.length(), false)
  for target in candidates {
    for other in candidates {
      if target != other {
        for gene in training.pairwise_markers[target][other] {
          chosen[gene] = true
        }
      }
    }
  }
  let output : Array[Int] = []
  for gene in 0..<chosen.length() {
    if chosen[gene] {
      output.push(gene)
    }
  }
  if output.length() < 2 {
    training.common_markers.copy()
  } else {
    output
  }
}

///|
fn single_r_adv_same_ints(left : Array[Int], right : Array[Int]) -> Bool {
  if left.length() != right.length() {
    return false
  }
  for index in 0..<left.length() {
    if left[index] != right[index] {
      return false
    }
  }
  true
}

///|
fn single_r_adv_fine_tune(
  data : SingleRAdvancedData,
  cell : Int,
  training : SingleRAdvancedTraining,
  initial_scores : Array[Double],
  config : SingleRAdvancedConfig,
) -> (Int, Double, Bool) raise SingleRAdvancedError {
  let all_labels : Array[Int] = []
  for label in 0..<training.label_names.length() {
    all_labels.push(label)
  }
  let (initial_best, _, initial_delta) = single_r_adv_best_score(
    all_labels, initial_scores,
  )
  if !config.fine_tune || all_labels.length() < 2 {
    return (initial_best, initial_delta, false)
  }
  let mut candidates = single_r_adv_candidates(
    all_labels,
    initial_scores,
    config.tune_threshold,
  )
  if candidates.length() < 2 {
    return (initial_best, initial_delta, false)
  }
  let mut best = initial_best
  let mut delta = initial_delta
  let mut iteration = 0
  while candidates.length() > 1 && iteration <= training.label_names.length() {
    let features = single_r_adv_candidate_markers(training, candidates)
    let scores = single_r_adv_score_labels(
      data,
      cell,
      training,
      features,
      candidates,
      config.quantile,
    )
    let (current_best, _, current_delta) = single_r_adv_best_score(
      candidates, scores,
    )
    best = current_best
    delta = current_delta
    let next = single_r_adv_candidates(
      candidates,
      scores,
      config.tune_threshold,
    )
    if single_r_adv_same_ints(candidates, next) {
      break
    }
    candidates = next
    iteration = iteration + 1
  }
  (best, delta, true)
}

///|
fn single_r_adv_row_median(values : Array[Double]) -> Double {
  single_r_adv_median(values)
}

///|
fn single_r_adv_mad(values : Array[Double]) -> Double {
  if values.length() == 0 {
    return 0.0
  }
  let center = single_r_adv_median(values)
  let deviations : Array[Double] = []
  for value in values {
    deviations.push((value - center).abs())
  }
  1.482602218505602 * single_r_adv_median(deviations)
}

///|
fn single_r_adv_apply_pruning(
  labels : Array[String],
  label_names : Array[String],
  delta_median : Array[Double],
  delta_next : Array[Double],
  config : SingleRAdvancedConfig,
) -> (Array[String?], Array[Double?]) {
  let output : Array[String?] = []
  for label in labels {
    output.push(Some(label))
  }
  let thresholds : Array[Double?] = Array::make(label_names.length(), None)
  if !config.prune {
    return (output, thresholds)
  }
  let keep = Array::make(labels.length(), true)
  for cell in 0..<labels.length() {
    if delta_median[cell] < config.minimum_delta_median ||
      delta_next[cell] < config.minimum_delta_next {
      keep[cell] = false
    }
  }
  for label in 0..<label_names.length() {
    let values : Array[Double] = []
    for cell in 0..<labels.length() {
      if keep[cell] && labels[cell] == label_names[label] {
        values.push(delta_median[cell])
      }
    }
    if values.length() > 0 {
      let threshold = single_r_adv_median(values) -
        config.nmads * single_r_adv_mad(values)
      thresholds[label] = Some(threshold)
      for cell in 0..<labels.length() {
        if keep[cell] &&
          labels[cell] == label_names[label] &&
          delta_median[cell] < threshold {
          keep[cell] = false
        }
      }
    }
  }
  for cell in 0..<labels.length() {
    if !keep[cell] {
      output[cell] = None
    }
  }
  (output, thresholds)
}

///|
pub fn single_r_advanced_classify(
  data : SingleRAdvancedData,
  training : SingleRAdvancedTraining,
  config? : SingleRAdvancedConfig = SingleRAdvancedConfig::default(),
) -> SingleRAdvancedResult raise SingleRAdvancedError {
  if !single_r_adv_same_strings(data.gene_names, training.test_gene_names) {
    raise SingleRAdvancedError(
      "test gene names must exactly match those used during training",
    )
  }
  let all_labels : Array[Int] = []
  for label in 0..<training.label_names.length() {
    all_labels.push(label)
  }
  let scores : Array[Array[Double]] = []
  let labels : Array[String] = []
  let first_labels : Array[String] = []
  let delta_next : Array[Double] = []
  let delta_median : Array[Double] = []
  let fine_tuned : Array[Bool] = []
  for cell in 0..<data.n_cells {
    let initial = single_r_adv_score_labels(
      data,
      cell,
      training,
      training.common_markers,
      all_labels,
      config.quantile,
    )
    scores.push(initial)
    let (first, _, _) = single_r_adv_best_score(all_labels, initial)
    first_labels.push(training.label_names[first])
    let (best, next_delta, tuned) = single_r_adv_fine_tune(
      data, cell, training, initial, config,
    )
    labels.push(training.label_names[best])
    delta_next.push(next_delta)
    fine_tuned.push(tuned)
    delta_median.push(initial[best] - single_r_adv_row_median(initial))
  }
  let (pruned_labels, pruning_thresholds) = single_r_adv_apply_pruning(
    labels,
    training.label_names,
    delta_median,
    delta_next,
    config,
  )
  SingleRAdvancedResult::{
    cell_names: data.cell_names.copy(),
    label_names: training.label_names.copy(),
    scores,
    labels,
    first_labels,
    delta_next,
    delta_median,
    pruned_labels,
    fine_tuned,
    pruning_thresholds,
    common_gene_names: training.common_gene_names.copy(),
    marker_gene_names: training.marker_gene_names(),
    reference_name: training.reference.reference_name,
    quantile: config.quantile,
  }
}

///|
pub fn single_r_advanced(
  data : SingleRAdvancedData,
  reference : SingleRAdvancedReference,
  restrict_genes? : Array[String] = [],
  config? : SingleRAdvancedConfig = SingleRAdvancedConfig::default(),
) -> SingleRAdvancedResult raise SingleRAdvancedError {
  let training = single_r_advanced_train(
    reference,
    data.gene_names,
    restrict_genes~,
    marker_count=config.marker_count,
  )
  single_r_advanced_classify(data, training, config~)
}

///|
pub fn SingleRAdvancedResult::n_cells(self : SingleRAdvancedResult) -> Int {
  self.cell_names.length()
}

///|
pub fn SingleRAdvancedResult::n_pruned(self : SingleRAdvancedResult) -> Int {
  let mut count = 0
  for value in self.pruned_labels {
    if value is None {
      count = count + 1
    }
  }
  count
}

///|
pub fn SingleRAdvancedResult::label_index(
  self : SingleRAdvancedResult,
  label : String,
) -> Int {
  single_r_adv_label_index(self.label_names, label)
}

///|
pub fn SingleRAdvancedResult::assigned_score(
  self : SingleRAdvancedResult,
  cell : Int,
) -> Double raise SingleRAdvancedError {
  if cell < 0 || cell >= self.n_cells() {
    raise SingleRAdvancedError("cell index is out of range")
  }
  let label = self.label_index(self.labels[cell])
  if label < 0 {
    raise SingleRAdvancedError("assigned label is missing from score columns")
  }
  self.scores[cell][label]
}

///|
pub fn SingleRAdvancedResult::summary(
  self : SingleRAdvancedResult,
) -> Map[String, Int] {
  let output : Map[String, Int] = Map([])
  for name in self.label_names {
    output[name] = 0
  }
  for label in self.labels {
    output[label] = output.get(label).unwrap_or(0) + 1
  }
  output["pruned"] = self.n_pruned()
  output["total"] = self.n_cells()
  output
}

///|
fn single_r_adv_aggregate_clusters(
  data : SingleRAdvancedData,
  clusters : Array[String],
) -> SingleRAdvancedData raise SingleRAdvancedError {
  if clusters.length() != data.n_cells {
    raise SingleRAdvancedError("clusters must contain one value per cell")
  }
  let cluster_names : Array[String] = []
  for cluster in clusters {
    if cluster.trim().length() == 0 {
      raise SingleRAdvancedError("cluster names cannot be empty")
    }
    if !cluster_names.contains(cluster) {
      cluster_names.push(cluster)
    }
  }
  let expression : Array[Array[Double]] = []
  for gene in 0..<data.n_genes {
    let row = Array::make(cluster_names.length(), 0.0)
    for cell in 0..<data.n_cells {
      let cluster = single_r_adv_label_index(cluster_names, clusters[cell])
      row[cluster] = row[cluster] + data.expression[gene][cell]
    }
    expression.push(row)
  }
  SingleRAdvancedData::create(
    expression,
    gene_names=data.gene_names,
    cell_names=cluster_names,
  )
}

///|
pub fn single_r_advanced_clusters(
  data : SingleRAdvancedData,
  clusters : Array[String],
  reference : SingleRAdvancedReference,
  restrict_genes? : Array[String] = [],
  config? : SingleRAdvancedConfig = SingleRAdvancedConfig::default(),
) -> SingleRAdvancedResult raise SingleRAdvancedError {
  let aggregated = single_r_adv_aggregate_clusters(data, clusters)
  single_r_advanced(aggregated, reference, restrict_genes~, config~)
}

///|
fn single_r_adv_marker_names_for_label(
  training : SingleRAdvancedTraining,
  label : Int,
  top : Int,
) -> Array[String] {
  let output : Array[String] = []
  for other in 0..<training.label_names.length() {
    let markers = training.pairwise_markers[label][other]
    let limit = if top < markers.length() { top } else { markers.length() }
    for position in 0..<limit {
      let name = training.common_gene_names[markers[position]]
      if !output.contains(name) {
        output.push(name)
      }
    }
  }
  output
}

///|
fn single_r_adv_global_common_genes(
  data : SingleRAdvancedData,
  references : Array[SingleRAdvancedReference],
) -> Array[String] {
  let reference_maps : Array[Map[String, Int]] = []
  for reference in references {
    reference_maps.push(single_r_adv_build_gene_index(reference.gene_names))
  }
  let output : Array[String] = []
  for name in data.gene_names {
    let mut present = true
    for index in reference_maps {
      if !index.contains(name) {
        present = false
      }
    }
    if present {
      output.push(name)
    }
  }
  output
}

///|
fn single_r_adv_named_correlation(
  data : SingleRAdvancedData,
  cell : Int,
  reference : SingleRAdvancedReference,
  sample : Int,
  genes : Array[String],
) -> Double raise SingleRAdvancedError {
  let data_index = single_r_adv_build_gene_index(data.gene_names)
  let reference_index = single_r_adv_build_gene_index(reference.gene_names)
  let left : Array[Double] = []
  let right : Array[Double] = []
  for gene in genes {
    if data_index.contains(gene) && reference_index.contains(gene) {
      left.push(data.expression[data_index[gene]][cell])
      right.push(reference.expression[reference_index[gene]][sample])
    }
  }
  single_r_advanced_spearman(left, right)
}

///|
fn single_r_adv_recomputed_score(
  data : SingleRAdvancedData,
  cell : Int,
  reference : SingleRAdvancedReference,
  label : String,
  genes : Array[String],
  quantile : Double,
) -> Double raise SingleRAdvancedError {
  let correlations : Array[Double] = []
  for sample in 0..<reference.n_samples {
    if reference.labels[sample] == label {
      correlations.push(
        single_r_adv_named_correlation(data, cell, reference, sample, genes),
      )
    }
  }
  single_r_adv_quantile(correlations, quantile)
}

///|
pub fn single_r_advanced_combine(
  data : SingleRAdvancedData,
  references : Array[SingleRAdvancedReference],
  restrict_genes? : Array[String] = [],
  config? : SingleRAdvancedConfig = SingleRAdvancedConfig::default(),
) -> SingleRAdvancedCombinedResult raise SingleRAdvancedError {
  if references.length() < 2 {
    raise SingleRAdvancedError(
      "multi-reference integration requires at least two references",
    )
  }
  let reference_names : Array[String] = []
  for reference in references {
    if reference_names.contains(reference.reference_name) {
      raise SingleRAdvancedError("reference names must be unique")
    }
    reference_names.push(reference.reference_name)
  }
  let global_common = single_r_adv_global_common_genes(data, references)
  if global_common.length() < 2 {
    raise SingleRAdvancedError(
      "all references and test data must share at least two genes",
    )
  }
  let effective_restriction : Array[String] = []
  if restrict_genes.length() == 0 {
    for name in global_common {
      effective_restriction.push(name)
    }
  } else {
    for name in global_common {
      if restrict_genes.contains(name) {
        effective_restriction.push(name)
      }
    }
  }
  if effective_restriction.length() < 2 {
    raise SingleRAdvancedError(
      "gene restriction leaves fewer than two globally shared genes",
    )
  }
  let trainings : Array[SingleRAdvancedTraining] = []
  let per_reference : Array[SingleRAdvancedResult] = []
  for reference in references {
    let training = single_r_advanced_train(
      reference,
      data.gene_names,
      restrict_genes=effective_restriction,
      marker_count=config.marker_count,
    )
    trainings.push(training)
    per_reference.push(single_r_advanced_classify(data, training, config~))
  }
  let combined_scores : Array[Array[Double]] = []
  let labels : Array[String] = []
  let pruned_labels : Array[String?] = []
  let selected_references : Array[String] = []
  let delta_next : Array[Double] = []
  let all_markers : Array[String] = []
  for cell in 0..<data.n_cells {
    let marker_names : Array[String] = []
    for reference in 0..<references.length() {
      let label = single_r_adv_label_index(
        trainings[reference].label_names,
        per_reference[reference].labels[cell],
      )
      let current = single_r_adv_marker_names_for_label(
        trainings[reference],
        label,
        config.combine_marker_count,
      )
      for name in current {
        if effective_restriction.contains(name) && !marker_names.contains(name) {
          marker_names.push(name)
        }
      }
    }
    if marker_names.length() < 2 {
      for name in effective_restriction {
        if !marker_names.contains(name) {
          marker_names.push(name)
        }
      }
    }
    for name in marker_names {
      if !all_markers.contains(name) {
        all_markers.push(name)
      }
    }
    let current_scores : Array[Double] = []
    for reference in 0..<references.length() {
      current_scores.push(
        single_r_adv_recomputed_score(
          data,
          cell,
          references[reference],
          per_reference[reference].labels[cell],
          marker_names,
          config.quantile,
        ),
      )
    }
    combined_scores.push(current_scores)
    let reference_indices : Array[Int] = []
    for index in 0..<references.length() {
      reference_indices.push(index)
    }
    let (best, _, delta) = single_r_adv_best_score(
      reference_indices, current_scores,
    )
    labels.push(per_reference[best].labels[cell])
    pruned_labels.push(per_reference[best].pruned_labels[cell])
    selected_references.push(references[best].reference_name)
    delta_next.push(delta)
  }
  SingleRAdvancedCombinedResult::{
    cell_names: data.cell_names.copy(),
    labels,
    pruned_labels,
    references: reference_names,
    selected_references,
    scores: combined_scores,
    delta_next,
    per_reference,
    marker_gene_names: all_markers,
  }
}

///|
pub fn SingleRAdvancedCombinedResult::n_cells(
  self : SingleRAdvancedCombinedResult,
) -> Int {
  self.cell_names.length()
}

///|
pub fn SingleRAdvancedCombinedResult::reference_summary(
  self : SingleRAdvancedCombinedResult,
) -> Map[String, Int] {
  let output : Map[String, Int] = Map([])
  for reference in self.references {
    output[reference] = 0
  }
  for reference in self.selected_references {
    output[reference] = output.get(reference).unwrap_or(0) + 1
  }
  output
}

///|
fn single_r_adv_copy_sce(source : SingleCellExperiment) -> SingleCellExperiment {
  let assays : Map[String, Array[Array[Double]]] = Map([])
  for key in source.assays.keys() {
    assays[key] = single_r_adv_copy_matrix(source.assays[key])
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
    reduced_dims[key] = single_r_adv_copy_matrix(source.reduced_dims[key])
  }
  let alternatives : Map[String, SingleCellExperiment] = Map([])
  for key in source.alternative_experiments.keys() {
    alternatives[key] = single_r_adv_copy_sce(
      source.alternative_experiments[key],
    )
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
fn single_r_adv_expand_cluster_strings(
  clusters : Array[String],
  cluster_names : Array[String],
  values : Array[String],
) -> Array[String] {
  let output : Array[String] = []
  for cluster in clusters {
    let index = single_r_adv_label_index(cluster_names, cluster)
    output.push(values[index])
  }
  output
}

///|
fn single_r_adv_expand_cluster_doubles(
  clusters : Array[String],
  cluster_names : Array[String],
  values : Array[Double],
) -> Array[Double] {
  let output : Array[Double] = []
  for cluster in clusters {
    let index = single_r_adv_label_index(cluster_names, cluster)
    output.push(values[index])
  }
  output
}

///|
fn single_r_adv_expand_cluster_options(
  clusters : Array[String],
  cluster_names : Array[String],
  values : Array[String?],
) -> Array[String?] {
  let output : Array[String?] = []
  for cluster in clusters {
    let index = single_r_adv_label_index(cluster_names, cluster)
    output.push(values[index])
  }
  output
}

///|
fn single_r_adv_option_string(value : String?) -> String {
  match value {
    Some(label) => label
    None => "NA"
  }
}

///|
pub fn single_r_advanced_sce(
  experiment : SingleCellExperiment,
  reference : SingleRAdvancedReference,
  assay_name? : String = "logcounts",
  cluster_column? : String = "",
  output_prefix? : String = "SingleR",
  restrict_genes? : Array[String] = [],
  config? : SingleRAdvancedConfig = SingleRAdvancedConfig::default(),
) -> SingleRAdvancedSCEOutput raise SingleRAdvancedError {
  if assay_name.trim().length() == 0 || output_prefix.trim().length() == 0 {
    raise SingleRAdvancedError("assay_name and output_prefix cannot be empty")
  }
  let assay = match experiment.assays.get(assay_name) {
    Some(value) => value
    None =>
      raise SingleRAdvancedError(
        "SingleCellExperiment assay not found: " + assay_name,
      )
  }
  let data = SingleRAdvancedData::create(
    assay,
    gene_names=experiment.row_names,
    cell_names=experiment.col_names,
  )
  let clusters = if cluster_column == "" {
    []
  } else {
    match experiment.col_data.get(cluster_column) {
      Some(values) => {
        if values.length() != data.n_cells {
          raise SingleRAdvancedError(
            "cluster column length must match the number of cells",
          )
        }
        values.copy()
      }
      None =>
        raise SingleRAdvancedError(
          "SingleCellExperiment cluster column not found: " + cluster_column,
        )
    }
  }
  let result = if clusters.length() == 0 {
    single_r_advanced(data, reference, restrict_genes~, config~)
  } else {
    single_r_advanced_clusters(
      data,
      clusters,
      reference,
      restrict_genes~,
      config~,
    )
  }
  let labels = if clusters.length() == 0 {
    result.labels.copy()
  } else {
    single_r_adv_expand_cluster_strings(
      clusters,
      result.cell_names,
      result.labels,
    )
  }
  let pruned = if clusters.length() == 0 {
    result.pruned_labels.copy()
  } else {
    single_r_adv_expand_cluster_options(
      clusters,
      result.cell_names,
      result.pruned_labels,
    )
  }
  let next_delta = if clusters.length() == 0 {
    result.delta_next.copy()
  } else {
    single_r_adv_expand_cluster_doubles(
      clusters,
      result.cell_names,
      result.delta_next,
    )
  }
  let median_delta = if clusters.length() == 0 {
    result.delta_median.copy()
  } else {
    single_r_adv_expand_cluster_doubles(
      clusters,
      result.cell_names,
      result.delta_median,
    )
  }
  let scores : Array[String] = []
  for cell in 0..<labels.length() {
    let result_cell = if clusters.length() == 0 {
      cell
    } else {
      single_r_adv_label_index(result.cell_names, clusters[cell])
    }
    scores.push(result.assigned_score(result_cell).to_string())
  }
  let enriched = single_r_adv_copy_sce(experiment)
  enriched.col_data[output_prefix + ".labels"] = labels
  enriched.col_data[output_prefix + ".pruned.labels"] = pruned.map(
    single_r_adv_option_string,
  )
  enriched.col_data[output_prefix + ".score"] = scores
  enriched.col_data[output_prefix + ".delta.next"] = next_delta.map(fn(value) {
    value.to_string()
  })
  enriched.col_data[output_prefix + ".delta.median"] = median_delta.map(fn(
    value,
  ) {
    value.to_string()
  })
  if clusters.length() > 0 {
    enriched.col_data[output_prefix + ".cluster"] = clusters
  }
  enriched.metadata[output_prefix + ".assay"] = assay_name
  enriched.metadata[output_prefix + ".reference"] = reference.reference_name
  enriched.metadata[output_prefix + ".quantile"] = config.quantile.to_string()
  enriched.metadata[output_prefix + ".fineTune"] = config.fine_tune.to_string()
  enriched.metadata[output_prefix + ".markers"] = result.marker_gene_names
    .length()
    .to_string()
  enriched.metadata[output_prefix + ".pruned"] = result.n_pruned().to_string()
  enriched.metadata[output_prefix + ".mode"] = if clusters.length() == 0 {
    "cell"
  } else {
    "cluster"
  }
  SingleRAdvancedSCEOutput::{ experiment: enriched, result }
}

///|
pub fn single_r_advanced_example() -> (
  SingleRAdvancedReference,
  SingleRAdvancedData,
) raise SingleRAdvancedError {
  let genes = [
    "CD3D", "CD3E", "TRAC", "MS4A1", "CD79A", "CD74", "LYZ", "S100A8", "FCGR3A",
    "NKG7", "KLRD1", "GNLY",
  ]
  let reference = SingleRAdvancedReference::create(
    [
      [8.0, 7.5, 0.5, 0.4, 0.5, 0.4, 0.3, 0.2],
      [7.5, 8.0, 0.4, 0.5, 0.4, 0.5, 0.3, 0.2],
      [7.0, 7.2, 0.8, 0.7, 0.6, 0.5, 0.5, 0.4],
      [0.4, 0.5, 8.0, 7.5, 0.4, 0.3, 0.2, 0.2],
      [0.5, 0.4, 7.4, 7.8, 0.5, 0.4, 0.2, 0.3],
      [0.7, 0.6, 7.0, 7.2, 0.8, 0.7, 0.3, 0.3],
      [0.5, 0.6, 0.8, 0.7, 8.0, 7.5, 0.5, 0.4],
      [0.3, 0.4, 0.3, 0.4, 7.5, 8.0, 0.3, 0.2],
      [0.4, 0.5, 0.4, 0.3, 6.8, 7.2, 0.7, 0.6],
      [0.8, 0.7, 0.3, 0.3, 0.5, 0.4, 8.0, 7.5],
      [0.7, 0.8, 0.2, 0.3, 0.4, 0.5, 7.5, 8.0],
      [0.6, 0.7, 0.2, 0.2, 0.3, 0.4, 7.2, 7.8],
    ],
    [
      "T cell", "T cell", "B cell", "B cell", "Monocyte", "Monocyte", "NK cell",
      "NK cell",
    ],
    gene_names=genes,
    sample_names=["T1", "T2", "B1", "B2", "M1", "M2", "NK1", "NK2"],
    reference_name="immune_reference",
  )
  let data = SingleRAdvancedData::create(
    [
      [7.8, 0.4, 0.5, 0.7, 4.0],
      [7.4, 0.5, 0.4, 0.6, 3.8],
      [7.0, 0.7, 0.5, 0.5, 3.6],
      [0.4, 7.8, 0.5, 0.3, 4.1],
      [0.5, 7.5, 0.4, 0.3, 3.9],
      [0.6, 7.0, 0.7, 0.3, 3.7],
      [0.5, 0.8, 7.8, 0.5, 0.7],
      [0.3, 0.4, 7.5, 0.3, 0.4],
      [0.4, 0.4, 7.0, 0.7, 0.5],
      [0.7, 0.3, 0.5, 7.8, 0.5],
      [0.8, 0.2, 0.4, 7.5, 0.4],
      [0.6, 0.2, 0.3, 7.2, 0.3],
    ],
    gene_names=genes,
    cell_names=["cell_T", "cell_B", "cell_M", "cell_NK", "cell_ambiguous"],
  )
  (reference, data)
}
