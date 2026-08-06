///|
/// Portable implementations of selected Bioconductor scuttle 1.23.1
/// utilities that are not covered by this package's scrapper module.
///
/// Matrices follow the Bioconductor convention: features in rows and cells in
/// columns. The module provides batch-aware MAD outlier detection, per-feature
/// quality-control metrics, feature-set aggregation, exact count
/// downsampling, batch coverage equalization and immutable
/// SingleCellExperiment wrappers.

///|
pub suberror ScuttleError {
  ScuttleError(String)
}

///|
pub enum ScuttleOutlierDirection {
  ScuttleBoth
  ScuttleLower
  ScuttleHigher
} derive(Eq, Debug)

///|
pub enum ScuttleBatchSummary {
  ScuttleBatchMedian
  ScuttleBatchMean
  ScuttleBatchGeometricMean
} derive(Eq, Debug)

///|
pub struct ScuttleOutlierConfig {
  nmads : Double
  direction : ScuttleOutlierDirection
  log_transform : Bool
  subset : Array[Int]?
  batches : Array[String]
  share_medians : Bool
  share_mads : Bool
  share_missing : Bool
  min_diff : Double?
}

///|
pub struct ScuttleOutlierResult {
  outliers : Array[Bool?]
  batch_names : Array[String]
  lower : Array[Double]
  higher : Array[Double]
  medians : Array[Double]
  mads : Array[Double]
  threshold_valid : Array[Bool]
}

///|
pub struct ScuttleNamedCellSubset {
  name : String
  indices : Array[Int]
}

///|
pub struct ScuttleFeatureQcMetrics {
  means : Array[Double]
  detected_percent : Array[Double]
  subset_names : Array[String]
  subset_means : Array[Array[Double]]
  subset_detected_percent : Array[Array[Double]]
  subset_ratios : Array[Array[Double]]
}

///|
pub struct ScuttleFeatureQcSceOutput {
  experiment : SingleCellExperiment
  metrics : ScuttleFeatureQcMetrics
}

///|
pub struct ScuttleFeatureSet {
  name : String
  indices : Array[Int]
}

///|
pub struct ScuttleFeatureAggregation {
  names : Array[String]
  values : Array[Array[Double]]
  sizes : Array[Int]
}

///|
pub struct ScuttleBatchDownsampleResult {
  matrix : Array[Array[Double]]
  group_names : Array[String]
  summaries : Array[Double]
  proportions : Array[Double]
  group_index : Array[Int]
}

///|
fn scuttle_fail(message : String) -> Unit raise ScuttleError {
  raise ScuttleError(message)
}

///|
fn scuttle_is_finite(value : Double) -> Bool {
  value == value && value.abs() <= 1.0e300
}

///|
fn scuttle_copy_matrix(matrix : Array[Array[Double]]) -> Array[Array[Double]] {
  let output : Array[Array[Double]] = []
  for row in matrix {
    output.push(row.copy())
  }
  output
}

///|
fn scuttle_matrix_shape(
  matrix : Array[Array[Double]],
  require_nonnegative : Bool,
) -> (Int, Int) raise ScuttleError {
  let rows = matrix.length()
  let columns = if rows > 0 { matrix[0].length() } else { 0 }
  for row in 0..<rows {
    if matrix[row].length() != columns {
      scuttle_fail("scuttle matrices must be rectangular")
    }
    for column in 0..<columns {
      let value = matrix[row][column]
      if !scuttle_is_finite(value) {
        scuttle_fail("scuttle matrices must contain finite values")
      }
      if require_nonnegative && value < 0.0 {
        scuttle_fail("scuttle count matrices must be non-negative")
      }
    }
  }
  (rows, columns)
}

///|
fn scuttle_validate_name(
  name : String,
  kind : String,
) -> Unit raise ScuttleError {
  if name.trim().length() == 0 {
    scuttle_fail("scuttle " + kind + " names must not be empty")
  }
}

///|
fn scuttle_validate_unique_names(
  names : Array[String],
  kind : String,
) -> Unit raise ScuttleError {
  for index in 0..<names.length() {
    scuttle_validate_name(names[index], kind)
    for prior in 0..<index {
      if names[index] == names[prior] {
        scuttle_fail("scuttle " + kind + " names must be unique")
      }
    }
  }
}

///|
fn scuttle_median(values : Array[Double]) -> Double {
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
fn scuttle_mad(values : Array[Double], center : Double) -> Double {
  if values.length() == 0 {
    return 0.0
  }
  let deviations : Array[Double] = []
  for value in values {
    if value <= -1.0e299 && center <= -1.0e299 {
      deviations.push(0.0)
    } else {
      deviations.push((value - center).abs())
    }
  }
  scuttle_median(deviations) * 1.4826
}

///|
fn scuttle_mean(values : Array[Double]) -> Double {
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
fn scuttle_string_index(values : Array[String], target : String) -> Int {
  for index in 0..<values.length() {
    if values[index] == target {
      return index
    }
  }
  -1
}

///|
fn scuttle_sorted_unique(values : Array[String]) -> Array[String] {
  let output : Array[String] = []
  for value in values {
    if scuttle_string_index(output, value) < 0 {
      output.push(value)
    }
  }
  output.sort_by(fn(left : String, right : String) -> Int {
    if left < right {
      -1
    } else if left > right {
      1
    } else {
      0
    }
  })
  output
}

///|
pub fn scuttle_outlier_both() -> ScuttleOutlierDirection {
  ScuttleBoth
}

///|
pub fn scuttle_outlier_lower() -> ScuttleOutlierDirection {
  ScuttleLower
}

///|
pub fn scuttle_outlier_higher() -> ScuttleOutlierDirection {
  ScuttleHigher
}

///|
pub fn scuttle_batch_median() -> ScuttleBatchSummary {
  ScuttleBatchMedian
}

///|
pub fn scuttle_batch_mean() -> ScuttleBatchSummary {
  ScuttleBatchMean
}

///|
pub fn scuttle_batch_geometric_mean() -> ScuttleBatchSummary {
  ScuttleBatchGeometricMean
}

///|
pub fn ScuttleOutlierConfig::create(
  nmads? : Double = 3.0,
  direction? : ScuttleOutlierDirection = ScuttleBoth,
  log_transform? : Bool = false,
  subset? : Array[Int]? = None,
  batches? : Array[String] = [],
  share_medians? : Bool = false,
  share_mads? : Bool = false,
  share_missing? : Bool = true,
  min_diff? : Double? = None,
) -> ScuttleOutlierConfig raise ScuttleError {
  if !scuttle_is_finite(nmads) || nmads < 0.0 {
    scuttle_fail("scuttle nmads must be finite and non-negative")
  }
  match min_diff {
    Some(value) =>
      if !scuttle_is_finite(value) || value < 0.0 {
        scuttle_fail("scuttle min_diff must be finite and non-negative")
      }
    None => ()
  }
  let copied_subset = match subset {
    Some(values) => Some(values.copy())
    None => None
  }
  ScuttleOutlierConfig::{
    nmads,
    direction,
    log_transform,
    subset: copied_subset,
    batches: batches.copy(),
    share_medians,
    share_mads,
    share_missing,
    min_diff,
  }
}

///|
pub fn ScuttleOutlierConfig::default() -> ScuttleOutlierConfig {
  ScuttleOutlierConfig::{
    nmads: 3.0,
    direction: ScuttleBoth,
    log_transform: false,
    subset: None,
    batches: [],
    share_medians: false,
    share_mads: false,
    share_missing: true,
    min_diff: None,
  }
}

///|
fn scuttle_log2_value(value : Double) -> Double {
  if value == 0.0 {
    -1.0e300
  } else {
    @math.log2(value)
  }
}

///|
fn scuttle_restore_log2_threshold(value : Double) -> Double {
  if value <= -1.0e299 {
    0.0
  } else if value >= 1.0e299 {
    1.0e300
  } else {
    @math.exp(value * 0.6931471805599453)
  }
}

///|
/// Identify median/MAD outliers with scuttle-compatible batch sharing.
///
/// Batch names are sorted as R factor levels. When a threshold cannot be
/// estimated, affected output entries are `None` and `threshold_valid` is
/// false for that batch.
pub fn scuttle_is_outlier(
  metric : Array[Double],
  config? : ScuttleOutlierConfig = ScuttleOutlierConfig::default(),
) -> ScuttleOutlierResult raise ScuttleError {
  for value in metric {
    if !scuttle_is_finite(value) {
      scuttle_fail("scuttle outlier metrics must contain finite values")
    }
    if config.log_transform && value < 0.0 {
      scuttle_fail(
        "scuttle log-transformed outlier metrics must be non-negative",
      )
    }
  }
  if config.batches.length() > 0 && config.batches.length() != metric.length() {
    scuttle_fail("scuttle batch length must equal metric length")
  }
  for batch in config.batches {
    scuttle_validate_name(batch, "batch")
  }
  let effective_batches = if config.batches.length() == 0 {
    Array::make(metric.length(), "1")
  } else {
    config.batches.copy()
  }
  let batch_names = if effective_batches.length() == 0 {
    ["1"]
  } else {
    scuttle_sorted_unique(effective_batches)
  }
  let selected : Array[Int] = match config.subset {
    Some(indices) => {
      for index in indices {
        if index < 0 || index >= metric.length() {
          scuttle_fail("scuttle outlier subset index is out of bounds")
        }
      }
      indices.copy()
    }
    None => {
      let indices : Array[Int] = []
      for index in 0..<metric.length() {
        indices.push(index)
      }
      indices
    }
  }
  let transformed : Array[Double] = []
  for value in metric {
    transformed.push(
      if config.log_transform {
        scuttle_log2_value(value)
      } else {
        value
      },
    )
  }
  let selected_by_batch : Array[Array[Double]] = []
  for _ in batch_names {
    selected_by_batch.push([])
  }
  for index in selected {
    let batch = scuttle_string_index(batch_names, effective_batches[index])
    selected_by_batch[batch].push(transformed[index])
  }
  let all_selected : Array[Double] = []
  for index in selected {
    all_selected.push(transformed[index])
  }
  let has_global = all_selected.length() > 0
  let global_median = if has_global {
    scuttle_median(all_selected)
  } else {
    0.0
  }
  let medians = Array::make(batch_names.length(), 0.0)
  let median_valid = Array::make(batch_names.length(), false)
  for batch in 0..<batch_names.length() {
    let has_batch = selected_by_batch[batch].length() > 0
    if !config.share_medians && has_batch {
      medians[batch] = scuttle_median(selected_by_batch[batch])
      median_valid[batch] = true
    }
    let use_shared = if config.share_medians && config.share_missing {
      true
    } else if config.share_missing {
      !has_batch
    } else if config.share_medians {
      has_batch
    } else {
      false
    }
    if use_shared && has_global {
      medians[batch] = global_median
      median_valid[batch] = true
    }
  }
  let mads = Array::make(batch_names.length(), 0.0)
  let mad_valid = Array::make(batch_names.length(), false)
  for batch in 0..<batch_names.length() {
    if selected_by_batch[batch].length() > 0 &&
      median_valid[batch] &&
      !config.share_mads {
      mads[batch] = scuttle_mad(selected_by_batch[batch], medians[batch])
      mad_valid[batch] = true
    }
  }
  let shared_deviations : Array[Double] = []
  if has_global {
    for index in selected {
      let batch = scuttle_string_index(batch_names, effective_batches[index])
      if median_valid[batch] {
        let value = transformed[index]
        if value <= -1.0e299 && medians[batch] <= -1.0e299 {
          shared_deviations.push(0.0)
        } else {
          shared_deviations.push((value - medians[batch]).abs())
        }
      }
    }
  }
  let has_shared_mad = shared_deviations.length() > 0
  let shared_mad = if has_shared_mad {
    scuttle_median(shared_deviations) * 1.4826
  } else {
    0.0
  }
  for batch in 0..<batch_names.length() {
    let has_batch = selected_by_batch[batch].length() > 0
    let use_shared = if config.share_mads && config.share_missing {
      true
    } else if config.share_missing {
      !has_batch
    } else if config.share_mads {
      has_batch
    } else {
      false
    }
    if use_shared && median_valid[batch] && has_shared_mad {
      mads[batch] = shared_mad
      mad_valid[batch] = true
    }
  }
  let lower = Array::make(batch_names.length(), 0.0)
  let higher = Array::make(batch_names.length(), 0.0)
  let threshold_valid = Array::make(batch_names.length(), false)
  for batch in 0..<batch_names.length() {
    if median_valid[batch] && mad_valid[batch] {
      let mad_difference = config.nmads * mads[batch]
      let difference = match config.min_diff {
        Some(value) => value.max(mad_difference)
        None => mad_difference
      }
      lower[batch] = match config.direction {
        ScuttleHigher => -1.0e300
        _ => medians[batch] - difference
      }
      higher[batch] = match config.direction {
        ScuttleLower => 1.0e300
        _ => medians[batch] + difference
      }
      threshold_valid[batch] = true
    }
  }
  let outliers : Array[Bool?] = []
  for index in 0..<metric.length() {
    let batch = scuttle_string_index(batch_names, effective_batches[index])
    if threshold_valid[batch] {
      outliers.push(
        Some(
          transformed[index] < lower[batch] ||
          transformed[index] > higher[batch],
        ),
      )
    } else {
      outliers.push(None)
    }
  }
  if config.log_transform {
    for batch in 0..<batch_names.length() {
      if threshold_valid[batch] {
        lower[batch] = scuttle_restore_log2_threshold(lower[batch])
        higher[batch] = scuttle_restore_log2_threshold(higher[batch])
        medians[batch] = scuttle_restore_log2_threshold(medians[batch])
      }
    }
  }
  ScuttleOutlierResult::{
    outliers,
    batch_names,
    lower,
    higher,
    medians,
    mads,
    threshold_valid,
  }
}

///|
pub fn ScuttleOutlierResult::discarded_count(
  self : ScuttleOutlierResult,
) -> Int {
  let mut count = 0
  for value in self.outliers {
    match value {
      Some(true) => count = count + 1
      _ => ()
    }
  }
  count
}

///|
pub fn ScuttleNamedCellSubset::create(
  name : String,
  indices : Array[Int],
) -> ScuttleNamedCellSubset raise ScuttleError {
  scuttle_validate_name(name, "cell subset")
  ScuttleNamedCellSubset::{ name, indices: indices.copy() }
}

///|
/// Compute per-feature means and detection percentages, including optional
/// subsets of cells and subset/global mean ratios.
pub fn scuttle_per_feature_qc(
  counts : Array[Array[Double]],
  subsets? : Array[ScuttleNamedCellSubset] = [],
  detection_limit? : Double = 0.0,
) -> ScuttleFeatureQcMetrics raise ScuttleError {
  let (features, cells) = scuttle_matrix_shape(counts, true)
  if !scuttle_is_finite(detection_limit) {
    scuttle_fail("scuttle detection_limit must be finite")
  }
  let subset_names : Array[String] = []
  for subset in subsets {
    subset_names.push(subset.name)
  }
  scuttle_validate_unique_names(subset_names, "cell subset")
  for subset in subsets {
    for index in subset.indices {
      if index < 0 || index >= cells {
        scuttle_fail("scuttle cell subset index is out of bounds")
      }
    }
  }
  let means = Array::make(features, 0.0)
  let detected_percent = Array::make(features, 0.0)
  for feature in 0..<features {
    let mut total = 0.0
    let mut detected = 0
    for cell in 0..<cells {
      let value = counts[feature][cell]
      total = total + value
      if value > detection_limit {
        detected = detected + 1
      }
    }
    if cells > 0 {
      means[feature] = total / cells.to_double()
      detected_percent[feature] = detected.to_double() *
        100.0 /
        cells.to_double()
    }
  }
  let subset_means : Array[Array[Double]] = []
  let subset_detected_percent : Array[Array[Double]] = []
  let subset_ratios : Array[Array[Double]] = []
  for subset in subsets {
    let current_means = Array::make(features, 0.0)
    let current_detected = Array::make(features, 0.0)
    let current_ratios = Array::make(features, 0.0)
    for feature in 0..<features {
      let mut total = 0.0
      let mut detected = 0
      for cell in subset.indices {
        let value = counts[feature][cell]
        total = total + value
        if value > detection_limit {
          detected = detected + 1
        }
      }
      if subset.indices.length() > 0 {
        current_means[feature] = total / subset.indices.length().to_double()
        current_detected[feature] = detected.to_double() *
          100.0 /
          subset.indices.length().to_double()
      }
      if means[feature] != 0.0 {
        current_ratios[feature] = current_means[feature] / means[feature]
      }
    }
    subset_means.push(current_means)
    subset_detected_percent.push(current_detected)
    subset_ratios.push(current_ratios)
  }
  ScuttleFeatureQcMetrics::{
    means,
    detected_percent,
    subset_names,
    subset_means,
    subset_detected_percent,
    subset_ratios,
  }
}

///|
pub fn ScuttleFeatureSet::create(
  name : String,
  indices : Array[Int],
) -> ScuttleFeatureSet raise ScuttleError {
  scuttle_validate_name(name, "feature set")
  ScuttleFeatureSet::{ name, indices: indices.copy() }
}

///|
fn scuttle_validate_feature_sets(
  feature_sets : Array[ScuttleFeatureSet],
  features : Int,
) -> Unit raise ScuttleError {
  let names : Array[String] = []
  for feature_set in feature_sets {
    names.push(feature_set.name)
  }
  scuttle_validate_unique_names(names, "feature set")
  for feature_set in feature_sets {
    for index in feature_set.indices {
      if index < 0 || index >= features {
        scuttle_fail("scuttle feature set index is out of bounds")
      }
    }
  }
}

///|
/// Aggregate arbitrary, potentially overlapping feature sets.
pub fn scuttle_aggregate_feature_sets(
  expression : Array[Array[Double]],
  feature_sets : Array[ScuttleFeatureSet],
  average? : Bool = false,
  detection_limit? : Double? = None,
) -> ScuttleFeatureAggregation raise ScuttleError {
  let (features, cells) = scuttle_matrix_shape(expression, false)
  scuttle_validate_feature_sets(feature_sets, features)
  match detection_limit {
    Some(value) =>
      if !scuttle_is_finite(value) {
        scuttle_fail("scuttle detection_limit must be finite")
      }
    None => ()
  }
  let names : Array[String] = []
  let values : Array[Array[Double]] = []
  let sizes : Array[Int] = []
  for feature_set in feature_sets {
    names.push(feature_set.name)
    sizes.push(feature_set.indices.length())
    let row = Array::make(cells, 0.0)
    for feature in feature_set.indices {
      for cell in 0..<cells {
        let contribution = match detection_limit {
          Some(limit) =>
            if expression[feature][cell] > limit {
              1.0
            } else {
              0.0
            }
          None => expression[feature][cell]
        }
        row[cell] = row[cell] + contribution
      }
    }
    if average && feature_set.indices.length() > 0 {
      for cell in 0..<cells {
        row[cell] = row[cell] / feature_set.indices.length().to_double()
      }
    }
    values.push(row)
  }
  ScuttleFeatureAggregation::{ names, values, sizes }
}

///|
/// Aggregate each feature according to a string identifier. Empty identifiers
/// are ignored, matching scuttle's treatment of NA factor entries.
pub fn scuttle_aggregate_features_by_ids(
  expression : Array[Array[Double]],
  ids : Array[String],
  average? : Bool = false,
  detection_limit? : Double? = None,
) -> ScuttleFeatureAggregation raise ScuttleError {
  let (features, _) = scuttle_matrix_shape(expression, false)
  if ids.length() != features {
    scuttle_fail("scuttle feature ids must match the number of rows")
  }
  let names = scuttle_sorted_unique(ids)
  let retained : Array[String] = []
  for name in names {
    if name.length() > 0 {
      retained.push(name)
    }
  }
  let feature_sets : Array[ScuttleFeatureSet] = []
  for name in retained {
    let indices : Array[Int] = []
    for feature in 0..<ids.length() {
      if ids[feature] == name {
        indices.push(feature)
      }
    }
    feature_sets.push(ScuttleFeatureSet::{ name, indices })
  }
  scuttle_aggregate_feature_sets(
    expression,
    feature_sets,
    average~,
    detection_limit~,
  )
}

///|
fn scuttle_sanitize_counts(
  counts : Array[Array[Double]],
) -> Array[Array[Int]] raise ScuttleError {
  let (features, cells) = scuttle_matrix_shape(counts, false)
  let output : Array[Array[Int]] = []
  let mut total = 0
  for feature in 0..<features {
    let row : Array[Int] = []
    for cell in 0..<cells {
      let rounded = counts[feature][cell].round()
      let value = if rounded < 0.0 { 0 } else { rounded.to_int() }
      if value > 2147483646 - total {
        scuttle_fail(
          "scuttle count total exceeds reproducible integer precision",
        )
      }
      total = total + value
      row.push(value)
    }
    output.push(row)
  }
  output
}

///|
fn scuttle_lcg(state : Int) -> (Int, Double) {
  let normalized = (state - 1) % 2147483646 + 1
  let quotient = normalized / 44488
  let remainder = normalized % 44488
  let candidate = 48271 * remainder - 3399 * quotient
  let next = if candidate > 0 { candidate } else { candidate + 2147483647 }
  (next, next.to_double() / 2147483647.0)
}

///|
fn scuttle_sample_categories(
  counts : Array[Int],
  required : Int,
  seed : Int,
) -> (Array[Int], Int) {
  let mut total = 0
  for value in counts {
    total = total + value
  }
  if required <= 0 {
    return (Array::make(counts.length(), 0), seed)
  }
  if required >= total {
    return (counts.copy(), seed)
  }
  let sample_kept = required <= total - required
  let draws = if sample_kept { required } else { total - required }
  let pool = counts.copy()
  let sampled = Array::make(counts.length(), 0)
  let mut remaining = total
  let mut state = seed
  for _ in 0..<draws {
    let (next, uniform) = scuttle_lcg(state)
    state = next
    let target = (uniform * remaining.to_double()).to_int()
    let mut cumulative = 0
    let mut chosen = pool.length() - 1
    for category in 0..<pool.length() {
      cumulative = cumulative + pool[category]
      if target < cumulative {
        chosen = category
        break
      }
    }
    sampled[chosen] = sampled[chosen] + 1
    pool[chosen] = pool[chosen] - 1
    remaining = remaining - 1
  }
  if sample_kept {
    (sampled, state)
  } else {
    let kept = Array::make(counts.length(), 0)
    for category in 0..<counts.length() {
      kept[category] = counts[category] - sampled[category]
    }
    (kept, state)
  }
}

///|
fn scuttle_validate_proportion(value : Double) -> Unit raise ScuttleError {
  if !scuttle_is_finite(value) || value < 0.0 || value > 1.0 {
    scuttle_fail("scuttle downsampling proportions must lie in [0, 1]")
  }
}

///|
/// Downsample each cell without replacement to an exact rounded total.
pub fn scuttle_downsample_columns(
  counts : Array[Array[Double]],
  proportions : Array[Double],
  seed? : Int = 1,
) -> Array[Array[Double]] raise ScuttleError {
  if seed <= 0 {
    scuttle_fail("scuttle random seed must be positive")
  }
  let sanitized = scuttle_sanitize_counts(counts)
  let features = sanitized.length()
  let cells = if features > 0 { sanitized[0].length() } else { 0 }
  if proportions.length() != cells {
    scuttle_fail("scuttle proportion count must equal the number of cells")
  }
  for proportion in proportions {
    scuttle_validate_proportion(proportion)
  }
  let output : Array[Array[Double]] = []
  for _ in 0..<features {
    output.push(Array::make(cells, 0.0))
  }
  let mut state = seed
  for cell in 0..<cells {
    let category_counts : Array[Int] = []
    let mut total = 0
    for feature in 0..<features {
      category_counts.push(sanitized[feature][cell])
      total = total + sanitized[feature][cell]
    }
    let required = (total.to_double() * proportions[cell]).round().to_int()
    let (sampled, next) = scuttle_sample_categories(
      category_counts, required, state,
    )
    state = next
    for feature in 0..<features {
      output[feature][cell] = sampled[feature].to_double()
    }
  }
  output
}

///|
/// Downsample a count matrix either per cell or over the entire matrix.
pub fn scuttle_downsample_matrix(
  counts : Array[Array[Double]],
  proportion : Double,
  by_column? : Bool = true,
  seed? : Int = 1,
) -> Array[Array[Double]] raise ScuttleError {
  if seed <= 0 {
    scuttle_fail("scuttle random seed must be positive")
  }
  scuttle_validate_proportion(proportion)
  let sanitized = scuttle_sanitize_counts(counts)
  let features = sanitized.length()
  let cells = if features > 0 { sanitized[0].length() } else { 0 }
  if by_column {
    return scuttle_downsample_columns(
      counts,
      Array::make(cells, proportion),
      seed~,
    )
  }
  let flattened : Array[Int] = []
  let mut total = 0
  for cell in 0..<cells {
    for feature in 0..<features {
      flattened.push(sanitized[feature][cell])
      total = total + sanitized[feature][cell]
    }
  }
  let required = (total.to_double() * proportion).round().to_int()
  let (sampled, _) = scuttle_sample_categories(flattened, required, seed)
  let output : Array[Array[Double]] = []
  for _ in 0..<features {
    output.push(Array::make(cells, 0.0))
  }
  let mut position = 0
  for cell in 0..<cells {
    for feature in 0..<features {
      output[feature][cell] = sampled[position].to_double()
      position = position + 1
    }
  }
  output
}

///|
fn scuttle_column_sums(matrix : Array[Array[Double]]) -> Array[Double] {
  let cells = if matrix.length() > 0 { matrix[0].length() } else { 0 }
  let sums = Array::make(cells, 0.0)
  for row in matrix {
    for cell in 0..<cells {
      sums[cell] = sums[cell] + row[cell]
    }
  }
  sums
}

///|
fn scuttle_batch_summary(
  values : Array[Double],
  summary_method : ScuttleBatchSummary,
) -> Double {
  match summary_method {
    ScuttleBatchMean => scuttle_mean(values)
    ScuttleBatchMedian => scuttle_median(values)
    ScuttleBatchGeometricMean => {
      if values.length() == 0 {
        return 0.0
      }
      let mut total = 0.0
      for value in values {
        total = total + @math.ln(value + 1.0)
      }
      @math.exp(total / values.length().to_double())
    }
  }
}

///|
/// Equalize average coverage across batches, optionally within blocks.
pub fn scuttle_downsample_batches(
  counts : Array[Array[Double]],
  batches : Array[String],
  blocks? : Array[String] = [],
  summary_method? : ScuttleBatchSummary = ScuttleBatchMedian,
  by_column? : Bool = true,
  seed? : Int = 1,
) -> ScuttleBatchDownsampleResult raise ScuttleError {
  if seed <= 0 {
    scuttle_fail("scuttle random seed must be positive")
  }
  let (_, cells) = scuttle_matrix_shape(counts, false)
  if batches.length() != cells {
    scuttle_fail("scuttle batch length must equal the number of cells")
  }
  if blocks.length() > 0 && blocks.length() != cells {
    scuttle_fail("scuttle block length must equal the number of cells")
  }
  for batch in batches {
    scuttle_validate_name(batch, "batch")
  }
  for block in blocks {
    scuttle_validate_name(block, "block")
  }
  let effective_blocks = if blocks.length() == 0 {
    Array::make(cells, "1")
  } else {
    blocks.copy()
  }
  let group_names : Array[String] = []
  let group_blocks : Array[String] = []
  let group_indices : Array[Array[Int]] = []
  let group_index = Array::make(cells, 0)
  for cell in 0..<cells {
    let group_name = if blocks.length() == 0 {
      batches[cell]
    } else {
      effective_blocks[cell] + "|" + batches[cell]
    }
    let mut group = scuttle_string_index(group_names, group_name)
    if group < 0 {
      group = group_names.length()
      group_names.push(group_name)
      group_blocks.push(effective_blocks[cell])
      group_indices.push([])
    }
    group_indices[group].push(cell)
    group_index[cell] = group
  }
  let libraries = scuttle_column_sums(counts)
  let summaries = Array::make(group_names.length(), 0.0)
  for group in 0..<group_names.length() {
    let values : Array[Double] = []
    for cell in group_indices[group] {
      values.push(libraries[cell])
    }
    summaries[group] = scuttle_batch_summary(values, summary_method)
  }
  let proportions = Array::make(group_names.length(), 1.0)
  for group in 0..<group_names.length() {
    let mut target = summaries[group]
    for other in 0..<group_names.length() {
      if group_blocks[other] == group_blocks[group] {
        target = target.min(summaries[other])
      }
    }
    proportions[group] = if summaries[group] > 0.0 {
      target / summaries[group]
    } else {
      1.0
    }
  }
  let sanitized = scuttle_sanitize_counts(counts)
  let features = sanitized.length()
  let output : Array[Array[Double]] = []
  for _ in 0..<features {
    output.push(Array::make(cells, 0.0))
  }
  for group in 0..<group_names.length() {
    let submatrix : Array[Array[Double]] = []
    for feature in 0..<features {
      let row : Array[Double] = []
      for cell in group_indices[group] {
        row.push(sanitized[feature][cell].to_double())
      }
      submatrix.push(row)
    }
    let group_seed = (seed + group * 7919 - 1) % 2147483646 + 1
    let downsampled = scuttle_downsample_matrix(
      submatrix,
      proportions[group],
      by_column~,
      seed=group_seed,
    )
    for feature in 0..<features {
      for position in 0..<group_indices[group].length() {
        output[feature][group_indices[group][position]] = downsampled[feature][position]
      }
    }
  }
  ScuttleBatchDownsampleResult::{
    matrix: output,
    group_names,
    summaries,
    proportions,
    group_index,
  }
}

///|
fn scuttle_copy_sce(source : SingleCellExperiment) -> SingleCellExperiment {
  let assays : Map[String, Array[Array[Double]]] = Map([])
  for key in source.assays.keys() {
    assays[key] = scuttle_copy_matrix(source.assays[key])
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
    reduced_dims[key] = scuttle_copy_matrix(source.reduced_dims[key])
  }
  let alternatives : Map[String, SingleCellExperiment] = Map([])
  for key in source.alternative_experiments.keys() {
    alternatives[key] = scuttle_copy_sce(source.alternative_experiments[key])
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
fn scuttle_doubles_to_strings(values : Array[Double]) -> Array[String] {
  let output : Array[String] = []
  for value in values {
    output.push(value.to_string())
  }
  output
}

///|
/// Compute per-feature QC metrics from an SCE assay and write them to a deep
/// copy of row_data.
pub fn scuttle_per_feature_qc_sce(
  experiment : SingleCellExperiment,
  assay_name? : String = "counts",
  subsets? : Array[ScuttleNamedCellSubset] = [],
  detection_limit? : Double = 0.0,
  output_prefix? : String = "scuttle",
) -> ScuttleFeatureQcSceOutput raise ScuttleError {
  scuttle_validate_name(output_prefix, "output prefix")
  let counts = match experiment.assays.get(assay_name) {
    Some(value) => value
    None => {
      scuttle_fail("scuttle assay '" + assay_name + "' was not found")
      []
    }
  }
  let (features, cells) = scuttle_matrix_shape(counts, true)
  if experiment.row_names.length() != features ||
    experiment.col_names.length() != cells {
    scuttle_fail("scuttle SCE dimensions do not match row or column names")
  }
  let metrics = scuttle_per_feature_qc(counts, subsets~, detection_limit~)
  let output = scuttle_copy_sce(experiment)
  output.row_data[output_prefix + ".mean"] = scuttle_doubles_to_strings(
    metrics.means,
  )
  output.row_data[output_prefix + ".detected"] = scuttle_doubles_to_strings(
    metrics.detected_percent,
  )
  for subset in 0..<metrics.subset_names.length() {
    let base = output_prefix + ".subsets." + metrics.subset_names[subset]
    output.row_data[base + ".mean"] = scuttle_doubles_to_strings(
      metrics.subset_means[subset],
    )
    output.row_data[base + ".detected"] = scuttle_doubles_to_strings(
      metrics.subset_detected_percent[subset],
    )
    output.row_data[base + ".ratio"] = scuttle_doubles_to_strings(
      metrics.subset_ratios[subset],
    )
  }
  output.metadata[output_prefix + ".assay"] = assay_name
  ScuttleFeatureQcSceOutput::{ experiment: output, metrics }
}

///|
/// Aggregate selected SCE assays across arbitrary feature sets. Cell metadata,
/// reduced dimensions, alternative experiments and metadata are deep-copied;
/// feature-level metadata are discarded because gene-level annotations are no
/// longer meaningful.
pub fn scuttle_aggregate_feature_sets_sce(
  experiment : SingleCellExperiment,
  feature_sets : Array[ScuttleFeatureSet],
  assay_names? : Array[String] = ["counts"],
  average? : Bool = false,
) -> SingleCellExperiment raise ScuttleError {
  if assay_names.length() == 0 {
    scuttle_fail("scuttle SCE aggregation requires at least one assay")
  }
  scuttle_validate_unique_names(assay_names, "assay")
  let assays : Map[String, Array[Array[Double]]] = Map([])
  let mut expected_features = -1
  let mut expected_cells = -1
  for assay_name in assay_names {
    let matrix = match experiment.assays.get(assay_name) {
      Some(value) => value
      None => {
        scuttle_fail("scuttle assay '" + assay_name + "' was not found")
        []
      }
    }
    let (features, cells) = scuttle_matrix_shape(matrix, false)
    if expected_features < 0 {
      expected_features = features
      expected_cells = cells
      scuttle_validate_feature_sets(feature_sets, features)
    } else if features != expected_features || cells != expected_cells {
      scuttle_fail("scuttle SCE assays must have identical dimensions")
    }
    assays[assay_name] = scuttle_aggregate_feature_sets(
        matrix,
        feature_sets,
        average~,
      ).values
  }
  if experiment.row_names.length() != expected_features ||
    experiment.col_names.length() != expected_cells {
    scuttle_fail("scuttle SCE dimensions do not match row or column names")
  }
  let col_data : Map[String, Array[String]] = Map([])
  for key in experiment.col_data.keys() {
    col_data[key] = experiment.col_data[key].copy()
  }
  let reduced_dims : Map[String, Array[Array[Double]]] = Map([])
  for key in experiment.reduced_dims.keys() {
    reduced_dims[key] = scuttle_copy_matrix(experiment.reduced_dims[key])
  }
  let alternatives : Map[String, SingleCellExperiment] = Map([])
  for key in experiment.alternative_experiments.keys() {
    alternatives[key] = scuttle_copy_sce(
      experiment.alternative_experiments[key],
    )
  }
  let metadata : Map[String, String] = Map([])
  for key in experiment.metadata.keys() {
    metadata[key] = experiment.metadata[key]
  }
  metadata["scuttle.aggregate.average"] = average.to_string()
  let names : Array[String] = []
  for feature_set in feature_sets {
    names.push(feature_set.name)
  }
  SingleCellExperiment::{
    assays,
    row_data: Map([]),
    col_data,
    reduced_dims,
    alternative_experiments: alternatives,
    metadata,
    row_names: names,
    col_names: experiment.col_names.copy(),
  }
}

///|
/// Add an exactly downsampled assay to a deep copy of an SCE.
pub fn scuttle_downsample_sce(
  experiment : SingleCellExperiment,
  proportion : Double,
  assay_name? : String = "counts",
  output_assay? : String = "downsampled",
  by_column? : Bool = true,
  seed? : Int = 1,
) -> SingleCellExperiment raise ScuttleError {
  scuttle_validate_name(output_assay, "output assay")
  let counts = match experiment.assays.get(assay_name) {
    Some(value) => value
    None => {
      scuttle_fail("scuttle assay '" + assay_name + "' was not found")
      []
    }
  }
  let (features, cells) = scuttle_matrix_shape(counts, false)
  if experiment.row_names.length() != features ||
    experiment.col_names.length() != cells {
    scuttle_fail("scuttle SCE dimensions do not match row or column names")
  }
  let downsampled = scuttle_downsample_matrix(
    counts,
    proportion,
    by_column~,
    seed~,
  )
  let output = scuttle_copy_sce(experiment)
  output.assays[output_assay] = downsampled
  output.metadata["scuttle.downsample.source"] = assay_name
  output.metadata["scuttle.downsample.proportion"] = proportion.to_string()
  output.metadata["scuttle.downsample.by_column"] = by_column.to_string()
  output.metadata["scuttle.downsample.seed"] = seed.to_string()
  output
}
