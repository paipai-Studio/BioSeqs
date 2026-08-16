///|
/// Core single-cell preprocessing utilities inspired by Bioconductor scrapper.
///
/// Matrices use the Bioconductor convention: features in rows and cells in
/// columns. The implementation covers RNA quality control, size-factor
/// normalization, mean-variance modelling, HVG selection and pseudo-bulk
/// aggregation.

///|
pub suberror ScrapperError {
  ScrapperError(String)
}

///|
pub struct ScrapperNamedSubset {
  name : String
  indices : Array[Int]
}

///|
pub struct ScrapperRnaQcMetrics {
  sums : Array[Double]
  detected : Array[Int]
  subset_names : Array[String]
  subset_proportions : Array[Array[Double]]
}

///|
pub struct ScrapperRnaQcThresholds {
  block_names : Array[String]
  sum_lower : Array[Double]
  detected_lower : Array[Double]
  subset_upper : Array[Array[Double]]
}

///|
pub struct ScrapperRnaQcResult {
  metrics : ScrapperRnaQcMetrics
  thresholds : ScrapperRnaQcThresholds
  keep : Array[Bool]
}

///|
pub enum ScrapperCenterMode {
  Lowest
  PerBlock
} derive(Eq, Debug)

///|
pub struct ScrapperVarianceTrend {
  fitted : Array[Double]
  residuals : Array[Double]
}

///|
pub struct ScrapperVarianceModel {
  means : Array[Double]
  variances : Array[Double]
  fitted : Array[Double]
  residuals : Array[Double]
}

///|
pub struct ScrapperFactor {
  name : String
  levels : Array[String]
}

///|
pub struct ScrapperAggregateResult {
  factor_names : Array[String]
  combinations : Array[Array[String]]
  group_names : Array[String]
  counts : Array[Int]
  index : Array[Int]
  sums : Array[Array[Double]]
  detected : Array[Array[Int]]
  medians : Array[Array[Double]]
}

///|
pub fn ScrapperNamedSubset::new(
  name : String,
  indices : Array[Int],
) -> ScrapperNamedSubset {
  ScrapperNamedSubset::{ name, indices: indices.copy() }
}

///|
pub fn ScrapperFactor::new(
  name : String,
  levels : Array[String],
) -> ScrapperFactor {
  ScrapperFactor::{ name, levels: levels.copy() }
}

///|
pub fn scrapper_center_lowest() -> ScrapperCenterMode {
  Lowest
}

///|
pub fn scrapper_center_per_block() -> ScrapperCenterMode {
  PerBlock
}

///|
fn scrapper_is_finite(value : Double) -> Bool {
  value == value && value.abs() < 1.0e300
}

///|
fn scrapper_matrix_shape(
  matrix : Array[Array[Double]],
  require_counts : Bool,
) -> (Int, Int) raise ScrapperError {
  let rows = matrix.length()
  let columns = if rows > 0 { matrix[0].length() } else { 0 }
  let mut row = 0
  while row < rows {
    if matrix[row].length() != columns {
      raise ScrapperError("scrapper matrices must be rectangular")
    }
    let mut column = 0
    while column < columns {
      let value = matrix[row][column]
      if !scrapper_is_finite(value) {
        raise ScrapperError("scrapper matrices must contain finite values")
      }
      if require_counts && value < 0.0 {
        raise ScrapperError("scrapper count matrices must be non-negative")
      }
      column = column + 1
    }
    row = row + 1
  }
  (rows, columns)
}

///|
fn scrapper_copy_matrix(matrix : Array[Array[Double]]) -> Array[Array[Double]] {
  let output : Array[Array[Double]] = []
  let mut row = 0
  while row < matrix.length() {
    output.push(matrix[row].copy())
    row = row + 1
  }
  output
}

///|
fn scrapper_median(values : Array[Double]) -> Double {
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
fn scrapper_mad(values : Array[Double], center : Double) -> Double {
  if values.length() == 0 {
    return 0.0
  }
  let deviations : Array[Double] = []
  let mut index = 0
  while index < values.length() {
    deviations.push((values[index] - center).abs())
    index = index + 1
  }
  scrapper_median(deviations) * 1.4826
}

///|
fn scrapper_mean(values : Array[Double]) -> Double {
  if values.length() == 0 {
    return 0.0
  }
  let mut total = 0.0
  let mut index = 0
  while index < values.length() {
    total = total + values[index]
    index = index + 1
  }
  total / values.length().to_double()
}

///|
fn scrapper_subset_masks(
  subsets : Array[ScrapperNamedSubset],
  features : Int,
) -> Array[Array[Bool]] raise ScrapperError {
  let masks : Array[Array[Bool]] = []
  let seen_names : Map[String, Bool] = Map([], capacity=subsets.length())
  let mut subset = 0
  while subset < subsets.length() {
    let current = subsets[subset]
    if current.name.length() == 0 {
      raise ScrapperError("scrapper subset names must not be empty")
    }
    if seen_names.contains(current.name) {
      raise ScrapperError("scrapper subset names must be unique")
    }
    seen_names[current.name] = true
    let mask = Array::make(features, false)
    let mut position = 0
    while position < current.indices.length() {
      let index = current.indices[position]
      if index < 0 || index >= features {
        raise ScrapperError("scrapper subset index is out of bounds")
      }
      mask[index] = true
      position = position + 1
    }
    masks.push(mask)
    subset = subset + 1
  }
  masks
}

///|
/// Compute per-cell RNA QC metrics.
///
/// `subset_proportions[s][c]` is the proportion of cell `c`'s total count
/// assigned to subset `s`.
pub fn scrapper_compute_rna_qc_metrics(
  counts : Array[Array[Double]],
  subsets? : Array[ScrapperNamedSubset] = [],
  detection_limit? : Double = 0.0,
) -> ScrapperRnaQcMetrics raise ScrapperError {
  let (features, cells) = scrapper_matrix_shape(counts, true)
  if !scrapper_is_finite(detection_limit) {
    raise ScrapperError("scrapper detection limit must be finite")
  }
  let masks = scrapper_subset_masks(subsets, features)
  let sums = Array::make(cells, 0.0)
  let detected = Array::make(cells, 0)
  let subset_sums : Array[Array[Double]] = []
  let mut subset = 0
  while subset < subsets.length() {
    subset_sums.push(Array::make(cells, 0.0))
    subset = subset + 1
  }
  let mut feature = 0
  while feature < features {
    let mut cell = 0
    while cell < cells {
      let value = counts[feature][cell]
      sums[cell] = sums[cell] + value
      if value > detection_limit {
        detected[cell] = detected[cell] + 1
      }
      subset = 0
      while subset < subsets.length() {
        if masks[subset][feature] {
          subset_sums[subset][cell] = subset_sums[subset][cell] + value
        }
        subset = subset + 1
      }
      cell = cell + 1
    }
    feature = feature + 1
  }
  let subset_proportions : Array[Array[Double]] = []
  subset = 0
  while subset < subsets.length() {
    let proportions = Array::make(cells, 0.0)
    let mut cell = 0
    while cell < cells {
      if sums[cell] > 0.0 {
        proportions[cell] = subset_sums[subset][cell] / sums[cell]
      }
      cell = cell + 1
    }
    subset_proportions.push(proportions)
    subset = subset + 1
  }
  let subset_names : Array[String] = []
  subset = 0
  while subset < subsets.length() {
    subset_names.push(subsets[subset].name)
    subset = subset + 1
  }
  ScrapperRnaQcMetrics::{ sums, detected, subset_names, subset_proportions }
}

///|
fn scrapper_block_layout(
  cells : Int,
  blocks : Array[String],
) -> (Array[String], Array[Int]) raise ScrapperError {
  if blocks.length() != 0 && blocks.length() != cells {
    raise ScrapperError(
      "scrapper block vector length must equal the number of cells",
    )
  }
  if blocks.length() == 0 {
    return (["all"], Array::make(cells, 0))
  }
  let names : Array[String] = []
  let mapping : Map[String, Int] = Map([], capacity=blocks.length())
  let index = Array::make(cells, 0)
  let mut cell = 0
  while cell < cells {
    let name = blocks[cell]
    if name.length() == 0 {
      raise ScrapperError("scrapper block names must not be empty")
    }
    match mapping.get(name) {
      Some(found) => index[cell] = found
      None => {
        let found = names.length()
        names.push(name)
        mapping[name] = found
        index[cell] = found
      }
    }
    cell = cell + 1
  }
  (names, index)
}

///|
fn scrapper_log_lower_threshold(
  values : Array[Double],
  nmads : Double,
) -> Double {
  let transformed : Array[Double] = []
  let mut index = 0
  while index < values.length() {
    transformed.push(@math.ln(values[index] + 1.0))
    index = index + 1
  }
  let center = scrapper_median(transformed)
  let spread = scrapper_mad(transformed, center)
  let result = @math.exp(center - nmads * spread) - 1.0
  if result > 0.0 {
    result
  } else {
    0.0
  }
}

///|
fn scrapper_upper_threshold(values : Array[Double], nmads : Double) -> Double {
  let center = scrapper_median(values)
  center + nmads * scrapper_mad(values, center)
}

///|
fn scrapper_threshold_tolerance(threshold : Double) -> Double {
  threshold.abs().max(1.0) * 1.0e-12
}

///|
/// Suggest batch-aware RNA QC thresholds using median absolute deviations.
pub fn scrapper_suggest_rna_qc_thresholds(
  metrics : ScrapperRnaQcMetrics,
  blocks? : Array[String] = [],
  sum_num_mads? : Double = 3.0,
  detected_num_mads? : Double = 3.0,
  subset_proportion_num_mads? : Double = 3.0,
) -> ScrapperRnaQcThresholds raise ScrapperError {
  let cells = metrics.sums.length()
  if metrics.detected.length() != cells {
    raise ScrapperError("scrapper RNA QC metrics have inconsistent lengths")
  }
  if sum_num_mads < 0.0 ||
    detected_num_mads < 0.0 ||
    subset_proportion_num_mads < 0.0 {
    raise ScrapperError("scrapper MAD multipliers must be non-negative")
  }
  let mut subset = 0
  while subset < metrics.subset_proportions.length() {
    if metrics.subset_proportions[subset].length() != cells {
      raise ScrapperError(
        "scrapper RNA QC subset metrics have inconsistent lengths",
      )
    }
    subset = subset + 1
  }
  if metrics.subset_names.length() != metrics.subset_proportions.length() {
    raise ScrapperError("scrapper RNA QC subset names and values differ")
  }
  let (block_names, block_index) = scrapper_block_layout(cells, blocks)
  let n_blocks = block_names.length()
  let sum_lower = Array::make(n_blocks, 0.0)
  let detected_lower = Array::make(n_blocks, 0.0)
  let subset_upper : Array[Array[Double]] = []
  subset = 0
  while subset < metrics.subset_proportions.length() {
    subset_upper.push(Array::make(n_blocks, 0.0))
    subset = subset + 1
  }
  let mut block = 0
  while block < n_blocks {
    let sums : Array[Double] = []
    let detected : Array[Double] = []
    let subset_values : Array[Array[Double]] = []
    subset = 0
    while subset < metrics.subset_proportions.length() {
      subset_values.push([])
      subset = subset + 1
    }
    let mut cell = 0
    while cell < cells {
      if block_index[cell] == block {
        sums.push(metrics.sums[cell])
        detected.push(metrics.detected[cell].to_double())
        subset = 0
        while subset < metrics.subset_proportions.length() {
          subset_values[subset].push(metrics.subset_proportions[subset][cell])
          subset = subset + 1
        }
      }
      cell = cell + 1
    }
    sum_lower[block] = scrapper_log_lower_threshold(sums, sum_num_mads)
    detected_lower[block] = scrapper_log_lower_threshold(
      detected, detected_num_mads,
    )
    subset = 0
    while subset < metrics.subset_proportions.length() {
      subset_upper[subset][block] = scrapper_upper_threshold(
        subset_values[subset],
        subset_proportion_num_mads,
      )
      subset = subset + 1
    }
    block = block + 1
  }
  ScrapperRnaQcThresholds::{
    block_names,
    sum_lower,
    detected_lower,
    subset_upper,
  }
}

///|
/// Apply RNA QC thresholds. `true` denotes a retained high-quality cell.
pub fn scrapper_filter_rna_qc_metrics(
  thresholds : ScrapperRnaQcThresholds,
  metrics : ScrapperRnaQcMetrics,
  blocks? : Array[String] = [],
) -> Array[Bool] raise ScrapperError {
  let cells = metrics.sums.length()
  let (block_names, block_index) = scrapper_block_layout(cells, blocks)
  if block_names != thresholds.block_names {
    raise ScrapperError("scrapper QC block levels do not match the thresholds")
  }
  if thresholds.sum_lower.length() != block_names.length() ||
    thresholds.detected_lower.length() != block_names.length() {
    raise ScrapperError("scrapper QC threshold lengths are inconsistent")
  }
  if thresholds.subset_upper.length() != metrics.subset_proportions.length() {
    raise ScrapperError("scrapper QC subset threshold count is inconsistent")
  }
  let keep = Array::make(cells, true)
  let mut cell = 0
  while cell < cells {
    let block = block_index[cell]
    let sum_lower = thresholds.sum_lower[block]
    let detected_lower = thresholds.detected_lower[block]
    if metrics.sums[cell] + scrapper_threshold_tolerance(sum_lower) < sum_lower ||
      metrics.detected[cell].to_double() +
      scrapper_threshold_tolerance(detected_lower) <
      detected_lower {
      keep[cell] = false
    }
    let mut subset = 0
    while subset < metrics.subset_proportions.length() {
      if thresholds.subset_upper[subset].length() != block_names.length() {
        raise ScrapperError(
          "scrapper QC subset threshold length is inconsistent",
        )
      }
      let subset_upper = thresholds.subset_upper[subset][block]
      if metrics.subset_proportions[subset][cell] >
        subset_upper + scrapper_threshold_tolerance(subset_upper) {
        keep[cell] = false
      }
      subset = subset + 1
    }
    cell = cell + 1
  }
  keep
}

///|
/// Run RNA QC metric calculation, threshold suggestion and filtering.
pub fn scrapper_quick_rna_qc(
  counts : Array[Array[Double]],
  subsets? : Array[ScrapperNamedSubset] = [],
  blocks? : Array[String] = [],
  sum_num_mads? : Double = 3.0,
  detected_num_mads? : Double = 3.0,
  subset_proportion_num_mads? : Double = 3.0,
) -> ScrapperRnaQcResult raise ScrapperError {
  let metrics = scrapper_compute_rna_qc_metrics(counts, subsets~)
  let thresholds = scrapper_suggest_rna_qc_thresholds(
    metrics,
    blocks~,
    sum_num_mads~,
    detected_num_mads~,
    subset_proportion_num_mads~,
  )
  let keep = scrapper_filter_rna_qc_metrics(thresholds, metrics, blocks~)
  ScrapperRnaQcResult::{ metrics, thresholds, keep }
}

///|
/// Replace invalid size factors with a positive fallback.
pub fn scrapper_sanitize_size_factors(
  size_factors : Array[Double],
  fallback? : Double = 1.0,
) -> Array[Double] raise ScrapperError {
  if !scrapper_is_finite(fallback) || fallback <= 0.0 {
    raise ScrapperError("scrapper size-factor fallback must be positive")
  }
  let output = size_factors.copy()
  let mut index = 0
  while index < output.length() {
    if !scrapper_is_finite(output[index]) || output[index] <= 0.0 {
      output[index] = fallback
    }
    index = index + 1
  }
  output
}

///|
/// Center positive size factors at unity, optionally within or across blocks.
pub fn scrapper_center_size_factors(
  size_factors : Array[Double],
  blocks? : Array[String] = [],
  mode? : ScrapperCenterMode = Lowest,
) -> Array[Double] raise ScrapperError {
  if size_factors.length() == 0 {
    return []
  }
  let mut index = 0
  while index < size_factors.length() {
    if !scrapper_is_finite(size_factors[index]) || size_factors[index] <= 0.0 {
      raise ScrapperError(
        "scrapper size factors must be finite and positive before centering",
      )
    }
    index = index + 1
  }
  let (block_names, block_index) = scrapper_block_layout(
    size_factors.length(),
    blocks,
  )
  let totals = Array::make(block_names.length(), 0.0)
  let counts = Array::make(block_names.length(), 0)
  index = 0
  while index < size_factors.length() {
    let block = block_index[index]
    totals[block] = totals[block] + size_factors[index]
    counts[block] = counts[block] + 1
    index = index + 1
  }
  let means = Array::make(block_names.length(), 0.0)
  let mut block = 0
  while block < block_names.length() {
    means[block] = totals[block] / counts[block].to_double()
    block = block + 1
  }
  let output = size_factors.copy()
  match mode {
    Lowest => {
      let mut lowest = means[0]
      block = 1
      while block < means.length() {
        if means[block] < lowest {
          lowest = means[block]
        }
        block = block + 1
      }
      index = 0
      while index < output.length() {
        output[index] = output[index] / lowest
        index = index + 1
      }
    }
    PerBlock => {
      index = 0
      while index < output.length() {
        output[index] = output[index] / means[block_index[index]]
        index = index + 1
      }
    }
  }
  output
}

///|
/// Compute library-size factors and center them at unity.
pub fn scrapper_library_size_factors(
  counts : Array[Array[Double]],
  blocks? : Array[String] = [],
  mode? : ScrapperCenterMode = Lowest,
) -> Array[Double] raise ScrapperError {
  let (features, cells) = scrapper_matrix_shape(counts, true)
  if cells == 0 {
    return []
  }
  let totals = Array::make(cells, 0.0)
  let mut feature = 0
  while feature < features {
    let mut cell = 0
    while cell < cells {
      totals[cell] = totals[cell] + counts[feature][cell]
      cell = cell + 1
    }
    feature = feature + 1
  }
  let sanitized = scrapper_sanitize_size_factors(totals)
  scrapper_center_size_factors(sanitized, blocks~, mode~)
}

///|
/// Apply scaling normalization and an optional logarithmic transformation.
pub fn scrapper_normalize_counts(
  counts : Array[Array[Double]],
  size_factors : Array[Double],
  log_transform? : Bool = true,
  pseudo_count? : Double = 1.0,
  log_base? : Double = 2.0,
) -> Array[Array[Double]] raise ScrapperError {
  let (features, cells) = scrapper_matrix_shape(counts, true)
  if size_factors.length() != cells {
    raise ScrapperError(
      "scrapper size-factor length must equal the number of cells",
    )
  }
  if log_transform && (!scrapper_is_finite(pseudo_count) || pseudo_count <= 0.0) {
    raise ScrapperError("scrapper pseudo-count must be finite and positive")
  }
  if log_transform &&
    (!scrapper_is_finite(log_base) || log_base <= 0.0 || log_base == 1.0) {
    raise ScrapperError("scrapper logarithm base must be positive and not one")
  }
  let mut cell = 0
  while cell < cells {
    if !scrapper_is_finite(size_factors[cell]) || size_factors[cell] <= 0.0 {
      raise ScrapperError("scrapper size factors must be finite and positive")
    }
    cell = cell + 1
  }
  let denominator = if log_transform { @math.ln(log_base) } else { 1.0 }
  let output : Array[Array[Double]] = []
  let mut feature = 0
  while feature < features {
    let row = Array::make(cells, 0.0)
    cell = 0
    while cell < cells {
      let normalized = counts[feature][cell] / size_factors[cell]
      row[cell] = if log_transform {
        @math.ln(normalized + pseudo_count) / denominator
      } else {
        normalized
      }
      cell = cell + 1
    }
    output.push(row)
    feature = feature + 1
  }
  output
}

///|
fn scrapper_lowess_at(
  x : Double,
  x_values : Array[Double],
  y_values : Array[Double],
  window : Int,
) -> Double {
  let distances : Array[(Double, Int)] = []
  let mut index = 0
  while index < x_values.length() {
    distances.push(((x_values[index] - x).abs(), index))
    index = index + 1
  }
  distances.sort_by(fn(left : (Double, Int), right : (Double, Int)) -> Int {
    if left.0 < right.0 {
      -1
    } else if left.0 > right.0 {
      1
    } else if left.1 < right.1 {
      -1
    } else if left.1 > right.1 {
      1
    } else {
      0
    }
  })
  let usable = if window < distances.length() {
    window
  } else {
    distances.length()
  }
  let radius = if usable > 0 { distances[usable - 1].0 } else { 0.0 }
  let mut weight_sum = 0.0
  let mut weighted_x = 0.0
  let mut weighted_y = 0.0
  index = 0
  while index < usable {
    let point = distances[index].1
    let ratio = if radius > 0.0 { distances[index].0 / radius } else { 0.0 }
    let one_minus_cube = 1.0 - ratio * ratio * ratio
    let weight = one_minus_cube * one_minus_cube * one_minus_cube
    weight_sum = weight_sum + weight
    weighted_x = weighted_x + weight * x_values[point]
    weighted_y = weighted_y + weight * y_values[point]
    index = index + 1
  }
  if weight_sum <= 0.0 {
    return y_values[distances[0].1]
  }
  let mean_x = weighted_x / weight_sum
  let mean_y = weighted_y / weight_sum
  let mut covariance = 0.0
  let mut variance_x = 0.0
  index = 0
  while index < usable {
    let point = distances[index].1
    let ratio = if radius > 0.0 { distances[index].0 / radius } else { 0.0 }
    let one_minus_cube = 1.0 - ratio * ratio * ratio
    let weight = one_minus_cube * one_minus_cube * one_minus_cube
    let delta_x = x_values[point] - mean_x
    covariance = covariance + weight * delta_x * (y_values[point] - mean_y)
    variance_x = variance_x + weight * delta_x * delta_x
    index = index + 1
  }
  if variance_x <= 1.0e-20 {
    mean_y
  } else {
    mean_y + covariance / variance_x * (x - mean_x)
  }
}

///|
/// Fit a LOWESS mean-variance trend with optional quarter-root stabilization.
pub fn scrapper_fit_variance_trend(
  means : Array[Double],
  variances : Array[Double],
  mean_filter? : Bool = true,
  min_mean? : Double = 0.1,
  transform? : Bool = true,
  span? : Double = 0.3,
  min_window_count? : Int = 5,
) -> ScrapperVarianceTrend raise ScrapperError {
  if means.length() != variances.length() {
    raise ScrapperError(
      "scrapper trend means and variances must have equal length",
    )
  }
  if span <= 0.0 || span > 1.0 {
    raise ScrapperError("scrapper LOWESS span must be in (0, 1]")
  }
  if min_window_count <= 0 {
    raise ScrapperError("scrapper LOWESS minimum window count must be positive")
  }
  let x_values : Array[Double] = []
  let y_values : Array[Double] = []
  let mut index = 0
  while index < means.length() {
    if !scrapper_is_finite(means[index]) ||
      !scrapper_is_finite(variances[index]) ||
      variances[index] < 0.0 {
      raise ScrapperError(
        "scrapper trend inputs must be finite with non-negative variances",
      )
    }
    if !mean_filter || means[index] >= min_mean {
      x_values.push(means[index])
      let variance = variances[index]
      y_values.push(if transform { variance.sqrt().sqrt() } else { variance })
    }
    index = index + 1
  }
  if x_values.length() == 0 {
    raise ScrapperError("scrapper trend has no genes above the mean filter")
  }
  let raw_window = (span * x_values.length().to_double()).ceil().to_int()
  let mut window = if raw_window > min_window_count {
    raw_window
  } else {
    min_window_count
  }
  if window > x_values.length() {
    window = x_values.length()
  }
  let mut left_x = x_values[0]
  index = 1
  while index < x_values.length() {
    if x_values[index] < left_x {
      left_x = x_values[index]
    }
    index = index + 1
  }
  let left_y = scrapper_lowess_at(left_x, x_values, y_values, window)
  let fitted = Array::make(means.length(), 0.0)
  let residuals = Array::make(means.length(), 0.0)
  index = 0
  while index < means.length() {
    let transformed_fit = if mean_filter && means[index] < left_x {
      if left_x > 0.0 {
        left_y * means[index].max(0.0) / left_x
      } else {
        left_y
      }
    } else {
      scrapper_lowess_at(means[index], x_values, y_values, window)
    }
    let non_negative = transformed_fit.max(0.0)
    fitted[index] = if transform {
      non_negative * non_negative * non_negative * non_negative
    } else {
      non_negative
    }
    residuals[index] = variances[index] - fitted[index]
    index = index + 1
  }
  ScrapperVarianceTrend::{ fitted, residuals }
}

///|
/// Compute per-gene means/variances and fit the mean-variance trend.
pub fn scrapper_model_gene_variances(
  expression : Array[Array[Double]],
  mean_filter? : Bool = true,
  min_mean? : Double = 0.1,
  transform? : Bool = true,
  span? : Double = 0.3,
  min_window_count? : Int = 5,
) -> ScrapperVarianceModel raise ScrapperError {
  let (features, cells) = scrapper_matrix_shape(expression, false)
  if features == 0 || cells == 0 {
    raise ScrapperError(
      "scrapper variance modelling requires a non-empty matrix",
    )
  }
  let means = Array::make(features, 0.0)
  let variances = Array::make(features, 0.0)
  let mut feature = 0
  while feature < features {
    means[feature] = scrapper_mean(expression[feature])
    if cells > 1 {
      let mut sum_squares = 0.0
      let mut cell = 0
      while cell < cells {
        let difference = expression[feature][cell] - means[feature]
        sum_squares = sum_squares + difference * difference
        cell = cell + 1
      }
      variances[feature] = sum_squares / (cells - 1).to_double()
    }
    feature = feature + 1
  }
  let trend = scrapper_fit_variance_trend(
    means,
    variances,
    mean_filter~,
    min_mean~,
    transform~,
    span~,
    min_window_count~,
  )
  ScrapperVarianceModel::{
    means,
    variances,
    fitted: trend.fitted,
    residuals: trend.residuals,
  }
}

///|
/// Select top genes from a variance-related statistic.
pub fn scrapper_choose_highly_variable_genes(
  statistics : Array[Double],
  top? : Int = 4000,
  larger? : Bool = true,
  keep_ties? : Bool = true,
  use_bound? : Bool = true,
  bound? : Double = 0.0,
) -> Array[Int] raise ScrapperError {
  if top < 0 {
    raise ScrapperError("scrapper HVG top count must be non-negative")
  }
  if !scrapper_is_finite(bound) {
    raise ScrapperError("scrapper HVG bound must be finite")
  }
  let candidates : Array[(Int, Double)] = []
  let mut index = 0
  while index < statistics.length() {
    let value = statistics[index]
    if !scrapper_is_finite(value) {
      raise ScrapperError("scrapper HVG statistics must be finite")
    }
    let accepted = !use_bound ||
      (if larger { value >= bound } else { value <= bound })
    if accepted {
      candidates.push((index, value))
    }
    index = index + 1
  }
  candidates.sort_by(fn(left : (Int, Double), right : (Int, Double)) -> Int {
    if left.1 == right.1 {
      if left.0 < right.0 {
        -1
      } else if left.0 > right.0 {
        1
      } else {
        0
      }
    } else if larger {
      if left.1 > right.1 {
        -1
      } else {
        1
      }
    } else if left.1 < right.1 {
      -1
    } else {
      1
    }
  })
  if top == 0 || candidates.length() == 0 {
    return []
  }
  let target = if top < candidates.length() { top } else { candidates.length() }
  let output : Array[Int] = []
  index = 0
  while index < target {
    output.push(candidates[index].0)
    index = index + 1
  }
  if keep_ties && target < candidates.length() {
    let cutoff = candidates[target - 1].1
    index = target
    while index < candidates.length() && candidates[index].1 == cutoff {
      output.push(candidates[index].0)
      index = index + 1
    }
  }
  output
}

///|
pub fn ScrapperVarianceModel::highly_variable_genes(
  self : ScrapperVarianceModel,
  top? : Int = 4000,
  keep_ties? : Bool = true,
  positive_only? : Bool = true,
) -> Array[Int] raise ScrapperError {
  scrapper_choose_highly_variable_genes(
    self.residuals,
    top~,
    keep_ties~,
    use_bound=positive_only,
    bound=0.0,
  )
}

///|
fn scrapper_combination_key(levels : Array[String]) -> String {
  let mut key = ""
  let mut index = 0
  while index < levels.length() {
    key = key + levels[index].length().to_string() + ":" + levels[index] + ";"
    index = index + 1
  }
  key
}

///|
fn scrapper_combination_name(
  factor_names : Array[String],
  levels : Array[String],
) -> String {
  let mut name = ""
  let mut index = 0
  while index < levels.length() {
    if index > 0 {
      name = name + ","
    }
    name = name + factor_names[index] + "=" + levels[index]
    index = index + 1
  }
  name
}

///|
/// Aggregate expression across one or more cell grouping factors.
pub fn scrapper_aggregate_across_cells(
  expression : Array[Array[Double]],
  factors : Array[ScrapperFactor],
  compute_sum? : Bool = true,
  compute_detected? : Bool = true,
  compute_median? : Bool = false,
  detection_limit? : Double = 0.0,
) -> ScrapperAggregateResult raise ScrapperError {
  let (features, cells) = scrapper_matrix_shape(expression, false)
  if factors.length() == 0 {
    raise ScrapperError("scrapper aggregation requires at least one factor")
  }
  if !scrapper_is_finite(detection_limit) {
    raise ScrapperError("scrapper aggregation detection limit must be finite")
  }
  let factor_names : Array[String] = []
  let seen_factor_names : Map[String, Bool] = Map([], capacity=factors.length())
  let mut factor = 0
  while factor < factors.length() {
    if factors[factor].name.length() == 0 {
      raise ScrapperError("scrapper factor names must not be empty")
    }
    if seen_factor_names.contains(factors[factor].name) {
      raise ScrapperError("scrapper factor names must be unique")
    }
    if factors[factor].levels.length() != cells {
      raise ScrapperError(
        "scrapper factor length must equal the number of cells",
      )
    }
    seen_factor_names[factors[factor].name] = true
    factor_names.push(factors[factor].name)
    factor = factor + 1
  }
  let combinations : Array[Array[String]] = []
  let group_names : Array[String] = []
  let mapping : Map[String, Int] = Map([], capacity=cells)
  let index = Array::make(cells, 0)
  let mut cell = 0
  while cell < cells {
    let levels : Array[String] = []
    factor = 0
    while factor < factors.length() {
      if factors[factor].levels[cell].length() == 0 {
        raise ScrapperError("scrapper factor levels must not be empty")
      }
      levels.push(factors[factor].levels[cell])
      factor = factor + 1
    }
    let key = scrapper_combination_key(levels)
    match mapping.get(key) {
      Some(found) => index[cell] = found
      None => {
        let found = combinations.length()
        combinations.push(levels)
        group_names.push(scrapper_combination_name(factor_names, levels))
        mapping[key] = found
        index[cell] = found
      }
    }
    cell = cell + 1
  }
  let groups = combinations.length()
  let counts = Array::make(groups, 0)
  cell = 0
  while cell < cells {
    counts[index[cell]] = counts[index[cell]] + 1
    cell = cell + 1
  }
  let sums : Array[Array[Double]] = []
  let detected : Array[Array[Int]] = []
  let medians : Array[Array[Double]] = []
  let mut feature = 0
  while feature < features {
    if compute_sum {
      sums.push(Array::make(groups, 0.0))
    }
    if compute_detected {
      detected.push(Array::make(groups, 0))
    }
    if compute_median {
      medians.push(Array::make(groups, 0.0))
    }
    feature = feature + 1
  }
  feature = 0
  while feature < features {
    cell = 0
    while cell < cells {
      let group = index[cell]
      let value = expression[feature][cell]
      if compute_sum {
        sums[feature][group] = sums[feature][group] + value
      }
      if compute_detected && value > detection_limit {
        detected[feature][group] = detected[feature][group] + 1
      }
      cell = cell + 1
    }
    if compute_median {
      let group = 0
      let mut current_group = group
      while current_group < groups {
        let values : Array[Double] = []
        cell = 0
        while cell < cells {
          if index[cell] == current_group {
            values.push(expression[feature][cell])
          }
          cell = cell + 1
        }
        medians[feature][current_group] = scrapper_median(values)
        current_group = current_group + 1
      }
    }
    feature = feature + 1
  }
  ScrapperAggregateResult::{
    factor_names,
    combinations,
    group_names,
    counts,
    index,
    sums,
    detected,
    medians,
  }
}

///|
pub fn ScrapperAggregateResult::means(
  self : ScrapperAggregateResult,
) -> Array[Array[Double]] raise ScrapperError {
  if self.sums.length() == 0 {
    raise ScrapperError("scrapper aggregate sums were not computed")
  }
  let output : Array[Array[Double]] = []
  let mut feature = 0
  while feature < self.sums.length() {
    let row = Array::make(self.counts.length(), 0.0)
    let mut group = 0
    while group < self.counts.length() {
      if self.counts[group] > 0 {
        row[group] = self.sums[feature][group] / self.counts[group].to_double()
      }
      group = group + 1
    }
    output.push(row)
    feature = feature + 1
  }
  output
}

///|
pub fn ScrapperAggregateResult::n_groups(self : ScrapperAggregateResult) -> Int {
  self.combinations.length()
}

///|
pub fn ScrapperAggregateResult::summary(
  self : ScrapperAggregateResult,
) -> String {
  "ScrapperAggregate(groups=" +
  self.n_groups().to_string() +
  ", cells=" +
  self.index.length().to_string() +
  ", features=" +
  self.sums
  .length()
  .max(self.detected.length())
  .max(self.medians.length())
  .to_string() +
  ")"
}

///|
fn scrapper_copy_string_map(
  source : Map[String, Array[String]],
) -> Map[String, Array[String]] {
  let keys = source.keys().collect()
  let output : Map[String, Array[String]] = Map([], capacity=keys.length())
  let mut index = 0
  while index < keys.length() {
    output[keys[index]] = source[keys[index]].copy()
    index = index + 1
  }
  output
}

///|
fn scrapper_copy_matrix_map(
  source : Map[String, Array[Array[Double]]],
) -> Map[String, Array[Array[Double]]] {
  let keys = source.keys().collect()
  let output : Map[String, Array[Array[Double]]] = Map(
    [],
    capacity=keys.length() + 1,
  )
  let mut index = 0
  while index < keys.length() {
    output[keys[index]] = scrapper_copy_matrix(source[keys[index]])
    index = index + 1
  }
  output
}

///|
fn scrapper_copy_metadata(source : Map[String, String]) -> Map[String, String] {
  let keys = source.keys().collect()
  let output : Map[String, String] = Map([], capacity=keys.length())
  let mut index = 0
  while index < keys.length() {
    output[keys[index]] = source[keys[index]]
    index = index + 1
  }
  output
}

///|
fn scrapper_copy_alt_experiments(
  source : Map[String, SingleCellExperiment],
) -> Map[String, SingleCellExperiment] {
  let keys = source.keys().collect()
  let output : Map[String, SingleCellExperiment] = Map(
    [],
    capacity=keys.length(),
  )
  let mut index = 0
  while index < keys.length() {
    output[keys[index]] = source[keys[index]]
    index = index + 1
  }
  output
}

///|
fn scrapper_clone_sce(
  sce : SingleCellExperiment,
  assays : Map[String, Array[Array[Double]]],
  col_data : Map[String, Array[String]],
) -> SingleCellExperiment {
  SingleCellExperiment::{
    assays,
    row_data: scrapper_copy_string_map(sce.row_data),
    col_data,
    reduced_dims: scrapper_copy_matrix_map(sce.reduced_dims),
    alternative_experiments: scrapper_copy_alt_experiments(
      sce.alternative_experiments,
    ),
    metadata: scrapper_copy_metadata(sce.metadata),
    row_names: sce.row_names.copy(),
    col_names: sce.col_names.copy(),
  }
}

///|
/// Normalize an SCE's count assay into a new assay without mutating the input.
pub fn scrapper_normalize_rna_counts_sce(
  sce : SingleCellExperiment,
  size_factors? : Array[Double] = [],
  assay_name? : String = "counts",
  output_assay? : String = "logcounts",
  blocks? : Array[String] = [],
  mode? : ScrapperCenterMode = Lowest,
  pseudo_count? : Double = 1.0,
) -> SingleCellExperiment raise ScrapperError {
  let counts = sce_get_assay(sce, assay_name)
  if counts.length() == 0 {
    raise ScrapperError("scrapper SCE input assay is missing or empty")
  }
  let factors = if size_factors.length() == 0 {
    scrapper_library_size_factors(counts, blocks~, mode~)
  } else {
    scrapper_center_size_factors(size_factors, blocks~, mode~)
  }
  let normalized = scrapper_normalize_counts(counts, factors, pseudo_count~)
  let assays = scrapper_copy_matrix_map(sce.assays)
  assays[output_assay] = normalized
  let col_data = scrapper_copy_string_map(sce.col_data)
  let factor_strings : Array[String] = []
  let mut index = 0
  while index < factors.length() {
    factor_strings.push(factors[index].to_string())
    index = index + 1
  }
  col_data["size_factor"] = factor_strings
  scrapper_clone_sce(sce, assays, col_data)
}

///|
/// Run RNA QC on an SCE and return a copied SCE with QC columns plus results.
pub fn scrapper_quick_rna_qc_sce(
  sce : SingleCellExperiment,
  subsets? : Array[ScrapperNamedSubset] = [],
  blocks? : Array[String] = [],
  assay_name? : String = "counts",
  sum_num_mads? : Double = 3.0,
  detected_num_mads? : Double = 3.0,
  subset_proportion_num_mads? : Double = 3.0,
) -> (SingleCellExperiment, ScrapperRnaQcResult) raise ScrapperError {
  let counts = sce_get_assay(sce, assay_name)
  if counts.length() == 0 {
    raise ScrapperError("scrapper SCE input assay is missing or empty")
  }
  let result = scrapper_quick_rna_qc(
    counts,
    subsets~,
    blocks~,
    sum_num_mads~,
    detected_num_mads~,
    subset_proportion_num_mads~,
  )
  let col_data = scrapper_copy_string_map(sce.col_data)
  let sums : Array[String] = []
  let detected : Array[String] = []
  let keep : Array[String] = []
  let mut cell = 0
  while cell < result.metrics.sums.length() {
    sums.push(result.metrics.sums[cell].to_string())
    detected.push(result.metrics.detected[cell].to_string())
    keep.push(result.keep[cell].to_string())
    cell = cell + 1
  }
  col_data["scrapper_sum"] = sums
  col_data["scrapper_detected"] = detected
  col_data["scrapper_keep"] = keep
  let mut subset = 0
  while subset < result.metrics.subset_names.length() {
    let values : Array[String] = []
    cell = 0
    while cell < result.metrics.sums.length() {
      values.push(result.metrics.subset_proportions[subset][cell].to_string())
      cell = cell + 1
    }
    col_data["scrapper_subset_" +
    result.metrics.subset_names[subset] +
    "_proportion"] = values
    subset = subset + 1
  }
  let assays = scrapper_copy_matrix_map(sce.assays)
  (scrapper_clone_sce(sce, assays, col_data), result)
}

///|
/// A deterministic feature-by-cell count matrix for examples and tests.
pub fn scrapper_sample_counts() -> Array[Array[Double]] {
  [
    [20.0, 18.0, 22.0, 60.0, 55.0, 65.0],
    [8.0, 7.0, 9.0, 3.0, 4.0, 2.0],
    [0.0, 1.0, 0.0, 15.0, 13.0, 17.0],
    [2.0, 2.0, 3.0, 1.0, 1.0, 2.0],
    [1.0, 1.0, 1.0, 12.0, 10.0, 11.0],
    [4.0, 5.0, 4.0, 6.0, 7.0, 6.0],
  ]
}

///|
pub fn ScrapperRnaQcResult::retained_count(self : ScrapperRnaQcResult) -> Int {
  let mut retained = 0
  let mut cell = 0
  while cell < self.keep.length() {
    if self.keep[cell] {
      retained = retained + 1
    }
    cell = cell + 1
  }
  retained
}

///|
pub fn ScrapperRnaQcResult::summary(self : ScrapperRnaQcResult) -> String {
  "ScrapperRNAQC(cells=" +
  self.keep.length().to_string() +
  ", retained=" +
  self.retained_count().to_string() +
  ", subsets=" +
  self.metrics.subset_names.length().to_string() +
  ")"
}
